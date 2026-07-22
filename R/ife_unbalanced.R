# ============================================================================
# ife_unbalanced.R — Unbalanced Panel IFE via EM Algorithm
#
# Model: Y_it = mu + alpha_i + xi_t + X_it' beta + lambda_i' F_t + u_it
#   Additive fixed effects alpha_i (unit) and xi_t (time) are controlled by
#   `force` ("none" | "unit" | "time" | "two-way"); they are estimated jointly
#   with the interactive factors (see .ife_fe_em_unb).  force = "none" is the
#   intercept-free interactive model.
#
# Theory:     estimation, inference, and bias correction follow Su, Wang and
#             Wang (2025); unbalanced extension of Bai (2009)
# Algorithm:  Bai (2009) Appendix B EM procedure (inner loop); missing-data
#             factor analysis of Bai and Ng (2021)
# Inference:  sandwich standard errors and analytical bias correction
#
# Functions (internal → public):
#   .ife_em_inner()           — inner EM loop: F and Lambda given beta (force="none")
#   .ife_fe_em_unb()          — inner EM for additive + interactive FE (force!="none")
#   .ife_nnr_unb()            — NNR soft-impute initial estimator
#   .ife_fit_unb()            — outer alternating loop: beta ↔ (F, Lambda [, alpha, xi])
#   .ife_delta_omega_unb()    — alternating LS for projected regressors
#   .ife_intermediates_unb()  — auxiliary r x r matrices Delta and Xi
#   .ife_se_unb()             — sandwich standard errors
#   .ife_bias_unb()           — analytical bias correction (terms b2-b6)
#   ife_unbalanced()          — public wrapper, returns class "ife_unb"
#   ife_select_r_unb()        — SVT factor selection
#   print.ife_unb()           — print method
# ============================================================================


# ----------------------------------------------------------------------------
# .ife_em_inner — Inner EM loop (Bai 2009 Appendix B)
#
# Given a fixed beta (already absorbed into W_mat), iterates between:
#   E-step: impute missing cells with current lambda_i' F_t
#   M-step: SVD/eigen of the completed matrix to update F and Lambda
#
# @param W_mat       T x N matrix: observed W_it = Y_it - X_it'beta for
#                    observed cells; 0 elsewhere (caller sets this)
# @param obs_mask    T x N logical matrix: TRUE = observed cell
# @param Lambda_prev N x r loading matrix from previous outer iteration
# @param F_prev      T x r factor matrix from previous outer iteration
# @param N           number of units
# @param TT          number of time periods (= max T_i)
# @param r           number of factors
# @param tol_em      convergence tolerance on max fitted-value change at
#                    observed cells (default 1e-7)
# @param max_iter_em maximum EM inner iterations (default 500L)
#
# @return list(F_hat, Lambda_hat, n_iter_em, converged_em)
# ----------------------------------------------------------------------------
.ife_em_inner <- function(W_mat, obs_mask, Lambda_prev, F_prev,
                          N, TT, r,
                          tol_em      = 1e-7,
                          max_iter_em = 500L) {

  W_fill <- W_mat

  if (any(r > 0L)) {
    fitted_prev <- F_prev %*% t(Lambda_prev)
    W_fill[!obs_mask] <- fitted_prev[!obs_mask]
  }

  F_new      <- F_prev
  Lambda_new <- Lambda_prev

  for (h in seq_len(max_iter_em)) {

    if (TT <= N) {
      M     <- tcrossprod(W_fill) / (N * TT)
      eig   <- eigen(M, symmetric = TRUE)
      F_new <- eig$vectors[, seq_len(r), drop = FALSE] * sqrt(TT)
    } else {
      M       <- crossprod(W_fill) / (N * TT)
      eig     <- eigen(M, symmetric = TRUE)
      Lam_tmp <- eig$vectors[, seq_len(r), drop = FALSE] * sqrt(N)
      F_raw   <- W_fill %*% Lam_tmp / N
      # Renormalise via SVD to guarantee F'F/TT = I_r (same as TT <= N branch)
      sv_F    <- svd(F_raw, nu = r, nv = 0L)
      F_new   <- sv_F$u[, seq_len(r), drop = FALSE] * sqrt(TT)     # TT x r
    }

    Lambda_new <- crossprod(W_fill, F_new) / TT   # recomputed with normalised F_new

    fitted_new  <- F_new  %*% t(Lambda_new)
    fitted_prev <- F_prev %*% t(Lambda_prev)
    delta <- max(abs((fitted_new - fitted_prev)[obs_mask]))

    F_prev      <- F_new
    Lambda_prev <- Lambda_new
    W_fill[!obs_mask] <- fitted_new[!obs_mask]

    if (delta < tol_em) {
      return(list(
        F_hat        = F_new,
        Lambda_hat   = Lambda_new,
        n_iter_em    = h,
        converged_em = TRUE
      ))
    }
  }

  list(
    F_hat        = F_new,
    Lambda_hat   = Lambda_new,
    n_iter_em    = max_iter_em,
    converged_em = FALSE
  )
}


# ----------------------------------------------------------------------------
# .ife_fe_em_unb — EM for ADDITIVE + INTERACTIVE fixed effects (force support)
#
# Estimates the structure  S_it = mu + alpha_i + xi_t + lambda_i' f_t  on an
# unbalanced panel by the standard EM: impute the missing cells with the current
# structure, then estimate (mu, alpha, xi, F, Lambda) on the COMPLETED (balanced)
# matrix via one-shot two-way demeaning + SVD.  Because all demeaning happens on
# the imputed balanced matrix, it is NOT selection-biased (cf. the old grand-mean
# centering, which demeaned observed cells only).  Mirrors fect::interFE's
# fe_ad_inter_iter.  force = "none" reduces to factor-only EM.
#
# @param Wm        TT x N residual matrix (Y - X beta), 0 at unobserved cells
# @param obs_mask  TT x N logical, TRUE = observed
# @param FE_init   TT x N warm-start structure (0 on first call)
# @param force     "none" | "unit" | "time" | "two-way"
# @return list(FE, F_hat, Lambda_hat, alpha, xi, mu, n_iter_em, converged_em)
# ----------------------------------------------------------------------------
.ife_fe_em_unb <- function(Wm, obs_mask, FE_init, force, r,
                           tol_em = 1e-7, max_iter_em = 500L) {
  TT <- nrow(Wm); N <- ncol(Wm)
  FE <- FE_init
  alpha <- numeric(N); xi <- numeric(TT); mu <- 0
  F_hat <- matrix(0, TT, r); Lambda_hat <- matrix(0, N, r)
  has_unit <- force %in% c("unit", "two-way")
  has_time <- force %in% c("time", "two-way")
  converged_em <- FALSE
  h <- 0L
  for (h in seq_len(max_iter_em)) {
    full <- Wm
    full[!obs_mask] <- FE[!obs_mask]          # E-step: impute with current structure

    # ---- additive FE on the completed balanced matrix (one-shot, unbiased) ----
    if (force == "none") {
      mu <- 0; alpha <- numeric(N); xi <- numeric(TT)
    } else {
      mu    <- mean(full)
      alpha <- if (has_unit) colMeans(full) - mu else numeric(N)
      xi    <- if (has_time) rowMeans(full) - mu else numeric(TT)
    }
    add <- matrix(alpha, TT, N, byrow = TRUE) + matrix(xi, TT, N) + mu
    D   <- full - add

    # ---- interactive factors via SVD (F'F/TT = I_r) ----
    sv         <- svd(D, nu = r, nv = 0L)
    F_hat      <- sv$u[, seq_len(r), drop = FALSE] * sqrt(TT)
    Lambda_hat <- crossprod(D, F_hat) / TT

    FE_new <- add + F_hat %*% t(Lambda_hat)
    if (max(abs((FE_new - FE)[obs_mask])) < tol_em) {
      FE <- FE_new; converged_em <- TRUE; break
    }
    FE <- FE_new
  }
  list(FE = FE, F_hat = F_hat, Lambda_hat = Lambda_hat,
       alpha = alpha, xi = xi, mu = mu,
       n_iter_em = h, converged_em = converged_em)
}


# ----------------------------------------------------------------------------
# .ife_nnr_unb — NNR soft-impute initial estimator
#
# Minimises  (1/2) sum_{it} d_it (Y_it - Theta_it - X_it' beta)^2 + nu ||Theta||_*
# via soft-impute iteration.  Cross-validates the penalty nu over
# nu_cands = c * sqrt(max(N, TT)) for c in c_grid.
#
# @param Y_long      n_obs outcome vector
# @param X_long      n_obs x p covariate matrix (may have 0 columns)
# @param obs_lin     n_obs integer vector: linear index into TT x N matrix
# @param N, TT       panel dimensions
# @param r           number of factors for factor extraction
# @param c_grid      grid of scaling constants; default c(0.01, 0.1, 1, 10)
# @param tol_nnr     convergence tolerance (default 1e-6)
# @param max_iter_nnr maximum iterations per candidate (default 200L)
#
# @return list(beta0, Theta0, F_hat_0, Lambda_hat_0, nu_used)
# ----------------------------------------------------------------------------
.ife_nnr_unb <- function(Y_long, X_long, obs_lin,
                          N, TT, r,
                          c_grid       = c(0.01, 0.1, 1, 10),
                          tol_nnr      = 1e-6,
                          max_iter_nnr = 200L) {

  p        <- ncol(X_long)
  n_obs    <- length(Y_long)
  nu_cands <- c_grid * sqrt(max(N, TT))

  best_loss   <- Inf
  best_nu     <- nu_cands[1L]
  best_Theta  <- matrix(0, TT, N)
  best_beta   <- if (p > 0L) numeric(p) else numeric(0L)

  for (nu in nu_cands) {

    # ---- Initialise ----
    if (p > 0L) {
      beta <- as.vector(solve(crossprod(X_long), crossprod(X_long, Y_long)))
    } else {
      beta <- numeric(0L)
    }
    Theta <- matrix(0, TT, N)

    # ---- Soft-impute iterations ----
    for (iter in seq_len(max_iter_nnr)) {
      Theta_old <- Theta

      # Residual at observed cells; fill W with Theta at missing
      if (p > 0L) {
        W_obs <- Y_long - as.vector(X_long %*% beta)
      } else {
        W_obs <- Y_long
      }
      W_fill          <- Theta
      W_fill[obs_lin] <- W_obs

      # Singular value soft-thresholding
      sv        <- svd(W_fill)
      d_shrunk  <- pmax(sv$d - nu, 0)
      Theta     <- sweep(sv$u, 2L, d_shrunk, `*`) %*% t(sv$v)  # TT x N

      # OLS update for beta given Theta
      if (p > 0L) {
        Theta_obs <- Theta[obs_lin]
        beta      <- as.vector(
          solve(crossprod(X_long), crossprod(X_long, Y_long - Theta_obs))
        )
      }

      # Convergence: max change at observed cells
      if (max(abs((Theta - Theta_old)[obs_lin])) < tol_nnr) break
    }

    # ---- Evaluate loss: rank-r approximation of Theta at observed cells ----
    r_eff    <- min(r, length(sv$d))
    d_r      <- sv$d[seq_len(r_eff)]
    Theta_r  <- sv$u[, seq_len(r_eff), drop = FALSE] %*%
                (d_r * t(sv$v[, seq_len(r_eff), drop = FALSE]))

    if (p > 0L) {
      resid <- Y_long - as.vector(X_long %*% beta) - Theta_r[obs_lin]
    } else {
      resid <- Y_long - Theta_r[obs_lin]
    }
    loss <- sum(resid^2) / n_obs

    if (loss < best_loss) {
      best_loss  <- loss
      best_nu    <- nu
      best_Theta <- Theta
      best_beta  <- beta
    }
  }

  # ---- Extract factors from best Theta using EM normalisation ----
  # F_hat_0 = sqrt(TT) * first r left singular vectors
  # Lambda_hat_0 = t(Theta0) %*% F_hat_0 / TT  (same as Lambda update in EM)
  sv0  <- svd(best_Theta, nu = r, nv = 0L)
  r_eff <- min(r, length(sv0$d))
  F_hat_0      <- sv0$u[, seq_len(r_eff), drop = FALSE] * sqrt(TT)   # TT x r
  Lambda_hat_0 <- crossprod(best_Theta, F_hat_0) / TT                 # N x r

  list(
    beta0        = best_beta,
    Theta0       = best_Theta,
    F_hat_0      = F_hat_0,
    Lambda_hat_0 = Lambda_hat_0,
    nu_used      = best_nu
  )
}


# ----------------------------------------------------------------------------
# .ife_fit_unb — Outer alternating loop for unbalanced IFE
#
# Alternates between:
#   (a) EM inner loop to update F and Lambda given beta
#   (b) Unit-specific FWL projection to update beta given F
#
# @param Y_long, X_long  long vectors/matrix (observed cells only)
# @param unit_idx        n_obs unit indices (1-based)
# @param time_idx        n_obs time indices (1-based)
# @param N, TT           panel dimensions
# @param r               number of factors
# @param beta_init       optional warm-start beta (NULL = OLS)
# @param F_init          optional warm-start F (NULL = zeros)
# @param Lambda_init     optional warm-start Lambda (NULL = zeros)
# @param tol, max_iter   outer-loop convergence
# @param tol_em, max_iter_em  inner EM convergence
#
# @return list: beta, F_hat, Lambda_hat, X_tilde_long, u_long,
#               obs_mask, obs_lin, unit_rows, unit_obs_t, n_iter, converged
# ----------------------------------------------------------------------------
.ife_fit_unb <- function(Y_long, X_long, unit_idx, time_idx,
                         N, TT, r,
                         force        = "none",
                         beta_init    = NULL,
                         F_init       = NULL,
                         Lambda_init  = NULL,
                         tol         = 1e-9,
                         max_iter    = 10000L,
                         tol_em      = 1e-7,
                         max_iter_em = 500L) {

  n_obs <- length(Y_long)
  p     <- ncol(X_long)

  obs_lin  <- (unit_idx - 1L) * TT + time_idx
  obs_mask <- matrix(FALSE, TT, N)
  obs_mask[obs_lin] <- TRUE

  unit_rows  <- lapply(seq_len(N), function(i) which(unit_idx == i))
  unit_obs_t <- lapply(unit_rows,  function(idx) time_idx[idx])

  # r = 0 special case
  if (r == 0L) {
    if (p > 0L) {
      beta <- as.vector(solve(crossprod(X_long), crossprod(X_long, Y_long)))
    } else {
      beta <- numeric(0)
    }
    F_hat        <- matrix(0, TT, 0)
    Lambda_hat   <- matrix(0, N,  0)
    X_tilde_long <- X_long
    u_long <- Y_long
    if (p > 0L) u_long <- u_long - as.vector(X_long %*% beta)
    return(list(
      beta         = beta,
      F_hat        = F_hat,
      Lambda_hat   = Lambda_hat,
      X_tilde_long = X_tilde_long,
      u_long       = u_long,
      obs_mask     = obs_mask,
      obs_lin      = obs_lin,
      unit_rows    = unit_rows,
      unit_obs_t   = unit_obs_t,
      n_iter       = 0L,
      converged    = TRUE
    ))
  }

  # ---- Initialisation (with optional warm start) ----
  if (!is.null(beta_init)) {
    beta <- beta_init
  } else if (p > 0L) {
    beta <- as.vector(solve(crossprod(X_long), crossprod(X_long, Y_long)))
  } else {
    beta <- numeric(0)
  }

  F_hat      <- if (!is.null(F_init))      F_init      else matrix(0, TT, r)
  Lambda_hat <- if (!is.null(Lambda_init)) Lambda_init else matrix(0, N,  r)

  # ==========================================================================
  # force != "none": joint ADDITIVE + INTERACTIVE estimation
  # (OLS-of-(Y - structure)-on-X with the structure estimated by EM on the
  #  imputed balanced residual; mirrors fect::interFE).  Needed because
  #  additive and interactive FE are NOT orthogonal, so they must be estimated
  #  jointly; the imputed-balanced demeaning is selection-robust.
  # ==========================================================================
  if (force != "none" && p > 0L) {
    Ym <- matrix(0, TT, N); Ym[obs_lin] <- Y_long
    Xc <- lapply(seq_len(p), function(k) {
      m <- matrix(0, TT, N); m[obs_lin] <- X_long[, k]; m
    })
    xxinv <- solve(crossprod(X_long))      # (X'X)^{-1} over observed cells
    FE    <- matrix(0, TT, N)              # warm-start full structure
    alpha <- numeric(N); xi <- numeric(TT); mu <- 0
    converged <- FALSE; n_iter <- 0L
    for (iter in seq_len(max_iter)) {
      n_iter   <- iter
      beta_old <- beta
      covar <- matrix(0, TT, N)
      for (k in seq_len(p)) covar <- covar + Xc[[k]] * beta[k]
      W_mat <- Ym - covar; W_mat[!obs_mask] <- 0
      em <- .ife_fe_em_unb(W_mat, obs_mask, FE, force, r,
                           tol_em = tol_em, max_iter_em = max_iter_em)
      FE <- em$FE; F_hat <- em$F_hat; Lambda_hat <- em$Lambda_hat
      alpha <- em$alpha; xi <- em$xi; mu <- em$mu
      beta_new <- as.vector(xxinv %*% crossprod(X_long, (Ym - FE)[obs_lin]))
      if (max(abs(beta_new - beta_old)) < tol) {
        beta <- beta_new; converged <- TRUE; break
      }
      beta <- beta_new
    }
    if (!converged)
      warning("ife_unbalanced (force = '", force, "'): outer loop did not ",
              "converge in ", max_iter, " iterations. Estimate may be ",
              "unreliable; consider force = 'none' or a larger panel.")
    covar <- matrix(0, TT, N)
    for (k in seq_len(p)) covar <- covar + Xc[[k]] * beta[k]
    u_long <- (Ym - covar - FE)[obs_lin]
    return(list(
      beta = beta, F_hat = F_hat, Lambda_hat = Lambda_hat,
      alpha = alpha, xi = xi, mu = mu, force = force,
      X_tilde_long = X_long, u_long = u_long,
      obs_mask = obs_mask, obs_lin = obs_lin,
      unit_rows = unit_rows, unit_obs_t = unit_obs_t,
      n_iter = n_iter, converged = converged
    ))
  }

  X_tilde_long <- X_long
  converged    <- FALSE
  n_iter       <- 0L

  for (iter in seq_len(max_iter)) {
    n_iter   <- iter
    beta_old <- beta

    if (p > 0L) {
      W_long <- Y_long - as.vector(X_long %*% beta)
    } else {
      W_long <- Y_long
    }

    W_mat          <- matrix(0, TT, N)
    W_mat[obs_lin] <- W_long

    em <- .ife_em_inner(
      W_mat       = W_mat,
      obs_mask    = obs_mask,
      Lambda_prev = Lambda_hat,
      F_prev      = F_hat,
      N = N, TT = TT, r = r,
      tol_em = tol_em, max_iter_em = max_iter_em
    )
    F_hat      <- em$F_hat
    Lambda_hat <- em$Lambda_hat

    if (!em$converged_em)
      warning("EM inner loop did not converge at outer iteration ", iter,
              ". Consider increasing max_iter_em or tol_em.")

    if (p > 0L) {
      A_mat <- matrix(0, p, p)
      b_vec <- numeric(p)

      for (i in seq_len(N)) {
        obs_i <- unit_rows[[i]]
        t_i   <- unit_obs_t[[i]]
        T_i   <- length(t_i)
        if (T_i == 0L) next

        F_i <- F_hat[t_i, , drop = FALSE]
        X_i <- X_long[obs_i, , drop = FALSE]
        Y_i <- Y_long[obs_i]

        FiFi_inv  <- solve(crossprod(F_i))
        FtX_i     <- crossprod(F_i, X_i)
        FtY_i     <- crossprod(F_i, Y_i)
        X_tilde_i <- X_i - F_i %*% (FiFi_inv %*% FtX_i)
        Y_tilde_i <- Y_i - as.vector(F_i %*% (FiFi_inv %*% FtY_i))

        A_mat <- A_mat + crossprod(X_tilde_i)
        b_vec <- b_vec + as.vector(crossprod(X_tilde_i, Y_tilde_i))

        X_tilde_long[obs_i, ] <- X_tilde_i
      }

      beta_new <- as.vector(solve(A_mat, b_vec))

      if (max(abs(beta_new - beta_old)) < tol) {
        beta      <- beta_new
        converged <- TRUE
        break
      }
      beta <- beta_new

    } else {
      X_tilde_long <- NULL
      converged    <- TRUE
      break
    }
  }

  fitted_lf <- rowSums(Lambda_hat[unit_idx, , drop = FALSE] *
                         F_hat[time_idx, , drop = FALSE])
  u_long <- Y_long - fitted_lf
  if (p > 0L) u_long <- u_long - as.vector(X_long %*% beta)

  list(
    beta         = beta,
    F_hat        = F_hat,
    Lambda_hat   = Lambda_hat,
    alpha        = numeric(N),
    xi           = numeric(TT),
    mu           = 0,
    force        = "none",
    X_tilde_long = X_tilde_long,
    u_long       = u_long,
    obs_mask     = obs_mask,
    obs_lin      = obs_lin,
    unit_rows    = unit_rows,
    unit_obs_t   = unit_obs_t,
    n_iter       = n_iter,
    converged    = converged
  )
}


# ----------------------------------------------------------------------------
# .ife_delta_omega_unb — Alternating LS for projected regressors
#
# For each regressor k = 1..p, solves:
#   min_{delta_ki, omega_kt} sum_{it} d_it (x_itk - delta_ki' f_t - omega_kt' lambda_i)^2
#
# Step A (per unit i):
#   delta_ki = [L_lam_inv]_i * sum_{t in T_i} f_t (x_itk - lambda_i' omega_kt)
#
# Step B (per time t):
#   omega_kt = [L_ff_inv]_t  * sum_{i in I_t} lambda_i (x_itk - f_t' delta_ki)
#
# where  [L_lam_inv]_i = (sum_{t in T_i} f_t f_t')^{-1}  (r x r)
#        [L_ff_inv]_t  = (sum_{i in I_t} lambda_i lambda_i')^{-1}  (r x r)
#
# @return list: delta_arr (list of p N x r matrices),
#               omega_arr (list of p TT x r matrices),
#               L_lam_inv, L_ff_inv, time_rows, time_obs_i
# ----------------------------------------------------------------------------
.ife_delta_omega_unb <- function(X_long, F_hat, Lambda_hat,
                                  unit_idx, time_idx,
                                  N, TT, r, p,
                                  unit_rows, unit_obs_t,
                                  tol_do      = 1e-6,
                                  max_iter_do = 200L) {

  # ---- Precompute per-time counterparts of unit_rows / unit_obs_t ----
  time_rows  <- lapply(seq_len(TT), function(t) which(time_idx == t))
  time_obs_i <- lapply(time_rows,   function(idx) unit_idx[idx])

  # ---- Precompute L_lam_inv[[i]] = (F_hat[T_i,]' F_hat[T_i,])^{-1} ----
  # NA-guard (mirrors L_ff_inv below): a unit observed fewer than r times gives a
  # rank-deficient crossprod(F_i); return NA so downstream anyNA() checks skip it
  # instead of solve() throwing.  Matters when force != "none" inflates r to r_eff.
  L_lam_inv <- lapply(seq_len(N), function(i) {
    t_i <- unit_obs_t[[i]]
    F_i <- F_hat[t_i, , drop = FALSE]    # T_i x r
    if (nrow(F_i) < r) return(matrix(NA_real_, r, r))
    solve(crossprod(F_i))
  })

  # ---- Precompute L_ff_inv[[t]] = (Lambda_hat[I_t,]' Lambda_hat[I_t,])^{-1} ----
  L_ff_inv <- lapply(seq_len(TT), function(t) {
    i_t <- time_obs_i[[t]]
    if (length(i_t) == 0L) return(matrix(NA_real_, r, r))
    Lam_t <- Lambda_hat[i_t, , drop = FALSE]   # |I_t| x r
    if (nrow(Lam_t) < r) return(matrix(NA_real_, r, r))
    solve(crossprod(Lam_t))
  })

  # ---- Alternating LS for each regressor k ----
  delta_arr <- vector("list", p)
  omega_arr <- vector("list", p)

  for (k in seq_len(p)) {
    x_k     <- X_long[, k]             # n_obs values
    delta_k <- matrix(0, N,  r)        # N x r
    omega_k <- matrix(0, TT, r)        # TT x r

    for (it in seq_len(max_iter_do)) {
      delta_old <- delta_k

      # ---- Step A: update delta_ki per unit i ----
      for (i in seq_len(N)) {
        obs_i    <- unit_rows[[i]]
        t_i      <- unit_obs_t[[i]]
        if (length(t_i) == 0L || anyNA(L_lam_inv[[i]])) next

        F_i      <- F_hat[t_i, , drop = FALSE]     # T_i x r
        lam_i    <- Lambda_hat[i, ]                 # r-vector
        omega_i  <- omega_k[t_i, , drop = FALSE]   # T_i x r
        x_ki     <- x_k[obs_i]                     # T_i

        # lambda_i' omega_kt for each t in T_i: (T_i x r) %*% (r x 1)
        lam_omega <- as.vector(omega_i %*% lam_i)  # T_i
        rhs_i     <- x_ki - lam_omega

        delta_k[i, ] <- as.vector(L_lam_inv[[i]] %*% crossprod(F_i, rhs_i))
      }

      # ---- Step B: update omega_kt per time t ----
      for (t in seq_len(TT)) {
        obs_t <- time_rows[[t]]
        i_t   <- time_obs_i[[t]]
        if (length(i_t) == 0L || anyNA(L_ff_inv[[t]])) next

        Lam_t   <- Lambda_hat[i_t, , drop = FALSE]   # |I_t| x r
        f_t     <- F_hat[t, ]                         # r-vector
        delta_t <- delta_k[i_t, , drop = FALSE]       # |I_t| x r
        x_kt    <- x_k[obs_t]                         # |I_t|

        # f_t' delta_ki for each i in I_t: (|I_t| x r) %*% (r x 1)
        f_delta <- as.vector(delta_t %*% f_t)         # |I_t|
        rhs_t   <- x_kt - f_delta

        omega_k[t, ] <- as.vector(L_ff_inv[[t]] %*% crossprod(Lam_t, rhs_t))
      }

      if (max(abs(delta_k - delta_old)) < tol_do) break
    }

    delta_arr[[k]] <- delta_k
    omega_arr[[k]] <- omega_k
  }

  list(
    delta_arr  = delta_arr,
    omega_arr  = omega_arr,
    L_lam_inv  = L_lam_inv,
    L_ff_inv   = L_ff_inv,
    time_rows  = time_rows,
    time_obs_i = time_obs_i
  )
}


# ----------------------------------------------------------------------------
# .ife_intermediates_unb — Auxiliary r x r matrices Delta_ki and Xi_kt
#
# Delta_ki = [L_lam_inv]_i (sum_{t in T_i} f_t omega_kt') [L_lam_inv]_i  (r x r)
# Xi_kt    = [L_ff_inv]_t  (sum_{i in I_t} lambda_i delta_ki') [L_ff_inv]_t  (r x r)
#
# @return list: Delta_arr[[k]][[i]] and Xi_arr[[k]][[t]]
# ----------------------------------------------------------------------------
.ife_intermediates_unb <- function(delta_arr, omega_arr,
                                    F_hat, Lambda_hat,
                                    L_ff_inv, L_lam_inv,
                                    unit_obs_t, time_obs_i,
                                    N, TT, r, p) {

  # Delta_arr[[k]][[i]]: r x r
  Delta_arr <- lapply(seq_len(p), function(k) {
    omega_k <- omega_arr[[k]]   # TT x r
    lapply(seq_len(N), function(i) {
      t_i     <- unit_obs_t[[i]]
      if (length(t_i) == 0L) return(matrix(0, r, r))
      F_i     <- F_hat[t_i, , drop = FALSE]       # T_i x r
      omega_i <- omega_k[t_i, , drop = FALSE]     # T_i x r
      # sum_t f_t omega_kt' = crossprod(F_i, omega_i)  (r x r)
      mid <- crossprod(F_i, omega_i)
      L   <- L_lam_inv[[i]]
      L %*% mid %*% L
    })
  })

  # Xi_arr[[k]][[t]]: r x r
  Xi_arr <- lapply(seq_len(p), function(k) {
    delta_k <- delta_arr[[k]]   # N x r
    lapply(seq_len(TT), function(t) {
      i_t <- time_obs_i[[t]]
      if (length(i_t) == 0L || anyNA(L_ff_inv[[t]])) return(matrix(0, r, r))
      Lam_t   <- Lambda_hat[i_t, , drop = FALSE]   # |I_t| x r
      delta_t <- delta_k[i_t, , drop = FALSE]       # |I_t| x r
      # sum_i lambda_i delta_ki' = crossprod(Lam_t, delta_t)  (r x r)
      mid <- crossprod(Lam_t, delta_t)
      L   <- L_ff_inv[[t]]
      L %*% mid %*% L
    })
  })

  list(Delta_arr = Delta_arr, Xi_arr = Xi_arr)
}


# ----------------------------------------------------------------------------
# .ife_se_unb — Sandwich standard errors
#
# Uses factor-projected regressors:
#   x_hat_itk = x_itk - delta_ki' f_t - omega_kt' lambda_i
#
# Sandwich: Var(beta_hat) = W_x^{-1} Omega_x W_x^{-1} / NT
#   W_x    = (1/NT) sum_{it} d_it x_hat_it x_hat_it'   (p x p)
#   Omega_x depends on se_type:
#     "standard" / "robust" — i.i.d. or HC1 (cross-sectionally independent)
#     "cluster"             — cluster-robust by unit
#     "hac"                 — HAC with Bartlett kernel (Newey-West 1987)
#
# @param beta         p-vector of estimated coefficients
# @param X_long       n_obs x p raw covariate matrix
# @param delta_arr    list of p N x r matrices (delta_ki)
# @param omega_arr    list of p TT x r matrices (omega_kt)
# @param F_hat        TT x r factor matrix
# @param Lambda_hat   N x r loading matrix
# @param u_long       n_obs residuals
# @param unit_idx, time_idx  n_obs index vectors
# @param N, TT, r     dimensions
# @param se_type      "standard" | "robust" | "cluster" | "hac"
# @param n_obs        number of observed cells
# @param unit_rows    list of N integer vectors — row indices per unit (for hac)
# @param unit_obs_t   list of N integer vectors — time positions per unit (for hac)
# @param L_T          Bartlett bandwidth (integer; used only for se_type = "hac")
#
# @return list(vcov_mat, df, W_x_inv, x_proj_long)
# ----------------------------------------------------------------------------
.ife_se_unb <- function(beta, X_long, delta_arr, omega_arr,
                         F_hat, Lambda_hat,
                         u_long, unit_idx, time_idx,
                         N, TT, r, se_type, n_obs,
                         unit_rows  = NULL,
                         unit_obs_t = NULL,
                         L_T        = NULL,
                         fe_dof     = 0L) {

  p <- length(beta)
  if (p == 0L) stop("Cannot compute SE with no covariates.")

  # ---- Build x_proj_long (n_obs x p): factor-projected regressors ----
  # x_hat_itk = x_itk - delta_ki' f_t - omega_kt' lambda_i
  x_proj_long <- X_long

  for (k in seq_len(p)) {
    dk <- delta_arr[[k]]    # N  x r
    wk <- omega_arr[[k]]    # TT x r

    # delta_ki' f_t for each obs l: rowSums( dk[i,] * F_hat[t,] )
    d_dot_f   <- rowSums(dk[unit_idx, , drop = FALSE] *
                         F_hat[time_idx, , drop = FALSE])

    # omega_kt' lambda_i for each obs l: rowSums( wk[t,] * Lambda_hat[i,] )
    w_dot_lam <- rowSums(wk[time_idx, , drop = FALSE] *
                         Lambda_hat[unit_idx, , drop = FALSE])

    x_proj_long[, k] <- X_long[, k] - d_dot_f - w_dot_lam
  }

  # ---- Degrees of freedom ----
  # r here is the TRUE number of interactive factors; fe_dof is the additive-FE
  # dof (0 / N-1 / T-1 / N+T-2).  The projection above uses the augmented
  # F_hat/delta/omega (which carry extra additive-FE columns), but the parameter
  # count must use the true interactive r plus the additive-FE dof.
  k_total <- p + r * (N + TT - r) + fe_dof
  df <- n_obs - k_total
  if (df <= 0L)
    stop("Degrees of freedom = ", df, " <= 0. Reduce r or use a larger panel.")

  # ---- A = X_proj' X_proj (p x p),  W_x = A / (N*TT) ----
  A     <- crossprod(x_proj_long)
  A_inv <- solve(A)

  # W_x^{-1} = N*TT * A^{-1}  — returned for bias correction
  W_x_inv <- A_inv * (N * TT)

  # ---- Variance estimator ----
  if (se_type == "standard") {

    sigma2   <- sum(u_long^2) / df
    vcov_mat <- sigma2 * A_inv

  } else if (se_type == "robust") {

    # HC1: B = X_proj' diag(u^2) X_proj  (vectorised)
    B        <- crossprod(x_proj_long * u_long)
    corr     <- n_obs / (n_obs - p)
    vcov_mat <- A_inv %*% B %*% A_inv * corr

  } else if (se_type == "cluster") {

    B <- matrix(0, p, p)
    for (i in seq_len(N)) {
      obs_i <- which(unit_idx == i)
      if (length(obs_i) == 0L) next
      Xp_i  <- x_proj_long[obs_i, , drop = FALSE]
      u_i   <- u_long[obs_i]
      psi_i <- as.vector(crossprod(Xp_i, u_i))
      B     <- B + tcrossprod(psi_i)
    }
    corr     <- (N / (N - 1L)) * ((n_obs - 1L) / (n_obs - p))
    vcov_mat <- A_inv %*% B %*% A_inv * corr

  } else {   # "hac" — Bartlett kernel (Newey-West 1987)

    # B_hac accumulates UNNORMALIZED sum (same convention as B in "robust"):
    #   B_hac = sum_i sum_{t,s in obs_i}
    #             Gamma(|t_val - s_val| / L_T) * v_it * v_is * xhat_it xhat_is'
    # Gamma(u) = (1 - u) * 1{u <= 1}  (Bartlett kernel)
    # gap = 0 → kern = 1, so diagonal (t=s) block = sum v_it^2 xhat_it xhat_it'
    #         → reduces to the "robust" B when L_T = 1 (only gap=0 survives).
    #
    # Normalization: V(beta_hat) = (1/NT) W_x^{-1} Omega_x W_x^{-1}
    #   with W_x = A/(NT) and Omega_x = B_hac/(NT) gives
    #   V(beta_hat) = A_inv %*% B_hac %*% A_inv  (same formula as "robust").
    #   Do NOT pre-divide B_hac by NT here.
    B_hac <- matrix(0, p, p)
    for (i in seq_len(N)) {
      rows_i <- unit_rows[[i]]
      obs_t  <- unit_obs_t[[i]]
      T_i    <- length(obs_t)
      for (idx_t in seq_len(T_i)) {
        l_t  <- rows_i[idx_t]
        v_t  <- u_long[l_t]
        xh_t <- x_proj_long[l_t, , drop = FALSE]    # 1 x p
        for (idx_s in seq_len(T_i)) {
          gap <- abs(obs_t[idx_s] - obs_t[idx_t])
          if (gap > L_T) next
          kern <- 1 - gap / L_T
          l_s  <- rows_i[idx_s]
          v_s  <- u_long[l_s]
          xh_s <- x_proj_long[l_s, , drop = FALSE]   # 1 x p
          B_hac <- B_hac + kern * v_t * v_s * crossprod(xh_t, xh_s)
        }
      }
    }
    # No division by NT — B_hac is unnormalized, parallel to B in "robust"
    vcov_mat <- A_inv %*% B_hac %*% A_inv

  }

  list(vcov_mat    = vcov_mat,
       df          = df,
       W_x_inv     = W_x_inv,
       x_proj_long = x_proj_long)
}


# ----------------------------------------------------------------------------
# .ife_bias_unb — Analytical bias correction
#
# Computes bias terms for two cases:
#
# Case A — strictly exogenous regressors, i.i.d. errors (exog = "strict"):
#   b3[k] = (1/sqrt(NT)) sum_{it} d_it v_it^2 delta_ki' [L_ff_inv]_t lambda_i
#   b4[k] = (1/sqrt(NT)) sum_{it} d_it v_it^2 omega_kt' [L_lam_inv]_i f_t
#   b5[k] = (1/sqrt(NT)) sum_{it} d_it v_it^2 lambda_i' Xi_kt lambda_i
#   b6[k] = (1/sqrt(NT)) sum_{it} d_it v_it^2 f_t' Delta_ki f_t
#   b2[k] = 0
#
# Case B — weakly exogenous regressors (exog = "weak"), adds:
#   b2[k] = (1/sqrt(NT)) sum_i sum_{t} sum_{s > t in obs_i}
#             Gamma((s-t)/L_T) * d_it v_it * d_is x_isk * f_t' [L_lam_inv]_i f_s
#           (one-sided kernel; uses RAW regressor x_isk, not projected)
#
# Case C — serial correlation (se_type = "hac"), replaces b4/b6 with:
#   b4[k] (HAC) = (1/sqrt(NT)) sum_i sum_{t,s}
#                   Gamma(|t-s|/L_T) * v_it * v_is * omega_kt' [L_lam_inv]_i f_s
#   b6[k] (HAC) = (1/sqrt(NT)) sum_i sum_{t,s}
#                   Gamma(|t-s|/L_T) * v_it * v_is * f_t' Delta_ki f_s
#
# Bias-corrected estimator:
#   beta_abc = beta_hat - (1/sqrt(NT)) W_x^{-1} b_hat
#   b_hat = b2 + b3 + b4 + b5 + b6
#
# @param beta         p-vector (the converged estimate)
# @param delta_arr    list of p N x r matrices
# @param omega_arr    list of p TT x r matrices
# @param F_hat        TT x r
# @param Lambda_hat   N x r
# @param u_long       n_obs residuals
# @param unit_idx, time_idx  n_obs index vectors
# @param L_ff_inv     list of TT  r x r matrices
# @param L_lam_inv    list of N   r x r matrices
# @param Delta_arr    list of p lists of N  r x r matrices
# @param Xi_arr       list of p lists of TT r x r matrices
# @param W_x_inv      p x p  (= NT * solve(X_proj' X_proj))
# @param N, TT, p     dimensions
# @param n_obs        observed cells
# @param X_long       n_obs x p raw covariate matrix (needed for b2)
# @param unit_rows    list of N integer vectors (needed for b2 and HAC b4/b6)
# @param unit_obs_t   list of N integer vectors (needed for b2 and HAC b4/b6)
# @param exog         "strict" | "weak"
# @param se_type      "standard" | "robust" | "cluster" | "hac"
# @param L_T          Bartlett bandwidth (used for b2 when exog="weak" and HAC terms)
#
# @return list(beta_abc, b_hat, b2, b3, b4, b5, b6)
# ----------------------------------------------------------------------------
.ife_bias_unb <- function(beta, delta_arr, omega_arr,
                           F_hat, Lambda_hat, u_long,
                           unit_idx, time_idx,
                           L_ff_inv, L_lam_inv,
                           Delta_arr, Xi_arr,
                           W_x_inv, N, TT, p, n_obs,
                           X_long     = NULL,
                           unit_rows  = NULL,
                           unit_obs_t = NULL,
                           exog       = "strict",
                           se_type    = "standard",
                           L_T        = NULL) {

  b3 <- numeric(p)
  b4 <- numeric(p)
  b5 <- numeric(p)
  b6 <- numeric(p)

  for (l in seq_len(n_obs)) {
    i   <- unit_idx[l]
    t   <- time_idx[l]
    v2  <- u_long[l]^2

    lam_i  <- Lambda_hat[i, ]     # r-vector
    f_t    <- F_hat[t, ]          # r-vector
    Lff_t  <- L_ff_inv[[t]]       # r x r
    Llam_i <- L_lam_inv[[i]]      # r x r

    if (anyNA(Lff_t) || anyNA(Llam_i)) next

    for (k in seq_len(p)) {
      delta_ki <- delta_arr[[k]][i, ]   # r-vector
      omega_kt <- omega_arr[[k]][t, ]   # r-vector

      # The bias-correction theory defines [Lbar_ff']_t = (-sum_i d_it lam_i lam_i')^{-1}
      # and [Lbar_lam']_i = (-sum_t d_it f_t f_t')^{-1}, i.e. the NEGATIVE of the
      # positive-definite inverses stored in Lff_t / Llam_i.  b3 and b4 use this
      # inverse ONCE, so they pick up the leading minus sign.  b5 and b6 use it
      # TWICE (inside Xi_kt and Delta_ki), so the two minus signs cancel and the
      # pre-computed Xi_arr / Delta_arr (built from the positive inverses) are
      # already correct.  Omitting this minus on b3/b4 makes the correction
      # double the b5/b6 magnitude instead of cancelling it (spurious overcorrection).

      # b3: delta_ki' [Lbar_ff']_t lambda_i = - delta_ki' (Lff_t) lambda_i
      b3[k] <- b3[k] - v2 * sum(delta_ki * as.vector(Lff_t %*% lam_i))

      # b4 (i.i.d.): omega_kt' [Lbar_lam']_i f_t = - omega_kt' (Llam_i) f_t
      b4[k] <- b4[k] - v2 * sum(omega_kt * as.vector(Llam_i %*% f_t))

      # b5: lambda_i' Xi_kt lambda_i  (Xi uses Lbar twice -> sign cancels)
      b5[k] <- b5[k] + v2 * sum(lam_i * as.vector(Xi_arr[[k]][[t]] %*% lam_i))

      # b6 (i.i.d.): f_t' Delta_ki f_t  (Delta uses Lbar twice -> sign cancels)
      b6[k] <- b6[k] + v2 * sum(f_t * as.vector(Delta_arr[[k]][[i]] %*% f_t))
    }
  }

  NT   <- N * TT
  sqNT <- sqrt(NT)
  b3   <- b3 / sqNT
  b4   <- b4 / sqNT
  b5   <- b5 / sqNT
  b6   <- b6 / sqNT

  # ---- b2: dynamic panel term (exog = "weak") ----
  # b2k = (1/sqrt(NT)) sum_i sum_{t < s, both in obs_i, s-t <= L_T}
  #         Gamma((s-t)/L_T) * v_it * x_isk * f_t' [L_lam_inv]_i f_s
  # One-sided kernel: s strictly AFTER t. Uses RAW regressor x_isk.
  b2 <- numeric(p)
  if (exog == "weak") {
    for (i in seq_len(N)) {
      rows_i <- unit_rows[[i]]
      obs_t  <- unit_obs_t[[i]]
      T_i    <- length(obs_t)
      Llam_i <- L_lam_inv[[i]]
      if (anyNA(Llam_i)) next
      for (idx_t in seq_len(T_i - 1L)) {
        t_val <- obs_t[idx_t]
        l_t   <- rows_i[idx_t]
        v_t   <- u_long[l_t]
        f_t   <- F_hat[t_val, ]
        for (idx_s in (idx_t + 1L):T_i) {
          s_val <- obs_t[idx_s]
          gap   <- s_val - t_val     # strictly positive
          if (gap > L_T) next
          kern <- 1 - gap / L_T
          l_s  <- rows_i[idx_s]
          f_s  <- F_hat[s_val, ]
          # f_t' [Lbar_lam']_i f_s = - f_t' (Llam_i) f_s  (single use -> minus sign)
          mid  <- -as.numeric(t(f_t) %*% Llam_i %*% f_s)
          for (k in seq_len(p)) {
            b2[k] <- b2[k] + kern * v_t * X_long[l_s, k] * mid
          }
        }
      }
    }
    b2 <- b2 / sqNT
  }

  # ---- HAC b4 and b6 — overwrite i.i.d. versions when se_type = "hac" ----
  # b4k (HAC) = (1/sqrt(NT)) sum_i sum_{t,s in obs_i, |t-s| <= L_T}
  #               Gamma(|t-s|/L_T) * v_it * v_is * omega_kt' [L_lam_inv]_i f_s
  # b6k (HAC) = (1/sqrt(NT)) sum_i sum_{t,s in obs_i, |t-s| <= L_T}
  #               Gamma(|t-s|/L_T) * v_it * v_is * f_t' Delta_ki f_s
  # Both are two-sided kernels (symmetric in t,s).
  if (se_type == "hac") {
    b4_hac <- numeric(p)
    b6_hac <- numeric(p)
    for (i in seq_len(N)) {
      rows_i <- unit_rows[[i]]
      obs_t  <- unit_obs_t[[i]]
      T_i    <- length(obs_t)
      Llam_i <- L_lam_inv[[i]]
      if (anyNA(Llam_i)) next
      for (idx_t in seq_len(T_i)) {
        t_val <- obs_t[idx_t]
        l_t   <- rows_i[idx_t]
        v_t   <- u_long[l_t]
        f_t   <- F_hat[t_val, ]
        for (idx_s in seq_len(T_i)) {
          gap <- abs(obs_t[idx_s] - obs_t[idx_t])
          if (gap > L_T) next
          kern  <- 1 - gap / L_T
          l_s   <- rows_i[idx_s]
          v_s   <- u_long[l_s]
          s_val <- obs_t[idx_s]
          f_s   <- F_hat[s_val, ]
          for (k in seq_len(p)) {
            omega_kt <- omega_arr[[k]][t_val, ]
            Delta_ki <- Delta_arr[[k]][[i]]
            # b4 uses [Lbar_lam']_i once -> minus sign (matches i.i.d. b4 above)
            b4_hac[k] <- b4_hac[k] - kern * v_t * v_s *
                         sum(omega_kt * as.vector(Llam_i %*% f_s))
            # b6 uses Delta_ki (Lbar twice -> sign cancels) -> no minus
            b6_hac[k] <- b6_hac[k] + kern * v_t * v_s *
                         sum(f_t * as.vector(Delta_ki %*% f_s))
          }
        }
      }
    }
    b4 <- b4_hac / sqNT
    b6 <- b6_hac / sqNT
  }

  b_hat    <- b2 + b3 + b4 + b5 + b6

  # beta_abc = beta - (1/sqrt(NT)) * W_x_inv %*% b_hat
  beta_abc <- beta - as.vector(W_x_inv %*% b_hat) / sqNT

  list(beta_abc = beta_abc, b_hat = b_hat,
       b2 = b2, b3 = b3, b4 = b4, b5 = b5, b6 = b6)
}


# ----------------------------------------------------------------------------
# ife_unbalanced — Public wrapper
#
#' Unbalanced Panel Interactive Fixed Effects Estimator
#'
#' Fits the interactive fixed effects model
#' \deqn{Y_{it} = \alpha_i + \xi_t + X_{it}'\beta + \lambda_i'F_t + u_{it}}
#' for unbalanced panels (units observed at different sets of time periods),
#' where the additive unit effects \eqn{\alpha_i} and time effects \eqn{\xi_t}
#' are controlled by \code{force}. The estimation and inference theory follows
#' Su, Wang and Wang (2025). Estimation uses an alternating outer loop
#' that updates \eqn{\hat\beta} and the structure \eqn{(\hat\alpha, \hat\xi,
#' \hat\lambda, \hat F)}, with an expectation-maximisation (EM) inner loop ---
#' in the spirit of Bai (2009, Appendix B) and the missing-data factor /
#' matrix-completion framework of Bai and Ng (2021) --- that imputes the
#' unobserved cells from the current structure and re-estimates the additive
#' and interactive components on the completed panel. An optional
#' nuclear-norm-regularised (soft-impute) warm start (Mazumder, Hastie and
#' Tibshirani 2010) is available via \code{init = "nnr"}.
#'
#' **Inference.** Standard errors use a sandwich estimator on factor-projected
#' regressors (the unbalanced analogue of the balanced formula of Bai 2009),
#' following Su, Wang and Wang (2025), with heteroskedasticity-robust,
#' cluster-robust (Arellano 1987; Cameron, Gelbach and Miller 2011) and HAC
#' (Newey and West 1987) variants. The optional analytical bias correction of
#' Su, Wang and Wang (2025) removes the leading incidental-parameter bias,
#' extending the corrections of Bai (2009) and Moon and Weidner (2017) to the
#' unbalanced and predetermined-regressor case.
#'
#' **Additive fixed effects.** With \code{force = "none"} (default) the model
#' is intercept-free and all heterogeneity is carried by the interactive
#' factors. Setting \code{force = "unit"}, \code{"time"} or \code{"two-way"}
#' estimates the additive effects jointly with the factors by demeaning the
#' \emph{imputed} (completed) panel inside the EM loop, which is robust to
#' informative missingness; the degrees-of-freedom adjustment is propagated to
#' the standard errors.
#'
#' @param formula R formula: \code{outcome ~ covariate1 + covariate2 + ...}
#' @param data    Data frame in long format (one row per observed unit-time
#'   pair).
#' @param index   Character vector of length 2: \code{c("unit_col", "time_col")}.
#' @param r       Positive integer. Number of interactive factors (default 1).
#' @param force   Additive fixed effects to remove jointly with the factors:
#'   \code{"none"} (default; intercept-free interactive model),
#'   \code{"unit"}, \code{"time"}, or \code{"two-way"}. Additive FE are
#'   estimated jointly with the factors via EM on the imputed (completed)
#'   panel, which is robust to informative (factor-correlated) missingness;
#'   use \code{force = "two-way"} for data with level/trend structure (matching
#'   the balanced \code{ife()} default). Note: with strong, near-collinear
#'   common trends, \code{"two-way"} can converge slowly.
#' @param se      SE type: \code{"standard"} (homoskedastic),
#'   \code{"robust"} (HC1), \code{"cluster"} (cluster-robust by unit), or
#'   \code{"hac"} (HAC with a Bartlett kernel, for serially correlated errors).
#'   Default \code{"standard"}.
#' @param init    Initialisation method: \code{"ols"} (default, grand-mean OLS)
#'   or \code{"nnr"} (nuclear-norm regularisation / soft-impute).
#' @param bias_corr Logical. Apply the analytical incidental-parameter bias
#'   correction. Supports both strictly and weakly exogenous regressors
#'   (controlled by \code{exog}). Default \code{FALSE}.
#' @param exog    Exogeneity assumption: \code{"strict"} (default, regressors
#'   uncorrelated with past and future errors) or \code{"weak"} (weakly
#'   exogenous, e.g., lagged dependent variable \eqn{x_{it} = y_{i,t-1}}).
#'   When \code{"weak"} and \code{bias_corr = TRUE}, an additional dynamic
#'   bias term \eqn{\hat{b}_2} is included.
#' @param L_T     Bartlett kernel bandwidth for HAC standard errors (\code{se
#'   = "hac"}) and the dynamic bias term \eqn{\hat{b}_2} (\code{exog =
#'   "weak"}, \code{bias_corr = TRUE}). If \code{NULL} (default), set to
#'   \eqn{\lfloor 2 T^{1/5} \rfloor} after the panel dimensions are known.
#' @param c_f     Singular-value-thresholding constant (default 0.6) used for
#'   factor-number selection. Used only when \code{init = "nnr"}.
#' @param nu_NT   NNR penalty grid. If \code{NULL} (default), cross-validates
#'   over \code{c * sqrt(max(N, TT))} for \code{c} in \code{c(0.01, 0.1, 1, 10)}.
#' @param tol     Outer-loop convergence tolerance on
#'   \eqn{\max|\hat\beta^{new} - \hat\beta^{old}|}. Default \code{1e-9}.
#' @param max_iter Maximum outer-loop iterations. Default \code{10000L}.
#' @param tol_em  Inner EM convergence tolerance. Default \code{1e-7}.
#' @param max_iter_em Maximum inner EM iterations per outer step.
#'   Default \code{500L}.
#'
#' @return An S3 object of class \code{"ife_unb"} with components:
#' \describe{
#'   \item{coef}{Named p-vector of estimated coefficients \eqn{\hat\beta}
#'     (bias-corrected when \code{bias_corr = TRUE}).}
#'   \item{coef_raw}{Named p-vector of uncorrected coefficients (only when
#'     \code{bias_corr = TRUE}).}
#'   \item{vcov}{p x p variance-covariance matrix.}
#'   \item{se}{Named p-vector of standard errors.}
#'   \item{tstat}{Named p-vector of t-statistics.}
#'   \item{pval}{Named p-vector of two-sided p-values.}
#'   \item{ci}{p x 2 matrix of 95 percent confidence intervals.}
#'   \item{table}{Data frame coefficient table.}
#'   \item{F_hat}{TT x r estimated factor matrix (normalised F'F/TT = I_r).}
#'   \item{Lambda_hat}{N x r estimated loading matrix.}
#'   \item{residuals}{n_obs numeric vector of full-model residuals at
#'     observed cells.}
#'   \item{sigma2}{Estimated error variance (\eqn{sum(u^2)/df}).}
#'   \item{df}{Residual degrees of freedom.}
#'   \item{n_obs}{Number of observed unit-time cells.}
#'   \item{n_iter}{Outer-loop iterations to convergence.}
#'   \item{converged}{Logical.}
#'   \item{N, TT, r, se_type}{Model dimensions and options.}
#'   \item{init, bias_corr, exog, L_T}{Options used.}
#'   \item{b_hat, b2, b3, b4, b5, b6}{Bias components (only when
#'     \code{bias_corr = TRUE}). \code{b2} is a zero vector when
#'     \code{exog = "strict"}.}
#'   \item{y_name, x_names, id_col, time_col}{Variable names.}
#'   \item{unit_vals, time_vals}{Unique unit and time identifiers.}
#'   \item{unit_idx, time_idx}{Integer index vectors for \code{residuals}.}
#'   \item{call}{Matched call.}
#' }
#'
#' @references
#' Su, L., Wang, F. and Wang, Y. (2025). Estimation and inference for
#' interactive fixed effects panel data models with unbalanced panels.
#' SSRN Working Paper No. 5177283. \doi{10.2139/ssrn.5177283}
#'
#' Bai, J. (2009). Panel data models with interactive fixed effects.
#' \emph{Econometrica}, 77(4), 1229--1279. \doi{10.3982/ECTA6135}
#'
#' Bai, J. and Ng, S. (2021). Matrix completion, counterfactuals, and factor
#' analysis of missing data. \emph{Journal of the American Statistical
#' Association}, 116(536), 1746--1763. \doi{10.1080/01621459.2021.1967163}
#'
#' Mazumder, R., Hastie, T. and Tibshirani, R. (2010). Spectral regularization
#' algorithms for learning large incomplete matrices. \emph{Journal of Machine
#' Learning Research}, 11, 2287--2322.
#'
#' Moon, H. R. and Weidner, M. (2017). Dynamic linear panel regression models
#' with interactive fixed effects. \emph{Econometric Theory}, 33, 158--195.
#' \doi{10.1017/S0266466615000328}
#'
#' @importFrom stats pt qt
#' @export
#'
#' @examples
#' data(cigar, package = "xtife")
#' # Drop ~10 % of rows to create an unbalanced panel
#' set.seed(1)
#' cigar_unb <- cigar[sample(nrow(cigar), 1200L), ]
#' fit <- ife_unbalanced(sales ~ price, data = cigar_unb,
#'                       index = c("state", "year"), r = 2L)
#' print(fit)
# ----------------------------------------------------------------------------
ife_unbalanced <- function(formula,
                            data,
                            index,
                            r           = 1L,
                            force       = "none",
                            se          = "standard",
                            init        = "ols",
                            bias_corr   = FALSE,
                            exog        = "strict",
                            L_T         = NULL,
                            c_f         = 0.6,
                            nu_NT       = NULL,
                            tol         = 1e-9,
                            max_iter    = 10000L,
                            tol_em      = 1e-7,
                            max_iter_em = 500L) {

  cl <- match.call()

  # ================================================================
  # Input validation
  # ================================================================
  if (!inherits(formula, "formula"))
    stop("'formula' must be an R formula object.")
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.")
  if (!is.character(index) || length(index) != 2L)
    stop("'index' must be a character vector of length 2: c('unit_col', 'time_col').")

  missing_idx <- setdiff(index, names(data))
  if (length(missing_idx) > 0L)
    stop("'index' column(s) not found in 'data': ",
         paste(missing_idx, collapse = ", "))

  if (!se %in% c("standard", "robust", "cluster", "hac"))
    stop("'se' must be one of: 'standard', 'robust', 'cluster', 'hac'.")

  if (!init %in% c("ols", "nnr"))
    stop("'init' must be 'ols' or 'nnr'.")

  if (!exog %in% c("strict", "weak"))
    stop("'exog' must be 'strict' (default) or 'weak' (dynamic/lagged dep. var.).")

  if (!force %in% c("none", "unit", "time", "two-way"))
    stop("'force' must be one of: 'none', 'unit', 'time', 'two-way'.")

  L_T_arg <- L_T   # store user input; resolved after TT is known

  r <- as.integer(r)
  if (r < 1L)
    stop("'r' must be a positive integer (>= 1). ",
         "For r = 0 (plain OLS), use lm() or plm().")

  tol         <- as.double(tol)
  tol_em      <- as.double(tol_em)
  max_iter    <- as.integer(max_iter)
  max_iter_em <- as.integer(max_iter_em)
  if (tol    <= 0) stop("'tol' must be a positive number.")
  if (tol_em <= 0) stop("'tol_em' must be a positive number.")
  if (max_iter    < 1L) stop("'max_iter' must be >= 1.")
  if (max_iter_em < 1L) stop("'max_iter_em' must be >= 1.")

  # Parse formula
  vars     <- all.vars(formula)
  y_name   <- vars[1L]
  x_names  <- vars[-1L]
  p        <- length(x_names)
  id_col   <- index[1L]
  time_col <- index[2L]

  if (force != "none" && p == 0L)
    stop("force = '", force, "' requires at least one covariate; the ",
         "no-covariate case is not supported (use force = 'none').")

  all_needed <- c(y_name, x_names, id_col, time_col)
  missing_v  <- setdiff(all_needed, names(data))
  if (length(missing_v) > 0L)
    stop("Variable(s) not found in 'data': ",
         paste(missing_v, collapse = ", "))

  dup_key <- paste(data[[id_col]], data[[time_col]], sep = "___IFE___")
  if (anyDuplicated(dup_key))
    stop("Duplicate (unit, time) pairs found. ",
         "Each (i, t) must appear at most once in 'data'.")

  for (v in c(y_name, x_names)) {
    if (anyNA(data[[v]]))
      stop("Missing values (NA) found in variable '", v, "'. ",
           "Structural missing observations should be represented by ",
           "absent rows, not NA values.")
  }

  # ================================================================
  # Prepare data
  # ================================================================
  data <- data[order(data[[id_col]], data[[time_col]]), ]

  unit_vals <- sort(unique(data[[id_col]]))
  time_vals <- sort(unique(data[[time_col]]))
  N  <- length(unit_vals)
  TT <- length(time_vals)

  unit_idx <- match(data[[id_col]],   unit_vals)
  time_idx <- match(data[[time_col]], time_vals)
  n_obs    <- nrow(data)

  obs_per_unit <- tabulate(unit_idx, nbins = N)
  if (any(obs_per_unit < r + 1L))
    stop("Some units have fewer than r + 1 = ", r + 1L,
         " observations. The FWL projection requires T_i > r for all units. ",
         "Reduce r or remove units with too few observations.")

  if (r > min(N, TT))
    stop("r = ", r, " exceeds min(N, TT) = min(", N, ", ", TT, ") = ",
         min(N, TT), ". Reduce r.")

  # ---- Resolve bandwidth L_T (needed for HAC SE and/or b2 bias term) ----
  if (is.null(L_T_arg)) {
    L_T <- floor(2L * TT^(1/5))
  } else {
    L_T <- as.integer(L_T_arg)
    if (L_T < 1L) stop("'L_T' must be a positive integer.")
  }

  Y_long <- as.double(data[[y_name]])
  X_long <- if (p > 0L) {
    as.matrix(data[, x_names, drop = FALSE])
  } else {
    matrix(0, n_obs, 0L)
  }
  storage.mode(X_long) <- "double"

  # ================================================================
  # No grand-mean centering.
  # the interactive-FE model has no intercept; the factor structure absorbs any
  # level.  Grand-mean centering subtracts the OBSERVED-cell mean, which is
  # selection-biased under informative (factor-correlated) missingness
  # (e.g. factor-correlated missingness, d_it ~ Phi(lambda_i'f_t)) and biases beta -- it
  # even sign-flips beta_1 in the dynamic Pattern-2 case.  Removing it makes
  # the estimator match the observed-cell objective (and fect::interFE).
  # ================================================================
  mu_Y     <- 0
  Y_long_c <- Y_long
  mu_X     <- if (p > 0L) numeric(p) else numeric(0L)
  X_long_c <- X_long

  # ================================================================
  # Initialisation: OLS or NNR
  # ================================================================
  beta_init   <- NULL
  F_init      <- NULL
  Lambda_init <- NULL
  nu_used     <- NA_real_

  if (init == "nnr" && p > 0L) {
    obs_lin_init <- (unit_idx - 1L) * TT + time_idx
    c_grid <- if (is.null(nu_NT)) c(0.01, 0.1, 1, 10) else nu_NT / sqrt(max(N, TT))

    nnr <- .ife_nnr_unb(
      Y_long  = Y_long_c,
      X_long  = X_long_c,
      obs_lin = obs_lin_init,
      N = N, TT = TT, r = r,
      c_grid  = c_grid
    )
    beta_init   <- nnr$beta0
    F_init      <- nnr$F_hat_0
    Lambda_init <- nnr$Lambda_hat_0
    nu_used     <- nnr$nu_used
  }

  # ================================================================
  # Core estimation
  # ================================================================
  fit <- .ife_fit_unb(
    Y_long      = Y_long_c,
    X_long      = X_long_c,
    unit_idx    = unit_idx,
    time_idx    = time_idx,
    N = N, TT = TT, r = r,
    force       = force,
    beta_init   = beta_init,
    F_init      = F_init,
    Lambda_init = Lambda_init,
    tol         = tol,
    max_iter    = max_iter,
    tol_em      = tol_em,
    max_iter_em = max_iter_em
  )

  if (!fit$converged)
    warning("Outer loop did not converge after ", max_iter,
            " iterations. Increase max_iter or relax tol.")

  # ================================================================
  # factor-projected regressors, exact SE, and bias correction
  # ================================================================
  coef_vec <- fit$beta
  coef_raw <- NULL      # pre-correction estimate (set below if bias_corr)
  b_hat    <- NULL
  b2 <- b3 <- b4 <- b5 <- b6 <- NULL

  vcov_mat <- matrix(NA_real_, p, p, dimnames = list(x_names, x_names))
  se_vec   <- rep(NA_real_, p); names(se_vec) <- x_names
  tstat    <- rep(NA_real_, p); names(tstat)  <- x_names
  pval     <- rep(NA_real_, p); names(pval)   <- x_names
  ci       <- matrix(NA_real_, p, 2L, dimnames = list(x_names, c("2.5 %", "97.5 %")))
  df_resid <- n_obs - p - r * (N + TT - r)

  if (p > 0L) {

    # ---- Additive-FE augmentation (force != "none") ----
    # Represent additive FE as constant-structure factor columns so the
    # projection machinery (delta/omega, SE, bias correction) removes the
    # additive FE together with the interactive factors:
    #   unit FE alpha_i : factor column 1_T, loading column alpha
    #   time FE xi_t    : factor column xi,  loading column 1_N
    # df is corrected to the TRUE additive-FE dof (not the augmented factor
    # count) further below.
    F_eff   <- fit$F_hat
    Lam_eff <- fit$Lambda_hat
    if (force %in% c("unit", "two-way")) {
      F_eff   <- cbind(F_eff,   rep(1, TT))
      Lam_eff <- cbind(Lam_eff, fit$alpha)
    }
    if (force %in% c("time", "two-way")) {
      F_eff   <- cbind(F_eff,   fit$xi)
      Lam_eff <- cbind(Lam_eff, rep(1, N))
    }
    # Re-orthonormalise the augmented factors so F_eff'F_eff/TT = I (    # Thm 4.2 normalisation), while preserving the fitted structure
    # F_eff %*% t(Lam_eff).  The SE projection is rotation-invariant, but the
    # bias-correction b-terms assume orthonormal factors; the raw cbind columns
    # (1_T, xi) are not orthonormal.  Skipped for force = "none" (fit$F_hat
    # already satisfies F'F/TT = I, leaving that path byte-identical).
    if (force != "none") {
      sv_eff  <- svd(F_eff)
      F_eff   <- sv_eff$u * sqrt(TT)
      Lam_eff <- Lam_eff %*% sv_eff$v %*% diag(sv_eff$d, length(sv_eff$d)) / sqrt(TT)
    }
    r_eff <- ncol(F_eff)

    fe_dof <- switch(force, none = 0L, unit = N - 1L,
                     time = TT - 1L, "two-way" = N + TT - 2L)

    # ---- Step 1: alternating LS for delta and omega ----
    do_res <- .ife_delta_omega_unb(
      X_long     = X_long_c,
      F_hat      = F_eff,
      Lambda_hat = Lam_eff,
      unit_idx   = unit_idx,
      time_idx   = time_idx,
      N = N, TT = TT, r = r_eff, p = p,
      unit_rows  = fit$unit_rows,
      unit_obs_t = fit$unit_obs_t
    )

    # ---- Step 2: auxiliary matrices Delta and Xi ----
    intermed <- .ife_intermediates_unb(
      delta_arr  = do_res$delta_arr,
      omega_arr  = do_res$omega_arr,
      F_hat      = F_eff,
      Lambda_hat = Lam_eff,
      L_ff_inv   = do_res$L_ff_inv,
      L_lam_inv  = do_res$L_lam_inv,
      unit_obs_t = fit$unit_obs_t,
      time_obs_i = do_res$time_obs_i,
      N = N, TT = TT, r = r_eff, p = p
    )

    # ---- Step 3: exact SE ----
    se_list  <- .ife_se_unb(
      beta        = fit$beta,
      X_long      = X_long_c,
      delta_arr   = do_res$delta_arr,
      omega_arr   = do_res$omega_arr,
      F_hat       = F_eff,
      Lambda_hat  = Lam_eff,
      u_long      = fit$u_long,
      unit_idx    = unit_idx,
      time_idx    = time_idx,
      N = N, TT = TT, r = r,
      se_type    = se,
      n_obs      = n_obs,
      unit_rows  = fit$unit_rows,
      unit_obs_t = fit$unit_obs_t,
      L_T        = L_T,
      fe_dof     = fe_dof
    )
    vcov_mat <- se_list$vcov_mat
    df_resid <- se_list$df

    # ---- Step 4: bias correction (optional) ----
    if (bias_corr) {
      bc <- .ife_bias_unb(
        beta       = fit$beta,
        delta_arr  = do_res$delta_arr,
        omega_arr  = do_res$omega_arr,
        F_hat      = F_eff,
        Lambda_hat = Lam_eff,
        u_long     = fit$u_long,
        unit_idx   = unit_idx,
        time_idx   = time_idx,
        L_ff_inv   = do_res$L_ff_inv,
        L_lam_inv  = do_res$L_lam_inv,
        Delta_arr  = intermed$Delta_arr,
        Xi_arr     = intermed$Xi_arr,
        W_x_inv    = se_list$W_x_inv,
        N = N, TT = TT, p = p, n_obs = n_obs,
        X_long     = X_long_c,
        unit_rows  = fit$unit_rows,
        unit_obs_t = fit$unit_obs_t,
        exog       = exog,
        se_type    = se,
        L_T        = L_T
      )
      coef_raw <- coef_vec
      coef_vec <- bc$beta_abc
      b_hat    <- bc$b_hat
      b2       <- bc$b2
      b3       <- bc$b3
      b4       <- bc$b4
      b5       <- bc$b5
      b6       <- bc$b6
    }

    names(coef_vec) <- x_names
    se_vec    <- sqrt(pmax(diag(vcov_mat), 0))
    names(se_vec) <- x_names
    tstat     <- coef_vec / se_vec
    pval      <- 2 * pt(-abs(tstat), df = df_resid)
    t_crit    <- qt(0.975, df = df_resid)
    ci[, 1L]  <- coef_vec - t_crit * se_vec
    ci[, 2L]  <- coef_vec + t_crit * se_vec
  }

  names(coef_vec) <- x_names
  sigma2 <- if (df_resid > 0L) sum(fit$u_long^2) / df_resid else NA_real_

  coef_table <- if (p > 0L) {
    data.frame(
      Estimate  = coef_vec,
      Std.Error = se_vec,
      t.value   = tstat,
      Pr.t      = pval,
      CI.lower  = ci[, 1L],
      CI.upper  = ci[, 2L],
      row.names = x_names,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  # ================================================================
  # Return
  # ================================================================
  out <- list(
    coef        = coef_vec,
    vcov        = vcov_mat,
    se          = se_vec,
    tstat       = tstat,
    pval        = pval,
    ci          = ci,
    table       = coef_table,
    F_hat       = fit$F_hat,
    Lambda_hat  = fit$Lambda_hat,
    residuals   = fit$u_long,
    sigma2      = sigma2,
    df          = df_resid,
    n_obs       = n_obs,
    n_iter      = fit$n_iter,
    converged   = fit$converged,
    N           = N,
    TT          = TT,
    r           = r,
    se_type     = se,
    init        = init,
    bias_corr   = bias_corr,
    exog        = exog,
    L_T         = L_T,
    y_name      = y_name,
    x_names     = x_names,
    id_col      = id_col,
    time_col    = time_col,
    unit_vals   = unit_vals,
    time_vals   = time_vals,
    unit_idx    = unit_idx,
    time_idx    = time_idx,
    call        = cl
  )

  # Append optional bias-correction fields
  if (bias_corr && p > 0L) {
    out$coef_raw <- coef_raw
    names(out$coef_raw) <- x_names
    out$b_hat <- b_hat
    out$b2    <- b2     # zero vector when exog = "strict"
    out$b3    <- b3
    out$b4    <- b4
    out$b5    <- b5
    out$b6    <- b6
  }
  if (init == "nnr") out$nu_used <- nu_used

  structure(out, class = "ife_unb")
}


# ----------------------------------------------------------------------------
# ife_select_r_unb — SVT factor selection for unbalanced panels
#
# Applies a singular value thresholding rule:
#   r_hat = #{s : sigma_s(Theta_hat / sqrt(NT)) >= c_f * sqrt(c_NT^{-1/4} * sigma_1)}
# where Theta_hat is the NNR soft-imputed matrix and c_NT = min(sqrt(N), sqrt(T)).
#
#' Factor Number Selection for Unbalanced Panel IFE via SVT
#'
#' Estimates the number of interactive factors in an unbalanced panel by
#' the singular value thresholding (SVT) rule of Su, Wang and Wang (2025),
#' applied to the nuclear-norm-regularised (soft-imputed) matrix --- a
#' missing-data counterpart of the information-criterion rules of Bai and Ng
#' (2002).
#'
#' @param formula R formula: \code{outcome ~ covariate1 + covariate2 + ...}
#' @param data    Data frame in long format.
#' @param index   Character vector of length 2: \code{c("unit_col", "time_col")}.
#' @param c_f     SVT threshold constant (default 0.6).
#' @param nu_NT   Optional scalar or vector of NNR penalty values. If
#'   \code{NULL} (default), cross-validates over
#'   \code{c(0.01, 0.1, 1, 10) * sqrt(max(N, TT))}.
#' @param verbose Logical; print result table. Default \code{TRUE}.
#'
#' @return Invisibly returns a list with components \code{r_hat}, \code{sv}
#'   (normalised singular values), \code{threshold}, \code{c_f},
#'   \code{c_NT}, and \code{nu_used}.
#'
#' @references
#' Su, L., Wang, F. and Wang, Y. (2025). Estimation and inference for
#' interactive fixed effects panel data models with unbalanced panels.
#' SSRN Working Paper No. 5177283. \doi{10.2139/ssrn.5177283}
#'
#' Bai, J. and Ng, S. (2002). Determining the number of factors in approximate
#' factor models. \emph{Econometrica}, 70(1), 191--221.
#' \doi{10.1111/1468-0262.00273}
#'
#' Bai, J. and Ng, S. (2021). Matrix completion, counterfactuals, and factor
#' analysis of missing data. \emph{Journal of the American Statistical
#' Association}, 116(536), 1746--1763. \doi{10.1080/01621459.2021.1967163}
#'
#' @export
#'
#' @examples
#' \donttest{
#'   data(cigar, package = "xtife")
#'   set.seed(42)
#'   cigar_unb <- cigar[sample(nrow(cigar), 1200L), ]
#'   ife_select_r_unb(sales ~ price, data = cigar_unb,
#'                    index = c("state", "year"))
#' }
# ----------------------------------------------------------------------------
ife_select_r_unb <- function(formula, data, index,
                              c_f     = 0.6,
                              nu_NT   = NULL,
                              verbose = TRUE) {

  # ---- Input validation (abbreviated) ----
  if (!inherits(formula, "formula")) stop("'formula' must be an R formula.")
  if (!is.data.frame(data))          stop("'data' must be a data.frame.")
  if (!is.character(index) || length(index) != 2L)
    stop("'index' must be character(2).")

  missing_idx <- setdiff(index, names(data))
  if (length(missing_idx) > 0L)
    stop("'index' columns not found: ", paste(missing_idx, collapse = ", "))

  vars     <- all.vars(formula)
  y_name   <- vars[1L]
  x_names  <- vars[-1L]
  p        <- length(x_names)
  id_col   <- index[1L]
  time_col <- index[2L]

  missing_v <- setdiff(c(y_name, x_names, id_col, time_col), names(data))
  if (length(missing_v) > 0L)
    stop("Variables not found in data: ", paste(missing_v, collapse = ", "))

  for (v in c(y_name, x_names)) {
    if (anyNA(data[[v]]))
      stop("Missing values (NA) in variable '", v, "'.")
  }

  dup_key <- paste(data[[id_col]], data[[time_col]], sep = "___")
  if (anyDuplicated(dup_key))
    stop("Duplicate (unit, time) pairs found.")

  # ---- Prepare ----
  data <- data[order(data[[id_col]], data[[time_col]]), ]

  unit_vals <- sort(unique(data[[id_col]]))
  time_vals <- sort(unique(data[[time_col]]))
  N  <- length(unit_vals)
  TT <- length(time_vals)

  unit_idx <- match(data[[id_col]],   unit_vals)
  time_idx <- match(data[[time_col]], time_vals)
  n_obs    <- nrow(data)

  Y_long <- as.double(data[[y_name]])
  X_long <- if (p > 0L) {
    as.matrix(data[, x_names, drop = FALSE])
  } else {
    matrix(0, n_obs, 0L)
  }
  storage.mode(X_long) <- "double"

  # No grand-mean centering (see note in the main wrapper above): the
  # observed-cell mean is selection-biased under informative missingness.
  mu_Y     <- 0
  Y_long_c <- Y_long
  mu_X     <- if (p > 0L) numeric(p) else numeric(0L)
  X_long_c <- X_long

  obs_lin <- (unit_idx - 1L) * TT + time_idx
  c_grid  <- if (is.null(nu_NT)) c(0.01, 0.1, 1, 10) else nu_NT / sqrt(max(N, TT))

  # ---- NNR: get Theta_hat_0 (no factor extraction; r = 1 here is a placeholder) ----
  nnr <- .ife_nnr_unb(
    Y_long  = Y_long_c,
    X_long  = X_long_c,
    obs_lin = obs_lin,
    N = N, TT = TT, r = 1L,
    c_grid  = c_grid
  )

  Theta0 <- nnr$Theta0   # TT x N

  # ---- SVT  ----
  sv_norm <- svd(Theta0 / sqrt(N * TT), nu = 0L, nv = 0L)$d
  c_NT    <- min(sqrt(N), sqrt(TT))
  thr     <- c_f * sqrt(c_NT^(-0.25) * sv_norm[1L])
  r_hat   <- max(1L, as.integer(sum(sv_norm >= thr)))

  if (verbose) {
    cat("\nSVT Factor Selection (singular value thresholding)\n")
    cat(strrep("-", 54L), "\n")
    cat(sprintf("N = %d  TT = %d  n_obs = %d  (%.1f%% fill)\n",
                N, TT, n_obs, 100 * n_obs / (N * TT)))
    cat(sprintf("c_f = %.2f  c_NT = %.4f  threshold = %.6f\n",
                c_f, c_NT, thr))
    cat(sprintf("NNR penalty used: %.4f\n\n", nnr$nu_used))

    k_show <- min(15L, length(sv_norm))
    sv_df  <- data.frame(
      s        = seq_len(k_show),
      sigma_s  = round(sv_norm[seq_len(k_show)], 6L),
      selected = sv_norm[seq_len(k_show)] >= thr
    )
    print(sv_df, row.names = FALSE)
    cat(strrep("-", 54L), "\n")
    cat(sprintf("Selected r_hat = %d\n\n", r_hat))
  }

  invisible(list(
    r_hat     = r_hat,
    sv        = sv_norm,
    threshold = thr,
    c_f       = c_f,
    c_NT      = c_NT,
    nu_used   = nnr$nu_used
  ))
}


# ----------------------------------------------------------------------------
# print.ife_unb — S3 print method
#
#' @export
# ----------------------------------------------------------------------------
print.ife_unb <- function(x, digits = 4L, ...) {

  cat("\n")
  cat("Unbalanced Panel Interactive Fixed Effects\n")
  cat(strrep("-", 58L), "\n")
  cat(sprintf("Panel    : N = %d units,  TT = %d periods (max)\n",
              x$N, x$TT))
  cat(sprintf("Observed : n_obs = %d cells  (%.1f%% of N x TT)\n",
              x$n_obs, 100 * x$n_obs / (x$N * x$TT)))
  cat(sprintf("Factors  : r = %d\n", x$r))

  se_label <- switch(x$se_type,
    "standard" = "standard (homoskedastic)",
    "robust"   = "robust (HC1)",
    "cluster"  = paste0("cluster-robust by unit (", x$id_col, ")"),
    "hac"      = paste0("HAC Bartlett (L_T = ", x$L_T, ")")
  )
  cat(sprintf("SE type  : %s\n", se_label))
  cat(sprintf("Exog     : %s\n",
              if (!is.null(x$exog) && x$exog == "weak")
                "weak (dynamic / lagged dep. var.)"
              else "strict (exogenous regressors)"))
  cat(sprintf("Init     : %s%s\n", x$init,
              if (x$init == "nnr" && !is.null(x$nu_used))
                sprintf("  (nu = %.4f)", x$nu_used)
              else ""))
  if (x$bias_corr)
    cat(sprintf("Bias corr: YES (analytical, %s)\n",
                if (!is.null(x$exog) && x$exog == "weak")
                  "weak exogeneity"
                else "strict exogeneity"))
  cat(sprintf("Outcome  : %s\n", x$y_name))
  cat(strrep("-", 58L), "\n")

  if (nrow(x$table) > 0L) {
    tbl   <- x$table
    stars <- ifelse(tbl$Pr.t < 0.01, "***",
              ifelse(tbl$Pr.t < 0.05, "**",
               ifelse(tbl$Pr.t < 0.10, "*", "")))
    out <- data.frame(
      Estimate   = formatC(tbl$Estimate,  digits = digits, format = "f"),
      Std.Error  = formatC(tbl$Std.Error, digits = digits, format = "f"),
      t.value    = formatC(tbl$t.value,   digits = digits, format = "f"),
      "Pr(>|t|)" = formatC(tbl$Pr.t,      digits = digits, format = "f"),
      "95% CI"   = paste0("[",
                     formatC(tbl$CI.lower, digits = digits, format = "f"),
                     ", ",
                     formatC(tbl$CI.upper, digits = digits, format = "f"),
                     "]"),
      " "        = stars,
      row.names  = rownames(tbl),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    print(out, quote = FALSE)
    cat("---\n")
    cat("Signif. codes: *** < 0.01  ** < 0.05  * < 0.10\n")

    if (x$bias_corr && !is.null(x$coef_raw)) {
      has_b2 <- !is.null(x$b2) && !is.null(x$exog) && x$exog == "weak"
      hdr    <- if (has_b2) "b2 + b3 + b4 + b5 + b6" else "b3 + b4 + b5 + b6"
      cat(sprintf("\nBias components (%s):\n", hdr))
      bc_df <- data.frame(
        raw  = formatC(x$coef_raw,           digits = digits, format = "f"),
        b3   = formatC(x$b3,                 digits = digits, format = "f"),
        b4   = formatC(x$b4,                 digits = digits, format = "f"),
        b5   = formatC(x$b5,                 digits = digits, format = "f"),
        b6   = formatC(x$b6,                 digits = digits, format = "f"),
        corr = formatC(x$coef - x$coef_raw,  digits = digits, format = "f"),
        row.names = x$x_names,
        stringsAsFactors = FALSE
      )
      col_names <- c("beta_raw", "b3", "b4", "b5", "b6", "correction")
      if (has_b2) {
        bc_df <- cbind(
          data.frame(b2 = formatC(x$b2, digits = digits, format = "f"),
                     stringsAsFactors = FALSE),
          bc_df
        )
        col_names <- c("beta_raw", "b2", "b3", "b4", "b5", "b6", "correction")
      }
      names(bc_df) <- col_names
      print(bc_df, quote = FALSE)
    }
  } else {
    cat("(No covariates specified)\n")
  }

  cat(strrep("-", 58L), "\n")
  cat(sprintf("sigma^2 = %.6g  |  df = %d\n", x$sigma2, x$df))
  cat(sprintf("Converged: %s  |  Outer iterations: %d\n",
              if (x$converged) "YES" else "NO (increase max_iter)",
              x$n_iter))
  cat("\n")
  invisible(x)
}
