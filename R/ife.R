# ==============================================================================
# 1.  DEMEANING FOR ADDITIVE FIXED EFFECTS
# ==============================================================================

#' @title Demean Panel Matrices for Additive Fixed Effects
#'
#' @description Removes additive unit and/or time fixed effects from a T x N
#' outcome matrix and a T x N x p covariate array using the within-group
#' transformation. Supports four specifications: no demeaning (`"none"`),
#' unit-only (`"unit"`), time-only (`"time"`), and two-way (`"two-way"`).
#'
#' @keywords internal
#'
#' @param Y_mat  T x N outcome matrix
#' @param X_arr  T x N x p covariate array  (NULL if p = 0)
#' @param force  character: "none" | "unit" | "time" | "two-way"
#' @return list with:
#'   Y_dm   T x N demeaned outcome
#'   X_dm   T x N x p demeaned covariate array (or NULL)
#'   mu_Y   grand mean of Y
#'   alpha_Y  N x 1 unit means of Y (after grand mean removal)
#'   xi_Y   T x 1 time means of Y  (after grand mean removal)
.demean_panel <- function(Y_mat, X_arr, force) {
  T <- nrow(Y_mat)
  N <- ncol(Y_mat)
  p <- if (is.null(X_arr)) 0L else dim(X_arr)[3]

  # grand means
  mu_Y <- mean(Y_mat)                      # scalar
  alpha_Y <- colMeans(Y_mat)               # N x 1 (unit means, before subtracting grand)
  xi_Y    <- rowMeans(Y_mat)               # T x 1 (time means, before subtracting grand)

  # demeaning (two-way: within-transformation, preserving balance)
  Y_dm <- switch(force,
    "none"    = Y_mat - mu_Y,
    "unit"    = Y_mat - matrix(alpha_Y, T, N, byrow = TRUE),
    "time"    = Y_mat - matrix(xi_Y,   T, N, byrow = FALSE),
    "two-way" = {
      Y_mat - matrix(alpha_Y, T, N, byrow = TRUE) -
               matrix(xi_Y,   T, N, byrow = FALSE) + mu_Y
    }
  )

  # demean covariates the same way
  X_dm <- NULL
  if (p > 0) {
    X_dm <- X_arr  # copy; will modify in-place below
    for (k in seq_len(p)) {
      Xk <- X_arr[, , k]
      mu_Xk    <- mean(Xk)
      alpha_Xk <- colMeans(Xk)
      xi_Xk    <- rowMeans(Xk)
      X_dm[, , k] <- switch(force,
        "none"    = Xk - mu_Xk,
        "unit"    = Xk - matrix(alpha_Xk, T, N, byrow = TRUE),
        "time"    = Xk - matrix(xi_Xk,   T, N, byrow = FALSE),
        "two-way" = {
          Xk - matrix(alpha_Xk, T, N, byrow = TRUE) -
               matrix(xi_Xk,   T, N, byrow = FALSE) + mu_Xk
        }
      )
    }
  }

  list(
    Y_dm   = Y_dm,
    X_dm   = X_dm,
    mu_Y   = mu_Y,
    alpha_Y = alpha_Y,   # unit means (of raw Y)
    xi_Y   = xi_Y        # time means (of raw Y)
  )
}


# ==============================================================================
# 2.  CORE SVD ALTERNATING PROJECTIONS
# ==============================================================================

#' @title Core IFE Estimation via SVD-Based Alternating Projections
#'
#' @description Estimates the Interactive Fixed Effects model by iterating
#' between an OLS step for the regression coefficients (given current factor
#' estimates) and an SVD step for the factor matrix (given current
#' coefficients). Convergence is declared when the maximum absolute change in
#' the coefficient vector falls below \code{tol}. Handles the degenerate case
#' \code{r = 0} (standard OLS on demeaned data) as a special case. Supports
#' both the strictly-exogenous (\code{"static"}, Bai 2009) and the
#' predetermined-regressor (\code{"dynamic"}, Moon and Weidner 2017) projection
#' schemes.
#'
#' @keywords internal
#'
#' @param Y_dm   T x N demeaned outcome
#' @param X_dm   T x N x p demeaned covariate array (or NULL if p = 0)
#' @param r      integer, number of factors (>= 0)
#' @param tol    convergence tolerance on max |beta_new - beta_old|
#' @param max_iter  maximum iterations
#' @param method  character: "static" (Bai 2009) or "dynamic" (Moon and Weidner 2017)
#' @return list:
#'   beta       p x 1 coefficient vector (numeric(0) if p = 0)
#'   F_hat      T x r factor matrix  (normalized F'F/T = I_r)
#'   Lambda_hat N x r loading matrix
#'   X_tilde    T x N x p factor-projected covariates (for SE computation)
#'   e_mat      T x N residual matrix E = Y_dm - X_dm * beta - F Lambda'
#'   n_iter     number of iterations used
#'   converged  logical
.ife_fit <- function(Y_dm, X_dm, r, tol = 1e-9, max_iter = 10000L,
                    method = "static") {
  T <- nrow(Y_dm)
  N <- ncol(Y_dm)
  p <- if (is.null(X_dm)) 0L else dim(X_dm)[3]

  ## ----- special case: r = 0 (standard TWFE OLS, no interactive FE) -----
  if (r == 0L) {
    if (p == 0L) {
      return(list(
        beta       = numeric(0),
        F_hat      = matrix(0, T, 0),
        Lambda_hat = matrix(0, N, 0),
        X_tilde    = X_dm,  # NULL
        e_mat      = Y_dm,
        n_iter     = 0L,
        converged  = TRUE
      ))
    }
    # stack covariates into (NT x p) design matrix, outcome into NT x 1
    Xlong <- matrix(0, T * N, p)
    for (k in seq_len(p)) Xlong[, k] <- as.vector(X_dm[, , k])
    Ylong <- as.vector(Y_dm)
    beta  <- as.vector(solve(crossprod(Xlong), crossprod(Xlong, Ylong)))
    # residuals (T x N)
    e_mat <- Y_dm
    for (k in seq_len(p)) e_mat <- e_mat - X_dm[, , k] * beta[k]
    return(list(
      beta       = beta,
      F_hat      = matrix(0, T, 0),
      Lambda_hat = matrix(0, N, 0),
      X_tilde    = X_dm,
      e_mat      = e_mat,
      n_iter     = 0L,
      converged  = TRUE
    ))
  }

  ## ----- helper: extract r factors from E (T x N matrix) via SVD/eigen -----
  # Bai (2009): F = sqrt(T) * V_r  where V_r are eigenvectors of EE'/(NT)
  # When T > N it is cheaper to eigen-decompose E'E/(NT) instead.
  extract_factors <- function(E) {
    if (T <= N) {
      # T x T moment matrix; eigenvectors are T-dimensional (= factors F)
      M    <- tcrossprod(E) / (N * T)          # T x T
      eig  <- eigen(M, symmetric = TRUE)
      # eigen() returns eigenvalues in DECREASING order (largest first)
      F_hat <- eig$vectors[, seq_len(r), drop = FALSE] * sqrt(T)  # T x r
    } else {
      # N x N moment matrix; eigenvectors are N-dimensional (= loadings Lambda)
      M      <- crossprod(E) / (N * T)         # N x N
      eig    <- eigen(M, symmetric = TRUE)
      Lambda <- eig$vectors[, seq_len(r), drop = FALSE] * sqrt(N)  # N x r
      F_hat  <- E %*% Lambda / N               # T x r
    }
    F_hat
  }

  ## ----- initialization: OLS beta ignoring factor structure -----
  if (p > 0L) {
    Xlong <- matrix(0, T * N, p)
    for (k in seq_len(p)) Xlong[, k] <- as.vector(X_dm[, , k])
    Ylong <- as.vector(Y_dm)
    beta  <- as.vector(solve(crossprod(Xlong), crossprod(Xlong, Ylong)))
  } else {
    beta <- numeric(0)
  }

  converged <- FALSE
  n_iter    <- 0L

  ## ----- main alternating projections loop -----
  for (iter in seq_len(max_iter)) {
    n_iter <- iter

    # (a) residual matrix E = Y_dm - sum_k X_dm[,,k] * beta[k]
    E <- Y_dm
    if (p > 0L) {
      for (k in seq_len(p)) E <- E - X_dm[, , k] * beta[k]
    }

    # (b) extract r factors from E
    F_hat <- extract_factors(E)   # T x r

    # (c) recover loadings: Lambda = E' F / T,  shape N x r
    Lambda_hat <- crossprod(E, F_hat) / T   # N x r

    # (d) project each covariate out of the factor space (FWL projection)
    # X_tilde_i = X_dm_i - F (F'F)^{-1} F' X_dm_i  for each unit i
    # Equivalently: X_tilde_k = X_dm[,,k] - F %*% solve(t(F)%*%F, t(F) %*% X_dm[,,k])
    # F'F / T ~= I_r by normalization, so (F'F)^{-1} ~= I/T
    FtF_inv <- solve(crossprod(F_hat))  # r x r (exact, not approximation)

    if (p > 0L) {
      # Build projected X: X_tilde is T x N x p
      X_tilde <- array(0, dim = c(T, N, p))
      for (k in seq_len(p)) {
        Xk <- X_dm[, , k]   # T x N
        # For each unit i (column), project out F:
        # X_tilde[,i,k] = Xk[,i] - F %*% (FtF_inv %*% (t(F) %*% Xk[,i]))
        # Vectorized over all units simultaneously:
        #   t(F) %*% Xk  gives  r x N
        #   F %*% FtF_inv %*% (t(F) %*% Xk)  gives  T x N
        FtXk            <- crossprod(F_hat, Xk)        # r x N
        X_tilde[, , k]  <- Xk - F_hat %*% (FtF_inv %*% FtXk)  # T x N
      }

      # (d2) M_Lambda additional projection for dynamic method (Moon & Weidner 2017)
      # Projects X_tilde columns (unit dimension) onto orthogonal complement of Lambda_hat.
      # Double projection: X_tilde -> M_Lambda * X_tilde  (applied as right multiplication
      # since T x N matrices have units as columns).
      # Also applied to Y_tilde for consistency of the OLS numerator.
      if (method == "dynamic") {
        LtL_inv_dyn <- solve(crossprod(Lambda_hat))    # r x r
        # For T x N matrix M: M_lambda_proj(M) = M - (M %*% L)(LtL)^{-1}L'
        LtLiLt <- LtL_inv_dyn %*% t(Lambda_hat)       # r x N
        for (k in seq_len(p)) {
          Xk <- X_tilde[, , k]  # T x N
          X_tilde[, , k] <- Xk - (Xk %*% Lambda_hat) %*% LtLiLt  # T x N
        }
      }

      # (e) update beta: OLS of X_tilde on Y_tilde (both projected)
      FtY      <- crossprod(F_hat, Y_dm)               # r x N
      Y_tilde  <- Y_dm - F_hat %*% (FtF_inv %*% FtY)  # T x N (M_F projection)

      # (e2) also apply M_Lambda to Y_tilde for dynamic method
      if (method == "dynamic") {
        Y_tilde <- Y_tilde - (Y_tilde %*% Lambda_hat) %*% LtLiLt
      }

      Xtlong   <- matrix(0, T * N, p)
      for (k in seq_len(p)) Xtlong[, k] <- as.vector(X_tilde[, , k])
      Ytlong   <- as.vector(Y_tilde)

      beta_new <- as.vector(solve(crossprod(Xtlong), crossprod(Xtlong, Ytlong)))

      # (f) convergence check
      if (max(abs(beta_new - beta)) < tol) {
        beta      <- beta_new
        converged <- TRUE
        break
      }
      beta <- beta_new

    } else {
      # p = 0: no covariates, nothing to iterate over
      X_tilde   <- NULL
      converged <- TRUE
      break
    }
  }  # end for iter

  # final residual matrix (full model)
  e_mat <- Y_dm
  if (p > 0L) {
    for (k in seq_len(p)) e_mat <- e_mat - X_dm[, , k] * beta[k]
  }
  if (r > 0L) {
    e_mat <- e_mat - F_hat %*% t(Lambda_hat)  # subtract F Lambda'
  }

  # if r > 0 but p = 0, X_tilde is still NULL (no beta to estimate)
  if (p == 0L) X_tilde <- NULL

  list(
    beta       = beta,
    F_hat      = F_hat,
    Lambda_hat = Lambda_hat,
    X_tilde    = X_tilde,
    e_mat      = e_mat,
    n_iter     = n_iter,
    converged  = converged
  )
}


# ==============================================================================
# 3.  STANDARD ERROR COMPUTATION
# ==============================================================================

#' @title Compute Standard Errors for the IFE Estimator
#'
#' @description Constructs the sandwich variance-covariance matrix for the IFE
#' coefficient vector using the Frisch-Waugh-Lovell principle. The effective
#' regressors are the factor-projected covariates \eqn{\tilde{X}_{it}} (i.e.,
#' demeaned X after removing the factor space). Three estimators are supported:
#' homoskedastic (`"standard"`), HC1 heteroskedasticity-robust (`"robust"`), and
#' cluster-robust by unit (`"cluster"`), following Cameron, Gelbach and Miller
#' (2011). Degrees of freedom account for regression coefficients, interactive FE
#' parameters, and additive FE parameters.
#'
#' @keywords internal
#'
#' @param beta       p x 1 coefficient vector
#' @param X_tilde    T x N x p projected covariate array
#' @param u_mat      T x N residual matrix from the full model
#' @param N,TT       panel dimensions
#' @param r          number of factors
#' @param force      additive FE specification (for df computation)
#' @param se_type    "standard" | "robust" | "cluster"
#' @return list:
#'   vcov_mat   p x p estimated variance-covariance matrix
#'   df         residual degrees of freedom
.ife_se <- function(beta, X_tilde, u_mat, N, TT, r, force, se_type) {
  p <- length(beta)
  if (p == 0L) stop("Cannot compute SE with no covariates.")

  # ---  degrees of freedom  ---
  # Parameter count: p (beta) + r(N+TT-r) (interactive FE, Bai 2009 Thm 1)
  # plus additive FE parameters
  fe_dof <- switch(force,
    "none"    = 0L,
    "unit"    = N - 1L,
    "time"    = TT - 1L,
    "two-way" = N + TT - 2L
  )
  # total parameters estimated:
  k_total <- p + r * (N + TT - r) + fe_dof
  df      <- N * TT - k_total
  if (df <= 0L) stop(
    "Degrees of freedom <= 0. Reduce r or use a larger panel."
  )

  # ---  A = X_tilde' X_tilde  (p x p) ---
  # Stack X_tilde into NT x p matrix
  Xt_long <- matrix(0, TT * N, p)
  for (k in seq_len(p)) Xt_long[, k] <- as.vector(X_tilde[, , k])
  A <- crossprod(Xt_long)   # p x p

  # residual vector (NT x 1), matching row order of X_tilde
  u_vec <- as.vector(u_mat)  # col-major: obs order is t=1..TT for unit 1, then unit 2, ...

  # ---  variance estimation  ---
  if (se_type == "standard") {
    # homoskedastic: Var(beta) = sigma^2 * A^{-1}
    sigma2   <- sum(u_vec^2) / df
    vcov_mat <- sigma2 * solve(A)

  } else if (se_type == "robust") {
    # HC1: Var(beta) = A^{-1} B A^{-1} * (NT / (NT - p))
    # B = sum_{i,t} u_it^2 * x_tilde_it x_tilde_it'
    # Vectorized: B = t(Xt_long) %*% diag(u^2) %*% Xt_long
    B        <- crossprod(Xt_long * u_vec)   # p x p;  u_vec broadcast row-wise
    A_inv    <- solve(A)
    corr     <- (N * TT) / (N * TT - p)     # HC1 finite-sample correction
    vcov_mat <- A_inv %*% B %*% A_inv * corr

  } else if (se_type == "cluster") {
    # Cluster by unit i.
    # Score for unit i:  psi_i = sum_{t=1}^{TT}  u_it * x_tilde_it  (p x 1)
    # B = sum_i psi_i psi_i'
    # Small-sample correction: [N/(N-1)] * [(NT-1)/(NT-p)]  (Cameron-Gelbach-Miller 2011)
    B <- matrix(0, p, p)
    for (i in seq_len(N)) {
      # rows in Xt_long for unit i: indices (i-1)*TT + 1 : i*TT
      # (col-major vectorization: units are columns, so unit i is column i)
      idx     <- (i - 1L) * TT + seq_len(TT) # TT indices for unit i
      Xt_i    <- Xt_long[idx, , drop = FALSE] # TT x p
      u_i     <- u_vec[idx]                   # TT x 1
      psi_i   <- crossprod(Xt_i, u_i)         # p x 1 (= t(Xt_i) %*% u_i)
      B       <- B + tcrossprod(psi_i)         # p x p outer product
    }
    A_inv    <- solve(A)
    corr     <- (N / (N - 1)) * ((N * TT - 1) / (N * TT - p))
    vcov_mat <- A_inv %*% B %*% A_inv * corr

  } else {
    stop("se_type must be 'standard', 'robust', or 'cluster'.")
  }

  list(vcov_mat = vcov_mat, df = df)
}


# ==============================================================================
# 4.  INFORMATION CRITERIA FOR FACTOR NUMBER SELECTION
# ==============================================================================

#' @title Compute Information Criteria for a Given Number of Factors (Internal)
#'
#' @description Evaluates five information criteria for a fitted IFE model with
#' `r` factors, given the mean squared residual `V_r`. Returns IC1, IC2, and IC3
#' from Bai and Ng (2002) Proposition 1 (ICp1/ICp2/ICp3), applied to IFE
#' residuals as suggested by Bai (2009) Section 9.4, plus a BIC-style criterion
#' (`IC_bic`) and a small-sample-corrected prediction criterion (`PC`) as
#' implemented in the \pkg{fect} package (C++ source, lines 196--224).
#' Called once per candidate `r` inside `ife_select_r()`.
#'
#' @keywords internal
#'
#' @param V_r    scalar: mean squared residual = mean(u_mat^2)  (not df-adjusted)
#' @param r      integer: number of factors
#' @param N,TT   panel dimensions
#' @param p      number of covariates
#' @param force  additive FE spec (for np calculation)
#' @return named list: IC1, IC2, IC3 (Bai & Ng 2002), IC_bic and PC (Bai 2009/fect)
.compute_ic <- function(V_r, r, N, TT, p, force) {
  CNT <- min(sqrt(N), sqrt(TT))

  # --- Bai and Ng (2002) Proposition 1 penalty functions (ICp1/ICp2/ICp3) ---
  # Applied here to IFE residuals per Bai (2009) Section 9.4.
  # g1: ICp1 penalty
  g1 <- (N + TT) / (N * TT) * log(N * TT / (N + TT))
  # g2: uses C_NT as the divergence rate
  g2 <- (N + TT) / (N * TT) * log(CNT^2)
  # g3: slower penalty, tends to select fewer factors
  g3 <- log(CNT^2) / (CNT^2)

  IC1 <- log(V_r) + r * g1
  IC2 <- log(V_r) + r * g2
  IC3 <- log(V_r) + r * g3

  # --- fect BIC-style IC (src/ife.cpp line 206-207) ---
  # np = total parameter count (interactive FE + covariates + additive FE)
  fe_dof <- switch(force,
    "none"    = 0L,
    "unit"    = N - 1L,
    "time"    = TT - 1L,
    "two-way" = N + TT - 2L
  )
  np      <- r * (N + TT - r) + p + 1L + fe_dof
  df_fect <- N * TT - np
  sigma2_adj <- if (df_fect > 0) V_r * N * TT / df_fect else NA_real_

  IC_bic <- if (!is.na(sigma2_adj))
    log(sigma2_adj) + np * log(N * TT) / (N * TT)
  else NA_real_

  # --- fect PC criterion (src/ife.cpp lines 209-224) ---
  # Small-sample inflation factor C (expands penalty when N < 60 or T < 60)
  m1 <- max(0L, 60L - N)
  m2 <- max(0L, 60L - TT)
  C_ss <- (N + m1) * (TT + m2) / (N * TT)

  PC <- if (!is.na(sigma2_adj))
    V_r + r * sigma2_adj * C_ss * (N + TT) / (N * TT) * log(N * TT / (N + TT))
  else NA_real_

  list(IC1 = IC1, IC2 = IC2, IC3 = IC3, IC_bic = IC_bic, PC = PC)
}


# ==============================================================================
# 4b.  BIAS CORRECTION (BAI 2009 SECTION 7)
# ==============================================================================

#' @title Compute Bias-Corrected IFE Coefficients (Bai 2009)
#'
#' @description Applies the two-term asymptotic bias correction from Bai (2009)
#' Theorems 7.1 and 7.2 to the raw IFE coefficient vector:
#' \deqn{\hat{\beta}^\dagger = \hat{\beta} - \hat{B}/N - \hat{C}/T}
#' where \eqn{\hat{B}} corrects for cross-sectional heteroskedasticity
#' (Equation 17) and \eqn{\hat{C}} corrects for time-varying heteroskedasticity
#' (Equation 19). Both terms require \eqn{T/N^2 \to 0} and \eqn{N/T^2 \to 0}
#' respectively (Theorem 7.2). For panels of the scale used in the package
#' examples (\eqn{N \approx 50}, \eqn{T \approx 30}) both conditions hold
#' approximately.
#'
#' @keywords internal
#'
#' @param beta       p-vector of uncorrected IFE coefficients
#' @param F_hat      T x r factor matrix (F'F/T = I_r enforced)
#' @param Lambda_hat N x r loading matrix
#' @param X_dm_arr   T x N x p array of demeaned covariates (after additive FE)
#' @param X_tilde    T x N x p array of factor-projected demeaned X
#' @param e_mat      T x N matrix of full-model residuals (from .ife_fit)
#' @param N,TT,p,r  panel dimensions
#' @return list with elements:
#'   beta_bc  p-vector: bias-corrected coefficients
#'   B_hat    p-vector: estimated cross-section bias term (pre-scaling by 1/N)
#'   C_hat    p-vector: estimated time-heteroskedasticity bias term (pre-scaling by 1/T)
.bias_correct <- function(beta, F_hat, Lambda_hat, X_dm_arr, X_tilde, e_mat,
                          N, TT, p, r) {

  # D0^{-1} = (A/NT)^{-1} where A = X_tilde'X_tilde  (same as in .ife_se)
  Xt_long <- matrix(0, TT * N, p)
  for (k in seq_len(p)) Xt_long[, k] <- as.vector(X_tilde[, , k])
  A      <- crossprod(Xt_long)           # p x p
  D0_inv <- solve(A / (N * TT))          # (N*TT) * solve(A)

  # G = (Lambda'Lambda/N)^{-1},  r x r
  G     <- solve(crossprod(Lambda_hat) / N)
  # A_mat = Lambda G Lambda',  N x N  (symmetric)
  A_mat <- Lambda_hat %*% G %*% t(Lambda_hat)

  # V_arr[t, i, k] = (1/N) * sum_j a_{ij} * X_dm[t, j, k]
  # Matrix form for covariate k: V_k = X_dm_arr[,,k] %*% A_mat / N  (TT x N)
  V_arr <- array(0, dim = c(TT, N, p))
  for (k in seq_len(p)) {
    V_arr[, , k] <- X_dm_arr[, , k] %*% A_mat / N
  }

  # unit error variances: sigma2_i = (1/T) sum_t eps_{it}^2  (N-vector)
  sigma2_i <- colMeans(e_mat^2)          # average over T rows for each of N units

  # --- B_hat (Eq. 17): bias from cross-section heteroskedasticity ---
  sum_B <- numeric(p)
  for (i in seq_len(N)) {
    diff_i  <- X_dm_arr[, i, ] - V_arr[, i, ]              # T x p
    # contribution: (p x T)(T x r)/T = p x r; then (p x r)(r x r)(r x 1)*scalar = p x 1
    # Parenthesise /T before next %*% to avoid R operator-precedence issue
    contrib <- (t(diff_i) %*% F_hat / TT) %*% G %*% Lambda_hat[i, ] * sigma2_i[i]
    sum_B   <- sum_B + as.vector(contrib)
  }
  B_hat <- -D0_inv %*% (sum_B / N)      # p x 1

  # --- C_hat (Eq. 19): bias from time-varying heteroskedasticity ---
  # omega_t = (1/N) sum_k eps_{kt}^2  -- T-vector (time-wise mean squared residual)
  omega <- rowMeans(e_mat^2)
  # M_F * Omega * F = Omega*F - (1/T)*F*(F'*Omega*F)   [T x r]
  # Parenthesise (F_hat / T) to avoid R precedence issue: F/T %*% ... parses wrong
  MF_Omega_F <- (omega * F_hat) -
                (F_hat / TT) %*% (t(F_hat) %*% (omega * F_hat))

  sum_C <- numeric(p)
  for (i in seq_len(N)) {
    X_i     <- X_dm_arr[, i, ]                             # T x p
    # (p x T)(T x r)(r x r)(r x 1) = p x 1
    contrib <- t(X_i) %*% MF_Omega_F %*% G %*% Lambda_hat[i, ]
    sum_C   <- sum_C + as.vector(contrib)
  }
  C_hat <- -D0_inv %*% (sum_C / N)      # p x 1

  beta_bc <- beta - as.vector(B_hat) / N - as.vector(C_hat) / TT

  list(beta_bc = beta_bc,
       B_hat   = as.vector(B_hat),
       C_hat   = as.vector(C_hat))
}


# ==============================================================================
# 4c.  BIAS CORRECTION -- MOON & WEIDNER (2017) DYNAMIC IFE
# ==============================================================================

#' @title Compute Bias-Corrected Coefficients for the Dynamic IFE Estimator
#'
#' @description Applies the three-term asymptotic bias correction from Moon and
#' Weidner (2017) Corollary 4.5 to the raw dynamic IFE coefficient vector:
#' \deqn{\hat{\beta}^* = \hat{\beta} + W^{-1}\!\left(\hat{B}_1/T + \hat{B}_2/N + \hat{B}_3/T\right)}
#' where \eqn{\hat{B}_1} corrects for the Nickell-type bias arising from
#' predetermined (lagged) regressors using a lag-truncation bandwidth `M1`,
#' and \eqn{\hat{B}_2}, \eqn{\hat{B}_3} correct for cross-sectional and
#' time-series heteroskedasticity respectively. The latter two terms are
#' algebraically equivalent to the Bai (2009) \eqn{\hat{B}} and \eqn{\hat{C}}
#' terms.
#'
#' @keywords internal
#'
#' @param beta       p-vector of uncorrected IFE coefficients
#' @param F_hat      T x r factor matrix
#' @param Lambda_hat N x r loading matrix
#' @param X_dm_arr   T x N x p array of demeaned covariates (after additive FE)
#' @param X_tilde    T x N x p double-projected covariates (M_Lambda M_F applied)
#' @param e_mat      T x N matrix of full-model residuals
#' @param N,TT,p,r   panel dimensions
#' @param M1         lag bandwidth for B1 (number of lags to include; default 1)
#' @return list:
#'   beta_bc  p-vector: bias-corrected coefficients
#'   B1_hat   p-vector: dynamic bias term (pre-multiplied by D0_inv; scaled by 1/(NT))
#'   B2_hat   p-vector: cross-section heteroscedasticity bias (same as Bai B_hat)
#'   B3_hat   p-vector: time heteroscedasticity bias (same as Bai C_hat)
.bias_correct_mw <- function(beta, F_hat, Lambda_hat, X_dm_arr, X_tilde,
                              e_mat, N, TT, p, r, M1 = 1L) {

  # --- B2 and B3: reuse Bai (2009) formulae (algebraically equivalent to MW) ---
  # Math-verifier confirmed: trunc(res*res',1,1) in MW B2 keeps diagonal only =
  # sigma2_i used in Bai B_hat.  trunc(res'*res,1,1) in MW B3 (M2=0) keeps
  # diagonal only = omega_t used in Bai C_hat.
  bc_bai <- .bias_correct(beta, F_hat, Lambda_hat, X_dm_arr, X_tilde,
                           e_mat, N, TT, p, r)

  # D0_inv = (A/NT)^{-1} where A = X_tilde'X_tilde (double-projected)
  Xt_long <- matrix(0, TT * N, p)
  for (k in seq_len(p)) Xt_long[, k] <- as.vector(X_tilde[, , k])
  D0_inv <- solve(crossprod(Xt_long) / (N * T))      # p x p

  # P_F = F (F'F)^{-1} F'  (T x T projection onto column space of F)
  FtF_inv <- solve(crossprod(F_hat))                  # r x r
  Pf      <- F_hat %*% (FtF_inv %*% t(F_hat))        # T x T (symmetric)

  # --- B1: dynamic bias from predetermined regressors (MW Eq. 4.4 estimator) ---
  # B1_hat[k] = trace(P_F * cross_trunc_k) / (N*T)
  # cross_k[t,s] = sum_i e_{it} * X_{k,is}  (T x T),
  # cross_trunc_k keeps only entries where 1 <= s-t <= M1 (MW: trunc(..., 0, M1+1))
  # Bias correction: beta_bc += D0_inv %*% B1_hat  (positive addition; math-verifier Q2)
  # The 1/(NT) factor absorbs both MW's 1/sqrt(NT) in B1 and 1/sqrt(NT) in bcorr1,
  # so NO additional /T is needed (math-verifier Q1c).
  B1_vec <- numeric(p)
  for (k in seq_len(p)) {
    X_k    <- X_dm_arr[, , k]                         # T x N
    # cross_k[t,s] = sum_i e_{it} * X_{k,is}  (T x T)
    cross_k <- e_mat %*% t(X_k)                       # (T x N)(N x T) = T x T
    # Keep only upper-triangle entries with lag s-t in {1, ..., M1}
    cross_trunc <- matrix(0, TT, TT)
    for (lag in seq_len(M1)) {
      t_idx <- seq_len(TT - lag)
      cross_trunc[cbind(t_idx, t_idx + lag)] <- cross_k[cbind(t_idx, t_idx + lag)]
    }
    # trace(Pf %*% cross_trunc) = sum(Pf * cross_trunc) since Pf symmetric
    B1_vec[k] <- sum(Pf * cross_trunc) / (N * TT)
  }

  # Apply all three corrections:
  #   B2 and B3 already incorporated in bc_bai$beta_bc (via beta - B_hat/N - C_hat/T)
  #   B1: add D0_inv %*% B1_vec (no /T; the 1/T is absorbed in the 1/(NT) normalization)
  beta_bc <- bc_bai$beta_bc + as.vector(D0_inv %*% B1_vec)

  list(beta_bc = beta_bc,
       B1_hat  = B1_vec,
       B2_hat  = bc_bai$B_hat,   # same as Bai B_hat (cross-section heteroscedasticity)
       B3_hat  = bc_bai$C_hat)   # same as Bai C_hat (time-series heteroscedasticity)
}


# ==============================================================================
# 5.  USER-FACING WRAPPER
# ==============================================================================

#' @title Estimate Interactive Fixed Effects Model (Bai 2009)
#'
#' @description Fits the panel model
#' \deqn{y_{it} = \alpha_i + \xi_t + X_{it}'\beta + \lambda_i'F_t + u_{it}}
#' for balanced panel data with analytical standard errors.
#'
#' @param formula  R formula: `outcome ~ covariate1 + covariate2 + ...`
#' @param data     data.frame in long format (one row per unit-time observation)
#' @param index    character(2): `c("unit_id_column", "time_id_column")`
#' @param r        integer >= 0, number of interactive factors (default 1)
#' @param force    additive FE specification: `"none"` | `"unit"` | `"time"` |
#'   `"two-way"` (default `"two-way"`)
#' @param se       SE type: `"standard"` | `"robust"` | `"cluster"`
#'   (default `"standard"`; `"cluster"` clusters by unit id)
#' @param bias_corr logical; if `TRUE` apply bias correction. For
#'   `method = "static"` uses the two-term Bai (2009) Sec. 7 correction
#'   (B/N + C/T). For `method = "dynamic"` uses the three-term Moon and
#'   Weidner (2017) correction (B1/T + B2/N + B3/T). Requires r > 0 and
#'   at least one covariate. (default `FALSE`)
#' @param method   `"static"` (default) for Bai (2009) strictly-exogenous
#'   regressors; `"dynamic"` for Moon and Weidner (2017) predetermined
#'   regressors (e.g. lagged dependent variable). The dynamic estimator uses
#'   double projection M_Lambda M_F on X in the SVD loop.
#' @param M1       integer; lag bandwidth for the B1 dynamic bias term
#'   (default `1L`). Only used when `method = "dynamic"` and
#'   `bias_corr = TRUE`.
#' @param tol      convergence tolerance (default `1e-9`)
#' @param max_iter maximum iterations (default `10000L`)
#'
#' @return An S3 object of class `"ife"` with the following components:
#'
#' * `coef` -- named p-vector of estimated coefficients
#' * `vcov` -- p x p variance-covariance matrix
#' * `se` -- named p-vector of standard errors
#' * `tstat` -- named p-vector of t-statistics
#' * `pval` -- named p-vector of two-sided p-values
#' * `ci` -- p x 2 matrix of 95% confidence intervals (CI.lower, CI.upper)
#' * `table` -- data.frame coefficient table (Estimate, Std.Error, t.value, Pr.t, CI.lower, CI.upper)
#' * `F_hat` -- T x r estimated factor matrix
#' * `Lambda_hat` -- N x r estimated loading matrix
#' * `residuals` -- T x N residual matrix (full model)
#' * `sigma2` -- estimated error variance
#' * `df` -- residual degrees of freedom
#' * `n_iter` -- iterations to convergence
#' * `converged` -- logical
#' * `N`, `T`, `r`, `force`, `se_type` -- model dimensions and options
#' * `call` -- matched call
#'
#' @references
#' Bai, J. (2009). Panel data models with interactive fixed effects.
#' *Econometrica*, 77(4), 1229--1279. \doi{10.3982/ECTA6135}
#'
#' Moon, H.R. and Weidner, M. (2017). Dynamic linear panel regression models
#' with interactive fixed effects. *Econometric Theory*, 33, 158--195.
#' \doi{10.1017/S0266466615000328}
#'
#' Bai, J. and Ng, S. (2002). Determining the number of factors in approximate
#' factor models. *Econometrica*, 70(1), 191--221. \doi{10.1111/1468-0262.00273}
#'
#' @importFrom stats lm.fit qt pt var
#' @export
#'
#' @examples
#' data(cigar, package = "xtife")
#' fit <- ife(sales ~ price, data = cigar, index = c("state", "year"),
#'            r = 2, force = "two-way", se = "standard")
#' print(fit)
ife <- function(formula,
                data,
                index,
                r          = 1L,
                force      = "two-way",
                se         = "standard",
                bias_corr  = FALSE,
                method     = "static",
                M1         = 1L,
                tol        = 1e-9,
                max_iter   = 10000L) {

  cl <- match.call()

  # ---- input validation ----
  if (!inherits(formula, "formula"))
    stop("'formula' must be an R formula object.")
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.")
  if (length(index) != 2 || !all(index %in% names(data)))
    stop("'index' must be a character vector of length 2 naming columns in 'data'.")
  if (!force %in% c("none", "unit", "time", "two-way"))
    stop("'force' must be one of: 'none', 'unit', 'time', 'two-way'.")
  if (!se %in% c("standard", "robust", "cluster"))
    stop("'se' must be one of: 'standard', 'robust', 'cluster'.")
  if (!method %in% c("static", "dynamic"))
    stop("'method' must be 'static' or 'dynamic'.")
  M1 <- as.integer(M1)
  if (M1 < 1L) stop("'M1' must be a positive integer.")
  r <- as.integer(r)
  if (r < 0L) stop("'r' must be a non-negative integer.")

  id_col   <- index[1]
  time_col <- index[2]

  # ---- parse formula ----
  vars     <- all.vars(formula)
  y_name   <- vars[1]
  x_names  <- vars[-1]
  p        <- length(x_names)

  all_needed <- c(y_name, x_names, id_col, time_col)
  missing_v  <- setdiff(all_needed, names(data))
  if (length(missing_v) > 0)
    stop("Variables not found in data: ", paste(missing_v, collapse = ", "))

  # ---- sort and check balance ----
  data <- data[order(data[[id_col]], data[[time_col]]), ]
  unit_vals <- unique(data[[id_col]])
  time_vals <- unique(data[[time_col]])
  N <- length(unit_vals)
  T <- length(time_vals)
  if (nrow(data) != N * T)
    stop("Panel is not balanced: nrow(data) = ", nrow(data),
         " but N*T = ", N * T, ". Only balanced panels are supported.")

  # ---- check for missing values ----
  for (v in c(y_name, x_names)) {
    if (any(is.na(data[[v]])))
      stop("Missing values in variable '", v,
           "'. Remove rows with NA before calling ife().")
  }

  # ---- covariate variation checks ----
  for (v in x_names) {
    unit_var <- tapply(data[[v]], data[[id_col]], var, na.rm = TRUE)
    if (all(unit_var == 0, na.rm = TRUE))
      stop("Variable '", v,
           "' has no within-unit variation. Remove it or use force = 'unit'.")
    time_var <- tapply(data[[v]], data[[time_col]], var, na.rm = TRUE)
    if (all(time_var == 0, na.rm = TRUE))
      stop("Variable '", v,
           "' has no within-time variation. Remove it or use force = 'time'.")
  }

  # check r does not exceed min(N,T)
  if (r > min(N, T))
    stop("r = ", r, " exceeds min(N, T) = ", min(N, T), ".")

  # ---- reshape to matrix / array form ----
  # Orientation: rows = time, columns = units (T x N)
  Y_mat <- matrix(data[[y_name]], nrow = T, ncol = N)
  colnames(Y_mat) <- as.character(unit_vals)
  rownames(Y_mat) <- as.character(time_vals)

  X_arr <- NULL
  if (p > 0L) {
    X_arr <- array(0, dim = c(T, N, p),
                   dimnames = list(as.character(time_vals),
                                   as.character(unit_vals),
                                   x_names))
    for (k in seq_len(p))
      X_arr[, , k] <- matrix(data[[x_names[k]]], nrow = T, ncol = N)
  }

  # ---- demean for additive FE ----
  dm  <- .demean_panel(Y_mat, X_arr, force)
  Y_dm <- dm$Y_dm
  X_dm <- dm$X_dm

  # ---- core estimation ----
  fit <- .ife_fit(Y_dm, X_dm, r, tol, max_iter, method = method)
  if (!fit$converged)
    warning("IFE algorithm did not converge after ", max_iter,
            " iterations. Results may be unreliable. Increase max_iter or relax tol.")

  # ---- recover full-model residuals (including additive FE) ----
  # Full model: y_it = alpha_i + xi_t + X_it' beta + lambda_i' F_t + u_it
  # Residual:   u_it = Y_it - alpha_i - xi_t - X_it' beta - lambda_i' F_t
  # Since Y_dm = Y - alpha - xi (plus grand mean adjustments), we have:
  #   u_mat (from .ife_fit) = Y_dm - X_dm beta - F Lambda' = full-model residuals
  # (The additive FE are absorbed into the demeaning step.)
  u_mat  <- fit$e_mat                  # T x N, full-model residuals
  sigma2 <- sum(u_mat^2) / (N * T)    # raw error variance (= V_r for IC)

  # ---- information criteria (Bai 2009 + fect) ----
  ic <- .compute_ic(V_r = sigma2, r = r, N = N, TT = T, p = p, force = force)

  # ---- standard errors ----
  se_list <- NULL
  if (p > 0L) {
    se_list <- .ife_se(
      beta    = fit$beta,
      X_tilde = fit$X_tilde,
      u_mat   = u_mat,
      N = N, TT = T, r = r,
      force   = force,
      se_type = se
    )
  }

  # ---- bias correction ----
  # Static (Bai 2009 Sec. 7):  beta_dagger = beta - B/N - C/T
  # Dynamic (Moon & Weidner 2017): beta* = beta + W^{-1}(B1 + B2/N + C/T)
  bias_applied <- FALSE
  coef_raw     <- NULL
  B_hat_out    <- NULL   # Bai B_hat (static) or B2_hat (dynamic)
  C_hat_out    <- NULL   # Bai C_hat (static) or B3_hat (dynamic)
  B1_hat_out   <- NULL   # dynamic only
  B2_hat_out   <- NULL   # dynamic only
  B3_hat_out   <- NULL   # dynamic only

  if (isTRUE(bias_corr) && r > 0L && p > 0L) {
    if (method == "static") {
      bc <- .bias_correct(
        beta       = fit$beta,
        F_hat      = fit$F_hat,
        Lambda_hat = fit$Lambda_hat,
        X_dm_arr   = X_dm,
        X_tilde    = fit$X_tilde,
        e_mat      = u_mat,
        N = N, TT = T, p = p, r = r
      )
      B_hat_out <- bc$B_hat;  names(B_hat_out) <- x_names
      C_hat_out <- bc$C_hat;  names(C_hat_out) <- x_names
    } else {
      # Dynamic: Moon & Weidner (2017) three-term correction
      bc <- .bias_correct_mw(
        beta       = fit$beta,
        F_hat      = fit$F_hat,
        Lambda_hat = fit$Lambda_hat,
        X_dm_arr   = X_dm,
        X_tilde    = fit$X_tilde,
        e_mat      = u_mat,
        N = N, TT = T, p = p, r = r, M1 = M1
      )
      B1_hat_out <- bc$B1_hat;  names(B1_hat_out) <- x_names
      B2_hat_out <- bc$B2_hat;  names(B2_hat_out) <- x_names
      B3_hat_out <- bc$B3_hat;  names(B3_hat_out) <- x_names
    }

    coef_raw  <- fit$beta
    names(coef_raw) <- x_names

    # recompute residuals: u_bc[t,i] = u[t,i] + X_dm[t,i,] %*% (beta - beta_bc)
    d_beta  <- fit$beta - bc$beta_bc               # p-vector of change
    delta_u <- matrix(0, T, N)
    for (k in seq_len(p)) delta_u <- delta_u + X_dm[, , k] * d_beta[k]
    u_mat_bc <- u_mat + delta_u

    # recompute SEs with bias-corrected residuals and coefficients
    se_list <- .ife_se(
      beta    = bc$beta_bc,
      X_tilde = fit$X_tilde,
      u_mat   = u_mat_bc,
      N = N, TT = T, r = r,
      force   = force,
      se_type = se
    )

    # overwrite quantities used to build coef_table below
    fit$beta <- bc$beta_bc
    u_mat    <- u_mat_bc
    sigma2   <- sum(u_mat^2) / (N * T)

    bias_applied <- TRUE
  }

  # ---- build output ----
  coef_vec <- fit$beta
  names(coef_vec) <- x_names

  if (!is.null(se_list)) {
    vcov_mat <- se_list$vcov_mat
    rownames(vcov_mat) <- colnames(vcov_mat) <- x_names
    df       <- se_list$df
    se_vec   <- sqrt(diag(vcov_mat))
    names(se_vec) <- x_names
    tstat    <- coef_vec / se_vec
    pval     <- 2 * pt(-abs(tstat), df = df)
    ci       <- cbind(
      coef_vec + qt(0.025, df = df) * se_vec,
      coef_vec + qt(0.975, df = df) * se_vec
    )
    colnames(ci) <- c("CI.lower", "CI.upper")
    rownames(ci) <- x_names

    coef_table <- data.frame(
      Estimate   = coef_vec,
      Std.Error  = se_vec,
      t.value    = tstat,
      Pr.t       = pval,
      CI.lower   = ci[, 1],
      CI.upper   = ci[, 2],
      row.names  = x_names,
      stringsAsFactors = FALSE
    )
  } else {
    vcov_mat   <- matrix(numeric(0), 0, 0)
    df         <- as.integer(N * T)
    se_vec     <- numeric(0)
    tstat      <- numeric(0)
    pval       <- numeric(0)
    ci         <- matrix(numeric(0), 0, 2)
    coef_table <- data.frame()
  }

  structure(
    list(
      coef       = coef_vec,
      vcov       = vcov_mat,
      se         = se_vec,
      tstat      = tstat,
      pval       = pval,
      ci         = ci,
      table      = coef_table,
      F_hat      = fit$F_hat,
      Lambda_hat = fit$Lambda_hat,
      residuals    = u_mat,
      sigma2       = sigma2,
      ic           = ic,          # list: IC1, IC2, IC3, IC_bic, PC
      df           = df,
      n_iter       = fit$n_iter,
      converged    = fit$converged,
      bias_applied = bias_applied,
      coef_raw     = coef_raw,    # uncorrected beta (NULL if bias_corr = FALSE)
      B_hat        = B_hat_out,   # static: cross-section bias; dynamic: NULL
      C_hat        = C_hat_out,   # static: time-series bias;   dynamic: NULL
      B1_hat       = B1_hat_out,  # dynamic only: predetermined-regressor bias
      B2_hat       = B2_hat_out,  # dynamic only: cross-section heteroscedasticity
      B3_hat       = B3_hat_out,  # dynamic only: time-series heteroscedasticity
      N            = N,
      T            = T,
      r            = r,
      force        = force,
      se_type      = se,
      method       = method,
      M1           = M1,
      y_name       = y_name,
      x_names      = x_names,
      id_col       = id_col,
      time_col     = time_col,
      call         = cl
    ),
    class = "ife"
  )
}


# ==============================================================================
# 6.  PRINT METHOD
# ==============================================================================

#' @title Print an IFE Model Summary
#'
#' @description Prints a formatted summary of an object of class `"ife"`,
#' including panel dimensions, number of factors, additive fixed effect
#' specification, SE type, and a coefficient table with standard errors,
#' t-statistics, p-values, and 95% confidence intervals. If bias correction
#' was applied, bias terms are also reported. Information criteria are printed
#' when the object contains them (i.e., when called from `ife_select_r()`).
#'
#' @param x      an object of class `"ife"`
#' @param digits number of significant digits (default 4)
#' @param ...    unused
#'
#' @return `x` invisibly.
#' @export
#'
#' @examples
#' data(cigar, package = "xtife")
#' fit <- ife(sales ~ price, data = cigar, index = c("state", "year"),
#'            r = 2, force = "two-way", se = "standard")
#' print(fit)
print.ife <- function(x, digits = 4, ...) {
  cat("\n")
  if (isTRUE(x$method == "dynamic"))
    cat("Interactive Fixed Effects -- Dynamic (Moon & Weidner 2017, ET)\n")
  else
    cat("Interactive Fixed Effects (Bai 2009, Econometrica)\n")
  cat(rep("-", 55), "\n", sep = "")

  # model summary line
  force_label <- switch(x$force,
    "none"    = "none",
    "unit"    = "unit",
    "time"    = "time",
    "two-way" = "two-way"
  )
  se_label <- switch(x$se_type,
    "standard" = "standard (homoskedastic)",
    "robust"   = "robust (HC1)",
    "cluster"  = paste0("cluster (by ", x$id_col, ")")
  )
  cat(sprintf("Panel    : N = %d units, T = %d periods\n", x$N, x$T))
  cat(sprintf("Factors  : r = %d\n", x$r))
  cat(sprintf("Force    : %s fixed effects\n", force_label))
  cat(sprintf("SE type  : %s\n", se_label))
  cat(sprintf("Outcome  : %s\n", x$y_name))
  cat(rep("-", 55), "\n", sep = "")

  if (nrow(x$table) > 0) {
    # format the table
    tbl <- x$table
    stars <- ifelse(tbl$Pr.t < 0.01,  "***",
              ifelse(tbl$Pr.t < 0.05,  "**",
               ifelse(tbl$Pr.t < 0.1,   "*", "")))
    print_tbl <- data.frame(
      Estimate   = formatC(tbl$Estimate,  digits = digits, format = "f"),
      Std.Error  = formatC(tbl$Std.Error, digits = digits, format = "f"),
      t.value    = formatC(tbl$t.value,   digits = digits, format = "f"),
      `Pr(>|t|)` = formatC(tbl$Pr.t,      digits = digits, format = "f"),
      `95% CI`   = paste0("[",
                           formatC(tbl$CI.lower, digits = digits, format = "f"),
                           ", ",
                           formatC(tbl$CI.upper, digits = digits, format = "f"),
                           "]"),
      ` ` = stars,
      row.names  = rownames(tbl),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    print(print_tbl, quote = FALSE)
    cat("---\n")
    cat("Signif. codes: *** <0.01  ** <0.05  * <0.1\n")
  } else {
    cat("(No covariates specified)\n")
  }

  cat(rep("-", 55), "\n", sep = "")
  cat(sprintf("sigma^2  : %.6f | df = %d\n", x$sigma2, x$df))
  cat(sprintf("Converged: %s | Iterations: %d\n",
              if (x$converged) "YES" else "NO (increase max_iter)",
              x$n_iter))

  # ---- bias correction summary ----
  if (isTRUE(x$bias_applied)) {
    cat(rep("-", 55), "\n", sep = "")
    if (isTRUE(x$method == "dynamic")) {
      cat("Bias correction (Moon & Weidner 2017): beta* = beta + W^{-1}(B1/T + B2/N + B3/T)\n")
      cat(sprintf("  Method: dynamic  N=%d  T=%d  M1=%d\n", x$N, x$T, x$M1))
      for (nm in x$x_names) {
        cat(sprintf(
          "  %-12s raw=%8.4f  B1/T=%9.6f  B2/N=%9.6f  B3/T=%9.6f  corrected=%8.4f\n",
          nm,
          x$coef_raw[nm],
          x$B1_hat[nm] / x$T,
          x$B2_hat[nm] / x$N,
          x$B3_hat[nm] / x$T,
          x$coef[nm]))
      }
    } else {
      cat("Bias correction (Bai 2009 Sec. 7): beta^ = beta_raw - B/N - C/T\n")
      cat(sprintf("  Conditions: T/N=%.3f  T/N^2=%.5f  N/T^2=%.5f\n",
                  x$T / x$N, x$T / x$N^2, x$N / x$T^2))
      for (nm in x$x_names) {
        cat(sprintf("  %-12s raw=%8.4f  B/N=%9.6f  C/T=%9.6f  corrected=%8.4f\n",
                    nm,
                    x$coef_raw[nm],
                    x$B_hat[nm] / x$N,
                    x$C_hat[nm] / x$T,
                    x$coef[nm]))
      }
    }
  }

  # ---- factor selection guidance ----
  cat(rep("-", 55), "\n", sep = "")
  cat(sprintf(
    "Factor selection criteria at r = %d [IC1-3: Bai & Ng 2002; IC_bic/PC: Bai 2009]:\n",
    x$r))
  cat(sprintf("  IC1 = %8.4f  |  IC2 = %8.4f  |  IC3 = %8.4f\n",
              x$ic$IC1, x$ic$IC2, x$ic$IC3))
  cat(sprintf("  PC  = %8.4f  |  IC (BIC-style) = %8.4f\n",
              x$ic$PC, x$ic$IC_bic))
  cat("  -> Run ife_select_r() to compare criteria across r = 0, 1, ..., r_max\n")
  cat("     and identify the IC-minimising number of factors.\n")

  cat("\n")
  invisible(x)
}


# ==============================================================================
# 7.  FACTOR NUMBER SELECTION
# ==============================================================================

#' @title Select the Number of Factors via Information Criteria
#'
#' @description Fits the IFE model for r = 0, 1, ..., `r_max` and evaluates
#' five information criteria at each value of r. Returns IC1, IC2, and IC3
#' from Bai and Ng (2002) Proposition 1, applied to IFE residuals per Bai
#' (2009) Section 9.4, plus a BIC-style penalty (`IC_bic`) and a
#' small-sample-corrected prediction criterion (`PC`) from Bai (2009).
#' The criterion-minimising r for each IC is flagged with `"*"` in the
#' printed table, and a data-driven recommendation (favouring `IC_bic` when
#' the Bai-Ng criteria decrease monotonically) is displayed.
#'
#' @param formula  R formula passed to `ife()`
#' @param data     long-format data.frame
#' @param index    character(2): `c("unit_id", "time_id")`
#' @param r_max    maximum r to consider (default: `min(8, floor(min(N,T)/2))`)
#' @param force    additive FE type (default `"two-way"`)
#' @param verbose  logical; if `TRUE` (default) print progress and results table
#'   to the console. Set to `FALSE` for silent operation.
#' @param tol      convergence tolerance (default `1e-9`)
#' @param max_iter maximum iterations (default `10000L`)
#'
#' @return (invisibly) a data.frame with columns `r`, `V_r`, `IC1`, `IC2`,
#'   `IC3`, `IC_bic`, `PC`, `converged`, and attribute `"suggested"` (named
#'   integer vector giving the IC-minimising r for each criterion).
#'
#' @references
#' Bai, J. (2009). Panel data models with interactive fixed effects.
#' *Econometrica*, 77(4), 1229--1279. \doi{10.3982/ECTA6135}
#'
#' Bai, J. and Ng, S. (2002). Determining the number of factors in approximate
#' factor models. *Econometrica*, 70(1), 191--221. \doi{10.1111/1468-0262.00273}
#'
#' @export
#'
#' @examples
#' \donttest{
#'   data(cigar, package = "xtife")
#'   sel <- ife_select_r(sales ~ price, data = cigar,
#'                       index = c("state", "year"), r_max = 4)
#' }
ife_select_r <- function(formula,
                         data,
                         index,
                         r_max    = NULL,
                         force    = "two-way",
                         verbose  = TRUE,
                         tol      = 1e-9,
                         max_iter = 10000L) {

  # ---- input validation (minimal; full validation done inside ife()) ----
  if (!force %in% c("none", "unit", "time", "two-way"))
    stop("'force' must be one of: 'none', 'unit', 'time', 'two-way'.")
  if (length(index) != 2 || !all(index %in% names(data)))
    stop("'index' must name two columns in 'data'.")

  id_col   <- index[1]
  time_col <- index[2]

  # ---- parse formula + basic data prep (mirrors ife() internals) ----
  vars    <- all.vars(formula)
  y_name  <- vars[1]
  x_names <- vars[-1]
  p       <- length(x_names)

  data <- data[order(data[[id_col]], data[[time_col]]), ]
  N    <- length(unique(data[[id_col]]))
  T    <- length(unique(data[[time_col]]))

  if (nrow(data) != N * T)
    stop("Panel is not balanced (nrow != N*T). Only balanced panels supported.")

  # default r_max: half of the smaller dimension, capped at 8
  if (is.null(r_max))
    r_max <- min(8L, floor(min(N, T) / 2L))
  r_max <- as.integer(r_max)
  if (r_max < 1L) stop("r_max must be at least 1.")
  if (r_max > min(N, T))
    stop("r_max = ", r_max, " exceeds min(N, T) = ", min(N, T), ".")

  # ---- reshape ----
  unit_vals <- unique(data[[id_col]])
  time_vals <- unique(data[[time_col]])
  Y_mat <- matrix(data[[y_name]], nrow = T, ncol = N)
  X_arr <- NULL
  if (p > 0L) {
    X_arr <- array(0, dim = c(T, N, p))
    for (k in seq_len(p))
      X_arr[, , k] <- matrix(data[[x_names[k]]], nrow = T, ncol = N)
  }

  # ---- demean once (same for all r) ----
  dm   <- .demean_panel(Y_mat, X_arr, force)
  Y_dm <- dm$Y_dm
  X_dm <- dm$X_dm

  # ---- fit for r = 0, 1, ..., r_max ----
  r_vals   <- 0L:r_max
  n_r      <- length(r_vals)
  V_r_vec  <- numeric(n_r)
  IC1_vec  <- numeric(n_r)
  IC2_vec  <- numeric(n_r)
  IC3_vec  <- numeric(n_r)
  ICb_vec  <- numeric(n_r)
  PC_vec   <- numeric(n_r)
  conv_vec <- logical(n_r)

  if (verbose)
    cat(sprintf("Fitting IFE for r = 0 to %d  (force = '%s') ...\n", r_max, force))

  for (i in seq_along(r_vals)) {
    r_i <- r_vals[i]
    fit_i <- .ife_fit(Y_dm, X_dm, r_i, tol, max_iter)
    if (!fit_i$converged)
      warning("r = ", r_i, ": algorithm did not converge.")
    conv_vec[i] <- fit_i$converged
    V_r_i       <- mean(fit_i$e_mat^2)
    V_r_vec[i]  <- V_r_i
    ic_i        <- .compute_ic(V_r = V_r_i, r = r_i, N = N, TT = T, p = p, force = force)
    IC1_vec[i]  <- ic_i$IC1
    IC2_vec[i]  <- ic_i$IC2
    IC3_vec[i]  <- ic_i$IC3
    ICb_vec[i]  <- ic_i$IC_bic
    PC_vec[i]   <- ic_i$PC
    if (verbose)
      cat(sprintf("  r = %d  V(r) = %.5f  IC1 = %8.4f  IC2 = %8.4f  IC3 = %8.4f\n",
                  r_i, V_r_i, ic_i$IC1, ic_i$IC2, ic_i$IC3))
  }

  # ---- IC-minimising r for each criterion ----
  r_IC1 <- r_vals[which.min(IC1_vec)]
  r_IC2 <- r_vals[which.min(IC2_vec)]
  r_IC3 <- r_vals[which.min(IC3_vec)]
  r_ICb <- r_vals[which.min(ICb_vec)]
  r_PC  <- r_vals[which.min(PC_vec)]

  # ---- print formatted table ----
  if (verbose) {
    cat("\n")
    cat(rep("=", 72), "\n", sep = "")
    cat("  Factor Number Selection (Bai 2009)  ---  Criterion-minimising r is starred\n")
    cat(rep("=", 72), "\n", sep = "")
    hdr <- sprintf("  %4s  %10s  %9s  %9s  %9s  %9s  %9s\n",
                   "r", "V(r)", "IC1", "IC2", "IC3", "IC(BIC)", "PC")
    cat(hdr)
    cat(rep("-", 72), "\n", sep = "")

    for (i in seq_along(r_vals)) {
      r_i <- r_vals[i]
      # mark minimum of each criterion
      star1 <- if (r_i == r_IC1) "*" else " "
      star2 <- if (r_i == r_IC2) "*" else " "
      star3 <- if (r_i == r_IC3) "*" else " "
      starb <- if (r_i == r_ICb) "*" else " "
      starP <- if (r_i == r_PC)  "*" else " "
      cat(sprintf("  %4d  %10.5f  %8.4f%s  %8.4f%s  %8.4f%s  %8.4f%s  %8.4f%s\n",
                  r_i, V_r_vec[i],
                  IC1_vec[i], star1, IC2_vec[i], star2, IC3_vec[i], star3,
                  ICb_vec[i], starb, PC_vec[i],  starP))
    }

    cat(rep("-", 72), "\n", sep = "")
    cat(sprintf("  Suggested r:  IC1 -> r=%d  |  IC2 -> r=%d  |  IC3 -> r=%d\n",
                r_IC1, r_IC2, r_IC3))
    cat(sprintf("                PC  -> r=%d  |  IC(BIC) -> r=%d\n", r_PC, r_ICb))
    cat("\n")
    cat("  Note on criteria (all: smaller = better):\n")
    cat("    IC1/IC2/IC3 : Bai & Ng (2002) Prop. 1, applied to IFE residuals per\n")
    cat("                  Bai (2009) Sec. 9.4. Consistent under sequential asymptotics\n")
    cat("                  (min(N,T) -> inf). May overselect in moderate-sized panels\n")
    cat("                  (rule of thumb: reliable when min(N,T) > 60).\n")
    cat("    IC(BIC)     : BIC-style with log(NT) penalty -- more conservative and\n")
    cat("                  generally more reliable for small-to-moderate panels.\n")
    cat("    PC          : Bai (2009) prediction criterion with small-sample correction.\n")
    cat("    Recommendation: use IC(BIC) when IC1/IC2/IC3 decrease monotonically\n")
    cat("    (a sign that the penalty is too weak relative to sample size).\n")
    cat(rep("=", 72), "\n", sep = "")
    cat("\n")
  }

  # ---- return invisible data.frame ----
  out <- data.frame(
    r      = r_vals,
    V_r    = V_r_vec,
    IC1    = IC1_vec,
    IC2    = IC2_vec,
    IC3    = IC3_vec,
    IC_bic = ICb_vec,
    PC     = PC_vec,
    converged = conv_vec,
    stringsAsFactors = FALSE
  )
  attr(out, "suggested") <- c(IC1 = r_IC1, IC2 = r_IC2, IC3 = r_IC3,
                               IC_bic = r_ICb, PC = r_PC)
  invisible(out)
}
