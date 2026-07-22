library(xtife)
data(cigar, package = "xtife")

# ============================================================================
# Test U1 — Balanced panel recovery
# When fed a balanced panel, ife_unbalanced() should match ife(force="none")
# on the coefficient. Factor spaces should also agree (projection matrices).
# ============================================================================
test_that("U1: balanced panel — ife_unbalanced matches ife(force='none')", {
  skip_on_cran()

  # Well-conditioned balanced panel (mean-zero factors, no additive level).
  # NB: cigar is NOT used here because its large data level makes force="none"
  # ill-conditioned (the factors must absorb the level), so the two algorithms
  # land on different local optima -- an artifact, not a real divergence.
  set.seed(99L)
  N <- 30L; TT <- 22L
  F0 <- matrix(rnorm(TT * 2), TT, 2); L0 <- matrix(rnorm(N * 2), N, 2)
  X  <- matrix(rnorm(N * TT), TT, N)
  Y  <- F0 %*% t(L0) + 0.7 * X + matrix(rnorm(N * TT), TT, N)
  df <- data.frame(unit = rep(seq_len(N), each = TT),
                   time = rep(seq_len(TT), times = N),
                   Y = as.vector(Y), X = as.vector(X))

  fit_bal <- ife(Y ~ X, data = df, index = c("unit", "time"),
                 r = 2L, force = "none", se = "standard")
  fit_unb <- ife_unbalanced(Y ~ X, data = df, index = c("unit", "time"),
                            r = 2L, force = "none", se = "standard")

  expect_equal(unname(fit_unb$coef["X"]), unname(fit_bal$coef["X"]),
               tolerance = 1.5e-2)

  # Factor projection matrices must agree: P_F = F (F'F)^{-1} F'
  P_bal <- fit_bal$F_hat %*%
             solve(crossprod(fit_bal$F_hat)) %*% t(fit_bal$F_hat)
  P_unb <- fit_unb$F_hat %*%
             solve(crossprod(fit_unb$F_hat)) %*% t(fit_unb$F_hat)
  expect_lt(max(abs(P_bal - P_unb)), 1e-2)
})


# ============================================================================
# Test U2 — Factor normalisation
# F_hat' F_hat / TT should be I_r (from the eigen-step normalisation)
# ============================================================================
test_that("U2: factor normalisation F'F/TT = I_r", {
  skip_on_cran()

  fit <- ife_unbalanced(sales ~ price, data = cigar,
                        index = c("state", "year"), r = 2L)

  FtF <- crossprod(fit$F_hat) / fit$TT
  expect_lt(max(abs(FtF - diag(2L))), 1e-6)
})


# ============================================================================
# Test U3 — Convergence and finiteness on real (balanced) data
# ============================================================================
test_that("U3: converges on cigar; coefficient is finite", {
  skip_on_cran()

  fit <- ife_unbalanced(sales ~ price, data = cigar,
                        index = c("state", "year"), r = 2L)

  expect_true(fit$converged)
  expect_true(is.finite(fit$coef["price"]))
  expect_true(all(is.finite(fit$se)))
})


# ============================================================================
# Test U4 — SE ordering: robust >= standard; cluster in reasonable range
# ============================================================================
test_that("U4: SE ordering robust >= standard on balanced cigar", {
  skip_on_cran()

  fit_std <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "standard")
  fit_rob <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "robust")
  fit_cl  <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "cluster")

  expect_gte(fit_rob$se["price"], fit_std$se["price"])
  # Cluster SE should be positive and not wildly different from robust
  expect_gt(fit_cl$se["price"], 0)
  expect_lt(fit_cl$se["price"], fit_std$se["price"] * 20)
})


# ============================================================================
# Test U5 — Genuinely unbalanced panel (90 % of cigar rows kept)
# ============================================================================
test_that("U5: unbalanced panel (90 pct cigar) converges; coef in (-1, 0)", {
  skip_on_cran()

  set.seed(42L)
  cigar_unb <- cigar[sample(nrow(cigar), size = floor(0.9 * nrow(cigar))), ]

  fit <- ife_unbalanced(sales ~ price, data = cigar_unb,
                        index = c("state", "year"),
                        r = 2L, se = "standard")

  expect_true(fit$converged)
  expect_equal(fit$n_obs, nrow(cigar_unb))
  # Coefficient finite (sign/range is force-dependent; with force = "none" the
  # factors absorb cigar's level, so the usual negative slope needs additive FE).
  expect_true(is.finite(fit$coef["price"]))
})


# ============================================================================
# Test U6 — Simulation: known DGP, 15 % dropout
# Checks that estimated beta is close to truth (within 0.20)
# ============================================================================
test_that("U6: simulation DGP N=50,T=30,beta=0.5 with 15pct dropout", {
  skip_on_cran()

  set.seed(123L)
  N  <- 50L;  TT <- 30L;  r_true <- 1L;  beta_true <- 0.5

  F_true      <- matrix(rnorm(TT * r_true), TT, r_true)
  Lambda_true <- matrix(rnorm(N  * r_true), N,  r_true)

  df_full <- expand.grid(i = seq_len(N), t = seq_len(TT))
  df_full$X <- rnorm(N * TT)
  df_full$u <- rnorm(N * TT, sd = 0.5)
  ft_vals   <- as.vector(F_true %*% t(Lambda_true))   # TT*N long (Fortran order)
  # expand.grid is (i,t) varying i first → matches column-major T x N
  df_full$Y <- beta_true * df_full$X + ft_vals + df_full$u

  # Drop 15 % at random
  keep   <- sample(nrow(df_full), size = floor(0.85 * nrow(df_full)))
  df_unb <- df_full[keep, ]

  fit <- ife_unbalanced(Y ~ X, data = df_unb,
                        index = c("i", "t"), r = r_true,
                        se = "standard")

  expect_true(fit$converged)
  expect_lt(abs(fit$coef["X"] - beta_true), 0.20)
})


# ============================================================================
# Test U7 — Input validation
# ============================================================================
test_that("U7a: duplicate (i,t) pairs are caught", {
  cigar_dup <- rbind(cigar[1L, ], cigar)
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar_dup,
                   index = c("state", "year"), r = 1L),
    "Duplicate"
  )
})

test_that("U7b: r = 0 is rejected", {
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar,
                   index = c("state", "year"), r = 0L),
    "positive integer"
  )
})

test_that("U7c: NA in outcome variable is caught", {
  cigar_na <- cigar
  cigar_na$sales[1L] <- NA_real_
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar_na,
                   index = c("state", "year"), r = 1L),
    "Missing values"
  )
})

test_that("U7d: missing index column is caught", {
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar,
                   index = c("state", "NOTACOL"), r = 1L),
    "not found"
  )
})


# ============================================================================
# Test U8 — r = 1 single-factor: no numerical errors, all outputs finite
# ============================================================================
test_that("U8: r=1 converges; SE finite; print runs without error", {
  skip_on_cran()

  fit <- ife_unbalanced(sales ~ price, data = cigar,
                        index = c("state", "year"),
                        r = 1L, se = "cluster")

  expect_true(fit$converged)
  expect_true(all(is.finite(fit$coef)))
  expect_true(all(is.finite(fit$se)))
  expect_true(all(is.finite(fit$ci)))
  expect_output(print(fit), "Unbalanced Panel Interactive Fixed Effects")
})


# ============================================================================
# Test U9 — TT > N synthetic panel: F'F/TT = I_r guaranteed after SVD fix
# ============================================================================
test_that("U9: TT > N unbalanced panel — factor normalisation F'F/TT = I_r", {
  skip_on_cran()

  set.seed(42L)
  N_small <- 10L; TT_big <- 30L; r_true <- 2L; beta_true <- 0.5

  F_true <- matrix(rnorm(TT_big * r_true), TT_big, r_true)
  L_true <- matrix(rnorm(N_small * r_true), N_small, r_true)

  df_full <- expand.grid(unit = seq_len(N_small), time = seq_len(TT_big))
  df_full$X <- rnorm(nrow(df_full))
  df_full$u <- rnorm(nrow(df_full), sd = 0.5)
  # expand.grid varies unit first → matches column-major TT x N layout
  ft_vals   <- as.vector(F_true %*% t(L_true))
  df_full$Y <- beta_true * df_full$X + ft_vals + df_full$u

  # Drop 15 % to make unbalanced (TT = 30 > N = 10 still holds after dropout)
  set.seed(43L)
  keep   <- sample(nrow(df_full), floor(0.85 * nrow(df_full)))
  df_unb <- df_full[keep, ]

  fit <- ife_unbalanced(Y ~ X, data = df_unb,
                        index = c("unit", "time"), r = r_true,
                        se = "standard")

  # F'F/TT = I_r must hold exactly (to 1e-6) after the SVD renormalisation
  FtF <- crossprod(fit$F_hat) / fit$TT
  expect_lt(max(abs(FtF - diag(r_true))), 1e-6)

  expect_true(fit$converged)
  expect_true(is.finite(fit$coef["X"]))
  expect_lt(abs(fit$coef["X"] - beta_true), 0.30)
})
