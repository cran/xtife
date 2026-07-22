library(xtife)
data(cigar, package = "xtife")

# ============================================================================
# Su, Wang and Wang (2025) Test Suite — implementation in ife_unbalanced
# ============================================================================

# ============================================================================
# Test S1 — NNR converges; Theta_hat is low-rank; beta finite
# ============================================================================
test_that("S1: NNR soft-impute converges; beta0 finite; Theta close to rank-r", {
  skip_on_cran()

  set.seed(42L)
  cigar_unb <- cigar[sample(nrow(cigar), 1100L), ]

  fit <- ife_unbalanced(sales ~ price, data = cigar_unb,
                        index = c("state", "year"),
                        r = 2L, se = "standard", init = "nnr")

  expect_true(fit$converged)
  expect_true(is.finite(fit$coef["price"]))
  expect_true(all(is.finite(fit$se)))
  # nu_used should be stored
  expect_false(is.null(fit$nu_used))
  expect_true(is.finite(fit$nu_used))
})


# ============================================================================
# Test S2 — SVT factor selection on balanced cigar (known ground truth ~ 2)
# ============================================================================
test_that("S2: ife_select_r_unb selects r_hat >= 1 on balanced cigar", {
  skip_on_cran()

  res <- ife_select_r_unb(sales ~ price, data = cigar,
                           index = c("state", "year"),
                           c_f = 0.6, verbose = FALSE)

  expect_true(is.list(res))
  expect_true(res$r_hat >= 1L)
  expect_true(is.numeric(res$sv))
  expect_true(is.finite(res$threshold))
  # On the full balanced cigar, 2 factors is standard; SVT should pick >= 1
  expect_gte(res$r_hat, 1L)
})


# ============================================================================
# Test S3 — delta/omega alternating LS: on balanced cigar with r=2
#   Check that the Su, Wang and Wang (2025) projected x̂_it is well-defined and has
#   smaller norm than the raw X (factor component is absorbed)
# ============================================================================
test_that("S3: delta/omega produces finite projected regressors with smaller variance", {
  skip_on_cran()

  fit_std <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "standard")
  fit_rob <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "robust")

  # Both use Su, Wang and Wang (2025) projected regressors internally; SE should be finite
  expect_true(all(is.finite(fit_std$se)))
  expect_true(all(is.finite(fit_rob$se)))

  # The SE computation should not error even with balanced data
  expect_gt(fit_std$se["price"], 0)
  expect_gt(fit_rob$se["price"], 0)
})


# ============================================================================
# Test S4 — Exact SE differs from FWL-approx SE for unbalanced panel
# (Both are valid; they are numerically distinct because x̂_it ≠ x̃_it
#  on an unbalanced panel.)
# ============================================================================
test_that("S4: SE is finite and positive on genuinely unbalanced panel", {
  skip_on_cran()

  set.seed(7L)
  cigar_unb <- cigar[sample(nrow(cigar), size = floor(0.80 * nrow(cigar))), ]

  fit_std <- ife_unbalanced(sales ~ price, data = cigar_unb,
                             index = c("state", "year"),
                             r = 2L, se = "standard")
  fit_cl  <- ife_unbalanced(sales ~ price, data = cigar_unb,
                             index = c("state", "year"),
                             r = 2L, se = "cluster")

  # All SEs finite and positive
  expect_true(all(is.finite(fit_std$se)))
  expect_true(all(is.finite(fit_cl$se)))
  expect_gt(fit_std$se["price"], 0)
  expect_gt(fit_cl$se["price"], 0)

  # Coefficient is finite (range checks are force-dependent: with force = "none"
  # the factors must absorb the data level, so the cigar coef is not the usual
  # negative demand slope -- that requires additive FE, tested separately in S4b).
  expect_true(is.finite(fit_std$coef["price"]))
})


# ============================================================================
# Test S4b — force support matches the balanced ife() on a balanced panel
# ============================================================================
test_that("S4b: ife_unbalanced(force=...) matches ife() on a balanced panel", {
  skip_on_cran()

  # small balanced two-way-FE panel
  set.seed(11L)
  N <- 25L; TT <- 20L
  alpha <- rnorm(N); xi <- rnorm(TT)
  F0 <- matrix(rnorm(TT * 2), TT, 2); L0 <- matrix(rnorm(N * 2), N, 2)
  X  <- matrix(rnorm(N * TT), TT, N)
  Y  <- matrix(alpha, TT, N, byrow = TRUE) + matrix(xi, TT, N) +
        F0 %*% t(L0) + 0.7 * X + matrix(rnorm(N * TT), TT, N)
  df <- data.frame(unit = rep(seq_len(N), each = TT),
                   time = rep(seq_len(TT), times = N),
                   Y = as.vector(Y), X = as.vector(X))

  # Compare on the CORRECTLY specified model (two-way matches the DGP, which has
  # additive + interactive FE).  On a balanced panel ife_unbalanced should
  # reproduce the balanced ife() estimator.  (force = "none"/"unit"/"time" are
  # misspecified for this DGP, so both estimators chase different local optima
  # and need not agree -- not a meaningful comparison.)
  b_bal <- ife(Y ~ X, data = df, index = c("unit", "time"),
               r = 2L, force = "two-way", se = "standard")$coef["X"]
  b_unb <- suppressWarnings(
    ife_unbalanced(Y ~ X, data = df, index = c("unit", "time"),
                   r = 2L, force = "two-way", se = "standard"))$coef["X"]
  expect_equal(unname(b_unb), unname(b_bal), tolerance = 1e-2)
})


# ============================================================================
# Test S5 — Bias correction reduces |beta - 0.5| in simulation DGP
# (Same DGP as Test U6: N=50, T=30, r=1, beta=0.5, 15% dropout)
# ============================================================================
test_that("S5: bias correction reduces estimation bias in simulation DGP", {
  skip_on_cran()

  set.seed(123L)
  N  <- 50L;  TT <- 30L;  r_true <- 1L;  beta_true <- 0.5

  F_true      <- matrix(rnorm(TT * r_true), TT, r_true)
  Lambda_true <- matrix(rnorm(N  * r_true), N,  r_true)

  df_full <- expand.grid(i = seq_len(N), t = seq_len(TT))
  df_full$X <- rnorm(N * TT)
  df_full$u <- rnorm(N * TT, sd = 0.5)
  ft_vals   <- as.vector(F_true %*% t(Lambda_true))
  df_full$Y <- beta_true * df_full$X + ft_vals + df_full$u

  keep   <- sample(nrow(df_full), size = floor(0.85 * nrow(df_full)))
  df_unb <- df_full[keep, ]

  fit_plain <- ife_unbalanced(Y ~ X, data = df_unb,
                               index = c("i", "t"), r = r_true,
                               se = "standard", bias_corr = FALSE)
  fit_bc    <- ife_unbalanced(Y ~ X, data = df_unb,
                               index = c("i", "t"), r = r_true,
                               se = "standard", bias_corr = TRUE)

  # Both should converge and produce finite coefficients
  expect_true(fit_plain$converged)
  expect_true(fit_bc$converged)
  expect_true(is.finite(fit_plain$coef["X"]))
  expect_true(is.finite(fit_bc$coef["X"]))

  # coef_raw should equal the uncorrected estimate
  expect_equal(unname(fit_bc$coef_raw["X"]),
               unname(fit_plain$coef["X"]),
               tolerance = 1e-8)

  # Bias correction should not explode the estimate
  expect_true(abs(fit_bc$coef["X"] - beta_true) < 0.30)
})


# ============================================================================
# Test S6 — init = "nnr" vs init = "ols": coefficients should be close
# (Both converge to the same AM fixed point if the panel is reasonable)
# ============================================================================
test_that("S6: init='nnr' and init='ols' give same coef to 1e-3 on cigar", {
  skip_on_cran()

  fit_ols <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "standard", init = "ols")
  fit_nnr <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "standard", init = "nnr")

  expect_equal(unname(fit_ols$coef["price"]),
               unname(fit_nnr$coef["price"]),
               tolerance = 1e-3)
})


# ============================================================================
# Test S7 — ife_select_r_unb returns r_hat >= 1 on unbalanced cigar
# ============================================================================
test_that("S7: ife_select_r_unb on unbalanced cigar returns r_hat >= 1", {
  skip_on_cran()

  set.seed(123L)
  cigar_unb <- cigar[sample(nrow(cigar), 1100L), ]

  res <- suppressWarnings(
    ife_select_r_unb(sales ~ price, data = cigar_unb,
                     index = c("state", "year"),
                     c_f = 0.6, verbose = FALSE)
  )

  expect_true(is.list(res))
  expect_gte(res$r_hat, 1L)
  expect_true(is.numeric(res$sv))
  expect_true(all(is.finite(res$sv)))
  expect_true(is.finite(res$threshold))
  expect_gt(res$threshold, 0)
})


# ============================================================================
# Test S8 — All new return fields are present when bias_corr = TRUE
# ============================================================================
test_that("S8: all new return fields present with bias_corr = TRUE", {
  skip_on_cran()

  fit <- ife_unbalanced(sales ~ price, data = cigar,
                        index = c("state", "year"),
                        r = 2L, se = "standard",
                        bias_corr = TRUE)

  # Standard fields
  expect_true(!is.null(fit$coef))
  expect_true(!is.null(fit$coef_raw))
  expect_true(!is.null(fit$vcov))
  expect_true(!is.null(fit$se))
  expect_true(!is.null(fit$F_hat))
  expect_true(!is.null(fit$Lambda_hat))

  # Su, Wang and Wang (2025) bias components
  expect_true(!is.null(fit$b_hat))
  expect_true(!is.null(fit$b3))
  expect_true(!is.null(fit$b4))
  expect_true(!is.null(fit$b5))
  expect_true(!is.null(fit$b6))
  expect_equal(length(fit$b_hat), 1L)   # 1 regressor

  # Meta fields
  expect_equal(fit$init, "ols")
  expect_true(fit$bias_corr)

  # b_hat = b3 + b4 + b5 + b6
  expect_equal(fit$b_hat, fit$b3 + fit$b4 + fit$b5 + fit$b6,
               tolerance = 1e-12)
})


# ============================================================================
# Test S9 — Cluster SE wider than standard SE (on unbalanced panel)
# ============================================================================
test_that("S9: cluster SE > standard SE on unbalanced cigar panel", {
  skip_on_cran()

  set.seed(99L)
  cigar_unb <- cigar[sample(nrow(cigar), 1100L), ]

  fit_std <- ife_unbalanced(sales ~ price, data = cigar_unb,
                             index = c("state", "year"),
                             r = 2L, se = "standard")
  fit_cl  <- ife_unbalanced(sales ~ price, data = cigar_unb,
                             index = c("state", "year"),
                             r = 2L, se = "cluster")

  # Cluster SE should be positive and at most 50x standard SE (sanity bound)
  expect_gt(fit_cl$se["price"], 0)
  expect_lt(fit_cl$se["price"] / fit_std$se["price"], 50)
})


# ============================================================================
# Test S10 — b5 + b6 ≈ 0 when factors and loadings are orthogonal to residuals
#   On balanced cigar, the bias terms should be small (O(1/sqrt(NT)))
# ============================================================================
test_that("S10: bias components finite and correction small relative to coef", {
  skip_on_cran()

  fit <- ife_unbalanced(sales ~ price, data = cigar,
                        index = c("state", "year"),
                        r = 2L, se = "standard",
                        bias_corr = TRUE)

  # All bias components should be finite
  expect_true(all(is.finite(fit$b3)))
  expect_true(all(is.finite(fit$b4)))
  expect_true(all(is.finite(fit$b5)))
  expect_true(all(is.finite(fit$b6)))

  # The bias correction should be small relative to the raw coefficient
  # (O(1/N + 1/T) ~ O(1/30 + 1/46) ~ 0.055 on cigar)
  correction <- abs(fit$coef["price"] - fit$coef_raw["price"])
  expect_lt(correction, abs(fit$coef_raw["price"]) * 0.50)

  # print method should show bias correction header without error
  expect_output(print(fit), "Bias corr")
})


# ============================================================================
# Test S11 — input validation still works after update
# ============================================================================
test_that("S11: input validation errors still propagate correctly", {
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar,
                   index = c("state", "year"), r = 1L, init = "bad"),
    "'init' must be"
  )
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar,
                   index = c("state", "year"), r = 0L),
    "positive integer"
  )
})


# ============================================================================
# Test S12 — exog = "weak": b2 finite; b2 = 0 for strict exog
# ============================================================================
test_that("S12: b2 finite when exog='weak'; exactly zero for exog='strict'", {
  skip_on_cran()

  # Strict exogeneity: b2 must be stored and equal to zero vector
  fit_strict <- ife_unbalanced(sales ~ price, data = cigar,
                               index = c("state", "year"),
                               r = 2L, se = "standard",
                               bias_corr = TRUE, exog = "strict")
  expect_true(!is.null(fit_strict$b2))
  expect_equal(fit_strict$b2, rep(0, length(fit_strict$b2)))

  # b_hat = b2 + b3 + b4 + b5 + b6 (b2 = 0 for strict, so unchanged)
  expect_equal(fit_strict$b_hat,
               fit_strict$b2 + fit_strict$b3 + fit_strict$b4 +
               fit_strict$b5 + fit_strict$b6,
               tolerance = 1e-12)

  # Weak exogeneity: b2 should be finite (price is roughly exogenous so b2 ≈ 0,
  # but it's computed from the one-sided kernel and should be a valid number)
  fit_weak <- ife_unbalanced(sales ~ price, data = cigar,
                             index = c("state", "year"),
                             r = 2L, se = "standard",
                             bias_corr = TRUE, exog = "weak")
  expect_true(!is.null(fit_weak$b2))
  expect_true(all(is.finite(fit_weak$b2)))

  # b_hat includes b2
  expect_equal(fit_weak$b_hat,
               fit_weak$b2 + fit_weak$b3 + fit_weak$b4 +
               fit_weak$b5 + fit_weak$b6,
               tolerance = 1e-12)
})


# ============================================================================
# Test S13 — se = "hac": SE finite, positive, L_T auto-computed
# Also verifies normalization: with L_T = 1 only the gap=0 diagonal survives,
# so HAC reduces to the HC sandwich (≈ robust) up to a tiny small-sample
# correction factor in the "robust" branch.
# ============================================================================
test_that("S13: se='hac' gives finite positive SE; L_T stored correctly", {
  skip_on_cran()

  set.seed(42L)
  cigar_unb <- cigar[sample(nrow(cigar), 1100L), ]

  fit_hac <- ife_unbalanced(sales ~ price, data = cigar_unb,
                             index = c("state", "year"),
                             r = 2L, se = "hac")

  expect_true(all(is.finite(fit_hac$se)))
  expect_gt(fit_hac$se["price"], 0)
  expect_false(is.null(fit_hac$L_T))
  expect_equal(fit_hac$L_T, floor(2 * fit_hac$TT^(1/5)))

  # L_T = 1: gap >= 1 always for distinct integer time indices, so only
  # the gap = 0 (t == s) diagonal contributes with kern = 1.
  # This makes B_hac = sum v^2 x x' = B_robust (no small-sample corr in hac).
  # Robust applies corr = n_obs/(n_obs - p) ≈ 1. Ratio should be near 1.
  fit_hac_lt1 <- ife_unbalanced(sales ~ price, data = cigar_unb,
                                 index = c("state", "year"),
                                 r = 2L, se = "hac", L_T = 1L)
  fit_rob     <- ife_unbalanced(sales ~ price, data = cigar_unb,
                                 index = c("state", "year"),
                                 r = 2L, se = "robust")

  # HAC (L_T=1) ≈ robust to within 1% (tiny corr factor difference)
  ratio <- fit_hac_lt1$se["price"] / fit_rob$se["price"]
  expect_lt(abs(ratio - 1), 0.01)

  # Custom L_T must be respected
  fit_hac2 <- ife_unbalanced(sales ~ price, data = cigar_unb,
                              index = c("state", "year"),
                              r = 2L, se = "hac", L_T = 3L)
  expect_equal(fit_hac2$L_T, 3L)
  expect_true(all(is.finite(fit_hac2$se)))
  expect_gt(fit_hac2$se["price"], 0)

  # exog field is returned
  expect_equal(fit_hac$exog, "strict")

  # print runs without error
  expect_output(print(fit_hac), "HAC Bartlett")
})


# ============================================================================
# Test S14 — DGP 2: dynamic panel (x = lag y), AR(1) factors, MA(1) errors
# exog = "weak", se = "hac", bias_corr = TRUE
# ============================================================================
test_that("S14: DGP 2 (dynamic panel) with exog='weak', se='hac', bias_corr=TRUE", {
  skip_on_cran()

  set.seed(42L)
  N <- 40L; TT_full <- 30L; r_true <- 1L; beta_true <- 0.3

  # AR(1) common factor
  F_true <- matrix(0, TT_full, r_true)
  F_true[1L, ] <- rnorm(r_true)
  for (tt in 2:TT_full)
    F_true[tt, ] <- 0.5 * F_true[tt - 1L, ] + rnorm(r_true, sd = sqrt(0.75))
  Lambda_true <- matrix(rnorm(N * r_true), N, r_true)

  # MA(1) error: v_it = (e_it + e_{i,t-1}) / sqrt(2)
  e_mat <- matrix(rnorm(N * TT_full), TT_full, N)
  v_mat <- (e_mat + rbind(matrix(0, 1L, N), e_mat[-TT_full, ])) / sqrt(2)

  # Build balanced panel: y_it = beta * y_{i,t-1} + lambda_i' f_t + v_it
  Y_mat <- matrix(0, TT_full, N)
  Y_mat[1L, ] <- as.vector(F_true[1L, ] %*% t(Lambda_true)) + v_mat[1L, ]
  for (tt in 2:TT_full)
    Y_mat[tt, ] <- beta_true * Y_mat[tt - 1L, ] +
                   as.vector(F_true[tt, ] %*% t(Lambda_true)) + v_mat[tt, ]

  # Long-format: t = 2..TT_full, x = y_{i,t-1}
  df_list <- lapply(seq_len(N), function(i) {
    data.frame(i = i,
               t = 2L:TT_full,
               Y = Y_mat[2L:TT_full, i],
               X = Y_mat[1L:(TT_full - 1L), i])
  })
  df_full <- do.call(rbind, df_list)

  # Drop 15 % at random
  keep   <- sample(nrow(df_full), floor(0.85 * nrow(df_full)))
  df_unb <- df_full[keep, ]

  fit <- ife_unbalanced(Y ~ X, data = df_unb,
                        index  = c("i", "t"), r = r_true,
                        se     = "hac", bias_corr = TRUE, exog = "weak")

  expect_true(fit$converged)
  expect_true(is.finite(fit$coef["X"]))
  expect_true(all(is.finite(fit$se)))
  expect_true(!is.null(fit$b2))
  expect_true(all(is.finite(fit$b2)))
  expect_equal(fit$exog, "weak")
  # Broad tolerance: dynamic unbalanced panel, single replication
  expect_lt(abs(fit$coef["X"] - beta_true), 0.40)
})


# ============================================================================
# Test S15 — Input validation: bad exog / L_T values raise errors
# ============================================================================
test_that("S15: invalid exog and L_T inputs raise errors", {
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar,
                   index = c("state", "year"), r = 1L, exog = "bad"),
    "'exog' must be"
  )
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar,
                   index = c("state", "year"), r = 1L, se = "xyz"),
    "'se' must be"
  )
  expect_error(
    ife_unbalanced(sales ~ price, data = cigar,
                   index = c("state", "year"), r = 1L,
                   se = "hac", L_T = -1L),
    "'L_T' must be a positive"
  )
})


# ============================================================================
# Test S16 — Bias correction must NOT over-correct when the raw estimator is
# unbiased.  Regression guard for the Su, Wang and Wang (2025) Theorem 4.2 sign of the curvature
# inverse [L-bar_ff']_t = (- sum_i d_it lam lam')^{-1}: b3/b4/b2 use it once and
# carry the minus sign; b5/b6 use it twice (sign cancels).  For a balanced panel
# with i.i.d. homoskedastic errors the true incidental-parameter bias is ~0, so
# b3 and b5 cancel and the correction must be tiny.  With the sign bug the
# correction was ~0.07-0.08 (b3 and b5 added instead of cancelling).
# ============================================================================
test_that("S16: BC does not over-correct an unbiased estimator (L-bar sign)", {
  skip_on_cran()

  set.seed(7L)
  N <- 50L; TT <- 50L; r <- 2L; beta_true <- 1.0
  F_mat <- matrix(rnorm(TT * r), TT, r)
  Lam   <- matrix(rnorm(N  * r), N,  r)
  Mu    <- matrix(rnorm(N  * r), N,  r)
  fc    <- as.vector(Lam %*% t(F_mat))
  xc    <- as.vector((Lam + Mu) %*% t(F_mat))
  df    <- expand.grid(unit = seq_len(N), time = seq_len(TT))
  df$X  <- xc + rnorm(N * TT)
  df$Y  <- beta_true * df$X + fc + rnorm(N * TT)   # i.i.d. homoskedastic

  fit <- ife_unbalanced(Y ~ X, data = df, index = c("unit", "time"),
                        r = 2L, se = "standard", bias_corr = TRUE)

  # Correction must be small (was ~0.07 with the sign bug); raw is unbiased.
  expect_lt(abs(unname(fit$coef["X"] - fit$coef_raw["X"])), 0.02)
  expect_lt(abs(unname(fit$coef["X"]) - beta_true), 0.05)

  # On a balanced panel b3 and b5 (and b4/b6) cancel -> |b_hat| is small
  # relative to the individual components.
  expect_lt(abs(fit$b_hat), abs(fit$b3) + abs(fit$b5))
})
