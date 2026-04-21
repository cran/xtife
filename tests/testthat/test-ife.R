## Tests for xtife package — converted from test_ife.R (12 tests)
## Data: cigar (46 US states x 30 years, 1963-1992)

data(cigar, package = "xtife")

# ==============================================================================
# Fixtures shared across tests
# ==============================================================================

fit2_std <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                r = 2, force = "two-way", se = "standard")

# ==============================================================================
# TEST 1: r=0 matches standard TWFE OLS to machine precision
# ==============================================================================

test_that("r=0 matches TWFE OLS to machine precision", {
  fit0 <- ife(sales ~ price, data = cigar, index = c("state", "year"),
              r = 0, force = "two-way", se = "standard")

  # two-way within transformation (manual demeaning)
  cigar2 <- cigar
  cigar2$y_dm <- cigar2$sales  - ave(cigar2$sales,  cigar2$state) -
                                   ave(cigar2$sales,  cigar2$year)  +
                                   mean(cigar2$sales)
  cigar2$x_dm <- cigar2$price  - ave(cigar2$price,  cigar2$state) -
                                   ave(cigar2$price,  cigar2$year)  +
                                   mean(cigar2$price)
  lm0 <- lm(y_dm ~ x_dm - 1, data = cigar2)

  expect_equal(unname(fit0$coef["price"]),
               unname(coef(lm0)["x_dm"]),
               tolerance = 1e-5)
})

# ==============================================================================
# TEST 2: r=2 coefficient in plausible range
# ==============================================================================

test_that("r=2 coefficient in plausible range (-0.65 to -0.30)", {
  expect_true(fit2_std$coef["price"] > -0.65 & fit2_std$coef["price"] < -0.30)
  expect_true(fit2_std$converged)
})

# ==============================================================================
# TEST 3: HC1 robust SE >= standard SE
# ==============================================================================

test_that("robust (HC1) SE >= standard SE", {
  fit2_rob <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                  r = 2, force = "two-way", se = "robust")
  expect_gte(fit2_rob$se["price"], fit2_std$se["price"])
})

# ==============================================================================
# TEST 4: SE ordering cluster >= robust >= standard
# ==============================================================================

test_that("SE ordering: cluster >= robust >= standard (80% tolerance)", {
  fit2_rob <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                  r = 2, force = "two-way", se = "robust")
  fit2_cl  <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                  r = 2, force = "two-way", se = "cluster")

  expect_gte(fit2_cl$se["price"],  fit2_rob$se["price"] * 0.8)
  expect_gte(fit2_rob$se["price"], fit2_std$se["price"] * 0.8)
})

# ==============================================================================
# TEST 5: Factor normalization F'F/T = I_r
# ==============================================================================

test_that("factor normalization F'F/T = I_r to machine precision", {
  FtF_over_T <- crossprod(fit2_std$F_hat) / fit2_std$T
  off_diag   <- max(abs(FtF_over_T - diag(nrow(FtF_over_T))))
  expect_lt(off_diag, 1e-8)
})

# ==============================================================================
# TEST 6: force options produce convergence
# ==============================================================================

test_that("all force options converge", {
  for (frc in c("none", "unit", "time", "two-way")) {
    fit_f <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                 r = 2, force = frc, se = "standard")
    expect_true(fit_f$converged,
                label = paste0("force='", frc, "' convergence"))
  }
})

# ==============================================================================
# TEST 7: Static bias correction (Bai 2009)
# ==============================================================================

test_that("static bias correction produces corrected coef ~-0.5309", {
  fit_bc <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                r = 2, force = "two-way", se = "standard", bias_corr = TRUE)

  expect_equal(unname(fit_bc$coef["price"]), -0.5309, tolerance = 0.005)
  expect_true(fit_bc$bias_applied)
  expect_false(is.null(fit_bc$B_hat))
  expect_false(is.null(fit_bc$C_hat))
})

# ==============================================================================
# TEST 8: Dynamic method, no bias correction
# ==============================================================================

test_that("dynamic method coef is in plausible range and converges", {
  fit_dyn <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                 r = 2, force = "two-way", se = "standard", method = "dynamic")

  expect_true(fit_dyn$coef["price"] > -0.70 & fit_dyn$coef["price"] < -0.30)
  expect_true(fit_dyn$converged)
  expect_identical(fit_dyn$method, "dynamic")
})

# ==============================================================================
# TEST 9: Dynamic bias correction (Moon & Weidner 2017)
# ==============================================================================

test_that("dynamic bias correction produces corrected coef ~-0.5343", {
  fit_dbc <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                 r = 2, force = "two-way", se = "standard",
                 method = "dynamic", bias_corr = TRUE, M1 = 1L)

  expect_equal(unname(fit_dbc$coef["price"]), -0.5343, tolerance = 0.005)
  expect_true(fit_dbc$bias_applied)
  expect_false(is.null(fit_dbc$B1_hat))
  expect_false(is.null(fit_dbc$B2_hat))
  expect_false(is.null(fit_dbc$B3_hat))
  # B1 near zero (price is approximately exogenous in cigar data)
  expect_lt(abs(fit_dbc$B1_hat["price"] / fit_dbc$T), 0.01)
})

# ==============================================================================
# TEST 10: Dynamic SE ordering
# ==============================================================================

test_that("dynamic SE ordering: cluster >= robust >= standard (80% tolerance)", {
  fit_dyn_cl  <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                     r = 2, force = "two-way", se = "cluster",  method = "dynamic")
  fit_dyn_rob <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                     r = 2, force = "two-way", se = "robust",   method = "dynamic")
  fit_dyn_std <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                     r = 2, force = "two-way", se = "standard", method = "dynamic")

  expect_gte(fit_dyn_cl$se["price"],  fit_dyn_rob$se["price"] * 0.8)
  expect_gte(fit_dyn_rob$se["price"], fit_dyn_std$se["price"] * 0.8)
})

# ==============================================================================
# TEST 11: ife_select_r() returns correct structure
# ==============================================================================

test_that("ife_select_r returns correct structure", {
  skip_on_cran()   # fitting r=0..4, takes ~15 s

  sel <- ife_select_r(sales ~ price, data = cigar,
                      index = c("state", "year"), r_max = 4, verbose = FALSE)

  expect_s3_class(sel, "data.frame")
  expect_equal(nrow(sel), 5L)
  expect_true(all(c("r", "IC1", "IC2", "IC3", "IC_bic", "PC") %in% names(sel)))
  sug <- attr(sel, "suggested")
  expect_true(!is.null(sug))
  expect_true(all(c("IC1", "IC2", "IC3", "IC_bic", "PC") %in% names(sug)))
})

# ==============================================================================
# TEST 12: method="static" explicit == default (regression guard)
# ==============================================================================

test_that("method='static' explicit equals default to machine precision", {
  fit_default  <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                      r = 2, force = "two-way", se = "standard")
  fit_explicit <- ife(sales ~ price, data = cigar, index = c("state", "year"),
                      r = 2, force = "two-way", se = "standard", method = "static")

  expect_equal(unname(fit_default$coef["price"]),
               unname(fit_explicit$coef["price"]),
               tolerance = 1e-10)
})
