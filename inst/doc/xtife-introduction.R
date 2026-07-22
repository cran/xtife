## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----data---------------------------------------------------------------------
library(xtife)
data(cigar)
dim(cigar)
head(cigar[, c("state", "year", "sales", "price")], 3)

## ----quickstart---------------------------------------------------------------
fit <- ife(sales ~ price,
           data  = cigar,
           index = c("state", "year"),
           r     = 2,
           force = "two-way",
           se    = "cluster")
print(fit)

## ----components---------------------------------------------------------------
fit$coef          # named coefficient vector
fit$se            # standard errors
fit$ci            # 95% confidence intervals (matrix: coef x 2)
fit$converged     # TRUE if outer loop converged
fit$n_iter        # number of outer iterations
dim(fit$F_hat)    # T x r estimated factors
dim(fit$Lambda_hat) # N x r estimated loadings

## ----se_compare---------------------------------------------------------------
fit_std <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "standard")
fit_rob <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "robust")
fit_cl  <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "cluster")

data.frame(
  se_type = c("standard", "robust (HC1)", "cluster"),
  coef    = round(c(fit_std$coef, fit_rob$coef, fit_cl$coef), 5),
  se      = round(c(fit_std$se,   fit_rob$se,   fit_cl$se),   5),
  t_stat  = round(c(fit_std$tstat, fit_rob$tstat, fit_cl$tstat), 3)
)

## ----twfe_compare-------------------------------------------------------------
fit0 <- ife(sales ~ price, data = cigar,
            index = c("state", "year"), r = 0, force = "two-way")

# Manual within-transformation
cigar$y_dm <- cigar$sales - ave(cigar$sales, cigar$state) -
                             ave(cigar$sales, cigar$year) + mean(cigar$sales)
cigar$x_dm <- cigar$price - ave(cigar$price, cigar$state) -
                             ave(cigar$price, cigar$year) + mean(cigar$price)
lm0 <- lm(y_dm ~ x_dm - 1, data = cigar)

cat(sprintf("ife (r=0): %.7f\n", fit0$coef["price"]))
cat(sprintf("lm TWFE:  %.7f\n", coef(lm0)["x_dm"]))
cat(sprintf("diff:     %.2e\n",
            abs(fit0$coef["price"] - coef(lm0)["x_dm"])))

## ----select_r, eval=FALSE-----------------------------------------------------
# # Not run during build (takes ~20 s); set verbose=FALSE in scripts
# sel <- ife_select_r(sales ~ price, data = cigar,
#                     index = c("state", "year"),
#                     r_max = 6, force = "two-way", verbose = FALSE)
# attr(sel, "suggested")   # named vector: r chosen by each IC

## ----bias_static--------------------------------------------------------------
fit_bc <- ife(sales ~ price, data = cigar,
              index = c("state", "year"), r = 2,
              se = "standard", bias_corr = TRUE)
cat(sprintf("Raw IFE:  %.5f\n", fit_bc$coef_raw["price"]))
cat(sprintf("Corrected:%.5f\n", fit_bc$coef["price"]))
cat(sprintf("B_hat:    %.5f  (B_hat/N = %.5f)\n",
            fit_bc$B_hat["price"], fit_bc$B_hat["price"] / fit_bc$N))
cat(sprintf("C_hat:    %.5f  (C_hat/T = %.5f)\n",
            fit_bc$C_hat["price"], fit_bc$C_hat["price"] / fit_bc$T))

## ----bias_dynamic-------------------------------------------------------------
fit_dyn <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2,
               se = "standard", method = "dynamic",
               bias_corr = TRUE, M1 = 1L)
cat(sprintf("Dynamic BC coef: %.5f\n", fit_dyn$coef["price"]))
cat(sprintf("B1/T (Nickell):  %.5f\n",
            fit_dyn$B1_hat["price"] / fit_dyn$T))

## ----unbalanced_setup---------------------------------------------------------
# Create a 10% randomly missing version of the cigar panel
set.seed(42)
cigar_unb <- cigar[sample(nrow(cigar), size = floor(0.9 * nrow(cigar))), ]
cat(sprintf("Observations: %d / %d  (%.0f%% fill)\n",
            nrow(cigar_unb), nrow(cigar), 100 * nrow(cigar_unb) / nrow(cigar)))

## ----unbalanced_fit-----------------------------------------------------------
fit_unb <- ife_unbalanced(sales ~ price,
                           data  = cigar_unb,
                           index = c("state", "year"),
                           r     = 2L,
                           se    = "cluster")
print(fit_unb)

## ----unbalanced_fields--------------------------------------------------------
fit_unb$coef       # estimated beta
fit_unb$se         # standard errors
fit_unb$ci         # 95% confidence intervals
fit_unb$converged  # outer-loop convergence
fit_unb$n_obs      # number of observed (i,t) pairs
fit_unb$N          # number of units
fit_unb$TT         # number of distinct time periods
dim(fit_unb$F_hat)       # TT x r
dim(fit_unb$Lambda_hat)  # N  x r

## ----init_compare-------------------------------------------------------------
# Both converge to the same coefficient for mild unbalancedness
fit_ols <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, init = "ols")
fit_nnr <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, init = "nnr")
cat(sprintf("OLS init: coef = %.5f,  %d outer iterations\n",
            fit_ols$coef["price"], fit_ols$n_iter))
cat(sprintf("NNR init: coef = %.5f,  %d outer iterations\n",
            fit_nnr$coef["price"], fit_nnr$n_iter))

## ----se_unb-------------------------------------------------------------------
fit_std <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, se = "standard")
fit_rob <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, se = "robust")
fit_cl  <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, se = "cluster")
fit_hac <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, se = "hac")

data.frame(
  se_type = c("standard", "robust", "cluster", "hac"),
  coef    = round(c(fit_std$coef, fit_rob$coef, fit_cl$coef, fit_hac$coef), 5),
  se      = round(c(fit_std$se,   fit_rob$se,   fit_cl$se,   fit_hac$se),   5)
)

## ----sel_unb, eval=FALSE------------------------------------------------------
# # Not run during build (NNR cross-validation takes ~10 s)
# sel_unb <- ife_select_r_unb(sales ~ price, data = cigar_unb,
#                               index = c("state", "year"), verbose = FALSE)
# cat(sprintf("SVT selects r_hat = %d\n", sel_unb$r_hat))

## ----bc_strict----------------------------------------------------------------
fit_bc_strict <- ife_unbalanced(sales ~ price, data = cigar_unb,
                                 index = c("state", "year"), r = 2,
                                 se = "standard", bias_corr = TRUE,
                                 exog = "strict")
cat(sprintf("Raw beta:       %.5f\n", fit_bc_strict$coef_raw["price"]))
cat(sprintf("Corrected beta: %.5f\n", fit_bc_strict$coef["price"]))
cat(sprintf("b_hat (b3+...+b6): %.5f\n", fit_bc_strict$b_hat["price"]))

## ----bc_weak, eval=FALSE------------------------------------------------------
# # Not run: exog="weak" adds the b2 dynamic term
# # (requires HAC SE for valid inference in dynamic models)
# fit_bc_weak <- ife_unbalanced(y ~ lag_y + x, data = df_dynamic,
#                                index = c("unit", "time"), r = 2,
#                                se = "hac", bias_corr = TRUE,
#                                exog = "weak")

## ----bal_vs_unbal-------------------------------------------------------------
set.seed(99)
N <- 30; TT <- 22
F0 <- matrix(rnorm(TT * 2), TT, 2); L0 <- matrix(rnorm(N * 2), N, 2)
X  <- matrix(rnorm(N * TT), TT, N)
Y  <- F0 %*% t(L0) + 0.7 * X + matrix(rnorm(N * TT), TT, N)
sim <- data.frame(unit = rep(seq_len(N), each = TT),
                  time = rep(seq_len(TT), times = N),
                  Y = as.vector(Y), X = as.vector(X))

fit_bal <- ife(Y ~ X, data = sim, index = c("unit", "time"),
               r = 2, force = "none", se = "standard")

fit_unb2 <- ife_unbalanced(Y ~ X, data = sim, index = c("unit", "time"),
                            r = 2, force = "none", se = "standard")

cat(sprintf("ife (balanced):     %.6f\n", fit_bal$coef["X"]))
cat(sprintf("ife_unbalanced:     %.6f\n", fit_unb2$coef["X"]))
cat(sprintf("Difference:         %.2e\n",
            abs(fit_bal$coef["X"] - fit_unb2$coef["X"])))

