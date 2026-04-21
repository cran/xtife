## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----quickstart---------------------------------------------------------------
library(xtife)
data(cigar)

# Fit IFE with r=2 factors, two-way FE, cluster-robust SE
fit <- ife(sales ~ price, data = cigar,
           index  = c("state", "year"),
           r      = 2,
           force  = "two-way",
           se     = "cluster")
print(fit)

## ----components---------------------------------------------------------------
# Estimated coefficients
fit$coef

# Standard errors
fit$se

# p-values
fit$pval

# 95% confidence intervals
fit$ci

# Estimated factors (T x r matrix)
dim(fit$F_hat)

# Estimated loadings (N x r matrix)
dim(fit$Lambda_hat)

## ----se_types-----------------------------------------------------------------
fit_std <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "standard")
fit_rob <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "robust")
fit_cl  <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "cluster")

# Compare standard errors
se_table <- data.frame(
  se_type  = c("standard", "robust (HC1)", "cluster"),
  coef     = c(fit_std$coef, fit_rob$coef, fit_cl$coef),
  se       = c(fit_std$se,   fit_rob$se,   fit_cl$se),
  t_stat   = c(fit_std$tstat, fit_rob$tstat, fit_cl$tstat)
)
print(se_table, digits = 4, row.names = FALSE)

## ----select_r, eval=FALSE-----------------------------------------------------
# # Not run during package build (takes ~30 s on cigar data)
# sel <- ife_select_r(sales ~ price, data = cigar,
#                     index = c("state", "year"),
#                     r_max = 6,
#                     force = "two-way")

## ----bias_static--------------------------------------------------------------
fit_bc <- ife(sales ~ price, data = cigar,
              index     = c("state", "year"),
              r         = 2,
              se        = "standard",
              bias_corr = TRUE)
print(fit_bc)

## ----bias_dynamic-------------------------------------------------------------
fit_dyn <- ife(sales ~ price, data = cigar,
               index     = c("state", "year"),
               r         = 2,
               se        = "standard",
               method    = "dynamic",
               bias_corr = TRUE,
               M1        = 1L)
print(fit_dyn)

## ----twfe_compare-------------------------------------------------------------
fit0 <- ife(sales ~ price, data = cigar,
            index = c("state", "year"), r = 0)

# Manual two-way demeaning
cigar$y_dm <- cigar$sales  - ave(cigar$sales,  cigar$state) -
                              ave(cigar$sales,  cigar$year)  + mean(cigar$sales)
cigar$x_dm <- cigar$price  - ave(cigar$price,  cigar$state) -
                              ave(cigar$price,  cigar$year)  + mean(cigar$price)
lm0 <- lm(y_dm ~ x_dm - 1, data = cigar)

cat(sprintf("ife (r=0): %.6f\n", fit0$coef["price"]))
cat(sprintf("lm TWFE:  %.6f\n", coef(lm0)["x_dm"]))
cat(sprintf("diff:     %.2e\n",
            abs(fit0$coef["price"] - coef(lm0)["x_dm"])))

