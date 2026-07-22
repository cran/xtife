# xtife: Interactive Fixed Effects Estimator for Panel Data

<!-- badges: start -->
[![R CMD Check](https://github.com/Rickchen0910/xtife/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Rickchen0910/xtife/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/xtife)](https://CRAN.R-project.org/package=xtife)
[![License: GPL v2/v3](https://img.shields.io/badge/License-GPL%20v2%2Fv3-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
<!-- badges: end -->

`xtife` provides a pure base-R implementation of the **Interactive Fixed Effects (IFE)** panel estimator for both **balanced** and **unbalanced** panels. It delivers full analytical standard errors, asymptotic bias correction, and factor number selection with no external dependencies beyond base R.

For a comprehensive review of interactive fixed effects, see [Ditzen & Karavias (2025)](https://doi.org/10.48550/arXiv.2507.19099).

---

## The Model

Standard two-way fixed effects (TWFE) assumes unobserved heterogeneity enters additively. IFE generalises this by allowing unobserved confounders to interact across units and time:

$$y_{it} = \alpha_i + \xi_t + X_{it}'\beta + \lambda_i'F_t + u_{it}$$

where $F_t \in \mathbb{R}^r$ are common factors and $\lambda_i \in \mathbb{R}^r$ are unit-specific loadings. Setting $r = 0$ reduces the model to standard TWFE. For **unbalanced panels**, `ife_unbalanced()` supports the same additive fixed effects via its `force` argument (`"none"`, `"unit"`, `"time"`, `"two-way"`), estimated jointly with the factors by EM on the imputed panel — robust to informative (factor-correlated) missingness. Its default is `force = "none"` (the intercept-free interactive model); use `force = "two-way"` for data with level or trend structure.

---

## Features

| Feature | Balanced (`ife`) | Unbalanced (`ife_unbalanced`) |
|---------|:---:|:---:|
| Estimator | [Bai (2009)](https://doi.org/10.3982/ECTA6135) SVD alternating projections | [Su, Wang & Wang (2025)](https://doi.org/10.2139/ssrn.5177283); EM with matrix completion ([Bai & Ng 2021](https://doi.org/10.1080/01621459.2021.1967163)); OLS or NNR init |
| Standard errors | Homoskedastic · HC1 robust · Cluster | Homoskedastic · HC1 robust · Cluster · HAC |
| Static bias correction | [Bai (2009)](https://doi.org/10.3982/ECTA6135) $\hat B/N + \hat C/T$ | Analytical incidental-parameter correction |
| Dynamic bias correction | [Moon & Weidner (2017)](https://doi.org/10.1017/S0266466615000328) | Analytical, incl. predetermined-regressor term |
| Factor number selection | IC1/2/3 · IC(BIC) · PC | SVT rule (singular value thresholding) |
| Initialisation | — | OLS · Nuclear-norm regularisation (NNR) |
| Dependencies | Base R only | Base R only |

---

## Installation

```r
# From CRAN
install.packages("xtife")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("Rickchen0910/xtife")
```

---

## Balanced Panel: Quick Start

```r
library(xtife)
data(cigar)   # 46 US states x 30 years cigarette panel (Baltagi 1995)

# Fit IFE with r = 2 factors, two-way FE, cluster-robust SE
fit <- ife(sales ~ price, data = cigar,
           index  = c("state", "year"),
           r      = 2,
           force  = "two-way",
           se     = "cluster")
print(fit)
```

```
Interactive Fixed Effects (Bai 2009, Econometrica)
-------------------------------------------------------
        Estimate Std.Error  t.value   Pr(>|t|) CI.lower CI.upper
price    -0.5242    0.0802  -6.5360     0.0000  -0.6814  -0.3670

Converged: TRUE  (10 iterations)
N = 46  T = 30  r = 2  force = two-way  se = cluster
```

### Standard error types

```r
fit_std <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "standard")
fit_rob <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "robust")
fit_cl  <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2, se = "cluster")
```

| `se =` | Assumption | Typical use |
|--------|-----------|-------------|
| `"standard"` | Homoskedasticity | Benchmark |
| `"robust"` | HC1 sandwich | Heteroskedasticity across cells |
| `"cluster"` | Cluster-robust by unit | Serial correlation within units (recommended) |

### Factor number selection

```r
sel <- ife_select_r(sales ~ price, data = cigar,
                    index = c("state", "year"),
                    r_max = 6, force = "two-way")
```

Prints a table of IC1, IC2, IC3 ([Bai & Ng 2002](https://doi.org/10.1111/1468-0262.00273)), IC(BIC), and PC ([Bai 2009](https://doi.org/10.3982/ECTA6135)) for each candidate $r$. **IC(BIC)** is recommended for panels with $\min(N, T) < 60$.

### Asymptotic bias correction

```r
# Static bias correction — Bai (2009)
fit_bc <- ife(sales ~ price, data = cigar,
              index = c("state", "year"), r = 2, bias_corr = TRUE)

# Dynamic bias correction — Moon & Weidner (2017)
# Use when regressors include lagged dependent variables
fit_dyn <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2,
               method = "dynamic", bias_corr = TRUE, M1 = 1L)
```

For the cigar panel ($N = 46$, $T = 30$, $T/N \approx 0.65$):

| Estimator | Price coefficient |
|-----------|------------------|
| TWFE ($r = 0$) | −0.3796 |
| IFE ($r = 2$) | −0.5242 |
| IFE + [Bai (2009)](https://doi.org/10.3982/ECTA6135) bias correction | −0.5309 |
| IFE dynamic + [Moon & Weidner (2017)](https://doi.org/10.1017/S0266466615000328) bias correction | −0.5343 |

### Comparison with TWFE

Setting `r = 0` recovers the standard two-way FE estimator, identical to `lm()` with unit and time dummies at machine precision:

```r
fit0 <- ife(sales ~ price, data = cigar,
            index = c("state", "year"), r = 0)
# Equivalent to plm(..., model = "within", effect = "twoways")
```

---

## Unbalanced Panel: Quick Start

`ife_unbalanced()` fits the IFE model on genuinely unbalanced panels following the estimation and inference theory of [Su, Wang & Wang (2025)](https://doi.org/10.2139/ssrn.5177283), which extends the interactive fixed effects estimator of [Bai (2009)](https://doi.org/10.3982/ECTA6135) to the unbalanced case using the missing-data factor analysis / matrix completion of [Bai & Ng (2021)](https://doi.org/10.1080/01621459.2021.1967163). The core algorithm is an alternating outer loop that updates β and the structure $(\hat\alpha, \hat\xi, \hat\lambda, \hat F)$ until convergence, with an expectation-maximisation inner loop that imputes the unobserved cells from the current structure. Initialisation is via plain OLS (`init = "ols"`, the default) or nuclear-norm-regularised soft-impute ([Mazumder, Hastie & Tibshirani 2010](https://www.jmlr.org/papers/v11/mazumder10a.html); `init = "nnr"`, recommended when the panel is severely unbalanced or $r \geq 3$). The resulting estimator is $\sqrt{NT}$-consistent and asymptotically normal.

```r
# Simulate a 10% randomly missing panel
set.seed(42)
cigar_unb <- cigar[sample(nrow(cigar), size = floor(0.9 * nrow(cigar))), ]

fit_unb <- ife_unbalanced(sales ~ price,
                           data  = cigar_unb,
                           index = c("state", "year"),
                           r     = 2L,
                           se    = "cluster")
print(fit_unb)
```

### Standard error types

| `se =` | Assumption | When to use |
|--------|-----------|-------------|
| `"standard"` | Homoskedastic i.i.d. errors | Benchmark |
| `"robust"` | HC1 | Cell-level heteroskedasticity |
| `"cluster"` | Cluster-robust by unit | Serial correlation within units |
| `"hac"` | Bartlett kernel (bandwidth $L_T = \lfloor 2T^{1/5}\rfloor$) | Serial correlation over time within units |

### Factor number selection

`ife_select_r_unb()` applies the singular value thresholding (SVT) rule of [Su, Wang & Wang (2025)](https://doi.org/10.2139/ssrn.5177283) to the nuclear-norm regularised matrix $\hat\Theta^{(0)}$ — a missing-data counterpart of the [Bai & Ng (2002)](https://doi.org/10.1111/1468-0262.00273) information criteria:

```r
sel_unb <- ife_select_r_unb(sales ~ price, data = cigar_unb,
                              index = c("state", "year"))
# Returns: r_hat, singular values, SVT threshold, nu_used (NNR penalty from cross-validation)
```

### Analytical bias correction

Setting `bias_corr = TRUE` applies the analytical incidental-parameter correction of [Su, Wang & Wang (2025)](https://doi.org/10.2139/ssrn.5177283),
$\hat\beta^{abc} = \hat\beta - (NT)^{-1/2}\hat W_X^{-1}\hat b$,
where the bias vector $\hat b$ depends on the exogeneity assumption:

```r
# Strictly exogenous regressors (default): b3 + b4 + b5 + b6
fit_bc <- ife_unbalanced(sales ~ price, data = cigar_unb,
                          index = c("state", "year"), r = 2,
                          se = "standard", bias_corr = TRUE,
                          exog = "strict")

# Weakly exogenous regressors (e.g. lagged dep. var.): b2 + b3 + ... + b6
fit_dyn_unb <- ife_unbalanced(y ~ lag_y,
                               data = df_dynamic, index = c("i", "t"),
                               r = 2, se = "hac", bias_corr = TRUE,
                               exog = "weak")
```

| `exog =` | Bias terms | Typical application |
|----------|-----------|---------------------|
| `"strict"` (default) | $\hat b_3 + \hat b_4 + \hat b_5 + \hat b_6$ | Standard panel regressors |
| `"weak"` | $\hat b_2 + \hat b_3 + \hat b_4 + \hat b_5 + \hat b_6$ | Lagged dependent variable |

### Initialisation options

```r
# Default OLS initialisation (fast; works well for fill >= 60%)
fit_ols <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, init = "ols")

# Nuclear-norm regularisation (recommended for fill < 60% or r >= 3)
fit_nnr <- ife_unbalanced(sales ~ price, data = cigar_unb,
                           index = c("state", "year"), r = 2, init = "nnr")
```

---

## Function Reference

| Function | Panel type | Description |
|----------|-----------|-------------|
| `ife()` | Balanced | Fit IFE model (Bai 2009); returns coefficients, SEs, factors, loadings |
| `print.ife()` | Balanced | Formatted coefficient table and model summary |
| `ife_select_r()` | Balanced | Fit IFE for $r = 0, \ldots, r_{\max}$; compare IC1/2/3, IC(BIC), PC |
| `ife_unbalanced()` | Unbalanced | Fit IFE via EM with matrix completion; OLS or NNR initialisation; analytical SE and bias correction |
| `print.ife_unb()` | Unbalanced | Formatted coefficient table with bias components |
| `ife_select_r_unb()` | Unbalanced | SVT factor selection (singular value thresholding) |

### Key `ife()` arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `formula` | — | `outcome ~ covariate1 + ...` |
| `data` | — | Long-format `data.frame` (balanced) |
| `index` | — | `c("unit_col", "time_col")` |
| `r` | `1` | Number of interactive factors |
| `force` | `"two-way"` | Additive FE: `"none"`, `"unit"`, `"time"`, `"two-way"` |
| `se` | `"standard"` | SE type: `"standard"`, `"robust"`, `"cluster"` |
| `bias_corr` | `FALSE` | Apply analytical bias correction |
| `method` | `"static"` | `"static"` ([Bai 2009](https://doi.org/10.3982/ECTA6135)) or `"dynamic"` ([Moon & Weidner 2017](https://doi.org/10.1017/S0266466615000328)) |
| `M1` | `1L` | Lag bandwidth for dynamic $\hat B_1$ bias term |

### Key `ife_unbalanced()` arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `formula` | — | `outcome ~ covariate1 + ...` |
| `data` | — | Long-format `data.frame` (balanced or unbalanced) |
| `index` | — | `c("unit_col", "time_col")` |
| `r` | `1L` | Number of interactive factors |
| `force` | `"none"` | Additive FE: `"none"`, `"unit"`, `"time"`, `"two-way"` (jointly estimated with the factors; default differs from `ife()` as the unbalanced model is intercept-free by default) |
| `se` | `"standard"` | SE type: `"standard"`, `"robust"`, `"cluster"`, `"hac"` |
| `init` | `"ols"` | Initialisation: `"ols"` or `"nnr"` (nuclear-norm) |
| `bias_corr` | `FALSE` | Apply analytical incidental-parameter bias correction |
| `exog` | `"strict"` | Exogeneity: `"strict"` or `"weak"` (dynamic regressors) |
| `L_T` | `NULL` | HAC bandwidth ($\lfloor 2T^{1/5}\rfloor$ if `NULL`) |

---

## About

### Author
Binzhi Chen (University of Essex)

Email: <Binzhi.Chen9@gmail.com>

Web: [https://rickchen0910.github.io/](https://rickchen0910.github.io/)

### Citation

Please cite as follows:

Chen, B. (2026). xtife: Interactive Fixed Effects Estimator for Panel Data.
R package version 0.1.4. https://CRAN.R-project.org/package=xtife.

```bibtex
@Manual{xtife,
  title  = {{xtife}: Interactive Fixed Effects Estimator for Panel Data},
  author = {Binzhi Chen},
  year   = {2026},
  note   = {R package version 0.1.4},
  url    = {https://CRAN.R-project.org/package=xtife},
}
```

---

## References

Bai, J. (2009). Panel data models with interactive fixed effects. *Econometrica*, 77(4), 1229–1279. [doi:10.3982/ECTA6135](https://doi.org/10.3982/ECTA6135)

Bai, J. and Ng, S. (2002). Determining the number of factors in approximate factor models. *Econometrica*, 70(1), 191–221. [doi:10.1111/1468-0262.00273](https://doi.org/10.1111/1468-0262.00273)

Bai, J. and Ng, S. (2021). Matrix completion, counterfactuals, and factor analysis of missing data. *Journal of the American Statistical Association*, 116(536), 1746–1763. [doi:10.1080/01621459.2021.1967163](https://doi.org/10.1080/01621459.2021.1967163)

Baltagi, B.H. (1995). *Econometric Analysis of Panel Data*. Wiley.

Cameron, A.C., Gelbach, J.B. and Miller, D.L. (2011). Robust inference with multiway clustering. *Journal of Business & Economic Statistics*, 29(2), 238–249. [doi:10.1198/jbes.2010.07136](https://doi.org/10.1198/jbes.2010.07136)

Ditzen, J. and Karavias, Y. (2025). Interactive, Grouped and Non-separable Fixed Effects: A Practitioner's Guide to the New Panel Data Econometrics. arXiv:2507.19099. [doi:10.48550/arXiv.2507.19099](https://doi.org/10.48550/arXiv.2507.19099)

Mazumder, R., Hastie, T. and Tibshirani, R. (2010). Spectral regularization algorithms for learning large incomplete matrices. *Journal of Machine Learning Research*, 11, 2287–2322.

Moon, H.R. and Weidner, M. (2017). Dynamic linear panel regression models with interactive fixed effects. *Econometric Theory*, 33, 158–195. [doi:10.1017/S0266466615000328](https://doi.org/10.1017/S0266466615000328)

Newey, W.K. and West, K.D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica*, 55(3), 703–708. [doi:10.2307/1913610](https://doi.org/10.2307/1913610)

Su, L., Wang, F. and Wang, Y. (2025). Estimation and inference for interactive fixed effects panel data models with unbalanced panels. *SSRN Working Paper No. 5177283*. [doi:10.2139/ssrn.5177283](https://doi.org/10.2139/ssrn.5177283)

---

## License

GPL-2 | GPL-3 © 2026 Binzhi Chen
