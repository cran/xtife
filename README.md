# xtife: Interactive Fixed Effects Estimator for Balanced Panel Data

<!-- badges: start -->
[![R CMD Check](https://github.com/Rickchen0910/xtife/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Rickchen0910/xtife/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/xtife)](https://CRAN.R-project.org/package=xtife)
[![License: GPL v2/v3](https://img.shields.io/badge/License-GPL%20v2%2Fv3-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
<!-- badges: end -->

`xtife` provides a pure base-R implementation of the **Interactive Fixed Effects (IFE)** panel estimator of Bai (2009) with full analytical standard errors, asymptotic bias correction, and information-criterion-based factor number selection. No external dependencies beyond base R are required. 
For a comprehensive review about interactive fixed effect, please refer to Ditzen, J., & Karavias, Y. (2025).

---

## The Model

Standard two-way fixed effects (TWFE) assumes unobserved heterogeneity enters additively. IFE generalises this by allowing unobserved confounders to interact across units and time:

$$y_{it} = \alpha_i + \xi_t + X_{it}'\beta + \lambda_i'F_t + u_{it}$$

where $F_t \in \mathbb{R}^r$ are common factors and $\lambda_i \in \mathbb{R}^r$ are unit-specific loadings. Setting $r = 0$ reduces the model to standard TWFE.

---

## Features

| Feature | Details |
|---------|---------|
| **Estimator** | Bai (2009) SVD-based alternating projections |
| **Standard errors** | Homoskedastic, HC1 robust, cluster-robust by unit |
| **Bias correction** | Bai (2009) static; Moon & Weidner (2017) dynamic |
| **Factor selection** | IC1, IC2, IC3 (Bai & Ng 2002); IC(BIC), PC (Bai 2009) |
| **Dynamic extension** | Predetermined regressors (Moon & Weidner 2017) |
| **Dependencies** | Base R only (`stats`) |
| **Panel type** | Balanced panels |

---

## Installation

```r
# From CRAN (once available)
install.packages("xtife")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("Rickchen0910/xtife")
```

---

## Quick Start

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

---

## Standard Error Types

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
| `"cluster"` | Cluster-robust by unit | Serial correlation within units |

---

## Factor Number Selection

```r
sel <- ife_select_r(sales ~ price, data = cigar,
                    index = c("state", "year"),
                    r_max = 6,
                    force = "two-way")
```

Prints a table of IC1, IC2, IC3 (Bai & Ng 2002), IC(BIC), and PC (Bai 2009) criteria for each candidate $r$. The recommended criterion for panels with $\min(N, T) < 60$ is **IC(BIC)**.

---

## Asymptotic Bias Correction

```r
# Static bias correction (Bai 2009)
fit_bc <- ife(sales ~ price, data = cigar,
              index = c("state", "year"), r = 2,
              bias_corr = TRUE)

# Dynamic bias correction (Moon & Weidner 2017)
# Use when regressors include lagged dependent variables
fit_dyn <- ife(sales ~ price, data = cigar,
               index = c("state", "year"), r = 2,
               method    = "dynamic",
               bias_corr = TRUE,
               M1        = 1L)
```

For the cigar panel ($N = 46$, $T = 30$, $T/N \approx 0.65$):

| Estimator | Price coefficient |
|-----------|------------------|
| IFE (r = 2) | −0.5242 |
| IFE + Bai (2009) bias correction | −0.5309 |
| IFE dynamic + Moon & Weidner (2017) bias correction | −0.5343 |

---

## Comparison with TWFE

Setting `r = 0` recovers the standard two-way FE estimator, identical to `lm()` with unit and time dummies at machine precision:

```r
fit0 <- ife(sales ~ price, data = cigar,
            index = c("state", "year"), r = 0)
# Equivalent to plm(..., model = "within", effect = "twoways")
```

---

## Function Reference

| Function | Description |
|----------|-------------|
| `ife()` | Fit IFE model; returns coefficients, SEs, factors, loadings |
| `print.ife()` | Print formatted coefficient table and model info |
| `ife_select_r()` | Fit IFE for r = 0, …, r_max and compare information criteria |

### Key `ife()` arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `formula` | — | `outcome ~ covariate1 + ...` |
| `data` | — | Long-format data.frame |
| `index` | — | `c("unit_col", "time_col")` |
| `r` | `1` | Number of interactive factors |
| `force` | `"two-way"` | Additive FE: `"none"`, `"unit"`, `"time"`, `"two-way"` |
| `se` | `"standard"` | SE type: `"standard"`, `"robust"`, `"cluster"` |
| `bias_corr` | `FALSE` | Apply analytical bias correction |
| `method` | `"static"` | `"static"` (Bai 2009) or `"dynamic"` (Moon & Weidner 2017) |
| `M1` | `1L` | Lag bandwidth for dynamic B1 bias term |

---

# About

### Author
Binzhi Chen (University of Essex)

Email: <Binzhi.Chen9@gmail.com>

Web: [https://sites.google.com/view/binzhichen/home](https://sites.google.com/view/binzhichen/home)


### Citation

Please cite as follows:

Chen, B. (2026). xtife: Interactive Fixed Effects Estimator for Balanced Panel Data.
R package version 0.1.0. https://github.com/Rickchen0910/xtife.

Or

@Manual{xtife,
  title  = {{xtife}: Interactive Fixed Effects Estimator for Balanced Panel Data},
  author = {Binzhi Chen},
  year   = {2026},
  note   = {R package version 0.1.0},
  url    = {https://github.com/Rickchen0910/xtife},
}


## References

Bai, J. (2009). Panel data models with interactive fixed effects. *Econometrica*, 77(4), 1229–1279. [doi:10.3982/ECTA6135](https://doi.org/10.3982/ECTA6135)

Bai, J. and Ng, S. (2002). Determining the number of factors in approximate factor models. *Econometrica*, 70(1), 191–221. [doi:10.1111/1468-0262.00273](https://doi.org/10.1111/1468-0262.00273)

Baltagi, B.H. (1995). *Econometric Analysis of Panel Data*. Wiley.

Ditzen, J., & Karavias, Y. (2025). Interactive, Grouped and Non-separable Fixed Effects: A Practitioner's Guide to the New Panel Data Econometrics. arXiv preprint arXiv:2507.19099. https://doi.org/10.48550/arXiv.2507.19099

Moon, H.R. and Weidner, M. (2017). Dynamic linear panel regression models with interactive fixed effects. *Econometric Theory*, 33, 158–195. [doi:10.1017/S0266466615000328](https://doi.org/10.1017/S0266466615000328)

---

## License

GPL-2 | GPL-3 © 2026 Binzhi Chen
