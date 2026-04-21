#' @title Dataset on US Cigarette Demand Panel
#'
#' @description Balanced panel of cigarette sales and prices across 46 US states
#' for 30 years (1963--1992). Originally used in Baltagi (1995) and
#' widely used as a benchmark dataset for panel estimators.
#'
#' @format A data frame with 1,380 rows and 9 variables:
#' \describe{
#'   \item{state}{US state identifier (integer, 1--46)}
#'   \item{year}{year (integer, 1963--1992)}
#'   \item{price}{cigarette price index}
#'   \item{pop}{state population}
#'   \item{pop16}{population aged 16 and over}
#'   \item{cpi}{consumer price index}
#'   \item{ndi}{per-capita disposable income}
#'   \item{sales}{per-capita cigarette sales (packs per person per year)}
#'   \item{pimin}{minimum cigarette price in adjoining states}
#' }
#'
#' @source Baltagi, B.H. (1995) \emph{Econometric Analysis of Panel Data}.
#'   Wiley. Distributed with the \pkg{plm} R package (Croissant and Millo 2008).
#'
#' @references
#' Baltagi, B.H. (1995). \emph{Econometric Analysis of Panel Data}. Wiley.
#'
#' Croissant, Y. and Millo, G. (2008). Panel data econometrics in R: the plm
#' package. \emph{Journal of Statistical Software}, 27(2), 1--43.
#' \doi{10.18637/jss.v027.i02}
"cigar"
