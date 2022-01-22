#' @title quickns: A package for faster natural spline
#' @description This package provides a function \code{quick_ns} that can generate natural cubic spline basis faster.
#' @details This package provides a faster natural cubic spline basis function, \code{\link{quick_ns}}.
#' 
#' @docType package
#' @useDynLib quickns
#' @import splines
#' @import Rcpp
#' @importFrom utils getFromNamespace
#' @importFrom stats quantile
#' @author Weidong Wang, \email{wang@umass.edu}
#' @references R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
"_PACKAGE"