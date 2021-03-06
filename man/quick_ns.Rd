% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quick_ns.R
\name{quick_ns}
\alias{quick_ns}
\title{quick_ns}
\usage{
quick_ns(
  x,
  df = NULL,
  knots = NULL,
  intercept = FALSE,
  Boundary.knots = range(x, na.rm = TRUE)
)
}
\arguments{
\item{x}{the predictor variable. Missing values are allowed.}

\item{df}{degrees of freedom. One can supply \code{df} rather than knots; \code{quick_ns()} then chooses \code{df - 1 - intercept} knots at suitably chosen quantiles of \code{x} (which will ignore missing values). The default, \code{df = NULL}, sets the number of inner knots as \code{length(knots)}.}

\item{knots}{breakpoints that define the spline. The default is no knots; together with the natural boundary conditions this results in a basis for linear regression on \code{x}. Typical values are the mean or median for one knot, quantiles for more knots. See also \code{Boundary.knots}.}

\item{intercept}{if \code{TRUE}, an intercept is included in the basis; default is \code{FALSE}.}

\item{Boundary.knots}{boundary points at which to impose the natural boundary conditions and anchor the B-spline basis (default the range of the data). If both \code{knots} and \code{Boundary.knots} are supplied, the basis parameters do not depend on \code{x}. Data can extend beyond \code{Boundary.knots}.}
}
\description{
quick_ns
}
\references{
R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
}
\seealso{
\code{\link[splines]{ns}} for more details.
}
\author{
Weidong Wang, \email{wang@umass.edu}
}
