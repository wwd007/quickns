#' quick_ns
#' 
#' @param x the predictor variable. Missing values are allowed.
#' @param df degrees of freedom. One can supply \code{df} rather than knots; \code{quick_ns()} then chooses \code{df - 1 - intercept} knots at suitably chosen quantiles of \code{x} (which will ignore missing values). The default, \code{df = NULL}, sets the number of inner knots as \code{length(knots)}.
#' @param knots breakpoints that define the spline. The default is no knots; together with the natural boundary conditions this results in a basis for linear regression on \code{x}. Typical values are the mean or median for one knot, quantiles for more knots. See also \code{Boundary.knots}.
#' @param intercept if \code{TRUE}, an intercept is included in the basis; default is \code{FALSE}.
#' @param Boundary.knots boundary points at which to impose the natural boundary conditions and anchor the B-spline basis (default the range of the data). If both \code{knots} and \code{Boundary.knots} are supplied, the basis parameters do not depend on \code{x}. Data can extend beyond \code{Boundary.knots}.
#' @seealso \code{\link[splines]{ns}} for more details.
#' @author Weidong Wang, \email{wang@umass.edu}
#' @references R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#' @useDynLib quickns spline_basis
#' @export
quick_ns <- function (x, df = NULL, knots = NULL, intercept = FALSE, 
                      Boundary.knots = range(x, na.rm = TRUE)) 
{
  splineDesign2 <- function (knots, x, ord = 4L, derivs = 0L) {
    #C_spline_basis <- utils::getFromNamespace("C_spline_basis", "splines")
    temp <- .Call(spline_basis, knots, ord, x, derivs)
    make_design2(nx = length(x), ncoef = length(knots) - ord, ord, temp)
  }
  splineDesign_no_intcpt <- function (knots, x, ord = 4L, derivs = 0L) { 
    #C_spline_basis <- utils::getFromNamespace("C_spline_basis", "splines")
    temp <- .Call(spline_basis, knots, ord, x, derivs)
    make_design_no_intcpt_cpp(nx = length(x), ncoef = length(knots) - ord,
                              ord = ord, temp)
  }
  ns_core <- function (x, Aknots, Boundary.knots) {
    basis <- splineDesign_no_intcpt(Aknots, x, ord = 4L)
    const <- splineDesign_no_intcpt(Aknots, Boundary.knots, ord = 4L, derivs = c(2L, 2L))
    get_basis(const, basis)[, -(1L:2L), drop = FALSE]
  }
  
  # Boundary.knots = range(x, na.rm = TRUE)
  ncore <- 4
  nx <- names(x)
  # x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  } else {
    if (length(x) == 1L) 
      Boundary.knots <- x * c(7, 9)/8
    FALSE
  }
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L
    knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                       2L)[-c(1L, nIknots + 2L)]
    knots <- quantile(x[!outside], knots)
  } else nIknots <- length(knots)
  Aknots <- sort(c(rep.int(Boundary.knots, 4L), knots))
  const <- splineDesign_no_intcpt(Aknots, Boundary.knots, ord = 4L, derivs = c(2L, 2L))
  
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splineDesign2(Aknots, rep(k.pivot, 2L), 4, 
                          c(0, 1))
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splineDesign2(Aknots, rep(k.pivot, 2L), 4, 
                          c(0, 1))
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside)) 
      basis[inside, ] <- splineDesign2(Aknots, x[inside], 
                                       4)
    if (T) { # !intercept
      basis <- get_basis(const, basis[, -1, drop = FALSE])[, -(1L:2L), drop = FALSE]
    }
  } else if (length(x) >= 3e7L & R.Version()$platform == "x86_64-pc-linux-gnu") {
    basis <- do.call(
      rbind,
      parallel::mclapply(parallel::splitIndices(length(x), ncore), 
                         function(i) ns_core(x[i], Aknots, Boundary.knots), 
                         mc.cores = ncore))
  } else {
    basis <- ns_core(x, Aknots, Boundary.knots)
  }
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = knots, 
            Boundary.knots = Boundary.knots, intercept = FALSE)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}
