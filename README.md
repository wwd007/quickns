
<!-- README.md is generated from README.Rmd. Please edit that file -->

# quickns

## Overview

quickns provides a useful function, `quick_ns()`, which generates a
B-spline basis matrix for a natural cubic spline faster. The usage is
totally the same as the basic `ns()` function in the base R `splines`
package.

The R code of `quick_ns()` is adapted from `ns()`. Some heavy
computation process in `ns()` originally written in R is rewritten in
C++ for better efficiency. Ideally, `quick_ns()` could save about 50%
time than `ns()` on large data.
