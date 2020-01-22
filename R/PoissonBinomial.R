#'@name Poisson-Binomial-package
#'
#'@title Exact and Approximate Implementations for Computing the Poisson Binomial Distribution
#'
#'@description
#'This package implements various algorithms for computing the probability mass
#'function the cumulative distribution function, quantiles and random numbers
#'of the Poisson binomial distribution.
#'
#'@docType package
#'@import Rcpp
#'@useDynLib PoissonBinomial, .registration = TRUE
#'
#'@section References:
#'Hong, Y. (2013). On computing the distribution function for the Poisson
#'    binomial distribution. \emph{Computational Statistics & Data Analysis},
#'    \strong{59}, pp. 41-51. doi:
#'    \href{https://doi.org/10.1016/j.csda.2012.10.006}{
#'    10.1016/j.csda.2012.10.006}
#'
#'Biscarri, W., Zhao, S. D. and Brunner, R. J. (2018) A simple and fast method
#'    for computing the Poisson binomial distribution.
#'    \emph{Computational Statistics and Data Analysis}, \strong{31}, pp.
#'    216–222. doi: \href{https://doi.org/10.1016/j.csda.2018.01.007}{
#'    10.1016/j.csda.2018.01.007}
#'    
#'@examples
#'set.seed(1)
#'pp <- c(0, 0, runif(995), 1, 1, 1)
#'dpbinom(NULL, pp)
#'ppbinom(450:550, pp, method = "DivideFFT")
#'qpbinom(pp, pp, method = "Convolve")
#'rpbinom(10, pp, method = "RefinedNormal")
"_PACKAGE"