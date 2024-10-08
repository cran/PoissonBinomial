#'@name PoissonBinomial-package
#'
#'@title Efficient Exact and Approximate Implementations for Computing Ordinary and Generalized Poisson Binomial Distributions
#'
#'@description
#'This package implements various algorithms for computing the probability mass
#'function, the cumulative distribution function, quantiles and random numbers
#'of both ordinary and generalized Poisson binomial distributions.
#'
#'@import Rcpp
#'@useDynLib PoissonBinomial, .registration = TRUE
#'
#'@section References:
#'Hong, Y. (2013). On computing the distribution function for the Poisson
#'    binomial distribution. \emph{Computational Statistics & Data Analysis},
#'    \strong{59}, pp. 41-51. \doi{10.1016/j.csda.2012.10.006}
#'
#'Biscarri, W., Zhao, S. D. and Brunner, R. J. (2018) A simple and fast method
#'    for computing the Poisson binomial distribution.
#'    \emph{Computational Statistics and Data Analysis}, \strong{31}, pp.
#'    216–222. \doi{10.1016/j.csda.2018.01.007}
#'    
#'Zhang, M., Hong, Y. and Balakrishnan, N. (2018). The generalized 
#'    Poisson-binomial distribution and the computation of its distribution
#'    function. \emph{Journal of Statistical Computational and Simulation},
#'    \strong{88}(8), pp. 1515-1527. \doi{10.1080/00949655.2018.1440294}
#'    
#'@examples
#'# Functions for ordinary Poisson binomial distributions
#'set.seed(1)
#'pp <- c(1, 0, runif(10), 1, 0, 1)
#'qq <- seq(0, 1, 0.01)
#'
#'dpbinom(NULL, pp)
#'ppbinom(7:10, pp, method = "DivideFFT")
#'qpbinom(qq, pp, method = "Convolve")
#'rpbinom(10, pp, method = "RefinedNormal")
#'
#'# Functions for generalized Poisson binomial distributions
#'va <- rep(5, length(pp))
#'vb <- 1:length(pp)
#'
#'dgpbinom(NULL, pp, va, vb, method = "Convolve")
#'pgpbinom(80:100, pp, va, vb, method = "Convolve")
#'qgpbinom(qq, pp, va, vb, method = "Convolve")
#'rgpbinom(100, pp, va, vb, method = "Convolve")
"_PACKAGE"