#'@name Poisson-Binomial-Distribution
#'
#'@importFrom stats dbinom pbinom runif
#'
#'@title The Poisson Binomial Distribution
#'
#'@description
#'Density, distribution function, quantile function and random generation for
#'the Poisson binomial distribution with probability vector \code{probs}.
#'
#'@param x           Either a vector of observed numbers of successes or NULL.
#'                   If NULL, probabilities of all possible observations are
#'                   returned.
#'@param p           Vector of probabilities for computation of quantiles.
#'@param n           Number of observations. If \code{length(n) > 1}, the
#'                   length is taken to be the number required.
#'@param probs       Vector of probabilities of success of each Bernoulli
#'                   trial.
#'@param method      Character string that specifies the method of computation
#'                   and must be one of \code{"DivideFFT"}, \code{"Convolve"},
#'                   \code{"Characteristic"}, \code{"Recursive"},
#'                   \code{"Mean"}, \code{"GeoMean"}, \code{"GeoMeanCounter"},
#'                   \code{"Poisson"}, \code{"Normal"} or
#'                   \code{"RefinedNormal"} (abbreviations are allowed).
#'@param wts         Vector of non-negative integer weights for the input
#'                   probabilities.
#'@param log,log.p   Logical value indicating if results are given as
#'                   logarithms.
#'@param lower.tail  Logical value indicating if results are \eqn{P[X \leq x]}
#'                   (if \code{TRUE}; default) or \eqn{P[X > x]} (if 
#'                   \code{FALSE}).
#'
#'@details
#'See the references for computational details. The \emph{Divide and Conquer}
#'(\code{"DivideFFT"}) and \emph{Direct Convolution} (\code{"Convolve"})
#'algorithms are derived and described in Biscarri et al. (2018). The 
#'\emph{Discrete Fourier Transformation of the Characteristic Function}
#'(\code{"Characteristic"}), the \emph{Recursive Formula} (\code{"Recursive"}),
#'the \emph{Poisson Approximation} (\code{"Poisson"}), the
#'\emph{Normal Approach} (\code{"Normal"}) and the
#'\emph{Refined Normal Approach} (\code{"RefinedNormal"}) are described in Hong
#'(2013). The calculation of the \emph{Recursive Formula} was modified to
#'overcome the excessive memory requirements of Hong's implementation.
#'
#'The \code{"Mean"} method is a naive binomial approach using the arithmetic
#'mean of the probabilities of success. Similarly, the \code{"GeoMean"} and
#'\code{"GeoMeanCounter"} procedures are binomial approximations, too, but
#'they form the geometric mean of the probabilities of success
#'(\code{"GeoMean"}) and their counter probabilities (\code{"GeoMeanCounter"}),
#'respectively.
#'
#'In some special cases regarding the values of \code{probs}, the \code{method}
#'parameter is ignored.
#'\describe{
#'  \item{All values are 0 or 1:}{The Distribution is not random can can only
#'        attain one value, which is the number of ones, \eqn{n_1}.}
#'  \item{Only one value is different from 0 or 1:}{The Distribution
#'        essentially is a Bernoulli Distribution and can only attain \eqn{n_1}
#'        or \eqn{n_1 + 1}.}
#'  \item{All values are equal:}{The Distribution is an ordinary Binomial
#'        Distribution.}
#'}
#'
#'@return
#'\code{dpbinom} gives the density, \code{ppbinom} gives the distribution
#'function, \code{qpbinom} gives the quantile function and \code{rpbinom}
#'generates random deviates.
#'
#'For \code{rpbinom}, the length of the result is determined by \code{n}, and
#'is the lengths of the numerical arguments for the other functions.
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
#'qq <- seq(0, 1, 0.01)
#'
#'dpbinom(NULL, pp, method = "DivideFFT")
#'ppbinom(450:550, pp, method = "DivideFFT")
#'qpbinom(qq, pp, method = "DivideFFT")
#'rpbinom(100, pp, method = "DivideFFT")
#'
#'dpbinom(NULL, pp, method = "Convolve")
#'ppbinom(450:550, pp, method = "Convolve")
#'qpbinom(qq, pp, method = "Convolve")
#'rpbinom(100, pp, method = "Convolve")
#'
#'dpbinom(NULL, pp, method = "Characteristic")
#'ppbinom(450:550, pp, method = "Characteristic")
#'qpbinom(qq, pp, method = "Characteristic")
#'rpbinom(100, pp, method = "Characteristic")
#'
#'dpbinom(NULL, pp, method = "Recursive")
#'ppbinom(450:550, pp, method = "Recursive")
#'qpbinom(qq, pp, method = "Recursive")
#'rpbinom(100, pp, method = "Recursive")
#'
#'dpbinom(NULL, pp, method = "Mean")
#'ppbinom(450:550, pp, method = "Mean")
#'qpbinom(qq, pp, method = "Mean")
#'rpbinom(100, pp, method = "Mean")
#'
#'dpbinom(NULL, pp, method = "GeoMean")
#'ppbinom(450:550, pp, method = "GeoMean")
#'qpbinom(qq, pp, method = "GeoMean")
#'rpbinom(100, pp, method = "GeoMean")
#'
#'dpbinom(NULL, pp, method = "GeoMeanCounter")
#'ppbinom(450:550, pp, method = "GeoMeanCounter")
#'qpbinom(qq, pp, method = "GeoMeanCounter")
#'rpbinom(100, pp, method = "GeoMeanCounter")
#'
#'dpbinom(NULL, pp, method = "Poisson")
#'ppbinom(450:550, pp, method = "Poisson")
#'qpbinom(qq, pp, method = "Poisson")
#'rpbinom(100, pp, method = "Poisson")
#'
#'dpbinom(NULL, pp, method = "Normal")
#'ppbinom(450:550, pp, method = "Normal")
#'qpbinom(qq, pp, method = "Normal")
#'rpbinom(100, pp, method = "Normal")
#'
#'dpbinom(NULL, pp, method = "RefinedNormal")
#'ppbinom(450:550, pp, method = "RefinedNormal")
#'qpbinom(qq, pp, method = "RefinedNormal")
#'rpbinom(100, pp, method = "RefinedNormal")
NULL

#'@rdname Poisson-Binomial-Distribution
#'@export
dpbinom <- function(x, probs, wts = NULL, method = "DivideFFT", log = FALSE){
  ## preliminary checks
  # number of probabilities
  n <- length(probs)
  
  # check if 'x' contains only integers
  if(!is.null(x) && any(abs(x - round(x)) > 1e-7))
    stop("'x' must contain integers only!")
  
  # check if 'probs' contains only probabilities
  if(is.null(probs) || any(is.na(probs) | probs < 0 | probs > 1))
    stop("'probs' must contain real numbers between 0 and 1!")
  
  # make sure that the value of 'method' matches one of the implemented procedures
  method <- match.arg(method, c("DivideFFT", "Convolve", "Characteristic", "Recursive", "Mean", "GeoMean", "GeoMeanCounter", "Poisson", "Normal", "RefinedNormal"))
  
  # check if 'wts' contains only integers (zeros are allowed)
  if(!is.null(wts) && any(is.na(wts) | wts < 0 | abs(wts - round(wts)) > 1e-07))
    stop("'wts' must contain non-negative integers!")
  
  if(!is.null(wts) && length(wts) != n)
    stop("'probs' and 'wts' (if not NULL) must have the same length!")
  
  ## expand 'probs' according to the counts in 'wts'
  # if 'wts' is NULL, set it to be a vector of ones
  if(is.null(wts))
    wts <- rep(1, n)
  
  # expand 'probs'
  probs <- rep(probs, wts)
  
  # re-compute length of 'probs' (= sum of 'wts')
  n <- sum(wts)
  
  # if x = NULL, return all possible probabilities
  if(is.null(x)) x <- 0:n
  
  # identify valid 'x' values (invalid ones will have 0-probability)
  idx.x <- which(x >= 0 & x <= n)
  
  # select valid observations
  y <- x[idx.x]
  
  ## compute probabilities
  # vector for storing the probabilities
  d <- double(length(x))
  
  # which probabilities are 0 or 1
  idx0 <- which(probs == 0)
  idx1 <- which(probs == 1)
  probs <- probs[probs > 0 & probs < 1]
  
  # number of zeros and ones
  n0 <- length(idx0)
  n1 <- length(idx1)
  np <- n - n0 - n1
  
  # relevant observations
  idx.y <- which(y %in% n1:(n - n0))
  z <- y[idx.y] - n1
  
  if(np == 0){
    # 'probs' contains only zeros and ones, i.e. only one possible observation
    if(length(idx.y)) d[idx.x][idx.y] <- 1
  }else if(np == 1){
    # 'probs' contains only one value that is not 0 or 1, i.e. a Bernoulli distribution
    if(length(idx.y)) d[idx.x][idx.y] <- c(1 - probs, probs)[z + 1]
  }else{
    if(all(probs == probs[1])){
      # all values of 'probs' are equal, i.e. a standard binomial distribution
      if(length(idx.y)) d[idx.x][idx.y] <- dbinom(z, np, probs)
    }else{
      # otherwise, compute distribution according to 'method'
      if(length(idx.y)) d[idx.x][idx.y] <- switch(method, DivideFFT = dpb_dc(z, probs),
                                                  Convolve = dpb_conv(z, probs),
                                                  Characteristic = dpb_dftcf(z, probs),
                                                  Recursive = dpb_rf(z, probs),
                                                  Mean = dpb_mean(z, probs),
                                                  GeoMean = dpb_gmba(z, probs, FALSE),
                                                  GeoMeanCounter = dpb_gmba(z, probs, TRUE),
                                                  Poisson = dpb_pa(z, probs),
                                                  Normal = dpb_na(z, probs, FALSE),
                                                  RefinedNormal = dpb_na(z, probs, TRUE))
    }
  }
  
  # return results
  return(d)
}

#'@rdname Poisson-Binomial-Distribution
#'@export
ppbinom <- function(x, probs, wts = NULL, method = "DivideFFT", lower.tail = TRUE, log.p = FALSE){
  ## preliminary checks
  # number of probabilities
  n <- length(probs)
  
  # check if 'x' contains only integers
  if(!is.null(x) && any(abs(x - round(x)) > 1e-7))
    stop("'x' must contain integers only!")
  
  # check if 'probs' contains only probabilities
  if(is.null(probs) || any(is.na(probs) | probs < 0 | probs > 1))
    stop("'probs' must contain real numbers between 0 and 1!")
  
  # make sure that the value of 'method' matches one of the implemented procedures
  method <- match.arg(method, c("DivideFFT", "Convolve", "Characteristic", "Recursive", "Mean", "GeoMean", "GeoMeanCounter", "Poisson", "Normal", "RefinedNormal"))
  
  # check if 'wts' contains only integers (zeros are allowed)
  if(!is.null(wts) && any(is.na(wts) | wts < 0 | abs(wts - round(wts)) > 1e-07))
    stop("'wts' must contain non-negative integers!")
  
  if(!is.null(wts) && length(wts) != n)
    stop("'probs' and 'wts' (if not NULL) must have the same length!")
  
  ## expand 'probs' according to the counts in 'wts'
  # if 'wts' is NULL, set it to be a vector of ones
  if(is.null(wts))
    wts <- rep(1, n)
  
  # expand 'probs'
  probs <- rep(probs, wts)
  
  # re-compute length of 'probs' (= sum of 'wts')
  n <- sum(wts)
  
  # if x = NULL, return all possible probabilities
  if(is.null(x)) x <- 0:n
  
  # identify valid 'x' values (invalid ones will have 0-probability)
  idx.x <- which(x >= 0 & x <= n)
  
  # select valid observations
  y <- x[idx.x]
  
  ## compute probabilities
  # vector for storing the probabilities
  d <- double(length(x))
  
  # which probabilities are 0 or 1
  idx0 <- which(probs == 0)
  idx1 <- which(probs == 1)
  probs <- probs[probs > 0 & probs < 1]
  
  # number of zeros and ones
  n0 <- length(idx0)
  n1 <- length(idx1)
  np <- n - n0 - n1
  
  # relevant observations
  idx.y <- which(y %in% n1:(n - n0))
  z <- y[idx.y] - n1
  idx.z <- if(length(idx.y)) which(y > n - n0) else integer(0)
  
  if(np == 0){
    # 'probs' contains only zeros and ones, i.e. only one possible observation
    if(length(idx.y)) d[idx.x][idx.y] <- 1
  }else if(np == 1){
    # 'probs' contains only one value that is not 0 or 1, i.e. a Bernoulli distribution
    if(length(idx.y)) d[idx.x][idx.y] <- c(1 - probs, 1)[z + 1]
  }else{
    if(all(probs == probs[1])){
      # all values of 'probs' are equal, i.e. a standard binomial distribution
      if(length(idx.y)) d[idx.x][idx.y] <- pbinom(z, np, probs[1])
    }else{
      # otherwise, compute distribution according tho 'method'
      if(length(idx.y)) d[idx.x][idx.y] <- switch(method, DivideFFT = ppb_dc(z, probs),
                                                  Convolve = ppb_conv(z, probs),
                                                  Characteristic = ppb_dftcf(z, probs),
                                                  Recursive = ppb_rf(z, probs),
                                                  Mean = ppb_mean(z, probs),
                                                  GeoMean = ppb_gmba(z, probs, FALSE),
                                                  GeoMeanCounter = ppb_gmba(z, probs, TRUE),
                                                  Poisson = ppb_pa(z, probs),
                                                  Normal = ppb_na(z, probs, FALSE),
                                                  RefinedNormal = ppb_na(z, probs, TRUE))
    }
  }
  # values above the relevant region must have a cumulative probability of 1
  if(length(idx.z)) d[idx.x][idx.z] <- 1
  
  # values above n must have a cumulative probability of 1
  d[x > n] <- 1
  
  # compute lower-tail counterparts, if necessary
  if(!lower.tail) d <- 1 - d
  
  # logarithm, if required
  if(log.p) d <- log(d)
  
  # return results
  return(d)
}

#'@rdname Poisson-Binomial-Distribution
#'@export
qpbinom <- function(p, probs, wts = NULL, method = "DivideFFT", lower.tail = TRUE, log.p = FALSE){
  ## preliminary checks
  # check if 'p' contains only probabilities
  if(!log.p){
    if(is.null(p) || any(is.na(p) | p < 0 | p > 1))
      stop("'p' must contain real numbers between 0 and 1!")
  }else{
    if(is.null(p) || any(is.na(p) | p > 0))
      stop("'p' must contain real numbers between -Inf and 0!")
  }
  
  # make sure that the value of 'method' matches one of the implemented procedures
  method <- match.arg(method, c("DivideFFT", "Convolve", "Characteristic", "Recursive", "Mean", "GeoMean", "GeoMeanCounter", "Poisson", "Normal", "RefinedNormal"))
  
  ## compute probabilities (does checking for the other variables)
  cdf <- ppbinom(NULL, probs, wts, method, lower.tail)
  # limit to first appearance of 1
  idx1 <- which(cdf == 1)
  idx0 <- which(cdf == 0)
  if(lower.tail && length(idx1)){
    cdf <- cdf[1:min(idx1)]
  }else if(!lower.tail && length(idx0)){
    cdf <- cdf[1:min(idx0)]
  }
  
  # size of distribution
  size <- length(probs)
  
  # length of cdf
  len <- length(cdf)
  
  # logarithm, if required
  if(log.p) p <- exp(p)
  
  ## compute quantiles
  # starting position in cdf
  pos <- 1L
  # reserve vector to store results
  res <- integer(length(p))
  
  if(lower.tail){
    # considering only the unique values in 'p' and sorting them saves time
    for(pr in sort(unique(p[p > 0 & p < 1]))){
      # identify positions of the current probability
      idx <- which(p == pr)
      # find the cdf value greater than current 'pr' value
      while(cdf[pos] <= pr && pos < len) pos <- pos + 1L
      # save respective observation (-1 because position k is for observation k - 1)
      res[idx] <- pos - 1L
    }
    idx <- which(p == 0)
    if(length(idx)) res[idx] <- length(which(probs == 1))
    idx <- which(p == 1)
    if(length(idx)){
      res[idx] <- size - length(which(probs == 0))
      if(method %in% c("Poisson", "Normal", "RefinedNormal")) res[idx] <- min(size, res[idx] + 1)
    }
  }else{
    # considering only the unique values in 'p' and sorting them saves time
    for(pr in sort(unique(p[p > 0 & p < 1]), decreasing = TRUE)){
      # identify positions of the current probability
      idx <- which(p == pr)
      # find the cdf value smaller than or equal to current 'pr' value
      while(cdf[pos] >= pr && pos < len) pos <- pos + 1L
      # save respective observation (-1 because position k is for observation k - 1)
      res[idx] <- pos - 1L
    }
    idx <- which(p == 1)
    if(length(idx)) res[idx] <- length(which(probs == 1))
    idx <- which(p == 0)
    if(length(idx)){
      res[idx] <- size - length(which(probs == 0))
      if(method %in% c("Poisson", "Normal", "RefinedNormal")) res[idx] <- min(size, res[idx] + 1)
    }
  }
  
  # return results
  return(res)
}

#'@rdname Poisson-Binomial-Distribution
#'@export
rpbinom <- function(n, probs, wts = NULL, method = "DivideFFT"){
  len <- length(n)
  if(len > 1) n <- len
  
  ## preliminary checks
  # check if 'n' is NULL
  if(is.null(n)) stop("'n' must not be NULL!")
  
  ## compute random numbers
  # generate random probabilities
  p <- runif(n)
  
  ## compute quantiles (does checking for the other variables)
  res <- qpbinom(p, probs, wts, method)
  
  # return results
  return(res)
}