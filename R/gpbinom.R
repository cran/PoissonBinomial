#'@name GenPoissonBinomial-Distribution
#'
#'@title The Generalized Poisson Binomial Distribution
#'
#'@description
#'Density, distribution function, quantile function and random generation for
#'the generalized Poisson binomial distribution with probability vector
#'\code{probs}.
#'
#'@param x           Either a vector of observed sums or NULL. If NULL,
#'                   probabilities of all possible observations are
#'                   returned.
#'@param p           Vector of probabilities for computation of quantiles.
#'@param n           Number of observations. If \code{length(n) > 1}, the
#'                   length is taken to be the number required.
#'@param probs       Vector of probabilities of success of each Bernoulli
#'                   trial.
#'@param val_p       Vector of values that each trial produces with probability
#'                   in \code{probs}.
#'@param val_q       Vector of values that each trial produces with probability
#'                   in \code{1 - probs}.
#'@param method      Character string that specifies the method of computation
#'                   and must be one of \code{"DivideFFT"}, \code{"Convolve"}, 
#'                   \code{"Characteristic"}, \code{"Normal"} or
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
#'algorithms are derived and described in Biscarri, Zhao & Brunner (2018). They
#'have been modified for use with the generalized Poisson binomial
#'distribution. The
#'\emph{Discrete Fourier Transformation of the Characteristic Function}
#'(\code{"Characteristic"}) is derived in Zhang, Hong & Balakrishnan (2018),
#'the \emph{Normal Approach} (\code{"Normal"}) and the
#'\emph{Refined Normal Approach} (\code{"RefinedNormal"}) are described in Hong
#'(2013). They were slightly adapted for the generalized Poisson binomial
#'distribution.
#'
#'In some special cases regarding the values of \code{probs}, the \code{method}
#'parameter is ignored (see Introduction vignette).
#'
#'@return
#'\code{dgpbinom} gives the density, \code{pgpbinom} computes the distribution
#'function, \code{qgpbinom} gives the quantile function and \code{rgpbinom}
#'generates random deviates.
#'
#'For \code{rgpbinom}, the length of the result is determined by \code{n}, and
#'is the lengths of the numerical arguments for the other functions.
#'
#'@section References:
#'Hong, Y. (2018). On computing the distribution function for the Poisson
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
#'Zhang, M., Hong, Y. and Balakrishnan, N. (2018). The generalized 
#'    Poisson-binomial distribution and the computation of its distribution
#'    function. \emph{Journal of Statistical Computational and Simulation},
#'    \strong{88}(8), pp. 1515-1527. doi:
#'    \href{https://doi.org/10.1080/00949655.2018.1440294}{
#'    10.1080/00949655.2018.1440294}
#'    
#'@examples
#'set.seed(1)
#'pp <- c(1, 0, runif(10), 1, 0, 1)
#'qq <- seq(0, 1, 0.01)
#'va <- rep(5, length(pp))
#'vb <- 1:length(pp)
#'
#'dgpbinom(NULL, pp, va, vb, method = "DivideFFT")
#'pgpbinom(75:100, pp, va, vb, method = "DivideFFT")
#'qgpbinom(qq, pp, va, vb, method = "DivideFFT")
#'rgpbinom(100, pp, va, vb, method = "DivideFFT")
#'
#'dgpbinom(NULL, pp, va, vb, method = "Convolve")
#'pgpbinom(75:100, pp, va, vb, method = "Convolve")
#'qgpbinom(qq, pp, va, vb, method = "Convolve")
#'rgpbinom(100, pp, va, vb, method = "Convolve")
#'
#'dgpbinom(NULL, pp, va, vb, method = "Characteristic")
#'pgpbinom(75:100, pp, va, vb, method = "Characteristic")
#'qgpbinom(qq, pp, va, vb, method = "Characteristic")
#'rgpbinom(100, pp, va, vb, method = "Characteristic")
#'
#'dgpbinom(NULL, pp, va, vb, method = "Normal")
#'pgpbinom(75:100, pp, va, vb, method = "Normal")
#'qgpbinom(qq, pp, va, vb, method = "Normal")
#'rgpbinom(100, pp, va, vb, method = "Normal")
#'
#'dgpbinom(NULL, pp, va, vb, method = "RefinedNormal")
#'pgpbinom(75:100, pp, va, vb, method = "RefinedNormal")
#'qgpbinom(qq, pp, va, vb, method = "RefinedNormal")
#'rgpbinom(100, pp, va, vb, method = "RefinedNormal")
#'
#'@export
dgpbinom <- function(x, probs, val_p, val_q, wts = NULL, method = "DivideFFT", log = FALSE){
  ## preliminary checks
  method <- check.args.GPB(x, probs, val_p, val_q, wts, method)
  
  if(all(val_p == 1) && all(val_q == 0)) return(dpbinom(x, probs, wts, method, log))
  if(all(val_p == 0) && all(val_q == 1)) return(dpbinom(x, 1 - probs, wts, method, log))
  if(all(val_p == 0 | val_p == 1) && all(val_q == 1 - val_p)){
    probs[val_p == 0] <- 1 - probs[val_p == 0]
    return(dpbinom(x, probs, wts, method, log))
  }
  
  ## transform input to relevant range
  transf <- transform.GPB(x, probs, val_p, val_q, wts)
  probs <- transf$probs
  val_p <- transf$val_p
  val_q <- transf$val_q
  n <- transf$n
  
  # if x = NULL, return all possible probabilities
  if(is.null(x)) x <- transf$compl.range
  
  # identify valid 'x' values (invalid ones will have 0-probability)
  idx.x <- which(x %in% transf$compl.range)
  
  # select valid observations in relevant range
  y <- x[idx.x]
  idx.y <- which(y %in% transf$inner.range)
  
  ## compute probabilities
  # vector for storing the probabilities
  d <- double(length(x))
  
  # if no input value is in relevant range, they are impossible (i.e. return 0-probabilities)
  if(!length(idx.y)) return(d)
  
  z <- y[idx.y] - transf$guaranteed
  
  if(n == 0){
    # 'probs' contains only zeros and ones, i.e. only one possible observation
    if(length(idx.y)) d[idx.x][idx.y] <- 1
  }else if(n == 1){
    # 'probs' contains only one value that is not 0 or 1, i.e. a Bernoulli distribution
    idx.z <- which(z %in% c(val_q, val_p))
    if(length(idx.y)) d[idx.x][idx.y][idx.z] <- c(1 - probs, probs)[match(z[idx.z], c(val_q, val_p))]
  }else{
    if(all(val_p == val_p[1]) && all(val_q == val_q[1])){
      # all values of 'probs' are equal, i.e. a standard binomial distribution
      u <- 0:n * val_p[1] + n:0 * val_q[1]
      idx.z <- which(z %in% u)
      idx.u <- which(u %in% z)
      if(length(idx.y) && length(idx.u) && length(idx.z)) d[idx.x][idx.y][idx.z] <- dpbinom(idx.u - 1, probs, method = method)
    }else{
      # compute distribution according to 'method'
      if(length(idx.y)) d[idx.x][idx.y] <- switch(method, DivideFFT = dgpb_dc(z, probs, val_p, val_q),
                                                  Convolve = dgpb_conv(z, probs, val_p, val_q),
                                                  Characteristic = dgpb_dftcf(z, probs, val_p, val_q),
                                                  Normal = dgpb_na(z, probs, val_p, val_q, FALSE),
                                                  RefinedNormal = dgpb_na(z, probs, val_p, val_q, TRUE))
    }
  }
  
  # logarithm, if required
  if(log) d <- log(d)
  
  # return results
  return(d)
}

#'@rdname GenPoissonBinomial-Distribution
#'@export
pgpbinom <- function(x, probs, val_p, val_q, wts = NULL, method = "DivideFFT", lower.tail = TRUE, log.p = FALSE){
  ## preliminary checks
  method <- check.args.GPB(x, probs, val_p, val_q, wts, method)
  
  if(all(val_p == 1) && all(val_q == 0)) return(ppbinom(x, probs, wts, method, log))
  if(all(val_p == 0) && all(val_q == 1)) return(ppbinom(x, 1 - probs, wts, method, log))
  if(all(val_p == 0 | val_p == 1) && all(val_q == 1 - val_p)){
    probs[val_p == 0] <- 1 - probs[val_p == 0]
    return(ppbinom(x, probs, wts, method, log))
  }
  
  ## transform input to relevant range
  transf <- transform.GPB(x, probs, val_p, val_q, wts)
  probs <- transf$probs
  val_p <- transf$val_p
  val_q <- transf$val_q
  n <- transf$n
  
  # if x = NULL, return all possible probabilities
  if(is.null(x)) x <- transf$compl.range
  
  # identify valid 'x' values (invalid ones will have 0-probability)
  idx.x <- which(x %in% transf$compl.range)
  
  # select valid observations in relevant range
  y <- x[idx.x]
  idx.y <- which(y %in% transf$inner.range)
  
  ## compute probabilities
  # vector for storing the probabilities
  d <- double(length(x))
  
  # if no input value is in relevant range, they are impossible (i.e. return 0-probabilities)
  if(!length(idx.y)) return(d)
  
  z <- y[idx.y] - transf$guaranteed
  idx.z <- if(length(idx.y)) which(y > max(transf$inner.range)) else integer(0)
  
  if(n == 0){
    # 'probs' contains only zeros and ones, i.e. only one possible observation
    if(length(idx.y)) d[idx.x][idx.y] <- 1
  }else if(n == 1){
    # 'probs' contains only one value that is not 0 or 1, i.e. a Bernoulli distribution
    idx.z <- which(z %in% c(val_q, val_p))
    if(length(idx.y)) d[idx.x][idx.y][idx.z] <- c(1 - probs, 1)[match(z[idx.z], c(val_q, val_p))]
  }else{
    if(all(val_p == val_p[1]) && all(val_q == val_q[1])){
      # all values of 'probs' are equal, i.e. a standard binomial distribution
      u <- 0:n * val_p[1] + n:0 * val_q[1]
      idx.z <- which(z %in% u)
      idx.u <- which(u %in% z)
      if(length(idx.y) && length(idx.u) && length(idx.z)) d[idx.x][idx.y][idx.z] <-ppbinom(idx.u - 1, probs, method = method)
    }else{
      # compute distribution according to 'method'
      if(length(idx.y)) d[idx.x][idx.y] <- switch(method, DivideFFT = pgpb_dc(z, probs, val_p, val_q),
                                                  Convolve = pgpb_conv(z, probs, val_p, val_q),
                                                  Characteristic = pgpb_dftcf(z, probs, val_p, val_q),
                                                  Normal = pgpb_na(z, probs, val_p, val_q, FALSE),
                                                  RefinedNormal = pgpb_na(z, probs, val_p, val_q, TRUE))
    }
  }
  # values above the relevant region must have a cumulative probability of 1
  if(length(idx.z)) d[idx.x][idx.z] <- 1
  
  # values above n must have a cumulative probability of 1
  d[x > max(transf$compl.range)] <- 1
  
  # compute lower-tail counterparts, if necessary
  if(!lower.tail) d <- 1 - d
  
  # logarithm, if required
  if(log.p) d <- log(d)
  
  # return results
  return(d)
}

#'@rdname GenPoissonBinomial-Distribution
#'@export
qgpbinom <- function(p, probs, val_p, val_q, wts = NULL, method = "DivideFFT", lower.tail = TRUE, log.p = FALSE){
  ## preliminary checks
  method <- check.args.GPB(NULL, probs, val_p, val_q, wts, method)
  
  # check if 'q' contains only probabilities
  if(!log.p){
    if(is.null(p) || any(is.na(p) | p < 0 | p > 1))
      stop("'p' must contain real numbers between 0 and 1!")
  }else{
    if(is.null(p) || any(is.na(p) | p > 0))
      stop("'p' must contain real numbers between -Inf and 0!")
  }
  
  ## transform input to relevant range
  transf <- transform.GPB(NULL, probs, val_p, val_q, wts)
  probs <- transf$probs
  val_p <- transf$val_p
  val_q <- transf$val_q
  
  ## compute probabilities (does checking for the other variables)
  cdf <- pgpbinom(NULL, probs, val_p, val_q, NULL, method, lower.tail)
  
  # bounds of relevant observations
  first <- min(transf$inner.range)
  last <- max(transf$inner.range)
  
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
      while(pos < len && cdf[pos] <= pr) pos <- pos + 1L
      # save respective observation (-1 because position k is for observation k - 1)
      res[idx] <- pos
    }
    idx <- which(p == 0)
    if(length(idx)) res[idx] <- as.numeric(transf$compl.range[1] == transf$inner.range[1])
    idx <- which(p == 1)
    if(length(idx)) res[idx] <- len
  }else{
    # considering only the unique values in 'p' and sorting them saves time
    for(pr in sort(unique(p[p > 0 & p < 1]), decreasing = TRUE)){
      # identify positions of the current probability
      idx <- which(p == pr)
      # find the cdf value smaller than or equal to current 'pr' value
      while(pos < len && cdf[pos] >= pr) pos <- pos + 1L
      # save respective observation (-1 because position k is for observation k - 1)
      res[idx] <- pos
    }
    idx <- which(p == 1)
    if(length(idx)) res[idx] <- as.numeric(transf$compl.range[1] == transf$inner.range[1])
    idx <- which(p == 0)
    if(length(idx)) res[idx] <- len
  }
  
  # return results
  return(res + first - 1)
}

#'@rdname GenPoissonBinomial-Distribution
#'@importFrom stats runif
#'@export
rgpbinom <- function(n, probs, val_p, val_q, wts = NULL, method = "DivideFFT"){
  len <- length(n)
  if(len > 1) n <- len
  
  ## preliminary checks
  # check if 'n' is NULL
  if(is.null(n)) stop("'n' must not be NULL!")
  
  ## compute random numbers
  # generate random probabilities
  p <- runif(n)
  
  ## compute quantiles (does checking for the other variables)
  res <- qgpbinom(p, probs, val_p, val_q, wts, method)
  
  # return results
  return(res)
}