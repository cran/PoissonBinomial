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
  
  # if x = NULL, return all possible probabilities
  if(is.null(x)) x <- transf$compl.range
  
  # identify valid 'x' values (invalid ones will have 0-probability)
  idx.x <- which(x %in% transf$compl.range)
  
  ## compute probabilities
  # vector for storing the probabilities
  d <- double(length(x))
  
  # no computation needed, if there are no valid observations in 'x'
  if(length(idx.x)){
    # select valid observations in relevant range
    y <- x[idx.x]
    
    # relevant observations
    idx.y <- which(y %in% transf$inner.range)
    
    # if no input value is in relevant range, they are impossible (i.e. return 0-probabilities)
    if(length(idx.y)){
      # transformed input parameters
      probs <- transf$probs
      val_p <- transf$val_p
      val_q <- transf$val_q
      n <- transf$n
      
      # select and rescale relevant observations
      z <- y[idx.y] - transf$guaranteed
    
      if(n == 0){
        # 'probs' contains only zeros and ones, i.e. only one possible observation
        d[idx.x][idx.y] <- 1
      }else if(n == 1){
        # 'probs' contains only one value that is not 0 or 1, i.e. a Bernoulli distribution
        idx.z <- which(z %in% c(val_q, val_p))
        if(length(idx.z)) d[idx.x][idx.y][idx.z] <- c(1 - probs, probs)[match(z[idx.z], c(val_q, val_p))]
      }else{
        if(all(val_p == val_p[1]) && all(val_q == val_q[1])){
          # all values of 'probs' are equal, i.e. a standard binomial distribution
          u <- 0:n * val_p[1] + n:0 * val_q[1]
          idx.u <- which(u %in% z)
          idx.z <- which(z %in% u)
          if(length(idx.u) && length(idx.z)) d[idx.x][idx.y][idx.z] <- dpbinom(idx.u - 1, probs, method = method)
        }else{
          # compute distribution according to 'method'
          d[idx.x][idx.y] <- switch(method,
                                    DivideFFT = dgpb_dc(z, probs, val_p, val_q),
                                    Convolve = dgpb_conv(z, probs, val_p, val_q),
                                    Characteristic = dgpb_dftcf(z, probs, val_p, val_q),
                                    Normal = dgpb_na(z, probs, val_p, val_q, FALSE),
                                    RefinedNormal = dgpb_na(z, probs, val_p, val_q, TRUE))
        }
      }
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
  
  if(all(val_p == 1) && all(val_q == 0)) return(ppbinom(x, probs, wts, method, lower.tail, log.p))
  if(all(val_p == 0) && all(val_q == 1)) return(ppbinom(x, 1 - probs, wts, method, lower.tail, log.p))
  if(all(val_p == 0 | val_p == 1) && all(val_q == 1 - val_p)){
    probs[val_p == 0] <- 1 - probs[val_p == 0]
    return(ppbinom(x, probs, wts, method, lower.tail, log.p))
  }
  
  ## transform input to relevant range
  transf <- transform.GPB(x, probs, val_p, val_q, wts)
  
  # if x = NULL, return all possible probabilities
  if(is.null(x)) x <- transf$compl.range
  
  # identify valid 'x' values (invalid ones will have 0-probability)
  idx.x <- which(x %in% transf$compl.range)
  
  ## compute probabilities
  # vector for storing the probabilities
  d <- rep(as.numeric(!lower.tail), length(x))
  
  # no computation needed, if there are no valid observations in 'x'
  if(length(idx.x)){
    # select valid observations in relevant range
    y <- x[idx.x]
    
    # relevant observations
    idx.y <- which(y %in% transf$inner.range)
    
    # which valid observations are outside relevant range
    idx.z <- which(y > max(transf$inner.range))
    
    if(length(idx.y)){
      # transformed input parameters
      probs <- transf$probs
      val_p <- transf$val_p
      val_q <- transf$val_q
      n <- transf$n
      
      # select and rescale relevant observations
      z <- y[idx.y] - transf$guaranteed
      
      if(n == 0){
        # 'probs' contains only zeros and ones, i.e. only one possible observation
        d[idx.x][idx.y] <- 1
      }else if(n == 1){
        # 'probs' contains only one value that is not 0 or 1, i.e. a Bernoulli distribution
        v <- min(val_q, val_p):max(val_q, val_p)
        idx.v <- which(z %in% v)
        if(length(idx.v)){
          pr <- numeric(length(v))
          if(lower.tail){
            pr[1] <- 1 - probs
            pr[2:length(v)] <- 1
          }else pr[1] <- probs
          d[idx.x][idx.y][idx.v] <- pr[match(z[idx.v], v)]
        }
      }else{
        if(all(val_p == val_p[1]) && all(val_q == val_q[1])){
          # all values of 'probs' are equal, i.e. a standard binomial distribution
          u <- 0:n * val_p[1] + n:0 * val_q[1]
          idx.v <- which(z %in% u)
          idx.u <- which(u %in% z)
          if(length(idx.u) && length(idx.v)) d[idx.x][idx.y][idx.v] <- ppbinom(idx.u - 1, probs, method = method, lower.tail = lower.tail)
        }else{
          # compute distribution according to 'method'
          d[idx.x][idx.y] <- switch(method,
                                    DivideFFT = pgpb_dc(z, probs, val_p, val_q, lower.tail),
                                    Convolve = pgpb_conv(z, probs, val_p, val_q, lower.tail),
                                    Characteristic = pgpb_dftcf(z, probs, val_p, val_q, lower.tail),
                                    Normal = pgpb_na(z, probs, val_p, val_q, FALSE, lower.tail),
                                    RefinedNormal = pgpb_na(z, probs, val_p, val_q, TRUE, lower.tail))
          
          # compute counter-probabilities, if necessary
          #if(!lower.tail && method < "Normal") d[idx.x][idx.y] <- 1 - d[idx.x][idx.y]
        }
      }
    }
    # fill cumulative probabilities of values above the relevant range
    if(length(idx.z)) d[idx.x][idx.z] <- as.double(lower.tail)
  }
  
  # fill cumulative probabilities of values above n
  d[x > max(transf$compl.range)] <- as.double(lower.tail)
  
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
  
  # order 'p' depending on if they are lower tail probabilities; save order
  ord <- order(p, decreasing = !lower.tail)
  p.s <- p[ord]
  # negative sign, if lower.tail = TRUE
  sign <- (-1)^as.numeric(!lower.tail)
  cdf <- sign*cdf
  # compute quantiles
  for(i in 1:length(p.s)){
    # handle 0's and 1's
    if(p.s[i] == as.numeric(!lower.tail)){
      res[i] <- ifelse(lower.tail, 0, first)#(1 - lower.tail) * max(0, (transf$complete.range[1] == transf$inner.range[1])) - first
      next
    }
    if(p.s[i] == as.numeric(lower.tail)){
      res[i] <- ifelse(lower.tail, last, max(transf$compl.range))#len - 1L
      next
    }
    # find the cdf value smaller than or equal to current 'p.s' value
    while(pos < len && cdf[pos] < sign*p.s[i]) pos <- pos + 1L
    # save respective observation (-1 because position k is for observation k - 1)
    res[i] <- pos - 1L + first
  }
  # arrange results in original order of 'p'
  res <- res[order(ord)]
  
  # return results
  return(res)
}

#'@rdname GenPoissonBinomial-Distribution
#'@importFrom stats runif
#'@export
rgpbinom <- function(n, probs, val_p, val_q, wts = NULL, method = "DivideFFT"){
  len <- length(n)
  if(len > 1) n <- len
  
  # check if 'n' is NULL
  if(is.null(n)) stop("'n' must not be NULL!")
  
  # generate random probabilities
  p <- runif(n)
  
  # compute quantiles (does checking for the other variables)
  res <- qgpbinom(p, probs, val_p, val_q, wts, method)
  
  # return results
  return(res)
}