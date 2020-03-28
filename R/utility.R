check.args.GPB <- function(x, probs, val_p, val_q, wts, method, log.p = FALSE){
  # check if 'x' contains only integers
  if(!is.null(x) && any(abs(x - round(x)) > 1e-7))
    stop("'x' must contain integers only!")
  
  # check if 'probs' contains only probabilities
  if(is.null(probs) || any(is.na(probs) | probs < 0 | probs > 1))
    stop("'probs' must contain real numbers between 0 and 1!")
  
  # number of probabilities
  n <- length(probs)
  
  # check if 'val_p' and 'val_q' have the same length as 'probs'
  if(length(val_p) != n || length(val_q) != n) stop("'probs', 'val_p' and 'val_q' must have the same length!")
  
  if(!is.null(wts) && length(wts) != n)
    stop("'probs' and 'wts' (if not NULL) must have the same length!")
  
  # check if 'val_p' contains only integers
  if(!is.null(val_p) && any(abs(val_p - round(val_p)) > 1e-7))
    stop("'val_p' must contain integers only!")
  
  # check if 'val_q' contains only integers
  if(!is.null(val_q) && any(abs(val_q - round(val_q)) > 1e-7))
    stop("'val_q' must contain integers only!")
  
  # check if 'wts' contains only integers (zeros are allowed)
  if(!is.null(wts) && any(is.na(wts) | wts < 0 | abs(wts - round(wts)) > 1e-07))
    stop("'wts' must contain non-negative integers!")
  
  # make sure that the value of 'method' matches one of the implemented procedures
  method <- match.arg(method, c("DivideFFT", "Convolve", "Characteristic", "Normal", "RefinedNormal"))
  
  # if all checks were successful, return matched 'method'
  return(method)
}


transform.GPB <- function(x, probs, val_p, val_q, wts){
  # number of probabilities
  n <- length(probs)
  
  ## expand 'probs', 'val_p' and 'val_q' according to the counts in 'wts'
  # if 'wts' is NULL, set it to be a vector of ones
  if(is.null(wts))
    wts <- rep(1, n)
  
  # expand 'probs', 'val_p', 'val_q'
  probs <- rep(probs, wts)
  val_p <- rep(val_p, wts)
  val_q <- rep(val_q, wts)
  
  # re-compute length of 'probs' (= sum of 'wts')
  n <- sum(wts)
  
  ## determine relevant range of observations
  # determine minimum and maximum possible observations
  sum_min <- sum(pmin(val_p, val_q))
  sum_max <- sum(pmax(val_p, val_q))
  
  # which probabilities are 0 or 1, which val_p and val_q are equal
  idx.0 <- which(probs == 0)
  idx.1 <- which(probs == 1)
  idx.v <- which(val_p == val_q & probs > 0 & probs < 1)
  idx.r <- setdiff(1:n, union(union(idx.0, idx.1), idx.v))
  
  # guaranteed
  val_p_sure <- val_p[idx.1]
  val_q_sure <- val_q[idx.0]
  vals_equal <- val_p[idx.v]# == val_q[idx.v]
  sum_sure <- sum(val_p_sure, val_q_sure, vals_equal)
  
  # limit 'probs', 'val_p' and 'val_q' to relevant range
  probs <- probs[idx.r]
  val_p <- val_p[idx.r]
  val_q <- val_q[idx.r]
  np <- length(idx.r)
  
  # bounds of relevant observations
  sum_min_in <- sum(pmin(val_p, val_q)) + sum_sure
  sum_max_in <- sum(pmax(val_p, val_q)) + sum_sure
  
  return(list(probs = probs, val_p = val_p, val_q = val_q, guaranteed = sum_sure, compl.range = sum_min:sum_max, inner.range = sum_min_in:sum_max_in, n = np))
}
