## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----ex-----------------------------------------------------------------------
library(PoissonBinomial)

# Case 1
dpbinom(NULL, rep(0.3, 7))
dbinom(0:7, 7, 0.3)
# equal results

# Case 2
dpbinom(NULL, c(0, 0, 0, 0, 0, 0, 0))
dpbinom(NULL, c(1, 1, 1, 1, 1, 1, 1))
dpbinom(NULL, c(0, 0, 0, 0, 1, 1, 1))

# Case 3
dpbinom(NULL, c(0, 0, 0.4, 0.2, 0.8, 0.1, 1))

## ----directconv---------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Convolve")
ppbinom(NULL, pp, wt, "Convolve")

## ----dividefft1---------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "DivideFFT")
ppbinom(NULL, pp, wt, "DivideFFT")

## ----dividefft2---------------------------------------------------------------
set.seed(1)
pp1 <- runif(751)
pp2 <- pp1[1:750]

sum(abs(dpbinom(NULL, pp2, method = "DivideFFT") - dpbinom(NULL, pp2, method = "Convolve")))
sum(abs(dpbinom(NULL, pp1, method = "DivideFFT") - dpbinom(NULL, pp1, method = "Convolve")))

## ----dividefft3---------------------------------------------------------------
set.seed(1)
pp1 <- runif(751)

d1 <- dpbinom(NULL, pp1, method = "DivideFFT")
d2 <- dpbinom(NULL, pp1, method = "Convolve")

min(d1[d1 > 0])
min(d2[d2 > 0])

## ----dftcf--------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Characteristic")
ppbinom(NULL, pp, wt, "Characteristic")

## ----rf1----------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Recursive")
ppbinom(NULL, pp, wt, "Recursive")

## ----rf2----------------------------------------------------------------------
set.seed(1)
pp <- runif(1000)
wt <- sample(1:10, 1000, TRUE)

sum(abs(dpbinom(NULL, pp, wt, "Convolve") - dpbinom(NULL, pp, wt, "Recursive")))

## ----benchmark----------------------------------------------------------------
library(microbenchmark)
set.seed(1)

f1 <- function() dpbinom(NULL, runif(4000), method = "DivideFFT")
f2 <- function() dpbinom(NULL, runif(4000), method = "Convolve")
f3 <- function() dpbinom(NULL, runif(4000), method = "Characteristic")
f4 <- function() dpbinom(NULL, runif(4000), method = "Recursive")

microbenchmark(f1(), f2(), f3(), f4())

