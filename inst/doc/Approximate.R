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
dpbinom(NULL, c(0, 0, 0.4, 0.2, 0.8, 0.1, 1), method = "RefinedNormal")

## ----pa1----------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Poisson")
ppbinom(NULL, pp, wt, "Poisson")

## ----pa2----------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
ppbinom(NULL, pp, method = "Poisson")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "Poisson") - ppbinom(NULL, pp))

# U(0, 0.01) random probabilities of success
pp <- runif(20, 0, 0.01)
ppbinom(NULL, pp, method = "Poisson")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "Poisson") - ppbinom(NULL, pp))

## ----am1----------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
mean(rep(pp, wt))

dpbinom(NULL, pp, wt, "Mean")
ppbinom(NULL, pp, wt, "Mean")

## ----am2----------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
ppbinom(NULL, pp, method = "Mean")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "Mean") - ppbinom(NULL, pp))

# U(0.4, 0.6) random probabilities of success
pp <- runif(20, 0.3, 0.5)
ppbinom(NULL, pp, method = "Mean")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "Mean") - ppbinom(NULL, pp))

# U(0.49, 0.51) random probabilities of success
pp <- runif(20, 0.39, 0.41)
ppbinom(NULL, pp, method = "Mean")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "Mean") - ppbinom(NULL, pp))

## ----gma1---------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
prod(rep(pp, wt))^(1/sum(wt))

dpbinom(NULL, pp, wt, "GeoMean")
ppbinom(NULL, pp, wt, "GeoMean")

## ----gma2---------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
ppbinom(NULL, pp, method = "GeoMean")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "GeoMean") - ppbinom(NULL, pp))

# U(0.4, 0.6) random probabilities of success
pp <- runif(20, 0.4, 0.6)
ppbinom(NULL, pp, method = "GeoMean")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "GeoMean") - ppbinom(NULL, pp))

# U(0.49, 0.51) random probabilities of success
pp <- runif(20, 0.49, 0.51)
ppbinom(NULL, pp, method = "GeoMean")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "GeoMean") - ppbinom(NULL, pp))

## ----gmb1---------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
1 - prod(1 - rep(pp, wt))^(1/sum(wt))

dpbinom(NULL, pp, wt, "GeoMeanCounter")
ppbinom(NULL, pp, wt, "GeoMeanCounter")

## ----gmb2---------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
ppbinom(NULL, pp, method = "GeoMeanCounter")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "GeoMeanCounter") - ppbinom(NULL, pp))

# U(0.4, 0.6) random probabilities of success
pp <- runif(20, 0.4, 0.6)
ppbinom(NULL, pp, method = "GeoMeanCounter")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "GeoMeanCounter") - ppbinom(NULL, pp))

# U(0.49, 0.51) random probabilities of success
pp <- runif(20, 0.49, 0.51)
ppbinom(NULL, pp, method = "GeoMeanCounter")
ppbinom(NULL, pp)
summary(ppbinom(NULL, pp, method = "GeoMeanCounter") - ppbinom(NULL, pp))

## ----na1----------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
mean(rep(pp, wt))

dpbinom(NULL, pp, wt, "Normal")
ppbinom(NULL, pp, wt, "Normal")

## ----na2----------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(10)
summary(ppbinom(NULL, pp, method = "Normal") - ppbinom(NULL, pp))

# U(0.4, 0.6) random probabilities of success
pp <- runif(1000)
summary(ppbinom(NULL, pp, method = "Normal") - ppbinom(NULL, pp))

# U(0.49, 0.51) random probabilities of success
pp <- runif(100000)
summary(ppbinom(NULL, pp, method = "Normal") - ppbinom(NULL, pp))

## ----rna1---------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
mean(rep(pp, wt))

dpbinom(NULL, pp, wt, "RefinedNormal")
ppbinom(NULL, pp, wt, "RefinedNormal")

## ----rna2---------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(10)
summary(ppbinom(NULL, pp, method = "RefinedNormal") - ppbinom(NULL, pp))

# U(0.4, 0.6) random probabilities of success
pp <- runif(1000)
summary(ppbinom(NULL, pp, method = "RefinedNormal") - ppbinom(NULL, pp))

# U(0.49, 0.51) random probabilities of success
pp <- runif(100000)
summary(ppbinom(NULL, pp, method = "RefinedNormal") - ppbinom(NULL, pp))

## ----benchmark----------------------------------------------------------------
library(microbenchmark)
set.seed(1)

f1 <- function() dpbinom(NULL, runif(4000), method = "Normal")
f2 <- function() dpbinom(NULL, runif(4000), method = "RefinedNormal")
f3 <- function() dpbinom(NULL, runif(4000), method = "Poisson")
f4 <- function() dpbinom(NULL, runif(4000), method = "Mean")
f5 <- function() dpbinom(NULL, runif(4000), method = "GeoMean")
f6 <- function() dpbinom(NULL, runif(4000), method = "GeoMeanCounter")
f7 <- function() dpbinom(NULL, runif(4000), method = "DivideFFT")

microbenchmark(f1(), f2(), f3(), f4(), f5(), f6(), f7())

