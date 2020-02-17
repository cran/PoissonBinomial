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
dpbinom(NULL, c(0, 0, 0.4, 0.2, 0.8, 0.1, 1), method = "Convolve")

