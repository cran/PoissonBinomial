# PoissonBinomial 1.1.1

* Fixed bugs in `ppbinom` and `pgpbinom` that caused incorrect calculation of
  logarithms and cumulative upper-tail probabilities.

# PoissonBinomial 1.1

* Added exact and approximate algorithms for the generalized Poisson binomial
  distribution described in Zhang, Hong & Balakrishnan (2018). The
  non-generalized distribution is now referred to as the 'ordinary' Poisson
  binomial distribution.
* Restructured vignettes. Added tables of content and fixed smaller issues.
* Minor bug fixes for `dbinom`, `ppbinom` and `qpbinom` functions.

# PoissonBinomial 1.0.2-1

* Fixes and improvements of the vignettes; no code changes.

# PoissonBinomial 1.0.2

* Improvements of C++ helper function "norm_dpb" to achieve better
  normalization.
* Bug fix of DFT-CF method ("Characteristic") so that negative probabilities
  are no longer possible.
* Reworked vignette structure.
* Added author acknowledgments to the Makevars.win file (original author was
  Geoff99 (https://github.com/Geoff99)).
  

# PoissonBinomial 1.0.1

* Fixed a bug in the C++ helper function "norm_dpb" that could cause infinite
  loops (the function is invisible to the user, since it is only used in the
  C++ domain).
  

# PoissonBinomial 1.0.0

* Initial release.
* Added a `NEWS.md` file to track changes to the package.
