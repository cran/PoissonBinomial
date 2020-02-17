# PoissonBinomial 1.0.2

* Improvements of C++ helper function "norm_dpb" to achieve better
  normalization.
* Bug fix of DFT-CF method ("Characteristic") so that negative probabilities
  are no longer possible.
* Reworked vignette structure.
* Added author acknowledgments to the Makevars.win file (original author was
  Geoff99 (https://github.com/Geoff99)).
  

# PoissonBinomial 1.0.1

* Fixed a bug in the C++ helper function "norm_dpb" that could cause infinite.
  loops (the function is invisible to the user, since it is only used in the
  C++ domain).
  

# PoissonBinomial 1.0.0

* Initial release.
* Added a `NEWS.md` file to track changes to the package.
