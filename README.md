# FastMix
R package for FastMix, a pipeline that combines the deconvolution step with the downstream analyses based on linear mixed eﬀects regression (LMER) model. 

To install this package, you can either: (a) download and unzip this package to a local directory, run R CMD build FastMix && R CMD INSTALL FastMix_<xxx>.tar.gz, or (b) use function `install_github()` (part of the `devtools` package) to install this package: `install_github("terrysun0302/FastMix", build_vignettes = TRUE)`

Once the package is installed, please use `library(FastMix)` to load this package and use `vignette("FastMix")` to view the package vignette.
