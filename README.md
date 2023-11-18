[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ETAS)](https://CRAN.R-project.org/package=ETAS)
[![CRAN_Download_Count](http://cranlogs.r-pkg.org/badges/ETAS)](https://CRAN.R-project.org/package=ETAS)
[![Build Status](https://travis-ci.org/jalilian/ETAS.svg?branch=master)](https://travis-ci.org/jalilian/ETAS)

# Temporal point process model to identify secondary crashes

The R code for the ETAS model obtained from Jalilian is adapted to identify the secondary crash events in the present work. The code was developed to analyze an earthquake catalog using the stochastic declustering approach. We modified the code for a temporal-only self-exciting point process for the current application. The input data contains the date and time of the crash events. We have classified potential secondary crashes from the dataset based on queue time and corresponding probability values obtained from the model.


## Installation

To install the package from [CRAN](https://CRAN.R-project.org/package=ETAS), run the following in R:
```R
install.packages('ETAS')
```

You can also install the current version of the package on GitHub by running:
```R
require(remotes)
install_github('jalilian/ETAS')
```

If [remotes](https://github.com/mangothecat/remotes) is not installed, you should first run:

```R
install.packages('remotes')
```
 
## Parallel computing

Computations of the conditional intensity function, the log-likelihood function, declustering probabilities and the Davidon-Fletcher-Powell algorithm for optimization are all written in C code. As of version 0.3, a new C++ code is implemented using the [Rcpp](http://www.rcpp.org/) package which allows multi-thread parallel computing on multi-core processors with OpenMP and [suported platforms](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#OpenMP-support). The argument `nthreads` in `etas` function determines the number of threads to be used in the parallel region of the code. If `nthreads = 1` (the default), then a serial version of the C++ code carries out the computations. The `detectCores` function in [parallel](http://stat.ethz.ch/R-manual/R-devel/library/parallel/html/parallel-package.html) package can be consulted to find out the overall number of available threads on a given machine:
```R
parallel::detectCores()
```
Parallel computing (`nthreads > 1`) reduces the computation time for large earthquake catalogs. However, resource usage and limitations should be considered when setting `nthreads`.
