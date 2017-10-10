designOptim
=============

[![Travis-CI Build Status](https://travis-ci.org/<USERNAME>/<REPO>.png?branch=master)](https://travis-ci.org/<USERNAME>/<REPO>)

An R package for building adaptive enrichment trials, simulating trials, and searching for optimal trials. The highest level user-facing functions are

* `buildTrial` for creating trials with controlled Type I error rates, and
* `optimizeTrial` for searching for trials

To install and see help pages:

```S
## if needed
install.packages("devtools")

## main package
library(devtools)

install_github('aaronjfisher/designOptim')

library(designOptim)

## to access help pages
help(package=designOptim)

```



