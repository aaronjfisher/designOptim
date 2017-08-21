designOptim
=============

[![Travis-CI Build Status](https://travis-ci.org/<USERNAME>/<REPO>.png?branch=master)](https://travis-ci.org/<USERNAME>/<REPO>)

An R package for building trials, simulating trials, and searching for optimal trials. The highest level user-facing functions are

* `buildTrial` for building trials, and
* `optimizeTrial` for searching for trials

To install and see help pages:

```S
## if needed
install.packages("devtools")

## main package
library(devtools)

install_github('aaronjfisher/designOptim', auth_token = auth_token)

library(designOptim)

## to access help pages
help(package=designOptim)

```



