
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdmTMB

[![Travis build
status](https://travis-ci.org/pbs-assess/sdmTMB.svg?branch=master)](https://travis-ci.org/pbs-assess/sdmTMB)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

The sdmTMB package implements spatiotemporal GLMM (Generalized Linear
Mixed Effects Model) predictive process models using Template Model
Builder ([TMB](https://github.com/kaskr/adcomp)),
[R-INLA](http://www.r-inla.org/), and Gaussian Markov random fields. One
common application is for spatial or spatiotemporal species distribution
models (SDMs).

## Installation

You can install sdmTMB with:

``` r
devtools::install_github("pbs-assess/sdmTMB")
```

## Functionality

sdmTMB:

  - Fits GLMMs with spatial, spatiotemporal, spatial and spatiotemporal,
    or AR1 spatiotemporal Gaussian Markov random fields with TMB.
  - Uses formula interfaces for fixed effects and any time-varying
    effects (dynamic regression) (e.g. `formula = y ~ 1 + x1,
    time_varying = ~ 0 + x2`).
  - Uses a `family(link)` format similar to `glm()` or `lme4::lmer()`.
    This includes Gaussian, Poisson, negative binomial, gamma, binomial,
    lognormal, Student-t, and Tweedie distributions with identity, log,
    inverse, and logit links. E.g. `family = tweedie(link = "log")`.
  - Has `predict()` and `residuals()` methods. The residuals are
    randomized-quantile residuals similar to those implemented in the
    [DHARMa](https://cran.r-project.org/package=DHARMa) package. The
    `predict()` function can take a `newdata` argument similar to `lm()`
    or `glm()` etc. The predictions are full predictive-process
    predictions (i.e. they make smooth pretty maps).
  - Includes functionality for estimating the centre of gravity or total
    biomass by time step for index standardization.
  - Implements multi-phase estimation for speed.
  - Can optionally allow for anisotropy in the random fields (spatial
    correlation that is directionally dependent).
  - Can generate an SPDE predictive-process mesh based on a clustering
    algorithm and R-INLA or can take any standard R-INLA mesh created
    externally as input.

## Example

The main function is `sdmTMB()`. See `?sdmTMB` and `?predict.sdmTMB` for
the most complete examples. There is also a simulation function `?sim`.
