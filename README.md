
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PostCard

<!-- badges: start -->

[![R-CMD-check](https://github.com/NNpackages/PostCard/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NNpackages/PostCard/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/NNpackages/PostCard/graph/badge.svg)](https://app.codecov.io/gh/NNpackages/PostCard)
<!-- badges: end -->

PostCard provides tools for accurately estimating marginal effects using
plug-in estimation with GLMs, including increasing precision using
prognostic covariate adjustment.

### Estimating marginal effects

`rctglm()` is used to estimate marginal effects. See introductory
examples of its usage in `vignette("PostCard")` and more details in
`vignette("model-fit")`.

*Marginal effects* are causal effects of the form $r(\Psi_1, \Psi_0)$,
where $\Psi_a=\mathbb{E}[Y(a)]$ are population mean outcomes under
exposure $a=0, 1$, respectively. These are sometimes referred to as
*counterfactual means*. The package uses *plug-in estimation for robust
estimation of any marginal effect estimand* as well as *influence
functions for robust estimation of the variance of the estimand*
(Rosenblum, M. and M. J. van der Laan, 2010: Simple, efficient
estimators of treatment effects in randomized trials using generalized
linear models to leverage baseline variables. The International Journal
of Biostatistics, 6, no. 1).

### Prognostic covariate adjustment

`rctglm_with_prognosticscore()` is used to estimate marginal effects
including the use of prognostic covariate adjustment. See introductory
examples of its usage in `vignette("PostCard")` and more details in
`vignette("model-fit")`.

*Prognostic covariate adjustment* involves training a *prognostic model*
on historical data to predict the response in that data. Assuming that
the historical data is representative of the comparator group in a
*“new”* data set, we can use the prognostic model to predict the
*comparator counterfactual outcome* for all observations (including the
ones in the comparator group). This prediction, which is called the
**prognostic score** is then used as an adjustment covariate in the GLM.

## Power approximation for linear models

Additionally, implementations of sample size/power approximation
formulas are available, enabling retrospective power analyses to be
performed using the package. Functionalities are available at the help
page `powerss()`, introductory examples in `vignette("PostCard")` and
more details in `vignette("retrospective-power")`.

## Installation

You can install the development version of PostCard from GitHub with:

``` r
pak::pak("NNpackages/PostCard")
```
