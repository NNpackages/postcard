
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PostCard

<!-- badges: start -->

[![R-CMD-check](https://github.com/NNpackages/PostCard/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NNpackages/PostCard/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/NNpackages/PostCard/graph/badge.svg)](https://app.codecov.io/gh/NNpackages/PostCard)
<!-- badges: end -->

PostCard provides tools for accurately estimating marginal effects using
GLMs, including increasing precision using prognostic covariate
adjustment (see introductory examples below, more details in `rctglm()`
and )

**Marginal effects** are causal effects of the form $r(\Psi_1, \Psi_0)$,
where $\Psi_a=\mathbb{E}[Y(a)]$ are population mean outcomes under
exposure $a=0, 1$, respectively. These are sometimes referred to as
*counterfactual means*.

The package uses **plug-in estimation for robust estimation of any
marginal effect estimand** as well as **influence functions for robust
estimation of the variance of the estimand** (Rosenblum, M. and M. J.
van der Laan, 2010: Simple, efficient estimators of treatment effects in
randomized trials using generalized linear models to leverage baseline
variables. The International Journal of Biostatistics, 6, no. 1).

## Installation

You can install the development version of PostCard from
[GitHub](https://github.com/) with:

``` r
pak::pak("NNpackages/PostCard")
```

Setup-chunk to load the package, set a seed and turn off verbosity for
the rendering of the README.

``` r
library(PostCard)
withr::local_seed(1395878)
withr::local_options(list(PostCard.verbose = 0))
```

# Plug-in estimation of marginal effects and variance estimation using influence functions

## Simulating data for exploratory analyses

First, we simulate some data to be able to enable showcasing of the
functionalities. For this we use the `glm_data()` function from the
package, where the user can specify an expression alongside variables
and a family of the response to then simulate a response from a GLM with
linear predictor given by the expression provided.

``` r
n <- 1000
b0 <- 1
b1 <- 3
b2 <- 2

# Simulate data with a non-linear effect
dat_treat <- glm_data(
  Y ~ b0+b1*sin(W)^2+b2*A,
  W = runif(n, min = -2, max = 2),
  A = rbinom(n, 1, prob = 1/2),
  family = gaussian() # Default value
)
```

### Fitting `rctglm()` without prognostic covariate adjustment

The `rctglm()` function estimates any specified estimand using plug-in
estimation for randomised clinical trials and estimates the variance
using the influence function of the marginal effect estimand.

The interface of `rctglm()` is similar to that of the `stats::glm()`
function but with an added mandatory specification of

- The randomisation variable in data, usually being the (name of) the
  treatment variable
- The randomisation ratio - the probability of being allocated to group
  1 (rather than 0)
  - As a default, a ratio of 1’s in data is used
- An estimand function
  - As a default, the function takes the average treatment effect (ATE)
    as the estimand

Thus, we can estimate the ATE by simply writing the below:

> Note that as a default, `verbose = 2`, meaning that information about
> the algorithm is printed to the console. However, here we suppress
> this behavior. See more in `vignette("model-fit")`.

``` r
ate <- rctglm(formula = Y ~ A * W,
              exposure_indicator = A,
              exposure_prob = 1/2,
              data = dat_treat,
              family = "gaussian") # Default value
```

This creates an `rctglm` object which prints as

``` r
ate
#> 
#> Object of class rctglm 
#> 
#> Call:  rctglm(formula = Y ~ A * W, exposure_indicator = A, exposure_prob = 1/2, 
#>     data = dat_treat, family = "gaussian")
#> 
#> Counterfactual control mean (psi_0=E[Y|X, A=0]) estimate: 2.776
#> Counterfactual control mean (psi_1=E[Y|X, A=1]) estimate: 4.867
#> Warning in body(fun): argument is not a function
#> Estimand function r: NULL
#> Estimand (r(psi_1, psi_0)) estimate (SE): 2.091 (0.09209)
```

The structure of such an `rctglm` object is broken down in the
**`Value`** section of the documentation in `rctglm()`.

Methods available are `estimand` (or the shorthand `est`) which prints a
`data.frame` with and estimate of the estimand and its standard error. A
method for `coef` is also available to extract coefficients from the
underlying `glm` fit.

``` r
est(ate)
#>   Estimate Std. Error
#> 1 2.091095 0.09209306
```

See more info in the documentation page `rctglm_methods()`.

### Using prognostic covariate adjustment

The `rctglm_with_prognosticscore()` function uses the
`fit_best_learner()` function to fit a prognostic model to historical
data and then uses the prognostic model to predict

for all observations in the current data set. These *prognostic scores*
are then used as a covariate in the GLM when running `rctglm()`.

Allowing the use of complex non-linear models to create such a
prognostic score allows utilising information from potentially many
variables, “catching” non-linear relationships and then using all this
information in the GLM model using a single covariate adjustment.

We simulate some historical data to showcase the use of this function as
well:

``` r
dat_notreat <- glm_data(
  Y ~ b0+b1*sin(W)^2,
  W = runif(n, min = -2, max = 2),
  family = gaussian # Default value
)
```

The call to `rctglm_with_prognosticscore()` is the same as to `rctglm()`
but with an added specification of

- (Historical) data to fit the prognostic model using
  `fit_best_learner()`
- A formula used when fitting the prognostic model
  - Default uses all covariates in the data.
- (Optionally) number folds in cross validation and a list of learners
  for fitting the best learner

Thus, a simple call which estimates the average treatment effect,
adjusting for a prognostic score, is seen below:

``` r
ate_prog <- rctglm_with_prognosticscore(
  formula = Y ~ A * W,
  exposure_indicator = A,
  exposure_prob = 1/2,
  data = dat_treat,
  family = gaussian(link = "identity"), # Default value
  data_hist = dat_notreat)
```

Quick results of the fit can be seen by printing the object:

``` r
ate_prog
#> 
#> Object of class rctglm_prog 
#> 
#> Call:  rctglm_with_prognosticscore(formula = Y ~ A * W, exposure_indicator = A, 
#>     exposure_prob = 1/2, data = dat_treat, family = gaussian(link = "identity"), 
#>     data_hist = dat_notreat)
#> 
#> Counterfactual control mean (psi_0=E[Y|X, A=0]) estimate: 2.827
#> Counterfactual control mean (psi_1=E[Y|X, A=1]) estimate: 4.821
#> Warning in body(fun): argument is not a function
#> Estimand function r: NULL
#> Estimand (r(psi_1, psi_0)) estimate (SE): 1.994 (0.06405)
```

It’s evident that in this case where there is a non-linear relationship
between the covariate we observe and the response, adjusting for the
prognostic score reduces the standard error of our estimand
approximation by quite a bit.

#### Investigating the prognostic model

Information on the prognostic model is available in the list element
`prognostic_info`, which the method `prog()` can be used to extract. A
breakdown of what this list includes, see the **`Value`** section of the
`rctglm_with_prognosticscore()` documentation.
