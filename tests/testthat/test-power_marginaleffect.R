test_that("`power_marginaleffect` snapshot tests", {
  withr::local_seed(42)
  exp_prob <- 1/2
  n <- 100
  dat <- glm_data(Y ~ X1 + 2*A,
                  X1 = rnorm(n),
                  A = rbinom(n, 1, exp_prob))
  mod <- glm(Y ~ X1 + A, data = dat)
  preds <- predict(mod, newdata = dat)

  pow <- power_marginaleffect(
    response = dat$Y,
    predictions = preds,
    target_effect = 2,
    exposure_prob = exp_prob
  )
  expect_snapshot(pow)

  spec_var_kappa <- power_marginaleffect(
    response = dat$Y,
    predictions = preds,
    var1 = function(var0) 2 * var0,
    kappa1_squared = 1.5,
    target_effect = 2,
    exposure_prob = exp_prob
  )
  expect_snapshot(spec_var_kappa)
})

test_that("`power_marginaleffect` gives errors", {
  Y <- rnorm(100)
  preds <- rnorm(50)

  expect_error(
    power_marginaleffect(
      response = Y,
      predictions = preds,
      target_effect = 2,
      exposure_prob = 2/3
    ),
    regexp = "Specify them with the same length"
  )

  expect_error(
    power_marginaleffect(
      margin = not_defined_value,
      response = Y,
      predictions = Y,
      target_effect = 2,
      exposure_prob = 2/3,
      verbose = 0
    ),
    regexp = "Specify `margin` explicitly as a `numeric`"
  )
})

#inverse
test_that("`inverse` works", {
  inv_sqrt <- inverse(sqrt, lower = 0.01, upper = 100)
  inv_sum <- inverse(function(x) x + 2, lower = -100, upper = 100)
  inv_div <- inverse(function(x) x / 3, lower = -100, upper = 100)

  expect_equal(inv_sqrt(2), 2^2)
  expect_equal(inv_sum(3), 1)
  expect_equal(inv_div(3), 9)
})

test_that("`inverse_val` works", {
  inv_sum_val <- inverse_val(function(x) x + 7, 1e2)
  inv_div_val <- inverse_val(function(x) x / 4, 5)
  expect_equal(inv_sum_val, 1e2-7)
  expect_equal(inv_div_val, 20)
})

# check_lower_upper
test_that("`check_lower_upper` works", {
  give_args <- check_lower_upper(lower = -100, upper = 100)
  expect_equal(give_args, list(lower = -100, upper = 100))

  negative_okay <- check_lower_upper(
    f = function(x) x / 4, f_arg = 5,
    default_lu_scale = 1e2
  )
  expect_equal(negative_okay, list(lower = -1e2, upper = 1e2))

  nonnegative_only <- check_lower_upper(
    sqrt, f_arg = 4,
    default_lu_scale = 1e5, default_lu_times = 2
  )
  expect_equal(nonnegative_only, list(lower = 1 / ((1+2)*1e5), upper = (1+2)*1e5))
})

# derive_check_psi1
test_that("`derive_check_psi1` works", {
  ate_fun <- function(psi1, psi0) psi1 - psi0
  ate_manual_inv <- derive_check_psi1(
    estimand_fun = ate_fun,
    inv_estimand_fun = function(psi0, target_effect) psi0 + target_effect,
    psi0 = 2, target_effect = 2
  )
  ate_auto_inv <- derive_check_psi1(
    estimand_fun = ate_fun,
    psi0 = 2, target_effect = 2
  )
  expect_equal(ate_manual_inv, ate_auto_inv)

  rate_ratio_fun <- function(psi1, psi0) psi1 / psi0
  rr_manual_inv <- derive_check_psi1(
    estimand_fun = rate_ratio_fun,
    inv_estimand_fun = function(psi0, target_effect) psi0 * target_effect,
    psi0 = 2, target_effect = 2
  )
  rr_auto_inv <- derive_check_psi1(
    estimand_fun = rate_ratio_fun,
    psi0 = 2, target_effect = 2
  )
  expect_equal(rr_manual_inv, rr_auto_inv)

  nonsense_fun <- function(psi1, psi0) psi1^2 * (psi0 - log(psi0))
  ns_manual_inv <- derive_check_psi1(
    estimand_fun = nonsense_fun,
    inv_estimand_fun = function(psi0, target_effect) sqrt(target_effect / (psi0 - log(psi0))),
    psi0 = 2, target_effect = 2,
    tolerance = 1e-2
  )
  ns_auto_inv <- derive_check_psi1(
    estimand_fun = nonsense_fun,
    psi0 = 2, target_effect = 2,
    tolerance = 1e-2
  )
  expect_equal(ns_manual_inv, ns_auto_inv, tolerance = 1e-2)
})

test_that("`derive_check_psi1` gives warning when inverse does not produce correct result", {
  nonsense_fun <- function(psi1, psi0) psi1^2 * (psi0 - log(psi0))
  expect_warning(
    derive_check_psi1(
      estimand_fun = nonsense_fun,
      psi0 = 2, target_effect = 2,
      tolerance = .Machine$double.eps
    ),
    regexp = "did not produce a good result"
  )
})
