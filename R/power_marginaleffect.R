#' Power approximation for estimating marginal effects in GLMs
#'
#' The functions implements the algorithm for power estimation described in
#' [Powering RCTs for marginal effects with GLMs using prognostic score adjustment](https://arxiv.org/abs/2503.22284)
#' by Højbjerre-Frandsen et. al (2025). See a bit more context in details and
#' all details in reference.
#'
#' @inheritParams rctglm
#' @inheritParams power_linear
#' @param response the response variable from comparator participants
#' @param predictions predictions of the response
#' @param target_effect a `numeric` minimum effect size that we should be able to detect.
#' See more in details.
#' @param exposure_prob a `numeric` with the probability of being in
#' "group 1" (rather than group 0). See more in details.
#' @param var1 a `numeric` variance of the potential outcome corresponding to group 1,
#' or a `function` with a single argument meant to obtain `var1` as a tranformation
#' of the variance of the potential outcome corresponding to group 0.
#' See more in details.
#' @param kappa1_squared a `numeric` mean-squared error from predicting potential
#' outcome corresponding to group 1, or a `function` with a single arguments meant
#' to obtain `kappa1_squared` as a transformation of the MSE in group 0.
#' See more in details.
#' @param inv_estimand_fun (optional) a `function` with arguments `psi0` and
#' `target_effect`, so `estimand_fun(psi1 = y, psi0 = x) = z` and
#' `inv_estimand_fun(psi0 = x, target_effect = z) = y` for all x, y, z.
#' If left as `NULL`, the inverse will be found automatically
#' @param margin a `numeric` superiority margin. As a default, the `estimand_fun`
#' is evaluated with the same counterfactual means `psi1` and `psi0`, corresponding
#' to a superiority margin assumin no difference (fx. 0 for ATE and 1 for rate ratio).
#' @param tolerance passed to [all.equal] when comparing calculated `target_effect`
#' from derivations and given `target_effect`.
#' @param ... arguments passed to `[stats::uniroot]`, which is used to find the
#' inverse of `estimand_fun`, when `inv_estimand_fun` is `NULL`.
#'
#' @returns a `numeric` with the estimated power.
#'
#' @details
#' The reference in the description shows in its "Prospective power" section a
#' derivation of a variance bound
#' \deqn{
#' v_\infty^{\uparrow 2} = r_0'^{\, 2}\sigma_0^2+
#' r_1'^{\, 2}\sigma_1^2+
#' \pi_0\pi_1\left(\frac{|r_0'|\kappa_0}{\pi_0} +
#' \frac{|r_1'|\kappa_1}{\pi_1} \right)^2
#' }
#'
#' where \eqn{r_a'} is the derivative of the `estimand_fun` with respect to
#' \eqn{\Psi_a}, \eqn{\sigma_a^2} is the variance of the potential outcome corresponding to
#' group \eqn{a}, \eqn{\pi_a} is the probablity of being assigned to group \eqn{a},
#' and \eqn{\kappa_a} is the expected mean-squared error when predicting the
#' potential outcome corresponding to group \eqn{a}.
#'
#' The variance bound is then used for calculating a lower bound of the power using
#' the distributions corresponding to the null and alternative hypotheses
#' \eqn{\mathcal{H}_0: \hat{\Psi} \sim F_0 = \mathcal{N}(\Delta ,v_\infty^{\uparrow 2} / n)}
#' and
#' \eqn{\mathcal{H}_1: \hat{\Psi} \sim F_1 = \mathcal{N}(\Psi,v_\infty^{\uparrow 2} / n)}.
#' See more details in the reference.
#'
#' ## Relating arguments to formulas
#' - `response`: Used to obtain both \eqn{\sigma_0^2} (by taking the sample
#' variance of the response) and \eqn{\kappa_0}.
#' - `predictions`: Used when calculating the MSE \eqn{\kappa_0}.
#' - `var1`: \eqn{\sigma_1^2}. As a default, chosen to be the same as
#' \eqn{\sigma_0^2}. Can specify differently through this argument fx. by
#'    - Inflating or deflating the value of \eqn{\sigma_0^2} by a scalar according
#' to prior beliefs. Fx. specify `var1 = function(x) 1.2 * x` to inflate
#' \eqn{\sigma_0^2} by 1.2.
#'    - If historical data is available for group 1, an estimate of the variance
#' from that data can be provided here.
#' - `kappa1_squared`: \eqn{\kappa_1}. Same as for `var1`, default is to assume
#' the same value as `kappa0_squared`, which is found by using the `response`
#' and `predictions` vectors. Adjust the value according to prior beliefs if
#' relevant.
#' - `target_effect`: \eqn{\Psi}.
#' - `exposure_prob`: \eqn{\pi_1}
#'
#' @examples
#' # Generate a data set to use as an example
#' n <- 100
#' exposure_prob <- .5
#'
#' dat_gaus <- glm_data(Y ~ 1+2*X1-X2+3*A+1.6*A*X2,
#'                 X1 = rnorm(n),
#'                 X2 = rgamma(n, shape = 2),
#'                 A = rbinom(n, size = 1, prob = exposure_prob),
#'                 family = gaussian())
#'
#' # ---------------------------------------------------------------------------
#' # Obtain out-of-sample (OOS) prediction using glm model
#' # ---------------------------------------------------------------------------
#' gaus1 <- dat_gaus[1:(n/2), ]
#' gaus2 <- dat_gaus[(n/2+1):n, ]
#'
#' glm1 <- glm(Y ~ X1 + X2 + A, data = gaus1)
#' glm2 <- glm(Y ~ X1 + X2 + A, data = gaus2)
#' preds_glm1 <- predict(glm2, newdata = gaus1, type = "response")
#' preds_glm2 <- predict(glm1, newdata = gaus2, type = "response")
#' preds_glm <- c(preds_glm1, preds_glm2)
#'
#' # Obtain power
#' power_marginaleffect(
#'   response = dat_gaus$Y,
#'   predictions = preds_glm,
#'   target_effect = 2,
#'   exposure_prob = exposure_prob
#' )
#'
#' # ---------------------------------------------------------------------------
#' # Get OOS predictions using discrete super learner and adjust variance
#' # ---------------------------------------------------------------------------
#' lrnr1 <- fit_best_learner(list(mod = Y ~ X1 + X2 + A), data = gaus1)
#' lrnr2 <- fit_best_learner(list(mod = Y ~ X1 + X2 + A), data = gaus2)
#' preds_lrnr1 <- dplyr::pull(predict(lrnr2, new_data = gaus1))
#' preds_lrnr2 <- dplyr::pull(predict(lrnr1, new_data = gaus2))
#' preds_lrnr <- c(preds_lrnr1, preds_lrnr2)
#'
#' # Estimate the power AND adjust the assumed variance in the "unknown"
#' # group 1 to be 20% larger than in group 0
#' power_marginaleffect(
#'   response = dat_gaus$Y,
#'   predictions = preds_lrnr,
#'   target_effect = 2,
#'   exposure_prob = exposure_prob,
#'   var1 = function(var0) 1.2 * var0
#' )
#'
#' @export
power_marginaleffect <- function(
    response,
    predictions,
    target_effect,
    exposure_prob,
    var1 = NULL,
    kappa1_squared = NULL,
    estimand_fun = "ate",
    estimand_fun_deriv0 = NULL, estimand_fun_deriv1 = NULL,
    inv_estimand_fun = NULL,
    margin = estimand_fun(1,1),
    alpha = 0.025,
    tolerance = 1e-2,
    verbose = options::opt("verbose"),
    ...
) {
  check_exposure_prob(exposure_prob = exposure_prob)

  n_resp <- length(response)
  n_pred <- length(predictions)
  if (n_resp != n_pred)
    cli::cli_abort("`response` is of length {n_resp} while `predictions` has
                   length {n_pred}. Specify them with the same length.")

  estimand_funs <- estimand_fun_setdefault_findderivs(
    estimand_fun = estimand_fun,
    estimand_fun_deriv0 = estimand_fun_deriv0,
    estimand_fun_deriv1 = estimand_fun_deriv1,
    verbose = verbose
  )
  estimand_fun <- estimand_funs$f
  psi0 <- mean(response)
  psi1 <- derive_check_psi1(
    psi0 = psi0,
    target_effect = target_effect,
    estimand_fun = estimand_fun,
    inv_estimand_fun = inv_estimand_fun,
    tolerance = tolerance
  )
  tryCatch(
    force(margin),
    error = function(e) {
      bullets <- c("The value of `margin` produced an error.",
                   i = "Specify `margin` explicitly as a `numeric`.")
      cli::cli_abort(bullets)
    }
  )
  d0_val <- estimand_funs$d0(psi1 = psi1, psi0 = psi0)
  d1_val <- estimand_funs$d1(psi1 = psi1, psi0 = psi0)

  var0 <- var(response)
  var1 <- setdefault_transform(val = var1, default = var0)

  kappa0_squared <- 1/n_resp * sum((response - predictions)^2)
  kappa1_squared <- setdefault_transform(val = kappa1_squared, default = kappa0_squared)

  v_bound <- var_bound_marginaleffect(
    var1 = var1, var0 = var0,
    d0 = d0_val, d1 = d1_val,
    kappa0_squared = kappa0_squared,
    kappa1_squared = kappa1_squared,
    exposure_prob = exposure_prob
  )

  sd <- sqrt(v_bound / n_resp)
  f0 <- qnorm(1 - alpha, mean = margin, sd = sd)
  f1 <- pnorm(f0, mean = target_effect, sd = sd)
  1 - f1
}

var_bound_marginaleffect <- function(
    var1, var0,
    d0, d1,
    kappa0_squared, kappa1_squared,
    exposure_prob) {
  d0^2 * var0 + d1^2 * var1 +
    (1 - exposure_prob) * exposure_prob *
    ((abs(d0) * sqrt(kappa0_squared) / (1 - exposure_prob) +
        abs(d1) * sqrt(kappa1_squared) / exposure_prob)^2)
}

setdefault_transform <- function(val, default) {
  if (is.null(val))
    return(default)
  if (is.function(val))
    val <- return(val(default))
  return(val)
}


# Check if lower_upper default of -1e3 and 1e3 gives an error and if so set
# lower to be positive
check_lower_upper <- function(f = NULL, f_arg = NULL, lower = NULL, upper = NULL,
                              default_lu_times = 0, default_lu_scale = 1e4,
                              ...) {
  lower_upper_notgiven <- is.null(lower) || is.null(upper)
  if (lower_upper_notgiven) {
    onegiven <- !is.null(lower) || !is.null(upper)
    if (onegiven)
      cli::cli_abort("You have specified `lower` or `upper` but not the other. Please specify both or none to get the default.")

    lu_range <- (1 + default_lu_times) * default_lu_scale
    lower <- -lu_range
    upper <- lu_range

    inv_f <- inverse(f, lower = lower, upper = upper, ...)
    lower_upper_bad <- tryCatch(
      inv_f(f_arg),
      error = function(e) TRUE,
      warning = function(w) TRUE
    )

    if (isTRUE(lower_upper_bad))
      lower <- 1 / lu_range
  }

  list(lower = lower, upper = upper)
}

# Find inverse function using uniroot
inverse <- function(f, lower, upper, ...){
  function(y){
    stats::uniroot(function(x){f(x) - y}, lower = lower, upper = upper, ...)$root
  }
}

# Check lower and upper, find the inverse and evaluate it
inverse_val <- function(
    f, f_arg, lower = NULL, upper = NULL,
    ...){
  args <- c(as.list(environment()), list(...))
  lower_upper <- do.call(check_lower_upper, args)
  inverse <- do.call(inverse, c(list(f = f), lower_upper))

  return(inverse(f_arg))
}

# Derive psi1 from psi0 and a target_effect using the inverse of the estimand_fun
# If inv_estimand_fun is not given, find it using functionalities from above
derive_check_psi1 <- function(
    psi0,
    target_effect,
    estimand_fun,
    inv_estimand_fun = NULL,
    tolerance = 1e-2,
    ...) {
  if (!is.null(inv_estimand_fun))
    psi1 <- inv_estimand_fun(psi0 = psi0, target_effect = target_effect)
  else {
    args_to_uniroot <- list(...)
    estimand_fun1 <- function(psi1) {
      estimand_fun(psi1 = psi1, psi0 = psi0)
    }

    psi1 <- do.call(
      inverse_val,
      c(list(f = estimand_fun1, f_arg = target_effect, default_lu_times = psi0),
        args_to_uniroot)
    )
  }

  calc_target_effect <- estimand_fun(psi0 = psi0, psi1 = psi1)
  diff_target_effect <- calc_target_effect - target_effect
  psi1_correct <- all.equal(diff_target_effect, 0, tolerance = tolerance)
  if (!isTRUE(psi1_correct)) {
    bullets <- c(
      "The inverse of the estimand fun `inv_estimand_fun` did not
      produce a good result. The calculated `target_effect` varies from the
      assumed/given `target_effect` by {.arg {diff_target_effect}}. Either",
      i = "specify `inv_estimand_fun` manually, or",
      i = "specify arguments `lower` and `upper` passed to `uniroot` that finds the inverse, or",
      i = "increase the value of the `tolerance` argument."
    )
    cli::cli_warn(bullets)
  }

  return(psi1)
}
