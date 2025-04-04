#' Power and sample size estimation for linear models
#'
#' @name power_linear
#'
#' @details
#'
#' This details section provides information about relation between arguments to
#' functions and the formulas described in sections below for each power
#' approximation formula.
#'
#' Note that all entities that carry the same name as an argument and in the formula
#' will not be mentioned below, as they are obviously linked (n, r, alpha)
#'
#' - `ate`: \eqn{\beta_1-\beta_0}
#' - `margin`: \eqn{\Delta_s}
#' - `variance`: \eqn{\widehat{\sigma}^2(1-\widehat{R}^2)}
#'
#' ## Finding the `variance` to use for approximation
#'
#' The `variance_ancova` function estimates \eqn{\sigma^2(1-R^2)} in data and
#' returns it as a `numeric` that can be passed directly as the `variance`
#' in `power_gs`. Corresponds to estimating the power from using an `lm` with
#' the same `formula` as specified in `variance_ancova`.
#'
#' The user can estimate the `variance` any way they see fit.
#'
#' @return
#' All functions return a `numeric`. `variance_ancova` returns a `numeric` with
#' a variance estimated from data to used for power estimation and sample size
#' estimation. `power_xx` and `samplesize_xx` functions return a `numeric` with
#' the power or sample size approximation.
#'
NULL

#' @rdname power_linear
#'
#' @description
#' `variance_ancova` provides a convenient function for estimating a
#' variance to use for power and sample size approximation.
#'
#' @inheritParams stats::model.frame
#' @param formula an object of class "formula" (or one that can be coerced to that class):
#' a symbolic description used in [stats::model.frame()] to create a `data.frame` with
#' response and covariates. This data.frame is used to estimate the \eqn{R^2}, which is
#' then used to find the variance. See more in details.
#' @param inflation a `numeric` to multiply the marginal variance of the response by.
#' Default is `1` which estimates the variance directly from data. Use values above `1`
#' to obtain a more conservative estimate of the marginal response variance.
#' @param deflation a `numeric` to multiply the \eqn{R^2} by.
#' Default is `1` which means the estimate of \eqn{R^2} is unchanged. Use values
#' below `1` to obtain a more conservative estimate of the coefficient of determination.
#' See details about how \eqn{R^2} related to the estimation.
#'
#' @export
variance_ancova <- function(formula, data, inflation = 1, deflation = 1) {
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  Y <- stats::model.response(mf)
  var_Y <- Y %>%
    var()
  var_Y_inf <- var_Y * inflation

  mt <- attr(mf, "terms")
  X <- stats::model.matrix(mt, mf)
  X_has_intercept <- attr(mt, "intercept")
  if (X_has_intercept)
    X <- X[, -1]
  if (ncol(as.matrix(X)) == 0)
    return(var_Y_inf)

  Sigma_X.I <- X %>%
    var() %>%
    chol() %>%
    chol2inv()
  sigma_XY <- stats::cov(X, Y)

  R2 <- as.numeric((t(sigma_XY) %*% Sigma_X.I %*% sigma_XY) / var_Y_inf)
  R2_def <- R2 * deflation

  var_bound <- var_Y_inf * (1 - R2_def)
  return(var_bound)
}

#' @rdname power_linear
#'
#' @description
#' The `power_gs` and `samplesize_gs` functions calculate the Guenther-Schouten
#' power approximation for ANOVA or ANCOVA.
#' The approximation is based in (Guenther WC. Sample Size Formulas for Normal Theory T Tests.
#' The American Statistician. 1981;35(4):243–244) and (Schouten HJA. Sample size formula with
#' a continuous outcome for unequal group sizes and unequal variances. Statistics in Medicine.
#' 1999;18(1):87–91).
#'
#' @param n               a `numeric` with number of participants in total.
#' From this number of participants in the treatment group is \eqn{n1=(r/(1+r))n}
#' and the control group is \eqn{n1=(1/(1+r))n}.
#' @param r               a `numeric` allocation ratio \eqn{r=n1/n0}. For one-to-one randomisation `r=1`.
#' @param variance  a `numeric` variance to use for the approximation. See more details
#' in documentation sections of each power approximating function.
#' @param ate             a `numeric` minimum effect size that we should be able to detect.
#' @param margin          a `numeric` superiority margin
#' (for non-inferiority margin, a negative value can be provided).
#' @param alpha           a `numeric` significance level. Due to regulatory guidelines when
#' using a one-sided test, half the specified significance level is used.
#' Thus, for standard significance level of 5%, the default is `alpha = 0.025`.
#'
#' @section Guenther-Schouten power approximation:
#' The estimation formula in the case of an ANCOVA model with multiple covariate adjustement is (see description for reference):
#'
#' \deqn{
#' n=\frac{(1+r)^2}{r}\frac{(z_{1-\alpha}+z_{1-\beta})^2\widehat{\sigma}^2(1-\widehat{R}^2)}{(\beta_1-\beta_0-\Delta_s)^2}+\frac{(z_{1-\alpha})^2}{2}
#' }
#'
#' where \eqn{\widehat{R}^2\coloneqq \frac{\widehat{\sigma}_{XY}^\top \widehat{\Sigma}_X^{-1}\widehat{\sigma}_{XY}}{\widehat{\sigma}^2}},
#' we denote by \eqn{\widehat{\sigma^2}} an estimate of the variance of the outcome,
#' \eqn{\widehat{\Sigma_X}} and estimate of the covariance matrix of the
#' covariates, and \eqn{\widehat{\sigma_{XY}}} a \eqn{p}-dimensional column vector consisting of
#' an estimate of the covariance
#' between the outcome variable and each covariate.
#' In the univariate case \eqn{R^2} is replaced by \eqn{\rho^2}
#'
#' @examples
#' # Generate a data set to use as an example
#' dat_gaus <- glm_data(Y ~ 1+2*X1-X2+3*A,
#'                 X1 = rnorm(100),
#'                 X2 = rgamma(100, shape = 2),
#'                 A = rbinom(100, size = 1, prob = 0.5),
#'                 family = gaussian())
#'
#' # Approximate the power using no adjustment covariates
#' va_nocov <- var(dat_gaus$Y)
#' power_gs(n = 200, variance = va_nocov, ate = 1)
#'
#' # Approximate the power with a model adjusting for both variables in the
#' # data generating process
#'
#' ## First estimate the variance sigma^2 * (1-R^2)
#' va_cov <- variance_ancova(Y ~ X1 + X2 + A, dat_gaus)
#' ## Then estimate the power using this variance
#' power_gs(n = 100, variance = va_cov, ate = 1.8, margin = 1, r = 2)
#'
#' # Approximate the sample size needed to obtain 90% power with same model as
#' # above
#' samplesize_gs(
#'   variance = va_cov, ate = 1.8, power = 0.9, margin = 1, r = 2
#' )
#'
#' @export
power_gs <- function(variance,
                     ate,
                     n,
                     r = 1,
                     margin = 0,
                     alpha = 0.025) {

  power <- stats::pnorm(
    sqrt(
      r/(1 + r)^2 *
        (ate - margin)^2/variance *
        (n - stats::qnorm(1 - alpha)^2/2)
    ) - stats::qnorm(1 - alpha)
  )
  return(power)
}


#' @rdname power_linear
#'
#' @param power a `numeric` giving the desired power when calculating the sample size
#'
#' @export
samplesize_gs <- function(variance,
                          ate,
                          r = 1,
                          margin = 0,
                          power = 0.9,
                          alpha = 0.025) {
  samplesize <- (1 + r)^2 / r *
    (stats::qnorm(1 - alpha) + stats::qnorm(power))^2 * variance /
    (ate - margin)^2 +
    stats::qnorm(1 - alpha)^2/2
  return(samplesize)
}

#' @rdname power_linear
#'
#' @description
#' The function `power_nc` calculates the power for ANOVA or ANCOVA based on the
#' non-centrality parameter and the exact t-distributions.
#'
#' See more details about each funtion in `Details` and in sections after `Value`.
#'
#' @param df a `numeric` degrees of freedom to use in the t-distribution.
#'
#' @section Power approximation using non-centrality parameter:
#' The prospective power estimations are based on (Kieser M. Methods and Applications of Sample Size Calculation and Recalculation in Clinical Trials. Springer; 2020).
#' The ANOVA power is calculated based on the non-centrality parameter given as
#'
#' \deqn{nc =\sqrt{\frac{r}{(1+r)^2}\cdot n}\cdot\frac{\beta_1-\beta_0-\Delta_s}{\sigma},}
#'
#' where we denote by \eqn{\sigma^2} the variance of the outcome, such that the power can be estimated as
#'
#' \deqn{1-\beta = 1 - F_{t,n-2,nc}\left(F_{t, n-2, 0}^{-1}(1-\alpha)\right).}
#'
#' The power of ANCOVA with univariate covariate adjustment and no interaction is calculated based on the non-centrality parameter given as
#'
#' \deqn{nc =\sqrt{\frac{rn}{(1+r)^2}}\frac{\beta_1-\beta_0-\Delta_s}{\sigma\sqrt{1-\rho^2}},}
#'
#' such that the power can be estimated as
#'
#' \deqn{1-\beta = 1 - F_{t,n-3,nc}\left(F_{t, n-3, 0}^{-1}(1-\alpha)\right).}
#'
#' The power of ANCOVA with either univariate covariate adjustment and interaction or multiple covariate adjustement with or without interaction is calculated based on the non-centrality parameter given as
#'
#' \deqn{nc =\frac{\beta_1-\beta_0-\Delta_s}{\sqrt{\left(\frac{1}{n_1}+\frac{1}{n_0} + X_d^\top\left((n-2)\Sigma_X\right)^{-1}X_d \right)\sigma^2\left(1-\widehat{R}^2\right)}}.}
#'
#' where \eqn{X_d \coloneqq \left(\overline{X}_1^1-\overline{X}_0^1, \ldots, \overline{X}_1^p-\overline{X}_0^p\right)^\top}, \eqn{\widehat{R}^2\coloneqq \frac{\widehat{\sigma}_{XY}^\top \widehat{\Sigma}_X^{-1}\widehat{\sigma}_{XY}}{\widehat{\sigma}^2}},
#' we denote by \eqn{\widehat{\sigma^2}} an estimate of the variance of the outcome,
#' \eqn{\widehat{\Sigma_X}} and estimate of the covariance matrix of the
#' covariates, and \eqn{\widehat{\sigma_{XY}}} a \eqn{p}-dimensional column vector consisting of
#' an estimate of the covariance
#' between the outcome variable and each covariate.
#' Since we are in the case of randomized trials the expected difference between the covariate
#' values between the to groups is 0. Furthermore, the elements of \eqn{\Sigma_X^{-1}} will be small, unless the variances are close to 0, or the covariates exhibit strong linear dependencies, so that the correlations are close to 1.
#' These scenarios are excluded since they could lead to potentially serious problems regarding inference either way. These arguments are used by Zimmermann et. al
#' (Zimmermann G, Kieser M, Bathke AC. Sample Size Calculation and Blinded Recalculation for Analysis of Covariance Models with Multiple Random Covariates. Journal of Biopharmaceutical Statistics. 2020;30(1):143–159.) to approximate
#' the non-centrality parameter as in the univariate case where \eqn{\rho^2} is replaced by \eqn{R^2}.
#'
#' Then the power for ANCOVA with `d` degrees of freedom can be estimated as
#'
#' \deqn{1-\beta = 1 - F_{t,d,nc}\left(F_{t, d,0), 0}^{-1}(1-\alpha)\right).}
#'
#' @export
#'
#' @examples
#' # No adjustment covariates
#' power_nc(n = 200, variance = va_nocov, df = 199, ate = 1)
#' # Adjusting for all covariates in data generating process
#' power_nc(n = 200, variance = va_cov, df = 196, ate = 1.8, margin = 1, r = 2)
#'
power_nc <- function(variance,
                     df,
                     ate,
                     n,
                     r = 1,
                     margin = 0,
                     alpha = 0.025){

  nc <- sqrt(r/(1 + r)^2*(n)) * (ate - margin)/sqrt(variance)
  power <- 1 - stats::pt(q = stats::qt(1 - alpha, df = df), df = df, ncp = nc)
  return(power)
}
