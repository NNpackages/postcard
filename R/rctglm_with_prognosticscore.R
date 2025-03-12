#' Use prognostic covariate adjustment when fitting an [rctglm]
#'
#' The procedure uses [fit_best_learner] to fit a prognostic model to historical data and uses
#' the model to produce counterfactual predictions as a prognostic score that is then adjusted
#' for as a covariate in the [rctglm] procedure.
#'
#' @inheritParams rctglm
#' @inheritParams fit_best_learner
#'
#' @param data_hist a `data.frame` with historical data on which to fit a prognostic model
#' @param prog_formula a `character` or `numeric` with the formula for fitting the prognostic
#' model on the historical data `data_hist`. Default models the response (assumed same as in
#' `formula`) using all columns in the `data_hist` data
#' @param cv_folds_prog a `numeric` with the number of cross-validation folds used when fitting and
#' evaluating models
#'
#' @details
#' More details on prognostic models and scores being predictions of counterfactual means
#' in control group.
#'
#' @return `rctglm_with_prognosticscore` returns an object of class `rctglm_prog`,
#' which inherits from [rctglm].
#'
#' Due to inheritance, methods described in the **Value** section of [rctglm] and
#' [rctglm_methods] also work for extracting information from this object.
#'
#' An `rctglm_prog` object is a list with the same components as an [rctglm] object,
#' but with an additional list element of:
#' - `prognostic_info`: List with information about the fitted prognostic model
#' on historical data. It has components:
#'    - `formula`: The `formula` with symbolic description of how the response
#'    is modelled as function of covariates in the models
#'    - `model_fit`: A trained `workflow` - the result of [fit_best_learner]
#'    - `learners`: A `list` of learners used for the discrete super learner
#'    - `cv_folds`: The amount of folds used for cross validation
#'    - `data`: The historical data used for cross validation when fitting and
#'    testing models
#' @export
#'
#' @examples
#' # Generate some data
#' n <- 100
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' W1 <- runif(n, min = -2, max = 2)
#'
#' dat_treat <- glm_data(
#'   b0+b1*abs(sin(W1))+b2*A,
#'   W1 = W1,
#'   A = rbinom (n, 1, .5)
#' )
#'
#' dat_notreat <- glm_data(
#'   b0+b1*abs(sin(W1)),
#'   W1 = W1
#' )
#'
#' ate <- rctglm_with_prognosticscore(
#'   formula = Y ~ .,
#'   group_indicator = A,
#'   data = dat_treat,
#'   family = gaussian(),
#'   estimand_fun = "ate",
#'   data_hist = dat_notreat)
#'
#' # Pull information on estimand
#' estimand(ate)
rctglm_with_prognosticscore <- function(
    formula,
    family,
    data,
    group_indicator,
    group_allocation_prob = NULL,
    estimand_fun = "ate",
    estimand_fun_deriv0 = NULL, estimand_fun_deriv1 = NULL,
    cv_variance = TRUE,
    cv_folds_variance = 5,
    ...,
    data_hist,
    prog_formula = NULL,
    cv_folds_prog = 5,
    learners = default_learners(),
    verbose = options::opt("verbose")) {

  call <- match.call()

  group_indicator <- rlang::enquo(group_indicator)
  named_args <- as.list(environment())
  extra_glm_args <- list(...)

  formula <- check_formula(formula)

  if (verbose >= 1) cli::cli_h2("Fitting prognostic model")
  if (is.null(prog_formula)) {
    is_response_in_data(formula, data = data_hist)
    prog_formula <- formula_everything(formula)
    if (verbose >= 1)
      cli::cli_alert_info("Created formula for fitting prognostic model as: {deparse(prog_formula)}")
  } else if (is.character(prog_formula)) {
    prog_formula <- formula(prog_formula)
  }

  data_hist <- do.call(
    model.frame,
    c(
      list(
        formula = prog_formula,
        data = data_hist
      ),
      extra_glm_args
    )
  )

  lrnr_fit <- fit_best_learner(formula = prog_formula,
                               data = data_hist,
                               cv_folds = cv_folds_prog,
                               learners = learners,
                               verbose = verbose)

  lrnr_pred <- predict(lrnr_fit, data) %>%
    dplyr::pull(.data$.pred)
  data %<>%
    dplyr::mutate(prog = lrnr_pred)

  if (verbose >= 2)
    cli::cli_alert_info("Investigate trained learners and fitted model in {.var prognostic_info} list element")

  formula_with_prognosticscore <- paste0(formula_to_str(formula), " + prog")

  rctglm_with_prognosticscore <- rctglm(
    formula = formula_with_prognosticscore,
    family = family,
    data = data,
    group_indicator = group_indicator,
    group_allocation_prob = group_allocation_prob,
    estimand_fun = estimand_fun,
    estimand_fun_deriv0 = estimand_fun_deriv0,
    estimand_fun_deriv1 = estimand_fun_deriv1,
    cv_variance = cv_variance,
    cv_variance_folds = cv_variance_folds,
    verbose = verbose,
    ...
  )

  prog_info <- list(
    prognostic_info = list(
      formula = prog_formula,
      model_fit = lrnr_fit,
      learners = learners,
      cv_folds = cv_folds_prog,
      data = data_hist
    ))

  list_with_prognostic_info <- c(rctglm_with_prognosticscore, prog_info)
  list_with_prognostic_info$call <- call

  out <- structure(list_with_prognostic_info,
                   class = c("rctglm_prog",
                             class(rctglm_with_prognosticscore)
                   )
  )

  return(out)
}
