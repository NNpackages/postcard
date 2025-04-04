# predict_counterfactual_mean
test_that("`predict_counterfactual_mean` predicts correctly", {
  treat_diff <- 10
  dat <- data.frame(
    Y = 1:(2*treat_diff),
    X = rep(1:treat_diff, 2),
    A = c(rep(0, treat_diff), rep(1, treat_diff)
    )
  )
  mod <- glm(Y ~ X + A, data = dat)
  pred0 <- predict_counterfactual_mean(
    model = mod,
    exposure_indicator_name = "A",
    group_val = 0)

  pred1 <- predict_counterfactual_mean(
    model = mod,
    exposure_indicator_name = "A",
    group_val = 1)

  expect_equal(pred0, pred1 - treat_diff)
})

test_that("`predict_counterfactual_mean` gives error when `exposure_indicator_name` not in model or data", {
  dat <- data.frame(
    Y = 1:10,
    X = 1:10,
    A = c(rep(0, 5), rep(1, 5)
    )
  )
  mod <- glm(Y ~ X + A, data = dat)

  expect_error(
    predict_counterfactual_mean(
      model = mod,
      exposure_indicator_name = "test",
      group_val = 0),
    regexp = "is not in"
  )
})

test_that("`predict_counterfactual_mean` works with and without data specification", {
  dat_fit <- data.frame(
    Y = 1:10,
    X = 1:10,
    A = c(rep(0, 5), rep(1, 5)
    )
  )
  mod <- glm(Y ~ X + A, data = dat_fit)

  pred_nodatspec <- predict_counterfactual_mean(
    model = mod,
    exposure_indicator_name = "A",
    group_val = 0)
  pred_datspec <- predict_counterfactual_mean(
    model = mod,
    exposure_indicator_name = "A",
    group_val = 0,
    data = dat_fit)

  expect_equal(pred_nodatspec, pred_datspec)

  dat_pred <- data.frame(
    Y = 1:11,
    X = -5:5
  )
  pred_newdata <- predict_counterfactual_mean(
    model = mod,
    exposure_indicator_name = "A",
    group_val = 0,
    data = dat_pred)

  expect_type(pred_newdata, "double")
})

# predict_counterfactual_means
test_that("`predict_counterfactual_mean` predicts correctly", {
  treat_diff <- 10
  dat <- data.frame(
    Y = 1:(2*treat_diff),
    X = rep(1:treat_diff, 2),
    A = c(rep(0, treat_diff), rep(1, treat_diff)
    )
  )
  mod <- glm(Y ~ X + A, data = dat)
  preds <- predict_counterfactual_means(
    model = mod,
    exposure_indicator_name = "A")

  expect_s3_class(preds, "data.frame")
  expect_equal(preds$psi0, preds$psi1 - treat_diff)
})

# default_estimand_funs
test_that("`default_estimand_funs` switches correctly", {
  default_ate <- default_estimand_funs()
  specify_ate <- default_estimand_funs("ate")
  specify_rr <- default_estimand_funs("rate_ratio")

  expect_identical(default_ate, specify_ate)
  expect_equal(default_ate(2, 1), 1)
  expect_equal(specify_rr(2, .5), 4)
})

test_that("`default_estimand_funs` error when giving non-legal default", {
  expect_error(default_estimand_funs("test"),
               "should be one of")
})

# oos_fitted.values_counterfactual
test_that("`oos_fitted.values_counterfactual` snapshot test", {
  dat <- data.frame(
    Y = 1:10,
    X = 1:10,
    A = c(rep(0, 5), rep(1, 5)
    )
  )

  args_glm <- list(
    formula = formula(Y ~ X + A)
  )

  oos <- oos_fitted.values_counterfactual(
    data = dat,
    exposure_indicator_name = "A",
    full_model.args_glm = args_glm
  )
  expect_named(oos, c("psi0", "psi1", "rowname"))
  expect_s3_class(oos, "data.frame")
  expect_snapshot(oos)
})

# extract_train_test
test_that("`extract_train_test` returns list of train and test data", {
  dat <- data.frame(
    Y = 1:10,
    X = 1:10,
    A = c(rep(0, 5), rep(1, 5)
    )
  )

  withr::local_seed(42)
  folds <- rsample::vfold_cv(dat, v = 2)
  single_fold <- folds$splits[[1]]
  list_of_train_test <- extract_train_test(single_fold)
  expect_type(list_of_train_test, "list")
  expect_named(list_of_train_test, c("train", "test"))
  expect_snapshot(list_of_train_test)
})
