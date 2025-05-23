% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_best_learner.R
\name{fit_best_learner}
\alias{fit_best_learner}
\title{Find the best learner in terms of RMSE among specified learners using cross validation}
\usage{
fit_best_learner(
  preproc,
  data,
  cv_folds = 5,
  learners = default_learners(),
  verbose = options::opt("verbose")
)
}
\arguments{
\item{preproc}{A list (preferably named) with preprocessing objects:
formulas, recipes, or \code{\link[workflows:add_variables]{workflows::workflow_variables()}}. Passed to
\code{\link[workflowsets:workflow_set]{workflowsets::workflow_set()}}.}

\item{data}{A data frame.}

\item{cv_folds}{a \code{numeric} with the number of cross-validation folds used when fitting and
evaluating models}

\item{learners}{a \code{list} (preferably named) containing named lists of elements
\code{model} and optionally \code{grid}. The \code{model} element should be a \code{parsnip}
model specification, which is passed to \link[workflowsets:workflow_set]{workflowsets::workflow_set} as the
\code{model} argument, while the \code{grid} element is passed as the \code{grid} argument
of \link[workflowsets:option_add]{workflowsets::option_add}}

\item{verbose}{\code{numeric} verbosity level. Higher values means more information is
printed in console. A value of 0 means nothing is printed to console during
execution (Defaults to \code{2}, overwritable using option 'postcard.verbose' or environment variable 'R_POSTCARD_VERBOSE')}
}
\value{
a trained \code{workflow}
}
\description{
Find the best learner in terms of RMSE among specified learners using cross validation
}
\details{
Ensure data compatibility with the learners.
}
\examples{
# Generate some synthetic 2-armed RCT data along with historical controls
n <- 100
dat_rct <- glm_data(
  Y ~ 1+2*x1+3*a,
  x1 = rnorm(n, 2),
  a = rbinom (n, 1, .5),
  family = gaussian()
)
dat_hist <- glm_data(
  Y ~ 1+2*x1,
  x1 = rnorm(n, 2),
  family = gaussian()
)

# Fit a learner to the historical control data
learners <- list(
  mars = list(
    model = parsnip::set_engine(
      parsnip::mars(
        mode = "regression", prod_degree = 3
      ),
      "earth"
    )
  )
)
fit <- fit_best_learner(
  preproc = list(mod = Y ~ .),
  data = dat_hist,
  learners = learners
)

# Use it fx. to predict the "control outcome" in the 2-armed RCT
predict(fit, new_data = dat_rct)

}
\seealso{
See \code{\link[=rctglm_with_prognosticscore]{rctglm_with_prognosticscore()}} for a function that utilises this
function to perform prognostic covariate adjustment.
}
