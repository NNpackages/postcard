% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rctglm_with_prognosticscore.R
\name{rctglm_with_prognosticscore}
\alias{rctglm_with_prognosticscore}
\title{Use prognostic covariate adjustment when fitting an \link{rctglm}}
\usage{
rctglm_with_prognosticscore(
  formula,
  exposure_indicator,
  exposure_prob,
  data,
  family = gaussian,
  estimand_fun = "ate",
  estimand_fun_deriv0 = NULL,
  estimand_fun_deriv1 = NULL,
  cv_variance = FALSE,
  cv_variance_folds = 10,
  ...,
  data_hist,
  prog_formula = NULL,
  cv_prog_folds = 5,
  learners = default_learners(),
  verbose = options::opt("verbose")
)
}
\arguments{
\item{formula}{an object of class "formula" (or one that can be coerced to that class):
a symbolic description of the model to be fitted. The details of model specification are
given under ‘Details’ in the \link{glm} documentation.}

\item{exposure_indicator}{(name of) the \emph{binary} variable in \code{data} that
identifies randomisation groups. The variable is required to be binary to
make the "orientation" of the \code{estimand_fun} clear.}

\item{exposure_prob}{a \code{numeric} with the probability of being in
"group 1" (rather than group 0) in groups defined by \code{exposure_indicator}.}

\item{data}{an optional data frame, list or environment (or object coercible
by as.data.frame to a data frame) containing the variables in the model. If
not found in data, the variables are taken from environment(formula), typically
the environment from which the function is called.}

\item{family}{a description of the error distribution and link
    function to be used in the model.  For \code{glm} this can be a
    character string naming a family function, a family function or the
    result of a call to a family function.  For \code{glm.fit} only the
    third option is supported.  (See \code{\link[stats]{family}} for details of
    family functions.)}

\item{estimand_fun}{a \code{function} with arguments \code{psi1} and \code{psi0} specifying
the estimand. Alternative, specify "ate" or "rate_ratio" as a \code{character}
to use one of the default estimand functions. See
more details in the "Estimand" section of \link{rctglm}.}

\item{estimand_fun_deriv0}{a \code{function} specifying the derivative of \code{estimand_fun} wrt. \code{psi0}. As a default
the algorithm will use symbolic differentiation to automatically find the derivative from \code{estimand_fun}}

\item{estimand_fun_deriv1}{a \code{function} specifying the derivative of \code{estimand_fun} wrt. \code{psi1}. As a default
the algorithm will use symbolic differentiation to automatically find the derivative from \code{estimand_fun}}

\item{cv_variance}{a \code{logical} determining whether to estimate the variance
using cross-validation (see details of \link{rctglm}).}

\item{cv_variance_folds}{a \code{numeric} with the number of folds to use for cross
validation if \code{cv_variance} is \code{TRUE}.}

\item{...}{Additional arguments passed to \code{\link[stats:glm]{stats::glm()}}}

\item{data_hist}{a \code{data.frame} with historical data on which to fit a prognostic model}

\item{prog_formula}{an object of class "formula" (or one that can be coerced to that class):
a symbolic description of the prognostic model to be fitted to \code{data_hist}. Passed in a
list as the \code{preproc} argument in \code{\link[=fit_best_learner]{fit_best_learner()}}. As a default,
the formula is created by modelling the response (assumed to have the same name as in
\code{formula}) using all columns in \code{data_hist}.}

\item{cv_prog_folds}{a \code{numeric} with the number of cross-validation folds used when fitting and
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
\code{rctglm_with_prognosticscore} returns an object of class \code{rctglm_prog},
which inherits from \link{rctglm}.

An \code{rctglm_prog} object is a list with the same components as an \link{rctglm} object
(see the \strong{\code{Value}} section of \link{rctglm} for a breakdown of the structure),
but with an additional list element of:
\itemize{
\item \code{prognostic_info}: List with information about the fitted prognostic model
on historical data. It has components:
\itemize{
\item \code{formula}: The \code{formula} with symbolic description of how the response
is modelled as function of covariates in the models
\item \code{model_fit}: A trained \code{workflow} - the result of \link{fit_best_learner}
\item \code{learners}: A \code{list} of learners used for the discrete super learner
\item \code{cv_folds}: The amount of folds used for cross validation
\item \code{data}: The historical data used for cross validation when fitting and
testing models
}
}
}
\description{
The procedure uses \link{fit_best_learner} to fit a prognostic model to historical data and uses
the model to produce counterfactual predictions as a prognostic score that is then adjusted
for as a covariate in the \link{rctglm} procedure. See
\href{https://arxiv.org/abs/2503.22284}{Powering RCTs for marginal effects with GLMs using prognostic score adjustment}
by Højbjerre-Frandsen et. al (2025) for more details.
}
\details{
Prognostic covariate adjustment involves training a prognostic model (using
\link{fit_best_learner}) on historical data (\code{data_hist}) to predict the response
in that data.

Assuming that the
historical data is representative of the comparator group in a “new” data
set (group 0 of the binary \code{exposure_indicator} in \code{data}), we can use the
prognostic model to predict the counterfactual
outcome of all observations (including the ones in the comparator group
for which the prediction of counterfactual outcome coincides with a
prediction of actual outcome).

This prediction, which is called the prognostic score, is then used as an
adjustment covariate in the GLM (the prognostic score is added to \code{formula}
before calling \link{rctglm} with the modified formula).

See much more details in the reference in the description.
}
\examples{
# Generate some data
n <- 100
b0 <- 1
b1 <- 1.5
b2 <- 2
W1 <- runif(n, min = -2, max = 2)
exp_prob <- .5

dat_treat <- glm_data(
  Y ~ b0+b1*abs(sin(W1))+b2*A,
  W1 = W1,
  A = rbinom (n, 1, exp_prob)
)

dat_notreat <- glm_data(
  Y ~ b0+b1*abs(sin(W1)),
  W1 = W1
)

learners <- list(
  mars = list(
    model = parsnip::set_engine(
      parsnip::mars(
        mode = "regression", prod_degree = 3
      ),
      "earth"
    )
 ),
    lm = list(
      model = parsnip::set_engine(
        parsnip::linear_reg(),
        "lm"
      )
    )
)
ate <- rctglm_with_prognosticscore(
  formula = Y ~ .,
  exposure_indicator = A,
  exposure_prob = exp_prob,
  data = dat_treat,
  family = gaussian(),
  estimand_fun = "ate",
  data_hist = dat_notreat,
  learners = learners)

# Pull information on estimand
estimand(ate)

}
\seealso{
Method to extract information of the prognostic model in \link{prog}. Function
used to fit the prognostic model is \code{\link[=fit_best_learner]{fit_best_learner()}}.

See \code{\link[=rctglm]{rctglm()}} for the function and class this inherits from.
}
