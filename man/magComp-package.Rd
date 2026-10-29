\name{magComp-package}
\alias{magComp-package}
\docType{package}
\title{magComp: model fitting, prediction, and summaries}
\description{
Bayesian magnitude-composition models for fixed-area stem counts and
point-sampling basal area. A negative-binomial model describes the total
count, and multinomial allocation distributes that count among categories.
Users choose one categorical variable or two in a conditional allocation;
neither their names nor their conditioning order is tied to a forestry response.
Models can include regression, group effects, spatial effects, and composition
factors, with smoothing across ordered categories. Fit either component alone
by setting the other formula to \code{NULL}; the existing nested composition
structure remains available without magnitude.
}
\section{Getting started}{
The \href{https://github.com/finleya/magComp#documentation}{website articles}
simulate data and illustrate magnitude, composition, and full models using
forest inventory examples. Their Quarto sources and HTML previews are in the
source repository; the articles are not installed with the package.

\describe{
  \item{Prepare data}{Use \code{\link{magCompData}} with a
  \code{\link{fixedArea}} or \code{\link{pointSampling}} design. Each sampling
  unit has total counts, counts by one category, or counts by two crossed
  categories. Standalone composition does not require a sampling design.}
  \item{Specify and fit a model}{Use \code{\link{magComp}}. For two categories,
  \code{\link{magCompComposition}} specifies the leading allocation and the
  conditional child allocation. Formula syntax is described in
  \link{magCompTerms}. Set priors with \code{\link{magCompPriors}},
  smoothing with \code{\link{magCompSmooth}}, and MCMC controls with
  \code{\link{magCompControl}}.}
  \item{Inspect settings}{Use \code{\link{magCompSettings}} for resolved priors,
  parameter constraints, and editable starting values. Use
  \code{\link{magCompStarting}} to extract starts from a previous fit.}
  \item{Check the chains}{Use \code{\link[=summary.magComp]{summary}} and
  \code{\link[=plot.magComp]{plot}}. Both use \code{\link{samples}} and support
  individual effects and combined coefficients through their \code{type} argument.
  \code{\link{spatialDiagnostics}} summarizes spatial parameters.
  \code{\link{samples}} extracts all retained parameter draws.}
  \item{Obtain fitted values and predictions}{Use
  \code{\link[=fitted.magComp]{fitted}} at observed units and
  \code{\link[=predict.magComp]{predict}} at new units. Both offer expected
  quantities or simulated observations. \link{magCompPrediction} explains
  the returned quantities and dimensions.}
  \item{Summarize and assess}{Use
  \code{\link[=summary.magCompPrediction]{summary}} on stored prediction draws
  for means and intervals. Use \code{\link{standSummary}} to aggregate stored
  expected or observation density draws with represented land-area weights, \code{\link{ppcSummary}} to compare observed and replicate
  counts, and \code{\link{magCompCV}} for cross-validation predictions.}
}
Call generic functions such as \code{summary(fit)}, \code{plot(fit)},
\code{fitted(fit)}, and
\code{predict(fit, newdata)};
R selects the appropriate method. \link{magCompNotation} and
\link{magCompPrediction} are explanatory help topics, not functions.
}
\section{Technical document}{
The \href{https://github.com/finleya/magComp/blob/master/dev/modeling_and_software_details/magComp.pdf}{Methods and Implementation Details}
guide gives the model equations, notation, priors, full conditional and
Metropolis updates, numerical solvers, and prediction calculations.
Its source and PDF are maintained in
\file{dev/modeling_and_software_details/} in the source repository.
Help pages refer to its section titles so those references remain useful
when equation and section numbers change.
}
\seealso{\code{\link{magComp}}, \link{magCompNotation}}
