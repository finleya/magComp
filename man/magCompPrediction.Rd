\name{magCompPrediction}
\alias{magCompPrediction}
\title{Results from fitted values and prediction}
\description{
This help topic describes the named list with class \code{magCompPrediction}, containing requested
\code{magnitude}, \code{proportion}, \code{density}, and \code{count}
components, followed by \code{metadata}. Predictions made with
\code{keep_effects = TRUE} also store a \code{latent} component for use by
\code{\link{samples}}; its internal indexing is not part of the public interface. \code{\link{standSummary}} also returns
this class, with a selected stand density, total, or proportion component.
Selecting one type still returns a
list: for example, \code{fitted(fit, type = "count")$count}.
}
\details{
The definitions below describe \code{fitted()} and \code{predict()} results.
See \code{\link{standSummary}} for the corresponding stand quantities and metadata.

Observed and new units each have two choices of target:
\tabular{lll}{
  Units \tab Conditional means \tab Simulated observations \cr
  Observed \tab \code{fitted(target = "expected")} \tab
    \code{fitted(target = "observation")} \cr
  New \tab \code{predict(target = "expected")} \tab
    \code{predict(target = "observation")} \cr
}
These calls also require the fitted object and, for prediction, \code{newdata}.
Defaults are \code{target = "expected"} and \code{summary = "mean"}.
With both model components, \code{type = NULL} selects all four quantities.
Standalone models return the quantities described below.

For fixed-area stems, let \eqn{M} be overall stem density, \eqn{p_j} the
probability of final category \eqn{j}, and \eqn{a} the sampling area.
For each retained posterior draw, expected quantities are constructed as
follows. Observation quantities additionally draw
\eqn{N \sim \mathrm{NB}(aM, r_M)} and then
\eqn{\mathbf{Y} \sim \mathrm{Multinomial}(N, \mathbf{p})}.
\tabular{lll}{
  Quantity \tab Expected \tab Observation \cr
  \code{magnitude} \tab \eqn{M} \tab \eqn{N/a} \cr
  \code{proportion} \tab \eqn{p_j} \tab \eqn{Y_j/N} \cr
  \code{density} \tab \eqn{M p_j} \tab \eqn{Y_j/a} \cr
  \code{count} \tab \eqn{a M p_j} \tab \eqn{Y_j} \cr
}
In a conditional two-category model, \eqn{p_j} is the product of the parent
probability and its conditional child probability. Size \eqn{r_M} uses the
fitted fixed value or its matching posterior draw when estimated.

Point-sampling basal-area models support both targets. Let \eqn{M^{BA}}
be basal-area density and \eqn{m} the in-tree tally. Expected total tally is
\eqn{M^{BA}/\mathrm{BAF}}. The observation target draws
\deqn{m \sim \mathrm{NB}(M^{BA}/\mathrm{BAF}, r_M^{BA}),}
using the fitted fixed size or its matching posterior draw, then allocates
\eqn{m} among categories when composition is fitted.

Expected magnitude and category density
are \eqn{M^{BA}} and \eqn{M^{BA}p_j}; observation values are the sampled
basal-area estimates \eqn{\mathrm{BAF}\,m} and \eqn{\mathrm{BAF}\,Y_j}.
Count is the expected or simulated in-tree tally, and proportion is \eqn{p_j}
or \eqn{Y_j/m}, respectively. The fixed-area expressions above therefore
apply with magnitude interpreted as basal-area density and exposure
\eqn{a = 1/\mathrm{BAF}}; this exposure is not a physical plot area.
BAF must be constant within each point but may vary between points.
BAF and basal-area density use matching units, such as square feet per acre.

For \code{summary = "mean"}, each component is a matrix with one row per
unit and one column per category, or one column named \code{"total"} for
magnitude. For \code{summary = "samples"}, category quantities have dimensions
unit by category by retained sample, plus a final chain dimension for multiple
chains. Magnitude omits the category dimension. For \code{summary = "quantile"},
category quantities have dimensions unit by category by probability, and
magnitude is a unit by probability matrix. Singleton dimensions are preserved;
means and quantiles pool chains. Row and category labels identify the fitted
units or prediction rows and final response cells.

Use \code{summary(result)} on a result retained with
\code{summary = "samples"} to obtain a table of means, standard deviations,
and quantiles by unit, quantity, and category. See
\code{\link[=summary.magCompPrediction]{summary}} for the returned columns.
This uses the stored draws and does not rerun prediction.

For magnitude-only models, the default quantities are \code{magnitude} and
\code{count}. Expected values are \eqn{M} and \eqn{aM}; observation values
are \eqn{N/a} and \eqn{N}, with a negative-binomial draw for \eqn{N}.
The count array retains a category dimension of length one named \code{"total"}.
\code{density} and \code{proportion} are unavailable because no composition
is fitted.

Composition-only models default to \code{proportion}. Expected probabilities
need neither area nor a supplied total. Conditional expected counts are
\eqn{N p_j}; replicate counts follow \eqn{\mathrm{Multinomial}(N,\mathbf p)}.
\code{fitted()} uses observed totals; \code{predict()} requires
\code{total_count}. Conditional densities divide these counts by exposure.
\code{magnitude} is unavailable. Expected probabilities remain defined at
zero totals, but observation proportions are \code{NA}; counts and densities
are zero. A supplied total is fixed across draws, not an estimate from a
separate abundance model. Nested composition retains the parent/child product.

These standalone-component rules apply to both sampling families. For
point sampling, conditional densities multiply tally counts by BAF, and
observation proportions are undefined when the tally is zero.

Metadata contains:
\describe{
  \item{\code{source}, \code{target}, \code{type}, \code{summary}, \code{family}}{
    The source function, target, requested quantities, summary, and response family.}
  \item{\code{units}, \code{unit_ids}, \code{exposure}}{
    Measurement units by quantity, row identifiers, and the exposures used.
    Exposure can be \code{NULL} for predictions requiring none.}
  \item{\code{components}}{Active model components.}
  \item{\code{count_conditioning}}{Count conditioning: \code{"observed_total"},
    \code{"supplied_total"}, or \code{"none"}. Composition-only CV identifies its
    assessment as conditional on observed totals.}
  \item{\code{total_count}}{Conditioning totals by unit when used, otherwise
    \code{NULL}.}
  \item{\code{n_samples}, \code{n_chains}, \code{n_draws}}{
    Retained draws per chain, number of chains, and their product.}
  \item{\code{sample_ids}, \code{chain_ids}, \code{probs}}{
    Draw positions within the original retained chain, chain labels, and requested probabilities
    for quantile summaries (otherwise \code{NULL}). Draws are aligned across
    quantities; observation simulation does not add another draw dimension.
    With \code{predict(n_draws = ...)}, positions identify the selected subset.}
  \item{\code{proportion_conditioning}}{
    \code{"positive_total"} for observation results, otherwise \code{"none"}.
    Observation proportions are undefined for zero totals. Their means and
    quantiles use positive-total draws only. Samples retain \code{NA} entries
    at zero totals, and all-zero units have \code{NA} proportion summaries.}
  \item{\code{zero_total_probability}, \code{n_positive_draws}}{
    For each unit, the empirical fraction of zero-total draws and number of
    positive-total draws across all chains. Returned for every observation
    result, including type subsets; \code{NULL} for expected results.}
}
All outputs within one result use the same draws. A zero simulated total
gives zero magnitude, density, and count, and undefined (\code{NA}) proportions.
Relations between draw-level quantities need not hold after separately taking
means or quantiles. Expected-target intervals describe uncertainty in conditional
means; observation-target intervals also include variation from the likelihood.
}
\seealso{\code{\link[=fitted.magComp]{fitted}},
\code{\link[=predict.magComp]{predict}};
the technical guide's Prediction and aggregation section in
\link{magComp-package}.}
