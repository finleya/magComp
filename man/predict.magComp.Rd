\name{predict.magComp}
\alias{predict.magComp}
\title{Expected values or observations at new sampling units}
\usage{
\method{predict}{magComp}(object, newdata,
                          type = NULL,
                          target = c("expected", "observation"),
                          summary = c("mean", "samples", "quantile"),
                          area = NULL, probs = c(0.025, 0.5, 0.975),
                          total_count = NULL, n_draws = NULL,
                          keep_effects = FALSE, ...)
}
\arguments{
  \item{object}{A fitted \code{magComp} object.}
  \item{newdata}{A data frame with one row per prediction unit, containing the
  covariates, group labels, and coordinates required by the fitted model.}
  \item{type}{Quantities to return; see \link{magCompPrediction} for supported
  choices. With \code{NULL}, defaults depend on the fitted components.}
  \item{target}{Use \code{"expected"} for conditional means based on posterior
  and new latent-effect draws, or \code{"observation"} to simulate observations
  from the fitted likelihood.}
  \item{summary}{Return means, aligned retained samples, or quantiles. Means
  and quantiles pool all retained draws across chains.}
  \item{area}{Prediction exposure: a positive finite scalar, a numeric vector
  with one value per row, or the name of a column in \code{newdata}. For
  \code{\link{fixedAreaStems}}, exposure is the sampling area for each predicted count.
  If omitted, the area column named in the fitted design must be in
  \code{newdata}. For \code{\link{pointBasalArea}}, supply exposure
  \code{1 / BAF}, or omit \code{area} and supply the fitted design's BAF column
  in \code{newdata}. Training exposures are never silently reused. Exposure
  is required for count predictions or observation simulation when magnitude
  is fitted. For composition-only models, it is required only for conditional
  density; probabilities and conditional counts need no exposure.}
  \item{probs}{One or more probabilities strictly between zero and one, used
  for \code{summary = "quantile"}.}
  \item{total_count}{For composition-only models, a numeric vector with one
  nonnegative integer total per prediction row, or the name of its column in
  \code{newdata}. Named vectors are matched to prediction unit IDs; unnamed
  vectors follow row order. Counts must not exceed R's integer limit. Required
  for conditional counts, conditional densities, or any observation target.
  Totals are treated as known and are not simulated. Not used when magnitude
  is fitted.}
  \item{n_draws}{Number of retained posterior draws to use per chain. The
  default \code{NULL} uses all retained draws. A positive integer selects
  approximately evenly spaced draws from each chain and must not exceed the
  number available. Selection preserves matching parameters and latent
  effects, and does not change the fitted object.}
  \item{keep_effects}{Retain the parameter and latent-effect draws used for
  prediction, for extraction with \code{\link{samples}}. Requires
  \code{summary = "samples"}. Defaults to \code{FALSE} to avoid storing
  additional arrays. Both expected and observation targets are supported.}
  \item{...}{Must be empty. Unknown arguments are rejected.}
}
\description{
Constructs posterior predictions at new units, including every fitted
regression, group, spatial, and composition-factor term. See
\link{magCompPrediction} for the returned quantities and the distinction
between expected values and simulated observations.
}
\details{
Prediction draws new latent effects conditional on each retained posterior
draw. Thus even \code{target = "expected"} can require random simulation.
With both components, \code{target = "observation"} additionally draws a
negative-binomial total and multinomial allocation. Magnitude-only models draw
only the total. Composition-only models instead allocate \code{total_count}
using a multinomial draw. Their expected category counts are
\eqn{N_0 p_{0j}}, conditional on the supplied total \eqn{N_0}; no magnitude
is estimated. Expected probabilities require neither totals nor areas. Use \code{set.seed()}
externally for reproducibility.

Use \code{n_draws} to reduce prediction time and memory while keeping all
fitting draws for parameter summaries and diagnostics. For example,
\code{n_draws = 100} with three fitted chains generates 300 predictive draws.
Selection occurs before prediction and retains every requested location in
each draw, including joint spatial effects. Selection itself uses no random
numbers. Output \code{metadata$sample_ids} identifies the positions selected
within the original retained chains. Fewer draws increase Monte Carlo
variability, especially in interval endpoints.

With \code{keep_effects = TRUE}, \code{samples(predictions, type = "effects")}
extracts the stored new-unit effects, and \code{type = "combined"} adds their
matching regression coefficients and loading-times-score contributions.
These are the same draws used to construct the response predictions; extraction
does not simulate again. All fitted latent terms are generated and retained,
even when \code{type} requests only magnitude. \code{standSummary()} aggregates
response draws and does not carry these location-level effects into its output.

Observed iid groups reuse posterior effect draws. Each unseen group receives
one new Gaussian effect per term and posterior draw, shared by all rows in
that group. Magnitude terms use their sampled variances; composition factor
scores use unit variance. For space-time terms, a location includes the final time coordinate: an exact
location-time match reuses its fitted effect, while a revisit at a new time
receives a conditional draw. Prediction inherits the fitted covariance family
and \code{st_scale} for neighbor searches.
Known spatial locations reuse fitted effects, and
repeated new coordinates share one sampled spatial effect per term and draw.
Dense-GP prediction draws new effects jointly from their conditional Gaussian
distribution. NNGP prediction orders new locations by their first coordinate,
using subsequent coordinates to break ties, and uses neighbors among observed and earlier new locations in the same call.
There is no conditionally independent spatial option. The NNGP approximation
depends on the requested location set; separate calls omit cross-call links
and do not reproduce a single joint call. Supply the complete desired support
in one call when its joint uncertainty matters.

Both targets support \code{\link{fixedAreaStems}} and
\code{\link{pointBasalArea}}.
For point sampling, the expected total tally is
\eqn{M_0^{BA}/\mathrm{BAF}_0}; multiplying simulated tallies by BAF gives
sampled basal-area estimates. Composition-only models condition on the supplied
tally total. For an observed location, prediction
reuses its stored spatial effect, but covariates and exposure still come from
\code{newdata} and \code{area}. Use \code{fitted()} to obtain results for all
original units with their original inputs.

\code{\link{standSummary}} aggregates stored density or magnitude samples
using explicitly supplied represented land areas. It preserves the selected
target. These weights serve a different purpose from the sampling exposures
used to generate counts here. Neither function determines areas from polygons.
}
\value{
A \code{\link{magCompPrediction}} named list with requested quantity components
and \code{metadata}, even for one requested quantity. All quantities share
posterior draw and chain identities. See that help page for dimensions and
metadata fields.
}
\seealso{\code{\link[=summary.magCompPrediction]{summary}} for summarizing stored draws;
\code{\link[=fitted.magComp]{fitted}}, \code{\link{standSummary}},
\code{\link{magCompCV}}}
