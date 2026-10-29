\name{fitted.magComp}
\alias{fitted.magComp}
\title{Expected values or replicate observations at fitted units}
\usage{
\method{fitted}{magComp}(object,
                         type = NULL,
                         target = c("expected", "observation"),
                         summary = c("mean", "samples", "quantile"),
                         probs = c(0.025, 0.5, 0.975), ...)
}
\arguments{
  \item{object}{A fitted \code{magComp} object.}
  \item{type}{Quantities to return; see \link{magCompPrediction} for supported
  choices. With \code{NULL}, defaults depend on the fitted components.}
  \item{target}{Use \code{"expected"} for conditional means or \code{"observation"}
  to simulate replicate observations from the fitted likelihood.}
  \item{summary}{Return posterior means, aligned retained samples, or quantiles.
  Means and quantiles pool all retained draws across chains.}
  \item{probs}{One or more probabilities strictly between zero and one, used
  for \code{summary = "quantile"}.}
  \item{...}{Must be empty. Unknown arguments are rejected.}
}
\description{
Returns expected quantities or replicate observations at the original fitted
units, using their original exposures and stored posterior latent effects.
Every fitted regression, group, spatial, and composition-factor term is included.
}
\details{
\code{target = "expected"} constructs quantities from the stored parameter
and latent-effect draws, without additional random simulation.
For the observation target, a full model draws an NB total followed by
multinomial allocation. With only magnitude, replicates draw just the total.
With only composition, replicates allocate the observed total among categories
with no abundance simulation. Expected composition-only category counts are
that observed total times the fitted probabilities. Conditional densities
additionally require the original area or BAF information.
\code{\link{ppcSummary}} compares replicates with observed counts.
Both targets support \code{\link{fixedAreaStems}} and
\code{\link{pointBasalArea}}. Point-sampling replicates are in-tree tallies;
observation magnitudes and category densities multiply those tallies by each
point's BAF to give sampled basal-area estimates.

See \link{magCompPrediction} for the four returned quantities, their array
dimensions, and the treatment of zero totals. New covariates or exposures belong
in \code{\link[=predict.magComp]{predict}}; \code{fitted()} always uses the
original fitted units. Observation simulation uses R's current random-number
state; use \code{set.seed()} externally for reproducibility.
}
\value{
A \code{\link{magCompPrediction}} named list with the requested quantity
components and \code{metadata}, even when only one quantity is requested.
Access counts with \code{result$count}. See that help page for array dimensions
and metadata fields.
}
\seealso{\code{\link[=summary.magCompPrediction]{summary}} for summarizing stored draws;
\code{\link[=predict.magComp]{predict}}, \code{\link{ppcSummary}},
\code{\link{standSummary}}}
