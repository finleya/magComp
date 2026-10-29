\name{standSummary}
\alias{standSummary}
\alias{standSummary.magCompPrediction}
\title{Summarize predictive densities for stands, strata, or forests}
\usage{
standSummary(object, ...)

\method{standSummary}{magCompPrediction}(object, data, unit, group = NULL, represented_area,
                            type = c("density", "total", "proportion"),
                            summary = c("mean", "samples", "quantile"),
                            probs = c(0.025, 0.5, 0.975), ...)
}
\arguments{
  \item{object}{A \code{\link{magCompPrediction}} result from \code{fitted()}
  or \code{predict()}, retained with \code{summary = "samples"} and containing
  \code{density} or \code{magnitude}. Both sampling families are supported.}
  \item{data}{A data frame or tibble with one row for each prediction
  location's contribution to a stand or other group. A cell crossing stand
  boundaries has a separate row and represented area for each stand.
  Rows may be in any order and may select a subset of prediction locations.}
  \item{unit}{The column in \code{data} containing prediction location IDs,
  supplied as a bare name or quoted column name. IDs must match the units
  stored in \code{object}; the column name can differ from the original ID column.}
  \item{group}{The column in \code{data} identifying stands or other groups,
  supplied as a bare name or quoted column name. Omission or \code{NULL}
  aggregates all rows into one group named \code{"all"}. Groups appear in
  order of first occurrence.}
  \item{represented_area}{The column in \code{data} containing positive finite
  land areas, supplied as a bare name or quoted column name. For cells crossing
  stand boundaries, use the cell's area within that row's stand. Use units
  consistent with the density, such as hectares for trees per hectare.}
  \item{type}{One stand quantity: area-weighted density (the default), total,
  or category proportion. Proportions require category density samples.}
  \item{summary}{Return posterior means, aligned retained samples, or quantiles.
  Means and quantiles pool retained draws across chains.}
  \item{probs}{Probabilities strictly between zero and one, used for
  \code{summary = "quantile"}.}
  \item{...}{Must be empty. Unknown arguments are rejected.}
}
\description{
Aggregates existing density draws using represented land areas. No new
predictions are generated: the supplied expected or observation target and
posterior draw alignment are preserved.
}
\details{
For each category and draw, total is the sum of density times represented area.
Density is that total divided by the group's summed area. Proportions divide
category totals by their sum across categories. They describe pooled stand
composition, not an area-weighted average of local proportions.

Full models supply category \code{density} draws. Magnitude-only models supply
\code{magnitude}, giving overall density or total under a single category named
\code{"total"}. When both quantities are present, category density is used.
Composition-only models can supply conditional category densities:
\code{fitted()} uses observed totals and \code{predict()} requires
\code{total_count}, with plot area or BAF supplied for the density conversion.
Those totals remain fixed across draws, and the conditioning is recorded in
metadata. Proportions alone lack the abundance information needed for pooled
stand composition and cannot be aggregated by this function.

For fixed-area sampling, stand densities are stems per unit area and totals
are stems. For point sampling, densities are basal area per unit area and
totals are basal area. Point-sampling \code{fitted()} and \code{predict()}
already multiply category tallies by BAF; \code{standSummary()} does not apply
BAF again. Point-sampling proportions describe basal-area composition.

Represented areas within a group define its summary area. For a cell grid,
intersect the cells with the target polygons and supply the intersection areas.
Each cell's same density draw contributes to every group it overlaps, treating
that density as constant within the cell. Multiple rows for the same location
and group add their represented areas, as for disconnected intersection pieces;
include each piece once. Omitting \code{group} pools all supplied areas.
For grid points, equal weights can be a stand's polygon area divided by its
number of points. These are aggregation weights, distinct from plot areas or
reciprocal BAF used to generate counts. No areas are inferred from polygons,
grid spacing, or sampling exposure. Fitted locations also require explicit
represented areas; sampled plots do not automatically represent a whole stand.

The same calculation applies to predictions from regression, stratum and
spatial models. Their effects are already included in the supplied draws.
Supply all locations for the desired support in one prediction call when
joint spatial uncertainty matters. Aggregation preserves dependence already
in those draws but cannot reconstruct dependence omitted by separate calls.
When model relationships differ by stratum, prediction units must represent
the appropriate stratum; aggregation weights do not change the fitted effects.

Internally, a sparse group-by-location area matrix multiplies the density
samples, preserving category, draw and chain alignment without constructing
a long table of draws. Density divides each group's totals by its summed area.

Totals and densities retain zero-valued draws. Proportions are \code{NA} for
zero-total draws; their means and quantiles condition on a positive total.
If no positive draws remain, proportion summaries are \code{NA}. Metadata
records the zero-total fraction and number of positive draws by stand.
}
\value{
A \code{\link{magCompPrediction}} named list containing the selected
\code{density}, \code{total}, or \code{proportion} component and \code{metadata}.
Means are stand-by-category matrices, samples are stand-by-category-by-draw
arrays with the input chain dimension retained, and quantiles are
stand-by-category-by-probability arrays. Magnitude inputs use a category
named \code{"total"}. Stand and category labels identify the rows and columns.
Use \code{summary(result)} on retained stand samples for a table of means,
standard deviations, quantiles, and contributing draw counts.

Metadata records the original target, family, units, count conditioning,
draw and chain identities, stand IDs in \code{unit_ids}, and summed
\code{represented_area} by stand. \code{prediction_source} identifies
\code{fitted} or \code{predict}. The \code{aggregation} list retains the
original prediction IDs, sampling exposures, and any conditioning totals.
Its \code{data} element records the supplied mapping with columns \code{unit},
\code{group}, and \code{represented_area}; these are location contributions,
not stand totals.
}
\examples{
\dontrun{
# Plot areas determine counts; intersection areas weight densities.
pred <- predict(fit, grid, area = "plot_area", target = "observation",
                type = "density", summary = "samples")
stands <- standSummary(pred, data = cell_stand_areas,
                       unit = prediction_id, group = stand_id,
                       represented_area = area_ha, summary = "samples")
summary(stands)
stands$density

# One row per cell is sufficient for a forest summary.
standSummary(pred, data = grid, unit = prediction_id,
             represented_area = cell_area_ha, summary = "quantile")

# The same interface handles basal-area densities from point sampling.
point_pred <- predict(point_fit, grid, target = "observation",
                      type = "density", summary = "samples")
standSummary(point_pred, data = grid, unit = prediction_id,
             group = stand_id, represented_area = expansion_area_acres,
             type = "total")$total

}
}
\seealso{\code{\link[=fitted.magComp]{fitted}},
\code{\link[=predict.magComp]{predict}},
\code{\link[=summary.magCompPrediction]{summary}};
the technical guide's Prediction and aggregation section in
\link{magComp-package}.}
