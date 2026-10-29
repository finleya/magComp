\name{magCompData}
\alias{magCompData}
\alias{fixedArea}
\alias{pointSampling}
\title{Prepare data for magnitude-composition models}
\usage{
magCompData(data, unit, category = NULL, categories = NULL,
            coords = NULL, tree_factor = NULL, tree_count = NULL,
            design = NULL, attributes = NULL, ordered = NULL)

fixedArea(area = NULL)

pointSampling(baf = NULL)
}
\arguments{
  \item{data}{Data frame with one row per observed tree or one row per category tally.}
  \item{unit}{Sampling-unit identifier column, supplied as a bare name or a string.}
  \item{category}{Optional single category column. Omit categories for total-count magnitude-only data.}
  \item{categories}{Optional one- or two-column category specification,
  supplied as \code{c(parent, child)} or \code{c("parent", "child")}. Use this instead of
  \code{category} for two-way composition tables.}
  \item{coords}{Optional coordinate columns, usually supplied as \code{c(x,y)}.}
  \item{tree_factor}{Optional tree expansion factor column retained for
  expanded attribute summaries and future weighted-response families. Current
  point-sampling basal-area fits use \code{pointSampling(baf = ...)} rather
  than a tree-factor likelihood.}
  \item{tree_count}{Optional integer tally-count column. When supplied,
  category count matrices are built from this column rather than assuming one
  observed tree per row. Zero-count rows can be used to retain sampling units
  with no observed trees in the modeled categories.}
  \item{design}{Sampling design helper, e.g., \code{fixedArea(area)}.}
  \item{attributes}{Optional tree-level attributes to expand by
  \code{tree_factor}. These summaries are stored with the prepared data for
  diagnostics and future continuous-response families.}
  \item{ordered}{Optional logical vector indicating which category columns are
  ordered. Names should match category column names. Ordered factor columns are
  also detected automatically. An ordered category must be a factor with levels
  in the intended order, or numeric (sorted numerically). Character labels
  require explicit factor levels to avoid accidental alphabetical ordering.}
  \item{area}{Fixed plot area, supplied as a scalar, a vector named by sampling-unit ID,
  or a column name constant within each unit. Non-scalar vectors must be named;
  values are matched by unit ID rather than position.}
  \item{baf}{Basal area factor, supplied as a scalar, a vector named by sampling-unit ID,
  or a column name constant within each unit. Non-scalar vectors must be named;
  values are matched by unit ID rather than position.}
}
\description{
Builds sampling-unit summaries, count matrices, allocation-block summaries,
coordinate matrices, design metadata, fixed-area exposure summaries, and
expanded tree-list attributes used by \code{magComp()}. A one-category response
creates one marginal allocation block. A two-category response creates a
leading marginal block and a conditional child block. Category ordering
metadata controls whether RW1/RW2 smoothing can be used for an allocation
block.
Coordinates must be constant within each sampling unit. When a model is
prepared, the same check is applied to its predictors and grouping columns.
Tree-level responses and attributes may vary within a unit.
}
\details{
Category columns can describe any categorical response dimensions; species
and diameter classes are examples. With two columns, their order in
\code{categories} determines the parent and child. Either can be ordered or
unordered. Supply any desired classes as category values before preparing
the data; the constructor does not create classes from continuous measurements.

Without category columns, counts are summed by sampling unit into one column
named \code{"total"}. Supply \code{tree_count} for tallied totals or use one
row per tree. Fit these data using \code{composition = NULL}.
Composition-only data still require integer category counts but can omit
\code{design}; provide area or BAF information if conditional density outputs
are needed later.

Column specifications accept bare names or strings. Sampling units retain their
order of first appearance. Factor levels determine category order; otherwise
numeric categories are sorted numerically and character categories alphabetically.
Use explicit factor levels when a particular category order is needed, and
ordered factors for categories to be smoothed. Unobserved category cells have
zero counts; unused factor levels are retained.

To retain an empty sampling unit, provide a row with \code{tree_count = 0}
and valid category labels, covariates, and design information. Without a tally
column, every row counts as one tree.

A fixed-area count model requires one common sampling area across all categories
within each unit. Nested subplots with different areas for different categories
must first be put on a common sampling support. Tree expansion factors do not
convert those unequal-area counts into a valid common-area likelihood.

Tree-level attributes are stored as sums of \code{tree_factor * attribute}
by unit and category. They are not responses in the implemented count models.
Basal-area fits use in-tree tallies with a
\code{pointSampling()} design and the \code{pointBasalArea()} family.
}
\value{
The data constructor returns a \code{magCompData} object.
Both design helpers return sampling-design objects.
}
\seealso{\code{\link{magComp}}, \code{\link{magCompComposition}},
\code{\link{fixedAreaStems}}, \code{\link{pointBasalArea}}}
\examples{
rows <- data.frame(plot = c("p1", "p1", "p2"),
                   dbh_class = c("small", "large", "small"),
                   count = c(5, 2, 0), area = c(0.05, 0.05, 0.08))
rows$dbh_class <- ordered(rows$dbh_class, levels = c("small", "large"))
dat <- magCompData(rows, unit = plot, category = dbh_class,
                   tree_count = count, design = fixedArea(area))
dat$counts
}
