\name{magCompTerms}
\alias{magCompTerms}
\alias{iid}
\alias{gp}
\alias{nngp}
\alias{factor_iid}
\alias{factor_gp}
\alias{factor_nngp}
\title{Formula terms for magnitude and composition models}
\usage{
iid(...)
gp(...)
nngp(...)
factor_iid(...)
factor_gp(...)
factor_nngp(...)
}
\arguments{
  \item{...}{Grouping or coordinate column names and term settings described below.}
}
\description{
Formula-only markers for magnitude iid, GP, and NNGP effects and composition
iid, GP, and NNGP factor terms. These functions are parsed inside model formulas.
}
\details{
Use \code{iid(group)} for magnitude group effects, \code{gp(x, y)} for a dense
Gaussian process, or \code{nngp(x, y, m = 15)} for a nearest-neighbor process.
Each magnitude process has its own estimated variance. Composition uses
\code{factor_iid(group, K = 1)}, \code{factor_gp(x, y, K = 1)}, or
\code{factor_nngp(x, y, K = 1, m = 15)}. Its factors have unit variance. Unconstrained loadings have normal priors;
positive anchor-diagonal loadings have half-normal priors.

For each spatial term, different sampling units must have distinct coordinate
vectors. Exact duplicate locations cause an error during model preparation
that identifies the term and affected unit IDs. Nearby but distinct locations
are allowed. Multiple tree records within one unit and the repeated unit rows
in conditional allocation blocks are allowed; those rows share the unit's
spatial effects. This restriction applies to fitting. Repeated prediction
coordinates share a spatial-effect draw as described in
\code{\link[=predict.magComp]{predict}}.

Spatial calculations use the specified correlation matrices without diagonal
additions or conditional-variance clipping. Failed covariance calculations
report an error during initialization, sampling, or prediction.

All terms accept \code{name}, a unique nonempty name within the component or
allocation block. Otherwise names such as \code{iid_1} and \code{gp_1} follow
formula order within each type. Use these names for overrides in
\code{\link{magCompPriors}} and \code{starting}; inspect the resolved names with
\code{\link{magCompSettings}}.

Spatial terms accept \code{cov_model}: \code{"exp"} (default), \code{"matern"},
\code{"sep_exp"}, \code{"multi_res_sep_exp"}, or \code{"gneiting"}.
For a space-time family, supply numeric time as the final coordinate, for example
\code{nngp(x, y, year, cov_model = "sep_exp", m = 15)}. Repeated sites at
\emph{different} times are then allowed; repeated full space-time coordinates
for different sampling units remain an error. A spatial-only term on revisited
sites is not supported.

NNGP neighbors use Euclidean spatial distance and absolute time separation,
combined as \eqn{\sqrt{h^2+(s u)^2}}, where \code{st_scale = s} defaults to one.
This positive multiplier has spatial units per time unit and affects only
neighbor selection. Covariance calculations use the original spatial and
temporal lags. Ordering and neighbors remain fixed during fitting. Prediction
inherits the same setting and conditions jointly on observed and earlier new
locations.

Covariance priors, starts, and tuning use each term's \code{covariance} list.
For example, \code{covariance = list(phi = uniform(0.1, 5))} sets a prior,
\code{covariance = list(phi = 1)} sets a start, and
\code{covariance = list(phi = fixed(1))} fixes that parameter.
See \link{magCompCovModels} for the kernels, parameter names, and defaults.

Factor terms accept \code{anchors} (coefficient-stick names or row indices)
and \code{factor_names}. Anchor rows define zeros above the anchor diagonal and
positive, estimated diagonal entries. Automatic anchors use variation in observed continuation
ratios. Default factor names are the anchor categories. Factor rank and anchor
constraints are model specifications; diagonal loadings are estimated under
half-normal priors controlled by \code{loading_variance}.

Colon interactions such as \code{height:iid(stand)} or
\code{height:factor_nngp(x, y, K = 1)} multiply the latent contribution by the
covariate. Include \code{height} separately for an average regression slope.
Multiple terms may be used in each component or allocation block.

Named arguments are limited to those described above. Unknown arguments
are rejected. Set priors and initial values through their dedicated interfaces.
}
\value{
These functions stop if evaluated outside a model formula. They specify terms
for \code{magComp()} to interpret; they do not generate effects themselves.
}
\seealso{\code{\link{magCompPriors}}, \code{\link{magCompStarting}},
\link{magCompNotation}; the technical guide's Model terms and Factor
identification sections in \link{magComp-package}.}
\examples{
# These are formula objects, to be supplied to magComp().
~ height + iid(stand, name = "stand") + nngp(x, y, m = 15)
~ height + factor_nngp(x, y, K = 2, name = "spatial")
}
