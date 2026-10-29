\name{magCompStarting}
\alias{magCompStarting}
\alias{magCompStarting.magComp}
\title{Inspect initial values or extract starts from posterior samples}
\usage{
magCompStarting(object, ...)

\method{magCompStarting}{magComp}(object,
  use = c("initial", "last", "mean"), chain = 1, ...)
}
\arguments{
  \item{object}{A prepared or fitted \code{magComp} model.}
  \item{use}{\code{"initial"} returns the resolved values used to initialize the
  model. \code{"last"} extracts the selected chain's last retained draw;
  \code{"mean"} extracts posterior means across all retained chains. The latter
  two choices require a fitted model.}
  \item{chain}{Chain index for \code{use = "last"}.}
  \item{...}{Currently unused.}
}
\description{
Returns a starting-value list organized by component, block, and term, as
in \code{\link{magCompSettings}}. Pass it to \code{magComp(starting = )}.
Extracting posterior values initializes a new run; it does not resume RNG state
or Metropolis adaptation from an earlier run.
}
\details{
Resolved model settings, starting values, and extracted samples contain only
active components. Supplied starting values or nondefault priors for a
disabled component are rejected by \code{magComp()}.

Estimated covariance-parameter starts must lie strictly between their uniform-prior bounds
because their Metropolis proposals use a bounded-logit transformation.
Starting values have the following structure:
\preformatted{
list(
  magnitude = list(
    coefficients = named_vector,
    size = positive_scalar,  # only when size is estimated
    terms = list(term_name = list(
      effects = named_vector, variance = positive_scalar,
      covariance = list(phi = value)  # spatial/space-time terms
    ))
  ),
  composition = list(blocks = list(block_name = list(
    coefficients = named_matrix,
    terms = list(term_name = list(
      loadings = named_matrix, scores = named_matrix,
      covariance = list(phi = value)  # one value per factor
    ))
  )))
)
}
Only the entries a user wants to override need to be supplied. Omitted entries
receive automatic starts. Obtain the correct term names, dimensions, and
dimnames from a prepared model via \code{magCompSettings(model)$starting}.
Use \code{name = "stand"} inside a formula term to give it a stable name.

Coefficient matrices have one row per coefficient stick and one column per
fixed-effect design column. Loading matrices have coefficient-stick rows and
factor columns; score matrices have group/location rows and factor columns.
Named inputs are matched to model names, including reordered vectors and
matrix rows/columns. Unnamed inputs must follow the template's order. Inputs
for coefficients, effects, scores, and loadings must have exactly the required
dimensions and finite values; scalar recycling is not performed for these arrays. Loading starts must have positive values on the anchor diagonal and
zeros above it. Supplied valid loading and score starts are preserved.

Automatic factor starts use residuals from regressions of transformed
empirical continuation ratios. Conditional child blocks have separate
regressions for each parent. Terms are initialized sequentially, with iid
terms before spatial terms and formula order within each group. Each term's
contribution is removed before initializing the next term by SVD.

Automatic loadings and scores are rotated and scaled to satisfy the anchor
constraints and give score columns unit mean square. This is an initialization
convention, not a constraint on sampled scores. Valid user-supplied loadings
or scores are preserved; an omitted matrix is initialized automatically.
The technical guide's Data preparation, aggregation, and initialization
section gives the details; see \link{magComp-package}.

Magnitude iid/spatial variance starts use the variance of initialized effects,
bounded below by \code{1e-4}, with a fallback of one if nonpositive or undefined.
Covariance starts have automatic coordinate-based defaults. Estimated values
must lie strictly inside their uniform-prior supports. To fix a covariance
parameter, wrap its starting value in \code{fixed()}, for example
\code{covariance = list(phi = fixed(2))}. Fixed markers are preserved when
extracting starts from a fitted model. Common scalars apply to all factors;
named vectors or lists specify individual factors. See \link{magCompCovModels}.

}
\value{
A named list with \code{magnitude} and \code{composition} branches. Initial
templates contain all fitted model parameter blocks, including fixed covariance values. Extracted starts require the same model structure and group/location
labels in the next fit.
}
\seealso{\link{magCompNotation} for the symbols corresponding to each entry;
\code{\link{samples}} for all retained draws in matching branches.}
