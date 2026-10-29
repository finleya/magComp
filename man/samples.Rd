\name{samples}
\alias{samples}
\alias{samples.magComp}
\alias{samples.magCompPrediction}
\alias{print.magCompSamples}
\title{Extract posterior draws and combined coefficients}
\usage{
samples(object, ...)

\method{samples}{magComp}(object, type = c("parameters", "effects", "combined"),
                              start = 1, thin = 1, chains = NULL,
                              combine_chains = FALSE, group = NULL, ...)

\method{samples}{magCompPrediction}(object,
                  type = c("parameters", "effects", "combined"),
                  start = 1, thin = 1, chains = NULL,
                  combine_chains = FALSE, group = NULL, ...)
}
\arguments{
  \item{object}{A fitted \code{\link{magComp}} model, or the result of
\code{predict(..., summary = "samples", keep_effects = TRUE)}.}
\item{type}{\code{"parameters"} returns all raw parameter arrays in the
component/block/term structure below. Its default summary reports scalar
parameters. \code{"effects"} selects individual magnitude effects and
composition factor scores, including them in summaries and traceplots.
\code{"combined"} constructs fixed coefficients plus their matching random
coefficients at fitted or predicted units. Only coefficients with random-effect terms
are returned.}
\item{start}{First saved draw to retain within each chain, counted after the
burn-in and thinning already applied by the fit. Defaults to 1.}
\item{thin}{Additional thinning interval among saved draws. Defaults to 1.}
\item{chains}{Chain numbers or names to select, in the requested order.
\code{NULL} selects all chains.}
\item{combine_chains}{Pool selected chains after constructing quantities.
Original draw and chain identities are retained for summaries and diagnostics.}
\item{group}{For \code{type = "combined"}, a column name in the fitted unit data or prediction
\code{newdata} for averaging within groups, or \code{".all"} for averaging over
all units. \code{NULL} returns each unit. Averages give each unit equal
weight and are calculated within each draw.}
\item{...}{Additional arguments, currently unused.}
}
\description{
Extracts retained posterior draws without running the sampler or generating
predictions. The result has \code{summary()} and \code{plot()} methods.
Draw selection is synchronized across every parameter and latent effect.
}
\details{
For prediction objects, the stored latent draws are those used to construct
the response predictions, before any negative-binomial or multinomial draws.
Extraction does not generate new effects. The parameter arrays contain the
matching subset selected by \code{predict(n_draws = ...)}; draw identities
and original MCMC iterations are preserved. \code{start} and \code{thin}
select positions within that subset. Magnitude effects are indexed by new
unit; spatial scores are indexed by new location, and iid composition scores
by allocation-block row (including parent context for conditional blocks).
Response draws themselves remain available in the prediction object's
\code{magnitude}, \code{proportion}, \code{density}, or \code{count} components.

Resolved model settings, starting values, and extracted samples contain only
active components. Supplied starting values or nondefault priors for a
disabled component are rejected by \code{magComp()}.

For \code{type = "parameters"}, the returned structure is:
\preformatted{
list(
  magnitude = list(
    coefficients = coefficient_draws,
    size = size_draws,                  # only when estimated
    terms = list(term_name = list(
      effects = effect_draws, variance = variance_draws,
      covariance = list(phi = phi_draws) # estimated parameters only
    ))
  ),
  composition = list(blocks = list(block_name = list(
    coefficients = coefficient_draws,
    terms = list(term_name = list(
      loadings = loading_draws, scores = score_draws,
      covariance = list(phi = phi_draws) # estimated parameters only
    ))
  )))
)
}
Magnitude coefficients and effects are parameter-by-sample matrices. Scalar
variance and size arrays have one parameter row; each covariance parameter has an array
with one row per factor for which it is estimated. Composition coefficient arrays are
stick-by-design-column-by-sample. Loading arrays are stick-by-factor-by-sample;
score arrays are group-or-location-by-factor-by-sample. Multiple chains add
a final chain dimension to every array. Parameter dimensions and singleton
dimensions are retained, and dimnames identify coefficients, sticks, groups,
locations, factors, and chains. Sample positions align across the result.

Fixed size and covariance values are omitted from samples. Fixed covariance
values are retained by \code{magCompStarting()} for restarting the same model.
Use \code{\link{magCompSettings}} to inspect these fixed values. Loading arrays retain
the fixed upper-anchor zeros and estimated positive diagonals. Factor variances and smoothing variances
are fixed controls and have no posterior sample arrays.

For example, \code{draws$magnitude$terms$stand$variance} contains the
posterior draws of the variance whose prior is specified by
\code{priors$magnitude$terms$stand$variance} and whose initial value is
\code{starting$magnitude$terms$stand$variance}. The prior entry holds
inverse-gamma shape and scale; the starting entry holds a single variance.

Use \code{\link[=fitted.magComp]{fitted}} and \code{\link[=predict.magComp]{predict}} for
draws of response quantities. Use \code{\link{magCompStarting}} with
\code{use = "last"} or \code{use = "mean"} to turn parameter samples into
a starting-value list.
}
\value{
An object of class \code{magCompSamples}. Raw parameters and effects are
nested lists of arrays. Combined coefficients are a named list of
saved-draw-by-chain matrices, or vectors when chains are pooled. Their
\code{info} attribute identifies component, block, coefficient, allocation
row, parent context, and fitted or predicted unit or group.

The \code{sample_info} attribute records type, retained draw positions,
original MCMC iteration numbers, chain identities, grouping, and pooling.
Pooled draws are concatenated one chain at a time in the selected order;
every chain uses the same selected draw positions.
}
\section{Combined coefficients}{
For a magnitude random intercept with both group and spatial terms, the
combined coefficient is \eqn{\beta_0 + \delta_{g(i)} + v(s_i)}. A varying
slope instead combines its fixed slope with the random coefficients attached
to that covariate. A missing fixed counterpart contributes zero.

Composition combinations sum the loading-times-score contributions across
all factors and matching terms, then add the corresponding fixed coefficient.
They support iid, GP, and NNGP terms and retain conditional parent contexts.
Loadings and scores are multiplied within matching draws and chains.
Covariates are not multiplied into these coefficients. Response means and
observation draws belong to \code{fitted()} and \code{predict()}.

Group averages summarize coefficients over the supplied units, not forest density
or abundance. Use \code{standSummary()} for area-weighted response summaries.
Where several latent terms contribute to one coefficient, all matching terms
are included. These are combined coefficients, not separately attributed
term effects. New groups or locations require \code{predict()} with
\code{keep_effects = TRUE}, followed by \code{samples()} on its result.
}
\seealso{\link{magCompNotation}, \code{\link{magCompSettings}},
\code{\link[=summary.magComp]{summary}}}
\examples{
rows <- expand.grid(plot = paste0("p", 1:4), cls = c("a", "b"))
rows$count <- rep(1:4, 2)
rows$area <- 1
rows$stand <- rep(c("s1", "s1", "s2", "s2"), 2)
dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                   design = fixedArea(area))
fit <- magComp(dat, magnitude = ~ iid(stand, name = "stand"),
  control = magCompControl(n_samples = 6, n_burn = 3,
    verbose = FALSE))
draws <- samples(fit)
draws$magnitude$coefficients
draws$magnitude$terms$stand$variance
summary(fit, type = "effects")
combined <- samples(fit, type = "combined", group = "stand")
summary(combined)
# Equivalent convenience call:
summary(fit, type = "combined", group = "stand")
}
