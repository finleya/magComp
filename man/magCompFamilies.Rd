\name{magCompFamilies}
\alias{magCompFamilies}
\alias{fixedAreaStems}
\alias{pointBasalArea}
\title{Response families for stem counts and basal area}
\usage{
fixedAreaStems()

pointBasalArea()
}
\description{
Constructs a response family for \code{\link{magComp}}. Both families model
negative-binomial totals and multinomial category allocation. Either component
can be fitted alone by setting the other formula to \code{NULL}.
}
\details{
\code{fixedAreaStems()} uses fixed-area plot counts. Its magnitude is stem
density \eqn{M_i}, so the expected plot total is \eqn{a_i M_i} for area
\eqn{a_i}. Supply a \code{\link{fixedArea}} design when preparing data.

\code{pointBasalArea()} uses horizontal point-sampling in-tree tallies.
Its magnitude is basal area per unit area \eqn{M_i^{BA}}, so the expected tally
is \eqn{M_i^{BA}/\mathrm{BAF}_i}. Supply a \code{\link{pointSampling}} design.
Tree-factor-expanded stem counts or other attributes are not fitted responses.

Negative-binomial size is positive and finite. Both families estimate it by
default under a Gamma(shape = 1, rate = 0.01) prior when the magnitude formula
has no latent process terms. With magnitude iid or spatial terms, size is
fixed at \code{control$nb_size} (default 1000). Use \code{\link{magCompPriors}}
to change the prior or explicitly fix size. Large sizes approximate a Poisson model;
the sampler does not fit an exact Poisson model.

Both families support expected fitted values and predictions, as well as
replicate or new observations from the fitted likelihood. Point-sampling count
outputs are in-tree tallies, with mean \eqn{M_i^{BA}/\mathrm{BAF}_i};
observation magnitude multiplies the tally by BAF to give a sampled basal-area
estimate. Composition-only models simulate allocations conditional on observed
or supplied totals under either family, including for posterior predictive checks. Plot area or BAF
is unnecessary for composition-only fitting and probabilities; it is needed
for conditional density. \code{\link{standSummary}} aggregates density or
magnitude samples under either family using represented land areas.
}
\value{
A \code{magCompFamily} object with \code{likelihood} equal to
\code{"negative_binomial_multinomial"}.
}
\seealso{\code{\link{magCompData}}, \code{\link{magCompControl}},
\code{\link{magCompPriors}}; the technical guide's Fixed-area count model
and Point-sampling basal area sections in \link{magComp-package}.}
