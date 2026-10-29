\name{magComp}
\alias{magComp}
\title{Fit a magnitude-composition model}
\usage{
magComp(data, family = fixedAreaStems(), magnitude = ~1,
        composition = ~1, smooth = magCompSmooth(),
        priors = magCompPriors(), starting = list(),
        control = magCompControl(), tuning = list(), ...)
}
\arguments{
  \item{data}{A \code{\link{magCompData}} object, or a data frame with its
  preparation arguments supplied through \code{...}.}
  \item{family}{A response family: \code{\link{fixedAreaStems}} (the default)
  or \code{\link{pointBasalArea}}.}
  \item{magnitude}{One-sided formula for the scalar magnitude component.
  Use \code{NULL} to fit composition only. See \link{magCompTerms} for group and spatial terms.}
  \item{composition}{A one-sided formula for a single-category response,
  or a \code{\link{magCompComposition}} specification with one formula per
  allocation block. Use \code{NULL} to fit magnitude only.}
  \item{smooth}{Composition smoothing from \code{\link{magCompSmooth}}.}
  \item{priors}{Prior settings from \code{\link{magCompPriors}}.}
  \item{starting}{Partial or complete starting-value list. Omitted values
  receive automatic starts. Use the \code{starting} template from
  \code{\link{magCompSettings}}; see \code{\link{magCompStarting}} for its structure.}
  \item{control}{MCMC settings from \code{\link{magCompControl}}. Sampling is the default.
  To inspect settings without sampling, use \code{run_sampler = FALSE}.}
  \item{tuning}{Optional covariance proposal standard deviations, organized
  like \code{starting} with a \code{covariance} list in each spatial term.
  Values are on the bounded-logit scale; default 0.1. See \link{magCompCovModels}.}
  \item{...}{Arguments passed to \code{\link{magCompData}} when \code{data}
  is a data frame.}
}
\description{
Fits a Bayesian model in which a scalar magnitude describes overall density
and a composition allocates that magnitude among categories. For fixed-area
stems, plot counts follow a negative-binomial total with multinomial allocation.
For point sampling, the count model describes in-tree tallies and the magnitude
is basal area per unit area.
}
\details{
Explicit \code{NULL} disables a component; omitting either formula still uses
its intercept-only default \code{~1}. At least one component is required.
With \code{composition = NULL}, fit total counts without category labels or
sum the supplied category counts. With \code{magnitude = NULL}, fit category
allocation conditional on observed totals. Composition may still have one
categorical response or the existing parent/child structure, with different covariates,
factors, and smoothing in each block. Composition-only fitting requires
integer category counts, not proportions, and does not require areas.

Only active components are initialized, sampled, and reported in parameter
settings, starts, samples, and diagnostics. Nondefault priors or starting
values for disabled components are rejected. Smoothing is unused and must
remain at its default when composition is disabled.
\code{\link[=fitted.magComp]{fitted}}, \code{\link[=predict.magComp]{predict}},
\code{\link{magCompCV}}, and \code{\link{ppcSummary}} support each configuration;
\code{\link{standSummary}} aggregates their density or magnitude samples using
explicit represented land areas.

Both components support regression effects. Magnitude formulas can add iid,
GP, or NNGP effects; composition formulas can add corresponding factor terms.
Colon interactions allow covariates to multiply latent effects. See
\link{magCompTerms} for syntax and identification constraints.

The \code{composition_link} control selects probit (the default) or logit. The selected link applies to all
allocation blocks, fitted values, and predictions. Ordered composition
coefficients can use RW1 or RW2 smoothing.

Each magnitude process term has an estimated variance. Composition factors
have unit variance, with positive estimated anchor-diagonal loadings and zeros
above the anchor diagonal. Covariance parameters are estimated unless their starting entry uses
\code{fixed(value)}; see \link{magCompCovModels} for family-specific defaults. Negative-binomial size is estimated by default
under a Gamma(shape = 1, rate = 0.01) prior when the magnitude formula has no
latent process terms. With magnitude iid or spatial terms, size is fixed at
\code{control$nb_size} (default 1000). See \code{\link{magCompPriors}} to change
the prior or explicitly fix size in a regression-only magnitude model.

Use \code{\link{magCompSettings}} to inspect the resolved model, priors, and
starts. The technical guide linked from \link{magComp-package} gives the model
and sampler details; \link{magCompNotation} connects its symbols to software
parameter names.
}
\value{
An object of class \code{magComp}. Use \code{\link{samples}} for
posterior parameter arrays, \code{\link[=summary.magComp]{summary}} for parameter summaries,
and \code{\link[=fitted.magComp]{fitted}} or
\code{\link[=predict.magComp]{predict}} for response quantities.
}
\seealso{\link{magComp-package}, \code{\link{magCompComposition}},
\code{\link{magCompControl}}}
\examples{
rows <- expand.grid(plot = paste0("p", 1:8),
                    dbh_class = c("small", "large"))
rows$count <- c(8, 5, 9, 6, 3, 7, 4, 5, 2, 4, 1, 3, 6, 2, 5, 4)
rows$area <- 0.05
rows$height <- rep(seq(-1, 1, length.out = 8), 2)
dat <- magCompData(rows, unit = plot, category = dbh_class,
                   tree_count = count, design = fixedArea(area))

# A short run to illustrate the interface, not assess convergence.
set.seed(1)
fit <- magComp(dat, magnitude = ~ height, composition = ~ height,
  control = magCompControl(n_samples = 60, n_burn = 20, verbose = FALSE))
summary(fit)
fitted(fit, type = "density")$density
predict(fit, data.frame(height = c(-0.5, 0.5)), area = 0.05)$count
}
