\name{magCompPriors}
\alias{magCompPriors}
\title{Set magnitude and composition priors}
\usage{
magCompPriors(magnitude = list(), composition = list())
}
\arguments{
  \item{magnitude}{Named list overriding magnitude prior defaults; see Details.}
  \item{composition}{Named list overriding composition prior defaults, with optional
  named \code{blocks} and their \code{terms}; see Details.}
}
\description{
Creates prior specifications for \code{\link{magComp}}.
Omitted settings retain automatic defaults. On a prepared model,
\code{\link{magCompSettings}} shows resolved priors and editable starts.
}
\details{
Resolved model settings, starting values, and extracted samples contain only
active components. Supplied starting values or nondefault priors for a
disabled component are rejected by \code{magComp()}.

Magnitude prior settings are:
\describe{
  \item{\code{coefficient_variance}}{Zero-mean normal variance for magnitude
  regression coefficients; default 100. It is independent of composition priors.}
  \item{\code{iid_variance}, \code{spatial_variance}}{Default inverse-gamma
  shape and scale for each magnitude iid or spatial term's variance;
  both default to \code{c(shape = 2, scale = 1)}.}
  \item{\code{size}}{The default \code{"auto"} estimates negative-binomial size
  with \code{c(shape = 1, rate = 0.01)} when the magnitude formula has no iid or
  spatial terms; otherwise size is fixed at \code{control$nb_size} (default 1000).
  Supply a positive Gamma shape-rate pair to change the prior, or \code{NULL}
  to fix size explicitly. Size estimation is not allowed with magnitude
  random-effect terms. Composition factor terms do not impose this restriction.}
  \item{\code{terms}}{List named by formula term names. Each entry can override
  \code{variance} (inverse-gamma shape and scale) and \code{covariance}, a named list of \code{uniform(lower, upper)} priors
  for covariance parameters.}
}
Composition defaults are:
\tabular{lr}{
  Setting \tab Value \cr
  \code{intercept_variance} \tab 9 \cr
  \code{slope_variance} \tab 100 \cr
  \code{intercept_smoothing_variance} \tab 1 \cr
  \code{slope_smoothing_variance} \tab 4 \cr
  \code{loading_variance} \tab 1 \cr
}
Coefficient and unconstrained loading priors are zero-mean normal.
Anchor-diagonal loadings are estimated under the same normal prior restricted
to \eqn{(0,\infty)}, a half-normal prior. \code{loading_variance} is the
underlying normal variance before truncation: one gives half-normal scale one;
four gives scale two. The half-normal variance itself is
\eqn{v(1-2/\pi)}, where \eqn{v} is \code{loading_variance}.
Fixed upper-anchor zeros have no prior.
Smoothing and normal prior variances are fixed hyperparameters, not sampled
variances. Composition factor scores retain unit scale.

The default Gamma size prior has mean 100 and a central 95\% interval of
approximately 2.5 to 369. It is a proper, adjustable prior, not a
noninformative one. Larger size means less extra count variation. Use
\code{magCompSettings(model)} to inspect the resolved prior and estimation
status. When size is estimated, \code{starting$magnitude$size} overrides its
automatic start at \code{control$nb_size}.

A named \code{composition$blocks} list overrides these defaults for individual
allocation blocks. Within each block, a named \code{terms} list can override
\code{loading_variance} and the \code{covariance} prior list for each factor term.
Within a factor term, the same \code{loading_variance} controls its half-normal
diagonals and its unconstrained normal loadings. The diagonal Gibbs updates
need no tuning parameter.
Covariance priors may be shared across factors or specified in named lists by
factor. Conditional block names use \code{"child|parent"}. Magnitude term
priors override magnitude defaults; composition block settings override
composition defaults, and factor-term settings override their block.
\link{magCompCovModels} describes covariance families and automatic supports.
The prior is uniform on each covariance parameter, not on an effective range.

For inverse-gamma priors, the density is proportional to
\eqn{v^{-a-1}\exp(-b/v)}, with shape \eqn{a} and scale \eqn{b}.
Magnitude variance updates are Gibbs steps and need no tuning.
Variance-prior parameters must be positive and finite. Covariance bounds must
be finite, increasing, and within their parameter domains. Unknown names,
unused covariance priors for fixed parameters or iid terms,
and incompatible model term/block names are rejected.
}
\value{
A \code{magCompPriors} object to pass as \code{magComp(priors = )}.
}
\seealso{\code{\link{magCompSmooth}}, \code{\link{magCompSettings}},
\link{magCompNotation}; the technical guide's Prior distributions section
in \link{magComp-package}.}
\examples{
# Estimate size with a different Gamma prior (no magnitude random effects).
magCompPriors(magnitude = list(size = c(shape = 2, rate = 0.1)))
# Fix size at the value supplied through magCompControl(nb_size = ...).
magCompPriors(magnitude = list(size = NULL))

magCompPriors(
  magnitude = list(coefficient_variance = 25,
    terms = list(stand = list(variance = c(3, 1)))),
  composition = list(slope_variance = 4,
    blocks = list(dbh = list(intercept_variance = 2,
      terms = list(shape = list(loading_variance = 4)))))
)
}
