\name{magCompCovModels}
\alias{magCompCovModels}
\alias{uniform}
\alias{fixed}
\title{Covariance models and their parameter settings}
\usage{
magCompCovModels()
uniform(lower, upper)
fixed(value)
}
\arguments{
  \item{lower, upper}{Finite bounds for a uniform covariance-parameter prior,
  with \code{lower < upper}. Bounds must respect the parameter's domain.}
  \item{value}{A fixed parameter value, or a named vector of values for factors.}
}
\description{
Lists the covariance families, parameter names, and domains supported by GP
and NNGP terms. The same covariance calculations are used for fitting and
prediction. \code{uniform()} specifies a prior; \code{fixed()} marks a starting
value as fixed throughout fitting.
}
\details{
Specify a family with \code{cov_model} inside \code{gp()}, \code{nngp()},
\code{factor_gp()}, or \code{factor_nngp()}. The default is \code{"exp"}.
For space-time families, the final coordinate is numeric time, for example
\code{nngp(easting, northing, year, cov_model = "sep_exp", m = 15)}.
Time gaps are used as supplied. Let \eqn{h} be spatial distance, \eqn{u} absolute
time separation, and \eqn{d} the number of spatial coordinates. Correlations are:
\describe{
  \item{\code{exp}}{\eqn{\exp(-\phi h)}, with \code{phi > 0}.}
  \item{\code{matern}}{\eqn{(\phi h)^\nu K_\nu(\phi h)/(2^{\nu-1}\Gamma(\nu))},
  with \code{phi > 0} and \code{nu > 0}. At \code{nu = 0.5} this is exponential.}
  \item{\code{sep_exp}}{\eqn{\exp(-\phi h-\lambda u)}, with positive
  \code{phi} and \code{lambda}.}
  \item{\code{multi_res_sep_exp}}{\eqn{\alpha\exp(-\phi_1h-\lambda_1u)
  +(1-\alpha)\exp(-\phi_2h-\lambda_2u)}. The mixing weight \code{alpha} is
  between zero and one; all four decays are positive.}
  \item{\code{gneiting}}{Write \eqn{T=1+a u^{2\alpha}}. Correlation is
  \eqn{T^{-\delta-d/2}\exp(-c h^{2\gamma}/T^{\beta\gamma})}.
  Parameters \code{a} and \code{c} are positive; \code{alpha} and \code{gamma}
  are in \eqn{(0,1]}; \code{beta} is in \eqn{[0,1]}; \code{delta} is nonnegative.}
}
Magnitude process variance is separate from these correlations and has an
inverse-gamma Gibbs update. Composition factors have unit marginal variance.

All estimated covariance parameters have independent bounded uniform priors.
Put overrides in a term's \code{covariance} list under \code{priors},
\code{starting}, or \code{tuning}. The same parameter names are used in
\code{samples()} and \code{summary()}. For composition, a common scalar setting
applies to all factors; a named vector of starts/tuning values or a named list
of priors can override individual factors. A named list of starting values can
mix \code{fixed(value)} and ordinary numeric values across factors.
Omitted entries retain defaults. \code{magCompSettings()} shows the resolved
supports and starts. Fixed parameters do not have priors or proposals.

By default, Matérn smoothness has prior \code{uniform(0.05, 2)}. Gneiting
\code{alpha = fixed(1)}, \code{gamma = fixed(0.5)}, and \code{delta = fixed(0)};
\code{a}, \code{c}, and \code{beta} are estimated. An ordinary starting value or
a prior override releases a default-fixed parameter.

Other default supports are calibrated to spatial distances and the time span.
For the mixture, the default supports of \code{phi_1} and \code{phi_2} meet at
\eqn{30/D}, where \eqn{D} is maximum spatial separation; component 1 describes
longer spatial dependence. Temporal supports are calibrated to dependence up
to twice the observed time span. User overrides are direct parameter bounds,
not range parameters. Overlapping mixture supports can permit label switching.
When a coordinate domain has no positive span, supply the corresponding
support or fix its parameter. The technical guide gives the default formulas.

Covariance parameters are proposed jointly within each magnitude process and
each composition factor, after bounded-logit transformation, with Jacobian
correction. The \code{tuning} entries are initial transformed-scale standard
deviations (default 0.1). \code{magCompControl(metropolis = )} controls adaptation.
}
\value{
\code{magCompCovModels()} returns a named list of registry entries.
\code{uniform()} returns a prior specification and \code{fixed()} a marked
numeric value for a covariance starting-value list.
}
\seealso{\link{magCompTerms}, \link{magCompPriors}, \link{magCompStarting},
\link{magCompSettings}, \link{magCompControl}}
\examples{
magCompCovModels()$sep_exp
~ nngp(easting, northing, year, cov_model = "sep_exp", name = "space_time")
magCompPriors(magnitude = list(terms = list(space_time = list(
  covariance = list(phi = uniform(0.1, 4), lambda = uniform(0.02, 2))))))
# A matching starting-value list for magComp(starting = ...).
list(magnitude = list(terms = list(space_time = list(
  covariance = list(phi = 1, lambda = fixed(0.3))))))
}
