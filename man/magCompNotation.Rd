\name{magCompNotation}
\alias{magCompNotation}
\title{Model notation and software parameter names}
\description{
Connects the mathematical notation used in the manuscript and technical
document to the package's prior settings, starting values, and sample labels.
This is a help topic, not a function. The technical guide is linked from
\link{magComp-package}.
}
\details{
The magnitude model uses regression coefficients \eqn{\beta}, iid group
effects \eqn{\delta_g}, and spatial effects \eqn{v(s)}. Composition uses
intercepts \eqn{\alpha_{0j}} and slope vectors \eqn{\alpha_j}. Iid factor
scores \eqn{u_h} have loadings \eqn{\gamma_j}; spatial factor scores
\eqn{w(s)} have loadings \eqn{\lambda_j}. Additional terms and allocation
blocks have their own parameter indices. Both factor types retain unit scale
and the specified loading-anchor constraints. The index \eqn{\ell} denotes
a retained posterior draw in the technical document.

The descriptive software names map to those symbols as follows:
\describe{
  \item{Magnitude coefficients \eqn{\beta}}{
  Within \code{starting$magnitude}, \code{coefficients} is a vector named
  by design columns, including \code{"(Intercept)"}. In the magnitude priors,
  \code{coefficient_variance} is \eqn{\sigma_\beta^2}.}
  \item{Composition coefficients \eqn{\alpha_{0j},\alpha_j}}{
  Each composition block has a \code{coefficients} matrix with one row
  per coefficient stick and one column per design column. Its intercept column
  contains \eqn{\alpha_{0j}}. Block prior settings \code{intercept_variance}
  and \code{slope_variance} are \eqn{\sigma_{\alpha,0}^2} and
  \eqn{\sigma_\alpha^2}, with defaults inherited from composition settings.}
  \item{Magnitude effects \eqn{\delta_{gt},v_t(s_i)}}{
  Within \code{starting$magnitude$terms[[term]]}, \code{effects} contains
  the group or spatial effects, \code{variance} is \eqn{\sigma_{\delta,t}^2}
  or \eqn{\sigma_{v,t}^2}, and an estimated spatial \code{covariance$phi} is
  \eqn{\phi_{v,t}}. The corresponding prior term uses \code{variance} for
  inverse-gamma shape and scale and \code{covariance$phi} for uniform bounds.}
  \item{Composition loadings and scores}{
  Within a composition block's \code{terms[[term]]} entry,
  \code{loadings} contains \eqn{\gamma} for iid terms or \eqn{\lambda} for
  spatial terms; \code{scores} contains \eqn{u} or \eqn{w}, respectively.
  An estimated spatial \code{covariance$phi} has one value \eqn{\phi_{w,t,k}} per
  factor column. For the term, \code{loading_variance} specifies
  \eqn{c_{\gamma,t}} or \eqn{c_{\lambda,t}}; \code{covariance$phi} gives uniform bounds.}
  \item{Negative-binomial size \eqn{r_M}}{
  The magnitude prior entry \code{size} gives Gamma shape and rate when size
  is estimated. The corresponding starting entry gives its initial value.
  By default size is estimated with Gamma(shape = 1, rate = 0.01) when there
  are no magnitude process terms. With those terms, or an explicit
  \code{size = NULL} prior setting, \code{control$nb_size} fixes its value.}
}
Space-time effects use \eqn{v(s,t)} and \eqn{w(s,t)}. Covariance entries
use the names listed by \code{\link{magCompCovModels}}, for example
\code{covariance$lambda} for temporal exponential decay.

Coefficient and loading prior variances and smoothing variances are fixed
hyperparameters. Magnitude process variances are estimated parameters.
For a composition coefficient vector \eqn{\boldsymbol{\theta}_q},
with \eqn{q=0} for intercepts and \eqn{q>0} for slopes, smoothing uses
\eqn{\tau_q^2}: \code{intercept_smoothing_variance} sets the intercept
value and \code{slope_smoothing_variance} sets the common slope value within
a block.

\code{magCompSmooth(intercept = , slope = )} controls composition intercept
and slope smoothing. \code{\link{samples}} returns posterior arrays
using the same descriptive branches and entry names as starting values.
Fixed size and covariance values are omitted.

\code{\link{magCompSettings}} and \code{\link[=summary.magComp]{summary}} include a
\code{symbol} column using portable ASCII notation, such as \code{beta},
\code{sigma_delta^2}, and \code{phi_w}. Summary columns separately identify
\code{component}, \code{block}, \code{term}, \code{parameter}, and
\code{index}. The index identifies a coefficient, stick/covariate pair,
stick/factor pair, or factor as appropriate. Parameter names are
\code{coefficients}, \code{size}, \code{variance}, \code{loadings}, and
covariance parameter names such as \code{phi}, \code{lambda}, and \code{nu}.
Traceplot labels join these identifying columns and display
the corresponding symbol.


}
\seealso{
\code{\link{magCompPriors}}, \code{\link{samples}},
\code{\link[=summary.magComp]{summary}}
}
