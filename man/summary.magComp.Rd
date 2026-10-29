\name{summary.magComp}
\alias{summary.magComp}
\alias{summary.magCompSamples}
\title{Posterior parameter summaries}
\usage{
\method{summary}{magComp}(object, probs = c(0.025, 0.5, 0.975), ...)

\method{summary}{magCompSamples}(object, probs = c(0.025, 0.5, 0.975), ...)
}
\arguments{
  \item{object}{A fitted \code{magComp} object or an object returned by \code{samples()}.}
  \item{probs}{Nonempty numeric vector of posterior quantile probabilities
    strictly between zero and one.}
  \item{...}{For a fitted model, arguments passed to \code{samples()}, including
  \code{type}, \code{start}, \code{thin}, \code{chains}, \code{combine_chains},
  and \code{group}. Unused for an already extracted sample object.}
}
\description{
Summarizes posterior draws. \code{summary(fit, ...)} calls
\code{samples(fit, ...)} and then summarizes that sample object. The default
reports scalar parameters; use \code{type = "effects"} or
\code{type = "combined"} for individual effects or combined coefficients. Multiple-chain fits include the classical unsplit Gelman-Rubin
\eqn{\widehat{R}} diagnostic; single-chain fits report \code{NA} for
\eqn{\widehat{R}}. Effective sample size is computed from within-chain
autocorrelation and summed across chains, without a between-chain correction.
The R-hat calculation is bounded below by one and does not split or rank-normalize
chains. See the technical guide's Diagnostics section in \link{magComp-package}
for the formulas.
}
\value{
A data frame of class \code{magCompSummary}, identifying \code{component}, \code{block}, \code{term},
\code{parameter}, \code{index}, and manuscript \code{symbol}, followed by
posterior means, standard deviations, requested quantiles,
\eqn{\widehat{R}}, effective sample size, retained samples per chain, and
number of chains.
}
\details{
Parameter names match the public sample entries: \code{coefficients},
\code{size}, \code{variance}, \code{loadings}, and covariance parameter
names such as \code{phi}, \code{lambda}, and \code{nu}. The
\code{index} column identifies the relevant design column, stick/covariate
pair, stick/factor pair, or factor. Component, block, and term identities stay
the same as in \code{\link{magCompSettings}} and \code{\link{samples}}.
Symbols use portable ASCII notation. Fixed size and covariance values are excluded;
positive anchor-diagonal loadings are included as estimated parameters.
Fixed upper-anchor zeros are omitted from the scalar summary. Latent effects
and scores are available through \code{summary(fit, type = "effects")}.
Combined coefficients also identify their fitted unit or group, coefficient,
allocation row, and parent context. R-hat and ESS use the selected draws
within their original chains, even when storage has been pooled. One-chain
selections report \code{NA} for R-hat.
}
\seealso{\link[=summary.magCompPrediction]{Fitted and prediction summaries}, \link{magCompNotation}, \code{\link[=plot.magComp]{plot}},
\code{\link{samples}}}
