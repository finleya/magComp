\name{ppcSummary}
\alias{ppcSummary}
\alias{ppcSummary.magComp}
\title{Posterior predictive check summaries}
\usage{
ppcSummary(object, ...)

\method{ppcSummary}{magComp}(object,
                             type = NULL,
                             group = NULL,
                             probs = c(0.025, 0.5, 0.975), ...)
}
\arguments{
  \item{object}{A fitted \code{magComp} object.}
  \item{type}{Summary scale: total count by sampling unit, total count by
  category, sampling-unit/category cell count, total count by group, or
  group/category count. The default is \code{NULL}: use \code{"unit_total"} when magnitude is fitted.
  Composition alone uses \code{"unit_category"}.}
  \item{group}{Grouping variable required for \code{"group_total"} or
  \code{"group_category"}. This can be a
  column name in the unit-level data or a vector with one value per sampling
  unit.}
  \item{probs}{Posterior predictive quantile probabilities.}
  \item{...}{Must be empty. Unknown arguments are rejected.}
}
\description{
Compares observed counts with posterior predictive count samples from a
fitted \code{magComp} model. Replicates are obtained from
\preformatted{fitted(object, target = "observation", type = "count",
       summary = "samples")$count
}
These use the original exposures when magnitude is fitted.
\code{fitted()} returns these simulated observations; \code{ppcSummary()}
compares replicate count summaries with the observed counts. Each call generates fresh replicates;
it does not reuse a previously created prediction object or refit the model.
}

\details{
Choose one supported \code{type}. Full models support all five checks.
Magnitude-only models support \code{"unit_total"} and \code{"group_total"},
using negative-binomial replicate totals. Composition-only models support
\code{"unit_category"}, \code{"category_total"}, and \code{"group_category"},
using multinomial allocations conditional on observed totals. Checking those
fixed totals would not assess fit, so total checks are rejected. The printed output states this conditioning. The
\code{count_conditioning} attribute also records it.
Both fixed-area stems and point-sampling basal-area models are supported.
For point sampling, checks compare observed and replicated in-tree tallies
using each point's BAF. Group summaries add the observed counts and replicate
counts from units in each group. These remain tally-count checks, including
when BAF varies between points; they are not expanded stand estimates.

For each check, \eqn{T} is the observed count or sum of counts selected by
\code{type}; \eqn{T^{rep}} is its replicate counterpart. The reported
\code{p_lower} and \code{p_upper} estimate
\deqn{
p_{lower} = \Pr(T^{rep} \le T \mid \mathrm{data}), \qquad
p_{upper} = \Pr(T^{rep} \ge T \mid \mathrm{data}).
}
They are computed as fractions of replicates across all retained draws and
chains. The \code{p_value} column is a two-sided posterior predictive p-value:
\deqn{p_{two} = \min\{1,\, 2\min(p_{lower}, p_{upper})\}.}
Ties contribute to both tails, so their probabilities need not sum to one.

A small \code{p_lower} indicates an unusually low observed count relative to
the replicates; a small \code{p_upper} indicates an unusually high count.
For example, \code{p_upper = 0.02} gives \code{p_value = 0.04} when the upper
tail is the smaller one. A \code{p_value} near one indicates a central observed
count. A reported zero means no simulated replicates reached that tail; it
is not evidence that its probability is exactly zero.

These are Bayesian posterior predictive checks, not probabilities that the
model is correct or incorrect. The same data inform the fit and the check;
posterior predictive p-values are not generally uniform under the fitted
model (Gelman, 2013), and a threshold such as 0.05 does not have the usual
hypothesis-test calibration. Use small values to flag discrepancies for
inspection, including their size and direction. Large values do not establish
overall model adequacy. These checks concern the fitting data, not
cross-validation performance.

The requested quantiles describe replicate observations, including sampling
variation. They are not credible intervals for expected counts. Separately
summarizing draws returned by \code{fitted()} gives the same kinds of
quantities, with Monte Carlo differences because separate calls simulate
separate replicates.

Replicate simulation uses R's current random-number state. Use
\code{set.seed()} externally for reproducibility.
}
\value{
A data frame with identification columns \code{type}, \code{group},
\code{unit}, and \code{category}, plus:
\describe{
  \item{observed}{The observed count or sum of counts.}
  \item{mean}{The posterior predictive replicate mean.}
  \item{q...}{Requested replicate quantiles; defaults are \code{q2_5},
  \code{q50}, and \code{q97_5}.}
  \item{p_lower, p_upper}{Fractions of replicates at or below, and at or
  above, the observed value.}
  \item{p_value}{Twice the smaller tail fraction, capped at one.}
  \item{n_draws}{The number of replicate draws pooled across chains.}
}
}
\references{
Gelman, A. (2013). Two simple examples for understanding posterior p-values
whose distributions are far from uniform. \emph{Electronic Journal of
Statistics}, 7, 2595--2602. \doi{10.1214/13-EJS854}.
}
\seealso{\code{\link[=fitted.magComp]{fitted}}, \code{\link{magCompCV}}}
