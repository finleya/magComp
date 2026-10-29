\name{summary.magCompPrediction}
\alias{summary.magCompPrediction}
\title{Summarize stored fitted, prediction, or stand draws}
\usage{
\method{summary}{magCompPrediction}(object, probs = c(0.025, 0.5, 0.975), ...)
}
\arguments{
  \item{object}{Stored draws returned by \code{fitted()}, \code{predict()},
    or \code{standSummary()} with \code{summary = "samples"}.}
  \item{probs}{Nonempty numeric vector of quantile probabilities strictly
    between zero and one. Defaults give the median and a 95\% interval.}
  \item{...}{Currently unused.}
}
\description{
Calculates means, standard deviations, and quantiles from stored draws,
pooling chains. No new observations or latent effects are simulated.
}
\value{
A data frame with one row per unit, quantity, and category. Its class is
\code{magCompPredictionSummary}. The identifying columns are \code{unit},
\code{quantity}, and \code{category}; magnitude uses category
\code{"total"}. Columns \code{mean}, \code{sd}, and the requested quantiles
(for example, \code{q2_5}, \code{q50}, and \code{q97_5}) summarize each row;
\code{n_draws} counts the draws used. Rows follow the requested quantity order,
then category order, with units in their original order within each category.
The \code{metadata} attribute preserves the source object's metadata,
including measurement units, target, and conditioning information. The
\code{probs} attribute records the summary's quantile probabilities.
}
\details{
Expected-target intervals describe uncertainty in conditional means.
Observation-target intervals also include variation from the observation model. Proportions from observation predictions or stand aggregation are undefined at zero totals:
their summaries use positive-total draws only. If none remain, the summary
statistics are \code{NA} and \code{n_draws} is zero. Counts, totals, densities, and
magnitudes retain zero-valued draws. Standard deviations are \code{NA} when
fewer than two draws contribute.

Objects containing only means or quantiles cannot supply these summaries.
Retain draws with \code{summary = "samples"} in the original call. For stand
results, units are stand IDs and quantities are density, total, or proportion;
the same calculations apply. This method
summarizes each unit and category separately; use \code{standSummary()} or
aggregate stored draws before summarizing to obtain area-level intervals.
}
\examples{
\dontrun{
draws <- predict(fit, newdata, summary = "samples")
summary(draws)
summary(draws, probs = c(0.05, 0.5, 0.95))

replicates <- fitted(fit, target = "observation", summary = "samples")
summary(replicates)
}
}
\seealso{\link{magCompPrediction}, \code{\link[=summary.magComp]{summary}}
for parameter summaries, \code{\link{standSummary}}, \code{\link{ppcSummary}}}
