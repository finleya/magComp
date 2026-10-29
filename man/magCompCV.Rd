\name{magCompCV}
\alias{magCompCV}
\alias{plot.magCompCV}
\title{Cross-validation at held-out sampling units}
\usage{
magCompCV(data, family = fixedAreaStems(), magnitude = ~1,
  composition = ~1, smooth = magCompSmooth(), priors = magCompPriors(),
  starting = list(), control = magCompControl(), tuning = list(), n_fold = 10,
  fold_id = NULL,
  type = NULL,
  target = c("expected", "observation"),
  summary = c("mean", "samples", "quantile"),
  probs = c(0.025, 0.5, 0.975), keep_fits = FALSE)

\method{plot}{magCompCV}(x, type = x$type[1L],
  level = c("cell", "total"), file = NULL, main = NULL, ...)
}
\arguments{
  \item{data}{A data object returned by \code{\link{magCompData}}.}
  \item{family}{A model family: \code{\link{fixedAreaStems}} or
  \code{\link{pointBasalArea}}. Both support expected and observation targets.}
  \item{magnitude}{Magnitude formula, or \code{NULL} to fit composition only.}
  \item{composition}{Composition model formula or
  \code{\link{magCompComposition}} specification passed to
  \code{\link{magComp}}, or \code{NULL} to fit magnitude only.}
  \item{smooth}{Composition smoothing specification.}
  \item{priors}{Prior specification.}
  \item{starting}{Partial or complete starting values specified for the full
  data. Supplied magnitude effects and composition scores are matched to the
  full-data template and subset to each training fold's groups or locations.
  Unspecified starts are computed using only that fold's training data.
  Supplied loadings must satisfy each fold's anchor constraints; use explicit
  anchors when reusing full-data loading starts.}
  \item{control}{MCMC control passed to each fold fit. Sampling is enabled.}
  \item{tuning}{Covariance proposal settings passed to each fold fit; see \link{magCompCovModels}.}
  \item{n_fold}{Number of folds when \code{fold_id} is not supplied. Use
  \code{n_fold = data$n_unit} for leave-one-out cross-validation.}
  \item{fold_id}{Optional vector of fold labels named by sampling unit,
  identifying at least two folds.}
  \item{type}{A nonempty, unique subset of \code{"magnitude"},
  \code{"proportion"}, \code{"density"}, and \code{"count"}, restricted to the
  quantities supported by the fitted components. Default \code{NULL} selects
  the same quantities as \code{\link[=predict.magComp]{predict}}. For plotting,
  select one returned quantity; the first requested quantity is the default.}
  \item{target}{\code{"expected"} returns expected quantities;
  \code{"observation"} simulates observations using the fitted likelihood.
  Meanings and defaults match \code{\link[=predict.magComp]{predict}}.}
  \item{summary}{Return means, aligned retained samples, or quantiles, with
  the same dimensions and behavior as \code{\link[=predict.magComp]{predict}}.}
  \item{probs}{One or more probabilities strictly between zero and one,
  used for \code{summary = "quantile"}. For example, \code{probs = 0.5}
  requests only the median. Also used for plot intervals when samples were
  retained.}
  \item{keep_fits}{Logical; if \code{TRUE}, retain each fold's fitted
  \code{magComp} object.}
  \item{x}{A \code{magCompCV} object.}
  \item{level}{Plot category values (\code{"cell"}) or within-unit totals
  across categories (\code{"total"}). Magnitude already represents overall
  density and has the same plot at either level. Proportions support only
  \code{"cell"}. Composition-only CV has no total-level plot because unit totals
  are conditioned on.}
  \item{file}{Optional PNG output path.}
  \item{main}{Optional plot title.}
  \item{...}{Additional graphical arguments passed to \code{plot()}.}
}
\description{
Fits a \code{magComp} model repeatedly while holding out complete sampling
units, retaining the original category levels. Each fold predicts all held-out
units together through \code{\link[=predict.magComp]{predict}}, using their covariates,
group labels, coordinates, and original sampling exposures. Prediction options
have the same defaults and meanings as \code{predict()}.
}
\details{
Standalone components are supported. Magnitude-only CV compares total counts
or overall densities. Composition-only CV predicts category probabilities;
conditional expected counts and simulated allocations use held-out observed
totals. These assess allocation given the total, not abundance prediction.
The result's \code{count_conditioning} is \code{"observed_total"} in this
case, otherwise \code{"none"}. Conditional count/density totals and total-level
plots are omitted because totals are fixed by construction. Composition-only
probability or count predictions need no areas; conditional densities do.

Folds are fitted sequentially. \code{control$n_cores} parallelizes chains
within a fit, not folds. With no supplied \code{fold_id}, units are randomly
assigned to approximately equal-sized folds.

Each fold uses one joint prediction call, including every fitted model term.
All requested quantities and within-unit totals use the same draws. Totals
are summed across categories before taking means or quantiles. Fold identities
are retained: posterior draws from separately fitted folds do not constitute
joint draws across all held-out units.

Observed counterparts are category counts, category counts divided by exposure
(density), total count divided by exposure (magnitude), and category count
fractions (proportion). Observed proportions are \code{NA} for zero observed
totals. For fixed-area stems, exposure is plot area; for point basal area it is
\eqn{1/\mathrm{BAF}}. CV supplies each held-out unit's original exposure;
there is no user-supplied replacement \code{area} or \code{newdata} argument.

The definitions and zero-total handling in \link{magCompPrediction} apply
to all returned predictions.

Expected-target intervals describe uncertainty in expected quantities.
Observation-target intervals also include observation variation. CV returns
predictions and observed values; it does not calculate predictive scores or
coverage statistics.

Plots use means for mean or sample results, and the requested quantile nearest
0.5 for quantile results, with the selected statistic labeled on the axis.
Intervals use the smallest and largest requested probabilities when at least
two are supplied. Mean-only results have no interval bars. Undefined observed
or predicted proportions do not produce plotted pairs.

Serial workflows can be reproduced with external \code{set.seed()}.
Parallel fitting uses \code{chain_seeds} in \code{magCompControl()}.
}
\value{
An object of class \code{magCompCV}, containing:
\describe{
  \item{\code{predictions}}{A list named by fold label, each entry a standard
    \code{\link{magCompPrediction}} result with the requested quantities,
    summary, row labels, and prediction metadata. A single requested type
    still returns that named result.}
  \item{\code{observed}}{A list named by requested type, with matrices in the
    original sampling-unit order. Category quantities have one column per
    category; magnitude has one column named \code{"total"}. Match prediction
    rows using each fold result's \code{metadata$unit_ids}.}
  \item{\code{totals}}{A list named by fold label, each containing within-unit
    category sums for requested \code{count} and \code{density} quantities.
    Each total uses magnitude-style dimensions: unit by sample (and chain)
    for samples, unit by one column for means, or unit by probability for
    quantiles. Magnitude is already overall density; no total proportion is
    returned. Entries are empty lists for composition-only models or if
    neither count nor density is requested.}
  \item{\code{count_conditioning}}{\code{"observed_total"} for composition-only
    assessment, otherwise \code{"none"}.}
  \item{\code{fold_id}}{Fold assignments named by original sampling unit.}
  \item{\code{fits}}{Fold models named by fold label, or \code{NULL} unless
    \code{keep_fits = TRUE}.}
  \item{\code{call}, \code{target}, \code{type}, \code{summary}, \code{probs}}{
    The call and selected prediction options.}
}
}
\examples{
rows <- expand.grid(plot = paste0("p", 1:6), cls = c("small", "large"))
rows$count <- rep(c(1, 3, 2, 4, 2, 1), 2)
rows$area <- rep(seq(0.05, 0.1, length.out = 6), 2)
dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                   design = fixedArea(area))
set.seed(42)
cv <- magCompCV(dat, n_fold = 2, target = "observation", type = "count",
                summary = "samples",
                control = magCompControl(n_samples = 6, n_burn = 2,
                                          verbose = FALSE))
cv$predictions[["1"]]$count
cv$observed$count
cv$totals[["1"]]$count
}
\seealso{\code{\link[=predict.magComp]{predict}}, \code{\link{magCompPrediction}}}
