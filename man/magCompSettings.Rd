\name{magCompSettings}
\alias{magCompSettings}
\title{Inspect resolved model parameters and editable settings}
\usage{
magCompSettings(object)
}
\arguments{
  \item{object}{A prepared or fitted \code{\link{magComp}} model.}
}
\description{
Shows which parameters are estimated, fixed, or constrained; their resolved
priors and initial values; and the mapping from term names to formula terms.
Returns editable prior, starting-value, and tuning templates without running the sampler.
}
\details{
Prepare a model with \code{magComp()} and
\code{control = magCompControl(run_sampler = FALSE)}. Inspect \code{magCompSettings(model)}, edit
the returned \code{priors} and/or \code{starting} lists, and pass them to a new
\code{magComp()} call with the same formulas and data.

The returned prior template explicitly records resolved block/term settings.
Editing a global default does not override an explicit block/term entry;
edit the relevant entry or remove it to inherit a default. The initial values
shown for fitted objects are their original starts, not posterior summaries.
Use \code{\link{magCompStarting}(fit, use = "last")} or \code{use = "mean"}
to extract starts from a previous fit.

Constraints on factor scales and loadings, fixed smoothing variances, and
the magnitude-size estimation restriction are displayed. Priors and starts
can be overridden selectively. Unsupported names and dimensions are rejected.
}
\value{
A \code{magCompSettings} object containing \code{parameters} (a data frame of
parameter settings, including a manuscript \code{symbol} column in ASCII),
\code{terms} (formula term names and labels),
\code{priors} (a resolved \code{magCompPriors} object), and \code{starting}
(a complete initial-value template), plus \code{tuning} (initial covariance
proposal standard deviations). Printing shows the two tables.
}
\seealso{\link{magCompNotation} for mathematical symbols and software names.}
\examples{
rows <- expand.grid(plot = paste0("p", 1:8), cls = c("a", "b"))
rows$count <- rep(1:8, 2)
rows$area <- 1
rows$stand <- rep(rep(c("s1", "s2"), each = 4), 2)
dat <- magCompData(rows, unit = plot, category = cls,
                   tree_count = count, design = fixedArea(area))
model <- magComp(dat, magnitude = ~ iid(stand, name = "stand"),
                 control = magCompControl(run_sampler = FALSE))
settings <- magCompSettings(model)
settings
settings$priors$magnitude$terms$stand$variance <- c(shape = 3, scale = 1)
settings$starting$magnitude$terms$stand$variance <- 0.5
fit <- magComp(dat, magnitude = ~ iid(stand, name = "stand"),
                priors = settings$priors, starting = settings$starting,
                control = magCompControl(n_samples = 10, n_burn = 5,
                  chain_seeds = 7, verbose = FALSE))
magCompStarting(fit, use = "last")
}
