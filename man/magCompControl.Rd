\name{magCompControl}
\alias{magCompControl}
\title{MCMC control for magComp}
\usage{
magCompControl(n_samples = 2000, n_burn = 1000, n_thin = 1,
               n_chains = 1, n_cores = 1,
               chain_seeds = NULL,
               n_report = 100, nb_size = 1000,
               nb_size_tuning = 0.25,
               metropolis = list(),
               composition_link = c("probit", "logit"),
               run_sampler = TRUE, verbose = TRUE)
}
\arguments{
  \item{n_samples}{Total number of MCMC iterations per chain, including
  burn-in. Must be an integer greater than \code{n_burn}.}
  \item{n_burn}{Number of initial iterations discarded from each chain.}
  \item{n_thin}{Thinning interval after burn-in. Must leave at least one
  retained draw per chain.}
  \item{n_chains}{Number of MCMC chains to run from the R side. Chains are
  run sequentially unless \code{n_cores > 1}.}
  \item{n_cores}{Number of R worker processes used for chain-level
  parallelism on Unix-like systems. Values greater than one fall back to
  serial execution on platforms without forked processes.}
  \item{chain_seeds}{Optional integer vector of length \code{n_chains}. When
  supplied, each chain is started after setting the corresponding R seed.
  With \code{NULL}, serial fitting uses R's current RNG state, so an external
  \code{set.seed()} can reproduce a fitting and prediction sequence. Supply
  explicit chain seeds to reproduce parallel fitting under R's default RNG
  settings.}
  \item{n_report}{Progress reporting interval.}
  \item{nb_size}{Fixed negative-binomial size, or its automatic starting value
  when size is estimated. By default size is estimated for magnitude formulas
  without iid or spatial terms and fixed otherwise. Set
  \code{magCompPriors(magnitude = list(size = NULL))} to fix size explicitly;
  setting \code{nb_size} alone does not change estimation status.
  Large fixed values approximate Poisson magnitude counts.}
  \item{nb_size_tuning}{Initial random-walk Metropolis tuning standard
  deviation on the log size scale when \code{nb_size} is estimated.}
  \item{metropolis}{List with \code{adapt} (default \code{TRUE}),
  \code{batch_length} (default 25), and \code{target_accept} (default zero,
  selecting 0.44 for scalar proposals and 0.234 for multivariate proposals).
  An explicit value between zero and one overrides the target. Covariance
  blocks use an initial scale search within burn-in, followed by diminishing
  adaptation of proposal covariance and scale, including after burn-in.
  The separate size update adapts its scalar scale only during burn-in.}
  \item{composition_link}{Composition stick-breaking link: \code{"probit"}
  (the default) or \code{"logit"}. The same link applies to all allocation
  blocks.}
  \item{run_sampler}{Logical; \code{TRUE} (the default) fits the model.
  \code{FALSE} prepares the model and resolves settings without sampling.}
  \item{verbose}{Logical; print sampler progress.}
}
\description{
Creates a control list for \code{magComp()}. Negative-binomial settings apply
only when magnitude is fitted; \code{composition_link} applies only when
composition is fitted.
}
\details{
Each chain runs exactly \code{n_samples} iterations. The number of retained
draws is \code{floor((n_samples - n_burn) / n_thin)}. The first saved draw is
at iteration \code{n_burn + n_thin}. Any final iterations after the last
thinning boundary still run but are not saved. For example, 25,000 total
iterations, 10,000 burn-in iterations, and thinning by 10 retain 1,500 draws
per chain. Output array dimensions and prediction metadata count retained
draws, not total iterations.

The initial covariance scale search uses at most 20 complete batches and never
exceeds burn-in. Each batch contains full sampler sweeps. After that search,
adaptation uses gain \eqn{(b+1)^{-0.8}} for adaptive batch number \eqn{b};
the empirical covariance uses the chain states, including rejected proposals.
There is one block per magnitude process and one per composition factor.
Fixed covariance parameters are excluded. Set \code{adapt = FALSE} to use the
initial proposal throughout; with no burn-in the initial search is skipped.

The composition link applies to all allocation blocks in a fit. Fitted values
and prediction use the matching inverse link when converting composition
linear predictors to stick-breaking probabilities.
}
\value{
A \code{magCompControl} object.
}
\examples{
magCompControl(n_samples = 25000, n_burn = 10000, n_thin = 10)
}
\seealso{\code{\link{magCompPriors}}, \code{\link[=summary.magComp]{summary}};
the technical guide's Sampler algorithm and Numerical implementation sections
in \link{magComp-package}.}
