\name{spatialDiagnostics}
\alias{spatialDiagnostics}
\alias{spatialDiagnostics.magComp}
\title{Spatial diagnostics for magComp fits}
\usage{
spatialDiagnostics(object, ...)

\method{spatialDiagnostics}{magComp}(object, probs = c(0.025, 0.5, 0.975), ...)
}
\arguments{
  \item{object}{A fitted \code{magComp} object.}
  \item{probs}{Posterior quantile probabilities.}
  \item{...}{Currently unused.}
}
\description{
Summarizes monitored spatial scalar parameters and simple posterior
dependence diagnostics from a fitted \code{magComp} object. The diagnostics
are intended to flag weak range identification and confounding between spatial
range parameters, process variances, and composition factor loading
magnitudes.
}
\details{
For each spatial or space-time term, the summary reports covariance parameters
under their registry names. Magnitude terms also report process variance;
composition factors report loading-column norms. Acceptance rates refer to
one joint covariance block per magnitude process or composition factor.
Posterior correlations compare each estimated covariance parameter with the
magnitude variance or composition loading norm. Fixed values are summarized
without acceptance or correlation rows for an entirely fixed block.

}
\value{
A list of data frames with elements \code{summary}, \code{acceptance}, and
\code{correlations}. The object has class \code{magCompSpatialDiagnostics}
and a compact print method.
}
\seealso{\link{magCompNotation} for the full parameter mapping.}
