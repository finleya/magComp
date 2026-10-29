\name{plot.magComp}
\alias{plot.magComp}
\alias{plot.magCompSamples}
\title{Traceplots for magComp fits}
\usage{
\method{plot}{magComp}(x, parameters = NULL, file = NULL,
                      n_col = 2, n_row = 3,
                      width = 8, height = 8, res = 150, ...)

\method{plot}{magCompSamples}(x, parameters = NULL, file = NULL,
                             n_col = 2, n_row = 3,
                             width = 8, height = 8, res = 150, ...)
}
\arguments{
  \item{x}{A fitted \code{magComp} object or extracted \code{magCompSamples} object.}
  \item{parameters}{Character vector selecting parameter names, indices, or full
  identifiers. Examples are \code{"variance"}, \code{"(Intercept)"}, and
  \code{"magnitude:stand:variance"}. Full identifiers join nonempty component,
  block, term, parameter, and index fields with colons.}
  \item{file}{Optional output file ending in \code{.pdf} or \code{.png}.}
  \item{n_col}{Number of plot columns per page.}
  \item{n_row}{Number of plot rows per page.}
  \item{width}{Graphics device width in inches.}
  \item{height}{Graphics device height in inches.}
  \item{res}{PNG resolution.}
  \item{...}{For a fitted model, arguments passed to \code{samples()}.
  Unused for a sample object.}
}
\description{
\code{plot(fit)} draws base-R traceplots for monitored scalar parameters in a fitted
\code{magComp} object. PDF output can contain multiple pages; PNG output is
limited to the first page.
Titles include the manuscript symbol alongside the descriptive identifier.
The horizontal axis uses original MCMC iteration numbers. Pooled sample
objects are separated into their original chains for plotting.
\code{plot(fit, ...)} delegates to \code{plot(samples(fit, ...))}; use
\code{type = "effects"} or \code{type = "combined"} to examine those quantities.
}
\value{
Invisibly returns the output file path when \code{file} is supplied, otherwise
the plotted parameter labels.
}
\seealso{\code{\link[=summary.magComp]{summary}} for the monitored parameters;
\code{\link{samples}} for selection and combined coefficients.}
