\name{continuousGaussian}
\alias{continuousGaussian}
\title{Unimplemented continuous response family}
\usage{
continuousGaussian(link = c("log", "identity"),
                   composition = c("logratio", "gaussian_stick"))
}
\arguments{
  \item{link}{Proposed magnitude link.}
  \item{composition}{Proposed composition representation.}
}
\description{
An unused development placeholder. This constructor records a proposed family,
but \code{magComp()} rejects it because no fitting method is implemented.
Use \code{\link{fixedAreaStems}} or \code{\link{pointBasalArea}} for supported
models.
}
\value{A \code{magCompFamily} object that cannot be fitted.}
\keyword{internal}
