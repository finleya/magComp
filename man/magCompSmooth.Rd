\name{magCompSmooth}
\alias{magCompSmooth}
\title{Smooth composition coefficients across ordered categories}
\usage{
magCompSmooth(intercept = c("none", "rw1", "rw2"),
              slope = c("none", "rw1", "rw2"),
              loading = c("none", "rw1", "rw2"),
              ordered = TRUE, ...)
}
\arguments{
  \item{intercept}{Smoothing for composition intercept coefficients:
  independent (\code{"none"}), first-order (\code{"rw1"}), or second-order
  (\code{"rw2"}) random-walk smoothing.}
  \item{slope}{Smoothing for each composition slope, with the same choices.}
  \item{loading}{Loading smoothing. Only \code{"none"} is implemented.}
  \item{ordered}{Ordering flag stored in this specification. When preparing
  a model, the category ordering recorded by \code{magCompData()} determines
  whether smoothing is permitted. This flag does not order the data.}
  \item{...}{Named block or response overrides. Each is a list containing
  \code{intercept}, \code{slope}, \code{loading}, and/or \code{ordered}, or a
  single method name applied to both intercept and slope. Conditional block
  names have the form \code{"child|parent"}.}
}
\description{
Specifies how regression coefficients vary across ordered composition
categories. Smoothing requires ordered categories in \code{\link{magCompData}}.
Unspecified methods default to independent coefficient priors.
}
\details{
RW1 penalizes adjacent coefficient differences; RW2 penalizes changes in
those differences. In conditional blocks, smoothing acts separately within
each parent's child-category sequence. It does not smooth across parents.

Global settings apply unless overridden by a named block or response entry.
An exact block name takes precedence over its response name. Smoothing
variances are fixed prior settings controlled by \code{\link{magCompPriors}};
they are not sampled parameters. The technical guide's section on regression
coefficients and ordered-category smoothing gives the precision matrices;
see \link{magComp-package}.
}
\value{A \code{magCompSmooth} object to pass as \code{magComp(smooth = )}.}
\seealso{\code{\link{magCompComposition}}, \code{\link{magCompPriors}}}
\examples{
magCompSmooth(intercept = "rw2", slope = "rw2")
magCompSmooth(dbh_class = list(intercept = "rw2", slope = "rw1"))
}
