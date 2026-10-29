\name{magCompComposition}
\alias{magCompComposition}
\title{Specify composition allocation blocks}
\usage{magCompComposition(...)}
\arguments{
  \item{...}{Two-sided formulas naming the response category on the left.
  Conditional blocks use \code{child | parent ~ ...}.}
}
\description{
Creates a composition specification with one formula per allocation block.
For two categories, the first formula allocates among parent categories and
the second allocates among child categories within each parent.
}
\details{
The category variables are chosen by the user and need not be species or
diameter classes. Either the parent or child can be ordered or unordered.
Smoothing uses the ordering recorded by \code{magCompData()}, not the
category's name or position.

For example, use \code{species ~ height} and
\code{dbh_class | species ~ height + cover}. The two blocks may use different
covariates and different factor terms. Their category names must match those
supplied to \code{\link{magCompData}}, in parent-then-child order.
To reverse the conditioning, reverse the columns in \code{categories} and
specify the corresponding formulas. This changes the model specification;
structured models and their priors need not give equivalent fits in both orders.

A leading block with \eqn{J_A} categories has \eqn{J_A - 1} coefficient rows.
A conditional block with \eqn{J_B} child categories has
\eqn{J_A(J_B - 1)} rows: a separate child allocation for each parent.
Within a conditional block, factor scores are shared across parents, with
separate loading rows for each parent/child stick. Each block can have its
own factor rank and anchors.

Right-hand sides use ordinary regression formulas and the markers described
in \link{magCompTerms}. For a single category, a two-sided formula such as
\code{dbh_class ~ height} is equivalent to supplying \code{composition = ~ height}
to \code{magComp()}.
}
\value{A \code{magCompComposition} object to pass to \code{magComp()}.}
\seealso{\code{\link{magComp}}, \code{\link{magCompSmooth}};
the technical guide's sections on composition in \link{magComp-package}.}
\examples{
magCompComposition(
  species ~ height + factor_nngp(x, y, K = 1),
  dbh_class | species ~ height + cover
)
}
