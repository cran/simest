\name{solve.pentadiag}
\title{Pentadiagonal Linear Solver}
\alias{solve_pentadiag}
%alias{solve.pentadiag} % --> ./simest-deprecated.Rd
\description{
  A function to solve pentadiagonal system of linear equations.
}
\usage{
solve_pentadiag(a, b)
## method  solve.pentadiag(a, b, ...)  is __deprecated__ now (Nov.18 2025)
}
\arguments{
  \item{a}{a numeric square matrix with pentadiagonal rows. The function does NOT check for pentadiagonal matrix.}
  \item{b}{a numeric vector of the same length as nrows(a). This argument cannot be a matrix.}
}
\details{
  This function is written mainly for use in this package. It may not be the most efficient code.

  Originally it was documented (and declared) as an S3 method for \code{\link{solve}()},
  and class \code{pentadiagonal} and as function named \code{solve.pentadiag},
  but such classed objects were never used.
}
\value{
  A vector containing the solution.
}
\author{Arun Kumar Kuchibhotla}
\examples{
A <- matrix(c(2,1,1,0,0,
              1,2,1,1,0,
              1,1,2,1,1,
              0,1,1,2,1,
              0,0,1,1,2), nrow = 5)
b <- 1:5
tools::assertWarning(tmp <- solve.pentadiag(A, b), verbose=TRUE) # deprecated
tmp # 0.5 0.75 -0.75 0.75 2.5
stopifnot( identical(tmp,   solve_pentadiag(A, b)))
}
\keyword{Pentadiagonal}
