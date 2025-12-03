\name{simest-deprecated}
\title{Deprecated Functions in Package \pkg{simest}}
\alias{solve.pentadiag} % deprecated, kept visible
%------ PLEASE: one \alias{.} for EACH ! (+ one \usage{} & \arguments{} for all)
\description{
  Deprecated functions in Package \pkg{simest}:
  \describe{
    \item{\code{solve.pentadiag()}: }{Solve pentadiagonal system of linear
      equations, replaced by \code{\link{solve_pentadiag}()}.}
  }
}
\usage{
solve.pentadiag(a, b, ...)  # is __deprecated__ now (Nov.18 2025)
}
\arguments{
  \item{a, b, ...}{see \code{\link{solve_pentadiag}()}; \code{\dots} is unused.}
}
\seealso{
  \code{\link{Deprecated}}
}
\keyword{internal}
