% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_funcs.R
\name{dual1d_Dmat}
\alias{dual1d_Dmat}
\title{Utility function for creating the penalty matrix D that corresponds to a chain graph of length n}
\usage{
dual1d_Dmat(m)
}
\arguments{
\item{n:}{length of the chain graph}
}
\value{
An (n-1)x(n) penalty matrix D
}
\description{
Utility function for creating the penalty matrix D that corresponds to a chain graph of length n
}
\examples{
n <- 10
D_1d <- dual1d_Dmat(n)
}
\keyword{internal}
