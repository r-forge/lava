\name{twinsim}
\alias{twinsim}
\title{Simulate twin data}
\usage{
  twinsim(n = 100, k1 = c(), k2 = 1, mu = 0,
    lambda = c(1, 1, 1), randomslope = NULL, type = "ace",
    binary = FALSE, ...)
}
\arguments{
  \item{n}{Number of monozygotic AND dizygotic twin-pairs
  to simulate (i.e. total of 4*n individuals).}

  \item{k1}{Number of covariates (labelled x1,x2,...) of
  type 1. One distinct covariate value for each
  twin/individual.}

  \item{k2}{Number of covariates (labelled g1,g2,...) of
  type 2. One covariate value for each twin pair.}

  \item{mu}{Intercept parameter.}

  \item{lambda}{Standard errors of variance components (in
  the order A,C/D,E)}

  \item{randomslope}{Logical indicating wether to include
  random slopes of the variance components w.r.t.
  x1,x2,...}

  \item{type}{Type of model (subset of "acde").}

  \item{binary}{If TRUE data is simulated from a
  probit-model (liability-model)}

  \item{...}{Additional arguments parsed on to lower-level
  functions}
}
\value{
  A \code{list} with the following elements
  \item{data}{Data in long format (one row for each
  individual.} \item{model}{The multigroup structural
  equation model.} \item{wide}{A list of data.frames (MZ
  data and DZ data) in wide format (one row for each pair)
  eady for a structural equation model.}
}
\description{
  Simulate twin data from a linear normal ACE/ADE/AE model.
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{twinlm}}
}
\keyword{models}
\keyword{regression}

