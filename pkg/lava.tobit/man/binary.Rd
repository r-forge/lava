\name{binary}
\alias{binary}
\alias{binary<-}
\alias{binary.lvm}
\alias{binary<-.lvm}

\title{Defines variables as binary}
\description{
  Defines variables of a latent variable model (\code{lvm}-object from
  the \code{lava} package) as binary.
}
\usage{
\method{binary}{lvm}(x, var=NULL, ...)
\method{binary}{lvm}(x, ...) <- value
}
\arguments{
  \item{x}{\code{lvm}-object}
  \item{var}{Formula or vector of variable names}
  \item{value}{Formula or vector of variable names}
  \item{\dots}{Additional arguments parsed to lower-level functions}
}
\details{
  This function sets the status of variables to 'binary'. For use with
  the \code{estimate} method a probit-link will be assumed. Used with
  the \code{sim} method normal distributed data will be simulated
  followed by thresholding at 0. To simulate data where the dichotomous
  variable has a direct effect on the outcome the \code{distribution}
  method can be used, e.g. \code{distribution(m,~x) <- binomial.lvm("probit")}.
}
\value{
  \code{lvm}-object (or vector of variable names if called without any arguments)
}
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
