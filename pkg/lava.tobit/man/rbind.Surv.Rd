\name{rbind.Surv}
\alias{rbind.Surv}
\title{Appending \code{Surv} objects}
\description{
  \code{rbind} method for \code{Surv} objects
}
\usage{
\method{rbind}{Surv}(...)
}
\arguments{
  \item{...}{\code{Surv} objects}
}
\value{
  \code{Surv} object
}
\examples{
y <- yl <- yr <- rnorm(10)
yl[1:5] <- NA; yr[6:10] <- NA
S1 <- Surv(yl,yr,type="interval2")
S2 <- Surv(y,y>0,type="right")
S3 <- Surv(y,y<0,type="left")

rbind(S1,S1)
rbind(S2,S2)
rbind(S3,S3)
}
\author{Klaus K. Holst}
\keyword{utilities}
