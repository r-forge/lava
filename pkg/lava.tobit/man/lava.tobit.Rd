\name{lava.tobit}
\alias{lava.tobit}
\alias{lava.tobit-package}
\docType{package}
\title{
  Estimation and simulation of probit and tobit latent variable models
}
\description{
  Framwork for estimating parameters and simulate data from 
  Latent Variable Models with binary and censored observations. Plugin
  for the \code{lava} package
}
\details{
\tabular{ll}{
Package: \tab lava.tobit \cr
Type: \tab Package \cr
Version: \tab 0.4-1 \cr
Date: \tab 2011-01-21 \cr
License: \tab GPL-3 \cr
LazyLoad: \tab yes \cr
}
}
\author{
Klaus K. Holst
Maintainer: <k.k.holst@biostat.ku.dk>
}
\keyword{ package }
\examples{
m <- lvm(list(c(y,z) ~ x, y~z))
## Simulate 200 observation from path analysis model
## with all slopes and residual variances set to 1 and intercepts 0:
d <- sim(m,200)
## Dichotomize y and introduce censoring on z
d <- transform(d, y=as.factor(y>0), z=Surv(z,z<2))
\dontrun{
e <- estimate(m,d,control=list(trace=1))
effects(e,y~x)
}
}