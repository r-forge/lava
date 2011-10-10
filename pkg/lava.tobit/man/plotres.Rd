\name{plotres}
\alias{plotres}
\title{Plot distribution of standardized residuals}
\description{
  Plot empirical (KM) and model-specific cumulative distribution function  of standardized residuals  
}
\usage{
plotres(x,var=endogenous(x),
                    ylab="Cumulative Distribution Function",
                    xlab="Standardized residuals",
                    main,
                    ...)}
\arguments{
  \item{x}{Model, \code{lvmfit} object}
  \item{var}{Character vector of (endogenous) variable names}
  \item{ylab}{Label of x-axis}
  \item{xlab}{Label of y-axis}
  \item{main}{Title of plot}
  \item{\dots}{Additional argument}
}
\examples{
\dontrun{

## Simulate data where (y01,y2)
## follows conditional bivariate normal distribution
## given covariate x. Instead of y01 we observe
## right censored version y2
n <- 200
m <- lvm(c(y01,y2) ~ x)
covariance(m) <- y01~y2
set.seed(1)
d <- sim(m,n)
d$cens1 <- rexp(n)
d$status1 <- with(d,y01<cens1)
d$y1 <- with(d, pmin(y01,cens1))

## Estimate model parameters
d$S1 <- with(d, Surv(y1,status1))
m <- lvm(c(S1,y2)~x); covariance(m) <- S1~y2
e <- estimate(m,d,control=list(trace=1))

## Plot cumulative distribution functions
par(mfrow=c(2,2)); plotres(e); plot(e)
}
}
\author{Klaus K. Holst}
\keyword{models}

