\name{twinlm}
\alias{twinlm}
\alias{print.twinlm}
\alias{summary.twinlm}
\title{Classic twin model for quantitative traits}
\description{
  Fits a classical twin model for quantitative traits.
}
\usage{
twinlm(formula, data, type=c("ace"), twinid="id", status="zyg", DZ, twinnum="twin",weight, keep=NULL,debug=FALSE,estimator="gaussian",...)

}
\arguments{
  \item{formula}{ Formula specifying effects of covariates on the
    response. }
  \item{data}{ \code{data.frame} with one observation pr row. In
    addition a column with the zygosity (DZ or MZ given as a factor) of each individual much be
    specified as well as a twin id variable giving a unique pair of
    numbers/factors to each twin pair.}
  \item{type}{ Character defining the type of analysis to be
    performed. Should be a subset of "aced" (additive genetic factors, common
    environmental factors, unique environmental factors, dominant
    genetic factors).}
  \item{twinid}{The name of the column in the dataset containing the twin-id variable.}
  \item{status}{The name of the column in the dataset containing the
    zygosity variable.}
  \item{DZ}{ Character defining the level in the status variable
    corresponding to the dyzogitic twins. If this argument is missing,
    the reference level (i.e. the first level) will be interpreted as
    the dyzogitic twins.}
  \item{twinnum}{The name of the column in the dataset numbering the
    twins (1,2). If it does not exist in \code{data} it will
    automatically be created.}
  \item{weight}{Weight matrix if needed by the chosen estimator}
  \item{keep}{Vector of variables from \code{data} that are not
specified in \code{formula}, to be added to data.frame}
  \item{debug}{Debug only...}
  \item{estimator}{Choice of estimator/model.}
  \item{\dots}{}
}
\value{
  Returns an object of class \code{twinlm}.
}
\author{Klaus K. Holst}
\seealso{ \code{\link{twinsim}} }
\examples{
## ACE model with confounders (x and z)
MZ<-lvm()
regression(MZ) <- y1 ~ f(a1,lambda[1]) + f(c1,lambda[2]) + f(e1,lambda[3]) + f(x1,beta[1]) + f(z1,beta[2])
regression(MZ) <- y2 ~ f(a2,lambda[1]) + f(c2,lambda[2]) + f(e2,lambda[3]) + f(x2,beta[1]) + f(z2,beta[2])
latent(MZ)<- ~a1+a2+c1+c2+e1+e2
covariance(MZ) <- c(y1,y2) ~ v(0)
intfix(MZ) <- c(y1,y2) ~ mu1
covfix(MZ, latent(MZ), var2=NULL) <- 1
covariance(MZ) <- c1 ~ f(c2,1)
covariance(MZ) <- a1 ~ f(a2,1)
DZ <- MZ
covariance(DZ) <- a1 ~ f(a2,0.5)

## Simulate data
n <- 200
mono <- sim(MZ,n,normal=TRUE)
di <- sim(DZ,n,normal=TRUE)
## Put data on long format (usual format with one row pr individual)
## and add zygosity status, twin-number (1,2), and twin-id variables to data.frame
d <- data.frame(y=c(mono$y1,mono$y2,di$y1,di$y2),
                x=c(mono$x1,mono$x2,di$x1,di$x2),
                z=c(mono$z1,mono$z2,di$z1,di$z2),
                zyg=rep(c("MZ","DZ"),each=2*n), twin=rep(c(1,2,1,2),each=n), id=rep(1:n,4))
d <- subset(d, select=c(y,x,z,zyg,id))


## E.g to fit the data to an ACE-model without any confounders we simply write
ace <- twinlm(y ~ 1, d)
ace
## An AE-model could be fitted as
ae <- twinlm(y ~ 1, d, type="ae")
## LRT:
compare(ae,ace)
## AIC
AIC(ae)-AIC(ace)
## To adjust for the covariates we simply alter the formula statement
ace2 <- twinlm(y~x+z,d)
## Summary/GOF
summary(ace2)
## An interaction could be analyzed as:
ace3 <- twinlm(y ~ x + z + x:I(z<0),d)
## Categorical variables are also supported
d2 <- transform(d,zcat=cut(z,3,labels=c("Low","Med","High")))
ace4 <- twinlm(y ~ x + zcat,d2)
## plot the model structure
plot(ace4)

}
\keyword{models}
\keyword{regression}
