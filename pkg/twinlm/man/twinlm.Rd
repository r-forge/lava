\name{twinlm}
\alias{twinlm}
\title{Classic twin model for quantitative traits}
\usage{
  twinlm(formula, data, id = "id", type = c("ace"),
    status = "zyg", DZ, twinnum = "twinnum", weight = NULL,
    binary = FALSE, probitscale = 1, keep = weight,
    estimator = "gaussian", ...)
}
\arguments{
  \item{formula}{Formula specifying effects of covariates
  on the response.}

  \item{data}{\code{data.frame} with one observation pr
  row. In addition a column with the zygosity (DZ or MZ
  given as a factor) of each individual much be specified
  as well as a twin id variable giving a unique pair of
  numbers/factors to each twin pair.}

  \item{id}{The name of the column in the dataset
  containing the twin-id variable.}

  \item{type}{Character defining the type of analysis to be
  performed. Should be a subset of "aced" (additive genetic
  factors, common environmental factors, unique
  environmental factors, dominant genetic factors).}

  \item{status}{The name of the column in the dataset
  containing the zygosity variable.}

  \item{DZ}{Character defining the level in the status
  variable corresponding to the dyzogitic twins. If this
  argument is missing, the reference level (i.e. the first
  level) will be interpreted as the dyzogitic twins.}

  \item{twinnum}{The name of the column in the dataset
  numbering the twins (1,2). If it does not exist in
  \code{data} it will automatically be created.}

  \item{weight}{Weight matrix if needed by the chosen
  estimator. For use with Inverse Probability Weights}

  \item{binary}{If \code{TRUE} a liability model is fitted}

  \item{probitscale}{Defines the scale in the
  liability-model}

  \item{keep}{Vector of variables from \code{data} that are
  not specified in \code{formula}, to be added to
  data.frame of the SEM}

  \item{estimator}{Choice of estimator/model.}

  \item{...}{Additional arguments parsed on to lower-level
  functions}
}
\value{
  Returns an object of class \code{twinlm}.
}
\description{
  Fits a classical twin model for quantitative traits.
}
\examples{
## Simulate data
set.seed(1)
d <- twinsim(n=200,k1=c(1,-1),k2=c(),lambda=c(1,1,1))$data
## E(y|z1,z2) = z1 - z2. var(A) = var(C) = var(E) = 1

## E.g to fit the data to an ACE-model without any confounders we simply write
ace <- twinlm(y ~ 1, d, DZ="DZ")
ace
## An AE-model could be fitted as
ae <- twinlm(y ~ 1, d, type="ae", DZ="DZ")
## LRT:
compare(ae,ace)
## AIC
AIC(ae)-AIC(ace)
## To adjust for the covariates we simply alter the formula statement
ace2 <- twinlm(y~z1+z2,d,DZ="DZ")
## Summary/GOF
summary(ace2)
## An interaction could be analyzed as:
ace3 <- twinlm(y ~ z1 + z2 + z1:I(z2<0),d, DZ="DZ")
## Categorical variables are also supported
d2 <- transform(d,z2cat=cut(z2,3,labels=c("Low","Med","High")))
ace4 <- twinlm(y ~ z1 + z2cat, d2, DZ="DZ")
## plot the model structure
\donttest{
plot(ace4)
}
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{twinsim}}
}
\keyword{models}
\keyword{regression}

