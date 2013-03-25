\name{binomial.twostage}
\alias{binomial.twostage}
\title{Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data
using marginals that are on logistic form.
If clusters contain more than two times, the algoritm uses a compososite likelihood
based on the pairwise bivariate models.}
\usage{
  binomial.twostage(margbin, data = sys.parent(),
    score.method = "nlminb", Nit = 60, detail = 0,
    clusters = NULL, silent = 1, weights = NULL,
    control = list(), theta = NULL, theta.des = NULL,
    var.link = 1, iid = 1, step = 0.5, notaylor = 1,
    model = "plackett", marginal.p = NULL, strata = NULL,
    max.clust = NULL, se.clusters = NULL, numDeriv = 0)
}
\arguments{
  \item{margbin}{Marginal binomial model}

  \item{data}{data frame}

  \item{score.method}{Scoring method}

  \item{Nit}{Number of iterations}

  \item{detail}{Detail}

  \item{clusters}{Cluster variable}

  \item{silent}{Debug information}

  \item{weights}{Weights for log-likelihood, can be used
  for each type of outcome in 2x2 tables.}

  \item{control}{Optimization arguments}

  \item{theta}{Starting values for variance components}

  \item{theta.des}{Variance component design}

  \item{var.link}{Link function for variance}

  \item{iid}{Calculate i.i.d. decomposition}

  \item{step}{Step size}

  \item{notaylor}{Taylor expansion}

  \item{model}{model}

  \item{marginal.p}{vector of marginal probabilities}

  \item{strata}{strata for fitting}

  \item{max.clust}{}

  \item{se.clusters}{clusters for iid decomposition for
  roubst standard errors}

  \item{numDeriv}{uses Fisher scoring aprox of second
  derivative if 0, otherwise numerical derivatives}
}
\description{
  The reported standard errors are based on the estimated
  information from the likelihood assuming that the
  marginals are known.
}
\examples{
data(twinstut)
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
margbin <- glm(stutter~factor(sex)+age,data=twinstut,family=binomial())
bin <- binomial.twostage(margbin,data=twinstut,
		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
	                 score.method="fisher.scoring")
summary(bin)

twinstut$cage <- scale(twinstut$age)
theta.des <- model.matrix( ~-1+factor(zyg)+cage,data=twinstut)
bina <- binomial.twostage(margbin,data=twinstut,
		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
	                 score.method="fisher.scoring")
summary(bina)

theta.des <- model.matrix( ~-1+factor(zyg)+factor(zyg)*cage,data=twinstut)
bina <- binomial.twostage(margbin,data=twinstut,
		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
	                 score.method="fisher.scoring")
summary(bina)

twinstut$binstut <- (twinstut$stutter=="yes")*1
out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut,response="binstut",id="tvparnr",
			      theta.formula=~-1+factor(zyg1),score.method="fisher.scoring")
summary(out)

desfs<-function(x,num1="zyg1",namesdes=c("mz","dz","os")) c(x[num1]=="dz",x[num1]=="mz",x[num1]=="os")*1
out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
      data=twinstut,response="binstut",id="tvparnr",type="zyg",
      score.method="fisher.scoring", theta.formula=desfs,desnames=c("mz","dz","os"))
summary(out3)
}
\author{
  Thomas Scheike
}
\references{
  Two-stage binomial modelling
}
\keyword{binomial}
\keyword{regression}
