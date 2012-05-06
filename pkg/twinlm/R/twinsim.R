##' Simulate twin data from a linear normal ACE/ADE/AE model.
##'
##' @title Simulate twin data
# ##' @return A \code{list} with the following elements
# ##'     \item{data}{Data in long format (one row for each individual.}
# ##'     \item{model}{The multigroup structural equation model.}
# ##'     \item{wide}{A list of data.frames (MZ data and DZ data) in wide format
# ##'     (one row for each pair) eady for a structural equation model.}
##' @author Klaus K. Holst
##' @export
##' @seealso \code{\link{twinlm}}
##' @keywords models
##' @keywords regression
##' @param nMZ Number of monozygotic twin pairs
##' @param nDZ Number of dizygotic twin pairs
##' @param b1 Effect of covariates (labelled x1,x2,...) of type 1. One
##'     distinct covariate value for each twin/individual.
##' @param b2 Effect of covariates (labelled g1,g2,...) of type 2. One
##'     covariate value for each twin pair.
##' @param mu Intercept parameter.
##' @param acde Variance of random effects (in the order A,C,D,E)
##' @param randomslope Logical indicating wether to include random slopes of
##'     the variance components w.r.t. x1,x2,...
##' @param threshold Treshold used to define binary outcome y0
##' @param cens Logical variable indicating whether to censor outcome
##' @param wide Logical indicating if wide data format should be returned
##' @param ... Additional arguments parsed on to lower-level functions
twinsim <- function(nMZ=100,nDZ=nMZ,b1=c(),b2=c(),mu=0,acde=c(1,1,0,1),randomslope=NULL,threshold=0,cens=FALSE,wide=FALSE,...) {
  n <- nMZ+nDZ
  sA <- acde[1]^0.5; sC <- acde[2]^0.5; sD <- acde[3]^0.5; sE <- acde[4]^0.5; 
  A.MZ <- rnorm(nMZ,sd=sA); A.MZ <- cbind(A.MZ,A.MZ)
  D.MZ <- rnorm(nMZ,sd=sD); D.MZ <- cbind(D.MZ,D.MZ)  
  S2 <- matrix(c(0,1,1,0),2)
  A.DZ <- sA*rmvnorm(nDZ,sigma=diag(2)+S2*0.5)
  D.DZ <- sD*rmvnorm(nDZ,sigma=diag(2)+S2*0.25)
  C.MZ <- rnorm(nMZ,sd=sC)
  C.DZ <- rnorm(nDZ,sd=sC)
  yMZ <- mu + A.MZ + cbind(C.MZ,C.MZ) + D.MZ + cbind(rnorm(nMZ,sd=sE),rnorm(nMZ,sd=sE))
  yDZ <- mu + A.DZ + cbind(C.DZ,C.DZ) + D.DZ + cbind(rnorm(nDZ,sd=sE),rnorm(nDZ,sd=sE))
  y <- rbind(yMZ,yDZ)
  if (length(b1)>0) {
    x1 <- rmvnorm(n,rep(0,length(b1)),diag(length(b1)))
    x2 <- rmvnorm(n,rep(0,length(b1)),diag(length(b1)))
    y <- y+cbind(x1%*%b1,x2%*%b1)
  }
  if (length(b2)>0) {
    g <- rmvnorm(n,rep(0,length(b2)),diag(length(b2)))
    ge <- g%*%b2
    y <- y+cbind(ge,ge)
  }
  Cens <- ifelse(cens,rep(Inf,nMZ+nDZ),rnorm(nMZ+nDZ,threshold+1))
  d <- data.frame(id=seq(n),y=y,zyg=c(rep("MZ",nMZ),rep("DZ",nDZ)),
                  cens=Cens)
  vary <- list(c("y1","y2"))
  colnames(d)[2:3] <- vary[[1]]
  if (length(b1)>0) {
    d <- cbind(d,x1=x1,x2=x2);
    vary <- c(vary,lapply(seq(length(b1)),FUN=function(x) paste(c("x1","x2"),x,sep="")))
  }
  if (length(b2)>0) { d <- cbind(d,g=g) }
  colnames(d) <- sub(".","",colnames(d),fixed=TRUE)
  if (wide) return(d)
  dd <- reshape(d,direction="long",varying=vary)
  dd <- transform(dd,y=(y1>threshold & y1<cens)*1,y0=(y1>threshold),status=y1<cens) 
  ## S.MZ <- diag(2)*vE+vC+vA
  ## S.DZ <- diag(2)*(vE+vA) + rho*S2*vA + vC
  ## Mu <- (threshold-mu)
  ## probs <- c(marginal=pnorm(mu,threshold,S.MZ[1]^0.5),concordance=pmvnorm(lower=c(threshold,threshold),mean=c(mu,mu),sigma=S.MZ))
  ## probs <- c(probs,casewise=probs[2]/probs[1])
  ## probs
  return(dd)
}


twinsimobs <- function(n=100,k1=c(),k2=1,mu=0,lambda=c(1,1,0,1),randomslope=NULL,type="acde",binary=FALSE,...) {

  mysep <- ""
  outcomes <- paste("y",mysep,1:2,sep="")
  covars <- covars2 <- NULL
  nk1 <- length(k1)
  nk2 <- length(k2)
  if (nk1>0) {
    covars <- paste("z",1:nk1,sep="")
    covars.1 <- paste(covars,mysep,"1",sep="")
    covars.2 <- paste(covars,mysep,"2",sep="")
  }
  if (nk2>0)
    covars2 <- paste("x",1:nk2 ,sep="")
  model1<-lvm(outcomes,silent=TRUE)
  f1 <- as.formula(paste(outcomes[1]," ~ ."))
  f2 <- as.formula(paste(outcomes[2]," ~ ."))
  regression(model1,silent=TRUE) <- update(f1, . ~ f(a1,lambda[1])+f(c1,lambda[2])+f(d1,lambda[3]) + f(e1,lambda[4]))
  regression(model1,silent=TRUE) <- update(f2, . ~ f(a1,lambda[1])+f(c1,lambda[2])+f(d1,lambda[3]) + f(e2,lambda[4]))    
  latent(model1) <- ~ a1+c1+d1+e1+e2
  if (!is.null(covars))
    for (i in 1:length(covars)) {
      regfix(model1, from=covars.1, to=outcomes[1],silent=TRUE) <- paste("beta[",i,"]",sep="")
      regfix(model1, from=covars.2, to=outcomes[2],silent=TRUE) <- paste("beta[",i,"]",sep="")
    }
  if (!is.null(covars2))
    for (i in 1:length(covars2)) {
      regfix(model1,from=covars2[i],to=outcomes[1],silent=TRUE) <- paste("beta[",i+length(covars2),"]",sep="")
      regfix(model1,from=covars2[i],to=outcomes[2],silent=TRUE) <- paste("beta[",i+length(covars2),"]",sep="")
    }
  ##    model1 <- regression(model1,to=endogenous(model1), from=covars)
  covariance(model1) <- update(f1, . ~  v(0))
  covariance(model1) <- update(f2, . ~  v(0))
  covfix(model1, latent(model1), var2=NULL) <- 1
  ##    covariance(model1) <- c1 ~ f(c2,1)
  ##    covariance(model1) <- a1 ~ f(a2,1)
  intfix(model1,outcomes) <- "mu1"
  model2 <- model1
  
  cancel(model2) <- update(f2, . ~ a1)
  cancel(model2) <- update(f2, . ~ d1)
  regression(model2,silent=TRUE) <- update(f2, . ~ f(a2,lambda[1]))
  regression(model2,silent=TRUE) <- update(f2, . ~ f(d2,lambda[3]))
  covariance(model2) <- a1 ~ f(a2,0.5)
  covariance(model2) <- d1 ~ f(d2,0.25)
  latent(model2) <- ~ a2+d2
  covariance(model2) <- c(a2,d2) ~ v(1)

  if (!is.null(randomslope)) {
    for (i in randomslope) {
      a1var <- paste("a1x",i,sep="")
      a2var <- paste("a2x",i,sep="")
      c1var <- paste("c1x",i,sep="")
      c2var <- paste("c2x",i,sep="")
      regfix(model1, to=a1var, from="a1",silent=TRUE) <- 1
      regfix(model1, to=a2var, from="a1",silent=TRUE) <- 1
      regfix(model1, to=c1var, from="c1",silent=TRUE) <- 1
      regfix(model1, to=c2var, from="c1",silent=TRUE) <- 1
      latent(model1) <- c(a1var,a2var,c1var,c2var)
      covfix(model1,c(a1var,a2var,c1var,c2var),NULL) <- 0
      regfix(model1, to=outcomes[1],from=a1var) <- covars.1[i]
      regfix(model1, to=outcomes[2],from=a2var) <- covars.2[i]
      regfix(model1, to=outcomes[1],from=c1var) <- covars.1[i]
      regfix(model1, to=outcomes[2],from=c2var) <- covars.2[i]
      model1 <- randomslope(model1, covar=covars.1[i],covars.2[i])
      
      regfix(model2, to=a1var, from="a1",silent=TRUE) <- 1
      regfix(model2, to=a2var, from="a2",silent=TRUE) <- 1
      regfix(model2, to=c1var, from="c1",silent=TRUE) <- 1
      regfix(model2, to=c2var, from="c1",silent=TRUE) <- 1
      latent(model2) <- c(a1var,a2var,c1var,c2var)
      covfix(model2,c(a1var,a2var,c1var,c2var),NULL) <- 0
      regfix(model2, to=outcomes[1],from=a1var) <- covars.1[i]
      regfix(model2, to=outcomes[2],from=a2var) <- covars.2[i]
      regfix(model2, to=outcomes[1],from=c1var) <- covars.1[i]
      regfix(model2, to=outcomes[2],from=c2var) <- covars.2[i]
      model2 <- randomslope(model2, covar=covars.1[i],covars.2[i])
      
##      randomslope.lvm(model2) <- c(covars.1[i],covars.2[i])
    }
  }

  
  full <- list(model1,model2)
  ## #######
  type <- tolower(type)
  isA <- length(grep("a",type))>0
  isC <- length(grep("c",type))>0
  isD <- length(grep("d",type))>0
  if (!isA) {
    kill(model1) <- ~ a1 + a2
    kill(model2) <- ~ a1 + a2
  }
  if (!isD) {
    kill(model1) <- ~ d1 + d2
    kill(model2) <- ~ d1 + d2
  }
  if (!isC) {
    kill(model1) <- ~ c1 + c2
    kill(model2) <- ~ c1 + c2
  }
  
  sim.model1 <- model1
##  intercept(sim.model1,outcomes) <- mu
  regfix(sim.model1,to=outcomes[1],from="a1") <- lambda[1]
  regfix(sim.model1,to=outcomes[2],from="a1") <- lambda[1]
  regfix(sim.model1,to=outcomes[1],from="c1") <- lambda[2]
  regfix(sim.model1,to=outcomes[2],from="c1") <- lambda[2]
  ## regfix(sim.model1,to=outcomes[1],from="d1") <- lambda[3]
  ## regfix(sim.model1,to=outcomes[2],from="d2") <- lambda[3]
  regfix(sim.model1,to=outcomes[1],from="e1") <- lambda[4]
  regfix(sim.model1,to=outcomes[2],from="e2") <- lambda[4]

  if (nk1>0) {
    for (i in 1:nk1) {
        regfix(sim.model1, from=covars.1[i], to=outcomes[1],silent=TRUE) <- k1[i]
        regfix(sim.model1, from=covars.2[i], to=outcomes[2],silent=TRUE) <- k1[i]
      }    
  }
  if (nk2>0) {
    for (i in 1:nk2) {
      regfix(sim.model1,from=covars2[i],to=outcomes[1],silent=TRUE) <- k2[i]
      regfix(sim.model1,from=covars2[i],to=outcomes[2],silent=TRUE) <- k2[i]
    }
  }
  intercept(sim.model1, latent(sim.model1)) <- 0
  
  sim.model2 <- model2
##  intercept(sim.model2,outcomes) <- mu
  regfix(sim.model2,to=outcomes[1],from="a1") <- lambda[1]
  regfix(sim.model2,to=outcomes[2],from="a2") <- lambda[1]
  regfix(sim.model2,to=outcomes[1],from="c1") <- lambda[2]
  regfix(sim.model2,to=outcomes[2],from="c1") <- lambda[2]
  ## regfix(sim.model2,to=outcomes[1],from="d1") <- lambda[3]
  ## regfix(sim.model2,to=outcomes[2],from="d2") <- lambda[3]
  regfix(sim.model2,to=outcomes[1],from="e1") <- lambda[4]
  regfix(sim.model2,to=outcomes[2],from="e2") <- lambda[4]
  if (nk1>0) {
    for (i in 1:nk1) {
        regfix(sim.model2, from=covars.1[i], to=outcomes[1],silent=TRUE) <- k1[i]
        regfix(sim.model2, from=covars.2[i], to=outcomes[2],silent=TRUE) <- k1[i]
      }    
  }
  if (nk2>0) {
    for (i in 1:nk2) {
      regfix(sim.model2,from=covars2[i],to=outcomes[1],silent=TRUE) <- k2[i]
      regfix(sim.model2,from=covars2[i],to=outcomes[2],silent=TRUE) <- k2[i]
    }
  }


  if (binary) {
    binary(sim.model1) <- outcomes
    binary(sim.model2) <- outcomes
  }

  d1 <- sim(sim.model1,n=n,...)
  d1$y1 <- d1$y1+mu
  d1$y2 <- d1$y2+mu
  d2 <- sim(sim.model2,n=n,...)
  d2$y1 <- d2$y1+mu
  d2$y2 <- d2$y2+mu
  

  d1$zyg <- "MZ"; d1$twinid <- 1:nrow(d1)
  d2$zyg <- "DZ"; d2$twinid <- 1:nrow(d2)  
  varying <- outcomes
  if(nk1>0)
    varying <- rbind(varying, cbind(covars.1,covars.2))

  dd1 <- reshape(d1[,c(manifest(model1),"zyg","twinid")], varying=varying,
              direction="long",
              timevar="twinnum",
              times=c(1,2),
              v.names=c("y",covars))
  dd2 <- reshape(d2[,c(manifest(model2),"zyg","twinid")], varying=varying,
                 direction="long",
                 timevar="twinnum",
                 times=c(1,2),
                 v.names=c("y",covars))
  long <- rbind(dd1,dd2)
  Wide <- rbind(subset(d1,select=c(manifest(model1),"zyg","twinid")),
                subset(d2,select=c(manifest(model2),"zyg","twinid")))
  
  return(list(data=long,model=list(model1,model2),wide=list(d1,d2),Wide=Wide))
}

