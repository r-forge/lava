
Dbvn <- function(p,design=function(p,...)
                 return(list(mu=cbind(p[1],p[1]),
                             dmu=cbind(1,1),
                             S=matrix(c(p[2],p[3],p[3],p[4]),ncol=2),
                             dS=rbind(c(1,0,0,0),c(0,1,1,0),c(0,0,0,1)))
                        ),
                        Y=cbind(0,0)) {
  mS <- design(p)
  U0 <- with(mS,.Call("biprobit0",
                          mu,
                          S,dS,Y,dmu,NULL,FALSE));
  return(c(U0,mS))
}

## reload <- function(...) {
##   myfun <- "bicif"; pkg <- "bptwin/src/bptwin.so"
##   if (is.loaded(myfun)) dyn.unload(pkg)
##   library(numDeriv)
##   library(mvtnorm)
##   dyn.load(pkg)
##   for (i in dir("bptwin/R/", pattern=glob2rx("*R"), full.names=TRUE)) { source(i) }
## };## reload()

## p <- c(0,c(1,0,1))
## Dbvn(p)
## fp <- function(p) Dbvn(p)$loglik[1]
## sp <- function(p) with(Dbvn(p),score*exp(loglik[1]))
## ##grad(fp,p)
## gp <- function(p) with(Dbvn(p),pmvnorm(upper=c(0,0),mean=mu[1,],sigma=S))[1]
## (gp(p)-gp(p+c(1,0,0,0)*1e-5))/1e-5
## (gp(p)-gp(p+c(0,1,0,0)*1e-5))/1e-5
## (gp(p)-gp(p+c(0,0,1,0)*1e-5))/1e-5
## (gp(p)-gp(p+c(0,0,0,1)*1e-5))/1e-5
## grad(gp,p,method="simple")

###{{{ bpACE/bptwin

bpACE <- bptwin <- function(formula, data, id, zyg, twinnum, DZ, weight=NULL,
                            truncweight=NULL,
                            entry=NULL, time,
                            B=NULL, degree=1, Bconstrain=TRUE,
                            control=list(trace=1),
                            type="ace",
                            eqmean=TRUE,
                            param=0,
                            robustvar=TRUE,                            
                            p, indiv=FALSE, debug=FALSE,...) {

###{{{ setup

  mycontrol <- list(trace=1,
                    method="nlminb",
                    mle=is.null(weight))
  if (length(control)>0)
    mycontrol[names(control)] <- control
  if (length(grep("flex",tolower(type)))>0) { type <- "u"; eqmean <- FALSE }

  Y <- suppressWarnings(model.frame(formula,data,na.action="na.pass")[,1])
  Blen <- 0
  Bord <- c()
  if (degree>0)
    if (is.Surv(Y) | !missing(time) | !is.null(B)) {
      if (is.Surv(Y)) {
        if (ncol(Y)==3) {
          entry <- Y[,1]
          time <- Y[,2]
          Y <- Y[,3]*1 # event
        } else {
          time <- Y[,1]
          Y <- Y[,2]*1
        }
      } else {
        if (is.character(time))
          time <- data[,time]
        if (is.character(entry))
          entry <- data[,entry]
      }
      if (is.null(B)) {
        B <- bs(time,degree=degree)
        Bord <- cbind(time,B)[order(time),,drop=FALSE]
        ##      B <- bs(data[,time],degree=degree)
        ##      Bord <- cbind(data[,time],B)[order(data[,time]),,drop=FALSE]
        colnames(B) <- paste("B","_",1:ncol(B),sep="")
        ##      data <- cbind(data,B)
        ##     formula <- update(formula, paste(".~.+",paste(colnames(B),collapse="+",sep="")))
        Blen <- ncol(B)   
      }
    }
  mycall <- match.call()
  idtab <- table(data[,id])
  ii0 <- order(data[as.character(data[,id])%in%names(idtab)[idtab==2],id])
  data0 <- data[ii0,]
  if (!missing(DZ)) {
    idx1 <- data0[,zyg]==DZ
    idx0 <- data0[,zyg]!=DZ
  } else {
    idx1 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[1] # DZ
    idx0 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[2] # MZ
  }

  N <- cbind(sum(idx0),sum(idx1)); colnames(N) <- c("MZ","DZ");rownames(N) <- ""

  B <- B[ii0,]
  ff <- as.formula(paste("~",paste(attributes(terms(formula))$term.labels,
                                   collapse="+"),"+1",sep=""))
  X <- cbind(model.matrix(ff,data0),B)
  ##  X <- cbind(model.matrix(formula,data0))
  nx <- ncol(X)
  yvar <- paste(deparse(formula[[2]]),collapse="")
  
  Y <- cbind(as.numeric(Y[ii0])-!is.numeric(Y[ii0]))
  ##cbind(as.numeric(data0[,yvar]))-(!is.numeric(data0[,yvar]))
  W0 <- W1 <- NULL
  if (!is.null(weight)) {
    W <- cbind(data0[,weight])
    W0 <- matrix(W[idx0,,drop=FALSE],ncol=2,byrow=TRUE)
    W1 <- matrix(W[idx1,,drop=FALSE],ncol=2,byrow=TRUE)
  }
  X0 <- X[idx0,,drop=FALSE]
  X1 <- X[idx1,,drop=FALSE]
  Y0 <- matrix(Y[idx0,,drop=FALSE],ncol=2,byrow=TRUE)
  Y1 <- matrix(Y[idx1,,drop=FALSE],ncol=2,byrow=TRUE)
  XX0 <- matrix(t(X0),ncol=nx*2,byrow=TRUE)
  XX1 <- matrix(t(X1),ncol=nx*2,byrow=TRUE)

  
  midx <- 1:nx
  dS0 <- rbind(rep(1,4),rep(1,4),rep(1,4)) ## MZ
  dS1 <- rbind(c(1,.5,.5,1),rep(1,4),c(1,.25,.25,1)) ## DZ
  ACDU <- sapply(c("a","c","d","e","u"),function(x) length(grep(x,tolower(type)))>0)
  if (ACDU["u"]) {
    dS0 <- rbind(rep(1,4))
    vidx0 <- 1
    vidx1 <- 2
    dS1 <- dS0
    nvar <- length(vidx0)+length(vidx1)
  } else {
    nvar <- sum(ACDU[1:3])
    vidx0 <- vidx1 <- 1:nvar
    dS0 <- dS0[ACDU[1:3],,drop=FALSE]
    dS1 <- dS1[ACDU[1:3],,drop=FALSE]
  }  
  if (eqmean) {
    midx0 <- midx1 <- midx    
  } else {
    midx0 <- 1:nx; midx1 <- midx0+nx
    nx <- 2*nx;
  }
  
  vidx0 <- vidx0+nx; vidx1 <- vidx1+nx
  vidx <- nx+seq_len(nvar)
  midx <- seq_len(nx)
  plen <- nx+nvar

  Am <- matrix(c(1,.5,.5,1),ncol=2)
  Dm <- matrix(c(1,.25,.25,1),ncol=2)
  Em <- diag(2)
  
  ##mytr <- function(x) x; dmytr <- function(x) 1
  ##  mytr <- function(x) x^2; dmytr <- function(x) 2*x
  ##mytr <- function(z) 1/(1+exp(-z)); dmytr <- function(z) exp(-z)/(1+exp(-z))^2
  mytr <- exp; dmytr <- exp
  Sigma <- function(p0) {    
    p0[vidx] <- mytr(p0[vidx])    
    if (ACDU["u"]) {     
      Sigma0 <- Em+p0[plen-1]; Sigma1 <- Em+p0[plen]
    } else {      
      pv <- ACDU*1;  pv[which(ACDU[1:3])] <- p0[vidx]
      Sigma0 <- Em*pv["e"] + pv["a"] + pv["c"] + pv["d"]
      Sigma1 <- Em*pv["e"] + pv["a"]*Am + pv["c"] + pv["d"]*Dm
    }
    return(list(Sigma0=Sigma0,Sigma1=Sigma1))
  }  

###}}} setupx
  
###{{{ U  

  U <- function(p,indiv=FALSE) {
    b0 <- cbind(p[midx0])
    b1 <- cbind(p[midx1])
    b00 <- b0; b11 <- b1
    if (Bconstrain) {
      b00 <- trMean(b0,Blen); b11 <- trMean(b1,Blen)
    }
    S <- Sigma(p)    
    mu0 <- X0%*%b00
    Mu0 <- matrix(mu0,ncol=2,byrow=TRUE)
    U0 <- .Call("biprobit0",
                Mu0,
                S$Sigma0,dS0,Y0,XX0,W0,!is.null(W0))         
    mu1 <- X1%*%b11
    Mu1 <- matrix(mu1,ncol=2,byrow=TRUE)      
    U1 <- .Call("biprobit0",
                Mu1,
                S$Sigma1,dS1,Y1,XX1,W1,!is.null(W1))    
    if (indiv) {
      val <- matrix(0,ncol=plen,nrow=nrow(U0$score)+nrow(U1$score))
      val[seq_len(nrow(U0$score)),c(midx0,vidx0)] <- U0$score
      val[nrow(U0$score)+seq_len(nrow(U1$score)),c(midx1,vidx1)] <- U1$score
      for (ii in vidx) {
        val[,ii] <- val[,ii]*dmytr(p[ii])
      }
      if (Bconstrain & Blen>0) {
        Bidx <- attributes(b00)$idx
        val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(b00)$D[Bidx,Bidx,drop=FALSE])
        if (!eqmean) {
          Bidx <- midx1[attributes(b11)$idx]
          val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(b11)$D[attributes(b11)$idx,attributes(b11)$idx])
        }
      }      
      attributes(val)$logLik <- c(U0$loglik,U1$loglik)
      return(val)
    }    
    val <- numeric(plen)
    val[c(midx0,vidx0)] <- colSums(U0$score)
    val[c(midx1,vidx1)] <- val[c(midx1,vidx1)]+colSums(U1$score)
    for (ii in vidx)
      val[ii] <- val[ii]*dmytr(p[ii])
    if (Bconstrain & Blen>0) {
      Bidx <- attributes(b00)$idx
      val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b00)$D[Bidx,Bidx])
      if (!eqmean) {
        Bidx <- midx1[attributes(b11)$idx]
        val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b11)$D[attributes(b11)$idx,attributes(b11)$idx])
      }
    }    
    attributes(val)$logLik <- sum(U0$loglik)+sum(U1$loglik)
    return(val)
  }

###}}} U

###{{{ Concordance model - Left Truncation/Delayed Entry

  U2 <- function(p,t,indiv=FALSE) {
    b0 <- cbind(p[midx0])
    b1 <- cbind(p[midx1])
    B0 <- trMean(b0,Blen); B1 <- trMean(b1,Blen)    
    S <- Sigma(p)
    X00 <- X0; X11 <- X1
    Y00 <- Y0; Y11 <- Y1
    XX00 <- XX0; XX11 <- XX1;
    if (!missing(t)) {
      timeidx <- which(t==Bord[,1])
      t0 <- data0[idx0,time]<=t
      t0. <- matrix(t0,ncol=2,byrow=TRUE)
      t1 <- data0[idx1,time]<=t
      t1. <- matrix(t1,ncol=2,byrow=TRUE)
      X00 <- X0[t0,,drop=FALSE];
      Y00 <- Y0[t0.[,1],,drop=FALSE];
      X11 <- X1[t1,,drop=FALSE];
      Y11 <- Y1[t1.[,1],,drop=FALSE];
      XX00 <- XX0[t0.[,1],,drop=FALSE]
      XX11 <- XX1[t1.[,1],,drop=FALSE]
    }
    
    mu0 <- X00%*%B0
    Mu0 <- matrix(mu0,ncol=2,byrow=TRUE)
    U0 <- .Call("biprobit0",
                Mu0,
                S$Sigma0,dS0,Y00,XX00,W0,FALSE)         
    mu1 <- X11%*%B1
    Mu1 <- matrix(mu1,ncol=2,byrow=TRUE)      
    U1 <- .Call("biprobit0",
                Mu1,
                S$Sigma1,dS1,Y11,XX11,W1,FALSE)
    
    if (indiv) {
      val <- matrix(0,ncol=plen,nrow=nrow(U0$score)+nrow(U1$score))
      val[seq_len(nrow(U0$score)),c(midx0,vidx0)] <- U0$score
      val[nrow(U0$score)+seq_len(nrow(U1$score)),c(midx1,vidx1)] <- U1$score
      for (ii in vidx) {
        val[,ii] <- val[,ii]*dmytr(p[ii])
      }
      if (Bconstrain & Blen>0) {
        Bidx <- attributes(B0)$idx
        val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(B0)$D[Bidx,Bidx,drop=FALSE])
        if (!eqmean) {
          Bidx <- midx1[attributes(B1)$idx]
          val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(B1)$D[attributes(B1)$idx,attributes(B1)$idx])
        }
      }      
      attributes(val)$logLik <- c(U0$loglik,U1$loglik)
      return(val)
    }
    
    val <- numeric(plen)
    val[c(midx0,vidx0)] <- colSums(U0$score)
    val[c(midx1,vidx1)] <- val[c(midx1,vidx1)]+colSums(U1$score)
    for (ii in vidx)
      val[ii] <- val[ii]*dmytr(p[ii])

    ##    browser()
    if (Bconstrain & Blen>0) {
      Bidx <- attributes(b00)$idx
      val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b00)$D[Bidx,Bidx])
      if (!eqmean) {
        Bidx <- midx1[attributes(b11)$idx]
        val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b11)$D[attributes(b11)$idx,attributes(b11)$idx])
      }
    }
    
    attributes(val)$logLik <- sum(U0$loglik)+sum(U1$loglik)
    return(val)
  }

###}}} Concordance model - Left Truncation/Delayed Entry
 
 
  p0 <- rep(0,plen); p0[vidx] <- 1
  if (!is.null(control$start)) {
    p0 <- control$start
    control$start <- NULL
  } else {
    g <- suppressWarnings(glm(Y~-1+X,family=binomial(probit)))
    p0[midx] <- coef(g)
    if (Blen>0) {
      pB <- p0[tail(midx,Blen)]
      pB[1] <- ifelse(pB[1]<0,-2,log(pB[1]))
      if (Blen>1) {
        pB[seq_len(Blen-1)+1] <- -2
      }
      p0[tail(midx,Blen)] <- pB
    }
  }

  if (!missing(p)) return(U(p,indiv=indiv))
  
  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))


  nlminbopt <- intersect(names(mycontrol),c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min"))
  ucminfopt <- intersect(names(mycontrol),c("trace","grtol","xtol","stepmax","maxeval","grad","gradstep","invhessian"))
  optimopt <- names(mycontrol)

  if (debug) browser()

  if (mycontrol$mle) {
    op <- switch(mycontrol$method,
                 ucminf=ucminf(p0,fn=f0,gr=g0,control=mycontrol[ucminfopt],hessian=F,...),
                 optim=optim(p0,fn=f0,gr=g0,control=mycontrol[ucminfopt],...),
                 nlminb(p0,f0,grad=g0,control=mycontrol[nlminbopt],...))
  } else {
    op <- switch(mycontrol$method,
                 ucminf=ucminf(p0,f,control=mycontrol[ucminfopt],hessian=F,...),
                 optim=optim(p0,f,control=mycontrol[ucminfopt],...),
                 nlminb(p0,f,control=mycontrol[nlminbopt],...))
    ##  op <- nlm(f,p0,print.level=2)
    ##  op <- spg(p0,f,control=control,...)
  }
  UU <- U(op$par,indiv=TRUE)
  I <- -numDeriv::jacobian(U,op$par)
  tol <- 1e-15
  V <- Inverse(I,tol)
  sqrteig <- attributes(V)$sqrteig
  J <- NULL
  if (robustvar) {
    J <- crossprod(UU)
    V <- V%*%J%*%V
  }
  if (any(sqrteig<tol)) warning("Near-singular covariance matrix (pseudo-inverse used)")
  
  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  rnames1 <- colnames(X0)
  vnames1 <- NULL
  trnam <- " "
  if (debug) browser()
  if (!eqmean) {
    rnames1 <- c(paste(rnames1,"MZ",sep=trnam),paste(rnames1,"DZ",sep=trnam))
  }
  if (ACDU["u"]) {
    rnames <- c(rnames1,paste(c("log(var(U))","log(var(U))"),c("MZ","DZ"),sep=trnam))
  } else {
    rnames <- c(rnames1,c("log(var(A))","log(var(C))","log(var(D))")[ACDU[1:3]])
  }
  rownames(cc) <- rnames
  S <- Sigma(op$par)
  val <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, Sigma0=S$Sigma0, Sigma1=S$Sigma1, dS0=dS0, dS1=dS1, call=mycall, N=N, data=data0, Blen=Blen, midx0=midx0, midx1=midx1, vidx0=vidx0, vidx1=vidx1, eqmean=eqmean, B=Bord,I=I,J=J)
  class(val) <- "bptwin"
  return(val)
}

###}}} bpACE

###{{{ bptwin methods

plot.bptwin <- function(x,n=50,rg=range(x$B[,1]),xlab="Time",ylab="Concordance",...) {
  require(mvtnorm)
  if (x$Blen>0) {
    ##    rg <- range(x$B[,1])
    t <- seq(rg[1],rg[2],length.out=n)
    B0 <- bs(t,degree=x$Blen)
    b0. <- coef(x)[x$midx0,1,drop=FALSE]
    b1. <- coef(x)[x$midx1,1,drop=FALSE]
    b0 <- trMean(b0.,x$Blen)
    b1 <- trMean(b1.,x$Blen)
##    t <- x$B[,1]
##    B0 <- x$B[,-1,drop=FALSE]    
    b00 <- tail(b0,x$Blen)
    b11 <- tail(b1,x$Blen)
##    b00 <- b0.[-1]
##    b11 <- b1.[-1]
    pr0 <- sapply(as.numeric(B0%*%b00+b0[1]), function(z)
                  pmvnorm(upper=rep(z,2),sigma=x$Sigma0))
    pr1 <- sapply(as.numeric(B0%*%b11+b1[1]), function(z)
                  pmvnorm(upper=rep(z,2),sigma=x$Sigma1))
    plot(pr0~t,type="l", xlab=xlab, ylab=ylab,...)
    lines(pr1~t,type="l",lty=2)
  }
  invisible(x)
}

sim <- function(x,...) UseMethod("sim")
sim.bptwin <- function(x,n=100,p,...) {
  browser()
  return(x)
}
score <- function(x,...) UseMethod("score")
score.bptwin <- function(x,indiv=FALSE,...) {
  if (indiv) { s <- x$score; attributes(s)$logLik <- NULL; return(s) }
  colSums(x$score)
}
logLik.bptwin <- function(object,indiv=FALSE,...) {
  if (indiv) return(object$logLik)
  n <- sum(object$N)/2
  p <- nrow(coef(object))
  loglik <- sum(object$logLik)
  attr(loglik, "nall") <- n
  attr(loglik, "nobs") <- n-p
  attr(loglik, "df") <- p
  class(loglik) <- "logLik"        
  return(loglik)
}
vcov.bptwin <- function(object,...) object$vcov
coef.bptwin <- function(object,...) object$coef
print.bptwin <- function(x,...) {
  printCoefmat(coef(x))
  S <- colSums(x$score);  names(S) <- rep("",length(S))
  cat("\n")
  print(x$N)
  cat("Score: "); cat(S);
  cat("\nlogLik: "); cat(sum(x$logLik),"\n");
}

summary.bptwin <- function(object,level=0.05,...) {
  trnam <- " "
  vcoef1 <- paste("log(var(",c("A","C","D"),"))",sep="")
  vcoef2 <- paste("log(var(",
                  c(paste("U))","MZ",sep=trnam),
                    paste("U))","DZ",sep=trnam)),sep="")
  idx1 <- na.omit(match(vcoef1,rownames(coef(object))))
  idx2 <- na.omit(match(vcoef2,rownames(coef(object))))
  if (length(idx2)>0) {
    idx <- idx2
    mz <- multinomlogit(coef(object)[idx2[1]]); names(mz) <- c("U","E")
    dz <- multinomlogit(coef(object)[idx2[2]]); names(dz) <- c("U","E")
    cc <- c(mz[1],dz[1]) ##,mz[2],dz[2])
    names(cc) <- c("Correlation MZ","Correlation DZ")
    corMZ <- mz[1]; corDZ <- dz[1]
    D <- (cbind(c(attributes(mz)$gradient[1],0),c(0,attributes(dz)$gradient[1])))
    h <- function(x) 2*(x[1]-x[2])
    dh <- function(x) c(2,-2)
    i1 <- 1:2
    corr <- NULL
  }
  if (length(idx1)>0) {
    idx <- idx1
    ACD <- match(rownames(coef(object))[idx1],vcoef1)
    nn <- c(c("A","C","D")[ACD],"E")
    dzsc <- c(1/2,1,1/4)[ACD]
    cc <- multinomlogit(coef(object)[idx1,1]); names(cc) <- nn
    D <- attributes(cc)$gradient
    K <- length(ACD)
    Ki <- seq_len(K)
    corMZ <- sum(cc[Ki]); corDZ <- sum(cc[Ki]*dzsc)
    i1 <- seq_len(length(dzsc))
    h <- function(x) 2*(sum(x[i1])-sum(x[i1]*dzsc))
    dh <- function(x) 2*(1-dzsc)    
  }
  V <- vcov(object)[idx,idx]
  Vc <- D%*%V%*%t(D)
  if (length(idx1)>0) {
    b <- cbind(rep(1,K))
    corMZ.sd <- (t(b)%*%Vc[Ki,Ki]%*%b)[1]^0.5
    corDZ.sd <- (t(dzsc)%*%Vc[Ki,Ki]%*%dzsc)[1]^0.5    
    corr <- rbind(c(corMZ,corMZ.sd),c(corDZ,corDZ.sd))
    rownames(corr) <- c("Correlation MZ","Correlation DZ")
  }
  newcoef <- rbind(cbind(cc,diag(Vc)^0.5),corr); colnames(newcoef) <- c("Estimate","Std.Err")  
  logit <- function(p) log(p/(1-p))
  tigol <- function(z) 1/(1+exp(-z))
  dlogit <- function(p) 1/(p*(1-p))
  logith <- function(x) logit(h(x))
  dlogith <- function(x) dlogit(h(x))*dh(x)
  Dlh <- dlogith(cc[i1])
  sdlh <- (t(Dlh)%*%Vc[i1,i1]%*%(Dlh))[1]^0.5
  H <- h(cc[i1])
  hstd <- t(dh(cc[i1]))%*%Vc[i1,i1]%*%dh(cc[i1])
  ci <- tigol(logith(cc[i1]) + qnorm(1-level/2)*c(-1,1)*sdlh)  
  
  concordance <-  conditional <- marg <- c()
  mu <- coef(object)[c(object$midx0[1],object$midx1[1]),1,drop=TRUE]
  Sigma <- list(object$Sigma0,object$Sigma1)
  for (i in 1:2) {    
    mu.cond <- function(x) mu+Sigma[[i]][1,2]/Sigma[[i]][2,2]*(x-mu[i])
    var.cond <- Sigma[[i]][1,1]-Sigma[[i]][1,2]^2/Sigma[[i]][2,2]    
    cc0 <- pmvnorm(upper=c(mu[i],mu[i]),sigma=Sigma[[i]])
    px <- pnorm(mu[i],sd=Sigma[[i]][2,2]^0.5)
    concordance <- c(concordance,cc0)
    marg <- c(marg,px)
    conditional <- c(conditional,cc0/px)
  }
  names(concordance) <- names(conditional) <- c("MZ","DZ")
  hval <- rbind(c(H,hstd^0.5,ci)); colnames(hval) <- c("Estimate","Std.Err",paste(100*c(level/2,1-level/2),"%",sep="")); rownames(hval) <- "Heritability"  
  res <- list(object=object, h=hval,
              coef=newcoef, concordance=concordance, conditional=conditional)
  class(res) <- "summary.bptwin"
  res
}

print.summary.bptwin <- function(x,...) {
  cat("\n")
  print(x$object)
  cat("\n")
  print(x$coef)
  cat("\nConcordance (MZ; DZ):\t\t", x$concordance,"\n")
  cat("Case-wise concordance (MZ; DZ):\t", x$conditional,"\n\n")
  print(x$h)
  cat("\n")
}

###}}} bptwin methods

###{{{ bptwin0

bptwin0 <- function(formula, data, id, zyg, twinnum, weight=NULL,
                    control=list(trace=1), intercepts=ifelse(unique,2,1),
                    epsilon=ifelse(unique,1,2), unique=TRUE, p,...) {

  mycall <- match.call()

  idtab <- table(data[,id])
  data0 <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  data0 <- data0[order(data0[,id]),]
  idx1 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[1]
  idx2 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[2]

  X <- model.matrix(formula,data0)
  nx <- ncol(X)
  yvar <- paste(deparse(formula[[2]]),collapse="")
  Y <- cbind(as.numeric(data0[,yvar]))-(!is.numeric(data0[,yvar]))
  W0 <- W1 <- NULL
  if (!is.null(weight)) {
    W <- cbind(data0[,weight])
    W0 <- matrix(W[idx1,,drop=FALSE],ncol=2,byrow=TRUE)
    W1 <- matrix(W[idx2,,drop=FALSE],ncol=2,byrow=TRUE)
  }
  X0 <- X[idx1,,drop=FALSE]
  X1 <- X[idx2,,drop=FALSE]
  Y0 <- matrix(Y[idx1,,drop=FALSE],ncol=2,byrow=TRUE)
  Y1 <- matrix(Y[idx2,,drop=FALSE],ncol=2,byrow=TRUE)
  XX0 <- matrix(t(X0),ncol=nx*2,byrow=TRUE)
  XX1 <- matrix(t(X1),ncol=nx*2,byrow=TRUE)
  midx <- 1:nx
  plen <- nx+1
  dS0 <- rbind(rep(1,4))
  dS1 <- dS0##rbind(dS0,c(1,0,0,1))
  
  mytr <- function(x) (x)  
  U <- function(p,group=1,indiv=FALSE) {
    if (group==1) {
      B <- cbind(p[midx])
      Sigma <- diag(2) + mytr(p[plen])
      mu <- X0%*%B
      Mu <- matrix(mu,ncol=2,byrow=TRUE)
      U1 <- .Call("biprobit0",
                  Mu,
                  Sigma,dS0,Y0,XX0,W0,!is.null(W0))     
      
    } else {
      B <- cbind(p[midx])
      Sigma <- diag(2) + mytr(p[plen])      
      mu <- X1%*%B
      Mu <- matrix(mu,ncol=2,byrow=TRUE)      
      U1 <- .Call("biprobit0",
                  Mu,
                  Sigma,dS1,Y1,XX1,W1,!is.null(W1))
    }
    if (indiv) {
      val <- U1$score
      attributes(val)$logLik <- U1$loglik
      return(val)
    }
    val <- colSums(U1$score)
    attributes(val)$logLik <- sum(U1$loglik)
    return(val)
  }

  mystart <- list(rep(0,plen),rep(0,plen))
  if (!is.null(control$start)) {
    mystart <- control$start
    control$start <- NULL
  }

  if (!missing(p)) return(U(p,indiv=FALSE))

  res <- list()
  for (i in 1:2) {
    p0 <- mystart[[i]]
    f <- function(p) crossprod(U(p,group=i))[1]
    f0 <- function(p) -sum(attributes(U(p,group=i))$logLik)
    if (!is.null(control$simple)) {
      control$simple <- NULL
      op <- nlminb(p0,f0,control=control,...)
    } else {
      op <- nlminb(p0,f,control=control,...)
    }

    
    UU <- U(op$par,group=i,indiv=TRUE)
    J <- crossprod(UU)
    U. <- function(p) U(p,group=i)
    I <- -numDeriv::jacobian(U.,op$par)
    iI <- solve(I)
    ##    I <- (I+t(I))/2
    V <- iI%*%J%*%iI
    ##    V <- iJ
    ##V <- solve(I)
    
    cc <- cbind(op$par,sqrt(diag(V)))
    cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
    colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
    rnames1 <- colnames(X0)
    paste(rnames1,i,sep=".")
    vnames1 <- NULL
    trnam <- ""
    rownames(cc) <- c(rnames1,paste(trnam,"var(U",i,")",sep=""))
    val <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall,tr=mytr)
    class(val) <- "bptwin"
    res <- c(res, list(val))
  }
  class(res) <- "bptwin.list"
  return(res)
}

###}}} bptwin0

###{{{ bptwin1

bptwin1 <- function(formula, data, id, zyg, twinnum, weight=NULL,
                    control=list(trace=1), intercepts=ifelse(unique,2,1),
                    epsilon=ifelse(unique,1,2), unique=TRUE, p,
                    
                    ...) {

  mycall <- match.call()
  idtab <- table(data[,id])
  data0 <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  data0 <- data0[order(data0[,id]),]
  idx1 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[1]
  idx2 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[2]

  X <- model.matrix(formula,data0)
  nx <- ncol(X)
  yvar <- paste(deparse(formula[[2]]),collapse="")
  Y <- cbind(as.numeric(data0[,yvar]))-(!is.numeric(data0[,yvar]))
  W0 <- W1 <- NULL
  if (!is.null(weight)) {
    W <- cbind(data0[,weight])
    W0 <- matrix(W[idx1,],ncol=2,byrow=TRUE)
    W1 <- matrix(W[idx2,],ncol=2,byrow=TRUE)
  }
  X0 <- X[idx1,,drop=FALSE]
  X1 <- X[idx2,,drop=FALSE]
  Y0 <- matrix(Y[idx1,,drop=FALSE],ncol=2,byrow=TRUE)
  Y1 <- matrix(Y[idx2,,drop=FALSE],ncol=2,byrow=TRUE)
  XX0 <- matrix(t(X0),ncol=nx*2,byrow=TRUE)
  XX1 <- matrix(t(X1),ncol=nx*2,byrow=TRUE)
  dS0 <- cbind(1,4)
  dS1 <- cbind(dS0,c(1,0,0,1))

  mytr <- function(x) (x)
  if (unique) {
    plen <- 2*nx+3
    midx1 <- seq_len(nx)
    midx2 <- midx1+nx
    if (intercepts==1) midx2 <- midx1
  } else {  
    plen <- nx+1+3
    midx1 <- c(1,2+seq_len(nx-1))
    midx2 <- c(2,2+seq_len(nx-1))
    if (intercepts==1) {
      plen <- nx+3
      midx2 <- midx1 <- seq_len(nx)    
    }
  }
  vidx1 <- c(plen-2)
  vidx2 <- c(plen-1,plen)
  if (epsilon==1) {
    plen <- plen-1
    vidx1 <- plen-1
    vidx2 <- plen
  }

  mytr <- function(x) x
  U <- function(p,indiv=FALSE) {
    B0 <- cbind(p[midx1])
    B1 <- cbind(p[midx2])
    if (epsilon==2) {
      Sigma0 <- diag(2)+(mytr(p[plen-2]))
      Sigma1 <- diag(rep(mytr(p[plen]),2)) + mytr(p[plen-1])
    } else {
      Sigma0 <- diag(2) + mytr(p[plen-1])
      Sigma1 <- diag(2) + mytr(p[plen])
    }

    ##    Sigma1 <- diag(2) +mytr(p[plen-1])
    ##    Sigma1 <- diag(2)+mytr(p[plen-1])
    ##    B0 <- cbind(p[1:2])
    ##    Sigma0 <- diag(2) + mytr(p[3])
    ##    B1 <- cbind(p[4:5])
    ##    Sigma1 <- diag(2) + mytr(p[6])
    ##    Sigma1 <- diag(rep(mytr(p[5]),2)) +mytr(p[4])
    
    mu0 <- X0%*%B0
    Mu0 <- matrix(mu0,ncol=2,byrow=TRUE)
    mu1 <- X1%*%B1
    Mu1 <- matrix(mu1,ncol=2,byrow=TRUE)
    U0 <- .Call("biprobit0",
                Mu0,
                Sigma0,dS0,Y0,XX0,W0,!is.null(W0),FALSE)
    U1 <- .Call("biprobit0",
                Mu1,
                Sigma1,dS1,Y1,XX1,W1,!is.null(W1),epsilon==2)

    ##    browser()
    if (indiv) {
      l <- c(U0$loglik,U1$loglik)
      val1 <- matrix(0,ncol=plen,nrow=(nrow(U0$score)))
      val1[,c(midx1,vidx1)] <- U0$score
      val2 <- matrix(0,ncol=length(p),nrow=(nrow(U1$score)))
      val2[,c(midx2,vidx2)] <- U1$score
      val <- rbind(val1,val2)
      attributes(val)$logLik <- l
      return(val)
    }
    val <- numeric(length(p))
    l <- 0
    val[c(midx1,vidx1)] <- colSums(U0$score)
    l <- l+sum(U0$loglik)
    val[c(midx2,vidx2)] <- val[c(midx2,vidx2)] + colSums(U1$score)
    ##    val[1:3] <- val[1:3] + colSums(U1$score)  
                                        #    val <- rep(0,length(p));    
                                        #    val[1:3] <- colSums(U0$score)
    ##    val[c(3,5,6)] <- val[c(3,5,6)] + colSums(U1$score)
                                        #    val[c(4:6)] <- val[4:6] + colSums(U1$score)
    l <- l+sum(U1$loglik)
    attributes(val)$logLik <- l
    return(val)
  }
  
  p0 <- rep(0,plen)  
  if (!is.null(control$start)) {    
    p0 <- control$start
    control$start <- NULL
  } 
  if (!missing(p)) return(U(p,indiv=FALSE))
  
  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  ##  browser()
  
  if (!is.null(control$simple)) {
    control$simple <- NULL
    op <- nlminb(p0,f0,control=control,...)
  } else {
    op <- nlminb(p0,f,control=control,...)
  }
  ##  UU <- U(op$par,indiv=FALSE)
  ##  res <- list(score=UU,logLik=sum(attributes(UU)$logLik),opt=op)
  ##  return(res)

  UU <- U(op$par,indiv=TRUE)
  J <- crossprod(UU)
  iJ <- solve(J)  
  I <- numDeriv::jacobian(U,op$par)
  Is <- (I+t(I))/2
  ## V <- J%*%Is%*%J
  V <- iJ
  
  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  rnames1 <- colnames(X0)
  if (unique) {
    rnames1 <- c(paste(rnames1,"1",sep="."),paste(rnames1,"2",sep="."))
  } else {
    rnames1 <- c(rep(rnames1[1],intercepts-1),rnames1)
    rnames1[seq_len(intercepts)] <- paste(rnames1[seq_len(intercepts)],seq_len(intercepts),sep=".")
  }
  vnames1 <- NULL
  trnam <- ""
  if (epsilon==2) vnames1 <- paste(trnam,"var(E2)",sep="")
  rownames(cc) <- c(rnames1,paste(trnam,"var(U1)",sep=""),paste(trnam,"var(U2)",sep=""),vnames1)
  res <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall,tr=mytr)
  class(res) <- "bptwin"
  return(res)
}

###}}} bptwin

###{{{ bptwin.list methods

summary.bptwin.list <- function(object,...) {
  concordance <- c()
  conditional <- c()
  for (x in object) {
    mu <- rep(coef(x)[1],1)
    Sigma <- diag(2) + x$tr(coef(x)[2])
    mu.cond <- function(x) mu+Sigma[1,2]/Sigma[2,2]*(x-mu)
    var.cond <- Sigma[1,1]-Sigma[1,2]^2/Sigma[2,2]    
    cp <- pnorm(mu[1])
    cc <- pmvnorm(lower=c(0,0),mean=c(mu,mu),sigma=Sigma)
    px <- 1-pnorm(0,mu,Sigma[2,2]^0.5)
    concordance <- c(concordance,cc)
    conditional <- c(conditional,px)
    print(x)
    cat("Concordance:", cc,"\n")
    cat("Conditional:", cc/px,"\n\n")
  }
  invisible(list(object,concordance=concordance,conditional=conditional))
}

print.bptwin.list <- function(x,...) summary(x,...)

###}}} bptwin.list methods

###{{{ biprobit

biprobit <- function(formula, data, id, weight=NULL,
                     control=list(trace=1), unique=TRUE, p,...) {

  mycall <- match.call()
  idtab <- table(data[,id])
  data0 <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  data0 <- data0[order(data0[,id]),]

  X0 <- model.matrix(formula,data0)
  nx <- ncol(X0)
  yvar <- paste(deparse(formula[[2]]),collapse="")
  Y <- cbind(as.numeric(data0[,yvar]))-(!is.numeric(data0[,yvar]))
  W0 <- NULL
  if (!is.null(weight)) {
    W <- cbind(data0[,weight])
    W0 <- matrix(W,ncol=2,byrow=TRUE)
  }
  Y0 <- matrix(Y,ncol=2,byrow=TRUE)
  XX0 <- matrix(t(X0),ncol=nx*2,byrow=TRUE)
  midx <- 1:nx
  plen <- nx+1
  dS0 <- rbind(rep(1,4))
  
  U <- function(p,indiv=FALSE) {
    B <- cbind(p[midx])
    Sigma <- diag(2) + p[plen]
    mu <- X0%*%B
    Mu <- matrix(mu,ncol=2,byrow=TRUE)
    U <- .Call("biprobit0",
               Mu,
               Sigma,dS0,Y0,XX0,W0,!is.null(W0))     
    if (indiv) {
      val <- U$score
      attributes(val)$logLik <- U$loglik
      return(val)
    }
    val <- colSums(U$score)
    attributes(val)$logLik <- sum(U$loglik)
    return(val)
  }

  p0 <- rep(0,plen)
  if (!is.null(control$start)) {
    p0 <- control$start
    control$start <- NULL
  }

  if (!missing(p)) return(U(p,indiv=FALSE))

  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  if (!is.null(control$simple)) {
    control$simple <- NULL
    op <- nlminb(p0,f0,control=control,...)
  } else {
    op <- nlminb(p0,f,control=control,...)
  }
  UU <- U(op$par,indiv=TRUE)
  J <- crossprod(UU)
  iJ <- solve(J)
  I <- -numDeriv::jacobian(U,op$par)    
  V <- iJ%*%I%*%iJ  
  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  rnames1 <- colnames(X0)
  rownames(cc) <- c(rnames1,"var(U)")
  val <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall)
  class(val) <- "bptwin"
  return(val)
}

###}}} biprobit

###{{{ utilities

trMean <- function(b,blen) {
##  mytr <- function(x) x^2; dmytr <- function(x) 2*x
  mytr <- dmytr <- exp
  if (blen==0) return(b) 
  k <- length(b)
  Bidx <- seq_len(blen)+(k-blen)
  b[Bidx[1]] <- mytr(b[Bidx[1]])
  D <- diag(nrow=k)
  D[Bidx[1]:k,Bidx[1]] <- b[Bidx[1]]
  for (i in Bidx[-1]) {
    D[i:k,i] <- dmytr(b[i])
    b[i] <- b[i-1]+mytr(b[i])
  }
  attributes(b)$D <- D
  attributes(b)$idx <- Bidx
  return(b)
}

multinomlogit <- function(x) {
  n <- length(x)
  ex <- exp(x)
  sx <- sum(ex)+1
  f <- c(ex,1)
  df <- c(ex,0)
  res <- f/sx
  dg <- -ex/sx^2   
  gradient <- matrix(ncol=n,nrow=n+1)
  I <- diag(n+1)
  for (i in 1:(n)) {
    gradient[,i] <- df[i]*I[i,]/sx+dg[i]*f
  }
  attributes(res)$gradient <- gradient
  return(res)
}

Inverse <- function(X,tol=1e-9) {
  n <- nrow(X)
  if (nrow(X)==1) {
    res <- 1/X
    if (det) attributes(res)$det <- X
    return(res)
  }
  svdX <- svd(X)
  id0 <- numeric(n)
  id0[svdX$d>tol] <- 1/svdX$d[svdX$d>tol]
  res <- with(svdX, v%*%diag(id0)%*%t(u))
  attributes(res)$sqrteig <- svdX$d
  return(res)
}

###}}} utilities
