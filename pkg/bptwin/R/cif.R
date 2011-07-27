###{{{ cif

cif <- function(t,theta0,knots,degree=1,theta,control=list(trace=1),
                ppar="probpar",vpar="varpar1",...) {
  PP <- prepcif(t,d=degree,knots=knots)
  loglik <- function(theta,d=degree,knots=knots,...) {
    with(PP, logLikcif(theta,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=TRUE,ppar=ppar,vpar=vpar,...))
  }
  if (PP$dim==1) {
    loglik <- function(theta,d=degree,knots=knots,...) {
      with(PP, logLikcif1(theta,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=TRUE,ppar=ppar,vpar=vpar,...))
    }
    if (missing(theta0)) {
      npar <- with(PP, ncol(B)*ncauses+ncauses-1)
      theta0 <- rep(0,npar)
    }
  } else {
    if (missing(theta0)) {
      npar <- with(PP, ncol(B)*ncauses+
        do.call(ppar,list(theta=NULL,ncauses=ncauses))+
          do.call(vpar,list(theta=NULL,ncauses=ncauses)))
      theta0 <- rep(0,npar)
    }
  }

  if (!missing(theta)) return(loglik(theta,...))
  opt <- nlminb(theta0,loglik,control=control)
  res <- c(list(coef=opt$par, opt=opt, ppar=ppar, vpar=vpar),PP)
  class(res) <- "cif"
  ll <- function(p,...) logLik(res,theta=p,...)
  ##  H <- -hessian(ll,coef(res))
  ##  res$information <- H
  l0 <- ll(coef(res),info=TRUE)
  res <- c(res, list(Sigma=l0$Sigma, b=l0$b, P=l0$P, P1=rowSums(as.matrix(l0$P))))
  class(res) <- "cif"
  return(res)

}
prepcif <- function(t,knots,d=1,...) {
  idx <- nn <- flip <- c()
  dim <- 2
  if (ncol(t)==4) {
    causes <- c(0,sort(setdiff(unique(t[,3:4,drop=TRUE]),0)))
    ncauses <- length(causes)-1
    for (i in seq_len(ncauses+1)) {
      for (j in seq(i+1,ncauses+1)) {
        ii <- which(t[,4]==causes[i] & t[,3]==causes[j])
        if (length(ii)>0)  {
          flip <- c(flip,ii)
          t[ii,] <- t[ii,c(2:1,4:3)]
        }
      }
    }
    for (i in seq_len(ncauses+1)) 
      for (j in seq(i,ncauses+1)) {
        ii <- which(t[,3]==causes[i] & t[,4]==causes[j])
        nn <- c(nn,paste(causes[i],causes[j]))
        idx <- c(idx, list(ii))
      }
    t. <- unique(sort(t[,1:2,drop=TRUE]))
    names(idx) <- nn    
  } else {
    causes <- c(0,sort(setdiff(unique(t[,2,drop=TRUE]),0)))
    ncauses <- length(causes)-1
    t. <- unique(sort(t[,1]))
    for (i in 1:(ncauses+1)) {
      nn <- c(nn,causes[i])
      idx <- c(idx, list(which(c(t[,2]==causes[i]))))
    }
    names(idx) <- nn
    dim <- 1
  }
  if (missing(knots))
    knots <- quantile(t.,c(.5))
  K <- c(rep(min(t.),d+1),knots,rep(max(t.),d+1)) #all knots(exterior+interior)
  B <- splinebasis(t.,degree=d,knots=knots,intercept=TRUE)
  dB <- attributes(B)$dB
  attributes(B)$dB <- NULL
  return(list(t=t,flip=flip,causes=causes,ncauses=ncauses,B=B,dB=dB,tB=t.,
              K=K,d=d,idx=idx,dim=dim))
}

###}}} cif

###{{{ cif methods

coef.cif <- function(object,...) object$opt$par
##vcov.cif <- function(object,...) object$H

print.cif <- function(x,...) {
  ##  ll <- logLik(x,info=TRUE)
  cat("\nVariance:\n"); print(x$Sigma)
  cat("\nCure-rates:\n"); print(x$P)
  cat("\nParameters:\n"); print(x$b)
  return(invisible(x))
}
summary.cif <- function(object,...) object

###}}} cif methods

###{{{ plot.cif
plot.cif <- function(x,col=c("seagreen","darkred","darkblue","goldenrod","mediumpurple"),curerate=TRUE,...) {
  t <- x$t
  if (length(x$flip)>0)
    t[x$flip,] <- t[x$flip,c(2,1,4,3)]
  plotcr(t,which=2,...)
  ll <- logLik(x,info=TRUE)
  mysd <- 1
  P <- as.matrix(ll$P)
  P1 <- rowSums(P)

  for (i in seq(x$ncauses)) {
    if (x$dim==2) mysd <- (ll$Sigma[2*(i-1)+1,2*(i-1)+1])^0.5
    y <- P1[i]*pnorm(ll$a[,i+1],sd=mysd)
    lines(ll$a[,1],y,col=col[i],lwd=2)
    if (curerate)
      abline(h=P1[i],col=Col(col[i],0.7),lty=2,lwd=0.5)
  }
  invisible(x)
}
###}}} plot.cif

###{{{ logLik.cif

logLik.cif <- function(object,theta=coef(object),info=FALSE,indiv=FALSE,...) {
  loglik <- function(p,d=degree,knots=knots) {
    with(object, logLikcif(p,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=FALSE,info=info,indiv=indiv,ppar=ppar,vpar=vpar,...))
  }
  if (object$dim==1) {
    loglik <- function(p,d=degree,knots=knots) {
      with(object, logLikcif1(p,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=FALSE,info=info,indiv=indiv,ppar=ppar,vpar=vpar,...))
    }
  }
  res <- loglik(theta,...)
  if (info) return(res)  
  n <- nrow(object$t)
  p <- length(coef(object))
  attr(res, "nall") <- n
  attr(res, "nobs") <- n-p
  attr(res, "df") <- p
  class(res) <- "logLik"        
  return(res)
}

score.cif <- function(x,theta=coef(x),...) {
  f <- function(p) logLik(x,theta=p)
  grad(f,theta)
}

###}}} logLik.cif

###{{{ logLikcif: univariate

logLikcif1 <-function(theta,indiv=FALSE,neg=TRUE,info=FALSE,##) {
                      ppar,vpar,
                      t,B,dB,K,causes,ncauses,tB,d,idx,...) {
  nb <- ncol(B)
  a1 <- da1 <- b <- db <- c()
  a <- tB
  b0 <- c()
  for (j in seq_len(ncauses)) {
    b0 <- theta[1+nb*(j-1)]
    if (d==1) {
      s <- 0
      for (i in seq_len(nb-1)+1) {
        b0 <- c(b0,exp(theta[i+nb*(j-1)])-s)
        s <- s+b0[i]
      }
      db0 <- b0[-1]
      da0 <- dB%*%db0
    } else {
      for (i in 1:(nb-1)) b0 <- c(b0,b0[i]+exp(theta[nb*(j-1)+i+1]))
      db0 <- diff(b0)*d/diff(K,lag=d)[2:length(b0)]
      da0 <- dB%*%db0    
    }
    b <- c(b, list(b0))
    db <- c(db, list(db0))
    a0 <- B%*%b0
    if (info)
      a <- cbind(a, a0)
    a1. <- fastapprox(tB,t[,1],a0)
    a1 <- c(a1, list(a1.$t))
    da1 <- c(da1, list(da0[a1.$pos+1]))
  }
  idx0 <- which(t[,2]==0) ## Right-censoring

  ## browser()
  ## Cure rate parameters
  pr <- do.call(ppar,list(theta=theta[(nb*ncauses+1):length(theta)],ncauses=ncauses,dim=1))
    
  sigma2 <- 1 
  res <- numeric(nrow(t))
  if (length(idx0)>0)
    S <- rep(0,length(idx0))
  for (i in seq_len(ncauses)) {
    idx1 <- which(t[,2]==causes[i+1])
    if (length(idx1)>0)
      res[idx1] <- log(pr[i]*da1[[i]][idx1]) + dnorm(a1[[i]][idx1],sd=sigma2^0.5,log=TRUE)
    if (length(idx0)>0)
      S <- S+pr[i]*pnorm(a1[[i]][idx0],sd=sigma2^0.5,lower.tail=FALSE,log=FALSE)
  }
  if (length(idx0)>0) {
    suppressWarnings(res[idx0] <- log(S))
    res[is.infinite(res)] <- -1e6
  }

  ll <- res  
  if (!indiv) ll <- sum(ll)
  res <- (ifelse(neg,-1,1)*ll)
  if (info) {
    res <- list(logLik=ll,a=a,Sigma=sigma2,b=b,P=pr)
    return(res)
  }
  return(res)
}

###}}} logLikcif1

###{{{ logLikcif: bivariate

logLikcif <- function(theta,indiv=FALSE,neg=TRUE,info=FALSE,
                      ppar,vpar,##) {
                      t,B,dB,K,causes,ncauses,tB,d,idx,...) {
  
  AD <- c()
  nb <- ncol(B)
  b <- db <- b0 <- c()

  a <- tB
  for (j in seq_len(ncauses)) {
    b0 <- theta[1+nb*(j-1)]
    if (d==0) {
      nb <- nb-1
      db0 <- 1; b0 <- c(theta[j],1)
    } else if (d==1) {
      s <- 0
      for (i in seq_len(nb-1)+1) {
        b0 <- c(b0,exp(theta[i+nb*(j-1)])-s)
        s <- s+b0[i]
      }
      db0 <- b0[-1] ## Derivative parameters
    } else {
      for (i in 1:(nb-1)) b0 <- c(b0,b0[i]+exp(theta[nb*(j-1)+i+1]))
      db0 <- diff(b0)*d/diff(K,lag=d)[2:length(b0)] ## Derivative parameters
    }
    da0 <- dB%*%db0 ## Derivative of parametrization of cause j
    b <- c(b, list(b0))
    db <- c(db, list(db0))
    a0 <- B%*%b0
    if (info) a <- cbind(a,a0)
    a1. <- fastapprox(tB,t[,1],a0)
    ##    af <- approxfun(t.,a)
    ##    a1.. <- af(t[,1])    
    a2. <- fastapprox(tB,t[,2],a0)
    ## nx4-matrix with columns a1(t1), a2(t2), a1'(t1), a2'(t2)
    AD <- c(AD, list(cbind(a1.$t,a2.$t,log(da0[a1.$pos+1]),log(da0[a2.$pos+1]))))
  }
  names(AD) <- causes[-1]

  theta1 <- theta[-seq(nb*ncauses)]
  Sigma <- do.call(vpar,list(theta=theta1,ncauses=ncauses))
  theta2 <- theta1[-seq(attributes(Sigma)$npar)]
  P <- do.call(ppar,list(theta=theta2,ncauses=ncauses))
  
  arglist <-  list(name="bicif",
                   n=nrow(t),ad=AD,ii=idx,causes=causes,m=0,S=Sigma,P=P,
                   DUP=FALSE) ##,PACKAGE="bptwin")
  res <- do.call(".Call",arglist)
  ll <- res$logLik
  ll[is.infinite(ll)] <- -1e6  
  if (!indiv) ll <- sum(ll)
  ll <- ifelse(neg,-1,1)*ll
  if (info) {
    res <- list(logLik=ll,a=a,Sigma=Sigma,b=b,P=P)
    return(res)
  }
  return(ll)
}

###}}} logLikcif

###{{{ Simulation

simcif <- function(n,theta=list(c(-10,1,0.1),c(-10,0.1,.2)),
                Sigma, pi,
                range=c(0,100),
                a=function(t,theta) theta[1]+theta[2]*t+theta[3]*exp(t),
                ia,
                cens
                ) {
  if (!is.list(theta)) theta <- list(theta)
  ncauses <- length(theta)
  if (missing(Sigma)) {
    Sigma <- diag(2*ncauses)+1
    for (i in seq_len(ncauses)) {
      Sigma[1:2+(i-1)*2,1:2+(i-1)*2] <- Sigma[1:2+(i-1)*2,1:2+(i-1)*2]+i
    }    
  }
  np <- (ncauses+1)*(ncauses)/2
  if (missing(pi)) {
    pr <- 1/2^seq_len(np)
  }
  P <- matrix(ncol=ncauses,nrow=ncauses)
  if (ncauses==1) {
    P[1,1] <- 1
  } else {
    P[upper.tri(P,diag=TRUE)] <- pr
    P[upper.tri(P)] <- P[lower.tri(P)] <- P[upper.tri(P)]/2
  }

  invexpl <- !missing(ia)
  t. <- seq(range[1],range[2],length.out=1000)
  pa <- matrix(ncol=1+ncauses,nrow=length(t.))
  pa[,1] <- t.

  simcauses <- rmultinom(1,n,as.vector(P))
  T <- Y <- matrix(ncol=4,nrow=n)
  N <- 0
  pos <- 0

  for (j in seq_len(ncauses)) {
    for (i in seq_len(ncauses)) {
      pos <- pos+1
      if (simcauses[pos]>0) {
        Nprev <- N+1
        N <- N + simcauses[pos]
        pseq <- seq(Nprev,N)
        T[pseq,3] <- Y[pseq,3] <- i
        T[pseq,4] <- Y[pseq,4] <- j
        idx <- c(1+(i-1)*2,2+(j-1)*2)
        S <- Sigma[idx,idx,drop=FALSE]
        val <- rnorm(simcauses[pos],sd=S[1,1]^0.5)
        cm <- CondMom(c(0,0),S,2,X=cbind(val))
        Y[pseq,1] <- val
        Y[pseq,2] <- with(cm, rnorm(simcauses[pos],mean=as.vector(mean),sd=var^0.5))
        
        theta0 <- theta[[i]];
        a. <- a(t.,theta0)
        if (!invexpl) {
          ia <- function(x,...) {
            fastapprox(a.,x,t.)$t
          }
        }
        T[pseq,1] <- ia(Y[pseq,1],theta0)
        theta0 <- theta[[j]]; # subject two
        a. <- a(t.,theta0)
        if (!invexpl) {
          ia <- function(x,...) {
            fastapprox(a.,x,t.)$t
          }
        }      
        T[pseq,2] <- ia(Y[pseq,2],theta0)
      }
    }
  }
##  browser()
  
  for (i in seq_len(ncauses)) { ## Marginals
    theta0 <- theta[[i]];
    a. <- a(t.,theta0)
    pa[,i+1] <- pnorm(a.,sd=Sigma[i*2,i*2]^0.5)
  }
  
  if (!missing(cens)) {    
    if (is.function(cens)) {
      cens <- cens(T[,1:2])
    }    
    T[T[,1]>=cens[,1],3] <- 0
    T[T[,2]>=cens[,2],4] <- 0
    T[T[,3]==0,1] <- cens[T[,3]==0,1]
    T[T[,4]==0,2] <- cens[T[,4]==0,2]
  }
  colnames(T) <- c("t1","t2","cause1","cause2")
  res <- list(data=T,prob=pa,var=Sigma,P=P)
  class(res) <- "simt"
  return(res)
}

print.simt <- function(x,...) {
  ##  cat("(t1,t2,cause1,cause2):\n")
  with(x, print(data))
  cat("Variance:\n")  
  with(x, print(var))
  cat("Prob. causes:\n")  
  with(x, print(P))
  invisible(x)
}

###}}} Simulation

###{{{ plotcr

plotcr <- function(x,col,legend=TRUE,which=1:2,
                   ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                   ...) {
  dots <- list(...)
  if ((!is.data.frame(x) | !is.matrix(x)) && ncol(x)<2) stop("Wrong type of data")
  if (ncol(x)==2) {
    co <- cuminc(x[,1],x[,2])
    if (any(x[,1]<0)) {
      for (i in seq(length(co))) {
        co[[i]]$time <- co[[i]]$time[-1]
        co[[i]]$est <- co[[i]]$est[-1]
      }
      if (is.null(dots$xlim)) dots$xlim <- range(x[,1])
    }
    do.call("plot", c(list(x=co), dots))
    return(invisible(co))
  }

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  t <- x[,1:2]; cause <- x[,3:4]
  causes <- sort(setdiff(unique(cause),0))
  if (missing(col)) col <- c("seagreen","darkred","darkblue","goldenrod","mediumpurple")
  if (1%in%which) {
    plot(t,type="n",...)
    count <- 1
    for (i in causes) {
      points(t[cause[,1]==causes[i],],col=Col(col[count],0.5),pch=2)
      points(t[cause[,2]==causes[i],],col=Col(col[count],0.5),pch=6)
      count <- count+1
    }
    points(t[cause[,1]==0 & cause[,2]==0,],col=Col("black",0.2),pch=1)  
    if (legend)
    legend("topleft", c("Subj 1, Cause 1", "Subj 2, Cause 1",
                        "Subj 1, Cause 2", "Subj 2, Cause 2",
                        "Double Censoring"), pch=c(2,6,2,6,1), col=c(rep(col[1:length(causes)],each=2),"black"))
  }
  if (2%in%which) {
    co1 <- cuminc(t[,1],cause[,1])
    co2 <- cuminc(t[,2],cause[,2])
    if (any(t[,1]<0)) {
      for (i in seq(length(co1))) {
        co1[[i]]$time <- co1[[i]]$time[-1]
        co1[[i]]$est <- co1[[i]]$est[-1]
      }
      for (i in seq(length(co2))) {
        co2[[i]]$time <- co1[[i]]$time[-1]
        co2[[i]]$est <- co1[[i]]$est[-1]
      }
    }
    if (is.null(dots$xlim)) dots$xlim <- range(t)
    do.call("plot", c(list(x=co1), dots))
    for (i in seq(length(co2)))
      with(co2[[i]], lines(time,est,type="s",lty=i,...))
  }
}

###}}} plotcr

###{{{ nonparcuminc
nonparcuminc <- function(t,status,cens=0) {
  ord <- order(t); t <- t[ord]; status <- status[ord]
  ud.cens<-survfit(Surv(t,status==cens)~1)
  Gfit<-cbind(ud.cens$time,ud.cens$surv)
  Gfit<-rbind(c(0,1),Gfit);
  causes <- setdiff(unique(status),cens)
  cweight<-fastapprox(Gfit[,1],t,Gfit[,2])[[1]];
  cc <- t
  for (i in 1:length(causes)) {
    c1 <- status==causes[i]
    cc <- cbind(cc,cumsum(c1/cweight)/length(c1))
  }
  return(cc)
}

###}}} nonparcuminc

###{{{ splinebasis

splinebasis <- function(t,knots=quantile(t,c(0.5)),
                        degree=1,intercept=TRUE) {
  if (degree==1) {
    B <- matrix(1,ncol=length(knots)+2,nrow=length(t))
    dB <- B[,-1,drop=FALSE]
    B[,2] <- t
    for (i in seq_len(length(knots))) {
      B[,i+2] <- pmax(t-knots[i],0)
      dB[B[,i+2]==0,i+1] <- 0
    }
  } else {
    B <- bs(t,degree=degree,knots=knots,intercept=intercept)
    dB <- bs(t,degree=degree-1,knots=knots,intercept=intercept)
  }
  if (!intercept) B <- B[,-1,drop=FALSE]
  attributes(B)$dB <- dB
  return(B)
}

###}}} splinebasis

###{{{ CondMom
CondMom <- function(mu,S,idx,X) {
  if (ncol(S)==length(idx)) {
    return(list(mean=mu,var=S))
  }  
  idxY <- idx  
  idxX <- setdiff(1:ncol(S),idxY)
  SXX <- S[idxX,idxX,drop=FALSE];
  SYY <- S[idxY,idxY,drop=FALSE]
  SYX <- S[idxY,idxX,drop=FALSE]
  iSXX <- solve(SXX)
  if (is.matrix(mu)) {
    muY <- mu[,idxY,drop=FALSE]
    muX <- mu[,idxX,drop=FALSE]
    Z <- t(X-muX)
  } else {
    muY <- mu[idxY]
    muX <- mu[idxX]
    Z <- apply(X,1,function(xx) xx-muX)
  }
  SZ  <- t(SYX%*%iSXX%*%Z)
  if (is.matrix(mu))
    condmean <- SZ+muY
  else
    condmean <- t(apply(SZ,1,function(x) muY+x))
  condvar <- SYY-SYX%*%iSXX%*%t(SYX)
  return(list(mean=condmean,var=condvar))
}

###}}} CondMom

###{{{ probpar

probpar <- function(theta,ncauses,dim=2,...) {
  if (is.null(theta)) return(ncauses*(ncauses-1)/2+ncauses-1)

  expit <- function(z,b=1) (b)/(1+exp(-z)) ## y\in(0,b)  

  if (dim==1) {
    pr <- rep(1,ncauses)
    if (ncauses>1 & length(theta)>0) {      
      npr <- ncauses
##    if (length(idx0)==0) npr <- npr-1
      A <- pos <- 0
      pr <- c()
      for (i in seq_len(npr-1)) {
        pos <- pos+1
        ##        pr0 <- expit(theta[nb*ncauses+i],b=1-A)
        pr0 <- expit(theta[pos],b=1-A)
        pr <- c(pr,pr0)
        A <- A+pr0
      }
      pr <- c(pr,1-A)
    }
    return(pr)
  }
  
  P <- matrix(1,ncol=ncauses,ncauses) ## Cure rate matrix 
  A <- 0
  if (ncauses>1 & length(theta)>0) {
    npr <- ncauses*(ncauses-1)/2
    pos <- 0
   for (i in seq_len(ncauses-1))
      for (j in seq_len(ncauses-i+1)) {        
        J <- j+i-1
        ##        cat(i,",",J,"\n")
        pos <- pos+1
        ##        pr0 <- expit(theta[nb*ncauses+pos+npvar],b=1-A)
        pr0 <- expit(theta[pos],b=1-A)
        if (i!=J) {
          P[i,J] <- P[J,i] <- pr0/2
       } else {
         P[i,i] <- pr0
       }        
        A <- A+pr0
      }
  } 
  P[ncauses,ncauses] <- 1-A
  return(P)
}

probpar0 <- function(theta,ncauses,...) {
  if (is.null(theta)) return(ncauses*(ncauses-1)/2+ncauses)

  P <- matrix(1,ncol=ncauses,ncauses) ## Cure rate matrix 
  expit <- function(z,b=1) (b)/(1+exp(-z)) ## y\in(0,b)  
  A <- 0
  ##  if (ncauses>1 & !(length(theta)<=nb*ncauses+npvar))
  if (ncauses>1 & length(theta)>0) {
    npr <- ncauses*(ncauses-1)/2
    pos <- 0
    for (i in seq_len(ncauses))
##      for (j in seq_len(ncauses-i+1)) {
      for (J in seq(i,ncauses)) {
##        J <- j+i-1
##        cat(i,",",J,"\n")
        pos <- pos+1
        ##  if (length(idx10)+length(idx20)==0) npr <- npr-1
        ##        for (i in seq_len(npr)) {
  ##      pr0 <- expit(theta[nb*ncauses+pos+npvar],b=1-A)
        pr0 <- expit(theta[pos],b=1-A)        
        if (i!=J) {
          P[i,J] <- P[J,i] <- pr0/2
       } else {
         P[i,i] <- pr0
       }        
        A <- A+pr0
      }
  } 
##  P[ncauses,ncauses] <- 1-A
  return(P)
}

###}}} probpar

###{{{ varpar

varpar1 <- function(theta,ncauses,...) {
  if (is.null(theta)) return(ncauses+1)
  Sigma <- diag(1,nrow=2*ncauses)
  for (i in seq_len(ncauses)) {
    Sigma[1:2+2*(i-1),1:2+2*(i-1)] <-  Sigma[1:2+2*(i-1),1:2+2*(i-1)]+
      exp(theta[i])
##        exp(theta[nb*ncauses+i])
  }
  npvar <- ncauses
  if (ncauses>1) {
    ##    browser()
    Sigma <- Sigma+exp(theta[ncauses+1])
    Sigma <- Sigma+theta[ncauses+1]
##    Sigma[3:4,1:2] <- Sigma[1:2,3:4] <- (theta[ncauses+1])
    npvar <- ncauses+1
  }
  attributes(Sigma)$npar <- npvar
  return(Sigma)
}  

###}}} varpar
