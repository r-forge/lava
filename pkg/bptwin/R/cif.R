###{{{ cif

cif <- function(t,X,start,knots,degree=1,theta,control=list(trace=1),
                ppar="prob",vpar="block",cpp="bicif",...) {
  PP <- prepcif(t,d=degree,knots=knots)
  ppar <- paste(ppar,"cif.ppar",sep="_")
  vpar <- paste(vpar,"cif.vpar",sep="_")
  loglik <- function(theta,d=degree,knots=knots,...) {
    with(PP, logLikcif(theta,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=TRUE,ppar=ppar,vpar=vpar,cpp=cpp,...))
  }
  if (PP$dim==1) {
    loglik <- function(theta,d=degree,knots=knots,...) {
      with(PP, logLikcif1(theta,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=TRUE,ppar=ppar,vpar=vpar,...))
    }
    if (missing(start)) {
      npar <- with(PP, ncol(B)*ncauses+ncauses-1)
      start <- rep(0,npar)
    }
  } else {
    if (missing(start)) {
      ## theta0 <-
      ##   with(PP, c(rep(0,ncol(B)*ncauses),
      ##              do.call(vpar,list(theta=NULL,ncauses=ncauses,
      ##                                start=TRUE)),
      ##              rep(0,do.call(ppar,list(theta=NULL,ncauses=ncauses)))),0,0,0)
      npar <- with(PP, ncol(B)*ncauses+
                   do.call(vpar,list(theta=NULL,ncauses=ncauses))+
                   do.call(ppar,list(theta=NULL,ncauses=ncauses)))
      start <- rep(0,npar)
      
    }
  }

  if (!missing(theta)) return(loglik(theta,...))
  opt <- nlminb(start,loglik,control=control)
  res <- c(list(coef=opt$par, opt=opt,
                ppar=ppar,vpar=vpar),
           PP)
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
              knots=knots,K=K,d=d,idx=idx,dim=dim))
}

###}}} cif

###{{{ cif methods

predict.cif <- function(object,cause=c(1,1),t,...) {
  l <- logLik(object,info=TRUE,...)
  if (missing(t)) t <- seq(min(l$a[,1]),max(l$a[,1]),length.out=100)

  B <- splinebasis(t,knots=object$knots,degree=object$d,Boundary.knots=object$K[c(1,length(object$K))])
  a <- cbind(t,B%*%l$b[[1]],B%*%l$b[[2]])
  
  if (object$dim & length(cause)==2) {
    idx <- c((cause[1]-1)*2+1,(cause[2]-1)*2+2)
    S <- object$Sigma[idx,idx]
    p <- object$P[cause[1],cause[2]]*apply(a[,cause+1],1,function(x) pmvnorm(upper=x,sigma=S))
##    p <- object$P[cause[1],cause[2]]*apply(l$a[,cause+1],1,function(x) pmvnorm(lower=x,sigma=S))
  } else {
    S <- object$Sigma[(cause[1]-1)*2+1,(cause[1]-1)*2+1]
    p <- sum(object$P[cause[1],])*pnorm(upper=x,sigma=S^0.5)
  }
  return(cbind(a[,1],p))  
}

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
  P[lower.tri(P)] <- P[upper.tri(P)]/2
  P[upper.tri(P)] <- P[lower.tri(P)]
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
    cat(".")
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
                      t,B,dB,K,causes,ncauses,tB,d,idx,cpp="bicif",...) {
  
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

  arglist <-  list(name=cpp,
                   n=nrow(t),ad=AD,ii=idx,causes=causes,m=0,S=Sigma,P=P,
                   DUP=FALSE) ##,PACKAGE="bptwin")
  res <- do.call(".Call",arglist)
  ll <- res$logLik
  ll[is.infinite(ll)] <- -1e6
  if (!indiv) ll <- sum(ll)
  if (!is.null(attributes(Sigma)$psd)) {
    ll <- ll - 1e4
  }

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
                Sigma, pr,
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
  if (missing(pr)) {
    pr <- 1/2^seq_len(np-1)
  }
  pr <- c(pr,1-sum(pr))
  P <- matrix(ncol=ncauses,nrow=ncauses)
  if (ncauses==1) {
    P[1,1] <- 1
  } else {    
    P[upper.tri(P,diag=TRUE)] <- pr
    P[lower.tri(P)] <- P[upper.tri(P)] <- P[upper.tri(P)]/2
  }

  invexpl <- !missing(ia)
  t. <- seq(range[1],range[2],length.out=1000)
  pa <- matrix(ncol=1+ncauses,nrow=length(t.))
  pa[,1] <- t.
  aa <- pa

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
    aa[,i+1] <- a.
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
  res <- list(data=T,prob=pa,var=Sigma,P=P,a=aa)
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

npc <- function(T,cause,same.cens=TRUE,sep=FALSE) {
  mtime <- apply(T[,1:2],1,max)
  ot <- order(mtime)
  mtime <- mtime[ot]
  T <- T[ot,]
  if (!sep) {
    time1 <- as.vector(T[,1:2]); status1 <- as.vector(T[,3:4])
    ud.cens1<-survfit(Surv(time1,status1==0)~+1);
    Gfit1<-cbind(ud.cens1$time,ud.cens1$surv)
    Gfit2 <- Gfit1<-rbind(c(0,1),Gfit1);
  } else {
    time1 <- as.vector(T[,1]); status1 <- as.vector(T[,3])
    ud.cens1<-survfit(Surv(time1,status1==0)~+1);
    time2 <- as.vector(T[,2]); status2 <- as.vector(T[,4])
    ud.cens2<-survfit(Surv(time2,status2==0)~+1);
    Gfit1<-cbind(ud.cens1$time,ud.cens1$surv)
    Gfit1<-rbind(c(0,1),Gfit1);
    Gfit2<-cbind(ud.cens2$time,ud.cens2$surv)
    Gfit2<-rbind(c(0,1),Gfit2);
  }
  cweights1<-fastapprox(Gfit1[,1],T[,1],Gfit1[,2])[[1]]
  cweights2<-fastapprox(Gfit2[,1],T[,2],Gfit2[,2])[[1]];
  weight11 <- apply(cbind(cweights1,cweights2),1,min)

  if (same.cens) {
    conc <- (T[,3]==cause[1])*(T[,4]==cause[2])/weight11
  } else {
    conc <-(T[,3]==cause[1])*(T[,4]==cause[2])/(cweights1*cweights2);
  }
  mtime <- mtime[!is.na(conc)]
  conc <- conc[!is.na(conc)]
  cbind(mtime,cumsum(conc)/length(conc))
}


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
                        degree=1,intercept=TRUE,...) {
  if (degree==1) {
    B <- matrix(1,ncol=length(knots)+2,nrow=length(t))
    dB <- B[,-1,drop=FALSE]
    B[,2] <- t
    for (i in seq_len(length(knots))) {
      B[,i+2] <- pmax(t-knots[i],0)
      dB[B[,i+2]==0,i+1] <- 0
    }
  } else {
    B <- bs(t,degree=degree,knots=knots,intercept=intercept,...)
    dB <- bs(t,degree=degree-1,knots=knots,intercept=intercept,...)
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

###{{{ ppar

prob_cif.ppar <- function(theta,ncauses,...) {
##  if (is.null(theta)) return(ncauses*(ncauses-1)/2+ncauses)
  if (is.null(theta)) return(ncauses*(ncauses-1)/2+ncauses)
  
  P <- matrix(0,ncol=ncauses,ncauses) ## Cure rate matrix 
  expit <- function(z,b=1) (b)/(1+exp(-z)) ## y\in(0,b)  
  A <- 0
  ##  if (ncauses>1 & !(length(theta)<=nb*ncauses+npvar))
  if (ncauses>1 & length(theta)>0) {
    npr <- ncauses*(ncauses-1)/2
    pos <- 0
    for (i in seq_len(ncauses))
      for (J in seq(i,ncauses)) {
##      for (J in seq_len(ncauses)) {
        pos <- pos+1
        pr0 <- expit(theta[pos],b=1-A)
        P[i,J] <- pr0
        A <- A+pr0
      }
  }
##  browser()
  return(P)
}

prob1_cif.ppar <- function(theta,ncauses,dim=2,...) {
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
  
  P <- matrix(0,ncol=ncauses,ncauses) ## Cure rate matrix 
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
          P[i,J]<- pr0
       } else {
         P[i,i] <- pr0
       }        
        A <- A+pr0
      }
  } 
  P[ncauses,ncauses] <- 1-A
  return(P)
}




probObs_cif.ppar <- function(theta,ncauses,dim=2,...) {
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

prob0Obs_cif.ppar <- function(theta,ncauses,...) {
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

###}}} ppar (Parametrization of cure rates)

###{{{ vpar (Parametrization of covariance parameters)

block2_cif.vpar <- function(theta,ncauses,start=FALSE,...) {
  npvar <- ncauses+1
  if (start) return(c(rep(1,npvar-1),0))
  if (is.null(theta)) return(npvar)
  Sigma <- diag(1,nrow=2*ncauses)  
  for (i in seq_len(ncauses)) {
    Sigma[1:2+2*(i-1),1:2+2*(i-1)] <-  Sigma[1:2+2*(i-1),1:2+2*(i-1)]+
      exp(theta[i])
  }
  npvar <- ncauses
  if (ncauses>1) {
##    Sigma <- Sigma+exp(theta[ncauses+1])
    Sigma[1:2,3:4] <- Sigma[3:4,1:2] <- theta[ncauses+1]
    npvar <- ncauses+1
  }
  ## ee <- eigen(Sigma)
  ## if (any(ee$values<0)) {
  ##   ee$values[ee$values<0] <- 0
  ##   Sigma <- with(ee,vectors%*%diag(values)%*%t(vectors))
  ##   attributes(Sigma)$psd <- FALSE
  ## }
  attributes(Sigma)$npar <- npvar
  return(Sigma)
}  

block0_cif.vpar <- function(theta,ncauses,start=FALSE,...) {
  res <- block0_cif.vpar(theta,ncauses,start=start,...)
  if (is.null(theta) | start) return(res)
  else res[1,3] <- res[3,1] <- res[2,4] <- res[4,2] <- 0
  return(res)
}


block_cif.vpar <- function(theta,ncauses,start=FALSE,...) {
  npvar <- ncauses+1
  if (start) return(rep(0,npvar))
  if (is.null(theta)) return(npvar)
  Sigma <- diag(1,nrow=2*ncauses)
  if (ncauses>1) {
    Sigma <- Sigma+tanh(theta[ncauses+1])
    npvar <- ncauses+1
  }
  diag(Sigma) <- 1 ##diag(1,nrow=2*ncauses)
  for (i in seq_len(ncauses)) {
    Sigma[1+2*(i-1),2+2*(i-1)] <- Sigma[2+2*(i-1),1+2*(i-1)] <- 
      tanh(theta[i])
  }
  attributes(Sigma)$npar <- npvar
  return(Sigma)
}



block1_cif.vpar <- function(theta,ncauses,start=FALSE,...) {
  npvar <- ncauses+1
  if (start) return(rep(0,npvar))
  if (is.null(theta)) return(npvar)
  Sigma <- diag(1,nrow=2*ncauses)
  for (i in seq_len(ncauses)) {
    Sigma[1:2+2*(i-1),1:2+2*(i-1)] <-  Sigma[1:2+2*(i-1),1:2+2*(i-1)]+
      exp(theta[i])
  }
  npvar <- ncauses
  if (ncauses>1) {
    Sigma <- Sigma+exp(theta[ncauses+1])
    npvar <- ncauses+1
  }
  attributes(Sigma)$npar <- npvar
  return(Sigma)
}

block3_cif.vpar <- function(theta,ncauses,...) {
  npvar <- ncauses + (ncauses-1)*ncauses/2  
  if (is.null(theta)) return(ncauses)
  Sigma <- diag(1,nrow=2*ncauses)
  for (i in seq_len(ncauses)) {
    Sigma[1:2+2*(i-1),1:2+2*(i-1)] <-  Sigma[1:2+2*(i-1),1:2+2*(i-1)]+
      exp(theta[i])
  }
  pos <- ncauses
  for (i in seq_len(ncauses-1)) {
    for (j in seq_len(ncauses-i)+i) {
      pos <- pos+1
      Sigma[cbind(2:1+2*(i-1),1:2+2*(j-1))] <-
        Sigma[cbind(1:2+2*(j-1),2:1+2*(i-1))] <- exp(theta[pos])
    }
  }
  attributes(Sigma)$npar <- npvar
  return(Sigma)
}

unstruct_cif.vpar <- function(theta,ncauses,...) {
  npvar <- (2*ncauses-1)*ncauses
  if (is.null(theta)) return(npvar)

  L <- diag(1,nrow=2*ncauses)
  pos <- 0
  for (i in seq_len(2*ncauses-1)) {
    for (j in seq_len(2*ncauses-i)+i) {
      pos <- pos+1
      L[j,i] <- theta[pos]
    }
  }
  Sigma <- crossprod(L) 
  attributes(Sigma)$npar <- npvar
  return(Sigma)
}  

###}}} varpar
