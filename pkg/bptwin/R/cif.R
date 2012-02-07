###{{{ cif

cif <- function(t,X=NULL,XX=NULL,Z=NULL,start,knots,degree=1,theta,control=list(trace=1),
                ppar="prob",vpar="block",hessian=FALSE,cpp="bicif",sym=TRUE,flex=FALSE,...) {
  PP <- prepcif(t,d=degree,knots=knots,sym=sym,flex=flex)
  ppar <- paste(ppar,"cif.ppar",sep="_")
  vpar <- paste(vpar,"cif.vpar",sep="_")
  nx <- ifelse(is.null(X),0,ncol(X))

  loglik <- function(theta,d=degree,knots=knots,...) {
    with(PP, logLikcif(theta,X1=X,X2=X,XX=XX,Z=Z,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=TRUE,ppar=ppar,vpar=vpar,cpp=cpp,sym=sym,...))
  }
  if (PP$dim==1) {
    loglik <- function(theta,d=degree,knots=knots,...) {
      with(PP, logLikcif1(theta,X1=X,X2=X,Z=Z,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=TRUE,ppar=ppar,vpar=vpar,...))
    }
    if (missing(start)) {
      npar <- with(PP, ncol(B)*ncauses+ncauses-1+nx*ncauses)
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
                   do.call(vpar,list(theta=NULL,ncauses=ncauses,sym=sym)) +
                   do.call(ppar,list(theta=NULL,ncauses=ncauses,sym=sym))) +
                     PP$ncauses*nx
      start <- rep(0,npar)
      
    }
  }

  if (!missing(theta)) return(loglik(theta,...))
  opt <- nlminb(start,loglik,control=control)
  res <- c(list(coef=opt$par, opt=opt,
                ppar=ppar,vpar=vpar,X=X,Z=Z),
           PP)
  class(res) <- "cif"
  ll <- function(p,...) logLik(res,theta=p,...)
  if (hessian) {
    H <- -hessian(ll,coef(res))
    res$information <- H
  }
  l0 <- ll(coef(res),info=TRUE)
  res <- c(res, list(Sigma=l0$Sigma, b=l0$b, b1=l0$b1, P=l0$P, P1=rowSums(as.matrix(l0$P)), sym=sym))
  class(res) <- "cif"
  return(res)

}

prepcif <- function(t,knots,d=1,sym=TRUE,flex=FALSE,...) {
  idx <- nn <- flip <- c()
  truecause <- NULL
  dim <- 2  
  if (ncol(t)==4) {
    causes <- c(0,sort(setdiff(unique(t[,3:4,drop=TRUE]),0)))
    ncauses <- length(causes)-1
    if (!sym) {
      for (j in seq_len(ncauses)+1) {
        t[t[,4]==causes[j],4] <- causes[j]+100 
      }
      causes <- c(causes,setdiff(causes,0)+100)
      ncauses <- ncauses*2
    }
    { ## Symm.
      nn <- c()
      for (i in seq_len(ncauses+1)) {
        for (j in seq(i+1,ncauses+1)) {
          ii <- which(t[,4]==causes[i] & t[,3]==causes[j])
          if (length(ii)>0)  {
            flip <- c(flip,ii)
            t[ii,] <- t[ii,c(2:1,4:3)]
          }
        }
      }

      ## if (flex) { ## Non-homogeneous case
      ##   count <- 0
      ##   truecause <- cbind(causes[-1],causes[-1])
      ##   for (i in seq_len(ncauses-1)) {
      ##     for (j in seq_len(ncauses-i)+i) {            
      ##       idx <- which(t[,3]==causes[i+1] & t[,4]==causes[j+1])
      ##       if (length(idx)>0) {
      ##         count <- count+1
      ##         t[idx,3:4] <- cbind(rep(causes[i+1]+10*count,length(idx)),
      ##                             rep(causes[j+1]+10*count,length(idx)))
      ##         truecause <- rbind(truecause, c(causes[i+1],causes[i+1]+10*count))
      ##         truecause <- rbind(truecause, c(causes[j+1],causes[j+1]+10*count))
      ##         ncauses <- ncauses+2
      ##         causes <- c(causes,causes[i+1]+10*count,causes[j+1]+10*count)
      ##         }
      ##     }
      ##   }
      ## }      
      
      for (i in seq_len(ncauses+1)) 
        for (j in seq(i,ncauses+1)) {
          ii <- which(t[,3]==causes[i] & t[,4]==causes[j])
          nn <- c(nn,paste(causes[i],causes[j]))
          idx <- c(idx, list(ii))
        }
    }
    
    t. <- unique(sort(t[,1:2,drop=TRUE]))
    names(idx) <- nn    
  } else { ## Univariate case
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
              knots=knots,K=K,d=d,idx=idx,dim=dim,sym=sym,truecause=truecause))
}

###}}} cif

###{{{ cif methods

predict.cif <- function(object,cause=c(1,1),t,...) {
  l <- logLik(object,info=TRUE,...)
  if (missing(t)) t <- seq(min(l$a[,1]),max(l$a[,1]),length.out=100)

  B <- splinebasis(t,knots=object$knots,degree=object$d,Boundary.knots=object$K[c(1,length(object$K))])
  a <- cbind(t,B%*%l$b[[cause[1]]],B%*%l$b[[cause[2]]])

  if (object$dim & length(cause)==2) {
    idx <- c((cause[1]-1)*2+1,(cause[2]-1)*2+2)
    S <- object$Sigma[idx,idx]
    p <- object$P[cause[1],cause[2]]*apply(a[,2:3],1,function(x) pmvnorm(upper=x,sigma=S))
##    p <- object$P[cause[1],cause[2]]*apply(l$a[,cause+1],1,function(x) pmvnorm(lower=x,sigma=S))
  } else {
    S <- object$Sigma[(cause[1]-1)*2+1,(cause[1]-1)*2+1]
    p <- sum(object$P[cause[1],])*pnorm(upper=x,sd=S^0.5)
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
  cat("\nMean parameters:\n"); print(x$b1)
  return(invisible(x))
}
summary.cif <- function(object,...) object

###}}} cif methods

###{{{ plot.cif

plot.cif <- function(x,col=c("seagreen","darkred","darkblue","goldenrod","mediumpurple"),curerate=TRUE,cause,nonpar=TRUE,add=FALSE,legend=!add,bounds=TRUE,...) {  
  dots <- list(...)
  if (is.null(dots$xlim)) dots$xlim <- range(x$t[,1:2])
  if (is.null(dots$ylim)) dots$ylim <- c(0,1)
  if (is.null(dots$xlab)) dots$xlab <- "Time"
  if (is.null(dots$ylab)) dots$ylab <- "Probability"
  myargs <- c(x=0,y=0,type="n",dots)
    
  if (!missing(cause) && (class(cause)=="logical" || ncol(cause)==2)) {
    if (class(cause)=="logical") {
      cause <- mycol <- c()
      for (i in seq(x$ncauses)) {        
        cause <- c(cause,rep(i,x$ncauses-i+1))
        mycol <- c(mycol,(i-1)+seq(x$ncauses-i+1))
      }
      cause <- cbind(cause,mycol)
    }
    cause <- rbind(cause)
    if (!add)
      do.call("plot",myargs)
    for (i in seq(nrow(cause))) {
      P <- predict(x,cause=cause[i,])
      lines(P,col=col[i],lwd=2)
      NP <- npc(x$t,cause[i,])
      pos <- unique(fastapprox(NP[,1],P[,1],NP[,2])$pos)+1
      lines(NP[pos,],col=col[i],lty=2,type="s")
      if (curerate)
        abline(h=x$P[cause[i,,drop=FALSE]],col=Col(col[i],0.7),lty=2,lwd=0.5)
      if (bounds) {
        whichcause <- match(cause[i,],setdiff(x$causes,0))
        ll <- logLik(x,info=TRUE)
        P <- as.matrix(ll$P)    
        P[lower.tri(P)] <- P[upper.tri(P)]/2
        P[upper.tri(P)] <- P[lower.tri(P)]
        P1 <- rowSums(P)
        y <- P1[whichcause[1]]*pnorm(ll$a[,whichcause[1]+1],sd=1)
        if (cause[i,1]==cause[i,2])
          lines(ll$a[,1],y,col=col[i],lwd=1,lty=3)
        y <- y*P1[whichcause[2]]*pnorm(ll$a[,whichcause[2]+1],sd=1)
        lines(ll$a[,1],y,col=col[i],lwd=1,lty=3)
      }
    }
    labs <- apply(cause,1,function(x) paste(x,collapse=","))
    if (legend)
      legend("topleft",labs,col=col[seq(nrow(cause))],lty=1,pch=-1,box.lwd=0)
    return(invisible(x))
  }

  t <- x$t
  if (length(x$flip)>0)
    t[x$flip,] <- t[x$flip,c(2,1,4,3)]
  if (nonpar) {
    plotcr(t,which=2,...)
  } else {
    if (!add)
      do.call("plot",myargs)
  }
  
  ll <- logLik(x,info=TRUE)
  mysd <- 1
  P1 <- ll$P
  if (x$dim==2) {
    P <- as.matrix(ll$P)    
    P[lower.tri(P)] <- P[upper.tri(P)]/2
    P[upper.tri(P)] <- P[lower.tri(P)]
    P1 <- rowSums(P)
  }
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
    with(object, logLikcif(p,X1=X,Z=Z,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=FALSE,info=info,indiv=indiv,ppar=ppar,vpar=vpar,sym=sym,...))
  }
  if (object$dim==1) {
    loglik <- function(p,d=degree,knots=knots) {
      with(object, logLikcif1(p,X1=X,Z=Z,t=t,B=B,dB=dB,K=K,causes=causes,ncauses=ncauses,tB=tB,d=d,idx=idx,neg=FALSE,info=info,indiv=indiv,ppar=ppar,vpar=vpar,...))
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
                      t,B,dB,X1=NULL,X2=X1,Z=NULL,K,causes,ncauses,tB,d,idx,...) {
  
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

  M <- c()
  curpos <- length(unlist(b))
  ##  curpos <- nb*ncauses;
  b1 <- c()##b2 <- c()
    for (i in seq_len(ncauses)) {
    if (!is.null(X1)) {      
      bcur <- theta[curpos+seq(ncol(X1))]
      b1 <- c(b1,list(bcur))
      ##      b2 <- c(b2,list(theta[curpos+ncol(X1)+ncol(X2)]))
      curpos <- curpos+ncol(X1)##+ncol(X2)
      M <- c(M,list(as.vector(X1%*%bcur)))
    } else {
      M <- c(M,list(rep(0,nrow(t))))
    }
  }
##  browser()
  
  ## Cure rate parameters
  pr <- do.call(ppar,list(theta=theta[(curpos+1):length(theta)],ncauses=ncauses,dim=1))
    
  sigma2 <- 1 
  res <- numeric(nrow(t))
  if (length(idx0)>0)
    S <- rep(0,length(idx0))
  for (i in seq_len(ncauses)) {
    idx1 <- which(t[,2]==causes[i+1])
    if (length(idx1)>0)
      res[idx1] <- log(pr[i]*da1[[i]][idx1]) + dnorm(a1[[i]][idx1]-M[[i]][idx1],sd=sigma2^0.5,log=TRUE)
    if (length(idx0)>0)
      S <- S+pr[i]*pnorm(a1[[i]][idx0]-M[[i]][idx0],sd=sigma2^0.5,lower.tail=FALSE,log=FALSE)
  }
  if (length(idx0)>0) {
    ##message(".")
    suppressWarnings(res[idx0] <- log(S))
    res[is.infinite(res)] <- -1e6
  }

  ll <- res  
  if (!indiv) ll <- sum(ll)
  res <- (ifelse(neg,-1,1)*ll)
  if (info) {
    res <- list(logLik=ll,a=a,Sigma=sigma2,b=b,P=pr,b1=b1)
    return(res)
  }
  return(res)
}

###}}} logLikcif1

###{{{ logLikcif: bivariate

logLikcif <- function(theta,indiv=FALSE,neg=TRUE,info=FALSE,
                      ppar,vpar,##) {
                      X1=NULL,X2=X1,Z=NULL,V=NULL,sym=TRUE,truecause=NULL,
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
    ## n x 4-matrix with columns a1(t1), a2(t2), a1'(t1), a2'(t2)
    AD <- c(AD, list(cbind(a1.$t,a2.$t,log(da0[a1.$pos+1]),log(da0[a2.$pos+1]))))
  }
  names(AD) <- causes[-1]  
  curpos <- nb*ncauses

  M1 <- b1 <- c()##b2 <- c()
    for (i in seq_len(ncauses)) {
    if (!is.null(X1)) {      
      bcur <- theta[curpos+seq(ncol(X1))]
      b1 <- c(b1,list(bcur))
      ##      b2 <- c(b2,list(theta[curpos+ncol(X1)+ncol(X2)]))
      curpos <- curpos+ncol(X1)##+ncol(X2)
      Mcur1 <- as.vector(X1%*%bcur)
      Mcur2 <- as.vector(X2%*%bcur)
      M1 <- c(M1,list(Mcur))      
      AD[[i]][,1] <- AD[[i]][,1]-Mcur
      AD[[i]][,2] <- AD[[i]][,2]-Mcur
    } else {
      M1 <- c(M1,list(rep(0,nrow(t))))
    }
  }
  
  theta1 <- theta[-seq(curpos)]
  Sigma <- do.call(vpar,list(theta=theta1,ncauses=ncauses,sym=sym))
  lambda <- eigen(Sigma)$values
  if (any(lambda>1e9)) stop("Variance matrix out of bounds")  
  theta2 <- theta1[-seq(attributes(Sigma)$npar)]
  P <- do.call(ppar,list(theta=theta2,ncauses=ncauses,sym=sym))

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
    res <- list(logLik=ll,a=a,Sigma=Sigma,b=b,db=db,b1=b1,P=P)
    return(res)
  }
  return(ll)
}

###}}} logLikcif

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

prob_cif.ppar <- function(theta,ncauses,dim=2,sym=TRUE,...) {
##  if (is.null(theta)) return(ncauses*(ncauses-1)/2+ncauses)
  if (is.null(theta)) {
    if (sym)
      return(ncauses*(ncauses-1)/2+ncauses)
    return((ncauses/2)^2)
  }
  
  P <- matrix(0,ncol=ncauses,ncauses) ## Cure rate matrix
  P[1,1] <- 1
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

  A <- 0
  P[] <- 0
  ##  if (ncauses>1 & !(length(theta)<=nb*ncauses+npvar))
  if (ncauses>1 & length(theta)>0) {
    if (!sym) {
      npr <- (ncauses/2)^2
      pos <- 0
      for (i in seq_len(ncauses/2)) {
        for (J in seq_len(ncauses/2)+ncauses/2) {
          pos <- pos+1
          pr0 <- expit(theta[pos],b=1-A)
          P[i,J] <- pr0
          A <- A+pr0          
        }
      }
    } else {
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
  }
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

block2_cif.vpar <- function(theta,ncauses,start=FALSE,sym=TRUE,...) {
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


block_cif.vpar <- function(theta,ncauses,start=FALSE,sym=TRUE,...) {
  npvar <- ifelse(sym,ncauses+ncauses*(ncauses-1)/2,ncauses/2+1)
  if (start) return(rep(0,npvar))
  if (is.null(theta)) return(npvar)
  Sigma <- diag(1,nrow=2*ncauses)
  if (sym) {
    pos <- ncauses
    for (i in seq_len(ncauses)) {
      
      if (i<ncauses) {
        for (j in seq_len(ncauses-i)+i) {
          pos <- pos+1
          Sigma[rbind(c(2*i-1,2*j),
                      c(2*i,2*j-1),
                      c(2*j-1,2*i),
                      c(2*j,2*i-1))] <- 
                        ##          Sigma[2*i+c(-1,0),2*j+c(-1,0)] <-
                        ##            Sigma[2*j+c(-1,0),2*i+c(-1,0)] <-
                        tanh(theta[pos])
        }
      }
      Sigma[1+2*(i-1),2+2*(i-1)] <- Sigma[2+2*(i-1),1+2*(i-1)] <- 
        tanh(theta[i])
    }
    
  } else {
    if ((ncauses/2)>1) {
  ##    Sigma <- Sigma+tanh(theta[ncauses/2+1])
    }
    diag(Sigma) <- 1 ##diag(1,nrow=2*ncauses)
##    for (i in seq_len(ncauses/2)) {
##      Sigma[1+2*(i-1),2+2*(i-1)] <- Sigma[2+2*(i-1),1+2*(i-1)] <- 
##        tanh(theta[i])
##    }
    Sigma[1,6] <- Sigma[6,1] <- tanh(theta[1])
    Sigma[3,8] <- Sigma[8,3] <- tanh(theta[2])
    Sigma[1,8] <- Sigma[8,1] <- tanh(theta[3])
    Sigma[3,6] <- Sigma[6,3] <- tanh(theta[4])
  }
##  browser()
  attributes(Sigma)$npar <- npvar
  return(Sigma)
}




flex_cif.ppar <- function(theta,ncauses,dim=2,sym=TRUE,...) {
  if (is.null(theta)) {
    return((ncauses-1))
  }  
  P <- matrix(0,ncol=ncauses,ncauses) ## Cure rate matrix
  expit <- function(z,b=1) (b)/(1+exp(-z)) ## y\in(0,b)

  A <- 0
  pos <- 0
  for (i in seq_len(ncauses/2)) {
    pos <- pos+1
    pr0 <- expit(theta[pos],b=1-A)
    P[i,i] <- pr0    
    A <- A+pr0
  }
  P[3,3] <- P[4,4] <- expit(theta[pos+1],b=1-A)/2
  return(P)
}

flex_cif.vpar <- function(theta,ncauses,start=FALSE,...) {
  npvar <- ncauses/2+1
  if (start) return(rep(0,npvar))
  if (is.null(theta)) return(npvar)
  Sigma <- diag(1,nrow=2*ncauses)
##  browser()
  for (i in seq_len(ncauses/2)) {
    revdiag(Sigma[1:2+2*(i-1),1:2+2*(i-1)]) <- tanh(theta[i])
  }
  revdiag(Sigma[seq_len(ncauses)+ncauses,seq_len(ncauses)+ncauses]) <-
    tanh(theta[npvar])
  
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

