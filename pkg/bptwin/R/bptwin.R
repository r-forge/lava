###{{{ biprobit

biprobit <- function(formula, data, id, zyg, twinnum, weight=NULL,
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
  dS1 <- cbind(dS1,c(1,0,0,1))

  
}

###}}} biprobit

###{{{ bptwin

bptwin <- function(formula, data, id, zyg, twinnum, weight=NULL,
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
  dS1 <- cbind(dS1,c(1,0,0,1))

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

###{{{ bptwin methods

score.bptwin <- function(x,indiv=FALSE,...) {
  x$score
}
logLik.bptwin <- function(object,indiv=FALSE,...) {
  if (indiv) object$logLik else sum(object$logLik)
}
vcov.bptwin <- function(object,...) object$vcov
coef.bptwin <- function(object,...) object$coef
print.bptwin <- function(x,...) {
  printCoefmat(coef(x))
  S <- colSums(x$score); names(S) <- rep("",length(S))
  cat("\nScore: "); cat(S);
  cat("\nlogLik: "); cat(sum(x$logLik),"\n");
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
    iJ <- solve(J)
    U. <- function(p) U(p,group=i)
    I <- -numDeriv::jacobian(U.,op$par)    
    ##    I <- (I+t(I))/2
    V <- iJ%*%I%*%iJ
##    V <- iJ

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

###{{{ bpACE

bpACE <- function(formula, data, id, zyg, twinnum, weight=NULL,
                  control=list(trace=1), intercepts=ifelse(unique,2,1),
                  unique=TRUE, p,...) {

  mycall <- match.call()
  idtab <- table(data[,id])
  data0 <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  data0 <- data0[order(data0[,id]),]
  idx2 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[1]
  idx1 <- as.factor(data0[,zyg])==levels(as.factor(data0[,zyg]))[2]

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
  plen <- nx+2
  dS0 <- rbind(rep(1,4),rep(1,4))
  dS1 <- rbind(c(1,.5,.5,1),rep(1,4))
  A <- matrix(dS1[1,],ncol=2)

  mytr <- function(x) x
  dmytr <- function(x) 1
  mytr <- exp
  dmytr <- exp
#  mytr <- function(x) x^2
#  dmytr <- function(x) 2*x
  
  U <- function(p,indiv=FALSE) {
    p1 <- p
    p[c(plen-1,plen)] <- mytr(p[c(plen-1,plen)])
    B <- cbind(p[midx])
    Sigma0 <- diag(2) + p[plen-1] + p[plen]
    mu0 <- X0%*%B
    Mu0 <- matrix(mu0,ncol=2,byrow=TRUE)
    U0 <- .Call("biprobit0",
                Mu0,
                Sigma0,dS0,Y0,XX0,W0,!is.null(W0))         
    Sigma1 <- diag(2) + p[plen-1]*A + p[plen]
    mu1 <- X1%*%B
    Mu1 <- matrix(mu1,ncol=2,byrow=TRUE)      
    U1 <- .Call("biprobit0",
                  Mu1,
                Sigma1,dS1,Y1,XX1,W1,!is.null(W1))
    if (indiv) {
      U0$score[,plen] <- U0$score[,plen]*dmytr(p1[plen])
      U0$score[,plen-1] <- U0$score[,plen-1]*dmytr(p1[plen-1])
      U1$score[,plen] <- U1$score[,plen]*dmytr(p1[plen])
      U1$score[,plen-1] <- U1$score[,plen-1]*dmytr(p1[plen-1])
      val <- rbind(U0$score,U1$score)
      attributes(val)$logLik <- c(U0$loglik,U1$loglik)
      return(val)
    }
    val <- colSums(U1$score)+colSums(U0$score)
    val[plen-1] <- val[plen-1]*dmytr(p1[plen-1])
    val[plen] <- val[plen]*dmytr(p1[plen])
    attributes(val)$logLik <- sum(U1$loglik)+sum(U0$loglik)
    return(val)
  }

  p0 <- rep(1,plen)
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
  vnames1 <- NULL
  trnam <- ""
  rownames(cc) <- c(rnames1,c("log(var(A))","log(var(C))"))
  val <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall)
  class(val) <- "bptwin"  
  return(val)
}

###}}} bpACE

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

