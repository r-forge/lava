###{{{ Dbvn

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

###}}} Dbvn

###{{{ bptwin

##--------------------------------------------------
bptwin <- function(formula, data, id, zyg, DZ, DZos,
                   weight=NULL,
                   biweight=function(x) 1/min(x),
                   strata=NULL,
                   messages=1,
                   control=list(trace=0),
                   type="ace",
                   eqmean=TRUE,
                   param=0,
                   pairsonly=FALSE,
                   stderr=TRUE,                  
                   robustvar=TRUE,                   
                   p, indiv=FALSE,
                   constrain,
                   samecens=TRUE,
                   allmarg=samecens&!is.null(weight),
                   bound=FALSE,
                   debug=FALSE,...) {

###{{{ setup

  mycall <- match.call()
  formulaId <- Specials(formula,"cluster")
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- paste("~.-cluster(",formulaId,")-strata(",paste(formulaStrata,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) strata <- formulaStrata
  mycall$formula <- formula
 
  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        mycall$data <- dd[[i]]
        eval(mycall)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("biprobit.strata","biprobit")
      res$coef <- unlist(lapply(res$model,coef))
      res$vcov <- blockdiag(lapply(res$model,vcov.biprobit))
      res$N <- length(dd)
      res$idx <- seq(length(coef(res$model[[1]])))
      rownames(res$vcov) <- colnames(res$vcov) <- names(res$coef)
      return(res)
    }
  }

##################################################
### No strata
    if (is.null(control$method)) {
    control$method <- "gradient"
    if (!samecens & !is.null(weight)) control$method <- "bhhh"
  }
  if (length(grep("flex",tolower(type)))>0) { type <- "u"; eqmean <- FALSE }

  yvar <- paste(deparse(formula[[2]]),collapse="")
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  if (pairsonly)
    data <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  if (is.logical(data[,yvar])) data[,yvar] <- data[,yvar]*1
  if (is.factor(data[,yvar])) data[,yvar] <- as.numeric(data[,yvar])-1  

  ## Y <- suppressWarnings(model.frame(formula,data,na.action="na.pass")[,1])    
  ## Blen <- 0
  ## Bord <- c()
  ## if (degree>0)
  ##   if (is.Surv(Y) | !missing(time) | !is.null(B)) {
  ##     if (is.Surv(Y)) {
  ##       if (ncol(Y)==3) {
  ##         entry <- Y[,1]
  ##         time <- Y[,2]
  ##         Y <- Y[,3]*1 # event
  ##       } else {
  ##         time <- Y[,1]
  ##         Y <- Y[,2]*1
  ##       }
  ##     } else {
  ##       if (is.character(time))
  ##         time <- data[,time]
  ##       if (is.character(entry))
  ##         entry <- data[,entry]
  ##     }
  ##     if (is.null(B)) {
  ##       B <- bs(time,degree=degree)
  ##       Bord <- cbind(time,B)[order(time),,drop=FALSE]
  ##       colnames(B) <- paste("B","_",1:ncol(B),sep="")
  ##       Blen <- ncol(B)   
  ##     }
  ##   }

  idtab <- table(data[,id])
  ##  ii0 <- which(as.character(data[,id])%in%names(idtab)[idtab==2])              
  ##  data0 <- data[ii0,]  
  idx2 <- NULL
  if (!missing(DZ)) {
    if (!missing(DZos))
      idx2 <- data[,zyg]==DZos
    idx1 <- data[,zyg]==DZ
    idx0 <- data[,zyg]!=DZ
    data[,zyg] <- (data[,zyg]!=DZ)*1
  } else {
    if (!missing(DZos))
      idx2 <- (as.factor(data[,zyg])==levels(as.factor(data[,zyg]))[3]) # DZos
    DZlev <- levels(as.factor(data[,zyg]))[1]
    idx1 <- (as.factor(data[,zyg])==DZlev) # DZ
    idx0 <- (as.factor(data[,zyg])==levels(as.factor(data[,zyg]))[2]) # MZ
    message("Using '",DZlev,"' as DZ",sep="")
    data[,zyg] <- (data[,zyg]!=levels(as.factor(data[,zyg]))[1])
  }
  
  
  time <- "time"
  while (time%in%names(data)) time <- paste(time,"_",sep="")
  data[,time] <- unlist(lapply(idtab,seq))
  
  ##  ff <- as.formula(paste("~",paste(attributes(terms(formula))$term.labels,
  ##                                   collapse="+"),"+1",sep=""))
  ff <- paste(as.character(formula)[3],"+",time,"+",id,"+",zyg)
  if (!is.null(weight))
    ff <- paste(weight,"+",ff)
  ff <- paste("~",yvar,"+",ff)
  formula0 <- as.formula(ff)
  Data <- model.matrix(formula0,data,na.action=na.pass)
  rnames1 <- setdiff(colnames(Data),c(yvar,time,id,weight,zyg))
  ##  X0 <- as.matrix(Data[,rnames1])
  nx <-length(rnames1) ##ncol(Data)-4 + !is.null(weight)
  if (nx==0) stop("Zero design not allowed")
  
  ##  X <- cbind(model.matrix(ff,data0)
  ##  nx <- ncol(X)
  ##  if (nx==0) stop("Zero design not allowed")
##################################################


  bidx0 <- seq(nx)
  midx0 <- bidx0; midx1 <- midx0+nx  
  dS0 <- rbind(rep(1,4),rep(1,4),rep(1,4)) ## MZ
  dS1 <- rbind(c(1,.5,.5,1),rep(1,4),c(1,.25,.25,1)) ## DZ

  ##mytr <- function(x) x; dmytr <- function(x) 1
  ##mytr <- function(x) x^2; dmytr <- function(x) 2*x
  ##mytr <- function(z) 1/(1+exp(-z)); dmytr <- function(z) exp(-z)/(1+exp(-z))^2
  mytr <- exp; dmytr <- exp; myinvtr <- log
  trname <- "exp"; invtrname <- "log"      
  ACDU <- sapply(c("a","c","d","e","u"),function(x) length(grep(x,tolower(type)))>0)
  if (ACDU["u"]) {
    ##      datanh <- function(r) 1/(1-r^2)
    dmytr <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
    mytr <- tanh;  myinvtr <- atanh
    trname <- "tanh"; invtrname <- "atanh"    
    dS0 <- rbind(c(0,1,1,0))
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
    bidx1 <- bidx0
##    midx0 <- midx1 <- midx    
  } else {
    bidx1 <- bidx0+nx
##    midx0 <- 1:nx; midx1 <- midx0+nx
    nx <- 2*nx;
  }
  
  vidx0 <- vidx0+nx; vidx1 <- vidx1+nx
  vidx <- nx+seq_len(nvar)
  midx <- seq_len(nx)
  plen <- nx+nvar

  Am <- matrix(c(1,.5,.5,1),ncol=2)
  Dm <- matrix(c(1,.25,.25,1),ncol=2)
  Em <- diag(2)

##################################################

  Wide <- reshape(as.data.frame(Data),idvar=c(id,zyg),timevar=time,direction="wide")
  yidx <- paste(yvar,1:2,sep=".")
  rmidx <- c(id,yidx,zyg)
  ## if (!is.null(weight)) {
  ##   W <- cbind(data[,weight])
  ##   widx <- paste(weight,1:2,sep=".")
  ##   WW <- as.matrix(Wide[,widx])
  ##   
  ## }
  W0 <- W1 <- W2 <- NULL
  if (!is.null(weight)) {
    widx <- paste(weight,1:2,sep=".")
    rmidx <- c(rmidx,widx)
    W0 <- as.matrix(Wide[Wide[,zyg]==1,widx,drop=FALSE])
    W1 <- as.matrix(Wide[Wide[,zyg]==0,widx,drop=FALSE])
  }
  XX <- as.matrix(Wide[,setdiff(colnames(Wide),rmidx)])
  XX[is.na(XX)] <- 0
  Y0 <- as.matrix(Wide[Wide[,zyg]==1,yidx,drop=FALSE])
  Y1 <- as.matrix(Wide[Wide[,zyg]==0,yidx,drop=FALSE])
  XX0 <- XX[Wide[,zyg]==1,,drop=FALSE]
  XX1 <- XX[Wide[,zyg]==0,,drop=FALSE]
  
##################################################

###}}} setup

###{{{ Mean/Var function

  ##Marginals etc.
  MyData0 <- ExMarg(Y0,XX0,W0,dS0,eqmarg=TRUE,allmarg=allmarg)
  MyData1 <- ExMarg(Y1,XX1,W1,dS1,eqmarg=TRUE,allmarg=allmarg) 
  N <- cbind(sum(idx0),sum(idx1),sum(idx2)); 
  if (missing(DZos)) N <- N[,-3,drop=FALSE]
  N <- cbind(N,
             2*nrow(MyData0$Y0)+nrow(MyData0$Y0_marg),
             2*nrow(MyData1$Y0)+nrow(MyData1$Y0_marg),
             nrow(MyData0$Y0),nrow(MyData1$Y0))

  if (samecens & !is.null(weight)) {
    MyData0$W0 <- cbind(apply(MyData0$W0,1,biweight))
    if (!is.null(MyData0$Y0_marg))
      MyData0$W0_marg <- cbind(apply(MyData0$W0_marg,1,biweight))
  }
  if (samecens & !is.null(weight)) {
    MyData1$W0 <- cbind(apply(MyData1$W0,1,biweight))
    if (!is.null(MyData1$Y0_marg))
      MyData1$W0_marg <- cbind(apply(MyData1$W0_marg,1,biweight))
  }
  rm(Y0,XX0,W0,Y1,XX1,W1)
  ##  suppressMessages(browser())
  ##  N <- rbind(##c("","","Complete","","Complete pairs",""),
  ##             rep(c("MZ","DZ"),ncol(N)/2), N)
  colnames(N) <- c("Total.MZ","Total.DZ","Complete.MZ","Complete.DZ","Complete pairs.MZ","Complete pairs.DZ")[seq(ncol(N))]
  rownames(N) <- rep("",nrow(N))
##  print(N,quote=FALSE)
  
  ## Mu <- function(p0) {
  ##   b0 <- cbind(p0[midx0])
  ##   b1 <- cbind(p0[midx1])
  ##   b00 <- b0; b11 <- b1
  ##   if (Bconstrain) {
  ##     b00 <- trMean(b0,Blen); b11 <- trMean(b1,Blen)
  ##   }
  ##     mu0 <- with(MyData0, X0%*%b00)
  ##     mu1 <- with(MyData1, X0%*%b11)
  ##   return(list(mu0=mu0,mu1=mu1))
  ## }
  Sigma <- function(p0) {    
    p0[vidx] <- mytr(p0[vidx])    
    if (ACDU["u"]) {     
      ##      Sigma0 <- Em+p0[plen-1]; Sigma1 <- Em+p0[plen]
      Sigma0 <- diag(2) + p0[plen-1]*matrix(c(0,1,1,0),2,2)
      Sigma1 <- diag(2) + p0[plen]*matrix(c(0,1,1,0),2,2)
    } else {    
      pv <- ACDU*1;  pv[which(ACDU[1:3])] <- p0[vidx]
      Sigma0 <- Em*pv["e"] + pv["a"] + pv["c"] + pv["d"]
      Sigma1 <- Em*pv["e"] + pv["a"]*Am + pv["c"] + pv["d"]*Dm
    }
    return(list(Sigma0=Sigma0,Sigma1=Sigma1))
  }

  ###}}} Mean/Var function
  
###{{{ U  

  U <- function(p,indiv=FALSE) {
    b0 <- cbind(p[bidx0])
    b1 <- cbind(p[bidx1])
    b00 <- b0; b11 <- b1
    ## if (Bconstrain) {
    ##   b00 <- trMean(b0,Blen); b11 <- trMean(b1,Blen)
    ## }
    if (bound) p[vidx] <- min(p[vidx],20)
    S <- Sigma(p)
    lambda <- eigen(S$Sigma0)$values
    if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
    ##    browser()
    ##mu0 <- with(MyData0, X0%*%b00)
    ##    Mu0 <- matrix(mu0,ncol=2,byrow=TRUE)
    ##        browser()
    
    Mu0 <- with(MyData0, cbind(XX0[,midx0,drop=FALSE]%*%b00,
                               XX0[,midx1,drop=FALSE]%*%b00))
    U0 <- with(MyData0, .Call("biprobit0",
                             Mu0,
                             S$Sigma0,dS0,Y0,XX0,W0,!is.null(W0),samecens))

    if (!is.null(MyData0$Y0_marg)) {
      mum <- with(MyData0, XX0_marg%*%b00)
       U_marg <- with(MyData0, .Call("uniprobit",
                                   mum,XX0_marg,
                                   S$Sigma0[1,1],t(dS0_marg),Y0_marg,
                                   W0_marg,!is.null(W0_marg),TRUE))
      U0$score <- rbind(U0$score,U_marg$score)
      U0$loglik <- c(U0$loglik,U_marg$loglik)
    }

##    mu1 <- with(MyData1, X0%*%b11) 
##    Mu1 <- matrix(mu1,ncol=2,byrow=TRUE)
    Mu1 <- with(MyData1, cbind(XX0[,midx0,drop=FALSE]%*%b11,
                               XX0[,midx1,drop=FALSE]%*%b11))

    U1 <- with(MyData1, .Call("biprobit0",
                             Mu1,
                             S$Sigma1,dS1,Y0,XX0,W0,!is.null(W0),samecens))
    if (!is.null(MyData1$Y0_marg)) {
      mum <- with(MyData1, XX0_marg%*%b11)
      U_marg <- with(MyData1, .Call("uniprobit",
                                    mum,XX0_marg,
                                    S$Sigma1[1,1],t(dS0_marg),Y0_marg,
                                    W0_marg,!is.null(W0_marg),TRUE))
      U1$score <- rbind(U1$score,U_marg$score)
      U1$loglik <- c(U1$loglik,U_marg$loglik)
    }

    if (indiv) {

      val0 <- U0$score[MyData0$id,,drop=FALSE]
      val1 <- U1$score[MyData1$id,,drop=FALSE]
      N0 <- length(MyData0$id)
      idxs0 <- seq_len(N0)
      for (i in seq_len(N0)) {
        idx0 <- which((MyData0$idmarg)==(MyData0$id[i]))+N0
        idxs0 <- c(idxs0,idx0)
        val0[i,] <- val0[i,]+colSums(U0$score[idx0,,drop=FALSE])
      }
      val0 <- rbind(val0, U0$score[-idxs0,,drop=FALSE])
      N1 <- length(MyData1$id)
      idxs1 <- seq_len(N1)
      for (i in seq_len(N1)) {
        idx1 <- which((MyData1$idmarg)==(MyData1$id[i]))+N1
        idxs1 <- c(idxs1,idx1)
        val1[i,] <- val1[i,]+colSums(U1$score[idx1,,drop=FALSE])
      }
      val1 <- rbind(val1, U1$score[-idxs1,,drop=FALSE])

      val <- matrix(0,ncol=plen,nrow=nrow(val0)+nrow(val1))
      val[seq_len(nrow(val0)),c(bidx0,vidx0)] <- val0
      val[nrow(val0)+seq_len(nrow(val1)),c(bidx1,vidx1)] <- val1
      for (ii in vidx) {
        val[,ii] <- val[,ii]*dmytr(p[ii])
      }
      attributes(val)$logLik <- c(U0$loglik,U1$loglik)
      return(val)
      
#########      
      ## val <- matrix(0,ncol=plen,nrow=nrow(U0$score)+nrow(U1$score))
      ## val[seq_len(nrow(U0$score)),c(bidx0,vidx0)] <- U0$score
      ## val[nrow(U0$score)+seq_len(nrow(U1$score)),c(bidx1,vidx1)] <- U1$score
      ## for (ii in vidx) {
      ##   val[,ii] <- val[,ii]*dmytr(p[ii])
      ## }
      ## ## if (Bconstrain & Blen>0) {
      ## ##   Bidx <- attributes(b00)$idx
      ## ##   val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(b00)$D[Bidx,Bidx,drop=FALSE])
      ## ##   if (!eqmean) {
      ## ##     Bidx <- bidx1[attributes(b11)$idx]
      ## ##     val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(b11)$D[attributes(b11)$idx,attributes(b11)$idx])
      ## ##   }
      ## ## }      
      ## attributes(val)$logLik <- c(U0$loglik,U1$loglik)
      ## return(val)
      
    }
    val <- numeric(plen)
    val[c(bidx0,vidx0)] <- colSums(U0$score)
    val[c(bidx1,vidx1)] <- val[c(bidx1,vidx1)]+colSums(U1$score)
    for (ii in vidx)
      val[ii] <- val[ii]*dmytr(p[ii])
    ## if (Bconstrain & Blen>0) {
    ##   Bidx <- attributes(b00)$idx
    ##   val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b00)$D[Bidx,Bidx])
    ##   if (!eqmean) {
    ##     Bidx <- bidx1[attributes(b11)$idx]
    ##     val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b11)$D[attributes(b11)$idx,attributes(b11)$idx])
    ##   }
    ## }    
    attributes(val)$logLik <- sum(U0$loglik)+sum(U1$loglik)
    return(val)
  }

###}}} U

###{{{ optim

  p0 <- rep(-1,plen); ##p0[vidx] <- 0
  if (type=="u")
    p0[vidx] <- 0.3
  if (!is.null(control$start)) {
    p0 <- control$start
    control$start <- NULL
  } else {
    X <- rbind(MyData0$XX0[,midx0,drop=FALSE],MyData0$XX0[,midx1,drop=FALSE])
    Y <- rbind(MyData0$Y0[,1,drop=FALSE],MyData0$Y0[,2,drop=FALSE])
    g <- suppressWarnings(glm(Y~-1+X,family=binomial(probit)))
    p0[midx] <- coef(g)
    ## if (Blen>0) {
    ##   pB <- p0[tail(midx,Blen)]
    ##   pB[1] <- ifelse(pB[1]<0,-2,log(pB[1]))
    ##   if (Blen>1) {
    ##     pB[seq_len(Blen-1)+1] <- -2
    ##   }
    ##   p0[tail(midx,Blen)] <- pB
    ## }
  }
 
  if (!missing(p)) return(U(p,indiv=indiv))


  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))
  h0 <- function(p) crossprod(U(p,indiv=TRUE))

  
  if (!missing(constrain)) {
    freeidx <- is.na(constrain)
    f <- function(p) {      
      p1 <- constrain; p1[freeidx] <- p
      res <- U(p1)[freeidx]
      crossprod(res)[1]
    }
    f0 <- function(p) {
      p1 <- constrain; p1[freeidx] <- p
      -sum(attributes(U(p1))$logLik)
    }
    g0 <- function(p) {
      p1 <- constrain; p1[freeidx] <- p
      -as.numeric(U(p1)[freeidx])
    }
    p0 <- p0[is.na(constrain)]    
  }


  controlstd <- list(hessian=0)
  controlstd[names(control)] <- control
  control <- controlstd
  
  nlminbopt <- intersect(names(control),c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min"))
  ucminfopt <- intersect(names(control),c("trace","grtol","xtol","stepmax","maxeval","grad","gradstep","invhessian.lt"))
  optimopt <- names(control) 

  if (debug) browser()
  op <- switch(tolower(control$method),
               nlminb=nlminb(p0,f0,grad=g0,control=control[nlminbopt],...),
               optim=optim(p0,fn=f0,gr=g0,control=control[ucminfopt],...),
               ucminf=,
               quasi=,
               gradient=ucminf(p0,fn=f0,gr=g0,control=control[ucminfopt],hessian=0,...),
               bhhh={
                 controlnr <- list(stabil=FALSE,
                                   gamma=0.1,
                                   gamma2=1,
                                   ngamma=5,
                                   iter.max=200,
                                   epsilon=1e-12,
                                   tol=1e-9,
                                   trace=1,
                                   stabil=FALSE)
                 controlnr[names(control)] <- control
                 lava:::NR(start=p0,NULL,g0, h0,control=controlnr)
               },
               ##                 op <- switch(mycontrol$method,
               ##                              ucminf=ucminf(p0,f,control=mycontrol[ucminfopt],hessian=F,...),
               ##                optim=optim(p0,f,control=mycontrol[ucminfopt],...),
                 nlminb(p0,f,control=control[nlminbopt],...))
  ##  op <- nlm(f,p0,print.level=2)
  ##  op <- spg(p0,f,control=control,...)
  

  if (stderr) {
    
    ## WW <- rbind(W0,W1)
    ## ff <- function(p) {
    ##   UU <- U(p,indiv=TRUE)
    ##   s <- apply(UU,2,function(x) x) ##*WW[,1,drop=FALSE])
    ##   l <- attributes(UU)$logLik##*WW[,1]
    ##   res <- structure(sum(l), grad=colSums(s))
    ##   res
    ## }
    ## browser()
    
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
  } else {
    UU <- matrix(NA,ncol=length(op$par),nrow=1)
    I <- J <- V <- matrix(NA,ncol=length(op$par),nrow=length(op$par))
  }

  ###}}} optim

###{{{ return

  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  vnames1 <- NULL
  trnam <- " "
  if (debug) browser()
  if (!eqmean) {
    rnames1 <- c(paste(rnames1,"MZ",sep=trnam),paste(rnames1,"DZ",sep=trnam))
  }

  if (ACDU["u"]) {
##    rnames <- c(rnames1,paste(c("log(var(U))","log(var(U))"),c("MZ","DZ"),sep=trnam))
    rnames <- c(rnames1,paste(c("atanh(rho)","atanh(rho)"),c("MZ","DZ"),sep=trnam))
  } else {
    rnames <- c(rnames1,c("log(var(A))","log(var(C))","log(var(D))")[ACDU[1:3]])
  }
  if (!missing(constrain)) rnames <- rnames[freeidx]
  rownames(cc) <- rnames
  rownames(V) <- colnames(V) <- rnames
  S <- Sigma(op$par)

  npar <- list(intercept=attributes(terms(formula))$intercept,
               pred=nrow(attributes(terms(formula))$factor)-1,
               var=sum(ACDU[-4]),
               ACDU=ACDU[-4]*1)
  
  npar[unlist(lapply(npar,length))==0] <- 0
##  npar$var <- nrow(cc)-sum(unlist(npar))
  
  val <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, Sigma0=S$Sigma0, Sigma1=S$Sigma1, dS0=dS0, dS1=dS1, N=N, midx0=midx0, midx1=midx1, vidx0=vidx0, vidx1=vidx1, eqmean=eqmean, I=I,J=J, robustvar=robustvar,
              transform=list(tr=mytr, invtr=myinvtr, dtr=dmytr,
                name=trname, invname=invtrname),
              SigmaFun=Sigma, ##MuFun=Mu,
              npar=npar
              )
  class(val) <- c("bptwin","biprobit")
  return(val)
}

###}}} return


##--------------------------------------------------

###}}} bptwin

###{{{ summary.bptwin

summary.bptwin <- function(object,level=0.05,...) {
  logit <- function(p) log(p/(1-p))
  tigol <- function(z) 1/(1+exp(-z))
  dlogit <- function(p) 1/(p*(1-p))
  trnam <- " "
  vcoef1 <- paste("log(var(",c("A","C","D"),"))",sep="")
  vcoef2 <- paste("atanh(",
                  c(paste("rho)","MZ",sep=trnam),
                    paste("rho)","DZ",sep=trnam)),sep="")
  idx1 <- na.omit(match(vcoef1,names(coef(object))))
  idx2 <- na.omit(match(vcoef2,names(coef(object))))
  CIs <- c()
  alpha <- level/2
  CIlab <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  V <- c()
  if (length(idx2)>0) {    
    idx <- idx2
    V <- vcov(object)[idx,idx]
    arho <- coef(object)[idx2[1:2]]
    mz <- multinomlogit(coef(object)[idx2[1]]); names(mz) <- c("U","E")
    dz <- multinomlogit(coef(object)[idx2[2]]); names(dz) <- c("U","E")
    cc <- tanh(arho)
  ##    cc <- c(mz[1],dz[1]) ##,mz[2],dz[2])
    names(cc) <- c("Correlation MZ","Correlation DZ")
##    corMZ <- mz[1]; corDZ <- dz[1]
    corMZ <- cc[1]; corDZ <- cc[2]
##    D <- (cbind(c(attributes(mz)$gradient[1],0),c(0,attributes(dz)$gradient[1])))
    D <- diag(object$tr$dtr(arho))
    h <- function(x) 2*(x[1]-x[2])
    dh <- function(x) c(2,-2)
    i1 <- 1:2
    corr <- NULL
  }
  if (length(idx1)>0) {
    idx <- idx1
    V <- vcov(object)[idx,idx]
    ACD <- match(names(coef(object))[idx1],vcoef1)
    nn <- c(c("A","C","D")[ACD],"E")
    dzsc <- c(1/2,1,1/4)[ACD]
    pp <- coef(object)[idx1]
    cc <- multinomlogit(pp); names(cc) <- nn
    D <- attributes(cc)$gradient  
    ##    p <- coef(object)[2:3,1]
    ##    F <- function(p) {
    ##      logit(multinomlogit(p))
    ##    }    
    cc2 <- logit(cc)
    D2 <- diag(dlogit(cc))
    DD <- D2%*%D
    Vc2 <- DD%*%V%*%t(DD)
    CIs <- tigol(cc2%x%cbind(1,1)+diag(Vc2)^0.5%x%cbind(-1,1)*qnorm(1-alpha))
    K <- length(ACD)
    Ki <- seq_len(K)
    corMZ <- sum(cc[Ki]); corDZ <- sum(cc[Ki]*dzsc)
    i1 <- seq_len(length(dzsc))
    h <- function(x) 2*(sum(x[i1])-sum(x[i1]*dzsc))
    dh <- function(x) 2*(1-dzsc)
    
  }
  Vc <- D%*%V%*%t(D)
  datanh <- function(r) 1/(1-r^2)
  if (length(idx1)>0) {
    pp <- coef(object)[idx]
    b <- cbind(rep(1,K))
    corMZ.sd <- (t(b)%*%Vc[Ki,Ki]%*%b)[1]^0.5
    corDZ.sd <- (t(dzsc)%*%Vc[Ki,Ki]%*%dzsc)[1]^0.5    
    corr <- rbind(c(corMZ,corMZ.sd),c(corDZ,corDZ.sd))
    zrho <- atanh(corr[,1])
    zrho.var <- datanh(corr[,1])^2*corr[,2]^2
    corr <- cbind(corr, tanh(zrho%x%cbind(1,1)+zrho.var^0.5%x%cbind(-1,1)*qnorm(1-alpha)))
    rownames(corr) <- c("Correlation MZ","Correlation DZ")    
  } else {   
    zrho <- atanh(cc)
    zrho.var <- datanh(cc)^2*diag(Vc)
    CIs <- tanh(zrho%x%cbind(1,1)+zrho.var^0.5%x%cbind(-1,1)*qnorm(1-alpha))
  }
  newcoef <- rbind(cbind(cc,diag(Vc)^0.5,CIs),corr);
  ##  CIs <- rbind(CIs,c(NA,NA),c(NA,NA))
  ##  newcoef <- cbind(newcoef,CIs)
  colnames(newcoef) <- c("Estimate","Std.Err",CIlab)
  logith <- function(x) logit(h(x))
  dlogith <- function(x) dlogit(h(x))*dh(x)
  Dlh <- dlogith(cc[i1])
  sdlh <- (t(Dlh)%*%Vc[i1,i1]%*%(Dlh))[1]^0.5
  H <- h(cc[i1])
  hstd <- t(dh(cc[i1]))%*%Vc[i1,i1]%*%dh(cc[i1])
  ci <- tigol(logith(cc[i1]) + qnorm(1-alpha)*c(-1,1)*sdlh)  
  
  concordance <-  conditional <- marg <- c()

  probs <- function(p,idx=1) {
    S <- (object$SigmaFun(p))[[idx]]
    m <- 0
    if((object$npar$intercept==1 & idx==1) | object$eqmean) m <- p[1]
    else m <- p[length(object$midx0)+1]
    mu.cond <- function(x) m+S[1,2]/S[2,2]*(x-m)
    var.cond <- S[1,1]-S[1,2]^2/S[2,2]    
    conc <- pmvnorm(upper=c(m,m),sigma=S,algorithm="Miwa")
    marg <- pnorm(m,sd=S[1,1]^0.5)
    cond <- conc/marg
    logit(c(conc,cond,marg))
  }

  mycoef <- coef(object)  
  formals(probs) <- alist(p=,idx=1)
  probMZ <- probs(mycoef)

  Dp0 <- jacobian(probs,mycoef)
  formals(probs) <- alist(p=,idx=2)
  probDZ <- probs(mycoef)
  Dp1 <- jacobian(probs,mycoef)
  sprobMZ <- diag((Dp0)%*%vcov(object)%*%t(Dp0))^0.5
  sprobDZ <- diag((Dp1)%*%vcov(object)%*%t(Dp1))^0.5
  probMZ <- tigol(cbind(probMZ,probMZ-qnorm(1-alpha)*sprobMZ,probMZ+qnorm(1-alpha)*sprobMZ))
  probDZ <- tigol(cbind(probDZ,probDZ-qnorm(1-alpha)*sprobDZ,probDZ+qnorm(1-alpha)*sprobDZ))
  rownames(probMZ) <- rownames(probDZ) <- c("Concordance","Conditional","Marginal")
  colnames(probMZ) <- colnames(probDZ) <- c("Estimate",CIlab)
 
  ## mu <- coef(object)[c(object$bidx0[1],object$bidx1[1])]
  ## Sigma <- list(object$Sigma0,object$Sigma1)
  ## for (i in 1:2) {
  ##   conc <- function()
  ##   mu.cond <- function(x) mu+Sigma[[i]][1,2]/Sigma[[i]][2,2]*(x-mu[i])
  ##   var.cond <- Sigma[[i]][1,1]-Sigma[[i]][1,2]^2/Sigma[[i]][2,2]    
  ##   cc0 <- pmvnorm(upper=c(mu[i],mu[i]),sigma=Sigma[[i]])
  ##   px <- pnorm(mu[i],sd=Sigma[[i]][2,2]^0.5)
  ##   concordance <- c(concordance,cc0)
  ##   marg <- c(marg,px)
  ##   conditional <- c(conditional,cc0/px)
  ## }
  ## names(concordance) <- names(conditional) <- c("MZ","DZ")

  hval <- rbind(c(H,hstd^0.5,ci)); colnames(hval) <- c("Estimate","Std.Err",CIlab); rownames(hval) <- "Heritability"

  Nstr <- object$N
  nN <- ncol(object$N)
  npos <- seq(nN/2)
  Nstr <- rbind(paste(Nstr[npos*2-1],Nstr[npos*2],sep="/"))
  rownames(Nstr) <- ""
  colnames(Nstr) <- unlist(lapply(strsplit(colnames(object$N)[npos*2-1],".",fixed=TRUE),
                                  function(x) paste(x[1], "MZ/DZ")))
  res <- list(object=object, h=hval,
              probMZ=probMZ, probDZ=probDZ, Nstr=Nstr,
              coef=newcoef) ##, concordance=concordance, conditional=conditional)
  class(res) <- "summary.bptwin"
  res
}

print.summary.bptwin <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  print(x$object,digits=digits,...)
  cat("\n")
  x$Nstr <- x$Nstr[,which((colnames(x$Nstr)!="Complete MZ/DZ")),drop=FALSE]
  print(x$Nstr,quote=FALSE)
  cat("\n")
  print(RoundMat(x$coef[,-2,drop=FALSE],digits=digits),quote=FALSE)
  cat("\nMZ:\n");
  print(RoundMat(x$probMZ,digits=digits),quote=FALSE)
  cat("DZ:\n")
  print(RoundMat(x$probDZ,digits=digits),quote=FALSE)
##  cat("\nConcordance (MZ; DZ):\t\t", x$concordance,"\n")
##  cat("Case-wise concordance (MZ; DZ):\t", x$conditional,"\n\n")
  cat("\n")
  print(RoundMat(x$h[,-2,drop=FALSE],digits=digits),quote=FALSE)
  cat("\n")
}

###}}} summary.bptwin

###{{{ bptwin/biprobit methods

do.biprobit.strata <- function(x,fun,print=FALSE,...) {
  res <- c()
  for (i in seq(length(x$model))) {    
    message(rep("-",60),sep="")
    message("Strata '",names(x$model)[i],"'",sep="")
    myargs <- c(list(x$model[[i]]),list(...))
    s <- do.call(fun, myargs)
    res <- c(res,list(s))
    if (print) print(s)
  }
  invisible(res)
}

plot.biprobit.strata <- function(x,...)
  suppressMessages(do.biprobit.strata(x,"plot",...))
print.biprobit.strata <- function(x,...)
  do.biprobit.strata(x,"print",...)
summary.biprobit.strata <- function(object,...)
  do.biprobit.strata(object,"summary",print=TRUE,...)
coef.biprobit.strata <- function(object,...) object$coef
logLik.biprobit.strata <- function(object,indiv=FALSE,list=FALSE,...) {
  ll <- lapply(object$model,function(x) logLik(x,indiv=indiv,...))
  if (list) return(ll)
  if (!indiv) {
    res <- structure(sum(unlist(ll)),df=0,nall=0)
    for (i in seq(length(ll))) {
      attributes(res)$nall <- attributes(res)$nall+attributes(ll[[i]])$nall
      attributes(res)$df <- attributes(res)$df+attributes(ll[[i]])$df
    }
    attributes(res)$nobs <- attributes(res)$nall-attributes(res)$df
    class(res) <- "logLik"
    return(res)
  }
  return(unlist(ll))
}
score.biprobit.strata <- function(object,...) {
  ss <- lapply(object$model,function(x) score(x,indiv=FALSE,...))
  return(unlist(ss))
}


print.biprobit <- function(x,...) {
  printCoefmat(x$coef,...)
  invisible(x)
}

print.summary.biprobit <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  print(x$object,digits=digits)
  S <- colSums(x$object$score);  names(S) <- rep("",length(S))
  cat("\n")
  print(x$object$N,quote=FALSE)
  ##  suppressMessages(browser())
  cat("Score: "); cat(formatC(S,...));
  cat("\nlogLik: "); cat(sum(x$object$logLik),"\n");
  if (!is.null(x$object$msg)) {
    cat(x$object$msg,"\n")
  }

  if (!is.null(x$varcomp)) {
    cat("\n")
    res <- x$varcomp
    if (!is.null(x$prob)) {
      res <- rbind(res,x$prob)
    }    
    print(RoundMat(res,digits=digits),quote=FALSE)
  }
  cat("\n")
}

summary.biprobit <- function(object,level=0.05,...) {
  alpha <- level/2
  varcomp <- object$coef[length(coef(object)),1:2]
  varcomp <- rbind(object$model$tr(c(varcomp[1],varcomp[1]%x%cbind(1,1) + qnorm(1-alpha)*varcomp[2]%x%cbind(-1,1))))
  colnames(varcomp)[2:3] <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  rownames(varcomp) <- ifelse(is.null(object$model$varcompname),"Variance component",object$model$varcompname)

  logit <- function(p) log(p/(1-p))
  tigol <- function(z) 1/(1+exp(-z))
  dlogit <- function(p) 1/(p*(1-p))
  probs <- function(p) {
    ##    S <- diag(2); S[1,2] <- S[2,1] <- exp(tail(p,1))
    S <- object$SigmaFun(p)
    m <- c(0,0)
    if (object$npar$intercept==1) m[1:2] <- p[1]
    if (object$npar$intercept==2) {
      m[1:2] <- p[c(1,1+object$npar$pred/2+1)]
    }
    mu.cond <- function(x) m[1]+S[1,2]/S[2,2]*(x-m[2])
    var.cond <- S[1,1]-S[1,2]^2/S[2,2]    
    conc <- pmvnorm(upper=m,sigma=S)
    marg <- pnorm(m[1],sd=S[1,1]^0.5)
    cond <- conc/marg
    logit(c(conc,cond,marg))
  }
  alpha <- level/2
  CIlab <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  mycoef <- coef(object)
  prob <- probs(mycoef)
  Dprob <- jacobian(probs,mycoef)
  sprob <- diag((Dprob)%*%vcov(object)%*%t(Dprob))^0.5
  pp <- tigol(cbind(prob,prob-qnorm(1-alpha)*sprob,prob+qnorm(1-alpha)*sprob))
  rownames(pp) <- c("Concordance","Case-wise/Conditional","Marginal")
  colnames(pp) <- c("Estimate",CIlab)
  
  res <- list(object=object,varcomp=varcomp,prob=pp)
  class(res) <- "summary.biprobit"
  res
}

plot.bptwin <- function(x,n=50,rg=range(x$B[,1]),xlab="Time",ylab="Concordance",...) {
  require(mvtnorm)
  if (x$Blen>0) {
    ##    rg <- range(x$B[,1])
    t <- seq(rg[1],rg[2],length.out=n)
    B0 <- bs(t,degree=x$Blen)
    b0. <- coef(x)[x$midx0]
    b1. <- coef(x)[x$midx1]
    b0 <- trMean(b0.,x$Blen)
    b1 <- trMean(b1.,x$Blen)
    b00 <- tail(b0,x$Blen)
    b11 <- tail(b1,x$Blen)
    pr0 <- sapply(as.numeric(B0%*%b00+b0[1]), function(z)
                  pmvnorm(upper=rep(z,2),sigma=x$Sigma0))
    pr1 <- sapply(as.numeric(B0%*%b11+b1[1]), function(z)
                  pmvnorm(upper=rep(z,2),sigma=x$Sigma1))
    plot(pr0~t,type="l", xlab=xlab, ylab=ylab,...)
    lines(pr1~t,type="l",lty=2)
  }
  invisible(x)
}

sim.bptwin <- function(x,n=100,p,...) {
  return(x)
}
score.biprobit <- function(x,indiv=FALSE,...) {
  if (indiv) { s <- x$score; attributes(s)$logLik <- NULL; return(s) }
  colSums(x$score)
}
logLik.biprobit <- function(object,indiv=FALSE,...) {
  if (indiv) return(object$logLik)
  n <- sum(object$N[1])
  p <- length(coef(object))
  loglik <- sum(object$logLik)
  attr(loglik, "nall") <- n
  attr(loglik, "nobs") <- n-p
  attr(loglik, "df") <- p
  class(loglik) <- "logLik"        
  return(loglik)
}

vcov.biprobit <- function(object,...) object$vcov
coef.biprobit <- function(object,matrix=FALSE,...) {
  if (matrix) return(object$coef)
  return(object$coef[,1])
}

###}}} bptwin/biprobit methods

###{{{ biprobit

uniprobit <- function(mu,dmu,S,dS,y,w=NULL,indiv=FALSE,...) {
  sigma <- S^0.5
  alpha <- alpha0 <- pnorm(mu,sd=sigma)
  alpha[y==0] <- 1-alpha[y==0]
  M <- -sigma*dnorm(mu/sigma)
  V <- S*alpha0+mu*M  
  U1 <- 0.5*t(dS)%*%(-alpha0+V/S)/S
  U <- matrix(0,ncol(dmu)+nrow(U1),ncol(U1))
  if (is.null(w)) w <- rep(1,ncol(U))
  for (i in seq(ncol(U))) {
    U[,i] <- c(- dmu[i,,drop=FALSE]*(1/S*M[i]),U1[,i])/alpha[i]*
       ifelse(y[i]==0,-1,1)*w[i]
  }  
  if (indiv) return(structure(t(U),logLik=log(alpha)))
  return(structure(rowSums(U),logLik=sum(log(alpha))))
}

biprobit <- function(formula, data, id, time, strata=NULL, eqmarg=TRUE,
                     indep=FALSE, weight=NULL,
                     biweight=function(x) 1/min(x),
                     samecens=TRUE, randomeffect=FALSE, vcov="robust",
                     pairsonly=FALSE,
                     allmarg=samecens&!is.null(weight),
                     control=list(trace=0,method="qausi"),
                     bound=FALSE,
                     messages=1,
                     p,...) {

  mycall <- match.call()
  formulaId <- Specials(formula,"cluster")
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- paste("~.-cluster(",formulaId,")-strata(",paste(formulaStrata,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) strata <- formulaStrata
  mycall$formula <- formula

  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    fit <- lapply(seq(length(dd)),function(i) {
      if (messages>0) message("Strata '",names(dd)[i],"'")
      mycall$data <- dd[[i]]
      eval(mycall)
    })
    res <- list(model=fit)
    res$strata <- names(res$model) <- names(dd)
    class(res) <- c("biprobit.strata","biprobit")
    res$coef <- unlist(lapply(res$model,coef))
    res$vcov <- blockdiag(lapply(res$model,vcov.biprobit))
    res$N <- length(dd)
    res$idx <- seq(length(coef(res$model[[1]])))
    rownames(res$vcov) <- colnames(res$vcov) <- names(res$coef)
    return(res)
  }
  
  if (missing(id)) {    
    if (!is.null(weight)) {
      weights <- data[,weight]
      return(glm(formula,data=data,family=binomial(probit),weights=weights,...))
    }
    return(glm(formula,data=data,family=binomial(probit),...))    
  }

  mycall <- match.call()
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  if (pairsonly) {
    data <- data[which(as.character(data[,id])%in%names(idtab)[idtab==2]),]
    idtab <- table(data[,id])
  }
  
  if (missing(time)) {
    time <- "time"
    while (time%in%names(data)) time <- paste(time,"_",sep="")
    data[,time] <- unlist(lapply(idtab,seq))
  }
  ff <- paste(as.character(formula)[3],"+",time,"+",id)
  yvar <- paste(deparse(formula[[2]]),collapse="")
  if (!is.null(weight))
    ff <- paste(weight,"+",ff)
  ff <- paste("~",yvar,"+",ff)

  if (is.logical(data[,yvar])) data[,yvar] <- data[,yvar]*1
  if (is.factor(data[,yvar])) data[,yvar] <- as.numeric(data[,yvar])-1  
##  Y <- cbind(as.numeric(data[,yvar]))-(!is.numeric(data[,yvar]))
  
  formula0 <- as.formula(ff)
  Data <- model.matrix(formula0,data,na.action=na.pass)
  rnames1 <- setdiff(colnames(Data),c(yvar,time,id,weight))
  X0 <- as.matrix(Data[,rnames1])
  
  nx <- length(rnames1)
  if (nx==0) stop("Zero design not allowed")
  
  ##  data0 <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  ##  data0 <- data
  ##  data0 <- data0[order(data0[,id]),]
  midx1 <- seq(nx)
  midx2 <- midx1+nx
  midx <- seq(2*nx)
  plen <- ifelse(eqmarg,nx+1,2*nx+1)

  Wide <- reshape(as.data.frame(Data),idvar=id,timevar=time,direction="wide")
  W0 <- NULL
  yidx <- paste(yvar,1:2,sep=".")
  rmidx <- c(id,yidx)
  if (!is.null(weight)) {
    W <- cbind(data[,weight])
    widx <- paste(weight,1:2,sep=".")
    W0 <- as.matrix(Wide[,widx])
    rmidx <- c(rmidx,widx)
  }
  Y0 <- as.matrix(Wide[,yidx])
  XX0 <- as.matrix(Wide[,setdiff(colnames(Wide),rmidx)])
  XX0[is.na(XX0)] <- 0

  datanh <- function(r) 1/(1-r^2)
  dtanh <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
  vartr <- tanh
  dvartr <- dtanh; varitr <- atanh
  trname <- "tanh"; itrname <- "atanh"    
  Sigma1 <- diag(2)  
  Sigma2 <- matrix(c(0,1,1,0),2,2)
  dS0 <- rbind(c(0,1,1,0))
  varcompname <- "Tetrachoric correlation"
  msg <- "Variance of latent residual term = 1 (standard probit link)"
  if (randomeffect) {
    dS0 <- rbind(rep(1,4))
    vartr <- dvartr <- exp; inv <- log
    trname <- "exp"; itrname <- "log"
    Sigma2 <- 1
    varcompname <- NULL
  }
  model <- list(tr=vartr,name=trname,inv=itrname,invname=itrname,deriv=dvartr,varcompname=varcompname,dS=dS0,eqmarg=eqmarg)

  MyData <- ExMarg(Y0,XX0,W0,dS0,midx1,midx2,eqmarg=eqmarg,allmarg=allmarg)
  if (samecens & !is.null(weight)) {
    MyData$W0 <- cbind(apply(MyData$W0,1,biweight))
    if (!is.null(MyData$Y0_marg)) {
      MyData$W0_marg <- cbind(apply(MyData$W0_marg,1,biweight))
    }
  }
  
  SigmaFun <- function(p,...) {
    Sigma <- Sigma1+Sigma2*vartr(p[plen])
    if (indep) Sigma <- diag(2)
    return(Sigma)
  }

  U <- function(p,indiv=FALSE) {
    if (bound) p[plen] <- min(p[plen],20)
    Sigma <- SigmaFun(p)
    lambda <- eigen(Sigma)$values
    if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
    Mu_marg <- NULL
    if (eqmarg) {
      B <- cbind(p[midx1])
      Mu <- with(MyData,
                 cbind(XX0[,midx1,drop=FALSE]%*%B,XX0[,midx2,drop=FALSE]%*%B))     
##      Mu <- with(MyData, matrix(X0%*%B,ncol=2,byrow=TRUE))
      if (!is.null(MyData$Y0_marg)) 
        Mu_marg <- with(MyData, XX0_marg%*%B)
    } else {
      B1 <- cbind(p[midx1])
      B2 <- cbind(p[midx2])
      Mu <- with(MyData,
                 cbind(XX0[,midx1,drop=FALSE]%*%B1,XX0[,midx2,drop=FALSE]%*%B2))
      if (!is.null(MyData$Y0_marg))
        Mu_marg <- with(MyData, rbind(X0_marg1%*%B1,X0_marg2%*%B2))
    }

    U <- with(MyData, .Call("biprobit2",
                             Mu,XX0,
                             Sigma,dS0*dvartr(p[plen]),Y0,W0,
                             !is.null(W0),TRUE,eqmarg))
    
    if (!is.null(MyData$Y0_marg)) {
      ## U_marg <- uniprobit(Mu_marg[,1],XX0_marg,
      ##                     Sigma[1,1],dS0_marg*dvartr(p[plen]),Y0_marg,
      ##                     W0_marg,indiv=TRUE)
      ## U$score <- rbind(U$score,U_marg)
      ## U$loglik <- c(U$loglik,attributes(U_marg)$logLik)
      ##      W0_marg <- rep(1,nrow(XX0_marg))
      ##      browser()
      U_marg <- with(MyData, .Call("uniprobit",
                                   Mu_marg,XX0_marg,
                                   Sigma[1,1],dS0_marg*dvartr(p[plen]),Y0_marg,
                                   W0_marg,!is.null(W0_marg),TRUE))
      U$score <- rbind(U$score,U_marg$score)
      U$loglik <- c(U$loglik,U_marg$loglik)
    }

    if (indiv) {
      val <- U$score[MyData$id,,drop=FALSE]
      N <- length(MyData$id)
      idxs <- seq_len(N)
      for (i in seq_len(N)) {
        idx <- which((MyData$idmarg)==(MyData$id[i]))+N
        idxs <- c(idxs,idx)
        val[i,] <- val[i,]+colSums(U$score[idx,,drop=FALSE])
      }
      val <- rbind(val, U$score[-idxs,,drop=FALSE])
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
##  browser()

  if (!missing(p)) return(U(p,indiv=FALSE))

  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))
  h0 <- function(p) crossprod(U(p,indiv=TRUE))

  if (is.null(control$method)) {
    ##    control$method <- ifelse(samecens & !is.null(weight), "bhhh","quasi")
    control$method <- "quasi"
  }
  control$method <- tolower(control$method)
  if (control$method=="score") {
    control$method <- NULL
    op <- nlminb(p0,f,control=control,...)
  } else if (control$method=="quasi") {
    control$method <- NULL
    op <- nlminb(p0,f0,gradient=g0,control=control,...)
  } else if (control$method=="bhhh") {
    controlnr <- list(stabil=FALSE,
                      gamma=0.1,
                      gamma2=1,
                      ngamma=5,
                      iter.max=200,
                      epsilon=1e-12,
                      tol=1e-9,
                      trace=1,
                      stabil=FALSE)
    controlnr[names(control)] <- control
    op <- lava:::NR(start=p0,NULL,g0, h0,control=controlnr)
  } else {
    control$method <- NULL
    op <- nlminb(p0,f0,control=control,...)
  }
  UU <- U(op$par,indiv=TRUE)
  J <- crossprod(UU)
  ##  iJ <- Inverse(J)
  iI <- Inverse(-numDeriv::jacobian(U,op$par))
  V <- switch(vcov,
              robust=,
              sandwich=iI%*%J%*%iI,##iJ%*%I%*%iJ,
              score=,
              outer=Inverse(J),
              hessian=iI              
              )
  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  rho <- ifelse(itrname=="log","U","rho")
  if (!eqmarg)
    rownames(cc) <- c(paste(rnames1,rep(c(1,2),each=length(rnames1)),sep="."),
                      paste(itrname,"(",rho,")",sep=""))
  else
    rownames(cc) <- c(rnames1,paste(itrname,"(",rho,")",sep=""))
  rownames(V) <- colnames(V) <- rownames(cc)
  npar <- list(intercept=attributes(terms(formula))$intercept,
              pred=nrow(attributes(terms(formula))$factor)-1)
  if (!eqmarg) npar <- lapply(npar,function(x) x*2)
  npar$var <- 1##nrow(cc)-sum(unlist(npar))
  N <- with(MyData, c(n=nrow(XX0)*2+length(margidx), pairs=nrow(XX0)))
  val <- list(coef=cc,N=N,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall, model=model,msg=msg,npar=npar,
              SigmaFun=SigmaFun)
  class(val) <- "biprobit"
  return(val)
}

###}}} biprobit

###{{{ utilities

ExMarg <- function(Y0,XX0,W0,dS0,midx1=seq(ncol(XX0)/2),midx2=seq(ncol(XX0)/2)+ncol(XX0)/2,eqmarg=TRUE,allmarg=FALSE) {  
  ii1 <- which(is.na(Y0[,2]) & !is.na(Y0[,1]))
  ii2 <- which(is.na(Y0[,1]) & !is.na(Y0[,2]))
  ii0 <- which(is.na(Y0[,1]) & is.na(Y0[,2]))
  margidx <- c(ii1,ii2)
  id1 <- id2 <-  NULL
  both <- setdiff(seq(nrow(Y0)),c(ii1,ii2,ii0))
  id <- seq_len(length(both))
  if (allmarg) {
    ##    id <- seq_len(length(both))
    id1 <- c(seq_len(length(ii1))+length(id), id)
    id2 <- c(seq_len(length(ii2))+length(id)+length(id1), id)
    ii1 <- c(ii1,both)
    ii2 <- c(ii2,both)
  }
  Y0_marg <- XX0_marg <- X0_marg1 <- X0_marg2 <- dS0_marg <- W0_marg <- NULL
  if (length(margidx)>0) {
    Y0_marg <- cbind(c(Y0[ii1,1],Y0[ii2,2]))
    X0_marg1 <- XX0[ii1,midx1,drop=FALSE]
    X0_marg2 <- XX0[ii2,midx2,drop=FALSE]
    dS0_marg <- dS0[,1,drop=FALSE]
    if (eqmarg) {
      XX0_marg <- rbind(X0_marg1,X0_marg2)
    } else {
      XX0_marg <- XX0[c(ii1,ii2),,drop=FALSE]
    }    
    if (!is.null(W0)) {
      W0_marg <- cbind(c(W0[ii1,1],W0[ii2,2]))
      W0 <- W0[-c(margidx,ii0),,drop=FALSE]
    }
    Y0 <- Y0[-c(margidx,ii0),,drop=FALSE]
    XX0 <- XX0[-c(margidx,ii0),,drop=FALSE]    
  }
  res <- list(Y0=Y0,XX0=XX0,W0=W0,
              Y0_marg=Y0_marg, XX0_marg=XX0_marg,
              X0_marg1=X0_marg1, X0_marg2=X0_marg2,
              dS0_marg=dS0_marg, W0_marg=W0_marg,
              id=id, idmarg=c(id1,id2),
              ii1=ii1,
              margidx=margidx)
}


RoundMat <- function(cc,digits = max(3, getOption("digits") - 2),...) format(round(cc,max(1,digits)),digits=digits)

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
  attributes(res)$warning <- any(svdX$d<tol)
  return(res)
}

blockdiag <- function(x,...,pad=0) {
  if (is.list(x)) xx <- x  else xx <- list(x,...)
  rows <- unlist(lapply(xx,nrow))
  crows <- c(0,cumsum(rows))
  cols <- unlist(lapply(xx,ncol))
  ccols <- c(0,cumsum(cols))
  res <- matrix(pad,nrow=sum(rows),ncol=sum(cols))
  for (i in 1:length(xx)) {
    idx1 <- 1:rows[i]+crows[i]; idx2 <- 1:cols[i]+ccols[i]
    res[idx1,idx2] <- xx[[i]]
  }
  colnames(res) <- unlist(lapply(xx,colnames)); rownames(res) <- unlist(lapply(xx,rownames))
  return(res)
}

Specials <- function(f,spec,split2="+",...) {
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
  if (is.null(pos)) return(NULL)
  x <- rownames(attributes(tt)$factors)[pos]
  st <- gsub(" ","",x)
  res <- unlist(strsplit(st,"[()]"))[2]
  if (is.null(split2)) return(res)
  unlist(strsplit(res,"+",fixed=TRUE))
}

###}}} utilities
