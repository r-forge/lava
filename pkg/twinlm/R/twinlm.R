
###{{{ twinlm

twinlm <- function(formula, data, type=c("ace"), twinid="id", status="zyg", DZ, twinnum="twinnum",weight=NULL,binary=FALSE,probitscale=1,keep=weight,estimator="gaussian",...) {
  type <- tolower(type)
  if ("u" %in% type) type <- c("ue")
  
  varnames <- all.vars(formula)
  latentnames <- c("a1","a2","c1","c2","d1","d2","e1","e2")
  if (any(latentnames%in%varnames))
    stop(paste(paste(latentnames,collapse=",")," reserved for names of latent variables.",sep=""))
  cl <- match.call()
  mf <- model.frame(formula,data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  formula <- update(formula, ~ . + 1)
  yvar <- getoutcome(formula)
  
  if (is.factor(data[,yvar]) | is.character(data[,yvar])) {
    data[,yvar] <- as.numeric(as.factor(data[,yvar]))-1
    binary <- TRUE
  }  
  
  opt <- options(na.action="na.pass")
  mm <- model.matrix(formula,mf)
  options(opt)
  
  covars <- colnames(mm)
  if (attr(terms(formula),"intercept")==1)
    covars <- covars[-1]
  if(length(covars)<1) covars <- NULL

  zygstat <- data[,status]
  if(!is.factor(zygstat)) {
    zygstat <- as.factor(zygstat)
  }
  zyglev <- levels(zygstat)
  if (length(zyglev)>2) stop("Only support for two zygosity levels")
  ## Get data on wide format and divide into two groups by zygosity
  if (!twinnum%in%names(data)) {
    mynum <- rep(1,nrow(data))
    mynum[zygstat==zyglev[1]][duplicated(data[zygstat==zyglev[1],twinid])] <- 2
    mynum[zygstat==zyglev[2]][duplicated(data[zygstat==zyglev[2],twinid])] <- 2
    data[,twinnum] <- mynum
  }
  
  cur <- cbind(data[,c(yvar,keep)],as.numeric(data[,twinnum]),as.numeric(data[,twinid]),as.numeric(zygstat));
  colnames(cur) <- c(yvar,keep,twinnum,twinid,status)
  mydata <- cbind(cur,mm)
  if (missing(DZ)) {
    warning("Using first level, `",zyglev[1],"', in status variable as indicator for 'dizygotic'", sep="")
    DZ <- zyglev[1]
  }
  myDZ <- which(levels(zygstat)==DZ)
  myMZ <- setdiff(1:2,myDZ)

  data1 <- mydata[mydata[,status]==myMZ,,drop=FALSE]
  data2 <- mydata[mydata[,status]==myDZ,,drop=FALSE]
    
  data1.1 <- data1[which(data1[,twinnum]==1),c(twinid,status,yvar,keep,covars)]; colnames(data1.1) <- c(twinid,status,paste(colnames(data1.1)[-c(1,2)],".1",sep=""))
  data1.2 <- data1[which(data1[,twinnum]==2),c(twinid,yvar,keep,covars),drop=FALSE]; colnames(data1.2) <- c(twinid, paste(colnames(data1.2)[-1],".2",sep=""))

  ##Missing data?
  id1 <- data1.1[,twinid]
  id2 <- data1.2[,twinid]
  d1.mis <- setdiff(id2,id1) # Id's not in data1.1
  d2.mis <- setdiff(id1,id2) # Id's not in data1.2
  if (length(d1.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data1.1), nrow=length(d1.mis))
    d.temp[,1] <- d1.mis; colnames(d.temp) <- colnames(data1.1)
    data1.1 <- rbind(data1.1,d.temp)
  }
  if (length(d2.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data1.2), nrow=length(d2.mis))
    d.temp[,1] <- d2.mis; colnames(d.temp) <- colnames(data1.2)
    data1.2 <- rbind(data1.2,d.temp)
  } 

  data2.1 <- data2[data2[,twinnum]==1,c(twinid,status,yvar,keep,covars)]; colnames(data2.1) <- c(twinid,status,paste(colnames(data2.1)[-c(1,2)],".1",sep=""))
  data2.2 <- data2[data2[,twinnum]==2,c(twinid,yvar,keep,covars),drop=FALSE]; colnames(data2.2) <- c(twinid, paste(colnames(data2.2)[-1],".2",sep=""))

  id1 <- data2.1[,twinid]
  id2 <- data2.2[,twinid]
  d1.mis <- setdiff(id2,id1) # Id's not in data1.1
  d2.mis <- setdiff(id1,id2) # Id's not in data1.2
  if (length(d1.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data2.1), nrow=length(d1.mis))
    d.temp[,1] <- d1.mis; colnames(d.temp) <- colnames(data2.1)
    data2.1 <- rbind(data2.1,d.temp)
  }
  if (length(d2.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data2.2), nrow=length(d2.mis))
    d.temp[,1] <- d2.mis; colnames(d.temp) <- colnames(data2.2)
    data2.2 <- rbind(data2.2,d.temp)
  }  

  wide1 <- merge(x=data1.1,y=data1.2, by=twinid); wide1[,status] <- myMZ
  wide2 <- merge(x=data2.1,y=data2.2, by=twinid); wide2[,status] <- myDZ
    
  ## ###### The SEM
  outcomes <- paste(yvar,".",1:2,sep="")
  model1<-lvm(outcomes,silent=TRUE)
  f1 <- as.formula(paste(outcomes[1]," ~ ."))
  f2 <- as.formula(paste(outcomes[2]," ~ ."))
##  parameter(model1) <- ~sdu1+sdu2
  regression(model1,silent=TRUE) <- update(f1, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e1,lambda[e]) + f(u,sdu1))
  regression(model1,silent=TRUE) <- update(f2, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e2,lambda[e]) + f(u,sdu1))
  latent(model1) <- ~ a1+c1+d1+e1+e2+u
  intercept(model1,latent(model1)) <- 0
  if (!is.null(covars))
    for (i in 1:length(covars)) {
      regfix(model1, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta[",i,"]",sep="")
      regfix(model1, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta[",i,"]",sep="")
    }
  covariance(model1) <- update(f1, . ~  v(0))
  covariance(model1) <- update(f2, . ~  v(0))
  covfix(model1, latent(model1), var2=NULL) <- 1
  intfix(model1,outcomes) <- "mu1"
  model2 <- model1
##  parameter(model2) <- ~sdu1+sdu2
  regression(model2) <- update(f1,.~f(u,sdu2))
  regression(model2) <- update(f2,.~f(u,sdu2))
  cancel(model2) <- update(f2, . ~ a1)
  cancel(model2) <- update(f2, . ~ d1)
  regression(model2,silent=TRUE) <- update(f2, . ~ f(a2,lambda[a]))
  regression(model2,silent=TRUE) <- update(f2, . ~ f(d2,lambda[d]))
  
  covariance(model2) <- a1 ~ f(a2,0.5)
  covariance(model2) <- d1 ~ f(d2,0.25)
  latent(model2) <- ~ a2+d2
  intercept(model2, ~ a2+d2) <- 0
  covariance(model2) <- c(a2,d2) ~ v(1)
  full <- list(model1,model2)
  ## #######
  isA <- length(grep("a",type))>0
  isC <- length(grep("c",type))>0
  isD <- length(grep("d",type))>0
  isE <- length(grep("e",type))>0
  isU <- length(grep("u",type))>0
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
  if (!isE) {
    kill(model1) <- ~ e1 + e2
    kill(model2) <- ~ e1 + e2
  }
  if (!isU) {
    kill(model1) <- ~ u##+sdu2
    kill(model2) <- ~ u##+sdu1
  }
  if (isU & isE) {
    regression(model1,outcomes[1],"e1") <- "lambda[e2]"
    regression(model1,outcomes[2],"e2") <- "lambda[e2]"
  }

  ## Full rank covariate/design matrix?
  for (i in covars) {
    myvars <- paste(i,c(1,2),sep=".")
    dif <- wide1[,myvars[1]]-wide1[,myvars[2]]   
    mykeep <- myvars
    if (all(na.omit(dif)==00)) {
      mykeep <- mykeep[-2]
    }   
    trash <- setdiff(myvars,mykeep)
    if (length(mykeep)==1) {
      regression(model1, to=outcomes[2], from=mykeep) <- regfix(model1)$label[trash,outcomes[2]]
      kill(model1) <- trash
    }

    dif <- wide2[,myvars[1]]-wide2[,myvars[2]]   
    mykeep <- myvars
    if (all(na.omit(dif)==00)) {
      mykeep <- mykeep[-2]
    }  
    trash <- setdiff(myvars,mykeep)
    if (length(mykeep)==1) {
      regression(model2, to=outcomes[2], from=mykeep) <- regfix(model2)$label[trash,outcomes[2]]
      kill(model2) <- trash
    }
  }

  if (!is.null(weight)) {
    weight <- paste(weight,1:2,sep=".")
    ## wide1[which(is.na(wide1[,weight[1]])),weight[1]] <- 0
    ## wide1[which(is.na(wide1[,weight[2]])),weight[2]] <- 0
    ## wide2[which(is.na(wide2[,weight[1]])),weight[1]] <- 0
    ## wide2[which(is.na(wide2[,weight[2]])),weight[2]] <- 0    
    ## wide1[which(is.na(wide1[,outcomes[1]])),c(outcomes[1],weight[1])] <- 0
    ## wide1[which(is.na(wide1[,outcomes[2]])),c(outcomes[2],weight[2])] <- 0
    ## wide2[which(is.na(wide2[,outcomes[1]])),c(outcomes[1],weight[1])] <- 0
    ## wide2[which(is.na(wide2[,outcomes[2]])),c(outcomes[2],weight[2])] <- 0
    if (!binary) estimator <- "weighted"
  }

  if (binary) {
    binary(model1) <- outcomes
    covariance(model1,outcomes) <- 0
    binary(model2) <- outcomes
    covariance(model2,outcomes) <- 0
    if (!is.null(probitscale))
      if (isE) {
        regression(model2,from=c("e1","e2"), to=outcomes) <-
          rep(probitscale,2)
        if (!isU)
          regression(model1,from=c("e1","e2"), to=outcomes) <-
            rep(probitscale,2) 
      } else { ## Testing:
        regression(model1,from=c("c1"), to=outcomes) <- probitscale
        regression(model2,from=c("c1"), to=outcomes) <- probitscale 
      }
  }

  ## Estimate
  newkeep <- unlist(sapply(keep, function(x) paste(x,1:2,sep=".")))
  suppressWarnings(mg <- multigroup(list(model1,model2), list(wide1,wide2), missing=TRUE,fix=FALSE,keep=newkeep))
  if (is.null(estimator)) return(mg)

  if (binary) {
    e <- estimate(mg,weight2=weight,estimator=estimator,fix=FALSE,...)
  } else {
    e <- estimate(mg,weight=weight,estimator=estimator,fix=FALSE,...)
  }
  res <- list(coefficients=e$opt$estimate, vcov=e$vcov, estimate=e, model=mg, full=full, call=cl, data=data, status=status, twinid=twinid, twinnum=twinnum, binary=binary, type=type, model.mz=model1, model.dz=model2, data.mz=wide1, data.dz=wide2,
              probitscale=probitscale)
  class(res) <- "twinlm"
  return(res)
}

###}}} twinlm

###{{{ print.twinlm

print.twinlm <- function(x,...) {
  print(summary(x,...))
  invisible(x)
}

###}}} print.twinlm

###{{{ summary.twinlm

summary.twinlm <- function(object,...) {
  e <- object$estimate
  zygtab <- with(object, table(data[,status]))
  theta <- pars(e)
  theta.sd <- sqrt(diag(e$vcov))
  myest <- cbind(theta,theta.sd,(Z <- theta/theta.sd),2*(1-pnorm(abs(Z))))
  colnames(myest) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
  
  if (length(grep("u",object$type))>0) {

    cc <- coef(e)
    pi <- seq_len(nrow(myest))
    pp <- modelPar(e$model,pi)$p
    nn <- rep(NA,nrow(myest))
    nn[pp[[1]]] <- coef(e$model$lvm[[1]])
    nn.na <- which(is.na(nn))
    i2 <- 0
    for (i in 1:length(pp)) {
      if (nn.na[1]%in%pp[[i]]) {
        i2 <- i
        break;
      }      
    }
    pp.i2 <- which(pp[[i2]]%in%nn.na)
    parnum <- pp[[i2]][pp.i2]
    nn[nn.na] <- coef(e$model$lvm[[i2]])[pp.i2]
    nn <- gsub(".1","",nn,fixed=TRUE)
    nn <- gsub(".2","",nn,fixed=TRUE)
    u.idx <- c(grep("<-u",nn))
    e.idx <- c(grep("<-e",nn))

    lastpos <- all(parnum>=u.idx) ## i2-coef positioined after coef. of model 1
    isMZ <- "sdu1"%in%parlabels(e$model$lvm[[i2]])
    if ((lastpos & isMZ)|(!lastpos & !isMZ)) {
      u.idx <- rev(u.idx); e.idx <- rev(e.idx)
    }
    nn <- c(nn); nn[u.idx] <- c("MZ:sd(U)","DZ:sd(U)")
    if (length(e.idx)==1) {
      nn[e.idx] <- "MZ:sd(E)"
    } else {
      nn[e.idx] <- c("MZ:sd(E)","DZ:sd(E)")
    }    
    rownames(myest) <- nn
    neword <- c(setdiff(seq_len(nrow(myest)),c(e.idx,u.idx)),e.idx,u.idx)

    MZcc <- DZcc <- NULL
    if (object$binary) {
      MZest <- myest[-match("DZ:sd(U)",rownames(myest)),1]
      DZest <- myest[-match(c("MZ:sd(E)","MZ:sd(U)"),rownames(myest)),1]
      M1 <- moments(object$model.mz,MZest,cond=FALSE)
      M2 <- moments(object$model.dz,DZest,cond=FALSE)
      myidx <- index(object$model.mz)$endo.obsidx
      MZcc <- with(M1, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
      myidx <- index(object$model.dz)$endo.obsidx
      DZcc <- with(M2, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
##      browser()
    }
    
    K <- ifelse (object$binary,object$probitscale,0)
    ##L <- binomial("logit")
    logit <- function(p) log(p/(1-p))
    tigol <- function(z) 1/(1+exp(-z))
    if (length(e.idx)==1) {
      corDZ <- function(x) (x[2]^2)/(x[2]^2+K)
      corMZ <- function(x) (x[1]^2)/(x[1]^2+x[3]^2)
    } else {
      corDZ <- function(x) (x[2]^2)/(x[2]^2+x[4]^2)
      corMZ <- function(x) (x[1]^2)/(x[1]^2+x[3]^2)
    }
    h <- function(x) 2*(corMZ(x)-corDZ(x))
    
    dh <- function(x) {
      x2 <- x^2;
      s2 <- ifelse(length(e.idx)==1,x2[2]+K^2,x2[2]+x2[4])
      s1 <- ifelse(length(e.idx)==1,x2[1]+x2[3],x2[1]+x2[3])
      if (length(e.idx)==1) {
        De <- -4*x[3]*x2[1]/s1^2
      } else {
        De <- c(-4*x[3]*x2[1]/s1^2,4*x[4]*x2[2]/s2^2)
      }      
      c(2*c(1/s1^2,-1/s2^2)*(2*x[1:2]*c(s1,s2)-x2[1:2]*2*x[1:2]),De)      
    }

    
    dlogit <- function(p) 1/(p*(1-p))
    logith <- function(x) logit(h(x))
    dlogith <- function(x) dlogit(h(x))*dh(x)

    V <- e$vcov[c(u.idx,e.idx),c(u.idx,e.idx)]
    b <- pars(e)[c(u.idx,e.idx)]
    Db <- dlogith(b)
    sd. <- (t(Db)%*%V%*%Db)^0.5
    hval <- c(h(b),(t(dh(b))%*%V%*%(dh(b)))^0.5)
    hci <- tigol(logith(b)+qnorm(0.975)*c(-1,1)*sd.)
    names(hci) <- c("2.5%","97.5%")
    res <- list(estimate=myest[neword,], zyg=zygtab,
                varEst=NULL, varSigma=NULL, heritability=hval, hci=hci,
                corMZ=corMZ(b), corDZ=corDZ(b),
                concMZ=MZcc, concDZ=DZcc)                
    class(res) <- "summary.twinlm"
    return(res)
  }

  rownames(myest) <- gsub(".1","",coef(Model(Model(e))[[1]],
                                       mean=e$meanstructure, silent=TRUE),
                          fixed=TRUE)
  rownames(myest) <- gsub(".2","",rownames(myest),fixed=TRUE)
##  rownames(myest) <- coef(Model(Model(e))[[1]],
##                                      mean=e$meanstructure, silent=TRUE)
                     
  lambda.idx <- sapply(c("<-a1","<-c1","<-d1","<-e1"),function(x) grep(x,rownames(myest)))
  lambda.w <- which(sapply(lambda.idx, function(x) length(x)>0))
  rownames(myest)[unlist(lambda.idx)] <- paste("sd(",c("A)","C)","D)","E)"),sep="")[lambda.w]


  varEst <- rep(0,4)
  varEst[lambda.w] <- myest[unlist(lambda.idx),1]
  varSigma <- matrix(0,4,4);
  varSigma[lambda.w,lambda.w] <- e$vcov[unlist(lambda.idx),unlist(lambda.idx)]

  L <- binomial("logit")
  varcomp <- c()
  genpos <- c()
  pos <- 0
  if ("e1"%in%latent(e) & object$binary) {
    if (length(object$probitscale)>0)
      varEst[4] <- object$probitscale
    if ("a1"%in%latent(e)) { varcomp <- "lambda[a]"; pos <- pos+1; genpos <- c(genpos,pos) }
    if ("c1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[c]"); pos <- pos+1 }
    if ("d1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[d]"); pos <- pos+1;
                             genpos <- c(genpos,pos) }
    f <- paste("h2~",paste(varcomp,collapse="+"))
    constrain(e, as.formula(f)) <- function(x) L$linkfun(sum(x[genpos]^2)/sum(c(x^2,varEst[4])))
    
  } else {
    varcomp <- c()
        if ("a1"%in%latent(e)) { varcomp <- "lambda[a]"; pos <- pos+1;
                                 genpos <- c(genpos,pos) }
    if ("c1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[c]"); pos <- pos+1 }
    if ("d1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[d]"); pos <- pos+1;
                             genpos <- c(genpos,pos) }
    if ("e1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[e]"); pos <- pos+1 }
    f <- paste("h2~",paste(varcomp,collapse="+"))

    constrain(e, as.formula(f)) <- function(x) L$linkfun(sum(x[genpos]^2)/sum(x^2))
  }
  ci.logit <- L$linkinv(constraints(e)["h2",5:6])
  
  h <- function(x) (x[1]^2)/(sum(x^2))
  dh <- function(x) {
    y <- x^2
    cbind(1/sum(y)^2*c(2*x[1]*(sum(y)-y[1]), -2*x[2:4]*y[1]))
  }
  ##  f(x1,x2,x3,x4) = x1/(x1+x2+x3+x4) = x1/s
  ## Quotient rule: (u/v)' = (u'v - uv')/v^2
  ##  f1(x1,x2,x3,x4) = (s - x1)/s^2 = (x2+x3+x4)/s^2
  ##  f2(x1,x2,x3,x4) = -1/(s)^2
  h2 <- function(x) (x[1]^2+x[3]^2)/sum(x^2)
  dh2 <- function(x) {
    y <- x^2
    cbind(1/sum(y)^2*c(
                       2*x[1]*(sum(y)-(y[1]+y[3])),
                       -2*x[2]*(y[1]+y[3]),
                       2*x[3]*(sum(y)-(y[1]+y[3])),
                       -2*x[4]*(y[1]+y[3])
                       ))
  }
  hval <- cbind(h(varEst), (t(dh(varEst))%*%varSigma%*%(dh(varEst)))^0.5)
  colnames(hval) <- c("Estimate", "Std.Err"); rownames(hval) <- "h squared"
  h2val <- cbind(h2(varEst), (t(dh2(varEst))%*%varSigma%*%(dh2(varEst)))^0.5)
  colnames(h2val) <- c("Estimate", "Std.Err"); rownames(h2val) <- "h squared"

  corMZ <- sum(varEst[1:3]^2)/sum(varEst^2)
  corDZ <- sum(varEst[1:3]^2*c(0.5,1,0.25))/sum(varEst^2)

  MZcc <- DZcc <- NULL
  if (object$binary) {    
    M1 <- moments(object$model.mz,coef(object),cond=FALSE)
    M2 <- moments(object$model.dz,coef(object),cond=FALSE)
    myidx <- index(object$model.mz)$endo.obsidx
    MZcc <- with(M1, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
    DZcc <- with(M2, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
##    browser()
  }

  
  
  res <- list(estimate=myest, zyg=zygtab, varEst=varEst, varSigma=varSigma, hval=hval, heritability=h2val, hci=ci.logit, corMZ=corMZ, corDZ=corDZ,
              concMZ=MZcc, concDZ=DZcc)                              
  class(res) <- "summary.twinlm"
  return(res)
}

###}}} summary.twinlm

###{{{ print.summary.twinlm

print.summary.twinlm <- function(x,signif.stars=FALSE,...) {
  printCoefmat(x$estimate,signif.stars=signif.stars,...)
  cat("\n")
  myzyg <- with(x,zyg)
  names(myzyg) <- c("Group 1 (DZ)", "Group 2 (MZ)")
  print(myzyg)
##   cat("\n")
##   mynames <- c("sigmaA","sigmaC","sigmaD","sigmaE")
##   for (i in 1:4) {
##     cat(mynames[i], "=", x$varEst[i], "\n")
##   }  
  cat("\n")
  ##  cat("hn2 = ", hn2, "\t hb2 = ", hb2, "\n\n")
##  cat("Narrow-sense heritability (additive genetic factors):\n")
##  print(x$hval)
##  cat("\n")  
  cat("heritability (total genetic factors):\n")
  h <- with(x, c(heritability,hci));
  names(h) <- c("Estimate","Std.Err",names(x$hci))
  h <- na.omit(h)
  print(h)  
  cat("\n")  
  cat("Correlation within MZ:", x$corMZ, "\n")
  cat("Correlation within DZ:", x$corDZ, "\n")
  cat("\n")
  if (!is.null(x$concMZ)) {
    cat("Concordance MZ:", x$concMZ, "\n")
    cat("Concordance DZ:", x$concDZ, "\n")
    cat("\n")
  }  
}

###}}} print.summary.twinlm

###{{{ compare.twinlm
compare.twinlm <- function(object,...) {
  objects <- list(object,...)
  if (length(objects)<2)
    return(summary(objects))
  res <- list()
  for (i in 1:(length(objects)-1)) {
    res <- c(res, list(compare(objects[[i]]$estimate,objects[[i+1]]$estimate)))
  }
  if (length(res)==1)
    return(res[[1]])
  return(res)
}
###}}} compare.twinlm

###{{{ plot.twinlm
plot.twinlm <- function(x,diag=TRUE,labels=TRUE,...) {
  op <- par(mfrow=c(2,1))
  plot(x$model,...)
  par(op)
}
###}}}

###{{{ vcov.twinlm
vcov.twinlm <- function(object,...) {
  return(object$vcov)
}
###}}} vcov.twinlm

###{{{ logLik.twinlm
logLik.twinlm <- function(object,...) logLik(object$estimate,...)
###}}} logLik.twinlm

###{{{ model.frame.twinlm
model.frame.twinlm <- function(formula,...) {
  return(formula$estimate$model$data)
}
###}}} model.frame.twinlm

###{{{ twinsim

twinsim <- function(n=100,k1=c(),k2=1,mu=0,lambda=c(1,1,1),randomslope=NULL,type="ace",binary=FALSE,...) {
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
  regfix(sim.model1,to=outcomes[1],from="e1") <- lambda[3]
  regfix(sim.model1,to=outcomes[2],from="e2") <- lambda[3]
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
  regfix(sim.model2,to=outcomes[1],from="e1") <- lambda[3]
  regfix(sim.model2,to=outcomes[2],from="e2") <- lambda[3]
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

###}}}

###{{{ acde

"acde" <- function(x,...) UseMethod("acde")
acde.twinlm <- function(x,...) {
  m <- x$estimate$model$lvm[[1]]
  lambdas <- c("lambda[a]","lambda[c]","lambda[d]","lambda[e]")
  ACDE <- lambdas%in%as.vector(m$par)
  lcur <- lambdas[ACDE]
  for (l in lcur) {
    pos <- which(lcur%in%l)
    par <- substr(strsplit(l,"[",fixed=TRUE)[[1]][2],1,1)
    f <- as.formula(paste(par,"~",paste(lcur,collapse="+")))
    myfun <- eval(parse(text=paste("function(x) x[",pos,"]^2/sum(x^2)")))
    constrain(x$estimate,f) <- myfun ##function(x) x[get("pos")]^2/sum(x^2)
  }
  constraints(x$estimate)
}

###}}} acde
