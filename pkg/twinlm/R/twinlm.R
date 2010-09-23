###{{{ twinlm

twinlm <- function(formula, data, type=c("ace"), twinid="id", status="zyg", DZ, twinnum="twin",weight, binary=FALSE,keep=NULL,debug=FALSE,estimator="gaussian",...) {

  if (!missing(weight)) {
    if (is.character(weight)) {
      weight <- data[,weight]
    }
  } else {
    weight <- NULL
  }
  
  varnames <- all.vars(formula)
  latentnames <- c("a1","a2","c1","c2","d1","d2","e1","e2")
  if (any(latentnames%in%varnames))
    stop(paste(paste(latentnames,collapse=",")," reserved for names of latent variables.",sep=""))
  cl <- match.call()
  mf <- model.frame(formula,data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  formula <- update(formula, ~ . + 1)
  yvar <- getoutcome(formula)

  require(survival)
  if (is.Surv(data[,yvar])) {
    weight <- 1-data[,yvar][,2]
    data[,yvar] <- data[,yvar][,1]
    estimator <- "tobit"
    data[,"_weight"] <- weight
    keep <- "_weight"
  }
  if (is.factor(data[,yvar]) | is.character(data[,yvar])) {
    data[,yvar] <- 1-as.numeric(as.factor(data[,yvar]))
    binary <- TRUE
  }  
  
  opt <- options(na.action="na.pass")
  ##  mm <- naresid(na.exclude, model.matrix(formula,data))
  mm <- model.matrix(formula,mf)
  options(opt)
  
  covars <- colnames(mm)
  if (attr(terms(formula),"intercept")==1)
    covars <- covars[-1]
  if(length(covars)<1) covars <- NULL
  
    zygstat <- data[,status]
    if(!is.factor(zygstat)) {
      ##      warning("Transforming zygosity status variable to a factor.")
      zygstat <- as.factor(zygstat)
    }
    zyglev <- levels(zygstat)
    if (length(zyglev)>2) stop("Only support for two zygosity levels")
  
    ## Get data on wide format and divide into two groups by zygosity
    if (!twinnum%in%names(data)) {
      mynum <- rep(2,nrow(data))
      firsttwin.Z1 <- sapply(unique(data[zygstat==zyglev[1],twinid]), function(x) which(data[zygstat==zyglev[1],twinid]==x)[1])
      firsttwin.Z2 <- sapply(unique(data[zygstat==zyglev[2],twinid]), function(x) which(data[zygstat==zyglev[2],twinid]==x)[1])
      mynum[zygstat==zyglev[1]][firsttwin.Z1] <- 1
      mynum[zygstat==zyglev[2]][firsttwin.Z2] <- 1
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
    
  data1.1 <- data1[data1[,twinnum]==1,c(twinid,status,yvar,keep,covars)]; colnames(data1.1) <- c(twinid,status,paste(colnames(data1.1)[-c(1,2)],".1",sep=""))
  data1.2 <- data1[data1[,twinnum]==2,c(twinid,yvar,keep,covars),drop=FALSE]; colnames(data1.2) <- c(twinid, paste(colnames(data1.2)[-1],".2",sep=""))

  ##Missing data?
  id1 <- data1.1[,twinid]
  id2 <- data1.2[,twinid]
  d1.mis <- setdiff(id2,id1) # Id's not in data1.1
  d2.mis <- setdiff(id1,id2) # Id's not in data1.2
  if (length(d1.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data1.1), nrow=length(d1.mis))
    d.temp[,1] <- d1.mis
    data1.1 <- rbind(data1.1,d.temp)
  }
  if (length(d2.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data1.2), nrow=length(d2.mis))
    d.temp[,1] <- d2.mis
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
    d.temp[,1] <- d1.mis
    data2.1 <- rbind(data2.1,d.temp)
  }
  if (length(d2.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data2.2), nrow=length(d2.mis))
    d.temp[,1] <- d2.mis
    data2.2 <- rbind(data2.2,d.temp)
  }  

  wide1 <- merge(x=data1.1,y=data1.2, by=twinid)
  wide2 <- merge(x=data2.1,y=data2.2, by=twinid)

  
##  Debug("Models...",debug)
    
    ## ###### The SEM
    outcomes <- paste(yvar,".",1:2,sep="")
    model1<-lvm(outcomes,silent=TRUE)
    f1 <- as.formula(paste(outcomes[1]," ~ ."))
    f2 <- as.formula(paste(outcomes[2]," ~ ."))
    regression(model1,silent=TRUE) <- update(f1, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e1,lambda[e]))
    regression(model1,silent=TRUE) <- update(f2, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e2,lambda[e]))    
    latent(model1) <- ~ a1+c1+d1+e1+e2
    if (!is.null(covars))
      for (i in 1:length(covars)) {
        regfix(model1, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta[",i,"]",sep="")
        regfix(model1, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta[",i,"]",sep="")
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
    regression(model2,silent=TRUE) <- update(f2, . ~ f(a2,lambda[a]))
    regression(model2,silent=TRUE) <- update(f2, . ~ f(d2,lambda[d]))
    covariance(model2) <- a1 ~ f(a2,0.5)
    covariance(model2) <- d1 ~ f(d2,0.25)
    latent(model2) <- ~ a2+d2
    covariance(model2) <- c(a2,d2) ~ v(1)
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

  ## Full rank covariate/design matrix?
  for (i in covars) {
    myvars <- paste(i,c(1,2),sep=".")
    dif <- wide1[,myvars[1]]-wide1[,myvars[2]]   
##    keep <- myvars[which(!duplicated(wide1[,myvars], MARGIN=2))]
    mykeep <- myvars
    if (all(na.omit(dif)==00)) {
      mykeep <- mykeep[-2]
    }   
    ##    keep <- myvars[which(!duplicated(wide1[,myvars], MARGIN=2))]    
    trash <- setdiff(myvars,mykeep)
    if (length(mykeep)==1) {
##      regression(model1, to=outcomes[2]) <- paste(mykeep,"@",regfix(model1)$label[trash,outcomes[2]],sep="")
      regression(model1, to=outcomes[2], from=mykeep) <- regfix(model1)$label[trash,outcomes[2]]
      kill(model1) <- trash
      ##      wide1 <- subset(wide1, select=-
      ##      warning(paste(" 'dizygotic'")      
    }

    ##    keep <- myvars[which(!duplicated(wide2[,myvars], MARGIN=2))]
    dif <- wide2[,myvars[1]]-wide2[,myvars[2]]   
##    keep <- myvars[which(!duplicated(wide1[,myvars], MARGIN=2))]
    mykeep <- myvars
    if (all(na.omit(dif)==00)) {
      mykeep <- mykeep[-2]
    }  
    trash <- setdiff(myvars,mykeep)
    if (length(mykeep)==1) {
      regfix(model2)
##      regression(model2, to=outcomes[2]) <- paste(mykeep,"@",regfix(model2)$label[trash,outcomes[2]],sep="")
      regression(model2, to=outcomes[2], from=mykeep) <- regfix(model2)$label[trash,outcomes[2]]
      kill(model2) <- trash
    }
  }


  myweights <- NULL
  if (!is.null(weight)) {
    Weight <- cbind(data[,twinid],weight); colnames(Weight)[1] <- twinid
    idx1 <- which(data[,twinid]%in%wide1[,twinid])
    idx2 <- which(data[,twinid]%in%wide2[,twinid])    
##    weight1 <- Weight[1:nrow(data1),]
    weight1 <- Weight[idx1,]
    weight1.1 <- weight1[data1[,twinnum]==1,]
    weight1.2 <- weight1[data1[,twinnum]==2,]
    ##    weight2 <- Weight[1:nrow(data2)+nrow(data1),]
    weight2 <- Weight[idx2,]
    weight2.1 <- weight2[data2[,twinnum]==1,]
    weight2.2 <- weight2[data2[,twinnum]==2,]
    weight1 <- merge(x=weight1.1,y=weight1.2,by=twinid)[,-1]
    weight2 <- merge(x=weight2.1,y=weight2.2,by=twinid)[,-1]
##    for (i in 1:length(exogenous(model1)))
##      weight1 <- cbind(weight1,1)
##    for (i in 1:length(exogenous(model2)))
##      weight2 <- cbind(weight2,1)
    colnames(weight1) <- outcomes
    colnames(weight2) <- outcomes
    myweights <- list(as.matrix(weight1),as.matrix(weight2))
  }
  
  if (binary) {
    binary(model1) <- outcomes
    regfix(model1,from=c("e1","e2"), to=outcomes) <- list(1,1) 
    binary(model2) <- outcomes
    regfix(model2,from=c("e1","e2"), to=outcomes) <- list(1,1) 
  }


###Estimate
  newkeep <- as.vector(sapply(keep, function(x) paste(x,1:2,sep=".")))
##  return(list(model=list(model1,model2), data=list(wide1,wide2)))
  mg <- multigroup(list(model1,model2), list(wide1,wide2),missing=TRUE,fix=FALSE,keep=newkeep)
  if (is.null(estimator)) return(mg)
  e <- estimate(mg,weight=myweights,debug=debug,estimator=estimator,...)
  res <- list(coefficients=e$opt$estimate, vcov=e$vcov, estimate=e, model=mg, full=full, call=cl, data=data, status=status, twinid=twinid, twinnum=twinnum, binary=binary)
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
  theta <- pars(e)
  theta.sd <- sqrt(diag(e$vcov))
  myest <- cbind(theta,theta.sd,(Z <- theta/theta.sd),2*(1-pnorm(abs(Z))))
  colnames(myest) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
  rownames(myest) <- coef(Model(Model(e))[[1]], mean=e$meanstructure, silent=TRUE)
  lambda.idx <- sapply(c("<-a1","<-c1","<-d1","<-e1"),function(x) grep(x,rownames(myest)))
  lambda.w <- which(sapply(lambda.idx, function(x) length(x)>0))
  rownames(myest)[unlist(lambda.idx)] <- paste("sigma",c("A","C","D","E"),sep="")[lambda.w]


  varEst <- rep(0,4)
  varEst[lambda.w] <- myest[unlist(lambda.idx),1]
  varSigma <- matrix(0,4,4);
  varSigma[lambda.w,lambda.w] <- e$vcov[unlist(lambda.idx),unlist(lambda.idx)]
  zygtab <- with(object, table(data[,status]))

  L <- binomial(logit)
  varcomp <- c()
  genpos <- c()
  pos <- 0
  if ("e1"%in%latent(e) & object$binary) {    
    varEst[4] <- 1
    if ("a1"%in%latent(e)) { varcomp <- "lambda[a]"; pos <- pos+1; genpos <- c(genpos,pos) }
    if ("c1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[c]"); pos <- pos+1 }
    if ("d1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[d]"); pos <- pos+1;
                             genpos <- c(genpos,pos) }
    f <- paste("h2~",paste(varcomp,collapse="+"))
    constrain(e, as.formula(f)) <- function(x) L$linkfun(sum(x[genpos])^2/sum(x^2+1))
  } else {
    varcomp <- c()
        if ("a1"%in%latent(e)) { varcomp <- "lambda[a]"; pos <- pos+1;
                                 genpos <- c(genpos,pos) }
    if ("c1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[c]"); pos <- pos+1 }
    if ("d1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[d]"); pos <- pos+1;
                             genpos <- c(genpos,pos) }
    if ("e1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[e]"); pos <- pos+1 }
    f <- paste("h2~",paste(varcomp,collapse="+"))
    constrain(e, as.formula(f)) <- function(x) L$linkfun(sum(x[genpos])^2/sum(x^2))
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
  
  res <- list(estimate=myest, zyg=zygtab, varEst=varEst, varSigma=varSigma, hval=hval, h2val=h2val, ci.logit=ci.logit)
  class(res) <- "summary.twinlm"
  return(res)
}

###}}} summary.twinlm

###{{{ print.summary.twinlm

print.summary.twinlm <- function(x,signif.stars=FALSE,...) {
  printCoefmat(x$estimate,signif.stars=signif.stars,...)
  cat("\n")
  myzyg <- with(x,zyg)
  names(myzyg) <- c("Group 1 (MZ)", "Group 2 (DZ)")
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
  h <- c(x$h2val,x$ci.logit);
  names(h) <- c("Estimate","Std.Err",names(x$ci.logit))
  print(h)
  cat("\n")
  cat("Correlation within MZ:", sum(x$varEst[1:3]^2)/sum(x$varEst^2), "\n")
  cat("Correlation within DZ:", sum(x$varEst[1:3]^2*c(0.5,1,0.25))/sum(x$varEst^2), "\n") 
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
      randomslope(model1) <- c(covars.1[i],covars.2[i])
      
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
      randomslope(model2) <- c(covars.1[i],covars.2[i])
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
  intfix(sim.model1,outcomes) <- mu
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
  
  sim.model2 <- model2
  intfix(sim.model2,outcomes) <- mu
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
  d2 <- sim(sim.model2,n=n,...)

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
