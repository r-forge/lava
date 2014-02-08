##' Ordinal regression models
##'
##' @title Univariate cumulative link regression models
##' @param formula formula
##' @param data data.frame
##' @param offset offset
##' @param family family (default proportional odds)
##' @param start optional starting values
##' @param fast If TRUE standard errors etc. will not be calculated
##' @param ... Additional arguments to lower level functions
##' @export
##' @author Klaus Kähler Holst
ordreg <- function(formula,data=parent.frame(),offset,family=binomial("logit"),start,fast=FALSE,...) {
    y <- ordered(model.frame(update(formula,.~0),data)[,1])
    lev <- levels(y)
    X <- model.matrix(update(formula,.~.+1),data=data)[,-1,drop=FALSE]
    up <- new.env()
    assign("h",family$linkinv,envir=up)
    assign("dh",family$mu.eta,envir=up)
    assign("y",as.numeric(y),envir=up)
    assign("X",X,envir=up)    
    assign("K",nlevels(y),envir=up)
    assign("n",length(y),envir=up)
    assign("p",NCOL(X),envir=up)    
    assign("threshold",
           function(theta,K) {
               a <- theta[1]
               if (K>2) a <- cumsum(c(a,exp(theta[seq(K-2)+1L])))
               return(a)
           }, envir=up)
    assign("dthreshold",
           function(theta,K) {
               Da <- matrix(0,K,K-1)
               Da[seq(K-1),1L] <- 1L
               for (i in seq_len(K-2)+1) Da[seq(i,K-1),i] <- exp(theta[i])
               Da
           },envir=up)
    ff <- function(theta) -ordreg_logL(theta,up)
    gg <- function(theta) -ordreg_score(theta,up)
    if (missing(start)) start <- with(up,c(rep(-1,K-1),rep(0,p)))
    op <- nlminb(start,ff,gg)
    cc <- op$par;
    if (fast) return(structure(cc,threshold=up$threshold(cc,up$K)))
    nn <- c(paste(lev[-length(lev)], lev[-1L], sep = "|"),
                   colnames(X))
    I <- ordreg_hessian(cc,up)
    names(cc) <- nn
    dimnames(I) <- list(nn,nn)
    res <- list(vcov=solve(I),coef=cc,call=match.call(),up=up)
    structure(res,class="ordreg")
}

##' @S3method print ordreg
print.ordreg <- function(x,...) {
    cat("Call:\n"); print(x$call)
    cat("\nParameter Estimates:\n")
    print(x$coef)
}

##' @S3method score ordreg
score.ordreg <- function(x,p=coef(x),indiv=FALSE,...) {
    ordreg_score(coef(x),x$up)
    if (!indiv) return(colSums(x$up$score))
    x$up$score
}

##' @S3method coef ordreg
logLik.ordreg <- function(object,p=coef(object),indiv=FALSE,...) {
    ordreg_logL(p,object$up)
    res <- log(object$up$pr)    
    if (!indiv) res <- sum(res)
    structure(res,nall=length(object$up$pr),nobs=object$up$pr,df=length(p),class="logLik")
}

##' @S3method coef ordreg
coef.ordreg <- function(object,...) object$coef

##' @S3method vcov ordreg
vcov.ordreg <- function(object,...) object$vcov

ordreg_logL <- function(theta,env,indiv=...) {
    if (length(theta)!=with(env,p+K-1)) stop("Wrong dimension")
    env$theta <- theta
    if (env$p>0) beta <- with(env,theta[seq(p)+K-1])
    alpha <- with(env, threshold(theta,K))
    env$alpha <- alpha
    env$beta <- beta
    if (env$p>0) eta <- env$X%*%beta else eta <- cbind(rep(0,env$n))
    env$lp <- kronecker(-eta,rbind(alpha),"+")
    F <- with(env,h(lp))
    Pr <- cbind(F,1)-cbind(0,F)
    pr <- Pr[with(env,cbind(seq(n),as.numeric(y)))]
    env$pr <- pr
    sum(log(pr))
}

ordreg_score <- function(theta,env,...) {
    if (!identical(theta,env$theta)) ordreg_logL(theta,env)
    Da <- with(env,dthreshold(theta,K))
    dF <- with(env, cbind(dh(lp),0))
    idx1 <- with(env,which(as.numeric(y)==1))
    S1 <- cbind(Da[as.numeric(env$y),,drop=FALSE],-env$X)
    S1 <- dF[with(env,cbind(seq(n),as.numeric(y)))]*S1
    y2 <- env$y-1; y2[idx1] <- env$K
    S2 <- cbind(Da[y2,,drop=FALSE],-env$X)
    S2 <- dF[cbind(seq(env$n),y2)]*S2
    env$score <- 1/env$pr*(S1-S2)
    colSums(env$score)
}
ordreg_hessian <- function(theta,env,...) {
    numDeriv::jacobian(function(p) ordreg_score(p,env,...),theta)
}


