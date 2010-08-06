lava.tobit.color.hook <- function(x,subset=vars(x),...) {
  return(list(vars=intersect(subset,binary(x)),col="indianred1"))
}

lava.tobit.estimate.hook <- function(x,data,weight,estimator,...) {
  dots <- list(...)
## Binary outcomes -> censored regression  
  if (length(binary(x))>0 & estimator%in%c("gaussian","tobit")) {
    if (is.null(weight)) {
      W <- data[,binary(x),drop=FALSE]; W[W==0] <- -1; colnames(W) <- binary(x)
      weight <- W
    } else {
      if (!all(binary(x)%in%colnames(data))) {
        W <- data[,binary(x),drop=FALSE]; W[W==0] <- -1; colnames(W) <- binary(x)
        weight[,binary(x)] <- W
      }
    }
    for (b in binary(x)) {
      data[!is.na(data[,b]),b] <- 0
    }
    ##    data[,binary(x)] <- 0
    estimator <- "tobit"
  }
##  if (!is.null(weight))
##  weight <- as.matrix(weight)
  
## Transform 'Surv' objects
  W <- mynames <- c()
  if (estimator%in%c("gaussian","tobit")) {
    for (i in setdiff(endogenous(x),binary(x))) {
      if (is.Surv(data[,i])) {
        estimator <- "tobit"
        S <- data[,i]
        y <- S[,1]
        if (attributes(S)$type=="left") 
          w <- S[,2]-1
        if (attributes(S)$type=="right") 
          w <- 1-S[,2]
        if (attributes(S)$type=="interval2") {
        w <- S[,3]; w[w==2] <- (-1)
      }
        mynames <- c(mynames,i)
        W <- cbind(W,w)
        data[,i] <- y
      }
    }
    if (length(W)>0) {
      colnames(W) <- mynames
      if (!is.null(weight)) {
        wW <- intersect(colnames(weight),colnames(W))
        if (length(wW)>0)
          weight[,wW] <- W[,wW]
        Wo <- setdiff(colnames(W),wW)
        if (length(Wo)>0)
        weight <- cbind(weight,W[,Wo,drop=FALSE])
      } else {
        weight <- W;
      }
    }
  }
  return(c(list(x=x,data=data,weight=weight,estimator=estimator),dots)) 
}
