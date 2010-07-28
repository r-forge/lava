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
  
  return(c(list(x=x,data=data,weight=weight,estimator=estimator),dots)) 
}
