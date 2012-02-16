
do.biprobit.strata <- function(x,fun,...) {
  res <- lapply(x$model,function(m) do.call(fun,c(list(m),list(...))))
  names(res) <- names(x$model)
  class(res) <- "do.biprobit.strata"
  return(res)
}
print.do.biprobit.strata <- function(x,...) {
  for (i in seq_len(length(x))) {    
    message(rep("-",60),sep="")
    message("Strata '",names(x)[i],"'",sep="")
    print(x[[i]])
  }
  return(invisible(x))
}

plot.biprobit.strata <- function(x,...)
  suppressMessages(do.biprobit.strata(x,"plot",...))
print.biprobit.strata <- function(x,...)
  print.do.biprobit.strata(x$model,...)
summary.biprobit.strata <- function(object,...)
  do.biprobit.strata(object,"summary",...)
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
score.biprobit.strata <- function(x,...) {
  ss <- lapply(x$model,function(m) score(m,indiv=FALSE,...))
  return(unlist(ss))
}
