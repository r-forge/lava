ipw <- function(formula,data,cluster,samecens=TRUE,obsonly=TRUE,weightname="w",
                pairs=TRUE) {
  require(rms)
  fit <- cph(formula,data=data,surv=TRUE,x=TRUE,y=TRUE)
  timevar <- as.character(terms(formula)[[2]][[2]])
  otimes <- data[,timevar]
  utimes <- sort(unique(otimes))
  delta <- min(diff(c(0,utimes)))/2 ## We want prediction just before event
  pr <- survest(fit,what="parallel",newdata=data,times=otimes-delta)
  if (samecens & !missing(cluster)) {
    noncens <- with(data,!eval(terms(formula)[[2]][[3]]))
    data[,weightname] <- pr
    data <- data[noncens,,drop=FALSE]
    data <- data[order(data[,cluster]),,drop=FALSE]
    id <-  table(data[,cluster])
    if (pairs) {
      gem <- data[,cluster]%in%(names(id)[id==2])
      id <- id[id==2]
      data <- data[gem,]
    }
    d0 <- subset(data,select=c(cluster,weightname))
    timevar <- paste("_",cluster,weightname,sep="")
    d0[,timevar] <- unlist(lapply(id,seq))
    W <- apply(reshape(d0,direction="wide",timevar=timevar,idvar=cluster),1,
               function(x) min(x,na.rm=TRUE))
    data[,weightname] <- 1/rep(W,id)
    return(data)
  }
  return(pr)
}
