ipw <- function(formula,data,cluster,samecens=TRUE,obsonly=TRUE,weight="w",
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
    data[,weight] <- pr
    data <- data[noncens,,drop=FALSE]
    data <- data[order(data[,cluster]),,drop=FALSE]
    id <-  table(data[,cluster])
    if (pairs) {
      gem <- data[,cluster]%in%(names(id)[id==2])
      data <- data[gem,]
    }
    id <- data[,cluster]
    for (i in unique(id)) {
      data[id==i,weight] <- min(data[id==i,weight])
    }
    data[,weight] <- 1/data[,weight]
    return(data)
  }
  return(pr)
}
