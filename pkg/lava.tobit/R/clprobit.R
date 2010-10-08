clprobit <- function(x,data,k=2,type=c("nearest","all"),pairlist,silent=TRUE,
                     ...) {
  bin <- binary(x)
  binpos <- which(colnames(data)%in%bin)
  if (length(bin)<(k+1)) stop("No need for composite likelihood analysis.")
  if (missing(pairlist)) {
    if (type[1]=="all") {
      mypar <- combn(length(bin),k) ## all pairs (or multiplets), k=2: k*(k-1)/2
    } else {
      mypar <- sapply(0:(length(bin)-k), function(x) x+1:k)
    }
  } else {
    mypar <- pairlist
  }
  mydata <- c()
  for (ii in 1:ncol(mypar)) {
    data0 <- data; data0[,binpos[-mypar[,ii]]] <- NA
    mydata <- rbind(mydata,data0)
  }  
  e0 <- estimate(x,data=mydata,missing=TRUE,silent=silent,
              ...)
  S <- score(e0,indiv=TRUE)
  Siid <- S[[1]]; Siid[is.na(Siid)] <- 0
  for (i in 2:length(S)) {
    S0 <- S[[i]]; S0[is.na(S0)] <- 0
    Siid <- Siid+S0
  }
  iI <- vcov(e0); J <- t(Siid)%*%(Siid)
  e0$vcov <- iI%*%J%*%iI
  cc <- e0$coef; cc[,2] <- sqrt(diag(e0$vcov))
  cc[,3] <- cc[,1]/cc[,2]; cc[,4] <- 2*(1-pnorm(cc[,3]))
  e0$coef <- cc
  class(e0) <- c("clprobit",class(e0))
  return(e0)
}
