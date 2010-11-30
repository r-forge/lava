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
  nblocks <- ncol(mypar)
  for (ii in 1:nblocks) {
    data0 <- data; data0[,binpos[-mypar[,ii]]] <- NA
    mydata <- rbind(mydata,data0)
  }
  suppressWarnings(e0 <- estimate(x,data=mydata,missing=TRUE,silent=silent,
              ...))
  S <- score(e0,indiv=TRUE)
  nd <- nrow(data)
  block1 <- which((1:nd)%in%(rownames(S)))
  blocks <- sapply(1:nblocks, function(x) 1:length(block1)+length(block1)*(x-1))
  Siid <- matrix(0,nrow=length(block1),ncol=ncol(S))
  for (j in 1:ncol(blocks)) {
    Siid <- Siid+S[blocks[,j],]
  }
  iI <- vcov(e0); J <- t(Siid)%*%(Siid)
  e0$iidscore <- Siid
  e0$blocks <- blocks
  e0$vcov <- iI%*%J%*%iI ## thetahat-theta0 :=(asymp) I^-1*S => var(thetahat) = iI*var(S)*iI 
  cc <- e0$coef; cc[,2] <- sqrt(diag(e0$vcov))
  cc[,3] <- cc[,1]/cc[,2]; cc[,4] <- 2*(1-pnorm(cc[,3]))
  e0$coef <- cc
  class(e0) <- c("clprobit",class(e0))
  return(e0)
}
score.clprobit <- function(x,indiv=FALSE,...) {
  if (!indiv)
    return(colSums(x$iidscore))
  x$iidscore
}
