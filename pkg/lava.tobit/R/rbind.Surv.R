rbind.Surv <- function(...) 
{
  dots <- list(...)
  type <- attributes(dots[[1]])$type
  ncol <- dim(dots[[1]])[2]
  nrow <- unlist(lapply(dots,nrow))
  cnrow <- c(0,cumsum(nrow))
  M <- matrix(ncol=ncol,nrow=sum(nrow))
  for (i in 1:length(dots)) {
    M[(cnrow[i]+1):cnrow[i+1],] <- dots[[i]]
  }
  x <- c(); for (i in 1:ncol(M)) x <- c(x,list(M[,i]))
  x <- c(x,list(type=type))
  do.call("Surv",x)
} 
