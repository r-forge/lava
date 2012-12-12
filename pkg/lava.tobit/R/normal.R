##' @export
normal_logLik.lvm <- function(object,p,data,data2,...) {
  res <- -normal_objective.lvm(x=object,p=p,data=data,data2=data2,...)
  return(res)
}

normal_method.lvm <- "nlminb0"

normal_gradient.lvm <- function(x,p,data,data2,...) {
  grad(function(p0) normal_objective.lvm(x,p=p0,data=data,data2=data2,...),p)
}

normal_objective.lvm <- function(x,p,data,data2,...) {
  set.seed(1)
  y.idx <- index(x)$endo.idx
  mu <- predict(x,data=data,p=p)
  S <- attributes(mu)$cond.var
  class(mu) <- "matrix"
  y <- endogenous(x)
  ord <- ordinal.lvm(x)
  status <- rep(0,length(y))
  status[match(binary(x),y)] <- 2
  status[match(ord,y)] <- 2  
  thres <- matrix(0,nrow=length(y),max(1,attributes(ord)$K-1)); rownames(thres) <- y
  for (i in seq_len(length(attributes(ord)$fix))) {
    nn <- names(attributes(ord)$idx)[i]
    ii <- attributes(ord)$idx[[nn]]
    val <- (attributes(mu)$e[ii])
    thres[nn,seq_len(length(val))] <-
      cumsum(c(val[1],exp(val[-1])))
  }  
  yl <- yu <- as.matrix(data[,y,drop=FALSE])
  if (!is.null(data2)) {
    yu[,colnames(data2)] <- data2
    status[match(colnames(data2),y)] <- 1
  }
##  suppressMessages(browser())
  l <- .Call("loglikMVN",
             yl=yl,
             yu=yu,
             status=status,
             mu=as.matrix(mu),dmu=NULL,s=S,ds=NULL,
             z=NULL,su=NULL,dsu=NULL,
             threshold=thres,
             dthreshold=NULL)

  return(-sum(l))  
}
