cumh <- function(formula,data,...,time,timestrata=seq(55,100,length=10),
               silent=FALSE) { 
  res <- list(); i <- 0
  ht <- c()
  outcome <- as.character(terms(formula)[[2]])
  y0 <- data[,outcome]
  for (i in seq(length(timestrata))) {
    t <- timestrata[i]
    data[,outcome] <- y0
    data[,outcome] <- data[,outcome]*(data[,time]<t)
    if (!silent) {
      message(t," ",sum(data[,outcome]))
    }
    res[[i]] <-bptwin(formula,data=data,...)
    ht <- rbind(ht,c(t,summary(res[[i]])$h))
  }
  rownames(ht) <- timestrata
  colnames(ht) <- c("time","Heritability","Std.Err","2.5%","97.5%")
  res <- (list(ht=ht,models=res))
  class(res) <- "cumh"
  res
}

summary.cumh <- function(object,...) object 
  
print.cumh <- function(x,...) {
  print(x$ht)
  invisible(x)
}

Col <- function (col, alpha = 0.2) {
    sapply(col, function(x) do.call(rgb, as.list(c(col2rgb(x)/255, 
        alpha))))
}

plot.cumh <- function(x,...,idx=seq(nrow(x$ht)),lwd=2,col,ylim=c(0,1),xlab="Time",ylab="Heritability",add=add) {

  if (!add) {
    plot(x$ht[idx,1:2,drop=FALSE],type="l",ylim=ylim,lwd=lwd,
         ylab=ylab,xlab=xlab,...)
  }
  if (missing(col)) col <- Col("darkblue")
  xx <- with(x, c(ht[idx,1],rev(ht[idx,1])))
  yy <- with(x, c(ht[idx,4],rev(ht[idx,5])))           
  polygon(xx,yy,col=col)
  lines(x$ht[idx,1:2,drop=FALSE],lwd=lwd,...)
  invisible(x)
}
