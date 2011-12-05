###{{{ plotcr

plotcr <- function(x,col,legend=TRUE,which=1:2,
                   ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                   ...) {
  dots <- list(...)
  if (is.null(dots$xlab)) dots$xlab <- "Time"
  if ((!is.data.frame(x) | !is.matrix(x)) && ncol(x)<2) stop("Wrong type of data")
  if (ncol(x)==2) {
    if (is.null(dots$curvlab)) {
      causes <- setdiff(unique(x[,2]),0)
      dots$curvlab <- seq(length(causes))
    }
    co <- cuminc(x[,1],x[,2])
    if (any(x[,1]<0)) {
      for (i in seq(length(co))) {
        co[[i]]$time <- co[[i]]$time[-1]
        co[[i]]$est <- co[[i]]$est[-1]
      }
      if (is.null(dots$xlim)) dots$xlim <- range(x[,1])
    }
    do.call("plot", c(list(x=co), dots))
    return(invisible(co))
  }

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  t <- x[,1:2]; cause <- x[,3:4]
  causes <- sort(setdiff(unique(cause),0))
  if (is.null(dots$curvlab))
    dots$curvlab <- seq(1:length(unique(causes)))
  if (missing(col)) col <- c("seagreen","darkred","darkblue","goldenrod","mediumpurple")
  if (1%in%which) {
    plot(t,type="n",...)
    count <- 1
    for (i in causes) {
      points(t[cause[,1]==causes[i],],col=Col(col[count],0.5),pch=2)
      points(t[cause[,2]==causes[i],],col=Col(col[count],0.5),pch=6)
      count <- count+1
    }
    points(t[cause[,1]==0 & cause[,2]==0,],col=Col("black",0.2),pch=1)  
    if (legend)
    legend("topleft", c("Subj 1, Cause 1", "Subj 2, Cause 1",
                        "Subj 1, Cause 2", "Subj 2, Cause 2",
                        "Double Censoring"), pch=c(2,6,2,6,1), col=c(rep(col[1:length(causes)],each=2),"black"))
  }
  if (2%in%which) {
    co1 <- cuminc(t[,1],cause[,1])
    co2 <- cuminc(t[,2],cause[,2])
    if (any(t[,1]<0)) {
      for (i in seq(length(co1))) {
        co1[[i]]$time <- co1[[i]]$time[-1]
        co1[[i]]$est <- co1[[i]]$est[-1]
      }
      for (i in seq(length(co2))) {
        co2[[i]]$time <- co1[[i]]$time[-1]
        co2[[i]]$est <- co1[[i]]$est[-1]
      }
    }
    if (is.null(dots$xlim)) dots$xlim <- range(t)
    do.call("plot", c(list(x=co1), dots))
    for (i in seq(length(co2)))
      with(co2[[i]], lines(time,est,type="s",lty=i,...))
  }
}

###}}} plotcr
