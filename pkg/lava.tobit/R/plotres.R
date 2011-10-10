plotres <- function(x,var=endogenous(x),
                    ylab="Cumulative Distribution Function",
                    xlab="Standardized residuals",
                    main,
                    ...) {
  require(survival)
  r <- residuals(x,std=TRUE)
  W <- Weight(x)
  
  for (v in var) {    
    if (v %in% colnames(W)) {
      S <- Surv(ifelse(W[,v]==-1,NA,r[,v]),
                ifelse(W[,v]==1,NA,r[,v]),
                type="interval2")      
    } else {
      S <- Surv(r[,v],rep(TRUE,length(r[,v])))
    }
    g <- survfit(S~1)
    mymain <- ifelse(!missing(main),main,v)    
    with(g,plot(1-surv~time,type="s",main=mymain,xlab=xlab,ylab=ylab))
    with(g,lines(1-upper~time,type="s",lty=2))
    with(g,lines(1-lower~time,type="s",lty=2))
    ro <- sort(r[,v]); 
    lines(ro,pnorm(ro),col="red",xlab=xlab,ylab=ylab)
    
   }
  invisible(x)
}
