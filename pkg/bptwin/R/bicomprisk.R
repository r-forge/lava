bicomprisk <- function(formula, data, cause=c(1,1), cens=0, causes, indiv, strata=NULL, id,num,sym=FALSE,ssym=FALSE,prodlim=FALSE,messages=TRUE,model,
                       return.data=0,...) {
  mycall <- match.call()
  formulaId <- Specials(formula,"id")
  formulaIndiv <- Specials(formula,"indiv")
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- paste("~.-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")",
                     "-indiv(",paste(formulaIndiv,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  indiv <- formulaIndiv
  if (!is.null(formulaIndiv)) {
    mycall$indiv <- indiv
  } 
  if (missing(id)) stop("Missing 'id' variable")
  
  timevar <- terms(formula)[[2]]
  if (is.call(timevar)) {
    causes <- timevar[[3]]
    timevar <- timevar[[2]]
  }  
  timevar <- as.character(timevar)
  causes <- as.character(causes)
  
  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq_len(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        mycall$data <- dd[[i]]
        eval(mycall)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("bicomprisk.strata","biprobit.strata")
      res$N <- length(dd)
      return(res)
    }
  }
  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  indiv2 <- covars2 <- NULL 
  
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  which(data[,id]%in%names(idtab==2))
  data <- data[which(data[,id]%in%names(idtab==2)),]
  if (missing(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,num] <- unlist(lapply(idtab,seq_len))
  }

  timevar2 <- paste(timevar,1:2,sep=".")
  causes2 <- paste(causes,1:2,sep=".")
  if (length(covars)>0)
    covars2 <- paste(covars,1,sep=".")
  for (i in seq_len(length(indiv)))
    indiv2 <- c(indiv2, paste(indiv[i],1:2,sep="."))
  
  ww0 <- reshape(data[,c(timevar,causes,covars,indiv,id,num)],
                 direction="wide",idvar=id,timevar=num)[,c(timevar2,causes2,indiv2,covars2)]

  if (sym & cause[1]!=cause[2]) {
    switchpos <- 1:4
    if (length(indiv)>0)
       switchpos <- c(switchpos, seq_len(2*length(indiv))+4)
     newpos <- switchpos
     for (i in seq_len(length(newpos)))   
       newpos[i] <- newpos[i] + ifelse(i%%2==1,1,-1)
    ww1 <- ww0[,newpos]; colnames(ww1) <- colnames(ww0)
    ww0 <- rbind(ww0,ww1)
    ## suppressMessages(browser())  
    ## switchers <- which(ww0[,causes2[1]]==cause[2] & ww0[,causes2[2]]==cause[1])
    ## switchpos <- 1:4
    ## if (length(indiv)>0)
    ##   switchpos <- c(switchpos, seq_len(2*length(indiv))+4)
    ## newpos <- switchpos
    ## for (i in seq_len(length(newpos)))   
    ##   newpos[i] <- newpos[i] + ifelse(i%%2==1,1,-1)
    ## ww0[switchers,switchpos] <- ww0[switchers,newpos]
  }
  ww0 <- na.omit(ww0)
 
  status <- rep(0,nrow(ww0))
  time <- ww0[,timevar2[2]]
  mycauses <- setdiff(unique(data[,causes]),0)

  time <- status <- rep(0,nrow(ww0))
  ##  time <- ww0[,"time.1"]

  ## cause = (i,j)
  primcond <- ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cause[2]
  if (ssym) primcond <- primcond | (ww0[,causes2[1]]==cause[2] & ww0[,causes2[2]]==cause[1]) ## (j,i)
  idx2 <- which(primcond)
  status[idx2] <- 1
  time[idx2] <- apply(ww0[idx2,timevar2[1:2]],1,max)
  ##  suppressMessages(browser())  
  
  ##(0,0), (0,j)
  cond <- ww0[,causes2[1]]==cens & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2])
  if (ssym) cond <- ww0[,causes2[1]]==cens & (ww0[,causes2[2]]%in%c(cens,cause)) ## (0,i)
  idx2 <- which(cond)
  status[idx2] <- 0
  time[idx2] <- ww0[idx2,timevar2[1]]

  ##(i,0)
  cond <- ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cens
  if (ssym) cond <- (ww0[,causes2[1]]%in%cause & ww0[,causes2[2]]==cens) ## (j,0)
  idx2 <- which(cond)
  status[idx2] <- 0
  time[idx2] <- ww0[idx2,timevar2[2]]
  
  ##(ic,jc)
  cond <- ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2])
  if (ssym) {
    cond <- (cond | ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[2] & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[1])) & !primcond
  }
  idx2 <- which(cond)
  status[idx2] <- 2
  time[idx2] <- apply(ww0[idx2,timevar2[1:2]],1,min)
  
  ##(ic,0), (ic,j) **
  cond <- ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2])
  if (ssym) cond <- (cond | (ww0[,causes2[1]]!=cens & ww0[,causes2[1]]==cause[2] & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[1]))) & !primcond  ##(jc,0), (jc,i)
  idx2 <- which(cond)
  status[idx2] <- 2
  time[idx2] <- ww0[idx2,timevar2[1]]

  ##(0,jc),(i,jc)
  cond <- (ww0[,causes2[1]]==cens | ww0[,causes2[1]]==cause[1]) & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2])
  if (ssym) cond <- (cond | (ww0[,causes2[1]]==cens | ww0[,causes2[1]]==cause[2]) & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[1])) & !primcond  ##(0,ic),(j,ic)
  idx2 <- which(cond)
  status[idx2] <- 2
  time[idx2] <- ww0[idx2,timevar2[2]]
  
  mydata0 <- mydata <- data.frame(time,status,ww0[,covars2],ww0[,indiv2])
  names(mydata) <- c(timevar,causes,covars,indiv2)


  
  if (!prodlim) {
    ff <- paste("Surv(",timevar,",",causes,"!=",cens,") ~ 1",sep="")
    if (length(c(covars,indiv))>0) {
      xx <- c(covars,indiv2)
      for (i in seq_len(length(xx)))
        xx[i] <- paste("const(",xx[i],")",sep="")
      ff <- paste(c(ff,xx),collapse="+")
      if (missing(model)) model <- "fg"      
    }
    if (missing(model)) model <- "additive"
    add<-comp.risk(as.formula(ff),data=mydata,
                   status,causeS=1,n.sim=0,resample.iid=1,model=model) ## time=100
    padd <- predict(add,X=1,se=1,uniform=1)
  } else {
    ff <- as.formula(paste("Hist(",timevar,",",causes,")~",paste(c("1",covars,indiv2),collapse="+")))
    padd <- prodlim(ff,data=mydata)
  }
  class(padd) <- c("bicomprisk",class(padd))
  if (return.data==1) return(list(comp.risk=padd,data=mydata)) else return(padd)}


plot.bicomprisk <- function(x,add=FALSE,...) {
  if ("predict.timereg"%in%class(x)) {    
    if (!add) { plot.predict.timereg(x,...) }
    else {
      with(x,lines(time,P1,...))
    }
    
  } else {
    plot(x,...)
  }
  return(invisible(x))
}


conc2case <- function(conc,marg) {
  out <- conc
  time <- conc$time
  margtime <- Cpred(cbind(marg$time,c(marg$P1)),time)[,2]
  out$P1 <- conc$P1/margtime
  if (!is.null(out$se.P1)) out$se.P1 <- conc$se.P1/margtime
  attr(out,"class") <- rev(attr(out,"class")) 
  return(out)
}

back2timereg <- function(obj){ 
  out <- obj
  attr(out,"class") <- rev(attr(out,"class")) 
  return(out)
}
