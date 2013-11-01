library(mets)
 data(prt);
###
ss <- sample(unique(prt$id),100)
prtid <-which( prt$id %in% ss)
prt <- prt[prtid,]
sum(table(prt$id)==1)

 ### marginal cumulative incidence of prostate cancer
library(mets)
data(prt)
times <- seq(60,100,by=10)
outm <- comp.risk(Surv(time,status==0)~+1,data=prt,
		  prt$status,causeS=2,max.clust=100)
###
 cifmz <- predict(outm,X=1,uniform=0,resample.iid=1)
 cifdz <- predict(outm,X=1,uniform=0,resample.iid=1)
 dim(cifmz$P1.iid)
 dim(cifdz$P1.iid)

 cc <- bicomprisk(Hist(time,status)~+1+id(id),data=prt,cause=c(2,2),max.clust=100)
 dim(cc$P1.iid)

 ### concordance for MZ and DZ twins
 cc <- bicomprisk(Hist(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2),max.clust=100)
 cdz <- cc$model$"DZ"
 cmz <- cc$model$"MZ"
 dim(cmz$P1.iid)
 dim(cdz$P1.iid)

 dim(cifmz$P1.iid)
 dim(cifdz$P1.iid)

 conc <- cmz
 marg <- cifmz
 casedz <- casewise.test(cdz,cifdz)
 casemz <- casewise.test(cmz,cifmz)
 casemz$casewise

 casemz <- casewise.test(cdz,cifmz,test="conc")
 casemz <- casewise.test(cmz,cifmz,test="case")
 par(mfrow=c(1,2))
 plot(casemz)
 plot(casedz)

 cifmziid <- predict(outm,X=1,uniform=0,resample.iid=1)
 conc <- cmz
 marg <- cifmziid
 casemziid <- casewise.test(cmz,cifmziid)
 casemziid.conc <- casewise.test(cmz,cifmziid,test="conc")
 casemziid.case <- casewise.test(cmz,cifmziid,test="case")

 par(mfrow=c(1,2))
 plot(casemz)
 plot(casedz)

 plot(casemz,ylim=c(0,0.5),xlim=c(60,100))
 par(new=TRUE)
 plot(casedz,ylim=c(0,0.5),xlim=c(60,100),col=3)

 ################# prodlim #####################################3
 ################################################################

 outm <- prodlim(Hist(time,status)~+1,data=prt)
      
 times <- seq(60,100,by=5)
      pcifmz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="MZ"))
      pcifdz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="DZ"))
           
 ### concordance for MZ and DZ twins
 pcc <- bicomprisk(Hist(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2),prodlim=TRUE)
 pcdz <- pcc$model$"DZ"
 pcmz <- pcc$model$"MZ"
		     
 pcdz <- casewise(pcdz,outm,cause.marg=2)
 pcmz <- casewise(pcmz,outm,cause.marg=2)
		          
 plot(pcmz,ci=NULL,ylim=c(0,0.5),xlim=c(60,100),legend=TRUE,col=c(3,2,1))
       par(new=TRUE)
     plot(pcdz,ci=NULL,ylim=c(0,0.5),xlim=c(60,100),legend=TRUE)

  summary(pcdz)
  summary(pcmz)
 cbind(casemz$casewise, casemz$se.casewise)


 plot(pcmz)
 lines(casemz$casewise)
 summary(pcmz)
 lines(casemz$se.casewise)


 ###############################################################3
 ###############################################################3
 ###############################################################3

library(mets)
data(prt)
###ss <- sample(unique(prt$id),100)
###prtid <-which( prt$id %in% ss)
###prt <- prt[prtid,]
###sum(table(prt$id)==1)
###
 ### marginal cumulative incidence of prostate cancer
times <- seq(60,100,by=1)
###times <- c(50,90)
outm <- comp.risk(Surv(time,status==0)~+1+cluster(id),data=prt,
		  prt$status,causeS=2,max.clust=100,conservative=1)
names(outm)
###
###
cifmz <- predict(outm,X=1,uniform=0,resample.iid=1)
cifdz <- predict(outm,X=1,uniform=0,resample.iid=1)
cifmz$se.P1
dim(cifmz$P1.iid)
head(cifmz$P1.iid)
###cifmz$clusters
###length(cifmz$clusters)
###length(unique(cifmz$clusters))

 ### concordance for MZ and DZ twins
head(prt)
 cc <- bicomprisk(Hist(time,status)~+1+id(id),data=prt,cause=c(2,2),
		  se.clusters=outm$clusters,conservative=1)
dim(cc$P1.iid)
plot(cc,ylim=c(0,0.05),se=1)
cc$se.P1
ccs <- bicomprisk(Hist(time,status)~+strata(zyg)+id(id),data=prt,cause=c(2,2),
		  se.clusters=outm$clusters,conservative=1)
### strata(zyg) går galt pga noget med data og strata
#############################################################
ccmz <- ccs$model$"MZ"
ccdz <- ccs$model$"DZ"
cc$P1
cc$time
dim(ccmz$P1.iid)
head(ccs$P1.iid)
length(unique(cc$clusters))

nn <- colnames(cc$P1.iid)
ii <- cifmz$P1.iid
dim(ii[,nn]) 
dim(cc$P1.iid)
cifmz$se.P1

ud <- casewise.test(ccmz,cifmz,test="conc")
ud$casewise

uddz <- casewise.test(ccdz,cifmz,test="case")
uddz$casewise

 dim(cmz$P1.iid)
 dim(cifmz$P1.iid)

 plot(ud)

 plot(ud)
 lines(ud$casewise[,1],ud$casewise[,2]-1.96*ud$casewise[,3],col=2)
 lines(ud$casewise[,1],ud$casewise[,2]+1.96*ud$casewise[,3],col=2)

################################################################################################
################################################################################################
################################################################################################


library(mets)
library(prodlim)

## {{{ simulation for gamma distributed cif model 

lap<-function(theta,t) { return( (1+t/theta)^(-theta)) }
ilap<-function(theta,t) {
itheta<-1/theta; return((t^(-itheta)-1)/(itheta)) }

F1clust<-function(t,rtheta=1,theta=1,lam0=0.5,beta=0.3,x=0) {
return(1-exp(-rtheta*ilap(theta,exp(-t*lam0-t*x*beta))))
}
F1clust(1); 
F1clust(0.1)
tt <- seq(0,3,by=0.01)


###################################################################
F1<-function(t,lam0=0.5,beta=0.3,x=0) # additive version
{ return( 1 - exp(-t*lam0-t*x*beta)) }

###F1<-function(t,lam0=0.5,beta=0.3,x=0) # proportional version 
###{ return( 1 - exp(-(t*lam0)*exp(x*beta))) }

sim.F1<-function(n,theta=1,lam0=0.5,beta=0.3,crate=2) 
{ ## {{{
x<-runif(n); tt<-seq(0,1,length=100)
F11x<-F1(1,x=x,beta=beta,lam0=lam0)
cause1<-rbinom(n,1,F11x)

stime<-rep(100,n); 
for (i in 1:n)
{
if (cause1[i]==1) {
myhazx<-F1(tt,x=x[i],beta=beta,lam0=lam0)/F11x[i]
stime[i]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
} 
}
ctime<-runif(n)*crate
time<-apply(cbind(ctime,stime),1,min)
status<-(stime<ctime); cause1[status==0]<-0; 
data<-data.frame(time=time,status=status,X=x,cause=cause1)
return(data)
} ## }}}

sim.F1clust<-function(n,theta=1,lam0=0.5,beta=0.3,alpha=0,crate=2,same.cens=FALSE,fix.cens=FALSE)
{ ## {{{ 
k<-n/2; tt<-seq(0,1,length=100)
rtheta<-rgamma(k,theta,scale=1/theta)
stime<-c();cause1<-c();id<-c();vtheta<-c(); X<-c()
cause1 <- c()

for (i in 1:k)
{ 
x<-runif(2); 
X<-c(X,x); 
F11x<-F1clust(1,rtheta=rtheta[i],theta=theta,x=x,beta=beta,lam0=lam0) 
cause<-rbinom(2,1,F11x); 
###cause1<-c(cause1,cause); 
id<-c(id,rep(i,2)); vtheta<-c(vtheta,rep(rtheta[i],2))

for (j in 1:2) {
if (cause[j]==1) {
myhazx<-F1clust(tt,x=x[j],rtheta=rtheta[i],theta=theta, beta=beta,lam0=lam0)/F11x[j]
stime<-c(stime,Cpred(cbind(myhazx,tt),runif(1))[1,2]+ runif(1,0,0.001))
cause1 <- c(cause1,1);
} else { stime<-c(stime,runif(1)); cause1 <- c(cause1,2);}
}
}

if (same.cens) ctime <- rep(ctime<-runif(k)*crate,each=2) else ctime<-runif(n)*crate;
if (fix.cens) ctime <- rep(ctime<-rep(crate,each=2))

time<-apply(cbind(ctime,stime),1,min)
status<-(stime<ctime); 
cause1[status==0]<-0; 
data<-data.frame(time=time,X=X,cause=cause1,stime=stime,ctime=ctime,
id=id,theta=vtheta)
return(data)
} ## }}}

###F1(tt)
## }}}

######################################################
### test ting  #######################################
######################################################

beta=0; n=200; theta=0.1; 
simdata<-sim.F1clust(n,theta=theta,beta=beta,crate=3.0,lam0=0.5)
table(simdata$cause)
simdata$tv <- rep(1:2,n/2)
###
tt <- seq(0,1,by=0.01)
F1(tt)
###
p11<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(1/theta, 2*ilap(1/theta,1-F1(tt)) )
###plot(tt,F1(tt),type="l")
###lines(tt,p11)
###lines(tt,(F1(tt)^2),col=2)

###beta=0; n=1000; theta=0.1; 
###gem <- c()
###for (i in 1:10) {
###simdata<-sim.F1clust(n,theta=theta,beta=beta,crate=3.0,lam0=0.5,same.cens=TRUE)
###table(simdata$cause)
###simdata$tv <- rep(1:2,n/2)
###simd <- reshape(simdata,direction="wide",idvar="id",timevar="tv")
###simd$cancer <- (simd$cause.1==1 ) + (simd$cause.2==1 )  
###tabd <- table(simd$cancer)
###casewise <- tabd["2"]/ ( tabd["2"] + tabd["1"]*0.5)  
###gem <- c(gem,casewise)
###}
###summary(gem)

### plot

casewiset <- p11/F1(tt)
plot(tt,casewiset,type="l")
F1t <- F1(tt)

###################################################################
###################################################################

### increase censoring by crate lower, dependence theta lower  

theta=1.0;  
p11t<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(theta, 2*ilap(theta,1-F1(tt)) )
F1t <- F1(tt)
Gc <- 1 - 2*tt
Gc[Gc<0] <- 0
mn11t <- cumsum(diff(p11t)* Gc[-1])
mndt <- cumsum(diff(F1t)* Gc[-1])
casewisem <- mn11t[100]/mndt[100]
casewiset <- p11t/F1(tt)

gem <- c()
for (i in 1:10) {
print(i)
beta=0; n=10000; 
theta=theta;  
simdata<-sim.F1clust(n,theta=theta,beta=beta,crate=0.5,lam0=0.5,same.cens=TRUE,fix.cens=FALSE)
table(simdata$cause)
simdata$tv <- rep(1:2,n/2)
simd <- reshape(simdata,direction="wide",idvar="id",timevar="tv")
simd$cancer <- (simd$cause.1==1 ) + (simd$cause.2==1 )  
tabd <- table(simd$cancer)
casewise <- tabd["2"]/ ( tabd["2"] + tabd["1"]*0.5)  
###
outm <- prodlim(Hist(time,cause)~+1,data=simdata)
par(mfrow=c(1,3))
plot(outm)
lines(tt,F1t,col=3,lwd=3)
legend("topleft",c("cif","sand"),lty=1,col=c(1,3))
###          
cc <- bicomprisk(Hist(time,cause)~+1+id(id),data=simdata,cause=c(1,1),prodlim=TRUE)
p11t<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(theta, 2*ilap(theta,1-F1(tt)) )
plot(cc)
lines(tt,p11t,col=3,lwd=3)
lines(tt,F1t^2,col=4,lwd=1)
legend("topleft",c("conc","sand","indep"),lty=1,col=c(1,3,4))
###
case <- casewise(cc,outm,cause.marg=1)
###
casewiset <- p11t/F1t
plot(case$casewise[,1],case$casewise[,2],type="l")
lines(tt,casewiset,col=3,lwd=3)
abline(h=casewise,col=2)
legend("topleft",c("est","sand","constant"),lty=1,col=c(1,3,2))
###
Gc <- 1 - 2*tt
Gc[Gc<0] <- 0
###
abline(h=casewisem,col=4)
gem <- c(gem,casewise)
}

abline(h=casewisem,col=4)
casewisem
mean(gem)

pdf("caswise-vs-surv.pdf")
par(mfrow=c(1,1))

dev.off()

############################################################################
############################################################################
############################################################################

library(doMC)
registerDoMC(45)
library(mets)

onesim <- function(i,n,theta=2,beta=0) { ## {{{ 
###	theta=2; beta=0; n=4000
	print(i)
simdata<-sim.F1clust(n,theta=theta,beta=beta,crate=0.5,lam0=0.5,same.cens=TRUE,fix.cens=FALSE)
simdata$tv <- rep(1:2,n/2)
simd <- reshape(simdata,direction="wide",idvar="id",timevar="tv")
simd$cancer <- (simd$cause.1==1 ) + (simd$cause.2==1 )  
tabd <- table(simd$cancer)
casewise <- tabd["2"]/ ( tabd["2"] + tabd["1"]*0.5)  
###
outm <- prodlim(Hist(time,cause)~+1,data=simdata)
###par(mfrow=c(1,3))
###plot(outm)
###lines(tt,F1t,col=3,lwd=3)
###legend("topleft",c("cif","sand"),lty=1,col=c(1,3))
###          
cc <- bicomprisk(Hist(time,cause)~+1+id(id),data=simdata,cause=c(1,1),prodlim=TRUE)
p11t<- 1 - (1-F1(tt)) - (1-F1(tt)) + lap(theta, 2*ilap(theta,1-F1(tt)) )
case <- casewise(cc,outm,cause.marg=1)
###case$concordance
###case$casewise
###simdata$time
times <- seq(0,0.4,by=0.01)
casepl <- Cpred(case$casewise,times)
concpl <- Cpred(case$concordance,times)

 outm <-comp.risk(Surv(time,cause==0)~+1+cluster(id),data=simdata,simdata$cause,
	                         causeS=1,times=times,n.sim=0,max.clust=NULL)
 cifmz <-predict(outm,X=1,uniform=0,resample.iid=1)
 cccr <-bicomprisk(Hist(time,cause)~+1+id(id),data=simdata,
                          cause=c(1,1),se.clusters=outm$clusters)
 casecr <- casewise.test(cccr,cifmz,test="case") ## test based on casewise
 casecr$test

ud <- list(prodlim=list(casewise=casepl,concordance=concpl),comprisk=casecr)
} ## }}} 
ud <- onesim(10,500)
###ud$prodlim$concordance
###ud$prodlim$case
###ud$comprisk$conc
###ud$comprisk$case

res <- foreach (i=1:1000) %dopar% onesim(i,1000)
gemcasepl <- gemconcpl <- gemcasecr <- gemconccr <- c()
gemcaseplse <- gemconcplse <- gemcasecrse <- gemconccrse <- c()
gemcasecrse2 <- gemconccrse2 <- c()
resm <- c()
for (i in 1:length(res)) {
	resm <- rbind(resm,c(res[[i]]$comprisk$test[4]))
	gemcasepl <- cbind(gemcasepl,res[[i]]$prodlim$case[,2])
	gemcaseplse <- cbind(gemcaseplse,res[[i]]$prodlim$case[,3])
	gemconcpl <- cbind(gemconcpl,res[[i]]$prodlim$concordance[,2])
	gemconcplse <- cbind(gemconcplse,res[[i]]$prodlim$concordance[,3])
	gemcasecr <- cbind(gemcasecr,res[[i]]$comprisk$case[,2])
	gemcasecrse <- cbind(gemcasecrse,res[[i]]$comprisk$case[,3])
	gemcasecrse2 <- cbind(gemcasecrse2,res[[i]]$comprisk$case[,4])
	gemconccr <- cbind(gemconccr,res[[i]]$comprisk$conc[,2])
	gemconccrse <- cbind(gemconccrse,res[[i]]$comprisk$conc[,3])
}
apply(gemcasepl,1,sd)/apply(gemcaseplse,1,mean)
apply(gemconcpl,1,sd)/apply(gemconcplse,1,mean)
apply(gemcasecr,1,sd)/apply(gemcasecrse,1,mean)
apply(gemcasecr,1,sd)/apply(gemcasecrse2,1,mean)
apply(gemconccr,1,sd)/apply(gemconccrse,1,mean)
mean(c(resm)<0.05)

