
### plackett simulation mini eks
########################################################
library(mets)
## {{{ 
library(numDeriv)
library(simecol)
n <- 100
x <- runif(n)
y <- pcu(x,alpha=10)
cor(x,y)
dats=list(time=c(x,y),id=rep(1:n,2),status=rep(1,2*n),x=runif(2*n))
margs <- coxph(Surv(time,status)~+prop(x),data=dats)
dats <- data.frame(dats)
dats <- dats[order(dats$id),]
datsl <- dats[1:4,]

## }}} 

#####################################################################################
##########################                     ######################################
#####################################################################################

library(mets)
library(simecol)
library(doMC)
registerDoMC()

## {{{ simulation functions set up 

onerun <- function(i,n=1000){ ## {{{ simple constant Clayton-Oakes model 
print(i) 
d <- subset(simClaytonOakes(n,2,0.5,0,stoptime=2,left=0.5),!truncated)
marg <- aalen(Surv(lefttime,time,status)~+1,data=d,n.sim=0,robust=0)
fit1<-two.stage(marg,data=d,clusters=d$cluster,var.link=1,theta=0.74)
fit2<-twostage(marg,data=d,clusters=d$cluster,model="clayton.oakes",theta=0.74) 
fit3<-twostage(marg,data=d,clusters=d$cluster,score.method="optimize",theta=2.4)

out <- c(fit1$theta,fit1$var.theta^.5,fit2$theta,fit2$var.theta^.5,fit3$theta,fit3$var.theta^.5)
names(out) <- c("log-t.s","t.s-SE","log-ts","ts-SE","log-plack","se-plack")
return(out)
} ## }}}
## log(2)=0.69
onerun(2,n=1000)

onerunpl <- function(i,pard=2.5,n=1000,trans=0,both=1,method="optimize",numD=1){ ## {{{
print(i)
### pard=2.5;n=1000;trans=0;both=TRUE; method="optimize"; numD=1; 
x <- runif(n)
y <- pcu(x,alpha=pard)
if (trans==1) {
yt <- y^.5
xt <- log(x)
}

d <- data.frame(time=c(x,y),id=rep(1:n,2),status=rep(1,2*n),x=runif(2*n))
d$cluster <- d$id

marg <- aalen(Surv(time,status)~+1,data=d,n.sim=0,robust=0)
###fit1  <- twostage(marg,data=d,score.method="fisher.scoring",clusters=d$id,numDeriv=1)
fit1  <- twostage(marg,data=d,score.method=method,clusters=d$id,numDeriv=numD)

if (both==1) { 
	fit2  <- twostage(marg,data=d,score.method=method,clusters=d$id,model="clayton.oakes",numDeriv=numD)
}

ud <- c(fit1$theta,fit1$var.theta^.5,fit1$robvar.theta^.5)
if (both) ud <- c(ud,fit2$theta,fit2$var.theta^.5,fit2$robvar.theta^.5)
return(ud)
} ## }}}
## log(2.5)=0.91
###for (i in 1:10) print( onerunpl(1,pard=1.5,n=2000,method="fisher.scoring"))
onerunpl(1,pard=1.5,n=2000,method="fisher.scoring")

onerunpc <- function(i,n=2000,pard=2,model="clayton.oakes",beta=0,stoptime=2,returndata=0,pc=1,ots=0, both=TRUE){ ## {{{
###n=2000;pard=2;model="clayton.oakes";beta=0;stoptime=2;returndata=0;pc=1;ots=0; both=TRUE
print(i)
if (model=="clayton.oakes")  {
d <- subset(simClaytonOakes(n,2,1/pard,beta,stoptime=stoptime,left=0),!truncated)
d$num <- rep(1:2,n) 
d$id <- d$cluster
} else { ## {{{ plackett model 
x <- runif(n)
y <- pcu(x,alpha=pard)
cor(x,y,method="spearman")
x <- 2*x
y <- 2*y
censt <- rbinom(n,1,0.9)+runif(2*n)
timet <- c(x,y)
status <- ifelse((timet<censt),1,0)
timet <- ifelse((timet<censt),timet,censt)
timett <- c(log(3+timet[1:n]),timet[(n+1):(2*n)]^.5 )
num <- c(rep(1,n),rep(2,n))

###d <- list(time=timet,id=rep(1:n,2),status=status,x=runif(2*n),num=c(rep(1,n),rep(2,n)))
d <- list(time=timet,timet=timett,id=rep(1:n,2),status=status,x=runif(2*n),num=num)
d <- data.frame(d)
d$cluster <- d$id
} ## }}} 

if (returndata==1) return(d) else {

if (pc==1) {
udp  <- piecewise.twostage(c(0,0.5,2),data=d,score.method="optimize",
		  id="cluster",timevar="time",status="status",num="num")
ud <-  summary(udp)

ud <- c(ud$estimates,ud$se,ud$cor,ud$se.cor)
names(ud) <- c(rep("pl",4),rep("pl-se",4),rep("sp-cor",4),rep("sp-cor-se",4))

if (both==TRUE)  {
udp2  <- piecewise.twostage(c(0,0.5,2),data=d,score.method="optimize",model="clayton.oakes",
		  id="cluster",timevar="time",status="status",num="num")
ud2 <-  summary(udp2)
ud2 <- c(ud2$estimates,ud2$se,ud2$cor,ud2$se.cor)
###names(ud2) <- c(rep("co",4),rep("se co",4))
names(ud2) <- c(rep("co",4),rep("co-se",4),rep("co-ken",4),rep("co-ken-se",4))
ud <- c(c(ud),c(ud2))
}

return(ud)
} else {
   margn <- aalen(Surv(time,status)~-1+factor(num),data=d,n.sim=0,robust=0)
   udp <- twostage(margn,data=d,score.method="optimize",clusters=d$id,model=model)
   if (ots==1) {
   udp2 <- two.stage(margn,data=d,clusters=d$id,var.link=1)
   ud <- c(udp$theta,udp$var.theta^.5,udp2$theta,udp2$var.theta^.5)
   names(ud) <- c("est","se","est","se"); 
   } else {
   ud <- c(udp$theta,udp$var.theta^.5)
   names(ud) <- c("est","se"); 
   }

return(ud)
}
}

} ## }}}
onerunpc(1,pard=0.7)

data.plack <- function(n,pard=2,censi=1)
{ ## {{{

x <- runif(n)
y <- pcu(x,pard)
if (censi==0) censt <- rep(1,n) else censt <- runif(2*n)*censi
timet <- c(x,y)
status <- ifelse((timet<censt),1,0)
timet <- ifelse((timet<censt),timet,censt)

d <- list(time=timet,id=rep(1:n,2),status=status,x=runif(2*n),num=c(rep(1,n),rep(2,n)))
d <- data.frame(d)
d$cluster <- d$id

return(d); 
} ## }}}

onerunpccc <- function(i,n=20000,kc=100,kk=100,pard=2,beta=0,model="clayton.oakes",returndata=0,pc=1){ ## {{{
###
	i <- 1;n=100;kc=5;kk=5;pard=2;beta=0;model="clayton.oakes";returndata=0;pc=1; returndata=2; 
print(i)
if (model=="clayton.oakes")  { ## {{{
d <- subset(simClaytonOakes(n,2,1/pard,beta,stoptime=2,left=0),!truncated)
num <-rep(1:2,n)  
od <- order(num,d$cluster)
d <- d[od,]
num <- c(rep(1,n),rep(2,n))
timet <- d$time
status <- d$status
id <- d$cluster 
d$id <- id ## }}}
} else { ## {{{ plackett model 
x <- runif(n)
y <- pcu(x,alpha=pard)
x <- 2*x
y <- 2*y
censt <- rbinom(n,1,0.9)+runif(2*n)
timet <- c(x,y)
status <- ifelse((timet<censt),1,0)
timet <- ifelse((timet<censt),timet,censt)
id <- rep(1:n,2)
num <- c(rep(1,n),rep(2,n))
} ## }}} 
tiln <- 1:n
###
stat1x <- sample(tiln[(status==1)[1:n] ],kc)
cens0 <- tiln[(status==0)[1:n]]
if (length(cens0)< kk) kk <- length(cens0); 
stat0x <- sample(cens0,kk)
###
id1 <- which(id %in% stat1x)
id0 <- which(id %in% stat0x)
###
###d <- list(time=timet,id=rep(1:n,2),status=status,x=runif(2*n),num=c(rep(1,n),rep(2,n)))
dd <- list(time=timet[c(id1,id0)],id=id[c(id1,id0)],status=status[c(id1,id0)],num=num[c(id1,id0)],
	  case=c(rep(1,2*kc),rep(0,2*kk)))
d <- data.frame(dd)
d$cluster <- dd$id

if (returndata==2) { 
d <- reshape(d[,c("time","status","id","num","case")],direction="wide",idvar="id",timevar="num")
}

if (returndata>=1) return(d) else {
margs <- aalen(Surv(time,status)~-1+factor(num),data=d)
plph<-twostage(margs,data=d,theta=0.1,detail=1,Nit=0,clusters=d$id,model="plackett",score.method="optimize")
ud <- c(plph$theta,plph$var.theta^.5)
return(ud)
}

} ## }}}

###d <- onerunpc(1,model="plackett",returndata=1)
###d <- onerunpc(1,model="plackett",pc=0,ots=1)
###d <- onerunpc(1,model="plackett",pc=1,ots=1)
###d <- onerunpccc(1,kc=400,kk=400,model="plackett",returndata=2)

###library(epitools)
###tt <- table(d$time.1<0.5,d$time.2<0.5)
###oddsratio(tt)
###
###gem <- c()
###for (i in 1:100) {
###print(i); 
###d <- onerunpccc(1,n=2000,kc=200,kk=200,model="plackett",returndata=2)
###tt <- table(d$time.1<1,d$time.2<1)
###gem <- c(gem,oddsratio(tt)$measure[2,1])
###}
###mean(gem)
###sd(gem)


## }}} 

########################################################################333

komud<-function(){  ## {{{ 


d <- onerunpccc(1,model="plack",kc=4,kk=10,returndata=1)

res <- c()
for (i in 1:100) {
print(i)
ud <- onerunpccc(1,model="plack",kc=300,kk=100,beta=0)
res <- rbind(res,ud)
}
apply(res,2,mean)
apply(res,2,sd)
log(2)

res <- c()
for (i in 1:100) {
print(i)
ud <- onerunpc(1,model="clayton.oakes",pc=0,ots=1)
res <- rbind(res,ud)
}
apply(res,2,mean)
apply(res,2,sd)
log(2)



n <- 1000

dd <- data.plack(10000,pard=5,censi=0)
dd2 <- data.frame(simple.reshape(dd,id="id"))
cor(dd2$time.1,dd2$time.2,method="spearman")
alpha2spear(log(5))

dd <- data.plack(10000,pard=5,censi=0)
dd2 <- data.frame(simple.reshape(dd,id="id"))
cor(dd2$time.1,dd2$time.2,method="spearman")
cor(dd2$time.1,dd2$time.2,method="kendall")
alpha2spear(log(5))
###
sa <- surv.boxarea(c(0.5,0.5),c(1,1),dd,id="id",timevar="time",status="status")
sa2 <- data.frame(simple.reshape(sa,id="id"))
cor(sa2$time.1,sa2$time.2,method="spearman")

d <- onerunpc(1,n=10000,model="plack",returndata=1)
sa <- surv.boxarea(c(0,0),c(0.5,0.5),d,id="cluster",timevar="time",status="status",num="num")
sa <- surv.boxarea(c(0.5,0.5),c(8,8),d,id="cluster",timevar="time",status="status")
sa <- surv.boxarea(c(0,0),c(0.5,0.5),d,id="cluster",timevar="time",status="status")
} ## }}}

##########################simulations   ##############################################333

komud<-function(){  ## {{{  for loops 

res <- c()
for (i in 1:500) {
ud <- onerunpl(i)
res <- rbind(res,ud)
}
pars <-  c(1,4)
separs <- c(2:3,5:6)
signif(apply(res[,pars,drop=F],2,mean),4)
signif(apply(res[,pars,drop=F],2,sd),4)
signif(apply(res[,separs],2,mean),4)
rep(apply(res[,pars],2,sd),each=2)/
apply(res[,separs],2,mean)


res <- c()
for (i in 1:500) {
ud <- onerunpc(i)
res <- rbind(res,ud)
}
colnames(res)
pars <-  c(1:4,9:12,17:20,25:28)
separs <- c(5:8,13:16,21:24,29:32)
signif(apply(res[,pars,drop=F],2,mean),4)
signif(apply(res[,pars,drop=F],2,sd),4)
signif(apply(res[,separs],2,mean),4)
signif(apply(res[,pars,drop=F],2,sd),4)/ signif(apply(res[,separs],2,mean),4)

} ## }}}

stop()


library(doMC)
registerDoMC()

nsim <- 1000
resm <- c()
res <- foreach (i=0:nsim) %dopar% onerunpl(i,n=2000,pard=2,both=1,method="fisher.scoring",numD=1)
for (i in 1:length(res)) resm <- rbind(resm,c(res[[i]]))
###pars <-  c(1)
###separs <- c(2:3)
pars <-  c(1,4)
separs <- c(2:3,5:6)
signif(apply(resm[,pars,drop=F],2,mean),4)
signif(apply(resm[,pars,drop=F],2,sd),4)
signif(apply(resm[,separs],2,mean),4)
rep(apply(resm[,pars],2,sd),each=2)/signif(apply(resm[,separs],2,mean),4)


nsim <- 1000
resm <- c()
res <- foreach (i=0:nsim) %dopar% onerunpc(i,n=10000,pard=2,model="plackett",both=TRUE)
###onerunpc(i,n=1000,pard=1,model="clayton.oakes");
for (i in 1:length(res)) resm <- rbind(resm,c(res[[i]]))
colnames(resm)
pars <-  c(1:4,9:12,17:20,25:28)
colnames(resm)[pars]
separs <- c(5:8,13:16,21:24,29:32)
colnames(resm)[separs]
signif(apply(resm[,pars],2,mean),4)
signif(apply(resm[,pars],2,sd),4)
signif(apply(resm[,separs],2,mean),4)
signif(apply(resm[,pars],2,sd),4)/signif(apply(resm[,separs],2,mean),4)

resm10000 <- resm
resm5000 <- resm

signif(apply(resm1000,2,mean),4)
apply(resm1000,2,sd)

signif(apply(resm10000,2,mean),4)
apply(resm10000,2,sd)



onerunpc(i,n=10000,pard=10,model="plackett")

nsim <- 100
resm <- c()
res <- foreach (i=0:nsim) %dopar% onerunpccc(i,kc=100,kk=200,pard=2.7,model="plackett")
for (i in 1:length(res)) resm <- rbind(resm,c(res[[i]]))
signif(apply(resm,2,mean),4)
apply(resm,2,sd)

nsim <- 100
resm <- c()
res <- foreach (i=0:nsim) %dopar% onerunpccc(i,kc=100,kk=200,pard=2,model="clayton.oakes")
for (i in 1:length(res)) resm <- rbind(resm,c(res[[i]]))
signif(apply(resm,2,mean),4)
apply(resm,2,sd)




install.packages("simecol")

ddco <- onerunpc(10000,returndata=1)

dd <- onerunpc(10000,returndata=1,model="pc")

names(dd)

table(dd$status)
table(ddco$status)

table(dd$cluster,dd$status)
ddd <- reshape(dd,direction="wide",idvar="id",varying=list(c("status")),timevar="num")

dd <- data.plack(100000,censi=100)
ddd <- reshape(dd,direction="wide",idvar="id",v.names=c("status","time"),timevar="num")
table(ddd$status.1,ddd$status.2)

ud <- twoby2( table(ddd$time.1>0.5,ddd$time.2>0.5))

library(mets)
library(simecol)
library(Epi)
ud <- twoby2( table(ddd$status.1,ddd$status.2))
ud$measures[2,1]


mantelhaen.test( table(ddd$id,ddd$status.1,ddd$status.2))




nsim <-  100
res <- foreach (i=0:nsim) %dopar% onerun(i,n=2000)
resm <- c()
for (i in 1:length(res)) resm <- rbind(resm,c(res[[i]]))
signif(apply(resm,2,mean),4)
apply(resm,2,sd)[c(1,3,5)]/apply(resm,2,mean)[c(2,4,6)]

nsim <-  500
res <- foreach (i=0:nsim) %dopar% onerunpl(i,n=2000)
###resm <- c()
log(2.5)
for (i in 1:length(res)) resm <- rbind(resm,c(res[[i]]))
signif(apply(resm,2,mean),4)
apply(resm,2,sd)

ud=cluster.index(c(1,2,1,2,1,2))

nsim <- 200 
res <- foreach (i=0:nsim) %dopar% onerunpc(i,n=2000)
###resm <- c()
for (i in 1:length(res)) resm <- rbind(resm,c(res[[i]]))
signif(apply(resm[,c(1:4,9:12)],2,mean),4)-c(rep(log(2),4),rep(0.5,4))
signif(apply(resm[,c(1:4,9:12)],2,sd),4)/ signif(apply(resm[,c(5:8,13:16)],2,mean),4) 

nsim <- 100 
res <- foreach (i=0:nsim) %dopar% onerunpc(i,n=2000,model="plackett")
resmp <- c()
for (i in 1:length(res)) resmp <- rbind(resmp,c(res[[i]]))
signif(apply(resmp[,c(1:4,9:12)],2,mean),4)
signif(apply(resmp[,c(1:4,9:12)],2,sd),4)/ signif(apply(resmp[,c(5:8,13:16)],2,mean),4) 


library(mets)

data(diabetes)

head(diabetes)
ud1=surv.boxarea(c(0,30),c(30,100),data=diabetes,id="id",timevar="time",status="status")
oi <- order(ud1$id)
ud1 <- ud1[oi,]
head(ud1)
table(ud1$num,ud1$lef)

library(simecol)
n <- 5000
x <- runif(n)
y <- pcu(x,2.5)
cor(x,y,method="spearman")
alpha2spear(2.5)
###
d <- list(time=c(x,y),id=rep(1:n,2),status=rep(1,2*n),x=runif(2*n))
d <- data.frame(d)
margs <- coxph(Surv(time,status)~+prop(x),data=d)
plph<-twostage(margs,data=d,theta=0.1,detail=1,Nit=0,clusters=d$id, model="plackett",score.method="nlminb")
summary(plph)
log(2.5)
#

udp  <- piecewise.twostage(c(0,0.5,1),data=d,score.method="optimize",id="id",timevar="time",status="status", model=model,silent=0)
summary(udp)



library(mets)
library(simecol)
n <- 5000
x <- runif(n)
y <- pcu(x,alpha=3)
cor(x,y,method="spearman")
alpha2spear(3,link=0)
d <- list(time=c(x,y),id=rep(1:n,2),status=rep(1,2*n),x=runif(2*n))
margs <- aalen(Surv(time,status)~+1,data=d,robust=0,n.sim=0)
plph<-twostage(margs,data=d,theta=log(3),detail=1,Nit=0,clusters=d$id, model="plackett",score.method="nlminb")
summary(plph)
log(3)

plph<-twostage(margs,data=d,theta=log(3),detail=1,Nit=20,clusters=d$id,model="plackett",
	       score.method="fisher.scoring",iid=1,step=0.2)
names(plph)


#

Cuvt <- function(x)
{
	u <- x[1]
	v <- x[2]
	theta <- x[3]
	S <- 1+(theta-1)*(u+v)
	cc <- (S - ( S^2 -4 *u*v*theta*(theta-1))^.5)/(2*(theta-1))
return(cc)
}

Cuv <- function(x)
{
	u <- x[1]
	v <- x[2]
	S <- 1+(theta-1)*(u+v)
	cc <- (S - ( S^2 -4 *u*v*theta*(theta-1))^.5)/(2*(theta-1))
return(cc)
}

theta <- 2

### conditional OR T> 0.5,S> 0.3

T<-runif(1)*0.5+0.5
S<-runif(1)*0.7+0.3
c(T,S)
psiti <-Cuv(c(0.5,0.3))-Cuv(c(T,0.3)) -Cuv(c(0.5,S))+Cuv(c(T,S)) 
ptsi <- Cuv(c(T,0.3))- Cuv(c(T,S))
psti <- Cuv(c(0.5,S))- Cuv(c(T,S))
c(S,Cuv(c(0.5,S))/Cuv(c(0.5,0.3)))
c(T,Cuv(c(T,0.3))/Cuv(c(0.5,0.3)))
c(psiti,ptsi,psti)
(Cuv(c(T,S))* psiti) / (ptsi*psti)
Cuv(c(0.5,0.3))
###
Cuv(c(T,S))*(1-T-S+ Cuv(c(T,S)))/((T- Cuv(c(T,S)))*(S- Cuv(c(T,S))))


library(numDeriv)

x <- c(0.9,0.5)
Cuv(x)

jacobian(Cuv,x)

hessian(Cuv,x)[1,2] * Cuv(x) / prod(jacobian(Cuv,x))

2*log(2)

hessian(Cuv,x)


xt <- c(0.9,0.5,2)
(hessian(Cuvt,xt+c(0,0,0.01))[1,2] - hessian(Cuvt,xt)[1,2] ) * jacobian(Cuvt,xt)[3] / (100 * prod( hessian(Cuvt,xt)[1:2,3]))


#########################################################

library(mets)
library(simecol)
library(numDeriv)

onerunpc <- function(i,n=2000,pard=2,model="clayton.oakes",beta=0,returndata=0,pc=1,ots=0){ ## {{{
print(i)
if (model=="clayton.oakes")  {
d <- subset(simClaytonOakes(n,2,1/pard,beta,stoptime=2,left=0),!truncated)
d$num <- rep(1:2,n) 
d$id <- d$cluster
} else { ## {{{ plackett model 
xv <- rbinom(n,1,0.5); 
yv <- rbinom(n,1,0.5); 
x <- runif(n)
y <- pcu(x,alpha=pard)
x <- -log(1-x)/exp(xv*beta)
y <- -log(1-y)/exp(yv*beta)
censt <- rbinom(n,1,0.9)+runif(2*n)
timet <- c(x,y)
status <- ifelse((timet<censt),1,0)
time <- ifelse((timet<censt),timet,censt)
###timett <- c(log(3+timet[1:n]),timet[(n+1):(2*n)]^.5 )
num <- c(rep(1,n),rep(2,n))

###d <- list(time=timet,id=rep(1:n,2),status=status,x=runif(2*n),num=c(rep(1,n),rep(2,n)))
d <- list(time=time,timet=timet,id=rep(1:n,2),status=status,x1=c(xv,yv),num=num)
d <- data.frame(d)
d$cluster <- d$id
} ## }}} 

if (returndata==1) return(d) else {

if (pc==1) {
udp  <- piecewise.twostage(c(0,0.5,2),data=d,score.method="optimize",id="cluster",timevar="time",status="status",num="num",model=model)
ud <-  summary(udp)
ud <- c(ud$estimates,ud$se,ud$cor,ud$se.cor)
names(ud) <- c(rep("cr",4),rep("se cr",4),rep("cor",4),rep("se cor",4))
return(ud)
} else {
   margn <- coxph(Surv(time,status)~-1+x1+strata(num),data=d)
###   margn <- cox.aalen(Surv(time,status)~num+prop(x1),data=d,n.sim=0,max.clust=NULL,max.timepoint.sim=10,silent=1)
###   margn <- aalen(Surv(time,status)~-1+factor(num),data=d,n.sim=0,robust=0)
   udp <- twostage(margn,data=d,score.method="optimize",clusters=d$id,model=model,theta=pard)
   if (ots==1) {
   udp2 <- two.stage(margn,data=d,clusters=d$id,var.link=1,theta=pard)
   ud <- c(udp$theta,udp$var.theta^.5,udp2$theta,udp2$var.theta^.5)
   names(ud) <- c("est","se","est","se"); 
   } else {
   ud <- c(udp$theta,udp$var.theta^.5)
   names(ud) <- c("est","se"); 
   }
return(ud)
}
}

} ## }}}

onerunpc(20,n=20000,pard=5,model="plackett",pc=1)
log(5)

	onerunpc(1,n=400,pard=2,model="plackett",beta=0.2,returndata=0,pc=1,ots=0)
ud <- onerunpc(1,n=30,pard=2,model="plackett",beta=0.2,returndata=1,pc=1,ots=0)
udp  <- piecewise.twostage(c(0,0.5,2),data=ud,score.method="optimize",id="cluster",timevar="time",status="status",
			   num="num",model="plackett",data.return=1,silent=-1)

cud <- cluster.index(ud$id)
udw <- cbind(ud[cud$idclust[,1]+1,], ud[cud$idclust[,2]+1,])
names(udw) <- c(paste(names(ud),".1",sep=""), paste(names(ud),".2",sep=""))
head(udw)

ud <- onerunpc(1,n=100,pard=2,model="plackett",beta=0.2,returndata=1,pc=1,ots=0)
udw <- fast.reshape(ud,id=ud$id)

loglikef <- function(par,ud=NULL,model="clayton.oakes")
{
###	par <- c(0.2,2)
   beta <- par[1]
   theta <- par[2]
   margn <- cox.aalen(Surv(time,status)~-1+prop(x1)+strata(num),data=ud,beta.fixed=1,beta=beta)
   udp <- twostage(margn,data=ud,score.method="fisher.scoring",clusters=ud$id,model=model,theta=theta,Nit=0)
   udp$score
###   udp$score1
###   udp$hess
###   udp$loglike
}

n <-100
ud <- onerunpc(1,n=n,pard=2,model=model,beta=0.2,returndata=1,pc=0,ots=0)


onerun <- function(i,n) {
print(i)
model <- "plackett"
ud <- onerunpc(1,n=n,pard=2,model=model,beta=0,returndata=1,pc=0,ots=0)
###hessian(loglikef,c(0.2,2),model=model)
jj <- jacobian(loglikef,c(0.0,2),model=model,ud=ud)
###
###model <- "clayton.oakes"
###ud <- onerunpc(1,n=n,pard=2,model=model,beta=0.2,returndata=1,pc=0,ots=0)
###hessian(loglikef,c(0.2,2))
###jjc <- jacobian(loglikef,c(0.2,2),ud=ud)
c(c(jj),c(jj))
}

library(doMC)
registerDoMC()

res <- list()
res <- foreach (i=0:100) %dopar% onerun(i,n=1000)
###resmp <- c()
for (i in 1:length(res)) resmp <- rbind(resmp,c(res[[i]]))
apply(resmp,2,mean)
apply(resmp,2,sd)/nrow(resmp)^.5

apply(resmp,2,mean)
apply(resmp,2,median)

res <- c()
for (i in 1:40) {
print(i)
model <- "plackett"
ud <- onerunpc(1,n=n,pard=2,model=model,beta=0.2,returndata=1,pc=0,ots=0)
###hessian(loglikef,c(0.2,2),model=model)
jj <- jacobian(loglikef,c(0.2,2),model=model)
###
###model <- "clayton.oakes"
###ud <- onerunpc(1,n=n,pard=2,model=model,beta=0.2,returndata=1,pc=0,ots=0)
###jjc <- jacobian(loglikef,c(0.2,2))
udc <- c(c(jj),c(jjc))
res <- rbind(res,udc)
}
apply(res,2,mean)
apply(res,2,sd)/nrow(res)^.5

dthetaloglike <- function(par) { ## {{{
x <- par[1]; y <- par[2]; z <- par[3]; 
###y <- exp(log(y)*exp(beta))
###z <- exp(log(z)*exp(beta))
### d/dx , d/dx d/dy, d/dx d/dz, d/dx d/dz d/dy
vec <- 
c( 
(y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*pow(-1 + x,2)),
  (1 + ((-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)),
  (1 + ((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)),
((-3*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(8.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),2.5)) + ((-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((2*pow(-1 + x,2) - 4*(-1 + x)*x)*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + (2*x)/Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*(-1 + x)) - (((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2))
  )
dutheta <- vec
dutheta
} ## }}}

ud <- onerunpc(1,n=30,pard=2,model="plackett",beta=0.2,returndata=1,pc=1,ots=0)

wolloglike <- function(par) { ## {{{
      x <- par[1]; y <- par[2]; z <- par[3]; beta <- par[4]
###  f,  d/dy f, d/dz f, d/dz d/dy f
y <- exp(log(y)*exp(beta))
z <- exp(log(z)*exp(beta))
vek <- c( (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ,
(-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)),
(-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) ,
(((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) )

return(vek); 
} ## }}}

wolloglike <- function(par) { ## {{{
 x <- par[1]; y <- par[2]; z <- par[3]; b <- par[4]
###  f,  d/dy f, d/dz f, d/dz d/dy f
y <- exp(log(y)*exp(b))
z <- exp(log(z)*exp(b))
vek <- c( 
 (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ,
(-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)),
(-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) ,
(((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) )


return(vek); 
} ## }}}

llp <- wolloglike(c(1.1,0.134387,0.893903,1.3))
ll <- dloglike(c(0.1,0.134387,0.893903,1.3))
ll
jj <- jacobian(dloglike,c(0.1,0.134387,0.893903,1.3))
jj
sum(jj[,4]*llp)


jj <- jacobian(dthetaloglike,c(1.1,0.134387,0.893903,1.3))
sum(jj[,4]*llp)



###udw <- reshape(ud[,c("id","time","x1","status","num")],
###	       direction="wide",idvar="id",varying=list(c("time","x1","status")),timevar="num") 


ddd <- reshape(dd,direction="wide",idvar="id",varying=list(c("status")),timevar="num")
	
	gem <- c()
	for (i in 1:500)
	gem <- rbind(gem, onerunpc(i,n=1000,pard=1,model="plackett",beta=0.2,returndata=0,pc=1,ots=0))
	mm <- apply(gem,2,mean)
	msd <- apply(gem,2,sd)
	mm[5:8]/msd[1:4]

	gemco <- c()
	for (i in 1:200)
	gemco <- rbind(gemco, onerunpc(i,n=1000,pard=2,model="clayton.oakes",beta=0.2,returndata=0,pc=1,ots=0))
	mm <- apply(gemco,2,mean)
	msd <- apply(gemco,2,sd)
	mm[5:8]/msd[1:4]


onerunpc(i,n=400,pard=2,model="clayton.oakes",beta=0.2,returndata=0,pc=0,ots=1)
onerunpc(i,n=400,pard=2,model="clayton.oakes",beta=0.2,returndata=0,pc=0,ots=1)


#########################################################################
###### simulation af betinget spearman                     ##############
#########################################################################
library(colorout)

oddsratioWald.proc <- function(tab,n00, n01, n10, n11, alpha = 0.05){## {{{ 	
n00 <- tab[1] 
n01 <- tab[2] 
n10 <- tab[3] 
n11 <- tab[4] 
	  #
	  #  Compute the odds ratio between two binary variables, x and y,
	  #  as defined by the four numbers nij:
	  #
	  #    n00 = number of cases where x = 0 and y = 0
	  #    n01 = number of cases where x = 0 and y = 1
	  #    n10 = number of cases where x = 1 and y = 0
	  #    n11 = number of cases where x = 1 and y = 1
	  #
	  OR <- (n00 * n11)/(n01 * n10)
  #
  #  Compute the Wald confidence intervals:
  #
  siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
    zalph <- qnorm(1 - alpha/2)
    logOR <- log(OR)
      loglo <- logOR - zalph * siglog
      loghi <- logOR + zalph * siglog
        #
        ORlo <- exp(loglo)
        ORhi <- exp(loghi)
	  #
	  oframe <- data.frame(LowerCI = ORlo, OR = OR, UpperCI = ORhi, alpha = alpha)
	  oframe
} ## }}} 

sgn <- function(x,y)
{
 (x<0)*(y<0) + (x>0)*(y>0) - (x>0)*(y<0) - (x<0)*(y>0)  
}

library(simecol)
or=3
rhoor <- alpha2rho(3)
n <- 200000
x <- runif(n)
y <- pcu(x,rho=rhoor)
cor(x,y,method="spearman")
rhoor
###
x1 <- runif(n)
y1 <- runif(n)
###
dat <- cbind(x,y,x1,y1)
minx <- apply(dat[,c(1,3)],1,min)
miny <- apply(dat[,c(2,4)],1,min)
###
s <- sgn((dat[,1]-dat[,3]),(dat[,2]-dat[,4]))
3*mean(sgn((dat[,1]-dat[,3]),(dat[,2]-dat[,4])))
datm <- dat[minx<0.3 & miny<0.7,]
datms <- dat[minx>0.3 & miny>0.7,]
###
cor(datm[,1],datm[,2],method="spearman")
3*mean(sgn((datm[,1]-datm[,3]),(datm[,2]-datm[,4])))

tt=table(dat[,1]<runif(1),dat[,2]<runif(1))
oddsratioWald.proc(c(tt))

ttm=table(datm[,1]<0.3,datm[,2]<0.7)
oddsratioWald.proc(c(ttm))
###
ttms=table(datms[,1]<0.3,datms[,2]<0.7)
oddsratioWald.proc(c(ttm))

censi1 <- runif(n)
censi2 <- runif(n)
datmi <- dat[minx<censi1 & miny<censi2,]
datmis <- dat[minx>censi1 & miny>censi2,]
censi11 <- censi1[minx<censi1 & miny<censi2]
censi21 <- censi1[minx<censi1 & miny<censi2]
ttm=table(datmi[,1]<censi11,datmi[,2]<censi21)
oddsratioWald.proc(c(ttm))
###
ttms=table(datmis[,1]<censi1,datmis[,2]<censi2)
oddsratioWald.proc(c(ttm))




