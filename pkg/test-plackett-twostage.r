
library(mets)
###example(twostage)
###example(binomial.twostage)

n=10
d <- subset(simClaytonOakes(n,2,0.5,0,stoptime=2,left=0),!truncated)
d$tv <- rep(1:2,n)
dd <- fast.reshape(d,id="cluster")
###
marg1 <- aalen(Surv(time,status)~+1,data=d,n.sim=0,robust=0)

cm1 <- coxph(Surv(time,status)~+x1,data=d)
cm1$linear.predictor
margsurv <- cm1

coxformula <- margsurv$formula
X <- model.matrix(coxformula,data=d)[,-1];
coef(margsurv)
 baseout <- basehaz(margsurv,centered=FALSE);
 baseout <- cbind(baseout$time,baseout$hazard)
cumh <-  Cpred(baseout,d$time)
 RR<-exp(X %*% coef(margsurv))
 psx<-exp(-baseout*RR)


cm1 <- phreg(Surv(time,status)~+1,data=d)

cm1$X







n=1000
d <- subset(simClaytonOakes(n,2,0.5,0,stoptime=2,left=0),!truncated)
d$tv <- rep(1:2,n)
dd <- fast.reshape(d,id="cluster")
###
marg1 <- aalen(Surv(time,status)~+1,data=d,n.sim=0,robust=0)
fitp1<-twostage(marg1,data=d,
	clusters=d$cluster,
	score.method="optimize",
###	model="plackett",
	model="clayton.oakes",
	detail=1,step=1.0,numDeriv=1)
###fitp1<-twostage(marg1,data=d,clusters=d$cluster,
###	score.method="optimize",model="clayton.oakes",detail=0,step=0.5,numDeriv=1)
fitp1$loglike
fitp1$score
fitp1$score1
fitp1$hess
fitp1$hessi
fitp1$Dscore
summary(fitp1)

fitp1<-twostage(marg1,data=d,
	clusters=d$cluster,
	score.method="optimize",
	model="plackett",
###	model="clayton.oakes",
	detail=1,step=1.0,numDeriv=1)
fitp1$hess
fitp1$Dscore
summary(fitp1)

fitp1<-twostage(marg1,data=d,
	clusters=d$cluster,
	score.method="fisher.scoring",
###	model="plackett",
	model="clayton.oakes",
	detail=1,step=1.0,numDeriv=1)
fitp1$hess
fitp1$Dscore
summary(fitp1)

fitp1<-twostage(marg1,data=d,
	clusters=d$cluster,
	score.method="nlminb",
###	model="plackett",
	model="clayton.oakes",
	detail=1,step=1.0,numDeriv=1)
fitp1$hess
fitp1$Dscore
summary(fitp1)


fitp1<-twostage(marg1,data=d,clusters=d$cluster,
score.method="optimize",model="plackett",detail=0,
step=0.5,numDeriv=1)
summary(fitp1)
fitp1$score
fitp1$score1
fitp1$hess
fitp1$hessi
fitp1$Dscore

fitp1<-twostage(marg1,data=d,clusters=d$cluster,
score.method="fisher.scoring",model="plackett",detail=0,
step=0.5,numDeriv=1,Nit=50,var.link=0)
summary(fitp1)
fitp1$score
fitp1$score1
fitp1$hess
fitp1$hessi
fitp1$Dscore


summary(fitp1)
fitp1$score
names(fitp1)
head(fitp1$theta.iid)

udss <- data.frame(fast.reshape(d,id="cluster"))
table(udss$status1,udss$status2)

uds <- piecewise.data(c(0,0.5,2),data=d,id="cluster",timevar="time",
status="status")
tdes <- model.matrix(~~-1+factor(strata),data=uds)
ud1 <- uds
marg1 <- aalen(Surv(boxtime,status)~-1+factor(strata):factor(num),data=ud1,n.sim=0,robust=0)
fitp1<-twostage(marg1,data=ud1,clusters=ud1$cluster,
		model="clayton.oakes",detail=0,step=0.5,
score.method="fisher.scoring",theta.des=tdes,strata=ud1$strata)
summary(fitp1)
fitp1$score

table(uds$strata,uds$status)
uds1 <- subset(uds,strata=="0.5-2,0.5-2")
udss <- fast.reshape(uds1,id="cluster")
table(udss$status1,udss$status2)


ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,
		 id="cluster",timevar="time",status="status",num="tv")
udss <- data.frame(fast.reshape(ud4,id="cluster",num="tv"))
table(udss$status1,udss$status2)
data=d; id="cluster";timevar="time";status="status"
num="tv"; covars=NULL; covars.pairs=NULL
left.trunc=c(0.5,0.5)
right.cens=c(2,2)
ww0 <- data.frame(ww0)
table(ww0$status1,ww0$status2)


udp <- piecewise.twostage(c(0,0.5,2),data=d,
	  score.method="fisher.scoring",id="cluster",timevar="time",
status="status",
###model="clayton.oakes",
silent=1,data.return=1)
summary(udp)

### Same model using the strata option, a bit slower
########################################################

ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")
ud1$strata <- 1; ud2$strata <- 2; ud3$strata <- 3; ud4$strata <- 4
ud <- rbind(ud1,ud2,ud3,ud4)

marg1 <- aalen(Surv(boxtime,status)~-1+factor(num),data=ud1,n.sim=0,robust=0)
fitp1<-twostage(marg1,data=ud1,clusters=ud1$cluster,score.method="optimize",model="clayton.oakes",
			   detail=0,step=0.5)
summary(fitp1)
fitp1$score
names(fitp1)

marg2 <- aalen(Surv(boxtime,status)~-1+factor(num):factor(strata),data=ud,n.sim=0,robust=0)
tdes <- model.matrix(~-1+factor(strata),data=ud)
fitp2<-twostage(marg2,data=ud,clusters=ud$cluster,score.method="fisher.scoring",model="clayton.oakes",
			   theta.des=tdes,step=1.0,detail=0,strata=ud$strata)
summary(fitp2)

### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
ud$stratas <- ud$strata; ud$stratas[ud$strata==3] <- 2;
tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
fitp3<-twostage(marg2,data=ud,clusters=ud$cluster,score.method="fisher.scoring",model="clayton.oakes",
				     theta.des=tdes2,step=1.0,detail=0,strata=ud$strata)
summary(fitp3)
     
     ### could also symmetry in marginal models








data(diabetes)
# Marginal Cox model  with treat as covariate
marg <- cox.aalen(Surv(time,status)~prop(treat)+cluster(id),data=diabetes,
		  robust=0,n.sim=0,max.clust=NULL)
fit<-two.stage(marg,data=diabetes,theta=1.0,detail=0,Nit=40)
summary(fit)
###
margph <- coxph(Surv(time,status)~treat,data=diabetes)
fitpho<-two.stage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,var.link=1)
summary(fitpho)

###source("plackettMLE.R")
margph <- coxph(Surv(time,status)~treat,data=diabetes)
fitph<-twostage(margph,data=diabetes,clusters=diabetes$id,theta=1.,
		iid=1,var.link=0,model="clayton.oakes")
coef(fitph)
c(fitph$score, fitph$score1)
fitph$Dscore
fitph$hess
names(fitph)
head(fitph$theta.iid)

diab <- diabetes[-1,]
margph2 <- coxph(Surv(time,status)~treat,data=diab)
fitph2<-twostage(margph2,data=diab,clusters=diab$id,theta=1,
	var.link=0,model="clayton.oakes")
coef(fitph2)
fitph2$Dscore
fitph2$hess
head(fitph2$theta.iid)


library(mets)
data(diabetes)
diabetes$idr <- diabetes$id
pfitph2<-piecewise.twostage(c(0,30),c(0,30),data=diabetes,silent=-1,
		   id="id", time="time",status="status", covars="idr")

pcdata<-piecewise.data(c(0,30,50),c(0,30,50),data=diabetes,silent=-1,
		    id="id", time="time",status="status", covars="idr")


stop()


fitphf<-twostage(margph,data=diabetes,clusters=diabetes$id,detail=1,theta=0.11,Nit=0,score.method="fisher.scoring",iid=1,var.link=0)

fitphf<-twostage(margph,data=diabetes,clusters=diabetes$id,theta=2.1,detail=1,Nit=0,score.method="fisher.scoring",iid=1,var.link=0,model="clayton.oakes")
c(fitphf$score, fitphf$score1)
fitph$Dscore
fitph$hess

fitphf<-twostage(margph,data=diabetes,clusters=diabetes$id,theta=0.1,detail=1,Nit=0,score.method="fisher.scoring",iid=1,var.link=0,model="clayton.oakes")
c(fitphf$score, fitphf$score1)
fitph$Dscore
fitph$hess

fitphf<-twostage(margph,data=diabetes,clusters=diabetes$id,theta=2.1,detail=1,Nit=0,score.method="fisher.scoring",iid=1,var.link=1,model="clayton.oakes")
c(fitphf$score, fitphf$score1)
fitph$Dscore
fitph$hess

fitphf<-twostage(margph,data=diabetes,clusters=diabetes$id,theta=0.1,detail=1,Nit=0,score.method="fisher.scoring",iid=1,var.link=1,model="clayton.oakes")
c(fitphf$score, fitphf$score1)
fitph$Dscore
fitph$hess



source("plackettMLE.R")
fitphf<-twostage(margph,data=diabetes,clusters=diabetes$id,theta=1.1,detail=1,Nit=40,score.method="fisher.scoring",iid=1,var.link=0)
fitph<-twostage(margph,data=diabetes,clusters=diabetes$id,theta=1.1,detail=1,Nit=0,score.method="nlminb",iid=1,var.link=0)
c(fitph$score, fitph$score1)
c(fitphf$score, fitphf$score1)
fitph$Dscore
fitph$hess

head(fitph$theta.iid)
head(fitph$margsurv)
head(diabetes)

wollike(c(exp(1),0.72,0.48))
dloglike(c(exp(1),0.72,0.48))

dloglike(c(exp(1),0.72,0.48))/wollike(c(exp(1),0.72,0.48))

dloglike(c(exp(1),0.73,0.57))/wollike(c(exp(1),0.73,0.57))




library(mets); 
###
d <- subset(simClaytonOakes(2000,2,2,1,stoptime=2,left=2),!truncated)
###
e <- ClaytonOakes(Surv(lefttime,time,status)~x1+cluster(~1,cluster),data=d)
summary(e)
###m <- coxph(Surv(lefttime,time,status)~+1,data=d)
m <- aalen(Surv(lefttime,time,status)~+1,data=d,n.sim=0,robust=0,max.clust=NULL)
summary(m)
class(m)
em <-two.stage(m,data=d,detail=0,clusters=d$cluster)
summary(em)

source("plackettMLE.R")
###
em <-twostage(m,data=d,detail=1,clusters=d$cluster,
	     model="plackett",score.method="nlminb",control=list(trace=TRUE),var.link=0)
summary(em)

em <-twostage(m,data=d,detail=1,clusters=d$cluster,theta=0.4,
	     model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=0)
summary(em)

em <-twostage(m,data=d,detail=1,clusters=d$cluster,
	     model="plackett",score.method="optimize",control=list(trace=TRUE),var.link=0)
summary(em)
names(em)
em$score
em$score1
em$Dscore
em$hess

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

source("plackettMLE.R")
plph<-twostage(margs,data=dats,detail=1,Nit=0,clusters=dats$id, model="plackett",score.method="nlminb",iid=1)
plph$theta
plph$score
plph$score1
plph$Dscore
plph$hess
plph$hessi

coef(plph)
summary(plph)
print(plph)
head(plph$theta.iid)
sum(plph$theta.iid)

source("plackettMLE.R")
datsl <- dats[1:4,]
plph1<-twostage(margs,data=datsl,theta=0.1,detail=1,Nit=0,clusters=datsl$id,model="plackett",score.method="fisher.scoring",var.link=1)
plph1<-twostage(margs,data=datsl,theta=0.1,detail=1,Nit=0,clusters=datsl$id,model="plackett",score.method="nlminb",var.link=1)
plph1$theta
c(plph1$score, plph1$score1)
c(plph1$Dscore, plph1$hess, plph$hessi)


plph1<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id,model="plackett",score.method="fisher.scoring",var.link=1)

summary(plph1)
plph1$theta
plph1$score
plph1$score1
plph1$Dscore
plph1$hess



###plph<-twostage(margph,data=diabetes,theta=0.5,detail=1,Nit=40,clusters=diabetes$id,
###	       control=list(trace=TRUE),model="plackett")
###plph$theta
###plph$score
###plph$hessi
###plph$var.theta
###coef(plph)

library(mets); 
library(simecol)
n <- 500
x <- runif(n)
y <- pcu(x,alpha=10)
cor(x,y,method="spearman")
alpha2spear(0.9,link=0)

alpha2spear(0.99)
source("plackettMLE.R")

dats=list(time=c(x,y),id=rep(1:n,2),status=rep(1,2*n),x=runif(2*n))
margs <- coxph(Surv(time,status)~+prop(x),data=dats)
plph<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id, model="plackett",score.method="nlminb")
plph$theta
coef(plph)
summary(plph)
print(plph)
###
plph1<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id,control=list(trace=TRUE),model="plackett",score.method="fisher.scoring",var.link=1)
summary(plph1)
plph1$score
###
plph0<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id,control=list(trace=TRUE),model="plackett",score.method="fisher.scoring",var.link=0)
###		score.method="nlminb",var.link=0)
plph0$theta
plph0$score
coef(plph0)
summary(plph0)
summary(plph1)

plack.cif2(0.5,0.5,-9.1)

plack2spear <- function(alpha) (alpha^2 - 1 - 2 * alpha * log(alpha))/(alpha - 1)^2
plack2spear(plph$theta)

plph2<-twostage(margs,data=dats,theta=0.1,detail=0,Nit=0,clusters=dats$id,control=list(trace=TRUE),model="clayton.oakes",score.method="nlminb")
plph2$theta
coef(plph2)
plph2ne<-twostage(margs,data=dats,theta=0.1,detail=0,Nit=0,clusters=dats$id,control=list(trace=TRUE),model="clayton.oakes",score.method="nlminb",var.link=0)
coef(plph2ne)

tph<-two.stage(margs,data=dats,theta=0.1,detail=0,Nit=40,clusters=dats$id)
tph$theta


## {{{ clayton oakes derivatives and 2 nd derivatives vs NumDeriv

sqr <- function(x) x^2
Sqrt <- function(x) sqrt(x) 
Exp <- function(x) exp(x)
pow <- function(x,p) x^p
Ln <- function(x) log(x)

co <- function(par) {
theta <- par[1]; S1 <- par[2]; S2 <- par[3]; 
x <- par[1]; y <- par[2]; z<- par[3]; 
###out <- ((S1^{-(1/theta)}+S2^{-(1/theta)})-1)^(-(theta))
Exp(-Ln(Exp(Ln(y)*(-1/x))+ Exp(Ln(z)*(-1/x))-1 )*x)
}
co(c(2,0.9,0.9))

###D(y) Exp(-Ln(Exp(Ln(y)*(-1/x))+ Exp(Ln(z)*(-1/x))-1 )*x)
###D(z) Exp(-Ln(Exp(Ln(y)*(-1/x))+ Exp(Ln(z)*(-1/x))-1 )*x)
###D(y) (Exp(-Ln(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1)*x)*x*Exp((-Ln(z))/x))/(z*x*(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1))

dco <- function(par)
{
x <- par[1]; y <- par[2]; z<- par[3]; 
(Exp(-Ln(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1)*x)*x*Exp((-Ln(z))/x))/(z*x*(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1))
}

d2co <- function(par)
{
x <- par[1]; y <- par[2]; z<- par[3]; 
((z*x*(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1)*Exp((-Ln(z))/x)*x*Exp(-Ln(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1)*x)*x*Exp((-Ln(y))/x))/(y*x*(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1))-(-Exp(-Ln(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1)*x)*x*Exp((-Ln(z))/x)*z*x*Exp((-Ln(y))/x))/(y*x))/(z*x*(Exp((-Ln(y))/x)+Exp((-Ln(z))/x)-1))^2
}


library(numDeriv)
jacobian(co,c(2,0.9,0.9))
dco(c(2,0.9,0.9))
###
hessian(co,c(2,0.9,0.9))
d2co(c(2,0.9,0.9))

## }}} 


library(mets)
pla <- function(par) {
theta <- par[1]; cif1 <- par[2]; cif2 <- par[3]; 
out <- plack.cif2(cif1,cif2,theta)
return(out)
}

pla2 <- function(par) { ## {{{
x <- exp(par[1]); y <- par[2]; z<- par[3]; 
(1+(y+z)*(x-1)-sqrt(sqr(1+(y+z)*(x-1))-4*x*(x-1)*y*z))/(2*(x-1))
} ## }}}

pla2 <- function(par) { ## {{{
x <- exp(par[1]); y <- par[2]; z<- par[3]; 
(1+(y+z)*(x-1)-sqrt(sqr(1+(y+z)*(x-1))-4*x*(x-1)*y*z))/(2*(x-1))
} ## }}}

pla3 <- function(par) { ## {{{
x <- exp(par[1]); y <- par[2]; z<- par[3]; 
(1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1))
} ## }}}


sqr <- function(x) x^2
Sqrt <- function(x) sqrt(x) 
Exp <- function(x) exp(x)
pow <- function(x,p) x^p
Ln <- function(x) log(x)
sqr(2)

###
###D(y) (1+(y+z)*(x-1)-Sqrt(Exp(Ln(1+(y+z)*(x-1))*2)-4*x*(x-1)*y*z))/(2*(x-1))
### (x-1-(Deriv(y)sqrt(sqr((y+z)*(x-1)+1)-4*x*(x-1)*y*z)))/(2*(x-1))
###D(z) (x-1-(sqrt(sqr((y+z)*(x-1)+1)-4*x*(x-1)*y*z)))/(2*(x-1))
### (-(sqrt(sqr((y+z)*(x-1)+1)-4*x*(x-1)*y*z)))/(2*(x-1))

dplayac <- function(par)
{
x <- exp(par[1]); y <- par[2]; z <- par[3]; 
###(x-1-(sqrt(sqr((y+z)*(x-1)+1)-4*x*(x-1)*y*z)))/(2*(x-1))
###(x-1-(sqrt(exp(2*log((y+z)*(x-1)+1))-4*x*(x-1)*y*z)))/(2*(x-1))
(x-1-((Exp(2*Ln((y+z)*(x-1)+1))*2*(x-1))/((y+z)*(x-1)+1)-4*x*(x-1)*z)/(2*Sqrt(Exp(2*Ln((y+z)*(x-1)+1))-4*x*(x-1)*y*z)))/(2*(x-1))
}

dwolframz<- function(par)
{
x <- exp(par[1]); y <- par[2]; z <- par[3]; 
###(x-1-(sqrt(sqr((y+z)*(x-1)+1)-4*x*(x-1)*y*z)))/(2*(x-1))
###(x-1-(sqrt(exp(2*log((y+z)*(x-1)+1))-4*x*(x-1)*y*z)))/(2*(x-1))
(-1 + y + x* y + z - x* z + Sqrt(-4* (-1 + x)* x* y* z + (1 + (-1 + x)* (y + z))^2))/(2 *Sqrt(-4* (-1 + x)* x* y* z + (1 + (-1 + x)* (y + z))^2))
}
x <- 1; y <- 0.5; z <- 0.5
dwolframz(c(1,0.3,0.5))
dplayac(c(1,0.3,0.5))

d2pla2 <- function(par) {
x <- exp(par[1]); y <- par[2]; z<- par[3]; 
###(-(sqrt(sqr((y+z)*(x-1)+1)-4*x*(x-1)*y*z)))/(2*(x-1))
###(-(2*Sqrt(Exp(2*Ln((y+z)*(x-1)+1))-4*x*(x-1)*y*z)*(((((y+z)*(x-1)+1)*(x-1)*2*Exp(2*Ln((y+z)*(x-1)+1))*2*(x-1))/((y+z)*(x-1)+1)-2*Exp(2*Ln((y+z)*(x-1)+1))* pow(x-1,2))/pow((y+z)*(x-1)+1,2)-4*x*(x-1))-(((2*Exp(2*Ln((y+z)*(x-1)+1))*(x-1))/((y+z)*(x-1)+1)-4*x*(x-1)*z)*2*((Exp(2*Ln((y+z)*(x-1)+1))*2*(x-1))/((y+z)*(x-1)+1)-4*x*(x-1)*y))/(2*Sqrt(Exp(2*Ln((y+z)*(x-1)+1))-4*x*(x-1)*y*z)))/pow(2*Sqrt(Exp(2*Ln((y+z)*(x-1)+1))-4*x*(x-1)*y*z),2))/(2*(x-1))
(x* (1 + (-1 + x)* z + y *(-1 + x + 2* z - 2* x* z)))/(-4* (-1 + x)* x* y* z + (1 + (-1 + x)* (y + z))^2)^(3/2)
}

### dz: (-1 + y + x* y + z - x* z + Sqrt(-4* (-1 + x)* x* y* z + (1 + (-1 + x)* (y + z))^2))/(2 *Sqrt(-4* (-1 + x)* x* y* z + (1 + (-1 + x)* (y + z))^2))
### dy: (-1 + y - x* y + z + x* z + Sqrt(-4* (-1 + x)* x* y* z + (1 + (-1 + x)* (y + z))^2))/(2* Sqrt(-4* (-1 + x)* x* y* z + (1 + (-1 + x)* (y + z))^2))
### dy dz : (x (1 + (-1 + x) z + y (-1 + x + 2 z - 2 x z)))/(-4 (-1 + x) x y z + (1 + (-1 + x) (y + z))^2)^(3/2)

### wolfram : (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1))


d2pla2(c(-2.302585,0.934387,0.893903))
###pla2(c(-2.302585,0.934387,0.893903))

hessian(pla,c(-2.302585,0.934387,0.893903))
hessian(pla2,c(-2.302585,0.934387,0.893903))
d2pla2(c(-2.302585,0.934387,0.893903))

hessian(pla2,c(2.302585,0.934387,0.893903))
d2pla2(c(2.302585,0.934387,0.893903))

sqr <- function(x) x^2
sqr(2)

d2pla <- function(par) { ## {{{
theta <- par[1]; cif1 <- par[2]; cif2 <- par[3]; 
theta <- exp(theta); 

cifs=cif1+cif2; S=1+cifs*(theta-1); S2=4*cif1*cif2*theta*(theta-1);
a=(1+(theta-1)*(cifs)); 

valn=theta*((theta-1)*(cif1-cif2)^2+1); 
val1=((S^2)-S2)^0.5; 
valr=valn/val1; 

return(valr)
} ## }}}

x <- seq(0,1,length=100)
f <- sapply(x,function(x) plack.cif2(-2,x,0.9))
plot(x,f)
f <- sapply(x,function(x) plack.cif2(-2,0.9,x))
plot(x,f)
x <- seq(-3,3,length=100)
f <- sapply(x,function(x) plack.cif2(x,0.9,0.9))
plot(x,f)

curve(pla2,-3,3)

co(c(2,0.9,0.9))
pla(c(2,0.9,0.9))
pla2(c(2,0.9,0.9))
pla3(c(2,0.9,0.9))
d2pla(c(2,0.9,0.9))

jacobian(pla,c(2,0.9,0.9))
jacobian(pla2,c(2,0.9,0.9))
dplayac(c(2,0.9,0.9))
dwolfram(c(2,0.9,0.9))

plack.cif2(0.934387,0.893903,-2.30)
pla2(c(-2.302585,0.934387,0.893903))
dplayac(c(-2.302585,0.934387,0.893903))
dwolframz(c(-2.302585,0.934387,0.893903))
jacobian(pla2,c(-2.302585,0.934387,0.893903))


library(numDeriv)
hessian(co,c(2,0.9,0.1))
hessian(pla,c(2,0.9,0.9))
d2pla(c(2,0.9,0.9))
d2pla2(c(2,0.9,0.9))


library(numDeriv)
hessian(co,c(2,0.9,0.1))
hessian(pla,c(2,0.9,0.9))
d2pla(c(2,0.9,0.9))
d2pla2(c(2,0.9,0.9))

## numerisk test for ortogonalitet i binomial OR model, 
binlike <- function(par) { 
beta <- par[1]; theta <- par[2]
cif1 <- exp(0.1+beta)/(1+exp(0.1+beta))
cif2 <- exp(0.9+beta)/(1+exp(0.9+beta))
p11 <- plack.cif2(cif1,cif2,theta)
p01 <- cif2-p11
p10 <- cif1-p11
p00 <- 1- cif1-cif2+p11 
vec <- c(p11,p01,p10,p00)
return(vec)
}
###
pp <- binlike(c(0.1,1))
jj <- jacobian(binlike,c(0.1,1))
sum(jj[,1]*jj[,2]/pp)


########################################################################################################
########################################################################################################
komud<-function(){  ## {{{ 

library(mets); 
###
d <- subset(simClaytonOakes(2000,2,2,1,stoptime=2,left=2),!truncated)
###
e <- ClaytonOakes(Surv(lefttime,time,status)~x1+cluster(~1,cluster),data=d)
summary(e)
m <- coxph(Surv(time,status)~x1,data=d)
summary(m)
em0 <-two.stage(m,data=d,detail=0,clusters=d$cluster,theta=1)
summary(em0)
em1 <-two.stage(m,data=d,detail=0,clusters=d$cluster,var.link=1)
summary(em1)

source("plackettMLE.R")
emn1 <-twostage(m,data=d,detail=1,clusters=d$cluster,model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=1)
emn1$score
emn0 <-twostage(m,data=d,detail=1,clusters=d$cluster,model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=0)
emn0$score


em <-twostage(m,data=d,detail=1,clusters=d$cluster,model="plackett",score.method="nlminb",control=list(trace=TRUE),var.link=1)
em <-twostage(m,data=d,detail=1,clusters=d$cluster,model="plackett",score.method="optimize",control=list(trace=TRUE),var.link=1)
em <-twostage(m,data=d,detail=1,clusters=d$cluster,model="plackett",score.method="nlminb",control=list(trace=TRUE),var.link=1)


dats=list(time=c(x,y),id=rep(1:n,2),status=rep(1,2*n),x=runif(2*n))
margs <- coxph(Surv(time,status)~+prop(x),data=dats)
plph<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id, model="plackett",score.method="nlminb")
plph$theta
coef(plph)
summary(plph)
plph$score
print(plph)

plph1<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id,
		control=list(trace=TRUE),model="plackett",score.method="fisher.scoring",var.link=1)
summary(plph1)


###
library(mets)
d <- subset(simClaytonOakes(2000,4,3,1,stoptime=2,left=2),!truncated)
###
e <- ClaytonOakes(Surv(time,status)~x1+cluster(~1,cluster),data=d)
summary(e)
m <- coxph(Surv(time,status)~x1,data=d)
summary(m)
em <-two.stage(m,data=d,detail=0,clusters=d$cluster)
source("plackettMLE.R")
em1 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=1)
em2 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=0)
em22 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=0,theta=exp(em1$theta))
summary(em)
summary(em1)
summary(em2)
summary(em22)
coef(em1)
coef(em2)
c(em1$score,em1$score1)
c(em1$Dscore,em1$hess)
em$theta.score
em2$score
em22$score
em1$score

em1 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="plackett",control=list(trace=TRUE),var.link=1)
em11 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="plackett",control=list(trace=TRUE),var.link=1,theta=log(em$theta))
em2 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="plackett",control=list(trace=TRUE),var.link=0,theta=0.5)
em2 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="plackett",control=list(trace=TRUE),var.link=0,theta=1)
em2 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="plackett",control=list(trace=TRUE),var.link=0,theta=2)
em1$score
em2$score
coef(em1)
coef(em2)

e <- ClaytonOakes(Surv(lefttime,time,status)~x1+cluster(~1,cluster),data=d)
summary(e)
m <- coxph(Surv(lefttime,time,status)~x1,data=d)
summary(m)
em <-two.stage(m,data=d,detail=0,clusters=d$cluster)
source("plackettMLE.R")
em1 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=1)
em2 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="clayton.oakes",score.method="nlminb",control=list(trace=TRUE),var.link=0)
summary(em)
summary(em1)
summary(em2)
coef(em1)
coef(em2)
em2$score
em1$score

em1 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="plackett",control=list(trace=TRUE),var.link=1)
em2 <-twostage(m,data=d,detail=0,clusters=d$cluster,model="plackett",control=list(trace=TRUE),var.link=0)
em1$score
em2$score
coef(em1)
coef(em2)



dats=list(time=c(x,y),id=rep(1:n,2),status=rep(1,2*n),x=runif(2*n))
margs <- coxph(Surv(time,status)~+prop(x),data=dats)
plph<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id, model="plackett",score.method="nlminb")
plph$theta
coef(plph)
summary(plph)
print(plph)

plph1<-twostage(margs,data=dats,theta=0.1,detail=1,Nit=0,clusters=dats$id,
		control=list(trace=TRUE),model="plackett",score.method="fisher.scoring",var.link=1)
summary(plph1)
} ## }}}

###########   Wofram afledte  Plackett

CForm[  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dy  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dz  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dy d/dz  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 

CForm[d/dx   (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dx d/dy  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dx d/dz  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dx d/dy d/dz  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 

wol <- function(par) {
      x <- par[1]; y <- par[2]; z <- par[3];
(1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1))  
}


CForm[  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 

CForm[d/dx  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dy  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 
CForm[d/dz  (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ] 

dwol <- function(par) { ## {{{
      x <- par[1]; y <- par[2]; z <- par[3];
###  d/dx, d/dy, d/dz
c(
(y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*pow(-1 + x,2)),
(-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)),
(-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) )
} ## }}}

wol <- function(par) { ## {{{
      x <- par[1]; y <- par[2]; z <- par[3];
###  f,  d/dy f, d/dz f, d/dz d/dy f
c( (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ,
(-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)),
(-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) ,
(((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) )
} ## }}}

wollike <- function(par) { ## {{{
      x <- par[1]; y <- par[2]; z <- par[3];
###  f,  d/dy f, d/dz f, d/dz d/dy f
c( (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ,
(-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)),
(-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) ,
(((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) )
} ## }}}

d2wol <- function(par) { ## {{{
      x <- par[1]; y <- par[2]; z <- par[3];
### d/dx d/dx, d/dx d/dy, d/dx d/dz
c( 
  (pow(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)),2)/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-8*y*z + 2*pow(y + z,2))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/pow(-1 + x,2) + (1 + (-1 + x)*(y + z) - Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/pow(-1 + x,3),
  (1 + ((-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)),

  (1 + ((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)))
} ## }}}

### d/dy d/dz

pow <- function(x,y) x^y
Sqrt <- function(x) sqrt(x)
library(numDeriv)

dCghosh <- function(par) {
x <- par[1]; y <- par[2]; z <- par[3];
0.5-( ((1+(y+z)*(x-1))/2)-2*x*y) /sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z)
}


dloglike <- function(par) { ## {{{
      x <- par[1]; y <- par[2]; z <- par[3];
### d/dx , d/dx d/dy, d/dx d/dz, d/dx d/dz d/dy
c( 
(y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*pow(-1 + x,2)),
  (1 + ((-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)),
  (1 + ((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)),
((-3*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(8.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),2.5)) + ((-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((2*pow(-1 + x,2) - 4*(-1 + x)*x)*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + (2*x)/Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*(-1 + x)) - (((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2))
  )
} ## }}}
dloglike(c(0.1,0.934387,0.893903))

jacobian(wol,c(0.1,0.134387,0.893903))
dloglike(c(0.1,0.134387,0.893903))
dCghosh(c(0.1,0.134387,0.893903))


### d/dx d/dy d/dz
((-3*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(8.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),2.5)) + ((-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((2*pow(-1 + x,2) - 4*(-1 + x)*x)*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + (2*x)/Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*(-1 + x)) - (((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2))

pow <- function(x,y) x^y
Sqrt <- function(x) sqrt(x)

library(numDeriv)
jacobian(wol,c(0.1,0.934387,0.893903))
dwol(c(0.1,0.934387,0.893903))

jacobian(wol,c(5.1,0.934387,0.893903))
dwol(c(5.1,0.934387,0.893903))
d2wol(c(5.1,0.934387,0.893903))

jacobian(dwol,c(0.1,0.934387,0.893903))
hessian(wol,c(0.1,0.934387,0.893903))
d2wol(c(0.1,0.934387,0.893903))

wollike(c(0.1,0.934387,0.893903))
jacobian(wollike,c(0.1,0.934387,0.893903))
dloglike(c(0.1,0.934387,0.893903))

wollike(c(3.17,.5,.6))
dloglike(c(3.17,.5,.6))

wollike <- function(par) { ## {{{
x <- par[1]; y <- par[2]; z <- par[3]; beta <- par[4]
y <- exp(log(y)*exp(beta))
z <- exp(log(z)*exp(beta))
###  f,  d/dy f, d/dz f, d/dz d/dy f
c( (1+(y+z)*(x-1)-sqrt((1+(y+z)*(x-1))^2-4*x*(x-1)*y*z))/(2*(x-1)) ,
(-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)),
(-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) ,
(((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) )
} ## }}}
dloglike <- function(par) { ## {{{
x <- par[1]; y <- par[2]; z <- par[3]; beta <- par[4]
y <- exp(log(y)*exp(beta))
z <- exp(log(z)*exp(beta))
### d/dx , d/dx d/dy, d/dx d/dz, d/dx d/dz d/dy
vec <- 
c( 
(y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*pow(-1 + x,2)),
  (1 + ((-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)),
  (1 + ((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2)),
((-3*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(8.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),2.5)) + ((-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((2*pow(-1 + x,2) - 4*(-1 + x)*x)*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + (2*x)/Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*(-1 + x)) - (((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*Sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2))
  )
} ## }}}
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

apply(dthetaloglike,c(1.1,0.134387,0.893903,1.3))


###########   Wofram afledte  Clayton-Oakes
co <- function(par) {
x <- par[1]; y <- par[2]; z<- par[3]; 
((y^(-1/x)+z^(-1/x))-1)^(-x)
###Exp(-Ln(Exp(Ln(y)*(-1/x))+ Exp(Ln(z)*(-1/x))-1 )*x)
###exp(-log(exp(log(y)*(-x))+ exp(log(z)*(-x))-1 )*(1/x));
}
co(c(2,0.5,0.6))

1/pow(exp(log(y)*(-1/x))+ exp(log(z)*(-1/x))-1 ,x)
sqr <- function(x) x^2
Sqrt <- function(x) sqrt(x) 
Exp <- function(x) exp(x)
pow <- function(x,p) x^p
Ln <- function(x) log(x)
sqr(2)


CForm[   ((1/y^(1/x)+1/z^(1/x))-1)^(-x) ]
CForm[d/dx   ((y^(-1/x)+z^(-1/x))-1)^(-x) ]
CForm[d/dz   ((y^(-1/x)+z^(-1/x))-1)^(-x) ]
CForm[d/dy   ((y^(-1/x)+z^(-1/x))-1)^(-x) ]
CForm[d/dz d/dy   ((y^(-1/x)+z^(-1/x))-1)^(-x) ]


colike <- function(par) {
x <- par[1]; y <- par[2]; z<- par[3]; 
val1= pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
val2=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
val3=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
val4= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
c(val1,val2,val3,val4)
}
colike(c(1,0.5,0.5))

dcolike <- function(par) {
x <- par[1]; y <- par[2]; z<- par[3]; 
dp1=(-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
dp2=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
dp3=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
dp4=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
c(dp1,dp2,dp3,dp4)
}

dcolike(c(1.1,0.8,0.5))
jacobian(colike,c(1.1,0.8,0.5))

komud<-function(){  ## {{{ 


if (status1==0 && status2==0) { // {{{
valr= pow(-1 + pow(y,-x) + pow(z,-x),-1/x);
dp(0)=(-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
} // }}}

if (status1==1 && status2==0) { // {{{
valr=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
dp(0)=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
} // }}}

if (status1==0 && status2==1) { // {{{
valr=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
dp(0)=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
} // }}}

if (status1==1 && status2==1) { // {{{
valr= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
dp(0)=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
} // }}}

} ## }}}


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



###########################################################################

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




