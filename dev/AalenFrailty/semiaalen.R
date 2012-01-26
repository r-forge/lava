library(RcppArmadillo)
library(inline)

con <- file("semiaalen.cpp"); src1 <- readLines(con); close(con)
semiaalen <- function(eventidx,time,X,Z,dt=NULL,gammas=NULL,allt=FALSE) csemiaalen(eventidx,time,X,Z,dt,gammas,allt)
csemiaalen <- cxxfunction(signature(ds="numeric", ts="numeric",
                                   Xs="numeric", Zs="numeric",
                                   dts="numeric", gammas="numeric",
                                   allt="boolean"),
                         paste(src1,collapse="\n"),
                         plugin="RcppArmadillo", verbose=TRUE)

t1 <- system.time(B1 <- semiaalen(event,t,X,Z,allt=FALSE))

t1 <- t[which(event==1)]
event1 <- rep(1,length(t1))
t2 <- system.time(B2 <- semiaalen(event1,t1,X,Z))



library(timereg)
data(sTRACE)
dd <- sTRACE
for (i in seq((NULL))) dd <- rbind(dd,dd)
## <- sTRACE
dd$time <- dd$time+runif(nrow(dd),0,1e-8)
dd <- dd[order(dd$time),]
X <- model.matrix(~1+chf,dd)
Z <- model.matrix(~-1+sex+age,dd)
t <- dd$time
event <- (dd$status==9)*1
dix <- which(dd$status==9)
dt0 <- diff(c(0,t[dix]))
t1 <- system.time(B <- semiaalen(event,t,X,Z,NULL,NULL))
system.time(out<-aalen(Surv(time,status==9)~1+chf+const(sex)+const(age),dd,n.sim=0,robust=0))


system.time(out<-aalen(Surv(time,status==9)~1+chf+const(sex)+const(age),dd,n.sim=0,robust=0))



dt <- diff(c(0,t))

system.time(B <- semiaalen(dix,t,X,Z,dt))


system.time(out<-aalen(Surv(time,status==9)~1+chf+const(sex)+const(age),dd,n.sim=0,robust=0))


system.time(out<-aalen(Surv(time,status==9)~1+chf+const(sex)+const(age),dd,n.sim=0,robust=0))#max.time=7.5))
(b <- coef(out))
out$intZHdN
out$intZHZ
system.time(B <- semiaalen(dix,t,X,Z,dt,coef(out)[,1]))
system.time(B <- semiaalen(dix,t,X,Z,dt,NULL))
idx <- 1
plot(out$cum[,c(1,idx+1)],type="s")
lines(cbind(t[dix],B$B[,idx]),type="s",col=Col("darkred",0.4),lwd=5)



library(timereg)
data(sTRACE)
dd <- sTRACE
for (i in seq((NULL))) dd <- rbind(dd,dd)
## <- sTRACE
dd$time <- dd$time+runif(nrow(dd),0,1e-8)
dd <- dd[order(dd$time),]
X <- model.matrix(~1+chf,dd)
Z <- model.matrix(~-1+sex+age,dd)
t <- dd$time
dix <- which(dd$status==9)
t1 <- system.time(B <- semiaalen(dix,t,X,Z,NULL,NULL))
t2 <- system.time(out<-aalen(Surv(time,status==9)~1+chf+const(sex)+const(age),dd,n.sim=0,robust=0))#max.time=7.5))
res <- c(t1[1],t2[1],B$gamma[,1],coef(out)[,1])



X <- model.matrix(~1+chf,dd)
Z <- model.matrix(~-1+sex+age,dd)
t <- dd$time
dix <- which(dd$status==9)
dt <- diff(c(0,t[dix]))
system.time(B <- semiaalen(dix,t,X,Z,dt))
B$gamma

system.time(out<-aalen(Surv(time,status==9)~1+chf+const(sex)+const(age),dd,n.sim=0,robust=0))#max.time=7.5))
(b <- coef(out))
out$intZHdN
out$intZHZ
system.time(B <- semiaalen(dix,t,X,Z,dt))
plot(out$cum,type="s")
lines(cbind(t[dix],B$B),type="s",col=Col("darkred",0.4),lwd=5)

surv <- with(dd,Surv(time+runif(nrow(dd),0,0.0001),status==9))
library(ahaz)
system.time(outss <- ahaz(surv,Z))


library(timereg)
library(bsplot)
data(sTRACE)
sTRACE <- sTRACE[order(sTRACE$time),]
system.time(out<-aalen(Surv(time,status==9)~1+const(sex)+const(age),sTRACE,n.sim=0,robust=0))
(b <- coef(out))
idx <- 1
plot(out,specific.comps=idx)
points(t[dix],B$B[,idx],type="s",col=Col(2),lwd=7)




##system.time(B <- EvilAalen(dix,X))

system.time(out<-aalen(Surv(time,status==9)~1,sTRACE,max.time=9,n.sim=0,robust=0))
idx <- 1
plot(out,specific.comps=idx)
points(t[dix],B$B[,idx],type="s",col=Col(2),lwd=7)


##################################################
### Benchmark
##################################################
f <- function(rep=NULL) {
  message(rep)
  library(timereg)
  data(sTRACE)
  dd <- sTRACE
  for (i in seq((rep))) dd <- rbind(dd,dd)
  ## <- sTRACE
  dd$time <- dd$time+runif(nrow(dd),0,1e-8)
  dd <- dd[order(dd$time),]
  X <- model.matrix(~1+chf,dd)
  Z <- model.matrix(~-1+sex+age,dd)
  t <- dd$time
  dix <- which(dd$status==9)
  t1 <- system.time(B <- semiaalen(dix,t,X,Z,NULL,NULL))
  res <- c()
  if (rep<9) {
    t2 <- system.time(out<-aalen(Surv(time,status==9)~1+chf+const(sex)+const(age),dd,n.sim=0,robust=0))#max.time=7.5))
    res <- c(n=nrow(dd),t.kkho=t1[1],t.ts=t2[1],B$gamma[,1],coef(out)[,1])
  } else {
    res <- c(n=nrow(dd),t.kkho=t1[1],t.ts=NA,B$gamma[,1],NA,NA)
  }
  return(res)
}

library(foreach)
res <- foreach(i=1:11) %do% f(i)

par(mfrow=c(1,2))
M <- matrix(unlist(res),byrow=TRUE,ncol=7)
plot(M[,1],M[,3],col="red",type="l")
lines(M[,1],M[,2])
legend("topleft",c("timereg::aalen","semiaalen"),lty=1,col=2:1)
plot(M[1:8,1],M[1:8,3],col="red",type="l")
lines(M[,1],M[,2])


dev.copy2pdf(file="test.pdf")

plot(M[,1],sqrt(M[,3]),col="red",type="l")
lines(M[,1],sqrt()(M[,2]))
