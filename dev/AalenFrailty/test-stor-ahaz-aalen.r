library(timereg)
library(MultiComp)
library(HaploSurvival)
data(TRACE)
ktrace <- rbind(TRACE,TRACE,TRACE,TRACE)
ktrace <- rbind(ktrace,ktrace,ktrace,ktrace)
ktrace <- rbind(ktrace,ktrace,ktrace,ktrace)
###ktrace <- rbind(ktrace,ktrace,ktrace,ktrace)
table(ktrace$status)

timent <- ktrace$time+runif(nrow(ktrace))*0.01

system.time(
out <- aalen(Surv(timent,status==9) ~ factor(vf)+factor(chf)+factor(sex),
		data=ktrace,robust=0)
)

plot(out)

X=model.matrix(~vf+chf+sex,data=ktrace)[,-1]
surv <- Surv(timent,ktrace$status==9)
library(ahaz)
system.time(
outss <- ahaz(surv,X)
)
outss$d
outss$D


system.time(
outs <- aalen(Surv(timent,status==9) ~ const(vf)+const(chf)+const(sex),data=ktrace,robust=0,n.sim=0)
)
c(outs$intZHdN)-outss$d
outs$intZHZ - outss$D

