.libPaths("~/.R/library.x86_64/")
library(mets)

set.seed(1)
d <- subset(simClaytonOakes(500,4,2,1,stoptime=2,left=2),!truncated)
e <- ClaytonOakes(Surv(lefttime,time,status)~x1+cluster(~1,cluster),cuts=c(0,0.5,1,2),data=d)


