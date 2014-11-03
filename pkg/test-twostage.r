
library(mets)
###example(twostage)
###example(easy.twostage)
###

library(mets)
x=1:2

.Call("claytonoakesR",1,1,1,0.5,0.3)
.Call("claytonoakesR",1,1,0,0.5,0.3)
.Call("claytonoakesR",1,0,1,0.5,0.3)
.Call("claytonoakesR",1,0,0,0.5,0.3)


theta <- runif(5)+0.1
s1 <- rbinom(5,1,0.5)
s2 <- rbinom(5,1,0.5)
p1 <- runif(5)
p2 <- runif(5)

.Call("claytonoakesR",theta,s1,s2,p1,p2)


