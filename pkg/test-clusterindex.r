
library(mets)

sessionInfo()

index <- c(1,1,2,2,1)

clusters <- cluster.index(index)

clusters <- cluster.index(index,Rindex=1)
ud <- familycluster.index(index)
ud

index <- c(5,5,2,2,10)
clusters <- cluster.index(index,Rindex=1)
ud <- familycluster.index(index)
ud

