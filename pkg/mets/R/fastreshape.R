##' Simple reshape/tranpose of data
##'
##' @title Fast reshape
##' @param data data.frame or matrix
##' @param id id-variable. If omitted then reshape from wide to long. 
##' @param varying Vector of prefix-names of the time varying
##' variables. Optional for long->wide reshaping.
##' @param num Optional number/time variable
##' @param sep String seperating prefix-name with number/time
##' @param all.numeric If TRUE (all variables numeric) avoid slow \code{cbind}
##' @param ... Optional additional arguments to the \code{reshape} function used in the wide->long reshape.
##' @author Thomas Scheike, Klaus K. Holst
##' @export
##' @examples
##' m <- lvm(c(y1,y2,y3,y4)~x)
##' d <- sim(m,1e1)
##' 
##' dd <- fast.reshape(d,var="y")
##' d1 <- fast.reshape(dd,"id")
##'
##' ## From wide-format
##' d1 <- fast.reshape(dd,"id")
##' d2 <- fast.reshape(dd,"id",var="y")
##' d3 <- fast.reshape(dd,"id",var="y",num="time")
##' 
##' ## From long-format
##' fast.reshape(d,var="y",idvar="a",timevar="b")
##' fast.reshape(d,var=list(c("y1","y2","y3","y4")),idvar="a",timevar="b")
fast.reshape <- function(data,id,varying,num,sep="",all.numeric=FALSE,...) {
  if (NCOL(data)==1) data <- cbind(data)
  
  if (missing(id)) {
    ## reshape from wide to long format. Fall-back to stats::reshape
    nn <- colnames(data)
    nsep <- nchar(sep)
    vnames <- NULL
    if (missing(varying)) stop("Prefix of time-varying variables needed")    
    ncvar <- sapply(varying,nchar)
    newlist <- c()
    if (!is.list(varying)) {
      for (i in seq_len(length(varying))) {
        ii <- which(varying[i]==substr(nn,1,ncvar[i]))
        tt <- as.numeric(substring(nn[ii],ncvar[i]+1+nsep))      
        newlist <- c(newlist,list(nn[ii[order(tt)]]))
      }
      vnames <- varying
      varying <- newlist
    }
    return(reshape(data,varying=varying,direction="long",v.names=vnames,...))
  }

  numvar <- idvar <- NULL 
  if (is.character(id)) {
    if (length(id)>1) stop("Expecting column name or vector of id's")
    idvar <- id
    id <- as.numeric(data[,id,drop=TRUE])
  } else {
    if (length(id)!=nrow(data)) stop("Length of ids and data-set does not agree")
  }    
  if (!missing(num)) {
    if (is.character(num)) {
      numvar <- num
      num <- as.numeric(data[,num,drop=TRUE])
    } else {
      if (length(num)!=nrow(data)) stop("Length of time and data-set does not agree")
    }
  } else {
    num <- NULL
  }    

  antpers <- nrow(data)
  unique.id <- unique(id)  
  clusters <- fast.approx(unique.id,id)$pos
  max.clust <- length(unique.id)
  nclust <- .C("nclusters", as.integer(antpers), as.integer(clusters), 
               as.integer(rep(0, antpers)), as.integer(0), as.integer(0), 
               package = "timereg")
  maxclust <- nclust[[5]]
  antclust <- nclust[[4]]
  cluster.size <- nclust[[3]][seq_len(antclust)]
  if (!is.null(num)) { ### different types in different columns
    mednum <- 1
    numnum <- numnum <- order(num)-1
  } else {
    numnum <- 0;
    mednum <- 0;
  }
  init <- -1
  clustud <- .C("clusterindex", as.integer(clusters), as.integer(antclust), 
                as.integer(antpers),
                as.integer(rep(init, antclust * maxclust)),
                as.integer(rep(0, antclust)), as.integer(mednum), 
                as.integer(numnum), package = "timereg")
  idclust <- matrix(clustud[[4]], antclust, maxclust)
  
  if (!is.null(numvar)) {
    ii <- which(colnames(data)==numvar)
    data <- data[,-ii,drop=FALSE]
  }
  if (missing(varying)) varying <- setdiff(colnames(data),c(idvar))
  vidx <- match(varying,colnames(data))
  N <- nrow(idclust)
  p <- length(varying)
  
  if (all(apply(data[1,],2,is.numeric))) {
  ## Everything numeric - we can work with matrices
    dataw <- matrix(NA, nrow = N, ncol = p * (maxclust-1) + ncol(data))
    for (i in seq_len(maxclust)) {
      if (i==1) {
        dataw[, seq(ncol(data))] <- as.matrix(data[idclust[, i] + 1,])
        mnames <- colnames(data);
        mnames[vidx] <- paste(mnames[vidx],i,sep=sep)
      } else {
        dataw[, seq(p) + (ncol(data)-p) + (i - 1) * p] <- as.matrix(data[idclust[, i] + 1,varying])
        mnames <- c(mnames,paste(varying,i,sep=sep))
      }
    }
##    colnames(dataw) <- mnames
    return(dataw)
  } ## Potentially slower with data.frame where we use cbind
  dataw <- c()  
  mnames <- c()  
  for (i in seq_len(maxclust)) {
     if (i==1) {
       dataw <- data[idclust[,i]+1,,drop=FALSE]
       mnames <- names(data);
       mnames[vidx] <- paste(mnames[vidx],sep,i,sep="")
     } else {
       dataw <- cbind(dataw,data[idclust[,i]+1,varying,drop=FALSE])
       mnames <- c(mnames,paste(varying,sep,i,sep=""))
     }
   }
  names(dataw) <- mnames
  
  return(dataw)
 } 


simple.reshape <- function (data, id = "id", num = NULL) {
    cud <- cluster.index(data[, c(id)], num = num, Rindex = 1)
    N <- nrow(cud$idclust)
    p <- ncol(data)
    dataw <- matrix(NA, nrow = N, ncol = p * cud$maxclust)
    for (i in seq_len(cud$maxclust)) {
           dataw[, seq(p) + (i - 1) * p] <- as.matrix(data[cud$idclust[, i] + 1, ])
    }
   colnames(dataw) <- paste(names(data), rep(seq_len(cud$maxclust), each = p), sep = ".")
   return(dataw)
}






###faster.reshape <- function(data,clusters,index.type=FALSE,num=NULL,Rindex=1)
###{ ## {{{
###data <- as.matrix(data)
###if (NCOL(data)==1) data <- cbind(data)
###
###antpers <- length(clusters)
###if (index.type==FALSE)  {
###	max.clust <- length(unique(clusters))
###	clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1 
###}
###
### nclust <- .C("nclusters",
###	as.integer(antpers), as.integer(clusters), as.integer(rep(0,antpers)), 
###	as.integer(0), as.integer(0), package="timereg")
###  maxclust <- nclust[[5]]
###  antclust <- nclust[[4]]
###  cluster.size <- nclust[[3]][1:antclust]
###
###if ((!is.null(num)) && (Rindex==1)) { ### different types in different columns
###   mednum <- 1
###   numnum <- as.integer(factor(num, labels = 1:maxclust)) -1
###} else { numnum <- 0; mednum <- 0; }
###
###data[is.na(data)] <- nan 
###p <- ncol(data); 
###init <- -1*Rindex;
###clustud <- .C("clusterindexdata",
###	        as.integer(clusters), as.integer(antclust),as.integer(antpers),
###                as.integer(rep(init,antclust*maxclust)),as.integer(rep(0,antclust)), as.integer(mednum), 
###		as.integer(numnum), as.double(c(data)), 
###		as.integer(p), as.double(rep(init*1.0,antclust*maxclust*p)), package="timereg")
###idclust <- matrix(clustud[[4]],antclust,maxclust)
###xny <- matrix(clustud[[10]],antclust,maxclust*p)
###if(Rindex==1) xny[idclust==-1] <- NA 
###if(Rindex==1) xny[idclust==-1] <- NA 
###if(Rindex==1) idclust[idclust==-1] <- NA 
###  mnames <- c()
###print(maxclust)
###  for (i in 1:maxclust) {
###     mnames <- c(mnames,paste(names(data),".",i,sep=""))
###  }
###  xny <- data.frame(xny)
###  names(xny) <- mnames
###out <- xny; 
###} ## }}}
###