`binary` <- function(x,...) UseMethod("binary")
"binary<-" <- function(x,...,value) UseMethod("binary<-")

"binary<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    return(binary(x,all.vars(value),...))
  }
  binary(x, value, ...)
}

`binary.lvm` <-
function(x,var=NULL, ...) {
  if (is.null(var)) {
    ## binidx <- tryCatch(unlist(nodeData(Graph(x), attr="binary")),error=function(e) NULL)
    binidx <- unlist(x$attributes$binary)
    if (length(binidx)>0)
      return(names(binidx)[binidx])
    else
      return(NULL)
  }
  ##  if (is.null(nodeDataDefaults(Graph(x))$binary)) {
  ##    nodeDataDefaults(Graph(x),"binary") <- FALSE
  ##  } 
  
  ##  x <- addattr(x,attr="shape",var=var,val="box")
  x$attributes$binary[var] <- TRUE
  x$attributes$normal[var] <- FALSE
  ## nodeData(Graph(x), var, attr="binary") <- TRUE
  ## nodeData(Graph(x), var, attr="normal") <- FALSE
  covfix(x,var,NULL) <- 1
  ##  distribution(x, var) <- probit.lvm
  return(x)
}

