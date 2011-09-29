bootstrap <- function (x, ...) UseMethod("bootstrap")
bootstrap.bptwin <- function(x,...) {
  data <- mycall$data
  
  mycall <- fitted$call
  mycall$constrain <- c(NA,0.423258,NA)
  mycall$stderr <- FALSE
  mycall$data=as.name("data")
  with(fitted, eval(mycall))  
}

profile.bptwin <- function(fitted,...) {
  mycall <- fitted$call
  mycall$constrain <- c(NA,0.423258,NA)
  mycall$stderr <- FALSE
  mycall$data=as.name("data")
  with(fitted, eval(mycall))
}
