###{{{ Misc
Col <- function (col, alpha = 0.2) 
{
    sapply(col, function(x) do.call(rgb, as.list(c(col2rgb(x)/255, 
        alpha))))
}

revdiag <- function(x) {
  n <- ncol(x)
  x[cbind(rev(seq(n)),seq(n))]
}
"revdiag<-" <- function(x,value,...) {
  n <- ncol(x)
  x[cbind(rev(seq(n)),seq(n))] <- value
  x
}


###}}} Misc

