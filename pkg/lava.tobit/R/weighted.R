tobitw_method.lvm <- "estfun"
tobitw_gradient.lvm <- tobit_gradient.lvm
tobitw_hessian.lvm <- tobit_hessian.lvm


atobitw_method.lvm <- "nlminb1"
atobitw_objective.lvm <- function(...) {
  S <- tobit_gradient.lvm(...)
  crossprod(S)[1]
}
atobitw_gradient.lvm <- function(p,...) {
  myfun <- function(p0) tobit_gradient.lvm(p=p0,...)
  H <- jacobian(myfun,p,method=lava.options()$Dmethod)
  S <- myfun(p)
  2*S%*%H
}

## tobitw_method.lvm <- "estfun"
## tobitw_objective.lvm <- function(...) {
##   S <- tobit_gradient.lvm(...)
##   crossprod(S)[1]
## }
## tobitw_gradient.lvm <- function(p,...) {
##   myfun <- function(p0) tobit_gradient.lvm(p=p0,...)
##   H <- jacobian(myfun,p,method=lava.options()$Dmethod)
##   S <- myfun(p)
##   2*S%*%H
## }
##obitw_hessian.lvm <- NULL

tobitw2_method.lvm <- "NR"
tobitw2_objective.lvm <- function(...) {
  S <- tobit_gradient.lvm(...)
  crossprod(S)[1]
}
tobitw2_gradient.lvm <- function(...) {
  tobit_gradient.lvm(...)
}
tobitw2_hessian.lvm <- function(p,...) {
  S <- tobit_gradient.lvm(p=p,...)
  myfun <- function(p0) tobitw_gradient.lvm(p=p0,...)
  H <- jacobian(myfun,p,method=lava.options()$Dmethod)
  attributes(H)$grad <- S
  H
}

