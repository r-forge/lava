tobitw_method.lvm <- "estfun"
tobitw_gradient.lvm <- tobit_gradient.lvm


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

p <- c(-0.0611481, -7.23247e-05, 0.532474, 0.000121360)
