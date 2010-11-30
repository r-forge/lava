'.onLoad' <- function(lib, pkg="lava.tobit") {
  addhook("lava.tobit.estimate.hook","estimate.hooks")
  addhook("lava.tobit.color.hook","color.hooks")
  addhook("lava.tobit.sim.hook","sim.hooks")
  addhook("lava.tobit.init.hook","init.hooks")

  lava.options(tobitAlgorithm=mvtnorm::GenzBretz(abseps=1e-5),
               tobitseed=1)
  desc <- utils::packageDescription(pkg)
  cat("Loading '", desc$Package, "' package...\n", sep="")
  cat("Version    : ", desc$Version, "\n", sep="")
}
