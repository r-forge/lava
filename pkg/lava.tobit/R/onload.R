'.onLoad' <- function(lib, pkg="lava.tobit") {
  addhook("lava.tobit.estimate.hook","estimate.hooks")
  addhook("lava.tobit.color.hook","color.hooks")
  desc <- packageDescription(pkg)
  cat("Loading '", desc$Package, "' package...\n", sep="")
  cat("Version    : ", desc$Version, "\n", sep="")
}
