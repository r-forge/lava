'.onLoad' <- function(lib, pkg="bptwin")
  {    
    desc <- packageDescription(pkg)
    cat("Loading '", desc$Package, "' package...\n", sep="")
    cat("Version:\t ", desc$Version, "\n", sep="")
    cat("Overview:\t help(package=", desc$Package, ")\n", sep="");   
  }
