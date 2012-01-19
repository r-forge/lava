##install(c("Rcpp","RcppArmadillo","inline"))
library(RcppArmadillo)
library(inline)

src <- 
'
Rcpp::NumericVector xa(a);
Rcpp::NumericVector xb(b);
int n_xa = xa.size(), n_xb = xb.size();
Rcpp::NumericVector xab(n_xa + n_xb - 1);
for (int i = 0; i < n_xa; i++)
for (int j = 0; j < n_xb; j++)
xab[i + j] += xa[i] * xb[j];
return xab;
'

## turn into a function that R can call
## compileargs redundant on Debian/Ubuntu as gsl headers are found anyway
funx <- cxxfunction(signature(a="numeric", b="numeric"),
                  src,
                  ##includes="#include <gsl/gsl_rng.h>",
                  ##,cppargs="-I/usr/include",libargs="-lgsl -lgslcblas"
                  plugin="RcppArmadillo")
funx(1:3,1:4)
