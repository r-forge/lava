##install(c("Rcpp","RcppArmadillo","inline"))
library(RcppArmadillo)
library(inline)


## now use Rcpp to pass down a parameter for the seed, and a vector size
src <- '
try {
  arma::colvec y = Rcpp::as<arma::colvec>(ys);    // direct from SEXP to arma::mat
  arma::mat X    = Rcpp::as<arma::mat>(Xs);
  int df = X.n_rows - X.n_cols;

  arma::colvec coef = arma::solve(X, y);          // fit model y ~ X
//  arma::colvec res  = y - X*coef;                 // residuals
//  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/df; // std.errors of coefficients
//  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec( arma::pinv(arma::trans(X)*X) ));
    return(Rcpp::List::create(Rcpp::Named("coefficients")=coef));
//  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
//                                  Rcpp::Named("stderr")       = std_err,
//                                  Rcpp::Named("df")           = df
//                                  );
  } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
  } catch(...) { 
       ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
'
ols <- cxxfunction(signature(ys="numeric", Xs="numeric"),
                  src,
                  ##includes="#include <gsl/gsl_rng.h>",
                  ##,cppargs="-I/usr/include",libargs="-lgsl -lgslcblas"
                  plugin="RcppArmadillo")


src0 <- '
using namespace arma;
try {
      colvec   event = Rcpp::as<colvec>(ds); 
      mat         X = Rcpp::as<mat>(Xs);

      mat b = zeros(event.n_rows,X.n_cols);
      mat Xc = zeros<mat>(X.n_cols,X.n_cols);
      unsigned start, ri;
      for (unsigned i=0; i<event.n_rows; i++) {
         ri = event.n_rows-i-1;
         start = event(ri)-1;
         mat Xi = X(span(start,X.n_rows-1), span::all);
         Xc = Xi.st()*Xi;
         mat U, V; vec s; 
         svd(U,s,V,Xc);
         mat Xci = U*diagmat(1/s)*V.st();
         b.row(ri) = trans(Xci*trans(X.row(start)));
      }
      mat B = cumsum(b);
      return(Rcpp::List::create(Rcpp::Named("B")=B));

  } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
  } catch(...) { 
       ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
'

src1 <- '
using namespace arma;
try {
      colvec   event = Rcpp::as<colvec>(ds); 
      mat         X = Rcpp::as<mat>(Xs);

      mat b = zeros(event.n_rows,X.n_cols);
      mat Xc = zeros<mat>(X.n_cols,X.n_cols);
      unsigned stop,start = X.n_rows;
      for (unsigned i=0; i<event.n_rows; i++) {
         unsigned ri = event.n_rows-i-1;
         stop  = start-1;
         start = event(ri)-1;
         mat Xi = X(span(start,stop), span::all);
         Xc = Xc + Xi.st()*Xi;
         mat U, V; vec s; 
         svd(U,s,V,Xc);
         mat Xci = U*diagmat(1/s)*V.st();
         b.row(ri) = trans(Xci*trans(X.row(start)));
      }
      mat B = cumsum(b);
      return(Rcpp::List::create(Rcpp::Named("B")=B));

  } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
  } catch(...) { 
       ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
'
EvilAalen0 <- cxxfunction(signature(ds="numeric", Xs="numeric"),
                         src0,plugin="RcppArmadillo", verbose=TRUE)

EvilAalen <- cxxfunction(signature(ds="numeric", Xs="numeric"),
                         src1,plugin="RcppArmadillo", verbose=TRUE)


data(sTRACE)
sTRACE <- sTRACE[order(sTRACE$time),]
system.time(out<-aalen(Surv(time,status==9)~age+sex,##+diabetes+chf+vf,
                       sTRACE,max.time=9,n.sim=0))

X <- model.matrix(~1+age+sex,sTRACE)
t <- sTRACE$time
dix <- which(sTRACE$status==9)
rdix <- rev(dix)
system.time(B <- EvilAalen(dix,X)$B)
system.time(B <- EvilAalen0(dix,X)$B)

plot(out,specific.comps=3)
lines(t[dix],B[,3],type="s",col="red")





data(iris)
d <- sim(lvm(Sepal.Length~Petal.Width+Petal.Length),1e5)
X <- model.matrix(Sepal.Length~Petal.Width+Petal.Length,d)
system.time(ols(d$Sepal.Length,X))
system.time(lm(d$Sepal.Length~X-1))


dim(sTRACE)
sTRACE1 <- sTRACE[sTRACE$status!=0,]

library(timereg)
library(bsplot)
data(sTRACE)
sTRACE <- sTRACE[order(sTRACE$time),]
X <- model.matrix(~1+age+sex,sTRACE)
system.time(out<-aalen(Surv(time,status==9)~age+sex,##+diabetes+chf+vf,
                       sTRACE,max.time=9,n.sim=0))

tt <- unique(sTRACE$time)
b <- matrix(0,ncol=ncol(X),nrow=length(tt))
for (i in seq(length(tt))) {
  cat(i,"\n")
  ss <- (sTRACE$time>=tt[i])  
  X0 <- apply(X,2,function(x) x*ss)
  Xi <- crossprod(X0)
  y0 <- (sTRACE$status==9)&(sTRACE$time==tt[i])
  b[i,] <- if (any(eigen(Xi)$values<1e-15)) {
    rep(0,ncol(X))
  } else {      
    tryCatch(ols(y0,X0)$coefficients[,1],
             error=function(...) rep(0,ncol(X)))
  }
}

par(mfrow=c(2,2))

plot(out,specific.comps=1)
lines(tt,apply(b,2,cumsum)[,1],type="s",col=Col(1,0.4),lwd=3)

plot(out,specific.comps=2)
lines(tt,apply(b,2,cumsum)[,2],type="s",col=Col(1,0.4),lwd=3)

plot(out,specific.comps=3)
lines(tt,apply(b,2,cumsum)[,3],type="s",col=Col(1,0.4),lwd=3)



f <- function() {

tidx <- rev(which(sTRACE$status==9))
b <- matrix(0,ncol=ncol(X),nrow=length(tidx))
Xc <- matrix(0,ncol(X0),ncol(X0))
idx <- nrow(sTRACE)+1
for (i in seq(length(tidx))) {
  idx <- seq(tidx[i],idx[1]-1)
  Xi <- X[idx,,drop=FALSE]
  ##  cat(i," ",idx," ",sTRACE$status[idx],"\n")
  Xc <- Xc+crossprod(Xi)
  DD <- svd(Xc)
  b[i,] <- if (any(abs(DD$d)<1e-15)) {
    rep(0,ncol(X))
  } else {      
    tryCatch((with(DD,u%*%diag(1/v)%*%t(v))%*%X[tidx[i],])[,1],
             error=function(...) rep(0,ncol(X)))
  }
}
return(b)
}

ts <- unique(sTRACE$time)
bs <- matrix(0,ncol=ncol(X),nrow=length(ts))
bs[rev(tidx),] <- b[rev(seq(nrow(b))),]

par(mfrow=c(1,1))
plot(out,specific.comps=1)
##abline(v=sTRACE$time[tidx])
points(sTRACE$time[rev(tidx)],cumsum(rev(b[,1])),col="red",type="s")

plot(out,specific.comps=2)
points(tt[rev(tidx)],cumsum(rev(b[,2])),col="red")


plot(out,specific.comps=2)
lines(tt[rev(tidx)],rev(apply(b,2,cumsum)[,3]),type="s",col=Col(1,0.4),lwd=3)

plot(out,specific.comps=2)
lines(tt,apply(b,2,cumsum)[,2],type="s",col=Col(1,0.4),lwd=3)

plot(out,specific.comps=3)
lines(tt,apply(b,2,cumsum)[,3],type="s",col=Col(1,0.4),lwd=3)





