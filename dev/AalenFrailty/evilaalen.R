##install(c("Rcpp","RcppArmadillo","inline"))
library(RcppArmadillo)
library(inline)

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

library(timereg)
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


par(mfrow=c(2,2))
for (i in 1:3) {
  plot(out,specific.comps=i); lines(t[dix],B[,i],type="s",col="red")
}
