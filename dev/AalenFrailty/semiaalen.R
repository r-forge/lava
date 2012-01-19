library(RcppArmadillo)
library(inline)

src1 <- '
#include <vector>
using namespace arma;
try {
      colvec   event = Rcpp::as<colvec>(ds);
      colvec      t = Rcpp::as<colvec>(ts);
      mat         X = Rcpp::as<mat>(Xs);
      mat         Z = Rcpp::as<mat>(Zs);
      mat b = zeros(event.n_rows,X.n_cols);
      mat Xc = zeros<mat>(X.n_cols,X.n_cols);
      colvec Ht(event.n_rows);
      colvec dt(event.n_rows);
      colvec zHz(event.n_rows);
      mat zH(event.n_rows,Z.n_cols);
      unsigned stop,start = X.n_rows;
      mat I = eye(Z.n_cols,Z.n_cols);
      mat zz;
      mat H;
      for (unsigned i=0; i<event.n_rows; i++) {
         unsigned ri = event.n_rows-i-1;
         stop  = start-1;
         start = event(ri)-1;
         mat Xi = X(span(start,stop), span::all);
         Xc = Xc + Xi.st()*Xi;
         mat U, V; vec s; 
         svd(U,s,V,Xc);
         mat Xci = U*diagmat(1/s)*V.st();
         mat Xpinv = Xci*trans(X.row(start));
         //         Xpinvs.push_back(Xpinv);
         b.row(ri) = trans(Xpinv);
//         mat XProj = Xc*Xci*trans(Xc);
//         cerr << XProj << endl;
         //         Xpinvs.push_back(Xpinv);
         Ht(ri) = 1-as_scalar(X.row(start)*Xpinv);
         for (unsigned k=0; k<s.n_elem; k++) if (s(k)<1.0e-15) Ht(ri)=0;
         zH.row(ri) = Z.row(start)*Ht(ri);
         zHz(ri) = as_scalar((zH.row(ri))*(Z.row(i).st()));
         zz = Z*Xc*Xci*trans(Xc)*trans(Z);
//         if (i==0) break;
     }
     

      colvec gamma(Z.n_cols);   
      mat B = cumsum(b);
      return(Rcpp::List::create(Rcpp::Named("zH")=zH,
                                Rcpp::Named("zHz")=zHz,
                                Rcpp::Named("zz")=zz,
                                Rcpp::Named("B")=B));
  } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
  } catch(...) { 
       ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
'

semiaalen <- cxxfunction(signature(ds="numeric", ts="numeric",
                                   Xs="numeric", Zs="numeric"),
                         src1,plugin="RcppArmadillo", verbose=TRUE)

X <- model.matrix(~1+age,sTRACE)
Z <- model.matrix(~-1+sex+chf,sTRACE)
t <- sTRACE$time
dix <- which(sTRACE$status==9)
dt <- diff(c(0,t[dix]))

system.time(B <- semiaalen(dix,t,X,Z))

colSums(B$zH)/sum(B$zHz*dt)

library(timereg)
library(bsplot)
data(sTRACE)
sTRACE <- sTRACE[order(sTRACE$time),]
system.time(out<-aalen(Surv(time,status==9)~1+const(sex)+const(chf)+age,sTRACE,max.time=9,n.sim=0,robust=0))
b <- coef(out)

system.time(B <- EvilAalen(dix,X))

system.time(out<-aalen(Surv(time,status==9)~1+age,sTRACE,max.time=9,n.sim=0,robust=0))
idx <- 2
plot(out,specific.comps=idx)
points(t[dix],B$B[,idx],type="s",col=Col(2),lwd=7)
