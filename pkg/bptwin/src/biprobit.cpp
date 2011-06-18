#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>

using namespace std;
using namespace Rcpp;
using namespace arma;

extern "C" double bvnd_(const double *dh, const double *dk, const double *r);

struct vecmat
{
  vec V;
  mat M;
};

RcppExport SEXP cdf(SEXP a, SEXP b, SEXP r) { 
  double u1=-Rcpp::as<double>(a);
  double u2=-Rcpp::as<double>(b);
  double cr=Rcpp::as<double>(r);  
  double val = bvnd_(&u1, &u2, &cr);
  NumericVector res(1); res[0] = val;
  return(res);
}

double dbvnorm(double y1, double y2, double R) {
  double detR = 1-R*R;
  // inv(R) = [1 -r; -r 1]/detR (prove by gauss elim.)
  double res = 1/(2*M_PI*sqrt(detR))*exp(-0.5/detR*(y1*y1+y2*y2-2*R*y1*y2));
  return(res);
}

vecmat Dbvn(double y1, double y2, double R) {
  vec DP(2);
  double R2 = R*R;
  DP(0) = Rf_dnorm4(y1,0.0,1.0,0)*Rf_pnorm5(y2-R*y1,0.0,sqrt(1-R2),1,0);
  DP(1) = Rf_dnorm4(y2,0.0,1.0,0)*Rf_pnorm5(y1-R*y2,0.0,sqrt(1-R2),1,0);
  mat HP(2,2);
  HP(1,0) = HP(0,1) = dbvnorm(y1,y2,R);
  HP(0,0) = -y1*DP(0) - R*HP(1,0);
  HP(1,1) = -y2*DP(1) - R*HP(1,0);  
  vecmat res;
  res.V = DP;
  res.M= HP;
  return(res);
}

double bvnd(double &u1, double &u2,double &r) {     
  double val = bvnd_(&u1, &u2, &r);
  return(val);
}

mat logLik(mat &y, mat &mu, mat &Sigma) {
  unsigned n = y.n_rows;
  colvec res(n);
  double l = Sigma(0,0);
  double R = Sigma(0,1)/l;
  for (unsigned i=0;i<n;i++) {
    double mysign = 1;
    rowvec lo = mu.row(i)/sqrt(l);
    if (y(i,0)==1)  lo(0) *= -1;
    if (y(i,1)==1) lo(1) *= -1;
    if (y(i,0)!=y(i,1)) mysign = -1;
    double R0 = R*mysign;
    res(i) = log(bvnd(lo(0),lo(1),R0)); //bvnd calculates lower tail (mult.by -1).
  }
  return(res);
}

vecmat score(mat &y, mat &x, mat &mu, mat &Sigma, mat &dS0, mat &w) {
  unsigned n = y.n_rows;
  mat iSigma = inv(Sigma);  
  double l = Sigma(0,0);
  vec ll = sqrt(diagvec(Sigma)); 
  mat Lambda = diagmat(ll);
  mat iLambda = diagmat(1/ll);
  mat R = iLambda*Sigma*iLambda;
  mat R0 = R;
  double r = Sigma(0,1)/l;
  
  mat LR = Lambda*R;
  mat LR0 = LR;
  mat S= Sigma;
  mat iS = iSigma;
  mat iSxiS1 = kron(iS,iS);
  iS(1,0) = iS(0,1) = -iSigma(1,0);
  mat iSxiS2 = kron(iS,iS);

  vec Lik(n);
  mat Score(n,x.n_cols/2+dS0.n_rows);
  mat U;
  for (unsigned i=0;i<n;i++) {
    int mysign = 1;
    double r0 = r;
    rowvec lo = mu.row(i)%(1/trans(ll));
    mat L = eye(2,2);
    if (y(i,0)==1) { lo(0) *= -1; L(0,0)= -1; } 
    if (y(i,1)==1) { lo(1) *= -1; L(1,1)= -1; }
    mat dS = dS0;    
    mat iSxiS = iSxiS1;    
    if (y(i,0)!=y(i,1)) { 
      dS.col(1) *= -1; dS.col(2) *= -1;
      iSxiS = iSxiS2;
      mysign = -1;
    } 
    S = Sigma;
    S(1,0) = mysign*Sigma(1,0); S(0,1) = S(1,0);
    iS = iSigma;
    iS(1,0) = mysign*iSigma(1,0); iS(0,1) = iS(1,0);
    LR = LR0;
    LR(0,1) = mysign*LR0(1,0); LR(1,0) = LR(0,1);
    r0 *= mysign;
    rowvec up = -lo;
    double alpha = bvnd(lo(0),lo(1),r0);    
    vecmat D = Dbvn(up(0),up(1),r0);
    mat M = -LR*D.V;    
    mat V = alpha*S + LR*D.M*trans(LR);
    if (!(w.is_empty())) {  
      mat W = diagmat(w.row(i));
      V = W*V;     
      iS = W*iS;
    }
    mat vecV = reshape(V,4,1);
    mat veciS = reshape(iS,4,1);
    mat dmu = L*trans(reshape(x.row(i),x.n_cols/2,2));     
    mat d1 = trans(dmu)*iS*M/alpha;
    mat d2 = -0.5*(dS)*(veciS - iSxiS*vecV/alpha);
    rowvec val = trans(join_cols(d1,d2));
    Score.row(i) = val;
    Lik(i) = alpha;
  }
  vecmat res; 
  res.V = Lik; res.M = Score;
  return(res);
}

RcppExport SEXP biprobit0(SEXP m, 
			  SEXP s, SEXP ds, SEXP y, SEXP x, SEXP w, SEXP std) {

  NumericMatrix mm(m); 
  NumericMatrix ss(s); 
  NumericMatrix dss(ds); 
  NumericMatrix xx(x);
  NumericMatrix yy(y);

  mat mu(mm.begin(), mm.nrow(), mm.ncol(), false);
  mat S(ss.begin(), ss.nrow(), ss.ncol(), false);
  mat dS(dss.begin(), dss.nrow(), dss.ncol(), false);
  mat X(xx.begin(), xx.nrow(), xx.ncol(), false);
  mat Y(yy.begin(), yy.nrow(), yy.ncol(), false);
  bool weights = Rcpp::as<bool>(std);
  //  bool efree = Rcpp::as<bool>(e); 

  mat W;
  if (weights) {
    NumericMatrix ww(w);
    mat W0(ww.begin(), ww.nrow(), ww.ncol(), false);
    W = W0;
   //    WW = trans(reshape(W1,2,N/2));    
  }
  vecmat U = score(Y,X,mu,S,dS,W);

  List res;
  res["loglik"] = log(U.V);
  res["score"] = U.M;
  return(res);
}


RcppExport SEXP FastApprox(const SEXP a,
			   const SEXP t,
			   const SEXP z) {
  NumericVector A(a);
  NumericVector T(t);
  NumericVector Z(z);
  vector<unsigned> idx(Z.size());
  vector<double> newT(Z.size());

  NumericVector::iterator it;  
  double lower,upper; int pos=200;
  for (int i=0; i<Z.size(); i++) {
    it = lower_bound(A.begin(), A.end(), Z[i]);
    if (it == A.begin()) { 
      pos = 0; 
      // upper = *it; // no smaller value  than val in vector
    } 
    else if (int(it-A.end())==0) {
      pos = A.size()-1;
      //lower = *(it-1); // no bigger value than val in vector
    } else {
      lower = *(it-1);
      upper = *it;
      pos = int(it- A.begin());
      if (abs(Z[i]-lower)<abs(Z[i]-upper)) {
	pos = int(it- A.begin());
      }
    }
    idx[i] = pos;
    newT[i] = T[pos];
  }
  List ans;
  ans["t"] = newT;
  ans["pos"] = idx;
  return(ans);
}
