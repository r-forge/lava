#include "mvn.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

RcppExport SEXP bicif(SEXP n, SEXP ad, SEXP ii, SEXP causes, 
		      SEXP m, SEXP S, SEXP P) {

  GenericVector AD(ad);
  NumericMatrix Pr(P);
  mat PP(Pr.begin(), Pr.nrow(), Pr.ncol(), false);
  NumericMatrix sigma(S);
  NumericVector Causes(causes);
  GenericVector II(ii);
  //mat T(tt.begin(), tt.nrow(), tt.ncol(), false);
  mat Sigma(sigma.begin(), sigma.nrow(), sigma.ncol(), false);
  //unsigned n = T.n_rows;
  unsigned ncauses = Causes.size()-1;

  //colvec ll(n);
  unsigned pos1=0, pos2=0;
  mat Sigma0(2,2);

  SEXP Idx0;
  CharacterVector save;

  vector<mat> TT;
  for (unsigned i=1;i<=ncauses;i++) {
    SEXP T0 = AD[numStr(Causes[i])];
    NumericMatrix T1(T0);
    mat T2(T1.begin(), T1.nrow(), T1.ncol(), false);
    TT.push_back(T2);
  }
  unsigned N = as<int>(n);//TT[1].n_rows;

  NumericVector ll(N); // log-Likelihood
  
  mat Pmarg1 = sum(PP,1);
  mat Pmarg2 = sum(PP,0);
    
  //  for (unsigned i=0; i<ncauses; i++)
  //    cerr << "P" << i << "=" << Pmarg(i) << endl;
  // cerr << endl;
  // cerr << "PP=" << PP << endl;
  // cerr << "P1=" << Pmarg1;
  // cerr << "P2=" << Pmarg2 << endl;

  // Censored observations: log(S(t1,t2))
  string causestr = "0 0";
  Idx0 = II(causestr);
  NumericVector Idx(Idx0);  
  if (Idx.size()>0) {
    save.push_back(causestr);          
    for (unsigned i=0; i<ncauses; i++) {
      pos1 = i*2;
      //      for (int k=0; k<Idx.size(); k++) ll[Idx[k]-1] = 0;
      double s1 = sqrt(Sigma(pos1,pos1));      
      
      for (int k=0; k<Idx.size(); k++) {
	unsigned idx = Idx[k]-1;
	double a0 = TT[i](idx,0);
	double a1 = TT[i](idx,1);
	//	double da0 = TT[i](idx,2);
	//	double da1 = TT[i](idx,3);
	double F1i = Pmarg1(i)*Rf_pnorm5(a0,0,s1,1,0); //Rcpp::stats::pnorm_0(a0,1,0); //Rf_pnorm5
	double F2i = Pmarg2(i)*Rf_pnorm5(a1,0,s1,1,0); //Rcpp::stats::pnorm_0(a1,1,0);
	ll[idx] = ll[idx] + F1i+F2i;
	// if (k==Idx.size()-1) {
	//   cerr << "Pmarg1[" << i << "]=" << Pmarg1(i) << endl;
	//   cerr << "Pmarg2[" << i << "]=" << Pmarg2(i) << endl;
	//   cerr << "a0=" << a0 <<endl;
	//   cerr << "a1=" << a1 <<endl;
	//   cerr << "s1=" << s1 << endl;
	//   cerr << "F1=" << F1i << endl;
	//   cerr << "F2=" << F2i << endl;
	//   cerr << "P1=" << Rf_pnorm5(a0,0,1,1,0) << endl;
	//   cerr << "P2=" << Rcpp::stats::pnorm_0(a1,1,0) << endl;
	// }
      }
      for (unsigned j=i; j<ncauses; j++) {
	//for (unsigned j=i; j<ncauses; j++) {
	pos2 = j*2+1;
	double s2 = sqrt(Sigma(pos2,pos2));
	double rho = Sigma(pos1,pos2)/(s1*s2);
	double pr = Pr(i,j);
	
	for (int k=0; k<Idx.size(); k++) {
	  unsigned idx = Idx[k]-1;
	  double a0 = (TT[i](idx,0))/s1;
	  double a1 = (TT[j](idx,1))/s2;
	  double Fij = pr*Fbvn(a0,a1,rho); // F(T1<a0,T2<a1)
	  ll[idx] -= Fij;
 	}
      }
    }
    for (int k=0; k<Idx.size(); k++) {
      unsigned idx = Idx[k]-1;
      ll[idx] = log(1-ll[idx]);
    }
  }

  // First censored, second observed
  for (unsigned i=0; i<ncauses; i++) {
    causestr = "0 "+numStr(Causes[i+1]);
    Idx0 = II(causestr);
    NumericVector Idx(Idx0);      
    if (Idx.size()>0) {
      save.push_back(causestr);            
      unsigned pos2 = 2*i+1;
      double s2 = Sigma(pos2,pos2); // Marginal variance of cause 
      //      double val0 = log(p2) -0.5*(log2pi+log(s2));      

      for (int k=0; k<Idx.size(); k++) {
	unsigned idx = Idx[k]-1;
	double z = TT[i](idx,1);
	double logdz = TT[i](idx,3);
	// += log(a2'(t2)) - log(f(a2(t2)))
	//	ll[idx] = dz + val0 - 0.5*z*z/s2;    
	ll[idx] = logdz + log(Pmarg2[i]) + Rf_dnorm4(z,0.0,sqrt(s2),1);
      
	double val = 0;
	for (unsigned j=0; j<ncauses; j++) {
	  double pos1 = 2*j;
	  double z1 = TT[j](idx,0);
	  double cov12 = Sigma(pos1,pos2);	  
	  double condvar = Sigma(pos1,pos1)-cov12*cov12/s2;
	  double condmean = cov12/s2*z;
	  val += Pr(j,i)/Pmarg2[i]*Rf_pnorm5(z1, condmean, sqrt(condvar), 1, 0); // lower tail, not log
	}
	ll[idx] += log(1-val);
      }
    }
  }

  // Uncensored observations
  for (unsigned i=1; i<=ncauses; i++) {
    for (unsigned j=i; j<=ncauses; j++) {    
      causestr = numStr(Causes[i])+" "+numStr(Causes[j]);
      Idx0 = II(causestr);
      NumericVector Idx1(Idx0);  

      if (Idx1.size()>0) {
	save.push_back(causestr);            
	pos1 = (i-1)*2;
	pos2 = (j-1)*2+1;
	Sigma0(0,0) = Sigma(pos1,pos1);
	Sigma0(1,1) = Sigma(pos2,pos2);
	Sigma0(0,1) = Sigma0(1,0) = Sigma(pos1,pos2);
	mat iSigma0 = inv(Sigma0);
	double detS; double sds;
	log_det(detS,sds,Sigma0);
	detS = -log2pi-0.5*detS;
	double logpr = log(Pr(i-1,j-1));

	colvec a(2);    
	for (int k=0; k<Idx1.size(); k++) {
	  unsigned idx = Idx1[k]-1;
	  a(0) = TT[i-1](idx,0);
	  a(1) = TT[j-1](idx,1);
	  double logd1 = TT[i-1](idx,2);
	  double logd2 = TT[j-1](idx,3);
	  ll[idx] = logd1+logd2+logpr +
	    detS-0.5*as_scalar(trans(a)*iSigma0*a);	
	}
      }
    }
  }

  List res;
  res["logLik"] = ll;
  res["save"] = save;
  res["Sigma"] = Sigma;
  res["causes"] = Causes;
  res["ncauses"] = ncauses;
  return(res);
}

