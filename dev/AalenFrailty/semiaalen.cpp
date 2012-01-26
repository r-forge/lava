#include <vector>
using namespace arma;

try {
  colvec event = Rcpp::as<colvec>(ds);
  colvec t = Rcpp::as<colvec>(ts);
  mat    X = Rcpp::as<mat>(Xs);
  mat    Z = Rcpp::as<mat>(Zs);
  unsigned Nevent = sum(event);  
  mat    b1 = zeros(Nevent,X.n_cols);
  mat    X2 = zeros<mat>(X.n_cols,X.n_cols);
  mat    Z2 = zeros<mat>(Z.n_cols,Z.n_cols);
  mat    zHdN(Nevent,Z.n_cols);
  mat    zHz = zeros<mat>(Z.n_cols,Z.n_cols);
  mat    ZX = zeros<mat>(Z.n_cols,X.n_cols);
  bool   all = Rcpp::as<bool>(allt);

  unsigned n = t.n_rows;
  colvec dt(n);  
  // if (!Rf_isNull(dts)) {
  //   dt = Rcpp::as<colvec>(dts);
  // } 
  double t0=0;
  for (unsigned i=0; i<n; i++) {
    //dt(i) = t(event(i))-t0;
    //t0 = t(event(i));
    dt(i) = t(i)-t0;
    t0 = t(i);
  }
  unsigned stop, start = n; 
  unsigned eventpos = Nevent-1;
  
  unsigned i=0;
  double dtlast = 0;
  t0 = t(n-1);
  while (i<n) {
    stop = start-1;
    start = stop;
    unsigned pos = n-i-1;
    double dt0 = dt(pos)+dtlast;
      
    if (!all) {
      bool newevent = (event(pos)==1);
      while (!newevent & pos>=0) {
	i++;
	pos--;
	dt0 += dt(pos);
	newevent = (event(pos)==1);
      }
      start = pos;
    }
    dtlast = dt(pos);
    //    dt0 = t0-t(start); t0 = t(start);
    cerr << start << "." << stop << endl;
    cerr << "dt=" << dt0 << endl;

    mat Xt = X(span(start,stop), span::all);
    mat Zt = Z(span(start,stop), span::all);
    //cerr << "Xt=" << Xt << endl;
    X2 += Xt.st()*Xt;
    Z2 += Zt.st()*Zt;
    ZX += Zt.st()*Xt;
    mat U, V; vec s; 
    svd(U,s,V,X2);
    mat X2i = U*diagmat(1/s)*V.st();
    bool nonSing = true;
    for (unsigned k=0; k<s.n_elem; k++) 
      if (s(k)<1.0e-15) {
	nonSing = false; break;
      }
    if (nonSing) {
      mat ZXp = ZX*X2i;
      mat zHznew = (Z2 - ZXp*ZX.st());
      zHz += zHznew*dt0;
      if (event(pos)==1) {
	//	cerr << "eventpos=" << eventpos << endl;
	b1.row(eventpos) = trans(X2i*trans(Xt.row(0)));
	zHdN.row(eventpos) = Zt.row(0) - trans(ZXp*trans(Xt.row(0)));
	eventpos--;
      }
    } else {
      if (pos==event(eventpos)) eventpos--;
    }
    i++;
  } 
  
  colvec gamma(Z.n_cols);
  if (!Rf_isNull(gammas)) {
    gamma = Rcpp::as<colvec>(gammas);
  } else {
    gamma = inv(zHz)*trans(sum(zHdN));
  }
  mat B = cumsum(b1);

  // mat b2 = zeros<mat>(b1.n_rows,b1.n_cols);
  // X2 *= 0; Z2 *= 0; ZX *= 0;
  // start = X.n_rows;
  // for (unsigned i=0; i<event.n_rows; i++) {
  //   unsigned ri = event.n_rows-i-1;
  //   stop  = start-1;
  //   start = event(ri)-1;
  //   mat Xt = X(span(start,stop), span::all);
  //   mat Zt = Z(span(start,stop), span::all);
  //   X2 += Xt.st()*Xt;
  //   Z2 += Zt.st()*Zt;
  //   ZX += Zt.st()*Xt;
  //   mat U, V; vec s; 
  //   svd(U,s,V,X2);
  //   mat X2i = U*diagmat(1/s)*V.st();   
  //   mat XpZ = (X2i*(ZX.st()));
  //   b2.row(ri) = trans(XpZ*gamma*dt(ri));
  // }
  // mat B = cumsum(b1-b2);

  return(Rcpp::List::create(Rcpp::Named("gamma")=gamma,
			    Rcpp::Named("zHdN")=zHdN,
			    Rcpp::Named("dt")=dt,
			    Rcpp::Named("zHz")=zHz,
			    Rcpp::Named("B")=B));


 } catch( std::exception &ex ) {
  forward_exception_to_r( ex );
 } catch(...) { 
  ::Rf_error( "C++ exception (unknown reason)" ); 
 }
return R_NilValue; // -Wall
