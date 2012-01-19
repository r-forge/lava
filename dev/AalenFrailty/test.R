  31 : using namespace arma;
  32 : try {
  33 :       colvec   tidx = Rcpp::as<colvec>(its); 
  34 :       mat         X = Rcpp::as<mat>(Xs);
  35 : 
  36 :       mat b(tidx.n_rows,X.n_cols);
  37 :       mat Xc = zeros<mat>(X.n_cols,X.n_cols);
  38 :       unsigned stop,start = X.n_rows;
  39 :       for (unsigned i=0; i<tidx.n_rows; i++) {
  40 :          stop  = start-1;
  41 :          start = tidx(i)-1;
  42 : //std::cerr << "start=" << start << endl;
  43 : //std::cerr << "stop=" << stop << endl;
  44 :          mat Xi = X(span(start,stop), span::all);
  45 :          Xc = Xc + Xi.st()*Xi;
  46 :          mat U, V; vec s; 
  47 :          svd(U,s,V,Xc);
  48 :          mat Xci = U*diagmat(1/s)*V.st();      
  49 :          b.row(tidx.n_rows-i-1) = trans(Xci*trans(X.row(tidx(i)-1)));
  50 :       }
  51 :       mat B = cumsum(b);
  52 :       return(Rcpp::List::create(Rcpp::Named("B")=B));
  53 : 
  54 : 
  55 :   } catch( std::exception &ex ) {
  56 :         forward_exception_to_r( ex );
  57 :   } catch(...) { 
  58 :        ::Rf_error( "c++ exception (unknown reason)" ); 
  59 :   }
  60 :   return R_NilValue; // -Wall
  61 : 
  62 : END_RCPP
  63 : }


[254,] -2.167428411 0.0429496545  1.786345e-01
[255,] -2.195046042 0.0433648290  1.935558e-01
[256,] -2.224391035 0.0443965262  1.604049e-01
[257,] -2.294743655 0.0460745324  1.280296e-01
[258,] -2.214247078 0.0456965319  7.042055e-02
[259,] -2.173538148 0.0450621968  1.206016e-01
