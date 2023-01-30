#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace R;

uword dyx (const uword &i, const uword &ii, const double &I)
{
   // indices used to store values in and call values from Y and Z
   // (vectorized) linear index from a lower triangular matrix position
   return( I*(I - 1)/2 - (I - ii)*(I - ii - 1)/2 + i - ii - 1 );
}

//[[Rcpp::export]]
rowvec get_inmat (const rowvec &Xi) 
{
   // vec lower triangular portion of inmat
   uword I = Xi.n_elem, i, ii;
   rowvec A(I*(I - 1)/2, fill::zeros);
   for (i = 1; i < I; ++i) {
      for (ii = 0; ii < i; ++ii) {
         if ( Xi(i) == Xi(ii) ) A( dyx(i, ii, I) ) = 1;
      }
   }
   return( A );
}

//[[Rcpp::export]]
rowvec get_inmat_hat (const mat &Xi_data)
{
  uword S = Xi_data.n_rows, I = Xi_data.n_cols, s;
  rowvec out(I*(I - 1)/2, fill::zeros);
  for (s = 0; s < S; ++s) out += get_inmat(Xi_data.row(s));
  return( out/(double)S );
}
