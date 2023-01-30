#ifndef samplers_H
#define samplers_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace R;

uword wsample (vec probs) 
{
     // weighted sampling according to probs
     // return indices fron 0 to n - 1, where n = probs.n_elem
     probs = probs/accu(probs);
     double u = R::runif(0.0, 1.0);
     if (u < probs(0)) {
          return 0;
     } else {
          uword n = probs.n_elem, i, out = 0;
          vec probsum = cumsum(probs);
          for (i = 1; i < n; ++i) {  
               if ( (probsum(i-1) <= u) && (u < probsum(i)) ) {
                    out = i;
                    goto endloop;
               }
          }
          endloop:
               return( out );
     }
}

double myround (const double &value, const uword &digits)
{
     // otherwise it will return 'nan' due to the log10() of zero
     if (value == 0.0) {
          return 0.0;
     } else {
          double factor = pow(10.0, digits - ceil(log10(fabs(value))));
          return( round(value * factor)/factor );
     }
}

char* mypaste0 (string path, string name)
{
     stringstream strname;
     strname << path << name << ".txt";
     string fullname = strname.str();
     string::iterator p = fullname.begin();
     char* chr = &(*p);
     return( chr );
}

vec mytable (const uvec &x)
{
     // x has continuous labels from 1 to N 
     uword N = max(x);
     vec out(N);
     for (uword n = 0; n < N; ++n) {
          uvec id = find(x == (n+1));
          out(n) = id.n_elem;
     }
     return( out );
}

vec myallelic (const uvec &x)
{
     uword I = x.n_rows;
     vec nk = mytable(x), out(I);
     for (uword i = 0; i < I; ++i) {
          uvec id = find(nk == (i+1));
          out(i) = id.n_elem;
     }
     return( out );
}

void ijx (uword &i, uword &j, const uword &ip, const uword &Ij0)
{
     // return record i and file j indices given a linear index ip
     if (ip < Ij0) {
          j = 0;
          i = ip;
     } else {
          j = 1;
          i = ip - Ij0;
     }
}

uword dyx (const uword &i, const uword &ii, const double &I)
{
     // indices used to store values in and call values from Y and Z
     // (vectorized) linear index from a lower triangular matrix position
     return( I*(I - 1)/2 - (I - ii)*(I - ii - 1)/2 + i - ii - 1 );
}

uword get_str_position(const std::string& str_to_find, const std::vector<std::string>& vector_to_search)
{
     std::vector<std::string>::const_iterator it = std::find(vector_to_search.begin(), vector_to_search.end(), str_to_find);
     if (it == vector_to_search.end()) {
          return -1;
     } else {
          return it - vector_to_search.begin();
     }
}
/*
double LDist (const string &s, const string &t) {
     // LevenshteinDistance
     // Number of elements
     double n = s.size();
     double m = t.size();
     Rcpp::IntegerMatrix d(n+1, m+1);
     //Rcpp::Rcout << "n:" << n << ", m:" << m << std::endl;
     if (n == 0){
          return m;
     }
     if (m == 0){
          return n;
     }
     for (int i = 0; i <= n; i++){
          d(i, 0) = i;
     }
     // No sense to revisit the (0,0) coordinate
     for (int j = 1; j <= m; j++){
          d(0, j) = j;
     }
     for (int j = 1; j <= m; j++){
          for (int i = 1; i <= n; i++){
               if (s[i - 1] == t[j - 1]){
                    d(i, j) = d(i - 1, j - 1);  // no operation
               } else {
                    d(i, j) = std::min(d(i - 1, j) + 1,    //a deletion
                      std::min(d(i, j - 1) + 1,    //an insertion
                               d(i - 1, j - 1) + 1));       //a substitution
               } // end if
          } // end inner for
     } // end outer for
     return d(n, m);
}
*/

double logH (const std::string& w, const vec &H, const vector<std::string> &SS)
{
     return( log(H( get_str_position(w, SS) )) );
}

double logecd (const std::string& w1, const std::string& w2, const double &nS, const vec &ecd, const vector<std::string> &SS)
{
     uword w1_idx = get_str_position(w1, SS);
     uword w2_idx = get_str_position(w2, SS);
     if (w2_idx == w1_idx) {
          return( 0.0 );
     } else {
          if (w2_idx < w1_idx)
               return( log(ecd( dyx(w1_idx, w2_idx, nS) )) );
          else
               return( log(ecd( dyx(w2_idx, w1_idx, nS) )) );
     }
}

double ll_iter_pro (const uvec &Ij, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const uvec &W, const uvec &La, const uvec &Xi, const vector<std::string> &Pi, const vector<std::string> &SS)
{
     // computes the log-likelihood in a given iteration of the MCMC
     double out = 0.0;
     for (uword ip = 0; ip < accu(Ij); ++ip) {
          if (W(ip) == 1) {
               out += exp( log(Alpha(La(ip) - 1)) + logecd(SS[La(ip) - 1], Pi[Xi(ip) - 1], nS, ecd, SS) + logH(Pi[Xi(ip) - 1], H, SS) );
          } else {
               if ( SS[La(ip) - 1] != Pi[Xi(ip) - 1] ) {
                    out = -datum::inf;
                    goto endloop;
               }
          }
     }
     endloop:
          return( out );
}

double ll_iter_pro_ij (const uword &ip, const uvec &Ij, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const uvec &W, const uvec &La, const uvec &Xi, const vector<std::string> &Pi, const vector<std::string> &SS)
{
     // computes the log-likelihood in a given iteration of the MCMC
     double out = 0.0;
     if (W(ip) == 1) {
          out += exp( log(Alpha(La(ip) - 1)) + logecd(SS[La(ip) - 1], Pi[Xi(ip) - 1], nS, ecd, SS) + logH(Pi[Xi(ip) - 1], H, SS) );
     } else {
          if ( SS[La(ip) - 1] != Pi[Xi(ip) - 1] ) {
               out = -datum::inf;
               goto endloop;
          }
     }
     endloop:
          return( out );
}

double ll_iter_net (const uvec &Ij, const vec &Beta, const mat &U, const uvec &Xi, const uvec &Y)
{
     // computes the log-likelihood in a given iteration of the MCMC
     uword j, i, ii, l, ll;
     uvec Dj(2, fill::zeros); Dj(0) = Ij(0)*(Ij(0) - 1)/2;
     double eta, out = 0.0;
     for (j = 0; j < 2; ++j) {
          l  = accu(Ij.rows(0, j)) - Ij(j);
          ll = accu(Dj.rows(0, j)) - Dj(j);
          for (i = 1; i < Ij(j); ++i) {  // lower triangualar matrix
               for (ii = 0; ii < i; ++ii) {
                    eta = Beta(j) - norm( U.row( Xi(i + l) - 1 ) - U.row( Xi(ii + l) - 1 ) );
                    out += eta * Y( dyx(i, ii, Ij(j)) + ll ) - log(1.0 + exp(eta));
               }
          }
     }
     return( out );
}

double ll_iter_net_ij (const uword &ip, const uvec &Ij, const vec &Beta, const mat &U, const uvec &Xi, const uvec &Y)
{
     // computes the contribution of the log-likelihood for a given i and j
     uword i, j;
     ijx(i, j, ip, Ij(0));
     uvec Dj(2, fill::zeros); Dj(0) = Ij(0)*(Ij(0) - 1)/2;
     uword l = accu(Ij.rows(0, j)) - Ij(j), ll = accu(Dj.rows(0, j)) - Dj(j), ii;
     double eta, out = 0.0;
     for (ii = 0; ii < i; ++ii) {
          eta = Beta(j) - norm( U.row( Xi(i + l) - 1 ) - U.row( Xi(ii + l) - 1 ) );
          out += eta * Y( dyx(i, ii, Ij(j)) + ll ) - log(1.0 + exp(eta));
     }
     for (ii = i+1; ii < Ij(j); ++ii) {
          eta = Beta(j) - norm( U.row( Xi(ii + l) - 1 ) - U.row( Xi(i + l) - 1 ) );
          out += eta * Y( dyx(ii, i, Ij(j)) + ll ) - log(1.0 + exp(eta));
     }
     return( out );
}

void tunning (double &delta, const double &mr0, const double &epsilon, const double &mix_rate, const double &s)
{
     // tunning paramter calibration
     // ntun = 100
     if ( ((uword)s % 100 == 0) && (abs(mix_rate - mr0) > epsilon) ) {
          double tmp = delta, cont = 1.0;
          do {
               tmp = delta + (0.1/cont) * (mix_rate - mr0);
               ++cont;
          } while ( (tmp <= 0.0) && (cont <= 100.0) );
          if (tmp > 0) delta = tmp;
     }
}

void relabel (uvec &Xi) 
{
     // re-label latent individuals from 1 to N
     // N is the number of unique elements in Xi
     uvec uXi = unique(Xi);  // unique elements of Xi, sorted in ascending order
     uword N = uXi.n_elem;
     if (N < uXi(N-1)) {
          for (uword n = 0; n < N; ++n) {
               uvec id = find(Xi == uXi(n));
               for (uword i = 0; i < id.n_elem; ++i) Xi(id(i)) = n+1;
          }
     }
}

void sample_Psi (const uvec &Ij, const double &J, const double &a_psi, const double &b_psi, vec &Psi, const uvec &W)
{
     double Suma;
     for (unsigned j = 0; j < J; ++j) {
          Suma = accu(W.rows( accu(Ij.rows(0, j)) - Ij(j), accu(Ij.rows(0, j)) - 1 ));
          Psi(j) = R::rbeta(a_psi + Suma, b_psi + Ij(j) - Suma);
     }
}

void sample_W_i (const uword ip, const uword &j, const uvec &Ij, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const vec &Psi, uvec &W, const uvec &La, const uvec &Xi, const vector<std::string> &Pi, const vector<std::string> &SS)
{
     W(ip) = 1;
     if ( SS[La(ip) - 1] == Pi[Xi(ip) - 1] ) {
          double logprob = log(Psi(j)) + log(Alpha(La(ip) - 1)) + logH(Pi[Xi(ip) - 1], H, SS) + logecd(SS[La(ip) - 1], Pi[Xi(ip) - 1], nS, ecd, SS);
          logprob -= log(exp(logprob) + 1.0 - Psi(j));
          W(ip) = R::rbinom(1.0, exp(logprob));
     }
}

void sample_W (const uvec &Ij, const double &J, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const vec &Psi, uvec &W, const uvec &La, const uvec &Xi, const vector<std::string> &Pi, const vector<std::string> &SS)
{
     
     uword i, j, ip;
     for (ip = 0; ip < accu(Ij); ++ip) {
          ijx(i, j, ip, Ij(0));  // indices
          W(ip) = 1;
          if ( SS[La(ip) - 1] == Pi[Xi(ip) - 1] ) {
               double logprob = log(Psi(j)) + log(Alpha(La(ip) - 1)) + logH(Pi[Xi(ip) - 1], H, SS) + logecd(SS[La(ip) - 1], Pi[Xi(ip) - 1], nS, ecd, SS);
               logprob -= log(exp(logprob) + 1.0 - Psi(j));
               W(ip) = R::rbinom(1.0, exp(logprob));
          }
     }
}

vec get_Phi ( const uvec &Rnj, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const uvec &W, const uvec &La, const vector<std::string> &SS)
{
     vec Phi(nS, fill::zeros);
     for (uword w = 0; w < nS; ++w) {
          double Suma = 0.0;
          for (uword r = 0; r < Rnj.n_elem; ++r) Suma += W(Rnj(r)) * ( log(H( w )) + logecd(SS[La(Rnj(r)) - 1], SS[w], nS, ecd, SS) );
          Phi(w) = exp( log(Alpha(w)) + Suma );
     }
     return( Phi );
}

void sample_Pi (const uvec &Ij, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const uvec &W, const uvec &La, const uvec &Xi, vector<std::string> &Pi, const vector<std::string> &SS)
{
     uword n, N = max(Xi);
     for (n = 0; n < N; ++n) {
          uvec Rnj = find(Xi == (n+1));
          uvec idx = find(W(Rnj) == 0);
          if (idx.n_elem > 0) {
               // If an undistorted X.ijl is linked to this Y.j'l, then Y.j'l=X.ijl
               Pi[n] = SS[La(Rnj(idx(0))) - 1];
          } else {
               // If all X.ijl linked to this Y.j'l are distorted, then we draw Y.j'l from a distribution
               // get empirical distribution
               vec Phi = get_Phi(Rnj, nS, H, ecd, Alpha, W, La, SS);
               Pi[n] = SS[wsample(Phi)];
          }
     }
}

double logcond_U (const rowvec &u, const uword &n, const uvec &id, const uvec &Ij, const double &K, const double &sigsq, const vec &Beta, mat U, const uvec &Xi, const uvec &Y) 
{
     U.row(n) = u;
     double out = 0.0;
     for (uword k = 0; k < id.n_elem; ++k) out += ll_iter_net_ij(id(k), Ij, Beta, U, Xi, Y);
     for (uword k = 0; k < K; ++k) out += R::dnorm(U(n, k), 0.0, sqrt(sigsq), 1);
     return( out );
}

void sample_U (const double &s, double &nstar_N, double &nstar_U, double &delta_U, const uvec &Ij, const double &K, const double &sigsq, const vec &Beta, mat &U, const uvec &Xi, const uvec &Y)
{
     // Xi = [\xi_1, ... , \xi_I] is indexed from 1 to N where N = max{Xi}
     // U  = [u*_{\xi_1}, ... , u*_{\xi_I}]^T : I x K matrix
     // U* = [u*_1, ... , u*_N]^T             : N x K matrix
     // Metropolis step
     uword N = max(Xi), n;
     rowvec u_p(K);
     double r;
     for (n = 0; n < N; ++n) {
          uvec id = find( Xi == (n+1) );
          u_p = U.row(n) + delta_U * randn<rowvec>(K);
          r = exp( logcond_U(u_p, n, id, Ij, K, sigsq, Beta, U, Xi, Y) - logcond_U(U.row(n), n, id, Ij, K, sigsq, Beta, U, Xi, Y) );
          if (R::runif(0.0, 1.0) < r) {
               U.row(n) = u_p;
               ++nstar_U;
          }
     }
     nstar_N += N;
     tunning(delta_U, 0.37, 0.025, nstar_U/nstar_N, s);
}

double logcond_Beta (const uvec &Ij, const double &omesq, const vec &Beta, const mat &U, const uvec &Xi, const uvec &Y)
{
     return( ll_iter_net(Ij, Beta, U, Xi, Y) + R::dnorm(Beta(0), 0.0, sqrt(omesq), 1) + R::dnorm(Beta(1), 0.0, sqrt(omesq), 1) );
}

void sample_Beta (const double &s, double &nstar_Beta, double &delta_Beta, const uvec &Ij, const double &omesq, vec &Beta, const mat &U, const uvec &Xi, const uvec &Y)
{
     // randrom draw from the fcd of Beta
     // Beta = [\Beta_1, \Beta_2]^T (2 x 1 vector)
     vec Beta_p = Beta + delta_Beta * randn<vec>(2); 
     double r = exp( logcond_Beta(Ij, omesq, Beta_p, U, Xi, Y) - logcond_Beta(Ij, omesq, Beta, U, Xi, Y) );
     if (R::runif(0.0, 1.0) < r) {
          Beta = Beta_p; 
          ++nstar_Beta;
     }
     tunning(delta_Beta, 0.37, 0.025, nstar_Beta/s, s);
}

void sample_sigsq (const double &K, const double &a_sig, const double &b_sig, double &sigsq, const mat &U, const uvec &Xi)
{
     double N = max(Xi);
     sigsq = 1.0/R::rgamma( a_sig + 0.5 * K * N, 1.0/(b_sig + 0.5 * accu(pow(U, 2.0))) ) ;
}

void get_val_two (const uword &i, const uword &j, const uword &ip, const uvec &Ij, const uword &newlab, uvec &val, const uvec &W, const uvec &La, const uvec &Xi, const vector<std::string> &Pi, const vector<std::string> &SS, vec &probs)
{
     // range of values for \Xi_{i,j}
     // takes from "val" all the elements of "ref" out
     // val = [\Xi_{1,j'},...,\Xi_{I_j',j'}]^T
     // ref = [\Xi_{1,j},...,\Xi_{i-1,j},\Xi_{i+1,j},...,\Xi_{I_j,j}]^T
     uvec ref = Xi.rows( accu(Ij.rows(0, j)) - Ij(j), accu(Ij.rows(0, j)) - 1 );
     ref.shed_row(i);
     uword nval0 = val.n_elem, c = 0;
     for (uword k = 0; k < nval0; ++k) {
          if ( any(ref == val(k - c)) ) {
               val.shed_row(k - c);
               probs.shed_row(k - c);
               ++c;
          }
     }
     // if w_ijl = 0 and x_ijl != y_cl, then Xi_ij = c is impossible
     nval0 = val.n_elem, c = 0;
     for (uword k = 0; k < nval0; ++k) {
          if (( W(ip) == 0 ) && ( SS[La(ip) - 1] != Pi[val(k - c) - 1] )) {
               val.shed_row(k - c);
               probs.shed_row(k - c);
               ++c;
          }
     }
     val.insert_rows(val.n_elem, 1);
     val(val.n_elem - 1) = newlab;
     probs.insert_rows(probs.n_elem, 1);
     probs(probs.n_elem - 1) = 1.0 - accu(probs);
}

void get_val_two_net (const uword &i, const uword &j, const uvec &Ij, const uword &newlab, uvec &val, const uvec &Xi, vec &probs)
{
     // range of values for \Xi_{i,j}
     // takes from "val" all the elements of "ref" out
     // val = [\Xi_{1,j'},...,\Xi_{I_j',j'}]^T
     // ref = [\Xi_{1,j},...,\Xi_{i-1,j},\Xi_{i+1,j},...,\Xi_{I_j,j}]^T
     uvec ref = Xi.rows( accu(Ij.rows(0, j)) - Ij(j), accu(Ij.rows(0, j)) - 1 );
     ref.shed_row(i);
     uword k, nval0 = val.n_elem, c = 0;
     for (k = 0; k < nval0; ++k) {
          if ( any(ref == val(k - c)) ) {
               val.shed_row(k - c);
               probs.shed_row(k - c);
               ++c;
          }
     }
     val.insert_rows(val.n_elem, 1);
     val(val.n_elem - 1) = newlab;
     probs.insert_rows(probs.n_elem, 1);
     probs(probs.n_elem - 1) = 1.0 - accu(probs);
}

void prior_int_pro (const uvec &Ij, const double &J, const double &a_psi, const double &b_psi, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, vec &Psi, uvec &W, const uvec &La, const uvec &Xi, vector<std::string> &Pi, const vector<std::string> &SS)
{
     // parameter initialization
     uword N = max(Xi), n, i, j, ip;
     vec probs(SS.size(), fill::ones);
     for (n = 0; n < N; ++n) Pi[n] = SS[wsample(probs)];
     for (j = 0; j < J; ++j) {
          Psi(j) = R::rbeta(a_psi, b_psi);
          for (i = 0; i < Ij(j); ++i) {
               if (j == 0) { ip = i; } else { ip = i + Ij(0); }  // linear index
               sample_W_i(ip, j, Ij, nS, H, ecd, Alpha, Psi, W, La, Xi, Pi, SS);
          }
     }
}

void prior_init_net (const double &J, const double &a_sig, const double &b_sig, double &sigsq, const double &omesq, vec &Beta, mat &U)
{
     sigsq = 1.0/R::rgamma(a_sig, 1.0/b_sig);
     for (uword j = 0; j < J; ++j) Beta(j) = R::rnorm(0.0, sqrt(omesq));
     for (uword i = 0; i < U.n_rows; ++i) {
          for (uword k = 0; k < U.n_cols; ++k) {
               U(i, k) = R::rnorm(0.0, sqrt(sigsq));
          }
     }
}

rowvec get_inmat (const uvec &Xi) 
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

#endif