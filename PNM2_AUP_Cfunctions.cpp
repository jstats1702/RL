#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "samplers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace R;

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const string &pi_new, const rowvec &u_new, const uword &ip, const uvec &Ij, const double &K, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const uvec &W, const uvec &La, vector<std::string> Pi, const vector<std::string> &SS, const double &sigsq, const vec &Beta, mat U, uvec Xi, const uvec &Y)
{
     // computes the log conditional of Xi
     Xi(ip) = Xi_i;
     relabel(Xi);
     if ((caso == 1) && (Xi_i != curlab)) {
          Pi.erase( Pi.begin() + curlab-1 );
          U.shed_row(curlab-1);
     }
     if (Xi_i == newlab) {
          Pi.push_back(pi_new);
          U = join_vert(U, u_new);
     }
     double out = ll_iter_pro_ij(ip, Ij, nS, H, ecd, Alpha, W, La, Xi, Pi, SS) + ll_iter_net_ij(ip, Ij, Beta, U, Xi, Y);
     out += log(Alpha( get_str_position(Pi[Xi(ip) - 1], SS) ));
     for (uword k = 0; k < K; ++k) out += R::dnorm(U(Xi(ip)-1, k), 0.0, sqrt(sigsq), 1);
     return( out );
}

void sample_Xi (const uvec &Ij, const double &K, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const vec &Psi, uvec &W, const uvec &La, vector<std::string> &Pi, const vector<std::string> &SS, const double &sigsq, const vec &Beta, mat &U, uvec &Xi, const uvec &Y, const mat &Simat, const uvec &Samip)
{
     // random draw Xi, modifying U accordingly
     uword ip, i, j, curlab, newlab, prolab, caso, newval;
     double acc;
     string pi_new;
     rowvec u_new(K);
     //for (ip = 0; ip < accu(Ij); ++ip) {
     for (uword n = 0; n < Samip.n_elem; ++n) {
          ip = Samip(n) - 1;
          ijx(i, j, ip, Ij(0));
          curlab = Xi(ip);       // current label
          newlab = max(Xi) + 1;  // new label
          pi_new = SS[wsample(Alpha)];
          u_new  = sqrt(sigsq) * randn<rowvec>(K);
          //caso 
          uvec id = find(Xi == curlab);  // id: actor indices pointing to curlab 
          caso = id.n_elem;  // i is currelntly unmatched (caso = 1) or matched (caso = 2)
          // support newval and probabilities
          uvec val = Xi.rows( accu(Ij.rows(0, 1-j)) - Ij(1-j), accu(Ij.rows(0, 1-j)) - 1 );
          vec probs(Ij(1-j));
          if (j == 0)
               probs = Simat.col(i);
          else 
               probs = Simat.row(i).t();
          probs = probs/accu(probs);
          get_val_two(i, j, ip, Ij, newlab, val, W, La, Xi, Pi, SS, probs);
          // Draw proposal
          prolab = val(wsample(probs));  // proposal
          // Compute acceptance rate          
          acc = exp( logcond_Xi(prolab, caso, curlab, newlab, pi_new, u_new, ip, Ij, K, nS, H, ecd, Alpha, W, La, Pi, SS, sigsq, Beta, U, Xi, Y) - logcond_Xi(curlab, caso, curlab, newlab, pi_new, u_new, ip, Ij, K, nS, H, ecd, Alpha, W, La, Pi, SS, sigsq, Beta, U, Xi, Y) );
          // Set new label
          if (R::runif(0.0, 1.0) < acc)
               newval = prolab;
          else
               newval = curlab;
          // modify Xi and U* accordingly
          Xi(ip) = newval;
          relabel(Xi);
          if (caso == 1 && newval != curlab) { 
               Pi.erase( Pi.begin() + curlab-1 );
               U.shed_row(curlab-1);
          }
          if (newval == newlab){ 
               Pi.push_back(pi_new);
               sample_W_i(ip, j, Ij, nS, H, ecd, Alpha, Psi, W, La, Xi, Pi, SS);
               U = join_vert(U, u_new);
          }
     }
}

//[[Rcpp::export]]
void MCMC (const double &a_psi, const double &b_psi, uvec Xi, const uvec &Ij, const vec &H, const vec &ecd, const vec &Alpha, const uvec &La, const vector<std::string> &SS, const double &K, const uvec &Y, const mat &Simat, const uword &nburn, const uword &nsams, const uword &nskip, const uword &ndisp, string path, const uvec &Samip)
{
     // some quantities
     double I  = accu(Ij);
     double J  = Ij.n_elem;
     double nS = SS.size();
     uvec Dj(J, fill::zeros); Dj(0) = Ij(0)*(Ij(0) - 1)/2;
     // hyperparameters
     double a_sig = 2.0 + pow(0.5, -2.0); // CV0 = 0.5
     double b_sig = (a_sig - 1.0) * ( sqrt(I)/(sqrt(I) - 2.0) ) * ( ( pow(datum::pi, K/2.0) / exp( lgamma(K/2.0 + 1.0) ) ) * pow(I, 2.0/K) );
     double omesq = 100.0;
     // parameter initialization
     vec  Psi(J);
     uvec W(I);
     vector<std::string> Pi(max(Xi));
     prior_int_pro (Ij, J, a_psi, b_psi, nS, H, ecd, Alpha, Psi, W, La, Xi, Pi, SS);
     double sigsq;
     vec Beta(2);
     mat U(max(Xi), K);
     prior_init_net(J, a_sig, b_sig, sigsq, omesq, Beta, U);
     // Metropolis
     double nstar_N = 0.0;
     double nstar_U = 0.0;
     double delta_U = 2.38/sqrt(K);
     double nstar_Beta = 0.0;
     double delta_Beta = 2.38/sqrt(2.0);
     // write samples: opening files
     char* full;
     string nam;
     nam = "inmat_out";   full = mypaste0(path, nam); ofstream inmat_out;   inmat_out.open(full);
     nam = "stats_out";   full = mypaste0(path, nam); ofstream stats_out;   stats_out.open(full);
     nam = "ll_hat_out";  full = mypaste0(path, nam); ofstream ll_hat_out;  ll_hat_out.open(full);
     nam = "ll_mean_out"; full = mypaste0(path, nam); ofstream ll_mean_out; ll_mean_out.open(full);
     nam = "lppd_out";    full = mypaste0(path, nam); ofstream lppd_out;    lppd_out.open(full);
     nam = "vlpd_out";    full = mypaste0(path, nam); ofstream vlpd_out;    vlpd_out.open(full);
     nam = "Xi_chain";    full = mypaste0(path, nam); ofstream Xi_chain;    Xi_chain.open(full);
     rowvec inmat(I*(I-1)/2, fill::zeros);
     // information criteria parameters
     double eta;
     double theta;
     double ll_mean = 0.0;
     vec lppd(Y.n_elem, fill::zeros);
     vec mell(Y.n_elem, fill::zeros);
     vec sqll(Y.n_elem, fill::zeros);
     vec That(Y.n_elem, fill::zeros);
     //chain
     double S = nburn + nskip * nsams, s; 
     double B = (double)nsams;
     uword j, i, ii, l1, l2, idx;
     for (s = 1; s <= S; ++s) {
          // update parameters
          sample_Xi    (Ij, K, nS, H, ecd, Alpha, Psi, W, La, Pi, SS, sigsq, Beta, U, Xi, Y, Simat, Samip);
          sample_Pi    (Ij, nS, H, ecd, Alpha, W, La, Xi, Pi, SS);          
          sample_W     (Ij, J, nS, H, ecd, Alpha, Psi, W, La, Xi, Pi, SS);
          sample_U     (s, nstar_N, nstar_U, delta_U, Ij, K, sigsq, Beta, U, Xi, Y);
          sample_Psi   (Ij, J, a_psi, b_psi, Psi, W);
          sample_Beta  (s, nstar_Beta, delta_Beta, Ij, omesq, Beta, U, Xi, Y);
          sample_sigsq (K, a_sig, b_sig, sigsq, U, Xi);
          // save samples
          if ((uword)s > nburn && (uword)s % nskip == 0) {
               for (i = 0; i < I; ++i) Xi_chain << Xi(i) << " ";
               Xi_chain << "\n";
               inmat += get_inmat(Xi)/B;
               vec r  = myallelic(Xi);
               stats_out << max(Xi) << " " << r(0) << " " << max(find(r)) + 1 << "\n";
               // compute information criteria inputs
               ll_mean += ll_iter_net(Ij, Beta, U, Xi, Y)/B;
               for (j = 0; j < J; ++j) {
                    l1 = accu(Ij.rows(0, j)) - Ij(j);
                    l2 = accu(Dj.rows(0, j)) - Dj(j);
                    for (i = 1; i < Ij(j); ++i) {  // lower triangualar matrix
                         for (ii = 0; ii < i; ++ii) {
                              idx   = dyx(i, ii, Ij(j)) + l2;
                              eta   = Beta(j) - norm( U.row( Xi(i + l1) - 1 ) - U.row( Xi(ii + l1) - 1 ) );
                              theta = 1.0/( 1.0 + exp(-eta) );
                              That(idx) += theta/B;
                              if (Y(idx) == 1) {
                                   lppd(idx) += theta/B;
                                   mell(idx) += log(theta)/B;
                                   sqll(idx) += pow(log(theta), 2.0);
                              } else {
                                   lppd(idx) += (1.0 - theta)/B;
                                   mell(idx) += log(1.0 - theta)/B;
                                   sqll(idx) += pow(log(1.0 - theta), 2.0);
                              }
                         }
                    }
               }
          }
          // display
          if ((uword)s % ndisp == 0) {
               vec r = myallelic(Xi);
               Rcout << " *** PROFIlE and NETWORK Model *** K = " << K << ", " << 100.0 * myround(s/S, 3) << "% done"
                     << ", mr_Beta = " << myround(nstar_Beta/s, 3)
                     << ", mr_U = "    << myround(nstar_U/nstar_N, 3)
                     << ", N = "       << max(Xi)
                     << endl;
          }
     }  // end chain
     // compute information criteria
     double ll_hat = 0.0;  // log p(y|theta hat)
     for (j = 0; j < J; ++j) {
          l2 = accu(Dj.rows(0, j)) - Dj(j);
          for (i = 1; i < Ij(j); ++i) {  // lower triangualar matrix
               for (ii = 0; ii < i; ++ii) {
                    idx = dyx(i, ii, Ij(j)) + l2;
                    if (Y(idx) == 1)
                         ll_hat += log(That(idx));
                    else
                         ll_hat += log(1.0 - That(idx));
               }
          }
     }
     // save information criteria and Theta hat
     ll_hat_out << ll_hat << "\n";
     ll_mean_out << ll_mean << "\n";
     lppd_out   << accu( log( lppd ) ) << "\n";
     vlpd_out   << accu( ( sqll - B * pow(mell, 2.0) )/(B - 1.0) ) << "\n";
     
     // save outputs
     for (i = 0; i < inmat.n_elem; ++i) inmat_out << inmat(i) << "\n";
     // close files
     stats_out.close();
     inmat_out.close();
     ll_hat_out.close();
     ll_mean_out.close();
     lppd_out.close();
     vlpd_out.close();
     Xi_chain.close();
     
}