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

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const rowvec &u_new, const uword &ip, const uvec &Ij, const double &K, const double &sigsq, const vec &Beta, mat U, uvec Xi, const uvec &Y)
{
     // computes the log conditional of Xi
     Xi(ip) = Xi_i;
     relabel(Xi);
     if ((caso == 1) && (Xi_i != curlab)) U.shed_row(curlab-1);
     if (Xi_i == newlab) U = join_vert(U, u_new);
     double out = ll_iter_net_ij(ip, Ij, Beta, U, Xi, Y);
     for (uword k = 0; k < K; ++k) out += R::dnorm(U(Xi(ip)-1, k), 0.0, sqrt(sigsq), 1);
     return( out );
}

void sample_Xi (const uvec &Ij, const double &K, const double &sigsq, const vec &Beta, mat &U, uvec &Xi, const uvec &Y, const mat &Simat)
{
     // random draw Xi, modifying U accordingly
     uword ip, i, j, curlab, newlab, prolab, caso, newval;
     double r;
     rowvec u_new(K);
     for (ip = 0; ip < accu(Ij); ++ip) {
          ijx(i, j, ip, Ij(0));
          curlab = Xi(ip);       // current label
          newlab = max(Xi) + 1;  // new label
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
          get_val_two_net(i, j, Ij, newlab, val, Xi, probs);
          // Draw proposal
          prolab = val(wsample(probs));  // proposal
          // Compute acceptance rate          
          r = exp( logcond_Xi(prolab, caso, curlab, newlab, u_new, ip, Ij, K, sigsq, Beta, U, Xi, Y) - logcond_Xi(curlab, caso, curlab, newlab, u_new, ip, Ij, K, sigsq, Beta, U, Xi, Y) );
          // Set new label
          if (R::runif(0.0, 1.0) < r)
               newval = prolab;
          else
               newval = curlab;
          // modify Xi and U* accordingly
          Xi(ip) = newval;
          relabel(Xi);
          if (caso == 1 && newval != curlab) U.shed_row(curlab-1);
          if (newval == newlab) U = join_vert(U, u_new);
     }
}

//[[Rcpp::export]]
void MCMC (uvec Xi, const uvec &Ij, const double &K, const uvec &Y, const mat &Simat, const uword &nburn, const uword &nsams, const uword &nskip, const uword &ndisp, string path)
{
     // some quantities
     double I = accu(Ij);
     double J = Ij.n_elem;
     // hyperparameters
     double a_sig = 2.0 + pow(0.5, -2.0); // CV0 = 0.5
     double b_sig = (a_sig - 1.0) * ( sqrt(I)/(sqrt(I) - 2.0) ) * ( ( pow(datum::pi, K/2.0) / exp( lgamma(K/2.0 + 1.0) ) ) * pow(I, 2.0/K) );
     double omesq = 100.0;
     // parameter initialization
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
     nam = "ll_chain"; full = mypaste0(path, nam); ofstream ll_chain; ll_chain.open(full);
     nam = "Xi_chain"; full = mypaste0(path, nam); ofstream Xi_chain; Xi_chain.open(full);
     //chain
     double S = nburn + nskip * nsams, s; 
     uword i; 
     for (s = 1; s <= S; ++s) {
          // update parameters
          sample_Xi    (Ij, K, sigsq, Beta, U, Xi, Y, Simat);
          sample_U     (s, nstar_N, nstar_U, delta_U, Ij, K, sigsq, Beta, U, Xi, Y);
          sample_Beta  (s, nstar_Beta, delta_Beta, Ij, omesq, Beta, U, Xi, Y);
          sample_sigsq (K, a_sig, b_sig, sigsq, U, Xi);
          // save samples
          if ((uword)s > nburn && (uword)s % nskip == 0) {
               for (i = 0; i < I; ++i) Xi_chain << Xi(i) << " ";
               Xi_chain << "\n";
               ll_chain << ll_iter_net(Ij, Beta, U, Xi, Y) << "\n";
          }
          // display
          if ((uword)s % ndisp == 0) {
               vec r = myallelic(Xi);
               Rcout << " * NETWORK Model * K = " << K << ", " << 100.0 * myround(s/S, 3) << "% done"
                     << ", mr_Beta = " << myround(nstar_Beta/s, 3)
                     << ", mr_U = "    << myround(nstar_U/nstar_N, 3)
                     << ", N = "       << max(Xi)
                     << ", SR = "      << myround(r(0)/I, 4) 
                     << ", Pairs = "   << r(1)
                     << endl;
          }
     }
     ll_chain.close();
     Xi_chain.close();
}