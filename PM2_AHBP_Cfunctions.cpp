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

double logp_Xi (const uvec &Ij, const double &vt, const uvec &Xi) 
{
     double I = accu(Ij), r2 = I - max(Xi), r1 = I - 2 * r2, Q2 = min(Ij), n;
     double out = 0.0;
     if (r2 == 0) {
          for (n = 2; n <= I; ++n) out += log(n);
     } else {
          if (r2 == Q2) {
               out += r2 * log(2.0);
               for (n = 2; n <= r2; ++n) out += log(n);
          } else {
               out += r2 * log(2.0); 
               for (n = 2; n <= r2; ++n) out += log(n);
               for (n = 2; n <= r1; ++n) out += log(n);
          }
     }
     out += R::dbinom(r2, Q2, vt, 1);
     return( out );
}

double logcond_Xi (const uword &Xi_i, const uword &caso, const uword &curlab, const uword &newlab, const string &pi_new, const uword &ip, const uvec &Ij, const double &vt, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const uvec &W, const uvec &La, vector<std::string> Pi, const vector<std::string> &SS, uvec Xi)
{
     // computes the log conditional of Xi
     Xi(ip) = Xi_i;
     relabel(Xi);
     if ((caso == 1) && (Xi_i != curlab)) {
          Pi.erase( Pi.begin() + curlab-1 );
     }
     if (Xi_i == newlab) {
          Pi.push_back(pi_new);
     }
     double out = ll_iter_pro_ij(ip, Ij, nS, H, ecd, Alpha, W, La, Xi, Pi, SS);
     out += log(Alpha( get_str_position(Pi[Xi(ip) - 1], SS) ));
     out += logp_Xi(Ij, vt, Xi);
     return( out );
}

void sample_Xi (const uvec &Ij, const double &vt, const double &nS, const vec &H, const vec &ecd, const vec &Alpha, const vec &Psi, uvec &W, const uvec &La, vector<std::string> &Pi, const vector<std::string> &SS, uvec &Xi, const mat &Simat)
{
     // random draw Xi, modifying U accordingly
     uword ip, i, j, curlab, newlab, prolab, caso, newval;
     double acc;
     string pi_new;
     for (ip = 0; ip < accu(Ij); ++ip) {
          ijx(i, j, ip, Ij(0));
          curlab = Xi(ip);       // current label
          newlab = max(Xi) + 1;  // new label
          pi_new = SS[wsample(Alpha)];
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
          acc = exp( logcond_Xi(prolab, caso, curlab, newlab, pi_new, ip, Ij, vt, nS, H, ecd, Alpha, W, La, Pi, SS, Xi) - logcond_Xi(curlab, caso, curlab, newlab, pi_new, ip, Ij, vt, nS, H, ecd, Alpha, W, La, Pi, SS, Xi) );
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
          }
          if (newval == newlab){ 
               Pi.push_back(pi_new);
               sample_W_i(ip, j, Ij, nS, H, ecd, Alpha, Psi, W, La, Xi, Pi, SS);
          }
     }
}

void sample_vt (const double &I, const double &a_vt, const double &b_vt, double &vt, const uvec &Xi)
{
     double r2 = I - max(Xi);
     vt = R::rbeta(a_vt + r2, b_vt + floor(I/2.0) - r2);
}

//[[Rcpp::export]]
void MCMC (const double &a_psi, const double &b_psi, uvec Xi, const double &a_vt, const double &b_vt, const uvec &Ij, const vec &H, const vec &ecd, const vec &Alpha, const uvec &La, const vector<std::string> &SS, const mat &Simat, const uword &nburn, const uword &nsams, const uword &nskip, const uword &ndisp, string path)
{
     // some quantities
     double I  = accu(Ij);    
     double J  = Ij.n_elem;
     double nS = SS.size();
     // parameter initialization
     double vt = R::rbeta(a_vt, b_vt);
     vec  Psi(J);
     uvec W(I);
     vector<std::string> Pi(max(Xi));
     prior_int_pro (Ij, J, a_psi, b_psi, nS, H, ecd, Alpha, Psi, W, La, Xi, Pi, SS);
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
          sample_Xi    (Ij, vt, nS, H, ecd, Alpha, Psi, W, La, Pi, SS, Xi, Simat);
          sample_Pi    (Ij, nS, H, ecd, Alpha, W, La, Xi, Pi, SS);          
          sample_W     (Ij, J, nS, H, ecd, Alpha, Psi, W, La, Xi, Pi, SS);
          sample_vt    (I, a_vt, b_vt, vt, Xi);
          sample_Psi   (Ij, J, a_psi, b_psi, Psi, W);
          // save samples
          if ((uword)s > nburn && (uword)s % nskip == 0) {
               for (i = 0; i < I; ++i) Xi_chain << Xi(i) << " ";
               Xi_chain << "\n";
               ll_chain << ll_iter_pro(Ij, nS, H, ecd, Alpha, W, La, Xi, Pi, SS) << "\n";
          }
          // display
          if ((uword)s % ndisp == 0) {
               vec r = myallelic(Xi);
               Rcout << " * PROFIlE Model * " << 100.0 * myround(s/S, 3) << "% complete"
                     << ", N = "       << max(Xi)
                     << ", SR = "      << myround(r(0)/I, 4) 
                     << ", Pairs = "    << r(1)
                     << endl;
          }
     }
     ll_chain.close();
     Xi_chain.close();
}