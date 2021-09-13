#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double tau2_update_pd(int m,
                      double a_tau2,
                      double b_tau2,
                      arma::vec eta,
                      arma::mat corr_inv){

double a_tau2_update = 0.50*m + 
                       a_tau2;

double b_tau2_update = 0.50*dot(eta, (corr_inv*eta)) + 
                       b_tau2;

double tau2 = 1.00/R::rgamma(a_tau2_update,
                             (1.00/b_tau2_update));

return(tau2);

}





