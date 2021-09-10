#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec eta_update(arma::mat vtv,
                     arma::mat v_trans,
                     int m,
                     arma::vec theta,
                     double sigma2_zeta,
                     double tau2_old,
                     arma::mat corr_inv){

arma::mat cov_eta = inv_sympd(vtv/sigma2_zeta + 
                              corr_inv/tau2_old);

arma::vec mean_eta = cov_eta*(v_trans*theta)/sigma2_zeta;

arma::mat ind_norms = arma::randn(1, 
                                  m);
arma::vec eta = mean_eta + 
                trans(ind_norms*arma::chol(cov_eta));

return(eta);

}






