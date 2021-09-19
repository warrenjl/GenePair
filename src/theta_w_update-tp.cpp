#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List theta_w_update_tp(arma::mat z,
                             arma::mat z_trans,
                             arma::mat zgtzg,
                             arma::mat v,
                             int n,
                             arma::vec w,
                             double sigma2_epsilon,
                             arma::vec theta_w_old,
                             double sigma2_zeta_w_old,
                             arma::vec eta_w_old,
                             arma::vec mu_w_old){

arma::mat cov_theta_w = inv_sympd(zgtzg/sigma2_epsilon + 
                                  eye(n,n)/sigma2_zeta_w_old);
  
arma::vec mean_theta_w = cov_theta_w*(z_trans*(w - mu_w_old + z*theta_w_old)/sigma2_epsilon + v*eta_w_old/sigma2_zeta_w_old);
  
arma::mat ind_norms = arma::randn(1, 
                                  n);
arma::vec theta_w = mean_theta_w + 
                    trans(ind_norms*arma::chol(cov_theta_w));
  
//Centering on the fly
theta_w = theta_w -
          mean(theta_w);

arma::vec mu_w = mu_w_old -
                 z*theta_w_old +
                 z*theta_w;

return Rcpp::List::create(Rcpp::Named("theta_w") = theta_w,
                          Rcpp::Named("mu_w") = mu_w);

}






