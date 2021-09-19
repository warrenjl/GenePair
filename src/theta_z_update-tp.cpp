#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List theta_z_update_tp(arma::mat z,
                             arma::mat z_trans,
                             int n,
                             arma::vec w_star,
                             arma::vec lambda,
                             arma::mat w_star_mat_theta,
                             arma::vec theta_z_old,
                             double sigma2_zeta_z_old,
                             arma::vec mu_z_old){

arma::mat cov_theta_z = inv_sympd(z_trans*(w_star_mat_theta%z) + 
                                  eye(n,n)/sigma2_zeta_z_old);
  
arma::vec mean_theta_z = cov_theta_z*(z_trans*(w_star%(lambda - mu_z_old + z*theta_z_old)));
  
arma::mat ind_norms = arma::randn(1, 
                                  n);
arma::vec theta_z = mean_theta_z + 
                    trans(ind_norms*arma::chol(cov_theta_z));
  
//Centering on the fly
theta_z = theta_z -
          mean(theta_z);

arma::vec mu_z = mu_z_old -
                 z*theta_z_old +
                 z*theta_z;

return Rcpp::List::create(Rcpp::Named("theta_z") = theta_z,
                          Rcpp::Named("mu_z") = mu_z);

}






