#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec theta_update_pd(arma::vec y,
                          arma::mat x_pair,
                          arma::mat x_ind,
                          arma::mat z_trans,
                          arma::mat ztz,
                          arma::mat v,
                          int n,
                          double sigma2_epsilon,
                          arma::vec beta,
                          arma::vec gamma,
                          double sigma2_zeta_old,
                          arma::vec eta_old){

arma::mat cov_theta = inv_sympd(ztz/sigma2_epsilon + 
                                eye(n, n)/sigma2_zeta_old);

arma::vec mean_theta = cov_theta*(z_trans*(y - x_pair*beta - x_ind*gamma)/sigma2_epsilon + v*eta_old/sigma2_zeta_old);

arma::mat ind_norms = arma::randn(1, 
                                  n);
arma::vec theta = mean_theta + 
                  trans(ind_norms*arma::chol(cov_theta));

//Centering on the fly
theta = theta -
        mean(theta);

return(theta);

}






