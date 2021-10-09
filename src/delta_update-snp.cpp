#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_update_snp(arma::mat x,
                            arma::mat x_trans,
                            arma::mat z,
                            int p_x,
                            int p_d,
                            arma::mat x_prior,
                            arma::vec omega,
                            arma::vec lambda,
                            arma::mat omega_mat_delta,
                            arma::vec theta_old){

arma::mat cov_delta = inv_sympd(x_trans*(omega_mat_delta%x) + 
                                x_prior);

arma::vec mean_delta = cov_delta*(x_trans*(omega%(lambda - z*theta_old)));

arma::mat ind_norms = arma::randn(1, 
                                  (p_x + p_d));
arma::vec delta = mean_delta + 
                  trans(ind_norms*arma::chol(cov_delta));

arma::vec beta = delta.subvec(0, (p_x - 1));
arma::vec gamma(p_d); gamma.fill(0.00);
if(p_d > 0){
  gamma = delta.subvec(p_x, (p_x + p_d - 1));
  }

return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("gamma") = gamma);

}



