#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_update_pd(arma::vec y,
                           arma::mat x_trans,
                           arma::mat xtx,
                           arma::mat z,
                           int p_x,
                           int p_d,
                           arma::mat x_prior,
                           double sigma2_epsilon,
                           arma::vec theta_old){

arma::mat cov_delta = inv_sympd(xtx/sigma2_epsilon + 
                                x_prior);

arma::vec mean_delta = cov_delta*(x_trans*(y - z*theta_old))/sigma2_epsilon;

arma::mat ind_norms = arma::randn(1, 
                                  (p_x + p_d));
arma::vec delta = mean_delta + 
                  trans(ind_norms*arma::chol(cov_delta));

arma::vec beta = delta.subvec(0, (p_x - 1));
arma::vec gamma = delta.subvec(p_x, (p_x + p_d - 1));

return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("gamma") = gamma);

}



