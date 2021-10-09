#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_w_update_tp(arma::mat x,
                             arma::mat x_trans,
                             arma::mat xtx,
                             int p_x,
                             int p_d,
                             arma::mat x_prior,
                             arma::vec w,
                             double sigma2_epsilon,
                             arma::vec beta_w_old,
                             arma::vec gamma_w_g_old,
                             arma::vec gamma_w_r_old,
                             arma::vec mu_w_old){

arma::vec delta_w_old = join_cols(beta_w_old,
                                  gamma_w_g_old,
                                  gamma_w_r_old);

arma::mat cov_delta_w = inv_sympd(xtx/sigma2_epsilon + 
                                  x_prior);

arma::vec mean_delta_w = cov_delta_w*(x_trans*(w - mu_w_old + x*delta_w_old))/sigma2_epsilon;

arma::mat ind_norms = arma::randn(1, 
                                  (p_x + 2*p_d));
arma::vec delta_w = mean_delta_w + 
                    trans(ind_norms*arma::chol(cov_delta_w));

arma::vec beta_w = delta_w.subvec(0, (p_x - 1));
arma::vec gamma_w_g(p_d); gamma_w_g.fill(0.00);
arma::vec gamma_w_r(p_d); gamma_w_r.fill(0.00);
if(p_d > 0){
  
  gamma_w_g = delta_w.subvec(p_x, (p_x + p_d - 1));
  gamma_w_r = delta_w.subvec((p_x + p_d), (p_x + 2*p_d - 1));
  
  }

arma::vec mu_w = mu_w_old -
                 x*delta_w_old +
                 x*delta_w;

return Rcpp::List::create(Rcpp::Named("beta_w") = beta_w,
                          Rcpp::Named("gamma_w_g") = gamma_w_g,
                          Rcpp::Named("gamma_w_r") = gamma_w_r,
                          Rcpp::Named("mu_w") = mu_w);

}



