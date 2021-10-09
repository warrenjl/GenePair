#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_z_update_tp(arma::mat x,
                             arma::mat x_trans,
                             int p_x,
                             int p_d,
                             arma::mat x_prior,
                             arma::vec w_star,
                             arma::vec lambda,
                             arma::mat w_star_mat_delta,
                             arma::vec beta_z_old,
                             arma::vec gamma_z_g_old,
                             arma::vec gamma_z_r_old,
                             arma::vec mu_z_old){

arma::vec delta_z_old = join_cols(beta_z_old,
                                  gamma_z_g_old,
                                  gamma_z_r_old);

arma::mat cov_delta_z = inv_sympd(x_trans*(w_star_mat_delta%x) + 
                                  x_prior);

arma::vec mean_delta_z = cov_delta_z*(x_trans*(w_star%(lambda - mu_z_old + x*delta_z_old)));

arma::mat ind_norms = arma::randn(1, 
                                  (p_x + 2*p_d));
arma::vec delta_z = mean_delta_z + 
                    trans(ind_norms*arma::chol(cov_delta_z));

arma::vec beta_z = delta_z.subvec(0, (p_x - 1));
arma::vec gamma_z_g(p_d); gamma_z_g.fill(0.00);
arma::vec gamma_z_r(p_d); gamma_z_r.fill(0.00);
if(p_d > 0){

  gamma_z_g = delta_z.subvec(p_x, (p_x + p_d - 1));
  gamma_z_r = delta_z.subvec((p_x + p_d), (p_x + 2*p_d - 1));
  
  }

arma::vec mu_z = mu_z_old -
                 x*delta_z_old +
                 x*delta_z;

return Rcpp::List::create(Rcpp::Named("beta_z") = beta_z,
                          Rcpp::Named("gamma_z_g") = gamma_z_g,
                          Rcpp::Named("gamma_z_r") = gamma_z_r,
                          Rcpp::Named("mu_z") = mu_z);

}



