#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List omega_update_snp(arma::vec y,
                            arma::mat x_pair,
                            arma::mat x_ind,
                            arma::mat z,
                            int n_star,
                            int n,
                            int p_x,
                            int p_d,
                            int r,
                            arma::vec beta_old,
                            arma::vec gamma_old,
                            arma::vec theta_old){

arma::vec mu = x_pair*beta_old +
               x_ind*gamma_old +
               z*theta_old;

arma::vec omega(n_star); omega.fill(0.00);
arma::vec mu_input(1); mu_input.fill(0.00);
for(int j = 0; j < n_star; ++j){
   
   mu_input.fill(mu(j));
   omega(j) = rcpp_pgdraw((r + y(j)),
                          mu_input)(0);
   
   }
arma::vec lambda = 0.50*(y - r)/omega;

arma::mat omega_mat_delta(n_star, (p_x + p_d));
for(int j = 0; j < (p_x + p_d); ++j){
   omega_mat_delta.col(j) = omega;
   }

arma::mat omega_mat_theta(n_star, n);
for(int j = 0; j < n; ++j){
  omega_mat_theta.col(j) = omega;
  }

return Rcpp::List::create(Rcpp::Named("omega") = omega,
                          Rcpp::Named("lambda") = lambda,
                          Rcpp::Named("omega_mat_delta") = omega_mat_delta,
                          Rcpp::Named("omega_mat_theta") = omega_mat_theta);

}
































































