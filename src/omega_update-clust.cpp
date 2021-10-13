#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List omega_update_clust(arma::vec y,
                              arma::mat x_pair,
                              arma::mat x_ind,
                              arma::mat z,
                              int n_star,
                              int n,
                              int p_x,
                              int p_d,
                              arma::vec beta_old,
                              arma::vec gamma_old,
                              arma::vec theta_old){

arma::vec mu = x_pair*beta_old +
               x_ind*gamma_old +
               z*theta_old;

arma::vec input(1); input.fill(1.00);
arma::vec omega = rcpp_pgdraw(input, 
                              mu);
   
arma::vec lambda = (y - 0.50)/omega;

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
































































