#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_star_update_tp(arma::vec y,
                            int n_star,
                            int n,
                            int p_x,
                            int p_d,
                            arma::vec mu_z){

arma::vec ones(n_star); ones.fill(1.00);
arma::vec w_star = rcpp_pgdraw(ones,
                               mu_z);

arma::uvec ids = find(y > 0.00);
arma::vec temp(n_star); temp.fill(0.00);
temp.elem(ids).fill(1.00);
temp = temp - 
       0.50;
  
arma::vec lambda = temp/w_star;

arma::mat w_star_mat_delta(n_star, (p_x + 2*p_d));
for(int j = 0; j < (p_x + 2*p_d); ++j){
   w_star_mat_delta.col(j) = w_star;
   }

arma::mat w_star_mat_theta(n_star, n);
for(int j = 0; j < n; ++j){
  w_star_mat_theta.col(j) = w_star;
  }

return Rcpp::List::create(Rcpp::Named("w_star") = w_star,
                          Rcpp::Named("lambda") = lambda,
                          Rcpp::Named("w_star_mat_delta") = w_star_mat_delta,
                          Rcpp::Named("w_star_mat_theta") = w_star_mat_theta);

}
































































