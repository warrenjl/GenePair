#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(arma::vec y,
                             arma::mat x_pair,
                             arma::mat x_ind,
                             arma::mat z,
                             int n_star,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec beta_old,
                             arma::vec gamma_old,
                             arma::vec theta_old){

double a_sigma2_epsilon_update = 0.50*n_star + 
                                 a_sigma2_epsilon;

double b_sigma2_epsilon_update = 0.50*dot((y - x_pair*beta_old - x_ind*gamma_old - z*theta_old), (y - x_pair*beta_old - x_ind*gamma_old - z*theta_old)) + 
                                 b_sigma2_epsilon;

double sigma2_epsilon = 1.00/R::rgamma(a_sigma2_epsilon_update,
                                       (1.00/b_sigma2_epsilon_update));

return(sigma2_epsilon);

}





