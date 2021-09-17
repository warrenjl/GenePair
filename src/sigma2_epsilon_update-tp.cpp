#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update_tp(int n_star,
                                double a_sigma2_epsilon,
                                double b_sigma2_epsilon,
                                arma::vec w,
                                arma::vec mu_w){

double a_sigma2_epsilon_update = 0.50*n_star + 
                                 a_sigma2_epsilon;

double b_sigma2_epsilon_update = 0.50*dot((w - mu_w), (w - mu_w)) + 
                                 b_sigma2_epsilon;

double sigma2_epsilon = 1.00/R::rgamma(a_sigma2_epsilon_update,
                                       (1.00/b_sigma2_epsilon_update));

return(sigma2_epsilon);

}





