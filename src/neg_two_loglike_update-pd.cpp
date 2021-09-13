#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update_pd(arma::vec y,
                                 arma::mat x_pair,
                                 arma::mat x_ind,
                                 arma::mat z, 
                                 int n_star,
                                 double sigma2_epsilon,
                                 arma::vec beta,
                                 arma::vec gamma,
                                 arma::vec theta){

arma::vec mu = x_pair*beta +
               x_ind*gamma +
               z*theta;

double dens = -0.50*n_star*log(2*datum::pi*sigma2_epsilon) -
              0.50*dot((y - mu), (y - mu))/sigma2_epsilon;

double neg_two_loglike = -2.00*dens;

return neg_two_loglike;

}

























































