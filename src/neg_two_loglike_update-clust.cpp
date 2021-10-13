#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update_clust(arma::vec y,
                                    arma::mat x_pair,
                                    arma::mat x_ind,
                                    arma::mat z, 
                                    arma::vec beta,
                                    arma::vec gamma,
                                    arma::vec theta){

arma::vec mu = x_pair*beta +
               x_ind*gamma +
               z*theta;
arma::vec prob = 1.00/(1.00 + exp(-mu));

double neg_two_loglike = -2*sum(y%log(prob) +
                                (1.00 - y)%log(1.00 - prob));
  
return neg_two_loglike;

}

























































