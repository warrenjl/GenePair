#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update_tp(arma::vec y,
                                 double sigma2_epsilon,
                                 arma::vec mu_z,
                                 arma::vec mu_w){
  
arma::vec probs_z = 1.00/(1.00 + exp(-mu_z));
  
arma::uvec ids1 = find(y > 0.00);
arma::vec y_temp = log(y.elem(ids1)/(1.00 - y.elem(ids1)));
arma::uvec ids2 = find(y == 0.00);

double loglike = sum(log(probs_z.elem(ids1))) +
                 sum(-log(y.elem(ids1)) - log(1.00 - y.elem(ids1))) +
                 -0.50*y_temp.size()*log(2*datum::pi*sigma2_epsilon) - (0.50/sigma2_epsilon)*dot((y_temp - mu_w.elem(ids1)), (y_temp - mu_w.elem(ids1))) +
                 sum(log(1.00 - probs_z.elem(ids2)));

double neg_two_loglike = -2.00*loglike;

return neg_two_loglike;

}

























































