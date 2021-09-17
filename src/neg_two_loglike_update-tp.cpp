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
  
arma::vec y_temp = log(y(y > 0.00)/(1.00 - y(y > 0.00)));

double neg_two_loglike = -2.00*sum(log(probs_z(y > 0.00)) -
                                   0.50*y_temp.size()*log(2*datum::pi*sigma2_epsilon) -
                                   (0.50/sigma2_epsilon)*dot((y_temp - mu_w(y > 0.00)), (y_temp - mu_w(y > 0.00)))) -
              2*sum(log(1.00 - probs_z(y == 0)));
           
return neg_two_loglike;

}

























































