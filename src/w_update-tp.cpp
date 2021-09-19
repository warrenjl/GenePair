#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec w_update_tp(arma::vec y,
                      int n_star,
                      double sigma2_epsilon_old,
                      arma::vec mu_w){

arma::vec w = mu_w +
              sqrt(sigma2_epsilon_old)*trans(arma::randn(1,
                                                         n_star)); 
  
arma::uvec ids = find(y > 0.00);
w.elem(ids) = log(y.elem(ids)/(1.00 - y.elem(ids)));

return w;

}

























































