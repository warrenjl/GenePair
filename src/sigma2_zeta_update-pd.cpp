#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_zeta_update_pd(arma::mat v,
                             int n,
                             double a_sigma2_zeta,
                             double b_sigma2_zeta,
                             arma::vec theta,
                             arma::vec eta_old){

double a_sigma2_zeta_update = 0.50*n + 
                              a_sigma2_zeta;

double b_sigma2_zeta_update = 0.50*dot((theta - v*eta_old), (theta - v*eta_old)) + 
                              b_sigma2_zeta;

double sigma2_zeta = 1.00/R::rgamma(a_sigma2_zeta_update,
                                    (1.00/b_sigma2_zeta_update));

return(sigma2_zeta);

}





