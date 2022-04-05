#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_nu_update(int n,
                        double a_sigma2_nu,
                        double b_sigma2_nu,
                        arma::vec nu){

double a_sigma2_nu_update = 0.50*n + 
                            a_sigma2_nu;

double b_sigma2_nu_update = 0.50*dot(nu, nu) + 
                            b_sigma2_nu;

double sigma2_nu = 1.00/R::rgamma(a_sigma2_nu_update,
                                  (1.00/b_sigma2_nu_update));

return(sigma2_nu);

}





