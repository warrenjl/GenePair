#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List spatial_corr_fun(int m,
                            arma::mat spatial_dists,
                            double phi){

double log_deter_inv = 0.00; 
double sign = 0.00;     

arma::mat spatial_corr = exp(-phi*spatial_dists);
arma::mat spatial_corr_inv = inv_sympd(spatial_corr);
log_det(log_deter_inv, sign, spatial_corr_inv);

return Rcpp::List::create(Rcpp::Named("spatial_corr_inv") = spatial_corr_inv,
                          Rcpp::Named("log_deter_inv") = log_deter_inv,
                          Rcpp::Named("spatial_corr") = spatial_corr);

}

