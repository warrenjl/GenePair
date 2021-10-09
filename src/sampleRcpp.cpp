#include "RcppArmadillo.h"
#include "RcppArmadilloExtensions/sample.h"
using namespace Rcpp ;

NumericVector sampleRcpp(Rcpp::NumericVector x,
                         int size,
                         bool replace,
                         Rcpp::NumericVector prob = Rcpp::NumericVector::create()) {
Rcpp::NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
return ret;

}

