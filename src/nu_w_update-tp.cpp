#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List nu_w_update_tp(arma::vec y,
                          int n,
                          arma::uvec index_nu,
                          double sigma2_epsilon,
                          arma::vec nu_w_old,
                          double sigma2_nu_w_old,
                          arma::vec mu_z_old,
                          arma::vec mu_w_old,
                          arma::vec metrop_var_nu_w,
                          arma::vec acctot_nu_w){

arma::vec nu_w = nu_w_old;
arma::vec nu_w_old_temp = vectorise(nu_w_old*trans(nu_w_old));
for(int j = 0; j < n; ++j){

   /*Second*/
   arma::vec nu_w_temp = vectorise(nu_w*trans(nu_w));
   arma::vec mu_w_temp = mu_w_old -
                         nu_w_old_temp.elem(index_nu) +
                         nu_w_temp.elem(index_nu);
   arma::vec probs_z = 1.00/(1.00 + exp(-mu_z_old));
   arma::uvec ids1 = find(y > 0.00);
   arma::vec y_temp = log(y.elem(ids1)/(1.00 - y.elem(ids1)));
   arma::uvec ids2 = find(y == 0.00);
   double loglike = sum(log(probs_z.elem(ids1))) +
                    sum(-log(y.elem(ids1)) - log(1.00 - y.elem(ids1))) +
                    -0.50*y_temp.size()*log(2*datum::pi*sigma2_epsilon) - (0.50/sigma2_epsilon)*dot((y_temp - mu_w_temp.elem(ids1)), (y_temp - mu_w_temp.elem(ids1))) +
                    sum(log(1.00 - probs_z.elem(ids2)));

   double second = loglike +
                   -nu_w(j)*nu_w(j)/(2.00*sigma2_nu_w_old);

   /*First*/
   nu_w(j) = R::rnorm(nu_w_old(j), 
                      sqrt(metrop_var_nu_w(j)));

   nu_w_temp = vectorise(nu_w*trans(nu_w));
   mu_w_temp = mu_w_old -
               nu_w_old_temp.elem(index_nu) +
               nu_w_temp.elem(index_nu);
   probs_z = 1.00/(1.00 + exp(-mu_z_old));
   ids1 = find(y > 0.00);
   y_temp = log(y.elem(ids1)/(1.00 - y.elem(ids1)));
   ids2 = find(y == 0.00);
   loglike = sum(log(probs_z.elem(ids1))) +
             sum(-log(y.elem(ids1)) - log(1.00 - y.elem(ids1))) +
             -0.50*y_temp.size()*log(2*datum::pi*sigma2_epsilon) - (0.50/sigma2_epsilon)*dot((y_temp - mu_w_temp.elem(ids1)), (y_temp - mu_w_temp.elem(ids1))) +
             sum(log(1.00 - probs_z.elem(ids2)));

   double first = loglike +
                  -nu_w(j)*nu_w(j)/(2.00*sigma2_nu_w_old);


   /*Decision*/
   double ratio = exp(first - second);   
   double acc = 1;
   if(ratio < R::runif(0.00, 1.00)){
  
     nu_w(j) = nu_w_old(j);
     acc = 0;
  
     }
   acctot_nu_w(j) = acctot_nu_w(j) + 
                    acc;

   }

arma::vec nu_w_temp = vectorise(nu_w*trans(nu_w));
arma::vec mu_w = mu_w_old -
                 nu_w_old_temp.elem(index_nu) +
                 nu_w_temp.elem(index_nu);

return Rcpp::List::create(Rcpp::Named("nu_w") = nu_w,
                          Rcpp::Named("acctot_nu_w") = acctot_nu_w,
                          Rcpp::Named("mu_w") = mu_w);

}



