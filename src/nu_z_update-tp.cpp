#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List nu_z_update_tp(arma::vec y,
                          int n,
                          arma::uvec index_nu,
                          arma::vec nu_z_old,
                          double sigma2_nu_z_old,
                          double sigma2_epsilon_old,
                          arma::vec mu_z_old,
                          arma::vec mu_w_old,
                          arma::vec metrop_var_nu_z,
                          arma::vec acctot_nu_z){

arma::vec nu_z = nu_z_old;
arma::vec nu_z_old_temp = vectorise(nu_z_old*trans(nu_z_old));
for(int j = 0; j < n; ++j){

   /*Second*/
   arma::vec nu_z_temp = vectorise(nu_z*trans(nu_z));
   arma::vec mu_z_temp = mu_z_old -
                         nu_z_old_temp.elem(index_nu) +
                         nu_z_temp.elem(index_nu);
   arma::vec probs_z = 1.00/(1.00 + exp(-mu_z_temp));
   arma::uvec ids1 = find(y > 0.00);
   arma::vec y_temp = log(y.elem(ids1)/(1.00 - y.elem(ids1)));
   arma::uvec ids2 = find(y == 0.00);
   double loglike = sum(log(probs_z.elem(ids1))) +
                    sum(-log(y.elem(ids1)) - log(1.00 - y.elem(ids1))) +
                    -0.50*y_temp.size()*log(2*datum::pi*sigma2_epsilon_old) - (0.50/sigma2_epsilon_old)*dot((y_temp - mu_w_old.elem(ids1)), (y_temp - mu_w_old.elem(ids1))) +
                    sum(log(1.00 - probs_z.elem(ids2)));

   double second = loglike +
                   -nu_z(j)*nu_z(j)/(2.00*sigma2_nu_z_old);

   /*First*/
   nu_z(j) = R::rnorm(nu_z_old(j), 
                      sqrt(metrop_var_nu_z(j)));

   nu_z_temp = vectorise(nu_z*trans(nu_z));
   mu_z_temp = mu_z_old -
               nu_z_old_temp.elem(index_nu) +
               nu_z_temp.elem(index_nu);
   probs_z = 1.00/(1.00 + exp(-mu_z_temp));
   ids1 = find(y > 0.00);
   y_temp = log(y.elem(ids1)/(1.00 - y.elem(ids1)));
   ids2 = find(y == 0.00);
   loglike = sum(log(probs_z.elem(ids1))) +
             sum(-log(y.elem(ids1)) - log(1.00 - y.elem(ids1))) +
             -0.50*y_temp.size()*log(2*datum::pi*sigma2_epsilon_old) - (0.50/sigma2_epsilon_old)*dot((y_temp - mu_w_old.elem(ids1)), (y_temp - mu_w_old.elem(ids1))) +
             sum(log(1.00 - probs_z.elem(ids2)));

   double first = loglike +
                  -nu_z(j)*nu_z(j)/(2.00*sigma2_nu_z_old);


   /*Decision*/
   double ratio = exp(first - second);   
   double acc = 1;
   if(ratio < R::runif(0.00, 1.00)){
  
     nu_z(j) = nu_z_old(j);
     acc = 0;
  
     }
   acctot_nu_z(j) = acctot_nu_z(j) + 
                    acc;

   }

arma::vec nu_z_temp = vectorise(nu_z*trans(nu_z));
arma::vec mu_z = mu_z_old -
                 nu_z_old_temp.elem(index_nu) +
                 nu_z_temp.elem(index_nu);

return Rcpp::List::create(Rcpp::Named("nu_z") = nu_z,
                          Rcpp::Named("acctot_nu_z") = acctot_nu_z,
                          Rcpp::Named("mu_z") = mu_z);

}



