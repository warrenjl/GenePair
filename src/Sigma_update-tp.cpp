#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List Sigma_update_tp(int m,
                           arma::mat Omega_Sigma_inv,
                           double nu_Sigma_inv,
                           arma::vec eta_z_g,
                           arma::vec eta_z_r,
                           arma::vec eta_w_g,
                           arma::vec eta_w_r,
                           arma::mat spatial_corr_inv){

arma::vec eta1_temp(4); eta1_temp.fill(0.00);
arma::vec eta2_temp(4); eta2_temp.fill(0.00);
arma::mat scale(4,4); scale.fill(0.00);
for(int k = 0; k < m; ++k){
   
   eta1_temp(0) = eta_z_g(k);
   eta1_temp(1) = eta_z_r(k);
   eta1_temp(2) = eta_w_g(k);
   eta1_temp(3) = eta_w_r(k);
   
   for(int j = 0; j < m; ++j){
     
      eta2_temp(0) = eta_z_g(j);
      eta2_temp(1) = eta_z_r(j);
      eta2_temp(2) = eta_w_g(j);
      eta2_temp(3) = eta_w_r(j);
      
      scale = scale +
              (eta1_temp*trans(eta2_temp))*spatial_corr_inv(j,k);
      
      }
   
   }
scale = inv_sympd(scale + Omega_Sigma_inv);

double df = m +
            nu_Sigma_inv;

//Bartlett Decomposition
arma::mat L = arma::chol(scale);
arma::mat A(4,4); A.fill(0.00);
A(1,0) = R::rnorm(0.00,
                  sqrt(1.00));
A(2,0) = R::rnorm(0.00,
                  sqrt(1.00));
A(2,1) = R::rnorm(0.00,
                  sqrt(1.00));
A(3,0) = R::rnorm(0.00,
                  sqrt(1.00));
A(3,1) = R::rnorm(0.00,
                  sqrt(1.00));
A(3,2) = R::rnorm(0.00,
                  sqrt(1.00));

A(0,0) = sqrt(R::rchisq(df));
A(1,1) = sqrt(R::rchisq(df - 1.00));
A(2,2) = sqrt(R::rchisq(df - 2.00));
A(3,3) = sqrt(R::rchisq(df - 3.00));

arma::mat Sigma_inv = L*A*trans(A)*trans(L);
arma::mat Sigma = inv_sympd(Sigma_inv);

return Rcpp::List::create(Rcpp::Named("Sigma") = Sigma,
                          Rcpp::Named("Sigma_inv") = Sigma_inv);

}