#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

int r_update_snp(arma::vec y,
                 arma::mat x_pair,
                 arma::mat x_ind,
                 arma::mat z,
                 int n_star,
                 double a_r,
                 double b_r,
                 arma::vec beta_old,
                 arma::vec gamma_old,
                 arma::vec theta_old){

arma::vec mu = x_pair*beta_old +
               x_ind*gamma_old +
               z*theta_old;
  
arma::vec prob = 1.00/(1.00 + exp(-mu));
  
arma::vec r_log_val(b_r - a_r + 1); r_log_val.fill(0.00);  
int counter = 0;
for(int j = (a_r - 1); j < b_r; ++j){

   for(int k = 0; k < n_star; ++k){
      r_log_val(counter) = r_log_val(counter) +
                           R::dnbinom(y(k),
                                      (j + 1),
                                      (1.00 - prob(k)),
                                      TRUE);
      }
   counter = counter +
             1;
  
   }
  
arma::vec r_prob(b_r - a_r + 1); r_prob.fill(0.00);
for(int j = 0; j < (b_r - a_r + 1); ++j){
   r_prob(j) = 1.00/sum(exp(r_log_val - r_log_val(j)));
   }

IntegerVector sample_set = seq(a_r, b_r);
int r = sampleRcpp(wrap(sample_set), 
                   1, 
                   TRUE, 
                   wrap(r_prob))(0);
    
return(r);

}





