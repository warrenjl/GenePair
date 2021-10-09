#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update_snp(arma::vec y,
                                  arma::mat x_pair,
                                  arma::mat x_ind,
                                  arma::mat z, 
                                  int n_star,
                                  int r,
                                  arma::vec beta,
                                  arma::vec gamma,
                                  arma::vec theta){

arma::vec mu = x_pair*beta +
               x_ind*gamma +
               z*theta;
arma::vec prob = 1.00/(1.00 + exp(-mu));
  
arma::vec dens(n_star); dens.fill(0.00);
for(int j = 0; j < n_star; ++j){
   dens(j) = R::dnbinom(y(j), 
                        r, 
                        (1.00 - prob(j)),        
                        TRUE);
   }

double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;

}

























































