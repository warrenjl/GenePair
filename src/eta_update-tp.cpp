#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec eta_update_tp(arma::mat v_trans,
                        arma::mat vtv,
                        int m,
                        arma::vec theta,
                        double sigma2_zeta,
                        arma::vec eta_other,
                        arma::mat Sigma11,
                        arma::mat Sigma12,
                        arma::mat Sigma22_inv){

arma::mat cov_piece = inv_sympd(Sigma11 - Sigma12*(Sigma22_inv*trans(Sigma12)));
arma::mat cov_eta = inv_sympd(vtv/sigma2_zeta + 
                              cov_piece);

arma::vec mean_eta = cov_eta*((v_trans*theta)/sigma2_zeta + cov_piece*(Sigma12*(Sigma22_inv*eta_other)));

arma::mat ind_norms = arma::randn(1, 
                                  m);
arma::vec eta = mean_eta + 
                trans(ind_norms*arma::chol(cov_eta));

return(eta);

}






