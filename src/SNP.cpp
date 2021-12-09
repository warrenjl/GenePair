#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List SNP(int mcmc_samples,
               arma::vec snp_distances,
               arma::mat x_pair,
               arma::mat x_ind,
               arma::mat z,
               arma::mat spatial_dists,
               arma::mat v,
               double metrop_var_phi_trans,
               Rcpp::Nullable<double> a_r_prior = R_NilValue,
               Rcpp::Nullable<double> b_r_prior = R_NilValue,
               Rcpp::Nullable<double> sigma2_regress_prior = R_NilValue,
               Rcpp::Nullable<double> a_sigma2_zeta_prior = R_NilValue,
               Rcpp::Nullable<double> b_sigma2_zeta_prior = R_NilValue,
               Rcpp::Nullable<double> a_tau2_prior = R_NilValue,
               Rcpp::Nullable<double> b_tau2_prior = R_NilValue,
               Rcpp::Nullable<double> a_phi_prior = R_NilValue,
               Rcpp::Nullable<double> b_phi_prior = R_NilValue,
               Rcpp::Nullable<double> r_init = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> gamma_init = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> theta_init = R_NilValue,
               Rcpp::Nullable<double> sigma2_zeta_init = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> eta_init = R_NilValue,
               Rcpp::Nullable<double> tau2_init = R_NilValue,
               Rcpp::Nullable<double> phi_init = R_NilValue){

//Defining Parameters and Quantities of Interest
arma::vec y = snp_distances;
int p_x = x_pair.n_cols;
int p_d = x_ind.n_cols;
int n = z.n_cols;
int m = v.n_cols;
int n_star = y.size();

arma::mat x = join_horiz(x_pair,
                         x_ind);
arma::mat x_trans = trans(x);
arma::mat xtx = x_trans*x;

arma::mat z_trans = trans(z);
arma::mat ztz = z_trans*z;

arma::mat v_trans = trans(v);
arma::mat vtv = v_trans*v;

arma::vec r(mcmc_samples); r.fill(0.00);
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
arma::mat gamma(p_d, mcmc_samples); gamma.fill(0.00);
arma::mat theta(n, mcmc_samples); theta.fill(0.00);
arma::vec sigma2_zeta(mcmc_samples); sigma2_zeta.fill(0.00);
arma::vec eta(m); eta.fill(0.00);
arma::vec tau2(mcmc_samples); tau2.fill(0.00);
arma::vec phi(mcmc_samples); phi.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
int a_r = 1;
if(a_r_prior.isNotNull()){
  a_r = Rcpp::as<int>(a_r_prior);
  }

int b_r = 50;
if(b_r_prior.isNotNull()){
  b_r = Rcpp::as<int>(b_r_prior);
  }

double sigma2_regress = 10000.00;
if(sigma2_regress_prior.isNotNull()){
  sigma2_regress = Rcpp::as<double>(sigma2_regress_prior);
  }
arma::mat x_prior = (1.00/sigma2_regress)*eye((p_x + p_d), (p_x + p_d));

double a_sigma2_zeta = 0.01;
if(a_sigma2_zeta_prior.isNotNull()){
  a_sigma2_zeta = Rcpp::as<double>(a_sigma2_zeta_prior);
  }
  
double b_sigma2_zeta = 0.01;
if(b_sigma2_zeta_prior.isNotNull()){
  b_sigma2_zeta = Rcpp::as<double>(b_sigma2_zeta_prior);
  }

double a_tau2 = 0.01;
if(a_tau2_prior.isNotNull()){
  a_tau2 = Rcpp::as<double>(a_tau2_prior);
  }

double b_tau2 = 0.01;
if(b_tau2_prior.isNotNull()){
  b_tau2 = Rcpp::as<double>(b_tau2_prior);
  }

double a_phi = 1.00;
if(a_phi_prior.isNotNull()){
  a_phi = Rcpp::as<double>(a_phi_prior);
  }
  
double b_phi = 1.00;
if(b_phi_prior.isNotNull()){
  b_phi = Rcpp::as<double>(b_phi_prior);
  }

//Initial Values
r(0) = b_r;
if(r_init.isNotNull()){
  r(0) = Rcpp::as<int>(r_init);
  }

beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

gamma.col(0).fill(0.00);
if(gamma_init.isNotNull()){
  gamma.col(0) = Rcpp::as<arma::vec>(gamma_init);
  }

theta.col(0).fill(0.00);
if(theta_init.isNotNull()){
  theta.col(0) = Rcpp::as<arma::vec>(theta_init);
  }

sigma2_zeta(0) = 1.00;
if(sigma2_zeta_init.isNotNull()){
  sigma2_zeta(0) = Rcpp::as<double>(sigma2_zeta_init);
  }

eta.fill(0.00);
if(eta_init.isNotNull()){
  eta = Rcpp::as<arma::vec>(eta_init);
  }

tau2(0) = 1.00;
if(tau2_init.isNotNull()){
  tau2(0) = Rcpp::as<double>(tau2_init);
  }

phi(0) = 1.00;
if(phi_init.isNotNull()){
  phi(0) = Rcpp::as<double>(phi_init);
  }

Rcpp::List spatial_corr_info = spatial_corr_fun(m, 
                                                spatial_dists,
                                                phi(0));
neg_two_loglike(0) = neg_two_loglike_update_snp(y,
                                                x_pair,
                                                x_ind,
                                                z, 
                                                n_star,
                                                r(0),
                                                beta.col(0),
                                                gamma.col(0),
                                                theta.col(0));

//Metropolis Settings
int acctot_phi_trans = 0;

//Main Sampling Loop
arma::vec omega(n_star); omega.fill(0.00);
arma::vec lambda(n_star); lambda.fill(0.00);
arma::mat omega_mat_delta(n_star, (p_x + p_d)); omega_mat_delta.fill(0.00);
arma::mat omega_mat_theta(n_star, n); omega_mat_theta.fill(0.00);
for(int j = 1; j < mcmc_samples; ++j){
  
   //r, omega Updates
   r(j) = r_update_snp(y,
                       x_pair,
                       x_ind,
                       z,
                       n_star,
                       a_r,
                       b_r,
                       beta.col(j-1),
                       gamma.col(j-1),
                       theta.col(j-1));
  
   Rcpp::List omega_output = omega_update_snp(y,
                                              x_pair,
                                              x_ind,
                                              z,
                                              n_star,
                                              n,
                                              p_x,
                                              p_d,
                                              r(j),
                                              beta.col(j-1),
                                              gamma.col(j-1),
                                              theta.col(j-1));
   omega = Rcpp::as<arma::vec>(omega_output[0]);
   lambda = Rcpp::as<arma::vec>(omega_output[1]);
   omega_mat_delta = Rcpp::as<arma::mat>(omega_output[2]);
   omega_mat_theta = Rcpp::as<arma::mat>(omega_output[3]);
    
   //beta, gamma Update
   Rcpp::List delta_output = delta_update_snp(x,
                                              x_trans,
                                              z,
                                              p_x,
                                              p_d,
                                              x_prior,
                                              omega,
                                              lambda,
                                              omega_mat_delta,
                                              theta.col(j-1));
  
   beta.col(j) = Rcpp::as<arma::vec>(delta_output[0]);
   gamma.col(j) = Rcpp::as<arma::vec>(delta_output[1]);
  
   //theta Update
   theta.col(j) = theta_update_snp(x_pair,
                                   x_ind,
                                   z,
                                   z_trans,
                                   v,
                                   n,
                                   omega,
                                   lambda,
                                   omega_mat_theta,
                                   beta.col(j),
                                   gamma.col(j),
                                   sigma2_zeta(j-1),
                                   eta);

   //sigma2_zeta Update
   sigma2_zeta(j) = sigma2_zeta_update(v,
                                       n,
                                       a_sigma2_zeta,
                                       b_sigma2_zeta,
                                       theta.col(j),
                                       eta);
  
   //eta Update
   eta = eta_update_pd(v_trans,
                       vtv,
                       m,
                       theta.col(j),
                       sigma2_zeta(j),
                       tau2(j-1),
                       spatial_corr_info[0]);
   
   //tau2 Update
   tau2(j) = tau2_update_pd(m,
                            a_tau2,
                            b_tau2,
                            eta,
                            spatial_corr_info[0]);
   
   //phi Update
   Rcpp::List phi_output = phi_update_pd(spatial_dists,
                                         m,
                                         a_phi,
                                         b_phi,
                                         eta,
                                         tau2(j),
                                         phi(j-1),
                                         spatial_corr_info,
                                         metrop_var_phi_trans,
                                         acctot_phi_trans);
 
   phi(j) = Rcpp::as<double>(phi_output[0]);
   acctot_phi_trans = phi_output[1];
   spatial_corr_info = phi_output[2];

   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update_snp(y,
                                                   x_pair,
                                                   x_ind,
                                                   z, 
                                                   n_star,
                                                   r(j),
                                                   beta.col(j),
                                                   gamma.col(j),
                                                   theta.col(j));
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     double accrate_phi_trans = round(100*(acctot_phi_trans/(double)j));
     Rcpp::Rcout << "phi Acceptance: " << accrate_phi_trans << "%" << std::endl;
     Rcpp::Rcout << "*******************" << std::endl;
     
     }
  
   }
                                  
return Rcpp::List::create(Rcpp::Named("r") = r,
                          Rcpp::Named("beta") = beta,
                          Rcpp::Named("gamma") = gamma,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("sigma2_zeta") = sigma2_zeta,
                          Rcpp::Named("tau2") = tau2,
                          Rcpp::Named("phi") = phi,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans);

}

