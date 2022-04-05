#include "RcppArmadillo.h"
#include "GenePair.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List Trans_Prob(int mcmc_samples,
                      arma::vec transmission_probabilities,
                      arma::mat x_pair,
                      arma::mat x_ind_g,
                      arma::mat x_ind_r,
                      arma::mat z_g,
                      arma::mat z_r,
                      arma::mat spatial_dists,
                      arma::mat v,
                      double metrop_var_phi_trans,
                      arma::vec metrop_var_nu_z,
                      arma::vec metrop_var_nu_w,
                      Rcpp::Nullable<double> sigma2_regress_prior = R_NilValue,  //Start of Priors
                      Rcpp::Nullable<double> a_sigma2_nu_z_prior = R_NilValue,
                      Rcpp::Nullable<double> b_sigma2_nu_z_prior = R_NilValue,
                      Rcpp::Nullable<double> a_sigma2_zeta_z_g_prior = R_NilValue,
                      Rcpp::Nullable<double> b_sigma2_zeta_z_g_prior = R_NilValue,
                      Rcpp::Nullable<double> a_sigma2_zeta_z_r_prior = R_NilValue,
                      Rcpp::Nullable<double> b_sigma2_zeta_z_r_prior = R_NilValue,
                      Rcpp::Nullable<double> a_sigma2_epsilon_prior = R_NilValue,
                      Rcpp::Nullable<double> b_sigma2_epsilon_prior = R_NilValue,
                      Rcpp::Nullable<double> a_sigma2_nu_w_prior = R_NilValue,
                      Rcpp::Nullable<double> b_sigma2_nu_w_prior = R_NilValue,
                      Rcpp::Nullable<double> a_sigma2_zeta_w_g_prior = R_NilValue,
                      Rcpp::Nullable<double> b_sigma2_zeta_w_g_prior = R_NilValue,
                      Rcpp::Nullable<double> a_sigma2_zeta_w_r_prior = R_NilValue,
                      Rcpp::Nullable<double> b_sigma2_zeta_w_r_prior = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericMatrix> Omega_Sigma_inv_prior = R_NilValue,
                      Rcpp::Nullable<double> nu_Sigma_inv_prior = R_NilValue,
                      Rcpp::Nullable<double> a_phi_prior = R_NilValue,
                      Rcpp::Nullable<double> b_phi_prior = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> beta_z_init = R_NilValue,  //Start of Initial Values
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_z_g_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_z_r_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> nu_z_init = R_NilValue,
                      Rcpp::Nullable<double> sigma2_nu_z_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_z_g_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_z_r_init = R_NilValue,
                      Rcpp::Nullable<double> sigma2_zeta_z_g_init = R_NilValue,
                      Rcpp::Nullable<double> sigma2_zeta_z_r_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_z_g_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_z_r_init = R_NilValue,
                      Rcpp::Nullable<double> sigma2_epsilon_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> beta_w_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_w_g_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_w_r_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> nu_w_init = R_NilValue,
                      Rcpp::Nullable<double> sigma2_nu_w_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_w_g_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_w_r_init = R_NilValue,
                      Rcpp::Nullable<double> sigma2_zeta_w_g_init = R_NilValue,
                      Rcpp::Nullable<double> sigma2_zeta_w_r_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_w_g_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_w_r_init = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericMatrix> Sigma_init = R_NilValue,
                      Rcpp::Nullable<double> phi_init = R_NilValue){

//Defining Parameters and Quantities of Interest
arma::vec y = transmission_probabilities;
int p_x = x_pair.n_cols;
int p_d = x_ind_g.n_cols;  //x_g and x_r have the Same # of Columns
int n = z_g.n_cols;  //z_g and z_r have the Same # of Columns
int m = v.n_cols;
int n_star = y.size();

arma::uvec index_nu(n*(n - 1)); index_nu.fill(0);
int counter1 = 0;
int counter2 = 0;
for(int j = 0; j < n; ++j){
   for(int k = 0; k < n; ++k){
      
      if(j != k){
         
        index_nu(counter2) = counter1;  //Needed for \nu Interactions
        ++counter2;
           
        }
      
      ++counter1;
         
      }
   }

arma::mat x = join_horiz(x_pair,
                         x_ind_g,
                         x_ind_r);
arma::mat x_trans = trans(x);
arma::mat xtx = x_trans*x;

arma::mat zg_trans = trans(z_g);
arma::mat zgtzg = zg_trans*z_g;

arma::mat zr_trans = trans(z_r);
arma::mat zrtzr = zr_trans*z_r;

arma::mat v_trans = trans(v);
arma::mat vtv = v_trans*v;

arma::mat beta_z(p_x, mcmc_samples); beta_z.fill(0.00);
arma::mat gamma_z_g(p_d, mcmc_samples); gamma_z_g.fill(0.00);
arma::mat gamma_z_r(p_d, mcmc_samples); gamma_z_r.fill(0.00);
arma::mat nu_z(n, mcmc_samples); nu_z.fill(0.00);
arma::vec sigma2_nu_z(mcmc_samples); sigma2_nu_z.fill(0.00);
arma::mat theta_z_g(n, mcmc_samples); theta_z_g.fill(0.00);
arma::mat theta_z_r(n, mcmc_samples); theta_z_r.fill(0.00);
arma::vec sigma2_zeta_z_g(mcmc_samples); sigma2_zeta_z_g.fill(0.00);
arma::vec sigma2_zeta_z_r(mcmc_samples); sigma2_zeta_z_r.fill(0.00);
arma::vec eta_z_g(m); eta_z_g.fill(0.00);
arma::vec eta_z_r(m); eta_z_r.fill(0.00);
arma::vec sigma2_epsilon(mcmc_samples); sigma2_epsilon.fill(0.00);
arma::mat beta_w(p_x, mcmc_samples); beta_w.fill(0.00);
arma::mat gamma_w_g(p_d, mcmc_samples); gamma_w_g.fill(0.00);
arma::mat gamma_w_r(p_d, mcmc_samples); gamma_w_r.fill(0.00);
arma::mat nu_w(n, mcmc_samples); nu_w.fill(0.00);
arma::vec sigma2_nu_w(mcmc_samples); sigma2_nu_w.fill(0.00);
arma::mat theta_w_g(n, mcmc_samples); theta_w_g.fill(0.00);
arma::mat theta_w_r(n, mcmc_samples); theta_w_r.fill(0.00);
arma::vec sigma2_zeta_w_g(mcmc_samples); sigma2_zeta_w_g.fill(0.00);
arma::vec sigma2_zeta_w_r(mcmc_samples); sigma2_zeta_w_r.fill(0.00);
arma::vec eta_w_g(m); eta_w_g.fill(0.00);
arma::vec eta_w_r(m); eta_w_r.fill(0.00);
Rcpp::List Sigma(mcmc_samples); Sigma.fill(0.00);
arma::vec phi(mcmc_samples); phi.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_regress = 10000.00;
if(sigma2_regress_prior.isNotNull()){
  sigma2_regress = Rcpp::as<double>(sigma2_regress_prior);
  }
arma::mat x_prior = (1.00/sigma2_regress)*eye((p_x + 2*p_d), (p_x + 2*p_d));

double a_sigma2_nu_z = 0.01;
if(a_sigma2_nu_z_prior.isNotNull()){
  a_sigma2_nu_z = Rcpp::as<double>(a_sigma2_nu_z_prior);
  }

double b_sigma2_nu_z = 0.01;
if(b_sigma2_nu_z_prior.isNotNull()){
  b_sigma2_nu_z = Rcpp::as<double>(b_sigma2_nu_z_prior);
  }

double a_sigma2_zeta_z_g = 0.01;
if(a_sigma2_zeta_z_g_prior.isNotNull()){
  a_sigma2_zeta_z_g = Rcpp::as<double>(a_sigma2_zeta_z_g_prior);
  }

double b_sigma2_zeta_z_g = 0.01;
if(b_sigma2_zeta_z_g_prior.isNotNull()){
  b_sigma2_zeta_z_g = Rcpp::as<double>(b_sigma2_zeta_z_g_prior);
  }

double a_sigma2_zeta_z_r = 0.01;
if(a_sigma2_zeta_z_r_prior.isNotNull()){
  a_sigma2_zeta_z_r = Rcpp::as<double>(a_sigma2_zeta_z_r_prior);
  }

double b_sigma2_zeta_z_r = 0.01;
if(b_sigma2_zeta_z_r_prior.isNotNull()){
  b_sigma2_zeta_z_r = Rcpp::as<double>(b_sigma2_zeta_z_r_prior);
  }

double a_sigma2_epsilon = 0.01;
if(a_sigma2_epsilon_prior.isNotNull()){
  a_sigma2_epsilon = Rcpp::as<double>(a_sigma2_epsilon_prior);
  }

double b_sigma2_epsilon = 0.01;
if(b_sigma2_epsilon_prior.isNotNull()){
  b_sigma2_epsilon = Rcpp::as<double>(b_sigma2_epsilon_prior);
  }

double a_sigma2_nu_w = 0.01;
if(a_sigma2_nu_w_prior.isNotNull()){
  a_sigma2_nu_w = Rcpp::as<double>(a_sigma2_nu_w_prior);
  }

double b_sigma2_nu_w = 0.01;
if(b_sigma2_nu_w_prior.isNotNull()){
  b_sigma2_nu_w = Rcpp::as<double>(b_sigma2_nu_w_prior);
  }

double a_sigma2_zeta_w_g = 0.01;
if(a_sigma2_zeta_w_g_prior.isNotNull()){
  a_sigma2_zeta_w_g = Rcpp::as<double>(a_sigma2_zeta_w_g_prior);
  }

double b_sigma2_zeta_w_g = 0.01;
if(b_sigma2_zeta_w_g_prior.isNotNull()){
  b_sigma2_zeta_w_g = Rcpp::as<double>(b_sigma2_zeta_w_g_prior);
  }

double a_sigma2_zeta_w_r = 0.01;
if(a_sigma2_zeta_w_r_prior.isNotNull()){
  a_sigma2_zeta_w_r = Rcpp::as<double>(a_sigma2_zeta_w_r_prior);
  }

double b_sigma2_zeta_w_r = 0.01;
if(b_sigma2_zeta_w_r_prior.isNotNull()){
  b_sigma2_zeta_w_r = Rcpp::as<double>(b_sigma2_zeta_w_r_prior);
  }

arma::mat Omega_Sigma_inv(4,4); Omega_Sigma_inv.eye();
if(Omega_Sigma_inv_prior.isNotNull()){
  Omega_Sigma_inv = Rcpp::as<arma::mat>(Omega_Sigma_inv_prior);
  }

double nu_Sigma_inv = 5.00;
if(nu_Sigma_inv_prior.isNotNull()){
  nu_Sigma_inv = Rcpp::as<double>(nu_Sigma_inv_prior);
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
beta_z.col(0).fill(0.00);
if(beta_z_init.isNotNull()){
  beta_z.col(0) = Rcpp::as<arma::vec>(beta_z_init);
  }

gamma_z_g.col(0).fill(0.00);
if(gamma_z_g_init.isNotNull()){
  gamma_z_g.col(0) = Rcpp::as<arma::vec>(gamma_z_g_init);
  }

gamma_z_r.col(0).fill(0.00);
if(gamma_z_r_init.isNotNull()){
  gamma_z_r.col(0) = Rcpp::as<arma::vec>(gamma_z_r_init);
  }

nu_z.col(0).fill(0.00);
if(nu_z_init.isNotNull()){
  nu_z.col(0) = Rcpp::as<arma::vec>(nu_z_init);
  }

sigma2_nu_z(0) = 1.00;
if(sigma2_nu_z_init.isNotNull()){
  sigma2_nu_z(0) = Rcpp::as<double>(sigma2_nu_z_init);
  }

theta_z_g.col(0).fill(0.00);
if(theta_z_g_init.isNotNull()){
  theta_z_g.col(0) = Rcpp::as<arma::vec>(theta_z_g_init);
  }

theta_z_r.col(0).fill(0.00);
if(theta_z_r_init.isNotNull()){
  theta_z_r.col(0) = Rcpp::as<arma::vec>(theta_z_r_init);
  }

sigma2_zeta_z_g(0) = 1.00;
if(sigma2_zeta_z_g_init.isNotNull()){
  sigma2_zeta_z_g(0) = Rcpp::as<double>(sigma2_zeta_z_g_init);
  }

sigma2_zeta_z_r(0) = 1.00;
if(sigma2_zeta_z_r_init.isNotNull()){
  sigma2_zeta_z_r(0) = Rcpp::as<double>(sigma2_zeta_z_r_init);
  }

eta_z_g.fill(0.00);
if(eta_z_g_init.isNotNull()){
  eta_z_g = Rcpp::as<arma::vec>(eta_z_g_init);
  }

eta_z_r.fill(0.00);
if(eta_z_r_init.isNotNull()){
  eta_z_r = Rcpp::as<arma::vec>(eta_z_r_init);
  }

sigma2_epsilon(0) = var(y);
if(sigma2_epsilon_init.isNotNull()){
  sigma2_epsilon(0) = Rcpp::as<double>(sigma2_epsilon_init);
  }

beta_w.col(0).fill(0.00);
if(beta_w_init.isNotNull()){
  beta_w.col(0) = Rcpp::as<arma::vec>(beta_w_init);
  }

gamma_w_g.col(0).fill(0.00);
if(gamma_w_g_init.isNotNull()){
  gamma_w_g.col(0) = Rcpp::as<arma::vec>(gamma_w_g_init);
  }

gamma_w_r.col(0).fill(0.00);
if(gamma_w_r_init.isNotNull()){
  gamma_w_r.col(0) = Rcpp::as<arma::vec>(gamma_w_r_init);
  }

nu_w.col(0).fill(0.00);
if(nu_w_init.isNotNull()){
  nu_w.col(0) = Rcpp::as<arma::vec>(nu_w_init);
  }

sigma2_nu_w(0) = 1.00;
if(sigma2_nu_w_init.isNotNull()){
  sigma2_nu_w(0) = Rcpp::as<double>(sigma2_nu_w_init);
  }

theta_w_g.col(0).fill(0.00);
if(theta_w_g_init.isNotNull()){
  theta_w_g.col(0) = Rcpp::as<arma::vec>(theta_w_g_init);
  }

theta_w_r.col(0).fill(0.00);
if(theta_w_r_init.isNotNull()){
  theta_w_r.col(0) = Rcpp::as<arma::vec>(theta_w_r_init);
  }

sigma2_zeta_w_g(0) = 1.00;
if(sigma2_zeta_w_g_init.isNotNull()){
  sigma2_zeta_w_g(0) = Rcpp::as<double>(sigma2_zeta_w_g_init);
  }

sigma2_zeta_w_r(0) = 1.00;
if(sigma2_zeta_w_r_init.isNotNull()){
  sigma2_zeta_w_r(0) = Rcpp::as<double>(sigma2_zeta_w_r_init);
  }

eta_w_g.fill(0.00);
if(eta_w_g_init.isNotNull()){
  eta_w_g = Rcpp::as<arma::vec>(eta_w_g_init);
  }

eta_w_r.fill(0.00);
if(eta_w_r_init.isNotNull()){
  eta_w_r = Rcpp::as<arma::vec>(eta_w_r_init);
  }

arma::mat Sigma_temp(4,4); Sigma_temp.eye();
Sigma[0] = Sigma_temp;
if(Sigma_init.isNotNull()){
  Sigma[0] = Rcpp::as<arma::mat>(Sigma_init);
  }

phi(0) = 1.00;
if(phi_init.isNotNull()){
  phi(0) = Rcpp::as<double>(phi_init);
  }

Rcpp::List spatial_corr_info = spatial_corr_fun(m, 
                                                spatial_dists,
                                                phi(0));

arma::vec nu_z_temp = vectorise(nu_z.col(0)*trans(nu_z.col(0)));
arma::vec mu_z = x_pair*beta_z.col(0) + 
                 x_ind_g*gamma_z_g.col(0) +
                 x_ind_r*gamma_z_r.col(0) +
                 z_g*theta_z_g.col(0) + 
                 z_r*theta_z_r.col(0) +
                 nu_z_temp.elem(index_nu);

arma::vec nu_w_temp = vectorise(nu_w.col(0)*trans(nu_w.col(0)));
arma::vec mu_w = x_pair*beta_w.col(0) + 
                 x_ind_g*gamma_w_g.col(0) +
                 x_ind_r*gamma_w_r.col(0) +
                 z_g*theta_w_g.col(0) + 
                 z_r*theta_w_r.col(0) +
                 nu_w_temp.elem(index_nu);

neg_two_loglike(0) = neg_two_loglike_update_tp(y,
                                               sigma2_epsilon(0),
                                               mu_z,
                                               mu_w);

//Metropolis Settings
int acctot_phi_trans = 0;
arma::vec acctot_nu_z(n); acctot_nu_z.fill(0);
arma::vec acctot_nu_w(n); acctot_nu_w.fill(0);

//Main Sampling Loop
arma::vec w_star(n_star); w_star.fill(0.00);
arma::vec lambda(n_star); lambda.fill(0.00);
arma::mat w_star_mat_delta(n_star, (p_x + 2*p_d)); w_star_mat_delta.fill(0.00);
arma::mat w_star_mat_theta(n_star, n); w_star_mat_theta.fill(0.00);
arma::vec w(n_star); w.fill(0.00);
arma::mat Sigma11(m,m); Sigma11.fill(0.00);
arma::mat Sigma12(m, (3*m)); Sigma12.fill(0.00);
arma::mat Sigma22_inv((3*m), (3*m)); Sigma22_inv.fill(0.00);
arma::uvec index(3); index.fill(0);
arma::vec eta_other(3*m); eta_other.fill(0.00);
arma::mat Sigma_inv(4,4); Sigma_inv.fill(0.00);
arma::vec eta_full(4*m); eta_full.fill(0.00);
for(int j = 1; j < mcmc_samples; ++j){
  
   //w_star Update
   Rcpp::List w_star_output = w_star_update_tp(y,
                                               n_star,
                                               n,
                                               p_x,
                                               p_d,
                                               mu_z);
   w_star = Rcpp::as<arma::vec>(w_star_output[0]);
   lambda = Rcpp::as<arma::vec>(w_star_output[1]);
   w_star_mat_delta = Rcpp::as<arma::mat>(w_star_output[2]);
   w_star_mat_theta = Rcpp::as<arma::mat>(w_star_output[3]);
   
   //beta_z, gamma_z_g, gamma_z_r Update
   Rcpp::List delta_z_output = delta_z_update_tp(x,
                                                 x_trans,
                                                 p_x,
                                                 p_d,
                                                 x_prior,
                                                 w_star,
                                                 lambda,
                                                 w_star_mat_delta,
                                                 beta_z.col(j-1),
                                                 gamma_z_g.col(j-1),
                                                 gamma_z_r.col(j-1),
                                                 mu_z);
   
   beta_z.col(j) = Rcpp::as<arma::vec>(delta_z_output[0]);
   gamma_z_g.col(j) = Rcpp::as<arma::vec>(delta_z_output[1]);
   gamma_z_r.col(j) = Rcpp::as<arma::vec>(delta_z_output[2]);
   mu_z = Rcpp::as<arma::vec>(delta_z_output[3]);
   
   //nu_z Update
   Rcpp::List nu_z_output = nu_z_update_tp(y,
                                           n,
                                           index_nu,
                                           nu_z.col(j-1),
                                           sigma2_nu_z(j-1),
                                           sigma2_epsilon(j-1),
                                           mu_z,
                                           mu_w,
                                           metrop_var_nu_z,
                                           acctot_nu_z);
   
   nu_z.col(j) = Rcpp::as<arma::vec>(nu_z_output[0]);
   acctot_nu_z = Rcpp::as<arma::vec>(nu_z_output[1]);
   mu_z = Rcpp::as<arma::vec>(nu_z_output[2]);
      
   //sigma2_nu_z Update
   sigma2_nu_z(j) = sigma2_nu_update(n,
                                     a_sigma2_nu_z,
                                     b_sigma2_nu_z,
                                     nu_z.col(j));
    
   //theta_z_g Update
   Rcpp::List theta_z_g_output = theta_z_update_tp(z_g,
                                                   zg_trans,
                                                   v,
                                                   n,
                                                   w_star,
                                                   lambda,
                                                   w_star_mat_theta,
                                                   theta_z_g.col(j-1),
                                                   sigma2_zeta_z_g(j-1),
                                                   eta_z_g,
                                                   mu_z);
   
   theta_z_g.col(j) = Rcpp::as<arma::vec>(theta_z_g_output[0]);
   mu_z = Rcpp::as<arma::vec>(theta_z_g_output[1]);
   
   //theta_z_r Update
   Rcpp::List theta_z_r_output = theta_z_update_tp(z_r,
                                                   zr_trans,
                                                   v,
                                                   n,
                                                   w_star,
                                                   lambda,
                                                   w_star_mat_theta,
                                                   theta_z_r.col(j-1),
                                                   sigma2_zeta_z_r(j-1),
                                                   eta_z_r,
                                                   mu_z);
   
   theta_z_r.col(j) = Rcpp::as<arma::vec>(theta_z_r_output[0]);
   mu_z = Rcpp::as<arma::vec>(theta_z_r_output[1]);
   
   //sigma2_zeta_z_g Update
   sigma2_zeta_z_g(j) = sigma2_zeta_update(v,
                                           n,
                                           a_sigma2_zeta_z_g,
                                           b_sigma2_zeta_z_g,
                                           theta_z_g.col(j),
                                           eta_z_g);
  
   //sigma2_zeta_z_r Update
   sigma2_zeta_z_r(j) = sigma2_zeta_update(v,
                                           n,
                                           a_sigma2_zeta_z_r,
                                           b_sigma2_zeta_z_r,
                                           theta_z_r.col(j),
                                           eta_z_r);
  
   //eta_z_g Update
   Sigma11 = Rcpp::as<arma::mat>(Sigma[j-1])(0,0)*Rcpp::as<arma::mat>(spatial_corr_info[2]);
   Sigma12 = join_horiz(Rcpp::as<arma::mat>(Sigma[j-1])(0,1)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(0,2)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(0,3)*Rcpp::as<arma::mat>(spatial_corr_info[2]));
  
   index(0) = 1; index(1) = 2; index(2) = 3;
   Sigma22_inv = kron(inv_sympd(Rcpp::as<arma::mat>(Sigma[j-1]).submat(index, index)), 
                      Rcpp::as<arma::mat>(spatial_corr_info[0]));
    
   eta_other = join_cols(eta_z_r,
                         eta_w_g,
                         eta_w_r);
    
   eta_z_g = eta_update_tp(v_trans,
                           vtv,
                           m,
                           theta_z_g.col(j),
                           sigma2_zeta_z_g(j),
                           eta_other,
                           Sigma11,
                           Sigma12,
                           Sigma22_inv);
  
   //eta_z_r Update
   Sigma11 = Rcpp::as<arma::mat>(Sigma[j-1])(1,1)*Rcpp::as<arma::mat>(spatial_corr_info[2]);
   Sigma12 = join_horiz(Rcpp::as<arma::mat>(Sigma[j-1])(1,0)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(1,2)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(1,3)*Rcpp::as<arma::mat>(spatial_corr_info[2]));
  
   index(0) = 0; index(1) = 2; index(2) = 3;
   Sigma22_inv = kron(inv_sympd(Rcpp::as<arma::mat>(Sigma[j-1]).submat(index, index)), 
                      Rcpp::as<arma::mat>(spatial_corr_info[0]));
  
   eta_other = join_cols(eta_z_g,
                         eta_w_g,
                         eta_w_r);
  
   eta_z_r = eta_update_tp(v_trans,
                           vtv,
                           m,
                           theta_z_r.col(j),
                           sigma2_zeta_z_r(j),
                           eta_other,
                           Sigma11,
                           Sigma12,
                           Sigma22_inv);
  
   //w Update
   w = w_update_tp(y,
                   n_star,
                   sigma2_epsilon(j-1),
                   mu_w);
  
   //sigma2_epsilon Update
   sigma2_epsilon(j) = sigma2_epsilon_update_tp(n_star,
                                                a_sigma2_epsilon,
                                                b_sigma2_epsilon,
                                                w,
                                                mu_w);
  
   //beta_w, gamma_w_g, gamma_w_r Update
   Rcpp::List delta_w_output = delta_w_update_tp(x,
                                                 x_trans,
                                                 xtx,
                                                 p_x,
                                                 p_d,
                                                 x_prior,
                                                 w,
                                                 sigma2_epsilon(j),
                                                 beta_w.col(j-1),
                                                 gamma_w_g.col(j-1),
                                                 gamma_w_r.col(j-1),
                                                 mu_w);
  
   beta_w.col(j) = Rcpp::as<arma::vec>(delta_w_output[0]);
   gamma_w_g.col(j) = Rcpp::as<arma::vec>(delta_w_output[1]);
   gamma_w_r.col(j) = Rcpp::as<arma::vec>(delta_w_output[2]);
   mu_w = Rcpp::as<arma::vec>(delta_w_output[3]);
   
   //nu_w Update
   Rcpp::List nu_w_output = nu_w_update_tp(y,
                                           n,
                                           index_nu,
                                           sigma2_epsilon(j),
                                           nu_w.col(j-1),
                                           sigma2_nu_w(j-1),
                                           mu_z,
                                           mu_w,
                                           metrop_var_nu_w,
                                           acctot_nu_w);
   
   nu_w.col(j) = Rcpp::as<arma::vec>(nu_w_output[0]);
   acctot_nu_w = Rcpp::as<arma::vec>(nu_w_output[1]);
   mu_w = Rcpp::as<arma::vec>(nu_w_output[2]);
   
   //sigma2_nu_w Update
   sigma2_nu_w(j) = sigma2_nu_update(n,
                                     a_sigma2_nu_w,
                                     b_sigma2_nu_w,
                                     nu_w.col(j));
  
   //theta_w_g Update
   Rcpp::List theta_w_g_output = theta_w_update_tp(z_g,
                                                   zg_trans,
                                                   zgtzg,
                                                   v,
                                                   n,
                                                   w,
                                                   sigma2_epsilon(j),
                                                   theta_w_g.col(j-1),
                                                   sigma2_zeta_w_g(j-1),
                                                   eta_w_g,
                                                   mu_w);
  
   theta_w_g.col(j) = Rcpp::as<arma::vec>(theta_w_g_output[0]);
   mu_w = Rcpp::as<arma::vec>(theta_w_g_output[1]);
  
   //theta_w_r Update
   Rcpp::List theta_w_r_output = theta_w_update_tp(z_r,
                                                   zr_trans,
                                                   zrtzr,
                                                   v,
                                                   n,
                                                   w,
                                                   sigma2_epsilon(j),
                                                   theta_w_r.col(j-1),
                                                   sigma2_zeta_w_r(j-1),
                                                   eta_w_r,
                                                   mu_w);
  
   theta_w_r.col(j) = Rcpp::as<arma::vec>(theta_w_r_output[0]);
   mu_w = Rcpp::as<arma::vec>(theta_w_r_output[1]);
   
   //sigma2_zeta_w_g Update
   sigma2_zeta_w_g(j) = sigma2_zeta_update(v,
                                           n,
                                           a_sigma2_zeta_w_g,
                                           b_sigma2_zeta_w_g,
                                           theta_w_g.col(j),
                                           eta_w_g);
   
   //sigma2_zeta_w_r Update
   sigma2_zeta_w_r(j) = sigma2_zeta_update(v,
                                           n,
                                           a_sigma2_zeta_w_r,
                                           b_sigma2_zeta_w_r,
                                           theta_w_r.col(j),
                                           eta_w_r);
   
   //eta_w_g Update
   Sigma11 = Rcpp::as<arma::mat>(Sigma[j-1])(2,2)*Rcpp::as<arma::mat>(spatial_corr_info[2]);
   Sigma12 = join_horiz(Rcpp::as<arma::mat>(Sigma[j-1])(2,0)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(2,1)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(2,3)*Rcpp::as<arma::mat>(spatial_corr_info[2]));
   
   index(0) = 0; index(1) = 1; index(2) = 3;
   Sigma22_inv = kron(inv_sympd(Rcpp::as<arma::mat>(Sigma[j-1]).submat(index, index)), 
                      Rcpp::as<arma::mat>(spatial_corr_info[0]));
   
   eta_other = join_cols(eta_z_g,
                         eta_z_r,
                         eta_w_r);
   
   eta_w_g = eta_update_tp(v_trans,
                           vtv,
                           m,
                           theta_w_g.col(j),
                           sigma2_zeta_w_g(j),
                           eta_other,
                           Sigma11,
                           Sigma12,
                           Sigma22_inv);
   
   //eta_w_r Update
   Sigma11 = Rcpp::as<arma::mat>(Sigma[j-1])(3,3)*Rcpp::as<arma::mat>(spatial_corr_info[2]);
   Sigma12 = join_horiz(Rcpp::as<arma::mat>(Sigma[j-1])(3,0)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(3,1)*Rcpp::as<arma::mat>(spatial_corr_info[2]),
                        Rcpp::as<arma::mat>(Sigma[j-1])(3,2)*Rcpp::as<arma::mat>(spatial_corr_info[2]));
   
   index(0) = 0; index(1) = 1; index(2) = 2;
   Sigma22_inv = kron(inv_sympd(Rcpp::as<arma::mat>(Sigma[j-1]).submat(index, index)), 
                      Rcpp::as<arma::mat>(spatial_corr_info[0]));
   
   eta_other = join_cols(eta_z_g,
                         eta_z_r,
                         eta_w_g);
   
   eta_w_r = eta_update_tp(v_trans,
                           vtv,
                           m,
                           theta_w_r.col(j),
                           sigma2_zeta_w_r(j),
                           eta_other,
                           Sigma11,
                           Sigma12,
                           Sigma22_inv);
   
   //Sigma Update
   Rcpp::List Sigma_output = Sigma_update_tp(m,
                                             Omega_Sigma_inv,
                                             nu_Sigma_inv,
                                             eta_z_g,
                                             eta_z_r,
                                             eta_w_g,
                                             eta_w_r,
                                             spatial_corr_info[0]);
  
   Sigma[j] = Rcpp::as<arma::mat>(Sigma_output[0]);
   Sigma_inv = Rcpp::as<arma::mat>(Sigma_output[1]);
  
   //phi Update
   eta_full = join_cols(eta_z_g,
                        eta_z_r,
                        eta_w_g,
                        eta_w_r);
     
   Rcpp::List phi_output = phi_update_tp(spatial_dists,
                                         m,
                                         a_phi,
                                         b_phi,
                                         eta_full,
                                         Sigma_inv,
                                         phi(j-1),
                                         spatial_corr_info,
                                         metrop_var_phi_trans,
                                         acctot_phi_trans);

   phi(j) = Rcpp::as<double>(phi_output[0]);
   acctot_phi_trans = phi_output[1];
   spatial_corr_info = phi_output[2];

   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update_tp(y,
                                                  sigma2_epsilon(j),
                                                  mu_z,
                                                  mu_w);
  
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   //Reminder Message About Sorting
   if(j == 1){
     
     Rcpp::Rcout << "###################################################################" << std::endl;
     Rcpp::Rcout << "Data Must Be Sorted As: (1,2), (1,3), ..., (2,1), (2,3),...(n,n-1)." << std::endl;
     Rcpp::Rcout << "###################################################################" << std::endl;
      
     }
   
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     
     double accrate_nu_z_min = round(100*(min(acctot_nu_z)/(double)j));
     Rcpp::Rcout << "nu_z Acceptance (min): " << accrate_nu_z_min << "%" << std::endl;
     
     double accrate_nu_z_max = round(100*(max(acctot_nu_z)/(double)j));
     Rcpp::Rcout << "nu_z Acceptance (max): " << accrate_nu_z_max << "%" << std::endl;
     
     double accrate_nu_w_min = round(100*(min(acctot_nu_w)/(double)j));
     Rcpp::Rcout << "nu_w Acceptance (min): " << accrate_nu_w_min << "%" << std::endl;
     
     double accrate_nu_w_max = round(100*(max(acctot_nu_w)/(double)j));
     Rcpp::Rcout << "nu_w Acceptance (max): " << accrate_nu_w_max << "%" << std::endl;
     
     double accrate_phi_trans = round(100*(acctot_phi_trans/(double)j));
     Rcpp::Rcout << "phi Acceptance: " << accrate_phi_trans << "%" << std::endl;
     Rcpp::Rcout << "**************************" << std::endl;
    
     }
  
   }
        
Rcpp::List regression_info = Rcpp::List::create(Rcpp::Named("beta_z") = beta_z,
                                                Rcpp::Named("gamma_z_g") = gamma_z_g,
                                                Rcpp::Named("gamma_z_r") = gamma_z_r,
                                                Rcpp::Named("beta_w") = beta_w,
                                                Rcpp::Named("gamma_w_g") = gamma_w_g,
                                                Rcpp::Named("gamma_w_r") = gamma_w_r);

Rcpp::List metropolis_info = Rcpp::List::create(Rcpp::Named("acctot_nu_z") = acctot_nu_z,
                                                Rcpp::Named("acctot_nu_w") = acctot_nu_w,
                                                Rcpp::Named("acctot_phi_trans") = acctot_phi_trans);

return Rcpp::List::create(Rcpp::Named("nu_z") = nu_z,
                          Rcpp::Named("sigma2_nu_z") = sigma2_nu_z,
                          Rcpp::Named("theta_z_g") = theta_z_g,
                          Rcpp::Named("theta_z_r") = theta_z_r,
                          Rcpp::Named("sigma2_zeta_z_g") = sigma2_zeta_z_g,
                          Rcpp::Named("sigma2_zeta_z_r") = sigma2_zeta_z_r,
                          Rcpp::Named("sigma2_epsilon") = sigma2_epsilon,
                          Rcpp::Named("nu_w") = nu_w,
                          Rcpp::Named("sigma2_nu_w") = sigma2_nu_w,
                          Rcpp::Named("theta_w_g") = theta_w_g,
                          Rcpp::Named("theta_w_r") = theta_w_r,
                          Rcpp::Named("sigma2_zeta_w_g") = sigma2_zeta_w_g,
                          Rcpp::Named("sigma2_zeta_w_r") = sigma2_zeta_w_r,
                          Rcpp::Named("Sigma") = Sigma,
                          Rcpp::Named("phi") = phi,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("regression_info") = regression_info,
                          Rcpp::Named("metropolis_info") = metropolis_info);

}

