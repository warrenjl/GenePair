#ifndef __GenePair__
#define __GenePair__

Rcpp::List spatial_corr_fun(int m,
                            arma::mat spatial_dists,
                            double phi);

double sigma2_epsilon_update(arma::vec y,
                             arma::mat x_pair,
                             arma::mat x_ind,
                             arma::mat z,
                             int n_star,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec beta_old,
                             arma::vec gamma_old,
                             arma::vec theta_old);

Rcpp::List delta_update(arma::vec y,
                        arma::mat xtx,
                        arma::mat x_trans,
                        arma::mat z,
                        int p_x,
                        int p_d,
                        arma::mat x_prior,
                        double sigma2_epsilon,
                        arma::vec theta_old);

arma::vec theta_update(arma::vec y,
                       arma::mat x_pair,
                       arma::mat x_ind,
                       arma::mat ztz,
                       arma::mat z_trans,
                       arma::mat v,
                       int n,
                       double sigma2_epsilon,
                       arma::vec beta,
                       arma::vec gamma,
                       double sigma2_zeta_old,
                       arma::vec eta_old);

double sigma2_zeta_update(arma::mat v,
                          int n,
                          double a_sigma2_zeta,
                          double b_sigma2_zeta,
                          arma::vec theta,
                          arma::vec eta_old);

arma::vec eta_update(arma::mat vtv,
                     arma::mat v_trans,
                     int m,
                     arma::vec theta,
                     double sigma2_zeta,
                     double tau2_old,
                     arma::mat corr_inv);

double tau2_update(int m,
                   double a_tau2,
                   double b_tau2,
                   arma::vec eta,
                   arma::mat corr_inv);

Rcpp::List phi_update(arma::mat spatial_dists,
                      int m,
                      double a_phi,
                      double b_phi,
                      Rcpp::List spatial_corr_info,
                      arma::vec eta,
                      double tau2,
                      double phi_old,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x_pair,
                              arma::mat x_ind,
                              arma::mat z, 
                              int n_star,
                              double sigma2_epsilon,
                              arma::vec beta,
                              arma::vec gamma,
                              arma::vec theta);

Rcpp::List Patristic(int mcmc_samples,
                     arma::vec log_patristic_distances,
                     arma::mat x_pair,
                     arma::mat x_ind,
                     arma::mat z,
                     arma::mat spatial_dists,
                     arma::mat v,
                     double metrop_var_phi_trans,
                     Rcpp::Nullable<double> a_sigma2_epsilon_prior,
                     Rcpp::Nullable<double> b_sigma2_epsilon_prior,
                     Rcpp::Nullable<double> sigma2_regress_prior,
                     Rcpp::Nullable<double> a_sigma2_zeta_prior,
                     Rcpp::Nullable<double> b_sigma2_zeta_prior,
                     Rcpp::Nullable<double> a_tau2_prior,
                     Rcpp::Nullable<double> b_tau2_prior,
                     Rcpp::Nullable<double> a_phi_prior,
                     Rcpp::Nullable<double> b_phi_prior,
                     Rcpp::Nullable<double> sigma2_epsilon_init,
                     Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                     Rcpp::Nullable<Rcpp::NumericVector> gamma_init,
                     Rcpp::Nullable<Rcpp::NumericVector> theta_init,
                     Rcpp::Nullable<double> sigma2_zeta_init,
                     Rcpp::Nullable<Rcpp::NumericVector> eta_init,
                     Rcpp::Nullable<double> tau2_init,
                     Rcpp::Nullable<double> phi_init); 

#endif // __GenePair__
