#ifndef __GenePair__
#define __GenePair__

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

Rcpp::List spatial_corr_fun(int m,
                            arma::mat spatial_dists,
                            double phi);

double sigma2_epsilon_update_pd(arma::vec y,
                                arma::mat x_pair,
                                arma::mat x_ind,
                                arma::mat z,
                                int n_star,
                                double a_sigma2_epsilon,
                                double b_sigma2_epsilon,
                                arma::vec beta_old,
                                arma::vec gamma_old,
                                arma::vec theta_old);

Rcpp::List delta_update_pd(arma::vec y,
                           arma::mat x_trans,
                           arma::mat xtx,
                           arma::mat z,
                           int p_x,
                           int p_d,
                           arma::mat x_prior,
                           double sigma2_epsilon,
                           arma::vec theta_old);

arma::vec theta_update_pd(arma::vec y,
                          arma::mat x_pair,
                          arma::mat x_ind,
                          arma::mat z_trans,
                          arma::mat ztz,
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

arma::vec eta_update_pd(arma::mat v_trans,
                        arma::mat vtv,
                        int m,
                        arma::vec theta,
                        double sigma2_zeta,
                        double tau2_old,
                        arma::mat corr_inv);

double tau2_update_pd(int m,
                      double a_tau2,
                      double b_tau2,
                      arma::vec eta,
                      arma::mat corr_inv);

Rcpp::List phi_update_pd(arma::mat spatial_dists,
                         int m,
                         double a_phi,
                         double b_phi,
                         arma::vec eta,
                         double tau2,
                         double phi_old,
                         Rcpp::List spatial_corr_info,
                         double metrop_var_phi_trans,
                         int acctot_phi_trans);

double neg_two_loglike_update_pd(arma::vec y,
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

Rcpp::List w_star_update_tp(arma::vec y,
                            int n_star,
                            int n,
                            int p_x,
                            int p_d,
                            arma::vec mu_z);

Rcpp::List delta_z_update_tp(arma::mat x,
                             arma::mat x_trans,
                             int p_x,
                             int p_d,
                             arma::mat x_prior,
                             arma::vec w_star,
                             arma::vec lambda,
                             arma::mat w_star_mat_delta,
                             arma::vec beta_z_old,
                             arma::vec gamma_z_g_old,
                             arma::vec gamma_z_r_old,
                             arma::vec mu_z_old);

Rcpp::List theta_z_update_tp(arma::mat z,
                             arma::mat z_trans,
                             int n,
                             arma::vec w_star,
                             arma::vec lambda,
                             arma::vec w_star_mat_theta,
                             arma::vec theta_z_old,
                             double sigma2_zeta_z_old,
                             arma::vec mu_z_old);

arma::vec eta_update_tp(arma::mat v_trans,
                        arma::mat vtv,
                        int m,
                        arma::vec theta,
                        double sigma2_zeta,
                        arma::vec eta_other,
                        arma::mat Sigma11,
                        arma::mat Sigma12,
                        arma::mat Sigma22_inv);

arma::vec w_update_tp(arma::vec y,
                      int n_star,
                      double sigma2_epsilon_old,
                      arma::vec mu_w);

double sigma2_epsilon_update_tp(int n_star,
                                double a_sigma2_epsilon,
                                double b_sigma2_epsilon,
                                arma::vec w,
                                arma::vec mu_w);

Rcpp::List delta_w_update_tp(arma::mat x,
                             arma::mat x_trans,
                             arma::mat xtx,
                             int p_x,
                             int p_d,
                             arma::mat x_prior,
                             arma::vec w,
                             double sigma2_epsilon,
                             arma::vec beta_w_old,
                             arma::vec gamma_w_g_old,
                             arma::vec gamma_w_r_old,
                             arma::vec mu_w_old);

Rcpp::List theta_w_update_tp(arma::mat z,
                             arma::mat z_trans,
                             arma::mat zgtzg,
                             arma::mat v,
                             int n,
                             arma::vec w,
                             arma::vec sigma2_epsilon,
                             arma::vec theta_w_old,
                             double sigma2_zeta_w_old,
                             arma::vec eta_w_old,
                             arma::vec mu_w_old);

Rcpp::List Sigma_update_tp(int m,
                           arma::mat Omega_Sigma_inv,
                           double nu_Sigma_inv,
                           arma::vec eta_z_g,
                           arma::vec eta_z_r,
                           arma::vec eta_w_g,
                           arma::vec eta_w_r,
                           arma::mat spatial_corr_inv);

Rcpp::List phi_update_tp(arma::mat spatial_dists,
                         int m,
                         double a_phi,
                         double b_phi,
                         arma::vec eta,
                         arma::mat Sigma_inv,
                         double phi_old,
                         Rcpp::List spatial_corr_info,
                         double metrop_var_phi_trans,
                         int acctot_phi_trans);

double neg_two_loglike_update_tp(arma::vec y,
                                 double sigma2_epsilon,
                                 arma::vec mu_z,
                                 arma::vec mu_w);

Rcpp::List Trans_Prob(int mcmc_samples,
                      arma::vec trans_probs,
                      arma::mat x_pair,
                      arma::mat x_ind_g,
                      arma::mat x_ind_r,
                      arma::mat z_g,
                      arma::mat z_r,
                      arma::mat spatial_dists,
                      arma::mat v,
                      double metrop_var_phi_trans,
                      Rcpp::Nullable<double> sigma2_regress_prior,  //Start of Priors
                      Rcpp::Nullable<double> a_sigma2_zeta_z_g_prior,
                      Rcpp::Nullable<double> b_sigma2_zeta_z_g_prior,
                      Rcpp::Nullable<double> a_sigma2_zeta_z_r_prior,
                      Rcpp::Nullable<double> b_sigma2_zeta_z_r_prior,
                      Rcpp::Nullable<double> a_sigma2_epsilon_prior,
                      Rcpp::Nullable<double> b_sigma2_epsilon_prior,
                      Rcpp::Nullable<double> a_sigma2_zeta_w_g_prior,
                      Rcpp::Nullable<double> b_sigma2_zeta_w_g_prior,
                      Rcpp::Nullable<double> a_sigma2_zeta_w_r_prior,
                      Rcpp::Nullable<double> b_sigma2_zeta_w_r_prior,
                      Rcpp::Nullable<Rcpp::NumericMatrix> Omega_Sigma_inv_prior,
                      Rcpp::Nullable<double> nu_Sigma_inv_prior,
                      Rcpp::Nullable<double> a_phi_prior,
                      Rcpp::Nullable<double> b_phi_prior,
                      Rcpp::Nullable<Rcpp::NumericVector> beta_z_init,  //Start of Initial Values
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_z_g_init,
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_z_r_init,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_z_g_init,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_z_r_init,
                      Rcpp::Nullable<double> sigma2_zeta_z_g_init,
                      Rcpp::Nullable<double> sigma2_zeta_z_r_init,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_z_g_init,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_z_r_init,
                      Rcpp::Nullable<double> sigma2_epsilon_init,
                      Rcpp::Nullable<Rcpp::NumericVector> beta_w_init,
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_w_g_init,
                      Rcpp::Nullable<Rcpp::NumericVector> gamma_w_r_init,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_w_g_init,
                      Rcpp::Nullable<Rcpp::NumericVector> theta_w_r_init,
                      Rcpp::Nullable<double> sigma2_zeta_w_g_init,
                      Rcpp::Nullable<double> sigma2_zeta_w_r_init,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_w_g_init,
                      Rcpp::Nullable<Rcpp::NumericVector> eta_w_r_init,
                      Rcpp::Nullable<Rcpp::NumericMatrix> Sigma_init,
                      Rcpp::Nullable<double> phi_init);

#endif // __GenePair__
