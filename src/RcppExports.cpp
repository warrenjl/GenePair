// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Patristic
Rcpp::List Patristic(int mcmc_samples, arma::vec log_patristic_distances, arma::mat x_pair, arma::mat x_ind, arma::mat z, arma::mat spatial_dists, arma::mat v, double metrop_var_phi_trans, Rcpp::Nullable<double> a_sigma2_epsilon_prior, Rcpp::Nullable<double> b_sigma2_epsilon_prior, Rcpp::Nullable<double> sigma2_regress_prior, Rcpp::Nullable<double> a_sigma2_zeta_prior, Rcpp::Nullable<double> b_sigma2_zeta_prior, Rcpp::Nullable<double> a_tau2_prior, Rcpp::Nullable<double> b_tau2_prior, Rcpp::Nullable<double> a_phi_prior, Rcpp::Nullable<double> b_phi_prior, Rcpp::Nullable<double> sigma2_epsilon_init, Rcpp::Nullable<Rcpp::NumericVector> beta_init, Rcpp::Nullable<Rcpp::NumericVector> gamma_init, Rcpp::Nullable<Rcpp::NumericVector> theta_init, Rcpp::Nullable<double> sigma2_zeta_init, Rcpp::Nullable<Rcpp::NumericVector> eta_init, Rcpp::Nullable<double> tau2_init, Rcpp::Nullable<double> phi_init);
RcppExport SEXP _GenePair_Patristic(SEXP mcmc_samplesSEXP, SEXP log_patristic_distancesSEXP, SEXP x_pairSEXP, SEXP x_indSEXP, SEXP zSEXP, SEXP spatial_distsSEXP, SEXP vSEXP, SEXP metrop_var_phi_transSEXP, SEXP a_sigma2_epsilon_priorSEXP, SEXP b_sigma2_epsilon_priorSEXP, SEXP sigma2_regress_priorSEXP, SEXP a_sigma2_zeta_priorSEXP, SEXP b_sigma2_zeta_priorSEXP, SEXP a_tau2_priorSEXP, SEXP b_tau2_priorSEXP, SEXP a_phi_priorSEXP, SEXP b_phi_priorSEXP, SEXP sigma2_epsilon_initSEXP, SEXP beta_initSEXP, SEXP gamma_initSEXP, SEXP theta_initSEXP, SEXP sigma2_zeta_initSEXP, SEXP eta_initSEXP, SEXP tau2_initSEXP, SEXP phi_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type mcmc_samples(mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type log_patristic_distances(log_patristic_distancesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_pair(x_pairSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_ind(x_indSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type spatial_dists(spatial_distsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi_trans(metrop_var_phi_transSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_sigma2_epsilon_prior(a_sigma2_epsilon_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_sigma2_epsilon_prior(b_sigma2_epsilon_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_regress_prior(sigma2_regress_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_sigma2_zeta_prior(a_sigma2_zeta_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_sigma2_zeta_prior(b_sigma2_zeta_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_tau2_prior(a_tau2_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_tau2_prior(b_tau2_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_phi_prior(a_phi_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_phi_prior(b_phi_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_epsilon_init(sigma2_epsilon_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type gamma_init(gamma_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_zeta_init(sigma2_zeta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type eta_init(eta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type tau2_init(tau2_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type phi_init(phi_initSEXP);
    rcpp_result_gen = Rcpp::wrap(Patristic(mcmc_samples, log_patristic_distances, x_pair, x_ind, z, spatial_dists, v, metrop_var_phi_trans, a_sigma2_epsilon_prior, b_sigma2_epsilon_prior, sigma2_regress_prior, a_sigma2_zeta_prior, b_sigma2_zeta_prior, a_tau2_prior, b_tau2_prior, a_phi_prior, b_phi_prior, sigma2_epsilon_init, beta_init, gamma_init, theta_init, sigma2_zeta_init, eta_init, tau2_init, phi_init));
    return rcpp_result_gen;
END_RCPP
}
// delta_update
Rcpp::List delta_update(arma::vec y, arma::mat xtx, arma::mat x_trans, arma::mat z, int p_x, int p_d, arma::mat x_prior, double sigma2_epsilon, arma::vec theta_old);
RcppExport SEXP _GenePair_delta_update(SEXP ySEXP, SEXP xtxSEXP, SEXP x_transSEXP, SEXP zSEXP, SEXP p_xSEXP, SEXP p_dSEXP, SEXP x_priorSEXP, SEXP sigma2_epsilonSEXP, SEXP theta_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xtx(xtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_trans(x_transSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< int >::type p_d(p_dSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_prior(x_priorSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_epsilon(sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_old(theta_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(delta_update(y, xtx, x_trans, z, p_x, p_d, x_prior, sigma2_epsilon, theta_old));
    return rcpp_result_gen;
END_RCPP
}
// eta_update
arma::vec eta_update(arma::mat vtv, arma::mat v_trans, int m, arma::vec theta, double sigma2_zeta, double tau2_old, arma::mat corr_inv);
RcppExport SEXP _GenePair_eta_update(SEXP vtvSEXP, SEXP v_transSEXP, SEXP mSEXP, SEXP thetaSEXP, SEXP sigma2_zetaSEXP, SEXP tau2_oldSEXP, SEXP corr_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type vtv(vtvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v_trans(v_transSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_zeta(sigma2_zetaSEXP);
    Rcpp::traits::input_parameter< double >::type tau2_old(tau2_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv(corr_invSEXP);
    rcpp_result_gen = Rcpp::wrap(eta_update(vtv, v_trans, m, theta, sigma2_zeta, tau2_old, corr_inv));
    return rcpp_result_gen;
END_RCPP
}
// neg_two_loglike_update
double neg_two_loglike_update(arma::vec y, arma::mat x_pair, arma::mat x_ind, arma::mat z, int n_star, double sigma2_epsilon, arma::vec beta, arma::vec gamma, arma::vec theta);
RcppExport SEXP _GenePair_neg_two_loglike_update(SEXP ySEXP, SEXP x_pairSEXP, SEXP x_indSEXP, SEXP zSEXP, SEXP n_starSEXP, SEXP sigma2_epsilonSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_pair(x_pairSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_ind(x_indSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n_star(n_starSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_epsilon(sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(neg_two_loglike_update(y, x_pair, x_ind, z, n_star, sigma2_epsilon, beta, gamma, theta));
    return rcpp_result_gen;
END_RCPP
}
// phi_update
Rcpp::List phi_update(arma::mat spatial_dists, int m, double a_phi, double b_phi, Rcpp::List spatial_corr_info, arma::vec eta, double tau2, double phi_old, double metrop_var_phi_trans, int acctot_phi_trans);
RcppExport SEXP _GenePair_phi_update(SEXP spatial_distsSEXP, SEXP mSEXP, SEXP a_phiSEXP, SEXP b_phiSEXP, SEXP spatial_corr_infoSEXP, SEXP etaSEXP, SEXP tau2SEXP, SEXP phi_oldSEXP, SEXP metrop_var_phi_transSEXP, SEXP acctot_phi_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type spatial_dists(spatial_distsSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type a_phi(a_phiSEXP);
    Rcpp::traits::input_parameter< double >::type b_phi(b_phiSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type spatial_corr_info(spatial_corr_infoSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< double >::type phi_old(phi_oldSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi_trans(metrop_var_phi_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_phi_trans(acctot_phi_transSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_update(spatial_dists, m, a_phi, b_phi, spatial_corr_info, eta, tau2, phi_old, metrop_var_phi_trans, acctot_phi_trans));
    return rcpp_result_gen;
END_RCPP
}
// sigma2_epsilon_update
double sigma2_epsilon_update(arma::vec y, arma::mat x_pair, arma::mat x_ind, arma::mat z, int n_star, double a_sigma2_epsilon, double b_sigma2_epsilon, arma::vec beta_old, arma::vec gamma_old, arma::vec theta_old);
RcppExport SEXP _GenePair_sigma2_epsilon_update(SEXP ySEXP, SEXP x_pairSEXP, SEXP x_indSEXP, SEXP zSEXP, SEXP n_starSEXP, SEXP a_sigma2_epsilonSEXP, SEXP b_sigma2_epsilonSEXP, SEXP beta_oldSEXP, SEXP gamma_oldSEXP, SEXP theta_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_pair(x_pairSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_ind(x_indSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n_star(n_starSEXP);
    Rcpp::traits::input_parameter< double >::type a_sigma2_epsilon(a_sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type b_sigma2_epsilon(b_sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_old(gamma_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_old(theta_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma2_epsilon_update(y, x_pair, x_ind, z, n_star, a_sigma2_epsilon, b_sigma2_epsilon, beta_old, gamma_old, theta_old));
    return rcpp_result_gen;
END_RCPP
}
// sigma2_zeta_update
double sigma2_zeta_update(arma::mat v, int n, double a_sigma2_zeta, double b_sigma2_zeta, arma::vec theta, arma::vec eta_old);
RcppExport SEXP _GenePair_sigma2_zeta_update(SEXP vSEXP, SEXP nSEXP, SEXP a_sigma2_zetaSEXP, SEXP b_sigma2_zetaSEXP, SEXP thetaSEXP, SEXP eta_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type a_sigma2_zeta(a_sigma2_zetaSEXP);
    Rcpp::traits::input_parameter< double >::type b_sigma2_zeta(b_sigma2_zetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_old(eta_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma2_zeta_update(v, n, a_sigma2_zeta, b_sigma2_zeta, theta, eta_old));
    return rcpp_result_gen;
END_RCPP
}
// spatial_corr_fun
Rcpp::List spatial_corr_fun(int m, arma::mat spatial_dists, double phi);
RcppExport SEXP _GenePair_spatial_corr_fun(SEXP mSEXP, SEXP spatial_distsSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type spatial_dists(spatial_distsSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(spatial_corr_fun(m, spatial_dists, phi));
    return rcpp_result_gen;
END_RCPP
}
// tau2_update
double tau2_update(int m, double a_tau2, double b_tau2, arma::vec eta, arma::mat corr_inv);
RcppExport SEXP _GenePair_tau2_update(SEXP mSEXP, SEXP a_tau2SEXP, SEXP b_tau2SEXP, SEXP etaSEXP, SEXP corr_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type a_tau2(a_tau2SEXP);
    Rcpp::traits::input_parameter< double >::type b_tau2(b_tau2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv(corr_invSEXP);
    rcpp_result_gen = Rcpp::wrap(tau2_update(m, a_tau2, b_tau2, eta, corr_inv));
    return rcpp_result_gen;
END_RCPP
}
// theta_update
arma::vec theta_update(arma::vec y, arma::mat x_pair, arma::mat x_ind, arma::mat ztz, arma::mat z_trans, arma::mat v, int n, double sigma2_epsilon, arma::vec beta, arma::vec gamma, double sigma2_zeta_old, arma::vec eta_old);
RcppExport SEXP _GenePair_theta_update(SEXP ySEXP, SEXP x_pairSEXP, SEXP x_indSEXP, SEXP ztzSEXP, SEXP z_transSEXP, SEXP vSEXP, SEXP nSEXP, SEXP sigma2_epsilonSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP sigma2_zeta_oldSEXP, SEXP eta_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_pair(x_pairSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_ind(x_indSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ztz(ztzSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z_trans(z_transSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_epsilon(sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_zeta_old(sigma2_zeta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_old(eta_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_update(y, x_pair, x_ind, ztz, z_trans, v, n, sigma2_epsilon, beta, gamma, sigma2_zeta_old, eta_old));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GenePair_Patristic", (DL_FUNC) &_GenePair_Patristic, 25},
    {"_GenePair_delta_update", (DL_FUNC) &_GenePair_delta_update, 9},
    {"_GenePair_eta_update", (DL_FUNC) &_GenePair_eta_update, 7},
    {"_GenePair_neg_two_loglike_update", (DL_FUNC) &_GenePair_neg_two_loglike_update, 9},
    {"_GenePair_phi_update", (DL_FUNC) &_GenePair_phi_update, 10},
    {"_GenePair_sigma2_epsilon_update", (DL_FUNC) &_GenePair_sigma2_epsilon_update, 10},
    {"_GenePair_sigma2_zeta_update", (DL_FUNC) &_GenePair_sigma2_zeta_update, 6},
    {"_GenePair_spatial_corr_fun", (DL_FUNC) &_GenePair_spatial_corr_fun, 3},
    {"_GenePair_tau2_update", (DL_FUNC) &_GenePair_tau2_update, 5},
    {"_GenePair_theta_update", (DL_FUNC) &_GenePair_theta_update, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_GenePair(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
