#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// // [[Rcpp::export]]
// inline double pnorm_cpp(double x) {
//   return 0.5 * std::erfc(-x / std::sqrt(2.0));
// }
// 
// // [[Rcpp::export]]
// arma::vec pnorm_vec(const arma::vec& x) {
//   arma::vec output = x;
//   output.transform( [](double val) { return pnorm_cpp(val); } );
//   return output;
// }
// 

/*
 * Negative log-likelihood function of the truncated normal distribution
 * based on the overall empirical null method.
 *
 * Parameters:
 *   arg         : two-dimensional vector containing the estimated mean (mu0) and standard deviation (sig0)
 *   z           : vector of Z-scores to be standardized
 *   trunc_lower : lower bound of the assumed truncated support
 *   trunc_upper : upper bound of the assumed truncated support
 *   p0          : probability that a cluster is 'null'
 *
 * Returns:
 *   double: negative log-likelihood value
 */
// [[Rcpp::export]]
double negloglik_overall(const arma::vec& arg, const arma::vec& z,
                         double trunc_lower, double trunc_upper, double p0) {
  double mu0 = arg[0], sig0 = arg[1];
  if (!std::isfinite(mu0) || !std::isfinite(sig0) || sig0 <= 0.0 ||
      !std::isfinite(p0) || p0 <= 0.0 || p0 >= 1.0) {
    return 1e100;
  }
  arma::uvec idx = arma::find((z >= trunc_lower) % (z <= trunc_upper));
  arma::vec z0 = z.elem(idx);
  double N0 = z0.n_elem, N1 = z.n_elem - N0;
  double Q = pnorm_cpp((trunc_upper - mu0) / sig0) - pnorm_cpp((trunc_lower - mu0) / sig0);
  double outside_prob = 1.0 - p0 * Q;
  if (!std::isfinite(outside_prob) || outside_prob <= 0.0) return 1e100;
  double loglik = N0*std::log(p0) - N0*std::log(sig0)
                 - arma::accu(arma::square(z0 - mu0)) / (2.0 * sig0 * sig0)
                 + N1*std::log(outside_prob);
  return -loglik;
}


/*
 * Negative log-likelihood function of the truncated normal distribution based
 * on the individualized empirical null method.
 *
 * Parameters:
 *  arg         : a single value for the variance component
 *  z           : a vector of Z-scores to be standardized
 *  n           : a vector of effective sample sizes corresponding to z
 *  trunc_lower : the lower bounds of the assumed truncated support
 *  trunc_upper : the upper bounds of the assumed truncated support
 *  p0          : the probability that a cluster is 'null'
 *
 * Returns:
 *  double      : negative log-likelihood value
 */
// [[Rcpp::export]]
double negloglik_individualized(double arg, arma::vec z, arma::vec n,
                                arma::vec trunc_lower, arma::vec trunc_upper,
                                double p0) {

  if (!std::isfinite(arg) || arg < 0.0 || !std::isfinite(p0) ||
      p0 <= 0.0 || p0 >= 1.0) return 1e100;

  arma::uvec idx_trunc = arma::find((z >= trunc_lower) && (z <= trunc_upper));

  arma::vec z0 = z.elem(idx_trunc);
  arma::vec n0 = n.elem(idx_trunc);
  double N0 = z0.n_elem;
  
  arma::vec sig0_square = 1.0 + n0 * arg;
  
  n.shed_rows(idx_trunc);
  arma::vec sig1_square = 1.0 + n * arg;

  if (trunc_upper.n_elem > 1) {
    trunc_lower.shed_rows(idx_trunc);
    trunc_upper.shed_rows(idx_trunc);
  }
  
  arma::vec Q = pnorm_vec(trunc_upper / arma::sqrt(sig1_square)) - 
                pnorm_vec(trunc_lower / arma::sqrt(sig1_square));

  arma::vec outside_prob = 1.0 - p0 * Q;
  if (!outside_prob.is_finite() || arma::any(outside_prob <= 0.0)) return 1e100;

  double loglik = N0 * std::log(p0)
                  - arma::accu(arma::log(sig0_square) / 2.0 + arma::square(z0) / (2.0 * sig0_square))
                  + arma::accu(arma::log(outside_prob));

  return -loglik;
}


/*
 * Negative log-likelihood function of the truncated normal distribution based
 * on the individualized empirical null method, adjusting for confounders.
 *
 * Parameters:
 *  arg         : a vector where the first elements are confounder coefficients (beta) 
 *      and the last element is the variance component.
 *  z           : a vector of Z-scores to be standardized.
 *  x_bar       : a matrix of confounders.
 *  n           : a vector of effective sample sizes corresponding to z.
 *  trunc_lower : the lower bounds of the assumed truncated support.
 *  trunc_upper : the upper bounds of the assumed truncated support.
 *  p0          : the probability that a cluster is 'null'.
 *
 * Returns:
 *  double      : negative log-likelihood value
 */
// [[Rcpp::export]]
double negloglik_confounder(const arma::vec& arg, const arma::vec& z, 
                            arma::mat x_bar, arma::vec n,
                            arma::vec trunc_lower, arma::vec trunc_upper,
                            double p0) {

  arma::uword arg_len = arg.n_elem;
  arma::vec beta = arg.subvec(0, arg_len - 2);
  double var_comp = arg(arg_len - 1);
  if (!arg.is_finite() || var_comp < 0.0 || !std::isfinite(p0) ||
      p0 <= 0.0 || p0 >= 1.0) return 1e100;

  arma::uvec idx_trunc = arma::find((z >= trunc_lower) && (z <= trunc_upper));

  arma::vec z0 = z.elem(idx_trunc);
  arma::vec n0 = n.elem(idx_trunc);
  arma::mat x_bar_0 = x_bar.rows(idx_trunc);
  double N0 = z0.n_elem;

  n.shed_rows(idx_trunc);
  x_bar.shed_rows(idx_trunc);
  
  arma::vec sig0_square = 1.0 + n0 * var_comp;
  arma::vec sig1_square = 1.0 + n * var_comp;

  if (trunc_upper.n_elem > 1) {
    trunc_lower.shed_rows(idx_trunc);
    trunc_upper.shed_rows(idx_trunc);
  }

  arma::vec mu_adj_0 = (x_bar_0 * beta) % arma::sqrt(n0);
  arma::vec mu_adj_1 = (x_bar * beta) % arma::sqrt(n);
  
  arma::vec term1 = (trunc_upper - mu_adj_1) / arma::sqrt(sig1_square);
  arma::vec term2 = (trunc_lower - mu_adj_1) / arma::sqrt(sig1_square);
  arma::vec Q = pnorm_vec(term1) - pnorm_vec(term2);
  arma::vec outside_prob = 1.0 - p0 * Q;
  if (!outside_prob.is_finite() || arma::any(outside_prob <= 0.0)) return 1e100;

  double loglik_part1 = N0 * std::log(p0);
  double loglik_part2 = arma::accu(
    arma::log(sig0_square) / 2.0 + 
    arma::square(z0 - mu_adj_0) / (2.0 * sig0_square)
  );
  double loglik_part3 = arma::accu(arma::log(outside_prob));

  double loglik = loglik_part1 - loglik_part2 + loglik_part3;

  return -loglik;
}



/*
 * Negative log-likelihood function of the truncated normal distribution based
 * on the individualized empirical null method.
 *
 * Parameters:
 *   arg         : a vector of coefficients for the confounders
 *   z           : a vector of Z-scores to be standardized
 *   x_bar       : a matrix of confounders
 *   n           : a vector of effective sample sizes corresponding to z
 *   trunc_lower : the lower bounds of the assumed truncated support
 *   trunc_upper : the upper bounds of the assumed truncated support
 *   p0          : the probability that a cluster is 'null'
 *   phi         : a scalar variance component
 *
 * Returns:
 *   double: negative log-likelihood value
 */
// [[Rcpp::export]]
double negloglik_confounder_gs(const arma::vec& arg, const arma::vec& z,
                               arma::mat x_bar, arma::vec n,
                               arma::vec trunc_lower, arma::vec trunc_upper,
                               double p0, double phi) {

  if (!arg.is_finite() || !std::isfinite(phi) || phi < 0.0 ||
      !std::isfinite(p0) || p0 <= 0.0 || p0 >= 1.0) return 1e100;

  arma::mat x_bar_sqrt_n = x_bar.each_col() % arma::sqrt(n);
  
  arma::uvec idx_trunc = arma::find((z >= trunc_lower) && (z <= trunc_upper));

  arma::vec z0 = z.elem(idx_trunc);
  arma::vec n0 = n.elem(idx_trunc);
  double N0 = z0.n_elem;
  arma::mat x_bar_sqrt_n_0 = x_bar_sqrt_n.rows(idx_trunc);
  arma::vec mu_all = x_bar_sqrt_n * arg;
  
  n.shed_rows(idx_trunc);
  mu_all.shed_rows(idx_trunc);
  
  arma::vec sig0_square = 1.0 + n0 * phi;
  arma::vec sig1_square = 1.0 + n * phi;

  if (trunc_upper.n_elem > 1) {
    trunc_lower.shed_rows(idx_trunc);
    trunc_upper.shed_rows(idx_trunc);
  }
  
  arma::vec Q = pnorm_vec((trunc_upper - mu_all) / arma::sqrt(sig1_square)) - 
                pnorm_vec((trunc_lower - mu_all) / arma::sqrt(sig1_square));
  arma::vec outside_prob = 1.0 - p0 * Q;
  if (!outside_prob.is_finite() || arma::any(outside_prob <= 0.0)) return 1e100;

  arma::vec mu_adj_0 = x_bar_sqrt_n_0 * arg;
  
  double loglik = N0 * std::log(p0)
                  - arma::accu(arma::log(sig0_square) / 2.0 + arma::square(z0 - mu_adj_0) / (2.0 * sig0_square))
                  + arma::accu(arma::log(outside_prob));

  return -loglik;
}


/*
 * Negative log-likelihood function of the truncated normal distribution based
 * on the individualized empirical null method.
 *
 * Parameters:
 * arg         : a single value for the variance component
 * z           : a vector of Z-scores to be standardized
 * n           : a vector of effective sample sizes corresponding to z
 * trunc_lower : the lower bounds of the assumed truncated support
 * trunc_upper : the upper bounds of the assumed truncated support
 * p0          : the probability that a cluster is 'null'
 *
 * Returns:
 * double      : negative log-likelihood value
 */
// [[Rcpp::export]]
double negloglik_confounder_CRE(double arg, const arma::vec& z, 
                                arma::vec n, arma::vec trunc_lower, 
                                arma::vec trunc_upper, double p0) {

  if (!std::isfinite(arg) || arg <= 0.0 || !std::isfinite(p0) ||
      p0 <= 0.0 || p0 >= 1.0) return 1e100;

  arma::uvec idx_trunc = arma::find((z >= trunc_lower) && (z <= trunc_upper));

  arma::vec z0 = z.elem(idx_trunc);
  arma::vec n0 = n.elem(idx_trunc);
  double N0 = z0.n_elem;

  n.shed_rows(idx_trunc);

  arma::vec sig0_square = n0 * arg;
  arma::vec sig1_square = n * arg;

  if (trunc_upper.n_elem > 1) {
    trunc_lower.shed_rows(idx_trunc);
    trunc_upper.shed_rows(idx_trunc);
  }
  
  arma::vec Q = pnorm_vec(trunc_upper / arma::sqrt(sig1_square)) - 
                pnorm_vec(trunc_lower / arma::sqrt(sig1_square));
  arma::vec outside_prob = 1.0 - p0 * Q;
  if (!outside_prob.is_finite() || arma::any(outside_prob <= 0.0)) return 1e100;
  
  double loglik = N0 * std::log(p0)
                  - arma::accu(arma::log(sig0_square) / 2.0 + arma::square(z0) / (2.0 * sig0_square))
                  + arma::accu(arma::log(outside_prob));

  return -loglik;
}



/*
 * Negative log-likelihood function of the truncated normal distribution based
 * on the Poisson Z-scores
 *
 * Parameters:
 *   theta       : parameter vector; first p entries are zeta (regression coefficients),
 *                 then sa2 = sigma_alpha^2, then se2 = sigma_epsilon^2
 *   z           : vector of observed Z-scores
 *   Xbar        : matrix (n_clusters × p) of cluster means
 *   E           : vector E_i = Σ_j exp(X_ij^T β)
 *   E2          : vector E_i^(2) = Σ_j exp(X_ij^T β)^2
 *   trunc_lower : scalar or vector lower truncation bounds
 *   trunc_upper : scalar or vector upper truncation bounds
 *   p0          : probability that a cluster is 'null'
 *
 * Returns:
 *   double      : negative log-likelihood value
 */
// [[Rcpp::export]]
double negloglik_pois2(const arma::vec& theta,
                       const arma::vec& z,
                       const arma::mat& Xbar,
                       const arma::vec& E,
                       const arma::vec& E2,
                       const arma::vec& trunc_lower,
                       const arma::vec& trunc_upper,
                       double p0)
{
  int p = (Xbar.n_cols == 0) ? 1 : static_cast<int>(Xbar.n_cols);
  arma::vec zeta = theta.subvec(0, p - 1);
  double sa2 = theta[p], se2 = theta[p + 1];
  if (!theta.is_finite() || sa2 < 0.0 || se2 < 0.0 ||
      !std::isfinite(p0) || p0 <= 0.0 || p0 >= 1.0 ||
      arma::any(E <= 0.0) || arma::any(E2 <= 0.0)) return 1e100;

  arma::vec eta = Xbar * zeta;

  double shift = 0.5 * sa2 + 0.5 * se2;
  arma::vec mu = (arma::exp(eta + shift) - 1.0) % arma::sqrt(E);

  arma::vec T1 = arma::exp(eta + shift);
  arma::vec T2 = arma::exp(2.0 * eta + sa2 + se2);
  arma::vec sigma2 = T1 + T2 % ( std::exp(sa2) * (std::exp(se2) - 1.0) * (E2 / E)
                               + (std::exp(sa2) - 1.0) * E );
  if (!sigma2.is_finite() || arma::any(sigma2 <= 0.0)) return 1e100;
  arma::vec sigma = arma::sqrt(sigma2);

  arma::uvec idx_in  = arma::find( (z >= trunc_lower) % (z <= trunc_upper) );
  arma::uvec idx_out = arma::find( (z <  trunc_lower) + (z >  trunc_upper) );

  arma::vec term_u = (trunc_upper - mu) / sigma;
  arma::vec term_l = (trunc_lower - mu) / sigma;
  arma::vec Q = pnorm_vec(term_u) - pnorm_vec(term_l);

  arma::vec term0 = std::log(p0) - arma::log(sigma) - arma::square(z - mu) / (2.0 * sigma2);
  double ll0 = arma::accu(term0.elem(idx_in));

  arma::vec outside_prob = 1.0 - p0 * Q.elem(idx_out);
  if (!outside_prob.is_finite() || arma::any(outside_prob <= 0.0)) return 1e100;
  double ll1 = arma::accu(arma::log(outside_prob));

  return -(ll0 + ll1);
}
