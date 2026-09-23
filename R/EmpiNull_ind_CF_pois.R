#' Estimate Individualized Empirical Null Parameters with Confounders (Poisson Setting)
#'
#' A variant of \code{\link{empirical_null_ind_CF}} for settings where the
#' outcome follows a Poisson distribution and the effective sample size is
#' characterised by the total expected count \eqn{E_i} rather than an explicit
#' sample size \eqn{n_i}.  The confounder coefficients (\code{zeta}), the
#' variance component (\code{phi}), and the null proportion (\code{p0}) are
#' estimated by a grid search over \code{p0}. At each grid point, \code{zeta}
#' and \code{phi} are optimized jointly using L-BFGS-B. The Gaussian working
#' null has mean \eqn{\sqrt{E_i}\,x_i^\top\zeta} and variance
#' \eqn{1 + \phi E_i}.
#'
#' The truncation interval for defining the null set is determined with the
#' following priority:
#' 1. If \code{range} is provided, a sample-size-dependent interval based on
#'    \eqn{E_i} is used.
#' 2. If \code{range} is missing but \code{sigma_est} is provided, a uniform
#'    interval is used.
#' 3. If both are missing, \code{range} is estimated from a robust regression
#'    of \code{z} on the columns of \code{x_bar} multiplied by \code{sqrt(E)}.
#'
#' @param z A numeric vector of Z-scores.
#' @param E A numeric vector of total expected counts
#'   \eqn{E_i = \sum_j \exp(X_{ij}^\top \beta)}, serving as the effective
#'   sample size for each facility.
#' @param E_2 A numeric vector of second-moment expected counts. It is validated
#'   for backward compatibility but is not used by this simplified
#'   \eqn{1+\varphi E_i} approximation. Use
#'   \code{\link{empirical_null_ind_pois}} when the second-moment term is needed.
#' @param x_bar A numeric matrix or vector of confounder means.  Each row
#'   corresponds to one facility.
#' @param p_grid A numeric vector of candidate null proportions to search over.
#'   Defaults to \code{seq(0.8, 0.999, 0.005)}.
#' @param phi_grid A numeric vector whose minimum and maximum define the
#'   nonnegative lower and upper bounds for the variance component. Defaults
#'   to \code{c(0, 5)}. Interior entries are not a grid of candidate values.
#' @param cutoff Central normal coverage used to compute per-facility
#'   truncation half-widths. Defaults to \code{0.95}.
#' @param sigma_est Optional.  A pre-estimated standard deviation defining a
#'   uniform truncation half-width for all facilities.
#' @param range Optional.  A pre-estimated range parameter for sample-size-
#'   dependent truncation.  Takes precedence over \code{sigma_est}.
#'
#' @return A list containing:
#'   \item{zeta}{Estimated confounder coefficients.}
#'   \item{phi}{Estimated variance component, so that the working null variance
#'     is \eqn{1+\phi E_i}.}
#'   \item{p0}{Estimated null proportion.}
#'   \item{diagnostics}{Convergence and search-boundary information, including
#'     \code{simplified_poisson = TRUE}.}
#'
#' @seealso \code{\link{empirical_null_ind_CF}} for the non-Poisson version.
#'
#' @importFrom MASS rlm
#' @importFrom stats qnorm median optim optimize
#' @export
empirical_null_ind_CF_pois <- function(z, E, E_2, x_bar,
                                       p_grid   = seq(0.8, 0.999, 0.005),
                                       phi_grid = c(0, 5),
                                       cutoff = 0.95, sigma_est, range) {
  .check_numeric_vector(E_2, "E_2", n = length(z), positive = TRUE)
  args <- list(z = z, n = E, x_bar = x_bar, p_grid = p_grid,
               phi_grid = phi_grid, cutoff = cutoff)
  if (!missing(sigma_est)) args$sigma_est <- sigma_est
  if (!missing(range)) args$range <- range
  fit <- do.call(empirical_null_ind_CF, args)
  fit$diagnostics$simplified_poisson <- TRUE
  fit
}
