#' Estimate Individualized Empirical Null Parameters (Correlated Random Effects)
#'
#' A variant of \code{\link{empirical_null_ind_RQ}} for the Correlated Random
#' Effects (CRE) model.  In the CRE model the null variance is proportional to
#' \eqn{n_i} alone (i.e., \eqn{\sigma_i^2 = n_i \cdot \phi}), rather than
#' \eqn{1 + n_i \cdot \phi} as in the standard random-effects model.  This is
#' appropriate when the within-facility baseline variance can be absorbed into
#' the random effect.
#'
#' Estimation is performed via a nested grid search over \code{phi} and
#' \code{p0} to minimise the negative log-likelihood of a truncated mixture
#' model.
#'
#' The truncation interval is determined with the following priority:
#' 1. If \code{range} is provided, a sample-size-dependent interval is used.
#' 2. If \code{range} is missing but \code{sigma_est} is provided, a uniform
#'    interval is used.
#' 3. If both are missing, \code{range} is estimated using robust regression.
#'
#' @param z A numeric vector of Z-scores.
#' @param n A numeric vector of effective sample sizes corresponding to each
#'   Z-score.
#' @param p_grid A numeric vector of candidate null proportions to search over.
#'   Defaults to \code{seq(0.8, 0.999, 0.005)}.
#' @param phi_grid A numeric vector of positive candidate variance components.
#' @param cutoff Central coverage used to compute provider-specific truncation
#'   half-widths. Defaults to \code{0.95}.
#' @param sigma_est Optional.  A pre-estimated standard deviation defining a
#'   uniform truncation half-width for all facilities.
#' @param range Optional.  A pre-estimated range parameter for sample-size-
#'   dependent truncation.  Takes precedence over \code{sigma_est}.
#'
#' @return A list containing:
#'   \item{phi}{Estimated variance component.}
#'   \item{p0}{Estimated null proportion.}
#'
#' @seealso \code{\link{empirical_null_ind_RQ}} for the standard random-effects
#'   model where the null variance is \eqn{1 + n_i \cdot \phi}.
#'
#' @importFrom MASS rlm
#' @importFrom stats qnorm median
#' @export
empirical_null_ind_CF_CRE <- function(z, n,
                                      p_grid   = seq(0.8, 0.999, 0.005),
                                      phi_grid = seq(0.01, 5, 0.01),
                                      cutoff = 0.95, sigma_est, range) {
  .check_numeric_vector(z, "z")
  .check_numeric_vector(n, "n", n = length(z), positive = TRUE)
  p_grid <- .check_grid(p_grid, "p_grid", 0, 1,
                        lower_closed = FALSE, upper_closed = FALSE)
  phi_grid <- .check_grid(phi_grid, "phi_grid", 0, Inf,
                          lower_closed = FALSE)
  c_val <- .central_z(cutoff)

  # Determine per-facility truncation half-widths.
  # Note: CRE variance is n*phi (no "+1"), so truncation uses sqrt(n * range_val).
  if (missing(sigma_est) & missing(range)) {
    rlm_est   <- MASS::rlm(z ~ 1)
    range_val <- max((rlm_est$s^2 - 1) / median(n), min(phi_grid))
    xlim <- c_val * sqrt(n * range_val)
  } else if (!missing(sigma_est) & missing(range)) {
    .check_numeric_vector(sigma_est, "sigma_est", n = 1L, positive = TRUE)
    xlim <- c_val * sigma_est
  } else {
    .check_numeric_vector(range, "range", n = 1L, positive = TRUE)
    xlim <- c_val * sqrt(n * range)
  }

  trunc_lower <- -xlim
  trunc_upper <-  xlim

  eval_res <- sapply(phi_grid, function(g) {
    eval_p_grid <- sapply(p_grid, function(p) {
      negloglik_confounder_CRE(arg = g, z, n, trunc_lower, trunc_upper, p0 = p)
    })
    c(p0 = p_grid[which.min(eval_p_grid)], neglik = min(eval_p_grid))
  })

  phi_selected <- phi_grid[which.min(eval_res["neglik", ])]
  p_selected   <- eval_res["p0", which.min(eval_res["neglik", ])]

  objective <- eval_res["neglik", which.min(eval_res["neglik", ])]
  return(list(
    phi = phi_selected,
    p0 = p_selected,
    diagnostics = .make_low_level_diagnostics(
      0L, objective, p_selected, p_grid,
      phi = phi_selected, phi_bounds = base::range(phi_grid),
      central_count = sum(z >= trunc_lower & z <= trunc_upper)
    )
  ))
}
