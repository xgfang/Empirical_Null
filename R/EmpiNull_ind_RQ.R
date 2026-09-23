#' Estimate Individualized Empirical Null Parameters with Random Effects
#'
#' Estimates a null proportion (p0) and a variance component (phi) for an
#' individualized empirical null model. The null distribution for facility
#' \eqn{i} is \eqn{N(\mu, 1 + n_i \cdot \phi)} under the null, where \eqn{\mu}
#' is either fixed at 0 or estimated as a shared intercept across all facilities.
#'
#' Estimation is performed via a nested grid search over \code{p0} and
#' \code{phi}.  When \code{common_intercept = TRUE}, a 1-D optimisation over
#' \eqn{\mu} is added as an inner step for each \code{(phi, p0)} grid cell,
#' using the robust-regression intercept as a starting point.
#'
#' The truncation interval for defining the null set is determined with the
#' following priority:
#' 1. If `range` is provided, a sample-size-dependent interval is used.
#' 2. If `range` is missing but `sigma_est` is provided, a uniform interval is used.
#' 3. If both are missing, `range` is estimated using robust regression.
#'
#' @param z A numeric vector of Z-scores.
#' @param n A numeric vector of effective sample sizes corresponding to each Z-score.
#' @param p_grid A numeric vector specifying the grid of null proportions (p0)
#'   to search. Defaults to a sequence from 0.8 to 0.999.
#' @param phi_grid A numeric vector specifying the grid of the variance
#'   component (phi) to search. Defaults to zero plus a sequence from 0.01 to 5
#'   in increments of 0.01.
#' @param cutoff Central coverage used for the initial truncation interval.
#'   Defaults to 0.95, corresponding to a standard-normal critical value of
#'   approximately 1.96.
#' @param common_intercept Logical. If \code{FALSE} (default), the null
#'   distribution is centered at 0. If \code{TRUE}, a shared null mean \eqn{\mu}
#'   is estimated jointly with \code{phi} and \code{p0}.
#' @param sigma_est Optional. A pre-estimated standard deviation. If provided,
#'   it defines a uniform truncation interval for all Z-scores.
#' @param range Optional. A parameter to define a sample-size-dependent
#'   truncation interval. Takes precedence over `sigma_est`.
#'
#' @return A list containing the estimated parameters:
#'   \item{phi}{The estimated variance component.}
#'   \item{p0}{The estimated null proportion.}
#'   \item{mu}{The estimated shared null mean. Returned only when
#'   \code{common_intercept = TRUE}.}
#'
#' @importFrom MASS rlm
#' @importFrom stats qnorm median optimize
#' @export
empirical_null_ind_RQ <- function(z, n, p_grid = seq(0.8, 0.999, 0.005),
                                  phi_grid = c(0, seq(0.01, 5, 0.01)),
                                  cutoff = 0.95,
                                  common_intercept = FALSE,
                                  sigma_est, range) {
  .check_numeric_vector(z, "z")
  .check_numeric_vector(n, "n", n = length(z), positive = TRUE)
  p_grid <- .check_grid(p_grid, "p_grid", 0, 1,
                        lower_closed = FALSE, upper_closed = FALSE)
  phi_grid <- .check_grid(phi_grid, "phi_grid", 0, Inf)
  c_val <- .central_z(cutoff)

  # Always fit rlm(z ~ 1): used for truncation-width estimation and, when
  # common_intercept = TRUE, for initialising the mu search interval.
  rlm_est <- MASS::rlm(z ~ 1)

  # Determine per-facility truncation half-widths
  if (missing(sigma_est) & missing(range)) {
    range_val <- max((rlm_est$s^2 - 1) / median(n), 0)
    xlim <- c_val * sqrt(1 + n * range_val)
  } else if (!missing(sigma_est) & missing(range)) {
    xlim <- c_val * sigma_est
  } else {
    .check_numeric_vector(range, "range", n = 1L, nonnegative = TRUE)
    xlim <- c_val * sqrt(1 + n * range)
  }

  trunc_lower <- -xlim
  trunc_upper <-  xlim

  if (!common_intercept) {
    # ------------------------------------------------------------------
    # Original model: null centered at 0
    # Grid search over (phi, p0); evaluate negloglik at each cell.
    # ------------------------------------------------------------------
    eval_res <- sapply(phi_grid, function(g) {
      eval_p_grid <- sapply(p_grid, function(p) {
        negloglik_individualized(arg = g, z, n, trunc_lower, trunc_upper, p0 = p)
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

  } else {
    # ------------------------------------------------------------------
    # Common-intercept model: null centered at shared mu
    # For each (phi, p0) cell, optimise over mu via 1-D search.
    # Passing (z - mu) to negloglik_individualized is equivalent to
    # fitting N(mu, 1 + n_i * phi) because the C++ function assumes a
    # zero-centred null and symmetric truncation bounds.
    # ------------------------------------------------------------------
    mu_init  <- as.numeric(rlm_est$coefficients)
    mu_range <- 3 * rlm_est$s

    eval_res <- sapply(phi_grid, function(g) {
      eval_p_grid <- sapply(p_grid, function(p) {
        opt <- optimize(
          f        = function(mu) {
            negloglik_individualized(arg = g, z - mu, n,
                                     trunc_lower, trunc_upper, p0 = p)
          },
          interval = c(mu_init - mu_range, mu_init + mu_range)
        )
        c(neglik = opt$objective, mu = opt$minimum)
      })

      best_idx <- which.min(eval_p_grid["neglik", ])
      c(p0     = unname(p_grid[best_idx]),
        neglik = unname(eval_p_grid["neglik", best_idx]),
        mu     = unname(eval_p_grid["mu",     best_idx]))
    })

    phi_selected <- phi_grid[which.min(eval_res["neglik", ])]
    p_selected   <- eval_res["p0",  which.min(eval_res["neglik", ])]
    mu_selected  <- eval_res["mu",  which.min(eval_res["neglik", ])]

    objective <- eval_res["neglik", which.min(eval_res["neglik", ])]
    return(list(
      mu = mu_selected,
      phi = phi_selected,
      p0 = p_selected,
      diagnostics = .make_low_level_diagnostics(
        0L, objective, p_selected, p_grid,
        phi = phi_selected, phi_bounds = base::range(phi_grid),
        central_count = sum(z >= trunc_lower & z <= trunc_upper)
      )
    ))
  }
}
