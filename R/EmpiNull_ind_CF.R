#' Estimate Individualized Empirical Null Parameters with Confounding Factors
#'
#' Estimates parameters for an individualized empirical null model that explicitly
#' accounts for confounding factors (\code{x_bar}).  The confounder coefficients
#' (\code{zeta}), a variance component (\code{phi}), and the null proportion
#' (\code{p0}) are estimated by a grid search over \code{p0}. At each grid point,
#' \code{zeta} and \code{phi} are optimized jointly using L-BFGS-B.
#'
#' The truncation interval for defining the null set is determined with the
#' following priority:
#' 1. If \code{range} is provided, a sample-size-dependent interval is used.
#' 2. If \code{range} is missing but \code{sigma_est} is provided, a uniform
#'    interval is used.
#' 3. If both are missing, \code{range} is estimated using robust regression.
#'
#' @param z A numeric vector of Z-scores.
#' @param n A numeric vector of effective sample sizes corresponding to each
#'   Z-score.
#' @param x_bar A numeric matrix or vector of confounders.  Each row corresponds
#'   to one observation in \code{z}.
#' @param p_grid A numeric vector of candidate null proportions to search over.
#'   Defaults to \code{seq(0.8, 0.999, 0.005)}.
#' @param phi_grid A length-two numeric vector giving the lower and upper bounds
#'   for the variance component. Defaults to \code{c(0, 5)}.
#' @param cutoff Central coverage used to compute provider-specific truncation
#'   half-widths. Defaults to \code{0.95}.
#' @param sigma_est Optional.  A pre-estimated standard deviation defining a
#'   uniform truncation half-width for all facilities.
#' @param range Optional.  A pre-estimated range parameter for sample-size-
#'   dependent truncation.  Takes precedence over \code{sigma_est}.
#'
#' @return A list containing:
#'   \item{zeta}{Estimated confounder coefficients.}
#'   \item{phi}{Estimated variance component.}
#'   \item{p0}{Estimated null proportion.}
#'
#' @seealso \code{\link{empirical_null_ind_RQ}} for the no-confounder version,
#'   \code{\link{empirical_null_ind_CF_pois}} for the Poisson variant.
#'
#' @importFrom MASS rlm
#' @importFrom stats qnorm median optim optimize
#' @export
empirical_null_ind_CF <- function(z, n, x_bar,
                                  p_grid   = seq(0.8, 0.999, 0.005),
                                  phi_grid = c(0, 5),
                                  cutoff = 0.95, sigma_est, range) {
  .check_numeric_vector(z, "z")
  .check_numeric_vector(n, "n", n = length(z), positive = TRUE)
  x_bar <- .as_design_matrix(x_bar, length(z), "x_bar")
  p_grid <- .check_grid(p_grid, "p_grid", 0, 1,
                        lower_closed = FALSE, upper_closed = FALSE)
  .check_numeric_vector(phi_grid, "phi_grid", nonnegative = TRUE)
  if (length(phi_grid) < 2L || max(phi_grid) <= min(phi_grid)) {
    stop("phi_grid must provide distinct lower and upper bounds.", call. = FALSE)
  }
  phi_bounds <- base::range(phi_grid)
  c_val <- .central_z(cutoff)

  rlm_est <- MASS::rlm(z ~ x_bar:sqrt(n) - 1)

  # Determine per-facility truncation half-widths
  if (missing(sigma_est) & missing(range)) {
    range_val <- max((rlm_est$s^2 - 1) / median(n), phi_bounds[1L])
    xlim <- c_val * sqrt(1 + n * range_val)
  } else if (!missing(sigma_est) & missing(range)) {
    .check_numeric_vector(sigma_est, "sigma_est", n = 1L, positive = TRUE)
    xlim <- c_val * sigma_est
    range_val <- max((rlm_est$s^2 - 1) / median(n), phi_bounds[1L])
  } else {
    .check_numeric_vector(range, "range", n = 1L, nonnegative = TRUE)
    range_val <- range
    xlim <- c_val * sqrt(1 + n * range)
  }

  initial_mean <- drop((x_bar * sqrt(n)) %*% rlm_est$coefficients)
  trunc_lower <- -xlim + initial_mean
  trunc_upper <-  xlim + initial_mean

  p <- ncol(x_bar)
  current <- c(as.numeric(rlm_est$coefficients),
               min(max(range_val, phi_bounds[1L] + 1e-8), phi_bounds[2L] - 1e-8))
  fits <- vector("list", length(p_grid))
  for (k in seq_along(p_grid)) {
    fits[[k]] <- stats::optim(
      par = current,
      fn = negloglik_confounder,
      z = z, x_bar = x_bar, n = n,
      trunc_lower = trunc_lower, trunc_upper = trunc_upper,
      p0 = p_grid[k], method = "L-BFGS-B",
      lower = c(rep(-Inf, p), phi_bounds[1L]),
      upper = c(rep(Inf, p), phi_bounds[2L])
    )
    if (all(is.finite(fits[[k]]$par))) current <- fits[[k]]$par
  }
  best_idx <- .best_optim_index(fits)
  best <- fits[[best_idx]]
  p_selected <- p_grid[best_idx]
  phi_selected <- unname(best$par[p + 1L])
  zeta <- best$par[seq_len(p)]
  names(zeta) <- colnames(x_bar)

  list(
    zeta = zeta,
    phi = phi_selected,
    p0 = p_selected,
    diagnostics = .make_low_level_diagnostics(
      best$convergence, best$value, p_selected, p_grid,
      phi = phi_selected, phi_bounds = phi_bounds,
      central_count = sum(z >= trunc_lower & z <= trunc_upper)
    )
  )
}
