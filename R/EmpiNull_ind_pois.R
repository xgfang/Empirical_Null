#' Estimate Empirical Null Parameters from Exact Poisson-Lognormal Moments
#'
#' Estimates the null proportion (\code{p0}) and parameters of a
#' moment-matched Gaussian working null whose mean and variance are the exact
#' Poisson-lognormal moments.
#' The null distribution accounts for between-facility (\eqn{\sigma_\alpha^2})
#' and within-facility (\eqn{\sigma_\varepsilon^2}) variance components, as
#' well as covariate effects encoded in \eqn{\zeta}.
#'
#' Estimation is performed via a grid search over \code{p0}.  At each grid
#' point the remaining parameters are optimised using L-BFGS-B.
#'
#' @param Z A numeric vector of Pearson provider scores
#'   \eqn{(O_i-E_i)/\sqrt{E_i}}.
#' @param Xbar A numeric matrix (n_facilities \eqn{\times} p) of facility-level
#'   covariate means.
#' @param E A numeric vector where \eqn{E_i = \sum_j \exp(X_{ij}^\top \beta)}.
#' @param E2 A numeric vector where
#'   \eqn{E_i^{(2)} = \sum_j \exp(2X_{ij}^\top \beta)}. This is the
#'   sum of squared expected contributions, not the square of \eqn{E_i}.
#' @param p_grid A numeric vector of candidate null proportions to search over.
#'   Defaults to \code{seq(0.5, 0.995, 0.005)}.
#' @param init Optional numeric vector of initial parameters
#'   \code{[zeta (length p), sigmaAlpha2, sigmaEps2]}.
#'   If \code{NULL} (default), initial values are constructed automatically via
#'   \code{init_trunc_pois()}.
#' @param cutoff Central normal coverage used to initialise truncation bounds
#'   in \code{init_trunc_pois}. Defaults to \code{0.95}.
#'
#' @return A named list with elements:
#'   \item{p0}{Estimated null proportion.}
#'   \item{negloglik}{Minimum negative log-likelihood.}
#'   \item{zeta}{Named vector of estimated null-mean coefficients.}
#'   \item{zeta01, zeta02, ...}{The same coefficient estimates as individual
#'     elements, retained for compatibility.}
#'   \item{sigmaAlpha2}{Estimated provider-level log-variance component
#'     \eqn{\sigma_\alpha^2}.}
#'   \item{sigmaEps2}{Estimated patient-level log-variance component
#'     \eqn{\sigma_\varepsilon^2}.}
#'   \item{diagnostics}{Convergence, selected objective value, search-boundary
#'     information, and the number of scores inside the central interval.}
#'
#' @importFrom MASS rlm psi.huber
#' @importFrom stats qnorm optim
#' @export
empirical_null_ind_pois <- function(Z, Xbar, E, E2,
                                    p_grid = seq(0.5, 0.995, 0.005),
                                    init   = NULL,
                                    cutoff = 0.95) {
  .check_numeric_vector(Z, "Z")
  .check_numeric_vector(E, "E", n = length(Z), positive = TRUE)
  .check_numeric_vector(E2, "E2", n = length(Z), positive = TRUE)
  Xbar <- .as_design_matrix(Xbar, length(Z), "Xbar")
  p_grid <- .check_grid(p_grid, "p_grid", 0, 1,
                        lower_closed = FALSE, upper_closed = FALSE)
  p <- ncol(Xbar)

  init_info <- init_trunc_pois(Z, Xbar, E, E2, cutoff)
  trunc_l   <- init_info$trunc_lower
  trunc_u   <- init_info$trunc_upper

  if (is.null(init)) {
    init <- c(
      init_info$zeta_init,
      sigmaAlpha2 = init_info$sa2_init,
      sigmaEps2   = init_info$se2_init
    )
  } else {
    .check_numeric_vector(init, "init", n = p + 2L)
    if (any(init[(p + 1L):(p + 2L)] < 0)) {
      stop("The variance components in init must be nonnegative.", call. = FALSE)
    }
  }

  results <- vector("list", length(p_grid))
  current <- init
  for (k in seq_along(p_grid)) {
    fit <- stats::optim(
      par     = current,
      fn      = negloglik_pois2,
      z       = Z,
      Xbar    = Xbar,
      E       = E,
      E2      = E2,
      trunc_lower = trunc_l,
      trunc_upper = trunc_u,
      p0      = p_grid[k],
      method  = "L-BFGS-B",
      lower   = c(rep(-Inf, p), 1e-6, 1e-6),
      upper   = c(rep( Inf, p),  Inf,  Inf)
    )
    results[[k]] <- list(p0 = p_grid[k], fit = fit)
    if (fit$convergence == 0L && all(is.finite(fit$par))) current <- fit$par
  }

  best_idx <- .best_optim_index(lapply(results, function(x) x$fit))
  best_fit <- results[[best_idx]]$fit
  zeta <- best_fit$par[seq_len(p)]
  names(zeta) <- colnames(Xbar)
  out <- list(
    p0 = p_grid[best_idx],
    negloglik = best_fit$value,
    zeta = zeta,
    sigmaAlpha2 = unname(best_fit$par[p + 1L]),
    sigmaEps2 = unname(best_fit$par[p + 2L]),
    diagnostics = .make_low_level_diagnostics(
      best_fit$convergence, best_fit$value, p_grid[best_idx], p_grid,
      central_count = sum(Z >= trunc_l & Z <= trunc_u)
    )
  )
  for (j in seq_len(p)) out[[sprintf("zeta%02d", j)]] <- unname(zeta[j])
  out
}


# Internal: compute initial parameter estimates and truncation bounds for the
# Poisson exact model.
init_trunc_pois <- function(Z, Xbar, E, E2, cutoff = 0.975) {
  c_val <- .central_z(cutoff)
  ratio <- pmax(Z / sqrt(E) + 1, 1e-6)
  sigma_est <- MASS::rlm(log(ratio) ~ Xbar - 1,
                         method = "M", scale.est = "MAD",
                         psi = MASS::psi.huber)
  zeta_init <- as.numeric(sigma_est$coefficients)
  eta0 <- drop(Xbar %*% zeta_init)
  base_var <- max(sigma_est$s^2, 1e-4)
  sa2_init <- base_var / 2
  se2_init <- base_var / 2
  moments <- poisson_null_moments(
    Xbar, E, E2, zeta_init, sa2_init, se2_init
  )
  lower <- moments$mean - c_val * moments$sd
  upper <- moments$mean + c_val * moments$sd

  list(
    zeta_init   = zeta_init,
    sa2_init    = sa2_init,
    se2_init    = se2_init,
    trunc_lower = lower,
    trunc_upper = upper
  )
}

#' Compute the Poisson-Lognormal Working-Null Moments
#'
#' Computes the exact mean and variance of the Pearson provider score under the
#' Poisson-lognormal model. The returned normal distribution is a working null,
#' not an exact distributional result when a shared provider effect is present.
#'
#' @param Xbar Provider-level design matrix.
#' @param E First expected-count summary.
#' @param E2 Second expected-count summary.
#' @param zeta Mean coefficient vector.
#' @param sigmaAlpha2 Provider-level log-variance component.
#' @param sigmaEps2 Patient-level log-variance component.
#'
#' @return A list with vectors \code{mean}, \code{variance}, and \code{sd}.
#' @export
poisson_null_moments <- function(Xbar, E, E2, zeta,
                                 sigmaAlpha2, sigmaEps2) {
  .check_numeric_vector(E, "E", positive = TRUE)
  .check_numeric_vector(E2, "E2", n = length(E), positive = TRUE)
  Xbar <- .as_design_matrix(Xbar, length(E), "Xbar")
  .check_numeric_vector(zeta, "zeta", n = ncol(Xbar))
  .check_numeric_vector(sigmaAlpha2, "sigmaAlpha2", n = 1L,
                        nonnegative = TRUE)
  .check_numeric_vector(sigmaEps2, "sigmaEps2", n = 1L,
                        nonnegative = TRUE)
  eta <- drop(Xbar %*% zeta)
  shift <- (sigmaAlpha2 + sigmaEps2) / 2
  t1 <- exp(eta + shift)
  t2 <- exp(2 * eta + sigmaAlpha2 + sigmaEps2)
  mu <- sqrt(E) * (t1 - 1)
  variance <- t1 + t2 * (
    exp(sigmaAlpha2) * (exp(sigmaEps2) - 1) * E2 / E +
      (exp(sigmaAlpha2) - 1) * E
  )
  list(mean = mu, variance = variance, sd = sqrt(variance))
}
