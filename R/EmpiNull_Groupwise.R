# --- Internal helper ----------------------------------------------------------

# Fit a robust linear model to the Z-scores within each group.
# Returns a data frame with columns: group, intercept, scale.
run_rlm_groupwise <- function(z, group, n_groups,
                              psi = MASS::psi.bisquare,
                              maxit = 1000, acc = 1e-8) {
  purrr::map_dfr(seq_len(n_groups), ~ {
    idx <- which(group == .x)
    res <- MASS::rlm(z[idx] ~ 1,
                     psi   = psi,
                     maxit = maxit,
                     acc   = acc)
    tibble::tibble(group = .x, intercept = res$coefficients, scale = res$s)
  })
}


# --- Exported functions -------------------------------------------------------

#' Compute Mid-P-Value-Based Z-Scores for Poisson Count Data
#'
#' Converts observed and expected Poisson counts into Z-scores using the
#' mid-p-value approach (equations 8 and 12 in Hartman et al., 2024).
#' These Z-scores can be passed directly to \code{\link{empirical_null_groupwise}}
#' or \code{\link{empirical_null_overall}}.
#'
#' @param obs A numeric vector of observed event counts.
#' @param exp A numeric vector of expected event counts (same length as \code{obs}).
#'
#' @return A numeric vector of Z-scores, one per facility.
#'
#' @importFrom stats ppois dpois qnorm
#' @export
cal_Z_htaz <- function(obs, exp) {
  .check_numeric_vector(obs, "obs", nonnegative = TRUE)
  .check_numeric_vector(exp, "exp", n = length(obs), positive = TRUE)
  P_min    <- 2 * ppois(obs, exp) - dpois(obs, exp)
  P_max    <- 2 * (1 - ppois(obs - 1, exp)) - dpois(obs, exp)
  pval_temp <- pmax(1e-6, pmin(P_min, P_max) / 2)
  Z_htaz_f  <- ifelse(P_min <= P_max, qnorm(pval_temp), -qnorm(pval_temp))
  return(Z_htaz_f)
}


#' Estimate Groupwise Empirical Null Parameters
#'
#' Estimates the empirical null distribution separately within exposure-defined
#' strata.  Facilities (or units) are partitioned into \code{n_groups} groups
#' based on quantiles of a size variable (e.g., years at risk or expected
#' counts).  Within each group, a robust linear model (Tukey bisquare) is
#' fitted to the Z-scores to estimate the null mean (intercept) and scale
#' (standard deviation).
#'
#' This approach is motivated by the fact that smaller facilities produce noisier
#' Z-scores.  Estimating the null separately per exposure stratum avoids the
#' large-facility null dominating inference for small facilities.
#'
#' @param z A numeric vector of Z-scores, one per facility.
#' @param size A numeric vector (same length as \code{z}) used to define
#'   exposure groups, e.g., total years at risk or total expected counts.
#'   Facilities are split into \code{n_groups} strata by quantiles of
#'   \code{size}.
#' @param n_groups A positive integer specifying the number of exposure groups.
#'   Defaults to \code{4}.
#' @param psi The influence function passed to \code{\link[MASS]{rlm}}. Defaults
#'   to \code{MASS::psi.bisquare}. Use \code{MASS::psi.huber} together with
#'   \code{maxit = 20} and \code{acc = 1e-4} to reproduce the legacy
#'   \code{rlm(z ~ 1, method = "M")} calls exactly.
#' @param maxit Maximum \code{rlm} iterations. Defaults to \code{1000}.
#' @param acc \code{rlm} convergence tolerance. Defaults to \code{1e-8}.
#'
#' @return A list containing:
#'   \item{intercept}{Numeric vector of length \code{n_groups}. Estimated null
#'   mean for each group.}
#'   \item{scale}{Numeric vector of length \code{n_groups}. Estimated null
#'   standard deviation for each group.}
#'   \item{group}{Integer vector of length \code{length(z)}. Group assignment
#'   (1 to \code{n_groups}) for each facility.}
#'
#' @importFrom MASS rlm psi.bisquare
#' @importFrom stats quantile
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @export
empirical_null_groupwise <- function(z, size, n_groups = 4,
                                     psi = MASS::psi.bisquare,
                                     maxit = 1000, acc = 1e-8) {
  .check_numeric_vector(z, "z")
  .check_numeric_vector(size, "size", n = length(z), positive = TRUE)
  if (length(n_groups) != 1L || !is.finite(n_groups) || n_groups < 2 ||
      n_groups != as.integer(n_groups))
    stop("n_groups must be an integer of at least 2.", call. = FALSE)

  # Assign each facility to a quantile-based stratum defined by size
  probs  <- seq(0, 1, length.out = n_groups + 1)[2:n_groups]
  breaks <- c(-Inf, quantile(size, probs), Inf)
  if (anyDuplicated(breaks)) {
    stop("size has too few distinct values for the requested groups.",
         call. = FALSE)
  }
  group  <- as.integer(cut(size, breaks = breaks, labels = seq_len(n_groups)))

  # Estimate null parameters within each stratum via robust LM
  params <- run_rlm_groupwise(z, group, n_groups, psi = psi,
                              maxit = maxit, acc = acc)

  return(list(
    intercept = params$intercept,
    scale     = params$scale,
    group     = group
  ))
}
