#' Estimate Parameters of an Overall Empirical Null Distribution
#'
#' This function estimates the parameters of an overall empirical null distribution
#' (mean, standard deviation, and null proportion p0) from a vector of Z-scores.
#' It uses a truncated normal distribution and optimizes the negative log-likelihood
#' over a specified grid of null proportions.
#'
#' @param z A numeric vector of Z-scores.
#' @param p_grid A numeric vector specifying the grid of null proportions (p0)
#'   to search. Defaults to a sequence from 0.8 to 0.999.
#' @param xlim An optional numeric vector of length two specifying the center
#'   and half-width (`c(center, half_width)`) of the truncation interval used
#'   to define the null distribution. If `NULL` (the default), the interval is
#'   calculated automatically based on the median and IQR of `z`.
#'
#' @return A list containing one element:
#'   \item{est}{A numeric vector of length three containing the estimated mean,
#'   standard deviation, and the selected null proportion (p0).}
#'
#' @export
empirical_null_overall <- function(z, p_grid = seq(0.8, 0.999, 0.005), xlim = NULL) {
  .check_numeric_vector(z, "z")
  p_grid <- .check_grid(p_grid, "p_grid", 0, 1,
                        lower_closed = FALSE, upper_closed = FALSE)
  N <- length(z)
  
  if (is.null(xlim)) {
    b <- ifelse(N > 500000, 1, 4.3 * exp(-0.26 * log(N, 10)))
    xlim <- c(stats::median(z),
              b * stats::IQR(z) / (2 * stats::qnorm(0.75)))
  }
  trunc_lower <- xlim[1] - xlim[2]
  trunc_upper <- xlim[1] + xlim[2]
  
  z0 <- z[which(z >= trunc_lower & z <= trunc_upper)]
  if (length(z0) < 3L) {
    stop("The truncation interval contains fewer than three providers.",
         call. = FALSE)
  }
  
  eval_p_grid <- sapply(p_grid, function(p) {
    optim(par = c(mean(z0), max(stats::sd(z0), 1e-4)),
          fn = negloglik_overall,
          z = z, trunc_lower = trunc_lower,
          trunc_upper = trunc_upper, p0 = p,
          hessian = FALSE, method = "L-BFGS-B",
          lower = c(-Inf, 1e-6), upper = c(Inf, Inf))$value
  })
  
  p0_selected <- p_grid[which.min(eval_p_grid)]
  res_optim <- optim(par = c(mean(z0), max(stats::sd(z0), 1e-4)),
                     fn = negloglik_overall,
                     z = z, trunc_lower = trunc_lower,
                     trunc_upper = trunc_upper, p0 = p0_selected,
                     hessian = FALSE, method = "L-BFGS-B",
                     lower = c(-Inf, 1e-6), upper = c(Inf, Inf))
  
  # Corrected to return the selected p0, not the last one from the grid.
  return(list(
    est = c(mu = res_optim$par[1L], sigma = res_optim$par[2L],
            p0 = p0_selected),
    diagnostics = .make_low_level_diagnostics(
      res_optim$convergence, res_optim$value, p0_selected, p_grid,
      central_count = length(z0)
    )
  ))
}
