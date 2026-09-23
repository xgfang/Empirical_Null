#' Estimate Empirical Null Parameters with Provider-Level Confounding Factors
#'
#' Estimates parameters of the empirical null distribution for facility
#' Z-scores that are confounded by provider-level covariates \code{W}.  The
#' null accounts for both a mean shift driven by \eqn{\nu^\top W_i} and an
#' overdispersion parameter \eqn{\varphi}.  A grid search over the null
#' proportion \code{p0} is used, with the remaining parameters optimised at
#' each grid point.
#'
#' @param Obs Numeric vector of observed event counts.
#' @param Exp Numeric vector of expected event counts.
#' @param W Numeric matrix of provider-level covariates (one row per facility).
#'   Columns are mean-centred internally.
#' @param p.grid Numeric vector of candidate null proportions to search over.
#'   Defaults to \code{seq(0.4, 0.999, 0.001)}.
#' @param cutoff Central coverage in \code{(0, 1)} used to form initial
#'   truncation bounds. Defaults to \code{0.95}.
#' @param family Character string specifying the working model for the mean and
#'   variance.  Either \code{"poisson"} (default) or \code{"approx poisson"}.
#'
#' @return A list with elements:
#'   \item{nu_hat}{Estimated covariate coefficient vector \eqn{\hat{\nu}}
#'   (length \code{ncol(W)}).}
#'   \item{sig2_hat}{Estimated overdispersion \eqn{\hat{\varphi}}.}
#'   \item{p0}{Selected null proportion (grid minimum).}
#'   \item{nu_var}{Estimated covariance matrix of \eqn{\hat{\nu}}.}
#'   \item{W}{The mean-centred design matrix used in fitting.}
#'   \item{out_idx}{Integer indices of facilities outside the truncation window.}
#'
#' @importFrom MASS rlm psi.huber
#' @importFrom stats pnorm dnorm qnorm optim
#' @export
empirical_null_CF_prov <- function(Obs, Exp, W,
                                   p.grid = seq(0.4, 0.999, 0.001),
                                   cutoff = 0.95,
                                   family = "poisson") {
  .check_numeric_vector(Obs, "Obs", nonnegative = TRUE)
  .check_numeric_vector(Exp, "Exp", n = length(Obs), positive = TRUE)
  W <- .as_design_matrix(W, length(Obs), "W")
  p.grid <- .check_grid(p.grid, "p.grid", 0, 1,
                        lower_closed = FALSE, upper_closed = FALSE)
  family <- match.arg(family, c("poisson", "approx poisson"))
  cval <- .central_z(cutoff)
  W_center <- colMeans(W)
  W <- sweep(W, 2L, W_center, "-")

  P     <- ncol(W)
  Z_FE  <- (Obs - Exp) / sqrt(Exp)
  niter <- length(p.grid)
  N     <- length(Z_FE)
  eval.p.grid <- rep(0, niter)
  W_pre <- W * sqrt(Exp)

  get_pdr <- function(Wp, arg) {
    (as.matrix(Wp) %*% arg[seq_len(P)])[, 1]
  }

  get_mean_var <- function(arg, ntilde, W, family) {
    if (family == "poisson") {
      temp <- exp(get_pdr(W, arg) + arg[P + 1] / 2)
      m <- sqrt(ntilde) * (temp - 1)
      v <- temp * (1 + temp * (exp(arg[P + 1]) - 1) * ntilde)
    } else if (family == "approx poisson") {
      pdrW    <- get_pdr(W, arg)
      pdrWpre <- get_pdr(W_pre, arg)
      m <- pdrWpre
      v <- 1 + pdrW + arg[P + 1] * ntilde
    }
    list(m, v)
  }

  # Initial estimates via robust regression
  sigma_est  <- MASS::rlm(Z_FE ~ 0 + W_pre,
                           method    = "M",
                           scale.est = "MAD",
                           psi       = MASS::psi.huber)
  nu_init    <- sigma_est$coefficients
  pdrW       <- get_pdr(W, nu_init)
  varphi_init <- max((sigma_est$s^2 - 1 - mean(pdrW)) / mean(Exp), 0)
  initial    <- c(nu_init, varphi_init)

  # Initial truncation bounds
  temp   <- exp(get_pdr(W, initial) + initial[P + 1] / 2)
  m_init <- sqrt(Exp) * (temp - 1)
  v_init <- temp * (1 + temp * (exp(initial[P + 1]) - 1) * Exp)
  aorig  <- m_init - cval * sqrt(v_init)
  borig  <- m_init + cval * sqrt(v_init)

  # Negative log-likelihood (closure over aorig, borig, Z_FE, Exp, W, W_pre)
  negloglik <- function(arg) {
    idx_in  <- which(Z_FE >= aorig & Z_FE <= borig)
    Z_FE0   <- Z_FE[idx_in]
    Exp0    <- Exp[idx_in]
    W0      <- W[idx_in, , drop = FALSE]

    null_mv <- get_mean_var(arg, Exp0, W0, family)
    m0      <- null_mv[[1]]
    v0      <- null_mv[[2]]
    N0      <- length(Exp0)

    idx_out <- setdiff(seq_along(Z_FE), idx_in)
    Exp1    <- Exp[idx_out]
    W1      <- W[idx_out, , drop = FALSE]
    aorig1  <- aorig[idx_out]
    borig1  <- borig[idx_out]

    out_mv <- get_mean_var(arg, Exp1, W1, family)
    m1     <- out_mv[[1]]
    v1     <- out_mv[[2]]

    v0 <- pmax(v0, 1e-4)
    v1 <- pmax(v1, 1e-4)

    Q <- pnorm(borig1, mean = m1, sd = sqrt(v1)) -
         pnorm(aorig1, mean = m1, sd = sqrt(v1))

    loglik <- N0 * log(p0) +
              sum(log(dnorm(Z_FE0, mean = m0, sd = sqrt(v0)))) +
              sum(log(1 - p0 * Q))

    return(-loglik)
  }

  # Grid search over p0
  for (i in seq_len(niter)) {
    p0 <- p.grid[i]
    eval.p.grid[i] <- optim(par = initial, fn = negloglik,
                            method = "L-BFGS-B",
                            lower = c(rep(-Inf, P), 0),
                            upper = rep(Inf, P + 1L))$value
  }

  p0            <- p.grid[which.min(eval.p.grid)]
  result.optim  <- optim(par = initial, fn = negloglik,
                         method = "L-BFGS-B",
                         lower = c(rep(-Inf, P), 0),
                         upper = rep(Inf, P + 1L))

  # Variance estimate for nu_hat
  idx_in_final <- which(Z_FE >= aorig & Z_FE <= borig)
  Exp0         <- Exp[idx_in_final]
  W0           <- W[idx_in_final, , drop = FALSE]
  W_pre0       <- W_pre[idx_in_final, , drop = FALSE]
  final_mv     <- get_mean_var(result.optim$par, Exp0, W0, family)
  Omega        <- diag(pmax(final_mv[[2]], 1e-8))

  XtX_inv  <- tryCatch(solve(t(W_pre0) %*% W_pre0),
                       error = function(e) MASS::ginv(t(W_pre0) %*% W_pre0))
  nu_var_est <- XtX_inv %*% t(W_pre0) %*% Omega %*% W_pre0 %*% XtX_inv

  return(list(
    nu_hat  = result.optim$par[seq_len(P)],
    sig2_hat = result.optim$par[P + 1],
    p0      = p0,
    nu_var  = nu_var_est,
    W       = W,
    W_center = W_center,
    out_idx = setdiff(seq_along(Z_FE), idx_in_final),
    diagnostics = .make_low_level_diagnostics(
      result.optim$convergence, result.optim$value, p0, p.grid,
      phi = result.optim$par[P + 1L], phi_bounds = c(0, Inf),
      central_count = length(idx_in_final)
    )
  ))
}
