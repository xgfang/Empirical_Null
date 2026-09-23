#' Fit an Assumption-Aware Empirical Null
#'
#' Provides the canonical package entry point intended for scripted workflows
#' and guided agents. It fits a supported empirical-null model, computes the
#' provider-specific null mean and standard deviation, and returns diagnostics
#' that indicate when the correction is weakly identified or numerically
#' unstable.
#'
#' @param z Provider-level test statistics.
#' @param size Effective provider sizes. For Poisson outcomes this is \eqn{E_i}.
#' @param covariates Optional provider-level covariates. An intercept is added by
#'   default when a covariate-dependent model is fitted.
#' @param second_moment Required for \code{family = "poisson"}. The vector
#'   \eqn{E_i^{(2)}}.
#' @param family Either \code{"gaussian"} or \code{"poisson"}.
#' @param model One of \code{"auto"}, \code{"overall"},
#'   \code{"individualized"}, \code{"correlated"}, \code{"cre"}, or
#'   \code{"moment"}. The last option is the Poisson exact-moment working null.
#' @param include_intercept Add an intercept when a design matrix has no constant
#'   column.
#' @param p_grid Candidate null proportions.
#' @param phi_bounds Lower and upper bounds for Gaussian variance components.
#' @param cutoff Central coverage used to initialize censored intervals.
#' @param alpha Two-sided flagging level stored in the returned object.
#'
#' @return An object of class \code{empinull_fit}.
#' @export
fit_empirical_null <- function(
    z, size, covariates = NULL, second_moment = NULL,
    family = c("gaussian", "poisson"),
    model = c("auto", "overall", "individualized", "correlated", "cre", "moment"),
    include_intercept = TRUE,
    p_grid = seq(0.5, 0.995, by = 0.005),
    phi_bounds = c(0, 5), cutoff = 0.95, alpha = 0.05) {
  family <- match.arg(family)
  model <- match.arg(model)
  .check_numeric_vector(z, "z")
  .check_numeric_vector(size, "size", n = length(z), positive = TRUE)
  p_grid <- .check_grid(p_grid, "p_grid", 0, 1,
                        lower_closed = FALSE, upper_closed = FALSE)
  .check_numeric_vector(alpha, "alpha", n = 1L, positive = TRUE)
  if (alpha >= 1) stop("alpha must be smaller than 1.", call. = FALSE)
  .check_numeric_vector(phi_bounds, "phi_bounds", nonnegative = TRUE)
  if (length(phi_bounds) != 2L || phi_bounds[2L] <= phi_bounds[1L]) {
    stop("phi_bounds must be two increasing nonnegative values.", call. = FALSE)
  }

  if (model == "auto") {
    if (family == "poisson") {
      model <- "moment"
    } else if (is.null(covariates)) {
      model <- "individualized"
    } else {
      model <- "correlated"
    }
  }
  if (family == "poisson" && model != "moment") {
    stop("Poisson agent fits currently require model = 'moment'.", call. = FALSE)
  }
  if (family == "gaussian" && model == "moment") {
    stop("model = 'moment' is only available for Poisson outcomes.", call. = FALSE)
  }

  design <- NULL
  if (model %in% c("correlated", "moment")) {
    design <- .prepare_agent_design(covariates, length(z), include_intercept)
    identification_design <- design * sqrt(size)
    if (qr(identification_design)$rank < ncol(identification_design)) {
      stop(
        "The provider design is rank deficient, so the requested null mean is not identified.",
        call. = FALSE
      )
    }
  }

  low_fit <- switch(
    model,
    overall = empirical_null_overall(z, p_grid = p_grid),
    individualized = {
      grid <- seq(phi_bounds[1L], phi_bounds[2L], length.out = 251L)
      empirical_null_ind_RQ(z, size, p_grid = p_grid, phi_grid = grid,
                            cutoff = cutoff)
    },
    correlated = empirical_null_ind_CF(
      z, size, design, p_grid = p_grid, phi_grid = phi_bounds,
      cutoff = cutoff
    ),
    cre = {
      lower <- max(phi_bounds[1L], 1e-6)
      grid <- seq(lower, phi_bounds[2L], length.out = 251L)
      empirical_null_ind_CF_CRE(
        z, size, p_grid = p_grid, phi_grid = grid, cutoff = cutoff
      )
    },
    moment = {
      if (is.null(second_moment)) {
        stop("second_moment is required for a Poisson moment fit.",
             call. = FALSE)
      }
      .check_numeric_vector(second_moment, "second_moment", n = length(z),
                            positive = TRUE)
      empirical_null_ind_pois(
        Z = z, Xbar = design, E = size, E2 = second_moment,
        p_grid = p_grid, cutoff = cutoff
      )
    }
  )

  if (model == "overall") {
    null_mean <- rep(unname(low_fit$est["mu"]), length(z))
    null_sd <- rep(unname(low_fit$est["sigma"]), length(z))
    parameters <- list(mu = low_fit$est["mu"], sigma = low_fit$est["sigma"],
                       p0 = low_fit$est["p0"])
  } else if (model == "individualized") {
    null_mean <- rep(if (!is.null(low_fit$mu)) low_fit$mu else 0, length(z))
    null_sd <- sqrt(1 + size * low_fit$phi)
    parameters <- list(phi = low_fit$phi, p0 = low_fit$p0)
    if (!is.null(low_fit$mu)) parameters$mu <- low_fit$mu
  } else if (model == "correlated") {
    null_mean <- sqrt(size) * drop(design %*% low_fit$zeta)
    null_sd <- sqrt(1 + size * low_fit$phi)
    parameters <- list(zeta = low_fit$zeta, phi = low_fit$phi,
                       p0 = low_fit$p0)
  } else if (model == "cre") {
    null_mean <- rep(0, length(z))
    null_sd <- sqrt(size * low_fit$phi)
    parameters <- list(phi = low_fit$phi, p0 = low_fit$p0)
  } else {
    moments <- poisson_null_moments(
      design, size, second_moment, low_fit$zeta,
      low_fit$sigmaAlpha2, low_fit$sigmaEps2
    )
    null_mean <- moments$mean
    null_sd <- moments$sd
    parameters <- list(
      zeta = low_fit$zeta,
      sigmaAlpha2 = low_fit$sigmaAlpha2,
      sigmaEps2 = low_fit$sigmaEps2,
      p0 = low_fit$p0
    )
  }

  standardized <- (z - null_mean) / null_sd
  diagnostics <- .agent_diagnostics(
    standardized = standardized, size = size, design = design,
    family = family, model = model, low = low_fit$diagnostics,
    p = if (is.null(design)) 0L else ncol(design),
    second_moment = second_moment, cutoff = cutoff
  )
  status <- if (length(diagnostics$warning_codes)) "warning" else "ok"

  structure(
    list(
      method = model,
      family = family,
      parameters = parameters,
      null_mean = null_mean,
      null_sd = null_sd,
      standardized = standardized,
      diagnostics = diagnostics,
      status = status,
      alpha = alpha,
      cutoff = cutoff,
      design = design,
      size = size,
      second_moment = second_moment,
      z = z,
      low_level_fit = low_fit,
      call = match.call()
    ),
    class = "empinull_fit"
  )
}

.agent_diagnostics <- function(standardized, size, design, family, model,
                               low, p, second_moment, cutoff) {
  id_design <- if (is.null(design)) NULL else design * sqrt(size)
  rank <- if (is.null(id_design)) NA_integer_ else qr(id_design)$rank
  condition <- .safe_condition_number(id_design)
  central_limit <- .central_z(cutoff)
  central <- standardized[abs(standardized) <= central_limit]
  warning_codes <- character()
  if (!isTRUE(low$converged)) warning_codes <- c(warning_codes, "OPTIMIZATION_FAILED")
  if (isTRUE(low$p0_boundary)) warning_codes <- c(warning_codes, "P0_AT_SEARCH_BOUNDARY")
  if (isTRUE(low$phi_boundary)) warning_codes <- c(warning_codes, "PHI_AT_SEARCH_BOUNDARY")
  if (!is.na(rank) && rank < p) warning_codes <- c(warning_codes, "RANK_DEFICIENT_DESIGN")
  if (is.finite(condition) && condition > 30) warning_codes <- c(warning_codes, "ILL_CONDITIONED_DESIGN")
  if (is.infinite(condition)) warning_codes <- c(warning_codes, "ILL_CONDITIONED_DESIGN")
  if (p > 0L && length(size) / p < 20) warning_codes <- c(warning_codes, "FEW_PROVIDERS_PER_PARAMETER")
  central_skew <- .skewness(central)
  if (is.finite(central_skew) && abs(central_skew) > 0.5) {
    warning_codes <- c(warning_codes, "CENTRAL_SKEWNESS")
  }
  if (length(central) < max(30L, 5L * max(p, 1L))) {
    warning_codes <- c(warning_codes, "SMALL_CENTRAL_SAMPLE")
  }
  variance_condition <- NA_real_
  if (family == "poisson") {
    variance_design <- cbind(second_moment / size, size)
    variance_condition <- .safe_condition_number(scale(variance_design))
    if (!is.finite(variance_condition) || variance_condition > 30) {
      warning_codes <- c(warning_codes, "POISSON_VARIANCE_COMPONENTS_WEAKLY_IDENTIFIED")
    }
  }
  size_group <- cut(rank(size, ties.method = "first"),
                    breaks = stats::quantile(seq_along(size),
                                             probs = seq(0, 1, length.out = 4)),
                    include.lowest = TRUE, labels = FALSE)
  subgroup <- do.call(rbind, lapply(split(standardized, size_group), function(x) {
    data.frame(n = length(x), mean = mean(x), sd = stats::sd(x))
  }))
  rownames(subgroup) <- paste0("size_group_", seq_len(nrow(subgroup)))
  list(
    n_providers = length(size),
    n_parameters = p,
    providers_per_parameter = if (p == 0L) Inf else length(size) / p,
    design_rank = rank,
    design_condition_number = condition,
    poisson_variance_condition_number = variance_condition,
    central_count = length(central),
    central_mean = if (length(central)) mean(central) else NA_real_,
    central_sd = if (length(central) > 1L) stats::sd(central) else NA_real_,
    central_skewness = central_skew,
    size_subgroup_calibration = subgroup,
    convergence = low,
    warning_codes = unique(warning_codes)
  )
}

#' Predict Provider-Specific Empirical-Null Parameters
#'
#' @param object An \code{empinull_fit} object.
#' @param ... Unused.
#' @return A data frame containing the fitted null mean and standard deviation.
#' @export
predict.empinull_fit <- function(object, ...) {
  data.frame(null_mean = object$null_mean, null_sd = object$null_sd)
}

#' Standardize Scores with a Fitted Empirical Null
#'
#' @param object An \code{empinull_fit} object.
#' @param z Optional replacement vector of scores with the same length as the
#'   fitted data.
#' @return A numeric vector of corrected scores.
#' @export
standardize_empirical_null <- function(object, z = NULL) {
  if (!inherits(object, "empinull_fit")) {
    stop("object must inherit from 'empinull_fit'.", call. = FALSE)
  }
  if (is.null(z)) return(object$standardized)
  .check_numeric_vector(z, "z", n = length(object$null_mean))
  (z - object$null_mean) / object$null_sd
}

#' Classify Providers with a Fitted Empirical Null
#'
#' @param object An \code{empinull_fit} object.
#' @param z Optional replacement scores.
#' @param alpha Two-sided flagging level.
#' @return A data frame with corrected scores, p-values, and directions.
#' @export
flag_empirical_null <- function(object, z = NULL, alpha = object$alpha) {
  .check_numeric_vector(alpha, "alpha", n = 1L, positive = TRUE)
  if (alpha >= 1) stop("alpha must be smaller than 1.", call. = FALSE)
  score <- standardize_empirical_null(object, z)
  critical <- stats::qnorm(1 - alpha / 2)
  direction <- ifelse(score < -critical, -1L,
                      ifelse(score > critical, 1L, 0L))
  data.frame(
    standardized = score,
    p_value = 2 * stats::pnorm(-abs(score)),
    direction = direction
  )
}

#' Produce an Audit Record for an Empirical-Null Fit
#'
#' @param object An \code{empinull_fit} object.
#' @return A structured list separating estimated quantities, diagnostics,
#'   assumptions, and interpretation limits.
#' @export
audit_empirical_null <- function(object) {
  if (!inherits(object, "empinull_fit")) {
    stop("object must inherit from 'empinull_fit'.", call. = FALSE)
  }
  list(
    status = object$status,
    method = object$method,
    family = object$family,
    estimates = object$parameters,
    diagnostics = object$diagnostics,
    assumptions = c(
      "Provider scores share the specified conditional null family.",
      "The central interval is dominated by null or near-null providers.",
      "The provider summaries contain enough between-provider variation.",
      if (object$family == "poisson")
        "The Gaussian distribution is a moment-matched working null."
    ),
    interpretation_limit = paste(
      "The fit estimates a decision-relevant provider-level null fingerprint.",
      "It does not identify the latent patient confounder or its causal decomposition."
    )
  )
}

#' @export
print.empinull_fit <- function(x, ...) {
  cat("EmpiNull fit\n")
  cat("  family:", x$family, "\n")
  cat("  model:", x$method, "\n")
  cat("  status:", x$status, "\n")
  cat("  providers:", x$diagnostics$n_providers, "\n")
  if (length(x$diagnostics$warning_codes)) {
    cat("  warnings:", paste(x$diagnostics$warning_codes, collapse = ", "), "\n")
  }
  invisible(x)
}
