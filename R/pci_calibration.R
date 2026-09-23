# Calibration stage used by fit_pci().  All scales are standard deviations.
# The existing fit_empirical_null() API is deliberately left unchanged.
.pci_calibrate <- function(
    z, calibration, size = NULL, covariates = NULL, second_moment = NULL,
    null_mean = 0, null_sd = 1, psi = "bisquare", groups = NULL,
    n_groups = 4, grouping = "quantile", common_mean = NULL,
    fit_mask = NULL, control = list()) {
  calibration <- match.arg(calibration, c(
    "theoretical", "fixed", "overall", "groupwise", "individualized",
    "correlated", "cre", "moment"
  ))
  if (!is.numeric(z) || !is.null(dim(z)) || !length(z)) {
    stop("z must be a non-empty numeric vector.", call. = FALSE)
  }
  n <- length(z)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       any(names(control) == "") || anyDuplicated(names(control))))) {
    stop("control must be a list with unique, non-empty names.", call. = FALSE)
  }
  estimator <- control$estimator
  if (is.null(estimator)) estimator <- "robust"
  estimator <- match.arg(estimator, c("robust", "likelihood"))
  robust <- calibration == "groupwise" ||
    (calibration == "overall" && estimator == "robust")
  if (calibration == "groupwise" && estimator != "robust") {
    stop("Groupwise calibration supports estimator = 'robust'.", call. = FALSE)
  }
  if (!is.null(control$estimator) &&
      !calibration %in% c("overall", "groupwise")) {
    stop("control$estimator applies only to overall or groupwise calibration.",
         call. = FALSE)
  }
  if (!is.null(covariates) && !calibration %in% c("correlated", "moment")) {
    stop("covariates apply only to correlated or moment calibration.",
         call. = FALSE)
  }
  if (!is.null(second_moment) && calibration != "moment") {
    stop("second_moment applies only to moment calibration.", call. = FALSE)
  }
  if (!is.null(groups) && calibration != "groupwise") {
    stop("groups apply only to groupwise calibration.", call. = FALSE)
  }
  if (!is.null(common_mean) && !robust) {
    stop("common_mean applies only to robust overall/groupwise calibration.",
         call. = FALSE)
  }
  if (is.null(fit_mask)) fit_mask <- rep(TRUE, n)
  if (!is.logical(fit_mask) || length(fit_mask) != n || anyNA(fit_mask)) {
    stop("fit_mask must be a non-missing logical vector with length(z) entries.",
         call. = FALSE)
  }
  if (!robust && !all(fit_mask)) {
    stop("fit_mask exclusions currently require robust overall/groupwise calibration.",
         call. = FALSE)
  }
  if (!is.null(size)) .check_numeric_vector(size, "size", n = n, positive = TRUE)
  means <- .pci_null_vector(null_mean, n, "null_mean")
  sds <- .pci_null_vector(null_sd, n, "null_sd", positive = TRUE)
  if (calibration != "fixed" && (any(means != 0) || any(sds != 1))) {
    stop("null_mean and null_sd are configurable only for fixed calibration.",
         call. = FALSE)
  }

  if (calibration %in% c("theoretical", "fixed")) {
    .pci_calibration_control(control, character())
    return(list(
      null_mean = means, null_sd = sds, group = rep(1L, n), fit = NULL,
      parameters = list(null_mean = null_mean, null_sd = null_sd),
      diagnostics = list(
        estimator = "specified", n_providers = n, n_fitted = 0L,
        n_nonfinite = sum(!is.finite(z)), parameter_uncertainty = "not propagated"
      ), warning_codes = character()
    ))
  }

  if (robust) {
    .pci_calibration_control(control, c(
      "estimator", "maxit", "acc", "tuning", "min_group_size", "small_group"
    ))
    return(.pci_robust_calibration(
      z, calibration, size, psi, groups, n_groups, grouping, common_mean,
      fit_mask, control
    ))
  }

  allowed <- c("p_grid", "phi_bounds", "cutoff")
  if (calibration == "overall") allowed <- c("estimator", "p_grid", "cutoff")
  if (calibration %in% c("correlated", "moment")) {
    allowed <- c(allowed, "include_intercept")
  }
  if (calibration == "moment") allowed <- setdiff(allowed, "phi_bounds")
  if (calibration == "individualized") allowed <- c(allowed, "common_intercept")
  .pci_calibration_control(control, allowed)
  .check_numeric_vector(z, "z")
  if (n < 3L) stop("Structured calibration requires at least three providers.",
                   call. = FALSE)
  if (is.null(size)) {
    if (calibration != "overall") {
      stop("size is required for ", calibration, " calibration.", call. = FALSE)
    }
    size <- rep(1, n)
  }
  if (calibration == "moment" && is.null(second_moment)) {
    stop("second_moment is required for moment calibration.", call. = FALSE)
  }
  if (!is.null(control$include_intercept) &&
      (!is.logical(control$include_intercept) ||
       length(control$include_intercept) != 1L || is.na(control$include_intercept))) {
    stop("control$include_intercept must be TRUE or FALSE.", call. = FALSE)
  }
  if (calibration %in% c("correlated", "moment") &&
      identical(control$include_intercept, FALSE) && is.null(covariates)) {
    stop("covariates are required when control$include_intercept = FALSE for ",
         calibration, " calibration.", call. = FALSE)
  }
  if (!is.null(control$common_intercept) &&
      (!is.logical(control$common_intercept) ||
       length(control$common_intercept) != 1L || is.na(control$common_intercept))) {
    stop("control$common_intercept must be TRUE or FALSE.", call. = FALSE)
  }
  args <- list(
    z = z, size = size, covariates = covariates, second_moment = second_moment,
    family = if (calibration == "moment") "poisson" else "gaussian",
    model = calibration
  )
  for (key in intersect(names(control), c(
    "p_grid", "phi_bounds", "cutoff", "include_intercept"
  ))) args[[key]] <- control[[key]]
  if (calibration == "individualized" && isTRUE(control$common_intercept)) {
    fit <- .pci_common_intercept_fit(args)
  } else {
    fit <- do.call(fit_empirical_null, args)
  }
  .pci_null_vector(fit$null_mean, n, "fitted null mean")
  .pci_null_vector(fit$null_sd, n, "fitted null standard deviation", positive = TRUE)
  diagnostics <- fit$diagnostics
  diagnostics$estimator <- "likelihood"
  diagnostics$n_fitted <- n
  diagnostics$parameter_uncertainty <- "not propagated"
  list(
    null_mean = fit$null_mean, null_sd = fit$null_sd, group = rep(1L, n),
    fit = fit, parameters = fit$parameters, diagnostics = diagnostics,
    warning_codes = fit$diagnostics$warning_codes
  )
}

.pci_calibration_control <- function(control, allowed) {
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    stop("Unsupported calibration control option(s): ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
}

.pci_null_vector <- function(x, n, name, positive = FALSE) {
  .check_numeric_vector(x, name, positive = positive)
  if (!length(x) %in% c(1L, n)) {
    stop(name, " must have length 1 or length(z).", call. = FALSE)
  }
  rep(as.numeric(x), length.out = n)
}

.pci_calibration_groups <- function(size, groups, n, n_groups, grouping) {
  if (!is.null(groups)) {
    if (!is.atomic(groups) || !is.null(dim(groups)) || length(groups) != n ||
        anyNA(groups)) {
      stop("groups must be a non-missing vector with length(z) entries.",
           call. = FALSE)
    }
    return(as.character(groups))
  }
  if (!is.numeric(n_groups) || length(n_groups) != 1L || !is.finite(n_groups) ||
      n_groups < 1 || n_groups > n || n_groups != floor(n_groups)) {
    stop("n_groups must be an integer between 1 and length(z).", call. = FALSE)
  }
  if (n_groups == 1L) return(rep(1L, n))
  if (is.null(size)) stop("size is required to construct groups.", call. = FALSE)
  grouping <- match.arg(grouping, c("quantile", "rank"))
  if (grouping == "rank") {
    ord <- order(size, seq_len(n))
    group <- integer(n)
    group[ord] <- ceiling(seq_len(n) / ceiling(n / n_groups))
    return(group)
  }
  breaks <- stats::quantile(size, seq(0, 1, length.out = n_groups + 1L),
                            names = FALSE, type = 7)
  if (anyDuplicated(breaks)) {
    stop(paste(
      "size has tied quantile breaks; use fewer groups, grouping = 'rank',",
      "or explicit groups."
    ), call. = FALSE)
  }
  as.integer(cut(size, breaks = breaks, include.lowest = TRUE, right = TRUE))
}

.pci_robust_calibration <- function(
    z, calibration, size, psi, groups, n_groups, grouping, common_mean,
    fit_mask, control) {
  psi <- match.arg(psi, c("bisquare", "huber"))
  maxit <- if (is.null(control$maxit)) 1000L else control$maxit
  acc <- if (is.null(control$acc)) 1e-8 else control$acc
  tuning <- if (is.null(control$tuning)) {
    if (psi == "huber") 1.345 else 4.685
  } else control$tuning
  min_n <- if (is.null(control$min_group_size)) 3L else control$min_group_size
  small <- if (is.null(control$small_group)) "error" else control$small_group
  small <- match.arg(small, c("error", "theoretical"))
  .check_numeric_vector(acc, "control$acc", n = 1L, positive = TRUE)
  .check_numeric_vector(tuning, "control$tuning", n = 1L, positive = TRUE)
  .check_numeric_vector(maxit, "control$maxit", n = 1L, positive = TRUE)
  .check_numeric_vector(min_n, "control$min_group_size", n = 1L, positive = TRUE)
  if (maxit != floor(maxit) || min_n < 2 || min_n != floor(min_n)) {
    stop("control$maxit must be an integer; min_group_size must be an integer >= 2.",
         call. = FALSE)
  }
  group <- if (calibration == "overall") rep(1L, length(z)) else
    .pci_calibration_groups(size, groups, length(z), n_groups, grouping)
  labels <- unique(group)
  fits <- stats::setNames(vector("list", length(labels)), as.character(labels))
  params <- data.frame(
    group = as.character(labels), fitted_mean = NA_real_, null_mean = NA_real_,
    null_sd = NA_real_, n_fitted = 0L, converged = FALSE,
    fallback = FALSE, stringsAsFactors = FALSE
  )
  warnings <- character()
  codes <- character()
  eligible <- fit_mask & is.finite(z)
  if (psi == "huber") {
    psi_function <- function(u, deriv = 0) MASS::psi.huber(u, k = tuning, deriv = deriv)
  } else {
    psi_function <- function(u, deriv = 0) MASS::psi.bisquare(u, c = tuning, deriv = deriv)
  }
  for (k in seq_along(labels)) {
    idx <- which(group == labels[k] & eligible)
    params$n_fitted[k] <- length(idx)
    if (length(idx) < min_n) {
      if (small == "error") {
        stop("Group '", labels[k], "' has ", length(idx),
             " eligible providers; at least ", min_n, " are required.", call. = FALSE)
      }
      params$fitted_mean[k] <- params$null_mean[k] <- 0
      params$null_sd[k] <- 1
      params$fallback[k] <- TRUE
      codes <- c(codes, "SMALL_GROUP_THEORETICAL_FALLBACK")
      next
    }
    values <- z[idx]
    fit <- withCallingHandlers(
      MASS::rlm(values ~ 1, psi = psi_function, method = "M", scale.est = "MAD",
                maxit = maxit, acc = acc),
      warning = function(w) {
        warnings <<- c(warnings, paste0("Group '", labels[k], "': ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
    mu <- unname(stats::coef(fit)[1L])
    if (!is.finite(mu) || !is.finite(fit$s) || fit$s <= 0) {
      stop("Group '", labels[k], "' has a nonpositive or nonfinite fitted scale/mean.",
           call. = FALSE)
    }
    fits[[k]] <- fit
    params$fitted_mean[k] <- params$null_mean[k] <- mu
    params$null_sd[k] <- fit$s
    params$converged[k] <- isTRUE(fit$converged)
    if (!isTRUE(fit$converged)) codes <- c(codes, "ROBUST_DID_NOT_CONVERGE")
  }
  # Means can be pooled or specified without refitting each group's scale.
  if (!is.null(common_mean)) {
    if (isTRUE(common_mean)) {
      if (!any(eligible)) stop("No eligible finite scores for a pooled mean.", call. = FALSE)
      params$null_mean[] <- mean(z[eligible])
    } else if (!identical(common_mean, FALSE)) {
      .check_numeric_vector(common_mean, "common_mean", n = 1L)
      params$null_mean[] <- common_mean
    }
  }
  which_group <- match(group, labels)
  list(
    null_mean = params$null_mean[which_group],
    null_sd = params$null_sd[which_group], group = group, fit = fits,
    parameters = params,
    diagnostics = list(
      estimator = "robust", psi = psi, tuning = tuning, maxit = maxit, acc = acc,
      grouping = if (calibration == "overall") "overall" else if (!is.null(groups))
        "supplied" else grouping,
      n_providers = length(z), n_fitted = sum(eligible),
      n_excluded = sum(!fit_mask), n_nonfinite = sum(!is.finite(z)),
      groups = params, common_mean = common_mean,
      pooled_mean_uses = "eligible finite providers",
      warnings = unique(warnings), warning_codes = unique(codes),
      parameter_uncertainty = "not propagated"
    ), warning_codes = unique(codes)
  )
}

# The legacy unified API does not expose the existing common-intercept option.
# Assemble its usual object directly while retaining the low-level estimator.
.pci_common_intercept_fit <- function(args) {
  p_grid <- if (is.null(args$p_grid)) seq(0.5, 0.995, by = 0.005) else args$p_grid
  bounds <- if (is.null(args$phi_bounds)) c(0, 5) else args$phi_bounds
  cutoff <- if (is.null(args$cutoff)) 0.95 else args$cutoff
  .check_numeric_vector(bounds, "control$phi_bounds", n = 2L, nonnegative = TRUE)
  if (bounds[2L] <= bounds[1L]) {
    stop("control$phi_bounds must be increasing.", call. = FALSE)
  }
  low <- empirical_null_ind_RQ(
    args$z, args$size, p_grid = p_grid,
    phi_grid = seq(bounds[1L], bounds[2L], length.out = 251L),
    cutoff = cutoff, common_intercept = TRUE
  )
  means <- rep(low$mu, length(args$z))
  sds <- sqrt(1 + args$size * low$phi)
  standardized <- (args$z - means) / sds
  diagnostics <- .agent_diagnostics(
    standardized, args$size, design = NULL, family = "gaussian",
    model = "individualized", low = low$diagnostics, p = 0L,
    second_moment = NULL, cutoff = cutoff
  )
  structure(list(
    method = "individualized", family = "gaussian",
    parameters = list(phi = low$phi, p0 = low$p0, mu = low$mu),
    null_mean = means, null_sd = sds, standardized = standardized,
    diagnostics = diagnostics,
    status = if (length(diagnostics$warning_codes)) "warning" else "ok",
    alpha = 0.05, cutoff = cutoff, design = NULL, size = args$size,
    second_moment = NULL, z = args$z, low_level_fit = low,
    call = match.call()
  ), class = "empinull_fit")
}
