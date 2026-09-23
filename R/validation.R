.central_z <- function(coverage) {
  if (length(coverage) != 1L || !is.finite(coverage) ||
      coverage <= 0 || coverage >= 1) {
    stop("coverage must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  stats::qnorm((1 + coverage) / 2)
}

.check_numeric_vector <- function(x, name, n = NULL, positive = FALSE,
                                  nonnegative = FALSE) {
  if (!is.numeric(x) || is.matrix(x) || length(x) == 0L) {
    stop(name, " must be a non-empty numeric vector.", call. = FALSE)
  }
  if (!is.null(n) && length(x) != n) {
    stop(name, " must have length ", n, ".", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (positive && any(x <= 0)) {
    stop(name, " must contain only positive values.", call. = FALSE)
  }
  if (nonnegative && any(x < 0)) {
    stop(name, " must contain only nonnegative values.", call. = FALSE)
  }
  invisible(x)
}

.check_grid <- function(x, name, lower, upper, lower_closed = TRUE,
                        upper_closed = TRUE) {
  .check_numeric_vector(x, name)
  lower_bad <- if (lower_closed) x < lower else x <= lower
  upper_bad <- if (upper_closed) x > upper else x >= upper
  if (any(lower_bad | upper_bad)) {
    left <- if (lower_closed) "[" else "("
    right <- if (upper_closed) "]" else ")"
    stop(name, " must lie in ", left, lower, ", ", upper, right, ".",
         call. = FALSE)
  }
  sort(unique(x))
}

.as_design_matrix <- function(x, n, name = "covariates") {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) != n || ncol(x) < 1L) {
    stop(name, " must have one row per provider and at least one column.",
         call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }
  x
}

.prepare_agent_design <- function(covariates, n, include_intercept = TRUE) {
  if (is.null(covariates)) {
    if (!include_intercept) return(NULL)
    out <- matrix(1, nrow = n, ncol = 1L,
                  dimnames = list(NULL, "(Intercept)"))
    return(out)
  }
  x <- .as_design_matrix(covariates, n)
  constant <- apply(x, 2L, function(v) diff(range(v)) < sqrt(.Machine$double.eps))
  if (include_intercept && !any(constant)) {
    x <- cbind(`(Intercept)` = 1, x)
  }
  x
}

.safe_condition_number <- function(x) {
  if (is.null(x) || ncol(x) == 0L) return(NA_real_)
  if (any(!is.finite(x))) return(Inf)
  if (qr(x)$rank < ncol(x)) return(Inf)
  kappa(x, exact = TRUE)
}

.skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(0)
  mean((x - mean(x))^3) / s^3
}

.boundary_hit <- function(value, grid, tolerance = sqrt(.Machine$double.eps)) {
  if (length(grid) < 2L || !is.finite(value)) return(FALSE)
  abs(value - min(grid)) <= tolerance * max(1, abs(min(grid))) ||
    abs(value - max(grid)) <= tolerance * max(1, abs(max(grid)))
}

.best_optim_index <- function(fits) {
  objectives <- vapply(fits, function(x) x$value, numeric(1L))
  converged <- vapply(
    fits,
    function(x) isTRUE(x$convergence == 0L) && is.finite(x$value),
    logical(1L)
  )
  eligible <- if (any(converged)) which(converged) else which(is.finite(objectives))
  if (!length(eligible)) {
    stop("All candidate optimizations failed to return a finite objective.",
         call. = FALSE)
  }
  eligible[which.min(objectives[eligible])]
}

.make_low_level_diagnostics <- function(convergence, objective, p0, p_grid,
                                        phi = NULL, phi_bounds = NULL,
                                        central_count = NULL) {
  list(
    convergence = as.integer(convergence),
    converged = isTRUE(convergence == 0L) && is.finite(objective),
    objective = as.numeric(objective),
    p0_boundary = .boundary_hit(p0, p_grid),
    phi_boundary = if (is.null(phi) || is.null(phi_bounds)) FALSE else
      phi <= phi_bounds[1L] + 1e-8 || phi >= phi_bounds[2L] - 1e-8,
    central_count = central_count
  )
}
