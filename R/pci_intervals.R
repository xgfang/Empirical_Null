# Conditional intervals: fitted null parameters are held fixed while the
# original measure/reference is varied. This is not a bootstrap of null fitting.
.pci_intervals <- function(stat, statistic, calibration, mu, sd, alpha,
                           alternative, ci, zero_method, mid_p,
                           probabilities = NULL, control = list()) {
  n <- length(stat$z)
  lower <- upper <- rep(NA_real_, n)
  scale <- rep(NA_character_, n)
  codes <- character()
  problems <- list()
  poisson_method <- statistic %in% c("poisson_midp", "poisson_exact", "poisson_score")
  binary_method <- statistic %in% c("poibin", "binary_score")
  selected <- ci
  if (selected == "auto") {
    selected <- if (statistic %in% c("wald", "log_ratio")) "wald" else
      if (statistic == "poisson_exact" && calibration == "theoretical") "poisson" else
        if (poisson_method || binary_method) "invert" else "none"
  }
  result <- function() list(lower = lower, upper = upper, scale = scale,
                           method = selected, warning_codes = unique(codes),
                           diagnostics = list(problems = problems,
                                              null_parameters_held_fixed = TRUE))
  if (selected == "none") return(result())
  if (selected == "wald" && !statistic %in% c("wald", "log_ratio")) {
    stop("Wald intervals require statistic='wald' or 'log_ratio' and estimate/SE inputs.", call. = FALSE)
  }
  if (selected == "poisson" && (!poisson_method || calibration != "theoretical")) {
    stop("Classical Poisson intervals require a Poisson statistic and theoretical calibration; use ci='invert' for calibrated intervals.", call. = FALSE)
  }
  if (selected == "invert" && !poisson_method && !binary_method) {
    stop("Test inversion supports Poisson, poibin, or binary_score inputs; simulation and score-only inputs have no automatic original-scale intervals.", call. = FALSE)
  }
  tail <- if (alternative == "two.sided") alpha / 2 else alpha
  critical <- stats::qnorm(tail, lower.tail = FALSE)
  low_z <- if (alternative == "greater") -Inf else -critical
  high_z <- if (alternative == "less") Inf else critical
  raw_low_z <- mu + sd * low_z
  raw_high_z <- mu + sd * high_z
  if (selected == "wald") {
    # A finite-df Wald statistic enters calibration as qnorm(pt(t)). Invert
    # that transform rather than scaling a t critical value approximately.
    to_t <- function(x, df) {
      ans <- x
      finite <- is.finite(df) & !is.na(x)
      pos <- finite & x >= 0
      neg <- finite & x < 0
      ans[pos] <- stats::qt(stats::pnorm(x[pos], lower.tail = FALSE, log.p = TRUE),
                            df[pos], lower.tail = FALSE, log.p = TRUE)
      ans[neg] <- stats::qt(stats::pnorm(x[neg], log.p = TRUE), df[neg], log.p = TRUE)
      ans
    }
    lo_t <- to_t(raw_low_z, stat$df)
    hi_t <- to_t(raw_high_z, stat$df)
    working_lower <- stat$working_estimate - stat$working_se * hi_t
    working_upper <- stat$working_estimate - stat$working_se * lo_t
    inverse <- switch(stat$transform, identity = identity, log = exp, logit = stats::plogis)
    lower <- inverse(working_lower)
    upper <- inverse(working_upper)
    scale[] <- if (stat$transform == "log") "ratio" else
      if (stat$transform == "logit") "rate" else "estimate"
    zero <- statistic == "log_ratio" & !is.na(stat$estimate) & stat$estimate == 0
    if (any(zero)) {
      lower[zero] <- upper[zero] <- NA_real_
      if (zero_method == "score") {
        lower[zero] <- 0
        upper[zero] <- if (alternative == "greater") Inf else critical * sd[zero] * stat$se[zero]
        codes <- c(codes, "ZERO_RATIO_SCORE_APPROXIMATION")
        problems$zero_ratio <- "The zero-ratio interval is the legacy linear score approximation, not an inversion of the log-ratio test."
      }
    }
    undefined <- is.na(stat$z)
    lower[undefined] <- upper[undefined] <- NA_real_
    return(result())
  }
  threshold <- .pci_ci_positive(control$normal_threshold, 100, "normal_threshold", allow_inf = TRUE)
  upper_cap <- .pci_ci_positive(control$upper_cap, 1e8, "upper_cap")
  tolerance <- .pci_ci_positive(control$root_tol, 1e-8, "root_tol")
  log_cap <- .pci_ci_positive(control$max_log_odds, 40, "max_log_odds")
  p_floor <- .pci_stat_floor(if (is.null(control$p_floor)) 1e-6 else control$p_floor)
  # Raw theoretical tails are not floored. Use a smaller computational floor
  # for their inversion so a stringent requested alpha does not inherit the
  # normal-score floor used only for empirical-null calibration.
  curve_floor <- if (calibration == "theoretical")
    max(.Machine$double.xmin * .Machine$double.eps, min(p_floor, alpha / 100)) else p_floor
  if (selected == "poisson") {
    if (statistic != "poisson_exact") codes <- c(codes, "INTERVAL_METHOD_DIFFERS_FROM_TEST")
    o <- stat$observed
    e <- stat$expected
    valid <- e > 0 & is.finite(e)
    lower <- rep(0, n)
    upper <- rep(Inf, n)
    exact <- e < threshold
    lower[exact & o > 0] <- stats::qchisq(tail, 2 * o[exact & o > 0]) / 2
    upper[exact] <- stats::qchisq(tail, 2 * (o[exact] + 1), lower.tail = FALSE) / 2
    byar <- !exact
    lower[byar & o > 0] <- o[byar & o > 0] *
      pmax(0, 1 - 1 / (9 * o[byar & o > 0]) - critical / (3 * sqrt(o[byar & o > 0])))^3
    upper[byar] <- (o[byar] + 1) *
      (1 - 1 / (9 * (o[byar] + 1)) + critical / (3 * sqrt(o[byar] + 1)))^3
    if (alternative == "greater") upper[] <- Inf
    if (alternative == "less") lower[] <- 0
    lower <- lower / e
    upper <- upper / e
    lower[!valid] <- upper[!valid] <- NA_real_
    scale[] <- "ratio"
    if (any(!valid)) codes <- c(codes, "ZERO_EXPECTED_RATIO_UNDEFINED")
    return(result())
  }
  if (is.null(probabilities)) probabilities <- stat$metadata$probabilities
  score_cap <- stats::qnorm(p_floor, lower.tail = FALSE)
  for (i in seq_len(n)) {
    if (is.na(stat$z[i])) next
    if (poisson_method) scale[i] <- "ratio" else scale[i] <- "log_odds_shift"
    if (poisson_method && stat$expected[i] <= 0) {
      codes <- c(codes, "ZERO_EXPECTED_RATIO_UNDEFINED")
      next
    }
    # Floors make extreme normal scores intentionally noninvertible beyond
    # their attainable range. Report that limitation rather than fabricating
    # a finite confidence bound.
    discrete_quantile <- statistic %in% c("poisson_midp", "poisson_exact", "poibin")
    needed <- c(raw_low_z[i], raw_high_z[i])
    if (discrete_quantile && calibration != "theoretical" &&
        any(abs(needed[is.finite(needed)]) >= score_cap)) {
      codes <- c(codes, "SCORE_FLOOR_LIMITS_INVERSION")
      problems[[paste0("provider_", i)]] <- "A requested calibrated cutoff exceeds the attainable floored normal-score range."
      next
    }
    bounds <- tryCatch({
      if (poisson_method) {
        poisson_curve <- function(mean) {
          if (statistic == "poisson_score") {
            if (mean == 0) return(if (stat$observed[i] == 0) 0 else Inf)
            return((stat$observed[i] - mean) / sqrt(mean))
          }
          .pci_poisson_z(stat$observed[i], mean,
                         mid_p = statistic == "poisson_midp", p_floor = curve_floor)
        }
        if (poisson_curve(0) < raw_low_z[i]) {
          .pci_empty_interval("The largest attainable Poisson score is below the calibrated acceptance interval.")
        }
        lo <- if (alternative == "less") 0 else
          .pci_count_root(poisson_curve, raw_high_z[i], stat$observed[i], stat$expected[i], upper_cap, tolerance)
        hi <- if (alternative == "greater") Inf else
          .pci_count_root(poisson_curve, raw_low_z[i], stat$observed[i], stat$expected[i], upper_cap, tolerance)
        c(lo, hi) / stat$expected[i]
      } else {
        probs <- probabilities[[i]]
        if (statistic == "binary_score" &&
            .pci_binary_score_nonmonotone(stat$observed[i], probs, log_cap)) {
          stop(structure(list(message = paste(
            "The heterogeneous binary score is nonmonotone; its inverted acceptance set may be disconnected.",
            "A single confidence interval is not reported. Use the Poisson-binomial test for monotone tail inversion."
          ), call = NULL), class = c("pci_nonmonotone_interval", "error", "condition")))
        }
        # Under theoretical calibration use the chosen one-sided binary tail,
        # which may be full-tail even when two-sided tests use mid-p.
        binary_midp <- mid_p
        if (calibration == "theoretical" && alternative != "two.sided") {
          binary_midp <- isTRUE(control$one_sided_midp)
        }
        binary_curve <- function(delta) {
          shifted <- if (delta == -Inf) as.numeric(probs == 1) else
            if (delta == Inf) as.numeric(probs > 0) else
              stats::plogis(stats::qlogis(probs) + delta)
          if (statistic == "binary_score") {
            v <- sum(shifted * (1 - shifted))
            d <- stat$observed[i] - sum(shifted)
            if (v == 0) return(if (d == 0) 0 else sign(d) * Inf)
            return(d / sqrt(v))
          }
          if (calibration == "theoretical" && alternative != "two.sided") {
            tails <- .pci_poibin_tails(stat$observed[i], shifted, binary_midp)
            # qnorm(F_less) and -qnorm(F_greater) agree for continuous
            # distributions; full discrete tails deliberately differ.
            if (alternative == "less") return(stats::qnorm(pmax(curve_floor, pmin(1 - curve_floor, tails$lower))))
            return(-stats::qnorm(pmax(curve_floor, pmin(1 - curve_floor, tails$upper))))
          }
          .pci_binary_z(stat$observed[i], shifted, mid_p, curve_floor)
        }
        if (binary_curve(-Inf) < raw_low_z[i] || binary_curve(Inf) > raw_high_z[i]) {
          .pci_empty_interval("No attainable binary score lies in the calibrated acceptance interval.")
        }
        lo <- if (alternative == "less") -Inf else
          .pci_binary_root(binary_curve, raw_high_z[i], log_cap, tolerance, "lower")
        hi <- if (alternative == "greater") Inf else
          .pci_binary_root(binary_curve, raw_low_z[i], log_cap, tolerance, "upper")
        c(lo, hi)
      }
    }, error = function(e) e)
    if (inherits(bounds, "error")) {
      codes <- c(codes, if (inherits(bounds, "pci_empty_interval"))
        "EMPTY_INVERTED_INTERVAL" else if (inherits(bounds, "pci_nonmonotone_interval"))
          "NONMONOTONE_SCORE_INVERSION" else "INTERVAL_INVERSION_FAILED")
      problems[[paste0("provider_", i)]] <- conditionMessage(bounds)
    } else if (bounds[1] > bounds[2]) {
      codes <- c(codes, "EMPTY_INVERTED_INTERVAL")
      problems[[paste0("provider_", i)]] <- "No interval satisfies the requested calibrated test cutoffs."
    } else {
      lower[i] <- bounds[1]
      upper[i] <- bounds[2]
    }
  }
  result()
}

.pci_ci_positive <- function(value, default, name, allow_inf = FALSE) {
  if (is.null(value)) return(default)
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0 ||
      (!allow_inf && !is.finite(value))) {
    stop("control$", name, " must be a positive number.", call. = FALSE)
  }
  value
}

.pci_count_root <- function(curve, target, observed, expected, cap, tolerance) {
  if (target == Inf) return(0)
  if (target == -Inf) return(Inf)
  at_zero <- curve(0) - target
  if (is.na(at_zero)) stop("Undefined score at the lower count boundary.")
  if (at_zero <= 0) return(0)
  upper <- min(cap, max(1, observed + 1, expected))
  while (curve(upper) > target && upper < cap) upper <- min(cap, upper * 2)
  if (curve(upper) > target) stop("Poisson root exceeds control$upper_cap.")
  stats::uniroot(function(x) curve(x) - target, c(0, upper),
                 tol = tolerance, check.conv = TRUE)$root
}

.pci_binary_root <- function(curve, target, cap, tolerance, side) {
  if (target == Inf) return(-Inf)
  if (target == -Inf) return(Inf)
  # Only asymptotic limits justify infinite endpoints. A finite search cap
  # failing to bracket a root is a numerical limitation, not unbounded support.
  left_limit <- curve(-Inf) - target
  right_limit <- curve(Inf) - target
  if (left_limit <= 0) {
    if (side == "lower") return(-Inf)
    .pci_empty_interval("No finite binary upper endpoint satisfies the calibrated test.")
  }
  if (right_limit >= 0) {
    if (side == "upper") return(Inf)
    .pci_empty_interval("No finite binary lower endpoint satisfies the calibrated test.")
  }
  a <- curve(-cap) - target
  b <- curve(cap) - target
  if (is.na(a) || is.na(b)) stop("Undefined binary score at an inversion boundary.")
  if (a <= 0) {
    stop("Binary root lies below the control$max_log_odds search range; increase that limit.")
  }
  if (b >= 0) {
    stop("Binary root lies above the control$max_log_odds search range; increase that limit.")
  }
  stats::uniroot(function(delta) curve(delta) - target,
                 c(-cap, cap), tol = tolerance, check.conv = TRUE)$root
}

.pci_empty_interval <- function(message) {
  stop(structure(list(message = message, call = NULL),
                 class = c("pci_empty_interval", "error", "condition")))
}

.pci_binary_score_nonmonotone <- function(observed, probs, cap) {
  # The derivative of (O-sum(p))/sqrt(sum(p*(1-p))) need not be negative
  # when Bernoulli probabilities differ. Check its analytic sign on a dense
  # grid augmented at the probability transitions before treating the score
  # as an invertible pivot. This is a numerical guard, not a general proof of
  # monotonicity for arbitrary probability vectors.
  finite_logits <- stats::qlogis(probs[probs > 0 & probs < 1])
  if (!length(finite_logits)) return(FALSE)
  transitions <- -as.numeric(stats::quantile(finite_logits,
    probs = seq(0, 1, length.out = min(41L, length(finite_logits))), names = FALSE))
  points <- sort(unique(c(seq(-cap, cap, length.out = max(161L, ceiling(8 * cap) + 1L)),
                          transitions[abs(transitions) <= cap])))
  for (delta in points) {
    p <- stats::plogis(stats::qlogis(probs) + delta)
    information <- sum(p * (1 - p))
    derivative_information <- sum(p * (1 - p) * (1 - 2 * p))
    if (information <= .Machine$double.eps) next
    numerator <- -information^2 - (observed - sum(p)) * derivative_information / 2
    if (numerator > 1e-10 * information^2) return(TRUE)
  }
  FALSE
}
