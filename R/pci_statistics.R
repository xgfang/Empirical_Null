# Raw, provider-level statistics used by fit_pci(). These helpers deliberately
# do not fit a patient-level model or estimate an empirical null.

.pci_stat_vector <- function(x, n, name, default = NULL,
                             infinite = FALSE, missing = FALSE) {
  if (is.null(x)) x <- default
  if (is.null(x) || !is.numeric(x) || is.matrix(x) ||
      !length(x) || !(length(x) %in% c(1L, n))) {
    stop(name, " must be numeric, with length 1 or one value per provider.",
         call. = FALSE)
  }
  if (length(x) == 1L) x <- rep(x, n)
  bad <- if (infinite) is.nan(x) else !is.finite(x)
  if (missing) bad[is.na(x)] <- FALSE
  if (any(bad) || (!missing && anyNA(x))) {
    stop(name, " contains invalid or missing values.", call. = FALSE)
  }
  as.numeric(x)
}

.pci_stat_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.pci_stat_floor <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || x >= 0.5) {
    stop("control$p_floor must be strictly between 0 and 0.5.", call. = FALSE)
  }
  x
}

.pci_stat_alternative <- function(alternative) {
  if (identical(alternative, "two_sided")) alternative <- "two.sided"
  match.arg(alternative, c("two.sided", "less", "greater"))
}

.pci_stat_normal_p <- function(z, alternative) {
  switch(alternative,
         two.sided = 2 * stats::pnorm(-abs(z)),
         less = stats::pnorm(z),
         greater = stats::pnorm(z, lower.tail = FALSE))
}

.pci_stat_tail_p <- function(lower, upper, alternative) {
  switch(alternative, two.sided = pmin(1, 2 * pmin(lower, upper)),
         less = lower, greater = upper)
}

.pci_stat_tail_z <- function(lower, upper, p_floor) {
  smaller <- pmin(0.5, pmax(p_floor, pmin(lower, upper)))
  magnitude <- stats::qnorm(smaller, lower.tail = FALSE)
  ifelse(lower <= upper, -magnitude, magnitude)
}

.pci_stat_count <- function(x, name, upper = Inf) {
  if (any(!is.finite(x)) || any(x < 0) ||
      any(abs(x - round(x)) > 1e-8) || any(x > upper)) {
    stop(name, " must contain nonnegative integer counts within their support.",
         call. = FALSE)
  }
  round(x)
}

.pci_stat_list <- function(x, n, name, default = NULL) {
  if (is.null(x)) x <- default
  if (is.numeric(x) && !is.matrix(x) && n == 1L) x <- list(x)
  if (!is.list(x) || length(x) != n) {
    stop(name, " must be a list with one numeric vector per provider",
         " (a numeric vector is allowed for one provider).", call. = FALSE)
  }
  for (i in seq_len(n)) {
    if (!is.numeric(x[[i]]) || is.matrix(x[[i]]) ||
        !length(x[[i]]) || any(!is.finite(x[[i]]))) {
      stop(name, "[[", i, "]] must be a nonempty finite numeric vector.",
           call. = FALSE)
    }
    x[[i]] <- as.numeric(x[[i]])
  }
  x
}

.pci_poisson_tails <- function(observed, expected, mid_p = TRUE) {
  mass <- stats::dpois(observed, expected)
  weight <- if (mid_p) 0.5 else 1
  list(lower = pmin(1, stats::ppois(observed - 1, expected) + weight * mass),
       upper = pmin(1, stats::ppois(observed, expected, lower.tail = FALSE) +
                      weight * mass),
       mass = mass)
}

# Return the signed normal quantile of a Poisson tail. The floor is applied
# to a ONE-sided tail, as in Python's poisson_midp_zscore(), not to the final
# two-sided p-value. A deterministic zero count under a zero mean has z = 0.
.pci_poisson_z <- function(observed, expected, mid_p = TRUE, p_floor = 1e-6) {
  n <- max(length(observed), length(expected))
  observed <- .pci_stat_vector(observed, n, "observed")
  expected <- .pci_stat_vector(expected, n, "expected")
  observed <- .pci_stat_count(observed, "observed")
  if (any(expected < 0)) stop("expected must be nonnegative.", call. = FALSE)
  mid_p <- .pci_stat_flag(mid_p, "mid_p")
  p_floor <- .pci_stat_floor(p_floor)
  tails <- .pci_poisson_tails(observed, expected, mid_p)
  if (mid_p) {
    z <- .pci_stat_tail_z(tails$lower, tails$upper, p_floor)
  } else {
    # Python's conventional Poisson test selects the tail by O versus E.
    side <- ifelse(observed > expected, tails$upper, tails$lower)
    p <- pmin(0.999, 2 * side)
    z <- ifelse(observed > expected, 1, -1) *
      stats::qnorm(pmax(p_floor, p / 2), lower.tail = FALSE)
  }
  z[expected == 0 & observed == 0] <- 0
  z
}

# Dynamic programming distribution of a sum of independent Bernoulli draws.
# Unlike probability clipping, retaining p = 0 and p = 1 preserves genuinely
# deterministic observations and impossible outcomes.
.pci_poibin_tails <- function(observed, probs, mid_p = TRUE) {
  probs <- .pci_stat_list(probs, 1L, "probs")[[1L]]
  if (any(probs < 0 | probs > 1)) {
    stop("probs must lie in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(observed) || length(observed) != 1L) {
    stop("observed must be one count.", call. = FALSE)
  }
  observed <- .pci_stat_count(observed, "observed", length(probs))
  mid_p <- .pci_stat_flag(mid_p, "mid_p")
  pmf <- 1
  for (prob in probs) {
    pmf <- c(pmf * (1 - prob), 0) + c(0, pmf * prob)
  }
  pmf <- pmf / sum(pmf)
  index <- observed + 1L
  mass <- pmf[index]
  weight <- if (mid_p) 0.5 else 1
  lower <- if (observed == 0) 0 else sum(pmf[seq_len(observed)])
  upper <- if (index == length(pmf)) 0 else
    sum(pmf[seq.int(index + 1L, length(pmf))])
  lower <- min(1, lower + weight * mass)
  upper <- min(1, upper + weight * mass)
  list(lower = lower, upper = upper,
       two.sided = min(1, 2 * min(lower, upper)),
       two_sided = min(1, 2 * min(lower, upper)),
       mass = mass, pmf = pmf)
}

.pci_binary_z <- function(observed, probs, mid_p = TRUE, p_floor = 1e-6) {
  tails <- .pci_poibin_tails(observed, probs, mid_p)
  .pci_stat_tail_z(tails$lower, tails$upper, .pci_stat_floor(p_floor))
}

.pci_stat_with_seed <- function(seed, code) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != round(seed)) {
    stop("control$seed must be a nonnegative integer seed.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed), kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
  force(code)
}

.pci_stat_sim_tails <- function(observed, simulated, mid_p) {
  mass <- mean(simulated == observed)
  weight <- if (mid_p) 0.5 else 1
  list(lower = mean(simulated < observed) + weight * mass,
       upper = mean(simulated > observed) + weight * mass,
       mass = mass)
}

.pci_stat_simulate <- function(probabilities, eta, re_mean, re_var,
                               null_effect, n_resample) {
  n_obs <- if (is.null(probabilities)) length(eta) else length(probabilities)
  result <- numeric(n_resample)
  # Bounded batches avoid allocating n_observations * n_resample at once.
  batch_size <- max(1L, min(n_resample, floor(1e6 / n_obs)))
  for (start in seq.int(1L, n_resample, by = batch_size)) {
    end <- min(n_resample, start + batch_size - 1L)
    count <- end - start + 1L
    if (is.null(probabilities)) {
      random_effect <- stats::rnorm(n_obs * count,
                                    mean = rep(re_mean, count),
                                    sd = rep(sqrt(re_var), count))
      probability <- stats::plogis(null_effect + rep(eta, count) + random_effect)
    } else {
      probability <- rep(probabilities, count)
    }
    draws <- matrix(stats::rbinom(n_obs * count, 1, probability), nrow = n_obs)
    result[seq.int(start, end)] <- colSums(draws)
  }
  result
}

.pci_compute_statistics <- function(n, statistic, z = NULL, estimate = NULL,
                                    se = NULL, observed = NULL, expected = NULL,
                                    p_value = NULL, direction = NULL,
                                    reference = NULL, transform = "identity",
                                    se_scale = "original", df = Inf,
                                    probabilities = NULL, eta = NULL,
                                    re_mean = NULL, re_var = NULL, mid_p = TRUE,
                                    zero_method = "score",
                                    alternative = "two.sided", control = list()) {
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
      n < 1 || n != round(n)) stop("n must be a positive integer.", call. = FALSE)
  n <- as.integer(n)
  if (identical(statistic, "poibin")) statistic <- "poibin_exact"
  statistic <- match.arg(statistic, c("z", "p_value", "wald", "log_ratio",
                                     "poisson_midp", "poisson_exact",
                                     "poisson_score", "poibin_exact",
                                     "binary_score", "resampling"))
  alternative <- .pci_stat_alternative(alternative)
  if (!is.list(control)) stop("control must be a list.", call. = FALSE)
  p_floor <- .pci_stat_floor(if (is.null(control$p_floor)) 1e-6 else control$p_floor)
  mid_p <- .pci_stat_flag(mid_p, "mid_p")
  one_sided_midp <- .pci_stat_flag(if (is.null(control$one_sided_midp)) FALSE else
                                    control$one_sided_midp, "control$one_sided_midp")
  transform <- match.arg(transform, c("identity", "log", "logit"))
  if (statistic == "log_ratio" && missing(se_scale)) se_scale <- "transformed"
  if (se_scale == "working") se_scale <- "transformed"
  se_scale <- match.arg(se_scale, c("original", "transformed"))
  zero_method <- match.arg(zero_method, c("score", "exclude"))
  df <- .pci_stat_vector(df, n, "df", infinite = TRUE)
  if (any(df <= 0)) stop("df must be positive or Inf.", call. = FALSE)
  blank <- rep(NA_real_, n)
  out <- list(z = blank, p = blank, estimate = blank, se = blank,
              reference = blank, working_estimate = blank, working_se = blank,
              working_reference = blank, df = df, transform = transform,
              observed = NULL, expected = NULL, statistic = statistic,
              test_statistic = blank,
              metadata = list(alternative = alternative, mid_p = mid_p,
                              p_floor = p_floor, se_scale = se_scale))

  if (statistic == "z") {
    out$z <- .pci_stat_vector(z, n, "z")
    out$p <- .pci_stat_normal_p(out$z, alternative)
    out$test_statistic <- out$z
    out$estimate <- out$working_estimate <- out$z
    out$se <- out$working_se <- rep(1, n)
    out$reference <- out$working_reference <- rep(0, n)
    return(out)
  }

  if (statistic == "p_value") {
    p_value <- .pci_stat_vector(p_value, n, "p_value")
    direction <- .pci_stat_vector(direction, n, "direction")
    if (any(p_value < 0 | p_value > 1)) {
      stop("p_value must lie in [0, 1].", call. = FALSE)
    }
    if (any(direction == 0 & p_value < 1)) {
      stop("direction must be nonzero when p_value is below 1.", call. = FALSE)
    }
    quantile_tail <- pmax(.Machine$double.xmin, p_value / 2)
    quantile_tail[p_value == 0] <- p_floor
    out$z <- sign(direction) * stats::qnorm(quantile_tail, lower.tail = FALSE)
    lower <- ifelse(direction < 0, p_value / 2, 1 - p_value / 2)
    upper <- ifelse(direction > 0, p_value / 2, 1 - p_value / 2)
    out$p <- if (alternative == "two.sided") p_value else
      if (alternative == "less") lower else upper
    out$test_statistic <- out$z
    out$metadata$input_p_value <- p_value
    out$metadata$direction <- sign(direction)
    return(out)
  }

  if (statistic %in% c("wald", "log_ratio")) {
    estimate <- .pci_stat_vector(estimate, n, "estimate")
    se <- .pci_stat_vector(se, n, "se")
    if (any(se <= 0)) stop("se must be strictly positive.", call. = FALSE)
    if (statistic == "log_ratio") {
      transform <- "log"
    }
    default_reference <- switch(transform, identity = 0, log = 1, logit = 0.5)
    reference <- .pci_stat_vector(reference, n, "reference", default_reference)
    if (transform == "identity") {
      working_estimate <- estimate
      working_reference <- reference
      working_se <- se
    } else if (transform == "log") {
      allow_zero <- statistic == "log_ratio"
      if (any(estimate < 0 | (!allow_zero & estimate == 0)) ||
          any(reference <= 0)) {
        stop("Log transformation requires positive estimates and references; ",
             "only log_ratio supports zero estimates.", call. = FALSE)
      }
      working_estimate <- log(estimate)
      working_reference <- log(reference)
      working_se <- if (se_scale == "original") se / estimate else se
    } else {
      if (any(estimate <= 0 | estimate >= 1) ||
          any(reference <= 0 | reference >= 1)) {
        stop("Logit transformation requires estimates and references in (0, 1).",
             call. = FALSE)
      }
      working_estimate <- stats::qlogis(estimate)
      working_reference <- stats::qlogis(reference)
      working_se <- if (se_scale == "original") se / (estimate * (1 - estimate)) else se
    }
    test_statistic <- (working_estimate - working_reference) / working_se
    if (statistic == "log_ratio" && any(estimate == 0)) {
      zeros <- estimate == 0
      if (se_scale != "transformed") {
        stop("Zero-ratio handling requires a transformed (log-scale) se.",
             call. = FALSE)
      }
      if (any(reference[zeros] != 1)) {
        stop("The zero-ratio score approximation requires reference = 1.",
             call. = FALSE)
      }
      test_statistic[zeros] <- if (zero_method == "score") -1 / se[zeros] else NA_real_
    }
    lower <- stats::pt(test_statistic, df = df)
    upper <- stats::pt(test_statistic, df = df, lower.tail = FALSE)
    out$p <- .pci_stat_tail_p(lower, upper, alternative)
    # Finite-df statistics must be calibrated on a normal-quantile scale;
    # normal Wald scores themselves need no artificial clipping.
    log_tail <- stats::pt(abs(test_statistic), df = df,
                          lower.tail = FALSE, log.p = TRUE)
    log_tail[is.infinite(log_tail) & log_tail < 0] <- log(p_floor)
    converted_z <- sign(test_statistic) *
      stats::qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)
    converted_z[is.infinite(df)] <- test_statistic[is.infinite(df)]
    out$z <- converted_z
    out$test_statistic <- test_statistic
    out$estimate <- estimate
    out$se <- se
    out$reference <- reference
    out$working_estimate <- working_estimate
    out$working_se <- working_se
    out$working_reference <- working_reference
    out$transform <- transform
    out$metadata$se_scale <- se_scale
    out$metadata$zero_method <- zero_method
    return(out)
  }

  if (statistic %in% c("poisson_midp", "poisson_exact", "poisson_score")) {
    observed <- .pci_stat_vector(observed, n, "observed")
    observed <- .pci_stat_count(observed, "observed")
    expected <- .pci_stat_vector(expected, n, "expected")
    if (any(expected < 0)) stop("expected must be nonnegative.", call. = FALSE)
    reference <- .pci_stat_vector(reference, n, "reference", 1)
    if (any(reference < 0)) stop("A Poisson ratio reference must be nonnegative.", call. = FALSE)
    null_expected <- expected * reference
    if (any(!is.finite(null_expected))) {
      stop("expected * reference must remain finite.", call. = FALSE)
    }
    deterministic <- null_expected == 0
    at_null <- deterministic & observed == 0
    pearson <- (observed - null_expected) / sqrt(null_expected)
    pearson[at_null] <- 0
    if (statistic == "poisson_score") {
      out$z <- pearson
      out$z[is.infinite(out$z)] <- stats::qnorm(p_floor, lower.tail = FALSE)
      out$p <- .pci_stat_normal_p(pearson, alternative)
    } else {
      use_midp <- statistic == "poisson_midp"
      # The explicitly named Poisson mid-p statistic uses half-mass tails
      # for every alternative. one_sided_midp controls the binary methods.
      tails <- .pci_poisson_tails(observed, null_expected, use_midp)
      if (statistic == "poisson_exact" && alternative == "two.sided") {
        out$p <- pmin(0.999, 2 * ifelse(observed > null_expected, tails$upper, tails$lower))
      } else {
        out$p <- .pci_stat_tail_p(tails$lower, tails$upper, alternative)
      }
      out$z <- .pci_poisson_z(observed, null_expected, use_midp, p_floor)
      out$metadata$mid_p <- use_midp
      out$metadata$one_sided_midp <- use_midp
    }
    out$p[at_null] <- if (alternative == "two.sided" ||
                           statistic != "poisson_midp") 1 else 0.5
    out$estimate <- observed / expected
    out$estimate[expected == 0 & observed == 0] <- NA_real_
    out$reference <- reference
    out$test_statistic <- if (statistic == "poisson_score") pearson else out$z
    out$observed <- observed
    out$expected <- expected
    out$metadata$null_expected <- null_expected
    return(out)
  }

  # Binary tests either take already-computed null probabilities, or a linear
  # predictor excluding the target effect together with posterior RE moments.
  observed <- .pci_stat_vector(observed, n, "observed")
  if (!is.null(probabilities) && !is.null(eta)) {
    stop("Supply probabilities or eta, not both.", call. = FALSE)
  }
  null_effect <- .pci_stat_vector(reference, n, "reference", 0)
  eta_input <- !is.null(eta)
  if (eta_input) {
    eta <- .pci_stat_list(eta, n, "eta")
    zero_list <- lapply(eta, function(x) rep(0, length(x)))
    re_mean <- .pci_stat_list(re_mean, n, "re_mean", zero_list)
    re_var <- .pci_stat_list(re_var, n, "re_var", zero_list)
    for (i in seq_len(n)) {
      if (length(re_mean[[i]]) != length(eta[[i]]) ||
          length(re_var[[i]]) != length(eta[[i]])) {
        stop("eta, re_mean, and re_var must have matching within-provider lengths.",
             call. = FALSE)
      }
      if (any(re_var[[i]] < 0)) stop("re_var must be nonnegative.", call. = FALSE)
    }
    probabilities <- lapply(seq_len(n), function(i)
      stats::plogis(null_effect[i] + eta[[i]] + re_mean[[i]]))
  } else {
    probabilities <- .pci_stat_list(probabilities, n, "probabilities")
    if (!is.null(re_mean) || !is.null(re_var)) {
      stop("re_mean and re_var require eta.", call. = FALSE)
    }
    if (any(null_effect != 0)) {
      stop("With probabilities, encode the null in those probabilities; ",
           "reference is only a log-odds effect for eta input.", call. = FALSE)
    }
  }
  if (any(vapply(probabilities, function(x) any(x < 0 | x > 1), logical(1)))) {
    stop("probabilities must lie in [0, 1].", call. = FALSE)
  }
  observed <- .pci_stat_count(observed, "observed", lengths(probabilities))
  expected <- vapply(probabilities, sum, numeric(1))
  variance <- vapply(probabilities, function(x) sum(x * (1 - x)), numeric(1))
  pearson <- (observed - expected) / sqrt(variance)
  pearson[variance == 0 & observed == expected] <- 0
  if (statistic == "binary_score") {
    out$z <- pearson
    infinite <- is.infinite(out$z)
    out$z[infinite] <- sign(out$z[infinite]) * stats::qnorm(p_floor, lower.tail = FALSE)
    out$p <- .pci_stat_normal_p(pearson, alternative)
    out$p[variance == 0 & observed == expected] <- 1
    out$test_statistic <- pearson
  } else {
    if (statistic == "resampling") {
      n_resample <- if (is.null(control$n_resample)) 10000 else control$n_resample
      if (!is.numeric(n_resample) || length(n_resample) != 1L ||
          !is.finite(n_resample) || n_resample < 1 ||
          n_resample != round(n_resample) || n_resample > .Machine$integer.max) {
        stop("control$n_resample must be a positive integer.", call. = FALSE)
      }
      seed <- if (is.null(control$seed)) 1 else control$seed
      simulations <- .pci_stat_with_seed(seed, lapply(seq_len(n), function(i)
        .pci_stat_simulate(if (eta_input) NULL else probabilities[[i]],
                          if (eta_input) eta[[i]] else NULL,
                          if (eta_input) re_mean[[i]] else NULL,
                          if (eta_input) re_var[[i]] else NULL,
                          null_effect[i], n_resample)))
      score_tails <- lapply(seq_len(n), function(i)
        .pci_stat_sim_tails(observed[i], simulations[[i]], mid_p))
      test_tails <- if (alternative == "two.sided" || one_sided_midp == mid_p)
        score_tails else lapply(seq_len(n), function(i)
          .pci_stat_sim_tails(observed[i], simulations[[i]], one_sided_midp))
      out$metadata$n_resample <- n_resample
      out$metadata$seed <- seed
      out$metadata$posterior_sampling <- eta_input
      out$metadata$simulation_mean <- vapply(simulations, mean, numeric(1))
    } else {
      score_tails <- lapply(seq_len(n), function(i)
        .pci_poibin_tails(observed[i], probabilities[[i]], mid_p))
      test_tails <- if (alternative == "two.sided" || one_sided_midp == mid_p)
        score_tails else lapply(seq_len(n), function(i)
          .pci_poibin_tails(observed[i], probabilities[[i]], one_sided_midp))
    }
    lower <- vapply(score_tails, `[[`, numeric(1), "lower")
    upper <- vapply(score_tails, `[[`, numeric(1), "upper")
    out$z <- .pci_stat_tail_z(lower, upper, p_floor)
    out$p <- .pci_stat_tail_p(vapply(test_tails, `[[`, numeric(1), "lower"),
                              vapply(test_tails, `[[`, numeric(1), "upper"),
                              alternative)
    out$test_statistic <- out$z
    out$metadata$one_sided_midp <- one_sided_midp
  }
  out$estimate <- observed
  out$se <- sqrt(variance)
  out$reference <- expected
  out$observed <- observed
  out$expected <- expected
  out$metadata$probabilities <- probabilities
  out$metadata$null_effect <- null_effect
  out$metadata$eta <- eta
  out$metadata$re_mean <- re_mean
  out$metadata$re_var <- re_var
  out
}
