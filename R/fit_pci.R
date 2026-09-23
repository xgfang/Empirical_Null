#' Provider Tests, Empirical-Null Calibration, and Confidence Intervals
#'
#' Construct a provider statistic, apply a specified null distribution, and
#' return provider-level p-values, directional flags, confidence intervals, and
#' diagnostics through one interface. Supply estimates and their
#' standard errors, observed/expected counts, null probabilities, or scores
#' computed from a fitted model. The function does not fit the outcome model.
#'
#' @param z Precomputed provider-level normal scores.
#' @param estimate,se Estimates and standard errors. For `statistic="log_ratio"`,
#'   `se` is on the log-ratio scale by default.
#' @param observed,expected Observed integer counts and nonnegative expected
#'   counts. Zero expected counts do not define an observed/expected ratio.
#' @param p_value,direction Supplied two-sided p-values and signed directions
#'   (negative below the reference, positive above).
#' @param statistic One of `z`, `wald`, `log_ratio`, `poisson_midp`,
#'   `poisson_exact`, `poisson_score`, `poibin`, `binary_score`, `resampling`,
#'   or `p_value`. `poibin_exact` is an alias for `poibin`.
#' @param calibration One of `theoretical`, `fixed`, `overall`, `groupwise`,
#'   `individualized`, `correlated`, `cre`, or `moment`.
#' @param reference Reference on the original estimate scale. The defaults are
#'   zero for identity, one for log, and 0.5 for logit. A Poisson reference is
#'   a nonnegative ratio multiplying `expected`, default one. For binary
#'   `eta` input this is an additive log-odds effect, default zero; when
#'   supplying `probabilities`, encode the null in those probabilities.
#' @param transform Working transformation for Wald inference: `identity`,
#'   `log`, or `logit`. Numeric reference values are transformed automatically.
#' @param se_scale Whether a Wald SE is on the `original` or `transformed` scale.
#' @param df Wald reference degrees of freedom. Infinite means normal; finite
#'   t statistics are converted to normal quantiles before calibration.
#' @param size Positive provider sizes for grouping or structured nulls. For
#'   Poisson moment calibration, this is the expected-count vector and defaults
#'   to `expected` when `statistic="poisson_score"`.
#' @param covariates Provider-level numeric covariates, with one row per
#'   provider; used by `correlated` and `moment` calibration.
#' @param second_moment For Poisson moment calibration, the sum of squared
#'   observation-level expected counts within each provider. This is not the
#'   square of the provider's total expected count. Requires Pearson scores.
#' @param null_mean,null_sd Supplied mean and SD for `calibration="fixed"`.
#'   Each may be a scalar or one value per provider. For example,
#'   `null_sd=1.81` specifies the SD of the score distribution; it does not
#'   set the outcome reference or supply a variance component.
#' @param psi Robust influence function: `bisquare` or `huber`.
#' @param groups Optional explicit group labels, one per provider.
#' @param n_groups Requested number of size groups for `groupwise` calibration;
#'   one places all providers in a single group. Ignored with explicit `groups`.
#' @param grouping `quantile` uses right-closed size quantiles; `rank` uses
#'   approximately equal-size blocks after stable sorting by size.
#' @param common_mean NULL keeps robust fitted means; TRUE substitutes the
#'   ordinary mean of finite scores among providers selected for fitting; a
#'   number supplies a common mean. Group scales are not refitted under this
#'   override. Applies only to robust `overall` or `groupwise` calibration.
#' @param fit_mask Logical selection of providers used for robust null fitting.
#'   Excluded providers remain in the result. Structured likelihood fits
#'   currently require all retained providers.
#' @param probabilities A list of null Bernoulli probabilities per provider;
#'   a numeric vector is accepted for a single provider.
#' @param eta,re_mean,re_var Lists of observation-level predictors and posterior
#'   random-effect means/variances. `eta` excludes the target effect;
#'   `re_mean` and `re_var` default to zero. Deterministic binary tests use
#'   `plogis(reference + eta + re_mean)` and require zero `re_var`.
#'   Resampling draws independent normal effects per observation with the
#'   supplied means and variances, then simulates Bernoulli outcomes.
#' @param alpha Significance level; also defines the confidence level `1-alpha`.
#' @param alternative `two.sided`, `greater`, or `less`; `two_sided` is an alias.
#' @param ci `auto`, `none`, `wald`, `poisson`, or `invert`. Auto supplies Wald
#'   intervals for estimate/SE input, classical Poisson intervals for an
#'   uncalibrated exact Poisson test, and test-inversion intervals for other
#'   Poisson or deterministic binary tests. Scores/p-values alone and simulated
#'   tests do not imply an original-scale confidence interval.
#' @param mid_p Use mid-p two-sided Poisson-binomial or resampling tails.
#'   One-sided binary tails are inclusive unless
#'   `control$one_sided_midp=TRUE`. The named Poisson statistics determine their
#'   own tail convention independently of this argument.
#' @param zero_method For zero log-ratios, use the `score` approximation
#'   (`-1/se`, requiring a log-scale SE and unit reference) or `exclude`.
#'   Score intervals at zero are approximations marked in the diagnostics.
#' @param id Optional unique provider identifiers.
#' @param na_action `fail` rejects missing inputs; `omit` omits incomplete
#'   providers from computation but preserves their rows with missing results.
#' @param control Named list of numerical and estimation controls; see Details.
#'
#' @details
#' ## Statistic and null distribution
#'
#' Choose `statistic` for the input and initial test, then `calibration` for
#' the score distribution. Under fixed or estimated normal nulls the adjusted
#' score is \eqn{(z_i-\mu_i)/\sigma_i}; its two-sided p-value is
#' \eqn{2\Phi(-|(z_i-\mu_i)/\sigma_i|)}. The theoretical mode preserves raw
#' tail probabilities, including Student-t and discrete tails. Thus theoretical
#' calibration and fixed `null_mean=0, null_sd=1` need not agree for discrete
#' tests whose score conversion reaches a numerical floor.
#'
#' For a count \eqn{O} under its null distribution, define the one-sided
#' mid-p tails \eqn{L=P(X<O)+P(X=O)/2} and
#' \eqn{U=P(X>O)+P(X=O)/2}. The two-sided mid-p is
#' \eqn{2\min(L,U)}. Equivalently, with \eqn{q_l=2L} and \eqn{q_u=2U},
#' it is \eqn{\min(q_l,q_u)}, without a further division by two.
#' Mid-p is a half-mass convention and need not provide conservative exact
#' test size. For `poisson_exact`, the two-sided convention instead doubles
#' the inclusive upper tail if `observed > expected * reference`, or the
#' inclusive lower tail otherwise, and caps the result at 0.999. It is not
#' a probability-ordering definition of a two-sided Poisson test.
#'
#' `overall` defaults to robust location/scale estimation. `groupwise` fits
#' that estimator separately by supplied groups or size groups, using
#' [MASS::rlm()] with MAD scale and the selected influence function.
#' `control=list(estimator="likelihood")` selects the overall likelihood
#' estimator. The structured `individualized`, `correlated`, `cre`, and
#' `moment` models use the estimators exposed by [fit_empirical_null()].
#' Moment calibration requires Poisson Pearson scores, unit reference, and
#' `second_moment`; for supplied scores, set `control$score_type="pearson"`.
#' Exact moment expressions do not make the working Gaussian tails exact
#' Poisson tests.
#'
#' ## Numerical and estimation controls
#'
#' Robust controls are `estimator`, `maxit` (1000), `acc` (1e-8), `tuning`
#' (1.345 for Huber, 4.685 for bisquare), `min_group_size` (3), and
#' `small_group` (`"error"`, or `"theoretical"` for a recorded fallback).
#' Likelihood controls are `p_grid`, `phi_bounds`, `cutoff`,
#' `include_intercept`, and individualized `common_intercept`, as applicable
#' to the selected model. Inapplicable calibration controls are rejected.
#'
#' Statistics controls are `p_floor` (one-sided normal-score floor, 1e-6),
#' `seed` (1), `n_resample` (10000), and `one_sided_midp` (FALSE).
#' Resampling preserves the caller's random-number state. Its tail proportions
#' are Monte Carlo estimates and may be zero; they are not exact probabilities.
#' Interval controls are `normal_threshold` (expected-count threshold for
#' Byar intervals, 100), `upper_cap` (Poisson root-search limit, 1e8),
#' `root_tol` (1e-8), and `max_log_odds` (binary inversion limit, 40).
#'
#' ## Intervals and interpretation
#'
#' Estimated/fixed normal nulls use signed normal scores; extreme discrete
#' tails can be limited by the chosen score floor. Normal-score
#' null parameters are held fixed during interval inversion; uncertainty in
#' estimated null parameters is not propagated. Binary inversion intervals
#' are on the log-odds-shift scale relative to supplied null probabilities,
#' not automatically on a standardized provider rate or ratio scale.
#' Numerical inversion failures, empty sets, and detected nonmonotone binary
#' score curves are reported in diagnostics and do not produce fabricated
#' interval endpoints. The binary monotonicity check is a numerical guard,
#' not a proof for every probability vector.
#'
#' Inference directions are -1 (lower), 0 (not flagged), +1 (higher), without
#' assigning clinical desirability. A supplied covariance/posterior SE may
#' come from a robust model fit; this function does not recompute covariance
#' matrices or refit models for likelihood-ratio tests. Such model-derived
#' test p-values can enter through `statistic="p_value"` with signed direction.
#' Converting such a supplied two-sided p-value to one-sided inference assumes
#' the symmetric-tail convention encoded by its direction. No multiplicity
#' adjustment is applied across providers.
#'
#' @return An object of class `fit_pci` containing:
#' \describe{
#'   \item{results}{One row per input provider, with raw and adjusted scores
#'     and p-values, null means/SDs, directions, interval endpoints and scale,
#'     and indicators `included` and `used_for_null`. `as.data.frame()` returns
#'     this table. `p_raw` and `p_value` use the requested alternative.}
#'   \item{null_fit, null_parameters}{Estimator output and the parameters used
#'     for calibration. Specified theoretical/fixed nulls have no fitted model.}
#'   \item{diagnostics}{Status, warning codes, omitted identifiers, and
#'     calibration/interval diagnostics. Inspect these before interpretation.}
#'   \item{method, statistic_details}{Selected methods, numerical settings,
#'     and statistic-specific input/calculation details.}
#'   \item{call, package_version}{Call and package-version metadata.}
#' }
#' @seealso [fit_empirical_null()], [empirical_null_groupwise()],
#'   `vignette("fit_pci", package = "EmpiNull")`
#' @export
#' @examples
#' # A specified score-distribution SD; this does not set an outcome benchmark.
#' fit_pci(z = c(-3, 0, 3), calibration = "fixed", null_sd = 1.81)
#'
#' # Poisson mid-p inference for observed/expected ratios.
#' fit_pci(observed = c(0, 13), expected = c(2, 7.385),
#'         statistic = "poisson_midp")
#'
#' # Test proportions against a 30% benchmark using original-scale SEs.
#' fit_pci(estimate = c(0.2, 0.4), se = c(0.04, 0.06),
#'         statistic = "wald", transform = "logit", reference = 0.3)
#'
#' # Group-specific robust location and scale from provider scores.
#' set.seed(17)
#' grouped <- fit_pci(z = rnorm(80, 0.1, 1.2), size = seq_len(80),
#'                    calibration = "groupwise", n_groups = 4, psi = "huber")
#' head(as.data.frame(grouped))
#' grouped$diagnostics$warning_codes
fit_pci <- function(
    z = NULL, estimate = NULL, se = NULL, observed = NULL, expected = NULL,
    p_value = NULL, direction = NULL,
    statistic = c("z", "wald", "log_ratio", "poisson_midp", "poisson_exact",
                  "poisson_score", "poibin", "binary_score", "resampling", "p_value"),
    calibration = c("theoretical", "fixed", "overall", "groupwise",
                    "individualized", "correlated", "cre", "moment"),
    reference = NULL, transform = c("identity", "log", "logit"),
    se_scale = c("original", "transformed"), df = Inf,
    size = NULL, covariates = NULL, second_moment = NULL,
    null_mean = 0, null_sd = 1, psi = c("bisquare", "huber"),
    groups = NULL, n_groups = 4, grouping = c("quantile", "rank"),
    common_mean = NULL, fit_mask = NULL,
    probabilities = NULL, eta = NULL, re_mean = NULL, re_var = NULL,
    alpha = 0.05, alternative = c("two.sided", "greater", "less"),
    ci = c("auto", "none", "wald", "poisson", "invert"), mid_p = TRUE,
    zero_method = c("score", "exclude"), id = NULL,
    na_action = c("fail", "omit"), control = list()) {
  original_call <- match.call()
  if (identical(statistic, "poibin_exact")) statistic <- "poibin"
  statistic <- match.arg(statistic)
  calibration <- match.arg(calibration)
  if (identical(alternative, "two_sided")) alternative <- "two.sided"
  alternative <- match.arg(alternative)
  ci <- match.arg(ci)
  transform <- if (statistic == "log_ratio" && missing(transform)) "log" else match.arg(transform)
  se_scale <- if (statistic == "log_ratio" && missing(se_scale)) "transformed" else match.arg(se_scale)
  psi <- match.arg(psi)
  grouping <- match.arg(grouping)
  zero_method <- match.arg(zero_method)
  na_action <- match.arg(na_action)
  .check_numeric_vector(alpha, "alpha", n = 1L, positive = TRUE)
  if (alpha >= 1) stop("alpha must be smaller than one.", call. = FALSE)
  if (!is.logical(mid_p) || length(mid_p) != 1L || is.na(mid_p)) {
    stop("mid_p must be TRUE or FALSE.", call. = FALSE)
  }
  .pci_check_control(control)
  raw_inputs <- list(z = z, estimate = estimate, se = se, observed = observed,
                    expected = expected, p_value = p_value, direction = direction,
                    reference = reference, probabilities = probabilities,
                    eta = eta, re_mean = re_mean, re_var = re_var)
  allowed_inputs <- switch(statistic,
    z = "z", p_value = c("p_value", "direction"),
    wald = c("estimate", "se", "reference"),
    log_ratio = c("estimate", "se", "reference"),
    poisson_midp = c("observed", "expected", "reference"),
    poisson_exact = c("observed", "expected", "reference"),
    poisson_score = c("observed", "expected", "reference"),
    c("observed", "probabilities", "eta", "re_mean", "re_var", "reference"))
  unused <- setdiff(names(raw_inputs)[!vapply(raw_inputs, is.null, logical(1L))], allowed_inputs)
  if (length(unused)) stop("Inputs not used by statistic='", statistic, "': ", paste(unused, collapse = ", "), ".", call. = FALSE)
  if (!statistic %in% c("wald", "log_ratio") && any(!is.na(df) & df != Inf)) {
    stop("df applies only to Wald or log-ratio statistics.", call. = FALSE)
  }
  if (statistic == "log_ratio" && transform != "log") {
    stop("log_ratio requires transform='log'.", call. = FALSE)
  }
  if (!statistic %in% c("wald", "log_ratio") && transform != "identity") {
    stop("transform applies only to Wald or log-ratio statistics.", call. = FALSE)
  }

  primary <- switch(statistic, z = z, p_value = p_value,
                    wald = estimate, log_ratio = estimate, observed)
  if (is.null(primary) || !is.numeric(primary) || !is.null(dim(primary)) || !length(primary)) {
    stop("Supply the provider-level input required by statistic.", call. = FALSE)
  }
  n <- length(primary)
  if (is.null(id)) id <- seq_len(n)
  if (!is.atomic(id) || !is.null(dim(id)) || length(id) != n || anyNA(id) || anyDuplicated(id)) {
    stop("id must contain one unique, nonmissing identifier per provider.", call. = FALSE)
  }
  if (is.null(fit_mask)) fit_mask <- rep(TRUE, n)
  if (!is.logical(fit_mask) || length(fit_mask) != n || anyNA(fit_mask)) {
    stop("fit_mask must be a nonmissing logical vector, one per provider.", call. = FALSE)
  }
  vectors <- list(z = z, estimate = estimate, se = se, observed = observed,
                  expected = expected, p_value = p_value, direction = direction,
                  reference = reference, df = df, size = size,
                  second_moment = second_moment, null_mean = null_mean, null_sd = null_sd)
  for (nm in names(vectors)) {
    x <- vectors[[nm]]
    if (!is.null(x) && (!is.numeric(x) || !is.null(dim(x)) ||
                        !length(x) || !length(x) %in% c(1L, n))) {
      stop(nm, " must be numeric with length one or the number of providers.", call. = FALSE)
    }
  }
  if (!is.null(groups) && (!is.atomic(groups) || length(groups) != n || !is.null(dim(groups)))) {
    stop("groups must have one label per provider.", call. = FALSE)
  }
  if (!is.null(covariates)) {
    if (is.null(dim(covariates))) covariates <- matrix(covariates, ncol = 1L)
    if (nrow(covariates) != n) stop("covariates must have one row per provider.", call. = FALSE)
  }
  keep <- rep(TRUE, n)
  for (x in vectors) if (!is.null(x)) keep <- keep & !is.na(rep_len(x, n))
  if (!is.null(groups)) keep <- keep & !is.na(groups)
  if (!is.null(covariates)) keep <- keep & stats::complete.cases(covariates)
  lists <- list(probabilities = probabilities, eta = eta, re_mean = re_mean, re_var = re_var)
  for (nm in names(lists)) {
    x <- lists[[nm]]
    if (is.null(x)) next
    if (!is.list(x)) {
      if (n != 1L || !is.numeric(x)) stop(nm, " must be a list with one vector per provider.", call. = FALSE)
      x <- list(x)
    }
    if (length(x) != n) stop(nm, " must have one entry per provider.", call. = FALSE)
    keep <- keep & !vapply(x, anyNA, logical(1L))
    lists[[nm]] <- x
  }
  if (any(!keep) && na_action == "fail") {
    stop("Missing provider inputs; use na_action='omit' to retain incomplete rows with NA results.", call. = FALSE)
  }
  if (!any(keep)) stop("No complete providers remain for inference.", call. = FALSE)
  subset_vector <- function(x) if (is.null(x)) NULL else rep_len(x, n)[keep]
  vv <- lapply(vectors, subset_vector)
  ll <- lapply(lists, function(x) if (is.null(x)) NULL else x[keep])
  if (calibration == "moment") {
    supplied_pearson <- statistic == "z" && identical(control$score_type, "pearson")
    if (statistic != "poisson_score" && !supplied_pearson) {
      stop("Moment calibration requires statistic='poisson_score', or z with control$score_type='pearson'; mid-p scores are not Pearson scores.", call. = FALSE)
    }
    if (is.null(vv$size)) vv$size <- vv$expected
    if (statistic == "poisson_score" && !is.null(vv$reference) && any(vv$reference != 1)) {
      stop("Moment calibration requires unit Poisson reference; supply the appropriate expected counts and second moments explicitly.", call. = FALSE)
    }
    if (statistic == "poisson_score" && !is.null(vv$size) &&
        !isTRUE(all.equal(vv$size, vv$expected, check.attributes = FALSE))) {
      stop("Moment calibration size must equal the expected-count vector.", call. = FALSE)
    }
  }

  stat <- .pci_compute_statistics(
    n = sum(keep), statistic = statistic, z = vv$z, estimate = vv$estimate,
    se = vv$se, observed = vv$observed, expected = vv$expected,
    p_value = vv$p_value, direction = vv$direction, reference = vv$reference,
    transform = transform, se_scale = se_scale, df = vv$df,
    probabilities = ll$probabilities, eta = ll$eta, re_mean = ll$re_mean,
    re_var = ll$re_var, mid_p = mid_p, zero_method = zero_method,
    alternative = alternative,
    control = control[intersect(names(control), .pci_stat_control())]
  )
  cal <- .pci_calibrate(
    z = stat$z, calibration = calibration, size = vv$size,
    covariates = if (is.null(covariates)) NULL else covariates[keep, , drop = FALSE],
    second_moment = vv$second_moment, null_mean = vv$null_mean,
    null_sd = vv$null_sd, psi = psi,
    groups = if (is.null(groups)) NULL else groups[keep],
    n_groups = n_groups, grouping = grouping, common_mean = common_mean,
    fit_mask = fit_mask[keep],
    control = control[intersect(names(control), .pci_cal_control())]
  )
  z_adj <- (stat$z - cal$null_mean) / cal$null_sd
  p <- if (calibration == "theoretical") stat$p else .pci_normal_tail(z_adj, alternative)
  flag <- ifelse(is.na(p), NA_integer_,
                 ifelse(p < alpha, ifelse(z_adj < 0, -1L, ifelse(z_adj > 0, 1L, 0L)), 0L))
  # For one-sided inference flag only the tested direction.
  if (alternative == "greater") flag[!is.na(flag) & flag < 0] <- 0L
  if (alternative == "less") flag[!is.na(flag) & flag > 0] <- 0L
  interval <- .pci_intervals(
    stat = stat, statistic = statistic, calibration = calibration,
    mu = cal$null_mean, sd = cal$null_sd, alpha = alpha,
    alternative = alternative, ci = ci, zero_method = zero_method,
    mid_p = mid_p, probabilities = ll$probabilities,
    control = control[intersect(names(control), .pci_ci_control())]
  )
  expand <- function(x, fill = NA_real_) {
    out <- rep(fill, n)
    if (!is.null(x)) out[keep] <- rep_len(x, sum(keep))
    out
  }
  result <- data.frame(
    provider_id = id, estimate = expand(stat$estimate), se = expand(stat$se),
    reference = expand(stat$reference), observed = expand(stat$observed),
    expected = expand(stat$expected), statistic_value = expand(stat$test_statistic),
    df = expand(stat$df), z_raw = expand(stat$z), p_raw = expand(stat$p),
    null_mean = expand(cal$null_mean), null_sd = expand(cal$null_sd),
    z_adjusted = expand(z_adj), p_value = expand(p), direction = expand(flag, NA_integer_),
    group = expand(as.character(cal$group), NA_character_),
    ci_lower = expand(interval$lower), ci_upper = expand(interval$upper),
    ci_scale = expand(interval$scale, NA_character_),
    included = keep & !is.na(expand(stat$z)),
    used_for_null = keep & fit_mask & is.finite(expand(stat$z)) &
      !calibration %in% c("theoretical", "fixed"),
    stringsAsFactors = FALSE
  )
  warning_codes <- unique(c(cal$warning_codes, interval$warning_codes,
                            if (any(!keep)) "INCOMPLETE_PROVIDERS_OMITTED",
                            if (anyNA(stat$z)) "UNDEFINED_SCORES_EXCLUDED"))
  floor_value <- if (is.null(control$p_floor)) 1e-6 else control$p_floor
  floored <- statistic %in% c("poisson_midp", "poisson_exact", "poibin", "resampling") &&
    any(abs(stat$z) >= stats::qnorm(floor_value, lower.tail = FALSE) - 1e-12, na.rm = TRUE)
  if (statistic == "p_value" && any(vv$p_value == 0)) floored <- TRUE
  if (floored) warning_codes <- unique(c(warning_codes, "NORMAL_SCORE_FLOOR_APPLIED"))
  diagnostics <- list(
    status = if (length(warning_codes)) "warning" else "ok",
    warning_codes = warning_codes, n_providers = n,
    n_included = sum(result$included), omitted_ids = id[!keep],
    calibration = cal$diagnostics, intervals = interval$diagnostics,
    parameter_uncertainty = "Null parameters held fixed for p-values and intervals; estimation uncertainty is not propagated.",
    interpretation = "Directions describe higher/lower outcomes, not better/worse performance. Structured null parameters do not identify latent patient-level confounders."
  )
  structure(list(
    results = result, null_fit = cal$fit, null_parameters = cal$parameters,
    diagnostics = diagnostics,
    method = list(statistic = statistic, calibration = calibration,
                  transform = stat$transform, alternative = alternative,
                  alpha = alpha, interval = interval$method,
                  mid_p = if (statistic == "poisson_midp") TRUE else
                    if (statistic == "poisson_exact") FALSE else
                      if (statistic %in% c("poibin", "resampling")) mid_p else NULL,
                  score_floor = if (is.null(control$p_floor)) 1e-6 else control$p_floor),
    statistic_details = stat$metadata, call = original_call,
    package_version = "1.2.0"
  ), class = "fit_pci")
}

.pci_stat_control <- function() c("p_floor", "seed", "n_resample", "one_sided_midp")
.pci_cal_control <- function() c("estimator", "maxit", "acc", "tuning", "min_group_size",
                                "small_group", "p_grid", "phi_bounds", "cutoff",
                                "include_intercept", "common_intercept")
.pci_ci_control <- function() c("normal_threshold", "upper_cap", "root_tol", "max_log_odds", "p_floor", "one_sided_midp")
.pci_check_control <- function(control) {
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || any(names(control) == "") || anyDuplicated(names(control))))) {
    stop("control must be a named list with unique names.", call. = FALSE)
  }
  allowed <- unique(c(.pci_stat_control(), .pci_cal_control(), .pci_ci_control(), "score_type"))
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) stop("Unknown control option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  invisible(control)
}

.pci_normal_tail <- function(z, alternative) {
  switch(alternative, two.sided = 2 * stats::pnorm(-abs(z)),
         greater = stats::pnorm(z, lower.tail = FALSE), less = stats::pnorm(z))
}

#' @export
print.fit_pci <- function(x, ...) {
  cat("Post-model provider inference\n")
  cat("  statistic:", x$method$statistic, " calibration:", x$method$calibration, "\n")
  cat("  providers:", nrow(x$results), " included:", x$diagnostics$n_included, "\n")
  cat("  alternative:", x$method$alternative, " alpha:", x$method$alpha, "\n")
  if (length(x$diagnostics$warning_codes)) {
    cat("  warnings:", paste(x$diagnostics$warning_codes, collapse = ", "), "\n")
  }
  print(utils::head(x$results), row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.fit_pci <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$results, row.names = row.names, optional = optional, ...)
}
