test_that("fixed 1.81 calibration reproduces the PPPW reference", {
  fit <- fit_pci(z = c(-3, 0, 3), calibration = "fixed", null_sd = 1.81,
                 id = c("low", "middle", "high"), ci = "none")
  expect_s3_class(fit, "fit_pci")
  expect_identical(fit$results$provider_id, c("low", "middle", "high"))
  expect_equal(fit$results$null_mean, rep(0, 3))
  expect_equal(fit$results$null_sd, rep(1.81, 3))
  expect_equal(fit$results$z_adjusted, c(-3, 0, 3) / 1.81)
  expect_equal(fit$results$p_value, 2 * pnorm(-abs(c(-3, 0, 3)) / 1.81))
  expect_gt(fit$results$p_value[3], 0.05)
  expect_lt(fit$results$p_raw[3], 0.05)
})

test_that("normal one-sided tests respect the signed statistic", {
  z <- c(-2, 0, 2)
  less <- fit_pci(z = z, alternative = "less", ci = "none")
  greater <- fit_pci(z = z, alternative = "greater", ci = "none")
  two <- fit_pci(z = z, alternative = "two_sided", ci = "none")
  expect_equal(less$results$p_value, pnorm(z))
  expect_equal(greater$results$p_value, pnorm(z, lower.tail = FALSE))
  expect_equal(two$results$p_value, 2 * pnorm(-abs(z)))
  expect_equal(less$results$p_value + greater$results$p_value, rep(1, 3))
})

test_that("supplied two-sided p-values recover signed normal scores", {
  p <- c(0.0027, 1, 0.05)
  fit <- fit_pci(p_value = p, direction = c(-1, 0, 1),
                 statistic = "p_value", ci = "none")
  expect_equal(fit$results$z_raw, c(-1, 0, 1) * qnorm(p / 2, lower.tail = FALSE))
  expect_equal(fit$results$p_value, p)
  expect_equal(fit$results$p_raw, p)
})

test_that("Wald inference uses original-scale references and calibrated intervals", {
  fit <- fit_pci(estimate = c(2, 4), se = c(0.5, 1), reference = 1,
                 statistic = "wald", calibration = "fixed",
                 null_mean = 0.4, null_sd = 1.3)
  crit <- qnorm(0.975)
  center <- c(2, 4) - 0.4 * c(0.5, 1)
  expect_equal(fit$results$z_raw, c(2, 3))
  expect_equal(fit$results$ci_lower, center - crit * 1.3 * c(0.5, 1))
  expect_equal(fit$results$ci_upper, center + crit * 1.3 * c(0.5, 1))
  excludes_reference <- fit$results$ci_lower > 1 | fit$results$ci_upper < 1
  expect_identical(fit$results$p_value < 0.05, excludes_reference)
})

test_that("theoretical Wald t-tests preserve Student reference and intervals", {
  fit <- fit_pci(estimate = 1, se = 0.5, reference = 0,
                 statistic = "wald", df = 8)
  expect_equal(fit$results$p_raw, 2 * pt(-2, df = 8))
  expect_equal(fit$results$p_value, 2 * pt(-2, df = 8))
  expect_equal(fit$results$ci_lower, 1 - qt(0.975, 8) * 0.5)
  expect_equal(fit$results$ci_upper, 1 + qt(0.975, 8) * 0.5)
})

test_that("calibrated t intervals invert normal-quantile calibration", {
  fit <- fit_pci(estimate = 2, se = 0.7, reference = 0, statistic = "wald", df = 5,
                 calibration = "fixed", null_mean = 0.4, null_sd = 1.3)
  endpoints <- c(fit$results$ci_lower, fit$results$ci_upper)
  endpoint_tests <- fit_pci(estimate = rep(2, 2), se = 0.7, reference = endpoints,
                            statistic = "wald", df = 5, calibration = "fixed",
                            null_mean = 0.4, null_sd = 1.3, ci = "none")
  expect_equal(endpoint_tests$results$p_value, rep(0.05, 2), tolerance = 1e-8)
})

test_that("log-ratio intervals use transformed standard errors and null centering", {
  fit <- fit_pci(estimate = c(0.8, 2), se = c(0.2, 0.3), reference = 1,
                 statistic = "log_ratio", se_scale = "transformed",
                 calibration = "fixed", null_mean = 0.35, null_sd = 1.2)
  center <- log(c(0.8, 2)) - 0.35 * c(0.2, 0.3)
  width <- qnorm(0.975) * 1.2 * c(0.2, 0.3)
  expect_equal(fit$results$z_raw, log(c(0.8, 2)) / c(0.2, 0.3))
  expect_equal(fit$results$ci_lower, exp(center - width))
  expect_equal(fit$results$ci_upper, exp(center + width))
  expect_identical(fit$results$p_value < 0.05,
                   fit$results$ci_lower > 1 | fit$results$ci_upper < 1)
})

test_that("an explicitly original-scale ratio SE receives the delta transformation", {
  original <- fit_pci(estimate = 2, se = 0.4, statistic = "log_ratio",
                      se_scale = "original", ci = "none")
  transformed <- fit_pci(estimate = 2, se = 0.2, statistic = "log_ratio",
                         se_scale = "transformed", ci = "none")
  expect_equal(original$results$z_raw, transformed$results$z_raw)
  expect_equal(original$results$p_value, transformed$results$p_value)
})

test_that("logit Wald references and intervals remain on the probability scale", {
  rate <- c(0.12, 0.35)
  se <- c(0.02, 0.04)
  reference <- 0.2
  fit <- fit_pci(estimate = rate, se = se, reference = reference,
                 statistic = "wald", transform = "logit",
                 calibration = "fixed", null_mean = 0.3, null_sd = 1.81)
  se_logit <- se / (rate * (1 - rate))
  expect_equal(fit$results$z_raw, (qlogis(rate) - qlogis(reference)) / se_logit)
  center <- qlogis(rate) - 0.3 * se_logit
  width <- qnorm(0.975) * 1.81 * se_logit
  expect_equal(fit$results$ci_lower, plogis(center - width))
  expect_equal(fit$results$ci_upper, plogis(center + width))
  expect_true(all(fit$results$ci_lower > 0 & fit$results$ci_upper < 1))
})

test_that("one-sided Wald intervals invert the selected alternative", {
  less <- fit_pci(estimate = 1, se = 0.2, statistic = "wald", alternative = "less")
  greater <- fit_pci(estimate = 1, se = 0.2, statistic = "wald", alternative = "greater")
  expect_equal(less$results$ci_lower, -Inf)
  expect_equal(less$results$ci_upper, 1 + qnorm(0.95) * 0.2)
  expect_equal(greater$results$ci_lower, 1 - qnorm(0.95) * 0.2)
  expect_equal(greater$results$ci_upper, Inf)
})

test_that("Poisson mid-p matches the hand-worked example and zero-count formula", {
  obs <- c(13, 0)
  expected <- c(7.385003064460942, 2)
  fit <- fit_pci(observed = obs, expected = expected, statistic = "poisson_midp", ci = "none")
  lower <- ppois(obs - 1, expected) + 0.5 * dpois(obs, expected)
  upper <- ppois(obs, expected, lower.tail = FALSE) + 0.5 * dpois(obs, expected)
  expect_equal(fit$results$p_raw, 2 * pmin(lower, upper))
  expect_equal(fit$results$p_value[1], 0.057842266323, tolerance = 1e-10)
  expect_equal(fit$results$p_value[2], exp(-2))
  expect_equal(fit$results$z_raw, cal_Z_htaz(obs, expected))
})

test_that("one-sided Poisson mid-p tails include half the observed mass", {
  obs <- c(0, 3, 8)
  expected <- rep(3, 3)
  lower <- ppois(obs - 1, expected) + 0.5 * dpois(obs, expected)
  upper <- ppois(obs, expected, lower.tail = FALSE) + 0.5 * dpois(obs, expected)
  less <- fit_pci(observed = obs, expected = expected, statistic = "poisson_midp",
                  alternative = "less", ci = "none")
  greater <- fit_pci(observed = obs, expected = expected, statistic = "poisson_midp",
                     alternative = "greater", ci = "none")
  expect_equal(less$results$p_value, lower)
  expect_equal(greater$results$p_value, upper)
})

test_that("exact Poisson inference returns full tails and Garwood ratio bounds", {
  obs <- c(0, 13)
  expected <- c(2, 7.385003064460942)
  fit <- fit_pci(observed = obs, expected = expected, statistic = "poisson_exact")
  p <- c(2 * ppois(0, 2), 2 * ppois(12, expected[2], lower.tail = FALSE))
  expect_equal(fit$results$p_value, p)
  expect_equal(fit$results$ci_lower, c(0, qchisq(0.025, 26) / 2 / expected[2]))
  expect_equal(fit$results$ci_upper, qchisq(0.975, 2 * (obs + 1)) / 2 / expected)
})

test_that("Poisson mid-p interval endpoints invert the calibrated test", {
  obs <- 13
  expected <- 7.385003064460942
  fit <- fit_pci(observed = obs, expected = expected, statistic = "poisson_midp",
                 calibration = "fixed", null_mean = 0.3, null_sd = 1.2, ci = "invert")
  endpoints <- c(fit$results$ci_lower, fit$results$ci_upper)
  expect_true(all(is.finite(endpoints)))
  expect_lt(endpoints[1], endpoints[2])
  endpoint_tests <- fit_pci(observed = rep(obs, 2), expected = rep(expected, 2),
                            reference = endpoints, statistic = "poisson_midp",
                            calibration = "fixed", null_mean = 0.3,
                            null_sd = 1.2, ci = "none")
  expect_equal(endpoint_tests$results$p_value, rep(0.05, 2), tolerance = 1e-6)
})

test_that("Poisson mid-p zero-count confidence limits remain finite", {
  fit <- fit_pci(observed = 0, expected = 2, statistic = "poisson_midp", ci = "invert")
  expect_equal(fit$results$ci_lower, 0)
  expect_equal(fit$results$ci_upper, -log(0.05) / 2, tolerance = 1e-6)
})

test_that("theoretical Poisson intervals retain tail accuracy below the score floor", {
  alpha <- 1e-8
  fit <- fit_pci(observed = 13, expected = 7.385, statistic = "poisson_midp",
                 alpha = alpha, ci = "invert")
  endpoints <- c(fit$results$ci_lower, fit$results$ci_upper)
  expect_true(all(is.finite(endpoints) & endpoints > 0))
  endpoint_tests <- fit_pci(observed = c(13, 13), expected = 7.385,
                            reference = endpoints, statistic = "poisson_midp", ci = "none")
  expect_equal(endpoint_tests$results$p_value / alpha, c(1, 1), tolerance = 1e-6)
})

test_that("Poisson-binomial mid-p reduces to binomial mid-p for equal probabilities", {
  observed <- c(0, 4, 8)
  probs <- rep(list(rep(0.3, 10)), 3)
  fit <- fit_pci(observed = observed, probabilities = probs, statistic = "poibin", ci = "none")
  lower <- pbinom(observed - 1, 10, 0.3) + 0.5 * dbinom(observed, 10, 0.3)
  upper <- pbinom(observed, 10, 0.3, lower.tail = FALSE) + 0.5 * dbinom(observed, 10, 0.3)
  expect_equal(fit$results$p_value, 2 * pmin(lower, upper), tolerance = 1e-12)
})

test_that("binary one-sided tails distinguish full tails from explicit mid-p", {
  full <- fit_pci(observed = 4, probabilities = rep(0.3, 10), statistic = "poibin",
                  alternative = "greater", ci = "none")
  mid <- fit_pci(observed = 4, probabilities = rep(0.3, 10), statistic = "poibin",
                 alternative = "greater", ci = "none",
                 control = list(one_sided_midp = TRUE))
  expect_equal(full$results$p_value, pbinom(3, 10, 0.3, lower.tail = FALSE))
  expect_equal(mid$results$p_value,
                pbinom(4, 10, 0.3, lower.tail = FALSE) + 0.5 * dbinom(4, 10, 0.3))
})

test_that("binary score tests use the Bernoulli variance", {
  probabilities <- list(c(0.1, 0.2, 0.4, 0.5), c(0.3, 0.4, 0.6))
  observed <- c(3, 0)
  fit <- fit_pci(observed = observed, probabilities = probabilities,
                 statistic = "binary_score", ci = "none")
  means <- vapply(probabilities, sum, numeric(1))
  vars <- vapply(probabilities, function(p) sum(p * (1 - p)), numeric(1))
  z <- (observed - means) / sqrt(vars)
  expect_equal(fit$results$z_raw, z)
  expect_equal(fit$results$p_value, 2 * pnorm(-abs(z)))
})

test_that("binary intervals invert tests on the log-odds-shift scale", {
  probabilities <- c(0.1, 0.3, 0.4, 0.8)
  for (statistic in c("poibin", "binary_score")) {
    fit <- fit_pci(observed = 2, probabilities = probabilities,
                   statistic = statistic, calibration = "fixed",
                   null_mean = 0.3, null_sd = 1.2, ci = "invert")
    endpoints <- c(fit$results$ci_lower, fit$results$ci_upper)
    expect_true(all(is.finite(endpoints)))
    expect_equal(fit$results$ci_scale, "log_odds_shift")
    shifted <- lapply(endpoints, function(delta) plogis(qlogis(probabilities) + delta))
    endpoint_tests <- fit_pci(observed = c(2, 2), probabilities = shifted,
                              statistic = statistic, calibration = "fixed",
                              null_mean = 0.3, null_sd = 1.2, ci = "none")
    expect_equal(endpoint_tests$results$p_value, rep(0.05, 2), tolerance = 1e-6)
  }
})

test_that("binary theoretical intervals match full and mid-p tail conventions", {
  probabilities <- c(0.1, 0.3, 0.4, 0.8)
  for (mid in c(TRUE, FALSE)) {
    for (alternative in c("two.sided", "less", "greater")) {
      fit <- fit_pci(observed = 2, probabilities = probabilities, statistic = "poibin",
                     mid_p = mid, alternative = alternative, ci = "invert",
                     control = list(one_sided_midp = mid))
      endpoints <- c(fit$results$ci_lower, fit$results$ci_upper)
      finite <- endpoints[is.finite(endpoints)]
      expect_length(finite, if (alternative == "two.sided") 2 else 1)
      shifted <- lapply(finite, function(delta) plogis(qlogis(probabilities) + delta))
      endpoint_tests <- fit_pci(observed = rep(2, length(finite)), probabilities = shifted,
                                statistic = "poibin", mid_p = mid,
                                alternative = alternative, ci = "none",
                                control = list(one_sided_midp = mid))
      expect_equal(endpoint_tests$results$p_value, rep(0.05, length(finite)), tolerance = 1e-6)
    }
  }
})

test_that("deterministic binary inputs give uninformative or unavailable intervals", {
  at_null <- fit_pci(observed = 2, probabilities = c(0, 0, 1, 1), statistic = "poibin")
  expect_equal(at_null$results$p_value, 1)
  expect_equal(at_null$results$ci_lower, -Inf)
  expect_equal(at_null$results$ci_upper, Inf)
  impossible <- fit_pci(observed = 1, probabilities = c(0, 0, 0), statistic = "poibin")
  expect_equal(impossible$results$p_value, 0)
  expect_true(is.na(impossible$results$ci_lower))
  expect_true(is.na(impossible$results$ci_upper))
  expect_true("EMPTY_INVERTED_INTERVAL" %in% impossible$diagnostics$warning_codes)
})

test_that("zero expected counts do not produce ratio confidence limits", {
  fit <- fit_pci(observed = c(0, 1), expected = 0, statistic = "poisson_exact")
  expect_equal(fit$results$p_value, c(1, 0))
  expect_true(all(is.na(fit$results$ci_lower)))
  expect_true(all(is.na(fit$results$ci_upper)))
  expect_true("ZERO_EXPECTED_RATIO_UNDEFINED" %in% fit$diagnostics$warning_codes)
})

test_that("score-floor limits produce unavailable calibrated intervals with a reason", {
  fit <- fit_pci(observed = 100, expected = 1, statistic = "poisson_midp",
                 calibration = "fixed", null_sd = 5, ci = "invert")
  expect_lt(fit$results$p_raw, 1e-6)
  expect_true(is.finite(fit$results$z_raw))
  expect_true(is.na(fit$results$ci_lower))
  expect_true(is.na(fit$results$ci_upper))
  expect_true("SCORE_FLOOR_LIMITS_INVERSION" %in% fit$diagnostics$warning_codes)
})

test_that("a short binary search range does not manufacture infinite confidence limits", {
  fit <- fit_pci(observed = 5, probabilities = rep(0.5, 10), statistic = "poibin",
                 control = list(max_log_odds = 0.1))
  expect_true(is.na(fit$results$ci_lower))
  expect_true(is.na(fit$results$ci_upper))
  expect_true("INTERVAL_INVERSION_FAILED" %in% fit$diagnostics$warning_codes)
})

test_that("nonmonotone heterogeneous binary scores do not yield a misleading interval", {
  probabilities <- c(0.5, 0.5, 1e-8)
  fit <- fit_pci(observed = 3, probabilities = probabilities,
                 statistic = "binary_score", calibration = "fixed", null_sd = 1.81)
  expect_true(is.na(fit$results$ci_lower))
  expect_true(is.na(fit$results$ci_upper))
  expect_true("NONMONOTONE_SCORE_INVERSION" %in% fit$diagnostics$warning_codes)
  # This point was incorrectly inside a single [-0.73, Inf) interval.
  shifted <- plogis(qlogis(probabilities) + 3)
  rejected_point <- fit_pci(observed = 3, probabilities = shifted,
                            statistic = "binary_score", calibration = "fixed",
                            null_sd = 1.81, ci = "none")
  expect_lt(rejected_point$results$p_value, 0.05)
})

test_that("an empty calibrated zero-count acceptance set is reported explicitly", {
  fit <- fit_pci(observed = 0, expected = 10, statistic = "poisson_midp",
                 calibration = "fixed", null_mean = 2, ci = "invert")
  expect_true(is.na(fit$results$ci_lower))
  expect_true(is.na(fit$results$ci_upper))
  expect_true("EMPTY_INVERTED_INTERVAL" %in% fit$diagnostics$warning_codes)
})

test_that("Poisson zero-count one-sided mid-p intervals match analytic limits", {
  fit <- fit_pci(observed = 0, expected = 2, statistic = "poisson_midp",
                 alternative = "less", ci = "invert")
  expect_equal(fit$results$ci_lower, 0)
  expect_equal(fit$results$ci_upper, -log(2 * 0.05) / 2, tolerance = 1e-6)
})

test_that("zero-ratio approximations and exclusions are marked explicitly", {
  score <- fit_pci(estimate = 0, se = 0.2, statistic = "log_ratio",
                   calibration = "fixed", null_mean = 0.3, null_sd = 1.2)
  expect_equal(score$results$ci_lower, 0)
  expect_equal(score$results$ci_upper, qnorm(0.975) * 1.2 * 0.2)
  expect_true("ZERO_RATIO_SCORE_APPROXIMATION" %in% score$diagnostics$warning_codes)
  exclude <- fit_pci(estimate = c(0, 2), se = 0.2, statistic = "log_ratio",
                     zero_method = "exclude", ci = "none")
  expect_true(is.na(exclude$results$p_value[1]))
  expect_identical(exclude$results$included, c(FALSE, TRUE))
})

test_that("Poisson moment calibration cannot accidentally consume mid-p scores", {
  expect_error(fit_pci(observed = c(1, 2, 3), expected = c(2, 3, 4),
                       statistic = "poisson_midp", calibration = "moment",
                       second_moment = c(2, 3, 4)), "Pearson|pearson|mid-p")
})

test_that("resampling is reproducible without changing the user's random stream", {
  set.seed(921)
  old_seed <- .Random.seed
  first <- fit_pci(observed = c(2, 4),
                   probabilities = list(rep(0.2, 10), rep(0.4, 10)),
                   statistic = "resampling", ci = "none",
                   control = list(seed = 192, n_resample = 4000))
  expect_identical(.Random.seed, old_seed)
  second <- fit_pci(observed = c(2, 4),
                    probabilities = list(rep(0.2, 10), rep(0.4, 10)),
                    statistic = "resampling", ci = "none",
                    control = list(seed = 192, n_resample = 4000))
  expect_equal(first$results$p_value, second$results$p_value)
})

test_that("groupwise robust calibration retains the existing MASS procedure", {
  z <- c(-2, -1, -0.5, 0.1, 0.8, 1.3, 1.9, 5,
         -1.8, -0.9, -0.1, 0.4, 0.9, 1.4, 2.5, 4.8)
  size <- seq_along(z)
  fit <- fit_pci(z = z, size = size, calibration = "groupwise", n_groups = 2,
                 psi = "huber", ci = "none")
  old <- empirical_null_groupwise(z, size, n_groups = 2, psi = MASS::psi.huber)
  expect_equal(fit$results$null_mean, unname(old$intercept[old$group]), tolerance = 1e-7)
  expect_equal(fit$results$null_sd, unname(old$scale[old$group]), tolerance = 1e-7)
  expect_equal(as.character(fit$results$group), as.character(old$group))
})

test_that("fit exclusions do not drop providers from the output", {
  z <- c(-1.2, -0.8, -0.2, 0.1, 0.7, 1.4, 10)
  fit <- fit_pci(z = z, calibration = "overall", psi = "huber",
                 fit_mask = c(rep(TRUE, 6), FALSE), ci = "none")
  expected <- MASS::rlm(z[1:6] ~ 1, psi = MASS::psi.huber, maxit = 1000, acc = 1e-8)
  expect_equal(nrow(fit$results), length(z))
  expect_equal(fit$results$null_mean, rep(unname(coef(expected)[1]), length(z)), tolerance = 1e-7)
  expect_equal(fit$results$null_sd, rep(expected$s, length(z)), tolerance = 1e-7)
  expect_true(is.finite(fit$results$p_value[7]))
})

test_that("custom labels fit separate nulls and small groups have an explicit policy", {
  z <- c(-1, -0.4, 0.3, 1.2, -2, -0.3, 0.4, 1.8)
  groups <- rep(c("small", "large"), each = 4)
  fit <- fit_pci(z = z, groups = groups, calibration = "groupwise", psi = "huber", ci = "none")
  expect_identical(as.character(fit$results$group), groups)
  expect_error(fit_pci(z = c(-1, 0, 1), groups = c("a", "b", "b"),
                       calibration = "groupwise", ci = "none"), "group|Group|few|size")
  fallback <- fit_pci(z = c(-1, 0, 1), groups = c("a", "b", "b"),
                      calibration = "groupwise", ci = "none",
                      control = list(small_group = "theoretical"))
  expect_equal(fallback$results$null_mean, c(0, 0, 0))
  expect_equal(fallback$results$null_sd, c(1, 1, 1))
})

test_that("degenerate robust nulls cannot silently produce infinite adjusted scores", {
  expect_error(fit_pci(z = rep(2, 8), calibration = "overall", ci = "none"),
               "scale|constant|degenerate|variation|zero")
})

test_that("groupwise common means override centers without changing fitted scales", {
  z <- c(-1, -0.4, 0.3, 1.2, 0.5, 1.2, 2.4, 3.8)
  groups <- rep(c("first", "second"), each = 4)
  original <- fit_pci(z = z, groups = groups, calibration = "groupwise", ci = "none")
  common <- fit_pci(z = z, groups = groups, calibration = "groupwise",
                    common_mean = 0, ci = "none")
  expect_equal(common$results$null_mean, rep(0, length(z)))
  expect_equal(common$results$null_sd, original$results$null_sd)
})

test_that("rank grouping is an explicit alternative for tied size quantiles", {
  z <- c(-1.4, -0.3, 0.2, 1.2, -1.1, -0.1, 0.5, 1.6)
  expect_error(fit_pci(z = z, size = rep(10, 8), calibration = "groupwise",
                       n_groups = 2, ci = "none"), "tie|distinct|quantile")
  rank_fit <- fit_pci(z = z, size = rep(10, 8), calibration = "groupwise",
                      n_groups = 2, grouping = "rank", ci = "none")
  expect_equal(as.character(rank_fit$results$group), rep(c("1", "2"), each = 4))
  expect_true(all(is.finite(rank_fit$results$p_value)))
})

test_that("omitting missing inputs retains alignment of provider identifiers", {
  fit <- fit_pci(z = c(-2, NA_real_, 2), id = c("A", "B", "C"),
                 na_action = "omit", ci = "none")
  expect_identical(fit$results$provider_id, c("A", "B", "C"))
  expect_equal(fit$results$z_raw, c(-2, NA_real_, 2))
  expect_equal(fit$results$p_value, c(2 * pnorm(-2), NA_real_, 2 * pnorm(-2)))
  expect_identical(fit$results$included, c(TRUE, FALSE, TRUE))
  expect_error(fit_pci(z = c(0, NA_real_), ci = "none"), "missing|finite|NA")
})

test_that("invalid statistical inputs receive useful errors", {
  expect_error(fit_pci(estimate = 1, se = 0, statistic = "wald"), "se|positive")
  expect_error(fit_pci(observed = 1.5, expected = 2, statistic = "poisson_midp"), "integer|count")
  expect_error(fit_pci(observed = 1, expected = -2, statistic = "poisson_midp"), "positive|expected")
  expect_error(fit_pci(observed = 2, probabilities = list(c(0.2, 1.2)), statistic = "poibin"), "probabilit|\\[0, ?1\\]")
  expect_error(fit_pci(z = 1, calibration = "fixed", null_sd = 0), "positive|sd|scale")
  expect_error(fit_pci(z = 1, alpha = 1), "alpha")
  expect_error(fit_pci(z = 1, observed = 1, expected = 2), "unused|Unused|incompatible|only|not used|not supported|does not")
})

test_that("the Python Poisson-binomial name is accepted as a compatibility alias", {
  fit <- fit_pci(observed = 2, probabilities = rep(0.2, 10),
                 statistic = "poibin_exact", ci = "none")
  canonical <- fit_pci(observed = 2, probabilities = rep(0.2, 10),
                       statistic = "poibin", ci = "none")
  expect_equal(fit$results$p_value, canonical$results$p_value)
  expect_equal(fit$method$statistic, "poibin")
})
