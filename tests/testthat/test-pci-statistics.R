test_that("Poisson mid-p uses two-sided half-mass tails and stable quantiles", {
  expected <- 7.385003064460942
  observed <- 13
  answer <- .pci_compute_statistics(1, "poisson_midp", observed = observed,
                                    expected = expected)
  lower <- ppois(observed - 1, expected) + dpois(observed, expected) / 2
  upper <- ppois(observed, expected, lower.tail = FALSE) + dpois(observed, expected) / 2
  expect_equal(answer$p, 0.057842266323, tolerance = 1e-11)
  expect_equal(answer$p, 2 * min(lower, upper))
  expect_equal(2 * pnorm(-abs(answer$z)), answer$p)
  expect_gt(answer$z, 0)
  lower_answer <- .pci_compute_statistics(1, "poisson_midp", observed = 0,
                                          expected = 4, alternative = "less")
  expect_equal(lower_answer$p, exp(-4) / 2)
  extreme <- .pci_compute_statistics(1, "poisson_midp", observed = 1000,
                                     expected = 1)
  expect_true(is.finite(extreme$z))
  expect_equal(extreme$z, qnorm(1e-6, lower.tail = FALSE))
  expect_equal(extreme$p, 0)
})

test_that("conventional Poisson tests reproduce the directional inclusive rule", {
  observed <- c(0, 10, 15)
  expected <- c(2, 10, 9)
  answer <- .pci_compute_statistics(3, "poisson_exact", observed = observed,
                                    expected = expected)
  manual <- pmin(.999, 2 * ifelse(observed > expected,
                                 ppois(observed - 1, expected, lower.tail = FALSE),
                                 ppois(observed, expected)))
  expect_equal(answer$p, manual)
  greater <- .pci_compute_statistics(3, "poisson_exact", observed = observed,
                                     expected = expected, alternative = "greater")
  expect_equal(greater$p, ppois(observed - 1, expected, lower.tail = FALSE))
  null <- .pci_compute_statistics(1, "poisson_exact", observed = 0, expected = 0)
  expect_equal(null$z, 0)
  expect_equal(null$p, 1)
})

test_that("Poisson score statistics distinguish deterministic and impossible counts", {
  answer <- .pci_compute_statistics(3, "poisson_score", observed = c(0, 1, 13),
                                    expected = c(0, 0, 7))
  expect_equal(answer$test_statistic[3], (13 - 7) / sqrt(7))
  expect_equal(answer$p[1:2], c(1, 0))
  expect_equal(answer$z[1], 0)
  expect_true(is.finite(answer$z[2]))
  expect_true(is.infinite(answer$test_statistic[2]))
})

test_that("Poisson references are ratios against the supplied baseline expected count", {
  answer <- .pci_compute_statistics(1, "poisson_midp", observed = 13,
                                    expected = 7, reference = 1.4)
  expect_equal(answer$expected, 7)
  expect_equal(answer$estimate, 13 / 7)
  expect_equal(answer$reference, 1.4)
  expect_equal(answer$metadata$null_expected, 9.8)
  expect_equal(answer$z, .pci_poisson_z(13, 9.8))
  expect_equal(answer$p, 2 * pnorm(-abs(answer$z)))
})

test_that("Poisson-binomial dynamic programming agrees with exhaustive enumeration", {
  probs <- c(.05, .2, .4, .75, .95)
  outcomes <- as.matrix(expand.grid(rep(list(c(0, 1)), length(probs))))
  weights <- apply(outcomes, 1, function(row) prod(ifelse(row == 1, probs, 1 - probs)))
  counts <- rowSums(outcomes)
  for (observed in 0:length(probs)) {
    answer <- .pci_poibin_tails(observed, probs)
    lower <- sum(weights[counts < observed]) + sum(weights[counts == observed]) / 2
    upper <- sum(weights[counts > observed]) + sum(weights[counts == observed]) / 2
    expect_equal(answer$lower, lower)
    expect_equal(answer$upper, upper)
    expect_equal(answer$two_sided, 2 * min(lower, upper))
    expect_equal(answer$pmf, vapply(0:length(probs), function(x) sum(weights[counts == x]), numeric(1)))
    expect_equal(2 * pnorm(-abs(.pci_binary_z(observed, probs))), answer$two_sided)
  }
})

test_that("binary one-sided tests default to inclusive rather than mid tails", {
  probs <- rep(.3, 8)
  args <- list(n = 1, statistic = "poibin_exact", observed = 1,
               probabilities = probs, alternative = "less")
  expect_equal(do.call(.pci_compute_statistics, args)$p, pbinom(1, 8, .3))
  args$control <- list(one_sided_midp = TRUE)
  expect_equal(do.call(.pci_compute_statistics, args)$p,
               pbinom(0, 8, .3) + dbinom(1, 8, .3) / 2)
  args$alternative <- "two.sided"
  expect_equal(do.call(.pci_compute_statistics, args)$p,
               2 * (pbinom(0, 8, .3) + dbinom(1, 8, .3) / 2))
  args$statistic <- "poibin"
  expect_equal(do.call(.pci_compute_statistics, args)$p,
               2 * (pbinom(0, 8, .3) + dbinom(1, 8, .3) / 2))
})

test_that("deterministic Bernoulli probabilities are preserved", {
  answer <- .pci_poibin_tails(2, c(0, 1, 1))
  expect_equal(answer$pmf, c(0, 0, 1, 0))
  expect_equal(answer$two_sided, 1)
  expect_equal(.pci_binary_z(2, c(0, 1, 1)), 0)
  expect_equal(.pci_poibin_tails(1, c(0, 1, 1))$two_sided, 0)
  score <- .pci_compute_statistics(1, "binary_score", observed = 2,
                                   probabilities = c(0, 1, 1), alternative = "greater")
  expect_equal(score$z, 0)
  expect_equal(score$p, 1)
})

test_that("binary score tests use the variance of unequal Bernoulli probabilities", {
  probs <- c(.1, .3, .7, .8)
  answer <- .pci_compute_statistics(1, "binary_score", observed = 3,
                                    probabilities = probs)
  manual <- (3 - sum(probs)) / sqrt(sum(probs * (1 - probs)))
  expect_equal(answer$z, manual)
  expect_equal(answer$p, 2 * pnorm(-abs(manual)))
  expect_equal(answer$expected, sum(probs))
})

test_that("Wald transformations transform the original-scale reference and SE", {
  answer <- .pci_compute_statistics(2, "wald", estimate = c(.3, .1),
                                    se = c(.04, .03), reference = .2,
                                    transform = "logit")
  transformed_se <- c(.04, .03) / (c(.3, .1) * (1 - c(.3, .1)))
  manual <- (qlogis(c(.3, .1)) - qlogis(.2)) / transformed_se
  expect_equal(answer$z, manual)
  expect_equal(answer$working_reference, rep(qlogis(.2), 2))
  expect_equal(answer$working_se, transformed_se)
  log_answer <- .pci_compute_statistics(1, "wald", estimate = 1.4, se = .2,
                                        reference = 1, transform = "log")
  expect_equal(log_answer$z, log(1.4) / (.2 / 1.4))
  transformed <- .pci_compute_statistics(1, "wald", estimate = 1.4, se = .2,
                                         reference = 1, transform = "log",
                                         se_scale = "transformed")
  expect_equal(transformed$z, log(1.4) / .2)
})

test_that("finite-df Wald scores preserve t tail probabilities on the normal scale", {
  answer <- .pci_compute_statistics(3, "wald", estimate = c(-2, 0, 6),
                                    se = .5, reference = 0, df = c(5, 8, 40))
  manual <- 2 * pt(-abs(c(-4, 0, 12)), df = c(5, 8, 40))
  expect_equal(answer$p, manual)
  expect_equal(2 * pnorm(-abs(answer$z)), manual, tolerance = 1e-12)
  expect_gt(answer$z[3], qnorm(1e-6, lower.tail = FALSE))
  one_sided <- .pci_compute_statistics(1, "wald", estimate = 2, se = .5,
                                       df = 5, alternative = "greater")
  expect_equal(one_sided$p, pt(4, 5, lower.tail = FALSE))
})

test_that("log-ratio zero handling is explicit", {
  answer <- .pci_compute_statistics(2, "log_ratio", estimate = c(0, 1.5),
                                    se = c(.4, .2))
  expect_equal(answer$z, c(-1 / .4, log(1.5) / .2))
  expect_equal(answer$transform, "log")
  excluded <- .pci_compute_statistics(1, "log_ratio", estimate = 0, se = .4,
                                      zero_method = "exclude")
  expect_true(is.na(excluded$z))
  expect_true(is.na(excluded$p))
  expect_error(.pci_compute_statistics(1, "log_ratio", estimate = 0, se = .4,
                                       reference = 2), "reference = 1")
  original <- .pci_compute_statistics(1, "log_ratio", estimate = 1.5, se = .3,
                                      se_scale = "original")
  expect_equal(original$z, log(1.5) / (.3 / 1.5))
  expect_error(.pci_compute_statistics(1, "log_ratio", estimate = 0, se = .4,
                                       se_scale = "original"), "log-scale")
})

test_that("supplied two-sided p-values have signed normal scores without clipping positives", {
  answer <- .pci_compute_statistics(4, "p_value", p_value = c(.05, 1, 1e-12, 0),
                                    direction = c(-1, 0, 1, -1))
  expect_equal(answer$p, c(.05, 1, 1e-12, 0))
  expect_equal(answer$z[1:3], c(qnorm(.025), 0, qnorm(5e-13, lower.tail = FALSE)))
  expect_true(is.finite(answer$z[4]))
  greater <- .pci_compute_statistics(2, "p_value", p_value = .1,
                                     direction = c(-1, 1), alternative = "greater")
  expect_equal(greater$p, c(.95, .05))
  expect_error(.pci_compute_statistics(1, "p_value", p_value = .5, direction = 0),
               "nonzero")
})

test_that("resampling is reproducible and restores an existing RNG state and kind", {
  set.seed(120)
  previous <- .Random.seed
  previous_kind <- RNGkind()
  args <- list(n = 2, statistic = "resampling", observed = c(2, 1),
               probabilities = list(c(.2, .3, .6), c(.1, .9)),
               control = list(seed = 22, n_resample = 10000))
  first <- do.call(.pci_compute_statistics, args)
  expect_identical(.Random.seed, previous)
  expect_identical(RNGkind(), previous_kind)
  second <- do.call(.pci_compute_statistics, args)
  expect_identical(first, second)
  exact <- .pci_compute_statistics(2, "poibin_exact", observed = args$observed,
                                    probabilities = args$probabilities)
  expect_equal(first$p, exact$p, tolerance = .025)
})

test_that("resampling restores the absence of an RNG state", {
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) saved <- get(".Random.seed", .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", saved, .GlobalEnv)
  }, add = TRUE)
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  invisible(.pci_compute_statistics(1, "resampling", observed = 1,
                                     probabilities = c(.2, .8),
                                     control = list(n_resample = 100)))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("posterior resampling includes the specified null effect and random effect", {
  args <- list(n = 1, statistic = "resampling", observed = 1,
               eta = c(-1, 1), re_mean = c(.2, -.1), re_var = c(.4, .7),
               reference = .3, control = list(seed = 81, n_resample = 30000))
  answer <- do.call(.pci_compute_statistics, args)
  expect_true(answer$metadata$posterior_sampling)
  expect_identical(answer, do.call(.pci_compute_statistics, args))
  # Independent quadrature integrates each logistic-normal success chance.
  probs <- vapply(1:2, function(i) integrate(function(x)
    plogis(args$reference + args$eta[i] + args$re_mean[i] + sqrt(args$re_var[i]) * x) * dnorm(x),
    lower = -10, upper = 10)$value, numeric(1))
  exact <- .pci_poibin_tails(1, probs)
  expect_equal(answer$p, exact$two_sided, tolerance = .025)
  expect_equal(answer$metadata$simulation_mean, sum(probs), tolerance = .025)
})

test_that("invalid raw statistics fail with targeted errors", {
  expect_error(.pci_compute_statistics(1, "poisson_midp", observed = 1.1,
                                       expected = 2), "integer counts")
  expect_error(.pci_compute_statistics(1, "poibin_exact", observed = 3,
                                       probabilities = c(.1, .2)), "support")
  expect_error(.pci_poibin_tails(0, c(-.1, .2)), "\\[0, 1\\]")
  expect_error(.pci_compute_statistics(2, "poibin_exact", observed = c(0, 1),
                                       probabilities = c(.1, .2)), "one numeric vector per")
  expect_error(.pci_compute_statistics(1, "wald", estimate = 1, se = 0),
               "strictly positive")
  expect_error(.pci_compute_statistics(1, "wald", estimate = 1, se = 1, df = -1),
               "df must be positive")
  expect_error(.pci_compute_statistics(1, "wald", estimate = 0, se = 1,
                                       transform = "logit"), "in \\(0, 1\\)")
  expect_error(.pci_compute_statistics(1, "resampling", observed = 1,
                                       eta = c(0, 1), re_var = c(0, -1)),
               "nonnegative")
  expect_error(.pci_compute_statistics(1, "resampling", observed = 1,
                                       eta = c(0, 1), re_mean = 0), "matching")
  expect_error(.pci_compute_statistics(1, "poibin_exact", observed = 1,
                                       probabilities = c(.1, .2), reference = .3),
               "encode the null")
})
