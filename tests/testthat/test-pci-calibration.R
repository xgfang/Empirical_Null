test_that("fixed and theoretical nulls preserve provider-specific parameters", {
  fit <- .pci_calibrate(c(-2, 0, 3), "fixed", null_sd = 1.81)
  expect_equal(fit$null_mean, rep(0, 3))
  expect_equal(fit$null_sd, rep(1.81, 3))
  expect_null(fit$fit)
  varying <- .pci_calibrate(c(-2, 0, 3), "fixed",
                            null_mean = c(0, 1, 2), null_sd = c(1, 2, 3))
  expect_equal(varying$null_mean, c(0, 1, 2))
  expect_equal(varying$null_sd, c(1, 2, 3))
  expect_equal(.pci_calibrate(c(0, NA_real_), "theoretical")$null_sd, c(1, 1))
  expect_equal(.pci_calibrate(1:3, "theoretical", null_mean = rep(0, 3),
                             null_sd = rep(1, 3))$null_sd, rep(1, 3))
  expect_error(.pci_calibrate(1:3, "fixed", null_sd = 0), "positive")
  expect_error(.pci_calibrate(1:3, "fixed", null_mean = c(1, 2)), "length")
  expect_error(.pci_calibrate(1:3, "theoretical", null_sd = 1.81), "only for fixed")
})

test_that("robust overall calibration matches direct MASS estimation", {
  z <- c(seq(-2.2, 2.1, length.out = 35), 8, 12)
  for (psi in c("huber", "bisquare")) {
    psi_fun <- if (psi == "huber") MASS::psi.huber else MASS::psi.bisquare
    direct <- MASS::rlm(z ~ 1, psi = psi_fun, maxit = 1000, acc = 1e-8)
    fit <- .pci_calibrate(z, "overall", psi = psi)
    expect_equal(fit$null_mean, rep(unname(coef(direct)[1]), length(z)))
    expect_equal(fit$null_sd, rep(direct$s, length(z)))
    expect_equal(fit$diagnostics$psi, psi)
    expect_equal(fit$diagnostics$parameter_uncertainty, "not propagated")
    expanded <- .pci_calibrate(z, "overall", psi = psi,
                               null_mean = rep(0, length(z)),
                               null_sd = rep(1, length(z)))
    expect_equal(expanded$null_mean, fit$null_mean)
    expect_equal(expanded$null_sd, fit$null_sd)
  }
})

test_that("quantile groups reproduce existing groupwise estimates", {
  size <- seq_len(48)
  z <- sin(size) * 2 + rep(c(0, 0.2, -0.3, 0.5), each = 12)
  old <- empirical_null_groupwise(z, size)
  fit <- .pci_calibrate(z, "groupwise", size = size)
  expect_equal(fit$group, old$group)
  expect_equal(fit$null_mean, unname(old$intercept[old$group]))
  expect_equal(fit$null_sd, unname(old$scale[old$group]))
})

test_that("rank groups use stable sorting and ceil-sized blocks", {
  z <- c(-2, 0, 1, 0, 1, -1, 1, 3, -1, 0)
  size <- c(2, 1, 1, 3, 4, 4, 5, 5, 6, 7)
  fit <- .pci_calibrate(
    z, "groupwise", size = size, n_groups = 3, grouping = "rank",
    control = list(min_group_size = 2)
  )
  expect_equal(fit$group, c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L))
  expect_equal(fit$parameters$n_fitted, c(4L, 4L, 2L))
})

test_that("explicit string labels and fitting exclusions retain row alignment", {
  z <- c(-1, 0.2, 1, 50, -2, 0, 2, NA_real_)
  groups <- rep(c("small providers", "large providers"), each = 4)
  mask <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  fit <- .pci_calibrate(z, "groupwise", groups = groups, fit_mask = mask)
  expect_identical(fit$group, groups)
  expect_equal(fit$parameters$n_fitted, c(3L, 3L))
  expect_equal(fit$diagnostics$n_nonfinite, 1)
  expect_equal(fit$diagnostics$n_excluded, 1)
  direct <- MASS::rlm(z[1:3] ~ 1, psi = MASS::psi.bisquare,
                     maxit = 1000, acc = 1e-8)
  expect_equal(fit$null_mean[4], unname(coef(direct)[1]))
  expect_equal(fit$null_sd[4], direct$s)
  expect_true(is.finite(fit$null_sd[8]))
})

test_that("pooled and specified means leave fitted group scales intact", {
  z <- c(-3, -1, 0.5, 1, 2, 3.5, 6, 100)
  group <- rep(c("a", "b"), each = 4)
  mask <- c(rep(TRUE, 7), FALSE)
  original <- .pci_calibrate(z, "groupwise", groups = group, fit_mask = mask)
  pooled <- .pci_calibrate(z, "groupwise", groups = group, fit_mask = mask,
                           common_mean = TRUE)
  zero <- .pci_calibrate(z, "groupwise", groups = group, fit_mask = mask,
                         common_mean = 0)
  expect_equal(pooled$null_mean, rep(mean(z[mask]), length(z)))
  expect_equal(pooled$null_sd, original$null_sd)
  expect_equal(zero$null_mean, rep(0, length(z)))
  expect_equal(zero$null_sd, original$null_sd)
})

test_that("one group is a robust overall fit and tied quantiles are explicit", {
  z <- c(-2, -1, 0, 1, 2, 4)
  one <- .pci_calibrate(z, "groupwise", n_groups = 1)
  overall <- .pci_calibrate(z, "overall")
  expect_equal(one$null_mean, overall$null_mean)
  expect_equal(one$null_sd, overall$null_sd)
  expect_error(.pci_calibrate(z, "groupwise", size = rep(10, 6), n_groups = 2),
               "tied quantile breaks")
})

test_that("small groups and degenerate scales cannot silently calibrate", {
  z <- c(-2, 0, 3, 1)
  g <- c("large", "large", "large", "small")
  expect_error(.pci_calibrate(z, "groupwise", groups = g), "eligible providers")
  fallback <- .pci_calibrate(z, "groupwise", groups = g,
                             control = list(small_group = "theoretical"))
  expect_equal(fallback$null_mean[4], 0)
  expect_equal(fallback$null_sd[4], 1)
  expect_true("SMALL_GROUP_THEORETICAL_FALLBACK" %in% fallback$warning_codes)
  expect_true(fallback$parameters$fallback[2])
  expect_error(.pci_calibrate(rep(0, 8), "overall"), "nonpositive")
  expect_error(.pci_calibrate(c(1, 2, 3), "overall", fit_mask = rep(FALSE, 3)),
               "eligible providers")
})

test_that("robust tuning and convergence controls are honored", {
  z <- c(seq(-2, 2, length.out = 20), 10)
  fit <- .pci_calibrate(z, "overall", psi = "huber",
                        control = list(tuning = 1, maxit = 100, acc = 1e-7))
  direct <- MASS::rlm(z ~ 1, psi = MASS::psi.huber, k = 1,
                     maxit = 100, acc = 1e-7)
  expect_equal(fit$null_mean, rep(unname(coef(direct)[1]), length(z)))
  expect_equal(fit$null_sd, rep(direct$s, length(z)))
  slow <- .pci_calibrate(z, "overall", psi = "huber",
                         control = list(maxit = 1, acc = 1e-15))
  expect_true("ROBUST_DID_NOT_CONVERGE" %in% slow$warning_codes)
  expect_true(length(slow$diagnostics$warnings) > 0)
})

test_that("incompatible options fail before fitting", {
  z <- c(-1, 0, 2, 3)
  expect_error(.pci_calibrate(z, "overall", covariates = z), "covariates")
  expect_error(.pci_calibrate(z, "overall", second_moment = z), "second_moment")
  expect_error(.pci_calibrate(z, "overall", groups = 1:4), "groups")
  expect_error(.pci_calibrate(z, "fixed", common_mean = 0), "common_mean")
  expect_error(.pci_calibrate(z, "overall", control = list(maxitt = 10)),
               "Unsupported calibration")
  expect_error(.pci_calibrate(z, "overall", control = list(10)), "unique")
  expect_error(.pci_calibrate(z, "individualized", size = 1:4,
                              fit_mask = c(FALSE, TRUE, TRUE, TRUE)), "fit_mask")
  expect_error(.pci_calibrate(z, "individualized", size = 1:4,
                              control = list(common_intercept = 1)), "TRUE or FALSE")
  expect_error(.pci_calibrate(z, "correlated", size = 1:4,
                              control = list(include_intercept = FALSE)),
               "covariates are required")
  expect_error(.pci_calibrate(z, "moment", size = 1:4, second_moment = 1:4,
                              control = list(include_intercept = FALSE)),
               "covariates are required")
})

test_that("structured routes return the established empirical-null fit", {
  skip_if_not(is.loaded("_EmpiNull_negloglik_individualized"),
              "Compiled EmpiNull likelihood is unavailable")
  set.seed(92)
  z <- rnorm(80)
  size <- seq_len(80)
  fit <- .pci_calibrate(z, "individualized", size = size,
                        control = list(p_grid = c(0.8, 0.95), phi_bounds = c(0, 0.1)))
  expect_s3_class(fit$fit, "empinull_fit")
  expect_equal(fit$null_sd, sqrt(1 + size * fit$parameters$phi))
  expect_equal(fit$null_mean, rep(0, 80))
})

test_that("individualized common intercept retains fit-object compatibility", {
  skip_if_not(is.loaded("_EmpiNull_negloglik_individualized"),
              "Compiled EmpiNull likelihood is unavailable")
  set.seed(93)
  z <- rnorm(80, 0.5, 1)
  fit <- .pci_calibrate(
    z, "individualized", size = rep(20, 80),
    control = list(common_intercept = TRUE, p_grid = c(0.8, 0.95),
                   phi_bounds = c(0, 0.1))
  )
  expect_s3_class(fit$fit, "empinull_fit")
  expect_equal(fit$null_mean, rep(fit$parameters$mu, 80))
  expect_equal(standardize_empirical_null(fit$fit), (z - fit$null_mean) / fit$null_sd)
})
