test_that("correlated Gaussian fit recovers the null fingerprint", {
  set.seed(105)
  m <- 300
  size <- sample(20:300, m, replace = TRUE)
  x <- matrix(rnorm(m), ncol = 1)
  design <- cbind(1, x)
  theta <- c(0.03, 0.20)
  phi <- 0.02
  z <- rnorm(m, sqrt(size) * drop(design %*% theta),
             sqrt(1 + size * phi))

  fit <- fit_empirical_null(
    z, size, covariates = x, model = "correlated",
    p_grid = seq(0.7, 0.99, by = 0.02), phi_bounds = c(0, 0.2)
  )
  expect_s3_class(fit, "empinull_fit")
  expect_lt(max(abs(unname(fit$parameters$zeta) - theta)), 0.08)
  expect_lt(abs(fit$parameters$phi - phi), 0.02)
  expect_lt(abs(mean(fit$standardized)), 0.2)
  expect_equal(sd(fit$standardized), 1, tolerance = 0.2)
  expect_equal(nrow(predict(fit)), m)
  expect_equal(nrow(flag_empirical_null(fit)), m)
})

test_that("phi always denotes a variance component", {
  set.seed(106)
  m <- 250
  size <- sample(20:200, m, replace = TRUE)
  phi <- 0.04
  z <- rnorm(m, 0, sqrt(1 + size * phi))
  fit <- empirical_null_ind_RQ(
    z, size, p_grid = seq(0.8, 0.98, by = 0.02),
    phi_grid = seq(0, 0.15, by = 0.002)
  )
  expect_lt(abs(fit$phi - phi), 0.02)
  expect_equal(sd(z / sqrt(1 + size * fit$phi)), 1, tolerance = 0.2)
})

test_that("rank-deficient designs are rejected before fitting", {
  set.seed(107)
  size <- rep(50, 100)
  x <- rnorm(100)
  z <- rnorm(100)
  expect_error(
    fit_empirical_null(
      z, size, covariates = cbind(x, x), model = "correlated",
      p_grid = c(0.8, 0.9, 0.98), phi_bounds = c(0, 0.2)
    ),
    "rank deficient"
  )
})

test_that("nearly collinear designs trigger an audit warning", {
  set.seed(109)
  size <- sample(30:100, 140, replace = TRUE)
  x <- rnorm(140)
  z <- rnorm(140)
  fit <- fit_empirical_null(
    z, size, covariates = cbind(x, x + rnorm(140, sd = 1e-4)),
    model = "correlated", p_grid = c(0.8, 0.9, 0.98),
    phi_bounds = c(0, 0.2)
  )
  expect_true("ILL_CONDITIONED_DESIGN" %in% fit$diagnostics$warning_codes)
  expect_equal(audit_empirical_null(fit)$status, "warning")
})

test_that("invalid inputs fail before optimization", {
  expect_error(fit_empirical_null(c(0, 1), c(10, 0)), "positive")
  expect_error(fit_empirical_null(c(0, NA), c(10, 20)), "finite")
  expect_error(
    fit_empirical_null(c(0, 1), c(10, 20), family = "poisson"),
    "second_moment"
  )
})
