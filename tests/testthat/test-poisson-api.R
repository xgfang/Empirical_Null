test_that("Poisson moment fit returns a usable agent object", {
  set.seed(108)
  m <- 180
  E <- runif(m, 20, 200)
  ratio <- runif(m, 1.2, 5)
  E2 <- E * ratio
  x <- matrix(rnorm(m, sd = 0.25), ncol = 1)
  design <- cbind(1, x)
  truth <- poisson_null_moments(design, E, E2, c(-0.02, 0.05), 0.01, 0.03)
  z <- rnorm(m, truth$mean, truth$sd)

  fit <- fit_empirical_null(
    z, E, covariates = x, second_moment = E2,
    family = "poisson", model = "moment",
    p_grid = seq(0.75, 0.99, by = 0.03)
  )
  expect_s3_class(fit, "empinull_fit")
  expect_true(all(is.finite(fit$null_mean)))
  expect_true(all(is.finite(fit$null_sd)))
  expect_true(all(fit$null_sd > 0))
  expect_true(fit$parameters$sigmaAlpha2 >= 0)
  expect_true(fit$parameters$sigmaEps2 >= 0)
})
