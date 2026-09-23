test_that("Poisson Pearson moments are correctly centered", {
  x <- cbind(`(Intercept)` = 1, x = c(-1, 1))
  ans <- poisson_null_moments(
    Xbar = x, E = c(20, 40), E2 = c(30, 70),
    zeta = c(0, 0), sigmaAlpha2 = 0, sigmaEps2 = 0
  )
  expect_equal(ans$mean, c(0, 0))
  expect_equal(ans$variance, c(1, 1))
  expect_equal(ans$sd, c(1, 1))
})

test_that("patient-level second moments affect Poisson variance", {
  x <- matrix(1, nrow = 2, ncol = 1)
  ans <- poisson_null_moments(
    Xbar = x, E = c(50, 50), E2 = c(60, 200),
    zeta = 0, sigmaAlpha2 = 0.02, sigmaEps2 = 0.1
  )
  expect_equal(ans$mean[1], ans$mean[2])
  expect_gt(ans$variance[2], ans$variance[1])
})
