test_that("Poisson inference agrees with the recorded Python reference", {
  ref <- read.csv(test_path("fixtures", "fit_pci_python_poisson.csv"))
  mid <- fit_pci(observed = ref$observed, expected = ref$expected,
                 statistic = "poisson_midp", ci = "none")$results
  exact <- fit_pci(observed = ref$observed, expected = ref$expected,
                   statistic = "poisson_exact")$results
  expect_equal(mid$z_raw, ref$z_midp, tolerance = 1e-10)
  expect_equal(mid$p_value, ref$p_midp, tolerance = 1e-10)
  expect_equal(exact$p_value, ref$p_exact, tolerance = 1e-10)
  expect_equal(exact$ci_lower, ref$lo_exact, tolerance = 1e-10)
  expect_equal(exact$ci_upper, ref$hi_exact, tolerance = 1e-10)
  expect_true(fit_pci(observed = 1, expected = 1,
                             statistic = "poisson_midp", ci = "none")$method$mid_p)
  expect_false(fit_pci(observed = 1, expected = 1,
                      statistic = "poisson_exact", ci = "none")$method$mid_p)
})

test_that("log-ratio inference agrees with the recorded Python reference", {
  ref <- read.csv(test_path("fixtures", "fit_pci_python_logratio.csv"))
  result <- fit_pci(estimate = ref$estimate, se = ref$se,
                    statistic = "log_ratio", calibration = "fixed",
                    null_mean = ref$mu, null_sd = ref$sd)$results
  expect_equal(result$z_raw, ref$z, tolerance = 1e-10)
  expect_equal(result$p_value, ref$p, tolerance = 1e-10)
  expect_equal(result$ci_lower, ref$lower, tolerance = 1e-10)
  expect_equal(result$ci_upper, ref$upper, tolerance = 1e-10)
})

test_that("Poisson-binomial inference agrees with the recorded Python reference", {
  ref <- read.csv(test_path("fixtures", "fit_pci_python_binary.csv"))
  probabilities <- lapply(strsplit(ref$probs, ";", fixed = TRUE), as.numeric)
  result <- fit_pci(observed = ref$observed, probabilities = probabilities,
                    statistic = "poibin", ci = "none")$results
  expect_equal(result$p_value, ref$p, tolerance = 1e-10)
})

test_that("Wald results retain the statistic and reference degrees of freedom", {
  result <- fit_pci(estimate = c(1, 2), se = c(0.5, 0.5),
                    statistic = "wald", df = c(5, Inf), ci = "none")$results
  expect_equal(result$statistic_value, c(2, 4))
  expect_equal(result$df, c(5, Inf))
})
