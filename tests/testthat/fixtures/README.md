# Python reference fixtures

These small numerical fixtures were generated on 2026-09-23 by executing
UM-KevinHe/pprof_py at commit
`f83cc0793e233354c18deb2c7862150752f02c6f`.

Source: https://github.com/UM-KevinHe/pprof_py/tree/f83cc0793e233354c18deb2c7862150752f02c6f

- `fit_pci_python_poisson.csv`: survival `poisson_midp_zscore`, the implied
  two-sided normal tail, and `poisson_exact_test` p-values and ratio intervals,
  at alpha 0.05 and the default expected-count approximation threshold 100.
- `fit_pci_python_logratio.csv`: survival `log_ratio_zscore`, normal null
  calibration using each row's `mu` and `sd`, and
  `log_ratio_confidence_intervals` at alpha 0.05, including the legacy zero-ratio
  score approximation. `se` is the log-ratio standard error.
- `fit_pci_python_binary.csv`: `poibin_exact_pvalue` with the two-sided mid-p
  convention. The semicolon-separated `probs` column stores each provider's
  null Bernoulli probabilities.

The regression tests compare these selected examples to a tolerance of 1e-10.
They do not claim bitwise agreement for robust fitting, all input combinations,
or Monte Carlo draws. EmpiNull also preserves some boundary cases differently,
including exact zero/one Bernoulli probabilities and untruncated theoretical
tail probabilities.
