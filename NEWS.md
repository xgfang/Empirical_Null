# EmpiNull 1.2.0

- Reworked the README, package help, and tutorials as self-contained R documentation with runnable synthetic examples and explicit statistical conventions.
- Fixed a named-scalar indexing error in individualized common-intercept fitting, exposed by the compiled integration tests.
- Added `fit_pci()` to combine post-model statistic construction, null calibration, p-values, flags, and supported confidence intervals.
- Added theoretical and fixed nulls, including a user-specified SD of 1.81; robust overall and groupwise Huber/bisquare calibration; explicit/rank/quantile grouping; fitting exclusions; and common-mean overrides.
- Connected existing likelihood-based overall, individualized, correlated, CRE, and Poisson moment nulls, including the individualized shared-mean option.
- Added identity/log/logit Wald inference with normal or Student reference, supplied-p conversion, Poisson mid-p/exact/Pearson tests, deterministic Poisson-binomial and binary-score tests, and reproducible Bernoulli/posterior resampling.
- Added one-sided inference and mean-adjusted transformed intervals, classical Poisson intervals and conditional test inversion. Undefined, empty, nonmonotone, or unbracketed intervals are reported explicitly rather than assigned misleading endpoints.
- Added provider-ID and omitted-row alignment, inference metadata, diagnostic codes, S3 display/data-frame conversion, and focused numerical and integration tests.
- Existing `fit_empirical_null()` calls retain their prior behavior. The new wrapper uses robust estimation by default for `calibration="overall"`; select `control=list(estimator="likelihood")` for the existing overall likelihood model.

# EmpiNull 1.1.0

- Added the canonical `fit_empirical_null()` interface for guided and automated workflows.
- Added provider-level predictions, score standardization, flagging, and structured audit output.
- Standardized `phi` to mean a variance component in every Gaussian function.
- Corrected the Poisson Pearson-score mean by restoring the required minus-one centering term.
- Corrected censored likelihood probabilities to use provider-specific candidate null means.
- Renamed the Poisson implementation conceptually as an exact-moment, moment-matched Gaussian working null.
- Added design-rank, condition-number, convergence, boundary, central-skewness, subgroup-calibration, and Poisson variance-identification diagnostics.
- Added unit tests for formula centering, second-moment sensitivity, parameter scale, recovery, validation, and agent output.
