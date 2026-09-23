# EmpiNull

**Empirical-null calibration and post-model inference for provider profiling in R.**

EmpiNull helps assess provider-level results when their null distribution may
have a shifted center or more variation than the theoretical reference allows.
It takes scores, estimates and standard errors, or count summaries from an
already fitted model and returns calibrated p-values, directional flags,
supported confidence intervals, and diagnostics.

The package provides two main entry points:

- **`fit_pci()`** combines statistic construction, null calibration, and inference
  in one call.
- **`fit_empirical_null()`** fits a structured null distribution directly when
  provider scores and effective sizes are already available.

EmpiNull is a standalone R package. It uses R and compiled C++ routines and does
not fit the underlying patient-level regression model. Supply the appropriate
model estimates, standard errors, expected counts, or null probabilities.

## Installation

EmpiNull requires **R 4.0 or later** and a C++ toolchain for installation from
source. On Windows, install the [Rtools version matching your R
version](https://cran.r-project.org/bin/windows/Rtools/). macOS requires the
Xcode command-line tools; Linux requires the standard R package build tools.

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("xgfang/Empirical_Null", upgrade = "never")
library(EmpiNull)
```

The repository is named `Empirical_Null`; the R package name is `EmpiNull`.
For installation from a downloaded or cloned copy, run this from its root:

```r
install.packages("remotes")
remotes::install_local(".", dependencies = TRUE, upgrade = "never")
library(EmpiNull)
```

Building the installed vignettes also requires
[Pandoc](https://pandoc.org/installing.html), which is included with RStudio:

```r
remotes::install_local(".", dependencies = TRUE, upgrade = "never",
                       build_vignettes = TRUE)
```

## Quick start

This complete example generates provider scores, estimates a robust null
within size groups, and reports adjusted inference.

```r
library(EmpiNull)

set.seed(42)
providers <- data.frame(
  id = sprintf("provider_%03d", 1:120),
  size = seq(20, 500, length.out = 120)
)
providers$z <- rnorm(120, mean = 0.15,
                     sd = sqrt(1 + 0.003 * providers$size))
providers$z[c(10, 100)] <- c(-5, 6)

pci <- fit_pci(
  z = providers$z,
  size = providers$size,
  id = providers$id,
  calibration = "groupwise",
  n_groups = 4,
  psi = "huber"
)

head(pci$results[, c("provider_id", "z_raw", "null_mean", "null_sd",
                     "z_adjusted", "p_value", "direction")])
pci$null_parameters
pci$diagnostics
```

For a raw normal score `z`, calibration uses
`z_adjusted = (z - null_mean) / null_sd`. The default test is two-sided at
`alpha = 0.05`. A direction of `-1` means significantly lower, `1` means
significantly higher, and `0` means not flagged. Lower and higher do not imply
better or worse care.

P-values are per-provider values. The package does not automatically adjust
for testing many providers; apply a separately chosen multiplicity procedure
when your analysis calls for one.

## Choose a statistic and a calibration

`statistic` defines the test constructed from your inputs. `calibration`
defines the null distribution used to interpret its normal-score
representation. The default is `statistic = "z"` with
`calibration = "theoretical"`.

| Input | `statistic` | What to supply |
|---|---|---|
| Normal scores | `"z"` | `z` |
| Estimates and standard errors | `"wald"` | `estimate`, `se`; optional `reference`, identity/log/logit `transform`, and `df` |
| Ratios with log-scale standard errors | `"log_ratio"` | `estimate`, `se`; default reference is 1 |
| Observed and expected counts | `"poisson_midp"`, `"poisson_exact"`, `"poisson_score"` | `observed`, `expected` |
| Binary outcomes | `"poibin"`, `"binary_score"` | Provider counts and a list of observation-level null probabilities |
| Simulated binary tests | `"resampling"` | Counts with null probabilities, or predictors and posterior effect summaries |
| Existing two-sided tests | `"p_value"` | `p_value` and signed `direction` |

| `calibration` | Null distribution or estimation method |
|---|---|
| `"theoretical"` | Preserve the original test's p-values, including discrete or t tails |
| `"fixed"` | User-specified normal mean and **standard deviation** |
| `"overall"` | Common robust location and scale; use `control = list(estimator = "likelihood")` for the likelihood estimator |
| `"groupwise"` | Robust location and scale within size groups or supplied groups |
| `"individualized"` | Gaussian null variance `1 + size * phi` |
| `"correlated"` | Covariate-dependent Gaussian mean and individualized variance |
| `"cre"` | Zero-centered Gaussian null with variance `size * phi` |
| `"moment"` | Poisson-lognormal moments used in a Gaussian working null; requires Pearson scores and `second_moment` |

Not every statistic/calibration combination is meaningful. In particular,
`"moment"` requires `"poisson_score"`, or explicitly identified precomputed
Pearson scores. See the [inference guide](vignettes/fit_pci.Rmd) for supported
combinations, controls, and interval definitions.

## Common workflows

### Observed versus expected events

```r
counts <- data.frame(
  id = c("A", "B", "C", "D"),
  observed = c(0, 8, 15, 30),
  expected = c(2, 10, 12, 18)
)

count_fit <- fit_pci(
  observed = counts$observed,
  expected = counts$expected,
  id = counts$id,
  statistic = "poisson_midp"
)
count_fit$results
```

The two-sided mid-p value gives half the observed probability mass to each
tail, then doubles the smaller tail. `"poisson_exact"` instead uses inclusive
tails with the package's documented directional convention. A mid-p test does
not have the same finite-sample conservativeness guarantee as an exact test.

### A specified null spread or a rate benchmark

```r
# An externally chosen null standard deviation.
fixed_fit <- fit_pci(z = c(-3, 0, 3), calibration = "fixed", null_sd = 1.81)
fixed_fit$results

# Test adjusted rates against a 20% benchmark.
rate_fit <- fit_pci(
  estimate = c(0.12, 0.19, 0.30),
  se = c(0.025, 0.03, 0.04),
  statistic = "wald",
  transform = "logit",
  reference = 0.20,
  alternative = "less"
)
rate_fit$results
```

`null_sd = 1.81` is an illustrative calibration choice, not a universal default
or a rate threshold. `reference` specifies the outcome benchmark on its
original scale. For Wald inference, standard errors are on the original scale
unless `se_scale = "transformed"`. The `"log_ratio"` shortcut instead defaults
to log-scale standard errors.

### Fit a structured empirical null directly

```r
# Reuse the synthetic provider scores and sizes from the quick start.
null_fit <- fit_empirical_null(
  z = providers$z,
  size = providers$size,
  family = "gaussian",
  model = "individualized",
  p_grid = c(0.8, 0.9, 0.95),
  phi_bounds = c(0, 0.05)
)

head(predict(null_fit))
head(flag_empirical_null(null_fit))
audit_empirical_null(null_fit)
```

The grid and bounds above are illustrative choices for these synthetic data.
The [empirical-null tutorial](vignettes/tutorial.Rmd) explains model assumptions,
the variance-component scale, effective sizes, and Poisson moment inputs.

## Results and interpretation

`fit_pci()` returns a `fit_pci` object. Its main components are:

| Component | Contents |
|---|---|
| `results` | Provider IDs, raw/adjusted scores and p-values, null parameters, flags, and available intervals |
| `null_fit`, `null_parameters` | Fitted calibration object and parameter summaries |
| `diagnostics` | Missing-input, convergence, identifiability, score-floor, and interval diagnostics as applicable |
| `method`, `call` | Selected methods and the original call |

Use `as.data.frame(pci)` for the results table. `na_action = "omit"` keeps
incomplete providers in their original rows with missing results. Inspect
diagnostics before reporting classifications.

Confidence intervals depend on the supplied information. Estimate/SE inputs
support transformed Wald intervals; counts support appropriate Poisson or
binary test inversion. Supplied Z/p-values and resampling alone do not define
an original-scale interval. Binary inversion reports a **log-odds shift** from
the supplied null probabilities. Estimated null parameters are held fixed in
the reported inference; their estimation uncertainty is not propagated.

Empirical-null parameters describe a provider-level reference distribution.
They do not identify an unmeasured patient-level confounder or establish a
causal explanation for differences between providers. Provider comparisons
also depend on the validity of the upstream risk-adjustment model.

## Documentation

- [Post-model inference guide](vignettes/fit_pci.Rmd): input formats, tests,
  calibration choices, confidence intervals, and diagnostics.
- [Empirical-null tutorial](vignettes/tutorial.Rmd): structured null models,
  assumptions, and reproducible examples.
- R help: `?fit_pci`, `?fit_empirical_null`, and `?poisson_null_moments`.
- Installed articles, when built: `browseVignettes("EmpiNull")`.
- [Release notes](NEWS.md) and [contributor instructions](CONTRIBUTING.md).

Report reproducible problems through the
[GitHub issue tracker](https://github.com/xgfang/Empirical_Null/issues).

## Citation and license

Use `citation("EmpiNull")` for the package citation. Background on individualized
empirical-null methods is provided by Hartman, N., Messana, J. M., Kang, J.,
Naik, A. S., Shearon, T. H., and He, K. (2024), *Composite scores for transplant
center evaluation: A new individualized empirical null method*, **The Annals
of Applied Statistics**, 18(1), 729–748.
[doi:10.1214/23-AOAS1809](https://doi.org/10.1214/23-AOAS1809).

EmpiNull is distributed under the **GNU General Public License, version 3
(GPL-3)**; see [the license text](LICENSE.md). Authorship and maintainer details
are recorded in [DESCRIPTION](DESCRIPTION).
