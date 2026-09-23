# EmpiNull agent contract

Use `fit_pci()` for unified post-model provider inference, or
`fit_empirical_null()` for direct empirical-null fitting. Do not call low-level
likelihood functions directly in guided analyses. Existing model-selection
rules below apply to the direct empirical-null fitting API.

## Post-model inference

- Establish what the supplied scores, estimates/SEs, or counts represent before choosing a statistic.
- Numeric estimate references are on the original scale. `null_sd=1.81` is a fixed null standard deviation, not a performance threshold.
- Keep Pearson Poisson scores separate from Poisson mid-p quantile scores. Moment calibration requires Pearson scores and second expected-count moments.
- Inspect `pci$diagnostics` and retain `pci$results`, `pci$method`, `pci$call`, and `pci$package_version` before reporting flags.
- Report missing confidence intervals and their diagnostic codes, including numerical search limits, empty acceptance sets, and detected nonmonotone score inversion. Do not replace them with guessed bounds.
- Numeric directions indicate higher/lower outcomes, not clinical desirability.
- Supplied Z/p-values and simulation tests alone do not determine original-scale confidence intervals. Model covariance, likelihood-ratio refits, and standardized-measure construction remain responsibilities of the fitted-model workflow.

## Required inputs

- `z`: one finite provider-level test statistic per provider
- `size`: one positive effective sample size per provider
- `family`: `"gaussian"` or `"poisson"`
- `covariates`: optional provider-level summaries with one row per provider
- `second_moment`: required for Poisson fits and equal to the provider-level second exposure moment

The agent must ask the user what each score, size measure, and provider summary represents. It must not invent, transform, or select covariates without reporting that choice.

## Model selection

- For Gaussian scores without provider summaries, use `model = "auto"`. This selects the individualized variance model.
- For Gaussian scores with provider summaries, use `model = "auto"`. This selects the correlated mean and individualized variance model.
- For Poisson scores, use `model = "auto"`. This selects the exact-moment Gaussian working null and requires `second_moment`.
- Use `model = "overall"` only when the user explicitly wants a common empirical null.

## Required outputs

Always return or store all of the following:

1. `predict(fit)` for provider-specific null means and standard deviations
2. `flag_empirical_null(fit)` for corrected scores, two-sided p-values, and flag directions
3. `audit_empirical_null(fit)` for assumptions, diagnostics, and interpretation limits
4. The original function call and package version for reproducibility

## Abstention and warnings

The agent must stop if fitting returns an error for invalid inputs or a rank-deficient design. It must surface every code in `fit$diagnostics$warning_codes`. A warning status is not permission to describe the correction as reliable.

The most important warning codes are:

- `OPTIMIZATION_FAILED`
- `P0_AT_SEARCH_BOUNDARY`
- `PHI_AT_SEARCH_BOUNDARY`
- `ILL_CONDITIONED_DESIGN`
- `FEW_PROVIDERS_PER_PARAMETER`
- `POISSON_VARIANCE_COMPONENTS_WEAKLY_IDENTIFIED`
- `CENTRAL_SKEWNESS`
- `SMALL_CENTRAL_SAMPLE`

## Interpretation boundary

The package estimates a provider-level null fingerprint that is useful for calibration and flagging. It does not identify the latent patient-level confounder, its distribution, or its causal decomposition. An agent must state this boundary whenever it explains results.
