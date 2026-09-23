#' Empirical-Null Calibration and Provider Inference
#'
#' EmpiNull estimates provider-level reference distributions and performs
#' post-model inference from supplied scores, estimates and standard errors,
#' or count summaries. Calibration can use a theoretical or specified null,
#' robust overall or groupwise estimates, or a structured empirical-null model.
#'
#' @details
#' Start with [fit_pci()] to construct a statistic, calibrate its null
#' distribution, and report p-values, directional flags, and supported
#' confidence intervals. Use [fit_empirical_null()] to fit a structured null
#' directly, then [standardize_empirical_null()], [flag_empirical_null()], and
#' [audit_empirical_null()] to inspect and use that fit.
#'
#' Models may allow the null mean to depend on provider summaries and the null
#' variance to depend on effective provider size. [poisson_null_moments()]
#' evaluates moments for the Poisson-lognormal working null. All parameters
#' named `phi` are variance components, whereas `null_sd` is a standard
#' deviation.
#'
#' The package uses outputs from an upstream risk-adjustment model. It does
#' not estimate patient-level regression coefficients or model covariance
#' matrices. Null parameters describe a provider-level distribution and do
#' not identify latent patient-level confounding or a causal decomposition.
#' Inspect diagnostics before interpreting provider flags.
#'
#' @references
#' Hartman, N., Messana, J. M., Kang, J., Naik, A. S., Shearon, T. H., and He, K.
#' (2024). Composite scores for transplant center evaluation: A new
#' individualized empirical null method. The Annals of Applied Statistics,
#' 18(1), 729--748. \doi{10.1214/23-AOAS1809}.
#'
#' @seealso [fit_pci()], [fit_empirical_null()], [poisson_null_moments()]
#' @keywords package
"_PACKAGE"
