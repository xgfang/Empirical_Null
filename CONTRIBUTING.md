# Contributing to EmpiNull

EmpiNull is an R package with native C++ likelihood routines. Development needs
R, an R-compatible C++ toolchain, and the package dependencies listed in
`DESCRIPTION`. Windows contributors should use the Rtools version matching
their R installation. Pandoc is needed to render the articles.

## Set up a checkout

```r
install.packages(c("remotes", "pkgbuild", "pkgload", "testthat",
                   "roxygen2", "rcmdcheck", "pkgdown"))
remotes::install_deps(dependencies = TRUE, upgrade = "never")
pkgbuild::check_build_tools(debug = TRUE)
```

## Test and document changes

Run these commands from the package root:

```r
pkgload::load_all()
testthat::test_local()
roxygen2::roxygenise()
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"),
                    error_on = "warning")
pkgdown::build_site()
```

The full checks load the compiled package. Sourcing R files alone is useful
for individual calculations but does not verify the compiled likelihood
routes. Do not describe skipped checks as passing.

Edit help text in `R/` and regenerate `man/` and `NAMESPACE` with roxygen2.
Edit articles in `vignettes/`. Generated HTML in `docs/` and `vignettes/`,
compiled libraries, personal session files, and local data should not be
committed. The GitHub check workflow builds documentation and runs tests from
the source package.

## Statistical changes

For a new statistic or calibration method, document its null hypothesis,
input scale, tail convention, assumptions, and interval interpretation. Add
numerical checks against independently derived formulas and relevant boundary
cases. Preserve provider IDs and row alignment, expose numerical failures in
diagnostics, and distinguish variance components from standard deviations.

Use synthetic data in examples and bug reports. Include the smallest
reproducible example, the observed and expected behavior, and `sessionInfo()`.
Keep comparisons with external implementations in validation tests with their
provenance; user-facing documentation should explain EmpiNull directly.
