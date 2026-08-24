# GOA rex sole bridge data (Model 25.1)

The 2025 Gulf of Alaska rex sole assessment (Model 25.1: two areas with
growth estimated separately in each, survey conditional age-at-length,
an ageing error matrix) parsed from its input files, with the
assessment's report quantities attached under `$ss3` and its maximum
likelihood estimates under `$mle`. Built by
`dev/rex_bridge/R/build_rex_data.R` from Carey McGilliard's goa_rex
repository. Used by `tests/testthat/helper-bridge_goa_rex.R` and the rex
sole regression test, which evaluates SPoRC at the assessment's own
estimate.

## Usage

``` r
mlt_rg_goa_rex_data
```

## Format

A named list of observation arrays, biological inputs, the assessment's
configuration, `$mle` (parameter values) and `$ss3` (report quantities:
numbers at age, spawning biomass, recruitment, total biomass, catch,
indices, growth tables, age-length keys, selectivity, one block of
conditional age-at-length expected values, likelihood components).

## Source

<https://github.com/careymcgilliard/goa_rex>
