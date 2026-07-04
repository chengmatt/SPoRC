# Set data indicators to unused for specified years

Sets `Use*` binary indicator arrays to `0` for the specified years
across one or more data types, and optionally removes conventional tag
cohorts released in those years. Used in MSE closed-loop simulations and
retrospective analyses to exclude future or withheld data from the
estimation model without modifying the underlying observation arrays.
Only years present in `data$years` are affected; out-of-range values in
`unused_years` are silently ignored.

## Usage

``` r
set_data_indicator_unused(
  data,
  unused_years,
  what = c("Catch", "Catch_pop", "Discard", "Discard_pop", "FishIdx", "FishIdx_pop",
    "FishAgeComps", "FishAgeComps_pop", "FishLenComps", "FishLenComps_pop",
    "FishAgeComps_discard", "FishAgeComps_discard_pop", "FishLenComps_discard",
    "FishLenComps_discard_pop", "SrvIdx", "SrvIdx_pop", "SrvAgeComps", "SrvAgeComps_pop",
    "SrvLenComps", "SrvLenComps_pop", "conv_tagging")
)
```

## Arguments

- data:

  Named list of model data as constructed by the `Setup_Mod_*` family of
  functions.

- unused_years:

  Integer vector. Year indices (relative to `data$years`) to mark as
  unused. Values not present in `1:length(data$years)` are dropped.

- what:

  Character vector. Data types to modify. Any combination of:

  `"Catch"`

  :   Sets `UseCatch[, unused_years, , ] <- 0`.

  `"Catch_pop"`

  :   Sets `UseCatch_pop[, , unused_years, , ] <- 0`.

  `"Discard"`

  :   Sets `UseDiscard[, unused_years, , ] <- 0`.

  `"Discard_pop"`

  :   Sets `UseDiscard_pop[, , unused_years, , ] <- 0`.

  `"FishIdx"`

  :   Sets `UseFishIdx[, unused_years, , ] <- 0`.

  `"FishIdx_pop"`

  :   Sets `UseFishIdx_pop[, , unused_years, , ] <- 0`.

  `"FishAgeComps"`

  :   Sets `UseFishAgeComps[, unused_years, , ] <- 0`.

  `"FishAgeComps_pop"`

  :   Sets `UseFishAgeComps_pop[, , unused_years, , ] <- 0`.

  `"FishLenComps"`

  :   Sets `UseFishLenComps[, unused_years, , ] <- 0`.

  `"FishLenComps_pop"`

  :   Sets `UseFishLenComps_pop[, , unused_years, , ] <- 0`.

  `"FishAgeComps_discard"`

  :   Sets `UseFishAgeComps_discard[, unused_years, , ] <- 0`.

  `"FishAgeComps_discard_pop"`

  :   Sets `UseFishAgeComps_discard_pop[, , unused_years, , ] <- 0`.

  `"FishLenComps_discard"`

  :   Sets `UseFishLenComps_discard[, unused_years, , ] <- 0`.

  `"FishLenComps_discard_pop"`

  :   Sets `UseFishLenComps_discard_pop[, , unused_years, , ] <- 0`.

  `"SrvIdx"`

  :   Sets `UseSrvIdx[, unused_years, , ] <- 0`.

  `"SrvIdx_pop"`

  :   Sets `UseSrvIdx_pop[, , unused_years, , ] <- 0`.

  `"SrvAgeComps"`

  :   Sets `UseSrvAgeComps[, unused_years, , ] <- 0`.

  `"SrvAgeComps_pop"`

  :   Sets `UseSrvAgeComps_pop[, , unused_years, , ] <- 0`.

  `"SrvLenComps"`

  :   Sets `UseSrvLenComps[, unused_years, , ] <- 0`.

  `"SrvLenComps_pop"`

  :   Sets `UseSrvLenComps_pop[, , unused_years, , ] <- 0`.

  `"conv_tagging"`

  :   Removes tag cohorts whose release year falls in `unused_years`
      from `conv_tagged_fish`, `obs_conv_tag_fish_recap`,
      `conv_tag_release_indicator`, and updates `n_conv_tag_cohorts`.
      Only applied when `any(data$use_conv_fish_tagging == 1)`.

## Value

The modified `data` list with `Use*` indicators set to `0` for
`unused_years` and, if applicable, tagging cohorts from those years
removed.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md)
