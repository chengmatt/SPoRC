# Set Data Indicators to Unused for Specified Years

Set Data Indicators to Unused for Specified Years

## Usage

``` r
set_data_indicator_unused(
  data,
  unused_years,
  what = c("Catch", "FishIdx", "FishAgeComps", "FishLenComps", "SrvIdx", "SrvAgeComps",
    "SrvLenComps", "Tagging")
)
```

## Arguments

- data:

  Data list for RTMB model

- unused_years:

  Integer vector specifying which years to mark as unused. Only years
  present in `data$years` are considered.

- what:

  Character vector specifying which data types to modify. Possible
  values include:

  "Catch"

  :   Catch data indicators.

  "FishIdx"

  :   Fishery index data indicators.

  "FishAgeComps"

  :   Fishery age composition data indicators.

  "FishLenComps"

  :   Fishery length composition data indicators.

  "SrvIdx"

  :   Survey index data indicators.

  "SrvAgeComps"

  :   Survey age composition data indicators.

  "SrvLenComps"

  :   Survey length composition data indicators.

  "Tagging"

  :   Tagging data and associated cohorts.

## Value

The modified `data` object, with indicators set to 0 for the specified
years and tagging cohorts removed if relevant.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/reference/rho_trans.md)
