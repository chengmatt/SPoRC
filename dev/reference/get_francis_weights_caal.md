# Francis weights for conditional age-at-length compositions

A conditional age-at-length row is the age composition of the fish aged
from one length bin, so a fleet's samples are the (year, season, length
bin) rows rather than the years alone. The Francis statistic is the same
standardized mean-age residual used for the marginal compositions,
pooled over every row a fleet carries, which gives one weight per fleet
(per region and sex where the composition type splits them). A weight
per length bin would rest on far fewer samples and would not be stable.

## Usage

``` r
get_francis_weights_caal(
  n_regions,
  n_sexes,
  n_fleets,
  n_years,
  n_seas,
  n_lens,
  Use,
  ISS,
  Pred_array,
  Obs_array,
  weights,
  bins,
  comp_type
)
```

## Arguments

- n_regions:

  Number of regions

- n_sexes:

  Number of sexes

- n_fleets:

  Number of fleets (fishery or survey)

- n_years:

  Number of years

- n_seas:

  Number of seasons

- n_lens:

  Number of length bins

- Use:

  Use flags `[region, year, season, len, fleet]`

- ISS:

  Input sample sizes `[region, year, season, len, sex, fleet]`

- Pred_array:

  Predicted proportions `[region, year, season, len, bin, sex, fleet]`

- Obs_array:

  Observed proportions, same shape as `Pred_array`

- weights:

  Array of weights to fill, shaped like `ISS`

- bins:

  Vector of age bins the compositions are over

- comp_type:

  Matrix of composition types `[year, fleet]`

## Value

List with the filled `weights` and a data frame of the observed and
expected mean ages behind them
