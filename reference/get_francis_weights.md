# Computes Francis weights, which is used internally by do_francis_reweighting

Computes Francis weights, which is used internally by
do_francis_reweighting

## Usage

``` r
get_francis_weights(
  n_regions,
  n_sexes,
  n_fleets,
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

- Use:

  Array from data list that specifies whether to use data that year

- ISS:

  Input sample size array

- Pred_array:

  Predicted values array dimensioned by n_regions, n_years, n_ages,
  n_sexes, n_fleets

- Obs_array:

  Observed values array dimensioned by n_regions, n_years, n_ages,
  n_sexes, n_fleets

- weights:

  Array of francis weights (NAs) to apply dimensioned by n_regions,
  n_years, n_sexes, n_fleets

- bins:

  Vector of bins used (age or length)

- comp_type:

  Matrix of composition structure types dimensioned by year and fleet

## Value

List of values for calculated francis weight, and a dataframe of
observed and expected means
