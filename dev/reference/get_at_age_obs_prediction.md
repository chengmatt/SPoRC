# Predicted values at the ages the observations are recorded on

An at-age observation counts fish by the age they were read as, so the
prediction at each observed age sums the predictions at every model age
read as it, weighted by the fleet's ageing error matrix, with units and
discard mortality applied at the model age first.

## Usage

``` r
get_at_age_obs_prediction(
  source,
  arrays,
  p_idx,
  r_idx,
  s_idx,
  y,
  seas,
  obs_ages,
  f,
  ageing_error = NULL
)
```

## Arguments

- source, arrays, p_idx, r_idx, s_idx, y, seas, f:

  See
  [`get_at_age_prediction`](https://chengmatt.github.io/SPoRC/dev/reference/get_at_age_prediction.md).

- obs_ages:

  Integer vector of the observed ages to predict.

- ageing_error:

  Array `[n_years, n_ages, n_obs_ages]` for this fleet, or `NULL` when
  the observed ages are the model ages.

## Value

A vector the length of `obs_ages`.
