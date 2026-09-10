# Lognormal index likelihoods for observations registered for OSA residuals

Evaluated after
[`RTMB::OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html) has
substituted the observation vector, so the objective reads the same
values the residuals are computed from. The observations arrive already
flattened and on the log scale.

## Usage

``` r
eval_index_osa_nLL(nLL_arr, obs_vec, obs_map, Pred, SD, seas_Type, const, pop)
```

## Arguments

- nLL_arr:

  Container, regional or population-specific.

- obs_vec:

  Flattened observations, log scale, as returned by
  [`RTMB::OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

- obs_map:

  Index of each observation, from `arrayInd`.

- Pred:

  Predicted index `[pop, region, year, season, fleet]`.

- SD:

  Index standard deviation.

- seas_Type:

  Integer vector, whether each fleet is fit as a season total.

- const:

  Constant added inside the log.

- pop:

  Logical, whether the data source is population-specific.

## Value

`nLL_arr` with the fitted cells filled.
