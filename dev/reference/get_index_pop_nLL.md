# Population-specific index likelihoods on the arithmetic scale

The population-specific counterpart of
[`get_index_regional_nLL`](https://chengmatt.github.io/SPoRC/dev/reference/get_index_regional_nLL.md),
for normal fleets only. A multivariate normal covariance describes the
regional series, so those fleets are fit on the log scale here and go
through
[`eval_index_osa_nLL`](https://chengmatt.github.io/SPoRC/dev/reference/eval_index_osa_nLL.md)
with the lognormal ones.

## Usage

``` r
get_index_pop_nLL(nLL_arr, Use, Obs, Pred, SD, LikeType, seas_Type, n_fleets)
```

## Arguments

- nLL_arr:

  Container `[pop, region, year, season, fleet]`.

- Use:

  Indicator array of the cells that are fit.

- Obs:

  Observed index, untransformed.

- Pred:

  Predicted index `[pop, region, year, season, fleet]`.

- SD:

  Index standard deviation.

- LikeType:

  Integer vector, one likelihood code per fleet.

- seas_Type:

  Integer vector, whether each fleet is fit as a season total.

- n_fleets:

  Number of fleets.

## Value

`nLL_arr` with the fitted cells filled.
