# Sum a population-specific prediction over the seasons a data source is fit against

The population-specific counterpart of
[`get_seas_pred`](https://chengmatt.github.io/SPoRC/dev/reference/get_seas_pred.md),
which keeps the population it is given rather than summing over them.

## Usage

``` r
get_seas_pred_pop(pred, p, r, y, seas, f, seas_agg)
```

## Arguments

- pred:

  Prediction array `[pop, region, year, season, fleet]`.

- p, r, y, seas, f:

  Population, region, year, season and fleet of the observation.

- seas_agg:

  Integer, `1` for a season total and `0` otherwise.

## Value

Scalar prediction.
