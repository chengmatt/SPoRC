# Sum a prediction over the seasons a data source is fit against

A data source set to `"aggSeas"` has one observation for the year, so
the prediction it is compared against is the whole year's total. One set
to `"spltSeas"` is compared against its own season alone. Populations
are summed over either way, because a regional observation does not see
them separately.

## Usage

``` r
get_seas_pred(pred, r, y, seas, f, seas_agg)
```

## Arguments

- pred:

  Prediction array `[pop, region, year, season, fleet]`.

- r, y, seas, f:

  Region, year, season and fleet of the observation.

- seas_agg:

  Integer, `1` for a season total and `0` otherwise.

## Value

Scalar prediction.
