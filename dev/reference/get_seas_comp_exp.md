# Predicted composition numbers for the seasons a data source is fit against

A composition set to `"aggSeas"` has one observation for the year, so
the numbers behind it are summed over every season before they are
turned into proportions. One set to `"spltSeas"` reads its own season
alone.

## Usage

``` r
get_seas_comp_exp(ExpArr, y, seas, f, seas_agg, p = NULL)
```

## Arguments

- ExpArr:

  Prediction array `[pop, region, year, season, bin, sex, fleet]`.

- y, seas, f:

  Year, season and fleet of the observation.

- seas_agg:

  Integer, `1` for a season total and `0` otherwise.

- p:

  Population index, or `NULL` for a regional data source, which sums
  over populations.

## Value

Array of predicted numbers with region, bin and sex left.
