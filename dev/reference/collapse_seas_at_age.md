# Sum an at-age or at-length array over seasons into season one

The counterpart of
[`collapse_seas_obs`](https://chengmatt.github.io/SPoRC/dev/reference/collapse_seas_obs.md)
for compositions, which are drawn from the numbers behind them rather
than from a total. Fleets left at `"spltSeas"` are returned untouched.

## Usage

``` r
collapse_seas_at_age(arr, seas_agg, y, sim, n_seas)
```

## Arguments

- arr:

  Array `[pop, region, year, season, bin, sex, fleet, sim]` of predicted
  numbers.

- seas_agg:

  Integer vector, one per fleet.

- y, sim:

  Year and replicate being drawn.

- n_seas:

  Number of seasons.

## Value

`arr` with each aggregated fleet's year summed into season one.
