# A deviation map run map_out over the operating model's years

Years past the fit are active in every column the fit varies, under
fresh levels that share what that column's last active fitted year
shares, as projection years are always active for movement. Years the
operating model does not run are cut.

## Usage

``` r
sim_map_over_years(map, n_yrs)
```

## Arguments

- map:

  Integer array `[pop, region, year, bin, sex]`, `NA` where a cell is
  fixed.

- n_yrs:

  Years the operating model runs.

## Value

The map with `n_yrs` years.
