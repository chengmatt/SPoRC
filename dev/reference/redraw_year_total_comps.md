# Redraw the compositions a fleet reports once a year

A composition set to `"aggSeas"` is one draw a year from the expected
numbers summed over every season, so a season landing more fish counts
for more, the way the estimation model builds the prediction. The draw
goes in the season the fit holds it in
([`seas_agg_season`](https://chengmatt.github.io/SPoRC/dev/reference/seas_agg_season.md)),
read at that season's sample size, and the seasonal draws the cycle made
for the fleet are cleared.

## Usage

``` r
redraw_year_total_comps(sim_env, platform, y, sim)
```

## Arguments

- sim_env:

  Simulation environment, after the year's seasonal draws.

- platform:

  `"fish"` or `"srv"`.

- y, sim:

  Year and replicate.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
