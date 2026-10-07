# Collapse a year's seasonal discards into one annual observation

The discard counterpart of
[`collapse_seas_obs`](https://chengmatt.github.io/SPoRC/dev/reference/collapse_seas_obs.md):
a fleet set to `"aggSeas"` gets one lognormal draw a year, from
[`year_discard_total`](https://chengmatt.github.io/SPoRC/dev/reference/year_discard_total.md),
written into the season the fit holds it in.

## Usage

``` r
collapse_seas_discards(sim_env, y, sim, bias_correct_oe = 0)
```

## Arguments

- sim_env:

  Simulation environment, after the year's discards are drawn season by
  season.

- y, sim:

  Year and replicate.

- bias_correct_oe:

  Whether lognormal draws sit at their mean.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
