# Write one replicate's redrawn selectivity at length into the operating model

The curve goes wherever the operating model reads selectivity at length:
the length compositions selected at length, and, when growth is rebuilt,
the curve `take_sim_growth_years` reads through each replicate's own
keys for its selectivity at age. Without a growth rebuild the
selectivity at age is formed here, the curve read through the keys the
operating model was given, as the estimation model forms it.

## Usage

``` r
sim_sel_at_length(sim_env, type, sel_l, yrs, sim)
```

## Arguments

- sim_env:

  Simulation environment.

- type:

  `"fish"`, `"ret"` or `"srv"`.

- sel_l:

  Selectivity at length `[region, year, len, sex, fleet]`.

- yrs:

  Fitted years to write.

- sim:

  Replicate.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
