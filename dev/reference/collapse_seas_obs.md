# Collapse a year's seasonal observations into one annual observation

A data source set to `"aggSeas"` reports once a year rather than once a
season, so the operating model has to draw one observation from the
year's total rather than eleven from its parts. The total is written
into season one by convention and the other seasons are left at zero,
which is where the estimation model's `Use` array should mark it.

## Usage

``` r
collapse_seas_obs(
  true_arr,
  obs_arr,
  se_arr,
  seas_agg,
  like_type,
  y,
  sim,
  n_seas,
  n_regions,
  n_fleets,
  pop = FALSE
)
```

## Arguments

- true_arr, obs_arr:

  Arrays `[region, year, season, fleet, sim]`, or with a leading
  population dim, of the true quantity and the observation.

- se_arr:

  Standard deviation array matching `true_arr` without the simulation
  dim.

- seas_agg:

  Integer vector, one per fleet.

- like_type:

  Integer vector, one per fleet, passed to
  [`draw_index_obs`](https://chengmatt.github.io/SPoRC/dev/reference/draw_index_obs.md).

- y, sim:

  Year and replicate being drawn.

- n_seas, n_regions, n_fleets:

  Dimension sizes.

- pop:

  Logical, whether the arrays have a leading population dim.

## Value

A list with the updated `true` and `obs` arrays.

## Details

The true quantity is summed first and the observation error is applied
once to that sum, so the error is on the annual total the way the
observation is reported, rather than on each season and then added up.
