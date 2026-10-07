# Collapse a year's seasonal observations into one annual observation

A data source set to `"aggSeas"` reports once a year rather than once a
season, so the operating model has to draw one observation from the
year's total rather than eleven from its parts. The total is written
into the season the fit holds it in, season one unless `slot` says
otherwise, and the other seasons are left at zero.

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
  pop = FALSE,
  bias_correct_oe = 0,
  slot = rep(1, n_fleets)
)
```

## Arguments

- true_arr, obs_arr:

  Arrays `[region, year, season, fleet, sim]`, or with a leading
  population dim, of the true quantity and the observation.

- se_arr:

  Standard deviation array matching `true_arr` without the simulation
  dim, read in the season the total is drawn into.

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

- bias_correct_oe:

  Whether lognormal draws sit at their mean.

- slot:

  Integer vector, one per fleet, of the season the total is drawn into.
  From
  [`seas_agg_season`](https://chengmatt.github.io/SPoRC/dev/reference/seas_agg_season.md).

## Value

A list with the updated `true` and `obs` arrays.

## Details

The true quantity is summed first and the observation error is applied
once to that sum, so the error is on the annual total the way the
observation is reported, rather than on each season and then added up.
