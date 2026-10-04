# Convert Catch Advice To Fishing Mortality In The Operating Model

Finds the F by region, season and fleet that makes the operating model's
retained catch in year `y` equal `target`, replaying the year one season
at a time with the operating model's own dynamics. Fleets fishing the
same region share total mortality, what earlier seasons caught is gone
before later ones are fished, fish move between regions as the operating
model moves them, and under `rec_lag = 0` this year's recruits come from
the spawning biomass the trial F leaves.

## Usage

``` r
catch_to_F_om(
  target,
  y,
  sim,
  sim_env,
  target_units = "biom",
  catch_f_max = 5,
  catch_tol = 1e-06,
  catch_max_iter = 100
)
```

## Arguments

- target:

  Numeric array `[n_regions, n_seas, n_fish_fleets]`. Retained catch for
  year `y`, in `target_units`. Zero means no fishing.

- y:

  Integer. Year being fished. Its start of year numbers at age must
  exist, so call after `run_annual_cycle(y - 1, ...)`.

- sim:

  Integer. Simulation replicate.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).
  Under cohort growth, year `y`'s growth is formed here, from the same
  numbers the annual cycle forms it from at the start of the year.

- target_units:

  Units of `target` by fleet, `"biom"` (default, through `WAA_fish`) or
  `"abd"` for numbers, or 1 and 0. One value is used for every fleet.
  Independent of the fleet's `catch_units`.

- catch_f_max, catch_tol, catch_max_iter:

  Upper bound on F, relative catch tolerance and iterations per solve,
  as in
  [`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md).

## Value

Named list with `Fmort` `[n_regions, n_seas, n_fish_fleets]`, to write
into `sim_env$Fmort[, y, , , sim]`, and `resid`, the relative catch miss
in each cell.

## See also

Other Closed Loop Simulations:
[`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/dev/reference/condition_closed_loop_simulations.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/dev/reference/get_closed_loop_reference_points.md)
