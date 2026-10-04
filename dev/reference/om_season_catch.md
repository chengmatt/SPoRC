# Operating Model Catch Over One Season At A Trial F

Runs one season of the operating model from the numbers at its start,
with the same mortality, movement and catch equation as `apply_pop_dy`
and `generate_fishery_catch_comp_idx`. Only the rec_lag 0 recruitment
step writes into `sim_env`, values the annual cycle sets again before
reading.

## Usage

``` r
om_season_catch(
  N,
  F_seas,
  y,
  seas,
  sim,
  sim_env,
  target_units,
  spawn_recruits = FALSE
)
```

## Arguments

- N:

  Numeric array `[n_pop, n_regions, n_ages, n_sexes]`. Numbers at age at
  the start of the season, with recruits already known added.

- F_seas:

  Numeric matrix `[n_regions, n_fish_fleets]`. Trial F.

- y, seas, sim:

  Year, season and replicate.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).

- target_units:

  Numeric `[n_fish_fleets]`. 1 sums catch in biomass through `WAA_fish`,
  0 in numbers.

- spawn_recruits:

  Logical. Whether this year's recruits arrive in this season from its
  own spawning biomass (`rec_lag = 0` at `spawn_seas`).

## Value

Named list with `catch` `[n_regions, n_fish_fleets]` in `target_units`,
`N_end`, the survivors at the end of the season before ageing, shaped
like `N`, and `Rec_y` `[n_pop, n_regions]`, this year's recruitment when
`spawn_recruits`, `NULL` otherwise.
