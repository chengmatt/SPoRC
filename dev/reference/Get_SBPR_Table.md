# Spawning biomass per recruit, by origin and destination region

The age-and-region projection behind
[`Get_Det_Recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Det_Recruitment.md)'s
Beverton-Holt and Ricker curves: one recruit per origin region is
projected through every age, season and the plus-group, unfished and
fished. This depends only on mortality, selectivity and movement, never
on `R0`, spawning biomass or the current year, so a caller that holds
those fixed across years (as the annual cycle does; see `det_rec_args`
in
[`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md))
should build this once and pass it back into `Get_Det_Recruitment` as
`sbpr_table`, rather than paying for this projection again on every year
it is otherwise identical.

## Usage

``` r
Get_SBPR_Table(
  rec_dd,
  rec_region_prop,
  rec_seas_prop,
  n_pop,
  n_regions,
  n_ages,
  n_fish_fleets,
  WAA,
  MatAA,
  natmort,
  Movement,
  sgl_seas_spawning_movement,
  do_recruits_move,
  t_spawn,
  init_F,
  dmr,
  fish_sel,
  ret_sel,
  n_seas,
  spawn_seas,
  seasdur,
  sexratio_f,
  Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0
)
```

## Arguments

- rec_dd:

  Integer. 0 = density dependence within each population or region, 1 =
  shared across regions, valid only when `n_pop = 1`.

- rec_region_prop:

  Matrix (`n_pop × n_regions`) of the proportion of recruitment
  allocated to each region.

- rec_seas_prop:

  Matrix (`n_pop × n_seas`) of seasonal recruitment proportions. Must be
  zero before `spawn_seas` when `rec_lag = 0`.

- n_pop:

  Number of populations.

- n_regions:

  Number of spatial regions.

- n_ages:

  Number of age classes (including the plus group).

- n_fish_fleets:

  Integer. Number of fishery fleets.

- WAA:

  Array (`n_pop × n_regions × n_seas × n_ages`) of weight-at-age.

- MatAA:

  Array (`n_pop × n_regions × n_seas × n_ages`) of maturity-at-age.

- natmort:

  Array (`n_pop × n_regions × n_seas × n_ages`) of natural mortality, a
  rate per year in each season.

- Movement:

  Array (`n_pop × origin × destination × n_seas × n_ages`) of seasonal
  movement probabilities.

- sgl_seas_spawning_movement:

  Array (`n_pop × origin × destination × n_ages`) of spawning movement
  when a single season is used and `n_pop > 1`.

- do_recruits_move:

  Indicator for whether recruits move in their first year.

- t_spawn:

  Fraction of the spawning season that occurs before spawning.

- init_F:

  Array (`n_regions × n_seas × n_fish_fleets`) of initial fishing
  mortality.

- dmr:

  Array (`n_regions × n_seas × n_fish_fleets`) of initial (first year)
  discard mortality.

- fish_sel:

  Array (`n_pop × n_regions × n_seas × n_ages × n_fish_fleets`) of total
  fishery selectivity.

- ret_sel:

  Array (`n_pop × n_regions × n_seas × n_ages × n_fish_fleets`) of
  retained fishery selectivity.

- n_seas:

  Number of seasons per year.

- spawn_seas:

  Season index in which spawning occurs.

- seasdur:

  Numeric vector (`n_seas`) of seasonal durations as fractions of a
  year.

- sexratio_f:

  Matrix (`n_pop × n_regions`) of female recruitment proportions.

## Value

List of `phi0` and `phiF`, the unfished and fished spawning biomass per
recruit at `R0 = 1`: a `[pop x region]` matrix under `rec_dd = 0`, a
scalar under `rec_dd = 1`. `Get_Det_Recruitment` recovers `S0`/`SF` by
scaling these by `R0`.
