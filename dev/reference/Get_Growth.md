# Growth module

Builds the size-age transition arrays and, when asked, weight at age and
maturity at age from the growth parameters, for every year. Every fleet
gets its own key and weight, read at that fleet's timing within the
season: each fishery fleet's at `t_fish`, each survey's at `t_srv`, and
the spawning weight at `t_spawn`.

## Usage

``` r
Get_Growth(
  ln_growth_pars,
  growth_A1,
  growth_A2,
  growth_L0,
  growth_len_lower,
  growth_cv_type,
  growth_sd_type,
  growth_dist,
  growth_plus_group,
  growth_L2_asymptote = 0,
  derive_waa,
  wt_len_pars,
  ages,
  seasdur,
  spawn_seas,
  t_spawn,
  n_pop,
  n_regions,
  n_yrs,
  n_seas,
  n_sexes,
  n_fish_fleets,
  n_srv_fleets,
  t_fish,
  t_srv,
  ln_growth_devs = NULL,
  growth_tv_model = NULL,
  growth_tv_link = 0,
  growth_par_bounds = NULL,
  growth_tv_type = 0,
  growth_cohort_styr = 1,
  years_eval = NULL,
  ln_growth_semipar_devs = NULL,
  growth_semipar = 0
)
```

## Arguments

- ln_growth_pars:

  Array `[pop, region, sex, n_gpars]` of log growth parameters in the
  order L1, L2, K, CV1, CV2 and, for the Richards form, rho.

- growth_A1, growth_A2:

  Reference ages; see
  [`get_laa_curve`](https://chengmatt.github.io/SPoRC/dev/reference/get_laa_curve.md).

- growth_L0:

  Length at age zero.

- growth_len_lower:

  Lower edges of the length bins.

- growth_cv_type, growth_sd_type, growth_dist:

  Integer switches passed to
  [`get_laa_curve`](https://chengmatt.github.io/SPoRC/dev/reference/get_laa_curve.md)
  and
  [`get_alk`](https://chengmatt.github.io/SPoRC/dev/reference/get_alk.md).

- growth_plus_group:

  Integer (0/1); apply
  [`plus_group_size`](https://chengmatt.github.io/SPoRC/dev/reference/plus_group_size.md).

- derive_waa:

  Integer (0/1); build weight at age from the key.

- wt_len_pars:

  Array `[pop, region, sex, 2]` of weight-length parameters \\a, b\\ in
  \\W = a L^b\\, read when `derive_waa = 1`.

- ages:

  Numeric vector of model ages.

- seasdur:

  Season durations as fractions of a year.

- spawn_seas, t_spawn:

  Spawning season and fraction of it elapsed at spawning.

- n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets, n_srv_fleets:

  Dimensions.

- t_fish:

  Array `[region, season, fish_fleet]` of fishery timings, the fraction
  of the season elapsed at which each fishery fleet's key and weight at
  age are read.

- t_srv:

  Array `[region, season, srv_fleet]` of survey timings, the fraction of
  the season elapsed when each survey is taken, at which that survey's
  key and weight at age are read.

- ln_growth_devs:

  Array `[pop, region, year, n_gpars, sex]` of time-varying deviations,
  or `NULL` for none.

- growth_tv_model:

  Integer vector `[n_gpars]`, 0 constant, 1 iid, 2 random walk, per
  parameter.

- growth_tv_link:

  Integer, 0 log link, 1 logit link within `growth_par_bounds`.

- growth_par_bounds:

  Matrix `[n_gpars x 2]` of bounds for the logit link.

- growth_tv_type:

  Integer, 0 each year on its own curve, 1 cohort propagation.

- growth_cohort_styr:

  Year index the cohort propagation starts from.

- years_eval:

  Integer vector of the years to evaluate; the default is every year.

- ln_growth_semipar_devs:

  Array `[pop, region, year, age, sex]` of log deviations on mean length
  at age, or `NULL` for none. Multiplies the parametric mean at age, so
  the curve stays the parametric part and the deviations have departures
  from it; the spread at age follows the deviated mean, leaving the
  coefficient of variation at age alone.

- growth_semipar:

  Integer process error code for those deviations, `0` for none. Only
  its being nonzero is read here; the structure is penalized in the
  objective.

## Value

List with `SizeAgeTrans_fish` and `SizeAgeTrans_srv`
`[pop, region, year, season, len, age, sex, fleet]`, one key per fleet
at that fleet's timing, the spawning key `SizeAgeTrans_spawn`
`[pop, region, year, len, age, sex]`, `mean_LAA_fish`, `sd_LAA_fish`,
`mean_LAA_srv`, `sd_LAA_srv`
`[pop, region, year, season, age, sex, fleet]`, `mean_LAA_spawn` and
`sd_LAA_spawn` `[pop, region, year, season, age, sex]`, `L_beg`
`[pop, region, year, age, sex]` (the start-of-year mean length), `Linf`
`[pop, region, year, sex]`, `growth_pars_y`
`[pop, region, year, 6, sex]`, and when `derive_waa = 1` also `WAA`,
`WAA_fish` and `WAA_srv`.

## Time variation

With no deviations the key is built once and broadcast over the years.
With deviations on any parameter and `growth_tv_type = 0` ("curve")
every year's sizes are read from that year's own curve. Under
`growth_tv_type = 1` ("cohort") the sizes are advanced cohort by cohort
from `growth_cohort_styr` on, which needs the numbers at age of each
year to blend the plus group; that form is run one year at a time from
the population dynamics through
[`Get_Growth_Year`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth_Year.md),
and this function only evaluates the years before the propagation
starts, which all sit on the first year's curve.

## Timing within the year

A season starts at the cumulative duration of the seasons before it ,
and a point inside a season is that start plus the fraction of the
season elapsed times the season's duration , so every evaluation is
anchored at the season start rather than compounded from the previous
evaluation. Mean length is the curve read at the real age, the integer
age plus that elapsed fraction of a year, and the plus group grows from
its adjusted size by the curve's increment over the same elapsed time.
Under constant growth the increment over a season is exactly the curve's
difference between the two real ages, so splitting the year into seasons
of any durations changes no mean length. The spawning weight is read at
the spawning fraction of the season, and each fleet's key and weight at
that fleet's own timing, `t_fish` or `t_srv`, so a composition and the
weight behind an index are formed at the point in the season the
observation is taken. Fleets that share a timing share one evaluation.
Seasonal multipliers on \\K\\ are not kept.
