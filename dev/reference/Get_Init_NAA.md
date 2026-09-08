# Initial Numbers-at-Age (NAA)

Numbers at age in the first model year, by population, region, age and
sex, at the equilibrium implied by constant recruitment, mortality and
movement. `init_age_strc` chooses how that equilibrium is solved, and
`ln_InitDevs` then deviates each age away from it.

## Usage

``` r
Get_Init_NAA(
  init_age_strc,
  init_iter,
  n_regions,
  n_pop,
  n_sexes,
  n_ages,
  n_seas,
  n_fish_fleets,
  seasdur,
  rec_seas_prop,
  natmort,
  natmort_annual = collapse_natmort_annual(natmort, seasdur, seas_dim = 3),
  init_F,
  dmr,
  fish_sel,
  ret_sel,
  R0_r,
  sexratio,
  Movement,
  do_recruits_move,
  ln_InitDevs,
  Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0
)
```

## Arguments

- init_age_strc:

  Integer, how the initial age structure is solved: `0` iterates the
  annual cycle `init_iter` times, `1` scalar geometric series with no
  movement at any age, `2` matrix geometric series with movement at
  every age, `3` movement below the plus group and a scalar series for
  the plus group, `4` no equilibrium at all, ages 2 and older are
  `exp(ln_InitDevs)` apportioned by sex ratio.

- init_iter:

  Integer, annual cycles run when `init_age_strc = 0`.

- n_regions, n_pop, n_sexes, n_ages, n_seas, n_fish_fleets:

  Integer dimensions. `n_ages` includes the plus group.

- seasdur:

  Numeric vector (`n_seas`) of each season's fraction of a year.

- rec_seas_prop:

  Matrix (`n_pop x n_seas`) of the share of annual recruitment entering
  in each season.

- natmort:

  Array (`n_pop x n_regions x n_seas x n_ages x n_sexes`) of natural
  mortality, a rate per year applied within each season.

- natmort_annual:

  Array (`n_pop x n_regions x n_ages x n_sexes`) of the annual total,
  the duration weighted sum over seasons, read by the steps that advance
  a whole year at once. Defaults to that sum.

- init_F:

  Numeric array (`n_regions x n_seas x n_fish_fleets`) of fully selected
  fishing mortality during initialization. Zero for an unfished
  population.

- dmr:

  Numeric array (`n_regions x n_seas x n_fish_fleets`) of the discard
  mortality rate during initialization.

- fish_sel, ret_sel:

  Arrays
  (`n_pop x n_regions x n_seas x n_ages x n_sexes x n_fish_fleets`) of
  total fishery selectivity at age and of the proportion of those fish
  retained.

- R0_r:

  Matrix (`n_pop x n_regions`) of unfished recruitment allocated to each
  region.

- sexratio:

  Array (`n_pop x n_regions x n_sexes`) of the proportion of recruits by
  sex.

- Movement:

  Array (`n_pop x n_regions x n_regions x n_seas x n_ages x n_sexes`) of
  seasonal movement probabilities, where `Movement[p,r,r2,,,]` is the
  fraction of the fish in region `r` that move to region `r2`.

- do_recruits_move:

  Integer, `0` recruits stay in their region for their first year, `1`
  recruits move with every other age.

- ln_InitDevs:

  Array (`n_pop x n_regions x (n_ages - 1) x n_sexes`) of log scale
  deviations for ages 2 and older. A 3-D array without the sex dimension
  is expanded across sexes as one shared curve. Under
  `init_age_strc = 4` these are the numbers themselves rather than
  multipliers on an equilibrium.

- Mrate:

  Array dimensioned like `Movement` of instantaneous movement rates.
  Required when `move_timing = 2`, ignored otherwise.

- move_timing:

  Integer ordering of movement and mortality within a season: `0`
  movement then mortality, `1` mortality then movement, `2` both at
  once. See
  [`build_seas_operator`](https://chengmatt.github.io/SPoRC/dev/reference/build_seas_operator.md).

- expm_nsub:

  Integer, how the matrix exponential is taken under `move_timing = 2`:
  `0` uses
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html),
  \\n \ge 1\\ the implicit backward Euler scheme. See
  [`mat_exp`](https://chengmatt.github.io/SPoRC/dev/reference/mat_exp.md).

## Value

Array (`n_pop x n_regions x n_ages x n_sexes`) of initial numbers at
age.
