# Cohort growth, one year at a time

The year-loop companion of
[`Get_Growth`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth.md)
under cohort propagation. For year `y` it builds every key, length,
weight and weight of the year from the start-of-year state, and carries
the state to the next year with the year's parameters: each propagated
age grows by the year's increment, the first propagated age is placed on
the year's curve, the ages in the linear phase keep the length at `A1`
their cohort was born with, and the plus group's size next year blends
the cohort entering it with the fish already there by their numbers at
the start of this year, \$\$\bar L\_{+,y+1} = \frac{(N\_{n-1} + 0.01)\\
g(L\_{n-1}) + (N\_{n} + 0.01)\\ g(L\_{+})}{N\_{n-1} + N_n + 0.02}\$\$
with \\g\\ the year's increment and \\N\\ the start-of-year numbers of
the stratum. The CV at age is held at the first year's sizes.

## Usage

``` r
Get_Growth_Year(
  growth,
  y,
  NAA_y,
  ln_growth_pars,
  ln_growth_devs,
  growth_tv_model,
  growth_tv_link,
  growth_par_bounds,
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
  n_seas,
  n_sexes,
  t_fish,
  t_srv,
  ln_growth_semipar_devs = NULL,
  growth_semipar = 0
)
```

## Arguments

- growth:

  The list
  [`Get_Growth`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth.md)
  returned, carrying the output containers and the years before the
  propagation started.

- y:

  Year index to evaluate.

- NAA_y:

  Numbers at age at the start of year `y`, array
  `[pop, region, age, sex]`, read for the plus-group blend.

- ln_growth_pars, ln_growth_devs, growth_tv_model, growth_tv_link,
  growth_par_bounds:

  As in
  [`Get_Growth`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth.md).

- growth_A1, growth_A2, growth_L0, growth_len_lower, growth_cv_type,
  growth_sd_type, growth_dist, growth_plus_group:

  As in
  [`Get_Growth`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth.md).

- derive_waa, wt_len_pars, ages, seasdur, spawn_seas, t_spawn, n_pop,
  n_regions, n_seas, n_sexes:

  As in
  [`Get_Growth`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth.md).

- t_fish, t_srv:

  Timings as in
  [`Get_Growth`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Growth.md).

## Value

The `growth` list with year `y` filled and, when a year follows,
`L_beg[,, y + 1,,]` set to the state the next year starts from.
