# Compute global SPR reference point for a spatially explicit model

Spatially integrated SPR from a per-recruit cohort tracked across every
region and season under movement, with one scalar \\F_x\\ applied across
regions and scaled by the region's fleet fractions and selectivity.
Returns the squared penalty \\100 (SPR - SPR_x)^2\\ for the optimizer.

## Usage

``` r
global_SPR(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters, holding `log_F_x`, the log-scale trial
  fishing mortality.

- data:

  Named list of RTMB data: the dimensions `n_pop`, `n_regions`, `n_ages`
  and `n_seas`; `seasdur` `[n_seas]`; `spawn_seas` and `t_spawn`;
  `F_fract_flt` `[n_regions, n_seas, n_fish_fleets]`, the fleet F
  fractions by region; `fish_sel` and `ret_sel`
  `[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]`, the female
  selectivity and the retained fraction of it; `dmr`
  `[n_regions, n_seas, n_fish_fleets]`, the fraction of discards that
  die; `natmort` `[n_pop, n_regions, n_ages]`; `WAA` and `MatAA`
  `[n_pop, n_regions, n_seas, n_ages]`; `Movement`
  `[n_pop, n_regions, n_regions, n_seas, n_ages]` and
  `sgl_seas_spawning_movement` `[n_pop, n_regions, n_regions, n_ages]`;
  `do_recruits_move`; `rec_region_prop` and `sex_ratio_f`
  `[n_pop, n_regions]`; `rec_seas_prop` `[n_pop, n_seas]`; `stray_rate`
  and `natal_region` `[n_pop]`; `n_pop_in_region` `[n_regions]`; and the
  target `SPR_x`.

## Value

Numeric scalar, the squared penalty \\(SPR - SPR_x)^2\\.

## Details

Single and multi-population models are both covered. When `n_pop > 1`,
effective SSB at each population's natal region accumulates straying
contributions from the others, and when `n_seas = 1` as well,
`sgl_seas_spawning_movement` redistributes fish to the natal grounds
before SSB is computed. The plus group is solved analytically through
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md).

Fishing mortality at age splits into retained,
`F_ret = F * selectivity * retention`, and dead discards,
`F_disc = F * selectivity * (1 - retention) * dmr`, so survival runs on
`Z = M + F_ret + F_disc`. Retained fish always die, only the fraction
`dmr` of discards die, and the surviving `(1 - dmr)` keeps ageing,
moving and spawning.
