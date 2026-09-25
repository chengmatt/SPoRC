# Compute SPR reference point for a single-region or non-spatial model

The spawning potential ratio at a trial fishing mortality \\F_x\\,
returned as the squared penalty \\100 (SPR - SPR_x)^2\\ the outer
optimizer minimizes to find \\F\_{SPR_x}\\. Several populations are
covered through stray rates, but there is no spatial movement.

## Usage

``` r
single_region_SPR(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters, holding `log_F_x`, the log-scale trial
  fishing mortality.

- data:

  Named list of RTMB data: the dimensions `n_pop`, `n_ages` and
  `n_seas`; `seasdur` `[n_seas]`; `spawn_seas` and `t_spawn`;
  `F_fract_flt` `[n_seas, n_fish_fleets]`; `fish_sel` and `ret_sel`
  `[n_pop, n_seas, n_ages, n_fish_fleets]`, the female selectivity and
  the retained fraction of it; `dmr` `[n_seas, n_fish_fleets]`, the
  fraction of discards that die; `natmort` `[n_pop, n_ages]`; `WAA` and
  `MatAA` `[n_pop, n_seas, n_ages]`; `sex_ratio_f` `[n_pop]`;
  `rec_seas_prop` `[n_pop, n_seas]`; `stray_rate` and `natal_region`
  `[n_pop]`; `n_pop_in_region` `[n_regions]`; and the target `SPR_x`.

## Value

Numeric scalar, the squared penalty \\(SPR - SPR_x)^2\\, zero at \\F =
F\_{SPR_x}\\.

## Details

Fishing mortality at age splits into retained,
`F_ret = F * selectivity * retention`, and dead discards,
`F_disc = F * selectivity * (1 - retention) * dmr`, so survival runs on
`Z = M + F_ret + F_disc`. Retained fish always die, only the fraction
`dmr` of discards die, and the surviving `(1 - dmr)` keeps ageing and
spawning.

## See also

[`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md)
