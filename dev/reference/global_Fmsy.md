# Compute Beverton-Holt Fmsy for a spatially explicit model

Equilibrium yield under a Beverton-Holt stock-recruit relationship,
built from spawning biomass per recruit, the equilibrium recruitment
formula and catch at age integrated over regions, seasons and movement.
Yield counts landings from fleets with `is_discard_fleet == 0` only; a
discard-only fleet's F stays in the \\Z\\ denominator, so the two
mortality sources compete correctly. Covers multi-region
single-population models with seasonal movement; straying needs
`single_region_Fmsy`.

## Usage

``` r
global_Fmsy(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters, holding `log_Fmsy`, the log-scale trial
  \\F\_{MSY}\\.

- data:

  Named list of RTMB data: the dimensions `n_regions`, `n_ages` and
  `n_seas`; `seasdur` `[n_seas]`; `spawn_seas` and `t_spawn`;
  `F_fract_flt` `[n_regions, n_seas, n_fish_fleets]`, the fleet F
  fractions by region; `fish_sel` and `ret_sel`
  `[1, n_regions, n_seas, n_ages, n_fish_fleets]`, the female
  selectivity and the retained fraction of it; `dmr`
  `[n_regions, n_seas, n_fish_fleets]`, the fraction of discards that
  die; `natmort` `[n_regions, n_ages]`; `WAA` and `MatAA`
  `[n_regions, n_seas, n_ages]`; `Movement`
  `[n_regions, n_regions, n_seas, n_ages]`; `rec_region_prop` and
  `sex_ratio_f` `[n_regions]`; `rec_seas_prop` `[n_seas]`; the steepness
  `h` and unfished recruitment `R0`; and `is_discard_fleet`
  `[n_fish_fleets]`, 1 for fleets whose catch is left out of landed
  yield while still contributing to \\Z\\.

## Value

Numeric scalar, the negative equilibrium yield, minimized to find
\\F\_{MSY}\\.

## Details

Fishing mortality at age splits into retained,
`F_ret = F * selectivity * retention`, and dead discards,
`F_disc = F * selectivity * (1 - retention) * dmr`, so survival runs on
`Z = M + F_ret + F_disc`. Retained fish always die, only the fraction
`dmr` of discards die, and the surviving `(1 - dmr)` keeps ageing,
moving and spawning. Landed yield comes from the Baranov equation on the
landed fraction alone.
