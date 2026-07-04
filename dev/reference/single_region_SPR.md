# Compute SPR reference point for a single-region or non-spatial model

Calculates the spawning potential ratio (SPR) as a function of a trial
fishing mortality \\F_x\\, then returns a squared penalty \\100 (SPR -
SPR_x)^2\\ that is minimised by the outer optimizer to find
\\F\_{SPR_x}\\. Supports multiple populations via stray rates but does
not include spatial movement.

## Usage

``` r
single_region_SPR(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters. Must contain:

  `log_F_x`

  :   Log-scale trial fishing mortality.

- data:

  Named list of RTMB data. Must contain:

  `n_pop`

  :   Integer. Number of populations.

  `n_ages`

  :   Integer. Number of age classes.

  `n_seas`

  :   Integer. Number of seasons.

  `seasdur`

  :   Numeric vector `[n_seas]`. Season durations.

  `spawn_seas`

  :   Integer. Index of the spawning season.

  `t_spawn`

  :   Numeric. Fraction of spawning season elapsed before spawning
      (mid-season mortality correction).

  `F_fract_flt`

  :   Numeric array `[n_seas, n_fish_fleets]`. Fleet F fractions.

  `fish_sel`

  :   Numeric array `[n_pop, n_seas, n_ages, n_fish_fleets]`. Fishery
      selectivity at age for females.

  `ret_sel`

  :   Numeric array `[n_pop, n_seas, n_ages, n_fish_fleets]`. Retention
      selectivity (fraction of selected fish retained).

  `dmr`

  :   Numeric array `[n_seas, n_fish_fleets]`. Discard mortality rate
      (fraction of discarded fish that die).

  `natmort`

  :   Numeric array `[n_pop, n_ages]`. Female natural mortality at age.

  `WAA`

  :   Numeric array `[n_pop, n_seas, n_ages]`. Female weight at age.

  `MatAA`

  :   Numeric array `[n_pop, n_seas, n_ages]`. Maturity at age.

  `sex_ratio_f`

  :   Numeric vector `[n_pop]`. Female sex ratio at recruitment.

  `rec_seas_prop`

  :   Numeric array `[n_pop, n_seas]`. Proportion of annual recruitment
      entering in each season.

  `stray_rate`

  :   Numeric vector `[n_pop]`. Per-population stray rate used to
      compute effective SSB across populations.

  `natal_region`

  :   Integer vector `[n_pop]`. Natal region index for each population.

  `n_pop_in_region`

  :   Integer vector `[n_regions]`. Number of populations per natal
      region.

  `SPR_x`

  :   Numeric. Target SPR fraction (e.g. 0.4).

## Value

Numeric scalar. Squared penalty \\(SPR - SPR_x)^2\\. Minimised to zero
at \\F = F\_{SPR_x}\\.

## Details

\*\*Fishing mortality decomposition\*\*

Fishing mortality at age is split into:

\- retained fishing mortality `F_ret = F * selectivity * retention`

\- discard fishing mortality (dead discards only)
`F_disc = F * selectivity * (1 - retention) * dmr`

where `dmr` is the discard mortality rate (fraction of discarded fish
that die). Only the dead fraction contributes to total instantaneous
mortality `Z`.

The total mortality used for survival is:

`Z = M + F_ret + F_disc`

This formulation assumes:

\- retained fish always die, - only a fraction `dmr` of discarded fish
die, - the surviving fraction `(1 - dmr)` of discards remains in the
population and continues aging and contributing to spawning biomass.

## See also

[`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md)
