# Compute global SPR reference point for a spatially explicit model

Calculates a single, spatially-integrated SPR using a per-recruit cohort
tracked across all regions and seasons under movement. A single scalar
\\F_x\\ is applied uniformly across regions (scaled by region-specific
fleet fractions and selectivity). Returns the squared penalty \\100
(SPR - SPR_x)^2\\ for optimisation.

## Usage

``` r
global_SPR(pars, data)
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

  `n_regions`

  :   Integer. Number of spatial regions.

  `n_ages`

  :   Integer. Number of age classes.

  `n_seas`

  :   Integer. Number of seasons.

  `seasdur`

  :   Numeric vector `[n_seas]`. Season durations.

  `spawn_seas`

  :   Integer. Index of the spawning season.

  `t_spawn`

  :   Numeric. Mid-season spawning timing correction.

  `F_fract_flt`

  :   Numeric array `[n_regions, n_seas, n_fish_fleets]`. Fleet F
      fractions by region.

  `fish_sel`

  :   Numeric array `[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]`.
      Female fishery selectivity.

  `ret_sel`

  :   Numeric array `[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]`.
      Retention selectivity (fraction of selected fish that are
      retained).

  `dmr`

  :   Numeric array `[n_regions, n_seas, n_fish_fleets]`. Discard
      mortality rate (fraction of discarded fish that die).

  `natmort`

  :   Numeric array `[n_pop, n_regions, n_ages]`. Female natural
      mortality.

  `WAA`

  :   Numeric array `[n_pop, n_regions, n_seas, n_ages]`. Female weight
      at age.

  `MatAA`

  :   Numeric array `[n_pop, n_regions, n_seas, n_ages]`. Maturity at
      age.

  `Movement`

  :   Numeric array `[n_pop, n_regions, n_regions, n_seas, n_ages]`.
      Seasonal movement transition matrices.

  `sgl_seas_spawning_movement`

  :   Numeric array `[n_pop, n_regions, n_regions, n_ages]`. Spawning
      movement for single-season natal homing models.

  `do_recruits_move`

  :   Integer (0/1). Whether age-1 recruits are subject to movement.

  `rec_region_prop`

  :   Numeric array `[n_pop, n_regions]`. Proportion of recruitment
      entering each region.

  `sex_ratio_f`

  :   Numeric array `[n_pop, n_regions]`. Female sex ratio at
      recruitment.

  `rec_seas_prop`

  :   Numeric array `[n_pop, n_seas]`. Seasonal recruitment proportions.

  `stray_rate`

  :   Numeric vector `[n_pop]`. Per-population stray rate.

  `natal_region`

  :   Integer vector `[n_pop]`. Natal region index for each population.

  `n_pop_in_region`

  :   Integer vector `[n_regions]`. Number of populations per natal
      region.

  `SPR_x`

  :   Numeric. Target SPR fraction.

## Value

Numeric scalar. Squared penalty \\(SPR - SPR_x)^2\\.

## Details

Supports single- and multi-population models. When `n_pop > 1`,
effective SSB at each population's natal region accumulates straying
contributions from other populations. When `n_seas = 1` and `n_pop > 1`,
`sgl_seas_spawning_movement` redistributes fish to natal grounds before
SSB is computed.

The plus-group is solved analytically using
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md).

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
population and continues aging, moving, and contributing to spawning
biomass.
