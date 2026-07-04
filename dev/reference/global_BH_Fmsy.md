# Compute Beverton-Holt Fmsy for a spatially explicit model

Calculates \\F\_{MSY}\\ by maximising equilibrium yield under a
Beverton-Holt stock-recruit relationship. Yield is computed from
spawning biomass per recruit (\\\phi_F\\), the BH equilibrium
recruitment formula, and catch-at-age integrated across all regions,
seasons, and movement transitions. Yield includes only landings from
fleets where `is_discard_fleet == 0`; discard-only fleets contribute to
total mortality but not to the yield being maximised.

## Usage

``` r
global_BH_Fmsy(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters. Must contain:

  `log_Fmsy`

  :   Log-scale trial \\F\_{MSY}\\.

- data:

  Named list of RTMB data. Must contain:

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

  :   Numeric array `[1, n_regions, n_seas, n_ages, n_fish_fleets]`.
      Fishery selectivity at age for females.

  `ret_sel`

  :   Numeric array `[1, n_regions, n_seas, n_ages, n_fish_fleets]`.
      Retention selectivity (fraction of selected fish retained).

  `dmr`

  :   Numeric array `[n_regions, n_seas, n_fish_fleets]`. Discard
      mortality rate (fraction of discarded fish that die).

  `natmort`

  :   Numeric array `[n_regions, n_ages]`. Female natural mortality at
      age.

  `WAA`

  :   Numeric array `[n_regions, n_seas, n_ages]`. Female weight at age.

  `MatAA`

  :   Numeric array `[n_regions, n_seas, n_ages]`. Maturity at age.

  `Movement`

  :   Numeric array `[n_regions, n_regions, n_seas, n_ages]`. Seasonal
      movement transition matrices.

  `rec_region_prop`

  :   Numeric vector `[n_regions]`. Proportion of recruitment entering
      each region.

  `sex_ratio_f`

  :   Numeric vector `[n_regions]`. Female sex ratio at recruitment.

  `rec_seas_prop`

  :   Numeric vector `[n_seas]`. Seasonal recruitment proportions.

  `h`

  :   Numeric. Beverton-Holt steepness.

  `R0`

  :   Numeric. Unfished equilibrium recruitment.

  `is_discard_fleet`

  :   Integer vector `[n_fish_fleets]`. Indicator for fleets whose catch
      is excluded from landed yield (0 = landing fleet, 1 = discard-only
      fleet). These fleets still contribute to total fishing mortality
      `Z` and affect population dynamics and spawning biomass.

## Value

Numeric scalar. Negative total equilibrium yield (minimised to find
\\F\_{MSY}\\).

## Details

Supports multi-region, single-population models with seasonal movement.
Straying is not included here (use `single_region_BH_Fmsy` for
multi-population non-spatial models).

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

Landed yield used in the objective function is computed via the Baranov
catch equation using only the landed fraction of fishing mortality
(excluding fleets where `is_discard_fleet == 1`). The discard fleet's F
remains in the Z denominator, so the partitioning correctly accounts for
competition between landing and discard mortality sources.
