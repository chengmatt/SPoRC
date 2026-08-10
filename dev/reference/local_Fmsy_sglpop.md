# Compute region-specific Beverton-Holt Fmsy for a spatially explicit single-population model

Computes the vector of regional \\F\_{MSY}\\ values that jointly
maximise total equilibrium yield across all regions under a spatially
explicit Beverton-Holt stock-recruit relationship. Unlike
[`global_Fmsy`](https://chengmatt.github.io/SPoRC/dev/reference/global_Fmsy.md),
which constrains all regions to share a single fishing mortality, this
function allows each region to have its own optimal \\F\\.

## Usage

``` r
local_Fmsy_sglpop(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters. Must contain:

  `log_Fmsy`

  :   Numeric vector `[n_regions]`. Log-scale trial \\F\_{MSY}\\ values,
      one per region.

- data:

  Named list of RTMB data. Must contain all spatial fields required by
  [`global_SPR`](https://chengmatt.github.io/SPoRC/dev/reference/global_SPR.md)
  (excluding `SPR_x`, `stray_rate`, and `natal_region`) plus:

  `h`

  :   Numeric vector `[n_regions]`. Beverton-Holt steepness by region.

  `R0`

  :   Numeric scalar. Total unfished equilibrium recruitment.

  `rec_region_prop`

  :   Numeric vector `[n_regions]`. Proportion of annual recruitment
      entering each region.

  `newton_steps`

  :   Integer. Number of Newton-Raphson iterations used to solve for
      equilibrium recruitment by origin region.

  `is_discard_fleet`

  :   Integer vector `[n_fish_fleets]`. Indicator for fleets whose catch
      is excluded from landed yield (0 = landing fleet, 1 = discard-only
      fleet). These fleets still contribute to total fishing mortality
      `Z` and affect population dynamics and spawning biomass.

## Value

Numeric scalar. Negative total equilibrium yield across all regions.
This is minimised to obtain the vector of regional \\F\_{MSY}\\ values.

## Details

Cohorts originating in each region are tracked separately through
seasonal movement, mortality, and ageing using an
`[origin, destination]` per-recruit accounting framework. Spawning
biomass per recruit is accumulated by origin and destination region, and
the plus group is solved analytically using
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md).

Equilibrium recruitment by origin region \\R\_{eq,o}\\ is solved using a
Newton-Raphson algorithm applied to the fixed-point condition that
recruitment produced at each destination region (via the BH relationship
applied to effective SSB) equals the recruitment attributed to that
origin. The Jacobian is derived analytically using the quotient rule and
the chain rule through the spatial redistribution of spawning biomass.

Yield is computed using only the landed fraction of fishing mortality,
excluding fleets flagged as discard-only via `is_discard_fleet`.
Discard-only fleets still contribute to total mortality `Z` and affect
population dynamics and spawning biomass.

Fishing mortality is decomposed into retained and discarded components:

- Retained fishing mortality: \$\$F^{\mathrm{ret}}\_{r,a,s,f} =
  F\_{MSY,r} \\ F\_{\mathrm{fract},r,s,f} \\ \mathrm{sel}\_{r,a,s,f} \\
  \mathrm{ret}\_{r,a,s,f}\$\$

- Discard fishing mortality (dead discards only):
  \$\$F^{\mathrm{disc}}\_{r,a,s,f} = F\_{MSY,r} \\
  F\_{\mathrm{fract},r,s,f} \\ \mathrm{sel}\_{r,a,s,f} \\ (1 -
  \mathrm{ret}\_{r,a,s,f}) \\ \mathrm{dmr}\_{r,s,f}\$\$

- Total instantaneous mortality: \$\$Z\_{r,a,s} = M\_{r,a} \\
  \mathrm{seasdur}\_s + F^{\mathrm{ret}}\_{r,a,s} +
  F^{\mathrm{disc}}\_{r,a,s}\$\$

Landed yield used in the objective function excludes catch from fleets
where `is_discard_fleet == 1`. The Baranov catch equation partitions
landed F out of total Z, so the discard fleet's contribution to
mortality is properly accounted for in the denominator.

Seasonal movement is applied using the
`Movement[origin, dest, seas, age]` array. Recruitment may move
immediately or only after age-1 depending on `do_recruits_move`.
Spawning biomass is accumulated at `spawn_seas` with fractional
mortality `t_spawn`.

The plus group is solved analytically using the transition matrices
produced by
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and the solver
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md).
