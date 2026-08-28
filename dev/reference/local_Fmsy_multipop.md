# Compute local Beverton-Holt Fmsy for a spatially explicit multi-population model

Multi-population extension of
[`local_Fmsy_sglpop`](https://chengmatt.github.io/SPoRC/dev/reference/local_Fmsy_sglpop.md).
Estimates a vector of region-specific \\F\_{MSY}\\ values that jointly
maximize total equilibrium yield when multiple populations, each with
distinct natal regions, movement schedules, and Beverton-Holt
parameters, co-occupy a shared spatial domain.

## Usage

``` r
local_Fmsy_multipop(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters. Must contain:

  `log_Fmsy`

  :   Log-scale trial \\F\_{MSY}\\ values, one per region (length
      `n_regions`).

- data:

  Named list of RTMB data. Must contain all fields required by
  [`global_SPR`](https://chengmatt.github.io/SPoRC/dev/reference/global_SPR.md)
  plus:

  `h`

  :   Numeric array `[n_pop, n_regions]`. Beverton-Holt steepness
      evaluated at each population's natal region.

  `R0`

  :   Numeric vector `[n_pop]`. Unfished equilibrium recruitment per
      population.

  `stray_rate`

  :   Numeric vector `[n_pop]`. Fraction of individuals contributing to
      non-natal spawning regions.

  `natal_region`

  :   Integer vector `[n_pop]`. Natal region index for each population.

  `n_pop_in_region`

  :   Integer vector `[n_regions]`. Number of populations sharing each
      natal region (used to normalize straying).

  `newton_steps`

  :   Integer. Number of Newton-Raphson iterations used to solve for
      equilibrium recruitment.

  `is_discard_fleet`

  :   Integer vector `[n_fish_fleets]`. Indicator for fleets whose catch
      is excluded from landed yield (0 = landing fleet, 1 = discard-only
      fleet). These fleets still contribute to total fishing mortality
      `Z` and affect population dynamics and spawning biomass.

## Value

Numeric scalar. Negative total equilibrium yield across all regions.
This objective is minimized to obtain the vector of regional
\\F\_{MSY}\\ values.

## Details

Cohorts are tracked using a per-recruit framework indexed by
`[population x origin x destination x age x season]`. Recruitment is
distributed across seasons (`rec_seas_prop`), regions
(`rec_region_prop`), and sex (`sex_ratio_f`). Initial recruits are
assigned in the first season at age-1, with additional seasonal
recruitment contributions added within the first age class prior to
movement and mortality.

Movement is applied at each seasonal step using region- and age-specific
transition matrices. Fishing mortality is decomposed into retained and
discarded components, and total mortality is applied continuously within
each season. Catch-at-age is accumulated across fleets, seasons, and
regions using only the landed fraction of fishing mortality (excluding
fleets flagged as discard-only via `is_discard_fleet`), while total
mortality `Z` includes all fleets.

Spawning biomass per recruit (SBPR) is computed at the spawning season
after applying movement and partial mortality up to the spawning time
(`t_spawn`). For single-season models with multiple populations,
`sgl_seas_spawning_movement` redistributes individuals to natal spawning
regions prior to SSB calculation.

The plus group is solved analytically using
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md),
ensuring a consistent equilibrium solution for the terminal age class.

Effective spawning biomass at each population's natal region includes
contributions from all populations via straying. Stray contributions are
scaled by `stray_rate` and normalized by `n_pop_in_region` to preserve
mass balance.

Equilibrium recruitment by population is obtained via a Newton-Raphson
algorithm that solves the coupled Beverton-Holt system. The Jacobian
accounts for cross-population dependence of spawning biomass induced by
straying.

Total equilibrium yield is computed by integrating catch-at-age over all
populations, regions, and seasons, scaled by equilibrium recruitment and
origin-region proportions.

Recruitment is implemented as a per-recruit process and later scaled by
equilibrium recruitment. Seasonal recruitment proportions are applied
within the first age class, allowing intra-annual timing of recruitment
before movement and mortality are applied.

The objective function is the negative of total yield, summed across all
populations and regions. Yield includes only landings from fleets where
`is_discard_fleet == 0`; discard-only fleets contribute to mortality but
not to the yield being maximized.
