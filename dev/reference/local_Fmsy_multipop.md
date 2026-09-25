# Compute local Beverton-Holt Fmsy for a spatially explicit multi-population model

The multi-population form of
[`local_Fmsy_sglpop`](https://chengmatt.github.io/SPoRC/dev/reference/local_Fmsy_sglpop.md):
region-specific \\F\_{MSY}\\ values that jointly maximize total
equilibrium yield when populations with distinct natal regions, movement
and Beverton-Holt parameters share a spatial domain. The objective is
the negative of that yield.

## Usage

``` r
local_Fmsy_multipop(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters, holding `log_Fmsy`, the log-scale trial
  \\F\_{MSY}\\ values, one per region.

- data:

  Named list of RTMB data, holding every field
  [`global_SPR`](https://chengmatt.github.io/SPoRC/dev/reference/global_SPR.md)
  needs plus `h` `[n_pop, n_regions]`, the steepness at each
  population's natal region; `R0` `[n_pop]`; `stray_rate` `[n_pop]`, the
  fraction contributing to non-natal spawning regions; `natal_region`
  `[n_pop]`; `n_pop_in_region` `[n_regions]`, which normalizes the
  straying; `newton_steps`; and `is_discard_fleet` `[n_fish_fleets]`, 1
  for fleets whose catch is left out of landed yield while still
  contributing to \\Z\\.

## Value

Numeric scalar, the negative total equilibrium yield across regions,
minimized to obtain the regional \\F\_{MSY}\\ vector.

## Details

Cohorts are tracked per recruit over
`[population x origin x destination x age x season]`, with recruitment
spread by `rec_seas_prop`, `rec_region_prop` and `sex_ratio_f`. Recruits
enter in the first season at age one, with later seasonal contributions
added within the first age class before movement and mortality. Movement
is applied each season through region and age-specific transition
matrices, fishing mortality splits into retained and discarded
components, and mortality acts continuously within a season. Catch at
age accumulates over fleets, seasons and regions on the landed fraction
alone, while `Z` includes every fleet.

Spawning biomass per recruit is taken at the spawning season, after
movement and partial mortality up to `t_spawn`; with one season and
several populations, `sgl_seas_spawning_movement` redistributes fish to
the natal regions first. The plus group is solved analytically through
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md).

Effective spawning biomass at each population's natal region takes
contributions from every population through straying, scaled by
`stray_rate` and normalized by `n_pop_in_region` to preserve mass
balance. Equilibrium recruitment per population comes from a
Newton-Raphson solve of the coupled Beverton-Holt system, whose Jacobian
accounts for the cross-population dependence straying induces. Total
yield is then catch at age integrated over populations, regions and
seasons, scaled by equilibrium recruitment and the origin-region
proportions.
