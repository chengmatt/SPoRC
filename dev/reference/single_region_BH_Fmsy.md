# Compute Beverton-Holt Fmsy for a single-region or non-spatial model

Finds \\F\_{MSY}\\ by maximising equilibrium yield under a Beverton-Holt
stock-recruit relationship. Yield is computed from spawning biomass per
recruit (\\\phi_F\\), the BH equilibrium recruitment formula, and
catch-at-age integrated across all seasons. Supports multiple
populations via stray rates but does not include spatial movement. Yield
includes only landings from fleets where `is_discard_fleet == 0`;
discard-only fleets contribute to total mortality but not to the yield
being maximised.

## Usage

``` r
single_region_BH_Fmsy(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters. Must contain:

  `log_Fmsy`

  :   Log-scale trial \\F\_{MSY}\\.

- data:

  Named list of RTMB data. Must contain all fields required by
  [`single_region_SPR`](https://chengmatt.github.io/SPoRC/dev/reference/single_region_SPR.md)
  plus:

  `h`

  :   Numeric vector `[n_pop]`. Beverton-Holt steepness.

  `R0`

  :   Numeric vector `[n_pop]`. Unfished equilibrium recruitment.

  `is_discard_fleet`

  :   Integer vector `[n_fish_fleets]`. Indicator for fleets whose catch
      is excluded from landed yield (0 = landing fleet, 1 = discard-only
      fleet). These fleets still contribute to total fishing mortality
      `Z` and affect population dynamics and spawning biomass.

## Value

Numeric scalar. Negative total equilibrium yield (minimised to find
\\F\_{MSY}\\).
