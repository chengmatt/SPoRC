# Compute region-specific Beverton-Holt Fmsy for a spatially explicit single-population model

The vector of regional \\F\_{MSY}\\ values that jointly maximize total
equilibrium yield, where
[`global_Fmsy`](https://chengmatt.github.io/SPoRC/dev/reference/global_Fmsy.md)
constrains every region to one fishing mortality. The objective is the
negative of that yield.

## Usage

``` r
local_Fmsy_sglpop(pars, data)
```

## Arguments

- pars:

  Named list of RTMB parameters, holding `log_Fmsy` `[n_regions]`, the
  log-scale trial values.

- data:

  Named list of RTMB data, holding every spatial field
  [`global_SPR`](https://chengmatt.github.io/SPoRC/dev/reference/global_SPR.md)
  needs apart from `SPR_x`, `stray_rate` and `natal_region`, plus `h`
  `[n_regions]`, the scalar `R0`, `rec_region_prop` `[n_regions]`,
  `newton_steps`, and `is_discard_fleet` `[n_fish_fleets]`, 1 for fleets
  whose catch is left out of landed yield while still contributing to
  \\Z\\.

## Value

Numeric scalar, the negative total equilibrium yield across regions,
minimized to obtain the regional \\F\_{MSY}\\ vector.

## Details

Cohorts from each region are tracked separately through seasonal
movement, mortality and ageing on an `[origin, destination]` per-recruit
accounting, with spawning biomass per recruit accumulated by origin and
destination and the plus group solved analytically through
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md).
Movement uses `Movement[origin, dest, seas, age]`, recruits move
immediately or from age one depending on `do_recruits_move`, and
spawning biomass accumulates at `spawn_seas` under the fractional
mortality `t_spawn`.

Equilibrium recruitment by origin region is solved by Newton-Raphson on
the fixed point where the recruitment each destination produces, through
the curve applied to effective SSB, equals the recruitment attributed to
that origin. The Jacobian is derived analytically with the quotient and
chain rules through the spatial redistribution of spawning biomass.

Fishing mortality splits into \$\$F^{\mathrm{ret}}\_{r,a,s,f} =
F\_{MSY,r} \\ F\_{\mathrm{fract},r,s,f} \\ \mathrm{sel}\_{r,a,s,f} \\
\mathrm{ret}\_{r,a,s,f}\$\$ and the dead discards
\$\$F^{\mathrm{disc}}\_{r,a,s,f} = F\_{MSY,r} \\
F\_{\mathrm{fract},r,s,f} \\ \mathrm{sel}\_{r,a,s,f} \\ (1 -
\mathrm{ret}\_{r,a,s,f}) \\ \mathrm{dmr}\_{r,s,f}\$\$ giving
\$\$Z\_{r,a,s} = M\_{r,a} \\ \mathrm{seasdur}\_s +
F^{\mathrm{ret}}\_{r,a,s} + F^{\mathrm{disc}}\_{r,a,s}\$\$ Landed yield
leaves out the catch of fleets with `is_discard_fleet == 1`, while the
Baranov equation keeps their F in the \\Z\\ denominator, so the two
mortality sources compete correctly.
