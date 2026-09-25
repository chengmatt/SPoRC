# Compute fishing and biological reference points from an assessment or simulation

Builds the data list, calls the inner objective through RTMB, and
returns the fishing and biological reference points a projection or
control rule needs. Covers single-region and spatially explicit models,
on SPR or Beverton-Holt MSY.

## Usage

``` r
Get_Reference_Points(
  data,
  rep,
  SPR_x = NULL,
  t_spawn = 0,
  sex_ratio_f = array(0.5, dim = c(data$n_pop, data$n_regions)),
  calc_rec_st_yr = 1,
  rec_age = 1,
  type,
  what,
  n_avg_yrs = 1,
  local_bh_msy_newton_steps = 6,
  is_discard_fleet = array(0, dim = data$n_fish_fleets)
)
```

## Arguments

- data:

  SPoRC data object holding the age structure, weight-at-age, maturity,
  natural mortality, seasons and spatial configuration.

- rep:

  SPoRC report object holding `Fmort`, `fish_sel`, `natmort`, `Rec`,
  `SSB`, `h_trans`, `R0`, `rec_region_prop`, `rec_seas_prop`, `Movement`
  and `stray_rate`.

- SPR_x:

  Target spawning potential ratio, required under `"SPR"`,
  `"independent_SPR"` and `"global_SPR"`.

- t_spawn:

  Fraction of the spawning season elapsed before spawning, for the
  mid-season mortality correction. Default 0.

- sex_ratio_f:

  Numeric array `[n_pop, n_regions]` of the female sex ratio at
  recruitment. Default 0.5.

- calc_rec_st_yr:

  First year included in the mean historical recruitment the biological
  reference points are scaled by. Default 1.

- rec_age:

  Recruitment lag in years, excluding the most recent years from that
  mean. Default 1.

- type:

  Spatial structure. `"single_region"` has no movement and takes `"SPR"`
  or `"MSY"`; `"multi_region"` takes `"independent_SPR"`,
  `"independent_MSY"`, `"global_SPR"`, `"global_MSY"` or `"local_MSY"`.

- what:

  Reference point method. `"SPR"` and `"MSY"` are the single-region
  \\F\_{SPR_x}\\ and Beverton-Holt \\F\_{MSY}\\. The `"independent_"`
  pair computes each region on its own without movement. `"global_SPR"`
  gives one shared \\F\_{SPR_x}\\ with movement, integrated across
  regions, and `"global_MSY"` the same for \\F\_{MSY}\\,
  single-population models only. `"local_MSY"` gives region-specific
  \\F\_{MSY}\\ values that jointly maximize total yield with movement,
  for single and multi-population models alike.

- n_avg_yrs:

  Terminal years the demographic rates (selectivity, natural mortality,
  weight, maturity, movement) are averaged over first. Default 1.

- local_bh_msy_newton_steps:

  Newton-Raphson iterations used to solve equilibrium recruitment by
  origin region under `what = "local_MSY"`. Raise it if convergence is
  suspect. Default 6.

- is_discard_fleet:

  Integer vector `[n_fish_fleets]`, 1 for fleets whose catch is left out
  of landed yield while still contributing to \\Z\\ and so to the
  population dynamics and spawning biomass. Read by the MSY methods
  only. Default all zeros.

## Value

A named list holding `f_ref_pt` `[n_regions]`, one value per region,
shared across regions under the global methods; `b_ref_pt`
`[n_pop, n_regions]`, the equilibrium spawning biomass at the reference
point (\\SBPR_F \times R\_{eq}\\ or \\SBPR_F \times \bar{R}\\);
`virgin_b_ref_pt`, its unfished counterpart (\\SBPR_0 \times R_0\\ or
\\SBPR_0 \times \bar{R}\\); and `pop_b_ref_pt` and
`virgin_pop_b_ref_pt`, the population-specific effective spawning
biomass at each population's natal region, with the stray contributions
from the other populations.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md)
