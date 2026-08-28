# Compute fishing and biological reference points from an assessment or simulation

Wrapper that constructs the appropriate data list, calls the relevant
inner objective function via RTMB, and returns fishing and biological
reference points for use in projections or harvest control rules.
Supports single-region and spatially explicit multi-region models, with
options for SPR-based or Beverton-Holt MSY-based reference points.

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

  List. SPoRC data object containing age structure, weight-at-age,
  maturity, natural mortality, seasons, and spatial configuration.

- rep:

  List. SPoRC report object from RTMB containing estimated or simulated
  quantities including `Fmort`, `fish_sel`, `natmort`, `Rec`, `SSB`,
  `h_trans`, `R0`, `rec_region_prop`, `rec_seas_prop`, `Movement`, and
  `stray_rate`.

- SPR_x:

  Numeric. Target spawning potential ratio fraction (e.g. 0.4). Required
  when `what` is `"SPR"`, `"independent_SPR"`, or `"global_SPR"`.

- t_spawn:

  Numeric. Fraction of the spawning season elapsed before spawning, used
  for the mid-season mortality correction. Default = 0.

- sex_ratio_f:

  Numeric array `[n_pop, n_regions]`. Female sex ratio at recruitment.
  Default = 0.5 everywhere.

- calc_rec_st_yr:

  Integer. First year included when computing mean historical
  recruitment for biological reference point scaling. Default = 1.

- rec_age:

  Integer. Recruitment lag in years, used to exclude the most recent
  years from the mean recruitment calculation. Default = 1.

- type:

  Character. Spatial structure of the model:

  `"single_region"`

  :   No spatial movement; supports `"SPR"` and `"MSY"`.

  `"multi_region"`

  :   Spatially explicit; supports `"independent_SPR"`,
      `"independent_MSY"`, `"global_SPR"`, `"global_MSY"`, and
      `"local_MSY"`.

- what:

  Character. Reference point method:

  `"SPR"`

  :   Single-region \\F\_{SPR_x}\\.

  `"MSY"`

  :   Single-region Beverton-Holt \\F\_{MSY}\\.

  `"independent_SPR"`

  :   Per-region \\F\_{SPR_x}\\ computed independently for each region
      without movement.

  `"independent_MSY"`

  :   Per-region \\F\_{MSY}\\ computed independently for each region
      without movement.

  `"global_SPR"`

  :   Single shared \\F\_{SPR_x}\\ with movement, integrated across all
      regions.

  `"global_MSY"`

  :   Single shared \\F\_{MSY}\\ with movement. Valid for
      single-population models only.

  `"local_MSY"`

  :   Region-specific \\F\_{MSY}\\ values that jointly maximize total
      yield with movement. Valid for both single- and multi-population
      models.

- n_avg_yrs:

  Integer. Number of terminal years over which demographic rates
  (selectivity, natural mortality, weight, maturity, movement) are
  averaged before computing reference points. Default = 1.

- local_bh_msy_newton_steps:

  Integer. Number of Newton-Raphson iterations used to solve for
  equilibrium recruitment by origin region when `what = "local_MSY"`.
  Increase if convergence is suspect. Default = 6.

- is_discard_fleet:

  Integer vector `[n_fish_fleets]`. Indicator for fleets whose catch
  should be excluded from landed yield when computing MSY-based
  reference points (0 = landing fleet, 1 = discard-only fleet). These
  fleets still contribute to total fishing mortality `Z` and affect
  population dynamics and spawning biomass. Only used by Beverton-Holt
  MSY methods (`"MSY"`, `"independent_MSY"`, `"global_MSY"`,
  `"local_MSY"`); ignored for SPR-based methods. Default is all zeros
  (all fleets are landing fleets).

## Value

A named list:

- `f_ref_pt`:

  Numeric vector `[n_regions]`. Fishing mortality reference point by
  region. All regions share the same value for global methods; regions
  have independent values for local or independent methods.

- `b_ref_pt`:

  Numeric array `[n_pop, n_regions]`. Equilibrium spawning biomass at
  the reference point by population and region (\\SBPR_F \times
  R\_{eq}\\ or \\SBPR_F \times \bar{R}\\).

- `virgin_b_ref_pt`:

  Numeric array `[n_pop, n_regions]`. Virgin (unfished) spawning biomass
  by population and region (\\SBPR_0 \times R_0\\ or \\SBPR_0 \times
  \bar{R}\\).

- `pop_b_ref_pt`:

  Numeric array `[n_pop, n_regions]`. Population-specific effective
  spawning biomass at the reference point, evaluated at each
  population's natal region and incorporating stray contributions from
  other populations.

- `virgin_pop_b_ref_pt`:

  Numeric array `[n_pop, n_regions]`. Population-specific effective
  virgin spawning biomass, evaluated at each population's natal region.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md)
