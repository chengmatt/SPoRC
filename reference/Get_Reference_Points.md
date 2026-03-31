# Wrapper function to get reference points

Wrapper function to compute fishing and biological reference points
given data and report objects from an assessment or simulation. Supports
both single-region and multi-region calculations with options for SPR or
Beverton–Holt MSY reference points.

## Usage

``` r
Get_Reference_Points(
  data,
  rep,
  SPR_x = NULL,
  t_spwn = 0,
  sex_ratio_f = rep(0.5, data$n_regions),
  calc_rec_st_yr = 1,
  rec_age = 1,
  type,
  what,
  n_avg_yrs = 1,
  local_bh_msy_newton_steps = 6
)
```

## Arguments

- data:

  List. Data object containing ages, years, weight-at-age, maturity,
  natural mortality, and other simulation/assessment info.

- rep:

  List. Report object from RTMB containing estimated parameters like
  Fmort, selectivity, recruitment, steepness.

- SPR_x:

  Numeric. Target Spawning Potential Ratio fraction. Required for
  SPR-based reference points.

- t_spwn:

  Numeric. Mortality time until spawning.

- sex_ratio_f:

  Numeric vector. Female sex ratio by region.

- calc_rec_st_yr:

  Integer. First year used to compute mean recruitment.

- rec_age:

  Integer. Age at recruitment.

- type:

  Character. "single_region" or "multi_region".

- what:

  Character. Type of reference point:

  SPR

  :   Single-region SPR reference point

  independent_SPR

  :   Multi-region SPR without movement

  global_SPR

  :   Multi-region SPR with movement

  BH_MSY

  :   Single-region Beverton–Holt MSY

  independent_BH_MSY

  :   Multi-region BH-MSY without movement

  global_BH_MSY

  :   Multi-region global BH-MSY with movement

  local_BH_MSY

  :   Multi-region local BH-MSY with movement

- n_avg_yrs:

  Integer. Number of years to average demographic rates when calculating
  reference points.

- local_bh_msy_newton_steps:

  Number of newton steps to take to solve for equilibrium recruitment in
  the origin region when local_BH_MSY is assumed.

## Value

A list with elements:

- f_ref_pt:

  Vector of fishing reference points for each region.

- b_ref_pt:

  Vector of biological reference points for each region.

- virgin_b_ref_pt:

  Vector of virgin biomass reference points for each region.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/reference/Do_Population_Projection.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/reference/get_key_quants.md)
