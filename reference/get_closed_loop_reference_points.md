# Get Closed Loop Reference Points

Computes fishery and biological reference points either using "true"
simulated values from the operating model or using assessment-derived
data and report objects. Supports single-region and multi-region
reference points.

## Usage

``` r
get_closed_loop_reference_points(
  use_true_values,
  sim_env,
  asmt_data = NULL,
  asmt_rep = NULL,
  y,
  sim,
  reference_points_opt = list(n_avg_yrs = 1, SPR_x = 0.4, calc_rec_st_yr = 1, rec_age =
    1, type = "single_region", what = "SPR"),
  n_proj_yrs
)
```

## Arguments

- use_true_values:

  Logical. If TRUE, uses values from the simulation environment
  (\`sim_env\`) for calculating reference points. If FALSE, uses
  \`asmt_data\` and \`asmt_rep\`.

- sim_env:

  Simulation environment

- asmt_data:

  Optional list. Assessment data object (from RTMB) if not using true
  values.

- asmt_rep:

  Optional list. Assessment report object (from RTMB) if not using true
  values.

- y:

  Integer. Number of years to include in calculations (usually the last
  year of the assessment or simulation).

- sim:

  Integer. Index of the simulation replicate in \`sim_env\`.

- reference_points_opt:

  List. Options for reference point calculations:

  n_avg_yrs

  :   Number of years to average over demographic rates. Default is 1.

  SPR_x

  :   Target SPR fraction for reference point calculations. Default is
      0.4.

  calc_rec_st_yr

  :   Year to start calculating mean recruitment. Default is 1.

  rec_age

  :   Age at recruitment. Default is 1.

  type

  :   Reference point type: "single_region" or "multi_region". Default
      is "single_region".

  what

  :   Method for reference point calculation. Options include "SPR",
      "BH_MSY", "independent_SPR", "independent_BH_MSY", "global_SPR",
      "global_BH_MSY". Default is "SPR".

- n_proj_yrs:

  Number of projection years

## Value

A list with elements:

- f_ref_pt:

  Array of fishing reference points by region and projection year.

- b_ref_pt:

  Array of biological reference points by region and projection year.

- virgin_b_ref_pt:

  Array of unfished biological reference points by region and projection
  year.

## See also

Other Closed Loop Simulations:
[`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_multifleet.md),
[`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_singlefleet.md),
[`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/reference/condition_closed_loop_simulations.md)
