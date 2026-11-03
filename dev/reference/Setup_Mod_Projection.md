# Set Up A Projection Inside The Estimation Model

Projects the population forward inside the model, so that `sdreport`
gives standard errors for projected spawning biomass, catch, fishing
mortality, numbers at age and stock status, with the uncertainty of the
fit and of future recruitment in them. The projection is the same one
[`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md)
runs after a fit, but fed the model's own estimates rather than copies
of them: the terminal year's numbers at age and F, the demographics
averaged over the terminal years, and mean recruitment. Projection year
1 is the terminal assessment year run again, so year 2 is the first year
after the assessment and the year catch advice is read from.

## Usage

``` r
Setup_Mod_Projection(
  input_list,
  n_proj_yrs = 0,
  f_opt = "terminal",
  f_input = NULL,
  catch_input = NULL,
  f_ref_pt = NULL,
  b_ref_pt = NULL,
  hcr_function = hcr_threshold,
  hcr_scope = "region",
  F_constraint = NULL,
  refpt_type = "single_region",
  refpt_what = "SPR",
  SPR_x = NULL,
  sex_ratio_f = NULL,
  is_discard_fleet = NULL,
  local_bh_msy_newton_steps = 6,
  refpt_newton_steps = 2,
  refpt_log_f_start = log(0.1),
  n_avg_yrs = 1,
  sel_opt = "average",
  move_opt = "average",
  growth_opt = "average",
  rec_opt = "propagate",
  naa_opt = "propagate",
  calc_rec_st_yr = 1,
  rec_age = 1,
  calc_rec_yrs = NULL
)
```

## Arguments

- input_list:

  List with `$data`, `$par` and `$map` from the earlier setup calls.

- n_proj_yrs:

  Integer. Number of years to project, counting the terminal assessment
  year, which is run again as projection year 1. 0 (default) leaves the
  projection off and the model as it was.

- f_opt:

  How fishing mortality is set in the projection years. `"terminal"`
  (default) fishes every year at the terminal year's F. `"input"` fishes
  at the F given in `f_input`. `"hcr"` runs a harvest control rule on
  reference points given in `f_ref_pt` and `b_ref_pt`, treated as known
  numbers. `"refpt"` runs the rule on reference points the model solves
  for itself from its own estimates, so the uncertainty in F40, B40 or
  MSY is part of the projection's standard errors.

- f_input:

  Array `[n_regions, n_proj_yrs]` of the F to fish at in each region and
  projection year, read under `f_opt = "input"`. Year 1 replays the
  terminal year, so column 1 is not read.

- catch_input:

  Retained catch to take in each region and projection year, an array
  `[n_regions, n_proj_yrs]`, or `[n_regions, n_proj_yrs, n_seas]` for
  catch by season, in the units of `proj_Catch`. The model finds the F
  that takes it, so the catch has no uncertainty and the uncertainty
  moves into that F and the projected biomass. A year left `NA` in every
  region is fished under `f_opt`, so known catch in the first years can
  be followed by the rule; a year is given in every region or in none.
  Year 1 replays the terminal year, so column 1 is not read.
  `F_constraint` applies only in the years left to `f_opt`. Default
  `NULL`, no catch given.

- f_ref_pt:

  The F reference point the rule reads under `f_opt = "hcr"`: the target
  F, one value per region, or `[n_regions, n_proj_yrs]` to vary it by
  year. A fixed number here, so its own uncertainty is not in the
  standard errors; `f_opt = "refpt"` is the way to include it.

- b_ref_pt:

  The biomass reference point the rule reads under `f_opt = "hcr"`: one
  value per population and region, or `[n_pop, n_regions, n_proj_yrs]`
  to vary it by year. A fixed number here, as `f_ref_pt` is.

- hcr_function:

  The harvest control rule: a function of `x`, the spawning biomass the
  rule reads, `frp`, the F reference point, and `brp`, the biomass
  reference point, returning the F to apply, as
  [`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md)
  takes it. Default
  [`hcr_threshold`](https://chengmatt.github.io/SPoRC/dev/reference/hcr_threshold.md):
  the target F above `brp`, no fishing below 0.05 times `brp`, and a
  straight line between. For another cutoff, pass
  `function(x, frp, brp) hcr_threshold(x, frp, brp, alpha = 0.1)`. A
  rule of your own runs inside the model, where the spawning biomass it
  is given is a quantity the model is differentiating, and such a
  quantity cannot be compared with `>` or `<`. Compare `ad_value(x)`,
  the plain number, to decide which part of the rule applies, and use
  `x` itself in the arithmetic that follows, as `hcr_threshold` does.

- hcr_scope:

  Which spawning biomass a region's rule reads. `"region"` (default),
  its own, against its own biomass reference point. `"global"`, the
  total over regions against the reference points summed over regions,
  so every region shares one stock status.

- F_constraint:

  A limit on the catch the rule's or input F would take.
  [`catch_cap`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md)
  lowers F wherever the catch would exceed a cap;
  [`catch_stability`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md)
  keeps the catch within a band of the previous year's, lowering or
  raising F; a list of them is applied in the order given, so put the
  cap last when nothing may exceed it. `NULL` (default) applies none.

- refpt_type:

  Under `f_opt = "refpt"`, `"single_region"` (default) or
  `"multi_region"`, named as in
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md).

- refpt_what:

  Under `f_opt = "refpt"`, which reference point the model solves, named
  as in
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md):
  `"SPR"` (default) or `"MSY"` for one region, and `"independent_SPR"`,
  `"independent_MSY"`, `"global_SPR"`, `"global_MSY"` or `"local_MSY"`
  for several. The MSY options need a model fitted with a stock-recruit
  curve.

- SPR_x:

  The target spawning potential ratio for the SPR reference points, for
  example 0.4 for F40 and B40.

- sex_ratio_f:

  Array `[n_pop, n_regions]`. Female share of recruits the per-recruit
  calculation behind the reference points starts from. Default `NULL`
  reads it from the model: in a two-sex model the recruitment sex ratio
  averaged over the last `n_avg_yrs` years, so an estimated sex ratio
  and its uncertainty reach the reference points; in a one-sex model
  0.5, since the model halves its spawning biomass. A value given here
  is used instead, for instance to match a post hoc
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md)
  run.

- is_discard_fleet:

  Integer vector `[n_fish_fleets]`, 1 for a fleet whose catch is
  discarded rather than landed, so the MSY options do not count it as
  yield. Default all zero.

- local_bh_msy_newton_steps:

  Integer. Steps the local MSY solve takes for the equilibrium
  recruitment at each trial F. Default 6, as in
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md).

- refpt_newton_steps:

  Integer. The reference points are solved numerically first, then
  refined with this many Newton steps that the model differentiates
  through, which is what gives a reference point its standard error.
  Default 2, which is enough. 0 leaves the reference points as fixed
  numbers with no standard error.

- refpt_log_f_start:

  Numeric. Where the numerical solve for a reference point starts, on
  the log scale. Default `log(0.1)`.

- n_avg_yrs:

  Integer. How many terminal years weight at age, maturity, natural
  mortality, sex ratio and straying are averaged over, for the
  projection and for the reference points alike. Default 1, the terminal
  year alone.

- sel_opt:

  Selectivity in the projection years. `"average"` (default) repeats its
  average over the last `n_avg_yrs`. `"propagate"` uses the years the
  model projected itself, which needs `n_proj_yrs_devs` in
  [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
  and a time-varying selectivity, so its projected deviations and their
  uncertainty come along.

- move_opt:

  Movement in the projection years, with the same two choices as
  `sel_opt`: `"average"` (default) over the last `n_avg_yrs`, or
  `"propagate"` from the years the model projected itself.

- growth_opt:

  Weight at age in the projection years, when a growth model derives it
  (`waa_model = "wt_len"` in
  [`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md)).
  `"average"` (default) repeats its average over the last `n_avg_yrs`.
  `"propagate"` uses the years the model projected growth itself, which
  needs `n_proj_yrs_devs` in
  [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
  time-varying growth, and `growth_tv_type = "curve"`, so the projected
  growth deviations and their uncertainty come along. Weight at age read
  from data has nothing to propagate and is always averaged.

- rec_opt:

  Recruitment in the projection years, always mean recruitment over the
  years `calc_rec_yrs` picks, times a deviation. `"propagate"` (default)
  uses the deviations the model estimated for the projection years,
  which needs `n_proj_yrs_devs` in
  [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md).
  With no data behind them they sit at zero, so the projection runs at
  the mean, but each has the recruitment process error as its standard
  error, so the projected standard errors widen as a forecast's should.
  `"mean"` applies no deviation: recruitment is the mean, and its
  uncertainty is only that of the mean, which is small. `"terminal"`
  applies the terminal year's deviation in every projection year, one
  value shared by all of them, so the projected cohorts rise or fall
  together rather than averaging out.

- naa_opt:

  The numbers at age state of a state-space model in the projection
  years. `"propagate"` (default) applies the state innovations the model
  estimated for the projection years, which needs `n_proj_yrs_devs` in
  [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md);
  they sit at their forecast, zero under `"iid"`, and bring the state's
  process error into the standard errors. `"terminal"` advances the
  terminal year's state with no further innovations. Ignored by a model
  without a state.

- calc_rec_st_yr:

  Integer. First model year of the window mean recruitment is averaged
  over, as in
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md);
  the window is `calc_rec_st_yr:(n_years - rec_age)`. Default 1, the
  first model year.

- rec_age:

  Integer. Years left off the end of that window, the most recent years
  whose recruitment is still poorly known. Default 1.

- calc_rec_yrs:

  Integer vector of model years to average recruitment over instead of
  the window, for a mean that runs to the terminal year or skips years.
  Model year indices, 1 for the first year, not calendar years. Default
  `NULL`, the window.

## Value

`input_list` with the settings stored in `$data` as `proj_*` fields.

## Details

[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
fits the model with the projection switched off, since the projection
adds nothing to the likelihood, then rebuilds it at the estimates with
the projection on.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`catch_cap()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md),
[`catch_stability()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md),
[`hcr_threshold()`](https://chengmatt.github.io/SPoRC/dev/reference/hcr_threshold.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# six years under the threshold rule on post hoc F40 and B40, capped at 25 kt
input_list <- Setup_Mod_Projection(input_list, n_proj_yrs = 6, f_opt = "hcr",
                                   f_ref_pt = refpts$f_ref_pt, b_ref_pt = refpts$b_ref_pt,
                                   F_constraint = catch_cap(25))
# 20 kt taken in each of the next two years, then the same rule
catch <- array(NA, dim = c(1, 6))
catch[1, 2:3] <- 20
input_list <- Setup_Mod_Projection(input_list, n_proj_yrs = 6, f_opt = "hcr",
                                   f_ref_pt = refpts$f_ref_pt, b_ref_pt = refpts$b_ref_pt,
                                   catch_input = catch)
fit <- fit_model(input_list$data, input_list$par, input_list$map, do_sdrep = TRUE)
fit$sd_rep # log_proj_SSB, log_proj_Catch, log_proj_F, log_proj_Total_Biom
} # }
```
