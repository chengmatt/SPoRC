# Generate Key Projection Quantities and Table Plot

Computes reference points and runs a short deterministic projection to
give terminal spawning biomass, catch advice and status ratios by model
and region, returned as a tidy data frame and a formatted table plot.

## Usage

``` r
get_key_quants(
  data,
  rep,
  reference_points_opt,
  proj_model_opt,
  model_names,
  pars = NULL
)
```

## Arguments

- data:

  List of length `n_models`, each a SPoRC data list holding the region,
  year, age, fleet and biological inputs.

- rep:

  List of length `n_models`, each a SPoRC report list from
  `obj$report()` after optimization, holding recruitment, selectivity,
  fishing mortality and numbers at age.

- reference_points_opt:

  Named list passed to
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
  holding `SPR_x` (the target spawning potential ratio, which may be
  `NULL` under `type = "bh_msy"`), `t_spawn` (the fraction of the year
  elapsed before spawning), `sex_ratio_f` `(n_pop, n_regions)`,
  `calc_rec_st_yr` (the first model year averaged over for recruitment),
  `rec_age`, `type` (e.g. `"multi_region"`) and `what` (e.g.
  `"global_SPR"`), and optionally `calc_rec_yrs`, the model years
  recruitment is averaged and projected from, which overrides
  `calc_rec_st_yr` and `rec_age`.

- proj_model_opt:

  Named list passed to
  [`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
  holding `n_proj_yrs`, `n_avg_yrs` (terminal years the demographic
  inputs are averaged over before being kept constant across the
  projection), `HCR_function` with signature
  `function(x, frp, brp, ...)`, `recruitment_opt` (`"mean_rec"`,
  `"bh_rec"`, `"zero_rec"` or `"inv_gauss"`), and `fmort_opt` (`"HCR"`
  to apply the control rule to each region's SSB, `"HCR_global"` to
  apply it to SSB summed across regions as global reference points need,
  or `"Input"` to fish at the F reference point every year without the
  control rule). It can also hold the choices
  [`Setup_Mod_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md)
  offers, so the advice here matches a projection inside the model:

  `rec_opt`

  :   Recruitment deviations in the projection years. `"mean"` (the
      default here) applies none, `"propagate"` applies the fit's
      projected deviations, and `"terminal"` repeats the terminal year's
      deviation. The last two need `pars`.

  `sel_opt`

  :   Fishery and retention selectivity. `"average"` (default) averages
      the last `n_avg_yrs` years, and `"propagate"` uses the years the
      fit projected selectivity, set by `n_proj_yrs_devs` in
      [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md).

  `move_opt`

  :   Movement. `"average"` (default) averages the last `n_avg_yrs`
      years, and `"propagate"` uses the years the fit projected
      movement.

  `growth_opt`

  :   Weight at age. `"average"` (default) averages the last `n_avg_yrs`
      years, and `"propagate"` uses the years the fit projected growth,
      which needs weight at age derived by a growth model.

  `naa_opt`

  :   Numbers at age under a state-space fit. `"terminal"` (default)
      projects from the terminal numbers with no innovations, and
      `"propagate"` applies the fit's projected innovations, which needs
      `pars`.

  `F_constraint`

  :   A catch cap or stability limit from
      [`catch_cap`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md)
      or
      [`catch_stability`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md),
      or a list of them applied in order. Default none.

- model_names:

  Character vector of length `n_models` of display names.

- pars:

  List of length `n_models`, each the fitted parameter list from
  `fit$env$parList()`. Needed only when `rec_opt` is `"propagate"` or
  `"terminal"`, or `naa_opt` is `"propagate"`, since those read the
  fit's deviations. Default `NULL`.

## Value

A list of two. The first is a data frame of key quantities by model and
region with columns `Model`, `Region`, `Terminal_SSB`, `Terminal_SSB0`,
`Terminal_F`, `Catch_Advice`, `B_Ref_Pt`, `F_Ref_Pt`, `B_over_B_Ref`,
`B_over_DynB_Ref` and `F_over_F_Ref`. The second is a `cowplot` `ggdraw`
object rendering the same quantities as a table.

## Details

Each model has its reference points computed, its demographic inputs
averaged over the last `n_avg_yrs` years, the population projected
forward `n_proj_yrs` years, and terminal SSB, dynamic unfished SSB,
catch advice (year 2 of the projection) and the status ratios extracted.

Under `recruitment_opt = "bh_rec"` the stock-recruit parameters reach
the projection through an internal `srr_opt` list built from year-1
demographics to approximate unfished SSB. Under `"inv_gauss"` a warning
is issued, since only one deterministic simulation is run and stochastic
recruitment belongs in a full MSE loop.

## Note

These quantities are approximate and are not catch advice. This is a
simplified projection interface for rapid model comparison and
diagnostic screening; stochastic recruitment, closed-loop feedback and
fleet-specific control rules need
[`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md)
directly.

## See also

[`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md)

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`Setup_Mod_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md),
[`catch_cap()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md),
[`catch_stability()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md),
[`hcr_threshold()`](https://chengmatt.github.io/SPoRC/dev/reference/hcr_threshold.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  reference_points_opt <- list(
    SPR_x          = 0.4,
    t_spawn        = 0,
    sex_ratio_f    = array(0.5, dim = c(n_pop, n_regions)),
    calc_rec_st_yr = 20,
    rec_age        = 2,
    type           = "multi_region",
    what           = "global_SPR"
  )

  proj_model_opt <- list(
    n_proj_yrs      = 2,
    n_avg_yrs       = 1,
    HCR_function    = function(x, frp, brp, alpha = 0.05) {
      stock_status <- x / brp
      if (stock_status >= 1)                            frp
      else if (stock_status > alpha) frp * (stock_status - alpha) / (1 - alpha)
      else                                              0
    },
    recruitment_opt = "mean_rec",
    fmort_opt       = "HCR_global"
  )

  out <- get_key_quants(
    data                  = list(mlt_rg_sable_data),
    rep                   = list(mlt_rg_sable_rep),
    reference_points_opt  = reference_points_opt,
    proj_model_opt        = proj_model_opt,
    model_names           = "Model 1"
  )

  out[[1]]  # key quantities data frame
  out[[2]]  # table plot
} # }
```
