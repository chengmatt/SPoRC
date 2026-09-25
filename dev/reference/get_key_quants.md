# Generate Key Projection Quantities and Table Plot

Computes reference points and runs a short deterministic projection to
give terminal spawning biomass, catch advice and status ratios by model
and region, returned as a tidy data frame and a formatted table plot.

## Usage

``` r
get_key_quants(data, rep, reference_points_opt, proj_model_opt, model_names)
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
  `"global_SPR"`).

- proj_model_opt:

  Named list passed to
  [`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
  holding `n_proj_yrs`, `n_avg_yrs` (terminal years the demographic
  inputs are averaged over before being held constant across the
  projection), `HCR_function` with signature
  `function(x, frp, brp, ...)`, `recruitment_opt` (`"mean_rec"`,
  `"bh_rec"`, `"zero_rec"` or `"inv_gauss"`), and `fmort_opt` (`"input"`
  to hold terminal F or `"HCR"`).

- model_names:

  Character vector of length `n_models` of display names.

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
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md)

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
    fmort_opt       = "HCR"
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
