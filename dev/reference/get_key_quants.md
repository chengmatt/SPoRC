# Generate Key Projection Quantities and Table Plot

Calculates biological and fishery reference points and performs
short-term population projections to estimate terminal spawning biomass,
catch advice, and reference point ratios by model and region. Returns
both a tidy data frame and a formatted table plot of the assembled
quantities.

## Usage

``` r
get_key_quants(data, rep, reference_points_opt, proj_model_opt, model_names)
```

## Arguments

- data:

  A list of length `n_models`, where each element is a SPoRC-formatted
  data list containing region, year, age, fleet, and biological inputs
  (e.g., weight-at-age, maturity, natural mortality).

- rep:

  A list of length `n_models`, where each element is a SPoRC-formatted
  report list (i.e., the output of `obj$report()` after optimization).
  Each element must include recruitment, selectivity, fishing mortality,
  and numbers-at-age arrays.

- reference_points_opt:

  A named list of options passed to
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md).
  Required elements:

  `SPR_x`

  :   Target spawning potential ratio (e.g., `0.4`) for SPR-based F
      reference points. May be `NULL` when `type = "bh_msy"`.

  `t_spawn`

  :   Fraction of the year elapsed before spawning occurs (e.g., `0` for
      start-of-year, `0.5` for mid-year).

  `sex_ratio_f`

  :   Array of dimension `(n_pop, n_regions)` giving the proportion of
      recruits that are female.

  `calc_rec_st_yr`

  :   Index of the first model year to include when averaging
      recruitment for reference point calculations.

  `rec_age`

  :   Age at recruitment (used to lag the recruitment series relative to
      terminal year).

  `type`

  :   Reference point calculation method; e.g., `"multi_region"` or
      `"single_region"`.

  `what`

  :   Output selector passed to `Get_Reference_Points`; e.g.,
      `"global_SPR"` or `"local_MSY"`.

- proj_model_opt:

  A named list of projection settings passed to
  [`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md).
  Required elements:

  `n_proj_yrs`

  :   Number of years to project forward.

  `n_avg_yrs`

  :   Number of terminal model years over which demographic inputs
      (selectivity, weight-at-age, maturity, natural mortality,
      movement) are averaged before being held constant across the
      projection period.

  `HCR_function`

  :   A harvest control rule function with signature
      `function(x, frp, brp, ...)`, where `x` is current biomass, `frp`
      is the F reference point, and `brp` is the biomass reference
      point.

  `recruitment_opt`

  :   Recruitment assumption for projection years. One of `"mean_rec"`,
      `"bh_rec"`, `"zero_rec"`, or `"inv_gauss"`. See Details.

  `fmort_opt`

  :   How fishing mortality is set during the projection. One of
      `"input"` (hold terminal F constant) or `"HCR"` (apply
      `HCR_function`).

- model_names:

  Character vector of length `n_models` giving display names for each
  model run (e.g., `c("Base", "Alt1")`).

## Value

A list of length 2:

- `[[1]]`:

  A data frame of key quantities by model and region, with columns
  `Model`, `Region`, `Terminal_SSB`, `Terminal_SSB0`, `Terminal_F`,
  `Catch_Advice`, `B_Ref_Pt`, `F_Ref_Pt`, `B_over_B_Ref`,
  `B_over_DynB_Ref`, and `F_over_F_Ref`.

- `[[2]]`:

  A `cowplot` `ggdraw` object rendering the same quantities as a
  formatted table, suitable for inclusion in a PDF report.

## Details

For each model, the function: (1) computes reference points via
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md);
(2) averages demographic inputs over the last `n_avg_yrs` model years;
(3) projects the population forward `n_proj_yrs` years via
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md);
and (4) extracts terminal SSB, dynamic unfished SSB, catch advice (year
2 of the projection), and status ratios.

When `recruitment_opt = "bh_rec"`, Beverton-Holt stock-recruit
parameters are passed to the projection via an internal `srr_opt` list
constructed from year-1 demographics to approximate unfished SSB. When
`recruitment_opt = "inv_gauss"`, a warning is issued because only a
single deterministic simulation is run; stochastic recruitment options
should be used within a full MSE loop rather than here.

## Note

The quantities returned by this function are **approximate and should
not be treated as official catch advice**. This wrapper provides only a
simplified projection interface; full projection capability, including
stochastic recruitment, closed-loop feedback, and fleet-specific harvest
control rules, requires calling
[`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md)
directly. Results here are intended for rapid model comparison and
diagnostic screening only.

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
