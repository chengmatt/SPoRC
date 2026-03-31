# Generate Key Projection Quantities and Table Plot

Calculates biological and fishery reference points and performs
population projections to estimate terminal spawning biomass, catch
advice, and reference point values by model and region. Also returns a
formatted table plot of key quantities.

## Usage

``` r
get_key_quants(data, rep, reference_points_opt, proj_model_opt, model_names)
```

## Arguments

- data:

  A list of model input data objects, one for each model (i.e., a list
  of SPoRC-formatted data lists). Each element should contain
  information on regions, years, ages, fleets, and biological inputs
  (e.g., weight-at-age, maturity, mortality).

- rep:

  A list of model output objects, one for each model (i.e., a list of
  SPoRC-formatted report lists). Each element must include recruitment,
  selectivity, mortality, and numbers-at-age.

- reference_points_opt:

  A named list specifying options for reference point calculations. See
  [`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/reference/Get_Reference_Points.md)
  for more details. Must include:

  SPR_x

  :   Spawning potential ratio (e.g., 0.4) for calculating F reference
      points. May be `NULL` if using `bh_rec`.

  t_spwn

  :   Fraction of year when spawning occurs (e.g., 0.5).

  sex_ratio_f

  :   Proportion of recruits that are female.

  calc_rec_st_yr

  :   Start year for averaging recruitment.

  rec_age

  :   Recruitment age.

  type

  :   Reference point calculation method (e.g., "multi_region").

  what

  :   Type of output requested from the reference point function.

- proj_model_opt:

  A named list of projection settings. See
  [`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/reference/Do_Population_Projection.md)
  for details. Must include:

  n_proj_yrs

  :   Number of years to project forward.

  HCR_function

  :   Harvest control rule function to use.

  recruitment_opt

  :   Recruitment assumption (e.g., "mean_rec", "bh_rec", "inv_gauss").

  fmort_opt

  :   Fishing mortality assumption (e.g., "input", "HCR").

  n_avg_yrs

  :   Number of years to average over for projection inputs.

- model_names:

  A character vector of model identifiers (e.g., c("Base", "Alt1",
  "Alt2")), one for each element in `data` and `rep`.

## Value

A list with two elements:

- `[[1]]`:

  A data.frame of key quantities by model and region, including terminal
  SSB, catch advice, and reference points.

- `[[2]]`:

  A cowplot table plot (ggdraw object) of the same key quantities.

## Details

This function checks input list completeness, calculates reference
points using
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/reference/Get_Reference_Points.md),
performs population projections with
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/reference/Do_Population_Projection.md),
and assembles both tabular and visual summaries.

If `recruitment_opt` is set to "inv_gauss", a warning is issued since
only a single simulation will be run. This is typically not appropriate
and an alternative assumption is recommended.

## See also

[`Get_Reference_Points`](https://chengmatt.github.io/SPoRC/reference/Get_Reference_Points.md),
[`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/reference/Do_Population_Projection.md)

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/reference/Do_Population_Projection.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/reference/Get_Reference_Points.md)

## Examples

``` r
if (FALSE) { # \dontrun{
reference_points_opt <- list(SPR_x = 0.4,
                             t_spwn = 0,
                             sex_ratio_f = 0.5,
                             calc_rec_st_yr = 20,
                             rec_age = 2,
                             type = "multi_region",
                             what = "global_SPR")

proj_model_opt <- list(
  n_proj_yrs = 2,
  n_avg_yrs = 1,
  HCR_function = function(x, frp, brp, alpha = 0.05) {
    stock_status <- x / brp
    if (stock_status >= 1) f <- frp
    if (stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
    if (stock_status < alpha) f <- 0
    return(f)
  },
  recruitment_opt = "mean_rec",
  fmort_opt = "HCR"
)

out <- get_key_quants(list(mlt_rg_sable_data),
                      list(mlt_rg_sable_rep),
                      reference_points_opt,
                      proj_model_opt,
                      "Model 1")
out[[1]]  # key quantities data.frame
out[[2]]  # table plot
} # }
```
