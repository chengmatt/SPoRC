# Get Francis Weights

Computes Francis composition reweighting factors for pooled and
population-specific fishery and survey age and length compositions. Used
inside
[`run_francis`](https://chengmatt.github.io/SPoRC/dev/reference/run_francis.md)
or directly in a user-defined reweighting loop.

## Usage

``` r
do_francis_reweighting(data, rep, age_labels, len_labels, year_labels)
```

## Arguments

- data:

  List of model data inputs containing observed compositions, usage
  flags, input sample sizes, weight arrays, and composition type
  matrices for both pooled and population-specific data sources.

- rep:

  Report list from a fitted RTMB model containing predicted compositions
  (`CAA`, `CAL`, `SrvIAA`, `SrvIAL`) with a leading population
  dimension.

- age_labels:

  Vector of observed age bin labels.

- len_labels:

  Vector of observed length bin labels.

- year_labels:

  Vector of assessment year labels.

## Value

A named list containing:

- `new_fish_age_wts`, `new_fish_len_wts`, `new_srv_age_wts`,
  `new_srv_len_wts`:

  Updated Francis weight arrays for pooled compositions, same dimensions
  as the corresponding `Wt_*` arrays in `data`.

- `new_fish_age_pop_wts`, `new_fish_len_pop_wts`, `new_srv_age_pop_wts`,
  `new_srv_len_pop_wts`:

  Updated Francis weight arrays for population-specific compositions,
  same dimensions as the corresponding `Wt_*_pop` arrays in `data`.
  Cells remain `NA` when the corresponding `Use*_pop` flag contains no
  ones.

- `mean_francis`:

  Long-format dataframe of observed and expected composition means
  across all data sources and populations, with columns `Region`,
  `Comp_Year`, `Comp_Seas`, `Sex`, `Fleet`, `obs`, `pred`, `Type`, and
  `Pop` (pooled rows have `Pop = NA`).

## See also

Other Francis Reweighting:
[`run_francis()`](https://chengmatt.github.io/SPoRC/dev/reference/run_francis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
for(j in 1:5) {
  if(j == 1) {
    data$Wt_FishAgeComps[] <- 1; data$Wt_FishLenComps[] <- 1
    data$Wt_SrvAgeComps[]  <- 1; data$Wt_SrvLenComps[]  <- 1
    data$Wt_FishAgeComps_pop[] <- 0; data$Wt_FishLenComps_pop[] <- 0
    data$Wt_SrvAgeComps_pop[]  <- 0; data$Wt_SrvLenComps_pop[]  <- 0
  } else {
    data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
    data$Wt_FishLenComps[] <- wts$new_fish_len_wts
    data$Wt_SrvAgeComps[]  <- wts$new_srv_age_wts
    data$Wt_SrvLenComps[]  <- wts$new_srv_len_wts
    if(any(data$UseFishAgeComps_pop == 1)) data$Wt_FishAgeComps_pop[] <- wts$new_fish_age_pop_wts
    if(any(data$UseFishLenComps_pop == 1)) data$Wt_FishLenComps_pop[] <- wts$new_fish_len_pop_wts
    if(any(data$UseSrvAgeComps_pop  == 1)) data$Wt_SrvAgeComps_pop[]  <- wts$new_srv_age_pop_wts
    if(any(data$UseSrvLenComps_pop  == 1)) data$Wt_SrvLenComps_pop[]  <- wts$new_srv_len_pop_wts
  }
  obj <- fit_model(data, parameters, mapping, random = NULL,
                   newton_loops = 3, silent = TRUE)
  rep <- obj$report(obj$env$last.par.best)
  wts <- do_francis_reweighting(data = data, rep = rep, age_labels = 2:31,
                                len_labels = seq(41, 99, 2), year_labels = 1960:2024)
}
} # }
```
