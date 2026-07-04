# Get Composition Proportions from RTMB Output

Extracts and standardizes age and length composition data for fishery
and survey fleets from RTMB model output. The function processes both
observed and expected compositions for pooled (all populations combined)
and population-specific data streams, returning results in both
long-format data frames and array formats suitable for diagnostics and
likelihood evaluation.

## Usage

``` r
get_comp_prop(data, rep, age_labels, len_labels, year_labels)
```

## Arguments

- data:

  List. RTMB input data containing observed compositions, usage flags,
  composition types, ageing error matrices, and fleet/region/season/
  population dimensions.

- rep:

  List. RTMB model output containing predicted compositions (`CAA`,
  `CAL`, `DAA`, `DAL`, `SrvIAA`, `SrvIAL`) with a leading population
  dimension that is summed to produce pooled outputs.

- age_labels:

  Vector. Age bin labels used in composition data.

- len_labels:

  Vector. Length bin labels used in composition data.

- year_labels:

  Vector. Labels for model years corresponding to projection or
  estimation periods.

## Value

A named list containing:

- Fishery_Ages, Fishery_Lens:

  Long-format data frames of observed and predicted pooled fishery age
  and length compositions, including metadata (Region, Year, Season,
  Age/Len, Sex, Fleet, Type).

- Fishery_Discard_Ages, Fishery_Discard_Lens:

  Long-format data frames for discard compositions.

- Survey_Ages, Survey_Lens:

  Long-format data frames for survey-based compositions.

- Pop_Fishery_Ages, Pop_Fishery_Lens, Pop_Survey_Ages, Pop_Survey_Lens:

  Population-specific long-format composition outputs with an additional
  Pop column.

- Obs\_\*\_mat:

  Arrays of observed compositions indexed by Region × Year × Season ×
  Bin × Sex × Fleet (and Pop where applicable).

- Pred\_\*\_mat:

  Arrays of model-predicted compositions with matching dimensions to
  observed arrays.

## Details

Compositions include retained catch, discard catch, and survey
observations, and are aligned across regions, years, seasons, fleets,
sexes, and age/length bins. Observed compositions are optionally
filtered by usage flags defined in the input data object.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
comp_props <- get_comp_prop(
  data = data,
  rep = rep,
  age_labels = 2:31,
  len_labels = seq(41, 99, 2),
  year_labels = 1960:2024
)
} # }
```
