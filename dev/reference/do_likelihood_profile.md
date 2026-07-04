# Run Likelihood Profile

Profiles the joint negative log-likelihood and all individual likelihood
components across a range of fixed values for a single parameter.
Supports both sequential and parallel execution.

## Usage

``` r
do_likelihood_profile(
  data,
  parameters,
  mapping,
  random = NULL,
  what,
  idx = NULL,
  min_val,
  max_val,
  inc = 0.05,
  do_par = FALSE,
  n_cores = NULL
)
```

## Arguments

- data:

  Data list from the fitted model.

- parameters:

  Parameter list from the fitted model.

- mapping:

  Mapping list from the fitted model.

- random:

  Character vector of random effects to estimate. Default `NULL`.

- what:

  Character string. Name of the parameter to profile.

- idx:

  Vector pointing to the specific elements to fix when
  `parameters[[what]]` is an array. `NULL` for scalar parameters.

- min_val:

  Numeric. Minimum value of the profile range.

- max_val:

  Numeric. Maximum value of the profile range.

- inc:

  Numeric. Increment between profile values. Default `0.05`.

- do_par:

  Logical. Whether to use parallel processing. Default `FALSE`.

- n_cores:

  Integer. Number of parallel workers. If `NULL` (default),
  `parallel::detectCores() - 1` is used.

## Value

A named list containing one dataframe per likelihood component, each
with a `prof_val` column indicating the profiled parameter value, plus
`agg_nLL` which aggregates all components across their respective
dimensions. Components include:

- Scalar penalties and priors:

  `jnLL_df`, `rec_nLL_df`, `M_nLL_df`, `sel_nLL_df`, `rec_prop_nLL_df`,
  `Movement_nLL_df`, `h_nLL_df`, `TagRep_nLL_df`, `Fmort_nLL_df`,
  `fish_q_nLL_df`, `srv_q_nLL_df`.

- Pooled data likelihoods:

  `Catch_nLL_df` `[Region × Year × Seas × Fleet]`, `FishIdx_nLL_df` and
  `SrvIdx_nLL_df` `[Region × Year × Seas × Fleet]`, `FishAge_nLL_df`,
  `FishLen_nLL_df`, `SrvAge_nLL_df`, `SrvLen_nLL_df`
  `[Region × Year × Seas × Sex × Fleet]`, `conv_fish_tag_nLL_df`
  `[Recap_Year × Recap_Seas × Tag_Cohort × Region × Fleet]`.

- Population-specific data likelihoods:

  `Catch_pop_nLL_df`, `FishIdx_pop_nLL_df`, `SrvIdx_pop_nLL_df`
  `[Pop × Region × Year × Seas × Fleet]`, `FishAge_pop_nLL_df`,
  `FishLen_pop_nLL_df`, `SrvAge_pop_nLL_df`, `SrvLen_pop_nLL_df`
  `[Pop × Region × Year × Seas × Sex × Fleet]`.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)
