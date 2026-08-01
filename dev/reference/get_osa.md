# Compute OSA residuals for composition data

Formats observed and expected composition data and calculates
one-step-ahead (OSA) residuals using multinomial, Dirichlet-multinomial,
or logistic-normal likelihoods. This function is the main interface for
residual diagnostics, internally calling \[run_external_comp_osa()\] to
perform the residual calculations.

## Usage

``` r
get_osa(
  obs_mat = NULL,
  exp_mat = NULL,
  N = NULL,
  DM_theta = NULL,
  LN_Sigma = NULL,
  years = NULL,
  seas = NULL,
  fleet = NULL,
  bins = NULL,
  comp_type = NULL,
  bin_label = NULL,
  comp_like = 0,
  addtocomp = 0,
  model = NULL,
  data = NULL,
  comp_source = NULL,
  index_source = NULL,
  family = "discrete",
  pop = FALSE,
  discard = FALSE,
  tag = FALSE,
  osa_method = NULL,
  parallel = FALSE
)
```

## Arguments

- obs_mat:

  Array of observed compositions, dimensioned by
  `[region, year, bin, sex, fleet]`. May contain `NA`s, which are
  removed when filtering by `years`.

- exp_mat:

  Array of expected compositions, dimensioned the same as `obs_mat`. May
  contain `NA`s, which are removed when filtering by `years`.

- N:

  Input (or effective if Multinomial) sample size. Dimensions depend on
  `comp_type`:

  - `comp_type = 0` (aggregated): vector of length `n_years`.

  - `comp_type = 1` (split by region and sex): array
    `[n_regions, n_years, n_sexes]`.

  - `comp_type = 2` (split by region, joint by sex): matrix
    `[n_regions, n_years]`.

  For years without data, users can simply input an NA or any abritary
  number (it gets filtered out within the function).

- DM_theta:

  Dirichlet-multinomial overdispersion parameter(s). Dimensions must
  match `N`:

  - aggregated: scalar

  - split by sex: matrix `[n_regions, n_sexes]`

  - joint by sex: vector of length `n_regions`

- LN_Sigma:

  Logistic-normal covariance matrix. Dimensions depend on `comp_type`:

  - aggregated: matrix `[n_bins, n_bins]`

  - split by region and sex: array
    `[n_regions, n_bins, n_bins, n_sexes]`

  - joint by sex: array `[n_regions, n_bins, n_bins]`

  Use \[get_logistN_Sigma()\] to help construct this input.

- years:

  Vector of years to filter to if composition type is aggregated (0).
  Otherwise, this expects a list where each list element is a vector of
  years for each region where compositions are available for use (split
  by region and sex, or split by region, joint by sex).

- seas:

  Season index

- fleet:

  Fleet identifier (character or numeric) to filter to.

- bins:

  Vector of age or length bin labels corresponding to the composition
  categories.

- comp_type:

  Integer specifying how compositions are structured:

  - 0 = aggregated across regions and sexes

  - 1 = split by region and sex

  - 2 = split by region, joint by sex

- bin_label:

  Character label describing whether bins represent ages or lengths.

- comp_like:

  Integer specifying the likelihood type (defaults to 0):

  - 0 = multinomial

  - 1 = Dirichlet-multinomial

  - 2-4 = logistic-normal variants

- addtocomp:

  Constant that is added to compositions

- model:

  A fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
  (built with `do_internal_comp_osa = TRUE` or
  `do_internal_conv_tag_osa = TRUE`). Supplying `model` switches
  `get_osa()` from the default external (post-hoc, compResidual-based)
  path to the internal path, which calls
  [`RTMB::oneStepPredict()`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
  directly on the model's internally tracked OSA vector. All `obs_mat`/
  `exp_mat`/`N`/`DM_theta`/`LN_Sigma`/`years`/ `comp_type`/`comp_like`
  arguments above are ignored in this mode.

- data:

  The model `data` list (e.g. `input_list$data`) used to build `model`.
  Required when `model` is supplied.

- comp_source:

  One of `"FishAge"`, `"FishLen"`, `"SrvAge"`, `"SrvLen"`, identifying
  which composition data source to pull internal OSA residuals for.
  Required when `model` is supplied and `index_source` is `NULL` and
  `tag = FALSE`.

- index_source:

  One of `"Catch"`, `"Discard"`, `"FishIdx"`, `"SrvIdx"`, identifying
  which continuous (log-normal) index-type data source to pull internal
  OSA residuals for. When supplied, takes precedence over
  `comp_source`/`tag`. Only used when `model` is supplied.

- family:

  Character, `"discrete"` or `"continuous"`; which of the two
  internally-tracked OSA vectors to use for `comp_source` (a source can
  have both, e.g. some fleets multinomial and others logistic-normal).
  Only used when `model` is supplied, `tag = FALSE`, and `index_source`
  is `NULL`.

- pop:

  Logical; population-specific composition or index source. Only used
  when `model` is supplied and `tag = FALSE`. Default `FALSE`.

- discard:

  Logical; discard composition source (only valid for
  `comp_source %in% c("FishAge","FishLen")`). Only used when `model` is
  supplied, `tag = FALSE`, and `index_source` is `NULL`. Default
  `FALSE`.

- tag:

  Logical; if `TRUE` (and `model` is supplied, and `index_source` is
  `NULL`), compute internal OSA residuals for conventional tag recapture
  data instead of composition data. Default `FALSE`.

- osa_method:

  Optional override for
  [`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)'s
  `method`, used only in internal mode. Must be one of
  `"oneStepGeneric"`, `"oneStepGaussianOffMode"`, or
  `"oneStepGaussian"`; the `"cdf"` method is not permitted (it is
  numerically fragile for the discrete likelihoods used here). Defaults
  to `"oneStepGeneric"` for discrete families/tags and
  `"oneStepGaussianOffMode"` for continuous (logistic-normal)
  composition families.

- parallel:

  Whether or not to parallelize OSA computation in internal mode.
  Defaults to `FALSE`.

## Value

A list with one element:

- res:

  Data frame of OSA residuals. Columns include: `fleet`, `index_label`,
  `year`, `index`, `resid`, `region`, `seas`, `sex`, and `comp_type`
  (composition sources, external or internal); `fleet`, `region`,
  `cohort`, `release_year`/`release_region`/`release_season`,
  `recovery_year`/`recovery_season`, `years_at_liberty`, `resid`, and
  `comp_type = "Tag"` (`tag = TRUE`); or `fleet`, `region`, `year`,
  `season`, `pop`, `resid`, and `idx_type` set to `index_source`
  (`index_source` supplied.

## Details

When computing OSA residuals for population-specific composition data,
slice the leading population dimension from the `obs_mat` and `exp_mat`
arrays before passing them to this function. For example, to compute
residuals for population `p`:


    get_osa(obs_mat = Obs_FishAge_pop_mat[p,,,,,,],
            exp_mat = Pred_FishAge_pop_mat[p,,,,,,],
            ...)

Population-specific composition arrays returned by
[`get_comp_prop`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md)
are dimensioned
`[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fleets]`.
Slicing on `p` yields a 6D array matching the expected input dimensions.

For internal (model-based) OSA residuals, fit the model with
`do_internal_comp_osa = TRUE` and/or `do_internal_conv_tag_osa = TRUE`
(set in
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)),
then call, e.g.:


    get_osa(model = fitted_obj, data = input_list$data, comp_source = "FishAge",
            family = "discrete", bins = input_list$data$ages, bin_label = "Age")
    get_osa(model = fitted_obj, data = input_list$data, tag = TRUE)
    get_osa(model = fitted_obj, data = input_list$data, index_source = "SrvIdx")

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)
