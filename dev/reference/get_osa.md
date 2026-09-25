# Compute OSA residuals for composition data

One-step-ahead residuals under the multinomial, Dirichlet-multinomial or
logistic-normal likelihoods. The external path formats the observed and
expected compositions and calls \[run_external_comp_osa()\]; supplying
`model` switches to the internal path, which calls
[`RTMB::oneStepPredict()`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
on the model's own tracked OSA vector.

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

- obs_mat, exp_mat:

  Observed and expected composition arrays
  `[region, year, bin, sex, fleet]`. `NA`s are removed when filtering by
  `years`.

- N:

  Input sample size, or the effective one under the multinomial, always
  at the model's full year dim: `years` selects the years used, the same
  way it selects them from `obs_mat`, so `N` is never pre-filtered. It
  is a vector of length `n_years` under `comp_type = 0`, an array
  `[n_regions, n_years, n_sexes]` under `1`, and a matrix
  `[n_regions, n_years]` under `2`. Years without data may hold `NA` or
  anything else, since they are filtered out. This is the `ISS_*Comps`
  array from the model's data list, so it usually passes straight
  through.

- DM_theta:

  Dirichlet-multinomial overdispersion, dimensioned to match `N`: a
  scalar when aggregated, a matrix `[n_regions, n_sexes]` when split by
  sex, and a vector of length `n_regions` when joint by sex.

- LN_Sigma:

  Logistic-normal covariance: a matrix `[n_bins, n_bins]` when
  aggregated, an array `[n_regions, n_bins, n_bins, n_sexes]` when split
  by region and sex, and `[n_regions, n_bins, n_bins]` when joint by
  sex. Use \[get_logistN_Sigma()\] to build it.

- years:

  Years with composition data, either a plain vector used for every
  region or a list with one vector per region. Both forms work for every
  composition type, and a region with no years is skipped.

- seas:

  Season index.

- fleet:

  Fleet identifier, character or numeric, to filter to.

- bins:

  Age or length bin labels for the composition categories.

- comp_type:

  Integer: 0 aggregated across regions and sexes, 1 split by region and
  sex, 2 split by region and joint by sex.

- bin_label:

  Character label saying whether the bins are ages or lengths.

- comp_like:

  Integer likelihood: 0 multinomial (default), 1 Dirichlet-multinomial,
  2-4 the logistic-normal forms.

- addtocomp:

  Constant added to the compositions.

- model:

  A fitted RTMB model from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
  built with `do_internal_comp_osa = TRUE` or
  `do_internal_conv_tag_osa = TRUE`. Supplying it switches to the
  internal path and every argument above is ignored.

- data:

  The model `data` list used to build `model`. Required when `model` is
  supplied.

- comp_source:

  Which composition data source to pull internal residuals for:
  `"FishAge"`, `"FishLen"`, `"SrvAge"` or `"SrvLen"`. Conditional
  age-at-length takes `"Fish_caal"` or `"Srv_caal"`, which return an
  extra `len` column giving the length bin each age composition was
  conditioned on and ignore `family`, CAAL having only the discrete
  likelihoods. Required when `model` is supplied, `index_source` is
  `NULL` and `tag = FALSE`.

- index_source:

  Which continuous index-type data source to pull internal residuals
  for: `"Catch"`, `"Discard"`, `"FishIdx"` or `"SrvIdx"`. Takes
  precedence over `comp_source` and `tag`.

- family:

  `"discrete"` or `"continuous"`, which of the two tracked OSA vectors
  to read for `comp_source`, since a source can have both. Read when
  `model` is supplied, `tag = FALSE` and `index_source` is `NULL`.

- pop:

  Logical, whether the source is population-specific. Read when `model`
  is supplied and `tag = FALSE`. Default `FALSE`.

- discard:

  Logical, whether the source is the discard compositions, valid for
  `comp_source %in% c("FishAge","FishLen")` only. Default `FALSE`.

- tag:

  Logical, `TRUE` to compute internal residuals for conventional tag
  recaptures instead of compositions. Default `FALSE`.

- osa_method:

  Optional override for
  [`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)'s
  `method` in internal mode: `"oneStepGeneric"`,
  `"oneStepGaussianOffMode"` or `"oneStepGaussian"`. The `"cdf"` method
  is not permitted, being numerically fragile for the discrete
  likelihoods. Defaults to `"oneStepGeneric"` for the discrete families
  and tags, and `"oneStepGaussianOffMode"` for the continuous ones.

- parallel:

  Whether to parallelize the internal computation. Default `FALSE`.

## Value

A list with one element, `res`, a data frame of residuals. A composition
source gives `fleet`, `index_label`, `year`, `index`, `resid`, `region`,
`seas`, `sex` and `comp_type`. `tag = TRUE` gives `fleet`, `region`,
`cohort`, the release and recovery year, region and season,
`years_at_liberty`, `resid` and `comp_type = "Tag"`. An `index_source`
gives `fleet`, `region`, `year`, `season`, `pop`, `resid` and
`idx_type`.

## Details

For population-specific compositions, slice the leading population dim
off `obs_mat` and `exp_mat` first. The arrays
[`get_comp_prop`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md)
returns are
`[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fleets]`,
so slicing on `p` gives the 6D array this function expects:


    get_osa(obs_mat = Obs_FishAge_pop_mat[p,,,,,,],
            exp_mat = Pred_FishAge_pop_mat[p,,,,,,],
            ...)

For internal residuals, fit with `do_internal_comp_osa = TRUE` or
`do_internal_conv_tag_osa = TRUE`, set in
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md):


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
