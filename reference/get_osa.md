# Compute OSA residuals for composition data

Formats observed and expected composition data and calculates
one-step-ahead (OSA) residuals using multinomial, Dirichlet-multinomial,
or logistic-normal likelihoods. This function is the main interface for
residual diagnostics, internally calling \[run_osa()\] to perform the
residual calculations.

## Usage

``` r
get_osa(
  obs_mat,
  exp_mat,
  N = NULL,
  DM_theta = NULL,
  LN_Sigma = NULL,
  years,
  fleet,
  bins,
  comp_type,
  bin_label,
  comp_like = 0,
  addtocomp = 0
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

  - 2–4 = logistic-normal variants

- addtocomp:

  Constant that is added to compositions

## Value

A list with one element:

- res:

  Data frame of OSA residuals. Columns include: `fleet`, `index_label`,
  `year`, `index`, `resid`, `region`, `sex`, and `comp_type`.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/reference/get_nLL_plot.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/reference/plot_resids.md)
