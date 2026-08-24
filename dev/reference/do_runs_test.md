# Runs Test for Residual Randomness

Performs a nonparametric runs test to evaluate whether a sequence of
residuals is randomly distributed around a reference mean. The function
also computes three-sigma control limits used to identify potential
residual outliers.

## Usage

``` r
do_runs_test(x, type = NULL, mixing = "two.sided")
```

## Arguments

- x:

  Numeric vector of residuals.

- type:

  Character string specifying the assumed mean of the residuals. If
  `"resid"` (default), the residual mean is assumed to be zero.
  Otherwise, the empirical mean of `x` is used.

- mixing:

  Character string specifying the alternative hypothesis for the runs
  test:

  - `"two.sided"`: tests for both positive and negative autocorrelation
    (default).

  - `"less"`: left-tailed test detecting positive autocorrelation.

## Value

A list containing:

- `sig3lim`: Numeric vector of length two giving the lower and upper
  three-sigma control limits for the residuals.

- `p.runs`: P-value from the runs test for randomness.

A small p-value (e.g., `< 0.05`) indicates evidence that the residual
sequence is not random and may exhibit autocorrelation or other
systematic patterns.

## Details

The runs test evaluates whether residuals exhibit non-random structure
(e.g., positive or negative autocorrelation).

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
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
