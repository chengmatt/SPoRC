# Construct logistic-normal covariance matrix

Helper function to generate the covariance matrix (\\\Sigma\\) used in
logistic-normal composition models. The structure depends on the
specification of `comp_like`:

- `comp_like = 2`: independent and identically distributed (iid) across
  categories (\\n\_\mathrm{categories}\\).

- `comp_like = 3`: first-order autoregressive (AR1) correlation across
  categories (\\n\_\mathrm{categories}\\).

- `comp_like = 4`: two-dimensional AR1 correlation across categories and
  sexes (\\n\_\mathrm{categories} \times n\_\mathrm{sexes}\\).

## Usage

``` r
get_logistN_Sigma(
  comp_like,
  n_bins,
  n_sexes,
  theta,
  corr_b = NULL,
  corr_s = NULL
)
```

## Arguments

- comp_like:

  Integer specifying the logistic-normal correlation structure:

  - 2 = iid across categories

  - 3 = AR1 across categories

  - 4 = AR1 across categories and sexes

- n_bins:

  Number of composition categories (e.g., age or length bins). For
  `comp_like = 2, 3`, the covariance matrix is dimensioned `n_bins`. For
  `comp_like = 4`, it is dimensioned `n_bins * n_sexes`.

- n_sexes:

  Number of sexes. Required when `comp_like = 4`.

- theta:

  Standard deviation parameter controlling the overall scale of the
  covariance.

- corr_b:

  Correlation parameter across categories, in the interval \\(-1, 1)\\.
  Used when `comp_like = 3` or `4`.

- corr_s:

  Correlation parameter across sexes, in the interval \\(-1, 1)\\. Used
  when `comp_like = 4`.

## Value

A covariance matrix \\\Sigma\\ with dimension:

- `n_bins` (`comp_like = 2, 3`)

- `n_bins * n_sexes` (`comp_like = 4`)

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/reference/fit_model.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/reference/set_data_indicator_unused.md)

## Examples

``` r
if (FALSE) { # \dontrun{
n_cat <- 5
n_sexes <- 2

# iid example (categories only)
get_logistN_Sigma(comp_like = 2, n_bins = n_cat, n_sexes = NULL, theta = 0.5)

# AR1 across categories
get_logistN_Sigma(comp_like = 3, n_bins = n_cat, n_sexes = NULL, theta = 0.5,
                  corr_b = 0.3)

# AR1 across categories and sexes
get_logistN_Sigma(comp_like = 4, n_bins = n_cat, n_sexes = n_sexes, theta = 0.5,
                  corr_b = 0.3, corr_s = 0.2)
} # }
```
