# Construct a logistic-normal covariance matrix

Builds the covariance matrix \\\Sigma\\ used in logistic-normal
composition likelihoods for a given correlation structure. Three
structures are supported, matching the `comp_like` codes used throughout
SPoRC: iid (`2`), AR(1) across bins (`3`), and AR(1) across bins with
constant correlation across sexes via a Kronecker product (`4`).

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

  Integer. Covariance structure: `2` = iid (diagonal \\\theta^2 I\\),
  `3` = AR(1) across bins (\\\theta^2 C\_{\text{AR1}}\\), `4` =
  Kronecker product of constant sex correlation and AR(1) bin
  correlation (\\\theta^2 (C\_{\text{sex}} \otimes C\_{\text{AR1}})\\).

- n_bins:

  Integer. Number of composition categories (ages or lengths). The
  resulting matrix has dimension `n_bins` for `comp_like` `2` and `3`,
  or `n_bins × n_sexes` for `comp_like = 4`.

- n_sexes:

  Integer. Number of sexes. Required for `comp_like = 4`; ignored
  otherwise.

- theta:

  Numeric. Marginal standard deviation \\\theta \> 0\\ controlling the
  overall scale of \\\Sigma\\.

- corr_b:

  Numeric. AR(1) correlation across bins in \\(-1, 1)\\. Required for
  `comp_like` `3` and `4`; ignored for `comp_like = 2`.

- corr_s:

  Numeric. Constant (exchangeable) correlation across sexes in \\(-1,
  1)\\. Required for `comp_like = 4`; ignored otherwise.

## Value

Numeric covariance matrix \\\Sigma\\ of dimension `n_bins × n_bins`
(`comp_like` `2`, `3`) or `(n_bins × n_sexes) × (n_bins × n_sexes)`
(`comp_like = 4`).

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# iid
get_logistN_Sigma(comp_like = 2, n_bins = 5, n_sexes = NULL, theta = 0.5)

# AR(1) across bins
get_logistN_Sigma(comp_like = 3, n_bins = 5, n_sexes = NULL,
                  theta = 0.5, corr_b = 0.3)

# AR(1) across bins x constant across sexes
get_logistN_Sigma(comp_like = 4, n_bins = 5, n_sexes = 2,
                  theta = 0.5, corr_b = 0.3, corr_s = 0.2)
} # }
```
