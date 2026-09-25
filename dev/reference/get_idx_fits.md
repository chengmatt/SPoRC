# Extract Index Fit Results

A tidy data frame of the observed and predicted survey and fishery
indices from a fitted model, with standard errors, confidence intervals,
a log-scale residual and the catchability blocks. Population-specific
indices are included when their `Use*_pop` flags hold any ones.

## Usage

``` r
get_idx_fits(data, rep, year_labs)
```

## Arguments

- data:

  Input data used in the RTMB model, holding `ObsSrvIdx`,
  `ObsSrvIdx_SE`, `ObsFishIdx`, `ObsFishIdx_SE`, `Wt_SrvIdx`,
  `Wt_FishIdx`, `srv_q_blocks`, `fish_q_blocks`, `UseFishIdx` and
  `UseSrvIdx`, plus the `_pop` counterparts for the population-specific
  indices.

- rep:

  RTMB report holding `PredSrvIdx` and `PredFishIdx`, both
  `[n_pop × n_regions × n_years × n_seas × n_fleets]`. The pooled
  indices are summed across populations; the population-specific ones
  read each slice directly.

- year_labs:

  Year labels for the year dim of the index arrays.

## Value

A data frame with `Region`, `Year`, `Seas`, `Fleet`, `Type` (`"Survey"`,
`"Fishery"`, `"Pop Survey"` or `"Pop Fishery"`), the observed value
`obs` and predicted `value`, the weight-adjusted `se`, the 95% lognormal
interval `lci` and `uci`, the `q_block`, the log-scale `resid`, and
`Category`, which combines the type, population, fleet, season and q
block.

## Details

The `resid` column is the plain log-scale residual \\\log(\text{obs}) -
\log(\text{predicted})\\, not a one-step-ahead residual. For
decorrelated OSA index residuals, with the QQ plots and SDNR from
[`plot_resids`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md),
use
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
with `index_source`. The fit plot built from this output is
[`get_idx_fits_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md).

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md),
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
idx_fits <- get_idx_fits(
  data = data,
  rep = rep,
  year_labs = seq(1960, 2024, 1)
)
} # }
```
