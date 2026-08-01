# Extract Index Fit Results

Generates a tidy dataframe of observed and predicted survey and fishery
indices from a fitted RTMB model, including standard errors, confidence
intervals, a raw log-scale (Pearson-style) residual, and catchability
blocks. Both pooled and population-specific indices are returned when
the corresponding `Use*_pop` flags contain any ones.

## Usage

``` r
get_idx_fits(data, rep, year_labs)
```

## Arguments

- data:

  List. Input data used in the RTMB model. Must contain `ObsSrvIdx`,
  `ObsSrvIdx_SE`, `ObsFishIdx`, `ObsFishIdx_SE`, `Wt_SrvIdx`,
  `Wt_FishIdx`, `srv_q_blocks`, `fish_q_blocks`, `UseFishIdx`, and
  `UseSrvIdx`. For population-specific indices, also requires
  `ObsSrvIdx_pop`, `ObsSrvIdx_pop_SE`, `ObsFishIdx_pop`,
  `ObsFishIdx_pop_SE`, `Wt_SrvIdx_pop`, `Wt_FishIdx_pop`,
  `UseSrvIdx_pop`, and `UseFishIdx_pop`.

- rep:

  List. RTMB report output containing `PredSrvIdx` and `PredFishIdx`,
  both dimensioned `[n_pop × n_regions × n_years × n_seas × n_fleets]`.
  Pooled indices are obtained by summing across the population
  dimension; population-specific indices use each population slice
  directly.

- year_labs:

  Vector. Year labels assigned to the year dimension of predicted and
  observed index arrays.

## Value

A dataframe containing pooled and, when active, population-specific
survey and fishery indices with the following columns:

- `Region`:

  Region label (prefixed with `"Region"`).

- `Year`:

  Year.

- `Seas`:

  Season.

- `Fleet`:

  Fleet identifier.

- `Type`:

  One of `"Survey"`, `"Fishery"`, `"Pop Survey"`, or `"Pop Fishery"`.

- `obs`:

  Observed index value.

- `value`:

  Predicted index value.

- `se`:

  Standard error of the observed index (weight-adjusted).

- `lci`, `uci`:

  95% log-normal confidence interval for the observed index.

- `q_block`:

  Catchability block identifier.

- `resid`:

  Log-scale residual (\\\log(\text{obs}) - \log(\text{predicted})\\).

- `Category`:

  Combined label of Type, population (for pop rows), Fleet, Season, and
  Q-block.

## Details

The `resid` column here is the simple log-scale residual
\\\log(\text{obs}) - \log(\text{predicted})\\, *not* a one-step-ahead
(OSA) residual. For properly decorrelated OSA index residuals (with
QQ-plots and SDNR diagnostics via
[`plot_resids`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)),
use
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
with `index_source = `. The observed-vs- predicted *fit* plot built from
this function's output is
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
