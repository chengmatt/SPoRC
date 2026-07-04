# Extract model report quantities from MCMC posterior samples

Discards warmup iterations, collapses all chains into a single matrix of
posterior draws, evaluates the RTMB report function at each draw in
parallel, and returns the requested report components as tidy
`data.table`s with a `posterior_sample` index column.

## Usage

``` r
get_model_rep_from_mcmc(rtmb_obj, mcmc_obj, what, n_cores)
```

## Arguments

- rtmb_obj:

  An RTMB `ADFun` object with a `$report()` method, as returned by
  `RTMB::MakeADFun`.

- mcmc_obj:

  An `adnuts` or `SparseNUTS` posterior object containing `$samples`
  (array `[n_iter × n_chain × n_param]`) and `$warmup` (number of warmup
  iterations to discard).

- what:

  Character vector. Names of components in the model report (i.e.,
  quantities passed to
  [`RTMB::REPORT`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)
  inside
  [`SPoRC_rtmb`](https://chengmatt.github.io/SPoRC/dev/reference/SPoRC_rtmb.md))
  to extract from each posterior draw.

- n_cores:

  Integer. Number of parallel workers to use via
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html).

## Value

Named list of `data.table`s, one per element of `what`. Each table is
the row-bound result of
[`reshape2::melt`](https://rdrr.io/pkg/reshape2/man/melt.html) applied
to the report component across all post-warmup draws, with an additional
integer column `posterior_sample` identifying the draw index (1 to
`n_iter × n_chain`).

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
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model_reports <- get_model_rep_from_mcmc(
  rtmb_obj, mcmc_obj,
  what = c("SSB", "Rec"),
  n_cores = 4
)
} # }
```
