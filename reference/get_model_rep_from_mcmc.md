# Extract model report from MCMC posterior samples

This function collapses MCMC chains from an RTMB/ADNUTS object,
generates model reports for each posterior draw, and extracts specified
components of the report.

## Usage

``` r
get_model_rep_from_mcmc(rtmb_obj, adnuts_obj, what, n_cores)
```

## Arguments

- rtmb_obj:

  An RTMB object created via \`ADFun\`.

- adnuts_obj:

  An \`adnuts\` object containing MCMC samples.

- what:

  Character vector specifying the names of components in the model
  report to extract.

- n_cores:

  Number of cores to use

## Value

A named list of \`data.table\`s, one for each element in \`what\`. Each
table contains the melted report component across all posterior samples,
with an additional column \`posterior_sample\` indicating the MCMC draw
index.

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
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model_reports <- get_model_rep_from_mcmc(rtmb_obj, adnuts_obj,
                                         what = c("SSB", "Rec"))
} # }
```
