# Get Retrospective Plot

Get Retrospective Plot

## Usage

``` r
get_retrospective_plot(retro_output, Rec_Age)
```

## Arguments

- retro_output:

  Dataframe generated from do_retrospective

- Rec_Age:

  Age in which recruitment occurs

## Value

A retrospective plot of recruitment and SSB in relative and absolute
scales, as well as a retrospective plot of recruitment by cohort (squid
plot)

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
[`get_osa()`](https://chengmatt.github.io/SPoRC/reference/get_osa.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# do retrospective
retro <- do_retrospective(n_retro = 7, # number of retro peels to run
data = data, # rtmb data
parameters = parameters, # rtmb parameters
mapping = mapping, # rtmb mapping
random = NULL, # if random effects are used
do_par = TRUE, # whether or not to parralleize
n_cores = 7, # if parallel, number of cores to use
do_francis = F, # if we want tod o Francis
n_francis_iter = NULL # Number of francis iterations to do
)
get_retrospective_plot(retro, Rec_Age = 2)
} # }
```
