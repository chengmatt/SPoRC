# Plots OSA residuals from outputs from get_osa. Much of this code is taken from the afscOM package, but with modificaitons to plot features.

Plots OSA residuals from outputs from get_osa. Much of this code is
taken from the afscOM package, but with modificaitons to plot features.

## Usage

``` r
plot_resids(osa_results)
```

## Arguments

- osa_results:

  List object obtained from get_osa, that contains a dataframe of
  residuals and aggregated fits.

## Value

A vareity of plots for OSA residuals (list)

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
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md)

## Examples

``` r
if (FALSE) { # \dontrun{
comp_props <- get_comp_prop(data = data, rep = sabie_rtmb_model$rep,
 age_labels = 2:31, len_labels = seq(41, 99, 2), year_labels = 1960:2024)
plot_resids(get_osa(obs_mat = comp_props$Obs_FishAge_mat,
                    exp_mat = comp_props$Pred_FishAge_mat,
                    N = rep(16.52215, length(1999:2023)),
                    years = which(1960:2024 %in% 1999:2023),
                    LN_Sigma = LN_Sigma,
                    fleet = 1,
                    bins = 2:31,
                    comp_type = 0,
                    comp_like = 0,
                    bin_label = "Age"))
osa_plot <- plot_resids(osa_results)
} # }
```
