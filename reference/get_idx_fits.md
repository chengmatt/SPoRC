# Gets index fits results

Gets index fits results

## Usage

``` r
get_idx_fits(data, rep, year_labs)
```

## Arguments

- data:

  Data list fed into RTMB

- rep:

  Report list output from RTMB

- year_labs:

  Year labels to use (vector)

## Value

Fits to indices as a dataframe

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/reference/get_comp_prop.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
idx_fits <- get_idx_fits(data = data, rep = rep,
 year_labs = seq(1960, 2024, 1))

idx_fits <- idx_fits %>%
  mutate(
    Idx = case_when(
      Type == "Fishery" & Year < 1995 ~
      "Japanese Fishery CPUE Index",
      Type == "Fishery" & Year >= 1995 ~
      "Domestic Fishery CPUE Index",
      Type == 'Survey' & Fleet == 1 ~
      "Domestic LL Survey Relative Population Numbers",
      Type == 'Survey' & Fleet == 2 ~
      "GOA Trawl Survey Biomass (kt)",
      Type == 'Survey' & Fleet == 3 ~
       'Japanese LL Survey Relative Population Numbers'
    )
  )
ggplot() +
  geom_line(idx_fits, mapping =
  aes(x = Year, y = value), lwd = 1.3, col = 'red') +
  geom_pointrange(idx_fits, mapping =
  aes(x = Year, y = obs, ymin = lci, ymax = uci), color = 'blue', pch = 1) +
  labs(x = "Year", y = 'Index') +
  theme_bw(base_size = 20) +
  facet_wrap(~Idx, scales = 'free', ncol = 2)

} # }
```
