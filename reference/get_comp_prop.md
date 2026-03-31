# Gets composition data proportions normalized according to the assessment specifications from RTMB

Gets composition data proportions normalized according to the assessment
specifications from RTMB

## Usage

``` r
get_comp_prop(data, rep, age_labels, len_labels, year_labels)
```

## Arguments

- data:

  list of data inputs

- rep:

  report file from RTMB

- age_labels:

  vector of observed age labels in assessment

- len_labels:

  vector of length labels in assessment

- year_labels:

  vector of years

## Value

List of fishery age, lengths, survey age, lengths dataframe as well as
in matrix form (dimensioned by region, year, bin, sex, fleet)

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_catch_fits_plot.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits.md),
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
comp_props <- get_comp_prop(data = data, rep = rep,
age_labels = 2:31, len_labels = seq(41, 99, 2), year_labels = 1960:2024)
comp_props$Fishery_Ages %>%
  filter(Fleet == 1, Sex == 1) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs)) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Year, ncol = 3)

  comp_props$Survey_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(lwr_obs = quantile(obs, 0.1),
              upr_obs = quantile(obs, 0.9),
              lwr_pred = quantile(pred, 0.1),
              upr_pred = quantile(pred, 0.9),
              obs = mean(obs),
              pred = mean(pred)) %>%
    ggplot() +
    geom_line(mapping = aes(x = Age, y = obs,
    color = 'Obs', lty = 'Obs'), lwd = 1.3) +
    geom_ribbon(mapping = aes(x = Age, y = obs,
    ymin = lwr_obs, ymax = upr_obs, fill = 'Obs'), alpha = 0.3) +
    geom_line(mapping = aes(x = Age, y = pred,
     color = 'Pred', lty = 'Pred'), lwd = 1.3) +
    geom_ribbon(mapping = aes(x = Age, y = pred,
     ymin = lwr_pred, ymax = upr_pred, fill = 'Pred'), alpha = 0.3) +
    facet_grid(Region~Fleet, labeller = labeller(
      Region = c('1' = "Region 1"),
      Fleet = c('1' = 'Domestic LL Survey', '3' = 'JP LL Survey')
    )) +
    labs(x = 'Age', y = 'Proportion', color = '', linetype = '', fill = '') +
    theme_bw(base_size = 20) +
    theme(legend.position = 'top')
} # }
```
