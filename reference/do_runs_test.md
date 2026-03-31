# Runs test function taken from SS3 diags.

Runs test function taken from SS3 diags.

## Usage

``` r
do_runs_test(x, type = NULL, mixing = "two.sided")
```

## Arguments

- x:

  Vector of residuals

- type:

  Whether to use mean 0 assumption of mean of residuals (default = use
  mean 0)

- mixing:

  Type of test to do, less = left tailed test that detects positive
  autocorrelation, two.sided = two sided test that tests whether there
  is positive and/or negative autocorrealtion. The null is that there
  isn't any, rejecting the null (\<0.05) indictes that there is some
  non-randomness.

## Value

List object with p value and limits for a three-sigma limit - (potential
data outlier, where residual is \> 3 standard deviations away from a
mean of 0)

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/reference/do_retrospective.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/reference/get_comp_prop.md),
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

  unique_idx <- unique(idx_fits$Idx)
  runs_all <- data.frame()
  for(i in 1:length(unique(idx_fits$Idx))) {
    tmp <- idx_fits %>% filter(Idx == unique_idx[i])
    runstest <- do_runs_test(x=as.numeric(tmp$resid),
    type="resid", mixing = "less")
    tmp_runs <- data.frame(p = runstest$p.runs,
     lwr = runstest$sig3lim[1], upr = runstest$sig3lim[2],
     Idx = unique_idx[i])
    runs_all <- rbind(runs_all, tmp_runs)
  } # end i

  ggplot() +
    geom_point(idx_fits, mapping = aes(x = Year, y = resid)) +
    geom_segment(idx_fits, mapping =
    aes(x = Year, xend = Year, y = 0, yend = resid)) +
    geom_smooth(idx_fits, mapping = aes(x = Year, y = resid), se = F) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_hline(runs_all, mapping = aes(yintercept = upr), lty = 2) +
    geom_hline(runs_all, mapping = aes(yintercept = lwr), lty = 2) +
    geom_text(data = runs_all, aes(x = -Inf, y = Inf,
    label = paste("p = ", round(p, 3))), hjust = -0.5,
    vjust = 8.2, size = 7)+
    labs(x = "Year", y = 'Residuals') +
    theme_bw(base_size = 20) +
    facet_wrap(~Idx, scales = 'free', ncol = 2)
} # }
```
