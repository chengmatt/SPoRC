# Get Retrospective Plot

Generates three retrospective diagnostic plots from a set of peeled
model runs: a relative-difference plot with Mohn's rho, an
absolute-scale trajectory plot, and a cohort-tracking squid plot for
recruitment. Together these plots diagnose systematic retrospective
patterns in SSB and recruitment estimates across sequential data
removals.

## Usage

``` r
get_retrospective_plot(retro_output, Rec_Age)
```

## Arguments

- retro_output:

  Data frame produced by `do_retrospective`, containing columns `Year`,
  `value`, `peel`, `Type` (e.g. `"SSB"`, `"Recruitment"`), `Pop`, and
  `Region`. Rows with `peel = 0` represent the full-data terminal run;
  rows with `peel > 0` represent sequential data removals. Rows with
  `value = 0` are excluded before plotting.

- Rec_Age:

  Integer giving the age at recruitment (e.g. `2` for age-2
  recruitment). Used to convert model year to cohort year in the squid
  plot (`cohort = Year - Rec_Age`).

## Value

A list of three `ggplot` objects:

- \[\[1\]\] retro_plot:

  Relative-difference plot. Each line shows the proportional deviation
  of a peeled run from the terminal-year estimate at each year, computed
  via `get_retrospective_relative_difference`. Terminal-year points
  (where `peel = max(Year) - Year`) are highlighted. Mohn's rho (mean
  relative difference across terminal-year points) is annotated on each
  facet panel. Faceted by Region × (Type + Population); lines and points
  colored by peel year on a continuous viridis scale.

- \[\[2\]\] abs_retro_plot:

  Absolute-scale trajectory plot. Peeled runs (`peel > 0`) are drawn as
  solid lines colored by peel; the full-data run (`peel = 0`) is
  overlaid as a dashed black line. Faceted by Region × (Type +
  Population) with free y-scales.

- \[\[3\]\] squid_plot:

  Cohort-tracking (squid) plot for recruitment. Restricted to the 10
  most recent cohorts. X-axis shows years since the cohort was last
  estimated (`terminal - Year - 1`), y-axis shows recruitment value, and
  lines are grouped and colored by cohort year. Faceted by Population ×
  Region with free y-scales.

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
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  retro <- do_retrospective(
    n_retro        = 7,
    data           = data,
    parameters     = parameters,
    mapping        = mapping,
    random         = NULL,
    do_par         = TRUE,
    n_cores        = 7,
    do_francis     = FALSE,
    n_francis_iter = NULL
  )
  plots <- get_retrospective_plot(retro, Rec_Age = 2)
  plots[[1]]  # relative difference + Mohn's rho
  plots[[2]]  # absolute trajectories
  plots[[3]]  # squid plot
} # }
```
