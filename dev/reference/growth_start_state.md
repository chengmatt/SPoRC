# Growth state at the start of a year

Evaluates the curve a stratum starts from and the quantities that stay
fixed over the years: the mean length at the start of every age in the
first year (the plus group adjusted), the CV at age and timing when it
is held at the first year's sizes, and the asymptote.

## Usage

``` r
growth_start_state(
  gp,
  ages,
  growth_A1,
  growth_A2,
  growth_L0,
  growth_cv_type,
  growth_sd_type,
  growth_plus_group,
  growth_L2_asymptote = 0
)
```
