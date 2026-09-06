# Specify the state-space numbers-at-age process for simulation

Turns on process error in the numbers at age for a simulated population.
The operating model applies the same centered state the estimation model
does: the deterministic mortality and ageing step is computed, then the
numbers are multiplied by \\\exp(\eta)\\ with \\\eta\\ drawn from the
covariance the arguments here describe.

## Usage

``` r
Setup_Sim_NAA_state(
  sim_list,
  NAA_re = "none",
  sigmaNAA = 0.3,
  rho_age = 0,
  rho_year = 0,
  rho_cohort = 0,
  NAA_re_pop = "iid",
  NAA_re_region = "iid",
  NAA_re_sex = "iid",
  NAA_re_season = "iid",
  pop_corr = 0,
  region_corr = 0,
  sex_corr = 0,
  season_corr = 0,
  NAA_re_ages = NULL,
  NAA_re_years = NULL,
  NAA_re_seasons = "annual"
)
```

## Arguments

- sim_list:

  Simulation list from
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

- NAA_re:

  Character. `"none"` (default) leaves the numbers at age deterministic
  past recruitment and the initial age structure. Otherwise one of
  `"iid"`, `"1dar1_a"`, `"1dar1_y"`, `"2dar1"`, `"3dcond"` or
  `"3dmarg"`.

- sigmaNAA:

  Numeric. Conditional standard deviation of the innovations, the same
  quantity `ln_sigmaNAA` holds in the estimation model. Under an
  autoregressive form the marginal standard deviation is larger by
  \\1/\sqrt{1 - \rho^2}\\ per correlated dim.

- rho_age, rho_year, rho_cohort:

  Numeric correlations in \\(-1, 1)\\ over the age, year and cohort
  dims. Only the ones the chosen form reads are used.

- NAA_re_pop, NAA_re_region, NAA_re_sex, NAA_re_season:

  Character, `"iid"` (default) or `"us"`, an unstructured correlation
  across that dim.

- pop_corr, region_corr, sex_corr, season_corr:

  Numeric vectors of length \\n(n-1)/2\\ giving the correlations for
  those dims, ordered as the strict lower triangle is filled by column.
  A single value is recycled.

- NAA_re_ages, NAA_re_years:

  Ages and year indices the state covers. `NULL` (default) uses
  everything from the second onward.

- NAA_re_seasons:

  Seasons the state covers. `"annual"` (default) puts a state at season
  one only, leaving the numbers deterministic between seasons; `"all"`
  puts one at the start of every season, and an integer vector selects
  specific seasons.

## Value

`sim_list` with the state-space settings attached.

## Details

Arguments mirror `Setup_Mod_Biologicals`'s state-space options so a
simulated population and the model fitted to it are written the same
way, which is what makes a self test a like-for-like comparison rather
than a translation.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
