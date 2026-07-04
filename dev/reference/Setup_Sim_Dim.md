# Initialise simulation dimension settings

Creates the foundational `sim_list` object used throughout a closed-loop
simulation or static simulation. All downstream setup functions
([`Setup_Sim_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
etc.) expect a `sim_list` produced by this function. Dimension scalars
are validated on input and stored alongside derived quantities such as
`init_iter` and the `natal_region` mapping.

## Usage

``` r
Setup_Sim_Dim(
  n_sims,
  n_yrs,
  n_pop = 1,
  n_seas = 1,
  n_regions,
  natal_region = NULL,
  n_ages,
  n_lens = NULL,
  n_obs_ages = n_ages,
  n_sexes,
  n_fish_fleets,
  n_srv_fleets,
  seasdur = if (n_seas == 1) 1 else rep(1/n_seas, n_seas),
  run_feedback = FALSE,
  feedback_start_yr = NULL
)
```

## Arguments

- n_sims:

  Positive integer. Number of simulation replicates.

- n_yrs:

  Positive integer. Number of projection years.

- n_pop:

  Positive integer. Number of distinct populations (default `1`).
  Populations share movement and selectivity schedules but can have
  independent stock-recruit relationships and natal regions.

- n_seas:

  Positive integer. Number of seasons within each year (default `1`).

- n_regions:

  Positive integer. Number of spatial regions.

- natal_region:

  Integer vector of length `n_pop` mapping each population to its natal
  region (1-indexed). Defaults are applied when `NULL`:

  - `n_regions == 1`: all populations assigned to region 1.

  - `n_pop == n_regions`: one-to-one mapping (`1:n_pop`).

  - `n_pop == 1`: assigned to region 1.

  Must be supplied explicitly when `n_pop > 1`, `n_pop != n_regions`,
  and `n_regions > 1`, as multiple populations would share a natal
  region and no sensible default exists.

- n_ages:

  Positive integer. Number of modelled age classes.

- n_lens:

  Positive integer. Number of length bins. Set to `NULL` (default) when
  length compositions are not simulated.

- n_obs_ages:

  Positive integer. Number of observed age bins in composition data. Can
  differ from `n_ages` when the plus group or youngest ages are pooled
  differently in observations. Defaults to `n_ages`.

- n_sexes:

  Integer. Number of sexes; must be either `1` (sex-aggregated) or `2`
  (sex-structured).

- n_fish_fleets:

  Positive integer. Number of fishery fleets.

- n_srv_fleets:

  Positive integer. Number of survey fleets.

- seasdur:

  Numeric vector of length `n_seas` giving the duration of each season
  as a fraction of a year. Values must sum to 1. Defaults to `1` when
  `n_seas == 1`, or `rep(1 / n_seas, n_seas)` for equal-length seasons
  otherwise.

- run_feedback:

  Logical. Whether to run a closed-loop feedback MSE in which an
  estimation model and harvest control rule update fishing mortality
  each year. Default `FALSE` (open-loop simulation).

- feedback_start_yr:

  Integer. First year in which feedback is applied. Required when
  `run_feedback = TRUE`; ignored otherwise.

## Value

A named list (`sim_list`) containing:

- `n_sims`, `n_yrs`, `n_pop`, `n_regions`, `n_seas`, `n_ages`,
  `n_obs_ages`, `n_lens`, `n_sexes`, `n_fish_fleets`, `n_srv_fleets`:

  Dimension scalars supplied as inputs.

- `natal_region`:

  Integer vector of length `n_pop` defining natal regions.

- `seasdur`:

  Season durations summing to 1.

- `init_iter`:

  Derived value `n_ages * 10` used to initialise equilibrium conditions.

- `run_feedback`, `feedback_start_yr`:

  Feedback control settings.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
