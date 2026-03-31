# Initialize Simulation Dimension Settings

Creates and returns a list of key dimension values used to set up a
simulation or management strategy evaluation (MSE). This list provides
structural information such as number of simulations, years, regions,
ages, fleets, and whether to include a feedback loop.

## Usage

``` r
Setup_Sim_Dim(
  n_sims,
  n_yrs,
  n_regions,
  n_ages,
  n_lens,
  n_obs_ages = n_ages,
  n_sexes,
  n_fish_fleets,
  n_srv_fleets,
  run_feedback = FALSE,
  feedback_start_yr = NULL
)
```

## Arguments

- n_sims:

  Integer. Number of simulation replicates.

- n_yrs:

  Integer. Number of years in the simulation.

- n_regions:

  Integer. Number of modeled regions.

- n_ages:

  Integer. Number of modeled age classes.

- n_lens:

  Integer. Number of modeled length bins.

- n_obs_ages:

  Integer. Number of observed age classes (can differ from `n_ages`,
  default = `n_ages`).

- n_sexes:

  Integer. Number of sexes.

- n_fish_fleets:

  Integer. Number of fishery fleets.

- n_srv_fleets:

  Integer. Number of survey fleets.

- run_feedback:

  Logical. Whether to include a feedback management loop (default =
  `FALSE`).

- feedback_start_yr:

  Integer. First year that feedback is applied (only used if
  `run_feedback = TRUE`).

## Value

A list containing the specified dimension values, with elements:

- `n_sims`, `n_yrs`, `n_regions`, `n_ages`, `n_lens`, `n_obs_ages`,
  `n_sexes`, `n_fish_fleets`, `n_srv_fleets`

- `init_iter` (set internally to `n_ages * 10`)

- `feedback_start_yr`, `run_feedback`

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
