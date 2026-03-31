# Set up recruitment dynamics for simulation

Set up recruitment dynamics for simulation

## Usage

``` r
Setup_Sim_Rec(
  do_recruits_move = 0,
  sexratio_input = array(if (sim_list$n_sexes == 1) 1 else 0.5, dim =
    c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_sims)),
  R0_input = array(10, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
  h_input = array(0.8, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
  ln_sigmaR = log(c(1, 1)),
  recruitment_opt = "bh_rec",
  rec_dd = "global",
  init_dd = "global",
  sim_list,
  init_age_strc = 2,
  t_spawn = 0,
  rec_lag = 1,
  Rec_input = NULL,
  ln_InitDevs_input = NULL
)
```

## Arguments

- do_recruits_move:

  Indicator for whether recruits move (default = 0):

  - `0`: No movement

  - `1`: Move

- sexratio_input:

  Sex ratio array \[n_regions × n_yrs × n_sexes × n_sims\] (default = 1
  if one sex, else 0.5 for each sex)

- R0_input:

  Unfished recruitment (R0) array \[n_regions × n_yrs × n_sims\]
  (default = 10)

- h_input:

  Steepness array \[n_regions × n_yrs × n_sims\] (default = 0.8)

- ln_sigmaR:

  Logarithmic standard deviation of recruitment \[2\]: 1st = sigma for
  initial devs, 2nd = sigma for latter devs (default = log(c(1, 1)))

- recruitment_opt:

  Recruitment type (default = "bh_rec"):

  - `"mean_rec"`: Mean recruitment

  - `"bh_rec"`: Beverton-Holt recruitment

  - `"resample_from_input"`: Resampling recruitment years from
    \`Rec_input\` and preserves covariance of recruitment among regions
    if spatially-explicit values are provided

- rec_dd:

  Recruitment density dependence (default = "global"):

  - `"global"`: Shared across regions

  - `"local"`: Region-specific

- init_dd:

  Initial age density dependence (default = "global"):

  - `"global"`: Shared across regions

  - `"local"`: Region-specific

- sim_list:

  Simulation list object from \`Setup_Sim_Dim()\`

- init_age_strc:

  Integer specifying the initialization method for the age structure: -
  0: Iterative solution to equilibrium - 1: Scalar geometric series
  solution w/o movement in any groups (no movement in all groups) - 2:
  Matrix geometric series solution (generalizes scalar solution with
  movement) - 3: Scalar geometric series solution w/o movement only in
  plus group (no movement in plus groups)

- t_spawn:

  Spawn timing fraction of the year (scalar, default = 0)

- rec_lag:

  Recruitment lag (default = 1)

- Rec_input:

  Recruitment array \[n_regions × n_yrs × n_sims\] (default = NULL)

- ln_InitDevs_input:

  Initial deviations \[n_regions × (n_ages-1) × n_sims\] (default =
  NULL)

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
