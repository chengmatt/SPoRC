# Set up simulation containers and inputs for biological parameters

Set up simulation containers and inputs for biological parameters

## Usage

``` r
Setup_Sim_Biologicals(
  natmort_input,
  WAA_input,
  WAA_fish_input,
  WAA_srv_input,
  MatAA_input,
  AgeingError_input = NULL,
  SizeAgeTrans_input = NULL,
  sim_list
)
```

## Arguments

- natmort_input:

  Natural mortality array \[n_regions × n_yrs × n_ages × n_sexes ×
  n_sims\]

- WAA_input:

  Spawning weight-at-age array \[n_regions × n_yrs × n_ages × n_sexes ×
  n_sims\]

- WAA_fish_input:

  Fishery weight-at-age array \[n_regions × n_yrs × n_ages × n_sexes ×
  n_sims\]

- WAA_srv_input:

  Survey weight-at-age array \[n_regions × n_yrs × n_ages × n_sexes ×
  n_sims\]

- MatAA_input:

  Maturity-at-age array \[n_regions × n_yrs × n_ages × n_sexes ×
  n_sims\]

- AgeingError_input:

  Ageing error matrix \[n_regions × n_model_ages × n_obs_ages × n_sims\]

- SizeAgeTrans_input:

  Size-age transition matrix \[n_regions × n_yrs × n_lens × n_ages ×
  n_sexes x n_sims\]

- sim_list:

  Simulation list object from \`Setup_Sim_Dim()\`

## See also

Other Simulation Setup:
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
