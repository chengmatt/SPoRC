# Set up biological parameter inputs for closed-loop simulation

Sets natural mortality, weight-at-age, maturity-at-age, ageing error and
an optional size-age transition for the operating model. All arrays are
validated against the dimensions in `sim_list`. Call after
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

## Usage

``` r
Setup_Sim_Biologicals(
  sim_list,
  natmort_input,
  WAA_input,
  WAA_fish_input,
  WAA_srv_input,
  MatAA_input,
  AgeingError_input = NULL,
  AgeingError_fish_input = NULL,
  AgeingError_srv_input = NULL,
  SizeAgeTrans_input = NULL,
  SizeAgeTrans_fish_input = NULL,
  SizeAgeTrans_srv_input = NULL
)
```

## Arguments

- sim_list:

  Simulation list returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

- natmort_input:

  Natural mortality array, either
  `[n_pop × n_regions × n_yrs × n_ages × n_sexes × n_sims]` or with
  `n_seas` between years and ages. Values are rates per year in both
  forms, so mortality within a season is the rate times `seasdur`.

- WAA_input:

  Spawning weight-at-age array
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]`.

- WAA_fish_input:

  Fishery weight-at-age array, `WAA_input` with an `n_fish_fleets` dim
  before the simulations.

- WAA_srv_input:

  Survey weight-at-age array, `WAA_input` with an `n_srv_fleets` dim
  before the simulations.

- MatAA_input:

  Maturity-at-age array dimensioned like `WAA_input`, in \\\[0, 1\]\\.
  Maturity at the first age must be exactly `0` under `rec_lag = 0`, set
  through
  [`Setup_Sim_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md).

- AgeingError_input:

  Ageing error array `[n_yrs × n_model_ages × n_obs_ages × n_sims]`,
  each slice row-stochastic. `NULL` (default) builds an identity matrix
  per year and simulation. For observed bins that are a subset of the
  model ages, supply a shifted identity such as
  `diag(1, n_model_ages)[, obs_age_index]` instead.

- AgeingError_fish_input:

  Optional per-fleet ageing error for the fishery fleets,
  `[n_yrs × n_ages × n_obs_ages × n_fish_fleets × n_sims]`, or `NULL`
  (default) to read `AgeingError_input`.

- AgeingError_srv_input:

  As `AgeingError_fish_input` with `n_srv_fleets` in place of
  `n_fish_fleets`.

- SizeAgeTrans_input:

  Size-age transition array
  `[n_pop × n_regions × n_yrs × n_seas × n_lens × n_ages × n_sexes × n_sims]`,
  column-stochastic over ages. Only needed when fitting length
  compositions. Default `NULL`.

- SizeAgeTrans_fish_input, SizeAgeTrans_srv_input:

  Optional per-fleet size-age arrays
  `[n_pop x n_regions x n_yrs x n_seas x n_lens x n_ages x n_sexes x n_fleets x n_sims]`,
  each read at that fleet's own timing and used in place of
  `SizeAgeTrans_input`. The self-test passes the fitted model's own keys
  here when growth was estimated.

## Value

`sim_list` with `$natmort`, `$WAA`, `$WAA_fish`, `$WAA_srv`, `$MatAA`,
`$AgeingError` (an identity matrix when none was supplied) and, when
supplied, `$SizeAgeTrans`.

## See also

Other Simulation Setup:
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
