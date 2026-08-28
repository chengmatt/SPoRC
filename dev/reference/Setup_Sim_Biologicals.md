# Set up biological parameter inputs for closed-loop simulation

Populates a simulation list (created by
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md))
with biological arrays needed to run the operating model: natural
mortality, weight-at-age (spawning, fishery, and survey),
maturity-at-age, ageing error, and an optional size-age transition
matrix for length compositions. All arrays must conform to the dimension
structure stored in `sim_list`.

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

  Simulation list object returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
  which defines the dimension sizes used to validate all input arrays.

- natmort_input:

  Natural mortality array with dimensions
  `[n_pop × n_regions × n_yrs × n_ages × n_sexes × n_sims]`. Note:
  natural mortality is not season-specific and therefore lacks an
  `n_seas` dimension.

- WAA_input:

  Spawning weight-at-age array with dimensions
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]`.
  Used to compute spawning stock biomass.

- WAA_fish_input:

  Fishery weight-at-age array with dimensions
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims]`.
  Used to compute fishery biomass and catch in weight.

- WAA_srv_input:

  Survey weight-at-age array with dimensions
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]`.
  Used to compute survey biomass indices.

- MatAA_input:

  Maturity-at-age array with dimensions
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]`.
  Values should be proportions in \\\[0, 1\]\\. When `rec_lag = 0`
  (age-0 recruitment, set via
  [`Setup_Sim_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)),
  maturity at the recruit age (the first age class) must be exactly `0`
  for all populations, regions, years, seasons, and sexes. An error is
  raised otherwise.

- AgeingError_input:

  Ageing error (age-length transition) array with dimensions
  `[n_yrs × n_model_ages × n_obs_ages × n_sims]`, where each
  `[n_model_ages × n_obs_ages]` slice is a row-stochastic matrix mapping
  true modeled ages to observed age bins. If `NULL` (default), an
  identity matrix is constructed for each year and simulation, which
  assumes that modeled and observed age bins are identical in number and
  alignment. **If observed age bins are a subset of modeled ages**
  (e.g., observed ages 2-10 vs. modeled ages 1-10), the default identity
  matrix will cause a dimensional mismatch. In that case, supply a
  shifted identity matrix such as
  `diag(1, n_model_ages)[, obs_age_index]` to correctly drop or collapse
  model ages into observed bins.

- AgeingError_fish_input:

  Optional fleet-specific ageing error for the simulated fishery fleets,
  dimensioned `[n_yrs × n_ages × n_obs_ages × n_fish_fleets × n_sims]`,
  or `NULL` (default) to give every fishery fleet `AgeingError_input`.

- AgeingError_srv_input:

  Optional fleet-specific ageing error for the simulated survey fleets,
  dimensioned `[n_yrs × n_ages × n_obs_ages × n_srv_fleets × n_sims]`,
  or `NULL` (default) to give every survey fleet `AgeingError_input`.

- SizeAgeTrans_input:

  Size-age transition matrix array with dimensions
  `[n_pop × n_regions × n_yrs × n_seas × n_lens × n_ages × n_sexes × n_sims]`.
  Each slice maps age classes to length bins and should be
  column-stochastic (columns sum to 1). Only required when fitting
  length compositions; defaults to `NULL`.

- SizeAgeTrans_fish_input, SizeAgeTrans_srv_input:

  Optional size-age transition arrays per fleet,
  `[n_pop x n_regions x n_yrs x n_seas x n_lens x n_ages x n_sexes x n_fleets x n_sims]`,
  each read at that fleet's own timing. When supplied they are used for
  that fleet type in place of `SizeAgeTrans_input`; the self-test passes
  the fitted model's own keys here when growth was estimated.

## Value

The input `sim_list` with the following fields added or updated:

- `$natmort`:

  Natural mortality array.

- `$WAA`:

  Spawning weight-at-age array.

- `$WAA_fish`:

  Fishery weight-at-age array.

- `$WAA_srv`:

  Survey weight-at-age array.

- `$MatAA`:

  Maturity-at-age array.

- `$AgeingError`:

  Ageing error array (identity matrix if not supplied).

- `$SizeAgeTrans`:

  Size-age transition array (only added if supplied).

## See also

Other Simulation Setup:
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
