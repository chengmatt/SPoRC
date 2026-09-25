# Set up recruitment dynamics for the operating model simulation

Sets the stock-recruit form and density dependence, \\R_0\\, steepness,
sex ratio, the recruitment and initial age deviations, seasonal
allocation, spawn timing and the equilibrium initialization. Call after
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

## Usage

``` r
Setup_Sim_Rec(
  sim_list,
  do_recruits_move = 0,
  sexratio_input = array(if (sim_list$n_sexes == 1) 1 else 0.5, dim = c(sim_list$n_pop,
    sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_sims)),
  R0_input = array(15, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_sims)),
  rinit_input = array(15, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims)),
  h_input = array(0.8, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_sims)),
  stray_rate_input = array(0, dim = c(sim_list$n_pop, sim_list$n_yrs, sim_list$n_sims)),
  ln_sigmaR = array(log(1), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
  rec_seas_prop_input = {
     rec_seas_prop = array(0, dim = c(sim_list$n_pop,
    sim_list$n_seas, sim_list$n_sims))
rec_seas_prop[, 1, ] <- 1
     rec_seas_prop

    },
  recruitment_opt = "bh_rec",
  rec_dd = "global",
  init_dd = "global",
  use_rinit = 0,
  init_age_strc = 2,
  spawn_seas = 1,
  t_spawn = 0,
  rec_lag = 1,
  SR_ref_yr = 1,
  Rec_input = NULL,
  ln_InitDevs_input = NULL,
  InitDevs_sex_spec = "est_shared_s",
  RecDevs_model = "iid",
  RecDevs_rho = array(0, dim = c(sim_list$n_pop, sim_list$n_regions)),
  rec_bias_correct = 1
)
```

## Arguments

- sim_list:

  Simulation list returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

- do_recruits_move:

  Integer flag. `0` (default) starts movement at age 2, `1` moves
  recruits from age 1.

- sexratio_input:

  Proportion of recruits per sex, array
  `[n_pop x n_regions x n_yrs x n_sexes x n_sims]`. Default `1` for one
  sex, `0.5` each for two.

- R0_input:

  Unfished equilibrium recruitment array
  `[n_pop x n_regions x n_yrs x n_sims]`. Default `15`.

- rinit_input:

  Equilibrium recruitment used for initialization under `use_rinit = 1`,
  array `[n_pop x n_regions x n_sims]`. Default `15`.

- h_input:

  Steepness array `[n_pop x n_regions x n_yrs x n_sims]`, values in
  \\(0.2, 1)\\. Default `0.8`.

- stray_rate_input:

  Natal-homing stray rate array `[n_pop x n_yrs x n_sims]`, the
  proportion straying from the natal region at spawning. Default `0`.

- ln_sigmaR:

  Log-scale sd of the recruitment deviations, array
  `[2 x n_pop x n_regions]`, index 1 for `ln_InitDevs` and 2 for
  `ln_RecDevs`. Default `log(1)`.

- rec_seas_prop_input:

  Seasonal allocation of annual recruitment, array
  `[n_pop x n_seas x n_sims]` summing to 1 across seasons. Default all
  in season 1. Must be zero before `spawn_seas` when `rec_lag = 0` and
  `spawn_seas > 1`.

- recruitment_opt:

  Recruitment model, default `"bh_rec"`. `0`/`"mean_rec"` has no
  stock-recruit relationship, `1`/`"bh_rec"` is Beverton-Holt and
  requires `rec_dd = "local"` when `n_pop > 1`, and
  `999`/`"resample_from_input"` resamples years from `Rec_input`, using
  historical years as they are and sampling projection years with
  replacement so spatial covariance within a year is kept.

- rec_dd:

  Density dependence for the stock-recruit relationship, default
  `"global"`. `0`/`"local"` gives each region its own \\R_0\\ and
  steepness, required when `n_pop > 1` under `"bh_rec"`; `1`/`"global"`
  pools across regions.

- init_dd:

  Density dependence for equilibrium initialization, same options as
  `rec_dd`. Default `"global"`.

- use_rinit:

  Integer (0/1). Whether `rinit_input` initializes the population
  separately from `R0_input`. Under `0` (default) `rinit_input` is
  ignored.

- init_age_strc:

  Equilibrium initialization method, default `2`. `0`/`"iterative"`
  iterates forward to approximate equilibrium, `1`/`"scalar_no_move"` is
  a scalar geometric series without movement, `2`/`"matrix"` is the
  matrix series with movement, and `3`/`"scalar_plus_only"` moves only
  the plus group. `4`/`"free"` projects no equilibrium: `ln_InitDevs`
  are the initial log numbers-at-age for ages 2 and above, apportioned
  by sex ratio, so the initial condition ignores `init_F_par` and
  `ln_rinit` and any penalty becomes a prior on initial abundance.

- spawn_seas:

  Integer index of the spawning season. Default `1`.

- t_spawn:

  Spawn timing as a fraction of `spawn_seas` elapsed before spawning.
  `0` (default) spawns before any mortality, `1` after all of it.

- rec_lag:

  Integer seasons between spawning and recruitment. `1` (default) uses
  SSB from that many seasons prior, in any season. `0` is age-0
  recruitment on the same year's SSB, so recruits may only enter in
  `spawn_seas` or later and `rec_seas_prop_input` must be zero before
  it.

- SR_ref_yr:

  Integer year index supplying the biological inputs to unfished
  spawning biomass per recruit, and so the curve's scale. Matches the
  estimation model's `SR_ref_yr`. Default `1`.

- Rec_input:

  External recruitment array `[n_pop x n_regions x n_yrs x n_sims]`,
  required under `recruitment_opt = "resample_from_input"`. Projection
  years beyond its length are resampled from historical years with
  replacement. Default `NULL`.

- ln_InitDevs_input:

  Optional log-scale initial age deviations, either
  `[n_pop x n_regions x (n_ages - 1) x n_sims]` for one shared curve or
  `[n_pop x n_regions x (n_ages - 1) x n_sexes x n_sims]` for one per
  sex. The `n_ages - 1` dim excludes the reference age. `NULL` (default)
  draws one shared curve per population and region; pass zeros to start
  in equilibrium.

- InitDevs_sex_spec:

  How initial age deviations are drawn across sexes when
  `ln_InitDevs_input` is not supplied. `"est_shared_s"` (default) draws
  one curve per population or region for every sex, `"est_all"` draws
  each sex its own. Names match
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

- RecDevs_model:

  Process error the deviations are drawn under. `"iid"` (default) is
  independent, `"rw"` a random walk, `"ar1"` reverts toward zero at rate
  `RecDevs_rho`. Year one is drawn at `ln_sigmaR` under the first two
  and from `ln_sigmaR / sqrt(1 - RecDevs_rho^2)` under `"ar1"`. Matches
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

- RecDevs_rho:

  Matrix `[n_pop x n_regions]` of AR1 correlations in \\(-1, 1)\\. Only
  read under `RecDevs_model = "ar1"`. Default zero.

- rec_bias_correct:

  Integer. `1` (default) draws the recruitment and initial age
  deviations as mean-one lognormal multipliers centered at
  \\-\sigma^2/2\\, matching an estimation model with the bias ramp on;
  `0` centers them at zero. A linked cell follows the same switch under
  the arrows.

## Value

`sim_list` with `$recruitment_opt`, `$rec_dd`, `$init_dd`, `$R0`, `$h`,
`$sexratio`, `$ln_sigmaR`, `$rec_seas_prop`, `$spawn_seas`, `$t_spawn`,
`$rec_lag`, `$init_age_strc`, `$do_recruits_move`, `$move_age`,
`$stray_rate`, and optionally `$Rec_input` and `$ln_InitDevs_input`.
Character codes are converted to integers before storage.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
