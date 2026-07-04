# Set up recruitment dynamics for the operating model simulation

Populates `sim_list` with all recruitment-related inputs needed by the
operating model: stock-recruit relationship type and density-dependence
structure, biological parameters (\\R_0\\, steepness, sex ratio),
recruitment and initial age-structure deviations, seasonal recruitment
allocation, spawn timing, and equilibrium initialisation method. Must be
called after
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
  Rec_input = NULL,
  ln_InitDevs_input = NULL
)
```

## Arguments

- sim_list:

  Simulation list returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

- do_recruits_move:

  Integer flag. `0` = age-1 fish do not move (default); movement begins
  at age 2. `1` = recruits participate in movement from age 1.

- sexratio_input:

  Proportion of recruits assigned to each sex, array
  `[n_pop x n_regions x n_yrs x n_sexes x n_sims]`. Default: `1` when
  `n_sexes = 1`; `0.5` per sex when `n_sexes = 2`.

- R0_input:

  Unfished equilibrium recruitment array
  `[n_pop x n_regions x n_yrs x n_sims]`. Default: `15` for all cells.

- rinit_input:

  Unfished equilibrium recruitment scalar used for population
  initialisation when `use_rinit = 1`, array
  `[n_pop x n_regions x n_sims]`. Ignored when `use_rinit = 0`. Default:
  `15` for all cells.

- h_input:

  Beverton-Holt steepness array `[n_pop x n_regions x n_yrs x n_sims]`.
  Values should be in \\(0.2, 1)\\. Default: `0.8`.

- stray_rate_input:

  Natal-homing stray rate array `[n_pop x n_yrs x n_sims]`. Proportion
  of individuals that stray from their natal region during spawning.
  Default: `0` (No individuals stray).

- ln_sigmaR:

  Log-scale standard deviation of recruitment deviations, array
  `[2 x n_pop x n_regions]`. The first element controls the SD for
  initial age-structure deviations (`ln_InitDevs`); the second controls
  the SD for annual recruitment deviations (`ln_RecDevs`). Default:
  `log(1)` for both.

- rec_seas_prop_input:

  Seasonal allocation of annual recruitment, array
  `[n_pop x n_seas x n_sims]`. Each population's values should sum to 1
  across seasons. Default: all recruitment assigned to season 1.

- recruitment_opt:

  Recruitment model. Default `"bh_rec"`. Options:

  `0`/`"mean_rec"`

  :   Mean recruitment; no stock-recruit relationship.

  `1`/`"bh_rec"`

  :   Beverton-Holt stock-recruit relationship. Requires
      `rec_dd = "local"` when `n_pop > 1`.

  `999`/`"resample_from_input"`

  :   Bootstrap recruitment by resampling years from `Rec_input`.
      Historical years are used as-is; projection years are sampled with
      replacement, preserving spatial covariance among regions within
      each resampled year. Requires `Rec_input` to be non-`NULL`.

- rec_dd:

  Density-dependence structure for the stock-recruit relationship.
  Default `"global"`. Options:

  `0`/`"local"`

  :   Region-specific spawner-recruit relationship; each region has its
      own \\R_0\\ and steepness. Required when `n_pop > 1` and
      `recruitment_opt = "bh_rec"`.

  `1`/`"global"`

  :   Single shared spawner-recruit relationship pooled across regions.

- init_dd:

  Density-dependence structure for equilibrium age-structure
  initialisation. Default `"global"`. Same options as `rec_dd`.

- use_rinit:

  Integer (0/1). Whether a separate initial recruitment scalar
  `rinit_input` is used to initialise the population independently of
  `R0_input`. When `0` (default), `rinit_input` is ignored and
  `R0_input` governs both initialisation and recruitment. When `1`,
  `rinit_input` is used exclusively for equilibrium initialisation and
  `R0_input` governs the recruitment relationship.

- init_age_strc:

  Integer specifying the equilibrium age-structure initialisation
  method. Default `2`. Options:

  `0`/`"iterative"`

  :   Iterates the population forward until approximate equilibrium.
      Slowest but most general.

  `1`/`"scalar_no_move"`

  :   Scalar geometric series solution assuming no movement in any age
      class.

  `2`/`"matrix"`

  :   Matrix geometric series solution that generalises the scalar
      approach to include movement. Recommended default for spatially
      explicit models.

  `3`/`"scalar_plus_only"`

  :   Scalar geometric series solution assuming no movement except in
      the plus group.

- spawn_seas:

  Integer index of the season in which spawning occurs. Default `1`.

- t_spawn:

  Spawn timing as a fraction of the season elapsed before spawning
  within `spawn_seas`. `0` (default) = spawning occurs before any
  mortality is applied in that season; `1` = spawning occurs after all
  mortality.

- rec_lag:

  Integer. Number of seasons between spawning and recruitment of age-1
  fish. Must be \\\geq 1\\; `0` is not permitted. Default `1`.

- Rec_input:

  External recruitment array `[n_pop x n_regions x n_yrs x n_sims]`.
  Required when `recruitment_opt = "resample_from_input"`; projection
  years beyond the length of `Rec_input` are filled by resampling
  historical years with replacement. Ignored for other recruitment
  options. Default `NULL`.

- ln_InitDevs_input:

  Optional log-scale initial age-structure deviations array
  `[n_pop x n_regions x (n_ages - 1) x n_sims]`. The `n_ages - 1`
  dimension excludes the reference age used during initialisation. If
  `NULL` (default), deviations are initialised to zero.

## Value

The input `sim_list` with recruitment-related fields appended:
`$recruitment_opt`, `$rec_dd`, `$init_dd`, `$R0`, `$h`, `$sexratio`,
`$ln_sigmaR`, `$rec_seas_prop`, `$spawn_seas`, `$t_spawn`, `$rec_lag`,
`$init_age_strc`, `$do_recruits_move`, `$move_age`, `$stray_rate`, and
optionally `$Rec_input` and `$ln_InitDevs_input`. Character-coded inputs
are converted to their integer equivalents before storage.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
