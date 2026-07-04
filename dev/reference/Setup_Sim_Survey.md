# Set up survey parameterisation for the operating model simulation

Populates `sim_list` with all survey-related inputs needed by the
operating model: catchability, selectivity, survey timing, index type,
and age/length composition likelihood settings including overdispersion
and correlation parameters. Must be called after
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

## Usage

``` r
Setup_Sim_Survey(
  sim_list,
  srv_sel_input,
  ObsSrvIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_srv_fleets)),
  ObsSrvIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets)),
  srv_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_srv_fleets, sim_list$n_sims)),
  t_srv = array(1, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets)),
  srv_idx_type = array(1, dim = c(sim_list$n_srv_fleets)),
  comp_srvage_like = rep(0, sim_list$n_srv_fleets),
  ISS_SrvAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
  ln_SrvAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets)),
  ln_SrvAge_theta_agg = rep(log(1), sim_list$n_srv_fleets),
  SrvAge_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
  SrvAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets, 2)),
  SrvAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
  comp_srvlen_like = rep(0, sim_list$n_srv_fleets),
  ISS_SrvLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
  ln_SrvLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets)),
  ln_SrvLen_theta_agg = rep(log(1), sim_list$n_srv_fleets),
  SrvLen_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
  SrvLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets, 2)),
  SrvLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
  comp_srvage_pop_like = rep(0, sim_list$n_srv_fleets),
  ISS_SrvAgeComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets,
    sim_list$n_sims)),
  ln_SrvAge_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_srv_fleets)),
  ln_SrvAge_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
  SrvAge_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
  SrvAge_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
  SrvAgeComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
  comp_srvlen_pop_like = rep(0, sim_list$n_srv_fleets),
  ISS_SrvLenComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets,
    sim_list$n_sims)),
  ln_SrvLen_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_srv_fleets)),
  ln_SrvLen_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
  SrvLen_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
  SrvLen_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
  SrvLenComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets))
)
```

## Arguments

- sim_list:

  Simulation list returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

- srv_sel_input:

  Survey selectivity array
  `[n_pop x n_regions x n_yrs x n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]`.
  No default; must be provided.

- ObsSrvIdx_SE:

  Lognormal observation error SD for survey index, array
  `[n_regions × n_yrs × n_seas × n_srv_fleets]`. Default: 0.2.

- ObsSrvIdx_pop_SE:

  As above, but for population-specific indices, array
  `[n_pop × n_regions × n_yrs × n_seas × n_srv_fleets]`.

- srv_q_input:

  Survey catchability array
  `[n_regions × n_yrs × n_srv_fleets × n_sims]`. Default: 1 for all
  cells.

- t_srv:

  Survey timing as fraction of year or season, array
  `[n_regions × n_seas × n_srv_fleets]`. Default: 1.

- srv_idx_type:

  Integer vector `[n_srv_fleets]` specifying survey index type. Default:
  all 1 (biomass). Options: 0/“abd” (abundance), 1/“biom” (biomass).

- comp_srvage_like:

  Integer or character vector `[n_srv_fleets]` specifying likelihood for
  survey age compositions. Default: all 0 (multinomial). Options:
  0/“Multinomial”, 1/“Dirichlet-Multinomial”, 2/“iid-Logistic-Normal”,
  3/“1d-Logistic-Normal”, 4/“2d-Logistic-Normal”.

- ISS_SrvAgeComps:

  Array `[n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]`
  of sample sizes or overdispersion for survey age compositions.
  Default: 100.

- ln_SrvAge_theta:

  Log-scale overdispersion array `[n_regions × n_sexes × n_srv_fleets]`.
  Used for likelihoods 1–4. Default: log(1).

- ln_SrvAge_theta_agg:

  Log-scale overdispersion for aggregated survey age compositions,
  vector `[n_srv_fleets]`. Default: log(1).

- SrvAge_corr_pars_agg:

  Vector `[n_srv_fleets]` for aggregated survey age correlations. Only
  for likelihood 3. Default: 0.01.

- SrvAge_corr_pars:

  Correlation parameters array
  `[n_regions × n_sexes × n_srv_fleets × 2]` (age AR1, sex). Only for
  likelihoods 3–4. Default: 0.01.

- SrvAgeComps_Type:

  Array `[n_yrs × n_srv_fleets]` specifying composition structure.
  Default: 2 (split by region, joint sexes). Options: 0/“agg”,
  1/“spltRspltS”, 2/“spltRjntS”, 999/“none”.

- comp_srvlen_like:

  Integer or character vector `[n_srv_fleets]` specifying likelihood for
  survey length compositions. Default: all 0.

- ISS_SrvLenComps:

  Array `[n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]`
  of sample sizes or overdispersion for survey length compositions.
  Default: 100.

- ln_SrvLen_theta:

  Log-scale overdispersion array `[n_regions × n_sexes × n_srv_fleets]`.
  Default: log(1).

- ln_SrvLen_theta_agg:

  Vector `[n_srv_fleets]` for aggregated length composition
  overdispersion. Default: log(1).

- SrvLen_corr_pars_agg:

  Vector `[n_srv_fleets]` for aggregated length composition
  correlations. Default: 0.01.

- SrvLen_corr_pars:

  Array `[n_regions × n_sexes × n_srv_fleets × 2]` correlation
  parameters for length comps. Default: 0.01.

- SrvLenComps_Type:

  Array `[n_yrs × n_srv_fleets]` specifying length composition
  structure. Default: 2.

- comp_srvage_pop_like:

  Integer or character vector `[n_srv_fleets]` specifying likelihood for
  population-specific survey age compositions. Default: all 0.

- ISS_SrvAgeComps_pop:

  Array
  `[n_pop × n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]`
  of population-specific sample sizes or overdispersion. Default: 100.

- ln_SrvAge_pop_theta:

  Log-scale overdispersion array
  `[n_pop × n_regions × n_sexes × n_srv_fleets]`. Default: log(1).

- ln_SrvAge_pop_theta_agg:

  Array `[n_pop × n_srv_fleets]` for aggregated population-specific
  overdispersion. Default: log(1).

- SrvAge_pop_corr_pars:

  Array `[n_pop × n_regions × n_sexes × n_srv_fleets × 2]` correlation
  parameters (age AR1, sex) for population-specific age compositions.
  Default: 0.01.

- SrvAge_pop_corr_pars_agg:

  Array `[n_pop × n_srv_fleets]` for aggregated population-specific age
  correlations. Default: 0.01.

- SrvAgeComps_pop_Type:

  Array `[n_yrs × n_srv_fleets]` specifying population-specific age
  composition structure. Default: 2.

- comp_srvlen_pop_like:

  Integer or character vector `[n_srv_fleets]` specifying likelihood for
  population-specific survey length compositions. Default: all 0.

- ISS_SrvLenComps_pop:

  Array
  `[n_pop × n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]`
  of population-specific sample sizes or overdispersion. Default: 100.

- ln_SrvLen_pop_theta:

  Array `[n_pop × n_regions × n_sexes × n_srv_fleets]` log-scale
  overdispersion for population-specific lengths. Default: log(1).

- ln_SrvLen_pop_theta_agg:

  Array `[n_pop × n_srv_fleets]` for aggregated population-specific
  length overdispersion. Default: log(1).

- SrvLen_pop_corr_pars:

  Array `[n_pop × n_regions × n_sexes × n_srv_fleets × 2]` correlation
  parameters for population-specific length comps. Default: 0.01.

- SrvLen_pop_corr_pars_agg:

  Array `[n_pop × n_srv_fleets]` for aggregated population-specific
  length correlations. Default: 0.01.

- SrvLenComps_pop_Type:

  Array `[n_yrs × n_srv_fleets]` specifying population-specific length
  composition structure. Default: 2.

## Value

The input `sim_list` with survey-related fields appended: `$srv_sel`,
`$srv_q`, `$ObsSrvIdx_SE`, `$ObsSrvIdx_pop_SE`, `$t_srv`,
`$srv_idx_type`, `$comp_srvage_like`, `$ISS_SrvAgeComps`,
`$ln_SrvAge_theta`, `$ln_SrvAge_theta_agg`, `$SrvAge_corr_pars_agg`,
`$SrvAge_corr_pars`, `$SrvAgeComps_Type`, `$comp_srvlen_like`,
`$ISS_SrvLenComps`, `$ln_SrvLen_theta`, `$ln_SrvLen_theta_agg`,
`$SrvLen_corr_pars_agg`, `$SrvLen_corr_pars`, `$SrvLenComps_Type`,
`$comp_srvage_pop_like`, `$ISS_SrvAgeComps_pop`, `$ln_SrvAge_pop_theta`,
`$ln_SrvAge_pop_theta_agg`, `$SrvAge_pop_corr_pars_agg`,
`$SrvAge_pop_corr_pars`, `$SrvAgeComps_pop_Type`,
`$comp_srvlen_pop_like`, `$ISS_SrvLenComps_pop`, `$ln_SrvLen_pop_theta`,
`$ln_SrvLen_pop_theta_agg`, `$SrvLen_pop_corr_pars_agg`,
`$SrvLen_pop_corr_pars`, `$SrvLenComps_pop_Type`. Character-coded inputs
are converted to integer equivalents before storage.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
