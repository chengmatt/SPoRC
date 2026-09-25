# Set up survey parameterization for the operating model simulation

Sets the survey catchability, selectivity, timing, index type and the
age and length composition settings, with their overdispersion and
correlation parameters. Call after
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

## Usage

``` r
Setup_Sim_Survey(
  sim_list,
  srv_sel_input,
  ObsSrvIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_srv_fleets)),
  ln_sigmaSrvIdxAA = array(log(0.2), dim = c(sim_list$n_obs_ages, sim_list$n_sexes,
    sim_list$n_srv_fleets)),
  UseSrvIdxAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets)),
  ObsSrvIdxAA_SE = NULL,
  SrvIdxAA_Type = "spltRaggS",
  SrvIdxAA_LikeType = "lognormal",
  SrvIdxAA_sigma_form = "none",
  use_srv_idx_aa = rep(0, sim_list$n_srv_fleets),
  ObsSrvIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets)),
  srv_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_srv_fleets, sim_list$n_sims)),
  t_srv = array(1, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets)),
  srv_idx_type = array(1, dim = c(sim_list$n_srv_fleets)),
  SrvIdx_LikeType = rep(0, sim_list$n_srv_fleets),
  SrvIdx_seas_Type = NULL,
  SrvIdx_pop_seas_Type = NULL,
  SrvAgeComps_seas_Type = NULL,
  SrvIdx_Cov = NULL,
  UseSrvIdx = NULL,
  comp_srv_caal_like = rep(999, sim_list$n_srv_fleets),
  ISS_Srv_caal = NULL,
  ln_Srv_caal_theta = NULL,
  ln_Srv_caal_theta_agg = NULL,
  Srv_caal_Type = array(999, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
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
  No default.

- ObsSrvIdx_SE, ObsSrvIdx_pop_SE:

  Lognormal observation error sd for the survey indices,
  `[n_regions × n_yrs × n_seas × n_srv_fleets]` with a leading `n_pop`
  for the second. Default 0.2.

- ln_sigmaSrvIdxAA:

  Log-scale observation error for the index at age, \`n_obs_ages x
  n_sexes x n_srv_fleets\`. The sex dim is required.

- UseSrvIdxAA:

  Integer array \`n_regions x n_yrs x n_seas x n_obs_ages x n_sexes x
  n_srv_fleets\`, \`1\` where a survey index at age is drawn, on the
  observed ages \`AgeingError_srv_input\` reads onto. The sex dim is
  required: a data source summed over sexes has its flag in sex slot
  one.

- ObsSrvIdxAA_SE:

  Reported standard errors shaped like \`UseSrvIdxAA\`, read only when
  \`SrvIdxAA_sigma_form\` asks for them.

- SrvIdxAA_Type:

  Which dims each fleet reports separately: \`"agg"\`, \`"spltRaggS"\`
  (default), \`"aggRspltS"\` or \`"spltRspltS"\`.

- SrvIdxAA_LikeType:

  \`"lognormal"\` (default) or \`"normal"\`, per fleet.

- SrvIdxAA_sigma_form:

  Where the observation error comes from: \`"none"\` (default),
  \`"data"\`, \`"est_additive"\` or \`"est_quadrature"\`.

- use_srv_idx_aa:

  Integer vector \`n_srv_fleets\`, \`1\` for fleets whose index at age
  is drawn.

- srv_q_input:

  Survey catchability array
  `[n_regions × n_yrs × n_srv_fleets × n_sims]`. Default 1.

- t_srv:

  Survey timing as a fraction of the year or season,
  `[n_regions × n_seas × n_srv_fleets]`. Default 1.

- srv_idx_type:

  Index type per fleet: 0/`"abd"` or 1/`"biom"` (default).

- SrvIdx_LikeType:

  Error structure each fleet's index is drawn under: `"lognormal"` (0,
  default), `"normal"` (1) or `"mvn"` (2), matching the estimation
  model. An mvn fleet draws from `SrvIdx_Cov` through a common-factor
  decomposition (see
  [`cov_to_factor`](https://chengmatt.github.io/SPoRC/dev/reference/cov_to_factor.md))
  instead of `ObsSrvIdx_SE`, and its population-specific data source
  stays lognormal.

- SrvIdx_seas_Type, SrvIdx_pop_seas_Type, SrvAgeComps_seas_Type:

  Whether the operating model reports a survey data source once a season
  (`"spltSeas"`, the default) or once a year as a season total
  (`"aggSeas"`), one value for every survey or one per survey. An annual
  total is written into season one with the other seasons left at zero
  and the observation error applied once to that total, so an estimation
  model reading it should mark season one in its `Use` array and set the
  matching argument in
  [`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md).

- SrvIdx_Cov:

  List with one element per fleet holding the fixed covariance over that
  fleet's fitted index observations, ordered by scanning `UseSrvIdx` in
  array order. Required for mvn fleets. Default `NULL`.

- UseSrvIdx:

  Fit flags `[n_regions x n_yrs x n_seas x n_srv_fleets]` from the
  estimation model, used to place each simulated cell in the covariance.
  Its year dim may be shorter than the simulation, in which case later
  years draw with the mean factor scale and loading. Required for mvn
  fleets. Default `NULL`.

- comp_srv_caal_like:

  Conditional age-at-length likelihood per fleet: \`"Multinomial"\` (0),
  \`"Dirichlet-Multinomial"\` (1) or \`"none"\` (999, default). The
  survey twin of \`comp_fish_caal_like\`, and only these two families
  exist for CAAL.

- ISS_Srv_caal:

  Number of fish aged within each length bin, \`n_regions x n_yrs x
  n_seas x n_lens x n_sexes x n_srv_fleets x n_sims\`. A bin whose
  sample size rounds to zero is skipped. \`NULL\` (default) draws no
  CAAL; supplying it alongside a likelihood other than \`"none"\`
  switches \`do_srv_caal\` on. Requires \`n_lens\`.

- ln_Srv_caal_theta:

  Log overdispersion for the Dirichlet-multinomial, \`n_regions x
  n_sexes x n_srv_fleets\`, read under the split types and ignored under
  the multinomial. Default log(1).

- ln_Srv_caal_theta_agg:

  The aggregated type's counterpart, length \`n_srv_fleets\`. Default
  log(1).

- Srv_caal_Type:

  Composition structure per year and fleet, \`n_yrs x n_srv_fleets\`,
  with the codes of \`Fish_caal_Type\`: \`"agg"\` (0), \`"spltRspltS"\`
  (1), \`"spltRjntS"\` (2) or \`"none"\` (999, default). The simulator
  takes the year by fleet array directly.

- comp_srvage_like, comp_srvlen_like, comp_srvage_pop_like,
  comp_srvlen_pop_like:

  Composition likelihood per fleet for the four survey composition data
  sources: 0/`"Multinomial"` (default), 1/`"Dirichlet-Multinomial"`,
  2/`"iid-Logistic-Normal"`, 3/`"1d-Logistic-Normal"` or
  4/`"2d-Logistic-Normal"`.

- ISS_SrvAgeComps, ISS_SrvLenComps:

  Input sample sizes
  `[n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]`.
  Default 100.

- ln_SrvAge_theta, ln_SrvLen_theta:

  Log-scale overdispersion `[n_regions × n_sexes × n_srv_fleets]`, read
  under likelihoods 1-4. Default log(1).

- ln_SrvAge_theta_agg, ln_SrvLen_theta_agg:

  The aggregated types' counterparts, length `n_srv_fleets`. Default
  log(1).

- SrvAge_corr_pars_agg, SrvLen_corr_pars_agg:

  The aggregated types' counterparts, length `n_srv_fleets`, read under
  likelihood 3. Default 0.01.

- SrvAge_corr_pars, SrvLen_corr_pars:

  Correlation parameters `[n_regions × n_sexes × n_srv_fleets × 2]`, the
  age AR1 and the sex correlation, read under likelihoods 3 and 4.
  Default 0.01.

- SrvAgeComps_Type, SrvLenComps_Type, SrvAgeComps_pop_Type,
  SrvLenComps_pop_Type:

  Composition structure `[n_yrs × n_srv_fleets]`: 0/`"agg"`,
  1/`"spltRspltS"`, 2/`"spltRjntS"` (default) or 999/`"none"`.

- ISS_SrvAgeComps_pop, ISS_SrvLenComps_pop:

  The population-specific counterparts, with a leading `n_pop` dim.
  Default 100.

- ln_SrvAge_pop_theta, ln_SrvLen_pop_theta:

  Log-scale overdispersion for the population-specific data sources
  `[n_pop × n_regions × n_sexes × n_srv_fleets]`. Default log(1).

- ln_SrvAge_pop_theta_agg, ln_SrvLen_pop_theta_agg:

  Their aggregated counterparts `[n_pop × n_srv_fleets]`. Default
  log(1).

- SrvAge_pop_corr_pars, SrvLen_pop_corr_pars:

  Correlation parameters for the population-specific data sources
  `[n_pop × n_regions × n_sexes × n_srv_fleets × 2]`. Default 0.01.

- SrvAge_pop_corr_pars_agg, SrvLen_pop_corr_pars_agg:

  Their aggregated counterparts `[n_pop × n_srv_fleets]`. Default 0.01.

## Value

`sim_list` with the survey fields appended: `$srv_sel`, `$srv_q`,
`$ObsSrvIdx_SE`, `$ObsSrvIdx_pop_SE`, `$t_srv`, `$srv_idx_type`, and,
for each of the four composition data sources, its likelihood, input
sample sizes, overdispersion, correlation parameters and composition
type. Character codes are converted to integers before storage.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
