# Setup Simulation Fishing Inputs

Sets and validates the fishing inputs of a \`sim_list\`: fishing
mortality, selectivity, catchability, observation error, and the age and
length composition settings for the aggregate and population-specific
data sources.

## Usage

``` r
Setup_Sim_Fishing(
  sim_list,
  ln_sigmaC = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_fish_fleets)),
  ln_sigmaC_pop = array(log(0.02), dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
  ln_sigmaCAA = array(log(0.2), dim = c(sim_list$n_obs_ages, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_sigmaDAA = array(log(0.2), dim = c(sim_list$n_obs_ages, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  UseCatchAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
  UseDiscardAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
  ObsCatchAA_SE = NULL,
  ObsDiscardAA_SE = NULL,
  CatchAA_Type = "spltRaggS",
  DiscardAA_Type = "spltRaggS",
  CatchAA_LikeType = "lognormal",
  DiscardAA_LikeType = "lognormal",
  CatchAA_sigma_form = "none",
  DiscardAA_sigma_form = "none",
  use_catch_aa = rep(0, sim_list$n_fish_fleets),
  use_discard_aa = rep(0, sim_list$n_fish_fleets),
  catch_units = array(1, dim = c(sim_list$n_fish_fleets)),
  init_F_val = array(0, dim = c(sim_list$n_regions, sim_list$n_seas,
    sim_list$n_fish_fleets)),
  Fmort_input = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_fish_fleets, sim_list$n_sims)),
  fish_sel_input,
  fish_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_fish_fleets, sim_list$n_sims)),
  ObsFishIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_fish_fleets)),
  ObsFishIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
  fish_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
  FishIdx_LikeType = rep(0, sim_list$n_fish_fleets),
  Catch_seas_Type = NULL,
  Catch_pop_seas_Type = NULL,
  FishIdx_seas_Type = NULL,
  FishIdx_pop_seas_Type = NULL,
  FishAgeComps_seas_Type = NULL,
  FishIdx_Cov = NULL,
  UseFishIdx = NULL,
  t_fish = array(0, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)),
  comp_fish_caal_like = rep(999, sim_list$n_fish_fleets),
  ISS_Fish_caal = NULL,
  ln_Fish_caal_theta = NULL,
  ln_Fish_caal_theta_agg = NULL,
  Fish_caal_Type = array(999, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  comp_fishage_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
  ln_FishAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_FishAge_theta_agg = rep(log(1), sim_list$n_fish_fleets),
  FishAge_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
  FishAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets, 2)),
  FishAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  comp_fishlen_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
  ln_FishLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_FishLen_theta_agg = rep(log(1), sim_list$n_fish_fleets),
  FishLen_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
  FishLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets, 2)),
  FishLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  comp_fishage_pop_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishAgeComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets,
    sim_list$n_sims)),
  ln_FishAge_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_fish_fleets)),
  ln_FishAge_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishAge_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
  FishAge_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishAgeComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  comp_fishlen_pop_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishLenComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets,
    sim_list$n_sims)),
  ln_FishLen_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_fish_fleets)),
  ln_FishLen_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishLen_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
  FishLen_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishLenComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  ret_sel_input = array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets,
    sim_list$n_sims)),
  dmr_input = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_fish_fleets, sim_list$n_sims)),
  discard_units = array(3, dim = c(sim_list$n_fish_fleets)),
  ln_sigmaD = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_fish_fleets)),
  ln_sigmaD_pop = array(log(0.02), dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
  comp_fishage_discard_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishAgeComps_discard = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
  ln_FishAge_discard_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_FishAge_discard_theta_agg = rep(log(1), sim_list$n_fish_fleets),
  FishAge_discard_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets, 2)),
  FishAge_discard_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
  FishAgeComps_discard_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  comp_fishlen_discard_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishLenComps_discard = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
  ln_FishLen_discard_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_FishLen_discard_theta_agg = rep(log(1), sim_list$n_fish_fleets),
  FishLen_discard_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets, 2)),
  FishLen_discard_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
  FishLenComps_discard_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  comp_fishage_discard_pop_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishAgeComps_discard_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets,
    sim_list$n_sims)),
  ln_FishAge_discard_pop_theta = array(log(1), dim = c(sim_list$n_pop,
    sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
  ln_FishAge_discard_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishAge_discard_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
  FishAge_discard_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishAgeComps_discard_pop_Type = array(2, dim = c(sim_list$n_yrs,
    sim_list$n_fish_fleets)),
  comp_fishlen_discard_pop_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishLenComps_discard_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets,
    sim_list$n_sims)),
  ln_FishLen_discard_pop_theta = array(log(1), dim = c(sim_list$n_pop,
    sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
  ln_FishLen_discard_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishLen_discard_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
  FishLen_discard_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop,
    sim_list$n_fish_fleets)),
  FishLenComps_discard_pop_Type = array(2, dim = c(sim_list$n_yrs,
    sim_list$n_fish_fleets))
)
```

## Arguments

- sim_list:

  Simulation list holding \`n_pop\`, \`n_regions\`, \`n_yrs\`,
  \`n_seas\`, \`n_ages\`, \`n_sexes\`, \`n_fish_fleets\` and \`n_sims\`.

- ln_sigmaC, ln_sigmaC_pop:

  Log-scale observation sd for total and population-specific catch,
  \`n_regions x n_yrs x n_seas x n_fish_fleets\` with a leading
  \`n_pop\` for the second. Default log(0.02).

- ln_sigmaCAA, ln_sigmaDAA:

  Log-scale observation error for the at-age data sources, \`n_obs_ages
  x n_sexes x n_fish_fleets\`. The sex dim is required.

- UseCatchAA, UseDiscardAA:

  Integer arrays \`n_regions x n_yrs x n_seas x n_obs_ages x n_sexes x
  n_fish_fleets\`, \`1\` where an at-age observation is drawn. The draws
  sit on the observed ages from \`Setup_Sim_Dim\`, read through
  \`AgeingError_fish_input\` the way the estimation model reads them.
  The sex dim is required: a data source summed over sexes has its flag
  in sex slot one.

- ObsCatchAA_SE, ObsDiscardAA_SE:

  Reported standard errors shaped like the use arrays, read only when
  the data source's \`sigma_form\` asks for them.

- CatchAA_Type, DiscardAA_Type:

  Which dims each fleet reports separately: \`"agg"\`, \`"spltRaggS"\`
  (default), \`"aggRspltS"\` or \`"spltRspltS"\`. A summed dim is drawn
  once, into slot one.

- CatchAA_LikeType, DiscardAA_LikeType:

  \`"lognormal"\` (default) or \`"normal"\`, per fleet.

- CatchAA_sigma_form, DiscardAA_sigma_form:

  Where the observation error comes from: \`"none"\` (default),
  \`"data"\`, \`"est_additive"\` or \`"est_quadrature"\`.

- use_catch_aa, use_discard_aa:

  Integer vectors \`n_fish_fleets\`, \`1\` for fleets whose at-age data
  sources are drawn.

- catch_units:

  Catch units per fleet, 0 = abundance, 1 = biomass (default).

- init_F_val:

  Initial fishing mortality, \`n_regions x n_seas x n_fish_fleets\`.
  Default 0.

- Fmort_input:

  Fishing mortality, \`n_regions x n_yrs x n_seas x n_fish_fleets x
  n_sims\`. Default 0.1.

- fish_sel_input:

  Fishery selectivity, \`n_pop x n_regions x n_yrs x n_seas x n_ages x
  n_sexes x n_fish_fleets x n_sims\`.

- fish_q_input:

  Catchability, \`n_regions x n_yrs x n_fish_fleets x n_sims\`. Default
  1.

- ObsFishIdx_SE, ObsFishIdx_pop_SE:

  Observation sd for the fishery indices, \`n_regions x n_yrs x n_seas x
  n_fish_fleets\` with a leading \`n_pop\` for the second. Default 0.2.

- fish_idx_type:

  Index type, 0 = abundance, 1 = biomass (default), \`n_regions x
  n_fish_fleets\`.

- FishIdx_LikeType:

  Error structure each fleet's index is drawn under: \`"lognormal"\` (0,
  default), \`"normal"\` (1) or \`"mvn"\` (2), matching the estimation
  model. An mvn fleet draws from \`FishIdx_Cov\` through a common-factor
  decomposition (see
  [`cov_to_factor`](https://chengmatt.github.io/SPoRC/dev/reference/cov_to_factor.md))
  instead of \`ObsFishIdx_SE\`, and its population-specific data source
  stays lognormal.

- Catch_seas_Type, Catch_pop_seas_Type, FishIdx_seas_Type,
  FishIdx_pop_seas_Type, FishAgeComps_seas_Type:

  Whether the operating model reports a data source once a season
  (\`"spltSeas"\`, the default) or once a year as a season total
  (\`"aggSeas"\`), one value for every fleet or one per fleet. An annual
  total is written into season one with the other seasons left at zero
  and the observation error applied once to that total, so an estimation
  model reading it should mark season one in its \`Use\` array and set
  the matching argument in
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md)
  or
  [`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md).

- FishIdx_Cov:

  List with one element per fleet holding the fixed covariance over that
  fleet's fitted index observations, ordered by scanning \`UseFishIdx\`
  in array order. Required for mvn fleets. Default \`NULL\`.

- UseFishIdx:

  Fit flags \`n_regions x n_yrs x n_seas x n_fish_fleets\` from the
  estimation model, used to place each simulated cell in the covariance.
  Its year dim may be shorter than the simulation, in which case later
  years draw with the mean factor scale and loading. Required for mvn
  fleets. Default \`NULL\`.

- t_fish:

  Fishery index timing, \`n_regions x n_seas x n_fish_fleets\`, the
  fraction of the season elapsed when the index is observed. Numbers at
  age are decayed by \`exp(-t_fish \* ZAA)\` before the index is formed,
  matching \`t_srv\` and the estimation model's own \`t_fish\`. Default
  0.

- comp_fish_caal_like:

  Conditional age-at-length likelihood per fleet: \`"Multinomial"\` (0),
  \`"Dirichlet-Multinomial"\` (1) or \`"none"\` (999, default). Only
  these two families exist for CAAL, since a CAAL row is the age
  composition of the otoliths from one length bin, usually a small and
  mostly zero sample.

- ISS_Fish_caal:

  Number of fish aged within each length bin, \`n_regions x n_yrs x
  n_seas x n_lens x n_sexes x n_fish_fleets x n_sims\`. A bin whose
  sample size rounds to zero is skipped. \`NULL\` (default) draws no
  CAAL; supplying it alongside a likelihood other than \`"none"\` is
  what switches \`do_fish_caal\` on. Requires \`n_lens\`.

- ln_Fish_caal_theta:

  Log overdispersion for the Dirichlet-multinomial, \`n_regions x
  n_sexes x n_fish_fleets\`. Read under the split types, \`\[r, s, f\]\`
  when sexes are split and \`\[r, 1, f\]\` when they are joint, and
  ignored under the multinomial. Default log(1).

- ln_Fish_caal_theta_agg:

  The aggregated type's counterpart, length \`n_fish_fleets\`. Default
  log(1).

- Fish_caal_Type:

  Composition structure per year and fleet, \`n_yrs x n_fish_fleets\`:
  \`"agg"\` (0) pools regions and sexes and is drawn once when the
  region loop reaches the last region, \`"spltRspltS"\` (1) draws each
  sex in a bin as its own sample, \`"spltRjntS"\` (2) draws one sample
  across the age by sex stack, and \`"none"\` (999, default) skips the
  fleet that year. The simulator takes the year by fleet array directly
  rather than the estimation model's \`"CompType_Year_x-y_Fleet_z"\`
  strings.

- comp_fishage_like, comp_fishlen_like, comp_fishage_pop_like,
  comp_fishlen_pop_like, comp_fishage_discard_like,
  comp_fishlen_discard_like, comp_fishage_discard_pop_like,
  comp_fishlen_discard_pop_like:

  Composition likelihood per fleet for the eight fishery composition
  data sources: 0 = multinomial (default), 1 = Dirichlet-multinomial,
  2-4 = logistic-normal, 999 = none.

- ISS_FishAgeComps, ISS_FishLenComps, ISS_FishAgeComps_discard,
  ISS_FishLenComps_discard:

  Input sample sizes, \`n_regions x n_yrs x n_seas x n_sexes x
  n_fish_fleets x n_sims\`. Default 100.

- ln_FishAge_theta, ln_FishLen_theta, ln_FishAge_discard_theta,
  ln_FishLen_discard_theta:

  Log-scale overdispersion, \`n_regions x n_sexes x n_fish_fleets\`.
  Default log(1).

- ln_FishAge_theta_agg, ln_FishLen_theta_agg,
  ln_FishAge_discard_theta_agg, ln_FishLen_discard_theta_agg:

  The aggregated types' counterparts, length \`n_fish_fleets\`. Default
  log(1).

- FishAge_corr_pars_agg, FishLen_corr_pars_agg,
  FishAge_discard_corr_pars_agg, FishLen_discard_corr_pars_agg:

  Their aggregated counterparts, length \`n_fish_fleets\`. Default 0.01.

- FishAge_corr_pars, FishLen_corr_pars, FishAge_discard_corr_pars,
  FishLen_discard_corr_pars:

  Correlation parameters, \`n_regions x n_sexes x n_fish_fleets x 2\`.
  Default 0.01.

- FishAgeComps_Type, FishLenComps_Type, FishAgeComps_pop_Type,
  FishLenComps_pop_Type, FishAgeComps_discard_Type,
  FishLenComps_discard_Type, FishAgeComps_discard_pop_Type,
  FishLenComps_discard_pop_Type:

  Composition structure per year and fleet, \`n_yrs x n_fish_fleets\`: 0
  = aggregated, 1 = split region and sex, 2 = split region joint sex
  (default), 999 = none.

- ISS_FishAgeComps_pop, ISS_FishLenComps_pop,
  ISS_FishAgeComps_discard_pop, ISS_FishLenComps_discard_pop:

  The population-specific counterparts, with a leading \`n_pop\` dim.
  Default 100.

- ln_FishAge_pop_theta, ln_FishLen_pop_theta,
  ln_FishAge_discard_pop_theta, ln_FishLen_discard_pop_theta:

  Log-scale overdispersion for the population-specific data sources,
  \`n_pop x n_regions x n_sexes x n_fish_fleets\`. Default log(1).

- ln_FishAge_pop_theta_agg, ln_FishLen_pop_theta_agg,
  ln_FishAge_discard_pop_theta_agg, ln_FishLen_discard_pop_theta_agg:

  Their aggregated counterparts, \`n_pop x n_fish_fleets\`. Default
  log(1).

- FishAge_pop_corr_pars, FishLen_pop_corr_pars,
  FishAge_discard_pop_corr_pars, FishLen_discard_pop_corr_pars:

  Correlation parameters for the population-specific data sources,
  \`n_pop x n_regions x n_sexes x n_fish_fleets x 2\`. Default 0.01.

- FishAge_pop_corr_pars_agg, FishLen_pop_corr_pars_agg,
  FishAge_discard_pop_corr_pars_agg, FishLen_discard_pop_corr_pars_agg:

  Their aggregated counterparts, \`n_pop x n_fish_fleets\`. Default
  0.01.

- ret_sel_input:

  Retained selectivity at age, \`n_pop x n_regions x n_yrs x n_seas x
  n_ages x n_sexes x n_fish_fleets x n_sims\`. Default 1.

- dmr_input:

  Discard mortality rate, \`n_regions x n_yrs x n_seas x n_fish_fleets x
  n_sims\`. Default 0.

- discard_units:

  Discard units per fleet: 0 = abundance, 1 = biomass, 2 = abundance
  fraction, 3 = biomass fraction (default).

- ln_sigmaD, ln_sigmaD_pop:

  Log-scale observation sd for discards, \`n_regions x n_yrs x n_seas x
  n_fish_fleets\` with a leading \`n_pop\` for the second. Default
  log(0.02).

## Value

A modified \`sim_list\` with validated fishing inputs.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
