# Setup Simulation Fishing Inputs

Initializes and validates fishing-related inputs for a simulation list
(\`sim_list\`). This includes fishing mortality, selectivity,
catchability, observation error, and age- and length-composition
parameters for both aggregate and population-specific data.

## Usage

``` r
Setup_Sim_Fishing(
  sim_list,
  ln_sigmaC = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_seas, sim_list$n_fish_fleets)),
  ln_sigmaC_pop = array(log(0.02), dim = c(sim_list$n_pop, sim_list$n_regions,
    sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
  ln_sigmaCAA = array(log(0.2), dim = c(sim_list$n_ages, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_sigmaDAA = array(log(0.2), dim = c(sim_list$n_ages, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  UseCatchAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
  UseDiscardAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
    sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
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

  A list containing simulation settings, including the number of
  populations (\`n_pop\`), regions (\`n_regions\`), years (\`n_yrs\`),
  seasons (\`n_seas\`), ages (\`n_ages\`), sexes (\`n_sexes\`), fishing
  fleets (\`n_fish_fleets\`), and simulations (\`n_sims\`).

- ln_sigmaC:

  Numeric array. Log-scale observation SD for total catch, dimensions
  \`n_regions x n_yrs x n_seas x n_fish_fleets\`. Default: log(0.02).

- ln_sigmaC_pop:

  Numeric array. Log-scale observation SD for population-specific catch,
  dimensions \`n_pop x n_regions x n_yrs x n_seas x n_fish_fleets\`.
  Default: log(0.02).

- ln_sigmaCAA, ln_sigmaDAA:

  Log-scale observation error for the at-age streams, \`n_ages x n_sexes
  x n_fish_fleets\`. An array without the sex margin is required.

- UseCatchAA, UseDiscardAA:

  Integer arrays \`n_regions x n_yrs x n_seas x n_ages x n_sexes x
  n_fish_fleets\`, \`1\` where an at-age observation is drawn. The sex
  margin is required: a stream summed over sexes carries its flag in sex
  slot one.

- ObsCatchAA_SE, ObsDiscardAA_SE:

  Reported standard errors shaped like the use arrays, read only when
  the stream's \`sigma_form\` asks for them.

- CatchAA_Type, DiscardAA_Type:

  Which margins each fleet reports separately: \`"agg"\`,
  \`"spltRaggS"\` (default), \`"aggRspltS"\` or \`"spltRspltS"\`. A
  summed margin is drawn once, into slot one.

- CatchAA_LikeType, DiscardAA_LikeType:

  \`"lognormal"\` (default) or \`"normal"\`, per fleet.

- CatchAA_sigma_form, DiscardAA_sigma_form:

  Where the observation error comes from: \`"none"\` (default),
  \`"data"\`, \`"est_additive"\` or \`"est_quadrature"\`.

- use_catch_aa, use_discard_aa:

  Integer vectors \`n_fish_fleets\`, \`1\` for fleets whose at-age
  streams are drawn.

- catch_units:

  Numeric vector. Catch units (0 = abundance, 1 = biomass), length
  \`n_fish_fleets\`. Default: 1.

- init_F_val:

  Numeric array. Initial fishing mortality, dimensions \`n_regions x
  n_seas x n_fish_fleets\`. Default: 0.

- Fmort_input:

  Numeric array. Fishing mortality, dimensions \`n_regions x n_yrs x
  n_seas x n_fish_fleets x n_sims\`. Default: 0.1.

- fish_sel_input:

  Numeric array. Fishery selectivity, dimensions \`n_pop x n_regions x
  n_yrs x n_seas x n_ages x n_sexes x n_fish_fleets x n_sims\`.

- fish_q_input:

  Numeric array. Catchability, dimensions \`n_regions x n_yrs x
  n_fish_fleets x n_sims\`. Default: 1.

- ObsFishIdx_SE:

  Numeric array. Observation SD for fishery indices, dimensions
  \`n_regions x n_yrs x n_seas x n_fish_fleets\`. Default: 0.2.

- ObsFishIdx_pop_SE:

  Numeric array. Observation SD for population-specific fishery indices,
  dimensions \`n_pop x n_regions x n_yrs x n_seas x n_fish_fleets\`.
  Default: 0.2.

- fish_idx_type:

  Numeric array. Index type (0 = abundance, 1 = biomass), dimensions
  \`n_regions x n_fish_fleets\`. Default: 1.

- FishIdx_LikeType:

  Character or numeric vector, length \`n_fish_fleets\`. Error structure
  each fleet's index is drawn under: `"lognormal"` (0), `"normal"` (1),
  or `"mvn"` (2), matching the estimation model's `FishIdx_LikeType`. An
  mvn fleet draws from `FishIdx_Cov` through a common-factor
  decomposition (see
  [`cov_to_factor`](https://chengmatt.github.io/SPoRC/dev/reference/cov_to_factor.md))
  instead of `ObsFishIdx_SE`, and its population-specific stream stays
  lognormal. Default: lognormal for every fleet.

- FishIdx_Cov:

  List with one element per fishery fleet holding the fixed covariance
  over that fleet's fitted index observations, ordered by scanning
  `UseFishIdx` in array order (region fastest, then year, then season).
  Required for mvn fleets. Default: `NULL`.

- UseFishIdx:

  Numeric array `[n_regions x n_yrs x n_seas x n_fish_fleets]` of fit
  flags from the estimation model, used to position each simulated cell
  in the covariance. Its year dimension may be shorter than the
  simulation, in which case later years draw with the mean factor scale
  and loading. Required for mvn fleets. Default: `NULL`.

- t_fish:

  Numeric array `[n_regions x n_seas x n_fish_fleets]` giving the
  fishery index timing, the fraction of the season elapsed when the
  index is observed. Numbers at age are decayed by `exp(-t_fish * ZAA)`
  before the index is formed, matching `t_srv` for surveys and the
  estimation model's own `t_fish`. Defaults to `0` (start of season).

- comp_fish_caal_like:

  Character or numeric vector \`n_fish_fleets\` giving the conditional
  age-at-length likelihood per fleet: \`"Multinomial"\` (0),
  \`"Dirichlet-Multinomial"\` (1), or \`"none"\` (999). Only these two
  families exist for CAAL: a CAAL row is the age composition of the
  otoliths taken from one length bin, usually a small and mostly zero
  sample, which the logistic-normal forms cannot support. Default:
  \`"none"\` for every fleet.

- ISS_Fish_caal:

  Numeric array. Number of fish aged within each length bin, dimensions
  \`n_regions x n_yrs x n_seas x n_lens x n_sexes x n_fish_fleets x
  n_sims\`. A bin whose sample size rounds to zero is skipped. \`NULL\`
  (the default) draws no CAAL; supplying it alongside a likelihood other
  than \`"none"\` is what switches \`do_fish_caal\` on. Requires
  \`n_lens\`.

- ln_Fish_caal_theta:

  Numeric array. Log overdispersion for the Dirichlet-multinomial,
  dimensions \`n_regions x n_sexes x n_fish_fleets\`. Read under the
  split types, \`\[r, s, f\]\` when sexes are split and \`\[r, 1, f\]\`
  when they are joint, and ignored under the multinomial. Default:
  log(1).

- ln_Fish_caal_theta_agg:

  Numeric vector \`n_fish_fleets\`. The aggregated type's counterpart to
  \`ln_Fish_caal_theta\`. Default: log(1).

- Fish_caal_Type:

  Numeric or character array giving the composition structure per year
  and fleet, dimensions \`n_yrs x n_fish_fleets\`: \`"agg"\` (0) pools
  regions and sexes and is drawn once when the region loop reaches the
  last region, \`"spltRspltS"\` (1) draws each sex in a bin as its own
  sample, \`"spltRjntS"\` (2) draws one sample across the age by sex
  stack, and \`"none"\` (999) skips the fleet in that year. Unlike the
  estimation model, which parses \`"CompType_Year_x-y_Fleet_z"\`
  strings, the simulator takes the year by fleet array directly.
  Default: \`"none"\` throughout.

- comp_fishage_like:

  Numeric vector. Likelihood for age composition (0 = Multinomial, 1 =
  Dirichlet-Multinomial, 2-4 = Logistic-Normal variants), length
  \`n_fish_fleets\`. Default: 0.

- ISS_FishAgeComps:

  Numeric array. Effective sample sizes for age compositions, dimensions
  \`n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims\`.
  Default: 100.

- ln_FishAge_theta:

  Numeric array. Log-scale overdispersion for fishery age compositions,
  dimensions \`n_regions x n_sexes x n_fish_fleets\`. Default: log(1).

- ln_FishAge_theta_agg:

  Numeric vector. Aggregated log-scale overdispersion for fishery age
  compositions, length \`n_fish_fleets\`. Default: log(1).

- FishAge_corr_pars_agg:

  Numeric vector. Aggregated correlation parameters for fishery age
  compositions, length \`n_fish_fleets\`. Default: 0.01.

- FishAge_corr_pars:

  Numeric array. Correlation parameters for fishery age compositions,
  dimensions \`n_regions x n_sexes x n_fish_fleets x 2\`. Default: 0.01.

- FishAgeComps_Type:

  Numeric array. Composition structure for fishery age compositions (0 =
  aggregated, 1 = split region/sex, 2 = split region joint sex, 999 =
  none), dimensions \`n_yrs x n_fish_fleets\`. Default: 2.

- comp_fishlen_like:

  Numeric vector. Likelihood for length composition (0 = Multinomial, 1
  = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants), length
  \`n_fish_fleets\`. Default: 0.

- ISS_FishLenComps:

  Numeric array. Effective sample sizes for length compositions,
  dimensions \`n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x
  n_sims\`. Default: 100.

- ln_FishLen_theta:

  Numeric array. Log-scale overdispersion for fishery length
  compositions, dimensions \`n_regions x n_sexes x n_fish_fleets\`.
  Default: log(1).

- ln_FishLen_theta_agg:

  Numeric vector. Aggregated log-scale overdispersion for fishery length
  compositions, length \`n_fish_fleets\`. Default: log(1).

- FishLen_corr_pars_agg:

  Numeric vector. Aggregated correlation parameters for fishery length
  compositions, length \`n_fish_fleets\`. Default: 0.01.

- FishLen_corr_pars:

  Numeric array. Correlation parameters for fishery length compositions,
  dimensions \`n_regions x n_sexes x n_fish_fleets x 2\`. Default: 0.01.

- FishLenComps_Type:

  Numeric array. Composition structure for fishery length compositions
  (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999
  = none), dimensions \`n_yrs x n_fish_fleets\`. Default: 2.

- comp_fishage_pop_like:

  Numeric vector. Likelihood for population-specific fishery age
  composition (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 =
  Logistic-Normal variants), length \`n_fish_fleets\`. Default: 0.

- ISS_FishAgeComps_pop:

  Numeric array. Effective sample sizes for population-specific fishery
  age compositions, dimensions \`n_pop x n_regions x n_yrs x n_seas x
  n_sexes x n_fish_fleets x n_sims\`. Default: 100.

- ln_FishAge_pop_theta:

  Numeric array. Log-scale overdispersion for population-specific
  fishery age compositions, dimensions \`n_pop x n_regions x n_sexes x
  n_fish_fleets\`. Default: log(1).

- ln_FishAge_pop_theta_agg:

  Numeric array. Aggregated log-scale overdispersion for
  population-specific fishery age compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: log(1).

- FishAge_pop_corr_pars:

  Numeric array. Correlation parameters for population-specific fishery
  age compositions, dimensions \`n_pop x n_regions x n_sexes x
  n_fish_fleets x 2\`. Default: 0.01.

- FishAge_pop_corr_pars_agg:

  Numeric array. Aggregated correlation parameters for
  population-specific fishery age compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: 0.01.

- FishAgeComps_pop_Type:

  Numeric array. Composition structure for population-specific fishery
  age compositions (0 = aggregated, 1 = split region/sex, 2 = split
  region joint sex, 999 = none), dimensions \`n_yrs x n_fish_fleets\`.
  Default: 2.

- comp_fishlen_pop_like:

  Numeric vector. Likelihood for population-specific fishery length
  composition (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 =
  Logistic-Normal variants), length \`n_fish_fleets\`. Default: 0.

- ISS_FishLenComps_pop:

  Numeric array. Effective sample sizes for population-specific fishery
  length compositions, dimensions \`n_pop x n_regions x n_yrs x n_seas x
  n_sexes x n_fish_fleets x n_sims\`. Default: 100.

- ln_FishLen_pop_theta:

  Numeric array. Log-scale overdispersion for population-specific
  fishery length compositions, dimensions \`n_pop x n_regions x n_sexes
  x n_fish_fleets\`. Default: log(1).

- ln_FishLen_pop_theta_agg:

  Numeric array. Aggregated log-scale overdispersion for
  population-specific fishery length compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: log(1).

- FishLen_pop_corr_pars:

  Numeric array. Correlation parameters for population-specific fishery
  length compositions, dimensions \`n_pop x n_regions x n_sexes x
  n_fish_fleets x 2\`. Default: 0.01.

- FishLen_pop_corr_pars_agg:

  Numeric array. Aggregated correlation parameters for
  population-specific fishery length compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: 0.01.

- FishLenComps_pop_Type:

  Numeric array. Composition structure for population-specific fishery
  length compositions (0 = aggregated, 1 = split region/sex, 2 = split
  region joint sex, 999 = none), dimensions \`n_yrs x n_fish_fleets\`.
  Default: 2.

- ret_sel_input:

  Numeric array. Retained selectivity at age, dimensions \`n_pop x
  n_regions x n_yrs x n_seas x n_ages x n_sexes x n_fish_fleets x
  n_sims\`. Default: 1.

- dmr_input:

  Numeric array. Discard mortality rate, dimensions \`n_regions x n_yrs
  x n_seas x n_fish_fleets x n_sims\`. Default: 0.

- discard_units:

  Numeric vector. Discard units (0 = abundance, 1 = biomass, 2 =
  abundance fraction, 3 = biomass fraction), length \`n_fish_fleets\`.
  Default: 3.

- ln_sigmaD:

  Numeric array. Log-scale observation SD for discards, dimensions
  \`n_regions x n_yrs x n_seas x n_fish_fleets\`. Default: log(0.02).

- ln_sigmaD_pop:

  Numeric array. Log-scale observation SD for population-specific
  discards, dimensions \`n_pop x n_regions x n_yrs x n_seas x
  n_fish_fleets\`. Default: log(0.02).

- comp_fishage_discard_like:

  Numeric vector. Likelihood for discard age composition (0 =
  Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal
  variants, 999 = none), length \`n_fish_fleets\`. Default: 0.

- ISS_FishAgeComps_discard:

  Numeric array. Effective sample sizes for discard age compositions,
  dimensions \`n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x
  n_sims\`. Default: 100.

- ln_FishAge_discard_theta:

  Numeric array. Log-scale overdispersion for discard age compositions,
  dimensions \`n_regions x n_sexes x n_fish_fleets\`. Default: log(1).

- ln_FishAge_discard_theta_agg:

  Numeric vector. Aggregated log-scale overdispersion for discard age
  compositions, length \`n_fish_fleets\`. Default: log(1).

- FishAge_discard_corr_pars:

  Numeric array. Correlation parameters for discard age compositions,
  dimensions \`n_regions x n_sexes x n_fish_fleets x 2\`. Default: 0.01.

- FishAge_discard_corr_pars_agg:

  Numeric vector. Aggregated correlation parameters for discard age
  compositions, length \`n_fish_fleets\`. Default: 0.01.

- FishAgeComps_discard_Type:

  Numeric array. Composition structure for discard age compositions (0 =
  aggregated, 1 = split region/sex, 2 = split region joint sex, 999 =
  none), dimensions \`n_yrs x n_fish_fleets\`. Default: 2.

- comp_fishlen_discard_like:

  Numeric vector. Likelihood for discard length composition (0 =
  Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal
  variants, 999 = none), length \`n_fish_fleets\`. Default: 0.

- ISS_FishLenComps_discard:

  Numeric array. Effective sample sizes for discard length compositions,
  dimensions \`n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x
  n_sims\`. Default: 100.

- ln_FishLen_discard_theta:

  Numeric array. Log-scale overdispersion for discard length
  compositions, dimensions \`n_regions x n_sexes x n_fish_fleets\`.
  Default: log(1).

- ln_FishLen_discard_theta_agg:

  Numeric vector. Aggregated log-scale overdispersion for discard length
  compositions, length \`n_fish_fleets\`. Default: log(1).

- FishLen_discard_corr_pars:

  Numeric array. Correlation parameters for discard length compositions,
  dimensions \`n_regions x n_sexes x n_fish_fleets x 2\`. Default: 0.01.

- FishLen_discard_corr_pars_agg:

  Numeric vector. Aggregated correlation parameters for discard length
  compositions, length \`n_fish_fleets\`. Default: 0.01.

- FishLenComps_discard_Type:

  Numeric array. Composition structure for discard length compositions
  (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999
  = none), dimensions \`n_yrs x n_fish_fleets\`. Default: 2.

- comp_fishage_discard_pop_like:

  Numeric vector. Likelihood for population-specific discard age
  composition (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 =
  Logistic-Normal variants, 999 = none), length \`n_fish_fleets\`.
  Default: 0.

- ISS_FishAgeComps_discard_pop:

  Numeric array. Effective sample sizes for population-specific discard
  age compositions, dimensions \`n_pop x n_regions x n_yrs x n_seas x
  n_sexes x n_fish_fleets x n_sims\`. Default: 100.

- ln_FishAge_discard_pop_theta:

  Numeric array. Log-scale overdispersion for population-specific
  discard age compositions, dimensions \`n_pop x n_regions x n_sexes x
  n_fish_fleets\`. Default: log(1).

- ln_FishAge_discard_pop_theta_agg:

  Numeric array. Aggregated log-scale overdispersion for
  population-specific discard age compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: log(1).

- FishAge_discard_pop_corr_pars:

  Numeric array. Correlation parameters for population-specific discard
  age compositions, dimensions \`n_pop x n_regions x n_sexes x
  n_fish_fleets x 2\`. Default: 0.01.

- FishAge_discard_pop_corr_pars_agg:

  Numeric array. Aggregated correlation parameters for
  population-specific discard age compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: 0.01.

- FishAgeComps_discard_pop_Type:

  Numeric array. Composition structure for population-specific discard
  age compositions (0 = aggregated, 1 = split region/sex, 2 = split
  region joint sex, 999 = none), dimensions \`n_yrs x n_fish_fleets\`.
  Default: 2.

- comp_fishlen_discard_pop_like:

  Numeric vector. Likelihood for population-specific discard length
  composition (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 =
  Logistic-Normal variants, 999 = none), length \`n_fish_fleets\`.
  Default: 0.

- ISS_FishLenComps_discard_pop:

  Numeric array. Effective sample sizes for population-specific discard
  length compositions, dimensions \`n_pop x n_regions x n_yrs x n_seas x
  n_sexes x n_fish_fleets x n_sims\`. Default: 100.

- ln_FishLen_discard_pop_theta:

  Numeric array. Log-scale overdispersion for population-specific
  discard length compositions, dimensions \`n_pop x n_regions x n_sexes
  x n_fish_fleets\`. Default: log(1).

- ln_FishLen_discard_pop_theta_agg:

  Numeric array. Aggregated log-scale overdispersion for
  population-specific discard length compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: log(1).

- FishLen_discard_pop_corr_pars:

  Numeric array. Correlation parameters for population-specific discard
  length compositions, dimensions \`n_pop x n_regions x n_sexes x
  n_fish_fleets x 2\`. Default: 0.01.

- FishLen_discard_pop_corr_pars_agg:

  Numeric array. Aggregated correlation parameters for
  population-specific discard length compositions, dimensions \`n_pop x
  n_fish_fleets\`. Default: 0.01.

- FishLenComps_discard_pop_Type:

  Numeric array. Composition structure for population-specific discard
  length compositions (0 = aggregated, 1 = split region/sex, 2 = split
  region joint sex, 999 = none), dimensions \`n_yrs x n_fish_fleets\`.
  Default: 2.

## Value

A modified \`sim_list\` with validated fishing-related inputs.

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
