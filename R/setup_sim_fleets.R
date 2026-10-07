# Stage 1 of 3: model setup
#
# Operating model fleet inputs. Setup_Sim_Fishing and Setup_Sim_Survey each populate a whole fleet at once,
# spanning every topic the setup_fishery_*.R and setup_survey_*.R files split apart on the estimation side.

#' Setup Simulation Fishing Inputs
#'
#' Sets and validates the fishing inputs of a `sim_list`: fishing mortality,
#' selectivity, catchability, observation error, and the age and length
#' composition settings for the aggregate and population-specific data sources.
#'
#' @param sim_list Simulation list holding `n_pop`, `n_regions`, `n_yrs`, `n_seas`,
#'   `n_ages`, `n_sexes`, `n_fish_fleets` and `n_sims`.
#' @param ln_sigmaC,ln_sigmaC_pop Log-scale observation sd for total and
#'   population-specific catch, `n_regions x n_yrs x n_seas x n_fish_fleets` with a
#'   leading `n_pop` for the second. Default log(0.02).
#' @param catch_units Catch units per fleet, 0 = abundance, 1 = biomass (default).
#' @param init_F_val Initial fishing mortality, `n_regions x n_seas x
#'   n_fish_fleets`. Default 0.
#' @param Fmort_input Fishing mortality, `n_regions x n_yrs x n_seas x
#'   n_fish_fleets x n_sims`. Default 0.1.
#' @param fish_sel_input Fishery selectivity, `n_pop x n_regions x n_yrs x n_seas x
#'   n_ages x n_sexes x n_fish_fleets x n_sims`.
#' @param fish_q_input Catchability, `n_regions x n_yrs x n_fish_fleets x n_sims`.
#'   Default 1.
#' @param ObsFishIdx_SE,ObsFishIdx_pop_SE Observation sd for the fishery indices,
#'   `n_regions x n_yrs x n_seas x n_fish_fleets` with a leading `n_pop` for the
#'   second. Default 0.2.
#' @param fish_idx_type Index type, 0 = abundance, 1 = biomass (default),
#'   `n_regions x n_fish_fleets`.
#' @param fish_idx_ages Ages counted in each fleet's index total, a 0/1 array
#'   `n_ages x n_fish_fleets`, the estimation model's `fish_idx_ages`. `NULL`
#'   (default) counts every age.
#' @param Catch_seas_Type,Catch_pop_seas_Type,Discard_seas_Type,Discard_pop_seas_Type,FishIdx_seas_Type,FishIdx_pop_seas_Type,FishAgeComps_seas_Type,FishLenComps_seas_Type,FishAgeComps_pop_seas_Type,FishLenComps_pop_seas_Type,FishAgeComps_discard_seas_Type,FishLenComps_discard_seas_Type,FishAgeComps_discard_pop_seas_Type,FishLenComps_discard_pop_seas_Type
#'   Whether the operating model reports a data source once a season
#'   (`"spltSeas"`, the default) or once a year as a season total (`"aggSeas"`),
#'   one value for every fleet or one per fleet. An annual total is written into
#'   the season `seas_agg_slot` names, season one by default, with the other
#'   seasons left at zero and the observation error applied once to that total.
#'   A discard fraction is the year's discards over the year's catch, and a
#'   composition is drawn from the numbers summed over the year. An estimation
#'   model reading it should mark the same season in its `Use` array and set the
#'   matching argument in \code{\link{Setup_Mod_Catch_and_F}} or
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#' @param seas_agg_slot Named list, by data source (`"Catch"`, `"FishIdx_pop"`,
#'   `"FishAgeComps_discard"` and so on), of the season each fleet's year total is
#'   written into, an integer matrix `n_yrs x n_fleets`. A data source left out
#'   writes into season one. \code{\link{simulation_self_test}} and
#'   \code{\link{condition_closed_loop_simulations}} fill it from the fit's `Use`
#'   arrays, so a total the fit holds in season two is drawn there. Default
#'   `NULL`.
#' @param FishIdx_LikeType Error structure each fleet's index is drawn under:
#'   `"lognormal"` (0, default), `"normal"` (1) or `"mvn"` (2), matching the
#'   estimation model. An mvn fleet draws from `FishIdx_Cov` through a
#'   common-factor decomposition (see \code{\link{cov_to_factor}}) instead of
#'   `ObsFishIdx_SE`, and its population-specific data source stays lognormal.
#' @param sigmaFishIdx_form,sigmaFishIdx_pop_form How the index sd combines the
#'   reported errors with an estimated part, as \code{sigmaFishIdx_spec} in
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}: `"fix"` (0, default) draws at
#'   `ObsFishIdx_SE`, `"est_additive"` (1) at `SE + sigma`, `"est_quadrature"` (2)
#'   at `sqrt(SE^2 + sigma^2)` and `"est_replace"` (3) at `sigma`. The reported
#'   errors stay in `ObsFishIdx_SE`, so a refit reading them estimates the same part.
#' @param ln_sigmaFishIdx,ln_sigmaFishIdx_pop Log of the estimated part, one per
#'   fleet, read under a form other than `"fix"`. Default `NULL`.
#' @param FishIdx_Cov List with one element per fleet holding the fixed covariance
#'   over that fleet's fitted index observations, ordered by scanning
#'   `UseFishIdx` in array order. Required for mvn fleets. Default `NULL`.
#' @param UseFishIdx Fit flags `n_regions x n_yrs x n_seas x n_fish_fleets` from the
#'   estimation model, used to place each simulated cell in the covariance. Its
#'   year dim may be shorter than the simulation, in which case later years draw
#'   with the mean factor scale and loading. Required for mvn fleets. Default
#'   `NULL`.
#' @param t_fish Fishery index timing, `n_regions x n_seas x n_fish_fleets`, the
#'   fraction of the season elapsed when the index is observed. Numbers at age are
#'   decayed by `exp(-t_fish * ZAA)` before the index is formed, matching `t_srv`
#'   and the estimation model's own `t_fish`. Default 0.
#' @param comp_fishage_like,comp_fishlen_like,comp_fishage_pop_like,comp_fishlen_pop_like,comp_fishage_discard_like,comp_fishlen_discard_like,comp_fishage_discard_pop_like,comp_fishlen_discard_pop_like
#'   Composition likelihood per fleet for the eight fishery composition data
#'   sources: 0 = multinomial (default), 1 = Dirichlet-multinomial, 2-4 =
#'   logistic-normal, 999 = none.
#' @param ISS_FishAgeComps,ISS_FishLenComps,ISS_FishAgeComps_discard,ISS_FishLenComps_discard
#'   Input sample sizes, `n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x
#'   n_sims`. Default 100.
#' @param ISS_FishAgeComps_pop,ISS_FishLenComps_pop,ISS_FishAgeComps_discard_pop,ISS_FishLenComps_discard_pop
#'   The population-specific counterparts, with a leading `n_pop` dim. Default 100.
#' @param ln_FishAge_theta,ln_FishLen_theta,ln_FishAge_discard_theta,ln_FishLen_discard_theta
#'   Log-scale overdispersion, `n_regions x n_sexes x n_fish_fleets`. Default
#'   log(1).
#' @param ln_FishAge_theta_agg,ln_FishLen_theta_agg,ln_FishAge_discard_theta_agg,ln_FishLen_discard_theta_agg
#'   The aggregated types' counterparts, length `n_fish_fleets`. Default log(1).
#' @param ln_FishAge_pop_theta,ln_FishLen_pop_theta,ln_FishAge_discard_pop_theta,ln_FishLen_discard_pop_theta
#'   Log-scale overdispersion for the population-specific data sources, `n_pop x
#'   n_regions x n_sexes x n_fish_fleets`. Default log(1).
#' @param ln_FishAge_pop_theta_agg,ln_FishLen_pop_theta_agg,ln_FishAge_discard_pop_theta_agg,ln_FishLen_discard_pop_theta_agg
#'   Their aggregated counterparts, `n_pop x n_fish_fleets`. Default log(1).
#' @param FishAge_corr_pars,FishLen_corr_pars,FishAge_discard_corr_pars,FishLen_discard_corr_pars
#'   Correlation parameters, `n_regions x n_sexes x n_fish_fleets x 2`. Default
#'   0.01.
#' @param FishAge_corr_pars_agg,FishLen_corr_pars_agg,FishAge_discard_corr_pars_agg,FishLen_discard_corr_pars_agg
#'   Their aggregated counterparts, length `n_fish_fleets`. Default 0.01.
#' @param FishAge_pop_corr_pars,FishLen_pop_corr_pars,FishAge_discard_pop_corr_pars,FishLen_discard_pop_corr_pars
#'   Correlation parameters for the population-specific data sources, `n_pop x
#'   n_regions x n_sexes x n_fish_fleets x 2`. Default 0.01.
#' @param FishAge_pop_corr_pars_agg,FishLen_pop_corr_pars_agg,FishAge_discard_pop_corr_pars_agg,FishLen_discard_pop_corr_pars_agg
#'   Their aggregated counterparts, `n_pop x n_fish_fleets`. Default 0.01.
#' @param FishAgeComps_Type,FishLenComps_Type,FishAgeComps_pop_Type,FishLenComps_pop_Type,FishAgeComps_discard_Type,FishLenComps_discard_Type,FishAgeComps_discard_pop_Type,FishLenComps_discard_pop_Type
#'   Composition structure per year and fleet, `n_yrs x n_fish_fleets`: 0 =
#'   aggregated, 1 = split region and sex, 2 = split region joint sex (default),
#'   999 = none.
#' @param ret_sel_input Retained selectivity at age, `n_pop x n_regions x n_yrs x
#'   n_seas x n_ages x n_sexes x n_fish_fleets x n_sims`. Default 1.
#' @param FishLenComps_sel Character vector `[n_fish_fleets]`, `"age"` (default) or
#'   `"length"`, as in [Setup_Mod_FishIdx_and_Comps()]. `"length"` spreads the fish
#'   available at each age over length and selects them length by length, so the
#'   length compositions read `fish_sel_l_input`.
#' @param fish_sel_l_input,ret_sel_l_input Fishery and retention selectivity at
#'   length, `n_regions x n_yrs x n_lens x n_sexes x n_fish_fleets x n_sims`, read
#'   under `FishLenComps_sel = "length"`. `ret_sel_l_input = NULL` (default) keeps
#'   retention at age.
#' @param dmr_input Discard mortality rate, `n_regions x n_yrs x n_seas x
#'   n_fish_fleets x n_sims`. Default 0.
#' @param discard_units Discard units per fleet: 0 = abundance, 1 = biomass, 2 =
#'   abundance fraction, 3 = biomass fraction (default).
#' @param ln_sigmaD,ln_sigmaD_pop Log-scale observation sd for discards, `n_regions
#'   x n_yrs x n_seas x n_fish_fleets` with a leading `n_pop` for the second.
#'   Default log(0.02).
#' @param UseCatchAA,UseDiscardAA Integer arrays `n_regions x n_yrs x n_seas x
#'   n_obs_ages x n_sexes x n_fish_fleets`, `1` where an at-age observation is
#'   drawn. The draws sit on the observed ages from `Setup_Sim_Dim`, read through
#'   `AgeingError_fish_input` the way the estimation model reads them. The sex dim
#'   is required: a data source summed over sexes has its flag in sex slot one.
#' @param use_catch_aa,use_discard_aa Integer vectors `n_fish_fleets`, `1` for
#'   fleets whose at-age data sources are drawn.
#' @param ln_sigmaCAA,ln_sigmaDAA Log-scale observation error for the at-age data
#'   sources, `n_obs_ages x n_sexes x n_fish_fleets`. The sex dim is required.
#' @param ObsCatchAA_SE,ObsDiscardAA_SE Reported standard errors shaped like the use
#'   arrays, read only when the data source's `sigma_form` asks for them.
#' @param CatchAA_Type,DiscardAA_Type Which dims each fleet reports separately:
#'   `"agg"`, `"spltRaggS"` (default), `"aggRspltS"` or `"spltRspltS"`. A summed dim
#'   is drawn once, into slot one.
#' @param CatchAA_LikeType,DiscardAA_LikeType `"lognormal"` (default) or
#'   `"normal"`, per fleet.
#' @param CatchAA_sigma_form,DiscardAA_sigma_form Where the observation error comes
#'   from: `"none"` (default), `"data"`, `"est_additive"` or `"est_quadrature"`.
#' @param AgeObsCorr_catch,AgeObsCorr_discard How each fleet's at-age residuals are
#'   correlated, as in [Setup_Mod_Catch_and_F()]: `"iid"` (default), `"1dar1"`
#'   across ages, `"us"` across ages, or `"2dar1"` across ages and the observed
#'   years. A correlated fleet's standardized residuals are drawn before the first
#'   year, so its years and ages are drawn together.
#' @param trans_rho_catch,trans_rho_catch_year,trans_rho_catch_us,trans_rho_discard,trans_rho_discard_year,trans_rho_discard_us
#'   Unconstrained correlations under the estimation model's names and shapes:
#'   across ages and across years `[n_regions, n_sexes, n_fish_fleets]`, and the
#'   unstructured parameters `[n_pairs, n_regions, n_sexes, n_fish_fleets]` over
#'   pairs of observed ages. `NULL` (default) is zero.
#' @param CatchAA_seas_Type,DiscardAA_seas_Type,CatchAA_pop_seas_Type,DiscardAA_pop_seas_Type
#'   Per fleet, `1` where the at-age observation is a year total, as `"aggSeas"`
#'   in [Setup_Mod_Catch_and_F()]: drawn once a year from every season summed,
#'   into the season its use flags name. `0` (default) draws each season.
#' @param UseCatchAA_pop,UseDiscardAA_pop,ln_sigmaCAA_pop,ln_sigmaDAA_pop,ObsCatchAA_pop_SE,ObsDiscardAA_pop_SE,CatchAA_pop_Type,DiscardAA_pop_Type,CatchAA_pop_LikeType,DiscardAA_pop_LikeType,CatchAA_pop_sigma_form,DiscardAA_pop_sigma_form
#'   The population-specific at-age data sources, as the aggregated ones with a
#'   leading `n_pop` dim on the arrays: each population is drawn on its own from
#'   its own numbers, never summed over populations. `NULL` draws none.
#' @param AgeObsCorr_catch_pop,AgeObsCorr_discard_pop,trans_rho_catch_pop,trans_rho_catch_pop_year,trans_rho_catch_pop_us,trans_rho_discard_pop,trans_rho_discard_pop_year,trans_rho_discard_pop_us
#'   Their correlation across ages, as the aggregated settings with a leading
#'   `n_pop` dim on the parameters.
#' @param comp_fish_caal_like Conditional age-at-length likelihood per fleet:
#'   `"Multinomial"` (0), `"Dirichlet-Multinomial"` (1) or `"none"` (999, default).
#'   Only these two families exist for CAAL, since a CAAL row is the age
#'   composition of the otoliths from one length bin, usually a small and mostly
#'   zero sample.
#' @param ISS_Fish_caal Number of fish aged within each length row, `n_regions x
#'   n_yrs x n_seas x n_caal_lens x n_sexes x n_fish_fleets x n_sims`. A bin whose sample
#'   size rounds to zero is skipped. `NULL` (default) draws no CAAL; supplying it
#'   alongside a likelihood other than `"none"` is what switches `do_fish_caal` on.
#'   Requires `n_lens`.
#' @param Fish_caal_Type Composition structure per year and fleet, `n_yrs x
#'   n_fish_fleets`: `"agg"` (0) pools regions and sexes and is drawn once when the
#'   region loop reaches the last region, `"spltRspltS"` (1) draws each sex in a bin
#'   as its own sample, `"spltRjntS"` (2) draws one sample across the age by sex
#'   stack, and `"none"` (999, default) skips the fleet that year. The simulator
#'   takes the year by fleet array directly rather than the estimation model's
#'   `"CompType_Year_x-y_Fleet_z"` strings.
#' @param ln_Fish_caal_theta Log overdispersion for the Dirichlet-multinomial,
#'   `n_regions x n_sexes x n_fish_fleets`. Read under the split types, `[r, s, f]`
#'   when sexes are split and `[r, 1, f]` when they are joint, and ignored under the
#'   multinomial. Default log(1).
#' @param ln_Fish_caal_theta_agg The aggregated type's counterpart, length
#'   `n_fish_fleets`. Default log(1).
#'
#' @return A modified `sim_list` with validated fishing inputs.
#'
#' @export Setup_Sim_Fishing
#' @family Simulation Setup
Setup_Sim_Fishing <- function(sim_list,

                              # Retained / total fishery dynamics
                              ln_sigmaC = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ln_sigmaC_pop = array(log(0.02), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ln_sigmaCAA = array(log(0.2), dim = c(sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_sigmaDAA = array(log(0.2), dim = c(sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              UseCatchAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              UseDiscardAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
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
                              AgeObsCorr_catch = "iid",
                              AgeObsCorr_discard = "iid",
                              trans_rho_catch = NULL,
                              trans_rho_catch_year = NULL,
                              trans_rho_catch_us = NULL,
                              trans_rho_discard = NULL,
                              trans_rho_discard_year = NULL,
                              trans_rho_discard_us = NULL,
                              CatchAA_seas_Type = 0,
                              DiscardAA_seas_Type = 0,
                              UseCatchAA_pop = NULL,
                              UseDiscardAA_pop = NULL,
                              ln_sigmaCAA_pop = NULL,
                              ln_sigmaDAA_pop = NULL,
                              ObsCatchAA_pop_SE = NULL,
                              ObsDiscardAA_pop_SE = NULL,
                              CatchAA_pop_Type = "spltRaggS",
                              DiscardAA_pop_Type = "spltRaggS",
                              CatchAA_pop_LikeType = "lognormal",
                              DiscardAA_pop_LikeType = "lognormal",
                              CatchAA_pop_sigma_form = "none",
                              DiscardAA_pop_sigma_form = "none",
                              CatchAA_pop_seas_Type = 0,
                              DiscardAA_pop_seas_Type = 0,
                              AgeObsCorr_catch_pop = "iid",
                              AgeObsCorr_discard_pop = "iid",
                              trans_rho_catch_pop = NULL,
                              trans_rho_catch_pop_year = NULL,
                              trans_rho_catch_pop_us = NULL,
                              trans_rho_discard_pop = NULL,
                              trans_rho_discard_pop_year = NULL,
                              trans_rho_discard_pop_us = NULL,
                              catch_units = array(1, dim = c(sim_list$n_fish_fleets)),
                              init_F_val = array(0, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)),
                              Fmort_input = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)),
                              fish_sel_input,
                              fish_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ObsFishIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ObsFishIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              fish_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
                              fish_idx_ages = NULL,
                              FishIdx_LikeType = rep(0, sim_list$n_fish_fleets),
                              sigmaFishIdx_form = "fix",
                              ln_sigmaFishIdx = NULL,
                              sigmaFishIdx_pop_form = "fix",
                              ln_sigmaFishIdx_pop = NULL,
                              Catch_seas_Type = NULL,
                              Catch_pop_seas_Type = NULL,
                              FishIdx_seas_Type = NULL,
                              FishIdx_pop_seas_Type = NULL,
                              FishAgeComps_seas_Type = NULL,
                              Discard_seas_Type = NULL,
                              Discard_pop_seas_Type = NULL,
                              FishLenComps_seas_Type = NULL,
                              FishAgeComps_pop_seas_Type = NULL,
                              FishLenComps_pop_seas_Type = NULL,
                              FishAgeComps_discard_seas_Type = NULL,
                              FishLenComps_discard_seas_Type = NULL,
                              FishAgeComps_discard_pop_seas_Type = NULL,
                              FishLenComps_discard_pop_seas_Type = NULL,
                              seas_agg_slot = NULL,
                              FishIdx_Cov = NULL,
                              UseFishIdx = NULL,
                              t_fish = array(0, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)),

                              # Conditional age-at-length. Off unless an ISS array is supplied; the
                              # thetas default to log(1) and are only read by the Dirichlet-multinomial
                              comp_fish_caal_like = rep(999, sim_list$n_fish_fleets),
                              ISS_Fish_caal = NULL,
                              ln_Fish_caal_theta = NULL,
                              ln_Fish_caal_theta_agg = NULL,
                              Fish_caal_Type = array(999, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retained age compositions
                              comp_fishage_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishAge_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retained length compositions
                              comp_fishlen_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishLen_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retained age compositions (population-specific)
                              comp_fishage_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAge_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAge_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAgeComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retained length compositions (population-specific)
                              comp_fishlen_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLen_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLen_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLenComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Length compositions selected at length
                              FishLenComps_sel = rep("age", sim_list$n_fish_fleets),
                              fish_sel_l_input = NULL,
                              ret_sel_l_input = NULL,

                              # Retention and discards
                              ret_sel_input = array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              dmr_input = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)),
                              discard_units = array(3, dim = c(sim_list$n_fish_fleets)),
                              ln_sigmaD = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ln_sigmaD_pop = array(log(0.02), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),

                              # Discard age compositions (non-population specific)
                              comp_fishage_discard_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps_discard = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_discard_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_discard_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishAge_discard_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAge_discard_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishAgeComps_discard_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Discard length compositions (non-population specific)
                              comp_fishlen_discard_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps_discard = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_discard_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_discard_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishLen_discard_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLen_discard_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishLenComps_discard_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Discard age compositions (population specific)
                              comp_fishage_discard_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps_discard_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_discard_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_discard_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAge_discard_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAge_discard_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAgeComps_discard_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Discard length compositions (population specific)
                              comp_fishlen_discard_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps_discard_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_discard_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_discard_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLen_discard_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLen_discard_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLenComps_discard_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets))

                              ) {

  # Convert Options to Codes ------------------------------------------------
  # Convert character inputs to numeric codes
  # the composition likelihoods as the estimation model codes them, the miss0 forms dropping empty bins
  comp_like_codes <- list(Multinomial = 0, `Dirichlet-Multinomial` = 1,
                          `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4,
                          `iid-Logistic-Normal-miss0` = 5, `1d-Logistic-Normal-miss0` = 6, `2d-Logistic-Normal-miss0` = 7,
                          none = 999)
  catch_units <- convert_to_numeric(catch_units,  list(abd = 0, biom = 1))
  fish_idx_type <- convert_to_numeric(fish_idx_type, list(abd = 0, biom = 1, none = 999))
  FishIdx_LikeType <- convert_to_numeric(FishIdx_LikeType, list(lognormal = 0, normal = 1, mvn = 2))
  sigmaFishIdx_form <- convert_to_numeric(sigmaFishIdx_form, list(fix = 0, est_additive = 1, est_quadrature = 2, est_replace = 3))
  sigmaFishIdx_pop_form <- convert_to_numeric(sigmaFishIdx_pop_form, list(fix = 0, est_additive = 1, est_quadrature = 2, est_replace = 3))
  comp_fishage_like <- convert_to_numeric(comp_fishage_like, comp_like_codes)
  comp_fishlen_like <- convert_to_numeric(comp_fishlen_like, comp_like_codes)
  comp_fishage_pop_like <- convert_to_numeric(comp_fishage_pop_like, comp_like_codes)
  comp_fishlen_pop_like <- convert_to_numeric(comp_fishlen_pop_like, comp_like_codes)
  FishAgeComps_Type <- convert_to_numeric(FishAgeComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_Type <- convert_to_numeric(FishLenComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishAgeComps_pop_Type <- convert_to_numeric(FishAgeComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_pop_Type <- convert_to_numeric(FishLenComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  comp_fishage_discard_like <- convert_to_numeric(comp_fishage_discard_like, comp_like_codes)
  comp_fishlen_discard_like <- convert_to_numeric(comp_fishlen_discard_like, comp_like_codes)
  FishAgeComps_discard_Type <- convert_to_numeric(FishAgeComps_discard_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_discard_Type <- convert_to_numeric(FishLenComps_discard_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  comp_fishage_discard_pop_like <- convert_to_numeric(comp_fishage_discard_pop_like, comp_like_codes)
  comp_fishlen_discard_pop_like <- convert_to_numeric(comp_fishlen_discard_pop_like, comp_like_codes)
  FishAgeComps_discard_pop_Type <- convert_to_numeric(FishAgeComps_discard_pop_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_discard_pop_Type <- convert_to_numeric(FishLenComps_discard_pop_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  discard_units <- convert_to_numeric(discard_units, list(abd = 0, biom = 1, abd_frac = 2, biom_frac = 3))

  # Input Validation --------------------------------------------------------
  # a 2d logistic normal correlates bins across sexes, so it needs the composition joint by sex
  check_sim_2d_comp(comp_fishage_like, FishAgeComps_Type, "comp_fishage_like", "FishAgeComps_Type")
  check_sim_2d_comp(comp_fishlen_like, FishLenComps_Type, "comp_fishlen_like", "FishLenComps_Type")
  check_sim_2d_comp(comp_fishage_pop_like, FishAgeComps_pop_Type, "comp_fishage_pop_like", "FishAgeComps_pop_Type")
  check_sim_2d_comp(comp_fishlen_pop_like, FishLenComps_pop_Type, "comp_fishlen_pop_like", "FishLenComps_pop_Type")
  check_sim_2d_comp(comp_fishage_discard_like, FishAgeComps_discard_Type, "comp_fishage_discard_like", "FishAgeComps_discard_Type")
  check_sim_2d_comp(comp_fishlen_discard_like, FishLenComps_discard_Type, "comp_fishlen_discard_like", "FishLenComps_discard_Type")
  check_sim_2d_comp(comp_fishage_discard_pop_like, FishAgeComps_discard_pop_Type, "comp_fishage_discard_pop_like", "FishAgeComps_discard_pop_Type")
  check_sim_2d_comp(comp_fishlen_discard_pop_like, FishLenComps_discard_pop_Type, "comp_fishlen_discard_pop_like", "FishLenComps_discard_pop_Type")

  # Validate dimensions of all input parameters
  check_sim_dimensions(
    ln_sigmaC,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_sigmaC"
  )
  check_sim_dimensions(
    ln_sigmaC_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_pop = sim_list$n_pop,
    what = "ln_sigmaC_pop"
  )
  check_sim_dimensions(
    Fmort_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "Fmort_input"
  )
  check_sim_dimensions(
    fish_sel_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_ages = sim_list$n_ages,
    n_sexes = sim_list$n_sexes,
    n_pop = sim_list$n_pop,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "fish_sel_input"
  )
  check_sim_dimensions(
    fish_q_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "fish_q_input"
  )
  check_sim_dimensions(
    ObsFishIdx_SE,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ObsFishIdx_SE"
  )
  check_sim_dimensions(
    ObsFishIdx_pop_SE,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ObsFishIdx_pop_SE"
  )
  check_sim_dimensions(FishIdx_LikeType, n_fish_fleets = sim_list$n_fish_fleets, what = "FishIdx_LikeType")

  # Multivariate normal index fleets draw from the supplied covariance rather than
  # the SE array, so the covariance is validated and factor-decomposed once here.
  fish_idx_mvn <- NULL
  if(any(FishIdx_LikeType == 2)) {
    if(is.null(UseFishIdx)) stop("UseFishIdx must be supplied when any FishIdx_LikeType is mvn, to position each observation in the covariance.")
    if(length(dim(UseFishIdx)) != 4 || any(dim(UseFishIdx)[c(1,3,4)] != c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)) || dim(UseFishIdx)[2] > sim_list$n_yrs)
      stop("UseFishIdx must be an n_regions x (at most n_yrs) x n_seas x n_fish_fleets array.")
    fish_idx_mvn <- build_idx_factor(FishIdx_Cov, FishIdx_LikeType, UseFishIdx, sim_list$n_fish_fleets, "FishIdx_Cov")
  }

  # Validate fishery age composition parameters
  check_sim_dimensions(comp_fishage_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_like")
  check_sim_dimensions(
    ISS_FishAgeComps,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishAgeComps"
  )
  check_sim_dimensions(
    ln_FishAge_theta,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishAge_theta"
  )
  check_sim_dimensions(ln_FishAge_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_theta_agg")
  check_sim_dimensions(FishAge_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_corr_pars_agg")
  check_sim_dimensions(
    FishAge_corr_pars,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAge_corr_pars"
  )
  check_sim_dimensions(
    FishAgeComps_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAgeComps_Type"
  )
  check_sim_dimensions(comp_fishage_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_pop_like")
  check_sim_dimensions(
    ISS_FishAgeComps_pop,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishAgeComps_pop"
  )
  check_sim_dimensions(
    ln_FishAge_pop_theta,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishAge_pop_theta"
  )
  check_sim_dimensions(
    ln_FishAge_pop_theta_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishAge_pop_theta_agg"
  )
  check_sim_dimensions(
    FishAge_pop_corr_pars_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAge_pop_corr_pars_agg"
  )
  check_sim_dimensions(
    FishAge_pop_corr_pars,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAge_pop_corr_pars"
  )
  check_sim_dimensions(
    FishAgeComps_pop_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAgeComps_pop_Type"
  )


  # Validate fishery length composition parameters
  check_sim_dimensions(comp_fishlen_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_like")
  check_sim_dimensions(
    ISS_FishLenComps,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishLenComps"
  )
  check_sim_dimensions(
    ln_FishLen_theta,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishLen_theta"
  )
  check_sim_dimensions(ln_FishLen_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_theta_agg")
  check_sim_dimensions(FishLen_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_corr_pars_agg")
  check_sim_dimensions(
    FishLen_corr_pars,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLen_corr_pars"
  )
  check_sim_dimensions(
    FishLenComps_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLenComps_Type"
  )
  check_sim_dimensions(comp_fishlen_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_pop_like")
  check_sim_dimensions(
    ISS_FishLenComps_pop,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishLenComps_pop"
  )
  check_sim_dimensions(
    ln_FishLen_pop_theta,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishLen_pop_theta"
  )
  check_sim_dimensions(
    ln_FishLen_pop_theta_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishLen_pop_theta_agg"
  )
  check_sim_dimensions(
    FishLen_pop_corr_pars_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLen_pop_corr_pars_agg"
  )
  check_sim_dimensions(
    FishLen_pop_corr_pars,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLen_pop_corr_pars"
  )
  check_sim_dimensions(
    FishLenComps_pop_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLenComps_pop_Type"
  )


  # Validate retention and discard inputs
  check_sim_dimensions(
    ret_sel_input,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_ages = sim_list$n_ages,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ret_sel_input"
  )

  # selectivity at length checks
  if(length(FishLenComps_sel) != sim_list$n_fish_fleets || !all(FishLenComps_sel %in% c("age", "length"))) stop("FishLenComps_sel must be one of age or length for each fishery fleet")
  fish_sel_l_dim <- c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)
  if(any(FishLenComps_sel == "length")) {
    if(is.null(fish_sel_l_input) || !identical(as.numeric(dim(fish_sel_l_input)), as.numeric(fish_sel_l_dim))) stop("FishLenComps_sel = 'length' selects the length compositions at length, so fish_sel_l_input must be n_regions x n_yrs x n_lens x n_sexes x n_fish_fleets x n_sims")
    if(!is.null(ret_sel_l_input) && !identical(as.numeric(dim(ret_sel_l_input)), as.numeric(fish_sel_l_dim))) stop("ret_sel_l_input must be n_regions x n_yrs x n_lens x n_sexes x n_fish_fleets x n_sims")
  }

  # index age checks
  if(!is.null(fish_idx_ages) && (!identical(as.numeric(dim(fish_idx_ages)), as.numeric(c(sim_list$n_ages, sim_list$n_fish_fleets))) || !all(fish_idx_ages %in% c(0, 1)))) {
    stop("fish_idx_ages must be a 0/1 array n_ages x n_fish_fleets, 1 for the ages each fleet's index counts")
  }

  check_sim_dimensions(
    dmr_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "dmr_input"
  )
  check_sim_dimensions(
    ln_sigmaD,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_sigmaD"
  )
  check_sim_dimensions(
    ln_sigmaD_pop,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_sigmaD_pop"
  )

  # Validate discard age composition parameters
  check_sim_dimensions(comp_fishage_discard_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_discard_like")
  check_sim_dimensions(
    ISS_FishAgeComps_discard,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishAgeComps_discard"
  )
  check_sim_dimensions(
    ln_FishAge_discard_theta,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishAge_discard_theta"
  )
  check_sim_dimensions(ln_FishAge_discard_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_discard_theta_agg")
  check_sim_dimensions(FishAge_discard_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_discard_corr_pars_agg")
  check_sim_dimensions(
    FishAge_discard_corr_pars,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAge_discard_corr_pars"
  )
  check_sim_dimensions(
    FishAgeComps_discard_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAgeComps_discard_Type"
  )

  # Validate discard length composition parameters
  check_sim_dimensions(comp_fishlen_discard_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_discard_like")
  check_sim_dimensions(
    ISS_FishLenComps_discard,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishLenComps_discard"
  )
  check_sim_dimensions(
    ln_FishLen_discard_theta,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishLen_discard_theta"
  )
  check_sim_dimensions(ln_FishLen_discard_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_discard_theta_agg")
  check_sim_dimensions(FishLen_discard_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_discard_corr_pars_agg")
  check_sim_dimensions(
    FishLen_discard_corr_pars,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLen_discard_corr_pars"
  )
  check_sim_dimensions(
    FishLenComps_discard_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLenComps_discard_Type"
  )

  # Validate population-specific discard age composition parameters
  check_sim_dimensions(comp_fishage_discard_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_discard_pop_like")
  check_sim_dimensions(
    ISS_FishAgeComps_discard_pop,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishAgeComps_discard_pop"
  )
  check_sim_dimensions(
    ln_FishAge_discard_pop_theta,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishAge_discard_pop_theta"
  )
  check_sim_dimensions(
    ln_FishAge_discard_pop_theta_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishAge_discard_pop_theta_agg"
  )
  check_sim_dimensions(
    FishAge_discard_pop_corr_pars_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAge_discard_pop_corr_pars_agg"
  )
  check_sim_dimensions(
    FishAge_discard_pop_corr_pars,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAge_discard_pop_corr_pars"
  )
  check_sim_dimensions(
    FishAgeComps_discard_pop_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishAgeComps_discard_pop_Type"
  )

  # Validate population-specific discard length composition parameters
  check_sim_dimensions(comp_fishlen_discard_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_discard_pop_like")
  check_sim_dimensions(
    ISS_FishLenComps_discard_pop,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_FishLenComps_discard_pop"
  )
  check_sim_dimensions(
    ln_FishLen_discard_pop_theta,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishLen_discard_pop_theta"
  )
  check_sim_dimensions(
    ln_FishLen_discard_pop_theta_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "ln_FishLen_discard_pop_theta_agg"
  )
  check_sim_dimensions(
    FishLen_discard_pop_corr_pars_agg,
    n_pop = sim_list$n_pop,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLen_discard_pop_corr_pars_agg"
  )
  check_sim_dimensions(
    FishLen_discard_pop_corr_pars,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLen_discard_pop_corr_pars"
  )
  check_sim_dimensions(
    FishLenComps_discard_pop_Type,
    n_years = sim_list$n_yrs,
    n_fish_fleets = sim_list$n_fish_fleets,
    what = "FishLenComps_discard_pop_Type"
  )

  # Populate Simulation List ------------------------------------------------
  # output variables into list
  sim_list$Fmort <- Fmort_input # input fishing mortality pattern
  sim_list$catch_units <- catch_units # catch units
  sim_list$ln_sigmaC <- ln_sigmaC # Observation sd for catch
  sim_list$ln_sigmaC_pop <- ln_sigmaC_pop
  aa_use_dim <- c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages,
                  sim_list$n_sexes, sim_list$n_fish_fleets)
  aa_sigma_dim <- c(sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets)
  for(data_name in c("ln_sigmaCAA", "ln_sigmaDAA")) check_at_age_shape(get(data_name), aa_sigma_dim, data_name)
  for(data_name in c("UseCatchAA", "UseDiscardAA")) check_at_age_shape(get(data_name), aa_use_dim, data_name)
  sim_list$ln_sigmaCAA <- ln_sigmaCAA
  sim_list$ln_sigmaDAA <- ln_sigmaDAA
  sim_list$UseCatchAA <- UseCatchAA
  sim_list$UseDiscardAA <- UseDiscardAA
  # a conditioned model built before these settings existed passes them through
  # as NULL, which means the default rather than an error
  or_default <- function(x, default) if(is.null(x)) default else x

  # the population-specific data sources have a leading population dim and nothing drawn by default
  aa_pop_use_dim <- c(sim_list$n_pop, aa_use_dim)
  aa_pop_sigma_dim <- c(sim_list$n_pop, aa_sigma_dim)
  sim_list$UseCatchAA_pop <- or_default(UseCatchAA_pop, array(0, dim = aa_pop_use_dim))
  sim_list$UseDiscardAA_pop <- or_default(UseDiscardAA_pop, array(0, dim = aa_pop_use_dim))
  sim_list$ln_sigmaCAA_pop <- or_default(ln_sigmaCAA_pop, array(log(0.2), dim = aa_pop_sigma_dim))
  sim_list$ln_sigmaDAA_pop <- or_default(ln_sigmaDAA_pop, array(log(0.2), dim = aa_pop_sigma_dim))
  for(data_name in c("UseCatchAA_pop", "UseDiscardAA_pop")) check_at_age_shape(sim_list[[data_name]], aa_pop_use_dim, data_name)
  for(data_name in c("ln_sigmaCAA_pop", "ln_sigmaDAA_pop")) check_at_age_shape(sim_list[[data_name]], aa_pop_sigma_dim, data_name)

  for(data_name in c("CatchAA", "DiscardAA", "CatchAA_pop", "DiscardAA_pop")) {
    se <- get(paste0("Obs", data_name, "_SE"))
    check_at_age_shape(se, if(grepl("_pop$", data_name)) aa_pop_use_dim else aa_use_dim, paste0("Obs", data_name, "_SE"))
    use_dim <- dim(sim_list[[paste0("Use", data_name)]])
    sim_list[[paste0("Obs", data_name, "_SE")]] <- if(is.null(se) && !is.null(use_dim)) array(0, dim = use_dim) else se
    # a year total is drawn once a year from every season summed
    sim_list[[paste0(data_name, "_seas_Type")]] <- parse_seas_agg_spec(get(paste0(data_name, "_seas_Type")), paste0(data_name, "_seas_Type"), sim_list$n_fish_fleets)
    sim_list[[paste0(data_name, "_Type")]] <- at_age_type_matrix(or_default(get(paste0(data_name, "_Type")), "spltRaggS"),
                                                          sim_list$n_fish_fleets, sim_list$n_yrs,
                                                          paste0(data_name, "_Type"))
    sim_list[[paste0(data_name, "_LikeType")]] <- rep_len(convert_to_numeric(or_default(get(paste0(data_name, "_LikeType")), "lognormal"),
                                                                      list(lognormal = 0, normal = 1)), sim_list$n_fish_fleets)
    sim_list[[paste0(data_name, "_sigma_form")]] <- rep_len(convert_to_numeric(or_default(get(paste0(data_name, "_sigma_form")), "none"),
                                                                        list(none = 0, data = 1, est_additive = 2, est_quadrature = 3)),
                                                     sim_list$n_fish_fleets)
  } # end data_name loop
  sim_list$use_catch_aa <- use_catch_aa
  sim_list$use_discard_aa <- use_discard_aa
  sim_list <- sim_at_age_corr_setup(sim_list, "catch", AgeObsCorr_catch, trans_rho_catch, trans_rho_catch_year, trans_rho_catch_us, sim_list$n_fish_fleets)
  sim_list <- sim_at_age_corr_setup(sim_list, "discard", AgeObsCorr_discard, trans_rho_discard, trans_rho_discard_year, trans_rho_discard_us, sim_list$n_fish_fleets)
  sim_list <- sim_at_age_corr_setup(sim_list, "catch_pop", AgeObsCorr_catch_pop, trans_rho_catch_pop, trans_rho_catch_pop_year, trans_rho_catch_pop_us,
                                    sim_list$n_fish_fleets, pop = TRUE)
  sim_list <- sim_at_age_corr_setup(sim_list, "discard_pop", AgeObsCorr_discard_pop, trans_rho_discard_pop, trans_rho_discard_pop_year, trans_rho_discard_pop_us,
                                    sim_list$n_fish_fleets, pop = TRUE)
  sim_list$init_F <- init_F_val # initial F value
  sim_list$fish_sel <- fish_sel_input # fishery selectivity
  sim_list$fish_q <- fish_q_input # fishery catchability
  sim_list$ObsFishIdx_SE <- ObsFishIdx_SE # fishery index SE
  sim_list$ObsFishIdx_pop_SE <- ObsFishIdx_pop_SE # fishery index SE pop-specific
  sim_list$fish_idx_type <- fish_idx_type # fishery index type
  sim_list$fish_idx_ages <- if(is.null(fish_idx_ages)) array(1, dim = c(sim_list$n_ages, sim_list$n_fish_fleets)) else fish_idx_ages # ages in the index total
  sim_list$FishIdx_LikeType <- FishIdx_LikeType # fishery index error structure
  sim_list <- store_idx_sigma(sim_list, "FishIdx", sigmaFishIdx_form, ln_sigmaFishIdx, sim_list$n_fish_fleets) # estimated part of the index sd
  sim_list <- store_idx_sigma(sim_list, "FishIdx_pop", sigmaFishIdx_pop_form, ln_sigmaFishIdx_pop, sim_list$n_fish_fleets)
  if(!is.null(fish_idx_mvn)) {
    sim_list$fish_idx_mvn <- fish_idx_mvn # factor parameters for mvn index fleets
    sim_list$fish_idx_u <- matrix(NA, sim_list$n_fish_fleets, sim_list$n_sims) # shared factor draw, filled per fleet and replicate
  }
  sim_list$t_fish <- t_fish # fishery index timing within the season

  # Fishery age compositions
  sim_list$comp_fishage_like <- comp_fishage_like
  sim_list$ISS_FishAgeComps <- ISS_FishAgeComps
  sim_list$ln_FishAge_theta <- ln_FishAge_theta
  sim_list$ln_FishAge_theta_agg <- ln_FishAge_theta_agg
  sim_list$FishAge_corr_pars_agg <- FishAge_corr_pars_agg
  sim_list$FishAge_corr_pars <- FishAge_corr_pars
  sim_list$FishAgeComps_Type <- FishAgeComps_Type

  # Fishery conditional age-at-length
  comp_fish_caal_like <- convert_to_numeric(comp_fish_caal_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, none = 999))
  Fish_caal_Type <- convert_to_numeric(Fish_caal_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  if(is.null(ln_Fish_caal_theta)) ln_Fish_caal_theta <- array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets))
  if(is.null(ln_Fish_caal_theta_agg)) ln_Fish_caal_theta_agg <- rep(log(1), sim_list$n_fish_fleets)
  sim_list$do_fish_caal <- !is.null(ISS_Fish_caal) && any(comp_fish_caal_like != 999)
  if(sim_list$do_fish_caal) {
    if(is.null(sim_list$n_lens)) stop("ISS_Fish_caal was supplied, but the simulation has no length bins (n_lens is NULL)")
    n_caal_lens <- if(is.null(sim_list$n_caal_lens)) sim_list$n_lens else sim_list$n_caal_lens # age-at-length rows
    if(length(dim(ISS_Fish_caal)) != 7 || !all(dim(ISS_Fish_caal) == c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, n_caal_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)))
      stop("Dimensions of ISS_Fish_caal are not correct. Should be n_regions, n_years, n_seas, n_caal_lens, n_sexes, n_fish_fleets, and n_sims")
  }
  sim_list$comp_fish_caal_like <- comp_fish_caal_like
  sim_list$ISS_Fish_caal <- ISS_Fish_caal
  sim_list$ln_Fish_caal_theta <- ln_Fish_caal_theta
  sim_list$ln_Fish_caal_theta_agg <- ln_Fish_caal_theta_agg
  sim_list$Fish_caal_Type <- Fish_caal_Type

  # Fishery length compositions
  sim_list$comp_fishlen_like <- comp_fishlen_like
  sim_list$ISS_FishLenComps <- ISS_FishLenComps
  sim_list$ln_FishLen_theta <- ln_FishLen_theta
  sim_list$ln_FishLen_theta_agg <- ln_FishLen_theta_agg
  sim_list$FishLen_corr_pars_agg <- FishLen_corr_pars_agg
  sim_list$FishLen_corr_pars <- FishLen_corr_pars
  sim_list$FishLenComps_Type <- FishLenComps_Type

  # Population-specific stuff
  sim_list$comp_fishage_pop_like <- comp_fishage_pop_like
  sim_list$ISS_FishAgeComps_pop <- ISS_FishAgeComps_pop
  sim_list$ln_FishAge_pop_theta <- ln_FishAge_pop_theta
  sim_list$ln_FishAge_pop_theta_agg <- ln_FishAge_pop_theta_agg
  sim_list$FishAge_pop_corr_pars <- FishAge_pop_corr_pars
  sim_list$FishAge_pop_corr_pars_agg <- FishAge_pop_corr_pars_agg
  sim_list$FishAgeComps_pop_Type <- FishAgeComps_pop_Type

  sim_list$comp_fishlen_pop_like <- comp_fishlen_pop_like
  sim_list$ISS_FishLenComps_pop <- ISS_FishLenComps_pop
  sim_list$ln_FishLen_pop_theta <- ln_FishLen_pop_theta
  sim_list$ln_FishLen_pop_theta_agg <- ln_FishLen_pop_theta_agg
  sim_list$FishLen_pop_corr_pars <- FishLen_pop_corr_pars
  sim_list$FishLen_pop_corr_pars_agg <- FishLen_pop_corr_pars_agg
  sim_list$FishLenComps_pop_Type <- FishLenComps_pop_Type

  # Length compositions selected at length
  sim_list$fish_len_comp_sel <- as.numeric(FishLenComps_sel == "length") # 1 where the length comps select at length
  sim_list$fish_sel_l <- fish_sel_l_input
  sim_list$ret_sel_l <- ret_sel_l_input # NULL keeps retention at age

  # Retention and discards
  sim_list$ret_sel <- ret_sel_input
  sim_list$dmr <- dmr_input
  sim_list$discard_units <- discard_units
  sim_list$ln_sigmaD <- ln_sigmaD
  sim_list$ln_sigmaD_pop <- ln_sigmaD_pop

  # Discard age compositions
  sim_list$comp_fishage_discard_like <- comp_fishage_discard_like
  sim_list$ISS_FishAgeComps_discard <- ISS_FishAgeComps_discard
  sim_list$ln_FishAge_discard_theta <- ln_FishAge_discard_theta
  sim_list$ln_FishAge_discard_theta_agg <- ln_FishAge_discard_theta_agg
  sim_list$FishAge_discard_corr_pars <- FishAge_discard_corr_pars
  sim_list$FishAge_discard_corr_pars_agg <- FishAge_discard_corr_pars_agg
  sim_list$FishAgeComps_discard_Type <- FishAgeComps_discard_Type

  # Discard length compositions
  sim_list$comp_fishlen_discard_like <- comp_fishlen_discard_like
  sim_list$ISS_FishLenComps_discard <- ISS_FishLenComps_discard
  sim_list$ln_FishLen_discard_theta <- ln_FishLen_discard_theta
  sim_list$ln_FishLen_discard_theta_agg <- ln_FishLen_discard_theta_agg
  sim_list$FishLen_discard_corr_pars <- FishLen_discard_corr_pars
  sim_list$FishLen_discard_corr_pars_agg <- FishLen_discard_corr_pars_agg
  sim_list$FishLenComps_discard_Type <- FishLenComps_discard_Type

  # Discard age compositions (population specific)
  sim_list$comp_fishage_discard_pop_like <- comp_fishage_discard_pop_like
  sim_list$ISS_FishAgeComps_discard_pop <- ISS_FishAgeComps_discard_pop
  sim_list$ln_FishAge_discard_pop_theta <- ln_FishAge_discard_pop_theta
  sim_list$ln_FishAge_discard_pop_theta_agg <- ln_FishAge_discard_pop_theta_agg
  sim_list$FishAge_discard_pop_corr_pars <- FishAge_discard_pop_corr_pars
  sim_list$FishAge_discard_pop_corr_pars_agg <- FishAge_discard_pop_corr_pars_agg
  sim_list$FishAgeComps_discard_pop_Type <- FishAgeComps_discard_pop_Type

  # Discard length compositions (population specific)
  sim_list$comp_fishlen_discard_pop_like <- comp_fishlen_discard_pop_like
  sim_list$ISS_FishLenComps_discard_pop <- ISS_FishLenComps_discard_pop
  sim_list$ln_FishLen_discard_pop_theta <- ln_FishLen_discard_pop_theta
  sim_list$ln_FishLen_discard_pop_theta_agg <- ln_FishLen_discard_pop_theta_agg
  sim_list$FishLen_discard_pop_corr_pars <- FishLen_discard_pop_corr_pars
  sim_list$FishLen_discard_pop_corr_pars_agg <- FishLen_discard_pop_corr_pars_agg
  sim_list$FishLenComps_discard_pop_Type <- FishLenComps_discard_pop_Type

  # whether a data source reports once a season or once a year, and the season an annual total lands in
  n_fish <- sim_list$n_fish_fleets
  sim_list$Catch_seas_Type <- parse_seas_agg_spec(Catch_seas_Type, "Catch_seas_Type", n_fish)
  sim_list$Catch_pop_seas_Type <- parse_seas_agg_spec(Catch_pop_seas_Type, "Catch_pop_seas_Type", n_fish)
  sim_list$FishIdx_seas_Type <- parse_seas_agg_spec(FishIdx_seas_Type, "FishIdx_seas_Type", n_fish)
  sim_list$FishIdx_pop_seas_Type <- parse_seas_agg_spec(FishIdx_pop_seas_Type, "FishIdx_pop_seas_Type", n_fish)
  sim_list$FishAgeComps_seas_Type <- parse_seas_agg_spec(FishAgeComps_seas_Type, "FishAgeComps_seas_Type", n_fish)
  sim_list$Discard_seas_Type <- parse_seas_agg_spec(Discard_seas_Type, "Discard_seas_Type", n_fish)
  sim_list$Discard_pop_seas_Type <- parse_seas_agg_spec(Discard_pop_seas_Type, "Discard_pop_seas_Type", n_fish)
  sim_list$FishLenComps_seas_Type <- parse_seas_agg_spec(FishLenComps_seas_Type, "FishLenComps_seas_Type", n_fish)
  sim_list$FishAgeComps_pop_seas_Type <- parse_seas_agg_spec(FishAgeComps_pop_seas_Type, "FishAgeComps_pop_seas_Type", n_fish)
  sim_list$FishLenComps_pop_seas_Type <- parse_seas_agg_spec(FishLenComps_pop_seas_Type, "FishLenComps_pop_seas_Type", n_fish)
  sim_list$FishAgeComps_discard_seas_Type <- parse_seas_agg_spec(FishAgeComps_discard_seas_Type, "FishAgeComps_discard_seas_Type", n_fish)
  sim_list$FishLenComps_discard_seas_Type <- parse_seas_agg_spec(FishLenComps_discard_seas_Type, "FishLenComps_discard_seas_Type", n_fish)
  sim_list$FishAgeComps_discard_pop_seas_Type <- parse_seas_agg_spec(FishAgeComps_discard_pop_seas_Type, "FishAgeComps_discard_pop_seas_Type", n_fish)
  sim_list$FishLenComps_discard_pop_seas_Type <- parse_seas_agg_spec(FishLenComps_discard_pop_seas_Type, "FishLenComps_discard_pop_seas_Type", n_fish)
  sim_list <- store_seas_agg_slot(sim_list, seas_agg_slot, n_fish)

  # an mvn index is one draw over a covariance the season layout defines, so it cannot be collapsed
  if(any(sim_list$FishIdx_seas_Type == 1 & FishIdx_LikeType == 2))
    stop("FishIdx_seas_Type is 'aggSeas' for a fleet whose FishIdx_LikeType is 'mvn'. A multivariate ",
         "normal index is drawn once over a covariance built from the observed cells, so its seasons ",
         "cannot be summed after the draw. Use a lognormal or normal index for an annual total.")

  return(sim_list)
}

#' Set up survey parameterization for the operating model simulation
#'
#' Sets the survey catchability, selectivity, timing, index type and the age and
#' length composition settings, with their overdispersion and correlation
#' parameters. Call after \code{\link{Setup_Sim_Dim}}.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}}.
#' @param srv_sel_input Survey selectivity array \code{[n_pop x n_regions x n_yrs x
#'   n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]}. No default.
#' @param SrvLenComps_sel Character vector \code{[n_srv_fleets]}, \code{"age"}
#'   (default) or \code{"length"}, as in \code{\link{Setup_Mod_SrvIdx_and_Comps}}.
#'   \code{"length"} spreads the numbers present at each age over length and selects
#'   them length by length, so the length compositions read \code{srv_sel_l_input}.
#' @param srv_sel_l_input Survey selectivity at length \code{[n_regions × n_yrs ×
#'   n_lens × n_sexes × n_srv_fleets × n_sims]}, read under
#'   \code{SrvLenComps_sel = "length"}.
#' @param srv_q_input Survey catchability array \code{[n_regions × n_yrs ×
#'   n_srv_fleets × n_sims]}. Default 1.
#' @param ObsSrvIdx_SE,ObsSrvIdx_pop_SE Lognormal observation error sd for the
#'   survey indices, \code{[n_regions × n_yrs × n_seas × n_srv_fleets]} with a
#'   leading \code{n_pop} for the second. Default 0.2.
#' @param t_srv Survey timing as a fraction of the year or season, \code{[n_regions
#'   × n_seas × n_srv_fleets]}. Default 1.
#' @param srv_idx_type Index type per fleet: 0/\code{"abd"} or 1/\code{"biom"}
#'   (default).
#' @param srv_idx_ages Ages counted in each fleet's index total, a 0/1 array
#'   \code{[n_ages × n_srv_fleets]}, the estimation model's \code{srv_idx_ages}.
#'   \code{NULL} (default) counts every age.
#' @param SrvIdx_seas_Type,SrvIdx_pop_seas_Type,SrvAgeComps_seas_Type,SrvLenComps_seas_Type,SrvAgeComps_pop_seas_Type,SrvLenComps_pop_seas_Type
#'   Whether the operating model reports a survey data source once a season
#'   (\code{"spltSeas"}, the default) or once a year as a season total
#'   (\code{"aggSeas"}), one value for every survey or one per survey. An annual
#'   total is written into the season \code{seas_agg_slot} names, season one by
#'   default, with the other seasons left at zero and the observation error
#'   applied once to that total. A composition is drawn from the numbers summed
#'   over the year. An estimation model reading it should mark the same season in
#'   its \code{Use} array and set the matching argument in
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}.
#' @param seas_agg_slot Named list, by survey data source (\code{"SrvIdx"},
#'   \code{"SrvAgeComps_pop"} and so on), of the season each survey's year total
#'   is written into, an integer matrix \code{[n_yrs x n_srv_fleets]}. See
#'   \code{\link{Setup_Sim_Fishing}}. Default \code{NULL}.
#' @param SrvIdx_LikeType Error structure each fleet's index is drawn under:
#'   \code{"lognormal"} (0, default), \code{"normal"} (1) or \code{"mvn"} (2),
#'   matching the estimation model. An mvn fleet draws from \code{SrvIdx_Cov}
#'   through a common-factor decomposition (see \code{\link{cov_to_factor}})
#'   instead of \code{ObsSrvIdx_SE}, and its population-specific data source stays
#'   lognormal.
#' @param sigmaSrvIdx_form,sigmaSrvIdx_pop_form How the index sd combines the
#'   reported errors with an estimated part, as \code{sigmaSrvIdx_spec} in
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}: \code{"fix"} (0, default),
#'   \code{"est_additive"} (1), \code{"est_quadrature"} (2) or
#'   \code{"est_replace"} (3). See \code{\link{Setup_Sim_Fishing}}.
#' @param ln_sigmaSrvIdx,ln_sigmaSrvIdx_pop Log of the estimated part, one per
#'   survey, read under a form other than \code{"fix"}. Default \code{NULL}.
#' @param SrvIdx_Cov List with one element per fleet holding the fixed covariance
#'   over that fleet's fitted index observations, ordered by scanning
#'   \code{UseSrvIdx} in array order. Required for mvn fleets. Default \code{NULL}.
#' @param UseSrvIdx Fit flags \code{[n_regions x n_yrs x n_seas x n_srv_fleets]}
#'   from the estimation model, used to place each simulated cell in the
#'   covariance. Its year dim may be shorter than the simulation, in which case
#'   later years draw with the mean factor scale and loading. Required for mvn
#'   fleets. Default \code{NULL}.
#' @param comp_srvage_like,comp_srvlen_like,comp_srvage_pop_like,comp_srvlen_pop_like
#'   Composition likelihood per fleet for the four survey composition data sources:
#'   0/\code{"Multinomial"} (default), 1/\code{"Dirichlet-Multinomial"},
#'   2/\code{"iid-Logistic-Normal"}, 3/\code{"1d-Logistic-Normal"} or
#'   4/\code{"2d-Logistic-Normal"}.
#' @param ISS_SrvAgeComps,ISS_SrvLenComps Input sample sizes \code{[n_regions ×
#'   n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]}. Default 100.
#' @param ISS_SrvAgeComps_pop,ISS_SrvLenComps_pop The population-specific
#'   counterparts, with a leading \code{n_pop} dim. Default 100.
#' @param ln_SrvAge_theta,ln_SrvLen_theta Log-scale overdispersion \code{[n_regions
#'   × n_sexes × n_srv_fleets]}, read under likelihoods 1-4. Default log(1).
#' @param ln_SrvAge_theta_agg,ln_SrvLen_theta_agg The aggregated types'
#'   counterparts, length \code{n_srv_fleets}. Default log(1).
#' @param ln_SrvAge_pop_theta,ln_SrvLen_pop_theta Log-scale overdispersion for the
#'   population-specific data sources \code{[n_pop × n_regions × n_sexes ×
#'   n_srv_fleets]}. Default log(1).
#' @param ln_SrvAge_pop_theta_agg,ln_SrvLen_pop_theta_agg Their aggregated
#'   counterparts \code{[n_pop × n_srv_fleets]}. Default log(1).
#' @param SrvAge_corr_pars,SrvLen_corr_pars Correlation parameters \code{[n_regions
#'   × n_sexes × n_srv_fleets × 2]}, the age AR1 and the sex correlation, read
#'   under likelihoods 3 and 4. Default 0.01.
#' @param SrvAge_corr_pars_agg,SrvLen_corr_pars_agg The aggregated types'
#'   counterparts, length \code{n_srv_fleets}, read under likelihood 3. Default
#'   0.01.
#' @param SrvAge_pop_corr_pars,SrvLen_pop_corr_pars Correlation parameters for the
#'   population-specific data sources \code{[n_pop × n_regions × n_sexes ×
#'   n_srv_fleets × 2]}. Default 0.01.
#' @param SrvAge_pop_corr_pars_agg,SrvLen_pop_corr_pars_agg Their aggregated
#'   counterparts \code{[n_pop × n_srv_fleets]}. Default 0.01.
#' @param SrvAgeComps_Type,SrvLenComps_Type,SrvAgeComps_pop_Type,SrvLenComps_pop_Type
#'   Composition structure \code{[n_yrs × n_srv_fleets]}: 0/\code{"agg"},
#'   1/\code{"spltRspltS"}, 2/\code{"spltRjntS"} (default) or 999/\code{"none"}.
#' @param UseSrvIdxAA Integer array `n_regions x n_yrs x n_seas x n_obs_ages x
#'   n_sexes x n_srv_fleets`, `1` where a survey index at age is drawn, on the
#'   observed ages `AgeingError_srv_input` reads onto. The sex dim is required: a
#'   data source summed over sexes has its flag in sex slot one.
#' @param use_srv_idx_aa Integer vector `n_srv_fleets`, `1` for fleets whose index
#'   at age is drawn.
#' @param ln_sigmaSrvIdxAA Log-scale observation error for the index at age,
#'   `n_obs_ages x n_sexes x n_srv_fleets`. The sex dim is required.
#' @param ObsSrvIdxAA_SE Reported standard errors shaped like `UseSrvIdxAA`, read
#'   only when `SrvIdxAA_sigma_form` asks for them.
#' @param SrvIdxAA_Type Which dims each fleet reports separately: `"agg"`,
#'   `"spltRaggS"` (default), `"aggRspltS"` or `"spltRspltS"`.
#' @param SrvIdxAA_LikeType `"lognormal"` (default) or `"normal"`, per fleet.
#' @param SrvIdxAA_sigma_form Where the observation error comes from: `"none"`
#'   (default), `"data"`, `"est_additive"` or `"est_quadrature"`.
#' @param AgeObsCorr_srv_idx How each fleet's index-at-age residuals are
#'   correlated, as in [Setup_Mod_SrvIdx_and_Comps()]: `"iid"` (default),
#'   `"1dar1"`, `"us"` or `"2dar1"`; see [Setup_Sim_Fishing()].
#' @param trans_rho_srv_idx,trans_rho_srv_idx_year,trans_rho_srv_idx_us
#'   Unconstrained correlations under the estimation model's names and shapes,
#'   as for the fishery's. `NULL` (default) is zero.
#' @param SrvIdxAA_seas_Type,SrvIdxAA_pop_seas_Type Per fleet, `1` where the index at
#'   age is a year total, drawn once a year from every season summed; see
#'   [Setup_Sim_Fishing()]. `0` (default) draws each season.
#' @param UseSrvIdxAA_pop,ln_sigmaSrvIdxAA_pop,ObsSrvIdxAA_pop_SE,SrvIdxAA_pop_Type,SrvIdxAA_pop_LikeType,SrvIdxAA_pop_sigma_form,AgeObsCorr_srv_idx_pop,trans_rho_srv_idx_pop,trans_rho_srv_idx_pop_year,trans_rho_srv_idx_pop_us
#'   The population-specific index at age, as the aggregated one with a leading
#'   `n_pop` dim on the arrays and parameters, each population drawn on its own.
#'   `NULL` draws none.
#' @param comp_srv_caal_like Conditional age-at-length likelihood per fleet:
#'   `"Multinomial"` (0), `"Dirichlet-Multinomial"` (1) or `"none"` (999, default).
#'   The survey twin of `comp_fish_caal_like`, and only these two families exist
#'   for CAAL.
#' @param ISS_Srv_caal Number of fish aged within each length row, `n_regions x
#'   n_yrs x n_seas x n_caal_lens x n_sexes x n_srv_fleets x n_sims`. A bin whose sample
#'   size rounds to zero is skipped. `NULL` (default) draws no CAAL; supplying it
#'   alongside a likelihood other than `"none"` switches `do_srv_caal` on. Requires
#'   `n_lens`.
#' @param Srv_caal_Type Composition structure per year and fleet, `n_yrs x
#'   n_srv_fleets`, with the codes of `Fish_caal_Type`: `"agg"` (0),
#'   `"spltRspltS"` (1), `"spltRjntS"` (2) or `"none"` (999, default). The
#'   simulator takes the year by fleet array directly.
#' @param ln_Srv_caal_theta Log overdispersion for the Dirichlet-multinomial,
#'   `n_regions x n_sexes x n_srv_fleets`, read under the split types and ignored
#'   under the multinomial. Default log(1).
#' @param ln_Srv_caal_theta_agg The aggregated type's counterpart, length
#'   `n_srv_fleets`. Default log(1).
#'
#' @return \code{sim_list} with the survey fields appended: \code{$srv_sel},
#'   \code{$srv_q}, \code{$ObsSrvIdx_SE}, \code{$ObsSrvIdx_pop_SE}, \code{$t_srv},
#'   \code{$srv_idx_type}, and, for each of the four composition data sources, its
#'   likelihood, input sample sizes, overdispersion, correlation parameters and
#'   composition type. Character codes are converted to integers before storage.
#'
#' @export Setup_Sim_Survey
#' @family Simulation Setup
Setup_Sim_Survey <- function(sim_list,
                             srv_sel_input,
                             ObsSrvIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,  sim_list$n_srv_fleets)),
                             ln_sigmaSrvIdxAA = array(log(0.2), dim = c(sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             UseSrvIdxAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ObsSrvIdxAA_SE = NULL,
                             SrvIdxAA_Type = "spltRaggS",
                             SrvIdxAA_LikeType = "lognormal",
                             SrvIdxAA_sigma_form = "none",
                             use_srv_idx_aa = rep(0, sim_list$n_srv_fleets),
                             AgeObsCorr_srv_idx = "iid",
                             trans_rho_srv_idx = NULL,
                             trans_rho_srv_idx_year = NULL,
                             trans_rho_srv_idx_us = NULL,
                             SrvIdxAA_seas_Type = 0,
                             UseSrvIdxAA_pop = NULL,
                             ln_sigmaSrvIdxAA_pop = NULL,
                             ObsSrvIdxAA_pop_SE = NULL,
                             SrvIdxAA_pop_Type = "spltRaggS",
                             SrvIdxAA_pop_LikeType = "lognormal",
                             SrvIdxAA_pop_sigma_form = "none",
                             SrvIdxAA_pop_seas_Type = 0,
                             AgeObsCorr_srv_idx_pop = "iid",
                             trans_rho_srv_idx_pop = NULL,
                             trans_rho_srv_idx_pop_year = NULL,
                             trans_rho_srv_idx_pop_us = NULL,
                             ObsSrvIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,  sim_list$n_srv_fleets)),
                             srv_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_srv_fleets, sim_list$n_sims)),
                             t_srv = array(1, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets)),
                             srv_idx_type = array(1, dim = c(sim_list$n_srv_fleets)),
                             srv_idx_ages = NULL,
                             SrvIdx_LikeType = rep(0, sim_list$n_srv_fleets),
                             sigmaSrvIdx_form = "fix",
                             ln_sigmaSrvIdx = NULL,
                             sigmaSrvIdx_pop_form = "fix",
                             ln_sigmaSrvIdx_pop = NULL,
                             SrvIdx_seas_Type = NULL,
                             SrvIdx_pop_seas_Type = NULL,
                             SrvAgeComps_seas_Type = NULL,
                             SrvLenComps_seas_Type = NULL,
                             SrvAgeComps_pop_seas_Type = NULL,
                             SrvLenComps_pop_seas_Type = NULL,
                             seas_agg_slot = NULL,
                             SrvIdx_Cov = NULL,
                             UseSrvIdx = NULL,
                             comp_srv_caal_like = rep(999, sim_list$n_srv_fleets),
                             ISS_Srv_caal = NULL,
                             ln_Srv_caal_theta = NULL,
                             ln_Srv_caal_theta_agg = NULL,
                             Srv_caal_Type = array(999, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvage_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvAge_theta_agg = rep(log(1), sim_list$n_srv_fleets),
                             SrvAge_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
                             SrvAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvlen_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvLen_theta_agg = rep(log(1), sim_list$n_srv_fleets),
                             SrvLen_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
                             SrvLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvage_pop_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvAgeComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvAge_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvAge_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvAge_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvAge_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvAgeComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvlen_pop_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvLenComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvLen_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvLen_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvLen_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvLen_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvLenComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),

                             # Length compositions selected at length
                             SrvLenComps_sel = rep("age", sim_list$n_srv_fleets),
                             srv_sel_l_input = NULL
                             ) {

  # Convert Options to Codes ------------------------------------------------
  # Convert character inputs to numeric codes
  # the composition likelihoods as the estimation model codes them, the miss0 forms dropping empty bins
  comp_like_codes <- list(Multinomial = 0, `Dirichlet-Multinomial` = 1,
                          `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4,
                          `iid-Logistic-Normal-miss0` = 5, `1d-Logistic-Normal-miss0` = 6, `2d-Logistic-Normal-miss0` = 7,
                          none = 999)
  srv_idx_type <- convert_to_numeric(srv_idx_type, list(abd = 0, biom = 1, none = 999))
  SrvIdx_LikeType <- convert_to_numeric(SrvIdx_LikeType, list(lognormal = 0, normal = 1, mvn = 2))
  sigmaSrvIdx_form <- convert_to_numeric(sigmaSrvIdx_form, list(fix = 0, est_additive = 1, est_quadrature = 2, est_replace = 3))
  sigmaSrvIdx_pop_form <- convert_to_numeric(sigmaSrvIdx_pop_form, list(fix = 0, est_additive = 1, est_quadrature = 2, est_replace = 3))
  comp_srvage_like <- convert_to_numeric(comp_srvage_like, comp_like_codes)
  comp_srvlen_like <- convert_to_numeric(comp_srvlen_like, comp_like_codes)
  SrvAgeComps_Type <- convert_to_numeric(SrvAgeComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  SrvLenComps_Type <- convert_to_numeric(SrvLenComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  SrvAgeComps_pop_Type <- convert_to_numeric(SrvAgeComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  SrvLenComps_pop_Type <- convert_to_numeric(SrvLenComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  comp_srvage_pop_like <- convert_to_numeric(comp_srvage_pop_like, comp_like_codes)
  comp_srvlen_pop_like <- convert_to_numeric(comp_srvlen_pop_like, comp_like_codes)

  # Input Validation --------------------------------------------------------
  # a 2d logistic normal correlates bins across sexes, so it needs the composition joint by sex
  check_sim_2d_comp(comp_srvage_like, SrvAgeComps_Type, "comp_srvage_like", "SrvAgeComps_Type")
  check_sim_2d_comp(comp_srvlen_like, SrvLenComps_Type, "comp_srvlen_like", "SrvLenComps_Type")
  check_sim_2d_comp(comp_srvage_pop_like, SrvAgeComps_pop_Type, "comp_srvage_pop_like", "SrvAgeComps_pop_Type")
  check_sim_2d_comp(comp_srvlen_pop_like, SrvLenComps_pop_Type, "comp_srvlen_pop_like", "SrvLenComps_pop_Type")

  # Validate dimensions of all input parameters
  check_sim_dimensions(
    srv_sel_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_ages = sim_list$n_ages,
    n_sexes = sim_list$n_sexes,
    n_pop = sim_list$n_pop,
    n_seas = sim_list$n_seas,
    n_srv_fleets = sim_list$n_srv_fleets,
    n_sims = sim_list$n_sims,
    what = "srv_sel_input"
  )

  # selex at length checks
  if(length(SrvLenComps_sel) != sim_list$n_srv_fleets || !all(SrvLenComps_sel %in% c("age", "length"))) stop("SrvLenComps_sel must be one of age or length for each survey fleet")
  srv_sel_l_dim <- c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)
  if(any(SrvLenComps_sel == "length") && (is.null(srv_sel_l_input) || !identical(as.numeric(dim(srv_sel_l_input)), as.numeric(srv_sel_l_dim)))) {
    stop("SrvLenComps_sel = 'length' selects the length compositions at length, so srv_sel_l_input must be n_regions x n_yrs x n_lens x n_sexes x n_srv_fleets x n_sims")
  }

  # index age checks
  if(!is.null(srv_idx_ages) && (!identical(as.numeric(dim(srv_idx_ages)), as.numeric(c(sim_list$n_ages, sim_list$n_srv_fleets))) || !all(srv_idx_ages %in% c(0, 1)))) {
    stop("srv_idx_ages must be a 0/1 array n_ages x n_srv_fleets, 1 for the ages each fleet's index counts")
  }

  check_sim_dimensions(
    srv_q_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_srv_fleets = sim_list$n_srv_fleets,
    n_sims = sim_list$n_sims,
    what = "srv_q_input"
  )
  check_sim_dimensions(
    ObsSrvIdx_SE,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ObsSrvIdx_SE"
  )
  check_sim_dimensions(
    ObsSrvIdx_pop_SE,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ObsSrvIdx_pop_SE"
  )
  check_sim_dimensions(
    t_srv,
    n_regions = sim_list$n_regions,
    n_seas = sim_list$n_seas,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "t_srv"
  )
  check_sim_dimensions(SrvIdx_LikeType, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvIdx_LikeType")

  # Multivariate normal index fleets draw from the supplied covariance rather than
  # the SE array, so the covariance is validated and factor-decomposed once here.
  srv_idx_mvn <- NULL
  if(any(SrvIdx_LikeType == 2)) {
    if(is.null(UseSrvIdx)) stop("UseSrvIdx must be supplied when any SrvIdx_LikeType is mvn, to position each observation in the covariance.")
    if(length(dim(UseSrvIdx)) != 4 || any(dim(UseSrvIdx)[c(1,3,4)] != c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets)) || dim(UseSrvIdx)[2] > sim_list$n_yrs)
      stop("UseSrvIdx must be an n_regions x (at most n_yrs) x n_seas x n_srv_fleets array.")
    srv_idx_mvn <- build_idx_factor(SrvIdx_Cov, SrvIdx_LikeType, UseSrvIdx, sim_list$n_srv_fleets, "SrvIdx_Cov")
  }

  # Validate survey age composition parameters
  check_sim_dimensions(comp_srvage_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvage_like")
  check_sim_dimensions(
    ISS_SrvAgeComps,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_SrvAgeComps"
  )
  check_sim_dimensions(
    ln_SrvAge_theta,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ln_SrvAge_theta"
  )
  check_sim_dimensions(ln_SrvAge_theta_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_theta_agg")
  check_sim_dimensions(SrvAge_corr_pars_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_corr_pars_agg")
  check_sim_dimensions(
    SrvAge_corr_pars,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvAge_corr_pars"
  )
  check_sim_dimensions(
    SrvAgeComps_Type,
    n_years = sim_list$n_yrs,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvAgeComps_Type"
  )
  check_sim_dimensions(comp_srvage_pop_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvage_pop_like")
  check_sim_dimensions(
    ISS_SrvAgeComps_pop,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_SrvAgeComps_pop"
  )
  check_sim_dimensions(
    ln_SrvAge_pop_theta,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ln_SrvAge_pop_theta"
  )
  check_sim_dimensions(
    ln_SrvAge_pop_theta_agg,
    n_pop = sim_list$n_pop,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ln_SrvAge_pop_theta_agg"
  )
  check_sim_dimensions(
    SrvAge_pop_corr_pars_agg,
    n_pop = sim_list$n_pop,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvAge_pop_corr_pars_agg"
  )
  check_sim_dimensions(
    SrvAge_pop_corr_pars,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvAge_pop_corr_pars"
  )
  check_sim_dimensions(
    SrvAgeComps_pop_Type,
    n_years = sim_list$n_yrs,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvAgeComps_pop_Type"
  )


  # Validate survey length composition parameters
  check_sim_dimensions(comp_srvlen_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvlen_like")
  check_sim_dimensions(
    ISS_SrvLenComps,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_SrvLenComps"
  )
  check_sim_dimensions(
    ln_SrvLen_theta,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ln_SrvLen_theta"
  )
  check_sim_dimensions(ln_SrvLen_theta_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_theta_agg")
  check_sim_dimensions(SrvLen_corr_pars_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_corr_pars_agg")
  check_sim_dimensions(
    SrvLen_corr_pars,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvLen_corr_pars"
  )
  check_sim_dimensions(
    SrvLenComps_Type,
    n_years = sim_list$n_yrs,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvLenComps_Type"
  )
  check_sim_dimensions(comp_srvlen_pop_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvlen_pop_like")
  check_sim_dimensions(
    ISS_SrvLenComps_pop,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_seas = sim_list$n_seas,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    n_sims = sim_list$n_sims,
    what = "ISS_SrvLenComps_pop"
  )
  check_sim_dimensions(
    ln_SrvLen_pop_theta,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ln_SrvLen_pop_theta"
  )
  check_sim_dimensions(
    ln_SrvLen_pop_theta_agg,
    n_pop = sim_list$n_pop,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "ln_SrvLen_pop_theta_agg"
  )
  check_sim_dimensions(
    SrvLen_pop_corr_pars_agg,
    n_pop = sim_list$n_pop,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvLen_pop_corr_pars_agg"
  )
  check_sim_dimensions(
    SrvLen_pop_corr_pars,
    n_pop = sim_list$n_pop,
    n_regions = sim_list$n_regions,
    n_sexes = sim_list$n_sexes,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvLen_pop_corr_pars"
  )
  check_sim_dimensions(
    SrvLenComps_pop_Type,
    n_years = sim_list$n_yrs,
    n_srv_fleets = sim_list$n_srv_fleets,
    what = "SrvLenComps_pop_Type"
  )

  # Populate Simulation List ------------------------------------------------
  # output into list
  sim_list$srv_sel <- srv_sel_input
  sim_list$srv_len_comp_sel <- as.numeric(SrvLenComps_sel == "length") # 1 where the length comps select at length
  sim_list$srv_sel_l <- srv_sel_l_input
  sim_list$srv_q <- srv_q_input
  sim_list$ObsSrvIdx_SE <- ObsSrvIdx_SE
  srv_aa_dim <- c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages,
                  sim_list$n_sexes, sim_list$n_srv_fleets) # on the observed ages the ageing error reads onto
  check_at_age_shape(ln_sigmaSrvIdxAA, c(sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets), "ln_sigmaSrvIdxAA")
  check_at_age_shape(UseSrvIdxAA, srv_aa_dim, "UseSrvIdxAA")
  check_at_age_shape(ObsSrvIdxAA_SE, srv_aa_dim, "ObsSrvIdxAA_SE")
  sim_list$ln_sigmaSrvIdxAA <- ln_sigmaSrvIdxAA
  sim_list$UseSrvIdxAA <- UseSrvIdxAA
  sim_list$ObsSrvIdxAA_SE <- if(is.null(ObsSrvIdxAA_SE) && !is.null(dim(sim_list$UseSrvIdxAA))) array(0, dim = dim(sim_list$UseSrvIdxAA))
                             else ObsSrvIdxAA_SE
  # a conditioned model built before these settings existed passes them through
  # as NULL, which means the default rather than an error
  if(is.null(SrvIdxAA_Type)) SrvIdxAA_Type <- "spltRaggS"
  if(is.null(SrvIdxAA_LikeType)) SrvIdxAA_LikeType <- "lognormal"
  if(is.null(SrvIdxAA_sigma_form)) SrvIdxAA_sigma_form <- "none"
  sim_list$SrvIdxAA_Type <- at_age_type_matrix(SrvIdxAA_Type, sim_list$n_srv_fleets, sim_list$n_yrs, "SrvIdxAA_Type")
  sim_list$SrvIdxAA_LikeType <- rep_len(convert_to_numeric(SrvIdxAA_LikeType, list(lognormal = 0, normal = 1)), sim_list$n_srv_fleets)
  sim_list$SrvIdxAA_sigma_form <- rep_len(convert_to_numeric(SrvIdxAA_sigma_form, list(none = 0, data = 1, est_additive = 2, est_quadrature = 3)), sim_list$n_srv_fleets)
  sim_list$use_srv_idx_aa <- use_srv_idx_aa
  sim_list <- sim_at_age_corr_setup(sim_list, "srv_idx", AgeObsCorr_srv_idx, trans_rho_srv_idx, trans_rho_srv_idx_year, trans_rho_srv_idx_us, sim_list$n_srv_fleets)
  sim_list$SrvIdxAA_seas_Type <- parse_seas_agg_spec(SrvIdxAA_seas_Type, "SrvIdxAA_seas_Type", sim_list$n_srv_fleets)

  # the population-specific index at age, with a leading population dim and nothing drawn by default
  srv_aa_pop_dim <- c(sim_list$n_pop, srv_aa_dim)
  srv_sigma_pop_dim <- c(sim_list$n_pop, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets)
  sim_list$UseSrvIdxAA_pop <- if(is.null(UseSrvIdxAA_pop)) array(0, dim = srv_aa_pop_dim) else UseSrvIdxAA_pop
  sim_list$ln_sigmaSrvIdxAA_pop <- if(is.null(ln_sigmaSrvIdxAA_pop)) array(log(0.2), dim = srv_sigma_pop_dim) else ln_sigmaSrvIdxAA_pop
  sim_list$ObsSrvIdxAA_pop_SE <- if(is.null(ObsSrvIdxAA_pop_SE)) array(0, dim = srv_aa_pop_dim) else ObsSrvIdxAA_pop_SE
  check_at_age_shape(sim_list$UseSrvIdxAA_pop, srv_aa_pop_dim, "UseSrvIdxAA_pop")
  check_at_age_shape(sim_list$ln_sigmaSrvIdxAA_pop, srv_sigma_pop_dim, "ln_sigmaSrvIdxAA_pop")
  check_at_age_shape(sim_list$ObsSrvIdxAA_pop_SE, srv_aa_pop_dim, "ObsSrvIdxAA_pop_SE")
  sim_list$SrvIdxAA_pop_Type <- at_age_type_matrix(if(is.null(SrvIdxAA_pop_Type)) "spltRaggS" else SrvIdxAA_pop_Type, sim_list$n_srv_fleets, sim_list$n_yrs, "SrvIdxAA_pop_Type")
  sim_list$SrvIdxAA_pop_LikeType <- rep_len(convert_to_numeric(if(is.null(SrvIdxAA_pop_LikeType)) "lognormal" else SrvIdxAA_pop_LikeType,
                                                               list(lognormal = 0, normal = 1)), sim_list$n_srv_fleets)
  sim_list$SrvIdxAA_pop_sigma_form <- rep_len(convert_to_numeric(if(is.null(SrvIdxAA_pop_sigma_form)) "none" else SrvIdxAA_pop_sigma_form,
                                                                 list(none = 0, data = 1, est_additive = 2, est_quadrature = 3)), sim_list$n_srv_fleets)
  sim_list$SrvIdxAA_pop_seas_Type <- parse_seas_agg_spec(SrvIdxAA_pop_seas_Type, "SrvIdxAA_pop_seas_Type", sim_list$n_srv_fleets)
  sim_list <- sim_at_age_corr_setup(sim_list, "srv_idx_pop", AgeObsCorr_srv_idx_pop, trans_rho_srv_idx_pop, trans_rho_srv_idx_pop_year, trans_rho_srv_idx_pop_us,
                                    sim_list$n_srv_fleets, pop = TRUE)
  sim_list$ObsSrvIdx_pop_SE <- ObsSrvIdx_pop_SE
  sim_list$t_srv <- t_srv
  sim_list$srv_idx_type <- srv_idx_type
  sim_list$srv_idx_ages <- if(is.null(srv_idx_ages)) array(1, dim = c(sim_list$n_ages, sim_list$n_srv_fleets)) else srv_idx_ages # ages in the index total
  sim_list$SrvIdx_LikeType <- SrvIdx_LikeType # survey index error structure
  sim_list <- store_idx_sigma(sim_list, "SrvIdx", sigmaSrvIdx_form, ln_sigmaSrvIdx, sim_list$n_srv_fleets) # estimated part of the index sd
  sim_list <- store_idx_sigma(sim_list, "SrvIdx_pop", sigmaSrvIdx_pop_form, ln_sigmaSrvIdx_pop, sim_list$n_srv_fleets)
  if(!is.null(srv_idx_mvn)) {
    sim_list$srv_idx_mvn <- srv_idx_mvn # factor parameters for mvn index fleets
    sim_list$srv_idx_u <- matrix(NA, sim_list$n_srv_fleets, sim_list$n_sims) # shared factor draw, filled per fleet and replicate
  }

  # Survey age compositions
  sim_list$comp_srvage_like <- comp_srvage_like
  sim_list$ISS_SrvAgeComps <- ISS_SrvAgeComps
  sim_list$ln_SrvAge_theta <- ln_SrvAge_theta
  sim_list$ln_SrvAge_theta_agg <- ln_SrvAge_theta_agg
  sim_list$SrvAge_corr_pars_agg <- SrvAge_corr_pars_agg
  sim_list$SrvAge_corr_pars <- SrvAge_corr_pars
  sim_list$SrvAgeComps_Type <- SrvAgeComps_Type

  # Survey conditional age-at-length
  comp_srv_caal_like <- convert_to_numeric(comp_srv_caal_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, none = 999))
  Srv_caal_Type <- convert_to_numeric(Srv_caal_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  if(is.null(ln_Srv_caal_theta)) ln_Srv_caal_theta <- array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets))
  if(is.null(ln_Srv_caal_theta_agg)) ln_Srv_caal_theta_agg <- rep(log(1), sim_list$n_srv_fleets)
  sim_list$do_srv_caal <- !is.null(ISS_Srv_caal) && any(comp_srv_caal_like != 999)
  if(sim_list$do_srv_caal) {
    if(is.null(sim_list$n_lens)) stop("ISS_Srv_caal was supplied, but the simulation has no length bins (n_lens is NULL)")
    n_caal_lens <- if(is.null(sim_list$n_caal_lens)) sim_list$n_lens else sim_list$n_caal_lens # age-at-length rows
    if(length(dim(ISS_Srv_caal)) != 7 || !all(dim(ISS_Srv_caal) == c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, n_caal_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)))
      stop("Dimensions of ISS_Srv_caal are not correct. Should be n_regions, n_years, n_seas, n_caal_lens, n_sexes, n_srv_fleets, and n_sims")
  }
  sim_list$comp_srv_caal_like <- comp_srv_caal_like
  sim_list$ISS_Srv_caal <- ISS_Srv_caal
  sim_list$ln_Srv_caal_theta <- ln_Srv_caal_theta
  sim_list$ln_Srv_caal_theta_agg <- ln_Srv_caal_theta_agg
  sim_list$Srv_caal_Type <- Srv_caal_Type

  # Survey length compositions
  sim_list$comp_srvlen_like <- comp_srvlen_like
  sim_list$ISS_SrvLenComps <- ISS_SrvLenComps
  sim_list$ln_SrvLen_theta <- ln_SrvLen_theta
  sim_list$ln_SrvLen_theta_agg <- ln_SrvLen_theta_agg
  sim_list$SrvLen_corr_pars_agg <- SrvLen_corr_pars_agg
  sim_list$SrvLen_corr_pars <- SrvLen_corr_pars
  sim_list$SrvLenComps_Type <- SrvLenComps_Type

  # Population-specific stuff
  sim_list$comp_srvage_pop_like <- comp_srvage_pop_like
  sim_list$ISS_SrvAgeComps_pop <- ISS_SrvAgeComps_pop
  sim_list$ln_SrvAge_pop_theta <- ln_SrvAge_pop_theta
  sim_list$ln_SrvAge_pop_theta_agg <- ln_SrvAge_pop_theta_agg
  sim_list$SrvAge_pop_corr_pars <- SrvAge_pop_corr_pars
  sim_list$SrvAge_pop_corr_pars_agg <- SrvAge_pop_corr_pars_agg
  sim_list$SrvAgeComps_pop_Type <- SrvAgeComps_pop_Type

  sim_list$comp_srvlen_pop_like <- comp_srvlen_pop_like
  sim_list$ISS_SrvLenComps_pop <- ISS_SrvLenComps_pop
  sim_list$ln_SrvLen_pop_theta <- ln_SrvLen_pop_theta
  sim_list$ln_SrvLen_pop_theta_agg <- ln_SrvLen_pop_theta_agg
  sim_list$SrvLen_pop_corr_pars <- SrvLen_pop_corr_pars
  sim_list$SrvLen_pop_corr_pars_agg <- SrvLen_pop_corr_pars_agg
  sim_list$SrvLenComps_pop_Type <- SrvLenComps_pop_Type

  # whether a data source reports once a season or once a year, and the season an annual total lands in
  n_srv <- sim_list$n_srv_fleets
  sim_list$SrvIdx_seas_Type <- parse_seas_agg_spec(SrvIdx_seas_Type, "SrvIdx_seas_Type", n_srv)
  sim_list$SrvIdx_pop_seas_Type <- parse_seas_agg_spec(SrvIdx_pop_seas_Type, "SrvIdx_pop_seas_Type", n_srv)
  sim_list$SrvAgeComps_seas_Type <- parse_seas_agg_spec(SrvAgeComps_seas_Type, "SrvAgeComps_seas_Type", n_srv)
  sim_list$SrvLenComps_seas_Type <- parse_seas_agg_spec(SrvLenComps_seas_Type, "SrvLenComps_seas_Type", n_srv)
  sim_list$SrvAgeComps_pop_seas_Type <- parse_seas_agg_spec(SrvAgeComps_pop_seas_Type, "SrvAgeComps_pop_seas_Type", n_srv)
  sim_list$SrvLenComps_pop_seas_Type <- parse_seas_agg_spec(SrvLenComps_pop_seas_Type, "SrvLenComps_pop_seas_Type", n_srv)
  sim_list <- store_seas_agg_slot(sim_list, seas_agg_slot, n_srv)

  # an mvn index is one draw over a covariance the season layout defines, so it cannot be collapsed
  if(any(sim_list$SrvIdx_seas_Type == 1 & SrvIdx_LikeType == 2))
    stop("SrvIdx_seas_Type is 'aggSeas' for a survey whose SrvIdx_LikeType is 'mvn'. A multivariate ",
         "normal index is drawn once over a covariance built from the observed cells, so its seasons ",
         "cannot be summed after the draw. Use a lognormal or normal index for an annual total.")

  return(sim_list)

} # end function
