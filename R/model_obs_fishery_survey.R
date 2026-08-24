# Stage 2 of 3: objective function
#
# Predicted observations from the population state: catch, discards, and index
# and composition at age or length, for the fishery and the survey. Reads the
# within season averages from model_transition.R so predictions land at the right
# point in the season.

#' Fishery observation model
#'
#' Converts realized fishing mortality and population state into predicted
#' catch-at-age/length, discards, and fishery indices. Called once from the
#' "Fishery Observation Model" section of \code{SPoRC_rtmb.R}. All array
#' arguments matching an output name (e.g. \code{fish_q}, \code{CAA}) are
#' passed in already dimensioned (typically all-zero) and returned fully
#' populated over \code{1:n_yrs}.
#'
#' @param n_pop,n_regions,n_yrs,n_seas,n_fish_fleets,n_sexes Dimension sizes.
#' @param fish_q_blocks Array \code{[region, year, fish_fleet]} of time-block
#'   catchability indices.
#' @param ln_fish_q Array \code{[region, block, fish_fleet]} of log fishery
#'   catchability.
#' @param fish_q Array \code{[region, year, fish_fleet]}, output container.
#' @param ret_FAA,disc_FAA Arrays \code{[pop, region, year, season, age, sex,
#'   fish_fleet]} of retained/discarded fishing mortality at age.
#' @param ZAA Array \code{[pop, region, year, season, age, sex]} of total
#'   mortality at age.
#' @param NAA Array \code{[pop, region, year, season, age, sex]} of numbers at
#'   age.
#' @param CAA,DAA Arrays \code{[pop, region, year, season, age, sex,
#'   fish_fleet]}, output containers for retained/discarded catch at age.
#' @param CAL,DAL Arrays \code{[pop, region, year, season, len, sex,
#'   fish_fleet]}, output containers for retained/discarded catch at length.
#' @param PredCatch,PredDiscard,PredFishIdx Arrays \code{[pop, region, year,
#'   season, fish_fleet]}, output containers.
#' @param fit_lengths Integer (0/1) switch for computing length compositions.
#' @param SizeAgeTrans Array \code{[pop, region, year, season, len, age, sex]}
#'   age-to-length transition matrix.
#' @param catch_units,discard_units Integer vectors \code{[fish_fleet]}
#'   selecting abundance/biomass/fraction units.
#' @param WAA_fish Array \code{[pop, region, year, season, age, sex,
#'   fish_fleet]} of fishery weight at age.
#' @param dmr Array \code{[region, year, season, fish_fleet]} of discard
#'   mortality rate.
#' @param t_fish Array \code{[region, season, fish_fleet]} of fishery index
#'   timing, the fraction of the season elapsed when the index is observed.
#'   Numbers at age are decayed by \code{exp(-t_fish * ZAA)} before the index is
#'   formed, mirroring \code{t_srv} for surveys. \code{NULL} (the default) skips
#'   the decay entirely and reproduces a start-of-season index.
#' @param fish_idx_type Integer vector \code{[fish_fleet]} selecting
#'   abundance/biomass fishery index type.
#' @param fish_sel,ret_sel Arrays \code{[pop, region, year, season, age, sex,
#'   fish_fleet]} of total/retained fishery selectivity.
#' @param do_caal Integer (0/1) switch for building the joint catch at length
#'   and age arrays. Off by default. 
#' @param Fish_caal,Fish_caal_discard Arrays \code{[pop, region, year, season, len, age, sex,
#'   fish_fleet]}, output containers for retained/discarded catch at length and
#'   age. Only written when \code{do_caal == 1} and \code{fit_lengths == 1}.
#' @param SizeAgeTrans_fish Array \code{[pop, region, year, season, len, age,
#'   sex, fish_fleet]}, each fleet's own key from the growth module, or
#'   \code{NULL} to read the shared data key.
#' @param fish_len_comp_sel Integer vector \code{[fish_fleet]}. \code{0}
#'   (default) selects the catch at age and spreads it over lengths afterwards,
#'   so every fish of an age is equally catchable. \code{1} selects at length
#'   instead: the fish at each age are spread over the lengths of the
#'   composition key and then selected length by length,
#'   \eqn{C_l = s(l) \sum_a P(l \mid a)\, N_a (1 - e^{-Z_a}) F / Z_a}, so the
#'   long fish of an age are taken more often. Requires
#'   \code{fish_selex_type == 1}.
#' @param fish_selex_type,ret_selex_type Integer (0 age, 1 length), the scale
#'   total and retention selectivity are defined on.
#' @param fish_sel_l,ret_sel_l Arrays \code{[region, year, len, sex,
#'   fish_fleet]} of selectivity at length, read under
#'   \code{fish_len_comp_sel == 1}.
#' @param Fmort Array \code{[region, year, season, fish_fleet]} of fishing
#'   mortality, read under \code{fish_len_comp_sel == 1}.
#'
#' @return List with elements \code{fish_q}, \code{CAA}, \code{DAA},
#'   \code{CAL}, \code{DAL}, \code{Fish_caal}, \code{Fish_caal_discard},
#'   \code{PredCatch}, \code{PredDiscard}, \code{PredFishIdx}.
#'
#' @keywords internal
#' @import RTMB
get_fishery_observation_model <- function(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets, n_sexes,
                                           fish_q_blocks, ln_fish_q, fish_q,
                                           ret_FAA, disc_FAA, ZAA, NAA,
                                           CAA, DAA, CAL, DAL, PredCatch, PredDiscard, PredFishIdx,
                                           fit_lengths, SizeAgeTrans,
                                           catch_units, discard_units, WAA_fish, dmr,
                                           fish_idx_type, fish_sel, ret_sel,
                                           Mrate = NULL, move_timing = 0, seasdur = rep(1, n_seas),
                                           NAA_int = NULL, t_fish = NULL, fish_idx_ages = NULL,
                                           fish_q_type = NULL, do_fish_q_cov = 0, fish_q_cov = NULL,
                                           fish_q_coeff = NULL, ObsFishIdx = NULL, UseFishIdx = NULL,
                                           do_caal = 0, Fish_caal = NULL, Fish_caal_discard = NULL,
                                           SizeAgeTrans_fish = NULL,
                                           fish_len_comp_sel = NULL, fish_selex_type = 0, ret_selex_type = 0,
                                           fish_sel_l = NULL, ret_sel_l = NULL, Fmort = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(is.null(fish_q_type)) fish_q_type <- rep(0, n_fish_fleets)
  if(is.null(fish_len_comp_sel)) fish_len_comp_sel <- rep(0, n_fish_fleets)

  n_lens <- dim(CAL)[5] # length bins, read off a container that is always allocated

  # Every age contributes to the index unless a fleet restricts it.
  if(is.null(fish_idx_ages)) fish_idx_ages <- array(1, dim = c(dim(NAA)[5], n_fish_fleets))

  if(move_timing == 2 && is.null(NAA_int)) {
    n_ages <- dim(NAA)[5]
    NAA_int <- array(0, dim = dim(NAA))
    for(p in 1:n_pop) {
      for(y in 1:n_yrs) {
        for(seas in 1:n_seas) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              NAA_int[p,,y,seas,a,s] <- integrate_seas_abundance(NAA[p,,y,seas,a,s], ZAA[p,,y,seas,a,s],
                                                                 Mrate[p,,,y,seas,a,s], seasdur[seas])
            } # end s loop
          } # end a loop
        } # end seas loop
      } # end y loop
    } # end p loop
  }

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(f in 1:n_fish_fleets) {

          fish_q_blk_idx <- fish_q_blocks[r,y,f] # get time-block catchability index
          # Analytically solved fleets get their catchability filled in after the
          # unscaled index exists, so hold them at 1 for now.
          if(fish_q_type[f] == 0) fish_q[r,y,f] <- exp(ln_fish_q[r,fish_q_blk_idx,f]) # Input into fishery catchability container
          else fish_q[r,y,f] <- 1
          if(do_fish_q_cov == 1 && fish_q_type[f] == 0) fish_q[r,y,f] <- fish_q[r,y,f] * exp(sum(fish_q_cov[r,y,f,] * fish_q_coeff[r,f,])) # adding covariate effects

          for(seas in 1:n_seas) {

            if(move_timing == 2) {
              # Spatial Baranov: F times the season-integrated abundance
              CAA[p,r,y,seas,,,f] <- ret_FAA[p,r,y,seas,,,f] * NAA_int[p,r,y,seas,,]
              DAA[p,r,y,seas,,,f] <- disc_FAA[p,r,y,seas,,,f] * NAA_int[p,r,y,seas,,]
            } else {
              # Retained Catch at Age
              CAA[p,r,y,seas,,,f] <- ret_FAA[p,r,y,seas,,,f] / ZAA[p,r,y,seas,,] * NAA[p,r,y,seas,,] * (1 - exp(-ZAA[p,r,y,seas,,]))
              # Dead Discarded Catch at Age
              DAA[p,r,y,seas,,,f] <- disc_FAA[p,r,y,seas,,,f] / ZAA[p,r,y,seas,,] * NAA[p,r,y,seas,,] * (1 - exp(-ZAA[p,r,y,seas,,]))
            }

            if(fit_lengths == 1) {
              for(s in 1:n_sexes) {
                # the fleet's own key when the growth module built one, the shared data key otherwise
                key_f <- if(is.null(SizeAgeTrans_fish)) SizeAgeTrans[p,r,y,seas,,,s] else SizeAgeTrans_fish[p,r,y,seas,,,s,f]
                if(fish_len_comp_sel[f] == 0) {
                  CAL[p,r,y,seas,,s,f] <- key_f %*% CAA[p,r,y,seas,,s,f] # Retained Catch at length
                  DAL[p,r,y,seas,,s,f] <- key_f %*% DAA[p,r,y,seas,,s,f] # Discarded Catch at length
                  # Joint catch at length and age. SizeAgeTrans holds P(len | age), so
                  # scaling each age column by the catch at that age gives an array whose
                  # length margin is CAL and whose age margin is CAA. Conditioning on a
                  # length bin is left to the reader, since a composition likelihood
                  # normalizes within whatever bins it is handed.
                  if(do_caal == 1) {
                    Fish_caal[p,r,y,seas,,,s,f] <- key_f * rep(CAA[p,r,y,seas,,s,f], each = n_lens) # Retained catch at length and age
                    Fish_caal_discard[p,r,y,seas,,,s,f] <- key_f * rep(DAA[p,r,y,seas,,s,f], each = n_lens) # Discarded catch at length and age
                  } # building caal
                } else {
                  # Selectivity applied at length: the fish available over the season at
                  # each age are spread over length by the composition key, then taken
                  # length by length by the length selectivity (and retention, at length
                  # when it is defined there, at age otherwise)
                  avail <- if(move_timing == 2) NAA_int[p,r,y,seas,,s] else NAA[p,r,y,seas,,s] * (1 - exp(-ZAA[p,r,y,seas,,s])) / ZAA[p,r,y,seas,,s]
                  avail <- avail * Fmort[r,y,seas,f]
                  if(ret_selex_type == 1) {
                    ret_l <- ret_sel_l[r,y,,s,f]; ret_a <- rep(1, length(avail))
                  } else {
                    ret_l <- rep(1, n_lens); ret_a <- ret_sel[p,r,y,seas,,s,f]
                  }
                  sel_at_l <- fish_sel_l[r,y,,s,f]
                  joint_ret <- (key_f * rep(avail * ret_a, each = n_lens)) * (sel_at_l * ret_l) # [len, age]
                  joint_disc <- (key_f * rep(avail * (1 - ret_a), each = n_lens)) * (sel_at_l * (1 - ret_l)) * dmr[r,y,seas,f]
                  CAL[p,r,y,seas,,s,f] <- rowSums(joint_ret)
                  DAL[p,r,y,seas,,s,f] <- rowSums(joint_disc)
                  if(do_caal == 1) {
                    Fish_caal[p,r,y,seas,,,s,f] <- joint_ret
                    Fish_caal_discard[p,r,y,seas,,,s,f] <- joint_disc
                  } # building caal
                }
              } # end s loop
            } # fitting lengths

            # Get catch
            if(catch_units[f] == 0) PredCatch[p,r,y,seas,f] <- sum(CAA[p,r,y,seas,,,f]) # abundance
            if(catch_units[f] == 1) PredCatch[p,r,y,seas,f] <- sum(CAA[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f]) # biomass

            # Get discards
            if(discard_units[f] == 0) PredDiscard[p,r,y,seas,f] <- sum(DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f])  # total discard abundance
            if(discard_units[f] == 1) PredDiscard[p,r,y,seas,f] <- sum(DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f] * WAA_fish[p,r,y,seas,,,f])  # total discard biomass
            if(discard_units[f] == 2) {
              total_catch <- CAA[p,r,y,seas,,,f] + DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f]
              PredDiscard[p,r,y,seas,f] <- 1 - sum(CAA[p,r,y,seas,,,f]) / sum(total_catch)
            } # abundance fraction
            if(discard_units[f] == 3) {
              total_catch <- CAA[p,r,y,seas,,,f] + DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f]
              PredDiscard[p,r,y,seas,f] <- 1 - sum(CAA[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f]) / sum(total_catch * WAA_fish[p,r,y,seas,,,f])
            } # biomass fraction

            # Get fishery index. t_fish is the fraction of the season elapsed when
            # the index is observed, so the numbers are decayed by that much total
            # mortality first; t_fish = 0 reproduces a start-of-season index, which
            # is what the model did before this was available.
            FishIdxN <- if(move_timing == 2) NAA_int[p,r,y,seas,,] else NAA[p,r,y,seas,,]
            if(!is.null(t_fish)) FishIdxN <- FishIdxN * exp(-t_fish[r,seas,f] * ZAA[p,r,y,seas,,])
            # fish_idx_ages restricts which ages enter the index total without
            # touching the selectivity the compositions share.
            fidx_ages <- array(fish_idx_ages[,f], dim = c(dim(NAA)[5], n_sexes))
            if(fish_idx_type[f] == 0) PredFishIdx[p,r,y,seas,f] <- fish_q[r,y,f] * sum(FishIdxN * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f] * fidx_ages) # retained abundance
            if(fish_idx_type[f] == 1) PredFishIdx[p,r,y,seas,f] <- fish_q[r,y,f] * sum(FishIdxN * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f] * fidx_ages) # retained biomass
          } # end seas loop

        } # end f loop
      } # end y loop
    } # end r loop
  } # end p loop

  # A fleet whose catchability is a scaling nuisance can be concentrated out of
  # the likelihood rather than estimated. The solve uses only the years with
  # observations, and is done on the population-summed index because that is
  # what the index likelihood compares against.
  if(any(fish_q_type != 0) && !is.null(ObsFishIdx) && !is.null(UseFishIdx)) {
    for(f in 1:n_fish_fleets) {
      if(fish_q_type[f] == 0) next
      for(r in 1:n_regions) {

        use_rf <- array(UseFishIdx[r,,,f], dim = c(n_yrs, n_seas))
        if(!any(use_rf == 1)) next
        obs_idx <- which(use_rf == 1)
        obs_map <- arrayInd(obs_idx, dim(use_rf))

        obs_vec <- rep(0, length(obs_idx))
        pred_vec <- rep(0, length(obs_idx))
        for(i in seq_along(obs_idx)) {
          y <- obs_map[i,1]
          seas <- obs_map[i,2]
          obs_vec[i] <- ObsFishIdx[r,y,seas,f]
          pred_vec[i] <- sum(PredFishIdx[,r,y,seas,f])
        } # end i loop

        if(fish_q_type[f] == 1) q_hat <- mean(obs_vec) / mean(pred_vec) # arithmetic scaling
        if(fish_q_type[f] == 2) q_hat <- exp(mean(log(obs_vec) - log(pred_vec))) # log-scale scaling

        fish_q[r,,f] <- q_hat
        PredFishIdx[,r,,,f] <- PredFishIdx[,r,,,f] * q_hat

      } # end r loop
    } # end f loop
  }

  return(list(fish_q = fish_q, CAA = CAA, DAA = DAA, CAL = CAL, DAL = DAL,
              Fish_caal = Fish_caal, Fish_caal_discard = Fish_caal_discard,
              PredCatch = PredCatch, PredDiscard = PredDiscard, PredFishIdx = PredFishIdx))
}

#' Survey observation model
#'
#' Converts population state into predicted survey index-at-age/length and
#' survey indices. Called once from the "Survey Observation Model" section of
#' \code{SPoRC_rtmb.R}. All array arguments matching an output name (e.g.
#' \code{srv_q}, \code{SrvIAA}) are passed in already dimensioned and returned
#' fully populated (or, for \code{srv_sel}, further updated) over
#' \code{1:n_yrs}.
#'
#' @param n_pop,n_regions,n_yrs,n_seas,n_srv_fleets,n_sexes Dimension sizes.
#' @param srv_q_blocks Array \code{[region, year, srv_fleet]} of time-block
#'   catchability indices.
#' @param ln_srv_q Array \code{[region, block, srv_fleet]} of log survey
#'   catchability.
#' @param srv_q Array \code{[region, year, srv_fleet]}, output container.
#' @param do_srv_q_cov Integer (0/1) switch for a catchability covariate
#'   effect.
#' @param srv_q_cov Array \code{[region, year, srv_fleet, covariate]} of
#'   covariate values.
#' @param srv_q_coeff Array \code{[region, srv_fleet, covariate]} of
#'   covariate coefficients.
#' @param srv_selex_type Integer (0 = age-based, 1 = length-based) switch.
#' @param srv_sel Array \code{[pop, region, year, season, age, sex,
#'   srv_fleet]} of survey selectivity at age; when \code{srv_selex_type == 1}
#'   this is derived here from \code{srv_sel_l} via \code{SizeAgeTrans} for
#'   \code{1:n_yrs}.
#' @param SizeAgeTrans_srv Array \code{[pop, region, year, season, len, age,
#'   sex, srv_fleet]}, each survey's own key from the growth module, or
#'   \code{NULL} to read the shared data key.
#' @param srv_len_comp_sel Integer vector \code{[srv_fleet]}. \code{0}
#'   (default) expands the index at age through the key; \code{1} applies the
#'   length selectivity at length to the numbers spread over the composition
#'   own key, \eqn{I_l = s(l) \sum_a P(l \mid a)\, N_a e^{-t Z_a}}. Needs
#'   \code{srv_selex_type == 1}.
#' @param srv_sel_l Array \code{[region, year, len, sex, srv_fleet]} of survey
#'   selectivity at length.
#' @param SizeAgeTrans Array \code{[pop, region, year, season, len, age, sex]}
#'   age-to-length transition matrix.
#' @param NAA Array \code{[pop, region, year, season, age, sex]} of numbers at
#'   age.
#' @param ZAA Array \code{[pop, region, year, season, age, sex]} of total
#'   mortality at age.
#' @param t_srv Array \code{[region, season, srv_fleet]} of survey timing
#'   (fraction of season elapsed).
#' @param SrvIAA Array \code{[pop, region, year, season, age, sex, srv_fleet]},
#'   output container for survey index at age.
#' @param fit_lengths Integer (0/1) switch for computing length compositions.
#' @param SrvIAL Array \code{[pop, region, year, season, len, sex, srv_fleet]},
#'   output container for survey index at length.
#' @param RecDev_anom Array \code{[pop, region, deviation]} of recruitment
#'   deviations measured from the centre their penalty asserts, or \code{NULL}
#'   when no fleet observes them. Read only by fleets with
#'   \code{srv_idx_type == 2}.
#' @param srv_idx_type Integer vector \code{[srv_fleet]} selecting
#'   abundance/biomass survey index type.
#' @param WAA_srv Array \code{[pop, region, year, season, age, sex,
#'   srv_fleet]} of survey weight at age.
#' @param PredSrvIdx Array \code{[pop, region, year, season, srv_fleet]},
#'   output container.
#' @param srv_idx_ages Array \code{[age, srv_fleet]} of 0/1 weights selecting
#'   which ages contribute to each fleet's index total, or \code{NULL} for all
#'   ages. Restricting to a single age turns that fleet into an index of that
#'   age alone, and the compositions keep using the full age range because the
#'   restriction is applied to the index sum rather than to selectivity.
#' @param srv_q_type Integer vector \code{[srv_fleet]} selecting how
#'   catchability is obtained: \code{0} estimated as \code{exp(ln_srv_q)},
#'   \code{1} solved analytically as the ratio of mean observed to mean
#'   predicted, \code{2} solved analytically on the log scale as
#'   \code{exp(mean(log(obs) - log(pred)))}. \code{NULL} means all estimated.
#' @param ObsSrvIdx,UseSrvIdx Arrays \code{[region, year, season, srv_fleet]}
#'   of observed index values and their use flags. Required only when a fleet
#'   solves catchability analytically.
#'
#' @param do_caal Integer (0/1) switch for building the joint survey index at
#'   length and age array. Off by default. 
#' @param Srv_caal Array \code{[pop, region, year, season, len, age, sex,
#'   srv_fleet]}, output container for the survey index at length and age. Only
#'   written when \code{do_caal == 1} and \code{fit_lengths == 1}.
#'
#' @return List with elements \code{srv_q}, \code{srv_sel}, \code{SrvIAA},
#'   \code{SrvIAL}, \code{Srv_caal}, \code{PredSrvIdx}.
#'
#' @keywords internal
#' @import RTMB
get_survey_observation_model <- function(n_pop, n_regions, n_yrs, n_seas, n_srv_fleets, n_sexes,
                                          srv_q_blocks, ln_srv_q, srv_q, do_srv_q_cov, srv_q_cov, srv_q_coeff,
                                          srv_selex_type, srv_sel, srv_sel_l, SizeAgeTrans,
                                          NAA, ZAA, t_srv, SrvIAA, fit_lengths, SrvIAL,
                                          srv_idx_type, WAA_srv, PredSrvIdx,
                                          Mrate = NULL, move_timing = 0, seasdur = rep(1, n_seas),
                                          srv_idx_ages = NULL, srv_q_type = NULL,
                                          ObsSrvIdx = NULL, UseSrvIdx = NULL,
                                          RecDev_anom = NULL,
                                          do_caal = 0, Srv_caal = NULL,
                                          SizeAgeTrans_srv = NULL,
                                          srv_len_comp_sel = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_ages <- dim(NAA)[5]
  n_lens <- dim(SrvIAL)[5] # length bins, read off a container that is always allocated

  # Every age contributes to the index unless a fleet restricts it, and every
  # fleet estimates its own catchability unless it solves one analytically.
  if(is.null(srv_idx_ages)) srv_idx_ages <- array(1, dim = c(n_ages, n_srv_fleets))
  if(is.null(srv_q_type)) srv_q_type <- rep(0, n_srv_fleets)

  # Numbers at survey timing
  if(move_timing == 2) {
    n_ages_srv <- dim(NAA)[5]
    SrvN <- array(0, dim = c(dim(NAA), n_srv_fleets))
    for(p in 1:n_pop) {
      for(y in 1:n_yrs) {
        for(seas in 1:n_seas) {
          for(sf in 1:n_srv_fleets) {
            for(a in 1:n_ages_srv) {
              for(s in 1:n_sexes) {
                SrvN[p,,y,seas,a,s,sf] <- survey_state(NAA[p,,y,seas,a,s], NULL, ZAA[p,,y,seas,a,s],
                                                       Mrate[p,,,y,seas,a,s], seasdur[seas],
                                                       t_srv[,seas,sf], move_timing)
              } # end s loop
            } # end a loop
          } # end sf loop
        } # end seas loop
      } # end y loop
    } # end p loop
  }

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(sf in 1:n_srv_fleets) {

          srv_q_blk_idx <- srv_q_blocks[r,y,sf] # get time-block catchability index
          # Analytically solved fleets get their catchability filled in after the
          # unscaled index exists, so hold them at 1 for now.
          if(srv_q_type[sf] == 0) srv_q[r,y,sf] <- exp(ln_srv_q[r,srv_q_blk_idx,sf]) # Input into survey catchability container
          else srv_q[r,y,sf] <- 1
          if(do_srv_q_cov == 1 && srv_q_type[sf] == 0) srv_q[r,y,sf] <- srv_q[r,y,sf] * exp(sum(srv_q_cov[r,y,sf,] * srv_q_coeff[r,sf,])) # adding covariate effects

          for(seas in 1:n_seas) {

            # Convert length-selex to age-selex
            # the survey's own key when the growth module built one, the shared data key otherwise
            key_sf <- function(s) if(is.null(SizeAgeTrans_srv)) SizeAgeTrans[p,r,y,seas,,,s] else SizeAgeTrans_srv[p,r,y,seas,,,s,sf]
            if(srv_selex_type == 1) for(s in 1:n_sexes) srv_sel[p,r,y,seas,,s,sf] <- srv_sel_l[r,y,,s,sf] %*% key_sf(s)

            if(move_timing == 2) {
              # Selectivity applies at the destination region, after propagation
              SrvIAA[p,r,y,seas,,,sf] <- SrvN[p,r,y,seas,,,sf] * srv_sel[p,r,y,seas,,,sf] # Survey index at age
            } else {
              SrvIAA[p,r,y,seas,,,sf] <- NAA[p,r,y,seas,,] * srv_sel[p,r,y,seas,,,sf] * exp(-t_srv[r,seas,sf] * ZAA[p,r,y,seas,,]) # Survey index at age
            }

            if(fit_lengths == 1) {
              for(s in 1:n_sexes) {
                if(is.null(srv_len_comp_sel) || srv_len_comp_sel[sf] == 0) {
                  SrvIAL[p,r,y,seas,,s,sf] <- key_sf(s) %*% SrvIAA[p,r,y,seas,,s,sf] # Survey index at length
                  # Joint survey index at length and age, built the same way as the
                  # fishery arrays: each age column of P(len | age) scaled by the index
                  # at that age, so the length margin is SrvIAL and the age margin SrvIAA.
                  if(do_caal == 1) Srv_caal[p,r,y,seas,,,s,sf] <- key_sf(s) * rep(SrvIAA[p,r,y,seas,,s,sf], each = n_lens) # Survey index at length and age
                } else {
                  # Selectivity applied at length: the numbers present at the survey's
                  # timing spread over length by the composition key, then selected
                  # length by length
                  present <- if(move_timing == 2) SrvN[p,r,y,seas,,s,sf] else NAA[p,r,y,seas,,s] * exp(-t_srv[r,seas,sf] * ZAA[p,r,y,seas,,s])
                  joint <- (key_sf(s) * rep(present, each = n_lens)) * srv_sel_l[r,y,,s,sf] # [len, age]
                  SrvIAL[p,r,y,seas,,s,sf] <- rowSums(joint)
                  if(do_caal == 1) Srv_caal[p,r,y,seas,,,s,sf] <- joint
                }
              } # end s loop
            } # fitting lengths

            # srv_idx_ages restricts which ages enter the index total without
            # touching the selectivity the compositions share. An index of a
            # single age is the limiting case.
            idx_ages <- array(srv_idx_ages[,sf], dim = c(n_ages, n_sexes))
            if(srv_idx_type[sf] == 0) PredSrvIdx[p,r,y,seas,sf] <- srv_q[r,y,sf] * sum(SrvIAA[p,r,y,seas,,,sf] * idx_ages) # abundance
            if(srv_idx_type[sf] == 1) PredSrvIdx[p,r,y,seas,sf] <- srv_q[r,y,sf] * sum(SrvIAA[p,r,y,seas,,,sf] * WAA_srv[p,r,y,seas,,,sf] * idx_ages) # biomass

            # Recruitment index as a survey fleet 
            if(srv_idx_type[sf] == 2 && !is.null(RecDev_anom) && y <= dim(RecDev_anom)[3]) {
              PredSrvIdx[p,r,y,seas,sf] <- srv_q[r,y,sf] * RecDev_anom[p,r,y]
            }

          } # end seas loop

        } # end sf loop
      } # end y loop
    } # end r loop
  } # end p loop

  ### Analytically solved catchability ---------------------------------------
  # A fleet's catchability can be a scaling nuisance and can be concentrated out of the likelihood
  if(any(srv_q_type != 0) && !is.null(ObsSrvIdx) && !is.null(UseSrvIdx)) {
    for(sf in 1:n_srv_fleets) {
      if(srv_q_type[sf] == 0) next
      for(r in 1:n_regions) {

        use_rf <- array(UseSrvIdx[r,,,sf], dim = c(n_yrs, n_seas))
        if(!any(use_rf == 1)) next
        obs_idx <- which(use_rf == 1)
        obs_map <- arrayInd(obs_idx, dim(use_rf))

        obs_vec <- rep(0, length(obs_idx))
        pred_vec <- rep(0, length(obs_idx))
        for(i in seq_along(obs_idx)) {
          y <- obs_map[i,1]
          seas <- obs_map[i,2]
          obs_vec[i] <- ObsSrvIdx[r,y,seas,sf]
          pred_vec[i] <- sum(PredSrvIdx[,r,y,seas,sf])
        } # end i loop

        if(srv_q_type[sf] == 1) q_hat <- mean(obs_vec) / mean(pred_vec) # arithmetic scaling
        if(srv_q_type[sf] == 2) q_hat <- exp(mean(log(obs_vec) - log(pred_vec))) # log-scale scaling

        srv_q[r,,sf] <- q_hat
        PredSrvIdx[,r,,,sf] <- PredSrvIdx[,r,,,sf] * q_hat

      } # end r loop
    } # end sf loop
  }

  return(list(srv_q = srv_q, srv_sel = srv_sel, SrvIAA = SrvIAA, SrvIAL = SrvIAL,
              Srv_caal = Srv_caal, PredSrvIdx = PredSrvIdx))
}
