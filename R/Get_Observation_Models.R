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
#' @param fish_idx_type Integer vector \code{[fish_fleet]} selecting
#'   abundance/biomass fishery index type.
#' @param fish_sel,ret_sel Arrays \code{[pop, region, year, season, age, sex,
#'   fish_fleet]} of total/retained fishery selectivity.
#'
#' @return List with elements \code{fish_q}, \code{CAA}, \code{DAA},
#'   \code{CAL}, \code{DAL}, \code{PredCatch}, \code{PredDiscard},
#'   \code{PredFishIdx}.
#'
#' @keywords internal
#' @import RTMB
get_fishery_observation_model <- function(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets, n_sexes,
                                           fish_q_blocks, ln_fish_q, fish_q,
                                           ret_FAA, disc_FAA, ZAA, NAA,
                                           CAA, DAA, CAL, DAL, PredCatch, PredDiscard, PredFishIdx,
                                           fit_lengths, SizeAgeTrans,
                                           catch_units, discard_units, WAA_fish, dmr,
                                           fish_idx_type, fish_sel, ret_sel) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(f in 1:n_fish_fleets) {

          fish_q_blk_idx <- fish_q_blocks[r,y,f] # get time-block catchability index
          fish_q[r,y,f] <- exp(ln_fish_q[r,fish_q_blk_idx,f]) # Input into fishery catchability container

          for(seas in 1:n_seas) {

            # Retained Catch at Age
            CAA[p,r,y,seas,,,f] <- ret_FAA[p,r,y,seas,,,f] / ZAA[p,r,y,seas,,] * NAA[p,r,y,seas,,] * (1 - exp(-ZAA[p,r,y,seas,,]))

            # Dead Discarded Catch at Age
            DAA[p,r,y,seas,,,f] <- disc_FAA[p,r,y,seas,,,f] / ZAA[p,r,y,seas,,] * NAA[p,r,y,seas,,] * (1 - exp(-ZAA[p,r,y,seas,,]))

            if(fit_lengths == 1) {
              for(s in 1:n_sexes) {
                CAL[p,r,y,seas,,s,f] <- SizeAgeTrans[p,r,y,seas,,,s] %*% CAA[p,r,y,seas,,s,f] # Retained Catch at length
                DAL[p,r,y,seas,,s,f] <- SizeAgeTrans[p,r,y,seas,,,s] %*% DAA[p,r,y,seas,,s,f] # Discarded Catch at length
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

            # Get fishery index
            if(fish_idx_type[f] == 0) PredFishIdx[p,r,y,seas,f] <- fish_q[r,y,f] * sum(NAA[p,r,y,seas,,] * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f]) # retained abundance
            if(fish_idx_type[f] == 1) PredFishIdx[p,r,y,seas,f] <- fish_q[r,y,f] * sum(NAA[p,r,y,seas,,] * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f]) # retained biomass
          } # end seas loop

        } # end f loop
      } # end y loop
    } # end r loop
  } # end p loop

  return(list(fish_q = fish_q, CAA = CAA, DAA = DAA, CAL = CAL, DAL = DAL,
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
#' @param srv_idx_type Integer vector \code{[srv_fleet]} selecting
#'   abundance/biomass survey index type.
#' @param WAA_srv Array \code{[pop, region, year, season, age, sex,
#'   srv_fleet]} of survey weight at age.
#' @param PredSrvIdx Array \code{[pop, region, year, season, srv_fleet]},
#'   output container.
#'
#' @return List with elements \code{srv_q}, \code{srv_sel}, \code{SrvIAA},
#'   \code{SrvIAL}, \code{PredSrvIdx}.
#'
#' @keywords internal
#' @import RTMB
get_survey_observation_model <- function(n_pop, n_regions, n_yrs, n_seas, n_srv_fleets, n_sexes,
                                          srv_q_blocks, ln_srv_q, srv_q, do_srv_q_cov, srv_q_cov, srv_q_coeff,
                                          srv_selex_type, srv_sel, srv_sel_l, SizeAgeTrans,
                                          NAA, ZAA, t_srv, SrvIAA, fit_lengths, SrvIAL,
                                          srv_idx_type, WAA_srv, PredSrvIdx) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(sf in 1:n_srv_fleets) {

          srv_q_blk_idx <- srv_q_blocks[r,y,sf] # get time-block catchability index
          srv_q[r,y,sf] <- exp(ln_srv_q[r,srv_q_blk_idx,sf]) # Input into survey catchability container
          if(do_srv_q_cov == 1) srv_q[r,y,sf] <- srv_q[r,y,sf] * exp(sum(srv_q_cov[r,y,sf,] * srv_q_coeff[r,sf,])) # adding covariate effects

          for(seas in 1:n_seas) {

            # Convert length-selex to age-selex
            if(srv_selex_type == 1) for(s in 1:n_sexes) srv_sel[p,r,y,seas,,s,sf] <- srv_sel_l[r,y,,s,sf] %*% SizeAgeTrans[p,r,y,seas,,,s]

            SrvIAA[p,r,y,seas,,,sf] <- NAA[p,r,y,seas,,] * srv_sel[p,r,y,seas,,,sf] * exp(-t_srv[r,seas,sf] * ZAA[p,r,y,seas,,]) # Survey index at age

            if(fit_lengths == 1) {
              for(s in 1:n_sexes) {
                SrvIAL[p,r,y,seas,,s,sf] <- SizeAgeTrans[p,r,y,seas,,,s] %*% SrvIAA[p,r,y,seas,,s,sf] # Survey index at length
              } # end s loop
            } # fitting lengths

            if(srv_idx_type[sf] == 0) PredSrvIdx[p,r,y,seas,sf] <- srv_q[r,y,sf] * sum(SrvIAA[p,r,y,seas,,,sf]) # abundance
            if(srv_idx_type[sf] == 1) PredSrvIdx[p,r,y,seas,sf] <- srv_q[r,y,sf] * sum(SrvIAA[p,r,y,seas,,,sf] * WAA_srv[p,r,y,seas,,,sf]) # biomass

          } # end seas loop

        } # end sf loop
      } # end y loop
    } # end r loop
  } # end p loop

  return(list(srv_q = srv_q, srv_sel = srv_sel, SrvIAA = SrvIAA, SrvIAL = SrvIAL, PredSrvIdx = PredSrvIdx))
}
