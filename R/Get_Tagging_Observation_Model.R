#' Tagging observation model
#'
#' Forward-simulates conventional tag cohort availability (release, movement,
#' mortality/ageing) and predicted recaptures across all recapture years,
#' seasons, and cohorts. Called once from the "Tagging Observation Model"
#' section of \code{SPoRC_rtmb.R}, guarded by \code{any(use_conv_fish_tagging
#' == 1)}. Depends on population state (\code{fish_sel}, \code{ret_sel},
#' \code{srv_sel}, \code{NAA_bef}, \code{Movement}, \code{natmort},
#' \code{Fmort}, \code{dmr}) already computed upstream by the population
#' projection and selectivity sections.
#'
#' \code{conv_tag_fish_reporting}, \code{conv_tag_fish_avail}, and
#' \code{pred_conv_tag_fish_recap} are passed in already dimensioned
#' (typically all-zero) and returned fully populated.
#'
#' @param n_fish_fleets,n_regions,n_conv_tag_cohorts,n_yrs,n_seas,n_pop,n_ages,n_sexes
#'   Dimension sizes.
#' @param conv_tag_fish_reporting_blocks Array \code{[region, block, fish_fleet]}
#'   mapping years to reporting-rate blocks.
#' @param conv_tag_fish_reporting_pars Array \code{[region, block, fish_fleet]}
#'   of tag reporting rate parameters on the logit scale.
#' @param conv_tag_fish_reporting Array \code{[region, year, fish_fleet]},
#'   output container for reporting rate on the natural scale.
#' @param conv_tag_release_indicator Matrix \code{[cohort, 3]} of
#'   release region/year/season per tag cohort.
#' @param conv_tag_max_liberty Integer maximum years at liberty tracked.
#' @param use_conv_fish_tagging Integer vector \code{[fish_fleet]} (0/1)
#'   flagging fleets with tagging data.
#' @param Fmort Array \code{[region, year, season, fish_fleet]} of fishing
#'   mortality.
#' @param fish_sel,ret_sel Arrays \code{[pop, region, year, season, age, sex,
#'   fish_fleet]} of total/retained fishery selectivity.
#' @param dmr Array \code{[region, year, season, fish_fleet]} of discard
#'   mortality rate.
#' @param natmort Array \code{[pop, region, year, age, sex]} of natural
#'   mortality at age.
#' @param seasdur Numeric vector \code{[season]} of season duration (fraction
#'   of year).
#' @param ln_conv_tag_shed Numeric vector \code{[cohort]} of log tag shedding
#'   rate.
#' @param conv_tag_t_tagging Numeric vector \code{[cohort]} of within-season
#'   timing of tag release (1 = start of season).
#' @param conv_tagged_fish Array \code{[cohort, age, sex, ...]} of raw tag
#'   release numbers, passed through to \code{release_conv_tag_attr}.
#' @param conv_fish_tag_attr,conv_tag_release_platform Arguments passed
#'   through to \code{release_conv_tag_attr} controlling how releases are
#'   apportioned across dimensions.
#' @param srv_sel Array \code{[pop, region, year, season, age, sex,
#'   srv_fleet]} of survey selectivity, used by \code{release_conv_tag_attr}
#'   when releases are platform-apportioned via survey gear.
#' @param NAA_bef Array \code{[pop, region, year, season, age, sex]} of
#'   numbers at age before movement, used by \code{release_conv_tag_attr} to
#'   apportion releases.
#' @param ln_init_conv_tag_mort Numeric vector \code{[cohort]} of log initial
#'   tag-induced mortality rate.
#' @param do_recruits_move Integer (0/1) switch for whether age-1 fish are
#'   subject to movement.
#' @param Movement Array \code{[pop, region_from, region_to, year, season,
#'   age, sex]} of movement rates.
#' @param conv_tag_fish_avail Array \code{[liberty+1, season, cohort, pop,
#'   region, age, sex]}, output container for tags available for recapture.
#' @param pred_conv_tag_fish_recap Array \code{[liberty, season, cohort, pop,
#'   region, age, sex, fish_fleet]}, output container for predicted
#'   recaptures.
#'
#' @return List with elements \code{conv_tag_fish_reporting},
#'   \code{conv_tag_fish_avail}, \code{pred_conv_tag_fish_recap}.
#'
#' @keywords internal
#' @import RTMB
get_tagging_observation_model <- function(n_fish_fleets, n_regions, n_conv_tag_cohorts, n_yrs, n_seas, n_pop, n_ages, n_sexes,
                                           conv_tag_fish_reporting_blocks, conv_tag_fish_reporting_pars, conv_tag_fish_reporting,
                                           conv_tag_release_indicator, conv_tag_max_liberty, use_conv_fish_tagging,
                                           Fmort, fish_sel, ret_sel, dmr, natmort, seasdur,
                                           ln_conv_tag_shed, conv_tag_t_tagging, conv_tagged_fish,
                                           conv_fish_tag_attr, conv_tag_release_platform, srv_sel, NAA_bef,
                                           ln_init_conv_tag_mort, do_recruits_move, Movement,
                                           conv_tag_fish_avail, pred_conv_tag_fish_recap) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Set up tag reporting rates
  for(f in 1:n_fish_fleets) {
    for(r in 1:n_regions) {
      conv_tagrep_blk_idx <- conv_tag_fish_reporting_blocks[r,,f]  # Get all blocks for this region
      conv_tag_fish_reporting[r,,f] <- RTMB::plogis(conv_tag_fish_reporting_pars[r,conv_tagrep_blk_idx,f])  # inverse logit transform
    } # end r loop
  } # end f loop

  for(tc in 1:n_conv_tag_cohorts) {

    tr <- conv_tag_release_indicator[tc,1] # extract tag release region
    ty <- conv_tag_release_indicator[tc,2] # extract tag release year
    tseas <- conv_tag_release_indicator[tc,3] # extract tag release season

    for(ry in 1:min(conv_tag_max_liberty, n_yrs - ty + 1)) {   # years
      y <- ty + ry - 1 # get real year
      for(rseas in 1:n_seas) { # seasons

        # get fishing mortality
        tmp_FAA <- array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
        tmp_ret_FAA <- array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
        tmp_disc_DAA <- array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
        for(p in 1:n_pop) for(f in 1:n_fish_fleets) {
          if(use_conv_fish_tagging[f] == 1) {
            tmp_ret_FAA[p,,1,,,f] <- Fmort[, y, rseas, f] * fish_sel[p,,y,rseas,,,f] * ret_sel[p,,y,rseas,,,f]  # Retained fishing mortality
            tmp_disc_DAA[p,,1,,,f] <- Fmort[, y, rseas, f] * fish_sel[p,,y,rseas,,,f] * (1 - ret_sel[p,,y,rseas,,,f]) * dmr[,y,rseas,f] # Dead discard fishing mortality
            tmp_FAA[p,,1,,,f] <- tmp_ret_FAA[p,,1,,,f] + tmp_disc_DAA[p,,1,,,f] # Total fishing mortality
          } # end if
        } # end p loop

        # get total mortality
        tmp_natmort <- array(natmort[,,y,,], dim = c(n_pop, n_regions, 1, n_ages, n_sexes))
        tmp_ZAA <- (tmp_natmort * seasdur[rseas]) + apply(tmp_FAA, 1:5, sum) + (exp(ln_conv_tag_shed[tc]) * seasdur[rseas])

        # Discount with tagging time (conv_tag_t_tagging) if it doesn't happen at the start of the season / year
        if(ry == 1 && rseas == tseas) {

          if(conv_tag_t_tagging[tc] != 1) tmp_ZAA <- tmp_ZAA * conv_tag_t_tagging[tc]

          # apportion tagged fish out to appropriate dimensions if necessary
          tmp_tagged_fish <- release_conv_tag_attr(array(conv_tagged_fish[tc, , , ], dim = c(n_pop, n_ages, n_sexes)),
                                                  conv_fish_tag_attr[tc],
                                                  conv_tag_release_platform[tc,],
                                                  srv_sel, fish_sel, NAA_bef,
                                                  ty, tseas, tr, n_pop,
                                                  n_ages, n_sexes)

          # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
          conv_tag_fish_avail[1, rseas, tc, , tr, , ] <- array(tmp_tagged_fish * exp(-exp(ln_init_conv_tag_mort[tc])), dim = c(n_pop, n_ages, n_sexes))
        }

        # get temporary survival value
        tmp_SAA <- exp(-tmp_ZAA)

        # Move tagged fish around (skip only in first release year + tagging season when tagging occurs mid-season)
        if(conv_tag_t_tagging[tc] == 1 || ry != 1 || rseas != tseas) {
          for(p in 1:n_pop) {
            # Movement of tag cohorts
            if(do_recruits_move == 0) {
              for(a in 2:n_ages) for(s in 1:n_sexes) {
                conv_tag_fish_avail[ry, rseas, tc, p, , a, s] <-
                  t(conv_tag_fish_avail[ry, rseas, tc, p, , a, s]) %*%
                  Movement[p, , , y, rseas, a, s]
              }
            } else { # if recruits move
              for(a in 1:n_ages) for(s in 1:n_sexes) {
                conv_tag_fish_avail[ry, rseas, tc, p, , a, s] <-
                  t(conv_tag_fish_avail[ry, rseas, tc, p, , a, s]) %*%
                  Movement[p, , , y, rseas, a, s]
              } # end s loop
            } # end else
          } # end p loop
        } # end if

        # Apply mortality and ageing to tagged fish
        if(rseas < n_seas) {

          # Season mortality within a given year, advance to next season same year/age
          conv_tag_fish_avail[ry, rseas + 1, tc, , , , ] <-
            conv_tag_fish_avail[ry, rseas, tc, , , , ] *
            tmp_SAA[,,1,,]

        } else {

          # End of year mortality and age advancement (end of season)
          conv_tag_fish_avail[ry + 1, 1, tc, , , 2:n_ages, ] <-
            conv_tag_fish_avail[ry, n_seas, tc, , , 1:(n_ages-1), ] *
            tmp_SAA[,,1,1:(n_ages - 1),]

          # Accumulate plus group
          conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, ] <-
            conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, ] +
            conv_tag_fish_avail[ry, n_seas, tc, , , n_ages, ] *
            tmp_SAA[,,1,n_ages,]
        }

        # # Apply Baranov's to get predicted recaptures
        for(f in 1:n_fish_fleets) {
          for(p in 1:n_pop) {
            pred_conv_tag_fish_recap[ry,rseas,tc,p,,,,f] <- conv_tag_fish_reporting[,y,f] *
              (tmp_ret_FAA[p,,1,,,f] / tmp_ZAA[p,,1,,]) *
              conv_tag_fish_avail[ry,rseas,tc,p,,,] *
              (1 - tmp_SAA[p,,1,,])
          } # end p loop
        } # end f loop


      }
    }
  }

  return(list(conv_tag_fish_reporting = conv_tag_fish_reporting,
              conv_tag_fish_avail = conv_tag_fish_avail,
              pred_conv_tag_fish_recap = pred_conv_tag_fish_recap))
}
