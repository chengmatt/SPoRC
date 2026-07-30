# Stage 2 of 3: objective function
#
# Predicted conventional tag recaptures. Follows each release cohort forward
# through movement, mortality, ageing, shedding and reporting.
# release_conv_tag_attr spreads releases that were recorded at coarser
# resolution than the population arrays out to full dimensions.

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
                                           conv_tag_fish_avail, pred_conv_tag_fish_recap,
                                           Mrate = NULL, move_timing = 0) {

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

        # Fraction of this season the tag cohort is actually at liberty for. Mid-season releases
        # only experience the remainder of the season.
        tag_frac <- if(ry == 1 && rseas == tseas) conv_tag_t_tagging[tc] else 1
        tag_dur <- seasdur[rseas] * tag_frac

        # Discount with tagging time (conv_tag_t_tagging) if it doesn't happen at the start of the
        # season / year. every mortality component is scaled, not just the total: Baranov's F/Z is
        # the fraction of deaths owing to fishing, so scaling Z while leaving F at full-season
        # scale turns it into F/(Z * t_tag), which exceeds 1 whenever t_tag < F/Z -- i.e. predicts
        # more recaptures than there are dead tags. Scaling F and Z together leaves the ratio at
        # F/Z and lets the (1 - exp(-Z * t_tag)) term carry the shorter exposure, which is the
        # standard partial-interval Baranov.
        if(tag_frac != 1) {
          tmp_ZAA      <- tmp_ZAA      * tag_frac
          tmp_FAA      <- tmp_FAA      * tag_frac
          tmp_ret_FAA  <- tmp_ret_FAA  * tag_frac
          tmp_disc_DAA <- tmp_disc_DAA * tag_frac
        }

        if(ry == 1 && rseas == tseas) {

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

        # Whether the discrete movement step applies this season. A mid-season release skips it,
        # since a full-season transition matrix cannot represent a partial interval.
        tag_moves <- (conv_tag_t_tagging[tc] == 1 || ry != 1 || rseas != tseas)

        # Continuous movement needs no such exemption, and applying one would be inconsistent:
        # the generator is scaled by tag_dur below, so a mid-season release diffuses for exactly
        # the fraction of the season it was at liberty for -- the same fraction its mortality is
        # already scaled by. Freezing movement while still discounting mortality partially would
        # have tags dying on a partial season but holding station for a whole one.
        tag_moves_seas <- if(move_timing == 2) TRUE else tag_moves

        # Move tagged fish around. Under move_timing == 0 movement is applied here, in place,
        # and mortality follows below; under timings 1 and 2 both are carried together by the
        # seasonal transition operator in the mortality step.
        if(move_timing == 0 && tag_moves) {
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

        # Post-season tag numbers, before the ageing shift. Under move_timing == 0 movement
        # was applied above so this is plain survival; under timings 1 and 2 the transition
        # operator carries movement and mortality together over tag_dur.
        if(move_timing == 0 || n_regions == 1) {
          tag_step <- array(conv_tag_fish_avail[ry, rseas, tc, , , , ] * tmp_SAA[,,1,,],
                            dim = c(n_pop, n_regions, n_ages, n_sexes))
        } else {
          tag_step <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
          for(p in 1:n_pop) {
            for(a in 1:n_ages) {
              moves <- tag_moves_seas && (do_recruits_move == 1 || a > 1)
              for(s in 1:n_sexes) {
                Mv <- if(moves) Movement[p,,,y,rseas,a,s] else diag(n_regions)
                Qv <- if(moves) Mrate[p,,,y,rseas,a,s] else matrix(0, n_regions, n_regions)
                tag_step[p,,a,s] <- advance_seas(conv_tag_fish_avail[ry,rseas,tc,p,,a,s], Mv,
                                                 tmp_ZAA[p,,1,a,s], Qv, tag_dur, move_timing)
              } # end s loop
            } # end a loop
          } # end p loop
        }

        # Apply mortality and ageing to tagged fish
        if(rseas < n_seas) {

          # Season mortality within a given year, advance to next season same year/age
          conv_tag_fish_avail[ry, rseas + 1, tc, , , , ] <- tag_step

        } else {

          # End of year mortality and age advancement (end of season)
          conv_tag_fish_avail[ry + 1, 1, tc, , , 2:n_ages, ] <- tag_step[,,1:(n_ages - 1),]

          # Accumulate plus group
          conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, ] <-
            conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, ] + tag_step[,,n_ages,]
        }

        # # Apply Baranov's to get predicted recaptures
        for(f in 1:n_fish_fleets) {
          for(p in 1:n_pop) {
            if(move_timing == 2) {
              # Spatial Baranov: tags redistribute among regions while being caught, so
              # recaptures use the season-integrated tag abundance rather than the
              # region-local (1 - exp(-Z)) / Z form.
              tag_int <- array(0, dim = c(n_regions, n_ages, n_sexes))
              for(a in 1:n_ages) {
                # must match the generator used for tag_step above, or the cohort's
                # dynamics and its recaptures would be built on different movement
                moves <- tag_moves_seas && (do_recruits_move == 1 || a > 1)
                for(s in 1:n_sexes) {
                  Qv <- if(moves) Mrate[p,,,y,rseas,a,s] else matrix(0, n_regions, n_regions)
                  tag_int[,a,s] <- integrate_seas_abundance(conv_tag_fish_avail[ry,rseas,tc,p,,a,s],
                                                            tmp_ZAA[p,,1,a,s], Qv, tag_dur)
                } # end s loop
              } # end a loop
              # array() guards against R dropping a length-1 sex (or age) dimension from the
              # F slice, which would leave it non-conformable with the 3-d tag_int
              tmp_ret_FAA_slice <- array(tmp_ret_FAA[p,,1,,,f], dim = c(n_regions, n_ages, n_sexes))
              pred_conv_tag_fish_recap[ry,rseas,tc,p,,,,f] <- conv_tag_fish_reporting[,y,f] *
                tmp_ret_FAA_slice * tag_int
            } else {
              pred_conv_tag_fish_recap[ry,rseas,tc,p,,,,f] <- conv_tag_fish_reporting[,y,f] *
                (tmp_ret_FAA[p,,1,,,f] / tmp_ZAA[p,,1,,]) *
                conv_tag_fish_avail[ry,rseas,tc,p,,,] *
                (1 - tmp_SAA[p,,1,,])
            }
          } # end p loop
        } # end f loop


      }
    }
  }

  return(list(conv_tag_fish_reporting = conv_tag_fish_reporting,
              conv_tag_fish_avail = conv_tag_fish_avail,
              pred_conv_tag_fish_recap = pred_conv_tag_fish_recap))
}
#' Distribute Tagged Fish Releases to Full Population Dimensions
#'
#' When tag release data are not recorded at full population resolution —
#' i.e. when one or more of the population, age, or sex dimensions are
#' unattended in \code{tag_attr} — this function distributes the known tag
#' totals to full \code{[n_pop, n_ages, n_sexes]} resolution using
#' apportionment weights derived from the release platform (population
#' abundance, fishery catch-at-age, or survey index-at-age). If all three
#' dimensions are attended (\code{tag_attr = "p_a_s"}), \code{tagged_fish}
#' is returned unchanged with no computation performed.
#'
#' Apportionment weights are constructed from numbers-at-age
#' (\code{platform = "population"}), numbers-at-age multiplied by fishery
#' selectivity (\code{platform = "fishery"}), or numbers-at-age multiplied
#' by survey selectivity (\code{platform = "survey"}), all evaluated at the
#' release region, year, and season. Weights are then normalised
#' conditionally on the attended dimensions: the denominator for cell
#' \code{[p, a, s]} is the sum of raw weights across all cells that share
#' the same indices in the attended dimensions. This ensures that the
#' marginal totals of \code{tagged_fish} are preserved exactly along every
#' attended dimension. For example, if only age is attended
#' (\code{tag_attr = "a"}), age-specific totals in \code{tagged_fish} are
#' preserved while tags are distributed across population and sex in
#' proportion to the platform weights.
#'
#' @param tagged_fish Numeric vector or array of released tagged fish for a
#'   single tag cohort. Unattended dimensions are expected to be collapsed
#'   to index 1. Reshaped internally to \code{[n_pop, n_ages, n_sexes]}.
#' @param tag_attr Character string specifying which population dimensions
#'   are attended in \code{tagged_fish}. Constructed from any combination of
#'   \code{"p"} (population), \code{"a"} (age), and \code{"s"} (sex),
#'   joined by underscores (e.g. \code{"p_a_s"}, \code{"a"}, \code{"p_a"}).
#'   See \code{\link{Setup_Mod_Tagging}} for the full set of valid strings.
#' @param tag_release_platform Character vector of length 2. Element 1 is
#'   the release platform: one of \code{"population"}, \code{"fishery"}, or
#'   \code{"survey"}. Element 2 is the fleet index as a character string
#'   (coerced to integer internally), or \code{NA} when
#'   \code{platform = "population"}.
#' @param srv_sel Numeric array of survey selectivity
#'   \code{[n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets]}. Used as the
#'   age-sex apportionment weight when \code{platform = "survey"}.
#' @param fish_sel Numeric array of fishery selectivity
#'   \code{[n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets]}. Used as the
#'   age-sex apportionment weight when \code{platform = "fishery"}.
#' @param NAA Numeric array of numbers-at-age prior to movement
#'   \code{[n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes]}. Used
#'   directly as weights when \code{platform = "population"}, and multiplied
#'   by selectivity for fishery and survey platforms.
#' @param ty Integer. Model year index of the tag release cohort.
#' @param tseas Integer. Season index of the tag release cohort.
#' @param tr Integer. Region index of the tag release cohort.
#' @param n_pop Integer. Number of populations.
#' @param n_ages Integer. Number of age classes.
#' @param n_sexes Integer. Number of sexes (1 or 2).
#'
#'
#' @keywords internal
release_conv_tag_attr <- function(tagged_fish,
                                  tag_attr,
                                  tag_release_platform,
                                  srv_sel,
                                  fish_sel,
                                  NAA,
                                  ty,
                                  tseas,
                                  tr,
                                  n_pop,
                                  n_ages,
                                  n_sexes
) {

  "c" = RTMB::ADoverload("c")
  "[<-" = RTMB::ADoverload("[<-")

  # Keep the original collapsed input for safe indexing into attended dims.
  tagged_fish_orig <- tagged_fish

  # get tagging attributes
  attr_parts = strsplit(tag_attr, "_")[[1]]
  attended_p = "p" %in% attr_parts
  attended_a = "a" %in% attr_parts
  attended_s = "s" %in% attr_parts

  # return tagged_fish (reshaped to full dims) if p_a_s — all dims attended,
  # no apportionment needed
  if(attended_p && attended_a && attended_s) {
    return(array(tagged_fish, dim = c(n_pop, n_ages, n_sexes)))
  }

  # get platform and fleet release
  platform = tag_release_platform[1]
  fleet = as.integer(tag_release_platform[2])

  # Get raw tag apportionment weights [n_pop, n_ages, n_sexes]
  weights = if(platform == "population") {
    array(NAA[, tr, ty, tseas, , ], dim = c(n_pop, n_ages, n_sexes))
  } else if(platform == "fishery") {
    NAA_slice = array(NAA[, tr, ty, tseas, , ], dim = c(n_pop, n_ages, n_sexes))
    wt = array(0, dim = c(n_pop, n_ages, n_sexes))
    for(p in 1:n_pop) {
      sel_slice = array(fish_sel[p, tr, ty, tseas, , , fleet], dim = c(n_ages, n_sexes))
      wt[p,,] = NAA_slice[p,,] * sel_slice
    }
    wt
  } else if(platform == "survey") {
    NAA_slice = array(NAA[, tr, ty, tseas, , ], dim = c(n_pop, n_ages, n_sexes))
    wt = array(0, dim = c(n_pop, n_ages, n_sexes))
    for(p in 1:n_pop) {
      sel_slice = array(srv_sel[p, tr, ty, tseas, , , fleet], dim = c(n_ages, n_sexes))
      wt[p,,] = NAA_slice[p,,] * sel_slice
    }
    wt
  }

  # Normalize weights within attended dimensions.
  # The denominator sums over unattended dims only, preserving totals
  # within each attended dim combination.
  norm_weights = array(0, dim = c(n_pop, n_ages, n_sexes))
  for(p in seq_len(n_pop)) {
    for(a in seq_len(n_ages)) {
      for(s in seq_len(n_sexes)) {
        denom = 0
        for(pp in seq_len(n_pop)) {
          for(aa in seq_len(n_ages)) {
            for(ss in seq_len(n_sexes)) {
              same_p = !attended_p || pp == p
              same_a = !attended_a || aa == a
              same_s = !attended_s || ss == s
              if(same_p && same_a && same_s) denom = denom + weights[pp, aa, ss]
            }
          }
        }
        norm_weights[p, a, s] = weights[p, a, s] / denom
      }
    }
  }

  # Build output by looking up the correct cell in the original collapsed input.
  # Use index 1 for any unattended dim (those dims are collapsed to size 1 in
  # tagged_fish_orig), and the actual loop index for attended dims.
  out = array(0, dim = c(n_pop, n_ages, n_sexes))
  for(p in seq_len(n_pop)) {
    for(a in seq_len(n_ages)) {
      for(s in seq_len(n_sexes)) {
        pi = if(attended_p) p else 1
        ai = if(attended_a) a else 1
        si = if(attended_s) s else 1
        out[p, a, s] = tagged_fish_orig[pi, ai, si] * norm_weights[p, a, s]
      }
    }
  }

  return(out)

}
