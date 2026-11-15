# Stage 2 of 3: objective function
#
# Predicted conventional tag recaptures, following each release cohort forward through movement,
# mortality, ageing, shedding and reporting.

#' Fishing and natural mortality by tagged fish in one year and season
#'
#' Builds the mortality-at-age components a tag cohort experiences over a single
#' season. These depend only on the year and season, not on which cohort is at
#' liberty, so both \code{\link{get_tagging_observation_model}} and the
#' simulation's \code{generate_fishery_conv_tags_recap} work them out once per
#' year and season and read them back for every cohort.
#'
#' Only fleets flagged in \code{use_conv_fish_tagging} contribute; the others are
#' left at zero, since tags are only ever returned by fleets that report them.
#'
#' \code{Z_before_shed} excludes tag shedding, which is a property of the release
#' cohort rather than of the year and season. Callers add
#' \code{exp(ln_conv_tag_shed[tc]) * seasdur[rseas]} to it.
#'
#' Every array argument is expected without a trailing simulation dimension, so
#' the operating model slices by \code{sim} before calling.
#'
#' @param y Integer year index.
#' @param rseas Integer season index.
#' @param n_pop,n_regions,n_ages,n_sexes,n_fish_fleets Dimension sizes.
#' @param use_conv_fish_tagging Integer vector \code{[fish_fleet]} (0/1)
#'   flagging fleets with tagging data.
#' @param Fmort Array \code{[region, year, season, fish_fleet]} of fishing
#'   mortality.
#' @param fish_sel,ret_sel Arrays \code{[pop, region, year, season, age, sex,
#'   fish_fleet]} of total/retained fishery selectivity.
#' @param dmr Array \code{[region, year, season, fish_fleet]} of discard
#'   mortality rate.
#' @param natmort Array \code{[pop, region, year, season, age, sex]} of natural
#'   mortality at age.
#' @param seasdur Numeric vector \code{[season]} of season duration (fraction
#'   of year).
#'
#' @return List with elements \code{FAA}, \code{ret_FAA}, \code{disc_DAA}, each
#'   \code{[pop, region, 1, age, sex, fish_fleet]}, and \code{Z_before_shed},
#'   \code{[pop, region, 1, age, sex]}.
#'
#' @keywords internal
#' @import RTMB
get_tag_mort <- function(y, rseas, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets,
                         use_conv_fish_tagging, Fmort, fish_sel, ret_sel, dmr, natmort, seasdur) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

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

  # get natural mortality in the recapture season
  tmp_natmort <- array(natmort[,,y,rseas,,], dim = c(n_pop, n_regions, 1, n_ages, n_sexes))

  return(list(
    FAA = tmp_FAA,
    ret_FAA = tmp_ret_FAA,
    disc_DAA = tmp_disc_DAA,
    Z_before_shed = (tmp_natmort * seasdur[rseas]) + apply(tmp_FAA, 1:5, sum)
  ))
}


# Tag Cohort Season Steps ---------------------------------------------------

# One season for the tag cohorts at liberty, in three steps both get_tagging_observation_model and the simulation's
# generate_fishery_conv_tags_recap call, so the two follow tagged fish the same way: each cohort's mortality and
# release, then every cohort's movement together, then each cohort's ageing and recaptures

#' Tag cohorts at liberty in a year and season
#'
#' The cohorts released at or before this season that are still within their
#' years at liberty.
#'
#' @param y,rseas Year and season.
#' @param conv_tag_release_indicator Matrix \code{[cohort, 3]} of release region,
#'   year and season.
#' @param conv_tag_max_liberty Integer, the most years a cohort is followed.
#'
#' @return Integer vector of cohort indices.
#' @keywords internal
get_tag_cohorts_at_liberty <- function(y, rseas, conv_tag_release_indicator, conv_tag_max_liberty) {
  tag_release_yr <- conv_tag_release_indicator[,2] # tag release year
  tag_release_seas <- conv_tag_release_indicator[,3] # tag release season
  released <- tag_release_yr < y | (tag_release_yr == y & tag_release_seas <= rseas) # released by this season
  within_liberty <- (y - tag_release_yr + 1) <= conv_tag_max_liberty # not past the years followed
  which(released & within_liberty)
}

#' Tag shedding groups
#'
#' Cohorts mapped to the same estimated shedding parameter, or fixed at the same
#' shedding value, are in one group, numbered by the first cohort in it. Under
#' \code{move_timing = 2} the cohorts of a group share one seasonal operator.
#'
#' @param conv_tag_shed_map Numeric vector \code{[cohort]}, the map of
#'   \code{ln_conv_tag_shed} (\code{NA} where it is fixed), or \code{NULL}, which
#'   makes every cohort its own group.
#' @param ln_conv_tag_shed Vector \code{[cohort]} of log shedding rates.
#'
#' @return Integer vector \code{[cohort]} of group numbers.
#' @keywords internal
get_tag_shed_group <- function(conv_tag_shed_map, ln_conv_tag_shed) {

  # setup tag shed groups
  n_cohorts <- length(ln_conv_tag_shed)
  tag_shed_group <- 1:n_cohorts
  if(is.null(conv_tag_shed_map)) return(tag_shed_group)

  shed_value <- ad_value(ln_conv_tag_shed) # make into plain numbers, only compared for fixed rates
  for(tc in 1:n_cohorts) {
    for(tc_first in 1:tc) {

      # figure out whether fixed, estimated, same fixed value, or same map
      both_fixed <- is.na(conv_tag_shed_map[tc]) && is.na(conv_tag_shed_map[tc_first])
      both_estimated <- !is.na(conv_tag_shed_map[tc]) && !is.na(conv_tag_shed_map[tc_first])
      same_fixed_value <- both_fixed && shed_value[tc] == shed_value[tc_first]
      same_parameter <- both_estimated && conv_tag_shed_map[tc] == conv_tag_shed_map[tc_first]

      # if same then part of the same group
      if(same_fixed_value || same_parameter) {
        tag_shed_group[tc] <- tag_shed_group[tc_first]
        break
      }
    } # end tc_first loop
  } # end tc loop

  tag_shed_group
}

#' Start of a season for one tag cohort: mortality and release
#'
#' Total mortality with the cohort's shedding, the share of the season it is at
#' liberty for, the release in its release season, and survival.
#'
#' @param avail_tc Array \code{[liberty + 1, season, pop, region, age, sex]} of
#'   the cohort's tags.
#' @param tagged_fish Array \code{[pop, age, sex]} of fish released, read in the
#'   release season only.
#' @param mort Output of \code{get_tag_mort} for this year and season.
#' @param ry,rseas Years at liberty and season.
#' @param tr,tseas Release region and season.
#' @param ln_conv_tag_shed_tc,conv_tag_t_tagging_tc,ln_init_conv_tag_mort_tc The
#'   cohort's log shedding rate, share of the release season remaining at release,
#'   and log tag-induced mortality.
#' @param seasdur Numeric vector \code{[season]} of season durations.
#' @param move_timing Where movement sits relative to mortality (0, 1 or 2).
#' @param n_pop,n_ages,n_sexes Dimension sizes.
#'
#' @return List with \code{avail_tc}, the cohort's tags with this season's
#'   release in, and \code{seas_vals}, what the rest of the season reads:
#'   \code{ry}, \code{tmp_ZAA}, \code{tmp_ret_FAA}, \code{tmp_SAA},
#'   \code{tag_frac}, \code{tag_dur}, \code{tag_moves} and \code{tag_moves_seas}.
#' @keywords internal
start_tag_cohort_season <- function(avail_tc,
                                    tagged_fish,
                                    mort,
                                    ry,
                                    rseas,
                                    tr,
                                    tseas,
                                    ln_conv_tag_shed_tc,
                                    conv_tag_t_tagging_tc,
                                    ln_init_conv_tag_mort_tc,
                                    seasdur,
                                    move_timing,
                                    n_pop,
                                    n_ages,
                                    n_sexes) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # get retained fishing mortality, and total mortality with this cohort's shedding rate
  tmp_ret_FAA <- mort$ret_FAA
  tmp_ZAA <- mort$Z_before_shed + (exp(ln_conv_tag_shed_tc) * seasdur[rseas])

  # Fraction of this season the tag cohort is actually at liberty for. Mid-season releases
  # only experience the remainder of the season.
  tag_frac <- if(ry == 1 && rseas == tseas) conv_tag_t_tagging_tc else 1
  tag_dur <- seasdur[rseas] * tag_frac

  # discount by tagging time when tagging is not at the start of the season. scale F and Z
  # together: scaling Z alone makes F/Z exceed one and overpredicts recaptures
  if(tag_frac != 1) {
    tmp_ZAA     <- tmp_ZAA     * tag_frac
    tmp_ret_FAA <- tmp_ret_FAA * tag_frac
  }

  # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced
  # mortality (exponential mortality process)
  if(ry == 1 && rseas == tseas) {
    avail_tc[1, rseas, , tr, , ] <- array(tagged_fish * exp(-exp(ln_init_conv_tag_mort_tc)), dim = c(n_pop, n_ages, n_sexes))
  }

  # get temporary survival value
  tmp_SAA <- exp(-tmp_ZAA)

  # Whether the discrete movement step applies this season. A mid-season release skips it,
  # since a full-season transition matrix cannot represent a partial interval.
  tag_moves <- (conv_tag_t_tagging_tc == 1 || ry != 1 || rseas != tseas)

  # continuous movement needs no such exemption: the generator is scaled by tag_dur, so a
  # mid-season release diffuses for the same fraction of the season its mortality is scaled by
  tag_moves_seas <- if(move_timing == 2) TRUE else tag_moves

  seas_vals <- list(ry = ry, # years at liberty
                    tmp_ZAA = tmp_ZAA, # total mortality, with shedding
                    tmp_ret_FAA = tmp_ret_FAA, # retained fishing mortality
                    tmp_SAA = tmp_SAA, # survival over the season
                    tag_frac = tag_frac, # share of the season at liberty
                    tag_dur = tag_dur, # time at liberty this season
                    tag_moves = tag_moves, # discrete movement applies
                    tag_moves_seas = tag_moves_seas) # movement applies under the seasonal operator
  list(avail_tc = avail_tc, seas_vals = seas_vals)
}

#' Move a set of tag vectors through one matrix
#'
#' Every vector is a cohort's tags by region, moved one cohort at a time,
#' \code{t(tag_vec) \%*\% move_mat} for each.
#'
#' @param tag_vecs List of numeric vectors \code{[region]}.
#' @param move_mat Matrix \code{[region_from, region_to]}.
#'
#' @return List of the moved vectors \code{[region_to]}, in the same order.
#' @keywords internal
move_tag_vectors <- function(tag_vecs, move_mat) {
  moved_vecs <- vector("list", length(tag_vecs))
  for(k in seq_along(tag_vecs)) moved_vecs[[k]] <- as.vector(t(tag_vecs[[k]]) %*% move_mat) # one product each
  moved_vecs
}

#' Movement for every tag cohort at liberty in a season
#'
#' Every cohort at liberty goes through the season's movement matrix, or its
#' seasonal operator under \code{move_timing = 2}, through
#' \code{move_tag_vectors}. Under \code{move_timing = 0} the cohorts' tags are
#' moved in place, and mortality follows in \code{end_tag_cohort_season}. Under
#' \code{move_timing = 1} each cohort is survived and then moved. Under
#' \code{move_timing = 2} the cohorts of a shedding group with the same time at
#' liberty share one seasonal operator and its integral.
#'
#' @param tag_avail_cohort List \code{[cohort]} of each cohort's tags.
#' @param cohort_seas_vals List \code{[cohort]} of \code{seas_vals} from
#'   \code{start_tag_cohort_season}.
#' @param cohorts_at_liberty Integer vector of the cohorts at liberty this season.
#' @param rseas Season.
#' @param move_seas Array \code{[pop, region_from, region_to, age, sex]}, the
#'   season's movement matrices.
#' @param rate_seas Array \code{[pop, region_from, region_to, age, sex]}, the
#'   season's movement rates under \code{move_timing = 2}, else \code{NULL}.
#' @param tag_shed_group Integer vector \code{[cohort]} from
#'   \code{get_tag_shed_group}.
#' @param move_timing,do_recruits_move,expm_nsub Movement settings.
#' @param n_pop,n_regions,n_ages,n_sexes Dimension sizes.
#'
#' @return List with \code{tag_avail_cohort} (moved under \code{move_timing = 0}),
#'   and \code{tag_step_cohort} and \code{tag_int_cohort}, each cohort's tags at the
#'   end of the season and integrated over it under \code{move_timing} 1 and 2.
#' @keywords internal
move_tag_cohorts <- function(tag_avail_cohort,
                             cohort_seas_vals,
                             cohorts_at_liberty,
                             rseas,
                             move_seas,
                             rate_seas,
                             tag_shed_group,
                             move_timing,
                             do_recruits_move,
                             expm_nsub,
                             n_pop,
                             n_regions,
                             n_ages,
                             n_sexes) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_cohorts <- length(tag_avail_cohort)
  tag_step_cohort <- vector("list", n_cohorts) # tags at the end of the season, by cohort (move_timing 1 and 2)
  tag_int_cohort <- vector("list", n_cohorts) # tags integrated over the season, by cohort (move_timing 2)

  # move tagged fish around. under move_timing == 0 movement is applied here and mortality follows
  if(move_timing == 0) {

    # setup dimensions to figure out how many cohorts to move
    tag_move_ages <- if(do_recruits_move == 0) 2:n_ages else 1:n_ages # ages that move
    moving_cohorts <- integer(0)
    for(tc in cohorts_at_liberty) if(cohort_seas_vals[[tc]]$tag_moves) moving_cohorts <- c(moving_cohorts, tc)
    n_moving_cohorts <- length(moving_cohorts)

    for(p in 1:n_pop) {

      # skip if no cohorts to move
      if(n_moving_cohorts == 0) next

      # input each moving cohort's tags this season, regions by ages by sexes
      moving_tags <- vector("list", n_moving_cohorts)
      for(k in 1:n_moving_cohorts) {
        tc <- moving_cohorts[k]
        moving_tags[[k]] <- array(tag_avail_cohort[[tc]][cohort_seas_vals[[tc]]$ry, rseas, p, , , ], dim = c(n_regions, n_ages, n_sexes))
      } # end k loop

      # every moving cohort at an age and sex through the same matrix
      for(a in tag_move_ages) {
        for(s in 1:n_sexes) {
          tag_vecs <- vector("list", n_moving_cohorts)
          for(k in 1:n_moving_cohorts) tag_vecs[[k]] <- moving_tags[[k]][, a, s]
          moved_vecs <- move_tag_vectors(tag_vecs = tag_vecs, move_mat = move_seas[p, , , a, s])
          for(k in 1:n_moving_cohorts) moving_tags[[k]][, a, s] <- moved_vecs[[k]]
        } # end s loop
      } # end a loop

      # input the moved tags back into each cohort's array
      for(k in 1:n_moving_cohorts) {
        tc <- moving_cohorts[k]
        tmp_avail_tc <- tag_avail_cohort[[tc]]
        tmp_avail_tc[cohort_seas_vals[[tc]]$ry, rseas, p, , , ] <- moving_tags[[k]]
        tag_avail_cohort[[tc]] <- tmp_avail_tc

      } # end k loop
    } # end p loop
  } # end if move_timing == 0

  # under move_timing 1 and 2 the seasonal transition operator does movement and mortality together
  if(move_timing %in% c(1, 2) && n_regions > 1) {
    for(tc in cohorts_at_liberty) {
      tag_step_cohort[[tc]] <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
      if(move_timing == 2) tag_int_cohort[[tc]] <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
    } # end tc loop

    # under 2, the tag cohorts that share a seasonal operator so do the same shedding group and the same time at liberty this season.
    operator_groups <- list() # container
    if(move_timing == 2) {
      ungrouped <- cohorts_at_liberty # cohorts not in a group yet
      while(length(ungrouped) > 0) {
        group_first_tc <- ungrouped[1]
        group_cohorts <- integer(0)
        for(tc in ungrouped) {
          same_shed_group <- tag_shed_group[tc] == tag_shed_group[group_first_tc]
          same_time_at_liberty <- cohort_seas_vals[[tc]]$tag_dur == cohort_seas_vals[[group_first_tc]]$tag_dur
          if(same_shed_group && same_time_at_liberty) group_cohorts <- c(group_cohorts, tc)
        } # end tc loop
        operator_groups[[length(operator_groups) + 1]] <- group_cohorts
        ungrouped <- setdiff(ungrouped, group_cohorts)
      } # end while ungrouped
    } # end if move_timing == 2

    for(p in 1:n_pop) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {

          # move_timing == 1: each cohort survived, then the ones that move this season moved
          if(move_timing == 1) {
            moving_cohorts <- integer(0)
            survived_vecs <- list()
            for(tc in cohorts_at_liberty) {
              tag_vec <- tag_avail_cohort[[tc]][cohort_seas_vals[[tc]]$ry, rseas, p, , a, s] # tags at the start of the season
              survived_vec <- tag_vec * exp(-cohort_seas_vals[[tc]]$tmp_ZAA[p,,1,a,s]) # survival over the season, with shedding
              tag_step_cohort[[tc]][p,,a,s] <- survived_vec # a cohort that does not move is only survived
              if(cohort_seas_vals[[tc]]$tag_moves_seas && (do_recruits_move == 1 || a > 1)) {
                moving_cohorts <- c(moving_cohorts, tc)
                survived_vecs[[length(moving_cohorts)]] <- survived_vec
              }
            } # end tc loop
            if(length(moving_cohorts) > 0) {
              moved_vecs <- move_tag_vectors(tag_vecs = survived_vecs, move_mat = move_seas[p, , , a, s])
              for(k in seq_along(moving_cohorts)) tag_step_cohort[[moving_cohorts[k]]][p,,a,s] <- moved_vecs[[k]]
            }
          } # end if move_timing == 1

          # move_timing == 2: continuous movement, the generator scaled by tag_dur, so a mid-season release diffuses
          # for the same fraction of the season its mortality is scaled by
          if(move_timing == 2) {
            moves <- do_recruits_move == 1 || a > 1
            Qv <- if(moves) rate_seas[p, , , a, s] else matrix(0, n_regions, n_regions)
            for(group_cohorts in operator_groups) {

              # the group's seasonal operator, from the first cohort's total mortality and time at liberty
              group_first_tc <- group_cohorts[1]
              tag_seas_operator <- seas_operator_and_integral(Z = cohort_seas_vals[[group_first_tc]]$tmp_ZAA[p,,1,a,s], # total mortality, with shedding
                                                              Q = Qv, # movement rates
                                                              dur = cohort_seas_vals[[group_first_tc]]$tag_dur, # time at liberty
                                                              expm_nsub = expm_nsub) # matrix exponential substeps

              # every cohort in the group through it: tags at the end of the season, and integrated over it
              tag_vecs <- vector("list", length(group_cohorts))
              for(k in seq_along(group_cohorts)) {
                tc <- group_cohorts[k]
                tag_vecs[[k]] <- tag_avail_cohort[[tc]][cohort_seas_vals[[tc]]$ry, rseas, p, , a, s]
              } # end k loop

              # get fish that moved and died
              stepped_vecs <- move_tag_vectors(tag_vecs = tag_vecs, move_mat = tag_seas_operator$T)
              integrated_vecs <- move_tag_vectors(tag_vecs = tag_vecs, move_mat = t(tag_seas_operator$Integral))

              # input back into array
              for(k in seq_along(group_cohorts)) {
                tag_step_cohort[[group_cohorts[k]]][p,,a,s] <- stepped_vecs[[k]]
                tag_int_cohort[[group_cohorts[k]]][p,,a,s] <- integrated_vecs[[k]]
              } # end k loop

            } # end group_cohorts loop
          } # end if move_timing == 2

        } # end s loop
      } # end a loop
    } # end p loop
  } # end if move_timing 1 or 2

  list(tag_avail_cohort = tag_avail_cohort, tag_step_cohort = tag_step_cohort, tag_int_cohort = tag_int_cohort)
}

#' End of a season for one tag cohort: ageing and recaptures
#'
#' Survival (under \code{move_timing = 0}, after movement), the shift to the next
#' season or the next year and age, the state-space numbers at age factor, and the
#' predicted recaptures. A season the cohort is at liberty for none of
#' (\code{conv_tag_t_tagging = 0} at release) has no recaptures.
#'
#' @param avail_tc,recap_tc Arrays of the cohort's tags
#'   \code{[liberty + 1, season, pop, region, age, sex]} and predicted recaptures
#'   \code{[liberty, season, pop, region, age, sex, fish_fleet]}.
#' @param seas_vals The cohort's \code{seas_vals} from
#'   \code{start_tag_cohort_season}.
#' @param tag_step,tag_int The cohort's tags at the end of the season and
#'   integrated over it, from \code{move_tag_cohorts} (\code{NULL} under
#'   \code{move_timing = 0}).
#' @param y,rseas,n_seas Year, season and number of seasons.
#' @param reporting_y Array \code{[region, fish_fleet]} of tag reporting rates in
#'   year \code{y}.
#' @param NAA_scalar Array \code{[pop, region, year, season, age, sex]} of the
#'   state-space numbers at age factor, or \code{NULL}.
#' @param move_timing Where movement sits relative to mortality (0, 1 or 2).
#' @param n_pop,n_regions,n_ages,n_sexes,n_fish_fleets Dimension sizes.
#'
#' @return List with \code{avail_tc} and \code{recap_tc}.
#' @keywords internal
end_tag_cohort_season <- function(avail_tc,
                                  recap_tc,
                                  seas_vals,
                                  tag_step,
                                  tag_int,
                                  y,
                                  rseas,
                                  n_seas,
                                  reporting_y,
                                  NAA_scalar,
                                  move_timing,
                                  n_pop,
                                  n_regions,
                                  n_ages,
                                  n_sexes,
                                  n_fish_fleets) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  ry <- seas_vals$ry
  tmp_ZAA <- seas_vals$tmp_ZAA
  tmp_ret_FAA <- seas_vals$tmp_ret_FAA
  tmp_SAA <- seas_vals$tmp_SAA

  # post-season tag numbers, before the ageing shift (from move_tag_cohorts under move_timing 1 and 2)
  if(move_timing == 0 || n_regions == 1) { # move then mortality
    tag_step <- array(avail_tc[ry, rseas, , , , ] * tmp_SAA[,,1,,], dim = c(n_pop, n_regions, n_ages, n_sexes))
  }

  # Apply mortality and ageing to tagged fish
  if(rseas < n_seas) {

    # Season mortality within a given year, advance to next season same year/age
    avail_tc[ry, rseas + 1, , , , ] <- tag_step

    # State-space numbers at age rescale tag cohorts at a within-year boundary too
    if(!is.null(NAA_scalar) && y <= dim(NAA_scalar)[3]) {
      avail_tc[ry, rseas + 1, , , , ] <-
        avail_tc[ry, rseas + 1, , , , ] * NAA_scalar[,,y,rseas + 1,,]
    }

  } else {

    # End of year mortality and age advancement (end of season)
    avail_tc[ry + 1, 1, , , 2:n_ages, ] <- tag_step[,,1:(n_ages - 1),]

    # Accumulate plus group
    avail_tc[ry + 1, 1, , , n_ages, ] <-
      avail_tc[ry + 1, 1, , , n_ages, ] + tag_step[,,n_ages,]

    # State-space numbers at age rescale tag cohorts as well if applied. The factor is
    # indexed by the calendar year the cohort has reached, not by its years at liberty.
    if(!is.null(NAA_scalar) && (y + 1) <= dim(NAA_scalar)[3]) {
      avail_tc[ry + 1, 1, , , 2:n_ages, ] <-
        avail_tc[ry + 1, 1, , , 2:n_ages, ] * NAA_scalar[,,y + 1,1,2:n_ages,]
    }
  }

  # Apply Baranov's to get predicted recaptures. a cohort at liberty for none of the season (released at its end)
  # catches nothing, where F / Z would be 0 / 0
  if(seas_vals$tag_frac > 0) {
    for(p in 1:n_pop) {

      # spatial Baranov: tags redistribute among regions while being caught, so recaptures use the season-integrated tag abundance
      for(f in 1:n_fish_fleets) {
        if(move_timing == 2) {
          tmp_ret_FAA_slice <- array(tmp_ret_FAA[p,,1,,,f], dim = c(n_regions, n_ages, n_sexes))
          tag_int_p <- array(tag_int[p,,,], dim = c(n_regions, n_ages, n_sexes))
          recap_tc[ry,rseas,p,,,,f] <- reporting_y[,f] * tmp_ret_FAA_slice * tag_int_p
        } else {
          recap_tc[ry,rseas,p,,,,f] <- reporting_y[,f] *
            (tmp_ret_FAA[p,,1,,,f] / tmp_ZAA[p,,1,,]) *
            avail_tc[ry,rseas,p,,,] *
            (1 - tmp_SAA[p,,1,,])
        }
      } # end f loop
    } # end p loop
  } # end if at liberty this season

  list(avail_tc = avail_tc, recap_tc = recap_tc)
}

# Tagging Observation Model -------------------------------------------------

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
#' Cohorts are followed calendar year by calendar year, so the cohorts at
#' liberty in a season share that season's movement matrix, and under
#' \code{move_timing = 2} a shedding group's seasonal operator is worked out
#' once for all of its cohorts with the same time at liberty.
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
#' @param natmort Array \code{[pop, region, year, season, age, sex]} of natural
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
#' @param NAA_scalar Array \code{[pop, region, year, season, age, sex]} of the
#'   factor the state-space numbers at age applied to the deterministic
#'   prediction, one wherever it did not apply. Tagged fish are part of the same
#'   population, and the innovation acts as mortality the model does not
#'   otherwise account for, so every cohort takes the same factor at each
#'   boundary, within a year as well as between years. \code{NULL} (the default)
#'   leaves them on the
#'   deterministic trajectory, which is correct only when the state is off.
#' @param conv_tag_shed_map Numeric vector \code{[cohort]}, the map of
#'   \code{ln_conv_tag_shed} (\code{NA} where it is fixed). Under
#'   \code{move_timing = 2}, cohorts with the same map level, or the same fixed
#'   value, share one seasonal operator. \code{NULL} gives every cohort its own.
#'
#' @return List with elements \code{conv_tag_fish_reporting},
#'   \code{conv_tag_fish_avail}, \code{pred_conv_tag_fish_recap}.
#'
#' @keywords internal
#' @import RTMB
get_tagging_observation_model <- function(
  n_fish_fleets,
  n_regions,
  n_conv_tag_cohorts,
  n_yrs,
  n_seas,
  n_pop,
  n_ages,
  n_sexes,
  conv_tag_fish_reporting_blocks,
  conv_tag_fish_reporting_pars,
  conv_tag_fish_reporting,
  conv_tag_release_indicator,
  conv_tag_max_liberty,
  use_conv_fish_tagging,
  Fmort,
  fish_sel,
  ret_sel,
  dmr,
  natmort,
  seasdur,
  ln_conv_tag_shed,
  conv_tag_t_tagging,
  conv_tagged_fish,
  conv_fish_tag_attr,
  conv_tag_release_platform,
  srv_sel,
  NAA_bef,
  ln_init_conv_tag_mort,
  do_recruits_move,
  Movement,
  conv_tag_fish_avail,
  pred_conv_tag_fish_recap,
  Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0,
  NAA_scalar = NULL,
  conv_tag_shed_map = NULL
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Set up tag reporting rates
  for(f in 1:n_fish_fleets) {
    for(r in 1:n_regions) {
      conv_tagrep_blk_idx <- conv_tag_fish_reporting_blocks[r,,f]  # Get all blocks for this region
      conv_tag_fish_reporting[r,,f] <- RTMB::plogis(conv_tag_fish_reporting_pars[r,conv_tagrep_blk_idx,f])  # inverse logit transform
    } # end r loop
  } # end f loop

  # Mortality at a given year and season is the same for every cohort at liberty
  # then, so it is worked out once per year and season and read back below.
  tag_mort <- vector("list", n_yrs * n_seas)
  for(tc in 1:n_conv_tag_cohorts) {
    ty <- conv_tag_release_indicator[tc,2]
    for(ry in 1:min(conv_tag_max_liberty, n_yrs - ty + 1)) {
      y <- ty + ry - 1
      for(rseas in 1:n_seas) {
        idx <- (y - 1) * n_seas + rseas
        if(is.null(tag_mort[[idx]]))
          tag_mort[[idx]] <- get_tag_mort(
            y = y,
            rseas = rseas,
            n_pop = n_pop,
            n_regions = n_regions,
            n_ages = n_ages,
            n_sexes = n_sexes,
            n_fish_fleets = n_fish_fleets,
            use_conv_fish_tagging = use_conv_fish_tagging,
            Fmort = Fmort,
            fish_sel = fish_sel,
            ret_sel = ret_sel,
            dmr = dmr,
            natmort = natmort,
            seasdur = seasdur
          )
      } # end rseas loop
    } # end ry loop
  } # end tc loop

  # Cohort specific containers
  tag_avail_cohort <- vector("list", n_conv_tag_cohorts)
  tag_recap_cohort <- vector("list", n_conv_tag_cohorts)
  for(tc in 1:n_conv_tag_cohorts) {
    tag_avail_cohort[[tc]] <- array(0, dim = c(conv_tag_max_liberty + 1, n_seas, n_pop, n_regions, n_ages, n_sexes))
    tag_recap_cohort[[tc]] <- array(0, dim = c(conv_tag_max_liberty, n_seas, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets))
  } # end tc loop

  # cohorts that share a shedding rate, and so a seasonal operator under move_timing == 2
  tag_shed_group <- get_tag_shed_group(conv_tag_shed_map, ln_conv_tag_shed)

  # Cohorts are followed calendar year by calendar year, each season in three steps the simulation shares
  # (generate_fishery_conv_tags_recap): each cohort's mortality and release, every cohort's movement together, then
  # each cohort's ageing and recaptures
  for(y in 1:n_yrs) {
    for(rseas in 1:n_seas) {

      cohorts_at_liberty <- get_tag_cohorts_at_liberty(y, rseas, conv_tag_release_indicator, conv_tag_max_liberty)
      if(length(cohorts_at_liberty) == 0) next

      # get fishing and natural mortality
      mort <- tag_mort[[(y - 1) * n_seas + rseas]]
      cohort_seas_vals <- vector("list", n_conv_tag_cohorts) # each cohort's values this season

      # Mortality and release, cohort by cohort
      for(tc in cohorts_at_liberty) {

        tr <- conv_tag_release_indicator[tc,1] # extract tag release region
        ty <- conv_tag_release_indicator[tc,2] # extract tag release year
        tseas <- conv_tag_release_indicator[tc,3] # extract tag release season
        ry <- y - ty + 1 # years at liberty

        # apportion tagged fish out to appropriate dimensions if necessary, in the release season
        tmp_tagged_fish <- NULL
        if(ry == 1 && rseas == tseas) {
          tmp_tagged_fish <- release_conv_tag_attr(array(conv_tagged_fish[tc, , , ], dim = c(n_pop, n_ages, n_sexes)),
                                                  conv_fish_tag_attr[tc],
                                                  conv_tag_release_platform[tc,],
                                                  srv_sel, fish_sel, NAA_bef,
                                                  ty, tseas, tr, n_pop,
                                                  n_ages, n_sexes)
        }

        season_start <- start_tag_cohort_season(avail_tc = tag_avail_cohort[[tc]], # this cohort's tags
                                                tagged_fish = tmp_tagged_fish, # released fish, release season only
                                                mort = mort, # fishing and natural mortality this season
                                                ry = ry, # years at liberty
                                                rseas = rseas, # season
                                                tr = tr, # release region
                                                tseas = tseas, # release season
                                                ln_conv_tag_shed_tc = ln_conv_tag_shed[tc], # shedding rate
                                                conv_tag_t_tagging_tc = conv_tag_t_tagging[tc], # share of the release season left
                                                ln_init_conv_tag_mort_tc = ln_init_conv_tag_mort[tc], # tag induced mortality
                                                seasdur = seasdur, # season durations
                                                move_timing = move_timing, # movement timing
                                                n_pop = n_pop,
                                                n_ages = n_ages,
                                                n_sexes = n_sexes)
        tag_avail_cohort[[tc]] <- season_start$avail_tc
        cohort_seas_vals[[tc]] <- season_start$seas_vals
      } # end tc loop

      # Movement, every cohort at liberty together
      move_seas <- array(Movement[,,,y,rseas,,], dim = c(n_pop, n_regions, n_regions, n_ages, n_sexes))
      rate_seas <- if(move_timing == 2) array(Mrate[,,,y,rseas,,], dim = c(n_pop, n_regions, n_regions, n_ages, n_sexes)) else NULL
      moved <- move_tag_cohorts(tag_avail_cohort = tag_avail_cohort, # each cohort's tags
                                cohort_seas_vals = cohort_seas_vals, # each cohort's values this season
                                cohorts_at_liberty = cohorts_at_liberty, # cohorts at liberty this season
                                rseas = rseas, # season
                                move_seas = move_seas, # movement matrices this season
                                rate_seas = rate_seas, # movement rates this season (move_timing 2)
                                tag_shed_group = tag_shed_group, # shedding groups
                                move_timing = move_timing, # movement timing
                                do_recruits_move = do_recruits_move, # whether age one moves
                                expm_nsub = expm_nsub, # matrix exponential substeps
                                n_pop = n_pop,
                                n_regions = n_regions,
                                n_ages = n_ages,
                                n_sexes = n_sexes)
      tag_avail_cohort <- moved$tag_avail_cohort

      # Survival, ageing and recaptures, cohort by cohort
      reporting_y <- array(conv_tag_fish_reporting[,y,], dim = c(n_regions, n_fish_fleets))
      for(tc in cohorts_at_liberty) {
        season_end <- end_tag_cohort_season(avail_tc = tag_avail_cohort[[tc]], # this cohort's tags
                                            recap_tc = tag_recap_cohort[[tc]], # this cohort's recaptures
                                            seas_vals = cohort_seas_vals[[tc]], # this cohort's values this season
                                            tag_step = moved$tag_step_cohort[[tc]], # tags at the end of the season
                                            tag_int = moved$tag_int_cohort[[tc]], # tags integrated over the season
                                            y = y, # year
                                            rseas = rseas, # season
                                            n_seas = n_seas, # number of seasons
                                            reporting_y = reporting_y, # reporting rates this year
                                            NAA_scalar = NAA_scalar, # state-space numbers at age factor
                                            move_timing = move_timing, # movement timing
                                            n_pop = n_pop,
                                            n_regions = n_regions,
                                            n_ages = n_ages,
                                            n_sexes = n_sexes,
                                            n_fish_fleets = n_fish_fleets)
        tag_avail_cohort[[tc]] <- season_end$avail_tc
        tag_recap_cohort[[tc]] <- season_end$recap_tc
      } # end tc loop

    } # end rseas loop
  } # end y loop

  # Store each cohort's tags and predicted recaptures
  for(tc in 1:n_conv_tag_cohorts) {
    conv_tag_fish_avail[,,tc,,,,] <- tag_avail_cohort[[tc]]
    pred_conv_tag_fish_recap[,,tc,,,,,] <- tag_recap_cohort[[tc]]
  } # end tc loop

  return(list(conv_tag_fish_reporting = conv_tag_fish_reporting,
              conv_tag_fish_avail = conv_tag_fish_avail,
              pred_conv_tag_fish_recap = pred_conv_tag_fish_recap))
}
#' Distribute Tagged Fish Releases to Full Population Dimensions
#'
#' Spreads tag totals recorded at less than full resolution out to
#' \code{[n_pop, n_ages, n_sexes]}, using weights from the release platform.
#' \code{tagged_fish} comes back unchanged when every dim is attended
#' (\code{tag_attr = "p_a_s"}).
#'
#' The weights are the numbers at age under \code{"population"}, or the numbers
#' at age times fishery or survey selectivity under the other two platforms, all
#' read at the release region, year and season. They are normalized conditionally
#' on the attended dims: the denominator for cell \code{[p, a, s]} sums the raw
#' weights over every cell sharing its indices in those dims, so the marginal
#' totals of \code{tagged_fish} are preserved exactly along each of them. With
#' only age attended, for instance, the age totals are kept and the tags are
#' spread over population and sex in proportion to the weights.
#'
#' @param tagged_fish Released tagged fish for one cohort, with the unattended
#'   dims collapsed to index 1 and reshaped internally to \code{[n_pop, n_ages,
#'   n_sexes]}.
#' @param tag_attr Which dims are attended, built from \code{"p"}, \code{"a"} and
#'   \code{"s"} joined by underscores. See \code{\link{Setup_Mod_Tagging}}.
#' @param tag_release_platform Character vector of length 2: the platform,
#'   \code{"population"}, \code{"fishery"} or \code{"survey"}, and the fleet index
#'   as a string, or \code{NA} under \code{"population"}.
#' @param srv_sel,fish_sel Survey and fishery selectivity \code{[n_regions,
#'   n_yrs, n_ages, n_sexes, n_fleets]}, used as the age and sex weight under
#'   their own platform.
#' @param NAA Numbers at age before movement \code{[n_pop, n_regions, n_yrs,
#'   n_seas, n_ages, n_sexes]}, the weights themselves under
#'   \code{"population"} and multiplied by selectivity otherwise.
#' @param ty,tseas,tr Year, season and region indices of the release cohort.
#' @param n_pop,n_ages,n_sexes Dimension sizes.
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

  # return tagged_fish (reshaped to full dims) if p_a_s, all dims attended,
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

  # normalize weights within attended dims. the denominator sums over unattended dims only,
  # which preserves totals within each attended dim combination
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

  # look up the matching cell in the original collapsed input: index 1 for any unattended dim,
  # since those are collapsed to size 1 in tagged_fish_orig, and the loop index otherwise
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
