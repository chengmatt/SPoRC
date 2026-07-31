# Stage 2 of 3: objective function
#
# Tag recapture likelihood assembly, with the same split between evaluation and
# one step ahead residual bookkeeping as the composition likelihoods.

#' Conventional Tagging Likelihood Block
#'
#' Computes the negative log-likelihood contribution from conventional
#' fish tagging data under Poisson, negative binomial, multinomial,
#' or Dirichlet–multinomial likelihoods.
#'
#' This function loops over tag cohorts, recapture years, seasons,
#' fleets, population pools, regions, age pools, and sex pools,
#' accumulating likelihood contributions into `conv_fish_tag_nLL`.
#'
#' @param n_conv_tag_cohorts Number of conventional tag cohorts.
#' @param conv_tag_release_indicator Matrix giving release region, year, season for each cohort.
#' @param conv_tag_max_liberty Maximum years at liberty to evaluate.
#' @param n_yrs Total number of modeled years.
#' @param n_seas Number of seasons per year.
#' @param conv_tag_mixing_period Minimum seasons at liberty before tags are modeled.
#' @param n_fish_fleets Number of fishing fleets.
#' @param use_conv_fish_tagging Vector indicating which fleets use conventional tagging.
#' @param n_conv_tag_pop_pool Number of population pooling groups.
#' @param n_regions Number of spatial regions.
#' @param n_conv_tag_age_pool Number of age pooling groups.
#' @param n_conv_tag_sex_pool Number of sex pooling groups.
#' @param conv_tag_pop_pool List of population index pools.
#' @param conv_tag_age_pool List of age index pools.
#' @param conv_tag_sex_pool List of sex index pools.
#' @param conv_fish_tag_like Likelihood type indicator (0–5).
#' @param conv_fish_tag_nLL Array of negative log-likelihood values to update.
#' @param obs_conv_tag_fish_recap Observed recapture array.
#' @param pred_conv_tag_fish_recap Predicted recapture array.
#' @param addtotag Small constant added to avoid zeros.
#' @param ln_conv_fish_tag_theta Log overdispersion parameter for NB/Dirichlet-multinomial.
#' @param conv_tagged_fish Array of numbers of tagged fish released.
#' @param zero_init Whether to zero the nLL array on entry.
#'
#' @returns Updated `conv_fish_tag_nLL` array.
#' @keywords internal
get_conv_tag_likelihoods <- function(n_conv_tag_cohorts,
                                     conv_tag_release_indicator,
                                     conv_tag_max_liberty,
                                     n_yrs,
                                     n_seas,
                                     conv_tag_mixing_period,
                                     n_fish_fleets,
                                     use_conv_fish_tagging,
                                     n_conv_tag_pop_pool,
                                     n_regions,
                                     n_conv_tag_age_pool,
                                     n_conv_tag_sex_pool,
                                     conv_tag_pop_pool,
                                     conv_tag_age_pool,
                                     conv_tag_sex_pool,
                                     conv_fish_tag_like,
                                     conv_fish_tag_nLL,
                                     obs_conv_tag_fish_recap,
                                     pred_conv_tag_fish_recap,
                                     addtotag,
                                     ln_conv_fish_tag_theta,
                                     conv_tagged_fish,
                                     zero_init = TRUE) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  d = dim(conv_fish_tag_nLL)
  if(zero_init) conv_fish_tag_nLL = RTMB::AD(as.vector(conv_fish_tag_nLL) * 0) else conv_fish_tag_nLL = RTMB::AD(as.vector(conv_fish_tag_nLL))
  dim(conv_fish_tag_nLL) = d

  # recapture arrays are [liberty, season, cohort, pop, region, age, sex, fleet],
  # so a single recovery event is the trailing five dimensions
  d_recap = dim(pred_conv_tag_fish_recap)[4:8]

  # number of pooled cells one event contributes, in (f, p, a, s, r) order to dimension vectors so they aren't 'growing'
  n_active_f = sum(use_conv_fish_tagging[1:n_fish_fleets] == 1)
  n_cells = n_active_f * n_conv_tag_pop_pool * n_conv_tag_age_pool * n_conv_tag_sex_pool * n_regions

  for(tc in 1:n_conv_tag_cohorts) {

    # set up tagging cohort indexing
    tr = conv_tag_release_indicator[tc,1] # extract tag release region
    ty = conv_tag_release_indicator[tc,2] # extract tag release year
    tseas = conv_tag_release_indicator[tc,3] # extract tag release season

    for(ry in 1:min(conv_tag_max_liberty, n_yrs - ty + 1)) { # loop through recapture years
      for(rseas in 1:n_seas) { # loop through recapture seasons

        # Dealing with tag mixing (not fitting to tags liberty < mixing period)
        # Skip seasons before release in the first year at liberty
        if(ry == 1 && rseas < tseas) next
        # Total seasonal time steps since release
        total_seas_at_liberty = (ry - 1) * n_seas + (rseas - tseas + 1)
        # Skip if within mixing period (in seasonal units)
        if(total_seas_at_liberty < conv_tag_mixing_period) next

        # This event's recaptures. Every option below pools these same
        # [pop, region, age, sex, fleet] numbers, so they are pulled out of the
        # cohort-spanning arrays once here to reduce overhead
        pred_ev = array(pred_conv_tag_fish_recap[ry,rseas,tc,,,,,], dim = d_recap)
        obs_ev  = array(obs_conv_tag_fish_recap[ry,rseas,tc,,,,,], dim = d_recap)

        # Poisson or Negative Binomial, one term per pooled cell
        if(conv_fish_tag_like %in% c(0, 1)) {
          for(f in 1:n_fish_fleets) {
            if(use_conv_fish_tagging[f] == 1) {
              for(r in 1:n_regions) {

                # every pool contributes a term to this region and fleet
                tmp_nLL = conv_fish_tag_nLL[ry,rseas,tc,r,f]

                for(p in 1:n_conv_tag_pop_pool) {
                  for(a in 1:n_conv_tag_age_pool) {
                    for(s in 1:n_conv_tag_sex_pool) {

                      pop_pool_idx = conv_tag_pop_pool[[p]] # extract movement pop pool indices
                      age_pool_idx = conv_tag_age_pool[[a]] # extract movement age pool indices
                      sex_pool_idx = conv_tag_sex_pool[[s]] # extract movement sex pool indices

                      # Poisson likelihood
                      if(conv_fish_tag_like == 0) {
                        tmp_nLL = tmp_nLL +
                          -dpois_noint(sum(obs_ev[pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag),
                                       sum(pred_ev[pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag),
                                       give_log = TRUE)
                      } # end if poisson likelihood

                      # Negative binomial likelihood
                      if(conv_fish_tag_like == 1) {
                        log_mu = log(sum(pred_ev[pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag)) # log mu
                        log_var_minus_mu = 2 * log_mu - ln_conv_fish_tag_theta # log var minus mu
                        tmp_nLL = tmp_nLL +
                          -dnbinom_robust_noint(x = sum(obs_ev[pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag),
                                                log_mu = log_mu, log_var_minus_mu = log_var_minus_mu, give_log = TRUE)
                      } # end if for negative binomial likelihood

                    } # end s loop
                  } # end a loop
                } # end p loop

                conv_fish_tag_nLL[ry,rseas,tc,r,f] = tmp_nLL

              } # end r loop
            } # end if
          } # end f loop
        } # end if count likelihood

        # # Release Conditioned for Multinomial or Dirichlet-Multinomial
        if(conv_fish_tag_like %in% c(2, 4)) {

          # Recaptured individuals, one cell per fleet, pool and region (containers).
          # Predicted cells are gathered in a list and joined in one go, which is
          # cheaper than writing them into a vector one at a time.
          tmp_pred_cells = vector("list", n_cells)
          tmp_obs_c_all = numeric(n_cells)
          ci = 0 # counter

          # number of tags released for a given tag cohort
          tmp_n_tags_released = sum(conv_tagged_fish[tc,,,] + addtotag)

          # Loop through age and sex pooling and combine vectors into the correct format
          for(f in 1:n_fish_fleets) {
            if(use_conv_fish_tagging[f] == 1) {
              for(p in 1:n_conv_tag_pop_pool) {
                for(a in 1:n_conv_tag_age_pool) {
                  for(s in 1:n_conv_tag_sex_pool) {

                    pop_pool_idx = conv_tag_pop_pool[[p]] # extract movement pop pool indices
                    age_pool_idx = conv_tag_age_pool[[a]] # extract movement age pool indices
                    sex_pool_idx = conv_tag_sex_pool[[s]] # extract movement sex pool indices

                    # Pool observed and expected if any pooling
                    for (r in 1:n_regions) {
                      ci = ci + 1
                      tmp_pred_cells[[ci]] = sum(pred_ev[pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                      tmp_obs_c_all[ci]    = sum(obs_ev[pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                    }

                  } # end s loop
                } # end a loop
              } # end p loop
            }
          }

          tmp_pred_c_all = do.call(c, tmp_pred_cells) / tmp_n_tags_released
          tmp_obs_c_all = tmp_obs_c_all / tmp_n_tags_released

          tmp_pred = c(tmp_pred_c_all, 1 - sum(tmp_pred_c_all))
          tmp_obs = c(tmp_obs_c_all, 1 - sum(tmp_obs_c_all))

          if(conv_fish_tag_like == 2)
            conv_fish_tag_nLL[ry,rseas,tc,1,1] = -tmp_n_tags_released * sum((tmp_obs) * log(tmp_pred))

          if(conv_fish_tag_like == 4)
            conv_fish_tag_nLL[ry,rseas,tc,1,1] =
            -1 * ddirmult(obs = tmp_obs, pred = tmp_pred,
                          Ntotal = tmp_n_tags_released,
                          ln_theta = ln_conv_fish_tag_theta, TRUE)

        } # end release conditioned

        # Recapture Conditioned (Multinomial or Dirichlet-Multinomial)
        if(conv_fish_tag_like %in% c(3,5)) {

          # setup containers
          tmp_pred_cells = vector("list", n_cells)
          tmp_obs_all = numeric(n_cells)
          ci = 0 # counter

          tmp_n_tags_recap = sum(obs_ev + addtotag)

          for(f in 1:n_fish_fleets) {
            if(use_conv_fish_tagging[f] == 1) {
              for(p in 1:n_conv_tag_pop_pool) {
                for(a in 1:n_conv_tag_age_pool) {
                  for(s in 1:n_conv_tag_sex_pool) {

                    pop_pool_idx = conv_tag_pop_pool[[p]]
                    age_pool_idx = conv_tag_age_pool[[a]]
                    sex_pool_idx = conv_tag_sex_pool[[s]]

                    for (r in 1:n_regions) {
                      ci = ci + 1
                      tmp_pred_cells[[ci]] = sum(pred_ev[pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                      tmp_obs_all[ci]      = sum(obs_ev[pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                    }

                  }
                }
              }
            }
          }

          tmp_pred_all = do.call(c, tmp_pred_cells)
          tmp_pred_all = tmp_pred_all / sum(tmp_pred_all)
          tmp_obs_all = tmp_obs_all / tmp_n_tags_recap

          if(conv_fish_tag_like == 3)
            conv_fish_tag_nLL[ry,rseas,tc,1,1] =
            -1 * tmp_n_tags_recap * sum(tmp_obs_all * log(tmp_pred_all))

          if(conv_fish_tag_like == 5)
            conv_fish_tag_nLL[ry,rseas,tc,1,1] =
            -1 * ddirmult(obs = tmp_obs_all, pred = tmp_pred_all,
                          Ntotal = tmp_n_tags_recap,
                          ln_theta = ln_conv_fish_tag_theta, TRUE)

        } # end recapture conditioned

      } # end rseas
    } # end ry

  } # end tc

  conv_fish_tag_nLL
}


#' Classify conventional tag likelihood family
#'
#' Maps a conventional tag likelihood type code to a coarse
#' family label used for OSA packing and evaluation.
#'
#' @param lt Integer likelihood type code:
#'   \itemize{
#'     \item 0, 1: count-based (Poisson / NB)
#'     \item 2, 3, 4, 5: composition-based (multinomial / Dirichlet-multinomial)
#'   }
#'
#' @return A character scalar: `"count"`, `"comp"`, or `NA_character_` if
#'   the code is not recognized.
#' @keywords internal
tag_fam_of = function(lt) {
  if(lt %in% c(0,1)) "count"
  else if(lt %in% c(2,3,4,5)) "comp"
  else NA_character_
}

#' Enumerate valid conventional-tag recovery events
#'
#' Enumerates all valid \code{(tc, ry, rseas)} recovery events applying the
#' exact skip logic used in the fitting loop (release season and mixing
#' period). Returns a \code{data.frame} with one row per event, in loop order.
#'
#' @param n_conv_tag_cohorts Number of conventional tag cohorts.
#' @param conv_tag_release_indicator Matrix giving release region, year, season for each cohort.
#' @param conv_tag_max_liberty Maximum years at liberty to evaluate.
#' @param n_yrs Total number of modeled years.
#' @param n_seas Number of seasons per year.
#' @param conv_tag_mixing_period Minimum seasons at liberty before tags are modeled.
#'
#' @return A \code{data.frame} with columns \code{tc}, \code{ry}, \code{rseas},
#'   \code{tr}, \code{ty}, and \code{tseas}, in the same order the fitting loop
#'   would visit them. If no valid events exist, returns an empty data.frame.
#' @keywords internal
tag_grid = function(n_conv_tag_cohorts, conv_tag_release_indicator,
                    conv_tag_max_liberty, n_yrs, n_seas, conv_tag_mixing_period) {
  rows = list()
  for(tc in 1:n_conv_tag_cohorts) {
    tr    = conv_tag_release_indicator[tc,1]
    ty    = conv_tag_release_indicator[tc,2]
    tseas = conv_tag_release_indicator[tc,3]
    for(ry in 1:min(conv_tag_max_liberty, n_yrs - ty + 1)) {
      for(rseas in 1:n_seas) {
        if(ry == 1 && rseas < tseas) next
        total_seas_at_liberty = (ry - 1) * n_seas + (rseas - tseas + 1)
        if(total_seas_at_liberty < conv_tag_mixing_period) next
        rows[[length(rows) + 1]] = c(tc = tc, ry = ry, rseas = rseas,
                                     tr = tr, ty = ty, tseas = tseas)
      }
    }
  }
  do.call(rbind, lapply(rows, function(z) as.data.frame(as.list(z))))
}

#' Pack conventional-tag observations for OSA
#'
#' Packs conventional-tag observations into a flat vector suitable for
#' one-step-ahead (OSA) analysis for a single likelihood family.
#'
#' For \code{family == "count"}, each valid event contributes one integer
#' observation per \code{[fleet, pop_pool, region, age_pool, sex_pool]} cell
#' with \code{use_fish_tagging[f] == 1}, in \code{(f, p, r, a, s)} loop order
#' within each event, and events in \code{tag_grid} order. This mirrors
#' \code{get_conv_tag_likelihoods()}'s exact accumulation structure -- a
#' separate Poisson/NB term is fit per pool/age/sex cell and their
#' -log-likelihoods summed -- so packing one count per \code{[region, fleet]}
#' pre-summed across pools would evaluate a different (non-equivalent)
#' likelihood whenever more than one pool is used.
#'
#' For \code{family == "comp"}, each valid event contributes one composition
#' vector:
#' \itemize{
#'   \item release-conditioned (\code{like_type} 2, 4): recap cells in
#'     \code{(f, p, a, s, r)} loop order plus a non-recapture tail; counts
#'     are \code{round(prop * n_tags_released)}.
#'   \item recapture-conditioned (\code{like_type} 3, 5): recap cells only,
#'     conditioned on total recaptures; counts are \code{round(prop * n_tags_recap)}.
#' }
#'
#' @param family Character, either \code{"count"} or \code{"comp"}.
#' @param like_type Integer likelihood type code (0–5).
#' @param obs_recap Observed recapture array.
#' @param pred_recap Predicted recapture array (used for scaling).
#' @param tagged_fish Array of numbers of tagged fish released.
#' @param conv_tag_release_indicator Matrix giving release region, year, season for each cohort.
#' @param conv_tag_max_liberty Maximum years at liberty to evaluate.
#' @param n_conv_tag_cohorts Number of conventional tag cohorts.
#' @param n_yrs Total number of modeled years.
#' @param n_seas Number of seasons per year.
#' @param n_regions Number of spatial regions.
#' @param n_fish_fleets Number of fishing fleets.
#' @param n_pop_pool Number of population pooling groups.
#' @param n_age_pool Number of age pooling groups.
#' @param n_sex_pool Number of sex pooling groups.
#' @param pop_pool List of population index pools.
#' @param age_pool List of age index pools.
#' @param sex_pool List of sex index pools.
#' @param use_fish_tagging Vector indicating which fleets use conventional tagging.
#' @param conv_tag_mixing_period Minimum seasons at liberty before tags are modeled.
#' @param addtotag Small constant added to avoid zeros.
#' @param return_labels Logical; if TRUE, also builds a per-element label
#'   data.frame identifying the origin (family, like_type, tag cohort/release
#'   region-year-season, recovery year/season, fleet, region, pop/age/sex pool,
#'   is_tail, last_in_group) of every entry in \code{vec}, in the same order.
#'   Intended for post-hoc relabeling of \code{TMB::oneStepPredict()} residuals
#'   (see [get_osa()]); left \code{FALSE} (default) inside the model itself to
#'   avoid the extra bookkeeping cost.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{vec}: flat numeric/AD vector of packed observations, or \code{NULL} if no events.
#'   \item \code{grp_end}: integer vector of end indices for each composition group
#'     (empty for \code{family == "count"}).
#'   \item \code{lengths}: integer vector of per-group lengths.
#'   \item \code{labels}: data.frame with one row per element of \code{vec}
#'     (\code{NULL} unless \code{return_labels = TRUE}).
#' }
#' @keywords internal
pack_tag_osa = function(family, like_type,
                        obs_recap, pred_recap, tagged_fish,
                        conv_tag_release_indicator, conv_tag_max_liberty,
                        n_conv_tag_cohorts, n_yrs, n_seas, n_regions, n_fish_fleets,
                        n_pop_pool, n_age_pool, n_sex_pool,
                        pop_pool, age_pool, sex_pool,
                        use_fish_tagging, conv_tag_mixing_period, addtotag,
                        return_labels = FALSE) {
  "c" <- RTMB::ADoverload("c")

  grid = tag_grid(n_conv_tag_cohorts, conv_tag_release_indicator,
                  conv_tag_max_liberty, n_yrs, n_seas, conv_tag_mixing_period)
  if(is.null(grid) || nrow(grid) == 0) return(list(vec = NULL, grp_end = integer(0), lengths = integer(0), labels = NULL))

  clean = list(); grp_end = integer(0); lengths = integer(0); pos = 0L
  label_rows = list()

  # one label row per element, mirroring the exact loop order used to build 'clean'
  tag_label_row = function(fleet, region, pop_pool_i, age_pool_i, sex_pool_i,
                           tc, ry, rseas, tr, ty, tseas, is_tail, last_in_group) {
    data.frame(family = family, like_type = like_type,
               tc = tc, ry = ry, rseas = rseas, tr = tr, ty = ty, tseas = tseas,
               fleet = fleet, region = region,
               pop_pool = pop_pool_i, age_pool = age_pool_i, sex_pool = sex_pool_i,
               is_tail = is_tail, last_in_group = last_in_group)
  }

  for(g in 1:nrow(grid)) {
    tc = grid$tc[g]; ry = grid$ry[g]; rseas = grid$rseas[g]
    tr = grid$tr[g]; ty = grid$ty[g]; tseas = grid$tseas[g]

    if(family == "count") {
      # one cell per [f,p,r,a,s] -- must mirror get_conv_tag_likelihoods()'s
      # exact accumulation structure: a separate Poisson/NB term is evaluated
      # per (pool, age_pool, sex_pool, region, fleet) cell and their
      # -log-likelihoods summed. Packing one combined count per [r,f] (summed
      # across pools first, as this used to do) is a DIFFERENT likelihood --
      # not equivalent even for Poisson -- whenever n_pop_pool/n_age_pool/
      # n_sex_pool > 1, and also destroys population identity before it ever
      # reaches the residual/label, so downstream plots can't facet by pop.
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(p in 1:n_pop_pool) {
          for(r in 1:n_regions) {
            for(a in 1:n_age_pool) {
              for(s in 1:n_sex_pool) {
                v = sum(obs_recap[ry, rseas, tc, pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
                clean[[length(clean)+1]] = round(v)   # integer count for cdf
                pos = pos + 1L
                if(return_labels) {
                  label_rows[[length(label_rows)+1]] = tag_label_row(
                    fleet = f, region = r, pop_pool_i = p, age_pool_i = a, sex_pool_i = s,
                    tc = tc, ry = ry, rseas = rseas, tr = tr, ty = ty, tseas = tseas,
                    is_tail = FALSE, last_in_group = NA
                  )
                }
              }
            }
          }
        }
      }
      lengths = c(lengths, 0L)  # count has no per-group determined bin

    } else {

      # composition (multinomial or dm): build recap-cell vector in fit-loop order (f,p,a,s,r). Preallocate to avoid 'growing'
      n_active_f = sum(use_fish_tagging[1:n_fish_fleets] == 1)
      n_cells = n_active_f * n_pop_pool * n_age_pool * n_sex_pool * n_regions
      obs_cells = numeric(n_cells)
      ci = 0 # counter

      if(return_labels) cell_labels = vector("list", n_cells)
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(p in 1:n_pop_pool) for(a in 1:n_age_pool) for(s in 1:n_sex_pool) {
          for(r in 1:n_regions) {
            ci = ci + 1
            obs_cells[ci] = sum(obs_recap[ry, rseas, tc, pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
            if(return_labels) {
              cell_labels[[ci]] = tag_label_row(
                fleet = f, region = r, pop_pool_i = p, age_pool_i = a, sex_pool_i = s,
                tc = tc, ry = ry, rseas = rseas, tr = tr, ty = ty, tseas = tseas,
                is_tail = FALSE, last_in_group = FALSE
              )
            }
          }
        }
      }

      if(like_type %in% c(2,4)) {
        # release-conditioned: normalize by tags released, append non-recap
        n_rel = sum(tagged_fish[tc,,,] + addtotag)
        prop  = obs_cells / n_rel
        tail  = 1 - sum(prop)
        if(tail < 0) tail = 0                     # guard to avoid rounding that pushes tail <0
        prop  = c(prop, tail)                     # non-recap tail (determined)
        prop  = prop / sum(prop)                  # renormalize after guard
        g_counts = round(prop * n_rel)
        if(return_labels) {
          cell_labels[[length(cell_labels)+1]] = tag_label_row(
            fleet = NA_integer_, region = NA_integer_, pop_pool_i = NA_integer_, age_pool_i = NA_integer_, sex_pool_i = NA_integer_,
            tc = tc, ry = ry, rseas = rseas, tr = tr, ty = ty, tseas = tseas,
            is_tail = TRUE, last_in_group = TRUE
          )
        }
      } else {
        # recapture-conditioned: normalize by total recaptures, no tail
        n_recap = sum(obs_recap[ry, rseas, tc, , , , , ] + addtotag)
        prop = obs_cells / n_recap
        g_counts = round(prop * n_recap)          # last cell determined
        if(return_labels) cell_labels[[n_cells]]$last_in_group = TRUE
      }

      L = length(g_counts)
      clean[[length(clean)+1]] = g_counts
      lengths = c(lengths, L)
      grp_end = c(grp_end, pos + L)               # determined bin = last element
      pos = pos + L
      if(return_labels) label_rows[[length(label_rows)+1]] = do.call(rbind, cell_labels)
    }
  }

  vec = do.call(c, clean)
  labels = if(return_labels) do.call(rbind, label_rows) else NULL
  list(vec = vec, grp_end = grp_end, lengths = lengths, labels = labels)
}

#' Evaluate conventional-tag OSA likelihood from a tracked vector
#'
#' Evaluates the one-step-ahead (OSA) negative log-likelihood for a single
#' conventional-tag likelihood family from a flat tracked observation vector.
#' Fills \code{nLL_arr} in the same grid order used by \code{pack_tag_osa()}.
#'
#' The Dirichlet-multinomial concentration matches the fitting likelihood: the
#' total \eqn{N} is the number of tags released for release-conditioned
#' families (\code{like_type} 2, 4) and the total recaptures for
#' recapture-conditioned families (\code{like_type} 3, 5). Both totals are taken
#' from raw data (\code{tagged_fish}, \code{obs_recap}) rather than from the
#' tracked observation, since summing an OSA-tagged S4 object is not permitted.
#'
#' @param nLL_arr Array of negative log-likelihood values to update.
#' @param tracked Flat tracked observation vector produced by \code{pack_tag_osa()}.
#' @param family Character, either \code{"count"} or \code{"comp"}.
#' @param like_type Integer likelihood type code (0-5).
#' @param pred_recap Predicted recapture array.
#' @param obs_recap Observed recapture array. Used to derive the recapture total
#'   \eqn{N} for recapture-conditioned Dirichlet-multinomial (\code{like_type} 5).
#' @param tagged_fish Array of numbers of tagged fish released.
#' @param conv_tag_release_indicator Matrix giving release region, year, season for each cohort.
#' @param conv_tag_max_liberty Maximum years at liberty to evaluate.
#' @param n_conv_tag_cohorts Number of conventional tag cohorts.
#' @param n_yrs Total number of modeled years.
#' @param n_seas Number of seasons per year.
#' @param n_regions Number of spatial regions.
#' @param n_fish_fleets Number of fishing fleets.
#' @param n_pop_pool Number of population pooling groups.
#' @param n_age_pool Number of age pooling groups.
#' @param n_sex_pool Number of sex pooling groups.
#' @param pop_pool List of population index pools.
#' @param age_pool List of age index pools.
#' @param sex_pool List of sex index pools.
#' @param use_fish_tagging Vector indicating which fleets use conventional tagging.
#' @param conv_tag_mixing_period Minimum seasons at liberty before tags are modeled.
#' @param addtotag Small constant added to avoid zeros.
#' @param ln_theta Log overdispersion parameter for the Dirichlet-multinomial
#'   composition family.
#' @param zero_init Logical; if \code{TRUE}, zero \code{nLL_arr} on entry.
#'
#' @return Updated \code{nLL_arr} array with OSA negative log-likelihood
#'   contributions filled in loop/grid order.
#' @keywords internal
eval_tag_osa = function(nLL_arr, tracked, family, like_type,
                        pred_recap, obs_recap, tagged_fish,
                        conv_tag_release_indicator, conv_tag_max_liberty,
                        n_conv_tag_cohorts, n_yrs, n_seas, n_regions, n_fish_fleets,
                        n_pop_pool, n_age_pool, n_sex_pool,
                        pop_pool, age_pool, sex_pool,
                        use_fish_tagging, conv_tag_mixing_period, addtotag,
                        ln_theta = 0, zero_init = TRUE) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  if(is.null(tracked)) return(nLL_arr)

  d = dim(nLL_arr)
  if(zero_init) nLL_arr = RTMB::AD(as.vector(nLL_arr) * 0) else nLL_arr = RTMB::AD(as.vector(nLL_arr))
  dim(nLL_arr) = d

  grid = tag_grid(n_conv_tag_cohorts, conv_tag_release_indicator,
                  conv_tag_max_liberty, n_yrs, n_seas, conv_tag_mixing_period)

  # recapture arrays are [liberty, season, cohort, pop, region, age, sex, fleet],
  # so a single recovery event is the trailing five dimensions
  d_recap = dim(pred_recap)[4:8]

  k = 1
  for(g in 1:nrow(grid)) {
    tc = grid$tc[g]; ry = grid$ry[g]; rseas = grid$rseas[g]

    # this event's predicted recaptures, pooled by both branches below
    pred_ev = array(pred_recap[ry,rseas,tc,,,,,], dim = d_recap)

    if(family == "count") {
      # mirrors pack_tag_osa()'s [f,p,r,a,s] cell order -- one Poisson/NB
      # term per (pool, age_pool, sex_pool, region, fleet) cell, accumulated
      # into the same [ry,rseas,tc,r,f] nLL_arr cell (matching
      # get_conv_tag_likelihoods()'s `+=` accumulation across pools).
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(p in 1:n_pop_pool) {
          for(r in 1:n_regions) {
            for(a in 1:n_age_pool) {
              for(s in 1:n_sex_pool) {
                mu = sum(pred_ev[pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
                if(like_type == 0) {          # Poisson
                  term = -RTMB::dpois(tracked[k], mu, log = TRUE)
                } else {                      # NB (robust): var = mu + mu^2/theta -> size = exp(ln_theta)
                  term = -RTMB::dnbinom_robust(x = tracked[k], log_mu = log(mu),
                                               log_var_minus_mu = 2 * log(mu) - ln_theta, log = TRUE)
                }
                nLL_arr[ry,rseas,tc,r,f] = nLL_arr[ry,rseas,tc,r,f] + term
                k = k + 1
              }
            }
          }
        }
      }

    } else {
      # composition: rebuild predicted proportions in the same order, then
      # dmultinom_osa / ddirmult_osa. build as a list rather than grow as a vec, to reduce overhead
      n_active_f = sum(use_fish_tagging[1:n_fish_fleets] == 1)
      n_cells = n_active_f * n_pop_pool * n_age_pool * n_sex_pool * n_regions
      cells = vector("list", n_cells)
      ci = 0
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(p in 1:n_pop_pool) for(a in 1:n_age_pool) for(s in 1:n_sex_pool) {
          for(r in 1:n_regions) {
            ci = ci + 1
            cells[[ci]] = sum(pred_ev[pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
          }
        }
      }
      pred_cells = do.call(c, cells)

      # DM total N matches fitting: released (2,4) vs total recaptures (3,5).
      # Both taken from RAW data -- never sum(tracked[idx]) (S4 not summable).
      if(like_type %in% c(2,4)) {
        Ntot_dm = sum(tagged_fish[tc, , , ] + addtotag)          # tags released
        pprop   = pred_cells / Ntot_dm
        pprop   = c(pprop, 1 - sum(pprop))                        # non-recap tail (determined)
      } else {
        Ntot_dm = sum(obs_recap[ry, rseas, tc, , , , , ] + addtotag)  # total recaptures
        pprop   = pred_cells / sum(pred_cells)
      }
      pprop = pprop / sum(pprop)                                  # ensure sums to 1

      L   = length(pprop)
      idx = k:(k + L - 1)

      if(like_type %in% c(2,3)) {                   # Multinomial (discrete)
        nLL_arr[ry,rseas,tc,1,1] = -dmultinom_osa(tracked[idx], pprop, log = TRUE)
      } else {                                      # Dirichlet-multinomial (4,5)
        nLL_arr[ry,rseas,tc,1,1] = -ddirmult_osa(tracked[idx], pprop * exp(ln_theta) * Ntot_dm, log = TRUE)
      }
      k = k + L
    }
  }
  nLL_arr
}
