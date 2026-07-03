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

        for(f in 1:n_fish_fleets) {
          if(use_conv_fish_tagging[f] == 1) {
            for(p in 1:n_conv_tag_pop_pool) {
              for(r in 1:n_regions) {
                for(a in 1:n_conv_tag_age_pool) {
                  for(s in 1:n_conv_tag_sex_pool) {

                    pop_pool_idx = conv_tag_pop_pool[[p]] # extract movement pop pool indices
                    age_pool_idx = conv_tag_age_pool[[a]] # extract movement age pool indices
                    sex_pool_idx = conv_tag_sex_pool[[s]] # extract movement sex pool indices

                    # Poisson likelihood
                    if(conv_fish_tag_like == 0) {
                      conv_fish_tag_nLL[ry,rseas,tc,r,f] = conv_fish_tag_nLL[ry,rseas,tc,r,f]  +
                        -dpois_noint(sum(obs_conv_tag_fish_recap[ry,rseas,tc,pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag),
                                     sum(pred_conv_tag_fish_recap[ry,rseas,tc,pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag),
                                     give_log = TRUE)
                    } # end if poisson likelihood

                    # Negative binomial likelihood
                    if(conv_fish_tag_like == 1) {
                      log_mu = log(sum(pred_conv_tag_fish_recap[ry,rseas,tc,pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag)) # log mu
                      log_var_minus_mu = 2 * log_mu - ln_conv_fish_tag_theta # log var minus mu
                      conv_fish_tag_nLL[ry,rseas,tc,r,f] = conv_fish_tag_nLL[ry,rseas,tc,r,f] +
                        -dnbinom_robust_noint(x = sum(obs_conv_tag_fish_recap[ry,rseas,tc,pop_pool_idx,r,age_pool_idx,sex_pool_idx,f] + addtotag),
                                              log_mu = log_mu, log_var_minus_mu = log_var_minus_mu, give_log = TRUE)
                    } # end if for negative binomial likelihood

                  } # end s loop
                } # end a loop
              } # end r loop
            } # end p loop
          } # end if
        } # end f loop

        # # Release Conditioned for Multinomial or Dirichlet-Multinomial
        if(conv_fish_tag_like %in% c(2, 4)) {

          # Temporary vectors for recaptured individuals
          tmp_pred_c_all = vector()
          tmp_obs_c_all = vector()

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
                      pred_val = sum(pred_conv_tag_fish_recap[ry, rseas, tc, pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                      obs_val  = sum(obs_conv_tag_fish_recap[ry, rseas, tc, pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                      tmp_pred_c_all = c(tmp_pred_c_all, pred_val)
                      tmp_obs_c_all  = c(tmp_obs_c_all,  obs_val)
                    }

                  } # end s loop
                } # end a loop
              } # end p loop
            }
          }

          tmp_pred_c_all = tmp_pred_c_all / tmp_n_tags_released
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

          tmp_pred_all = vector()
          tmp_obs_all = vector()

          tmp_n_tags_recap = sum(obs_conv_tag_fish_recap[ry,rseas,tc,,,,,] + addtotag)

          for(f in 1:n_fish_fleets) {
            if(use_conv_fish_tagging[f] == 1) {
              for(p in 1:n_conv_tag_pop_pool) {
                for(a in 1:n_conv_tag_age_pool) {
                  for(s in 1:n_conv_tag_sex_pool) {

                    pop_pool_idx = conv_tag_pop_pool[[p]]
                    age_pool_idx = conv_tag_age_pool[[a]]
                    sex_pool_idx = conv_tag_sex_pool[[s]]

                    for (r in 1:n_regions) {
                      pred_val = sum(pred_conv_tag_fish_recap[ry, rseas, tc, pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                      obs_val  = sum(obs_conv_tag_fish_recap[ry, rseas, tc, pop_pool_idx, r, age_pool_idx, sex_pool_idx, f] + addtotag)
                      tmp_pred_all = c(tmp_pred_all, pred_val)
                      tmp_obs_all  = c(tmp_obs_all,  obs_val)
                    }

                  }
                }
              }
            }
          }

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
#' observation per \code{[region, fleet]} with \code{use_fish_tagging[f] == 1},
#' in region-fastest then fleet order within each event, and events in
#' \code{tag_grid} order.
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
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{vec}: flat numeric/AD vector of packed observations, or \code{NULL} if no events.
#'   \item \code{grp_end}: integer vector of end indices for each composition group
#'     (empty for \code{family == "count"}).
#'   \item \code{lengths}: integer vector of per-group lengths.
#' }
#' @keywords internal
pack_tag_osa = function(family, like_type,
                        obs_recap, pred_recap, tagged_fish,
                        conv_tag_release_indicator, conv_tag_max_liberty,
                        n_conv_tag_cohorts, n_yrs, n_seas, n_regions, n_fish_fleets,
                        n_pop_pool, n_age_pool, n_sex_pool,
                        pop_pool, age_pool, sex_pool,
                        use_fish_tagging, conv_tag_mixing_period, addtotag) {
  "c" <- RTMB::ADoverload("c")

  grid = tag_grid(n_conv_tag_cohorts, conv_tag_release_indicator,
                  conv_tag_max_liberty, n_yrs, n_seas, conv_tag_mixing_period)
  if(is.null(grid) || nrow(grid) == 0) return(list(vec = NULL, grp_end = integer(0), lengths = integer(0)))

  clean = list(); grp_end = integer(0); lengths = integer(0); pos = 0L

  for(g in 1:nrow(grid)) {
    tc = grid$tc[g]; ry = grid$ry[g]; rseas = grid$rseas[g]

    if(family == "count") {
      # one cell per [r,f] with active fleet
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(r in 1:n_regions) {
          v = 0
          for(p in 1:n_pop_pool) for(a in 1:n_age_pool) for(s in 1:n_sex_pool) {
            v = v + sum(obs_recap[ry, rseas, tc, pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
          }
          clean[[length(clean)+1]] = round(v)   # integer count for cdf
          pos = pos + 1L
        }
      }
      lengths = c(lengths, 0L)  # count has no per-group determined bin

    } else {
      # composition (comp or dm): build recap-cell vector in fit-loop order (f,p,a,s,r).
      # obs are DATA (numeric). Preallocate to avoid growing from a NULL seed (the
      # ADoverloaded c() cannot advector() NULL).
      n_active_f = sum(use_fish_tagging[1:n_fish_fleets] == 1)
      n_cells = n_active_f * n_pop_pool * n_age_pool * n_sex_pool * n_regions
      obs_cells = numeric(n_cells)
      ci = 0
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(p in 1:n_pop_pool) for(a in 1:n_age_pool) for(s in 1:n_sex_pool) {
          for(r in 1:n_regions) {
            ci = ci + 1
            obs_cells[ci] = sum(obs_recap[ry, rseas, tc, pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
          }
        }
      }

      if(like_type %in% c(2,4)) {
        # release-conditioned: normalize by tags released, append non-recap
        n_rel = sum(tagged_fish[tc,,,] + addtotag)
        prop  = obs_cells / n_rel
        tail  = 1 - sum(prop)
        if(tail < 0) tail = 0                     # guard: rounding can push tail <0
        prop  = c(prop, tail)                     # non-recap tail (determined)
        prop  = prop / sum(prop)                  # renormalize after guard
        g_counts = round(prop * n_rel)
      } else {
        # recapture-conditioned: normalize by total recaptures, no tail
        n_recap = sum(obs_recap[ry, rseas, tc, , , , , ] + addtotag)
        prop = obs_cells / n_recap
        g_counts = round(prop * n_recap)          # last cell determined
      }

      L = length(g_counts)
      clean[[length(clean)+1]] = g_counts
      lengths = c(lengths, L)
      grp_end = c(grp_end, pos + L)               # determined bin = last element
      pos = pos + L
    }
  }

  vec = do.call(c, clean)
  list(vec = vec, grp_end = grp_end, lengths = lengths)
}

#' Evaluate conventional-tag OSA likelihood from a tracked vector
#'
#' Evaluates the one-step-ahead (OSA) negative log-likelihood for a single
#' conventional-tag likelihood family from a flat tracked observation vector.
#' Fills \code{nLL_arr} in the same grid order used by \code{pack_tag_osa()}.
#'
#' @param nLL_arr Array of negative log-likelihood values to update.
#' @param tracked Flat tracked observation vector produced by \code{pack_tag_osa()}.
#' @param family Character, either \code{"count"} or \code{"comp"}.
#' @param like_type Integer likelihood type code (0–5).
#' @param pred_recap Predicted recapture array.
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
#' @param ln_theta Log overdispersion parameter for NB/Dirichlet-multinomial (composition family).
#' @param zero_init Logical; if \code{TRUE}, zero \code{nLL_arr} on entry.
#'
#' @return Updated \code{nLL_arr} array with OSA negative log-likelihood
#'   contributions filled in loop/grid order.
#' @importFrom RTMBdist ddirmult
#' @keywords internal
eval_tag_osa = function(nLL_arr, tracked, family, like_type,
                        pred_recap, tagged_fish,
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

  k = 1
  for(g in 1:nrow(grid)) {
    tc = grid$tc[g]; ry = grid$ry[g]; rseas = grid$rseas[g]

    if(family == "count") {
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(r in 1:n_regions) {
          mu = 0
          for(p in 1:n_pop_pool) for(a in 1:n_age_pool) for(s in 1:n_sex_pool) {
            mu = mu + sum(pred_recap[ry, rseas, tc, pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
          }
          if(like_type == 0) {          # Poisson
            nLL_arr[ry,rseas,tc,r,f] = -RTMB::dpois(tracked[k], mu, log = TRUE)
          } else {                      # NB (robust): var = mu + mu^2/theta -> size = exp(ln_theta)
            nLL_arr[ry,rseas,tc,r,f] = -RTMB::dnbinom_robust(x = tracked[k], log_mu = log(mu),
                                                             log_var_minus_mu = 2 * log(mu) - ln_theta, log = TRUE)
          }
          k = k + 1
        }
      }

    } else {
      # composition: rebuild predicted proportions in the same order, then dmultinom / ddirmult.
      # Preallocate an AD vector (don't grow from empty c() -- the ADoverloaded c()
      # cannot advector() a NULL seed).
      n_active_f = sum(use_fish_tagging[1:n_fish_fleets] == 1)
      n_cells = n_active_f * n_pop_pool * n_age_pool * n_sex_pool * n_regions
      pred_cells = RTMB::AD(numeric(n_cells))
      ci = 0
      for(f in 1:n_fish_fleets) {
        if(use_fish_tagging[f] != 1) next
        for(p in 1:n_pop_pool) for(a in 1:n_age_pool) for(s in 1:n_sex_pool) {
          for(r in 1:n_regions) {
            ci = ci + 1
            pred_cells[ci] = sum(pred_recap[ry, rseas, tc, pop_pool[[p]], r, age_pool[[a]], sex_pool[[s]], f] + addtotag)
          }
        }
      }

      if(like_type %in% c(2,4)) {
        n_rel = sum(tagged_fish[tc,,,] + addtotag)
        pprop = pred_cells / n_rel
        pprop = c(pprop, 1 - sum(pprop))            # non-recap tail
      } else {
        pprop = pred_cells / sum(pred_cells)
      }
      pprop = pprop / sum(pprop)                    # ensure sums to 1

      L = length(pprop)
      idx = k:(k + L - 1)

      if(like_type %in% c(2,3)) {                   # Multinomial (discrete)
        nLL_arr[ry,rseas,tc,1,1] = -RTMB::dmultinom(tracked[idx], prob = pprop, log = TRUE)
      } else {                                      # DM (4,5): fit works, OSA gated externally
        nLL_arr[ry,rseas,tc,1,1] = -RTMBdist::ddirmult(tracked[idx], sum(tracked[idx]), pprop * exp(ln_theta) * sum(tracked[idx]), log = TRUE)
      }
      k = k + L
    }
  }
  nLL_arr
}
