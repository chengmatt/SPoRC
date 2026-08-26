# Stage 2 of 3: objective function
#
# Conditional age-at-length likelihoods. A CAAL observation is an age composition
# within one length bin, so everything here loops length bins and hands each row to
# the composition functions in model_lik_comps.R rather than restating the
# multinomial and Dirichlet-multinomial algebra. Get_CAAL_Likelihoods evaluates the
# likelihood; pack_caal_osa and eval_caal_osa carry the one step ahead residual
# bookkeeping.
#

#' Conditional Age-at-Length Likelihood
#'
#' Computes the negative log-likelihood contribution for conditional
#' age-at-length (CAAL) data for a single year, season and fleet. A CAAL
#' observation is the age composition of the fish sampled from one length bin, so
#' each length bin is treated as an independent age composition and evaluated
#' through \code{\link{Get_Comp_Likelihoods}}.
#'
#'
#' @param Exp Expected joint numbers at length and age (a \code{Fish_caal},
#'   \code{Fish_caal_discard} or \code{Srv_caal} slice), indexed by
#'   \eqn{[region \times len \times model\_age \times sex]}.
#' @param Obs Observed CAAL counts or proportions indexed by
#'   \eqn{[region \times len \times observed\_age \times sex]}.
#' @param ISS Input sample size indexed by \eqn{[region \times len \times sex]}.
#'   This is the number aged within the length bin, not the number measured.
#' @param Wt_Mltnml Multinomial weighting indexed by
#'   \eqn{[region \times len \times sex]}.
#' @param ln_theta Log overdispersion for the Dirichlet-multinomial, indexed by
#'   \eqn{[region \times sex]}. One value per region and sex is shared across
#'   length bins, since the bins come from one length-stratified sample rather
#'   than from independent surveys.
#' @param ln_theta_agg Log overdispersion used when \code{Comp_Type = 0}.
#' @param Comp_Type Integer specifying the composition parameterization, as in
#'   \code{\link{Get_Comp_Likelihoods}} (0 aggregated, 1 split by region and sex,
#'   2 joint across sexes and split by region).
#' @param Likelihood_Type Integer specifying the likelihood family. Only
#'   \code{0} (multinomial) and \code{1} (Dirichlet-multinomial) are supported.
#' @param n_regions Number of regions modeled.
#' @param n_lens Number of length bins.
#' @param n_model_bins Number of age bins used internally in the model.
#' @param n_obs_bins Number of observed age bins.
#' @param n_sexes Number of sexes modeled.
#' @param AgeingError Ageing error matrix mapping model ages to observed ages.
#' @param use Integer matrix \eqn{[region \times len]} indicating which cells
#'   have observations (\code{1} = use, \code{0} = ignore). A length bin with no
#'   aged fish in any region is skipped entirely.
#' @param addtocomp Small constant added to compositions to avoid numerical
#'   issues when zeros are present.
#' @param comp_bins Integer vector of age bins the composition is fitted over, or
#'   \code{NULL} (default) for all of them.
#' @param comp_const_obs Integer (0 or 1). Whether \code{addtocomp} is added to
#'   the observed proportions that weight the multinomial.
#'
#' @return Array \eqn{[region \times len \times sex]} of negative
#'   log-likelihood contributions. Length bins with no observations stay at zero.
#'
#' @keywords internal
#' @import RTMB
Get_CAAL_Likelihoods = function(Exp,
                                Obs,
                                ISS,
                                Wt_Mltnml,
                                ln_theta,
                                ln_theta_agg,
                                Comp_Type,
                                Likelihood_Type,
                                n_regions,
                                n_lens,
                                n_model_bins,
                                n_obs_bins,
                                n_sexes,
                                AgeingError,
                                use,
                                addtocomp,
                                comp_bins = NULL,
                                comp_const_obs = 1
                                ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  caal_nLL = array(0, dim = c(n_regions, n_lens, n_sexes)) # initialize nLL here

  # Coerce to full shape
  Exp = array(Exp, dim = c(n_regions, n_lens, n_model_bins, n_sexes))
  Obs = array(Obs, dim = c(n_regions, n_lens, n_obs_bins, n_sexes))
  ISS = array(ISS, dim = c(n_regions, n_lens, n_sexes))
  Wt_Mltnml = array(Wt_Mltnml, dim = c(n_regions, n_lens, n_sexes))
  use = array(use, dim = c(n_regions, n_lens))

  for(l in 1:n_lens) {

    # a length bin with no aged fish anywhere contributes nothing and must be
    # skipped, since the composition machinery normalizes by a row sum - could cause an Inf
    if(sum(use[,l]) < 1) next

    caal_nLL[,l,] = Get_Comp_Likelihoods(
      Exp = array(Exp[,l,,], dim = c(n_regions, n_model_bins, n_sexes)),
      Obs = array(Obs[,l,,], dim = c(n_regions, n_obs_bins, n_sexes)),
      ISS = array(ISS[,l,], dim = c(n_regions, n_sexes)),
      Wt_Mltnml = array(Wt_Mltnml[,l,], dim = c(n_regions, n_sexes)),
      ln_theta = ln_theta,
      ln_theta_agg = ln_theta_agg,
      Comp_Type = Comp_Type,
      Likelihood_Type = Likelihood_Type,
      n_regions = n_regions,
      n_model_bins = n_model_bins,
      n_obs_bins = n_obs_bins,
      n_sexes = n_sexes,
      age_or_len = 0, # a CAAL row is an age composition, so ageing error applies
      AgeingError = AgeingError,
      use = use[,l],
      addtocomp = addtocomp,
      comp_bins = comp_bins,
      comp_const_obs = comp_const_obs
    )

  } # end l loop

  return(caal_nLL)
}


#' Pack observed CAAL data into a single flat OBS vector (OSA)
#'
#' Produces the flat tracked OBS vector required by \code{RTMB::oneStepPredict}
#' for conditional age-at-length data. Mirrors \code{\link{pack_comp_osa}} with a
#' length-bin loop inserted inside the season loop, so each length bin becomes its
#' own tracked group of age bins. Only the discrete families are handled, since
#' those are the only ones the CAAL likelihood supports.
#'
#' Counts are formed the same way as for the marginal compositions:
#' \itemize{
#'   \item Multinomial (0): counts = round(prop x ISS x Wt)
#'   \item Dirichlet-multinomial (1): counts = round(prop x ISS)
#' }
#'
#' The order is year, fleet, season, length bin, then region-fastest within the
#' group, and \code{\link{eval_caal_osa}} walks the vector in exactly that order.
#'
#' @param ObsArr Observed CAAL array \eqn{[region \times year \times season
#'   \times len \times age \times sex \times fleet]}.
#' @param ISSArr Input sample sizes \eqn{[region \times year \times season \times
#'   len \times sex \times fleet]}.
#' @param WtArr Multinomial weights, same shape as \code{ISSArr}.
#' @param UseArr Use flags \eqn{[region \times year \times season \times len
#'   \times fleet]}.
#' @param TypeMat Composition type matrix \eqn{[year \times fleet]}.
#' @param LikeTypeVec Likelihood type per fleet.
#' @param n_yrs Number of model years.
#' @param n_seas Number of seasons per year.
#' @param n_lens Number of length bins.
#' @param n_fleets Number of fleets.
#' @param n_sexes Number of sexes.
#' @param addtocomp Small constant added to proportions before normalization.
#' @param return_labels Logical; if \code{TRUE}, also builds a per-element label
#'   data.frame identifying the origin of every entry, in the same order, for
#'   post-hoc relabeling of \code{TMB::oneStepPredict()} residuals.
#'
#' @return If \code{return_labels = FALSE} (default), the flat OBS vector, or
#'   \code{NULL} when no fleet carries CAAL data. If \code{return_labels = TRUE},
#'   a list with \code{vec} and \code{labels}.
#'
#' @keywords internal
#' @import RTMB
pack_caal_osa = function(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                         n_yrs, n_seas, n_lens, n_fleets, n_sexes, addtocomp,
                         return_labels = FALSE) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_obs_bins = dim(ObsArr)[5]

  # the label rows for one accepted (y, seas, f, l) group, in the same element
  # order as the packed values
  make_labels = function(used, n_ru, ct, lt, y, seas, f, l) {
    if(ct == 0) {
      region = rep(used[1], n_obs_bins); sex = rep(1L, n_obs_bins); bin = 1:n_obs_bins
      last_in_group = (bin == n_obs_bins)
    } else {
      region = rep(used, times = n_obs_bins * n_sexes)
      bin = rep(rep(1:n_obs_bins, each = n_ru), times = n_sexes)
      sex = rep(1:n_sexes, each = n_ru * n_obs_bins)
      # split by sex: each (region, sex) is its own multinomial, so each has its own
      # determined bin. Joint by sex: the whole stack per region is one multinomial,
      # so only the single last cell is determined.
      last_in_group = if(ct == 1) (bin == n_obs_bins) else (bin == n_obs_bins) & (sex == n_sexes)
    }
    data.frame(region = region, year = y, season = seas, len = l, fleet = f,
               sex = sex, bin = bin, comp_type = ct, likelihood_type = lt,
               family = "discrete", last_in_group = last_in_group)
  }

  clean = list()
  label_rows = list()

  for(y in 1:n_yrs) {
    for(f in 1:n_fleets) {
      if(!LikeTypeVec[f] %in% c(0,1)) next # CAAL supports the discrete families only
      for(seas in 1:n_seas) {
        for(l in 1:n_lens) {

          use_vec = UseArr[, y, seas, l, f]
          if(sum(use_vec) < 1) next

          ct = TypeMat[y, f]
          lt = LikeTypeVec[f]
          used = which(use_vec == 1)
          n_ru = length(used)

          obs_slice = ObsArr[used, y, seas, l, , , f, drop = FALSE]
          dim(obs_slice) = c(n_ru, n_obs_bins, n_sexes)

          if(ct == 2) {
            # Joint by sex: the true sample is one joint draw across sexes
            for(rr in 1:n_ru) {
              r_orig = used[rr]
              iss = ISSArr[r_orig, y, seas, l, 1, f]
              wt = WtArr[r_orig, y, seas, l, 1, f]
              v = as.vector(obs_slice[rr, , ]) # bin-fastest-then-sex
              pr = (v + addtocomp) / sum(v + addtocomp)
              if(lt == 0) v = round(pr * iss * wt) # multinomial
              if(lt == 1) v = round(pr * iss) # DM (no Wt)
              for(s in 1:n_sexes) obs_slice[rr, , s] = v[((s - 1) * n_obs_bins + 1):(s * n_obs_bins)]
            } # end rr loop
          } else {
            # Aggregated (ct 0) or split by region and sex (ct 1): each region and sex
            # cell is its own independent sample
            for(rr in 1:n_ru) {
              r_orig = used[rr]
              for(s in 1:n_sexes) {
                iss = ISSArr[r_orig, y, seas, l, s, f]
                wt = WtArr[r_orig, y, seas, l, s, f]
                pr = (obs_slice[rr, , s] + addtocomp) / sum(obs_slice[rr, , s] + addtocomp)
                if(lt == 0) obs_slice[rr, , s] = round(pr * iss * wt) # multinomial
                if(lt == 1) obs_slice[rr, , s] = round(pr * iss) # DM (no Wt)
              } # end s loop
            } # end rr loop
          } # end composition type

          g = if(ct == 0) as.vector(obs_slice[1, , 1]) else as.vector(obs_slice)

          clean[[length(clean) + 1]] = g
          if(return_labels) label_rows[[length(label_rows) + 1]] = make_labels(used, n_ru, ct, lt, y, seas, f, l)

        } # end l loop
      } # end seas loop
    } # end f loop
  } # end y loop

  if(length(clean) == 0) return(NULL)

  vec = unlist(clean)
  if(!return_labels) return(vec)
  return(list(vec = vec, labels = do.call(rbind, label_rows)))
}


#' Evaluate CAAL likelihoods from a tracked OBS vector (OSA)
#'
#' Walks the flat vector built by \code{\link{pack_caal_osa}} in the same order
#' and evaluates each length bin's age composition through
#' \code{\link{Get_Comp_Likelihoods_OSA}}, so the observations stay on the tape as
#' the tracked quantity \code{RTMB::oneStepPredict} needs.
#'
#' @param nLL_arr Array \eqn{[region \times year \times season \times len \times
#'   sex \times fleet]} receiving negative log-likelihood contributions.
#' @param tracked Flat tracked OBS vector from \code{\link{pack_caal_osa}}.
#' @param ExpArrFn Function of \code{(y, seas, l, f)} returning the expected joint
#'   numbers at length and age for that group, indexed
#'   \eqn{[region \times age \times sex]}.
#' @param UseArr Use flags \eqn{[region \times year \times season \times len
#'   \times fleet]}.
#' @param TypeMat Composition type matrix \eqn{[year \times fleet]}.
#' @param LikeTypeVec Likelihood type per fleet.
#' @param ISSArr Input sample sizes \eqn{[region \times year \times season \times
#'   len \times sex \times fleet]}.
#' @param lnThetaArr Log overdispersion \eqn{[region \times sex \times fleet]}.
#' @param lnThetaAggVec Aggregated log overdispersion, one per fleet.
#' @param n_regions,n_yrs,n_seas,n_lens,n_fleets,n_sexes Dimension sizes.
#' @param n_model_bins Number of model age bins.
#' @param n_obs_bins Number of observed age bins.
#' @param AgeingErrorFn Function returning the ageing error matrix for a year.
#' @param addtocomp Small constant added to proportions before normalization.
#'
#' @return Updated \code{nLL_arr}.
#'
#' @keywords internal
#' @import RTMB
eval_caal_osa = function(nLL_arr, tracked, ExpArrFn,
                         UseArr, TypeMat, LikeTypeVec,
                         ISSArr, lnThetaArr, lnThetaAggVec,
                         n_regions, n_yrs, n_seas, n_lens, n_fleets, n_sexes,
                         n_model_bins, n_obs_bins,
                         AgeingErrorFn, addtocomp) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(is.null(tracked)) return(nLL_arr)

  d = dim(nLL_arr)
  nLL_arr = RTMB::AD(as.vector(nLL_arr) * 0)
  dim(nLL_arr) = d

  k = 1
  for(y in 1:n_yrs) {
    for(f in 1:n_fleets) {
      if(!LikeTypeVec[f] %in% c(0,1)) next
      for(seas in 1:n_seas) {
        for(l in 1:n_lens) {

          use_vec = UseArr[, y, seas, l, f]
          if(sum(use_vec) < 1) next

          ct = TypeMat[y, f]
          lt = LikeTypeVec[f]
          n_ru = sum(use_vec == 1)

          slice_length = if(ct == 0) n_obs_bins else n_ru * n_obs_bins * n_sexes
          active_obs_slice = tracked[k:(k + slice_length - 1)]

          nLL_arr[, y, seas, l, , f] = Get_Comp_Likelihoods_OSA(
            Exp = ExpArrFn(y, seas, l, f), Obs = active_obs_slice,
            ISS = array(ISSArr[, y, seas, l, , f], dim = c(n_regions, n_sexes)),
            ln_theta = lnThetaArr[, , f], ln_theta_agg = lnThetaAggVec[f],
            Comp_Type = ct, Likelihood_Type = lt,
            n_regions = n_regions, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
            n_sexes = n_sexes, age_or_len = 0,
            AgeingError = if(is.null(AgeingErrorFn)) NA else AgeingErrorFn(y),
            use = use_vec, addtocomp = addtocomp
          )

          k = k + slice_length

        } # end l loop
      } # end seas loop
    } # end f loop
  } # end y loop

  return(nLL_arr)
}


#' Sum a conditional age-at-length array across populations
#'
#' The conditional age-at-length arrays carry a population dimension the
#' likelihood does not use, so it is summed away before the comparison. A single
#' population needs only a reshape, which avoids an apply over a degenerate
#' margin.
#'
#' @param arr Array indexed population, region, year, season, length, age, sex, fleet.
#' @param y,seas,f Year, season and fleet to extract.
#' @param n_pop,n_regions,n_lens,n_ages,n_sexes Model dimensions.
#'
#' @return An array indexed region, length, age, sex.
#'
#' @keywords internal
caal_sum_pop = function(arr, y, seas, f, n_pop, n_regions, n_lens, n_ages, n_sexes) {
  extracted = arr[,,y,seas,,,,f, drop = FALSE]
  if(n_pop == 1) dim(extracted) = c(n_regions, n_lens, n_ages, n_sexes)
  else extracted = apply(extracted, c(2,5,6,7), sum)
  return(extracted)
}

#' Sum a conditional age-at-length array across populations, for one length bin
#'
#' As \code{\link{caal_sum_pop}}, for a single length bin.
#'
#' @param arr Array indexed population, region, year, season, length, age, sex, fleet.
#' @param y,seas,l,f Year, season, length bin and fleet to extract.
#' @param n_pop,n_regions,n_ages,n_sexes Model dimensions.
#'
#' @return An array indexed region, age, sex.
#'
#' @keywords internal
caal_sum_pop_len = function(arr, y, seas, l, f, n_pop, n_regions, n_ages, n_sexes) {
  extracted = arr[,,y,seas,l,,,f, drop = FALSE]
  if(n_pop == 1) dim(extracted) = c(n_regions, n_ages, n_sexes)
  else extracted = apply(extracted, c(2,6,7), sum)
  return(extracted)
}
