# At-age observation likelihoods, shared by eight streams: retained catch,
# discards, fishery index and survey index, each aggregated and
# population-specific. Streams differ only in the prediction array they read and
# the parameter supplying their observation error.
#
# Every stream is stored region by year by season by age by sex by fleet, with a
# leading population dimension for the population-specific form, so an
# observation can be split or summed over regions and sexes independently. The
# margins a fleet sums over are named by its Type code, and the observation
# occupies slot one of every summed margin.

#' Decode an at-age aggregation type into its split margins
#'
#' The Type codes follow the composition vocabulary: \code{"agg"} sums over both
#' regions and sexes, \code{"spltRaggS"} keeps regions apart and sums over sexes,
#' \code{"aggRspltS"} does the reverse, and \code{"spltRspltS"} keeps both apart.
#'
#' @param code Integer, \code{0} to \code{3} in the order above.
#'
#' @return A list with logical \code{region} and \code{sex}, \code{TRUE} where
#'   that margin is split.
#'
#' @keywords internal
at_age_split = function(code) {
  list(region = code %in% c(1, 3), sex = code %in% c(2, 3))
}

#' Standard deviation for one at-age observation
#'
#' An at-age observation may carry its own reported standard error, an estimated
#' component, or both, matching what the aggregated index streams allow. The
#' parameter alone is the default and is what a stream with no reported errors
#' means.
#'
#' @param se Reported standard errors for the ages in one cell.
#' @param extra Estimated component for the same ages, on the natural scale.
#' @param form Integer. \code{0} the parameter alone, \code{1} the reported
#'   errors alone, \code{2} additive, \code{3} in quadrature.
#'
#' @return A vector of standard deviations the length of \code{extra}.
#'
#' @keywords internal
at_age_obs_sd = function(se, extra, form) {
  if(form == 1) return(se)                    # reported error alone
  if(form == 2) return(se + extra)            # additive
  if(form == 3) return(sqrt(se^2 + extra^2))  # independent variances
  return(extra)                               # the parameter alone
}

#' Transform at-age observations onto the scale their likelihood is written on
#'
#' A lognormal stream is fit on the log scale and a normal stream on the natural
#' scale, and the choice is per fleet, so the transformation is applied cell by
#' cell before \code{\link[RTMB]{OBS}} registration. Registration must happen
#' against the name \code{getAll} supplied, so the caller does it: a vector
#' registered under a local name does not link to the data element, and the
#' objective then diverges from the reported likelihood.
#'
#' @param obs Observation array, shaped like \code{use}.
#' @param use Integer array flagging which cells are fit.
#' @param like_type Integer per fleet. \code{0} lognormal, \code{1} normal.
#' @param const Small constant added inside the log of a lognormal cell.
#'
#' @return A numeric vector, one element per flagged cell in \code{which()}
#'   order, on the scale its fleet's likelihood uses.
#'
#' @keywords internal
prep_at_age_obs = function(obs, use, like_type, const = 0) {

  fit_cells = which(use == 1)
  if(length(fit_cells) == 0) return(numeric(0))

  d = dim(use)
  cells_per_fleet = prod(d[-length(d)])
  fleet = 1L + (fit_cells - 1L) %/% cells_per_fleet   # fleet is the last dimension

  x = as.numeric(obs)[fit_cells]
  lognormal = like_type[fleet] == 0
  x[lognormal] = log(x[lognormal] + const)

  return(x)
}

#' Predicted value for one age-disaggregated observation
#'
#' Catch and discards are already at age and only need their units applied. The
#' discards are dead discards, so they are raised by the discard mortality rate
#' to the total the observation counts, exactly as the aggregated stream does.
#' The indices apply an age-specific catchability to the numbers available to
#' that fleet.
#'
#' Population, region and sex arrive as index vectors rather than single
#' indices. A margin the fleet splits over is a single index, and a margin it
#' sums over is the whole extent, so one expression covers every aggregation.
#'
#' @param source Character, one of \code{"catch"}, \code{"discard"},
#'   \code{"fish_index"} or \code{"srv_index"}.
#' @param arrays Named list of the model arrays the prediction reads:
#'   \code{CAA}, \code{DAA}, \code{SrvIAA}, \code{FishIAA}, \code{WAA_fish},
#'   \code{dmr}, \code{catch_units} and \code{discard_units}. The two index
#'   arrays already carry their fleet's selectivity, timing and movement
#'   treatment, and the age shape of catchability lives in that selectivity: a
#'   fleet fit age by age uses the \code{"nonparfree"} selectivity form, whose
#'   values carry the height of the curve as well as its shape.
#' @param p_idx,r_idx,s_idx Population, region and sex indices, each either one
#'   index or the whole extent of that margin.
#' @param y,seas,a,f Year, season, age and fleet indices.
#'
#' @return The predicted observation, a scalar.
#'
#' @keywords internal
get_at_age_prediction = function(source, arrays, p_idx, r_idx, s_idx, y, seas, a, f) {

  if(source == "catch") {
    numbers = arrays$CAA[p_idx,r_idx,y,seas,a,s_idx,f]
    if(arrays$catch_units[f] == 0) return(sum(numbers))                            # abundance
    return(sum(numbers * arrays$WAA_fish[p_idx,r_idx,y,seas,a,s_idx,f]))           # biomass
  }

  if(source == "discard") { # dead discards raised to the total discarded
    total = 0
    for(rr in r_idx) { # the mortality rate is region specific, so raise region by region
      numbers = arrays$DAA[p_idx,rr,y,seas,a,s_idx,f] / arrays$dmr[rr,y,seas,f]
      if(arrays$discard_units[f] == 0) total = total + sum(numbers)                # abundance
      else total = total + sum(numbers * arrays$WAA_fish[p_idx,rr,y,seas,a,s_idx,f]) # biomass
    } # end rr loop
    return(total)
  }

  if(source == "fish_index") { # index numbers at age, selectivity already applied
    return(sum(arrays$FishIAA[p_idx,r_idx,y,seas,a,s_idx,f]))
  }

  return(sum(arrays$SrvIAA[p_idx,r_idx,y,seas,a,s_idx,f])) # survey available numbers
}

#' Evaluate one age-disaggregated observation stream
#'
#' Computes the at-age negative log likelihood for every fleet in one stream.
#' Observations arrive already transformed by \code{\link{prep_at_age_obs}} and
#' registered through \code{\link[RTMB]{OBS}}.
#'
#' Everything that can differ between fleets does: the margins summed over, the
#' error structure, whether reported standard errors enter, and whether the
#' density is lognormal or normal. Ages within a cell may be independent, an
#' AR(1) across ages, or an unstructured correlation matrix; a fleet may instead
#' correlate over both age and year through a separable AR(1), which needs the
#' age by year block it is given to be complete.
#'
#' @param obs_t Registered observations for this stream, one element per cell
#'   flagged in \code{use}, in \code{which()} order, on the scale its fleet's
#'   likelihood uses.
#' @param use Integer array flagging which cells are fit, dimensioned region by
#'   year by season by age by sex by fleet, with a leading population dimension
#'   when \code{pop} is \code{TRUE}.
#' @param ln_sigma Log-scale observation error, over age by sex by fleet, with a
#'   leading population dimension when \code{pop} is \code{TRUE}.
#' @param source,arrays Passed to \code{\link{get_at_age_prediction}}.
#' @param pop Logical. \code{TRUE} for the population-specific stream, whose
#'   arrays carry a leading population dimension and whose observations are
#'   never summed over populations.
#' @param obs_se Reported standard errors shaped like \code{use}, read only by
#'   fleets whose \code{sd_form} asks for them.
#' @param sd_form Integer per fleet, see \code{\link{at_age_obs_sd}}.
#' @param like_type Integer per fleet. \code{0} lognormal, \code{1} normal.
#' @param const Small constant added inside the log of a lognormal cell,
#'   matching the aggregated stream's convention.
#' @param corr_type Integer per fleet. \code{0} \code{"iid"}, \code{1}
#'   \code{"1dar1"}, \code{2} \code{"us"}, \code{3} \code{"2dar1"}.
#' @param trans_rho Unconstrained correlation across ages, over sex by fleet.
#' @param trans_rho_year Unconstrained correlation across years, over sex by
#'   fleet, read under \code{"2dar1"}.
#' @param us_pars Unconstrained correlation parameters, over pair by sex by
#'   fleet, read under \code{"us"}.
#' @param aa_type Integer per fleet naming the split margins, see
#'   \code{\link{at_age_split}}.
#'
#' @return A list with \code{nLL} and \code{pred}, both arrays shaped like
#'   \code{use} and zero wherever nothing is fit. The predictions are returned so
#'   they can be reported and plotted directly rather than reconstructed.
#'
#' @keywords internal
get_at_age_stream_nLL = function(obs_t, use, ln_sigma, source, pop, arrays,
                                 obs_se = NULL, sd_form = 0, like_type = 0,
                                 const = 0, corr_type = 0, trans_rho = 0, trans_rho_year = 0,
                                 us_pars = NULL, aa_type = 1) {

  "[<-" <- RTMB::ADoverload("[<-")

  d = dim(use)
  stream_nLL = array(0, dim = d)  # zero wherever nothing is fit
  stream_pred = array(0, dim = d)
  if(!any(use == 1)) return(list(nLL = stream_nLL, pred = stream_pred))

  nd = length(d)
  n_fleets = d[nd]
  # every per-fleet setting is read by fleet, so a single setting stands for all
  sd_form = rep_len(sd_form, n_fleets)
  like_type = rep_len(like_type, n_fleets)
  corr_type = rep_len(corr_type, n_fleets)
  aa_type = rep_len(aa_type, n_fleets)

  i_r = if(pop) 2 else 1        # dimension positions within one fleet's slice
  i_y = i_r + 1; i_seas = i_y + 1; i_a = i_seas + 1; i_s = i_a + 1

  # the prediction array supplies the full extent of every margin a fleet sums
  # over, so an aggregated margin reads as the whole dimension
  pred_arr = switch(source, catch = arrays$CAA, discard = arrays$DAA,
                    fish_index = arrays$FishIAA, arrays$SrvIAA)
  all_pop = seq_len(dim(pred_arr)[1])
  all_reg = seq_len(dim(pred_arr)[2])
  all_sex = seq_len(dim(pred_arr)[6])

  fit_cells = which(use == 1)               # linear positions of fitted observations
  obs_slot = array(NA_integer_, dim = d)
  obs_slot[fit_cells] = seq_along(fit_cells) # array position -> slot in obs_t

  df = d[-nd]                               # one fleet's slice, fleet being last
  n_cell = prod(df)
  strides = c(1, cumprod(df)[-length(df)])  # linear offsets within that slice

  for(f in seq_len(n_fleets)) {

    fleet_off = (f - 1) * n_cell
    use_f = array(use[fleet_off + seq_len(n_cell)], dim = df)
    if(!any(use_f == 1)) next

    split = at_age_split(aa_type[f])
    r_all = if(split$region) NULL else all_reg  # NULL means "take the cell's own index"
    s_all = if(split$sex) NULL else all_sex

    # an unstructured correlation is one matrix per sex and fleet, built once and
    # subset to whichever ages a cell observes
    us_corr = NULL
    if(corr_type[f] == 2) {
      us_corr = vector("list", df[i_s])
      for(s in seq_len(df[i_s])) us_corr[[s]] = build_us_corr(us_pars[,s,f], df[i_a])
    } # end s loop

    # a separable correlation runs over years as well as ages, so year leaves the
    # cell definition and the block of years by ages is evaluated at once
    free_dims = if(corr_type[f] == 3) c(i_y, i_a) else i_a
    margins = setdiff(seq_along(df), free_dims)
    active_cells = which(apply(use_f == 1, margins, any))
    cell_index = arrayInd(active_cells, df[margins])

    for(cell in seq_len(nrow(cell_index))) {

      idx = integer(length(df))
      idx[margins] = cell_index[cell,]
      p = if(pop) idx[1] else 1
      s = idx[i_s]
      p_idx = if(pop) p else all_pop
      r_idx = if(is.null(r_all)) idx[i_r] else r_all
      s_idx = if(is.null(s_all)) s else s_all

      # linear position of this cell with every free dimension at its first slot
      base = 1 + sum((pmax(idx, 1) - 1) * strides)
      age_step = (seq_len(df[i_a]) - 1) * strides[i_a]

      if(corr_type[f] == 3) {

        yr_step = (seq_len(df[i_y]) - 1) * strides[i_y]
        block = base::matrix(use_f[as.vector(base + outer(yr_step, age_step, "+"))], nrow = df[i_y])
        obs_yrs = which(rowSums(block) > 0)
        obs_ages = which(colSums(block) > 0)
        if(!all(block[obs_yrs, obs_ages] == 1)) {
          stop("A fleet fitting at-age observations as '2dar1' must observe a complete ",
               "block of ages by years, since a separable correlation is defined over the ",
               "whole grid. Fleet ", f, " has gaps in that block. Use '1dar1' or 'us', ",
               "which are defined over whatever ages a cell observes.")
        }

        lin = as.vector(fleet_off + base + outer(yr_step[obs_yrs], age_step[obs_ages], "+"))
        slot = obs_slot[lin]
        extra = if(pop) exp(ln_sigma[p,obs_ages,s,f]) else exp(ln_sigma[obs_ages,s,f])

        pred = rep(0, length(obs_yrs) * length(obs_ages))
        k = 1
        for(ay in seq_along(obs_ages)) {   # column major, matching the slot matrix
          for(yy in seq_along(obs_yrs)) {
            pred[k] = get_at_age_prediction(source, arrays, p_idx, r_idx, s_idx,
                                            obs_yrs[yy], idx[i_seas], obs_ages[ay], f)
            k = k + 1
          } # end yy loop
        } # end ay loop

        pred_t = if(like_type[f] == 0) log(pred + const) else pred
        sd_vec = extra[rep(seq_along(obs_ages), each = length(obs_yrs))]
        if(sd_form[f] != 0) sd_vec = at_age_obs_sd(as.numeric(obs_se)[lin], sd_vec, sd_form[f])

        resid = matrix(obs_t[slot] - pred_t, nrow = length(obs_yrs))
        scale = matrix(sd_vec, nrow = length(obs_yrs))
        cell_nLL = rep(0, length(pred))
        cell_nLL[1] = get_at_age_2dar1_nLL(resid, scale, trans_rho[s,f], trans_rho_year[s,f])

        stream_nLL[lin] = cell_nLL
        stream_pred[lin] = pred
        next
      } # end 2dar1

      obs_ages = which(use_f[base + age_step] == 1)
      if(length(obs_ages) == 0) next

      lin = fleet_off + base + age_step[obs_ages]
      # subset the registered vector directly; copying element by element loses the OBS tagging
      slot = obs_slot[lin]
      extra = if(pop) exp(ln_sigma[p,obs_ages,s,f]) else exp(ln_sigma[obs_ages,s,f])
      sd_vec = if(sd_form[f] == 0) extra else at_age_obs_sd(as.numeric(obs_se)[lin], extra, sd_form[f])

      pred = rep(0, length(obs_ages))
      for(age_slot in seq_along(obs_ages)) { # prediction for each observed age
        pred[age_slot] = get_at_age_prediction(source, arrays, p_idx, r_idx, s_idx,
                                               idx[i_y], idx[i_seas], obs_ages[age_slot], f)
      } # end age_slot loop

      pred_t = if(like_type[f] == 0) log(pred + const) else pred

      # iid returns one value per age; a correlated cell puts its density on the first
      cell_nLL = get_at_age_nLL(obs_t[slot], pred_t, sd_vec, corr_type[f], rho_trans(trans_rho[s,f]),
                                ages = obs_ages,
                                corr_mat = if(corr_type[f] == 2) us_corr[[s]][obs_ages,obs_ages] else NULL)

      stream_nLL[lin] = cell_nLL
      stream_pred[lin] = pred
    } # end cell loop
  } # end f loop

  return(list(nLL = stream_nLL, pred = stream_pred))
}
