# At-age observation likelihoods for six data sources: retained catch, discards and survey index, each
# aggregated and population-specific. Stored region x year x season x observed age x sex x fleet.

#' Decode an at-age aggregation type into its split dims
#'
#' The Type codes follow the composition vocabulary: \code{"agg"} sums over both
#' regions and sexes, \code{"spltRaggS"} keeps regions apart and sums over sexes,
#' \code{"aggRspltS"} does the reverse, and \code{"spltRspltS"} keeps both apart.
#'
#' @param code Integer, \code{0} to \code{3} in the order above.
#'
#' @return A list with logical \code{region} and \code{sex}, \code{TRUE} where
#'   that dim is split.
#'
#' @keywords internal
at_age_split = function(code) {
  list(region = code %in% c(1, 3), sex = code %in% c(2, 3))
}

#' Standard deviation for one at-age observation
#'
#' An at-age observation may have its own reported standard error, an estimated
#' component, or both, matching what the aggregated index data sources allow. The
#' parameter alone is the default and is what a data source with no reported errors
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

# At-Age Observation Helpers ------------------------------------------------

#' Transform at-age observations onto the scale their likelihood is written on
#'
#' A lognormal data source is fit on the log scale and a normal data source on the natural
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
  fleet = 1 + (fit_cells - 1) %/% cells_per_fleet   # fleet is the last dimension

  x = as.numeric(obs)[fit_cells]
  lognormal = like_type[fleet] == 0
  x[lognormal] = log(x[lognormal] + const)

  return(x)
}

#' Predicted value for one age-disaggregated observation
#'
#' Catch and discards are already at age and only need their units applied. The
#' discards are dead discards, so they are raised by the discard mortality rate
#' to the total the observation counts, exactly as the aggregated data source does.
#' The survey index applies an age-specific catchability to the numbers
#' available to that fleet.
#'
#' Population, region and sex arrive as index vectors rather than single
#' indices. A dim the fleet splits over is a single index, and a dim it
#' sums over is the whole extent, so one expression covers every aggregation.
#'
#' @param source Character, one of \code{"catch"}, \code{"discard"} or
#'   \code{"srv_index"}.
#' @param arrays Named list of the model arrays the prediction reads:
#'   \code{CAA}, \code{DAA}, \code{SrvIAA}, \code{WAA_fish},
#'   \code{dmr}, \code{catch_units} and \code{discard_units}. The two index
#'   arrays already have their fleet's selectivity, timing and movement
#'   treatment, and the age shape of catchability lives in that selectivity: a
#'   fleet fit age by age uses the \code{"nonparfree"} selectivity form, whose
#'   values hold the height of the curve as well as its shape.
#' @param p_idx,r_idx,s_idx Population, region and sex indices, each either one
#'   index or the whole extent of that dim.
#' @param y,seas,a,f Year, season, model age and fleet indices. \code{seas} is one
#'   season, or every season of the year for a data source reported as a season
#'   total.
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

  return(sum(arrays$SrvIAA[p_idx,r_idx,y,seas,a,s_idx,f])) # survey available numbers
}

#' Is one fleet's ageing error the identity in every year?
#'
#' An identity matrix reads every age as itself, so a fleet with one is
#' predicted on the model's own ages without passing through the map.
#'
#' @param ageing_error Array \code{[n_years, n_ages, n_obs_ages]} for one fleet.
#'
#' @return \code{TRUE} when every year is the identity matrix.
#'
#' @keywords internal
is_identity_ageing_error = function(ageing_error) {

  d = dim(ageing_error)
  if(d[2] != d[3]) return(FALSE)

  identity_mat = base::diag(d[2])
  for(y in seq_len(d[1])) {
    if(!all(ageing_error[y,,] == identity_mat)) return(FALSE)
  } # end y loop

  return(TRUE)
}

#' Predicted values at the ages the observations are recorded on
#'
#' An at-age observation counts fish by the age they were read as, so the prediction at each
#' observed age sums the predictions at every model age read as it, weighted by the fleet's
#' ageing error matrix, with units and discard mortality applied at the model age first.
#'
#' @param source,arrays,p_idx,r_idx,s_idx,y,seas,f See
#'   \code{\link{get_at_age_prediction}}.
#' @param obs_ages Integer vector of the observed ages to predict.
#' @param ageing_error Array \code{[n_years, n_ages, n_obs_ages]} for this
#'   fleet, or \code{NULL} when the observed ages are the model ages.
#'
#' @return A vector the length of \code{obs_ages}.
#'
#' @keywords internal
get_at_age_obs_prediction = function(source, arrays, p_idx, r_idx, s_idx, y, seas, obs_ages, f, ageing_error = NULL) {

  "[<-" <- RTMB::ADoverload("[<-")

  pred = rep(0, length(obs_ages))

  # observed ages are the model ages
  if(is.null(ageing_error)) {
    for(k in seq_along(obs_ages)) {
      pred[k] = get_at_age_prediction(source, arrays, p_idx, r_idx, s_idx, y, seas, obs_ages[k], f)
    } # end k loop
    return(pred)
  }

  ae = base::matrix(ageing_error[y,,obs_ages], ncol = length(obs_ages)) # model age by observed age
  from_ages = which(base::rowSums(ae != 0) > 0) # model ages read as any of these observed ages

  # prediction at each of those model ages
  pred_model = rep(0, length(from_ages))
  for(j in seq_along(from_ages)) {
    pred_model[j] = get_at_age_prediction(source, arrays, p_idx, r_idx, s_idx, y, seas, from_ages[j], f)
  } # end j loop

  # each observed age collects the model ages read as it
  for(k in seq_along(obs_ages)) {
    read_as = which(ae[from_ages,k] != 0)
    pred[k] = sum(pred_model[read_as] * ae[from_ages[read_as],k])
  } # end k loop

  return(pred)
}

#' Evaluate age-disaggregated data source
#'
#' Computes the at-age negative log likelihood for every fleet in one data source.
#' Observations arrive already transformed by \code{\link{prep_at_age_obs}} and
#' registered through \code{\link[RTMB]{OBS}}.
#'
#' Everything that can differ between fleets does: the dims summed over, the
#' error structure, whether reported standard errors enter, and whether the
#' density is lognormal or normal. Ages within a cell may be independent, an
#' AR(1) across ages, or an unstructured correlation matrix; a fleet may instead
#' correlate over both age and year through a separable AR(1), which needs the
#' age by year block it is given to be complete.
#'
#' @param obs_t Registered observations for this data source, one element per cell
#'   flagged in \code{use}, in \code{which()} order, on the scale its fleet's
#'   likelihood uses.
#' @param use Integer array flagging which cells are fit, dimensioned region by
#'   year by season by observed age by sex by fleet, with a leading population
#'   dimension when \code{pop} is \code{TRUE}.
#' @param ln_sigma Log-scale observation error, over observed age by sex by
#'   fleet, with a leading population dimension when \code{pop} is \code{TRUE}.
#' @param source,arrays Passed to \code{\link{get_at_age_prediction}}.
#' @param seas_agg Integer vector, one per fleet. \code{1} compares the
#'   observation against every season of the year summed together, \code{0}
#'   against the season it sits in.
#' @param pop Logical. \code{TRUE} for the population-specific data source, whose
#'   arrays have a leading population dimension and whose observations are
#'   never summed over populations.
#' @param obs_se Reported standard errors shaped like \code{use}, read only by
#'   fleets whose \code{sd_form} asks for them.
#' @param sd_form Integer per fleet, see \code{\link{at_age_obs_sd}}.
#' @param like_type Integer per fleet. \code{0} lognormal, \code{1} normal.
#' @param const Small constant added inside the log of a lognormal cell,
#'   matching the aggregated data source's convention.
#' @param corr_type Integer per fleet. \code{0} \code{"iid"}, \code{1}
#'   \code{"1dar1"}, \code{2} \code{"us"}, \code{3} \code{"2dar1"}.
#' @param trans_rho Unconstrained correlation across ages, over region by sex by
#'   fleet, with a leading population dim when \code{pop} is \code{TRUE}.
#' @param trans_rho_year Unconstrained correlation across years, shaped like
#'   \code{trans_rho}, read under \code{"2dar1"}.
#' @param us_pars Unconstrained correlation parameters, over pair by the dims
#'   of \code{trans_rho}, read under \code{"us"}.
#' @param aa_type Integer codes naming the split dims, as a matrix over year
#'   by fleet or a vector per fleet standing for every year, see
#'   \code{\link{at_age_split}}.
#' @param ageing_error Array \code{[n_years, n_ages, n_obs_ages, n_fleets]}
#'   reading model ages as observed ages, the fishery or survey ageing error, or
#'   \code{NULL} when the observed ages are the model ages. Its observed ages
#'   must be the age dim of \code{use}.
#'
#' @return A list with \code{nLL} and \code{pred}, both arrays shaped like
#'   \code{use} and zero wherever nothing is fit. The predictions are returned so
#'   they can be reported and plotted directly rather than reconstructed.
#'
#' @keywords internal
get_at_age_source_nLL = function(
  obs_t,
  use,
  ln_sigma,
  source,
  pop,
  arrays,
  obs_se = NULL,
  sd_form = 0,
  like_type = 0,
  const = 0,
  corr_type = 0,
  trans_rho = 0,
  trans_rho_year = 0,
  us_pars = NULL,
  aa_type = 1,
  seas_agg = 0,
  ageing_error = NULL,
  bias_correct_oe = 0
) {

  "[<-" <- RTMB::ADoverload("[<-")

  obs_dim = dim(use)                    # the shape the caller gets back

  # refuse data that are not region by year by season by age by sex by fleet
  if(length(obs_dim) != (if(pop) 7 else 6)) {
    stop("An at-age data source is dimensioned region by year by season by observed age by sex ",
         "by fleet, with a leading population dim for a population-specific source. This one ",
         "arrived with ", length(obs_dim), " dims.")
  }

  source_nLL = array(0, dim = obs_dim)  # zero wherever nothing is fit
  source_pred = array(0, dim = obs_dim)
  if(!any(use == 1)) return(list(nLL = source_nLL, pred = source_pred))

  # obs_t is a flat list of the observations, so record which of them each array position holds
  fit_cells = which(use == 1)
  obs_slot = array(NA, dim = obs_dim)
  obs_slot[fit_cells] = seq_along(fit_cells)

  # data not split by population are one dim short. add a population of one, so the code
  # below indexes both kinds of data the same way
  if(!pop) {
    dim(use) = c(1, dim(use))
    dim(obs_slot) = c(1, dim(obs_slot))
    dim(source_nLL) = c(1, dim(source_nLL))
    dim(source_pred) = c(1, dim(source_pred))
    dim(ln_sigma) = c(1, dim(ln_sigma))
    dim(trans_rho) = c(1, dim(trans_rho))
    dim(trans_rho_year) = c(1, dim(trans_rho_year))
    if(!is.null(obs_se)) dim(obs_se) = c(1, dim(obs_se))
    if(!is.null(us_pars)) dim(us_pars) = c(dim(us_pars)[1], 1, dim(us_pars)[-1])
  }

  n_pop = dim(use)[1]
  n_regions = dim(use)[2]
  n_years = dim(use)[3]
  n_seas = dim(use)[4]
  n_obs_ages = dim(use)[5]
  n_sexes = dim(use)[6]
  n_fleets = dim(use)[7]

  # one setting given for the whole data source is repeated out to one per fleet
  sd_form = rep_len(sd_form, n_fleets)
  like_type = rep_len(like_type, n_fleets)
  corr_type = rep_len(corr_type, n_fleets)
  seas_agg = rep_len(seas_agg, n_fleets)

  # the data must be on the same ages the ageing error reads model ages onto
  if(!is.null(ageing_error) && dim(ageing_error)[3] != n_obs_ages) {
    stop("The at-age observations are on ", n_obs_ages, " ages, but the ageing error reads model ages onto ",
         dim(ageing_error)[3], " observed ages. At-age data are recorded on the observed ages of ",
         "AgeingError, so rebuild the input list through its Setup_Mod_ functions.")
  }

  # a fleet can report its ages one way in some years and another way in others, so this
  # setting is held as year by fleet whatever shape it arrived in
  if(is.null(dim(aa_type)))
    aa_type = base::matrix(
      rep_len(aa_type, n_fleets),
      nrow = n_years,
      ncol = n_fleets,
      byrow = TRUE
    )

  # the full list of populations, regions and sexes, used when a fleet reports them together
  pred_arr = switch(source, catch = arrays$CAA, discard = arrays$DAA, arrays$SrvIAA)
  all_pop = seq_len(dim(pred_arr)[1])
  all_reg = seq_len(dim(pred_arr)[2])
  all_sex = seq_len(dim(pred_arr)[6])

  for(f in seq_len(n_fleets)) {

    if(!any(use[,,,,,,f] == 1)) next

    # 2dar1 fits all of a fleet's years at once, so the fleet cannot have reported its
    # regions or sexes one way in some of those years and another way in others
    if(corr_type[f] == 3 && length(unique(aa_type[,f])) > 1) {
      stop("Fleet ", f, " fits at-age observations as '2dar1', whose correlation runs ",
           "over the whole block of years by ages, but its aggregation changes between ",
           "years. Hold the aggregation constant for this fleet, or give each period its ",
           "own fleet so that each block is its own observation.")
    }

    # this fleet's ageing error, NULL when every age is read as itself
    ae_f = NULL
    if(!is.null(ageing_error)) {
      ae_f = array(ageing_error[,,,f], dim = dim(ageing_error)[1:3])
      if(is_identity_ageing_error(ae_f)) ae_f = NULL
    }

    for(p in seq_len(n_pop)) {
      for(r in seq_len(n_regions)) {
        for(seas in seq_len(n_seas)) {
          for(s in seq_len(n_sexes)) {

            if(!any(use[p,r,,seas,,s,f] == 1)) next

            # a fleet reporting once a year is compared against all seasons added together
            pred_seas = if(seas_agg[f] == 1) seq_len(n_seas) else seas

            # one correlation matrix across ages, shared by every year of this region and sex
            us_corr = if(corr_type[f] == 2) build_us_corr(us_pars[,p,r,s,f], n_obs_ages) else NULL

            if(corr_type[f] == 3) {

              # 2dar1 treats every year and age together as a single observation, so the
              # fleet must have aged its catch in every year of the block
              block = base::matrix(use[p,r,,seas,,s,f], nrow = n_years, ncol = n_obs_ages)
              obs_years = which(base::rowSums(block) > 0)
              obs_ages = which(base::colSums(block) > 0)
              if(!all(block[obs_years,obs_ages] == 1)) {
                stop("A fleet fitting at-age observations as '2dar1' must observe a complete ",
                     "block of ages by years, since a separable correlation is defined over the ",
                     "whole grid. Fleet ", f, " has gaps in that block. Use '1dar1' or 'us', ",
                     "which are defined over whatever ages a cell observes.")
              }

              n_block_years = length(obs_years)
              n_block_ages = length(obs_ages)
              slot = as.vector(obs_slot[p,r,obs_years,seas,obs_ages,s,f]) # year runs fastest

              # add the prediction up over whatever regions and sexes the data were reported over
              split = at_age_split(aa_type[1,f])

              pred = rep(0, n_block_years * n_block_ages)
              k = 1
              for(a in seq_len(n_block_ages)) {   # age by age, so year still runs fastest
                for(y in seq_len(n_block_years)) {
                  pred[k] = get_at_age_obs_prediction(
                    source = source,
                    arrays = arrays,
                    p_idx = if(pop) p else all_pop,
                    r_idx = if(split$region) r else all_reg,
                    s_idx = if(split$sex) s else all_sex,
                    y = obs_years[y],
                    seas = pred_seas,
                    obs_ages = obs_ages[a],
                    f = f,
                    ageing_error = ae_f
                  )
                  k = k + 1
                } # end y loop
              } # end a loop

              sigma = exp(ln_sigma[p,obs_ages,s,f])[rep(seq_len(n_block_ages), each = n_block_years)]
              if(sd_form[f] != 0) {
                sigma = at_age_obs_sd(as.numeric(obs_se[p,r,obs_years,seas,obs_ages,s,f]), sigma, sd_form[f])
              }

              # do lognormal bias correction here
              oe = if(bias_correct_oe == 1) 0.5 * sigma^2 else 0
              pred_t = if(like_type[f] == 0) log(pred + const) - oe else pred

              # one density covers the whole block, stored on its first year and age
              block_nLL = rep(0, length(pred))
              block_nLL[1] = get_at_age_2dar1_nLL(
                matrix(obs_t[slot] - pred_t, nrow = n_block_years),
                matrix(sigma, nrow = n_block_years),
                trans_rho[p,r,s,f],
                trans_rho_year[p,r,s,f]
              )

              source_nLL[p,r,obs_years,seas,obs_ages,s,f] = block_nLL
              source_pred[p,r,obs_years,seas,obs_ages,s,f] = pred

            } else {

              for(y in seq_len(n_years)) {

                obs_ages = which(use[p,r,y,seas,,s,f] == 1)
                if(length(obs_ages) == 0) next

                # subset the vector directly b/c copying element by element loses the OBS tagging
                slot = obs_slot[p,r,y,seas,obs_ages,s,f]
                sigma = exp(ln_sigma[p,obs_ages,s,f])
                if(sd_form[f] != 0) {
                  sigma = at_age_obs_sd(as.numeric(obs_se[p,r,y,seas,obs_ages,s,f]), sigma, sd_form[f])
                }

                # add the prediction up over whatever regions and sexes the data were reported over
                split = at_age_split(aa_type[y,f])
                pred = get_at_age_obs_prediction(
                  source = source,
                  arrays = arrays,
                  p_idx = if(pop) p else all_pop,
                  r_idx = if(split$region) r else all_reg,
                  s_idx = if(split$sex) s else all_sex,
                  y = y,
                  seas = pred_seas,
                  obs_ages = obs_ages,
                  f = f,
                  ageing_error = ae_f
                )

                oe = if(bias_correct_oe == 1) 0.5 * sigma^2 else 0
                pred_t = if(like_type[f] == 0) log(pred + const) - oe else pred

                # iid returns one value per age; a correlated cell puts its density on the first
                cell_nLL = get_at_age_nLL(
                  obs_t[slot],
                  pred_t,
                  sigma,
                  corr_type[f],
                  rho_trans(trans_rho[p,r,s,f]),
                  ages = obs_ages,
                  corr_mat = if(is.null(us_corr)) NULL else us_corr[obs_ages,obs_ages]
                )

                source_nLL[p,r,y,seas,obs_ages,s,f] = cell_nLL
                source_pred[p,r,y,seas,obs_ages,s,f] = pred

              } # end y loop
            } # end correlation form
          } # end s loop
        } # end seas loop
      } # end r loop
    } # end p loop
  } # end f loop

  # coerce shapes back
  dim(source_nLL) = obs_dim
  dim(source_pred) = obs_dim

  return(list(nLL = source_nLL, pred = source_pred))
}
