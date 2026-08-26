# Lognormal at-age observation likelihoods, shared by eight streams: retained
# catch, discards, fishery index and survey index, each aggregated and
# population-specific. Streams differ only in prediction array and sigma.

#' Predicted value for one age-disaggregated observation
#'
#' Catch and discards are already at age and only need their units applied.
#' The indices apply an age-specific catchability to the numbers available to
#' that fleet.
#'
#' @param source Character, one of \code{"catch"}, \code{"discard"},
#'   \code{"fish_index"} or \code{"srv_index"}.
#' @param pop Logical. \code{TRUE} for the population-specific stream, in which
#'   case \code{p} indexes a single population rather than summing over all.
#' @param p,r,y,seas,a,f Population, region, year, season, age and fleet indices.
#' @param arrays Named list of the model arrays the prediction reads:
#'   \code{CAA}, \code{DAA}, \code{SrvIAA}, \code{FishIAA}, \code{WAA_fish}
#'   and \code{catch_units}. The two index arrays already carry their fleet's
#'   selectivity, timing and movement treatment, and the age shape of
#'   catchability lives in that selectivity: a fleet fit age by age uses the
#'   \code{"nonparfree"} selectivity form, whose values carry the height of the
#'   curve as well as its shape.
#'
#' @return The predicted observation, a scalar.
#'
#' @keywords internal
get_at_age_prediction = function(source, pop, p, r, y, seas, a, f, arrays) {

  if(source == "catch" || source == "discard") {
    at_age = if(source == "catch") arrays$CAA else arrays$DAA
    numbers = if(pop) at_age[p,r,y,seas,a,,f] else at_age[,r,y,seas,a,,f] # sum over sex, and pop unless pop specific
    if(arrays$catch_units[f] == 0) return(sum(numbers))                   # abundance
    weight = if(pop) arrays$WAA_fish[p,r,y,seas,a,,f] else arrays$WAA_fish[,r,y,seas,a,,f]
    return(sum(numbers * weight))                                          # biomass
  }

  if(source == "fish_index") { # index numbers at age, selectivity already applied
    numbers = if(pop) arrays$FishIAA[p,r,y,seas,a,,f] else arrays$FishIAA[,r,y,seas,a,,f]
    return(sum(numbers))
  }

  numbers = if(pop) arrays$SrvIAA[p,r,y,seas,a,,f] else arrays$SrvIAA[,r,y,seas,a,,f] # survey available numbers
  return(sum(numbers))
}

#' Evaluate one age-disaggregated observation stream
#'
#' Walks the observations a stream fits and returns their negative log
#' likelihood, shaped like the stream's use array so it can be weighted and
#' reported alongside the aggregated streams.
#'
#' Observations arrive log-scale and already registered through
#' \code{RTMB::OBS}. Registration must happen against the name \code{getAll}
#' supplied, so the caller does it: a vector registered under a local name does
#' not link to the data element, and the objective then diverges from the
#' reported likelihood.
#'
#' @param obs_log Registered log-scale observations for this stream, one element
#'   per cell flagged in \code{use}, in \code{which()} order.
#' @param use Integer array flagging which cells are fit, dimensioned region by
#'   year by season by age by fleet, with a leading population dimension when
#'   \code{pop} is \code{TRUE}.
#' @param ln_sigma Log-scale observation error, over age and fleet, with a
#'   leading population dimension when \code{pop} is \code{TRUE}.
#' @param source,pop,arrays Passed to \code{\link{get_at_age_prediction}}.
#' @param const Small constant added inside the log, matching the aggregated
#'   stream's convention.
#' @param corr_type Integer. \code{0} is \code{"iid"}, \code{1} is
#'   \code{"1dar1"}, an AR(1) across ages within a cell.
#' @param rho Correlation per fleet, used when \code{corr_type} is \code{1}.
#'
#' @return A list with \code{nLL} and \code{pred}, both arrays shaped like
#'   \code{use} and zero wherever nothing is fit. The predictions are returned so
#'   they can be reported and plotted directly rather than reconstructed.
#'
#' @keywords internal
get_at_age_stream_nLL = function(obs_log, use, ln_sigma, source, pop, arrays,
                                 const = 0, corr_type = 0, rho = 0) {

  "[<-" <- RTMB::ADoverload("[<-")

  stream_nLL = array(0, dim = dim(use))  # zero wherever nothing is fit
  stream_pred = array(0, dim = dim(use))
  if(!any(use == 1)) return(list(nLL = stream_nLL, pred = stream_pred))

  fit_cells = which(use == 1)             # linear positions of fitted observations
  obs_slot = array(NA_integer_, dim = dim(use))
  obs_slot[fit_cells] = seq_along(fit_cells) # array position -> slot in obs_log

  cell_margins = if(pop) c(1,2,3,4,6) else c(1,2,3,5) # everything but age
  active_cells = which(apply(use == 1, cell_margins, any))
  cell_index = arrayInd(active_cells, dim(use)[cell_margins])

  for(cell in seq_len(nrow(cell_index))) {

    if(pop) { # population stream carries a leading pop index
      p = cell_index[cell,1]; r = cell_index[cell,2]; y = cell_index[cell,3]
      seas = cell_index[cell,4]; f = cell_index[cell,5]
      observed_ages = which(use[p,r,y,seas,,f] == 1)
    } else {
      p = 0; r = cell_index[cell,1]; y = cell_index[cell,2]
      seas = cell_index[cell,3]; f = cell_index[cell,4]
      observed_ages = which(use[r,y,seas,,f] == 1)
    }
    if(length(observed_ages) == 0) next

    # subset the registered vector directly; copying element by element loses the OBS tagging
    slot = if(pop) obs_slot[p,r,y,seas,observed_ages,f] else obs_slot[r,y,seas,observed_ages,f]
    obs_sd = if(pop) exp(ln_sigma[p,observed_ages,f]) else exp(ln_sigma[observed_ages,f])

    pred = rep(0, length(observed_ages))
    for(age_slot in seq_along(observed_ages)) { # prediction for each observed age
      pred[age_slot] = get_at_age_prediction(source, pop, p, r, y, seas,
                                             observed_ages[age_slot], f, arrays)
    } # end age_slot loop

    # iid returns one value per age; 1dar1 puts the cell's density on the first
    cell_nLL = get_at_age_nLL(obs_log[slot], log(pred + const), obs_sd, corr_type, rho[f])

    for(age_slot in seq_along(observed_ages)) { # scatter back to the arrays
      a = observed_ages[age_slot]
      if(pop) {
        stream_nLL[p,r,y,seas,a,f] = cell_nLL[age_slot]
        stream_pred[p,r,y,seas,a,f] = pred[age_slot]
      } else {
        stream_nLL[r,y,seas,a,f] = cell_nLL[age_slot]
        stream_pred[r,y,seas,a,f] = pred[age_slot]
      }
    } # end age_slot loop
  } # end cell loop

  return(list(nLL = stream_nLL, pred = stream_pred))
}
