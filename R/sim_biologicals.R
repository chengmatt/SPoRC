# Operating model
#
# Growth and movement rebuilt from deviation arrays the operating model holds, by the fit's own functions, so a
# drawn series (a dsem's today) reaches weight at age, the keys, selectivity at age and the movement matrix.

#' Arguments of a model function read from the data and parameter lists
#'
#' Matches each argument of \code{fn} by name: a value given in \code{...}
#' first, then the parameter list, then the data list. Arguments found nowhere
#' keep the function's own default.
#'
#' @param fn A model function such as \code{Get_Growth}.
#' @param data Data list of the fit.
#' @param pars Parameter list at the fitted values.
#' @param ... Values for arguments the lists lack or hold under another name.
#'
#' @return Named list ready for \code{do.call(fn, ...)}.
#'
#' @keywords internal
match_model_args <- function(fn, data, pars, ...) {

  args <- list(...)
  for(nm in setdiff(names(formals(fn)), names(args))) {
    if(!is.null(pars[[nm]])) args[[nm]] <- pars[[nm]]
    else if(!is.null(data[[nm]])) args[[nm]] <- data[[nm]]
  } # end nm loop

  return(args)

} # end function

#' Rebuild growth for every replicate from its deviation arrays
#'
#' Runs the fit's own \code{Get_Growth} at each replicate's
#' \code{ln_growth_devs} and \code{ln_growth_semipar_devs} and writes what it
#' returns over that replicate's slices through \code{take_sim_growth_years},
#' so a drawn growth series reaches the population the way it does in the fit.
#' Under cohort growth only the years before \code{growth_cohort_styr} are
#' built here; the rest come one year at a time from
#' \code{advance_sim_growth_year} inside the annual cycle, and the growth
#' list each replicate advances from is kept in \code{growth_state}.
#'
#' @param sim_env Simulation environment holding \code{dsem_growth_args} from
#'   \code{Setup_Sim_DSEM} and the deviation arrays with the replicate dim last.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
derive_sim_growth <- function(sim_env) {

  args <- sim_env$dsem_growth_args
  cohort <- isTRUE(args$growth_tv_type == 1)
  if(cohort) sim_env$growth_state <- vector("list", sim_env$n_sims)
  # the fit builds the years before the propagation starts up front and advances the rest
  built_yrs <- if(cohort) seq_len(max(1, args$growth_cohort_styr - 1)) else seq_len(sim_env$n_yrs)

  for(sim in seq_len(sim_env$n_sims)) {
    args$ln_growth_devs <- sim_growth_devs(sim_env, "ln_growth_devs", sim)
    args$ln_growth_semipar_devs <- sim_growth_devs(sim_env, "ln_growth_semipar_devs", sim)
    growth <- do.call(Get_Growth, args)
    take_sim_growth_years(sim_env, sim, growth, built_yrs)
    if(cohort) sim_env$growth_state[[sim]] <- growth
  } # end sim loop

  return(invisible(NULL))

} # end function

#' One replicate's growth deviations without the replicate dim
#'
#' @param sim_env Simulation environment.
#' @param par_name \code{"ln_growth_devs"} or \code{"ln_growth_semipar_devs"}.
#' @param sim Replicate.
#'
#' @return The five-dim array \code{Get_Growth} takes, or \code{NULL} when the
#'   operating model holds no such array.
#'
#' @keywords internal
sim_growth_devs <- function(sim_env, par_name, sim) {
  devs <- sim_env[[par_name]]
  if(is.null(devs)) return(NULL)
  array(devs[,,,,,sim], dim = dim(devs)[-6])
}

#' Advance one replicate's cohort growth by a year
#'
#' The fit's \code{Get_Growth_Year} at this replicate's deviations and its own
#' start of year numbers at age, which blend the plus group, called from
#' \code{run_annual_cycle} before anything else in the year is formed, as the
#' fit's population loop does.
#'
#' @param y Year.
#' @param sim Replicate.
#' @param sim_env Simulation environment holding \code{growth_state} from
#'   \code{derive_sim_growth}.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
advance_sim_growth_year <- function(y, sim, sim_env) {

  args <- sim_env$dsem_growth_args
  args <- args[names(args) %in% names(formals(Get_Growth_Year))]
  args$ln_growth_devs <- sim_growth_devs(sim_env, "ln_growth_devs", sim)
  args$ln_growth_semipar_devs <- sim_growth_devs(sim_env, "ln_growth_semipar_devs", sim)
  NAA_y <- array(sim_env$NAA[,,y,1,,,sim], dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_ages, sim_env$n_sexes))

  growth <- do.call(Get_Growth_Year, c(args, list(growth = sim_env$growth_state[[sim]], y = y, NAA_y = NAA_y)))
  sim_env$growth_state[[sim]] <- growth
  take_sim_growth_years(sim_env, sim, growth, y)

  return(invisible(NULL))

} # end function

#' Write years of one replicate's growth into the operating model
#'
#' Weight at age, the size-age keys, and then what the fit forms from them in
#' \code{compute_mortality_year}: a fleet asking for it gets the weight of
#' what it selects, and with selectivity at length each fleet's selectivity at
#' age is its curve read through the key. Weight at age is left alone when the
#' fit takes it as data (\code{derive_waa = 0}).
#'
#' @param sim_env Simulation environment holding \code{dsem_length_sel} from
#'   \code{Setup_Sim_DSEM}.
#' @param sim Replicate.
#' @param growth Growth list from \code{Get_Growth} or \code{Get_Growth_Year}.
#' @param yrs Years to write.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
take_sim_growth_years <- function(sim_env, sim, growth, yrs) {

  args <- sim_env$dsem_growth_args
  sel <- sim_env$dsem_length_sel
  len_mid <- growth_len_mid(args$growth_len_lower)
  n_pop <- sim_env$n_pop
  n_regions <- sim_env$n_regions
  n_seas <- sim_env$n_seas
  n_sexes <- sim_env$n_sexes

  for(y in yrs) {

    if(!is.null(sim_env$SizeAgeTrans_fish)) sim_env$SizeAgeTrans_fish[,,y,,,,,,sim] <- growth$SizeAgeTrans_fish[,,y,,,,,]
    if(!is.null(sim_env$SizeAgeTrans_srv)) sim_env$SizeAgeTrans_srv[,,y,,,,,,sim] <- growth$SizeAgeTrans_srv[,,y,,,,,]

    if(!is.null(growth$WAA)) {
      WAA_fish <- growth$WAA_fish
      WAA_srv <- growth$WAA_srv
      if(sel$fish_selex_type == 1 && any(sel$fish_waa_selected == 1)) WAA_fish <- growth_selected_waa_year(WAA_fish, growth$SizeAgeTrans_fish, sel$fish_sel_l, args$wt_len_pars, len_mid, sel$fish_waa_selected, y, n_pop, n_regions, n_seas, n_sexes)
      if(sel$srv_selex_type == 1 && any(sel$srv_waa_selected == 1)) WAA_srv <- growth_selected_waa_year(WAA_srv, growth$SizeAgeTrans_srv, sel$srv_sel_l, args$wt_len_pars, len_mid, sel$srv_waa_selected, y, n_pop, n_regions, n_seas, n_sexes)
      sim_env$WAA[,,y,,,,sim] <- growth$WAA[,,y,,,]
      sim_env$WAA_fish[,,y,,,,,sim] <- WAA_fish[,,y,,,,]
      sim_env$WAA_srv[,,y,,,,,sim] <- WAA_srv[,,y,,,,]
    }

    # selectivity at length read through this replicate's key, one curve per fleet, sex, season and origin
    for(p in 1:n_pop) for(r in 1:n_regions) for(seas in 1:n_seas) for(s in 1:n_sexes) {
      if(sel$fish_selex_type == 1 || sel$ret_selex_type == 1) for(f in seq_len(dim(growth$SizeAgeTrans_fish)[8])) {
        key <- growth$SizeAgeTrans_fish[p,r,y,seas,,,s,f]
        if(sel$fish_selex_type == 1) sim_env$fish_sel[p,r,y,seas,,s,f,sim] <- sel$fish_sel_l[r,y,,s,f] %*% key
        if(sel$ret_selex_type == 1) sim_env$ret_sel[p,r,y,seas,,s,f,sim] <- sel$ret_sel_l[r,y,,s,f] %*% key
      } # end f loop
      if(sel$srv_selex_type == 1) for(sf in seq_len(dim(growth$SizeAgeTrans_srv)[8])) {
        sim_env$srv_sel[p,r,y,seas,,s,sf,sim] <- sel$srv_sel_l[r,y,,s,sf] %*% growth$SizeAgeTrans_srv[p,r,y,seas,,,s,sf]
      } # end sf loop
    } # end p, r, seas, s loop

  } # end y loop

  return(invisible(NULL))

} # end function

#' Rebuild movement for every replicate from its deviation array
#'
#' Runs the fit's own \code{Get_Movement} at each replicate's \code{move_devs}
#' and writes the movement matrix (and the rate matrix under continuous
#' movement) over that replicate's slice.
#'
#' @param sim_env Simulation environment holding \code{dsem_move_args} from
#'   \code{Setup_Sim_DSEM} and \code{move_devs} with the replicate dim last.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
derive_sim_movement <- function(sim_env) {

  args <- sim_env$dsem_move_args

  for(sim in seq_len(sim_env$n_sims)) {
    args$move_devs <- array(sim_env$move_devs[,,,,,,,sim], dim = dim(sim_env$move_devs)[-8])
    moved <- do.call(Get_Movement, args)
    sim_env$Movement[,,,,,,,sim] <- moved$Movement
    if(!is.null(sim_env$Mrate) && !is.null(moved$Mrate)) sim_env$Mrate[,,,,,,,sim] <- moved$Mrate
  } # end sim loop

  return(invisible(NULL))

} # end function
