# Operating model
#
# Growth and movement deviations drawn from the processes the fit penalizes, and both rebuilt from the
# deviation arrays by the fit's own functions, so a draw reaches weight at age, selectivity and movement.

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

  matched <- list(...)
  for(name in setdiff(names(formals(fn)), names(matched))) { # get function arg names
    if(!is.null(pars[[name]])) matched[[name]] <- pars[[name]] # if matching pars then input
    else if(!is.null(data[[name]])) matched[[name]] <- data[[name]] # if matching data then input
  } # end name loop

  return(matched)

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
#' @param sim_env Simulation environment holding \code{growth_args} from
#'   \code{Setup_Sim_Growth_RE} or \code{dsem_growth_args} from
#'   \code{Setup_Sim_DSEM}, and the deviation arrays with the replicate dim last.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
derive_sim_growth <- function(sim_env) {

  growth_args <- sim_growth_args(sim_env)
  cohort_growth <- isTRUE(growth_args$growth_tv_type == 1)
  if(cohort_growth) sim_env$growth_state <- vector("list", sim_env$n_sims)

  # the fit builds the years before the propagation starts up front and advances the rest
  built_yrs <- if(cohort_growth) seq_len(max(1, growth_args$growth_cohort_styr - 1)) else seq_len(sim_env$n_yrs)

  for(sim in seq_len(sim_env$n_sims)) {
    growth_args$ln_growth_devs <- sim_growth_devs(sim_env, "ln_growth_devs", sim)
    growth_args$ln_growth_semipar_devs <- sim_growth_devs(sim_env, "ln_growth_semipar_devs", sim)
    growth <- do.call(Get_Growth, growth_args)
    take_sim_growth_years(sim_env, sim, growth, built_yrs)
    if(cohort_growth) sim_env$growth_state[[sim]] <- growth
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

  growth_args <- sim_growth_args(sim_env)
  growth_args <- growth_args[names(growth_args) %in% names(formals(Get_Growth_Year))]
  growth_args$ln_growth_devs <- sim_growth_devs(sim_env, "ln_growth_devs", sim)
  growth_args$ln_growth_semipar_devs <- sim_growth_devs(sim_env, "ln_growth_semipar_devs", sim)
  NAA_y <- array(sim_env$NAA[,,y,1,,,sim], dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_ages, sim_env$n_sexes))

  growth <- do.call(Get_Growth_Year, c(growth_args, list(growth = sim_env$growth_state[[sim]], y = y, NAA_y = NAA_y)))
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
#' @param sim_env Simulation environment holding \code{growth_length_sel} from
#'   \code{Setup_Sim_Growth_RE} or \code{Setup_Sim_DSEM}.
#' @param sim Replicate.
#' @param growth Growth list from \code{Get_Growth} or \code{Get_Growth_Year}.
#' @param yrs Years to write.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
take_sim_growth_years <- function(sim_env, sim, growth, yrs) {

  growth_args <- sim_growth_args(sim_env)
  length_sel <- sim_env$growth_length_sel
  len_mid <- growth_len_mid(growth_args$growth_len_lower)
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
      if(length_sel$fish_selex_type == 1 && any(length_sel$fish_waa_selected == 1)) WAA_fish <- growth_selected_waa_year(WAA_fish, growth$SizeAgeTrans_fish, length_sel$fish_sel_l, growth_args$wt_len_pars, len_mid, length_sel$fish_waa_selected, y, n_pop, n_regions, n_seas, n_sexes)
      if(length_sel$srv_selex_type == 1 && any(length_sel$srv_waa_selected == 1)) WAA_srv <- growth_selected_waa_year(WAA_srv, growth$SizeAgeTrans_srv, length_sel$srv_sel_l, growth_args$wt_len_pars, len_mid, length_sel$srv_waa_selected, y, n_pop, n_regions, n_seas, n_sexes)
      sim_env$WAA[,,y,,,,sim] <- growth$WAA[,,y,,,]
      sim_env$WAA_fish[,,y,,,,,sim] <- WAA_fish[,,y,,,,]
      sim_env$WAA_srv[,,y,,,,,sim] <- WAA_srv[,,y,,,,]
    }

    # selectivity at length read through this replicate's key, one curve per fleet, sex, season and origin
    for(p in 1:n_pop) for(r in 1:n_regions) for(seas in 1:n_seas) for(s in 1:n_sexes) {
      if(length_sel$fish_selex_type == 1 || length_sel$ret_selex_type == 1) for(f in seq_len(dim(growth$SizeAgeTrans_fish)[8])) {
        key <- growth$SizeAgeTrans_fish[p,r,y,seas,,,s,f]
        if(length_sel$fish_selex_type == 1) sim_env$fish_sel[p,r,y,seas,,s,f,sim] <- length_sel$fish_sel_l[r,y,,s,f] %*% key
        if(length_sel$ret_selex_type == 1) sim_env$ret_sel[p,r,y,seas,,s,f,sim] <- length_sel$ret_sel_l[r,y,,s,f] %*% key
      } # end f loop
      if(length_sel$srv_selex_type == 1) for(sf in seq_len(dim(growth$SizeAgeTrans_srv)[8])) {
        sim_env$srv_sel[p,r,y,seas,,s,sf,sim] <- length_sel$srv_sel_l[r,y,,s,sf] %*% growth$SizeAgeTrans_srv[p,r,y,seas,,,s,sf]
      } # end sf loop
    } # end p, r, seas, s loop

  } # end y loop

  return(invisible(NULL))

} # end function

#' The growth arguments the operating model rebuilds growth with
#'
#' \code{growth_args} from \code{Setup_Sim_Growth_RE} when the fit penalizes
#' its own growth deviations, else \code{dsem_growth_args} from
#' \code{Setup_Sim_DSEM}.
#'
#' @param sim_env Simulation environment.
#'
#' @return Named list of \code{Get_Growth} arguments, or \code{NULL}.
#'
#' @keywords internal
sim_growth_args <- function(sim_env) {
  if(!is.null(sim_env$growth_args)) sim_env$growth_args else sim_env$dsem_growth_args
}

#' Selectivity at length the operating model reads the rebuilt keys through
#'
#' What the fit forms after growth needs the selectivity at length the report
#' holds, run out to the operating model's years the way the closed loop extends
#' its other inputs, and which fleets weigh their catch by what they select.
#'
#' @param data Data list of the fit.
#' @param rep Report of the fit, needed only under selectivity at length.
#' @param n_sim_yrs Years the operating model runs.
#'
#' @return List with the three selectivity type flags, the selected weight
#'   flags and the fishery, retention and survey selectivity at length.
#'
#' @keywords internal
sim_growth_length_sel <- function(data, rep, n_sim_yrs) {

  # set sel at length true
  sel_at_length <- c(fish = isTRUE(data$fish_selex_type == 1),
                     ret = isTRUE(data$ret_selex_type == 1),
                     srv = isTRUE(data$srv_selex_type == 1))

  if(any(sel_at_length) && is.null(rep))
    stop("The fit has selectivity at length, so rebuilding growth in the operating model needs rep, ",
         "the fit's report, for the selectivity at length the rebuilt keys are read through.")

  sel_length_arrays <- list()
  for(sel_name in names(sel_at_length)[sel_at_length]) {
    sel_array <- rep[[paste0(sel_name, "_sel_l")]] # [region, year, len, sex, fleet]
    n_rep_yrs <- dim(sel_array)[2]
    if(n_rep_yrs > n_sim_yrs) sel_array <- sel_array[,seq_len(n_sim_yrs),,,,drop = FALSE]
    if(n_rep_yrs < n_sim_yrs) sel_array <- extend_years(sel_array, n_sim_yrs - n_rep_yrs, 2, "last")
    sel_length_arrays[[sel_name]] <- sel_array
  } # end sel_name loop

  # return arguments
  list(fish_selex_type = as.integer(sel_at_length[["fish"]]),
       ret_selex_type = as.integer(sel_at_length[["ret"]]),
       srv_selex_type = as.integer(sel_at_length[["srv"]]),
       fish_waa_selected = if(is.null(data$fish_waa_selected)) 0 else data$fish_waa_selected,
       srv_waa_selected = if(is.null(data$srv_waa_selected)) 0 else data$srv_waa_selected,
       fish_sel_l = sel_length_arrays$fish,
       ret_sel_l = sel_length_arrays$ret,
       srv_sel_l = sel_length_arrays$srv)

} # end function

#' Growth deviations in the operating model
#'
#' Stores what each replicate needs to rebuild growth from its own deviations:
#' the fit's time variation codes, the semi-parametric form and its ages, the
#' deviation maps, the process error parameters, the \code{Get_Growth}
#' arguments and, under selectivity at length, the report's selectivity at
#' length the rebuilt keys are read through. With the fit's own penalty on the
#' deviations, \code{\link{draw_sim_growth_devs}} draws them in
#' \code{Setup_sim_env}. With a dsem holding or linking any growth deviation,
#' \code{\link{Setup_Sim_DSEM}} builds both arrays itself and nothing is stored
#' here, so the deviations of a parameter outside the arrows keep the fit's
#' values. Nothing is stored either when growth is data or has no deviations.
#'
#' @param sim_list Simulation list with \code{n_yrs} and \code{n_sims}.
#' @param data Data list of the fit.
#' @param pars Parameter list at the fitted values.
#' @param rep Report of the fit, needed under selectivity at length.
#'
#' @return \code{sim_list} with \code{growth_args} and the fields above added.
#'
#' @export Setup_Sim_Growth_RE
Setup_Sim_Growth_RE <- function(sim_list, data, pars, rep = NULL) {

  # return if no growth model
  if(is.null(data$growth_model) || data$growth_model == 0) return(sim_list)

  # return if no tv devs
  tv_devs_on <- !is.null(data$growth_tv_model) && any(data$growth_tv_model > 0)
  semipar_devs_on <- isTRUE(data$growth_semipar > 0)
  if(!tv_devs_on && !semipar_devs_on) return(sim_list)

  # return if using dsem
  dsem_holds <- isTRUE(any(data$growth_tv_dsem == 1)) || isTRUE(data$growth_semipar_dsem == 1) ||
    isTRUE(any(data$dsem_link_par %in% c("ln_growth_devs", "ln_growth_semipar_devs")))
  if(dsem_holds) return(sim_list)

  # return if not using dsem to get devs
  for(name in c("growth_tv_model", "growth_semipar", "growth_semipar_bins", "map_ln_growth_devs", "map_ln_growth_semipar_devs")) sim_list[[name]] <- data[[name]]
  sim_list$growth_pe_pars <- pars$growth_pe_pars
  sim_list$growth_args <- match_model_args(Get_Growth, data, pars, n_yrs = sim_list$n_yrs)
  sim_list$growth_length_sel <- sim_growth_length_sel(data, rep, sim_list$n_yrs)
  return(sim_list)

} # end function

#' A deviation map run map_out over the operating model's years
#'
#' Years past the fit are active in every column the fit varies, under fresh
#' levels that share what that column's last active fitted year shares, as
#' projection years are always active for movement. Years the operating model
#' does not run are cut.
#'
#' @param map Integer array \code{[pop, region, year, bin, sex]}, \code{NA}
#'   where a cell is fixed.
#' @param n_yrs Years the operating model runs.
#'
#' @return The map with \code{n_yrs} years.
#'
#' @keywords internal
sim_map_over_years <- function(map, n_yrs) {

  map_dim <- dim(map)
  if(map_dim[3] >= n_yrs) return(map[,,seq_len(n_yrs),,,drop = FALSE])
  map_out <- array(NA_real_, dim = c(map_dim[1:2], n_yrs, map_dim[4:5]))
  map_out[,,seq_len(map_dim[3]),,] <- map
  top_level <- if(all(is.na(map))) 0 else max(map, na.rm = TRUE)
  for(bin in seq_len(map_dim[4])) {
    active_yrs <- which(apply(map[,,,bin,,drop = FALSE], 3, function(m) any(!is.na(m))))
    if(length(active_yrs) == 0) next
    for(extra_yr in seq_len(n_yrs - map_dim[3])) map_out[,,map_dim[3] + extra_yr,bin,] <- map[,,max(active_yrs),bin,] + top_level * extra_yr
  } # end bin loop
  return(map_out)

} # end function

#' One draw of a deviation surface from the process the estimation model penalizes
#'
#' The reverse of \code{\link{Get_PE_loglik}}. A shared level is drawn once and
#' written wherever it appears, at the sd of the slot the penalty reads it at.
#' Under iid every level is its own normal. Under the random walk each level
#' steps from the year before; the first year starts at the walk's own sd, as
#' the recruitment walk does in the operating model, since the diffuse start
#' the estimation model gives it only leaves the level free. The 3D GMRF and
#' the separable AR1 draw each population, region and sex's whole surface over
#' years and \code{bins} from the form's precision, conditional on the cells
#' the map fixes at zero, which is the density the penalty evaluates at them.
#'
#' @param PE_model Integer process error code, as \code{Get_PE_loglik} reads it.
#' @param map Integer array \code{[pop, region, year, bin, sex]} of the levels,
#'   \code{NA} where a cell is fixed.
#' @param pe_pars Array \code{[pop, region, slot, sex]}, the half of
#'   \code{growth_pe_pars} the form reads.
#' @param bins Bins the correlated forms run over.
#'
#' @return Array shaped as \code{map}, zero where the map is \code{NA}.
#'
#' @keywords internal
draw_growth_pe_surface <- function(PE_model, map, pe_pars, bins) {

  map_dim <- dim(map)
  devs <- array(0, dim = map_dim)
  levels <- sort(unique(as.vector(map))) # sort drops the fixed cells
  if(length(levels) == 0) return(devs)

  if(PE_model %in% c(1, 2)) {
    for(level in levels) {
      level_cells <- which(map == level)
      first_cell <- arrayInd(level_cells[1], map_dim) # the slot the penalty reads this level at
      p <- first_cell[1]; r <- first_cell[2]; y <- first_cell[3]; bin <- first_cell[4]; s <- first_cell[5]
      sigma <- exp(pe_pars[p,r,bin,s])
      draw <- if(PE_model == 1 || y == 1) stats::rnorm(1, 0, sigma) else devs[p,r,y - 1,bin,s] + stats::rnorm(1, 0, sigma)
      devs[level_cells] <- draw
    } # end level loop
    return(devs)
  } # end iid or random walk

  # one surface per population, region and sex, following the specified form's precision matrix
  for(p in seq_len(map_dim[1])) for(r in seq_len(map_dim[2])) for(s in seq_len(map_dim[5])) {
    map_slice <- array(map[p,r,,bins,s], dim = c(map_dim[3], length(bins))) # [year, bin]
    if(all(is.na(map_slice))) next
    if(PE_model == 5) {
      rho_bin <- rho_trans(pe_pars[p,r,1,s])
      rho_year <- rho_trans(pe_pars[p,r,2,s])
      marginal_sd <- exp(pe_pars[p,r,4,s]) / sqrt(1 - rho_year^2) / sqrt(1 - rho_bin^2) # marginal sd from the conditional one
      precision <- kronecker(ar1_precision(length(bins), rho_bin), ar1_precision(map_dim[3], rho_year)) / marginal_sd^2 # year fastest
      is_free <- !is.na(as.vector(map_slice))
    } else {
      precision <- Get_3d_precision(length(bins), map_dim[3], pe_pars[p,r,1,s], pe_pars[p,r,2,s], pe_pars[p,r,3,s], pe_pars[p,r,4,s],
                                    Var_Type = if(PE_model == 3) 0 else 1) # bin fastest
      precision <- as.matrix(precision)
      precision <- (precision + t(precision)) / 2
      is_free <- !is.na(as.vector(t(map_slice)))
    }
    # the free cells given the fixed ones at zero have the precision's free block
    nodes <- numeric(length(is_free))
    nodes[is_free] <- backsolve(chol(precision[is_free, is_free, drop = FALSE]), stats::rnorm(sum(is_free)))
    devs[p,r,,bins,s] <- if(PE_model == 5) matrix(nodes, map_dim[3], length(bins)) else t(matrix(nodes, length(bins), map_dim[3]))
  } # end p, r, s loop

  return(devs)

} # end function

#' Draw growth deviations for every replicate and rebuild growth
#'
#' Each varying growth parameter's series and the semi-parametric surface are
#' drawn by \code{\link{draw_growth_pe_surface}} from the form the estimation
#' model penalizes them under, at the fitted process error parameters: a
#' parameter's sd in the time-varying half of \code{growth_pe_pars} and the
#' surface's in the semi-parametric half, as \code{Get_PE_loglik} reads them.
#' The fit's maps decide what is drawn: a cell the map fixes stays at zero, a
#' shared level takes one draw, and years past the fit are active. Growth is
#' then rebuilt through \code{derive_sim_growth}; under cohort growth that
#' builds only the years before the propagation starts, and the annual cycle
#' advances the rest from each replicate's own numbers at age.
#'
#' @param sim_env Simulation environment holding what \code{Setup_Sim_Growth_RE}
#'   stored, \code{n_yrs} and \code{n_sims}.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
draw_sim_growth_devs <- function(sim_env) {

  # get dimensions
  growth_args <- sim_env$growth_args
  n_yrs <- sim_env$n_yrs
  n_sims <- sim_env$n_sims
  pe_pars <- sim_env$growth_pe_pars # [pop, region, slot, sex, half]
  growth_tv_model <- if(is.null(sim_env$growth_tv_model)) rep(0, dim(growth_args$ln_growth_devs)[4]) else sim_env$growth_tv_model
  growth_semipar <- if(is.null(sim_env$growth_semipar)) 0 else sim_env$growth_semipar

  # get mapping stuff
  map_tv_devs <- sim_map_over_years(sim_env$map_ln_growth_devs, n_yrs)
  map_semipar_devs <- sim_map_over_years(sim_env$map_ln_growth_semipar_devs, n_yrs)
  bins <- if(is.null(sim_env$growth_semipar_bins)) seq_len(dim(map_semipar_devs)[4]) else sim_env$growth_semipar_bins
  tv_dim <- dim(map_tv_devs)
  semipar_dim <- dim(map_semipar_devs)

  # do draws of growth devs
  tv_devs <- array(0, dim = c(tv_dim, n_sims))
  semipar_devs <- array(0, dim = c(semipar_dim, n_sims))
  for(sim in seq_len(n_sims)) {
    # one column wide surface per varying parameter, its sd in that parameter's time-varying slot
    for(par_idx in which(growth_tv_model > 0)) {
      pe_par <- array(pe_pars[,,par_idx,,1], dim = c(tv_dim[1:2], 1, tv_dim[5]))
      tv_devs[,,,par_idx,,sim] <- draw_growth_pe_surface(growth_tv_model[par_idx], map_tv_devs[,,,par_idx,,drop = FALSE], pe_par, 1)
    } # end par_idx loop
    if(growth_semipar > 0) semipar_devs[,,,,,sim] <- draw_growth_pe_surface(growth_semipar, map_semipar_devs, array(pe_pars[,,,,2], dim = dim(pe_pars)[1:4]), bins)
  } # end sim loop

  # return stuff
  sim_env$ln_growth_devs <- tv_devs
  sim_env$ln_growth_semipar_devs <- semipar_devs
  derive_sim_growth(sim_env)
  return(invisible(NULL))

} # end function

#' Rebuild movement for every replicate from its deviation array
#'
#' Runs the fit's own \code{Get_Movement} at each replicate's \code{move_devs}
#' and writes the movement matrix (and the rate matrix under continuous
#' movement) over that replicate's slice.
#'
#' @param sim_env Simulation environment holding \code{move_args} from
#'   \code{Setup_Sim_Movement} or \code{dsem_move_args} from
#'   \code{Setup_Sim_DSEM}, and \code{move_devs} with the replicate dim last.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
derive_sim_movement <- function(sim_env) {

  move_args <- if(!is.null(sim_env$move_args)) sim_env$move_args else sim_env$dsem_move_args

  for(sim in seq_len(sim_env$n_sims)) {
    move_args$move_devs <- array(sim_env$move_devs[,,,,,,,sim], dim = dim(sim_env$move_devs)[-8])
    movement <- do.call(Get_Movement, move_args)
    sim_env$Movement[,,,,,,,sim] <- movement$Movement
    if(!is.null(sim_env$Mrate) && !is.null(movement$Mrate)) sim_env$Mrate[,,,,,,,sim] <- movement$Mrate
  } # end sim loop

  return(invisible(NULL))

} # end function

#' Movement deviations in the operating model
#'
#' Stores what each replicate needs to rebuild movement from its own deviations:
#' the fit's switches, deviation map, surfaces, active years, seasons and ages,
#' process error and correlation parameters, and the \code{Get_Movement}
#' arguments, and the fit's own deviations. With the fit's own penalty on the
#' deviations, \code{\link{draw_sim_move_devs}} keeps the fit's values over
#' the simulation list's \code{n_cond_yrs} and draws the years after them in
#' \code{Setup_sim_env}; with a dsem holding or linking them,
#' \code{\link{Setup_Sim_DSEM}} draws them and the arguments stored here
#' rebuild movement. Nothing is stored when the fit has no movement deviations
#' or movement is fixed.
#'
#' @param sim_list Simulation list with \code{n_yrs}, \code{n_sims} and
#'   \code{expm_nsub}.
#' @param data Data list of the fit.
#' @param pars Parameter list at the fitted values.
#'
#' @return \code{sim_list} with \code{move_args} and the fields above added.
#'
#' @export Setup_Sim_Movement
Setup_Sim_Movement <- function(sim_list, data, pars) {

  if(is.null(data$map_move_devs) || all(is.na(data$map_move_devs)) || isTRUE(data$use_fixed_movement == 1)) return(sim_list) # no deviations to draw

  n_fit_yrs <- length(data$years)
  n_sim_yrs <- sim_list$n_yrs
  for(name in c("move_pop_re", "move_year_re", "move_seas_re", "move_age_re", "move_sex_re", "move_dsem", "move_pairs", "map_move_devs",
              "move_pe_block", "move_pop_block", "move_year_block", "move_seas_block", "move_age_block", "move_sex_block")) sim_list[[name]] <- data[[name]]
  for(name in c("move_pe_pars", "move_pop_corr_pars", "move_seas_corr_pars", "move_sex_corr_pars")) sim_list[[name]] <- pars[[name]]
  sim_list$move_devs_fit <- pars$move_devs # the conditioned years reproduce these

  # movement keeps the fit's year count so covariates stop at the data, as the dsem route does
  sim_list$move_args <- match_model_args(Get_Movement, data, pars,
                                         n_yrs = min(n_fit_yrs, n_sim_yrs),
                                         n_proj_yrs_devs = max(0, n_sim_yrs - n_fit_yrs),
                                         n_ages = length(data$ages),
                                         expm_nsub = if(is.null(data$move_expm_nsub)) 0 else data$move_expm_nsub)
  return(sim_list)

} # end function

#' Draw movement deviations for every replicate and rebuild movement
#'
#' The first \code{n_cond_yrs} years of every replicate take the fit's own
#' deviations, and the years after them are drawn from the process the
#' estimation model penalizes: a stationary AR1 over ages where that dim is
#' \code{"ar1"}, an unstructured correlation across populations, seasons or
#' sexes where a dim is \code{"us"}, independent otherwise, with the
#' conditional sd \code{exp(move_pe_pars[..., 1])} raised to the marginal by
#' \code{1 / sqrt(1 - rho^2)} for each \code{"ar1"} dim, as the estimation
#' model reads it. An AR1 over years continues from the last conditioned year
#' rather than restarting, a shared year deviation (\code{"none"}) continues
#' unchanged, and with no conditioned years the whole series is drawn from its
#' stationary distribution. One value is drawn per block of each dim and
#' written into every level of that block, as the estimation model's map ties
#' them; inactive levels stay at zero, projection years are always active, and a
#' cell the map holds as \code{NA} stays at zero. Movement is then rebuilt
#' through \code{derive_sim_movement}.
#'
#' @param sim_env Simulation environment holding what \code{Setup_Sim_Movement}
#'   stored, \code{n_yrs}, \code{n_sims} and \code{n_cond_yrs} (\code{NULL} or
#'   \code{0} draws every year).
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
draw_sim_move_devs <- function(sim_env) {

  fit_dims <- dim(sim_env$move_args$move_devs) # [pop, from, to, year, season, age, sex] of the fit
  n_years_fit <- fit_dims[4]
  n_years_sim <- sim_env$n_yrs
  n_sims <- sim_env$n_sims
  move_pairs <- sim_env$move_pairs
  n_cond <- min(n_years_sim, if(is.null(sim_env$n_cond_yrs)) 0 else sim_env$n_cond_yrs) # years that take the fit's deviations
  drawn_years <- if(n_cond < n_years_sim) (n_cond + 1):n_years_sim else integer(0) # years drawn after them
  year_ar1 <- sim_env$move_year_re == 2
  year_shared <- sim_env$move_year_re == 0

  # the conditioned years are the fit's, for every replicate
  move_devs <- array(0, dim = c(fit_dims[1:3], n_years_sim, fit_dims[5:7], n_sims))
  for(sim in seq_len(n_sims)) move_devs[,,,seq_len(n_cond),,,,sim] <- sim_env$move_devs_fit[,,,seq_len(n_cond),,,]
  if(length(drawn_years) == 0) {
    sim_env$move_devs <- move_devs
    derive_sim_movement(sim_env)
    return(invisible(NULL))
  }

  # the block of every level, with the fit's year blocks extended so projection years stay active
  year_block <- sim_env$move_year_block[seq_len(min(n_years_fit, n_years_sim))]
  if(n_years_sim > n_years_fit) {
    extra <- (n_years_fit + 1):n_years_sim
    year_block <- c(year_block, if(year_shared) rep(1, length(extra)) else max(year_block, na.rm = TRUE) + seq_along(extra))
  }
  blocks <- list(pop = sim_env$move_pop_block, season = sim_env$move_seas_block, age = sim_env$move_age_block, sex = sim_env$move_sex_block)
  active <- lapply(blocks, function(block_of) which(!is.na(block_of))) # the levels that hold a deviation
  take <- lapply(blocks, function(block_of) block_of[!is.na(block_of)]) # the drawn block each active level reads
  n_blocks <- sapply(blocks, function(block_of) max(c(0, block_of), na.rm = TRUE))
  first_of_block <- function(block_of) match(seq_len(max(c(0, block_of), na.rm = TRUE)), block_of)
  reads <- lapply(blocks, first_of_block) # the level the map is read at for each block
  drawn_year_blocks <- unique(year_block[drawn_years][!is.na(year_block[drawn_years])]) # the year blocks drawn, in order
  active_years <- drawn_years[!is.na(year_block[drawn_years])]
  take_year <- match(year_block[active_years], drawn_year_blocks)
  n_drawn_years <- length(drawn_year_blocks)
  read_year <- pmin(match(drawn_year_blocks, year_block), n_years_fit) # the map has the fit's years only
  last_cond_year <- if(n_cond > 0) max(which(!is.na(year_block[seq_len(n_cond)])), 0) else 0 # the year an ar1 continues from

  # the factors that color each dim: ar1 over ages, unstructured over populations, seasons and sexes, identity otherwise
  ar1_chol <- function(n, rho) if(n == 1 || rho == 0) diag(n) else t(chol(rho^abs(outer(1:n, 1:n, "-"))))
  us_or_identity <- function(code, pars, n) if(code == 2) build_us_chol(pars, n) else diag(n)
  chol_pop <- us_or_identity(sim_env$move_pop_re, sim_env$move_pop_corr_pars, n_blocks[["pop"]])
  chol_season <- us_or_identity(sim_env$move_seas_re, sim_env$move_seas_corr_pars, n_blocks[["season"]])
  chol_sex <- us_or_identity(sim_env$move_sex_re, sim_env$move_sex_corr_pars, n_blocks[["sex"]])

  for(b in seq_len(max(sim_env$move_pe_block))) {

    # the pairs in this block, each drawn on its own, and which drawn cells the map estimates
    pairs_drawn <- which(sim_env$move_pe_block == b)
    n_pairs_drawn <- length(pairs_drawn)
    n_drawn <- c(n_pairs_drawn, n_blocks[["pop"]], n_drawn_years, n_blocks[["season"]], n_blocks[["age"]], n_blocks[["sex"]]) # [pair, pop, year, season, age, sex]
    estimated <- array(FALSE, dim = n_drawn)
    last_fit <- array(0, dim = n_drawn[-3]) # the last conditioned year's deviation, per drawn cell
    for(j in 1:n_pairs_drawn) {
      region_from <- move_pairs[pairs_drawn[j], 1]
      region_to <- move_pairs[pairs_drawn[j], 2]
      estimated[j,,,,,] <- !is.na(sim_env$map_move_devs[reads$pop, region_from, region_to, read_year, reads$season, reads$age, reads$sex, drop = FALSE])
      if(last_cond_year > 0) last_fit[j,,,,] <- sim_env$move_devs_fit[reads$pop, region_from, region_to, last_cond_year, reads$season, reads$age, reads$sex, drop = FALSE]
    } # end j loop
    if(!any(estimated)) next # nothing estimated in this block

    # this block's correlations and marginal sd
    region_from <- move_pairs[pairs_drawn[1], 1]
    region_to <- move_pairs[pairs_drawn[1], 2]
    rho_year <- if(year_ar1) rho_trans(sim_env$move_pe_pars[region_from, region_to, 3]) else 0
    rho_age <- if(sim_env$move_age_re == 2) rho_trans(sim_env$move_pe_pars[region_from, region_to, 2]) else 0
    sd_marginal <- exp(sim_env$move_pe_pars[region_from, region_to, 1]) / sqrt(1 - rho_year^2) / sqrt(1 - rho_age^2)
    chol_age <- ar1_chol(n_blocks[["age"]], rho_age)
    continues <- year_ar1 && last_cond_year > 0 # the ar1 runs on from the fit's last year rather than restarting
    chol_year <- if(year_ar1 && !continues) ar1_chol(n_drawn_years, rho_year) else diag(n_drawn_years)

    for(sim in 1:n_sims) {
      # color a standard normal array along each dim, in the drawn order [pair, pop, year, season, age, sex]
      draw <- array(stats::rnorm(prod(n_drawn)), dim = n_drawn)
      draw <- color_naa_dim(draw, chol_pop, 2)
      draw <- color_naa_dim(draw, chol_year, 3)
      draw <- color_naa_dim(draw, chol_season, 4)
      draw <- color_naa_dim(draw, chol_age, 5)
      draw <- color_naa_dim(draw, chol_sex, 6)
      draw <- sd_marginal * draw

      # an ar1 continuing from the fit: each year is rho times the year before plus an innovation of sd sqrt(1 - rho^2)
      if(continues) {
        previous <- last_fit
        for(y in 1:n_drawn_years) {
          draw[,,y,,,] <- rho_year * previous + sqrt(1 - rho_year^2) * draw[,,y,,,]
          previous <- array(draw[,,y,,,], dim = n_drawn[-3])
        } # end y loop
      }
      if(year_shared && last_cond_year > 0) for(y in 1:n_drawn_years) draw[,,y,,,] <- last_fit # one deviation shared across years continues unchanged
      draw[!estimated] <- 0 # a cell the map fixes stays at zero

      # every pair takes its own draw, and every active level reads its block's
      for(j in 1:n_pairs_drawn) {
        move_devs[active$pop, move_pairs[pairs_drawn[j], 1], move_pairs[pairs_drawn[j], 2], active_years, active$season, active$age, active$sex, sim] <-
          draw[j, take$pop, take_year, take$season, take$age, take$sex, drop = FALSE]
      } # end j loop
    } # end sim loop

  } # end b loop

  sim_env$move_devs <- move_devs
  derive_sim_movement(sim_env)
  return(invisible(NULL))

} # end function
