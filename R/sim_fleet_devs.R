# Operating model
#
# Fishing mortality, discard mortality and selectivity deviations drawn from the processes the fit
# penalizes, and the fleet arrays rebuilt from them through the estimation model's own functions.

#' The selectivity a fleet type draws deviations for, by prefix
#'
#' @keywords internal
sim_sel_prefixes <- c(fish = "fishsel", ret = "retsel", srv = "srvsel")

#' Arguments to Get_Selex_Array for one fleet type
#'
#' The same call the objective makes for fishery, retention or survey
#' selectivity, read off the fit's data and a parameter list.
#'
#' @param type \code{"fish"}, \code{"ret"} or \code{"srv"}.
#' @param data Data list of the fit.
#' @param pars Parameter list.
#'
#' @return Named list of \code{Get_Selex_Array} arguments.
#'
#' @keywords internal
sim_selex_args <- function(type, data, pars) {

  stem <- sim_sel_prefixes[[type]] # fishsel, retsel or srvsel
  selex_type <- data[[paste0(type, "_selex_type")]]
  list(selex_type = selex_type,
       bins = if(isTRUE(selex_type == 1)) data$lens else data$ages, # selectivity at length runs over length bins
       sel_blocks = data[[paste0(type, "_sel_blocks")]],
       sel_model = data[[paste0(type, "_sel_model")]],
       fixed_sel_pars = pars[[paste0(type, "_fixed_sel_pars")]],
       cont_tv_sel = data[[paste0("cont_tv_", type, "_sel")]],
       ln_seldevs = pars[[paste0("ln_", stem, "_devs")]],
       use_fixed_sel = data[[paste0("use_fixed_", type, "_sel")]],
       bin_devs = pars[[paste0("ln_", stem, "_bin_devs")]],
       bin_dev_bins = data[[paste0(type, "_sel_bin_dev_bins")]],
       sel_norm_bins = data[[paste0(type, "_sel_norm_bins")]],
       sex_par_offset = data[[paste0(stem, "_sex_par_offset")]],
       sex_scale_offset = data[[paste0(stem, "_sex_scale_offset")]],
       sex_apical_offset = data[[paste0(stem, "_sex_apical_offset")]],
       sex_scale = pars[[paste0("ln_", stem, "_sex_scale")]],
       nselbins = data[[paste0(type, "_sel_bicubic_nselbins")]],
       dbnrml_raw = data[[paste0(type, "_dbnrml_raw")]],
       dbnrml_startbin = data[[paste0(type, "_dbnrml_startbin")]],
       sel_input = data[[paste0(type, "_sel_input")]],
       bicubic_Wbin = data[[paste0(type, "_sel_bicubic_Wbin")]],
       bicubic_Wyr = data[[paste0(type, "_sel_bicubic_Wyr")]],
       bicubic_binnodes = data[[paste0(type, "_sel_bicubic_binnodes")]],
       bicubic_yrnodes = data[[paste0(type, "_sel_bicubic_yrnodes")]],
       n_pop = data$n_pop,
       n_regions = data$n_regions,
       n_yrs = length(data$years),
       n_proj_yrs_devs = if(is.null(data$n_proj_yrs_devs)) 0 else data$n_proj_yrs_devs,
       n_seas = data$n_seas,
       n_ages = length(data$ages),
       n_lens = length(data$lens),
       n_sexes = data$n_sexes,
       n_fleets = if(type == "srv") data$n_srv_fleets else data$n_fish_fleets)

} # end function

#' Fishing mortality, discard mortality and selectivity deviations in the operating model
#'
#' Stores what each replicate needs to draw these deviations and rebuild the
#' fleet arrays from them, and \code{\link{draw_sim_fleet_devs}} does so in
#' \code{Setup_sim_env}, keeping the fit's deviations over the simulation list's
#' \code{n_cond_yrs} and drawing the fitted years after them.
#'
#' Selectivity deviations are drawn wherever the fit varies selectivity, under
#' the process \code{\link{Get_PE_loglik}} penalizes, as growth and movement
#' deviations are. F and discard mortality deviations are drawn by default only
#' when the fit integrates them out (\code{random}): as fixed effects their
#' penalty is usually a loose constraint, \code{sigmaF = 1} by default, and
#' drawing from it would replace the fitted F with noise. \code{F_devs_draw =
#' "all"} draws them as fixed effects too, from the same penalty. Either way they
#' are kept at the fit when the penalty is off or centered on the deviations' own
#' mean, neither of which is a density to draw from.
#'
#' A penalty weighted by \code{Wt_F}, \code{Wt_D} or \code{*_pe_wt} is the
#' density of the same process with its variance divided by the weight, and
#' that is what is drawn; a self test keeps the weight in its refits.
#'
#' Selectivity at length (\code{*_selex_type = 1}) is drawn and rebuilt at
#' length, and its selectivity at age is the curve read through each
#' replicate's own size-age transition matrix, by the growth rebuild when the
#' operating model rebuilds growth and directly otherwise. A penalty weight
#' \code{w} (\code{*_pe_wt}) scales the penalty to the density of a process
#' whose variance is divided by \code{w}, which is what is drawn. When the
#' operating model runs past the fit, as a closed loop does, selectivity
#' deviations are drawn in those years too, on from the fitted ones (see
#' \code{\link{draw_sim_fleet_devs}}). F and discard mortality never are, since a
#' closed loop sets F in its projection years by its control rule.
#'
#' @param sim_list Simulation list with \code{n_yrs}, \code{n_sims},
#'   \code{Fmort}, \code{dmr} and the selectivity arrays already set.
#' @param data Data list of the fit.
#' @param pars Parameter list at the fitted values.
#' @param random Character vector of the fit's random effects. Default \code{NULL}.
#' @param pars_by_sim Optional list of one parameter list per replicate, for
#'   replicates that each run on their own parameter draw (\code{sim_type =
#'   "joint"} in \code{\link{simulation_self_test}}). \code{NULL} (default) gives
#'   every replicate \code{pars}.
#' @param F_devs_draw \code{"random"} (default) draws F and discard mortality
#'   deviations only when \code{random} integrates them out; \code{"all"} also
#'   draws them when the fit holds them as fixed effects, from their penalty.
#'
#' @return \code{sim_list} with \code{fleet_devs_on} (which of \code{F},
#'   \code{dmr}, \code{fish}, \code{ret} and \code{srv} are drawn),
#'   \code{fleet_dev_data} and \code{fleet_dev_pars} added, or unchanged when
#'   nothing is drawn.
#'
#' @export Setup_Sim_Fleet_Devs
#' @family Simulation Setup
Setup_Sim_Fleet_Devs <- function(sim_list, data, pars, random = NULL, pars_by_sim = NULL, F_devs_draw = c("random", "all")) {

  F_devs_draw <- match.arg(F_devs_draw)
  is_est <- function(map) !is.null(map) && any(!is.na(map)) # whether any deviation is estimated

  # F and discard mortality deviations, when integrated out or asked for as fixed effects, and only
  # where the penalty is a density
  F_pen_off <- !isTRUE(data$Use_F_pen == 1) || isTRUE(data$Wt_F == 0)
  F_on <- ("ln_F_devs" %in% random || F_devs_draw == "all") && is_est(data$map_ln_F_devs)
  if(F_on && (F_pen_off || isTRUE(data$Fdev_pen_center == 1))) {
    message("Setup_Sim_Fleet_Devs: ln_F_devs would be drawn, but its penalty is ", if(F_pen_off) "off" else "centered on the deviations' own mean",
            ", which is not a density to draw from, so F keeps the fit's deviations.")
    F_on <- FALSE
  }
  dmr_on <- ("logit_dmr_devs" %in% random || F_devs_draw == "all") && is_est(data$map_logit_dmr_devs)
  if(dmr_on && (!isTRUE(data$Use_dmr_pen == 1) || isTRUE(data$Wt_D == 0))) {
    message("Setup_Sim_Fleet_Devs: logit_dmr_devs would be drawn, but its penalty is off, so discard mortality keeps the fit's deviations.")
    dmr_on <- FALSE
  }

  # selectivity deviations wherever the fit varies selectivity
  sel_on <- c(fish = FALSE, ret = FALSE, srv = FALSE)
  for(type in names(sel_on)) {
    stem <- sim_sel_prefixes[[type]]
    varies <- any(data[[paste0("cont_tv_", type, "_sel")]] > 0) && is_est(data[[paste0("map_ln_", stem, "_devs")]])
    bin_varies <- any(data[[paste0("cont_tv_", stem, "_bin_devs")]] > 0) && is_est(data[[paste0("map_ln_", stem, "_bin_devs")]])
    sel_on[[type]] <- varies || bin_varies
  } # end type loop

  if(!F_on && !dmr_on && !any(sel_on)) return(sim_list)

  # the data each draw reads, and the parameters, one list per replicate under a joint self test
  sim_list$fleet_devs_on <- c(F = F_on, dmr = dmr_on, sel_on)
  sim_list$fleet_dev_data <- data
  par_list <- function(p) p[intersect(names(p), c("ln_F_mean", "ln_F_devs", "ln_sigmaF", "Fdev_rho", "logit_dmr_mean", "logit_dmr_devs", "ln_sigma_dmr",
                                                  unlist(lapply(names(sim_sel_prefixes), function(type) {
                                                    stem <- sim_sel_prefixes[[type]]
                                                    c(paste0(type, "_fixed_sel_pars"), paste0("ln_", stem, c("_devs", "_bin_devs", "_sex_scale")),
                                                      paste0(stem, c("_pe_pars", "_bin_devs_pe_pars")))
                                                  }))))]
  sim_list$fleet_dev_pars <- if(is.null(pars_by_sim)) list(par_list(pars)) else lapply(pars_by_sim, par_list)
  return(sim_list)

} # end function

#' One replicate's parameter list for the fleet deviations
#'
#' @param sim_env Simulation environment.
#' @param sim Replicate.
#'
#' @return The replicate's own under a joint self test, else the one every replicate shares.
#'
#' @keywords internal
sim_fleet_dev_pars <- function(sim_env, sim) {
  pars <- sim_env$fleet_dev_pars
  if(length(pars) == 1) pars[[1]] else pars[[sim]]
}

#' One draw of an F deviation series from the process the estimation model penalizes
#'
#' The reverse of \code{\link{Get_Fdev_PE_loglik}} for one region, season and
#' fleet. Only estimated years are stepped through, and a gap between two of
#' them is bridged at the elapsed number of years. A random walk's first year
#' keeps the replicate's own value, since the estimation model gives it a
#' diffuse \eqn{N(0, 5)} and leaves the F level to the data, as the recruitment
#' walk does in the operating model; an AR1's starts at its stationary sd.
#'
#' @param devs Numeric vector over years, the fit's deviations.
#' @param is_est Logical vector over years, where a deviation is estimated.
#' @param PE_model \code{1} iid, \code{2} random walk, \code{3} AR1.
#' @param sigma Process sd.
#' @param rho AR1 correlation on the natural scale.
#' @param n_cond Leading years that keep the fit's deviations.
#'
#' @return \code{devs} with the estimated years after \code{n_cond} drawn.
#'
#' @keywords internal
draw_F_dev_series <- function(devs, is_est, PE_model, sigma, rho, n_cond) {

  last_y <- NA # the previous estimated year
  for(y in which(is_est)) {
    if(y > n_cond) {
      if(PE_model == 1) devs[y] <- stats::rnorm(1, 0, sigma)
      if(PE_model == 2 && !is.na(last_y)) devs[y] <- stats::rnorm(1, devs[last_y], sigma * sqrt(y - last_y)) # the diffuse first year keeps its value
      if(PE_model == 3) {
        if(is.na(last_y)) devs[y] <- stats::rnorm(1, 0, sigma / sqrt(1 - rho^2)) # stationary
        else {
          gap <- y - last_y # elapsed years
          devs[y] <- stats::rnorm(1, rho^gap * devs[last_y], sigma * sqrt((1 - rho^(2 * gap)) / (1 - rho^2)))
        }
      } # end ar1
    } # end if drawn
    last_y <- y
  } # end y loop
  devs

} # end function

#' Process error parameters a weighted penalty draws at
#'
#' A penalty multiplied by \code{wt} is, up to a constant, the density of the
#' same process with its variance divided by \code{wt}: each log sd is lowered by
#' \code{log(wt) / 2} under iid, a random walk or the separable AR1, and the
#' 3D GMRF's log variance by \code{log(wt)}. The correlations are unchanged.
#'
#' @param pe_pars Array \code{[1, 1, slot, sex]} for one region and fleet.
#' @param PE_model Integer process code, as \code{Get_PE_loglik} reads it.
#' @param wt Penalty weight.
#'
#' @return \code{pe_pars} at that weight.
#'
#' @keywords internal
pe_pars_at_weight <- function(pe_pars, PE_model, wt) {
  if(wt == 1) return(pe_pars)
  if(PE_model %in% c(1, 2)) pe_pars[] <- pe_pars - 0.5 * log(wt) # every slot is a log sd
  if(PE_model == 5) pe_pars[,,4,] <- pe_pars[,,4,] - 0.5 * log(wt) # log sd
  if(PE_model %in% c(3, 4)) pe_pars[,,4,] <- pe_pars[,,4,] - log(wt) # log variance
  pe_pars
}

#' One draw of a fleet type's selectivity deviations
#'
#' Each region and fleet's surface is drawn by \code{\link{draw_growth_pe_surface}}
#' under that cell's process, the region in the place it reads as population. A
#' level shared over regions, sexes, bins or fleets then takes the value drawn
#' at its first cell, where the penalty reads its sd, which also fills the bins
#' the correlated forms leave out of their shared groups. A level shared over
#' units whose sds differ is drawn at the first unit's.
#'
#' @param PE_codes Integer array \code{[n_regions, n_fleets]} of process codes,
#'   or a vector \code{[n_fleets]} shared by every region.
#' @param map Array \code{[n_regions, n_yrs, n_bins, n_sexes, n_fleets]} of
#'   levels, \code{NA} where fixed.
#' @param pe_pars Array \code{[n_regions, slot, n_sexes, n_fleets]}.
#' @param fit_devs The fit's deviations, shaped as \code{map}.
#' @param bins Bins the correlated forms run over.
#' @param n_cond Leading years that keep the fit's deviations.
#' @param pe_wt Penalty weight by fleet; a fleet at zero has no penalty and keeps the fit's.
#' @param rw_init_sigma The sd the estimation model gives a walk's first year, by
#'   fleet, \code{NA} for the walk's own. \code{NULL} (default) is \code{NA} for every fleet.
#'
#' @return Array shaped as \code{map}.
#'
#' @keywords internal
draw_sel_dev_surface <- function(PE_codes, map, pe_pars, fit_devs, bins, n_cond, pe_wt = NULL, rw_init_sigma = NULL) {

  map_dim <- dim(map)
  n_fleets <- map_dim[5]
  if(is.null(dim(PE_codes))) PE_codes <- matrix(PE_codes, map_dim[1], n_fleets, byrow = TRUE) # one code per fleet, every region
  if(is.null(pe_wt)) pe_wt <- rep(1, n_fleets)
  devs <- fit_devs

  for(f in seq_len(n_fleets)) {
    if(pe_wt[f] == 0) next # no penalty, nothing to draw from
    for(r in seq_len(map_dim[1])) {
      if(PE_codes[r,f] == 0 || all(is.na(map[r,,,,f]))) next
      pe_rf <- pe_pars_at_weight(array(pe_pars[r,,,f], dim = c(1, 1, dim(pe_pars)[2:3])), PE_codes[r,f], pe_wt[f])
      drawn <- draw_growth_pe_surface(PE_model = PE_codes[r,f],
                                      map = array(map[r,,,,f], dim = c(1, 1, map_dim[2:4])),
                                      pe_pars = pe_rf,
                                      bins = bins,
                                      fit_devs = array(fit_devs[r,,,,f], dim = c(1, 1, map_dim[2:4])),
                                      n_cond = n_cond,
                                      rw_init_sigma = if(is.null(rw_init_sigma)) NA else rw_init_sigma[f])
      devs[r,,,,f] <- drawn[1,1,,,]
    } # end r loop
  } # end f loop

  # a shared level takes its first cell's value everywhere it appears, across fleets too under est_shared_f_x
  for(level in sort(unique(as.vector(map)))) {
    cells <- which(map == level)
    devs[cells] <- devs[cells[1]]
  } # end level loop

  devs

} # end function

#' Write one replicate's redrawn selectivity at length into the operating model
#'
#' The curve goes wherever the operating model reads selectivity at length: the
#' length compositions selected at length, and, when growth is rebuilt, the
#' curve \code{take_sim_growth_years} reads through each replicate's own keys
#' for its selectivity at age. Without a growth rebuild the selectivity at age
#' is formed here, the curve read through the keys the operating model was
#' given, as the estimation model forms it.
#'
#' @param sim_env Simulation environment.
#' @param type \code{"fish"}, \code{"ret"} or \code{"srv"}.
#' @param sel_l Selectivity at length \code{[region, year, len, sex, fleet]}.
#' @param yrs Fitted years to write.
#' @param sim Replicate.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
sim_sel_at_length <- function(sim_env, type, sel_l, yrs, sim) {

  sel_l_name <- paste0(type, "_sel_l")
  if(!is.null(sim_env[[sel_l_name]])) sim_env[[sel_l_name]][,yrs,,,,sim] <- sel_l[,yrs,,,,drop = FALSE] # length comps selected at length

  # a growth rebuild reads the curve through each replicate's own keys
  if(!is.null(sim_env$growth_args) || !is.null(sim_env$dsem_growth_args)) {
    if(is.null(sim_env$growth_length_sel_by_sim)) sim_env$growth_length_sel_by_sim <- rep(list(sim_env$growth_length_sel), sim_env$n_sims)
    sim_env$growth_length_sel_by_sim[[sim]][[sel_l_name]][,yrs,,,] <- sel_l[,yrs,,,,drop = FALSE]
    return(invisible(NULL))
  }

  # otherwise through the keys the operating model was given, the fleet's own or the shared one
  key_name <- if(type == "srv") "SizeAgeTrans_srv" else "SizeAgeTrans_fish"
  sel_name <- paste0(type, "_sel")
  for(p in seq_len(sim_env$n_pop)) for(r in seq_len(sim_env$n_regions)) for(y in yrs) for(seas in seq_len(sim_env$n_seas)) {
    for(s in seq_len(sim_env$n_sexes)) for(f in seq_len(dim(sel_l)[5])) {
      key <- if(!is.null(sim_env[[key_name]])) sim_env[[key_name]][p,r,y,seas,,,s,f,sim] else sim_env$SizeAgeTrans[p,r,y,seas,,,s,sim] # [len, age]
      sim_env[[sel_name]][p,r,y,seas,,s,f,sim] <- sel_l[r,y,,s,f] %*% key
    } # end s, f loop
  } # end p, r, y, seas loop

  return(invisible(NULL))

} # end function

#' Draw fleet deviations for every replicate and rebuild the fleet arrays
#'
#' F deviations under \code{Get_Fdev_PE_loglik}'s process, discard mortality
#' deviations iid, and selectivity deviations under \code{Get_PE_loglik}'s, each
#' kept at the fit over \code{n_cond_yrs} and drawn in the fitted years after
#' them. F and discard mortality are rebuilt in the cells whose deviation was
#' drawn, as \code{exp(ln_F_mean + ln_F_devs)} and
#' \code{plogis(logit_dmr_mean + logit_dmr_devs)}, and selectivity through
#' \code{Get_Selex_Array} over every fitted year. The draws are kept as
#' \code{ln_F_devs}, \code{logit_dmr_devs} and \code{ln_<type>sel_devs} with
#' \code{ln_<type>sel_bin_devs}, the replicate dim last.
#'
#' Years past the fit, a closed loop's projection, draw selectivity deviations on
#' from the last fitted year under the same process (\code{\link{extend_devs_past_fit}}),
#' and the projection years' selectivity is rebuilt over the longer span and
#' scaled so the fitted years keep the fit's standardization. Those draws are
#' kept as \code{ln_<type>sel_devs_proj}. A bicubic surface has no year weights
#' past the fit and keeps its last fitted year's selectivity.
#'
#' @param sim_env Simulation environment holding what
#'   \code{\link{Setup_Sim_Fleet_Devs}} stored, \code{n_yrs}, \code{n_sims} and
#'   \code{n_cond_yrs} (\code{NULL} or \code{0} draws every fitted year).
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
draw_sim_fleet_devs <- function(sim_env) {

  data <- sim_env$fleet_dev_data
  on <- sim_env$fleet_devs_on
  n_sims <- sim_env$n_sims
  n_yrs <- min(sim_env$n_yrs, length(data$years)) # fitted years the operating model runs
  n_cond <- min(n_yrs, if(is.null(sim_env$n_cond_yrs)) 0 else sim_env$n_cond_yrs) # years that keep the fit's deviations
  yrs <- seq_len(n_yrs)

  for(sim in seq_len(n_sims)) {

    pars <- sim_fleet_dev_pars(sim_env, sim)

    # F deviations, each region, season and fleet's estimated years in sequence
    if(on[["F"]]) {
      if(sim == 1) sim_env$ln_F_devs <- array(0, dim = c(dim(pars$ln_F_devs), n_sims))
      F_devs <- pars$ln_F_devs
      is_est <- !is.na(data$map_ln_F_devs)
      F_wt <- if(is.null(data$Wt_F)) 1 else data$Wt_F # the penalty's weight, a density at the sd over its root
      for(r in seq_len(dim(F_devs)[1])) for(seas in seq_len(dim(F_devs)[3])) for(f in seq_len(dim(F_devs)[4])) {
        rho <- if(data$Fdev_model == 3) rho_trans(pars$Fdev_rho[r,seas,f]) else 0
        F_devs[r,yrs,seas,f] <- draw_F_dev_series(F_devs[r,yrs,seas,f], is_est[r,yrs,seas,f], data$Fdev_model, exp(pars$ln_sigmaF[r,seas,f]) / sqrt(F_wt), rho, n_cond)
        for(y in which(is_est[r,yrs,seas,f] & yrs > n_cond)) sim_env$Fmort[r,y,seas,f,sim] <- exp(pars$ln_F_mean[r,seas,f] + F_devs[r,y,seas,f])
      } # end r, seas, f loop
      sim_env$ln_F_devs[,,,,sim] <- F_devs
    } # end F deviations

    # discard mortality deviations, iid
    if(on[["dmr"]]) {
      if(sim == 1) sim_env$logit_dmr_devs <- array(0, dim = c(dim(pars$logit_dmr_devs), n_sims))
      dmr_devs <- pars$logit_dmr_devs
      dmr_wt <- if(is.null(data$Wt_D)) 1 else data$Wt_D
      drawn <- !is.na(data$map_logit_dmr_devs) & slice.index(dmr_devs, 2) > n_cond & slice.index(dmr_devs, 2) <= n_yrs
      for(cell in which(drawn)) {
        idx <- arrayInd(cell, dim(dmr_devs)) # region, year, season, fleet
        dmr_devs[cell] <- stats::rnorm(1, 0, exp(pars$ln_sigma_dmr[idx[1],idx[3],idx[4]]) / sqrt(dmr_wt))
        sim_env$dmr[idx[1],idx[2],idx[3],idx[4],sim] <- stats::plogis(pars$logit_dmr_mean[idx[1],idx[3],idx[4]] + dmr_devs[cell])
      } # end cell loop
      sim_env$logit_dmr_devs[,,,,sim] <- dmr_devs
    } # end discard mortality deviations

    # selectivity, each fleet type's deviations then the array rebuilt from them
    for(type in names(sim_sel_prefixes)) {
      if(!on[[type]]) next
      stem <- sim_sel_prefixes[[type]]
      dev_name <- paste0("ln_", stem, "_devs")
      bin_name <- paste0("ln_", stem, "_bin_devs")
      if(sim == 1) {
        sim_env[[dev_name]] <- array(0, dim = c(dim(pars[[dev_name]]), n_sims))
        sim_env[[bin_name]] <- array(0, dim = c(dim(pars[[bin_name]]), n_sims))
      }
      args <- sim_selex_args(type, data, pars)
      args$ln_seldevs <- draw_sel_dev_surface(PE_codes = data[[paste0("cont_tv_", type, "_sel")]],
                                              map = data[[paste0("map_", dev_name)]],
                                              pe_pars = pars[[paste0(stem, "_pe_pars")]],
                                              fit_devs = pars[[dev_name]],
                                              bins = data[[paste0(stem, "_devs_min_shared_bins")]],
                                              n_cond = n_cond,
                                              pe_wt = data[[paste0(stem, "_pe_wt")]],
                                              rw_init_sigma = data[[paste0(stem, "_rw_init_sigma")]])
      args$bin_devs <- draw_sel_dev_surface(PE_codes = data[[paste0("cont_tv_", stem, "_bin_devs")]],
                                            map = data[[paste0("map_", bin_name)]],
                                            pe_pars = pars[[paste0(stem, "_bin_devs_pe_pars")]],
                                            fit_devs = pars[[bin_name]],
                                            bins = seq_len(dim(pars[[bin_name]])[3]),
                                            n_cond = n_cond,
                                            rw_init_sigma = data[[paste0(stem, "_bin_devs_rw_init_sigma")]])
      selex <- do.call(Get_Selex_Array, args)
      if(isTRUE(data[[paste0(type, "_selex_type")]] == 1)) sim_sel_at_length(sim_env, type, selex$sel_l, yrs, sim)
      else sim_env[[paste0(type, "_sel")]][,,yrs,,,,,sim] <- selex$sel[,,yrs,,,,,drop = FALSE]
      sim_env[[dev_name]][,,,,,sim] <- args$ln_seldevs
      sim_env[[bin_name]][,,,,,sim] <- args$bin_devs

      # years past the fit, drawn on from the fitted ones; a bicubic surface has no year weights there
      n_proj_om <- sim_env$n_yrs - length(data$years)
      if(n_proj_om > 0 && !any(data[[paste0(type, "_sel_model")]] == 8)) {
        draw_sel_devs_past_fit(sim_env, type, data, pars, args, selex, n_proj_om, sim)
      }
    } # end type loop

  } # end sim loop

  return(invisible(NULL))

} # end function

#' A fleet's deviations and their levels extended past the fit
#'
#' The fitted years keep the deviations given and their levels. Each year past
#' the fit starts at zero and takes the sharing of the last fitted year that
#' estimates any deviation, under levels of its own, so
#' \code{\link{draw_sel_dev_surface}} draws it fresh rather than copying a fitted
#' year's value.
#'
#' @param devs Deviations \code{[region, year, bin, sex, fleet]} over at least the
#'   fitted years.
#' @param map Their levels, shaped as \code{devs}, \code{NA} where fixed.
#' @param n_fit Fitted years.
#' @param n_total Years the operating model runs.
#'
#' @return List of \code{devs} and \code{map} over \code{n_total} years.
#'
#' @keywords internal
extend_devs_past_fit <- function(devs, map, n_fit, n_total) {

  fit_yrs <- seq_len(n_fit)
  ext_dim <- dim(devs)
  ext_dim[2] <- n_total
  ext_devs <- array(0, dim = ext_dim)
  ext_map <- array(NA, dim = ext_dim)
  ext_devs[,fit_yrs,,,] <- devs[,fit_yrs,,,,drop = FALSE]
  if(is.null(map) || all(is.na(map))) return(list(devs = ext_devs, map = ext_map)) # nothing estimated, nothing drawn
  ext_map[,fit_yrs,,,] <- map[,fit_yrs,,,,drop = FALSE]

  # the last fitted year with an estimated deviation sets the sharing every projection year copies
  est_yrs <- which(apply(!is.na(map[,fit_yrs,,,,drop = FALSE]), 2, any))
  if(!length(est_yrs)) return(list(devs = ext_devs, map = ext_map))
  last_pattern <- map[,max(est_yrs),,,,drop = FALSE]
  level_step <- max(map, na.rm = TRUE) # keeps each projection year's levels apart from the fit's and from each other's
  for(y in n_fit + seq_len(n_total - n_fit)) ext_map[,y,,,] <- last_pattern + level_step * (y - n_fit)

  list(devs = ext_devs, map = ext_map)

} # end extend_devs_past_fit

#' Draw one replicate's selectivity deviations past the fit and rebuild those years
#'
#' The deviations are drawn on from the replicate's fitted ones under the fitted
#' process and rebuilt with \code{Get_Selex_Array} over the fitted and
#' projection years, the projection years reading the terminal year's blocks and
#' model as the estimation model's projection does. A form centered over all its
#' years would put the projection years in the mean and shift every fitted year,
#' so each region, sex and fleet's projection years are scaled by what takes the
#' longer build's terminal fitted year back to the fit's.
#'
#' @param sim_env Simulation environment.
#' @param type \code{"fish"}, \code{"ret"} or \code{"srv"}.
#' @param data,pars The fit's data list and this replicate's parameters.
#' @param args The \code{Get_Selex_Array} arguments the fitted years were rebuilt
#'   with, their deviations the replicate's.
#' @param selex That rebuild's result.
#' @param n_proj_om Years the operating model runs past the fit.
#' @param sim Replicate.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
draw_sel_devs_past_fit <- function(sim_env, type, data, pars, args, selex, n_proj_om, sim) {

  stem <- sim_sel_prefixes[[type]]
  dev_name <- paste0("ln_", stem, "_devs")
  bin_name <- paste0("ln_", stem, "_bin_devs")
  n_fit <- length(data$years)
  n_total <- n_fit + n_proj_om
  proj_yrs <- n_fit + seq_len(n_proj_om)

  # the deviations, on from the replicate's fitted ones
  sel_ext <- extend_devs_past_fit(args$ln_seldevs, data[[paste0("map_", dev_name)]], n_fit, n_total)
  proj_args <- args
  proj_args$n_proj_yrs_devs <- n_proj_om # Get_Selex_Array reads the terminal year's blocks past the fit
  proj_args$ln_seldevs <- draw_sel_dev_surface(PE_codes = data[[paste0("cont_tv_", type, "_sel")]],
                                               map = sel_ext$map,
                                               pe_pars = pars[[paste0(stem, "_pe_pars")]],
                                               fit_devs = sel_ext$devs,
                                               bins = data[[paste0(stem, "_devs_min_shared_bins")]],
                                               n_cond = n_fit,
                                               pe_wt = data[[paste0(stem, "_pe_wt")]])
  if(!is.null(args$bin_devs)) {
    bin_ext <- extend_devs_past_fit(args$bin_devs, data[[paste0("map_", bin_name)]], n_fit, n_total)
    proj_args$bin_devs <- draw_sel_dev_surface(PE_codes = data[[paste0("cont_tv_", stem, "_bin_devs")]],
                                               map = bin_ext$map,
                                               pe_pars = pars[[paste0(stem, "_bin_devs_pe_pars")]],
                                               fit_devs = bin_ext$devs,
                                               bins = seq_len(dim(args$bin_devs)[3]),
                                               n_cond = n_fit)
  }
  proj <- do.call(Get_Selex_Array, proj_args)

  # keep the projection draws
  proj_name <- paste0(dev_name, "_proj")
  dev_dim <- dim(proj_args$ln_seldevs) # region, year, bin, sex, fleet
  if(is.null(sim_env[[proj_name]])) sim_env[[proj_name]] <- array(0, dim = c(dev_dim[1], n_proj_om, dev_dim[3:5], sim_env$n_sims))
  sim_env[[proj_name]][,,,,,sim] <- proj_args$ln_seldevs[,proj_yrs,,,,drop = FALSE]

  # scale each region, sex and fleet's projection years to the fit's standardization, then write them
  if(isTRUE(data[[paste0(type, "_selex_type")]] == 1)) {
    sel_l <- proj$sel_l # [region, year, len, sex, fleet]
    for(r in seq_len(dim(sel_l)[1])) for(s in seq_len(dim(sel_l)[4])) for(f in seq_len(dim(sel_l)[5])) {
      to_fit <- sum(selex$sel_l[r,n_fit,,s,f]) / sum(sel_l[r,n_fit,,s,f]) # 1 unless the form is centered over its years
      if(is.finite(to_fit)) sel_l[r,proj_yrs,,s,f] <- sel_l[r,proj_yrs,,s,f] * to_fit
    } # end r, s, f loop
    sim_sel_at_length(sim_env, type, sel_l, proj_yrs, sim)
  } else {
    sel <- proj$sel # [pop, region, year, season, age, sex, fleet]
    for(p in seq_len(dim(sel)[1])) for(r in seq_len(dim(sel)[2])) for(seas in seq_len(dim(sel)[4])) {
      for(s in seq_len(dim(sel)[6])) for(f in seq_len(dim(sel)[7])) {
        to_fit <- sum(selex$sel[p,r,n_fit,seas,,s,f]) / sum(sel[p,r,n_fit,seas,,s,f]) # 1 unless the form is centered over its years
        if(is.finite(to_fit)) sel[p,r,proj_yrs,seas,,s,f] <- sel[p,r,proj_yrs,seas,,s,f] * to_fit
      } # end s, f loop
    } # end p, r, seas loop
    sim_env[[paste0(type, "_sel")]][,,proj_yrs,,,,,sim] <- sel[,,proj_yrs,,,,,drop = FALSE]
  }

  return(invisible(NULL))

} # end draw_sel_devs_past_fit
