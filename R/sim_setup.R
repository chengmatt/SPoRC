# Operating model
#
# Builds the environment the operating model runs inside. Keeping simulation state in its own
# environment lets the annual cycle helpers update shared arrays without threading them through.

#' Construct and populate a simulation execution environment
#'
#' Creates a new R environment populated with all objects from \code{sim_list}
#' and binds the SPoRC simulation functions required by
#' \code{\link{run_annual_cycle}}. Isolating the simulation state in a
#' dedicated environment prevents name collisions with the calling frame and
#' allows \code{with()} / \code{<<-} assignment patterns used internally by
#' the annual-cycle helpers to modify shared state without polluting the
#' global workspace.
#'
#' @param sim_list Named list returned by \code{\link{Setup_Sim_Rec}} (or the
#'   last upstream setup function called). All elements are copied into the
#'   new environment via \code{list2env}.
#'
#' @return A new environment (parent = calling frame) containing every element
#'   of \code{sim_list} as a named object, plus bound references to the
#'   following SPoRC simulation functions: \code{generate_initial_age_structure},
#'   \code{generate_recruitment}, \code{apply_pop_dy}, \code{compute_biom_y_sim},
#'   \code{generate_fishery_catch_comp_idx}, \code{generate_survey_comp_idx},
#'   \code{release_conv_tags}, \code{generate_fishery_conv_tags_recap},
#'   \code{Get_Det_Recruitment}, \code{Get_Init_NAA},
#'   \code{predict_sim_fish_iss_fmort}, \code{rho_trans},
#'   \code{simulate_comps}, \code{simulate_conv_tag_fish_recaptures},
#'   \code{draw_index_obs}, \code{resolve_idx_factor}.
#'
#'
#' @export Setup_sim_env
#' @family Simulation Setup
#'
#' @examples
#' \dontrun{
#' sim_env <- Setup_sim_env(sim_list)
#' }
Setup_sim_env <- function(sim_list) {

  # Guard movement stuff
  if(is.null(sim_list$move_timing)) sim_list$move_timing <- 0
  if(is.null(sim_list$expm_nsub)) sim_list$expm_nsub <- 0

  sim_env <- new.env(parent = parent.frame()) # define new environment for simulation

  # Get SPoRC functions in simulation environment
  sim_env$generate_initial_age_structure <- generate_initial_age_structure
  sim_env$generate_recruitment <- generate_recruitment
  sim_env$apply_pop_dy <- apply_pop_dy
  sim_env$compute_biom_y_sim <- compute_biom_y_sim
  sim_env$generate_fishery_catch_comp_idx <- generate_fishery_catch_comp_idx
  sim_env$generate_survey_comp_idx <- generate_survey_comp_idx
  sim_env$release_conv_tags <- release_conv_tags
  sim_env$generate_fishery_conv_tags_recap <- generate_fishery_conv_tags_recap
  sim_env$Get_Det_Recruitment <- Get_Det_Recruitment
  sim_env$Get_Init_NAA <- Get_Init_NAA
  sim_env$predict_sim_fish_iss_fmort <- predict_sim_fish_iss_fmort
  sim_env$rho_trans <- rho_trans
  sim_env$simulate_comps <- simulate_comps
  sim_env$simulate_caal <- simulate_caal
  sim_env$simulate_conv_tag_fish_recaptures <- simulate_conv_tag_fish_recaptures
  # bound explicitly like the helpers above: the annual-cycle with() blocks resolve functions
  # through this environment's parent chain, which only reaches the namespace inside a package call
  sim_env$draw_index_obs <- draw_index_obs
  sim_env$resolve_idx_factor <- resolve_idx_factor
  sim_env$draw_naa_innovations <- draw_naa_innovations
  sim_env$color_naa_dim <- color_naa_dim
  sim_env$Get_3d_precision <- Get_3d_precision

  # state-space numbers at age; lists built before the option existed have none of it
  if(is.null(sim_list$NAA_re)) sim_list$NAA_re <- 0
  if(is.null(sim_list$sigmaNAA)) sim_list$sigmaNAA <- 0
  if(is.null(sim_list$naa_rho)) sim_list$naa_rho <- c(age = 0, year = 0, cohort = 0)
  for(opt_name in c("NAA_re_pop", "NAA_re_region", "NAA_re_sex", "NAA_re_season")) if(is.null(sim_list[[opt_name]])) sim_list[[opt_name]] <- 0
  for(opt_name in c("naa_pop_corr", "naa_region_corr", "naa_sex_corr", "naa_season_corr")) if(is.null(sim_list[[opt_name]])) sim_list[[opt_name]] <- 0
  if(is.null(sim_list$naa_re_ages)) sim_list$naa_re_ages <- if(isTRUE(sim_list$n_ages > 1)) 2:sim_list$n_ages else integer(0)
  if(is.null(sim_list$naa_re_yrs)) sim_list$naa_re_yrs <- if(isTRUE(sim_list$n_yrs > 1)) 2:sim_list$n_yrs else integer(0)
  if(is.null(sim_list$naa_re_seas)) sim_list$naa_re_seas <- 1L

  # recruitment deviation process error; lists built before the option existed drew independently
  if(is.null(sim_list$RecDevs_model)) sim_list$RecDevs_model <- 1
  if(is.null(sim_list$rec_bias_correct)) sim_list$rec_bias_correct <- 1
  if(is.null(sim_list$RecDevs_rho)) sim_list$RecDevs_rho <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions))

  # output into simulation environment
  list2env(sim_list, envir = sim_env)

  # State-space containers live here rather than in Simulate_Pop_Static
  naa_dims <- c(sim_env$n_pop, sim_env$n_regions, sim_env$n_yrs, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_sims)
  if(length(naa_dims) == 7 && all(is.finite(naa_dims))) {
    sim_env$naa_eta_all <- array(0, dim = naa_dims)
    sim_env$NAA_pred <- array(0, dim = naa_dims)
  }

  # catchability deviations; lists built before the option existed have none, so every fleet is "none"
  for(prefix in c("fish", "srv")) {
    n_fleets <- sim_env[[paste0("n_", prefix, "_fleets")]]
    if(is.null(n_fleets) || is.null(sim_env[[paste0(prefix, "_q")]])) next
    if(is.null(sim_env[[paste0(prefix, "_q_model")]])) sim_env[[paste0(prefix, "_q_model")]] <- rep(1L, n_fleets)
    if(is.null(sim_env[[paste0("sigma_", prefix, "_q")]])) sim_env[[paste0("sigma_", prefix, "_q")]] <- array(0, dim = c(sim_env$n_regions, n_fleets))
    if(is.null(sim_env[[paste0(prefix, "_q_rho")]])) sim_env[[paste0(prefix, "_q_rho")]] <- array(0, dim = c(sim_env$n_regions, n_fleets))

    # a self test or closed loop supplies the fit's own deviations here; everything else starts at zero
    devs_name <- paste0("ln_", prefix, "_q_devs")
    q_dim <- dim(sim_env[[paste0(prefix, "_q")]])
    if(is.null(sim_env[[devs_name]])) sim_env[[devs_name]] <- array(0, dim = q_dim)
    else if(!identical(dim(sim_env[[devs_name]]), q_dim))
      stop(devs_name, " is [", paste(dim(sim_env[[devs_name]]), collapse = ", "), "] and ", prefix, "_q is [",
           paste(q_dim, collapse = ", "), "]. The deviations must sit on the same region, year, fleet and ",
           "replicate grid as the catchability they scale.")

  } # end prefix loop

  # dsem stuff
  sim_env$dsem_drawn <- list()
  for(entry in dsem_process_table()) {
    if(is.null(sim_env[[entry$sim_par]])) next
    dims <- dim(sim_env[[entry$sim_par]])
    sim_env$dsem_drawn[[entry$sim_par]] <- array(FALSE, dim = dims[-length(dims)]) # every dim but the replicate
  } # end entry loop

  # the conditioning years reproduce the fit's reported catchability, so a dsem does not rewrite them
  q_devs_cond <- stats::setNames(lapply(q_dev_par_names(), function(nm) sim_env[[nm]]), q_dev_par_names())
  if(!is.null(sim_env$dsem_model)) draw_dsem_sim(sim_env)
  n_cond_q <- if(is.null(sim_env$n_cond_yrs)) 0L else as.integer(sim_env$n_cond_yrs)
  for(nm in names(q_devs_cond)) {
    if(is.null(q_devs_cond[[nm]]) || n_cond_q < 1) next
    cond_yr <- seq_len(min(n_cond_q, dim(q_devs_cond[[nm]])[2]))
    sim_env[[nm]][,cond_yr,,] <- q_devs_cond[[nm]][,cond_yr,,,drop = FALSE]
  } # end nm loop

  check_q_dsem_drawable(sim_env) # check to see if dsem can do draws

  # movement stuff
  if(is.null(sim_env$Mrate)) sim_env$Mrate <- NULL
  if(sim_env$move_timing == 2 && is.null(sim_env$Mrate)) stop("move_timing == 2 (continuous movement) requires Mrate, the instantaneous rate matrix, in the simulation list.")

  return(sim_env)
}
