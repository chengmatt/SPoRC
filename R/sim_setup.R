# Operating model
#
# Builds the environment the operating model runs inside. Holding simulation
# state in a dedicated environment is what lets the annual cycle helpers update
# shared arrays without threading them through every call.

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
  # Bound explicitly like the helpers above: the annual-cycle with() blocks resolve
  # functions through this environment's parent chain, which only reaches the package
  # namespace when the environment was built inside a package function.
  sim_env$draw_index_obs <- draw_index_obs
  sim_env$resolve_idx_factor <- resolve_idx_factor
  sim_env$draw_naa_innovations <- draw_naa_innovations
  sim_env$color_naa_margin <- color_naa_margin
  sim_env$Get_3d_precision <- Get_3d_precision

  # state-space numbers at age; lists built before the option existed carry none of it
  if(is.null(sim_list$NAA_re)) sim_list$NAA_re <- 0
  if(is.null(sim_list$sigmaNAA)) sim_list$sigmaNAA <- 0
  if(is.null(sim_list$naa_rho)) sim_list$naa_rho <- c(age = 0, year = 0, cohort = 0)
  for(nm in c("NAA_re_pop", "NAA_re_region", "NAA_re_sex")) if(is.null(sim_list[[nm]])) sim_list[[nm]] <- 0
  for(nm in c("naa_pop_corr", "naa_region_corr", "naa_sex_corr")) if(is.null(sim_list[[nm]])) sim_list[[nm]] <- 0
  # isTRUE rather than a bare comparison: a minimal list assembled for one component carries no
  # dimensions, and NULL > 1 is logical(0), which if() rejects rather than treating as false
  if(is.null(sim_list$naa_re_ages)) sim_list$naa_re_ages <- if(isTRUE(sim_list$n_ages > 1)) 2:sim_list$n_ages else integer(0)
  if(is.null(sim_list$naa_re_yrs)) sim_list$naa_re_yrs <- if(isTRUE(sim_list$n_yrs > 1)) 2:sim_list$n_yrs else integer(0)

  # output into simulation environment
  list2env(sim_list, envir = sim_env)

  # State-space containers live here rather than in Simulate_Pop_Static 
  naa_dims <- c(sim_env$n_pop, sim_env$n_regions, sim_env$n_yrs, sim_env$n_ages, sim_env$n_sexes, sim_env$n_sims)
  if(length(naa_dims) == 6 && all(is.finite(naa_dims))) {
    sim_env$naa_eta_all <- array(0, dim = naa_dims)
    sim_env$NAA_pred <- array(0, dim = naa_dims)
  }

  # A NULL element is dropped by list2env, but the movement-timing branches reference
  # Mrate unconditionally. Bind it explicitly so it resolves to NULL (which the
  # move_timing 0 and 1 operators ignore) instead of escaping to the enclosing frame.
  if(is.null(sim_env$Mrate)) sim_env$Mrate <- NULL
  if(sim_env$move_timing == 2 && is.null(sim_env$Mrate))
    stop("move_timing == 2 (continuous movement) requires Mrate, the instantaneous rate matrix, in the simulation list.")

  return(sim_env)
}
