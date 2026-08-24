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

  # Movement / mortality sequencing. Simulation lists assembled directly through the
  # Setup_Sim_* helpers never set this, so default it here rather than at every use.
  if(is.null(sim_list$move_timing)) sim_list$move_timing <- 0

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

  # output into simulation environment
  list2env(sim_list, envir = sim_env)

  # A NULL element is dropped by list2env, but the movement-timing branches reference
  # Mrate unconditionally. Bind it explicitly so it resolves to NULL (which the
  # move_timing 0 and 1 operators ignore) instead of escaping to the enclosing frame.
  if(is.null(sim_env$Mrate)) sim_env$Mrate <- NULL
  if(sim_env$move_timing == 2 && is.null(sim_env$Mrate))
    stop("move_timing == 2 (continuous movement) requires Mrate, the instantaneous rate matrix, in the simulation list.")

  return(sim_env)
}
