#' Initialize Simulation Dimension Settings
#'
#' Creates and returns a list of key dimension values used to set up a
#' simulation or management strategy evaluation (MSE). This list provides
#' structural information such as number of simulations, years, regions,
#' ages, fleets, and whether to include a feedback loop.
#'
#' @param n_sims Integer. Number of simulation replicates.
#' @param n_yrs Integer. Number of years in the simulation.
#' @param n_regions Integer. Number of modeled regions.
#' @param n_ages Integer. Number of modeled age classes.
#' @param n_lens Integer. Number of modeled length bins.
#' @param n_obs_ages Integer. Number of observed age classes (can differ from \code{n_ages}, default = \code{n_ages}).
#' @param n_sexes Integer. Number of sexes.
#' @param n_fish_fleets Integer. Number of fishery fleets.
#' @param n_srv_fleets Integer. Number of survey fleets.
#' @param run_feedback Logical. Whether to include a feedback management loop (default = \code{FALSE}).
#' @param feedback_start_yr Integer. First year that feedback is applied (only used if \code{run_feedback = TRUE}).
#' @param n_seas Integer. Number of seasons
#' @param seasdur Duration of season. Default is 1 if 1 season, or 1 / n_seas if n_seas > 1
#'
#' @return
#' A list containing the specified dimension values, with elements:
#' \itemize{
#'   \item \code{n_sims}, \code{n_yrs}, \code{n_regions}, \code{n_ages}, \code{n_lens},
#'   \code{n_obs_ages}, \code{n_sexes}, \code{n_fish_fleets}, \code{n_srv_fleets}
#'   \item \code{init_iter} (set internally to \code{n_ages * 10})
#'   \item \code{feedback_start_yr}, \code{run_feedback}
#' }
#'
#' @export Setup_Sim_Dim
#' @family Simulation Setup
Setup_Sim_Dim <- function(n_sims,
                          n_yrs,
                          n_seas = 1,
                          n_regions,
                          n_ages,
                          n_lens,
                          n_obs_ages = n_ages,
                          n_sexes,
                          n_fish_fleets,
                          n_srv_fleets,
                          seasdur = if(n_seas == 1) 1 else rep(1 / n_seas, n_seas),
                          run_feedback = FALSE,
                          feedback_start_yr = NULL
                          ) {

  sim_list <- list() # setup empty list

  # output dimensions into list
  sim_list$n_sims <- n_sims
  sim_list$n_yrs <- n_yrs
  sim_list$n_seas <- n_seas
  sim_list$n_regions <- n_regions
  sim_list$n_ages <- n_ages
  sim_list$n_lens <- n_lens
  sim_list$n_obs_ages <- n_obs_ages
  sim_list$n_sexes <- n_sexes
  sim_list$n_fish_fleets <- n_fish_fleets
  sim_list$n_srv_fleets <- n_srv_fleets
  sim_list$init_iter <- n_ages * 10
  sim_list$seasdur <- seasdur
  sim_list$feedback_start_yr <- feedback_start_yr
  sim_list$run_feedback <- run_feedback

  return(sim_list)

}

#' Set up model dimensions
#'
#' @param n_regions Integer specifying the number of spatial regions.
#' @param ages Numeric vector of age classes.
#' @param n_sexes Integer specifying the number of sexes.
#' @param n_fish_fleets Integer specifying the number of fishery fleets.
#' @param n_srv_fleets Integer specifying the number of survey fleets.
#' @param years Numeric vector of years.
#' @param lens Numeric vector of length bins; can be set to \code{1} if length data are not modeled.
#' @param verbose Logical flag indicating whether to print progress messages (default \code{FALSE}).
#' @param n_proj_yrs_devs Number of projection years for deviation parameters (ln_RecDevs, move_devs, ln_fishsel_devs, ln_srvsel_devs)
#'
#' @returns A list containing three named elements:
#' \describe{
#'   \item{\code{data}}{List of data inputs dimensioned by the model dimensions.}
#'   \item{\code{parameters}}{List of model parameters initialized according to dimensions.}
#'   \item{\code{map}}{List of parameter mappings for model fitting.}
#' }
#' @export Setup_Mod_Dim
#' @family Model Setup
Setup_Mod_Dim <- function(years,
                          ages,
                          lens,
                          n_regions,
                          n_sexes,
                          n_fish_fleets,
                          n_srv_fleets,
                          n_proj_yrs_devs = 0,
                          verbose = FALSE
                          ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Create empty list
  input_list <- list(data = list(), par = list(), map = list())

  # ouput variables into list
  input_list$data$years <- years
  input_list$data$n_regions <- n_regions
  input_list$data$ages <- ages
  input_list$data$lens <- if(is.null(lens)) 1 else lens
  input_list$data$n_sexes <- n_sexes
  input_list$data$n_fish_fleets <- n_fish_fleets
  input_list$data$n_srv_fleets <- n_srv_fleets
  input_list$data$n_proj_yrs_devs <- n_proj_yrs_devs
  input_list$verbose <- verbose

  collect_message("Number of Years: ", length(years))
  collect_message("Number of Projection Years for Dev Pars: ", n_proj_yrs_devs)
  collect_message("Number of Regions: ", n_regions)
  collect_message("Number of Age Bins: ", length(ages))
  collect_message("Number of Length Bins: ", length(lens))
  collect_message("Number of Sexes: ", n_sexes)
  collect_message("Number of Fishery Fleets: ", n_fish_fleets)
  collect_message("Number of Survey Fleets: ", n_srv_fleets)

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}
