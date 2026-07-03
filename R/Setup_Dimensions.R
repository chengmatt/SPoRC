#' Initialise simulation dimension settings
#'
#' Creates the foundational \code{sim_list} object used throughout a closed-loop
#' simulation or static simulation. All downstream setup
#' functions (\code{\link{Setup_Sim_Biologicals}}, \code{\link{Setup_Sim_Containers}},
#' etc.) expect a \code{sim_list} produced by this function. Dimension scalars are
#' validated on input and stored alongside derived quantities such as
#' \code{init_iter} and the \code{natal_region} mapping.
#'
#' @param n_sims Positive integer. Number of simulation replicates.
#' @param n_yrs Positive integer. Number of projection years.
#' @param n_pop Positive integer. Number of distinct populations (default \code{1}).
#'   Populations share movement and selectivity schedules but can have independent
#'   stock-recruit relationships and natal regions.
#' @param n_regions Positive integer. Number of spatial regions.
#' @param natal_region Integer vector of length \code{n_pop} mapping each population
#'   to its natal region (1-indexed). Defaults are applied when \code{NULL}:
#'   \itemize{
#'     \item \code{n_regions == 1}: all populations assigned to region 1.
#'     \item \code{n_pop == n_regions}: one-to-one mapping (\code{1:n_pop}).
#'     \item \code{n_pop == 1}: assigned to region 1.
#'   }
#'   Must be supplied explicitly when \code{n_pop > 1}, \code{n_pop != n_regions},
#'   and \code{n_regions > 1}, as multiple populations would share a natal region
#'   and no sensible default exists.
#' @param n_seas Positive integer. Number of seasons within each year (default \code{1}).
#' @param seasdur Numeric vector of length \code{n_seas} giving the duration of
#'   each season as a fraction of a year. Values must sum to 1. Defaults to
#'   \code{1} when \code{n_seas == 1}, or \code{rep(1 / n_seas, n_seas)} for
#'   equal-length seasons otherwise.
#' @param n_ages Positive integer. Number of modelled age classes.
#' @param n_obs_ages Positive integer. Number of observed age bins in composition
#'   data. Can differ from \code{n_ages} when the plus group or youngest ages are
#'   pooled differently in observations. Defaults to \code{n_ages}.
#' @param n_lens Positive integer. Number of length bins. Set to \code{NULL}
#'   (default) when length compositions are not simulated.
#' @param n_sexes Integer. Number of sexes; must be either \code{1}
#'   (sex-aggregated) or \code{2} (sex-structured).
#' @param n_fish_fleets Positive integer. Number of fishery fleets.
#' @param n_srv_fleets Positive integer. Number of survey fleets.
#' @param run_feedback Logical. Whether to run a closed-loop feedback MSE in
#'   which an estimation model and harvest control rule update fishing mortality
#'   each year. Default \code{FALSE} (open-loop simulation).
#' @param feedback_start_yr Integer. First year in which feedback is applied.
#'   Required when \code{run_feedback = TRUE}; ignored otherwise.
#'
#' @return A named list (\code{sim_list}) containing:
#'   \describe{
#'     \item{\code{n_sims}, \code{n_yrs}, \code{n_pop}, \code{n_regions},
#'       \code{n_seas}, \code{n_ages}, \code{n_obs_ages}, \code{n_lens},
#'       \code{n_sexes}, \code{n_fish_fleets}, \code{n_srv_fleets}}{
#'       Dimension scalars supplied as inputs.}
#'     \item{\code{natal_region}}{Integer vector of length \code{n_pop}
#'       defining natal regions.}
#'     \item{\code{seasdur}}{Season durations summing to 1.}
#'     \item{\code{init_iter}}{Derived value \code{n_ages * 10} used to initialise
#'       equilibrium conditions.}
#'     \item{\code{run_feedback}, \code{feedback_start_yr}}{
#'       Feedback control settings.}
#'   }
#'
#'
#' @export Setup_Sim_Dim
#' @family Simulation Setup
Setup_Sim_Dim <- function(n_sims,
                          n_yrs,
                          n_pop = 1,
                          n_seas = 1,
                          n_regions,
                          natal_region = NULL,
                          n_ages,
                          n_lens = NULL,
                          n_obs_ages = n_ages,
                          n_sexes,
                          n_fish_fleets,
                          n_srv_fleets,
                          seasdur = if(n_seas == 1) 1 else rep(1 / n_seas, n_seas),
                          run_feedback = FALSE,
                          feedback_start_yr = NULL
                          ) {

  sim_list <- list() # setup empty list1

  if(n_sexes > 2) stop("The number of sexes modelled cannot be larger than 2!")

  if(is.null(natal_region)) {
    if(n_regions == 1) {
      natal_region <- rep(1, n_pop) # all home to region 1
    } else if(n_pop == n_regions) {
      natal_region <- seq_len(n_pop)       # 1:1 mapping
    } else if(n_pop == 1) {
      natal_region <- rep(1, n_pop)
    } else {
      stop("natal_region must be specified when n_pop != n_regions and n_regions > 1")
    }
  }

  # output dimensions into list
  sim_list$n_sims <- n_sims
  sim_list$n_yrs <- n_yrs
  sim_list$n_pop <- n_pop
  sim_list$natal_region <- natal_region
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

#' Initialise model dimension settings
#'
#' Creates the foundational \code{input_list} object used throughout the estimation
#' model setup. All downstream configuration functions
#' (\code{\link{Setup_Mod_Biologicals}}, etc.) expect an \code{input_list} produced by
#' this function. Dimension vectors and scalars are stored in \code{$data}, with
#' empty \code{$par} and \code{$map} sublists ready to be populated by subsequent
#' setup calls. Optionally, a \code{$config} sublist can be stored to retain
#' model configuration options.
#'
#' @param years Numeric vector of calendar years included in the model (e.g.,
#'   \code{1990:2024}). The length of this vector determines \code{n_years}
#'   throughout the model.
#' @param ages Numeric vector of modelled age classes (e.g., \code{2:31} for a
#'   model spanning ages 2–31). The final element is treated as a plus-group.
#' @param lens Numeric vector of length bin midpoints. Set to \code{NULL} when
#'   length data are not modelled; a scalar placeholder of \code{1} is stored
#'   internally in that case.
#' @param n_pop Positive integer. Number of distinct populations (default \code{1}).
#'   Populations can have independent stock-recruit relationships and natal regions
#'   but share the spatial domain defined by \code{n_regions}.
#' @param n_regions Positive integer. Number of spatial regions.
#' @param natal_region Integer vector of length \code{n_pop} mapping each population
#'   to its natal region (1-indexed). Defaults are applied when \code{NULL}:
#'   \itemize{
#'     \item \code{n_regions == 1}: all populations assigned to region 1.
#'     \item \code{n_pop == n_regions}: one-to-one mapping (\code{1:n_pop}).
#'     \item \code{n_pop == 1}: assigned to region 1.
#'   }
#'   Must be supplied explicitly when \code{n_pop > 1}, \code{n_pop != n_regions},
#'   and \code{n_regions > 1}, as no sensible default exists when populations
#'   must share a natal region.
#' @param n_seas Positive integer. Number of seasons within each year (default
#'   \code{1}).
#' @param seasdur Numeric vector of length \code{n_seas} giving the duration of
#'   each season as a fraction of a year (values should sum to 1). Defaults to
#'   \code{1} for a single season or \code{rep(1 / n_seas, n_seas)} for
#'   equal-length seasons when \code{n_seas > 1}.
#' @param n_sexes Integer. Number of sexes; must be \code{1} (sex-aggregated) or
#'   \code{2} (sex-structured).
#' @param n_fish_fleets Positive integer. Number of fishery fleets.
#' @param n_srv_fleets Positive integer. Number of survey fleets.
#' @param n_proj_yrs_devs Non-negative integer. Number of projection years for
#'   which deviation parameters (\code{ln_RecDevs}, \code{move_devs},
#'   \code{ln_fishsel_devs}, \code{ln_srvsel_devs}) are allocated. Set to
#'   \code{0} (default) when projections beyond the assessment period are not
#'   required.
#' @param verbose Logical. If \code{TRUE}, prints a summary of all dimension
#'   settings to the console via \code{message()} after setup. Default
#'   \code{FALSE}.
#' @param store_config Logical. If \code{TRUE}, stores configuration of the model. Default is \code{FALSE}
#' @param do_internal_comp_osa Logical. If \code{TRUE}, allows OSA residuals for composition datasets.
#' Default \code{FALSE}.
#' @param do_internal_conv_tag_osa Logical. If \code{TRUE}, allows OSA residuals for tagging datasets.
#' Default \code{FALSE}.
#'
#' @return A named list (\code{input_list}) with three sublists:
#'   \describe{
#'     \item{\code{$data}}{Dimension scalars and vectors stored for use by the
#'       TMB/RTMB objective function: \code{years}, \code{ages}, \code{lens},
#'       \code{n_pop}, \code{n_regions}, \code{natal_region}, \code{n_seas},
#'       \code{seasdur}, \code{n_sexes}, \code{n_fish_fleets},
#'       \code{n_srv_fleets}, \code{n_proj_yrs_devs}.}
#'     \item{\code{$par}}{Empty list to be populated with parameter starting
#'       values by downstream setup functions.}
#'     \item{\code{$map}}{Empty list to be populated with TMB/RTMB factor maps
#'       by downstream setup functions.}
#'   }
#'   The top-level \code{$verbose} flag is also set and respected by all
#'   subsequent setup functions.
#'
#'
#' @export Setup_Mod_Dim
#' @family Model Setup
Setup_Mod_Dim <- function(years,
                          ages,
                          lens,
                          n_pop = 1,
                          natal_region = NULL,
                          n_seas = 1,
                          seasdur = if(n_seas == 1) 1 else rep(1 / n_seas, n_seas),
                          n_regions,
                          n_sexes,
                          n_fish_fleets,
                          n_srv_fleets,
                          n_proj_yrs_devs = 0,
                          do_internal_comp_osa = FALSE,
                          do_internal_conv_tag_osa = FALSE,
                          verbose = FALSE,
                          store_config = FALSE
                          ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Create empty list
  input_list <- list(data = list(), par = list(), map = list())
  if(store_config) {
    input_list$config <- list()
    input_list$config$Setup_Mod_Dim <- mget(names(formals()))
  }

  if(is.null(natal_region)) {
    if(n_regions == 1) {
      natal_region <- rep(1, n_pop) # all home to region 1
    } else if(n_pop == n_regions) {
      natal_region <- seq_len(n_pop)       # 1:1 mapping
    } else if(n_pop == 1) {
      natal_region <- rep(1, n_pop)
    } else {
      stop("natal_region must be specified when n_pop != n_regions and n_regions > 1")
    }
  }

  # ouput variables into list
  input_list$data$years <- years
  input_list$data$n_pop <- n_pop
  input_list$data$n_regions <- n_regions
  input_list$data$natal_region <- natal_region
  input_list$data$ages <- ages
  input_list$data$lens <- if(is.null(lens)) 1 else lens
  input_list$data$n_sexes <- n_sexes
  input_list$data$n_fish_fleets <- n_fish_fleets
  input_list$data$n_srv_fleets <- n_srv_fleets
  input_list$data$n_proj_yrs_devs <- n_proj_yrs_devs
  input_list$data$n_seas <- n_seas
  input_list$data$seasdur <- seasdur
  input_list$verbose <- verbose
  input_list$store_config <- store_config
  input_list$data$do_internal_comp_osa <- do_internal_comp_osa
  input_list$data$do_internal_conv_tag_osa <- do_internal_conv_tag_osa

  collect_message("Number of Years: ", length(years))
  collect_message("Number of Seasons: ", n_seas)
  for(i in 1:n_seas) collect_message("Duration of season ", i, ": ", seasdur[i])
  collect_message("Number of Projection Years for Dev Pars: ", n_proj_yrs_devs)
  collect_message("Number of Regions: ", n_regions)
  collect_message("Number of Populations: ", n_pop)
  collect_message("Number of Age Bins: ", length(ages))
  collect_message("Number of Length Bins: ", length(lens))
  collect_message("Number of Sexes: ", n_sexes)
  collect_message("Number of Fishery Fleets: ", n_fish_fleets)
  collect_message("Number of Survey Fleets: ", n_srv_fleets)

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}
