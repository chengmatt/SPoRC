#' Set up recruitment dynamics for simulation
#'
#' @param sim_list Simulation list object from `Setup_Sim_Dim()`
#' @param do_recruits_move Indicator for whether recruits move (default = 0):
#'   \itemize{
#'     \item \code{0}: No movement
#'     \item \code{1}: Move
#'   }
#' @param recruitment_opt Recruitment type (default = "bh_rec"):
#'   \itemize{
#'     \item \code{0} or \code{"mean_rec"}: Mean recruitment
#'     \item \code{1} or \code{"bh_rec"}: Beverton-Holt recruitment
#'     \item \code{999} or \code{"resample_from_input"}: Resampling recruitment years from `Rec_input` and preserves covariance of recruitment among regions if spatially-explicit values are provided
#'   }
#' @param rec_dd Recruitment density dependence (default = "global"):
#'   \itemize{
#'     \item \code{0} or \code{"local"}: Region-specific
#'     \item \code{1} or \code{"global"}: Shared across regions
#'   }
#' @param init_dd Initial age density dependence (default = "global"):
#'   \itemize{
#'     \item \code{0} or \code{"local"}: Region-specific
#'     \item \code{1} or \code{"global"}: Shared across regions
#'   }
#' @param rec_lag Recruitment lag (default = 1)
#' @param sexratio_input Sex ratio array [n_regions × n_yrs × n_sexes × n_sims]
#'   (default = 1 if one sex, else 0.5 for each sex)
#' @param R0_input Unfished recruitment (R0) array [n_regions × n_yrs × n_sims]
#'   (default = 10)
#' @param h_input Steepness array [n_regions × n_yrs × n_sims]
#'   (default = 0.8)
#' @param ln_sigmaR Logarithmic standard deviation of recruitment [2]:
#'   1st = sigma for initial devs, 2nd = sigma for latter devs
#'   (default = log(c(1, 1)))
#' @param Rec_input Recruitment array [n_regions × n_yrs × n_sims] (default = NULL)
#' @param ln_InitDevs_input Initial deviations [n_regions × (n_ages-1) × n_sims] (default = NULL)
#' @param init_age_strc Integer specifying the initialization method for the age structure (default = 2):
#'   \itemize{
#'     \item \code{0} or \code{"iterative"}: Iterative solution to equilibrium
#'     \item \code{1} or \code{"scalar_no_move"}: Scalar geometric series solution w/o movement in any groups
#'     \item \code{2} or \code{"matrix"}: Matrix geometric series solution (generalizes scalar solution with movement)
#'     \item \code{3} or \code{"scalar_plus_only"}: Scalar geometric series solution w/o movement only in plus group
#'   }
#' @param t_spawn Spawn timing fraction within the year / season (scalar, default = 0, where spawning happens before mortality processes)
#' @param spawn_seas Season in which spawning occurs
#'
#' @export Setup_Sim_Rec
#' @family Simulation Setup
Setup_Sim_Rec <- function(
    do_recruits_move = 0,
    sexratio_input = array(if(sim_list$n_sexes == 1) 1 else 0.5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_sims)),
    R0_input = array(10, dim = c(sim_list$n_pop,sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    h_input = array(0.8, dim = c(sim_list$n_pop,sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    ln_sigmaR = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions)),
    recruitment_opt = 'bh_rec',
    rec_dd = 'global',
    init_dd = 'global',
    sim_list,
    init_age_strc = 2,
    spawn_seas = 1,
    t_spawn = 0,
    rec_lag = 1,
    Rec_input = NULL,
    ln_InitDevs_input = NULL
    ) {

  if(rec_dd == 'global' && sim_list$n_pop > 1 && recruitment_opt == 'bh_rec') stop("Invalid recruitment density-dependence option! When n_pop > 1 and recruitment_opt == 'bh_rec', rec_dd must be local (0).")

  # Convert character inputs to numeric codes
  recruitment_opt <- convert_to_numeric(recruitment_opt, list(mean_rec = 0, bh_rec = 1, resample_from_input = 999))
  rec_dd <- convert_to_numeric(rec_dd, list(local = 0, global = 1))
  init_dd <- convert_to_numeric(init_dd, list(local = 0, global = 1))
  init_age_strc <- convert_to_numeric(init_age_strc, list(iterative = 0, scalar_no_move = 1, matrix = 2, scalar_plus_only = 3))

  check_sim_dimensions(sexratio_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_sexes = sim_list$n_sexes, n_sims = sim_list$n_sims, n_pop = sim_list$n_pop, what = "sexratio_input")
  check_sim_dimensions(R0_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_sims = sim_list$n_sims, n_pop = sim_list$n_pop, what = "R0_input")
  check_sim_dimensions(h_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_sims  = sim_list$n_sims, n_pop = sim_list$n_pop, what = "h_input")
  if(!is.null(ln_InitDevs_input)) check_sim_dimensions(ln_InitDevs_input, n_regions = sim_list$n_regions, n_ages = sim_list$n_ages, n_sims = sim_list$n_sims, n_pop = sim_list$n_pop, what = "ln_InitDevs_input")

  # Recruitment options
  sim_list$do_recruits_move <- do_recruits_move
  if(sim_list$do_recruits_move == 0) sim_list$move_age <- 2 else sim_list$move_age <- 1 # what age to start movement of individuals

  if(recruitment_opt == 999) {
    if(is.null(Rec_input)) stop("Recruitment input is NULL, but future recruitment is specified to be resampled!")
    rec_input_yrs <- dim(Rec_input)[2] # get years from Rec_input
    tmp_Rec_input <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
    # loop through simulations to resample years
    for(i in 1:sim_list$n_sims) {
      tmp_Rec_input[,1:rec_input_yrs,i] <- Rec_input[,,i]
      resampled_years <- sample(1:rec_input_yrs, length(tmp_Rec_input[1,-c(1:rec_input_yrs),i]), TRUE)
      tmp_Rec_input[,-c(1:rec_input_yrs),i] <- Rec_input[,resampled_years,i]
    } # end i loop
    Rec_input <- tmp_Rec_input # overwrite
  } # resampling

  # Output recruitment stuff into environment
  sim_list$recruitment_opt <- recruitment_opt
  sim_list$rec_dd <- rec_dd
  sim_list$init_dd <- init_dd
  sim_list$h <- h_input
  sim_list$R0 <- R0_input
  sim_list$sexratio <- sexratio_input
  sim_list$rec_lag <- rec_lag
  sim_list$ln_sigmaR <- ln_sigmaR
  sim_list$t_spawn <- t_spawn
  sim_list$init_age_strc <- init_age_strc
  sim_list$spawn_seas <- spawn_seas
  if(!is.null(Rec_input)) sim_list$Rec_input <- Rec_input
  if(!is.null(ln_InitDevs_input)) sim_list$ln_InitDevs_input <- ln_InitDevs_input

  return(sim_list)

}

#' Title Helper function to handle sigmaR mapping
#'
#' @param input_list Input list
#' @param sigmaR_spec Character vector for specifying sigmaR mapping
#'
#' @returns Input list with mapping modified
#' @keywords internal
do_sigmaR_mapping <- function(input_list, sigmaR_spec) {

  # Define valid sigmaR options
  valid_options <- c("est_all", "est_shared", "fix_early_est_late", "fix")

  # Checking to see if valid options
  if (!is.null(sigmaR_spec)) {
    if (!sigmaR_spec %in% valid_options) {
      stop("Invalid sigmaR_spec. Must be one of: ", paste(valid_options, collapse = ", "))
    }

    # Switch for defining sigmaR mapping
    input_list$map$ln_sigmaR <- switch(
      sigmaR_spec,
      est_all = factor(c(1,2)),
      est_shared = factor(c(1, 1)),
      fix_early_est_late = factor(c(NA, 1)),
      fix = factor(c(NA, NA))
    )

    collect_message("Recruitment Variability is specified as: ", sigmaR_spec)
  } else {
    collect_message("Recruitment Variability is estimated for both early and late periods")
  }
  return(input_list)
}

#' Helper function to handle mapping for ln_InitDevs
#'
#' @param input_list Input list
#' @param InitDevs_spec Character vector for specifying InitDevs mapping
#' @param rec_dd Recruitment density dependence indicator (global vs. local)
#' @keywords internal
do_InitDevs_mapping <- function(input_list, InitDevs_spec, rec_dd) {

  # Initial age deviations (equilibrium)
  if(input_list$data$equil_init_age_strc == 0) {
    input_list$par$ln_InitDevs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$ages) - 1)) # override starting values if previously specified
    input_list$map$ln_InitDevs <- factor(rep(NA, length(input_list$par$ln_InitDevs))) # set mapping
    collect_message("Initial Age Structure is specified to be in equilibrium. No initial age deviations are estimated.")
  }

  # Initial age deviations (stochastic for all ages, including plus group)
  if(!is.null(InitDevs_spec)) {

    # Validate options
    if(rec_dd == 'global' && InitDevs_spec != "est_shared_r" && input_list$data$n_regions > 1) {
      stop("Please specify a valid initial age deviations option for global recruitment density dependence (should be est_shared_r or leave as NULL)!")
    }

    if(!InitDevs_spec %in% c("est_shared_r", "fix")) stop("Please specify a valid initial deviations option. These include: fix, est_shared_r. Conversely, leave at NULL to estimate all initial deviations.")
    else collect_message("Initial Deviations is stochastic and specified as: ", InitDevs_spec)

    # set up mapping for initial age deviations
    map_InitDevs <- input_list$par$ln_InitDevs

    # Fix all initial deviations
    if(InitDevs_spec == "fix") input_list$map$ln_InitDevs <- factor(rep(NA, prod(dim(map_InitDevs))))

    # Share across regions and estimate
    if(InitDevs_spec == "est_shared_r") {

      # share parameters, but no stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 1) {
        for(r in 1:input_list$data$n_regions) {
          map_InitDevs[r,-dim(input_list$par$ln_InitDevs)[2]] <- 1:length(map_InitDevs[1,-dim(input_list$par$ln_InitDevs)[2]]) # share parameters across regions (but don't estimate plus group)
          map_InitDevs[r,dim(input_list$par$ln_InitDevs)[2]] <- NA # NA for plus group
          input_list$par$ln_InitDevs[r,dim(input_list$par$ln_InitDevs)[2]] <- 0 # reset plus group starting value to 0
        } # end r loop
        collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
      }

      # share parameters across regions, with stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 2) {
        for(r in 1:input_list$data$n_regions) {
          map_InitDevs[r,] <- 1:length(map_InitDevs[1,])
        }
        collect_message("Initial age deviations are stochastic and estimated for all ages, including the plus group")
      }
      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
    } # end if

  } else { # If NULL, then estimating age deviations across all dimensions
    map_InitDevs <- input_list$par$ln_InitDevs # set up mapping for initial age deviations

    if(input_list$data$equil_init_age_strc == 1) { # estimating all deviations across all dimensions, except for plus group
      map_InitDevs[,-dim(input_list$par$ln_InitDevs)[2]] <- 1:length(map_InitDevs[,-dim(input_list$par$ln_InitDevs)[2]]) # don't estimate plus group
      map_InitDevs[,dim(input_list$par$ln_InitDevs)[2]] <- NA # NA for plus group
      input_list$par$ln_InitDevs[,dim(input_list$par$ln_InitDevs)[2]] <- 0 # reset plus group starting value to 0
      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
      collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
    }

    # Plus group and estimating deviations for all dimensions
    if(input_list$data$equil_init_age_strc == 2) {
      input_list$map_ln_InitDevs <- factor(1:length(map_InitDevs)) # input into map
      collect_message("Initial Age Deviations is estimated for all dimensions. They are are stochastic and estimated for all ages, including the plus group")
    }
  }
  return(input_list)
}

#' Helper function to handle RecDevs mapping
#'
#' @param input_list Input list
#' @param RecDevs_spec Character vector for specifying RecDevs mapping
#' @param rec_dd Recruitment density dependence indicator (global vs. local)
#' @keywords internal
do_RecDevs_mapping <- function(input_list, RecDevs_spec, rec_dd) {

  map_RecDevs <- input_list$par$ln_RecDevs # set up mapping for recruitment deviations

  # Recruitment deviations
  if(!is.null(RecDevs_spec)) {

    # Validate options
    if(rec_dd == 'global' && RecDevs_spec != "est_shared_r" && input_list$data$n_regions > 1) stop("Please specify a valid recruitment deviations option for global recruitment density dependence (should be est_shared_r)!")
    if(!RecDevs_spec %in% c("est_shared_r", "fix"))  stop("Please specify a valid recruitment deviations option. These include: fix, est_shared_r. Conversely, leave at NULL to estimate all recruitment deviations.")

    # Share across regions and estimate
    if(RecDevs_spec == "est_shared_r") {
      for(r in 1:input_list$data$n_regions) map_RecDevs[r,] <- 1:length(map_RecDevs[1,]) # share parameters across regions
      input_list$map$ln_RecDevs <- factor(map_RecDevs)
    } # end if

    # Fix all recruitment deviations
    if(RecDevs_spec == "fix") input_list$map$ln_RecDevs <- factor(rep(NA, prod(dim(map_RecDevs))))

    # print message
    collect_message("Recruitment Deviations is specified as: ", RecDevs_spec)

  } else { # if NULL, estimating all dimensions
    input_list$map$ln_RecDevs <- factor(1:length(map_RecDevs)) # input into mapping
    collect_message("Recruitment Deviations is estimated for all dimensions")
  }
  return(input_list)
}

#' Helper function to setup steepness mapping
#'
#' @param input_list Input list
#' @param h_spec Character vector for specifying steepness mapping
#' @param rec_dd Recruitment density dependence indicator (global vs. local)
#' @keywords internal
do_h_mapping <- function(input_list, h_spec, rec_dd) {

  # Mean recruitment
  if(input_list$data$rec_model == 0) {
    input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))
  } else if(!is.null(h_spec)) {

    # Validate options
    if(rec_dd == 'global' && !h_spec %in% c("est_shared_r", "fix") && input_list$data$n_regions > 1) stop("Please specify a valid steepness option for global recruitment density dependence (should be est_shared_r or fix)!")
    if(!h_spec %in% c("est_shared_r", "fix"))  stop("Please specify a valid steepness option. These include: fix, est_shared_r. Conversely, leave at NULL to estimate all steepness values.")

    # Share across regions and estimate
    if(h_spec == "est_shared_r") input_list$map$steepness_h <- factor(rep(1, length(input_list$par$steepness_h)))

    # Fix all steepness values
    if(h_spec == "fix") input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))

    collect_message("Steepness is specified as: ", h_spec) # output message
  } else {
    # if beverton holt and estimating all steepness parameters
    if(input_list$data$rec_model == 1) {
      # Validate options
      if(rec_dd == 'global' && input_list$data$n_regions > 1) stop("Please specify a valid steepness option for global recruitment density dependence (should be est_shared_r)!")
      input_list$map$steepness_h <- factor(c(1:input_list$data$n_regions)) # estimating all steepness parameters
    }
    collect_message("Steepness is estimated for all dimensions")
  }
  return(input_list)
}

#' Helper function to set up sex ratio parameters
#'
#' @param input_list Input list
#' @param sexratio_spec Charcacter specifying sex ratio parameterization
#' @keywords internal
do_sexratio_pars_mapping <- function(input_list, sexratio_spec) {

  # Initialize arrays and counters
  map_sexratio <- input_list$par$sexratio_pars
  map_sexratio[] <- NA
  sexratio_counter <- 1

  # Validate inputs here
  if(!sexratio_spec %in% c("est_all", "est_shared_r", "fix")) stop("Sex Ratio Specificaiton is not correctly specified. Needs to be fix, est_all, or est_shared_r")

  # if we want to fix
  if(sexratio_spec == 'fix') map_sexratio[] <- NA

  for(r in 1:input_list$data$n_regions) {

    # Get number of sex ratio rate blocks
    sexratio_blocks_tmp <- unique(as.vector(input_list$data$sexratio_blocks[r,]))

    for(b in 1:length(sexratio_blocks_tmp)) {

      # Estimate for all regions
      if(sexratio_spec == 'est_all') {
        map_sexratio[r,b] <- sexratio_counter
        sexratio_counter <- sexratio_counter + 1
      }

      # Estimate but share sex ratio across regions
      if(sexratio_spec == 'est_shared_r' && r == 1) {
        for(rr in 1:input_list$data$n_regions) {
          # only assign if this value exists for this region
          if(sexratio_blocks_tmp[b] %in% input_list$data$sexratio_blocks[rr,]) {
            map_sexratio[rr, b] <- sexratio_counter
          } # end if
        } # end rr loop
        sexratio_counter <- sexratio_counter + 1
      }

    } # end b loop
  } # end r loop

  collect_message("Sex ratio is specified as: ", sexratio_spec)

  # input sex ratio rates into mapping list
  input_list$map$sexratio_pars <- factor(map_sexratio) # sex ratio rates

  return(input_list)
}

#' Helper function to map recruitment proportions
#'
#' @param input_list Input list
#' @keywords internal
do_Rec_prop_mapping <- function(input_list) {
  # map off parameters if single region
  if(input_list$data$n_regions == 1) input_list$map$Rec_prop <- factor(rep(NA, length(input_list$par$Rec_prop)))
  return(input_list)
}

#' Setup model objects for specifying recruitment module and associated processes
#'
#' @param input_list List containing data, parameters, and map lists used by the model.
#' @param rec_model Character string specifying the recruitment model. Options are:
#' \itemize{
#'   \item \code{"mean_rec"}: Recruitment is a fixed mean value.
#'   \item \code{"bh_rec"}: Beverton-Holt recruitment with steepness parameter.
#' }
#' @param rec_lag Integer specifying the recruitment lag duration relative to spawning stock biomass (SSB).
#' @param Use_h_prior Integer flag (0 or 1) indicating whether to apply a prior on steepness \code{h}.
#' @param do_rec_bias_ramp Integer flag (0 or 1) indicating whether to apply a recruitment bias correction ramp.
#' @param bias_year Numeric vector of length 4 defining the recruitment bias ramp periods:
#'   \itemize{
#'     \item Element 1: End year of no bias correction period.
#'     \item Element 2: End year of ascending bias ramp period.
#'     \item Element 3: End year of full bias correction period.
#'     \item Element 4: Start year of final no bias correction period.
#'   }
#'   For example, with 65 years total, \code{c(21, 31, 60, 64)} means:
#'   \itemize{
#'     \item Years 1–21: No bias correction.
#'     \item Years 22–31: Ascending bias correction.
#'     \item Years 32–60: Full bias correction.
#'     \item Years 61–63: Descending bias ramp.
#'     \item Years 64–65: No bias correction.
#'   }
#' @param sigmaR_switch Integer year indicating when \code{sigmaR} switches from early to late values (0 disables switching).
#' @param init_age_strc Initialization method for initial age structure (default = 2):
#'   \itemize{
#'     \item \code{0} or \code{"iterative"}: Initialize by iteration.
#'     \item \code{1} or \code{"scalar_no_move"}: Scalar geometric series w/o any movement in all groups.
#'     \item \code{2} or \code{"matrix"}: Matrix geometric series (accounts for movement).
#'     \item \code{3} or \code{"scalar_plus_only"}: Scalar geometric series w/o movement only in plus groups.
#'   }
#' @param init_F_prop Numeric vector (n_seas) specifying the initial fishing mortality proportion relative to mean fishing mortality for initializing age structure.
#' @param sigmaR_spec Character string specifying estimation of recruitment variability (\code{sigmaR}):
#' \itemize{
#'   \item \code{NULL or "est_all"}: Estimate separate \code{sigmaR} for early and late periods.
#'   \item \code{"est_shared"}: Estimate one \code{sigmaR} shared across periods.
#'   \item \code{"fix"}: Fix both \code{sigmaR} values.
#'   \item \code{"fix_early_est_late"}: Fix early \code{sigmaR}, estimate late \code{sigmaR}.
#' }
#' @param InitDevs_spec Character string specifying estimation of initial age deviations:
#' \itemize{
#'   \item \code{NULL}: Estimate deviations for all ages and regions.
#'   \item \code{"est_shared_r"}: Estimate deviations shared across regions.
#'   \item \code{"fix"}: Fix all deviations.
#' }
#' @param RecDevs_spec Character string specifying recruitment deviation estimation:
#' \itemize{
#'   \item \code{NULL}: Estimate deviations for all regions and years.
#'   \item \code{"est_shared_r"}: Estimate deviations shared across regions (global recruitment deviations).
#'   \item \code{"fix"}: Fix all recruitment deviations.
#' }
#' @param dont_est_recdev_last Integer specifying how many of the most recent recruitment deviations to not estimate. Default is 0.
#' @param h_spec Character string specifying steepness estimation:
#' \itemize{
#'   \item \code{NULL}: Estimate steepness for all regions if \code{rec_model == "bh_rec"}.
#'   \item \code{"est_shared_r"}: Estimate steepness shared across regions.
#'   \item \code{"fix"}: Fix steepness values.
#' }
#'   If \code{rec_model == "mean_rec"}, steepness is fixed.
#' @param rec_dd Character string specifying recruitment density dependence, options:
#' \code{"local"}, \code{"global"}, or \code{NULL}.
#' @param t_spawn Spawn timing fraction within the year / season (scalar, default = 0, where spawning happens before mortality processes)
#' @param ... Additional arguments specifying starting values for recruitment parameters such as \code{ln_global_R0}, \code{Rec_prop}, \code{h}, \code{ln_InitDevs}, \code{ln_RecDevs}, and \code{ln_sigmaR}.
#' @param equil_init_age_strc How initial age structure deviations should be initialized (default = 1):
#'   \itemize{
#'     \item \code{0} or \code{"equil"}: Equilibrium initial age structure.
#'     \item \code{1} or \code{"stoch_no_plus"}: Stochastic for all ages except plus group (follows equilibrium).
#'     \item \code{2} or \code{"stoch_all"}: Stochastic initial age structure for all ages.
#'   }
#' @param max_bias_ramp_fct Numeric specifying the maximum bias correction to apply to the recruitment bias ramp (should be between 0 and 1)
#' @param h_prior Data frame specifying beta prior distributions for the `h_trans` parameters.
#'   Must include the following columns:
#'   - `region`: Integer region index corresponding to the element in `h_trans` being penalized.
#'   - `mu`: Mean of the prior in normal space (used to calculate the corresponding beta distribution).
#'   - `sd`: Standard deviation of the prior in normal space.
#'   For each row, a beta distribution is scaled to the interval [0.2, 1], and the corresponding element of `h_trans` is transformed to that scale and penalized using the log-density from the beta distribution.
#' @param Use_Rec_prop_Prior Integer flag (0 or 1) indicating whether to apply a prior on recruitment proportions.
#' @param Rec_prop_prior Scalar or array specifying prior values for recruitment proportion parameters.
#'   If scalar, a constant uniform prior is applied across all dimensions.
#' @param sexratio_blocks Character vector specifying blocks of years and regions for sex ratio. Format examples:
#'   \itemize{
#'     \item \code{"Block_1_Year_1-15_Region_1"}
#'     \item \code{"Block_2_Year_16-terminal_Region_2"}
#'     \item \code{"none_Region_3"} (means no block, constant for that region; this is the default option)
#'   }
#' @param sexratio_spec Character string specifying sex ratio estimation scheme:
#'   \itemize{
#'     \item \code{"est_all"} estimates sex ratio for all blocks and regions independently
#'     \item \code{"est_shared_r"} estimates sex ratio shared across regions but varying by block
#'     \item \code{"fix"} fixes all sex ratio (no estimation)
#'     }
#' @param spawn_seas Season in which spawning occurs. Default is season 1.
#' @export Setup_Mod_Rec
#' @family Model Setup
Setup_Mod_Rec <- function(input_list,
                          rec_model,
                          rec_dd = NULL,
                          rec_lag = 1,
                          Use_h_prior = 0,
                          h_prior = NULL,
                          Use_Rec_prop_Prior = 0,
                          Rec_prop_prior = NULL,
                          do_rec_bias_ramp = 0,
                          bias_year = NA,
                          max_bias_ramp_fct = 1,
                          sigmaR_switch = 1,
                          dont_est_recdev_last = 0,
                          init_age_strc = 2,
                          equil_init_age_strc = 1,
                          init_F_prop = rep(0, input_list$data$n_seas),
                          sigmaR_spec = NULL,
                          InitDevs_spec = NULL,
                          RecDevs_spec = NULL,
                          h_spec = NULL,
                          t_spawn = 0,
                          spawn_seas = 1,
                          sexratio_spec = 'fix',
                          sexratio_blocks = c(
                            paste("none_Region_", c(1:input_list$data$n_regions), sep = '')
                          ),
                          ...
                          ) {

  messages_list <<- character(0)
  starting_values <- list(...)

  # Convert character inputs to numeric codes for init_age_strc and equil_init_age_strc
  init_age_strc <- convert_to_numeric(init_age_strc, list(iterative = 0, scalar_no_move = 1, matrix = 2, scalar_plus_only = 3))
  equil_init_age_strc <- convert_to_numeric(equil_init_age_strc, list(equil = 0, stoch_no_plus = 1, stoch_all = 2))

  # Recruitment Model Type and Options --------------------------------------

  # Recruitment model
  rec_model_map <- list(mean_rec = 0, bh_rec = 1)
  if (!rec_model %in% names(rec_model_map)) stop("Invalid recruitment model. Use 'mean_rec' or 'bh_rec'")
  rec_model_val <- rec_model_map[[rec_model]]
  collect_message("Recruitment is specified as: ", rec_model)

  # Recruitment density dependence
  if (!is.null(rec_dd)) {
    rec_dd_map <- list(local = 0, global = 1)
    if (!rec_dd %in% names(rec_dd_map)) stop("Invalid rec_dd. Use 'local' or 'global'")
    rec_dd_val <- rec_dd_map[[rec_dd]]
    collect_message("Recruitment Density Dependence is specified as: ", rec_dd)
  } else {
    rec_dd_val <- ifelse(rec_model == "mean_rec", 999, 1)
  }

  # Recruitment lag
  if (rec_model != "mean_rec") collect_message("Recruitment and SSB lag is specified as: ", rec_lag)

  # Recruitment proportion prior
  if(!Use_Rec_prop_Prior %in% c(0,1)) stop("Use_Rec_prop_Prior must be 0 or 1")
  if(Use_Rec_prop_Prior == 1 && input_list$data$n_regions == 1) stop("Priors should not be applied to recruitment proportions when n_regions = 1.")
  collect_message("Recruitment proportion priors are: ", ifelse(Use_Rec_prop_Prior == 1, "Used", "Not Used"))

  # Steepness Settings ------------------------------------------------------

  if (rec_model == "bh_rec") {
    if (!Use_h_prior %in% c(0, 1)) stop("Use_h_prior must be 0 or 1")
    if (Use_h_prior == 1) {
      required_cols <- c("region", "mu", "sd")
      missing_cols <- setdiff(required_cols, names(h_prior))
      if (length(missing_cols) > 0) stop("h_prior is missing columns: ", paste(missing_cols, collapse = ", "))
    }
    collect_message("Steepness priors are: ", ifelse(Use_h_prior == 1, "Used", "Not Used"))
  }


  # Sex Ratio Options ---------------------------------------------
  sexratio_blocks_mat <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years)))

  if(!is.null(sexratio_blocks)) {
    for(i in 1:length(sexratio_blocks)) {

      # Extract out components from list
      tmp <- sexratio_blocks[i]
      tmp_vec <- unlist(strsplit(tmp, "_"))

      if(!tmp_vec[1] %in% c("none", "Block")) stop("Sex Ratio Blocks not correctly specified. This should be either none_Region_x or Block_x_Year_x-y_Region_x")

      # extract out fleets if constant
      if(tmp_vec[1] == "none") {
        region <- as.numeric(tmp_vec[3]) # get region index
        sexratio_blocks_mat[region,] <- 1 # input sex ratio time block
      }

      if(tmp_vec[1] == "Block") {
        block_val <- as.numeric(tmp_vec[2]) # get block value
        region <- as.numeric(tmp_vec[6]) # get region value

        # get year ranges
        if(!str_detect(tmp, "terminal")) { # if not terminal year
          year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
          years <- year_range[1]:year_range[2] # get sequence of years
        } else { # if terminal year
          year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
          years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
        }

        sexratio_blocks_mat[region,years] <- block_val # input sex ratio time block
      }

    } # end i loop
  }

  for(r in 1:input_list$data$n_regions) collect_message("Sex Ratios estimated with ", length(unique(sexratio_blocks_mat[r,])), " block for region ", r)

  # Input Validation --------------------------------------------------------

  # Helper function
  check_in <- function(x, valid, name) {
    if (!x %in% valid) stop(name, " must be one of: ", paste(valid, collapse = ", "))
  }

  # Validation
  check_in(do_rec_bias_ramp, 0:1, "do_rec_bias_ramp")
  check_in(init_age_strc, 0:3, "init_age_strc")
  if(!is.numeric(sigmaR_switch)) stop("sigmaR_switch must be numeric")
  if(max_bias_ramp_fct > 1 || max_bias_ramp_fct < 0) stop("max_bias_ramp_fct must be between 0 and 1")

  # print messages
  collect_message("Recruitment Bias Ramp is: ", ifelse(do_rec_bias_ramp == 0, "Off", 'On'))
  init_age_methods <- c("Iterated", "No Movement and Scalar Geometric Series", "Movement and Matrix Geometric Series", "Movement but Scalar Geometric Series for plus group")
  collect_message("Initial Age Structure is: ", init_age_methods[init_age_strc + 1])
  if(sigmaR_switch > 1) collect_message("Sigma R switches from an early period value to a late period value at year: ", sigmaR_switch)
  collect_message("Recruitment deviations for ", ifelse(dont_est_recdev_last == 0, "every year are estimated", paste("terminal year not estimated -", dont_est_recdev_last)))
  if(dont_est_recdev_last != 0 && input_list$data$n_proj_yrs_devs != 0) {
    collect_message(
      "Recruitment deviations were specified to not be estimated for the last ",
      dont_est_recdev_last,
      " years, but n_proj_yrs_devs != 0. Because projected deviations are still computed (penalized toward the mean), those `unestimated` years are stil effectively estimated. Setting dont_est_recdev_last to 0."
    )
    dont_est_recdev_last <- 0 # overwrite at 0
  }

  # Populate Data List ------------------------------------------------------

  # input variables into data list
  collect_message("Spawning season occurs in season ", spawn_seas)
  input_list$data$spawn_seas <- spawn_seas
  input_list$data$rec_model <- rec_model_val
  input_list$data$rec_dd <- rec_dd_val
  input_list$data$rec_lag <- rec_lag
  input_list$data$Use_h_prior <- Use_h_prior
  input_list$data$h_prior <- h_prior
  input_list$data$do_rec_bias_ramp <- do_rec_bias_ramp
  input_list$data$bias_year <- bias_year
  input_list$data$sigmaR_switch <- sigmaR_switch
  input_list$data$init_age_strc <- init_age_strc
  input_list$data$init_F_prop <- init_F_prop
  input_list$data$t_spawn <- t_spawn
  input_list$data$equil_init_age_strc <- equil_init_age_strc
  input_list$data$max_bias_ramp_fct <- max_bias_ramp_fct
  input_list$data$Use_Rec_prop_Prior <- Use_Rec_prop_Prior
  Rec_prior_vals = ifelse(is.null(Rec_prop_prior), rep(1, input_list$data$n_regions), Rec_prop_prior)
  input_list$data$Rec_prop_prior <- array(Rec_prior_vals, dim = c(input_list$data$n_regions))
  input_list$data$sexratio_blocks <- sexratio_blocks_mat

  # Populate Parameter List -------------------------------------------------

  # Global R0
  if("ln_global_R0" %in% names(starting_values)) input_list$par$ln_global_R0 <- starting_values$ln_global_R0
  else input_list$par$ln_global_R0 <- log(15)

  # R0 proportion
  if("Rec_prop" %in% names(starting_values)) input_list$par$Rec_prop <- starting_values$Rec_prop
  else input_list$par$Rec_prop <- array(rep(0, input_list$data$n_regions - 1))

  # Steepness in bounded logit space (0.2 and 1)
  if("steepness_h" %in% names(starting_values)) input_list$par$steepness_h <- starting_values$steepness_h
  else input_list$par$steepness_h <- rep(0, input_list$data$n_regions)

  # Initial age deviations
  if("ln_InitDevs" %in% names(starting_values)) input_list$par$ln_InitDevs <- starting_values$ln_InitDevs
  else input_list$par$ln_InitDevs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$ages) - 1))

  # Recruitment deviations
  if("ln_RecDevs" %in% names(starting_values)) input_list$par$ln_RecDevs <- starting_values$ln_RecDevs
  else input_list$par$ln_RecDevs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) - dont_est_recdev_last + input_list$data$n_proj_yrs_devs))

  # Recruitment variability
  if("ln_sigmaR" %in% names(starting_values)) input_list$par$ln_sigmaR <- starting_values$ln_sigmaR
  else input_list$par$ln_sigmaR <- c(0,0)  # (early period 1st element, late period 2nd element)

  # sexratio parameters
  max_sexratio_blks <- max(apply(input_list$data$sexratio_blocks, 1, FUN = function(x) length(unique(x)))) # figure out maximum number ofsex ratio blocks for each region
  if("sexratio_pars" %in% names(starting_values)) input_list$par$sexratio_pars <- starting_values$sexratio_pars
  else input_list$par$sexratio_pars <- array(0, dim = c(input_list$data$n_regions, max_sexratio_blks)) # specified at 0.5 in inverse logit space

  # Mapping Options -----------------------------------------------------------

  input_list <- do_sigmaR_mapping(input_list, sigmaR_spec) # sigmaR mapping
  input_list <- do_InitDevs_mapping(input_list, InitDevs_spec, rec_dd) # InitDevs mapping
  input_list <- do_RecDevs_mapping(input_list, RecDevs_spec, rec_dd) # RevDevs mapping
  input_list <- do_h_mapping(input_list, h_spec, rec_dd) # steepness mapping
  input_list <- do_Rec_prop_mapping(input_list) # Recruitment proportion mapping
  input_list <- do_sexratio_pars_mapping(input_list, sexratio_spec) # sex ratio parameters

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}

