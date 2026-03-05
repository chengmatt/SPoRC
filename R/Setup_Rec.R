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
#' @param sexratio_input Sex ratio array [n_pop x n_regions × n_yrs × n_sexes × n_sims]
#'   (default = 1 if one sex, else 0.5 for each sex)
#' @param R0_input Unfished recruitment (R0) array [n_pop x n_regions × n_yrs × n_sims]
#'   (default = 10)
#' @param h_input Steepness array [n_pop x n_regions × n_yrs × n_sims]
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
#' @param stray_rate_input straying rate array [n_pop × n_yrs × n_sims]
#' @param rec_seas_prop_input Seasonal recruitment allocation. Vector dimensioned by [n_pop x n_seas x n_sims].
#'
#' @export Setup_Sim_Rec
#' @family Simulation Setup
Setup_Sim_Rec <- function(
    do_recruits_move = 0,
    sexratio_input = array(if(sim_list$n_sexes == 1) 1 else 0.5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_sims)),
    R0_input = array(10, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    h_input = array(0.8, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    stray_rate_input = array(1, dim = c(sim_list$n_pop, sim_list$n_yrs, sim_list$n_sims)),
    ln_sigmaR = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions)),
    rec_seas_prop_input = {
      rec_seas_prop = array(0, dim = c(sim_list$n_pop, sim_list$n_seas, sim_list$n_sims))
      rec_seas_prop[, 1, ] <- 1
      rec_seas_prop
    },
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
  check_sim_dimensions(stray_rate_input, n_years = sim_list$n_yrs, n_sims  = sim_list$n_sims, n_pop = sim_list$n_pop, what = "stray_rate_input")
  check_sim_dimensions(rec_seas_prop_input, n_seas = sim_list$n_seas, n_sims  = sim_list$n_sims, n_pop = sim_list$n_pop, what = "rec_seas_prop_input")
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
  sim_list$rec_seas_prop <- rec_seas_prop_input
  sim_list$init_age_strc <- init_age_strc
  sim_list$spawn_seas <- spawn_seas
  sim_list$stray_rate <- stray_rate_input
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

  map_sigmaR = input_list$par$ln_sigmaR
  map_sigmaR[] = NA

  # Define valid sigmaR options
  valid_options <- c("est_all", "est_shared_all", "fix_early_est_late", "fix")

  # Checking to see if valid options
  if (!sigmaR_spec %in% valid_options) stop("Invalid sigmaR_spec. Must be one of: ", paste(valid_options, collapse = ", "))

    # SigmaR specifications
    sigmaR_counter <- 1
    for(i in 1:2) {
      for(p in 1:input_list$data$n_pop) {

        if(sigmaR_spec == 'est_all') {
          if(input_list$data$n_pop == 1 && input_list$data$rec_dd == 0) {
            # single pop + local DD: each region gets its own sigmaR
            for(r in 1:input_list$data$n_regions) {
              map_sigmaR[i,p,r] <- sigmaR_counter
              sigmaR_counter <- sigmaR_counter + 1
            }
          } else {
            # multi-pop or global DD: each population gets its own sigmaR, shared across regions
            for(r in 1:input_list$data$n_regions) map_sigmaR[i,p,r] <- sigmaR_counter
            sigmaR_counter <- sigmaR_counter + 1
          }
        }

        if(sigmaR_spec == 'fix_early_est_late' && i == 2) {
          if(input_list$data$n_pop == 1 && input_list$data$rec_dd == 0) {
            # single pop + local DD: each region gets its own sigmaR
            for(r in 1:input_list$data$n_regions) {
              map_sigmaR[i,p,r] <- sigmaR_counter
              sigmaR_counter <- sigmaR_counter + 1
            }
          } else {
            # multi-pop or global DD: each population gets its own sigmaR, shared across regions
            for(r in 1:input_list$data$n_regions) map_sigmaR[i,p,r] <- sigmaR_counter
            sigmaR_counter <- sigmaR_counter + 1
          }
        }
      } # end p loop
    } # end i loop

    if(sigmaR_spec == 'est_shared_all') map_sigmaR[] <- 1
    if(sigmaR_spec == 'fix') map_sigmaR[] <- NA
    input_list$map$ln_sigmaR <- factor(map_sigmaR)
    collect_message("Recruitment Variability is specified as: ", sigmaR_spec)

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
    input_list$par$ln_InitDevs <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$ages) - 1)) # override starting values if previously specified
    input_list$map$ln_InitDevs <- factor(rep(NA, length(input_list$par$ln_InitDevs))) # set mapping
    collect_message("Initial Age Structure is specified to be in equilibrium. No initial age deviations are estimated.")
  }

  # Initial age deviations (stochastic for all ages, including plus group)
  if(!is.null(InitDevs_spec)) {

    # Validate options
    if(!is.null(rec_dd) && rec_dd == 'global' && !InitDevs_spec %in% c("est_shared_r", "est_shared_pop_r") && input_list$data$n_regions > 1) {
      stop("Please specify a valid initial age deviations option for global recruitment density dependence (should be est_shared_r or est_shared_pop_r)!")
    }

    if(!InitDevs_spec %in% c("est_shared_pop_r", "est_shared_r", "fix")) stop("Please specify a valid initial deviations option. These include: fix, est_shared_r, est_shared_pop_r. Conversely, leave at NULL to estimate all initial deviations.")
    else collect_message("Initial Deviations is stochastic and specified as: ", InitDevs_spec)

    # set up mapping for initial age deviations
    map_InitDevs <- input_list$par$ln_InitDevs

    # Fix all initial deviations
    if(InitDevs_spec == "fix") input_list$map$ln_InitDevs <- factor(rep(NA, prod(dim(map_InitDevs))))

    # Share across regions and populations
    if(InitDevs_spec == 'est_shared_pop_r') {

      # share parameters, but no stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 1) {

        # get indices
        n_ages <- dim(input_list$par$ln_InitDevs)[3] - 1
        # each populaiton and reigon has the same indices
        map_InitDevs[,,-dim(input_list$par$ln_InitDevs)[3]] <- rep(1:n_ages, each = input_list$data$n_regions * input_list$data$n_pop)
        map_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- NA  # NA for plus group
        input_list$par$ln_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- 0
        collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
      }

      # share parameters across regions, with stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 2) {

        # get indices
        n_ages_all <- dim(input_list$par$ln_InitDevs)[3]
        map_InitDevs[] <- rep(1:n_ages_all, each = input_list$data$n_regions * input_list$data$n_pop)

        collect_message("Initial age deviations are stochastic and estimated for all ages, including the plus group")
      }

      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
    }

    # Share across regions and estimate for each population
    if(InitDevs_spec == "est_shared_r") {

      # share parameters, but no stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 1) {

        # get indices
        n_ages <- dim(input_list$par$ln_InitDevs)[3] - 1
        n_region <- dim(input_list$par$ln_InitDevs)[2]

        for(p in 1:input_list$data$n_pop) {
          age_indices <- (1:n_ages) + (p-1) * n_ages # get age indices
          # each region gets the same index for a given age (repeat each index across regions)
          map_InitDevs[p,,-dim(input_list$par$ln_InitDevs)[3]] <- matrix(rep(age_indices, each = n_region), nrow = n_region)
          map_InitDevs[p,,dim(input_list$par$ln_InitDevs)[3]] <- NA  # NA for plus group
          input_list$par$ln_InitDevs[p,,dim(input_list$par$ln_InitDevs)[3]] <- 0
        } # end p loop

        collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
      }

      # share parameters across regions, with stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 2) {

        # get indices
        n_ages_all <- dim(input_list$par$ln_InitDevs)[3]
        n_region <- dim(input_list$par$ln_InitDevs)[2]

        for(p in 1:input_list$data$n_pop) {
          age_indices <- (1:n_ages_all) + (p-1) * n_ages_all
          # each region gets the same index for a given age (repeat each index across regions)
          map_InitDevs[p,,] <- matrix(rep(age_indices, each = n_region), nrow = n_region)
        } # end p loop

        collect_message("Initial age deviations are stochastic and estimated for all ages, including the plus group")
      }

      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
    } # end if

  } else { # If NULL, then estimating age deviations across all dimensions

    if(input_list$data$n_pop > 1 && input_list$data$rec_region_prop_spec == 1)
      stop("Can't estimate initial age deviations for all populations and regions if no recruitment dispersal is occuring within a given region! Please specify est_shared_r or est_shared_pop_r instead!")

    map_InitDevs <- input_list$par$ln_InitDevs # set up mapping for initial age deviations

    if(input_list$data$equil_init_age_strc == 1) { # estimating all deviations across all dimensions, except for plus group
      map_InitDevs[,,-dim(input_list$par$ln_InitDevs)[3]] <- 1:length(map_InitDevs[,,-dim(input_list$par$ln_InitDevs)[3]]) # don't estimate plus group
      map_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- NA # NA for plus group
      input_list$par$ln_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- 0 # reset plus group starting value to 0
      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
      collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
    }

    # Plus group and estimating deviations for all dimensions
    if(input_list$data$equil_init_age_strc == 2) {
      input_list$map$ln_InitDevs <- factor(1:length(map_InitDevs)) # input into map
      collect_message("Initial Age Deviations is estimated for all dimensions. They are are stochastic and estimated for all ages, including the plus group")
    }
  }

  # When no_dispersal, non-natal regions have no recruitment so their
  # deviations are structurally zero. Fix them regardless of InitDevs_spec.
  if(input_list$data$rec_region_prop_spec == 1 && input_list$data$n_pop > 1) {

    # extract mapping
    map_tmp <- as.integer(input_list$map$ln_InitDevs)
    dim(map_tmp) <- dim(input_list$par$ln_InitDevs)

    for(p in seq_len(input_list$data$n_pop)) {
      for(r in seq_len(input_list$data$n_regions)) {
        if(r != input_list$data$natal_region[p]) {
          input_list$par$ln_InitDevs[p, r, ] <- 0  # fix starting value
          map_tmp[p, r, ] <- NA                     # turn off estimation
        }
      }
    }

    # Re-index non-NA values sequentially (1, 2, 3, ...)
    non_na <- !is.na(map_tmp)
    map_tmp[non_na] <- as.integer(factor(map_tmp[non_na]))
    input_list$map$ln_InitDevs <- factor(map_tmp)
    collect_message("No dispersal: initial age deviations for non-natal regions fixed to 0 and not estimated.")
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
    if(!is.null(rec_dd) && rec_dd == 'global' && !RecDevs_spec %in% c("est_shared_r", "est_shared_pop_r") && input_list$data$n_regions > 1) stop("Please specify a valid recruitment deviations option for global recruitment density dependence (should be est_shared_r or est_shared_pop_r)!")
    if(!RecDevs_spec %in% c("est_shared_pop_r", "est_shared_r", "fix"))  stop("Please specify a valid recruitment deviations option. These include: fix, est_shared_r. Conversely, leave at NULL to estimate all recruitment deviations.")

    # Share across regions and estimate by population
    if(RecDevs_spec == "est_shared_r") {

      # get indices
      n_yrs <- dim(input_list$par$ln_RecDevs)[3]
      n_region <- dim(input_list$par$ln_RecDevs)[2]

      for(p in 1:input_list$data$n_pop) {
        yr_indices <- (1:n_yrs) + (p-1) * n_yrs # get age indices
        # each region gets the same index for a given age (repeat each index across regions)
        map_RecDevs[p,,] <- matrix(rep(yr_indices, each = n_region), nrow = n_region)
      } # end p loop

      input_list$map$ln_RecDevs <- factor(map_RecDevs)
    } # end if

    # Share across regions and populations
    if(RecDevs_spec == 'est_shared_pop_r') {
      # get indices
      n_yrs <- dim(input_list$par$ln_RecDevs)[3]
      map_RecDevs[] <- rep(1:n_yrs, each = input_list$data$n_regions * input_list$data$n_pop)
      input_list$map$ln_RecDevs <- factor(map_RecDevs)
    }

    # Fix all recruitment deviations
    if(RecDevs_spec == "fix") input_list$map$ln_RecDevs <- factor(rep(NA, prod(dim(map_RecDevs))))

    # print message
    collect_message("Recruitment Deviations is specified as: ", RecDevs_spec)

  } else { # if NULL, estimating all dimensions

    if(input_list$data$n_pop > 1 && input_list$data$rec_region_prop_spec == 1)
      stop("Can't estimate recruitment eviations for all populations and regions if no recruitment dispersal is occuring within a given region! Please specify est_shared_r or est_shared_pop_r instead!")

    input_list$map$ln_RecDevs <- factor(1:length(map_RecDevs)) # input into mapping
    collect_message("Recruitment Deviations is estimated for all dimensions")
  }

  # When no_dispersal, non-natal regions have no recruitment so their
  # deviations are structurally zero. Fix them regardless of RecDevs_spec.
  if(input_list$data$rec_region_prop_spec == 1 && input_list$data$n_pop > 1) {

    # extract mapping
    map_tmp <- as.integer(input_list$map$ln_RecDevs)
    dim(map_tmp) <- dim(input_list$par$ln_RecDevs)

    for(p in seq_len(input_list$data$n_pop)) {
      for(r in seq_len(input_list$data$n_regions)) {
        if(r != input_list$data$natal_region[p]) {
          input_list$par$ln_RecDevs[p, r, ] <- 0  # fix starting value
          map_tmp[p, r, ] <- NA                     # turn off estimation
        }
      }
    }

    # Re-index non-NA values sequentially (1, 2, 3, ...)
    non_na <- !is.na(map_tmp)
    map_tmp[non_na] <- as.integer(factor(map_tmp[non_na]))
    input_list$map$ln_RecDevs <- factor(map_tmp)
    collect_message("No dispersal: Recruitment deviations for non-natal regions fixed to 0 and not estimated.")
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

  # Validate h_spec given rec_dd context
  if(input_list$data$rec_model != 0 && !is.null(rec_dd) && rec_dd == "global") {
    if(is.null(h_spec)) {
      stop("When rec_dd == `global` (global density dependence), h_spec cannot be NULL. ",
           "Steepness must be shared across the global SR relationship: use 'est_shared_pop_r' or 'est_shared_r'.")
    }
    if(!h_spec %in% c("est_shared_pop_r", "est_shared_r", "fix")) {
      stop("When rec_dd == `global` (global density dependence), h_spec must be ",
           "'est_shared_pop_r', 'est_shared_r', or 'fix'.")
    }
  }

  # Mean recruitment
  if(input_list$data$rec_model == 0) {
    input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))
  } else if(!is.null(h_spec)) {

    # Validate options
    if(!h_spec %in% c("est_shared_pop_r", "est_shared_r", "fix")) stop("Please specify a valid steepness option. These include: est_shared_pop_r, fix, est_shared_r. Conversely, leave at NULL to estimate all steepness values.")

    # Share across populations and regions and estimate
    if(h_spec == "est_shared_pop_r") input_list$map$steepness_h <- factor(rep(1, length(input_list$par$steepness_h)))

    # Share across regions but estimate for each population
    if(h_spec == "est_shared_r") input_list$map$steepness_h <- factor(rep(1:input_list$data$n_pop, times = input_list$data$n_regions))

    # Fix all steepness values
    if(h_spec == "fix") input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))

    collect_message("Steepness is specified as: ", h_spec) # output message

  } else {
    # if beverton holt and estimating all steepness parameters
    if(input_list$data$rec_model == 1) {
      # estimate steepness for all populations
      if(input_list$data$n_pop > 1) input_list$map$steepness_h <- factor(rep(1:input_list$data$n_pop, times = input_list$data$n_regions)) # estimating all steepness parameters by  population
      if(input_list$data$n_pop == 1) input_list$map$steepness_h <- factor(1:input_list$data$n_regions) # estimating all steepness parameters by region
    }
    collect_message("Steepness is estimated for all relavant dimensions")
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
  if(!sexratio_spec %in% c("est_all", "est_shared_pop_r", "est_shared_r", "fix")) stop("Sex Ratio Specificaiton is not correctly specified. Needs to be fix, est_all, est_shared_pop_r, or est_shared_r")
  if(input_list$data$n_sexes == 1 && sexratio_spec != 'fix') stop('Sex Ratio is being estiamted, but there is only 1 sex!')

  # Validate whether blocking structure is appropriate
  if(sexratio_spec == 'est_shared_pop_r') {
    ref_blocks <- input_list$data$sexratio_blocks[1,1,]
    for(pp in 1:input_list$data$n_pop) {
      for(rr in 1:input_list$data$n_regions) {
        if(!identical(as.vector(input_list$data$sexratio_blocks[pp,rr,]), as.vector(ref_blocks))) {
          stop("est_shared_pop_r requires consistent sex ratio block structure across all populations and regions.")
        }
      } # end rr loop
    } # end pp loop
  } # end if

  # if we want to fix
  if(sexratio_spec == 'fix') map_sexratio[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(r in 1:input_list$data$n_regions) {

      # Get number of sex ratio rate blocks
      sexratio_blocks_tmp <- unique(as.vector(input_list$data$sexratio_blocks[p,r,]))

      for(b in 1:length(sexratio_blocks_tmp)) {

        # Estimate for all regions
        if(sexratio_spec == 'est_all') {
          if(input_list$data$n_pop > 1 && input_list$data$rec_region_prop_spec == 1)
            stop("Can't estimate recruitment sex ratio for all populations and regions if no recruitment dispersal is occuring within a given region! Please specify est_shared_r or est_shared_pop_r instead!")
          map_sexratio[p,r,b] <- sexratio_counter
          sexratio_counter <- sexratio_counter + 1
        }

        # Estimate but share sex ratio across regions
        if(sexratio_spec == 'est_shared_r' && r == 1) {
          for(rr in 1:input_list$data$n_regions) {
            # only assign if this value exists for this region
            if(sexratio_blocks_tmp[b] %in% input_list$data$sexratio_blocks[p,rr,]) {
              map_sexratio[p,rr, b] <- sexratio_counter
            } # end if
          } # end rr loop
          sexratio_counter <- sexratio_counter + 1
        }

        if(sexratio_spec == 'est_shared_pop_r' && p == 1 && r == 1) {
          for(pp in 1:input_list$data$n_pop) {
            for(rr in 1:input_list$data$n_regions) {
              # only assign if this block exists for this pop/region combo
              if(sexratio_blocks_tmp[b] %in% input_list$data$sexratio_blocks[pp,rr,]) {
                map_sexratio[pp,rr,b] <- sexratio_counter
              }
            } # end rr loop
          } # end pp loop
          sexratio_counter <- sexratio_counter + 1
        }

      } # end b loop
    } # end r loop
  } # end p loop

  collect_message("Sex ratio is specified as: ", sexratio_spec)

  # input sex ratio rates into mapping list
  input_list$map$sexratio_pars <- factor(map_sexratio) # sex ratio rates

  return(input_list)
}

#' Helper function to map recruitment proportions by region
#'
#' @param input_list Input list
#' @keywords internal
do_rec_region_prop_mapping <- function(input_list, rec_region_prop_spec) {

  # Validate spec options
  valid_specs <- c("no_dispersal")
  if(!is.null(rec_region_prop_spec) && !rec_region_prop_spec %in% valid_specs) {
    stop("Invalid rec_region_prop_spec: '", rec_region_prop_spec, "'. Valid options are: ", paste(valid_specs, collapse=", "), ", or NULL to estimate all.")
  }

  # no_dispersal only makes sense with multiple populations
  if(!is.null(rec_region_prop_spec) && rec_region_prop_spec == "no_dispersal" && input_list$data$n_pop == 1 && input_list$data$n_regions == 1) stop("'no_dispersal' is only valid when n_pop > 1 and n_regions > 1.")

  # par is [n_pop, n_regions-1]
  if(!is.null(rec_region_prop_spec) && rec_region_prop_spec == 'no_dispersal') {
    par_mat <- matrix(-20, nrow = input_list$data$n_pop, ncol = input_list$data$n_regions - 1)
    natal_region <- input_list$data$natal_region
    for(p in seq_len(input_list$data$n_pop)) {
      if(natal_region[p] > 1) par_mat[p, natal_region[p] - 1] <- 20
    }
    input_list$par$rec_region_prop_pars <- par_mat # fix values
    input_list$map$rec_region_prop_pars <- factor(rep(NA, length(par_mat)))  # fix all
  }

  # estimate all recruitment propostions if n_regions > 1
  if(is.null(rec_region_prop_spec) && input_list$data$n_regions > 1) input_list$map$rec_region_prop_pars <- factor(1:length(input_list$par$rec_region_prop_pars))

  # single region - not even a parameter
  if(input_list$data$n_regions == 1) {
    input_list$par$rec_region_prop_pars <- NULL
    input_list$map$rec_region_prop_pars <- NULL
  }

  return(input_list)
}

#' Helper function to map recruitment proportions by seasons
#'
#' @param input_list Input list
#' @keywords internal
do_rec_seas_prop_mapping <- function(input_list, rec_seas_prop_spec) {

  # Validate spec options
  valid_specs <- c("fix", "est_shared_p")
  if(!is.null(rec_seas_prop_spec) && !rec_seas_prop_spec %in% valid_specs) {
    stop("Invalid rec_seas_prop_spec: '", rec_seas_prop_spec, "'. Valid options are: ", paste(valid_specs, collapse=", "), ", or NULL to estimate all.")
  }

  if((is.null(rec_seas_prop_spec) || rec_seas_prop_spec == 'est_shared_p') && input_list$data$use_fixed_rec_seas_prop == 1) {
    input_list$data$use_fixed_rec_seas_prop <- 0
    warning("Recruitment seasonal apportionment is specified as estimated, but use_fixed_rec_seas_prop == 1 (fixed). Changing to use_fixed_rec_seas_prop == 0.")
  }

  # estimating recruitment seasonal apporitonment is only valid for seasonal models
  if(!is.null(rec_seas_prop_spec) && rec_seas_prop_spec == 'est_shared_p' && input_list$data$n_seas == 1)
    stop("Estimating recruitment seasonal apportionment is only applicable for seasonal models. ")

  # estimate all recruitment seasonal proportions if n_seas > 1
  if(is.null(rec_seas_prop_spec)) {
    input_list$map$rec_seas_prop_pars <- factor(1:length(input_list$par$rec_seas_prop_pars))
  } else if(rec_seas_prop_spec == 'est_shared_p') { # estimate recruitment seasonal proportions but share across populations
    counter <- 1
    tmp_map = input_list$par$rec_seas_prop_pars
    for(seas in 1:(input_list$data$n_seas - 1)) {
      tmp_map[,seas] <- counter
      counter + 1
    }
    input_list$map$rec_seas_prop_pars <- factor(tmp_map)
  } else if(rec_seas_prop_spec == 'fix') input_list$map$rec_seas_prop_pars <- factor(rep(NA, length(input_list$par$rec_seas_prop_pars)))

  # single seas - not even a parameter
  if(input_list$data$n_seas == 1) {
    input_list$par$rec_seas_prop_pars <- NULL
    input_list$map$rec_seas_prop_pars <- NULL
  }

  return(input_list)
}

#' Setup Recruitment Module and Associated Processes
#'
#' Configures all recruitment-related components of the model, including
#' the recruitment function, density dependence structure, recruitment lag,
#' dispersal processes, steepness estimation and priors, recruitment
#' variability (\eqn{\sigma_R}), recruitment deviations, initial age
#' structure, spawning movement, bias ramp options, and sex ratio dynamics.
#'
#' This function performs the following tasks:
#' \itemize{
#'   \item Validates structural compatibility among recruitment options
#'   \item Initializes recruitment-related parameters in \code{input_list$par}
#'   \item Populates recruitment-related entries in \code{input_list$data}
#'   \item Constructs parameter mapping objects used to control estimation
#' }
#'
#' All recruitment processes must be configured before model compilation.
#'
#' @param input_list Model input list created during earlier setup steps.
#'   Must contain \code{input_list$data} with population, region, age,
#'   year, and season dimensions defined.
#'
#' @section Recruitment Model:
#'
#' @param rec_model Character string specifying the recruitment model:
#' \itemize{
#'   \item \code{"mean_rec"} — Fixed mean recruitment (no stock–recruit relationship)
#'   \item \code{"bh_rec"} — Beverton–Holt stock–recruit relationship
#' }
#' When \code{rec_model = "mean_rec"}, steepness is fixed and not estimated.
#'
#' @param rec_dd Character string specifying the density dependence structure:
#' \itemize{
#'   \item \code{"local"} — Independent stock–recruit relationship per population
#'   \item \code{"global"} — Pooled spawning biomass drives a single SR relationship
#' }
#'
#' If \code{n_pop > 1}, density dependence must be \code{"local"}.
#'
#' @param rec_lag Integer. Lag (in years) between spawning biomass and
#' resulting recruitment.
#'
#' @section Recruitment Variability:
#'
#' @param sigmaR_spec Character controlling estimation of recruitment
#' variability \eqn{\sigma_R}. Options include:
#' \itemize{
#'   \item \code{"est_all"}
#'   \item \code{"est_shared_all"}
#'   \item \code{"fix"}
#'   \item \code{"fix_early_est_late"}
#' }
#'
#' Parameter dimension:
#' \preformatted{
#' ln_sigmaR: [2, n_pop, n_regions]
#' }
#'
#' The first element represents the early-period value and the second
#' represents the late-period value when \code{sigmaR_switch > 1}.
#'
#' @param sigmaR_switch Numeric. Year index at which \eqn{\sigma_R}
#' switches from the early-period value to the late-period value.
#' If \code{<= 1}, only a single value is used.
#'
#' @param RecDevs_spec Character controlling recruitment deviation
#' estimation structure:
#' \itemize{
#'   \item \code{NULL}
#'   \item \code{"est_shared_r"}
#'   \item \code{"est_shared_pop_r"}
#'   \item \code{"fix"}
#' }
#'
#' Parameter dimension:
#' \preformatted{
#' ln_RecDevs: [n_pop, n_regions, n_years]
#' }
#'
#' @param dont_est_recdev_last Integer. Number of terminal years for which
#' recruitment deviations are not estimated.
#'
#' @section Initial Age Structure:
#'
#' @param init_age_strc Method used to initialize the age structure:
#' \itemize{
#'   \item \code{0} — Iterative equilibrium
#'   \item \code{1} — Scalar geometric series (no movement)
#'   \item \code{2} — Matrix geometric series (movement allowed)
#'   \item \code{3} — Scalar geometric series with movement only for the plus group
#' }
#'
#' @param equil_init_age_strc Plus-group treatment for stochastic
#' initialization:
#' \itemize{
#'   \item \code{0} — Deterministic equilibrium
#'   \item \code{1} — Stochastic without plus-group deviations
#'   \item \code{2} — Fully stochastic
#' }
#'
#' @param InitDevs_spec Character controlling estimation structure of
#' initial age deviations:
#' \itemize{
#'   \item \code{NULL}
#'   \item \code{"est_shared_r"}
#'   \item \code{"est_shared_pop_r"}
#'   \item \code{"fix"}
#' }
#'
#' Parameter dimension:
#' \preformatted{
#' ln_InitDevs: [n_pop, n_regions, n_ages - 1]
#' }
#'
#' @section Recruitment Spatial Structure:
#'
#' @param rec_region_prop_spec Recruitment dispersal specification:
#' \itemize{
#'   \item \code{NULL} — Regional recruitment proportions are estimated
#'   \item \code{"no_dispersal"} — Recruitment fixed to natal regions
#' }
#'
#' Parameter dimension:
#' \preformatted{
#' rec_region_prop_pars: [n_pop, n_regions - 1]
#' }
#'
#' @param use_rec_region_prop_prior Integer (0/1). Whether Dirichlet priors
#' are applied to recruitment regional proportions.
#'
#' @param rec_region_prop_prior Optional data frame specifying Dirichlet
#' prior concentration parameters. Must contain columns:
#' \code{pop} and \code{alpha}, where \code{alpha} is a list-column
#' containing length-\code{n_regions} vectors.
#'
#' @section Recruitment Seasonal Structure:
#'
#' @param rec_seas_prop_spec Character controlling estimation of seasonal
#' recruitment proportions.
#'
#' @param use_rec_seas_prop_prior Integer (0/1). Whether priors are applied
#' to seasonal recruitment proportions.
#'
#' @param rec_seas_prop_prior Optional data frame defining Dirichlet prior
#' concentration parameters for seasonal recruitment proportions.
#'
#' @param use_fixed_rec_seas_prop Integer (0/1). Whether fixed seasonal
#' recruitment proportions are used.
#'
#' @param fixed_rec_seas_prop Array specifying fixed seasonal recruitment
#' proportions with dimension:
#' \preformatted{
#' [n_pop, n_seas]
#' }
#'
#' @section Steepness:
#'
#' @param h_spec Character controlling steepness parameter mapping:
#' \itemize{
#'   \item \code{NULL}
#'   \item \code{"est_shared_r"}
#'   \item \code{"est_shared_pop_r"}
#'   \item \code{"fix"}
#' }
#'
#' @param Use_h_prior Integer (0/1). Whether steepness priors are applied.
#'
#' @param h_prior Data frame specifying steepness prior information.
#' Required columns:
#' \code{pop}, \code{region}, \code{mu}, \code{sd}.
#'
#' Steepness parameter dimension:
#' \preformatted{
#' steepness_h: [n_pop, n_regions]
#' }
#'
#' @section Spawning Processes:
#'
#' @param spawn_seas Integer. Season index in which spawning occurs.
#'
#' @param t_spawn Numeric. Fraction of the year at which spawning occurs.
#'
#' @param sgl_seas_spawning_movement Optional spawning movement array:
#' \preformatted{
#' [n_pop, n_regions, n_regions, n_years, n_ages, n_sexes]
#' }
#'
#' If \code{NA}, 100\% natal homing is assumed.
#'
#' @param stray_rate Array of stray rates with dimension:
#' \preformatted{
#' [n_pop, n_years]
#' }
#'
#' Values must lie between 0 and 1.
#'
#' @section Sex Ratio Dynamics:
#'
#' @param sexratio_spec Character controlling sex ratio estimation:
#' \itemize{
#'   \item \code{"est_all"}
#'   \item \code{"est_shared_r"}
#'   \item \code{"est_shared_pop_r"}
#'   \item \code{"fix"}
#' }
#'
#' @param sexratio_blocks Character vector specifying temporal block
#' structure. Valid formats include:
#' \itemize{
#'   \item \code{"none_Pop_x_Region_x"}
#'   \item \code{"Block_k_Year_a-b_Pop_x_Region_x"}
#' }
#'
#' Parameter dimension:
#' \preformatted{
#' sexratio_pars: [n_pop, n_regions, n_blocks]
#' }
#'
#' @section Bias Ramp:
#'
#' @param do_rec_bias_ramp Integer (0/1). Whether a recruitment bias ramp
#' is applied.
#'
#' @param bias_year Numeric. Year at which the bias ramp begins.
#'
#' @param max_bias_ramp_fct Numeric in [0,1]. Maximum bias correction factor.
#'
#' @section Additional Inputs:
#'
#' @param init_F_prop Numeric vector of length \code{n_seas} specifying
#' the seasonal distribution of initial fishing mortality.
#'
#' @param ... Optional named starting values for parameters. Supported names:
#' \itemize{
#'   \item \code{ln_global_R0}
#'   \item \code{rec_region_prop_pars}
#'   \item \code{rec_seas_prop_pars}
#'   \item \code{steepness_h}
#'   \item \code{ln_InitDevs}
#'   \item \code{ln_RecDevs}
#'   \item \code{ln_sigmaR}
#'   \item \code{sexratio_pars}
#' }
#'
#' @section Structural Compatibility Rules:
#'
#' \strong{Multiple populations:}
#' If \code{n_pop > 1}, recruitment density dependence must be local.
#'
#' \strong{Global density dependence:}
#' When \code{rec_dd = "global"}, steepness and deviation parameters must
#' be shared or fixed across spatial units.
#'
#' \strong{No dispersal constraint:}
#' When \code{rec_region_prop_spec = "no_dispersal"} and
#' \code{n_pop > 1}, recruitment and initial deviations cannot be
#' fully independent across populations and regions.
#'
#' @family Model Setup
#' @export Setup_Mod_Rec
Setup_Mod_Rec <- function(input_list,
                          rec_model,
                          rec_dd = "global",
                          rec_lag = 1,
                          Use_h_prior = 0,
                          h_prior = NULL,
                          rec_region_prop_spec = NULL,
                          use_rec_region_prop_prior = 0,
                          rec_region_prop_prior = NULL,
                          rec_seas_prop_spec = "fix",
                          use_rec_seas_prop_prior = 0,
                          rec_seas_prop_prior = NULL,
                          use_fixed_rec_seas_prop = 1,
                          fixed_rec_seas_prop = array(rep(c(1, rep(0, input_list$data$n_seas - 1)),
                                                          each = input_list$data$n_seas), dim = c(input_list$data$n_pop, input_list$data$n_seas)),
                          do_rec_bias_ramp = 0,
                          bias_year = NA,
                          max_bias_ramp_fct = 1,
                          sigmaR_switch = 1,
                          dont_est_recdev_last = 0,
                          init_age_strc = 2,
                          equil_init_age_strc = 1,
                          init_F_prop = rep(0, input_list$data$n_seas),
                          sigmaR_spec = "est_all",
                          InitDevs_spec = NULL,
                          RecDevs_spec = NULL,
                          h_spec = NULL,
                          sgl_seas_spawning_movement = NA,
                          t_spawn = 0,
                          stray_rate = array(1, dim = c(input_list$data$n_pop, length(input_list$data$years))),
                          spawn_seas = 1,
                          sexratio_spec = 'fix',
                          sexratio_blocks = {
                            grid <- expand.grid(region = 1:input_list$data$n_regions, pop = 1:input_list$data$n_pop)
                            blks <- paste0("none_Pop_", grid$pop, "_Region_", grid$region)
                            blks
                          },
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

  # Recruitment regional proportion prior
  if(!use_rec_region_prop_prior %in% c(0,1)) stop("use_rec_region_prop_prior must be 0 or 1")
  if(use_rec_region_prop_prior == 1 && input_list$data$n_regions == 1) stop("Priors should not be applied to recruitment regional proportions when n_regions = 1.")
  if(use_rec_region_prop_prior == 1) {
    required_cols <- c("pop", "alpha")
    missing_cols <- setdiff(required_cols, names(rec_region_prop_prior))
    if (length(missing_cols) > 0) stop("rec_region_prop_prior is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  collect_message("Recruitment regional proportion priors are: ", ifelse(use_rec_region_prop_prior == 1, "Used", "Not Used"))

  # recruitment seasonal priors
  if(!use_rec_seas_prop_prior %in% c(0,1)) stop("use_rec_seas_prop_prior must be 0 or 1")
  if(use_rec_seas_prop_prior == 1 && input_list$data$n_seas == 1) stop("Priors should not be applied to recruitment seasonal proportions when n_seass = 1.")
  if(use_rec_seas_prop_prior == 1) {
    required_cols <- c("pop", "alpha")
    missing_cols <- setdiff(required_cols, names(rec_seas_prop_prior))
    if (length(missing_cols) > 0) stop("rec_seas_prop_prior is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  collect_message("Recruitment seasonal proportion priors are: ", ifelse(use_rec_seas_prop_prior == 1, "Used", "Not Used"))
  collect_message("Recruitment seasonal proportions is: ", ifelse(is.null(rec_seas_prop_spec), "estimated for all dimensions", rec_seas_prop_spec))

  # Checking that rec_dd is local when n_pop > 1
  if(input_list$data$n_pop > 1 && rec_dd != 'local') stop("When n_pop > 1, rec_dd must be local!")


  # Spawning Movement -------------------------------------------------------
  if(is.na(sum(sgl_seas_spawning_movement))) {
    natal_region <- input_list$data$natal_region
    arr <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))
    for(p in seq_len(input_list$data$n_pop)) arr[p, , natal_region[p], , , ] <- 1
    tmp_sgl_seas_spawning_movement <- arr
    if(input_list$data$n_pop > 1 && input_list$data$n_seas == 1 && rec_model_val == 1) collect_message("Using 100% natal homing rate.")
  } else {
    check_data_dimensions(sgl_seas_spawning_movement, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years),
                          n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'sgl_seas_spawning_movement')
    tmp_sgl_seas_spawning_movement <- sgl_seas_spawning_movement
    if(input_list$data$n_pop > 1 && input_list$data$n_seas == 1 && rec_model_val == 1) collect_message("Using user input natal homing rate.")
  }


  # Straying Rates ----------------------------------------------------------
  check_data_dimensions(stray_rate, n_pop = input_list$data$n_pop, n_years = length(input_list$data$years),  what = 'stray_rate')


  # Steepness Settings ------------------------------------------------------
  if (rec_model == "bh_rec") {
    if (!Use_h_prior %in% c(0, 1)) stop("Use_h_prior must be 0 or 1")
    if (Use_h_prior == 1) {
      required_cols <- c("pop", "region", "mu", "sd")
      missing_cols <- setdiff(required_cols, names(h_prior))
      if (length(missing_cols) > 0) stop("h_prior is missing columns: ", paste(missing_cols, collapse = ", "))
    }
    collect_message("Steepness priors are: ", ifelse(Use_h_prior == 1, "Used", "Not Used"))
  }


  # Sex Ratio Options ---------------------------------------------
  sexratio_blocks_mat <- array(NA, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years)))

  for(i in 1:length(sexratio_blocks)) {

    # Extract out components from list
    tmp <- sexratio_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    if(!tmp_vec[1] %in% c("none", "Block")) stop("Sex Ratio Blocks not correctly specified. This should be either none_Pop_x_Region_x or Block_x_Year_x-y_Pop_x_Region_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      pop <- as.numeric(tmp_vec[3]) # get pop index
      region <- as.numeric(tmp_vec[5]) # get region index
      sexratio_blocks_mat[pop, region,] <- 1 # input sex ratio time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      pop <- as.numeric(tmp_vec[6]) # get pop value
      region <- as.numeric(tmp_vec[8]) # get region value

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      sexratio_blocks_mat[pop,region,years] <- block_val # input sex ratio time block
    }

  } # end i loop

  for(p in 1:input_list$data$n_pop) for(r in 1:input_list$data$n_regions)
    collect_message("Sex Ratios specified with ", length(unique(sexratio_blocks_mat[p,r,])), " block for population ", p, " and region ", r)

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

  # # input variables into data list
  collect_message("Spawning season occurs in season ", spawn_seas)
  input_list$data$spawn_seas <- spawn_seas
  input_list$data$sgl_seas_spawning_movement <- tmp_sgl_seas_spawning_movement
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
  input_list$data$use_rec_region_prop_prior <- use_rec_region_prop_prior
  input_list$data$rec_region_prop_spec <- ifelse(is.null(rec_region_prop_spec), 0, 1) # 0 = Full dispersal, 1 = no dispersal
  input_list$data$rec_region_prop_prior <- rec_region_prop_prior
  input_list$data$stray_rate <- stray_rate
  input_list$data$sexratio_blocks <- sexratio_blocks_mat
  input_list$data$use_fixed_rec_seas_prop <- use_fixed_rec_seas_prop
  input_list$data$fixed_rec_seas_prop <- fixed_rec_seas_prop
  input_list$data$use_rec_seas_prop_prior <- use_rec_seas_prop_prior
  input_list$data$rec_seas_prop_prior <- rec_seas_prop_prior

  # Populate Parameter List -------------------------------------------------

  # Global R0
  if("ln_global_R0" %in% names(starting_values)) input_list$par$ln_global_R0 <- starting_values$ln_global_R0
  else input_list$par$ln_global_R0 <- array(log(15), dim = c(input_list$data$n_pop))

  # R0 regional proportion (not availiable when n_regions == 1; altered in do_rec_region_prop_mapping)
  if("rec_region_prop_pars" %in% names(starting_values)) input_list$par$rec_region_prop_pars <- starting_values$rec_region_prop_pars
  else input_list$par$rec_region_prop_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions - 1))

  # R0 seasonal proportion (not availiable when n_seas == 1; altered in do_rec_seas_prop_mapping)
  if("rec_seas_prop_pars" %in% names(starting_values)) input_list$par$rec_seas_prop_pars <- starting_values$rec_seas_prop_pars
  else input_list$par$rec_seas_prop_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_seas - 1))

  # Steepness in bounded logit space (0.2 and 1)
  if("steepness_h" %in% names(starting_values)) input_list$par$steepness_h <- starting_values$steepness_h
  else input_list$par$steepness_h <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions))

  # Initial age deviations
  if("ln_InitDevs" %in% names(starting_values)) input_list$par$ln_InitDevs <- starting_values$ln_InitDevs
  else input_list$par$ln_InitDevs <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$ages) - 1))

  # Recruitment deviations
  if("ln_RecDevs" %in% names(starting_values)) input_list$par$ln_RecDevs <- starting_values$ln_RecDevs
  else input_list$par$ln_RecDevs <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years) - dont_est_recdev_last + input_list$data$n_proj_yrs_devs))

  # Recruitment variability
  if("ln_sigmaR" %in% names(starting_values)) input_list$par$ln_sigmaR <- starting_values$ln_sigmaR
  else input_list$par$ln_sigmaR <- array(log(1), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)) # (early period 1st element, late period 2nd element)

  # sexratio parameters
  max_sexratio_blks <- max(apply(input_list$data$sexratio_blocks, c(1,2), FUN = function(x) length(unique(x)))) # figure out maximum number ofsex ratio blocks for each region
  if("sexratio_pars" %in% names(starting_values)) input_list$par$sexratio_pars <- starting_values$sexratio_pars
  else input_list$par$sexratio_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, max_sexratio_blks)) # specified at 0.5 in inverse logit space

  # Mapping Options -----------------------------------------------------------

  input_list <- do_rec_region_prop_mapping(input_list, rec_region_prop_spec) # Recruitment regional proportion mapping
  input_list <- do_rec_seas_prop_mapping(input_list, rec_seas_prop_spec) # Recruitment seasonal proportion mapping
  input_list <- do_sigmaR_mapping(input_list, sigmaR_spec) # sigmaR mapping
  input_list <- do_InitDevs_mapping(input_list, InitDevs_spec, rec_dd) # InitDevs mapping
  input_list <- do_RecDevs_mapping(input_list, RecDevs_spec, rec_dd) # RevDevs mapping
  input_list <- do_h_mapping(input_list, h_spec, rec_dd) # steepness mapping
  input_list <- do_sexratio_pars_mapping(input_list, sexratio_spec) # sex ratio parameters

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}
