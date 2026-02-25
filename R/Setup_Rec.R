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

    if(input_list$data$n_pop > 1 && input_list$data$Rec_prop_spec == 1)
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

    if(input_list$data$n_pop > 1 && input_list$data$Rec_prop_spec == 1)
      stop("Can't estimate recruitment eviations for all populations and regions if no recruitment dispersal is occuring within a given region! Please specify est_shared_r or est_shared_pop_r instead!")

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
          if(input_list$data$n_pop > 1 && input_list$data$Rec_prop_spec == 1)
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

#' Helper function to map recruitment proportions
#'
#' @param input_list Input list
#' @keywords internal
do_Rec_prop_mapping <- function(input_list, Rec_prop_spec) {

  # Validate spec options
  valid_specs <- c("no_dispersal")
  if(!is.null(Rec_prop_spec) && !Rec_prop_spec %in% valid_specs) {
    stop("Invalid Rec_prop_spec: '", Rec_prop_spec, "'. Valid options are: ", paste(valid_specs, collapse=", "), ", or NULL to estimate all.")
  }

  # no_dispersal only makes sense with multiple populations
  if(!is.null(Rec_prop_spec) && Rec_prop_spec == "no_dispersal" && input_list$data$n_pop == 1) stop("'no_dispersal' is only valid when n_pop > 1.")

  # single region - map off entirely, transform returns 1
  if(input_list$data$n_regions == 1) {
    input_list$map$Rec_prop <- factor(rep(NA, length(input_list$par$Rec_prop)))
  }

  # par is [n_pop, n_regions-1]
  if(!is.null(Rec_prop_spec) && Rec_prop_spec == "no_dispersal" && input_list$data$n_regions == 1)
    stop("'no_dispersal' requires n_regions > 1.")
  if(!is.null(Rec_prop_spec) && Rec_prop_spec == 'no_dispersal') {
    par_mat <- matrix(-20, nrow = input_list$data$n_pop, ncol = input_list$data$n_regions - 1)
    for(p in seq_len(input_list$data$n_pop)) if(p > 1) par_mat[p, p - 1] <- 20  # concentrate weight on natal region p
    input_list$par$Rec_prop <- par_mat # fix values
    input_list$map$Rec_prop <- factor(rep(NA, length(par_mat)))  # fix all
    collect_message("No dispersal: recruitment fixed to natal regions.")
  }

  # estimate all recruitment propostions
  if(is.null(Rec_prop_spec)) input_list$map$Rec_prop <- factor(1:length(input_list$par$Rec_prop))

  return(input_list)
}

#' Setup Recruitment Module and Associated Processes
#'
#' Configures all recruitment-related components of the model, including
#' recruitment form, density dependence structure, recruitment deviations,
#' steepness, recruitment dispersal, initial age structure, recruitment
#' variability, spawning movement, and sex ratio dynamics.
#'
#' This function initializes parameter arrays, validates option compatibility,
#' and constructs mapping objects for parameter estimation.
#'
#' @param rec_model Character string specifying the recruitment model:
#' \itemize{
#'   \item \code{"mean_rec"}: Recruitment is a fixed mean value (no stock–recruit relationship).
#'   \item \code{"bh_rec"}: Beverton–Holt recruitment with steepness parameter \code{h}.
#' }
#' If \code{rec_model = "mean_rec"}, steepness is fixed automatically.
#'
#' @param rec_dd Character string specifying recruitment density dependence:
#' \itemize{
#'   \item \code{"local"}: Separate stock-recruit relationship per population/region.
#'   \item \code{"global"}: Single pooled SSB drives one global SR relationship.
#' }
#' When \code{n_pop > 1}, \code{rec_dd} must be \code{"local"}.
#' When \code{rec_dd = "global"} and \code{n_regions > 1}, all shared deviation and
#' steepness specs (\code{RecDevs_spec}, \code{InitDevs_spec}, \code{h_spec}) must
#' use \code{"est_shared_r"}, \code{"est_shared_pop_r"}, or \code{"fix"} —
#' \code{NULL} is not permitted for any of these.
#'
#' @param rec_lag Integer specifying lag between spawning biomass and recruitment.
#'
#' @param Rec_prop_spec Character string controlling recruitment dispersal:
#' \itemize{
#'   \item \code{NULL}: Estimate all recruitment proportions if appropriate.
#'   \item \code{"no_dispersal"}: Fix recruitment to natal regions only.
#' }
#' Recruitment proportion parameter dimension:
#' \preformatted{
#' Rec_prop: [n_pop, n_regions - 1]
#' }
#' When \code{n_regions == 1}, recruitment proportions are fixed.
#' The option \code{"no_dispersal"} is only valid when \code{n_pop > 1}.
#'
#' @param sigmaR_spec Character string controlling mapping of recruitment
#' variability (\eqn{\sigma_R}). Valid options are:
#' \itemize{
#'   \item \code{"est_all"}: Estimate separate \eqn{\sigma_R} values for each
#'   period (early/late) and for each population. If
#'   \code{n_pop == 1} and recruitment density dependence is local,
#'   \eqn{\sigma_R} may be estimated separately by region. Otherwise,
#'   values are shared across regions within population. Default.
#'
#'   \item \code{"est_shared_all"}: Estimate a single \eqn{\sigma_R}
#'   shared across early and late periods, populations, and regions.
#'
#'   \item \code{"fix"}: Fix all \eqn{\sigma_R} values at their initial values.
#'
#'   \item \code{"fix_early_est_late"}: Fix the early-period \eqn{\sigma_R}
#'   and estimate the late-period \eqn{\sigma_R}. The definition of
#'   "early" and "late" is controlled by \code{sigmaR_switch}.
#' }
#'
#' The internal parameter array has dimension:
#' \preformatted{
#' ln_sigmaR: [2, n_pop, n_regions]
#' }
#'
#' The first dimension indexes:
#' \itemize{
#'   \item 1 = early period
#'   \item 2 = late period
#' }
#'
#' When \code{sigmaR_switch <= 1}, only a single period is effectively used,
#' but the array retains two elements for structural consistency.
#'
#' If \code{rec_dd = "global"} and \code{n_regions > 1},
#' \eqn{\sigma_R} is shared across regions within population.
#'
#' @param InitDevs_spec Character string controlling mapping of initial age deviations:
#' \itemize{
#'   \item \code{NULL}: Estimate deviations for all populations, regions, and ages.
#'   \item \code{"est_shared_r"}: Share deviations across regions within population.
#'   \item \code{"est_shared_pop_r"}: Share deviations across populations and regions.
#'   \item \code{"fix"}: Fix all deviations.
#' }
#'
#' Initial age deviation parameter dimension:
#' \preformatted{
#' ln_InitDevs: [n_pop, n_regions, n_ages - 1]
#' }
#'
#' If \code{equil_init_age_strc = 1}, the plus group follows equilibrium
#' and the final age index is fixed (not estimated).
#'
#' If \code{rec_dd = "global"} and \code{n_regions > 1},
#' only shared specifications (\code{"est_shared_r"} or
#' \code{"est_shared_pop_r"}) are valid.
#'
#' If \code{n_pop > 1} and recruitment dispersal is fixed to natal regions
#' (\code{Rec_prop_spec = "no_dispersal"}), full estimation across all
#' populations and regions (\code{InitDevs_spec = NULL}) is not allowed.
#' In this case, deviations must be shared using
#' \code{"est_shared_r"} or \code{"est_shared_pop_r"}.
#'
#' @param RecDevs_spec Character string controlling recruitment deviation mapping:
#' \itemize{
#'   \item \code{NULL}: Estimate recruitment deviations for all populations,
#'   regions, and years.
#'   \item \code{"est_shared_r"}: Share deviations across regions within
#'   population (dimension reduces to [n_pop, 1, n_years]).
#'   \item \code{"est_shared_pop_r"}: Share deviations across populations
#'   and regions (dimension reduces to [1, 1, n_years]).
#'   \item \code{"fix"}: Fix all recruitment deviations to zero.
#' }
#'
#' Recruitment deviation parameter dimension:
#' \preformatted{
#' ln_RecDevs: [n_pop, n_regions, n_years]
#' }
#'
#' If \code{n_pop > 1} and recruitment dispersal is fixed to natal regions
#' (\code{Rec_prop_spec = "no_dispersal"}), full estimation across all
#' populations and regions (\code{RecDevs_spec = NULL}) is not allowed.
#' In this case, deviations must be shared using
#' \code{"est_shared_r"} or \code{"est_shared_pop_r"}.
#'
#' Under \code{rec_dd = "global"} with multiple regions, recruitment
#' deviations should be shared across regions to maintain consistency
#' with the density dependence structure.
#'
#'
#' @param h_spec Character string controlling steepness estimation:
#' \itemize{
#'   \item \code{NULL}: Estimate steepness for all relevant dimensions. When
#'     \code{n_pop > 1}, shares across regions and estimates per population
#'     (equivalent to \code{"est_shared_r"}). When \code{n_pop == 1}, estimates
#'     separately by region. \strong{Not valid when \code{rec_dd = "global"}}.
#'   \item \code{"est_shared_r"}: Share steepness across regions, estimate per population.
#'   \item \code{"est_shared_pop_r"}: Single shared steepness across all populations and regions.
#'   \item \code{"fix"}: Fix steepness to starting values (not estimated).
#' }
#' When \code{rec_model = "mean_rec"}, steepness is fixed automatically regardless of \code{h_spec}.
#'
#' Valid \code{h_spec} by context:
#' \tabular{lll}{
#'   \strong{rec_dd} \tab \strong{n_pop} \tab \strong{Valid h_spec} \cr
#'   \code{"global"} \tab any \tab \code{"est_shared_pop_r"}, \code{"est_shared_r"}, \code{"fix"} \cr
#'   \code{"local"} \tab 1 \tab any, including \code{NULL} \cr
#'   \code{"local"} \tab > 1 \tab any, including \code{NULL} \cr
#' }
#'
#' Steepness parameter dimension:
#' \preformatted{steepness_h: [n_pop, n_regions]}
#'
#' @param sgl_seas_spawning_movement Array specifying spawning movement
#' in single-season multi-population Beverton–Holt models.
#'
#' Dimension:
#' \preformatted{
#' [n_pop, n_regions, n_regions, n_years, n_ages, n_sexes]
#' }
#'
#' The second dimension corresponds to origin region and the third
#' dimension corresponds to spawning region.
#'
#' Values represent the proportion of individuals moving from origin
#' region \eqn{r} to spawning region \eqn{r'} in year \eqn{y}.
#'
#' If not supplied, the default assumption is 100% natal homing, such that
#' population \eqn{i} spawns only in region \eqn{i}. Under this default,
#' the effective spawning structure reduces to population-specific
#' self-recruitment only.
#'
#' If recruitment dispersal is also fixed to natal regions
#' (\code{Rec_prop_spec = "no_dispersal"}), spawning movement and
#' recruitment assignment are structurally aligned.
#'
#'
#' @param sexratio_spec Character string specifying sex ratio estimation:
#' \itemize{
#'   \item \code{"est_all"}: Estimate sex ratio for all populations,
#'   regions, and blocks.
#'   \item \code{"est_shared_r"}: Share sex ratio across regions within
#'   population.
#'   \item \code{"est_shared_pop_r"}: Share sex ratio across populations
#'   and regions.
#'   \item \code{"fix"}: Fix sex ratio to input values (not estimated).
#' }
#'
#' Sex ratio parameter dimension:
#' \preformatted{
#' sexratio_pars: [n_pop, n_regions, n_blocks]
#' }
#'
#' If sex ratio is shared across populations and regions
#' (\code{"est_shared_pop_r"}), the block structure
#' (\code{n_blocks} and block timing) must be identical across all
#' populations and regions.
#'
#' Sharing across regions is generally recommended when recruitment
#' or density dependence is specified at a global scale.
#'
#' @description
#' \strong{Parameter compatibility quick reference:}
#'
#' When \code{rec_dd = "global"}, \code{h_spec}, \code{RecDevs_spec}, and
#' \code{InitDevs_spec} must all be shared or fixed. \code{NULL} is not allowed
#' for any of these because a single pooled SR relationship is inconsistent with
#' separately estimated regional/population parameters.
#'
#' \tabular{llll}{
#'   \strong{rec_dd} \tab \strong{n_pop} \tab \strong{h_spec} \tab \strong{RecDevs/InitDevs_spec} \cr
#'   \code{"global"} \tab 1 \tab shared or fix only \tab shared or fix only \cr
#'   \code{"local"} \tab 1 \tab any incl. NULL \tab any incl. NULL \cr
#'   \code{"local"} \tab > 1 \tab any incl. NULL \tab NULL blocked if \code{Rec_prop_spec = "no_dispersal"} \cr
#' }
#'
#' @export Setup_Mod_Rec
#' @family Model Setup
Setup_Mod_Rec <- function(input_list,
                          rec_model,
                          rec_dd = "global",
                          rec_lag = 1,
                          Use_h_prior = 0,
                          h_prior = NULL,
                          Rec_prop_spec = NULL,
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
                          sigmaR_spec = "est_all",
                          InitDevs_spec = NULL,
                          RecDevs_spec = NULL,
                          h_spec = NULL,
                          sgl_seas_spawning_movement = NA,
                          t_spawn = 0,
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

  # Recruitment proportion prior
  if(!Use_Rec_prop_Prior %in% c(0,1)) stop("Use_Rec_prop_Prior must be 0 or 1")
  if(Use_Rec_prop_Prior == 1 && input_list$data$n_regions == 1) stop("Priors should not be applied to recruitment proportions when n_regions = 1.")
  collect_message("Recruitment proportion priors are: ", ifelse(Use_Rec_prop_Prior == 1, "Used", "Not Used"))

  # Checking that rec_dd is local when n_pop > 1
  if(input_list$data$n_pop > 1 && rec_dd != 'local') stop("When n_pop > 1, rec_dd must be local!")


  # Spawning Movement -------------------------------------------------------
  if(is.na(sum(sgl_seas_spawning_movement))) {
    arr <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))
    for (p in seq_len(input_list$data$n_pop)) arr[p, , p, , , ] <- 1 # natal homing with 100% probability
    tmp_sgl_seas_spawning_movement <- arr
    if(input_list$data$n_pop > 1 && input_list$data$n_seas == 1 && rec_model_val == 1) collect_message("Using 100% natal homing rate.")
  } else {
    check_data_dimensions(sgl_seas_spawning_movement, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years),
                          n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'sgl_seas_spawning_movement')
    tmp_sgl_seas_spawning_movement <- sgl_seas_spawning_movement
    if(input_list$data$n_pop > 1 && input_list$data$n_seas == 1 && rec_model_val == 1) collect_message("Using user input natal homing rate.")
  }

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
    collect_message("Sex Ratios estimated with ", length(unique(sexratio_blocks_mat[p,r,])), " block for population ", p, " and region ", r)

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
  input_list$data$Use_Rec_prop_Prior <- Use_Rec_prop_Prior
  input_list$data$Rec_prop_spec <- ifelse(is.null(Rec_prop_spec), 0, 1) # 0 = Full dispersal, 1 = no dispersal
  Rec_prior_vals = ifelse(is.null(Rec_prop_prior), array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions)), Rec_prop_prior)
  input_list$data$Rec_prop_prior <- array(Rec_prior_vals, dim = c(input_list$data$n_pop, input_list$data$n_regions))
  input_list$data$sexratio_blocks <- sexratio_blocks_mat

  # Populate Parameter List -------------------------------------------------

  # Global R0
  if("ln_global_R0" %in% names(starting_values)) input_list$par$ln_global_R0 <- starting_values$ln_global_R0
  else input_list$par$ln_global_R0 <- array(log(15), dim = c(input_list$data$n_pop))

  # R0 proportion
  if("Rec_prop" %in% names(starting_values)) input_list$par$Rec_prop <- starting_values$Rec_prop
  else input_list$par$Rec_prop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions - 1))

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

  input_list <- do_Rec_prop_mapping(input_list, Rec_prop_spec) # Recruitment proportion mapping
  input_list <- do_sigmaR_mapping(input_list, sigmaR_spec) # sigmaR mapping
  input_list <- do_InitDevs_mapping(input_list, InitDevs_spec, rec_dd) # InitDevs mapping
  input_list <- do_RecDevs_mapping(input_list, RecDevs_spec, rec_dd) # RevDevs mapping
  input_list <- do_h_mapping(input_list, h_spec, rec_dd) # steepness mapping
  input_list <- do_sexratio_pars_mapping(input_list, sexratio_spec) # sex ratio parameters

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}

