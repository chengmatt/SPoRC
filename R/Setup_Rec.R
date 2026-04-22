#' Set up recruitment dynamics for the operating model simulation
#'
#' Populates \code{sim_list} with all recruitment-related inputs needed by the
#' operating model: stock-recruit relationship type and density-dependence
#' structure, biological parameters (\eqn{R_0}, steepness, sex ratio),
#' recruitment and initial age-structure deviations, seasonal recruitment
#' allocation, spawn timing, and equilibrium initialisation method. Must be
#' called after \code{\link{Setup_Sim_Dim}}.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}}.
#' @param recruitment_opt Recruitment model. Default \code{"bh_rec"}. Options:
#'   \describe{
#'     \item{\code{0}/\code{"mean_rec"}}{Mean recruitment; no stock-recruit
#'       relationship.}
#'     \item{\code{1}/\code{"bh_rec"}}{Beverton-Holt stock-recruit relationship.
#'       Requires \code{rec_dd = "local"} when \code{n_pop > 1}.}
#'     \item{\code{999}/\code{"resample_from_input"}}{Bootstrap recruitment by
#'       resampling years from \code{Rec_input}. Historical years are used
#'       as-is; projection years are sampled with replacement, preserving
#'       spatial covariance among regions within each resampled year.
#'       Requires \code{Rec_input} to be non-\code{NULL}.}
#'   }
#' @param rec_dd Density-dependence structure for the stock-recruit
#'   relationship. Default \code{"global"}. Options:
#'   \describe{
#'     \item{\code{0}/\code{"local"}}{Region-specific spawner-recruit
#'       relationship; each region has its own \eqn{R_0} and steepness.
#'       Required when \code{n_pop > 1} and \code{recruitment_opt = "bh_rec"}.}
#'     \item{\code{1}/\code{"global"}}{Single shared spawner-recruit
#'       relationship pooled across regions.}
#'   }
#' @param init_dd Density-dependence structure for equilibrium age-structure
#'   initialisation. Default \code{"global"}. Same options as \code{rec_dd}.
#' @param R0_input Unfished equilibrium recruitment array
#'   \code{[n_pop x n_regions x n_yrs x n_sims]}. Default: \code{10} for all
#'   cells.
#' @param h_input Beverton-Holt steepness array
#'   \code{[n_pop x n_regions x n_yrs x n_sims]}. Values should be in
#'   \eqn{(0.2, 1)}. Default: \code{0.8}.
#' @param sexratio_input Proportion of recruits assigned to each sex, array
#'   \code{[n_pop x n_regions x n_yrs x n_sexes x n_sims]}. Default: \code{1}
#'   when \code{n_sexes = 1}; \code{0.5} per sex when \code{n_sexes = 2}.
#' @param ln_sigmaR Log-scale standard deviation of recruitment deviations,
#'   array \code{[n_pop x n_regions]}. The first element controls the SD for
#'   initial age-structure deviations (\code{ln_InitDevs}); the second controls
#'   the SD for annual recruitment deviations (\code{ln_RecDevs}). Default:
#'   \code{log(1)} for both.
#' @param rec_seas_prop_input Seasonal allocation of annual recruitment, array
#'   \code{[n_pop x n_seas x n_sims]}. Each population's values should sum to
#'   1 across seasons. Default: all recruitment assigned to season 1.
#' @param spawn_seas Integer index of the season in which spawning occurs.
#'   Default \code{1}.
#' @param t_spawn Spawn timing as a fraction of the season elapsed before
#'   spawning within \code{spawn_seas}. \code{0} (default) = spawning occurs
#'   before any mortality is applied in that season; \code{1} = spawning occurs
#'   after all mortality.
#' @param rec_lag Integer. Number of seasons between spawning and recruitment
#'   of age-1 fish. Must be \eqn{\geq 1}; \code{0} is not permitted. Default
#'   \code{1}.
#' @param init_age_strc Integer specifying the equilibrium age-structure
#'   initialisation method. Default \code{2}. Options:
#'   \describe{
#'     \item{\code{0}/\code{"iterative"}}{Iterates the population forward until
#'       approximate equilibrium. Slowest but most general.}
#'     \item{\code{1}/\code{"scalar_no_move"}}{Scalar geometric series solution
#'       assuming no movement in any age class.}
#'     \item{\code{2}/\code{"matrix"}}{Matrix geometric series solution that
#'       generalises the scalar approach to include movement. Recommended
#'       default for spatially explicit models.}
#'     \item{\code{3}/\code{"scalar_plus_only"}}{Scalar geometric series
#'       solution assuming no movement except in the plus group.}
#'   }
#' @param do_recruits_move Integer flag. \code{0} = age-1 fish do not move
#'   (default); movement begins at age 2. \code{1} = recruits participate in
#'   movement from age 1.
#' @param stray_rate_input Natal-homing stray rate array
#'   \code{[n_pop x n_yrs x n_sims]}. Proportion of individuals that stray
#'   from their natal region during spawning. Default: \code{0} (No individuals stray).
#' @param Rec_input External recruitment array
#'   \code{[n_pop x n_regions x n_yrs x n_sims]}. Required when
#'   \code{recruitment_opt = "resample_from_input"}; projection years beyond
#'   the length of \code{Rec_input} are filled by resampling historical years
#'   with replacement. Ignored for other recruitment options. Default
#'   \code{NULL}.
#' @param ln_InitDevs_input Optional log-scale initial age-structure deviations
#'   array \code{[n_pop x n_regions x (n_ages - 1) x n_sims]}. The
#'   \code{n_ages - 1} dimension excludes the reference age used during
#'   initialisation. If \code{NULL} (default), deviations are initialised to
#'   zero.
#'
#' @return The input \code{sim_list} with recruitment-related fields appended:
#'   \code{$recruitment_opt}, \code{$rec_dd}, \code{$init_dd}, \code{$R0},
#'   \code{$h}, \code{$sexratio}, \code{$ln_sigmaR}, \code{$rec_seas_prop},
#'   \code{$spawn_seas}, \code{$t_spawn}, \code{$rec_lag},
#'   \code{$init_age_strc}, \code{$do_recruits_move}, \code{$move_age},
#'   \code{$stray_rate}, and optionally \code{$Rec_input} and
#'   \code{$ln_InitDevs_input}. Character-coded inputs are converted to their
#'   integer equivalents before storage.
#'
#'
#' @export Setup_Sim_Rec
#' @family Simulation Setup
Setup_Sim_Rec <- function(
    do_recruits_move = 0,
    sexratio_input = array(if(sim_list$n_sexes == 1) 1 else 0.5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_sims)),
    R0_input = array(10, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    h_input = array(0.8, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    stray_rate_input = array(0, dim = c(sim_list$n_pop, sim_list$n_yrs, sim_list$n_sims)),
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
  if(rec_lag == 0) stop("rec_lag cannot be 0!")

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
    rec_input_yrs <- dim(Rec_input)[3] # get years from Rec_input
    tmp_Rec_input <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
    # loop through simulations to resample years
    for(i in 1:sim_list$n_sims) {
      tmp_Rec_input[,,1:rec_input_yrs,i] <- Rec_input[,,,i]
      resampled_years <- sample(1:rec_input_yrs, length(tmp_Rec_input[1,1,-c(1:rec_input_yrs),i]), TRUE)
      tmp_Rec_input[,,-c(1:rec_input_yrs),i] <- Rec_input[,,resampled_years,i]
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

#' Map recruitment variability (sigma_R) parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct
#' the TMB/RTMB factor map for \code{ln_sigmaR}, the log-scale standard
#' deviation of recruitment deviations. The \code{ln_sigmaR} array has
#' dimensions \code{[2 x n_pop x n_regions]}, where the first index
#' distinguishes initial age-structure deviations (\code{i = 1}) from annual
#' recruitment deviations (\code{i = 2}).
#'
#' Sharing behaviour adapts to the density-dependence structure: under local
#' density dependence (\code{rec_dd = 0}) with a single population, each
#' region receives its own \code{ln_sigmaR} estimate; under global density
#' dependence or with multiple populations, a single \code{ln_sigmaR} is
#' shared across regions within each population.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigmaR_spec Character string specifying the estimation structure for
#'   \code{ln_sigmaR}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Estimate \code{ln_sigmaR} separately for both
#'       the initial deviation period (\code{i = 1}) and the annual deviation
#'       period (\code{i = 2}), respecting the density-dependence sharing rules
#'       described above.}
#'     \item{\code{"est_shared_all"}}{Single \code{ln_sigmaR} shared across
#'       both deviation periods, all populations, and all regions.}
#'     \item{\code{"fix_early_est_late"}}{Fix \code{ln_sigmaR} for the initial
#'       deviation period (\code{i = 1}) at its starting value; estimate
#'       \code{ln_sigmaR} for the annual deviation period (\code{i = 2}).
#'       Useful when initial age-structure uncertainty is assumed known.}
#'     \item{\code{"fix"}}{Fix all \code{ln_sigmaR} parameters at their
#'       starting values (all mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigmaR} set to a
#'   factor vector of length \code{prod(dim(par$ln_sigmaR))}. Active parameters
#'   receive sequential integer indices; fixed parameters are \code{NA}.
#'
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

#' Map initial age-structure deviation parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct
#' the TMB/RTMB factor map for \code{ln_InitDevs}, the log-scale deviations
#' from the equilibrium initial age structure. The \code{ln_InitDevs} array has
#' dimensions \code{[n_pop x n_regions x (n_ages - 1)]}, where the age
#' dimension excludes the plus group by default (see \code{equil_init_age_strc}
#' below).
#'
#' Mapping behaviour is governed by three interacting considerations:
#' \enumerate{
#'   \item \strong{Equilibrium initialisation} (\code{equil_init_age_strc}):
#'     if \code{0}, all deviations are fixed at zero (no stochastic initial
#'     structure). If \code{1}, plus-group deviations are fixed and the
#'     remaining ages are estimated or shared. If \code{2}, all ages including
#'     the plus group receive stochastic deviations.
#'   \item \strong{Sharing specification} (\code{InitDevs_spec}): controls
#'     whether deviations are shared across regions and/or populations.
#'   \item \strong{No-dispersal constraint}: when \code{rec_region_prop_spec = 1}
#'     and \code{n_pop > 1}, non-natal regions receive no recruitment and their
#'     initial age deviations are structurally zero; these are automatically
#'     fixed to \code{NA} regardless of \code{InitDevs_spec}, and the remaining
#'     indices are re-numbered sequentially.
#' }
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$equil_init_age_strc},
#'   \code{$data$rec_region_prop_spec}, \code{$data$natal_region}, and
#'   \code{$data$rec_dd} to be set by upstream setup functions.
#' @param InitDevs_spec Character string specifying the sharing structure for
#'   \code{ln_InitDevs}, or \code{NULL} to estimate all deviations independently
#'   across all dimensions. Options when non-\code{NULL}:
#'   \describe{
#'     \item{\code{"est_shared_pop_r"}}{A single set of age-specific deviations
#'       shared across all populations and regions. Each age class receives one
#'       estimated parameter regardless of how many populations or regions
#'       are modelled. Required when \code{rec_dd = "global"} and
#'       \code{n_regions > 1}.}
#'     \item{\code{"est_shared_r"}}{Separate age-specific deviations per
#'       population, shared across regions within each population. Regions
#'       within the same population are constrained to identical initial
#'       age structure. Also valid under global density dependence.}
#'     \item{\code{"fix"}}{All \code{ln_InitDevs} parameters fixed at zero
#'       (mapped to \code{NA}). Equivalent to assuming a fully deterministic
#'       initial age structure.}
#'     \item{\code{NULL}}{Estimate all deviations independently across
#'       populations, regions, and ages. Not permitted when
#'       \code{rec_region_prop_spec = 1} and \code{n_pop > 1}, as non-natal
#'       regions have no recruitment.}
#'   }
#' @param rec_dd Recruitment density-dependence structure inherited from
#'   \code{\link{Setup_Mod_Rec}}. \code{"global"} restricts valid
#'   \code{InitDevs_spec} choices to \code{"est_shared_r"} or
#'   \code{"est_shared_pop_r"} when \code{n_regions > 1}.
#'
#' @return The input \code{input_list} with \code{$map$ln_InitDevs} set to a
#'   factor vector of length \code{prod(dim(par$ln_InitDevs))}. Active
#'   parameters receive sequential integer indices; plus-group slots (when
#'   \code{equil_init_age_strc = 1}) and non-natal region slots (when
#'   \code{rec_region_prop_spec = 1}) are \code{NA}. Starting values in
#'   \code{$par$ln_InitDevs} are also reset to \code{0} for any fixed cells.
#'
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

#' Map annual recruitment deviation parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct
#' the TMB/RTMB factor map for \code{ln_RecDevs}, the log-scale annual
#' recruitment deviations. The \code{ln_RecDevs} array has dimensions
#' \code{[n_pop x n_regions x n_years]}.
#'
#' Mapping behaviour is governed by two interacting considerations:
#' \enumerate{
#'   \item \strong{Sharing specification} (\code{RecDevs_spec}): controls
#'     whether deviations are shared across regions and/or populations, or
#'     estimated independently.
#'   \item \strong{No-dispersal constraint}: when \code{rec_region_prop_spec = 1}
#'     and \code{n_pop > 1}, non-natal regions receive no recruitment and their
#'     deviations are structurally zero; these are automatically fixed to
#'     \code{NA} regardless of \code{RecDevs_spec}, and the remaining indices
#'     are re-numbered sequentially.
#' }
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$rec_region_prop_spec},
#'   \code{$data$natal_region}, \code{$data$rec_dd}, and \code{$data$n_pop}
#'   to be set by upstream setup functions.
#' @param RecDevs_spec Character string specifying the sharing structure for
#'   \code{ln_RecDevs}, or \code{NULL} to estimate all deviations independently
#'   across all dimensions. Options when non-\code{NULL}:
#'   \describe{
#'     \item{\code{"est_shared_r"}}{Separate year-specific deviations per
#'       population, shared across regions within each population. All regions
#'       of a given population follow the same annual deviation time series.
#'       Valid under both local and global density dependence.}
#'     \item{\code{"est_shared_pop_r"}}{A single set of year-specific
#'       deviations shared across all populations and regions. Each year
#'       receives one estimated parameter regardless of how many populations
#'       or regions are modelled. Required when \code{rec_dd = "global"} and
#'       \code{n_regions > 1}.}
#'     \item{\code{"fix"}}{All \code{ln_RecDevs} parameters fixed at zero
#'       (mapped to \code{NA}). Equivalent to deterministic recruitment with
#'       no interannual variability.}
#'     \item{\code{NULL}}{Estimate all deviations independently across
#'       populations, regions, and years. Not permitted when
#'       \code{rec_region_prop_spec = 1} and \code{n_pop > 1}, as non-natal
#'       regions have no recruitment.}
#'   }
#' @param rec_dd Recruitment density-dependence structure inherited from
#'   \code{\link{Setup_Mod_Rec}}. \code{"global"} restricts valid
#'   \code{RecDevs_spec} choices to \code{"est_shared_r"} or
#'   \code{"est_shared_pop_r"} when \code{n_regions > 1}.
#'
#' @return The input \code{input_list} with \code{$map$ln_RecDevs} set to a
#'   factor vector of length \code{prod(dim(par$ln_RecDevs))}. Active
#'   parameters receive sequential integer indices; non-natal region slots
#'   (when \code{rec_region_prop_spec = 1}) and fixed deviations are
#'   \code{NA}. Starting values in \code{$par$ln_RecDevs} are reset to
#'   \code{0} for any fixed cells.
#'
#' @seealso \code{\link{do_InitDevs_mapping}} for the analogous initial
#'   age-structure deviation mapping, which shares the same sharing options
#'   and no-dispersal constraint logic.
#'
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

#' Map Beverton-Holt steepness parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct
#' the TMB/RTMB factor map for \code{steepness_h}, the Beverton-Holt
#' steepness parameter. The \code{steepness_h} array has dimensions
#' \code{[n_pop x n_regions]}.
#'
#' When \code{rec_model = 0} (mean recruitment), steepness has no role in the
#' stock-recruit relationship and all elements are mapped to \code{NA}
#' regardless of \code{h_spec}. For Beverton-Holt recruitment
#' (\code{rec_model = 1}), mapping follows \code{h_spec} subject to the
#' density-dependence constraint: global density dependence requires steepness
#' to be shared across regions (\code{"est_shared_r"} or
#' \code{"est_shared_pop_r"}), since a single pooled spawner-recruit
#' relationship cannot support region-specific steepness values.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$rec_model}, \code{$data$n_pop},
#'   \code{$data$n_regions}, and \code{$data$rec_dd} to be set by upstream
#'   setup functions.
#' @param h_spec Character string specifying the sharing structure for
#'   \code{steepness_h}, or \code{NULL} to estimate steepness independently
#'   across all relevant dimensions. Options when non-\code{NULL}:
#'   \describe{
#'     \item{\code{"est_shared_pop_r"}}{Single steepness value shared across
#'       all populations and regions. All elements of \code{steepness_h} share
#'       factor level \code{1}. Required when \code{rec_dd = "global"} and
#'       \code{n_regions > 1}.}
#'     \item{\code{"est_shared_r"}}{Separate steepness per population, shared
#'       across regions within each population. Produces \code{n_pop}
#'       estimated parameters. Also valid under global density dependence.}
#'     \item{\code{"fix"}}{All \code{steepness_h} parameters fixed at their
#'       starting values (mapped to \code{NA}).}
#'     \item{\code{NULL}}{Estimate steepness independently: by population when
#'       \code{n_pop > 1} (shared across regions within each population), or
#'       by region when \code{n_pop = 1}. Not permitted when
#'       \code{rec_dd = "global"} and \code{n_regions > 1}.}
#'   }
#' @param rec_dd Recruitment density-dependence structure inherited from
#'   \code{\link{Setup_Mod_Rec}}. \code{"global"} restricts valid
#'   \code{h_spec} values to \code{"est_shared_r"}, \code{"est_shared_pop_r"},
#'   or \code{"fix"}, and prohibits \code{NULL}.
#'
#' @return The input \code{input_list} with \code{$map$steepness_h} set to a
#'   factor vector of length \code{prod(dim(par$steepness_h))}. Active
#'   parameters receive sequential integer indices; unused parameters
#'   (mean recruitment model or fixed steepness) are \code{NA}.
#'
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

#' Map sex ratio parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct
#' the TMB/RTMB factor map for \code{sexratio_pars}, the proportion of
#' recruits assigned to the first sex. The \code{sexratio_pars} array has
#' dimensions \code{[n_pop x n_regions x n_sexratio_blocks]}, where
#' \code{n_sexratio_blocks} is the maximum number of time blocks across all
#' population-region combinations as defined in \code{$data$sexratio_blocks}.
#'
#' When \code{n_sexes = 1}, estimation is meaningless and \code{sexratio_spec}
#' must be \code{"fix"}. Under the no-dispersal constraint
#' (\code{rec_region_prop_spec = 1} with \code{n_pop > 1}), \code{"est_all"}
#' is prohibited because non-natal regions receive no recruitment and cannot
#' support independent sex ratio estimates. The \code{"est_shared_pop_r"}
#' option additionally requires that all population-region combinations share
#' the same block structure; a mismatch raises an error.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_sexes}, \code{$data$n_pop},
#'   \code{$data$n_regions}, \code{$data$sexratio_blocks}, and
#'   \code{$data$rec_region_prop_spec} to be set by upstream setup functions.
#' @param sexratio_spec Character string specifying the estimation structure
#'   for \code{sexratio_pars}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate sex ratio parameter per population x
#'       region x block. Not permitted when \code{rec_region_prop_spec = 1}
#'       and \code{n_pop > 1}.}
#'     \item{\code{"est_shared_r"}}{Separate sex ratio per population x block,
#'       shared across regions within each population. Block membership is
#'       checked per region to ensure only valid blocks are assigned.}
#'     \item{\code{"est_shared_pop_r"}}{Single sex ratio per block, shared
#'       across all populations and regions. Requires identical block
#'       structures across all population-region combinations.}
#'     \item{\code{"fix"}}{All \code{sexratio_pars} fixed at their starting
#'       values (mapped to \code{NA}). Required when \code{n_sexes = 1}.}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$sexratio_pars} set to
#'   a factor vector of length \code{prod(dim(par$sexratio_pars))}. Active
#'   parameters receive sequential integer indices; fixed or invalid cells
#'   are \code{NA}.
#'
#'
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

#' Map recruitment regional apportionment parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct the
#' TMB/RTMB factor map for \code{rec_region_prop_pars}, the logit-scale
#' parameters controlling the proportion of recruits assigned to each region.
#' The array has dimensions \code{[n_pop x (n_regions - 1)]}, using a
#' sum-to-one soft-max parameterisation with one reference region omitted.
#'
#' When \code{n_regions = 1}, the parameter is structurally irrelevant and
#' both \code{$par$rec_region_prop_pars} and \code{$map$rec_region_prop_pars}
#' are set to \code{NULL}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_pop}, \code{$data$n_regions}, and
#'   \code{$data$natal_region} to be set by upstream setup functions.
#' @param rec_region_prop_spec Character string specifying the dispersal
#'   structure, or \code{NULL} to estimate all regional proportions freely.
#'   Options when non-\code{NULL}:
#'   \describe{
#'     \item{\code{"no_dispersal"}}{Recruits are assigned entirely to their
#'       natal region. Starting values for \code{rec_region_prop_pars} are
#'       overwritten with large-magnitude values (\code{-20} for non-natal
#'       regions, \code{+20} for the natal region when
#'       \code{natal_region > 1}), and all elements are mapped to \code{NA}
#'       so the parameters are not estimated. Requires \code{n_pop > 1} and
#'       \code{n_regions > 1}.}
#'     \item{\code{NULL}}{All \code{rec_region_prop_pars} are estimated
#'       independently. Only available when \code{n_regions > 1}.}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$rec_region_prop_pars}
#'   set to a factor vector of length \code{n_pop * (n_regions - 1)}, or
#'   \code{NULL} when \code{n_regions = 1}. Under \code{"no_dispersal"},
#'   \code{$par$rec_region_prop_pars} starting values are also overwritten.
#'
#'
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

#' Map stray rate parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct the
#' TMB/RTMB factor map for \code{stray_rate_pars}, the logit-scale stray rate
#' parameters. The array has dimensions \code{[n_pop x max_stray_blocks]},
#' where \code{max_stray_blocks} is the maximum number of time blocks across
#' all populations as defined by \code{$data$stray_rate_blocks}.
#'
#' When \code{n_pop = 1}, straying is not applicable and all parameters are
#' automatically fixed to \code{NA}. When \code{use_fixed_stray_rate = 1},
#' the objective function reads from \code{fixed_stray_rate} directly and
#' \code{stray_rate_pars} are not used -- all elements are fixed regardless
#' of \code{stray_rate_spec}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_pop}, \code{$data$stray_rate_blocks},
#'   and \code{$data$use_fixed_stray_rate} to be set by upstream functions.
#' @param stray_rate_spec Character string specifying the estimation structure.
#'   One of:
#'   \describe{
#'     \item{\code{"fix"}}{All parameters fixed at starting values (mapped to
#'       \code{NA}).}
#'     \item{\code{"est_all"}}{Estimate independently per population x block.
#'       Produces \code{n_pop x n_unique_blocks} estimated parameters.}
#'     \item{\code{"est_shared_pop"}}{Single parameter per block, shared across
#'       all populations. Requires identical block structures across all
#'       populations -- an error is raised if block indices differ.}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$stray_rate_pars} set to
#'   a factor vector of length \code{prod(dim(par$stray_rate_pars))}. Active
#'   parameters receive sequential integer indices; fixed parameters are
#'   \code{NA}.
#'
#' @keywords internal
do_stray_rate_mapping <- function(input_list, stray_rate_spec) {

  map_stray      <- input_list$par$stray_rate_pars
  map_stray[]    <- NA
  stray_counter  <- 1

  valid_specs <- c("fix", "est_all", "est_shared_pop")
  if (!stray_rate_spec %in% valid_specs)
    stop("Invalid stray_rate_spec. Must be one of: ", paste(valid_specs, collapse = ", "))

  # Not applicable for single population
  if (input_list$data$n_pop == 1) {
    input_list$map$stray_rate_pars <- factor(map_stray)
    collect_message("Stray rates fixed (n_pop == 1, straying not applicable).")
    return(input_list)
  }

  # Externally fixed -- pars not used by objective function
  if (input_list$data$use_fixed_stray_rate == 1) {
    input_list$map$stray_rate_pars <- factor(map_stray)
    collect_message("Stray rates are externally fixed (use_fixed_stray_rate == 1).")
    return(input_list)
  }

  if (stray_rate_spec == "fix") {
    input_list$map$stray_rate_pars <- factor(map_stray)
    collect_message("Stray rates fixed at starting values.")
    return(input_list)
  }

  # Validate block consistency for shared estimation
  if (stray_rate_spec == "est_shared_pop") {
    ref_blks <- sort(unique(as.vector(input_list$data$stray_rate_blocks[1, ])))
    for (p in seq_len(input_list$data$n_pop)) {
      if (!identical(sort(unique(as.vector(input_list$data$stray_rate_blocks[p, ]))), ref_blks))
        stop("est_shared_pop requires identical stray rate block structure across all populations.")
    }
  }

  for (p in seq_len(input_list$data$n_pop)) {

    blks <- unique(as.vector(input_list$data$stray_rate_blocks[p, ]))

    for (b in blks) {

      if (stray_rate_spec == "est_all") {
        map_stray[p, b] <- stray_counter
        stray_counter    <- stray_counter + 1
      }

      if (stray_rate_spec == "est_shared_pop" && p == 1) {
        for (pp in seq_len(input_list$data$n_pop)) map_stray[pp, b] <- stray_counter
        stray_counter <- stray_counter + 1
      }
    }
  }

  input_list$map$stray_rate_pars <- factor(map_stray)
  collect_message("Stray rates specified as: ", stray_rate_spec)
  return(input_list)
}

#' Map recruitment seasonal apportionment parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Rec}} to construct the
#' TMB/RTMB factor map for \code{rec_seas_prop_pars}, the logit-scale
#' parameters controlling the proportion of annual recruitment assigned to
#' each season. The array has dimensions \code{[n_pop x (n_seas - 1)]},
#' using a sum-to-one soft-max parameterisation with one reference season
#' omitted.
#'
#' When \code{n_seas = 1}, the parameter is structurally irrelevant and both
#' \code{$par$rec_seas_prop_pars} and \code{$map$rec_seas_prop_pars} are set
#' to \code{NULL}. If estimation is requested but
#' \code{use_fixed_rec_seas_prop = 1}, a warning is issued and
#' \code{$data$use_fixed_rec_seas_prop} is automatically reset to \code{0}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_pop}, \code{$data$n_seas}, and
#'   \code{$data$use_fixed_rec_seas_prop} to be set by upstream setup
#'   functions.
#' @param rec_seas_prop_spec Character string specifying the seasonal
#'   apportionment structure, or \code{NULL} to estimate all proportions
#'   independently across populations and seasons. Options when
#'   non-\code{NULL}:
#'   \describe{
#'     \item{\code{"est_shared_pop"}}{Estimate seasonal proportions but share
#'       them across populations, so a single set of \code{n_seas - 1}
#'       parameters applies to all populations. Only valid for seasonal models
#'       (\code{n_seas > 1}). Also resets \code{use_fixed_rec_seas_prop} to
#'       \code{0} if it was previously \code{1}.}
#'     \item{\code{"fix"}}{All \code{rec_seas_prop_pars} fixed at their
#'       starting values (mapped to \code{NA}).}
#'     \item{\code{NULL}}{Estimate all \code{n_pop x (n_seas - 1)} parameters
#'       independently. Also resets \code{use_fixed_rec_seas_prop} to
#'       \code{0} if it was previously \code{1}.}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$rec_seas_prop_pars}
#'   set to a factor vector of length \code{n_pop * (n_seas - 1)}, or
#'   \code{NULL} when \code{n_seas = 1}. \code{$data$use_fixed_rec_seas_prop}
#'   may be modified as a side effect when estimation is requested alongside
#'   a previously fixed seasonal proportion flag.
#'
#'
#' @keywords internal
do_rec_seas_prop_mapping <- function(input_list, rec_seas_prop_spec) {

  # Validate spec options
  valid_specs <- c("fix", "est_shared_pop")
  if(!is.null(rec_seas_prop_spec) && !rec_seas_prop_spec %in% valid_specs) {
    stop("Invalid rec_seas_prop_spec: '", rec_seas_prop_spec, "'. Valid options are: ", paste(valid_specs, collapse=", "), ", or NULL to estimate all.")
  }

  if((is.null(rec_seas_prop_spec) || rec_seas_prop_spec == 'est_shared_pop') && input_list$data$use_fixed_rec_seas_prop == 1) {
    input_list$data$use_fixed_rec_seas_prop <- 0
    warning("Recruitment seasonal apportionment is specified as estimated, but use_fixed_rec_seas_prop == 1 (fixed). Changing to use_fixed_rec_seas_prop == 0.")
  }

  # estimating recruitment seasonal apporitonment is only valid for seasonal models
  if(!is.null(rec_seas_prop_spec) && rec_seas_prop_spec == 'est_shared_pop' && input_list$data$n_seas == 1)
    stop("Estimating recruitment seasonal apportionment is only applicable for seasonal models. ")

  # estimate all recruitment seasonal proportions if n_seas > 1
  if(is.null(rec_seas_prop_spec)) {
    input_list$map$rec_seas_prop_pars <- factor(1:length(input_list$par$rec_seas_prop_pars))
  } else if(rec_seas_prop_spec == 'est_shared_pop') { # estimate recruitment seasonal proportions but share across populations
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

#' Set up the recruitment module and associated processes
#'
#' Configures all recruitment-related components of the estimation model:
#' stock-recruit relationship type and density-dependence structure,
#' Beverton-Holt steepness and priors, recruitment variability
#' (\eqn{\sigma_R}), annual and initial age-structure deviations, regional
#' and seasonal recruitment apportionment, spawning movement and stray rates,
#' sex ratio dynamics, equilibrium initialisation method, and the recruitment
#' bias ramp. Delegates parameter mapping to a family of internal helpers
#' (\code{\link{do_sigmaR_mapping}}, \code{\link{do_RecDevs_mapping}},
#' \code{\link{do_InitDevs_mapping}}, \code{\link{do_h_mapping}},
#' \code{\link{do_sexratio_pars_mapping}},
#' \code{\link{do_rec_region_prop_mapping}},
#' \code{\link{do_rec_seas_prop_mapping}}). Must be called after
#' \code{\link{Setup_Mod_Dim}} and \code{\link{Setup_Mod_Biologicals}}, and
#' before model compilation.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#'   Population, region, age, year, and season dimensions must already be
#'   defined in \code{$data}.
#'
#' @param rec_model Character string (required). Stock-recruit relationship:
#'   \describe{
#'     \item{\code{"mean_rec"}}{Fixed mean recruitment; no stock-recruit
#'       relationship. Steepness is automatically fixed and not estimated.}
#'     \item{\code{"bh_rec"}}{Beverton-Holt stock-recruit relationship.}
#'   }
#' @param rec_dd Density-dependence structure. Default \code{"global"}.
#'   \describe{
#'     \item{\code{"local"}}{Independent stock-recruit relationship per
#'       population. Required when \code{n_pop > 1}.}
#'     \item{\code{"global"}}{Single pooled spawner-recruit relationship
#'       across all regions. Constrains \code{h_spec},
#'       \code{RecDevs_spec}, and \code{InitDevs_spec} to shared or fixed
#'       options when \code{n_regions > 1}.}
#'   }
#' @param rec_lag Integer. Lag between spawning biomass and age-1
#'   recruitment (in seasons). Must be \eqn{\geq 1}. Default \code{1}.
#'
#' @param sigmaR_spec Character. Estimation structure for \eqn{\sigma_R},
#'   stored in \code{ln_sigmaR} \code{[2 x n_pop x n_regions]}, where index
#'   1 = initial deviation period and index 2 = annual deviation period.
#'   Default \code{"est_all"}. See \code{\link{do_sigmaR_mapping}} for full
#'   option descriptions.
#' @param sigmaR_switch Integer. Year index at which \eqn{\sigma_R} switches
#'   from the early-period value (index 1) to the late-period value (index 2).
#'   If \eqn{\leq 1}, a single \eqn{\sigma_R} is applied throughout. Default
#'   \code{1}.
#' @param RecDevs_spec Character or \code{NULL}. Sharing structure for annual
#'   recruitment deviations \code{ln_RecDevs} \code{[n_pop x n_regions x
#'   n_years]}. Default \code{NULL} (estimate all independently). See
#'   \code{\link{do_RecDevs_mapping}} for full option descriptions.
#' @param dont_est_recdev_last Non-negative integer. Number of terminal years
#'   for which recruitment deviations are not estimated. Automatically
#'   overridden to \code{0} if \code{n_proj_yrs_devs > 0}, since projected
#'   deviation years are penalised toward the mean and are effectively
#'   estimated regardless. Default \code{0}.
#'
#' @param init_age_strc Equilibrium initialisation method. Default \code{2}.
#'   \describe{
#'     \item{\code{0}/\code{"iterative"}}{Iterates the population to
#'       approximate equilibrium. Slowest but most general.}
#'     \item{\code{1}/\code{"scalar_no_move"}}{Scalar geometric series
#'       assuming no movement.}
#'     \item{\code{2}/\code{"matrix"}}{Matrix geometric series incorporating
#'       movement. Recommended default for spatial models.}
#'     \item{\code{3}/\code{"scalar_plus_only"}}{Scalar geometric series
#'       with movement only in the plus group.}
#'   }
#' @param equil_init_age_strc Plus-group treatment during stochastic
#'   initialisation. Default \code{1}.
#'   \describe{
#'     \item{\code{0}/\code{"equil"}}{Deterministic equilibrium; no
#'       \code{ln_InitDevs} are estimated.}
#'     \item{\code{1}/\code{"stoch_no_plus"}}{Stochastic deviations for all
#'       ages except the plus group.}
#'     \item{\code{2}/\code{"stoch_all"}}{Stochastic deviations for all ages
#'       including the plus group.}
#'   }
#' @param InitDevs_spec Character or \code{NULL}. Sharing structure for
#'   initial age-structure deviations \code{ln_InitDevs} \code{[n_pop x
#'   n_regions x (n_ages - 1)]}. Default \code{NULL} (estimate all
#'   independently). See \code{\link{do_InitDevs_mapping}} for full option
#'   descriptions.
#' @param init_F_prop Numeric vector \code{[n_seas]}. Seasonal distribution of
#'   fishing mortality applied during equilibrium initialisation. Default:
#'   zero for all seasons.
#'
#' @param rec_region_prop_spec Character or \code{NULL}. Regional recruitment
#'   dispersal structure. Default \code{NULL} (estimate all proportions
#'   freely). See \code{\link{do_rec_region_prop_mapping}} for full option
#'   descriptions including \code{"no_dispersal"}. Stored as
#'   \code{$data$rec_region_prop_spec}: \code{0} = full dispersal,
#'   \code{1} = no dispersal.
#' @param use_rec_region_prop_prior Integer (0/1). Whether Dirichlet priors
#'   are applied to regional recruitment proportions. Not valid when
#'   \code{n_regions = 1}. Default \code{0}.
#' @param rec_region_prop_prior Data frame of Dirichlet prior concentration
#'   parameters. Required columns: \code{pop} and \code{alpha}, where
#'   \code{alpha} is a list-column of length-\code{n_regions} vectors.
#'   Ignored when \code{use_rec_region_prop_prior = 0}. Default \code{NULL}.
#'
#'
#' @param rec_seas_prop_spec Character or \code{NULL}. Seasonal recruitment
#'   apportionment structure. Default \code{"fix"}. See
#'   \code{\link{do_rec_seas_prop_mapping}} for full option descriptions
#'   including \code{"est_shared_pop"}.
#' @param use_fixed_rec_seas_prop Integer (0/1). Whether fixed (non-estimated)
#'   seasonal proportions from \code{fixed_rec_seas_prop} are used. Automatically
#'   reset to \code{0} with a warning if \code{rec_seas_prop_spec} requests
#'   estimation. Default \code{1}.
#' @param fixed_rec_seas_prop Array \code{[n_pop x n_seas]}. Fixed seasonal
#'   recruitment proportions used when \code{use_fixed_rec_seas_prop = 1}.
#'   Default: all recruitment assigned to season 1.
#' @param use_rec_seas_prop_prior Integer (0/1). Whether Dirichlet priors are
#'   applied to seasonal recruitment proportions. Not valid when
#'   \code{n_seas = 1}. Default \code{0}.
#' @param rec_seas_prop_prior Data frame of Dirichlet prior concentration
#'   parameters for seasonal proportions. Required columns: \code{pop} and
#'   \code{alpha}. Ignored when \code{use_rec_seas_prop_prior = 0}. Default
#'   \code{NULL}.
#'
#' @param h_spec Character or \code{NULL}. Sharing structure for
#'   Beverton-Holt steepness \code{steepness_h} \code{[n_pop x n_regions]},
#'   parameterised in bounded logit space \eqn{(0.2, 1)}. Default \code{NULL}
#'   (estimate by population when \code{n_pop > 1}, by region when
#'   \code{n_pop = 1}). Ignored when \code{rec_model = "mean_rec"}. See
#'   \code{\link{do_h_mapping}} for full option descriptions.
#' @param Use_h_prior Integer (0/1). Whether normal priors on steepness are
#'   applied. Only relevant for \code{rec_model = "bh_rec"}. Default \code{0}.
#' @param h_prior Data frame of steepness prior parameters. Required columns:
#'   \code{pop}, \code{region}, \code{mu}, \code{sd}. Ignored when
#'   \code{Use_h_prior = 0}. Default \code{NULL}.
#'
#' @param spawn_seas Integer. Season index in which spawning occurs. Default
#'   \code{1}.
#' @param t_spawn Numeric. Spawn timing as a fraction of the season elapsed
#'   before spawning. \code{0} (default) = spawning before any mortality;
#'   \code{1} = spawning after all mortality.
#' @param sgl_seas_spawning_movement Spawning movement array
#'   \code{[n_pop x n_regions x n_regions x n_years x n_ages x n_sexes]}.
#'   Each \code{[p, , r, y, a, s]} slice is a row-stochastic movement matrix
#'   giving the probability of fish from each origin region spawning in region
#'   \code{r}. If \code{NA} (default), 100\% natal homing is assumed and the
#'   array is constructed internally.
#'
#' @param use_fixed_stray_rate Integer (0/1). Whether stray rates are supplied
#'   as a fixed external array (\code{fixed_stray_rate}) rather than estimated
#'   as model parameters. Default \code{1} (fixed), preserving existing
#'   behaviour. Set to \code{0} to estimate stray rates via
#'   \code{stray_rate_pars}.
#' @param fixed_stray_rate Array \code{[n_pop x n_years]}. Fixed stray rate
#'   values used when \code{use_fixed_stray_rate = 1}. Values should be in
#'   \eqn{[0, 1]}. Default: \code{0} (no straying) for all populations and
#'   years. Ignored when \code{use_fixed_stray_rate = 0}.
#' @param stray_rate_spec Character string. Estimation structure for
#'   \code{stray_rate_pars} \code{[n_pop x max_stray_blocks]}, parameterised
#'   on the logit scale. Ignored when \code{use_fixed_stray_rate = 1} or
#'   \code{n_pop = 1}. Default \code{"fix"}. Options:
#'   \describe{
#'     \item{\code{"fix"}}{All parameters fixed at starting values (mapped to
#'       \code{NA}). Use this alongside \code{use_fixed_stray_rate = 0} to
#'       hold stray rates at a specified value without estimating.}
#'     \item{\code{"est_all"}}{Estimate independently per population x block.
#'       Produces one parameter per population per unique block.}
#'     \item{\code{"est_shared_pop"}}{Single parameter per block, shared across
#'       all populations. Requires identical block structures across all
#'       populations -- an error is raised if block indices differ.}
#'   }
#' @param stray_rate_blocks Character vector of length \code{n_pop} defining
#'   the temporal block structure for stray rate parameters. Valid formats:
#'   \describe{
#'     \item{\code{"none_Pop_x"}}{Constant stray rate for population \code{x}
#'       across all years (single block).}
#'     \item{\code{"Block_k_Year_a-b_Pop_x"}}{Block \code{k} applies to years
#'       \code{a} through \code{b} for population \code{x}. Use
#'       \code{"terminal"} in place of the end year to extend through the
#'       final model year.}
#'   }
#'   Default: a single constant block for every population.
#'   \strong{Note:} stray rate is generally unidentifiable from fisheries data
#'   alone. Time-blocking is provided for completeness but regularisation via
#'   \code{use_stray_rate_prior} in the penalty setup is strongly recommended
#'   whenever \code{stray_rate_spec != "fix"}.
#' @param use_stray_rate_prior Integer (0/1). Whether Beta priors are applied
#'   to estimated stray rate parameters. Only relevant when
#'   \code{use_fixed_stray_rate = 0} and \code{n_pop > 1}. An error is raised
#'   if \code{use_stray_rate_prior = 1} alongside \code{use_fixed_stray_rate = 1}
#'   since \code{stray_rate_pars} would not be estimated. Default \code{0}.
#' @param stray_rate_prior Data frame of Beta prior parameters for stray rates.
#'   Required columns: \code{pop} (population index), \code{block} (block
#'   index matching \code{stray_rate_blocks}), \code{mu} (prior mean, in
#'   \eqn{(0,1)}), \code{sd} (prior standard deviation). One row per
#'   population x block combination to penalise. Ignored when
#'   \code{use_stray_rate_prior = 0}. Default \code{NULL}.
#'
#' @param sexratio_spec Character. Estimation structure for sex ratio
#'   parameters \code{sexratio_pars} \code{[n_pop x n_regions x n_blocks]}.
#'   Default \code{"fix"}. See \code{\link{do_sexratio_pars_mapping}} for
#'   full option descriptions. Must be \code{"fix"} when \code{n_sexes = 1}.
#' @param sexratio_blocks Character vector defining temporal block structure
#'   for sex ratio parameters. One entry per population-region combination.
#'   Valid formats:
#'   \describe{
#'     \item{\code{"none_Pop_x_Region_x"}}{Constant sex ratio for population
#'       \code{x} and region \code{x} (single block across all years).}
#'     \item{\code{"Block_k_Year_a-b_Pop_x_Region_x"}}{Block \code{k}
#'       applies to years \code{a} through \code{b}. Use \code{"terminal"}
#'       in place of the end year to extend through the final model year.}
#'   }
#'   Default: a single constant block for every population-region combination.
#' @param do_rec_bias_ramp Integer (0/1). Whether a recruitment bias ramp is
#'   applied to \code{ln_RecDevs} to account for reduced information in early
#'   and terminal years. Default \code{0}.
#' @param bias_year Numeric. Calendar year at which the bias ramp reaches its
#'   maximum correction. Only used when \code{do_rec_bias_ramp = 1}. Default
#'   \code{NA}.
#' @param max_bias_ramp_fct Numeric in \eqn{[0, 1]}. Maximum bias correction
#'   factor applied at \code{bias_year}. Default \code{1}.
#'
#' @param ... Optional named starting values for parameters. Any of:
#'   \code{ln_global_R0} \code{[n_pop]},
#'   \code{rec_region_prop_pars} \code{[n_pop x (n_regions - 1)]},
#'   \code{rec_seas_prop_pars} \code{[n_pop x (n_seas - 1)]},
#'   \code{steepness_h} \code{[n_pop x n_regions]} (bounded logit scale),
#'   \code{ln_InitDevs} \code{[n_pop x n_regions x (n_ages - 1)]},
#'   \code{ln_RecDevs} \code{[n_pop x n_regions x n_years]},
#'   \code{ln_sigmaR} \code{[2 x n_pop x n_regions]},
#'   \code{sexratio_pars} \code{[n_pop x n_regions x n_blocks]}.
#'   Unspecified parameters use internal defaults.
#'
#' @return The input \code{input_list} with all recruitment-related fields
#'   populated in \code{$data} and \code{$par}, and factor maps constructed
#'   in \code{$map} for: \code{rec_region_prop_pars}, \code{rec_seas_prop_pars},
#'   \code{ln_sigmaR}, \code{ln_InitDevs}, \code{ln_RecDevs},
#'   \code{steepness_h}, \code{sexratio_pars}, and \code{stray_rate_pars}. Character-coded inputs
#'   for \code{init_age_strc} and \code{equil_init_age_strc} are converted to
#'   integer codes before storage.
#'
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
                          fixed_rec_seas_prop = {
                            rec_seas_prop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_seas))
                            rec_seas_prop[, 1] <- 1
                            rec_seas_prop
                          },
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
                          stray_rate_spec = "fix",
                          stray_rate_blocks = paste0("none_Pop_", seq_len(input_list$data$n_pop)),
                          use_fixed_stray_rate = if(stray_rate_spec != 'fix') 0 else 1,
                          fixed_stray_rate = array(0, dim = c(input_list$data$n_pop, length(input_list$data$years))),
                          use_stray_rate_prior = 0,
                          stray_rate_prior = NULL,
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
  if(input_list$store_config) input_list$config$Setup_Mod_Rec <- mget(names(formals()))[-1]

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
  if(rec_model != "mean_rec") collect_message("Recruitment and SSB lag is specified as: ", rec_lag)
  if(rec_lag == 0) stop("rec_lag cannot be 0!")

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
  stray_rate_blocks_mat <- array(NA, dim = c(input_list$data$n_pop,
                                             length(input_list$data$years)))

  for (i in seq_along(stray_rate_blocks)) {

    # parse
    tmp     <- stray_rate_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    if (!tmp_vec[1] %in% c("none", "Block"))
      stop("stray_rate_blocks not correctly specified. ",
           "Use 'none_Pop_x' or 'Block_k_Year_a-b_Pop_x'.")

    # if none
    if (tmp_vec[1] == "none") {
      pop <- as.numeric(tmp_vec[3])
      stray_rate_blocks_mat[pop, ] <- 1
    }

    # if blocks
    if (tmp_vec[1] == "Block") {
      block_val <- as.numeric(tmp_vec[2])
      pop        <- as.numeric(tmp_vec[6])
      if (!str_detect(tmp, "terminal")) {
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        yrs        <- year_range[1]:year_range[2]
      } else {
        yrs <- as.numeric(unlist(strsplit(tmp_vec[4], "-"))[1]):length(input_list$data$years)
      }
      stray_rate_blocks_mat[pop, yrs] <- block_val
    }
  }

  for (p in seq_len(input_list$data$n_pop))
    collect_message("Stray rates for population ", p, " specified with ", length(unique(stray_rate_blocks_mat[p, ])), " block(s).")

  # validate priors
  if (!use_stray_rate_prior %in% c(0, 1)) stop("use_stray_rate_prior must be 0 or 1")
  if (use_stray_rate_prior == 1 && input_list$data$n_pop == 1)
    stop("Stray rate priors are not applicable when n_pop == 1.")
  if (use_stray_rate_prior == 1 && use_fixed_stray_rate == 1)
    stop("use_stray_rate_prior == 1 but use_fixed_stray_rate == 1 - stray_rate_pars are not estimated so a prior has no effect.")
  if (use_stray_rate_prior == 1) {
    required_cols <- c("pop", "block", "mu", "sd")
    missing_cols  <- setdiff(required_cols, names(stray_rate_prior))
    if (length(missing_cols) > 0)
      stop("stray_rate_prior is missing columns: ", paste(missing_cols, collapse = ", "))
    if (any(stray_rate_prior$mu <= 0 | stray_rate_prior$mu >= 1))
      stop("stray_rate_prior$mu must be in (0, 1).")
  }
  collect_message("Stray rate prior is: ", ifelse(use_stray_rate_prior == 1, "Used", "Not Used"))

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
  input_list$data$use_fixed_stray_rate <- use_fixed_stray_rate
  input_list$data$fixed_stray_rate     <- fixed_stray_rate
  input_list$data$stray_rate_blocks    <- stray_rate_blocks_mat
  input_list$data$sexratio_blocks <- sexratio_blocks_mat
  input_list$data$use_fixed_rec_seas_prop <- use_fixed_rec_seas_prop
  input_list$data$fixed_rec_seas_prop <- fixed_rec_seas_prop
  input_list$data$use_rec_seas_prop_prior <- use_rec_seas_prop_prior
  input_list$data$rec_seas_prop_prior <- rec_seas_prop_prior
  input_list$data$use_stray_rate_prior <- use_stray_rate_prior
  input_list$data$stray_rate_prior     <- stray_rate_prior

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

  # Stray rate parameters (logit scale)
  max_stray_blks <- if (input_list$data$n_pop > 1) max(apply(stray_rate_blocks_mat, 1, function(x) length(unique(x))))
  else 1

  max_stray_blks <- if (input_list$data$n_pop > 1)  max(apply(stray_rate_blocks_mat, 1, function(x) length(unique(x)))) else 1
  if ("stray_rate_pars" %in% names(starting_values)) input_list$par$stray_rate_pars <- starting_values$stray_rate_pars
  else input_list$par$stray_rate_pars <- array(0,  dim = c(input_list$data$n_pop, max_stray_blks))

  # Mapping Options -----------------------------------------------------------

  input_list <- do_rec_region_prop_mapping(input_list, rec_region_prop_spec) # Recruitment regional proportion mapping
  input_list <- do_rec_seas_prop_mapping(input_list, rec_seas_prop_spec) # Recruitment seasonal proportion mapping
  input_list <- do_sigmaR_mapping(input_list, sigmaR_spec) # sigmaR mapping
  input_list <- do_InitDevs_mapping(input_list, InitDevs_spec, rec_dd) # InitDevs mapping
  input_list <- do_RecDevs_mapping(input_list, RecDevs_spec, rec_dd) # RevDevs mapping
  input_list <- do_h_mapping(input_list, h_spec, rec_dd) # steepness mapping
  input_list <- do_sexratio_pars_mapping(input_list, sexratio_spec) # sex ratio parameters
  input_list <- do_stray_rate_mapping(input_list, stray_rate_spec) # stray rates

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}
