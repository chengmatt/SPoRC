#' Initialise output containers for the operating model simulation
#'
#' Allocates and appends zero-filled arrays to \code{sim_list} for all biological,
#' fishery, and survey quantities tracked during a closed-loop simulation. This
#' function should be called after \code{\link{Setup_Sim_Dim}} and before any
#' operating model dynamics are run. All arrays are pre-filled with zeros and
#' populated in subsequent simulation steps.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}},
#'   containing the following named dimension scalars used to size the containers:
#'   \code{n_pop}, \code{n_regions}, \code{n_yrs}, \code{n_seas}, \code{n_ages},
#'   \code{n_sexes}, \code{n_sims}, \code{n_fish_fleets}, \code{n_srv_fleets},
#'   \code{n_obs_ages}, \code{n_lens}.
#'
#' @return The input \code{sim_list} with the following zero-filled arrays appended:
#'
#'   **Biological containers**
#'   \describe{
#'     \item{\code{$NAA}}{Numbers-at-age \code{[n_pop × n_regions × (n_yrs+1) × n_seas × n_ages × n_sexes × n_sims]}.
#'       The \code{+1} year dimension accommodates the plus-group projection year.}
#'     \item{\code{$NAA_bef}, \code{$NAA_aft}}{Numbers-at-age immediately before and
#'       after fishing mortality is applied, same dimensions as \code{$NAA}.}
#'     \item{\code{$NAA0}}{Unfished numbers-at-age, same dimensions as \code{$NAA}.
#'       Used to compute dynamic \eqn{B_0} reference quantities.}
#'     \item{\code{$ZAA}}{Total instantaneous mortality-at-age
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]}.}
#'     \item{\code{$Rec}}{Recruitment \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$SSB}}{Spawning stock biomass \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$eff_SSB}}{Effective (population-aggregated) spawning stock biomass
#'       \code{[n_pop × n_yrs × n_sims]}, summed across regions for stock-recruit calculations.}
#'     \item{\code{$Dynamic_SSB0}}{Dynamic unfished spawning stock biomass
#'       \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$Total_Biom}}{Total biomass \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$ln_RecDevs}}{Log-scale recruitment deviations
#'       \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$ln_InitDevs}}{Log-scale initial age-structure deviations
#'       \code{[n_pop × n_regions × (n_ages - 1) × n_sims]}.
#'       The minus-one accounts for the reference age fixed during initialisation.}
#'   }
#'
#'   **Fishery containers**
#'   \describe{
#'     \item{\code{$ObsCatch}}{Observed (noisy) catch in weight
#'       \code{[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]}.}
#'     \item{\code{$TrueCatch}}{True (error-free) catch in weight, same dimensions as
#'       \code{$ObsCatch}.}
#'     \item{\code{$ObsFishIdx}}{Observed fishery CPUE index, same dimensions as
#'       \code{$ObsCatch}.}
#'     \item{\code{$TrueFishIdx}}{True fishery index, same dimensions as \code{$ObsCatch}.}
#'     \item{\code{$CAA}}{Catch-at-age (true)
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims]}.}
#'     \item{\code{$CAL}}{Catch-at-length (true)
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]}.}
#'     \item{\code{$ObsFishAgeComps}}{Observed fishery age compositions
#'       \code{[n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims]}.
#'       Uses \code{n_obs_ages} rather than \code{n_ages} to reflect the observed (binned) age structure.}
#'     \item{\code{$ObsFishLenComps}}{Observed fishery length compositions
#'       \code{[n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]}.}
#'   }
#'
#'   **Survey containers**
#'   \describe{
#'     \item{\code{$ObsSrvIdx}}{Observed survey abundance or biomass index
#'       \code{[n_regions × n_yrs × n_seas × n_srv_fleets × n_sims]}.}
#'     \item{\code{$TrueSrvIdx}}{True survey index, same dimensions as \code{$ObsSrvIdx}.}
#'     \item{\code{$SrvIAA}}{Survey index-at-age (true)
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]}.}
#'     \item{\code{$SrvIAL}}{Survey index-at-length (true), same dimensions as
#'       \code{$SrvIAA} but with \code{n_lens} replacing \code{n_ages}.}
#'     \item{\code{$ObsSrvAgeComps}}{Observed survey age compositions
#'       \code{[n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_srv_fleets × n_sims]}.}
#'     \item{\code{$ObsSrvLenComps}}{Observed survey length compositions
#'       \code{[n_regions × n_yrs × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims]}.}
#'   }
#'
#' @export Setup_Sim_Containers
#' @family Simulation Setup
Setup_Sim_Containers <- function(sim_list) {

  # Biological Containers
  sim_list$NAA_aft = sim_list$NAA_bef = sim_list$NAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs+1, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # numbers at age
  sim_list$NAA0 <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs+1,  sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # unfished numbers at age
  sim_list$ZAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # total mortality at age
  sim_list$Rec <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # recruitment
  sim_list$SSB <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # spawning biomass
  sim_list$eff_SSB <- array(0, dim = c(sim_list$n_pop, sim_list$n_yrs, sim_list$n_sims)) # effective spawning biomass by population
  sim_list$Dynamic_SSB0 <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,  sim_list$n_sims)) # Dynamic unfished spawning biomass
  sim_list$Total_Biom <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # total biomass
  sim_list$ln_RecDevs <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # rec devs
  sim_list$ln_InitDevs <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_ages - 1, sim_list$n_sims)) # initial age devs

  # Fishery Containers
  sim_list$ObsCatch <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # observed catch
  sim_list$ObsFishAgeComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # observed fishery age comps
  sim_list$ObsFishLenComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # observed fishery length comps
  sim_list$CAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # catch at age
  sim_list$CAL <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # catch at length
  sim_list$TrueCatch <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # true catches
  sim_list$ObsFishIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # observed fishery index
  sim_list$TrueFishIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # true fishery index

  # Survey Containers
  sim_list$ObsSrvIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims)) # observed survey index
  sim_list$TrueSrvIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims)) # true survey index
  sim_list$ObsSrvAgeComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # observed survey age comps
  sim_list$ObsSrvLenComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # observed survey length comps
  sim_list$SrvIAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # survey index at age
  sim_list$SrvIAL <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # survey index at length

  return(sim_list)
}
