# Stage 1 of 3: model setup
#
# Allocates the arrays the operating model fills in. Called early in the
# Setup_Sim_* chain so the later stages have somewhere to write.

#' Initialise output containers for the operating model simulation
#'
#' Allocates and appends zero-initialised arrays to \code{sim_list} for all
#' biological, fishery, and survey quantities tracked during simulation. This
#' function should be called after \code{\link{Setup_Sim_Dim}} and before any
#' operating model dynamics are run. All arrays are pre-allocated and populated
#' in subsequent simulation steps.
#'
#' @param sim_list A simulation list returned by \code{\link{Setup_Sim_Dim}}.
#'   Dimension elements (e.g., \code{n_pop}, \code{n_regions}, \code{n_yrs},
#'   \code{n_seas}, \code{n_ages}, \code{n_sexes}, \code{n_sims},
#'   \code{n_fish_fleets}, \code{n_srv_fleets}, \code{n_obs_ages},
#'   \code{n_lens}) are used to size all output containers.
#'
#'
#' @return The input \code{sim_list} with the following zero-initialised arrays added:
#'
#'   **Biological containers**
#'   \describe{
#'     \item{\code{$NAA}}{Numbers-at-age
#'       \code{[n_pop × n_regions × (n_yrs+1) × n_seas × n_ages × n_sexes × n_sims]}.
#'       The \code{+1} year dimension stores initial conditions and propagates
#'       population state through the final year.}
#'     \item{\code{$NAA_bef}, \code{$NAA_aft}}{Numbers-at-age immediately before and
#'       after fishing mortality is applied; same dimensions as \code{$NAA}.}
#'     \item{\code{$NAA0}}{Unfished numbers-at-age; same dimensions as \code{$NAA}.
#'       Used to compute dynamic \eqn{B_0} reference quantities.}
#'     \item{\code{$ZAA}}{Total instantaneous mortality-at-age
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]}.}
#'     \item{\code{$Rec}}{Recruitment
#'       \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$SSB}}{Spawning stock biomass
#'       \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$eff_SSB}}{Effective (population-aggregated) spawning stock biomass
#'       \code{[n_pop × n_yrs × n_sims]}.}
#'     \item{\code{$Dynamic_SSB0}}{Dynamic unfished spawning stock biomass
#'       \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$Total_Biom}}{Total biomass
#'       \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$ln_RecDevs}}{Log-scale recruitment deviations
#'       \code{[n_pop × n_regions × n_yrs × n_sims]}.}
#'     \item{\code{$ln_InitDevs}}{Log-scale initial age-structure deviations
#'       \code{[n_pop × n_regions × (n_ages - 1) × n_sexes × n_sims]}.}
#'   }
#'
#'   **Fishery containers**
#'
#'   *Retained catch and indices (aggregated)*
#'   \describe{
#'     \item{\code{$ObsCatch}, \code{$TrueCatch}}{Observed and true catch
#'       \code{[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]}.}
#'     \item{\code{$ObsFishIdx}, \code{$TrueFishIdx}}{Observed and true fishery CPUE index;
#'       same dimensions as \code{$ObsCatch}.}
#'     \item{\code{$ObsFishAgeComps}}{Observed fishery age compositions
#'       \code{[n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims]}.}
#'     \item{\code{$ObsFishLenComps}}{Observed fishery length compositions
#'       \code{[n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]}.}
#'   }
#'
#'   *Discards (aggregated)*
#'   \describe{
#'     \item{\code{$ObsDiscard}, \code{$TrueDiscard}}{Observed and true discard
#'       \code{[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]}.}
#'     \item{\code{$ObsFishAgeComps_discard}}{Observed discard age compositions;
#'       same dimensions as \code{$ObsFishAgeComps}.}
#'     \item{\code{$ObsFishLenComps_discard}}{Observed discard length compositions;
#'       same dimensions as \code{$ObsFishLenComps}.}
#'   }
#'
#'   *Population-specific quantities*
#'   \describe{
#'     \item{\code{$ObsCatch_pop}, \code{$TrueCatch_pop}}{Observed and true catch
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]}.}
#'     \item{\code{$ObsFishIdx_pop}, \code{$TrueFishIdx_pop}}{Observed and true fishery index;
#'       same dimensions as \code{$ObsCatch_pop}.}
#'     \item{\code{$ObsFishAgeComps_pop}}{Observed fishery age compositions
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims]}.}
#'     \item{\code{$ObsFishLenComps_pop}}{Observed fishery length compositions
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]}.}
#'   }
#'
#'   *Discards (population-specific)*
#'   \describe{
#'     \item{\code{$ObsDiscard_pop}, \code{$TrueDiscard_pop}}{Observed and true discard;
#'       same dimensions as \code{$ObsCatch_pop}.}
#'     \item{\code{$ObsFishAgeComps_discard_pop}}{Observed discard age compositions;
#'       same dimensions as \code{$ObsFishAgeComps_pop}.}
#'     \item{\code{$ObsFishLenComps_discard_pop}}{Observed discard length compositions;
#'       same dimensions as \code{$ObsFishLenComps_pop}.}
#'   }
#'
#'   *True catch/discard at age and length*
#'   \describe{
#'     \item{\code{$CAA}, \code{$DAA}}{Catch- and discard-at-age (true)
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims]}.}
#'     \item{\code{$CAL}, \code{$DAL}}{Catch- and discard-at-length (true)
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]}.}
#'   }
#'
#'   **Survey containers**
#'
#'   *Aggregated*
#'   \describe{
#'     \item{\code{$ObsSrvIdx}, \code{$TrueSrvIdx}}{Observed and true survey index
#'       \code{[n_regions × n_yrs × n_seas × n_srv_fleets × n_sims]}.}
#'     \item{\code{$ObsSrvAgeComps}}{Observed survey age compositions
#'       \code{[n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_srv_fleets × n_sims]}.}
#'     \item{\code{$ObsSrvLenComps}}{Observed survey length compositions
#'       \code{[n_regions × n_yrs × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims]}.}
#'   }
#'
#'   *Population-specific*
#'   \describe{
#'     \item{\code{$ObsSrvIdx_pop}, \code{$TrueSrvIdx_pop}}{Observed and true survey index
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_srv_fleets × n_sims]}.}
#'     \item{\code{$ObsSrvAgeComps_pop}}{Observed survey age compositions
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_srv_fleets × n_sims]}.}
#'     \item{\code{$ObsSrvLenComps_pop}}{Observed survey length compositions
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims]}.}
#'   }
#'
#'   *True index at age and length*
#'   \describe{
#'     \item{\code{$SrvIAA}}{Survey index-at-age (true)
#'       \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]}.}
#'     \item{\code{$SrvIAL}}{Survey index-at-length (true), same structure as
#'       \code{$SrvIAA} with \code{n_lens} replacing \code{n_ages}.}
#'   }
#'
#' @export Setup_Sim_Containers
#' @family Simulation Setup
Setup_Sim_Containers <- function(sim_list) {

  # Biological Containers
  sim_list$NAA_aft = sim_list$NAA_bef = sim_list$NAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs+1, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  sim_list$NAA0 <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs+1, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  sim_list$ZAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  sim_list$Rec <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
  sim_list$SSB <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
  sim_list$eff_SSB <- array(0, dim = c(sim_list$n_pop, sim_list$n_yrs, sim_list$n_sims))
  sim_list$Dynamic_SSB0 <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
  sim_list$Total_Biom <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
  sim_list$ln_RecDevs <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
  sim_list$ln_InitDevs <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_ages - 1, sim_list$n_sexes, sim_list$n_sims))

  # Fishery Containers
  # Aggregated (pooled across populations)
  sim_list$ObsCatch <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$TrueCatch <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$TrueFishIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishAgeComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishLenComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  # Conditional age-at-length, one age composition per length bin
  sim_list$ObsFish_caal <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsDiscard <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$TrueDiscard <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishAgeComps_discard <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishLenComps_discard <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))

  # Population-specific
  sim_list$ObsCatch_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$TrueCatch_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishIdx_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$TrueFishIdx_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishAgeComps_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishLenComps_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsDiscard_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$TrueDiscard_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishAgeComps_discard_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$ObsFishLenComps_discard_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))

  # True catch-at-age/length (always pop-resolved)
  sim_list$CAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$CAL <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$DAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$DAL <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))

  # Survey Containers
  # Aggregated (pooled across populations)
  sim_list$ObsSrvIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$TrueSrvIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$ObsSrvAgeComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$ObsSrvLenComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$ObsSrv_caal <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))

  # Population-specific
  sim_list$ObsSrvIdx_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$TrueSrvIdx_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$ObsSrvAgeComps_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$ObsSrvLenComps_pop <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))

  # True index-at-age/length (always pop-resolved)
  sim_list$SrvIAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$SrvIAL <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))

  return(sim_list)
}
