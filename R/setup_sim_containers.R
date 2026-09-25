# Stage 1 of 3: model setup
#
# Allocates the arrays the operating model fills in. Called early in the Setup_Sim_* chain so the
# later stages have somewhere to write.

#' Initialize output containers for the operating model simulation
#'
#' Allocates the zero-filled arrays the simulation writes into, sized off the
#' dimensions in \code{sim_list}. Call after \code{\link{Setup_Sim_Dim}} and before
#' any operating model dynamics.
#'
#' @param sim_list A simulation list returned by \code{\link{Setup_Sim_Dim}}, whose
#'   \code{n_pop}, \code{n_regions}, \code{n_yrs}, \code{n_seas}, \code{n_ages},
#'   \code{n_sexes}, \code{n_sims}, \code{n_fish_fleets}, \code{n_srv_fleets},
#'   \code{n_obs_ages} and \code{n_lens} size every container.
#'
#' @return \code{sim_list} with the containers added.
#'
#'   Biological: \code{$NAA} \code{[n_pop × n_regions × (n_yrs+1) × n_seas × n_ages
#'   × n_sexes × n_sims]}, whose extra year holds the initial conditions and
#'   advances the population through the final year; \code{$NAA_bef} and
#'   \code{$NAA_aft}, the numbers before and after fishing mortality; \code{$NAA0},
#'   the unfished numbers, for dynamic \eqn{B_0}; \code{$ZAA}, total mortality at
#'   age over \code{n_yrs}; and \code{$Rec}, \code{$SSB}, \code{$Dynamic_SSB0} and
#'   \code{$Total_Biom} \code{[n_pop × n_regions × n_yrs × n_sims]}, with
#'   \code{$eff_SSB} \code{[n_pop × n_yrs × n_sims]}, \code{$ln_RecDevs} on the SSB
#'   dims and \code{$ln_InitDevs} \code{[n_pop × n_regions × (n_ages - 1) × n_sexes
#'   × n_sims]}.
#'
#'   Fishery: \code{$ObsCatch} and \code{$TrueCatch} \code{[n_regions × n_yrs ×
#'   n_seas × n_fish_fleets × n_sims]}, with \code{$ObsFishIdx},
#'   \code{$TrueFishIdx}, \code{$ObsDiscard} and \code{$TrueDiscard} on the same
#'   dims; \code{$ObsFishAgeComps} and \code{$ObsFishAgeComps_discard} with
#'   \code{n_obs_ages × n_sexes} before the fleet dim, and
#'   \code{$ObsFishLenComps} and \code{$ObsFishLenComps_discard} with
#'   \code{n_lens × n_sexes}. Each has a \code{_pop} counterpart with a leading
#'   \code{n_pop}. The true catch and discards at age and length are \code{$CAA},
#'   \code{$DAA} \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes ×
#'   n_fish_fleets × n_sims]} and \code{$CAL}, \code{$DAL} with \code{n_lens} in
#'   place of \code{n_ages}.
#'
#'   Survey: \code{$ObsSrvIdx} and \code{$TrueSrvIdx} \code{[n_regions × n_yrs ×
#'   n_seas × n_srv_fleets × n_sims]}, \code{$ObsSrvAgeComps} and
#'   \code{$ObsSrvLenComps} with the bin and sex dims before the fleet dim, each
#'   with a \code{_pop} counterpart, and the true \code{$SrvIAA} \code{[n_pop ×
#'   n_regions × n_yrs × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]} and
#'   \code{$SrvIAL} with \code{n_lens} in place of \code{n_ages}.
#'
#' @export Setup_Sim_Containers
#' @family Simulation Setup
Setup_Sim_Containers <- function(sim_list) {

  # Biological Containers
  sim_list$NAA_aft = sim_list$NAA_bef = sim_list$NAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs + 1, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  sim_list$NAA0 <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs + 1, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
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
  aa_fish <- c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages,
               sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)
  aa_srv <- c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages,
              sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)
  sim_list$TrueCatchAA <- array(0, dim = aa_fish)
  sim_list$ObsCatchAA <- array(0, dim = aa_fish)
  sim_list$TrueDiscardAA <- array(0, dim = aa_fish)
  sim_list$ObsDiscardAA <- array(0, dim = aa_fish)
  sim_list$TrueSrvIdxAA <- array(0, dim = aa_srv)
  sim_list$ObsSrvIdxAA <- array(0, dim = aa_srv)
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
