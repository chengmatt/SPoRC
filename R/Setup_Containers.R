#' Setup containers for simulation and output
#'
#' @param sim_list List from `Setup_Sim_Dim()` containing core simulation dimensions
#' (`n_regions`, `n_yrs`, `n_seas` `n_ages`, `n_sexes`, `n_sims`, `n_fish_fleets`, `n_srv_fleets`, `n_obs_ages`, `n_lens`).
#' The function appends container arrays for biological, fishery, and survey quantities.
#'
#' @export Setup_Sim_Containers
#' @family Simulation Setup
Setup_Sim_Containers <- function(sim_list) {

  # Biological Containers
  sim_list$NAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs+1, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # numbers at age
  sim_list$NAA0 <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs+1,  sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # unfished numbers at age
  sim_list$ZAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # total mortality at age
  sim_list$Rec <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # recruitment
  sim_list$SSB <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # spawning biomass
  sim_list$Dynamic_SSB0 <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs,  sim_list$n_sims)) # Dynamic unfished spawning biomass
  sim_list$Total_Biom <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # total biomass
  sim_list$ln_RecDevs <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # rec devs
  sim_list$ln_InitDevs <- array(0, dim = c(sim_list$n_regions, sim_list$n_ages - 1, sim_list$n_sims)) # initial age devs

  # Fishery Containers
  sim_list$ObsCatch <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # observed catch
  sim_list$ObsFishAgeComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # observed fishery age comps
  sim_list$ObsFishLenComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # observed fishery length comps
  sim_list$CAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # catch at age
  sim_list$CAL <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)) # catch at length
  sim_list$TrueCatch <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # true catches
  sim_list$ObsFishIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # observed fishery index
  sim_list$TrueFishIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)) # true fishery index

  # Survey Containers
  sim_list$ObsSrvIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims)) # observed survey index
  sim_list$TrueSrvIdx <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets, sim_list$n_sims)) # true survey index
  sim_list$ObsSrvAgeComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_obs_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # observed survey age comps
  sim_list$ObsSrvLenComps <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # observed survey length comps
  sim_list$SrvIAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # survey index at age
  sim_list$SrvIAL <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)) # survey index at length

  return(sim_list)
}
