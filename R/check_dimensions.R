#' Helper function used to check data dimensions to ensure they are correct
#'
#' @param x Object to evaluate
#' @param n_regions Number of regions
#' @param n_years Number of years
#' @param n_ages Number of ages
#' @param n_lens Number of lengths
#' @param n_sexes Number of sexes
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#' @param what charcter specifying what to be evaluated
#' @param n_seas Number of seasons
#' @param max_tag_liberty Maximum tag liberty to track a given cohort
#' @param n_tag_cohorts Number of tag cohorts to fit to
#' @param n_pop Number of populations
#'
#' @keywords internal
check_data_dimensions <- function(x,
                                  n_pop = NULL,
                                  n_regions = NULL,
                                  n_years = NULL,
                                  n_ages = NULL,
                                  n_seas = NULL,
                                  n_lens = NULL,
                                  n_sexes = NULL,
                                  n_fish_fleets = NULL,
                                  n_srv_fleets = NULL,
                                  max_tag_liberty = NULL,
                                  n_tag_cohorts = NULL,
                                  what
                                  ) {

# Biologicals -------------------------------------------------------------

  # Weight at age (spawning), maturity
  if(what %in% c('WAA', 'MatAA')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop,, n_regions, n_years, n_seas, n_ages, and n_sexes"))
  }

  if(what %in% c("Fixed_natmort")) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_ages, n_sexes)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_ages, and n_sexes"))
  }

  # weight at age for the fishery
  if(what == 'WAA_fish') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets)) != 7)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, and n_fish_fleets"))
  }

  # weight at age for the survey
  if(what == 'WAA_srv') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets)) != 7)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, and n_srv_fleets"))
  }

  if(what == 'AgeingError') { # Not checking the age dimension
    if(sum(dim(x)[1] == n_ages) != 1)
      stop("Dimensions of AgeingError are not correct. Should be n_ages, number of observed composition ages")
  }

  if(what == 'AgeingError_t') { # Not checking the age dimension
    if(sum(dim(x)[1:2] == c(n_years, n_ages)) != 2)
      stop("Dimensions of AgeingError are not correct. Should be n_years, n_ages, number of observed composition ages")
  }

  if(what == 'SizeAgeTrans') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_lens, n_ages, n_sexes)) != 7)
       stop("Dimensions of SizeAgeTrans are not correct. Should be n_pop, n_regions, n_years, n_seas, n_lens, n_ages, and n_sexes")
  }

  if(what == 'Fixed_Movement') {
    if(sum(dim(x) == c(n_pop, n_regions, n_regions, n_years, n_seas, n_ages, n_sexes)) != 7)
      stop("Fixed Movement Matrix does not have the correct dimensions. This should be n_pop, n_regions, n_regions, n_years, n_seas, n_ages, n_sexes")
  }

  if(what == 'sgl_seas_spawning_movement') {
    if(sum(dim(x) == c(n_pop, n_regions, n_regions, n_years, n_ages, n_sexes)) != 6)
      stop("Fixed Movement Matrix does not have the correct dimensions. This should be n_pop, n_regions, n_regions, n_years, n_ages, n_sexes")
  }


# Fishery Stuff -----------------------------------------------------------------

  if(what %in% c('ObsCatch', "UseCatch", 'ObsFishIdx', 'ObsFishIdx_SE', 'UseFishIdx', 'UseFishAgeComps', 'UseFishLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_fish_fleets)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_seas, n_fish_fleets"))
  }


  if(what %in% c("FishAgeComps_LikeType", "FishLenComps_LikeType")) {
    if(length(x) != n_fish_fleets)
      stop(paste(what, "needs to have a length of n_fish_fleets"))
  }

  if(what == "ObsFishAgeComps") { # Not checking the age dimension
    if(sum(dim(x)[-4] == c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)) != 5)
      stop(paste("ObsFishAgeComps is not the correct dimension. Should be n_regions, n_years, n_seas, number of observed composition ages, n_sexes, n_fish_fleets"))
  }

  if(what == "ObsFishLenComps") {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets)) != 6)
      stop(paste("ObsFishLenComps is not the correct dimension. Should be n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('ISS_FishLenComps', 'ISS_FishAgeComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)) != 5)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_seas, n_sexes, n_fish_fleets"))
  }


# Survey Stuff ------------------------------------------------------------

  if(what %in% c('ObsSrvIdx', 'ObsSrvIdx_SE', 'UseSrvIdx', 'UseSrvAgeComps', 'UseSrvLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_srv_fleets)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_seas, n_srv_fleets"))
  }

  if(what %in% c("SrvAgeComps_LikeType", "SrvLenComps_LikeType")) {
    if(length(x) != n_srv_fleets)
      stop(paste(what, "needs to have a length of n_srv_fleets"))
  }

  if(what == "ObsSrvAgeComps") { # Not checking the age dimension
    if(sum(dim(x)[-4] == c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets)) != 5)
      stop(paste("ObsSrvAgeComps is not the correct dimension. Should be n_regions, n_years, n_seas, number of observed composition ages, n_sexes, n_srv_fleets"))
  }

  if(what == "ObsSrvLenComps") {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets)) != 6)
      stop(paste("ObsSrvLenComps is not the correct dimension. Should be n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets"))
  }

  if(what %in% c('ISS_SrvLenComps', 'ISS_SrvAgeComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets)) != 5)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_seas, n_sexes, n_srv_fleets"))
  }


  # Tagging Stuff -----------------------------------------------------------
  if(what %in% c("conv_tagged_fish")) {
    if(sum(dim(x) == c(n_tag_cohorts, n_pop, n_ages, n_sexes)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_tag_cohorts, n_pop, n_ages, n_sexes"))
  }

  if(what %in% c("obs_conv_tag_fish_recap")) {
    if(sum(dim(x) == c(max_tag_liberty, n_seas, n_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets)) != 8)
      stop(paste(what, " is not the correct dimension. Should be max_tag_liberty, n_seas, n_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets"))
  }

}

#' Helper function used to check simulation object dimensions to ensure they are correct
#'
#' @param x Object to evaluate
#' @param n_regions Number of regions
#' @param n_years Number of years
#' @param n_ages Number of ages
#' @param n_lens Number of lengths
#' @param n_sexes Number of sexes
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#' @param n_sims Number of simulations
#' @param what character specifying what to be evaluated
#' @param n_seas Number of seasons
#' @param n_pop
#'
#' @keywords internal
check_sim_dimensions <- function(x,
                                 n_pop = NULL,
                                 n_regions = NULL,
                                 n_years = NULL,
                                 n_ages = NULL,
                                 n_lens = NULL,
                                 n_sexes = NULL,
                                 n_seas = NULL,
                                 n_fish_fleets = NULL,
                                 n_srv_fleets = NULL,
                                 n_sims = NULL,
                                 what
                                 ) {

  # Biologicals -------------------------------------------------------------
  if(what %in% c('WAA_input', 'MatAA_input')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_sims)) != 7)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_sims"))
  }

  if(what %in% c('natmort_input')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_ages, n_sexes, n_sims)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_ages, n_sexes, n_sims"))
  }

  if(what == 'WAA_fish_input') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets, n_sims)) != 8)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets, n_sims"))
  }

  if(what == 'WAA_srv_input') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets, n_sims)) != 8)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets, n_sims"))
  }

  if(what == 'SizeAgeTrans_input') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_lens, n_ages, n_sexes, n_sims)) != 8)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_lens, n_ages, n_sexes n_sims"))
  }

  # Fishing Stuff  -------------------------
  if(what %in% c('Fmort_input')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_fish_fleets, n_sims)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_fish_fleets, n_sims"))
  }

  if(what %in% c('fish_q_input')) {
    if(sum(dim(x) == c(n_regions, n_years, n_fish_fleets, n_sims)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_fish_fleets, n_sims"))
  }

  if(what == 'fish_sel_input') {
    if(sum(dim(x) == c(n_regions, n_years, n_ages, n_sexes, n_fish_fleets, n_sims)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_ages, n_sexes, n_fish_fleets, n_sims"))
  }

  if(what %in% c('ln_sigmaC', 'ObsFishIdx_SE')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_fish_fleets)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_fish_fleets"))
  }

  if(what %in% c('comp_fishage_like', 'ln_FishAge_theta_agg', 'FishAge_corr_pars_agg',
                 'comp_fishlen_like', 'ln_FishLen_theta_agg', 'FishLen_corr_pars_agg')) {
    if(length(x) != n_fish_fleets)
      stop(paste(what, "needs to have a length of n_fish_fleets"))
  }

  if(what %in% c('ISS_FishAgeComps', 'ISS_FishLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims"))
  }

  if(what %in% c('ln_FishAge_theta', 'ln_FishLen_theta')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_fish_fleets)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('FishAge_corr_pars', 'FishLen_corr_pars')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_fish_fleets, 2)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_fish_fleets, 2"))
  }

  if(what %in% c('FishAgeComps_Type', 'FishLenComps_Type')) {
    if(sum(dim(x) == c(n_years, n_fish_fleets)) != 2)
      stop(paste("Dimensions of", what, "are not correct. Should be n_years, n_fish_fleets"))
  }

  # Survey Stuff  -------------------------
  if(what %in% c('srv_q_input')) {
    if(sum(dim(x) == c(n_regions, n_years, n_srv_fleets, n_sims)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_srv_fleets, n_sims"))
  }

  if(what == 'srv_sel_input') {
    if(sum(dim(x) == c(n_regions, n_years, n_ages, n_sexes, n_srv_fleets, n_sims)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_ages, n_sexes, n_srv_fleets, n_sims"))
  }

  if(what %in% c('ObsSrvIdx_SE')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_srv_fleets)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_srv_fleets"))
  }

  if(what %in% c('comp_srvage_like', 'ln_SrvAge_theta_agg', 'SrvAge_corr_pars_agg',
                 'comp_srvlen_like', 'ln_SrvLen_theta_agg', 'SrvLen_corr_pars_agg')) {
    if(length(x) != n_srv_fleets)
      stop(paste(what, "needs to have a length of n_srv_fleets"))
  }

  if(what %in% c('ISS_SrvAgeComps', 'ISS_SrvLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims"))
  }

  if(what %in% c('ln_SrvAge_theta', 'ln_SrvLen_theta')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_srv_fleets)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_srv_fleets"))
  }

  if(what %in% c('SrvAge_corr_pars', 'SrvLen_corr_pars')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_srv_fleets, 2)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_srv_fleets, 2"))
  }

  if(what %in% c('SrvAgeComps_Type', 'SrvLenComps_Type')) {
    if(sum(dim(x) == c(n_years, n_srv_fleets)) != 2)
      stop(paste("Dimensions of", what, "are not correct. Should be n_years, n_srv_fleets"))
  }


  if(what %in% c('t_srv')) {
    if(sum(dim(x) == c(n_regions, n_seas, n_srv_fleets)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_seas, n_srv_fleets"))
  }

  # Recruitment Stuff  -------------------------
  if(what == 'sexratio_input') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_sexes, n_sims)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_sexes, n_sims"))
  }

  if(what %in% c("R0_input", "h_input")) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_sims)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_sims"))
  }

  if(what == 'ln_InitDevs_input') {
    if(sum(dim(x) == c(n_pop, n_regions, n_ages - 1, n_sims)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_ages - 1, n_sims"))
  }


  # Tagging Stuff -----------------------------------------------------------
  if(what == 'conv_tag_fish_reporting_input') {
    if(sum(dim(x) == c(n_regions, n_years, n_fish_fleets, n_sims)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_fish_fleets, n_sims"))
  }

}
