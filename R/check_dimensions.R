#' Validate Dimensions of Model Input Objects
#'
#' Internal utility that verifies the dimensional structure of model input
#' objects. The function checks that the supplied object \code{x} has the
#' expected dimensions for the specified object type \code{what}. Expected
#' dimension ordering depends on the data structure used in the population
#' model (e.g., population, region, year, season, age, sex, fleet).
#'
#' The argument \code{what} determines which dimension template is applied.
#' If the dimensions of \code{x} do not match the required structure, the
#' function stops with an informative error message.
#'
#' @param x Object to evaluate. Typically a numeric array (or vector) whose
#'   dimensions must match the expected structure associated with \code{what}.
#' @param n_pop Integer. Number of populations.
#' @param n_regions Integer. Number of spatial regions.
#' @param n_years Integer. Number of model years.
#' @param n_ages Integer. Number of age classes.
#' @param n_seas Integer. Number of seasons.
#' @param n_lens Integer. Number of length bins.
#' @param n_sexes Integer. Number of sexes.
#' @param n_fish_fleets Integer. Number of fishery fleets.
#' @param n_srv_fleets Integer. Number of survey fleets.
#' @param conv_tag_max_liberty Integer. Maximum tag liberty for conventional
#'   tagging data (number of seasons between release and recapture).
#' @param n_conv_tag_cohorts Integer. Number of conventional tagging cohorts.
#' @param what Character string identifying the type of object being validated.
#'   This determines the expected dimension ordering (e.g., biological inputs,
#'   fishery observations, survey observations, or tagging data).
#'
#' @details
#' The function supports validation of several model input classes including:
#'
#' \itemize{
#'   \item Biological quantities (e.g., weight-at-age, maturity, natural mortality)
#'   \item Movement and spatial processes
#'   \item Fishery observations (catch, indices, composition data)
#'   \item Survey observations (indices and compositions)
#'   \item Conventional tagging data and recaptures
#' }
#'
#' Some objects (e.g., ageing-error matrices or age compositions) may contain
#' additional dimensions representing observed ages or length bins. These
#' dimensions are not always validated explicitly because they depend on the
#' structure of the observed data.
#'
#' @return
#' No return value. The function is used for validation and will terminate
#' execution with an error if the dimensions of \code{x} do not match the
#' expected structure.
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
                                  conv_tag_max_liberty = NULL,
                                  n_conv_tag_cohorts = NULL,
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

  if(what == 'stray_rate') {
    if(sum(dim(x) == c(n_pop, n_years)) != 2)
      stop("stray_ratedoes not have the correct dimensions. This should be n_pop, n_years")
  }

# Fishery Stuff -----------------------------------------------------------------

  if(what %in% c('ObsCatch', "UseCatch", 'ObsDiscard', "UseDiscard",
                 'ObsFishIdx', 'ObsFishIdx_SE', 'UseFishIdx',
                 'UseFishAgeComps', 'UseFishLenComps',
                 'UseFishAgeComps_discard', 'UseFishLenComps_discard')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_fish_fleets)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_seas, n_fish_fleets"))
  }

  if(what %in% c('ObsCatch_pop', "UseCatch_pop", 'ObsDiscard_pop', "UseDiscard_pop",
                 'ObsFishIdx_pop', 'ObsFishIdx_pop_SE', 'UseFishIdx_pop',
                 'UseFishAgeComps_pop', 'UseFishLenComps_pop',
                 'UseFishAgeComps_discard_pop', 'UseFishLenComps_discard_pop')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_fish_fleets)) != 5)
      stop(paste(what, " is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_fish_fleets"))
  }

  if(what %in% c("FishAgeComps_LikeType", "FishLenComps_LikeType",
                 "pop_FishAgeComps_LikeType", "pop_FishLenComps_LikeType",
                 "FishAgeComps_discard_LikeType", "FishLenComps_discard_LikeType",
                 "pop_FishAgeComps_discard_LikeType", "pop_FishLenComps_discard_LikeType")) {
    if(length(x) != n_fish_fleets)
      stop(paste(what, "needs to have a length of n_fish_fleets"))
  }

  if(what %in% c("ObsFishAgeComps", "ObsFishAgeComps_discard")) { # not checking age dimension
    if(sum(dim(x)[-4] == c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)) != 5)
      stop(paste(what, "is not the correct dimension. Should be n_regions, n_years, n_seas, number of observed composition ages, n_sexes, n_fish_fleets"))
  }

  if(what %in% c("ObsFishLenComps", "ObsFishLenComps_discard")) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets)) != 6)
      stop(paste(what, "is not the correct dimension. Should be n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets"))
  }

  if(what %in% c("ObsFishAgeComps_pop", "ObsFishAgeComps_discard_pop")) { # not checking age dimension
    if(sum(dim(x)[-5] == c(n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets)) != 6)
      stop(paste(what, "is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, number of observed composition ages, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('fish_sel_input_age', 'ret_sel_input_age')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets)) != 7)
      stop(paste(what, "is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets"))
  }

  if(what %in% c("ObsFishLenComps_pop", "ObsFishLenComps_discard_pop", 'fish_sel_input_len', 'ret_sel_input_len')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets)) != 7)
      stop(paste(what, "is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('ISS_FishLenComps', 'ISS_FishAgeComps',
                 'ISS_FishLenComps_discard', 'ISS_FishAgeComps_discard')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)) != 5)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_seas, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('ISS_FishLenComps_pop', 'ISS_FishAgeComps_pop',
                 'ISS_FishLenComps_discard_pop', 'ISS_FishAgeComps_discard_pop')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets)) != 6)
      stop(paste(what, " is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets"))
  }


# Survey Stuff ------------------------------------------------------------

  if(what %in% c('ObsSrvIdx', 'ObsSrvIdx_SE', 'UseSrvIdx', 'UseSrvAgeComps', 'UseSrvLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_srv_fleets)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_seas, n_srv_fleets"))
  }

  if(what %in% c('ObsSrvIdx_pop', 'ObsSrvIdx_pop_SE', 'UseSrvIdx_pop', 'UseSrvAgeComps_pop', 'UseSrvLenComps_pop')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_srv_fleets)) != 5)
      stop(paste(what, " is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_srv_fleets"))
  }

  if(what %in% c("SrvAgeComps_LikeType", "SrvLenComps_LikeType",
                 "pop_SrvAgeComps_LikeType", "pop_SrvLenComps_LikeType")) {
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

  if(what %in% c('ISS_SrvLenComps_pop', 'ISS_SrvAgeComps_pop')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets)) != 6)
      stop(paste(what, " is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets"))
  }

  if(what %in% c("ObsSrvAgeComps_pop")) { # not checking age dimension
    if(sum(dim(x)[-5] == c(n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets)) != 6)
      stop(paste("ObsSrvAgeComps_pop is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, number of observed composition ages, n_sexes, n_srv_fleets"))
  }

  if(what %in% c('srv_sel_input_age')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets)) != 7)
      stop(paste(what, "is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets"))
  }

  if(what %in% c("ObsSrvLenComps_pop", 'srv_sel_input_len')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets)) != 7)
      stop(paste(what, "is not the correct dimension. Should be n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets"))
  }



  # Tagging Stuff -----------------------------------------------------------
  if(what %in% c("conv_tagged_fish")) {
    if(sum(dim(x) == c(n_conv_tag_cohorts, n_pop, n_ages, n_sexes)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_conv_tag_cohorts, n_pop, n_ages, n_sexes"))
  }

  if(what %in% c("obs_conv_tag_fish_recap")) {
    if(sum(dim(x) == c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets)) != 8)
      stop(paste(what, " is not the correct dimension. Should be conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets"))
  }

}


#' Validate Dimensions of Simulation Input Objects
#'
#' Internal utility that verifies the dimensional structure of simulation
#' input objects used in the operating model. The function checks that the
#' supplied object \code{x} has the expected dimensions for the object type
#' specified by \code{what}. If the dimensions do not match the required
#' structure, the function stops with an informative error.
#'
#' Simulation inputs generally include an additional dimension representing
#' the number of stochastic simulations (\code{n_sims}). The position of this
#' dimension depends on the object type but is typically the final dimension
#' of the array.
#'
#' @param x Object to evaluate. Typically a numeric array (or vector) whose
#'   dimensions must match the expected structure associated with
#'   \code{what}.
#' @param n_pop Integer. Number of populations.
#' @param n_regions Integer. Number of spatial regions.
#' @param n_years Integer. Number of model years.
#' @param n_ages Integer. Number of age classes.
#' @param n_seas Integer. Number of seasons.
#' @param n_lens Integer. Number of length bins.
#' @param n_sexes Integer. Number of sexes.
#' @param n_fish_fleets Integer. Number of fishery fleets.
#' @param n_srv_fleets Integer. Number of survey fleets.
#' @param n_sims Integer. Number of stochastic simulations.
#' @param what Character string identifying the type of object to validate.
#'   This determines the expected dimension ordering and structure.
#'
#' @details
#' The function validates the dimensional structure of several classes of
#' simulation inputs used in the operating model:
#'
#' \itemize{
#'   \item Biological processes (e.g., weight-at-age, maturity, natural mortality)
#'   \item Fishery processes (e.g., fishing mortality, selectivity, catchability)
#'   \item Survey processes (e.g., selectivity, catchability)
#'   \item Observation models for indices and composition data
#'   \item Recruitment and demographic processes
#'   \item Conventional tagging parameters
#' }
#'
#' Arrays generally follow the dimension ordering
#' \preformatted{
#' population × region × year × season × age × sex × fleet × simulation
#' }
#' although only the dimensions relevant to a given object are included.
#'
#' @return
#' No return value. The function is used for validation and terminates with an
#' error if the dimensions of \code{x} do not match the expected structure.
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
  if(what %in% c('Fmort_input', 'dmr_input')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_fish_fleets, n_sims)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_fish_fleets, n_sims"))
  }

  if(what %in% c('fish_q_input')) {
    if(sum(dim(x) == c(n_regions, n_years, n_fish_fleets, n_sims)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_fish_fleets, n_sims"))
  }

  if(what %in% c('fish_sel_input', 'ret_sel_input')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets, n_sims)) != 8)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets, n_sims"))
  }

  if(what %in% c('ln_sigmaC', 'ObsFishIdx_SE', 'ln_sigmaD')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_fish_fleets)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_fish_fleets"))
  }

  if(what %in% c('ln_sigmaC_pop', 'ObsFishIdx_pop_SE', 'ln_sigmaD_pop')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_fish_fleets)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_fish_fleets"))
  }

  if(what %in% c('comp_fishage_like', 'ln_FishAge_theta_agg', 'FishAge_corr_pars_agg',
                 'comp_fishlen_like', 'ln_FishLen_theta_agg', 'FishLen_corr_pars_agg',
                 'pop_comp_fishage_like', 'pop_comp_fishlen_like',
                 'comp_fishage_discard_like', 'ln_FishAge_discard_theta_agg', 'FishAge_discard_corr_pars_agg',
                 'comp_fishlen_discard_like', 'ln_FishLen_discard_theta_agg', 'FishLen_discard_corr_pars_agg',
                 'pop_comp_fishage_discard_like', 'pop_comp_fishlen_discard_like')) {
    if(length(x) != n_fish_fleets)
      stop(paste(what, "needs to have a length of n_fish_fleets"))
  }

  if(what %in% c('FishAge_pop_corr_pars_agg', 'FishLen_pop_corr_pars_agg', 'ln_FishAge_pop_theta_agg', 'ln_FishLen_pop_theta_agg',
                 'FishAge_discard_pop_corr_pars_agg', 'FishLen_discard_pop_corr_pars_agg', 'ln_FishAge_discard_pop_theta_agg', 'ln_FishLen_discard_pop_theta_agg')) {
    if(sum(dim(x) == c(n_pop, n_fish_fleets)) != 2)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_fish_fleets"))
  }

  if(what %in% c('ISS_FishAgeComps', 'ISS_FishLenComps', 'ISS_FishAgeComps_discard', 'ISS_FishLenComps_discard')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims"))
  }

  if(what %in% c('ISS_FishAgeComps_pop', 'ISS_FishLenComps_pop', 'ISS_FishAgeComps_discard_pop', 'ISS_FishLenComps_discard_pop')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims)) != 7)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims"))
  }

  if(what %in% c('ln_FishAge_theta', 'ln_FishLen_theta', 'ln_FishAge_discard_theta', 'ln_FishLen_discard_theta')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_fish_fleets)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('ln_FishAge_pop_theta', 'ln_FishLen_pop_theta', 'ln_FishAge_discard_pop_theta', 'ln_FishLen_discard_pop_theta')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_sexes, n_fish_fleets)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('FishAge_corr_pars', 'FishLen_corr_pars', 'FishAge_discard_corr_pars', 'FishLen_discard_corr_pars')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_fish_fleets, 2)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_fish_fleets, 2"))
  }

  if(what %in% c('FishAge_pop_corr_pars', 'FishLen_pop_corr_pars', 'FishAge_discard_pop_corr_pars', 'FishLen_discard_pop_corr_pars')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_sexes, n_fish_fleets, 2)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_sexes, n_fish_fleets, 2"))
  }

  if(what %in% c('FishAgeComps_Type', 'FishLenComps_Type', 'pop_FishAgeComps_Type', 'pop_FishLenComps_Type',
                 'FishAgeComps_discard_Type', 'FishLenComps_discard_Type', 'pop_FishAgeComps_discard_Type', 'pop_FishLenComps_discard_Type')) {
    if(sum(dim(x) == c(n_years, n_fish_fleets)) != 2)
      stop(paste("Dimensions of", what, "are not correct. Should be n_years, n_fish_fleets"))
  }

  # Survey Stuff  -------------------------
  if(what %in% c('srv_q_input')) {
    if(sum(dim(x) == c(n_regions, n_years, n_srv_fleets, n_sims)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_srv_fleets, n_sims"))
  }

  if(what == 'srv_sel_input') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets, n_sims)) != 8)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets, n_sims"))
  }

  if(what %in% c('ObsSrvIdx_SE')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_srv_fleets)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_srv_fleets"))
  }

  if(what %in% c('ObsSrvIdx_pop_SE')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_srv_fleets)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_srv_fleets"))
  }

  if(what %in% c('comp_srvage_like', 'ln_SrvAge_theta_agg', 'SrvAge_corr_pars_agg',
                 'comp_srvlen_like', 'ln_SrvLen_theta_agg', 'SrvLen_corr_pars_agg',
                 "pop_comp_srvage_like", "pop_comp_srvlen_like")) {
    if(length(x) != n_srv_fleets)
      stop(paste(what, "needs to have a length of n_srv_fleets"))
  }

  if(what %in% c('ISS_SrvAgeComps', 'ISS_SrvLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims)) != 6)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims"))
  }

  if(what %in% c('ISS_SrvAgeComps_pop', 'ISS_SrvLenComps_pop')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims)) != 7)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims"))
  }

  if(what %in% c('ln_SrvAge_theta', 'ln_SrvLen_theta')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_srv_fleets)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_srv_fleets"))
  }

  if(what %in% c('ln_SrvAge_pop_theta', 'ln_SrvLen_pop_theta')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_sexes, n_srv_fleets)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_sexes, n_srv_fleets"))
  }

  if(what %in% c('SrvAge_corr_pars', 'SrvLen_corr_pars')) {
    if(sum(dim(x) == c(n_regions, n_sexes, n_srv_fleets, 2)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_sexes, n_srv_fleets, 2"))
  }

  if(what %in% c("SrvAge_pop_corr_pars_agg", "SrvLen_pop_corr_pars_agg", "ln_SrvAge_pop_theta_agg", "ln_SrvLen_pop_theta_agg")) {
    if(sum(dim(x) == c(n_pop, n_srv_fleets)) != 2)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_srv_fleets"))
  }

  if(what %in% c('SrvAgeComps_Type', 'SrvLenComps_Type')) {
    if(sum(dim(x) == c(n_years, n_srv_fleets)) != 2)
      stop(paste("Dimensions of", what, "are not correct. Should be n_years, n_srv_fleets"))
  }


  if(what %in% c('t_srv')) {
    if(sum(dim(x) == c(n_regions, n_seas, n_srv_fleets)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_seas, n_srv_fleets"))
  }

  if(what %in% c('SrvAge_pop_corr_pars', 'SrvLen_pop_corr_pars')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_sexes, n_srv_fleets, 2)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_sexes, n_srv_fleets, 2"))
  }

  if(what %in% c('SrvAge_pop_corr_pars', 'SrvLen_pop_corr_pars')) {
    if(sum(dim(x) == c(n_pop, n_regions, n_sexes, n_srv_fleets, 2)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_sexes, n_srv_fleets, 2"))
  }

  if(what %in% c('SrvAgeComps_Type', 'SrvLenComps_Type', "pop_SrvAgeComps_Type", "pop_SrvLenComps_Type")) {
    if(sum(dim(x) == c(n_years, n_srv_fleets)) != 2)
      stop(paste("Dimensions of", what, "are not correct. Should be n_years, n_srv_fleets"))
  }

  # Recruitment Stuff  -------------------------
  if(what == 'sexratio_input') {
    if(sum(dim(x) == c(n_pop, n_regions, n_years, n_sexes, n_sims)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_regions, n_years, n_sexes, n_sims"))
  }

  if(what == 'stray_rate_input') {
    if(sum(dim(x) == c(n_pop, n_years, n_sims)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_years, n_sims"))
  }

  if(what == 'rec_seas_prop') {
    if(sum(dim(x) == c(n_pop, n_seas, n_sims)) != 3)
      stop(paste("Dimensions of", what, "are not correct. Should be n_pop, n_seas, n_sims"))
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
