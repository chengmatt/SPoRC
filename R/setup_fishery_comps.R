# Stage 1 of 3: model setup
#
# Fishery index and composition inputs, for both the retained stream
# (Setup_Mod_FishIdx_and_Comps) and the discarded stream
# (Setup_Mod_Discard_Comps). Composition likelihood choice and its
# overdispersion and correlation parameters are set here.

#' Set up discard age and length composition inputs
#'
#' @param input_list Main model input list containing data, parameters, and mapping structures.
#'
#' @param ObsFishAgeComps_discard 5D array of observed discard age compositions:
#'   \code{[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param UseFishAgeComps_discard Matrix indicating whether discard age comps are used:
#'   \code{[n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishAgeComps_discard Optional ISS (effective sample size) array for discard age comps:
#'   \code{[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param ObsFishLenComps_discard 6D array of observed discard length compositions:
#'   \code{[n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]}
#'
#' @param UseFishLenComps_discard Matrix indicating whether discard length comps are used:
#'   \code{[n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishLenComps_discard Optional ISS array for discard length comps:
#'   \code{[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param FishAgeComps_discard_LikeType Character vector (length n_fish_fleets) specifying likelihood type:
#'   one of \code{c("none","Multinomial","Dirichlet-Multinomial","iid-Logistic-Normal","1d-Logistic-Normal","2d-Logistic-Normal")}
#'
#' @param FishLenComps_discard_LikeType Character vector (length n_fish_fleets) specifying likelihood type
#'
#' @param FishAgeComps_discard_Type List/encoded character strings defining composition structure by year and fleet.
#'
#' @param FishLenComps_discard_Type List/encoded character strings defining composition structure by year and fleet.
#'
#' @param ObsFishAgeComps_discard_pop 6D array of population-specific discard age comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param UseFishAgeComps_discard_pop 5D array indicating use of population age comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishAgeComps_discard_pop Optional ISS array for population discard age comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param ObsFishLenComps_discard_pop 7D array of population-specific discard length comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]}
#'
#' @param UseFishLenComps_discard_pop 5D array indicating use of population length comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishLenComps_discard_pop Optional ISS array for population length comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param FishAgeComps_discard_pop_LikeType Character vector (length n_fish_fleets) for population age likelihood types.
#'
#' @param FishLenComps_discard_pop_LikeType Character vector (length n_fish_fleets) for population length likelihood types.
#'
#' @param FishAgeComps_discard_pop_Type Encoded structure definitions for population age comps by year/fleet.
#'
#' @param FishLenComps_discard_pop_Type Encoded structure definitions for population length comps by year/fleet.
#'
#' @param ... Optional starting values for parameters (e.g., dispersion, correlation)
#'
#' @return The input \code{input_list} updated with:
#' \itemize{
#'   \item discard age and length composition data structures
#'   \item ISS (effective sample size) arrays (computed or supplied)
#'   \item likelihood type mappings (integer-coded)
#'   \item composition type matrices by year and fleet
#'   \item population-specific discard composition structures
#'   \item parameter arrays for dispersion and correlation
#'   \item mapping configurations for estimation
#' }
#'
#' @details
#' All composition arrays follow consistent indexing conventions:
#' \itemize{
#'   \item Age compositions: \code{[region, year, season, sex, fleet]}
#'   \item Length compositions: \code{[region, year, season, length, sex, fleet]}
#'   \item Population age comps: \code{[pop, region, year, season, sex, fleet]}
#'   \item Population length comps: \code{[pop, region, year, season, length, sex, fleet]}
#' }
#'
#'
#' @keywords internal
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Discard_Comps     <- function(input_list,
                                        ObsFishAgeComps_discard,
                                        UseFishAgeComps_discard,
                                        ISS_FishAgeComps_discard,
                                        ObsFishLenComps_discard,
                                        UseFishLenComps_discard,
                                        ISS_FishLenComps_discard,
                                        FishAgeComps_discard_LikeType,
                                        FishLenComps_discard_LikeType,
                                        FishAgeComps_discard_Type,
                                        FishLenComps_discard_Type,
                                        ObsFishAgeComps_discard_pop,
                                        UseFishAgeComps_discard_pop,
                                        ISS_FishAgeComps_discard_pop,
                                        ObsFishLenComps_discard_pop,
                                        UseFishLenComps_discard_pop,
                                        ISS_FishLenComps_discard_pop,
                                        FishAgeComps_discard_pop_LikeType ,
                                        FishLenComps_discard_pop_LikeType,
                                        FishAgeComps_discard_pop_Type,
                                        FishLenComps_discard_pop_Type,
                                        ...
) {

  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_FishIdx_and_Comps <- c(input_list$config$Setup_Mod_FishIdx_and_Comps, mget(names(formals()))[-1])

  # Input Validation ---------------------------------------------------------
  # Discard Fishery compositions
  check_data_dimensions(ObsFishAgeComps_discard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps_discard')
  check_data_dimensions(UseFishAgeComps_discard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps_discard')
  check_data_dimensions(UseFishLenComps_discard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps_discard')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsFishLenComps_discard, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = obs_len_bins(input_list), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps_discard')
  if(!is.null(ISS_FishAgeComps_discard)) check_data_dimensions(ISS_FishAgeComps_discard, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps_discard')
  if(!is.null(ISS_FishLenComps_discard)) check_data_dimensions(ISS_FishLenComps_discard, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps_discard')
  check_data_dimensions(FishAgeComps_discard_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_discard_LikeType')
  check_data_dimensions(FishLenComps_discard_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_discard_LikeType')
  if(!all(FishAgeComps_discard_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_discard_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_discard_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_discard_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # Discard Fishery compositions (population-specific)
  if(any(UseFishAgeComps_discard_pop == 1)) check_data_dimensions(ObsFishAgeComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps_discard_pop')
  check_data_dimensions(UseFishAgeComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps_discard_pop')
  check_data_dimensions(UseFishLenComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps_discard_pop')
  if(input_list$data$fit_lengths == 1 && any(UseFishLenComps_discard_pop == 1)) check_data_dimensions(ObsFishLenComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = obs_len_bins(input_list), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps_discard_pop')
  if(!is.null(ISS_FishAgeComps_discard_pop)) check_data_dimensions(ISS_FishAgeComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps_discard_pop')
  if(!is.null(ISS_FishLenComps_discard_pop)) check_data_dimensions(ISS_FishLenComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps_discard_pop')
  check_data_dimensions(FishAgeComps_discard_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_discard_pop_LikeType')
  check_data_dimensions(FishLenComps_discard_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_discard_pop_LikeType')
  if(!all(FishAgeComps_discard_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_discard_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_discard_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_discard_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # checking to make sure defaults are not applied
  if(any(UseFishAgeComps_discard_pop == 1)) {
    if(is.null(ObsFishAgeComps_discard_pop)) stop("ObsFishAgeComps_discard_pop is NULL, but UseFishAgeComps_discard_pop contains 1s!")
    if(any(str_detect(FishAgeComps_discard_pop_LikeType, "none"))) warning("FishAgeComps_discard_pop_LikeType has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
    if(any(str_detect(FishAgeComps_discard_pop_Type, "none"))) warning("FishAgeComps_discard_pop_Type has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
  }

  if(any(UseFishLenComps_discard_pop == 1)) {
    if(is.null(ObsFishLenComps_discard_pop)) stop("ObsFishLenComps_discard_pop is NULL, but UseFishAgeComps_discard_pop contains 1s!")
    if(any(str_detect(FishLenComps_discard_pop_LikeType, "none"))) warning("FishLenComps_discard_pop_LikeType has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
    if(any(str_detect(FishLenComps_discard_pop_Type, "none"))) warning("FishLenComps_discard_pop_Type has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
  }

  # Discard Fishery Age Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishage_discard_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_discard_LikeType[f] == 'none') comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 999)
    if(FishAgeComps_discard_LikeType[f] == "Multinomial") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 0)
    if(FishAgeComps_discard_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 1)
    if(FishAgeComps_discard_LikeType[f] == "iid-Logistic-Normal") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 2)
    if(FishAgeComps_discard_LikeType[f] == "1d-Logistic-Normal") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 3)
    if(FishAgeComps_discard_LikeType[f] == "2d-Logistic-Normal") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 4)
    collect_message(paste("Discard Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_discard_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_discard_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_discard_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_discard_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_discard_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_discard_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_discard_like_vals[fleet] == 4) stop("Discard Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_discard_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_discard_Type_Mat))) stop("FishAgeComps_discard_Type is returning an NA. Did you update the year range of FishAgeComps_discard_Type?")

  # Specifying composition likelihood for population-specific data
  comp_fishage_discard_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_discard_pop_LikeType[f] == 'none') comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 999)
    if(FishAgeComps_discard_pop_LikeType[f] == "Multinomial") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 0)
    if(FishAgeComps_discard_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 1)
    if(FishAgeComps_discard_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 2)
    if(FishAgeComps_discard_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 3)
    if(FishAgeComps_discard_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 4)
    collect_message(paste("Discard Population Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_discard_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_discard_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_discard_pop_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_discard_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_discard_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_discard_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_discard_pop_like_vals[fleet] == 4) stop("Discard Population Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_discard_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_discard_pop_Type_Mat))) stop("FishAgeComps_discard_pop_Type is returning an NA. Did you update the year range of FishAgeComps_discard_pop_Type?")

  # Fishery Length Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishlen_discard_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_discard_LikeType[f] == 'none') comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 999)
    if(FishLenComps_discard_LikeType[f] == "Multinomial") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 0)
    if(FishLenComps_discard_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 1)
    if(FishLenComps_discard_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 2)
    if(FishLenComps_discard_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 3)
    if(FishLenComps_discard_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 4)
    collect_message(paste("Discard Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_discard_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_discard_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_discard_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_discard_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # define composition types
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_discard_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_discard_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_discard_like_vals[fleet] == 4) stop("Discard Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_discard_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_discard_Type_Mat))) stop("FishLenComps_discard_Type_Mat is returning an NA. Did you update the year range of FishLenComps_discard_Type_Mat?")


  # Specifying composition likelihood for population-specific data
  comp_fishlen_discard_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_discard_pop_LikeType[f] == 'none') comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 999)
    if(FishLenComps_discard_pop_LikeType[f] == "Multinomial") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 0)
    if(FishLenComps_discard_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 1)
    if(FishLenComps_discard_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 2)
    if(FishLenComps_discard_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 3)
    if(FishLenComps_discard_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 4)
    collect_message(paste("Discard Population Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_discard_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_discard_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_discard_pop_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_discard_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_discard_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_discard_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_discard_pop_like_vals[fleet] == 4) stop("Discard Population Len composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_discard_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_discard_pop_Type_Mat))) stop("FishLenComps_discard_pop_Type is returning an NA. Did you update the year range of FishLenComps_discard_pop_Type?")

  # ISS Munging -------------------------------------------------------------

  # Discard Fishery Ages
  if(is.null(ISS_FishAgeComps_discard)) {
    collect_message("No ISS is specified for FishAgeComps_discard. ISS weighting is calculated by summing up values from ObsFishAgeComps_discard each year")
    ISS_FishAgeComps_discard <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0) or joint across sexes
          if(FishAgeComps_discard_Type_Mat[y,f] == 0) ISS_FishAgeComps_discard[1,y,seas,1,f] <- sum(ObsFishAgeComps_discard[,y,seas,,,f])
          # if split by region and sex
          if(FishAgeComps_discard_Type_Mat[y,f] == 1) ISS_FishAgeComps_discard[,y,seas,,f] <- apply(ObsFishAgeComps_discard[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishAgeComps_discard_Type_Mat[y,f] == 2) ISS_FishAgeComps_discard[,y,seas,1,f] <- apply(ObsFishAgeComps_discard[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps_discard)) {
    collect_message("No ISS is specified for FishLenComps_discard. ISS weighting is calculated by summing up values from ObsFishLenComps_discard each year")
    ISS_FishLenComps_discard <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(FishLenComps_discard_Type_Mat[y,f] == 0) ISS_FishLenComps_discard[1,y,seas,1,f] <- sum(ObsFishLenComps_discard[,y,seas,,,f])
          # if split by region and sex
          if(FishLenComps_discard_Type_Mat[y,f] == 1) ISS_FishLenComps_discard[,y,seas,,f] <- apply(ObsFishLenComps_discard[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishLenComps_discard_Type_Mat[y,f] == 2) ISS_FishLenComps_discard[,y,seas,1,f] <- apply(ObsFishLenComps_discard[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Ages
  if(is.null(ISS_FishAgeComps_discard_pop)) {
    collect_message("No ISS is specified for pop_FishAgeComps_discard. ISS weighting is calculated by summing up values from ObsFishAgeComps_discard_pop each year")
    ISS_FishAgeComps_discard_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0) or joint across sexes
            if(FishAgeComps_discard_pop_Type_Mat[y,f] == 0) ISS_FishAgeComps_discard_pop[p,1,y,seas,1,f] <- sum(ObsFishAgeComps_discard_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishAgeComps_discard_pop_Type_Mat[y,f] == 1) ISS_FishAgeComps_discard_pop[p,,y,seas,,f] <- apply(ObsFishAgeComps_discard_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishAgeComps_discard_pop_Type_Mat[y,f] == 2) ISS_FishAgeComps_discard_pop[p,,y,seas,1,f] <- apply(ObsFishAgeComps_discard_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps_discard_pop)) {
    collect_message("No ISS is specified for pop_FishLenComps_discard. ISS weighting is calculated by summing up values from ObsFishLenComps_discard_pop each year")
    ISS_FishLenComps_discard_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0)
            if(FishLenComps_discard_pop_Type_Mat[y,f] == 0) ISS_FishLenComps_discard_pop[p,1,y,seas,1,f] <- sum(ObsFishLenComps_discard_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishLenComps_discard_pop_Type_Mat[y,f] == 1) ISS_FishLenComps_discard_pop[p,,y,seas,,f] <- apply(ObsFishLenComps_discard_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishLenComps_discard_pop_Type_Mat[y,f] == 2) ISS_FishLenComps_discard_pop[p,,y,seas,1,f] <- apply(ObsFishLenComps_discard_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_FishAgeComps_discard <- ISS_FishAgeComps_discard
  input_list$data$ISS_FishLenComps_discard <- ISS_FishLenComps_discard
  input_list$data$ISS_FishAgeComps_discard_pop <- ISS_FishAgeComps_discard_pop
  input_list$data$ISS_FishLenComps_discard_pop <- ISS_FishLenComps_discard_pop
  input_list$data$ObsFishAgeComps_discard <- ObsFishAgeComps_discard
  input_list$data$UseFishAgeComps_discard <- UseFishAgeComps_discard
  input_list$data$ObsFishLenComps_discard <- ObsFishLenComps_discard
  input_list$data$UseFishLenComps_discard <- UseFishLenComps_discard
  input_list$data$ObsFishAgeComps_discard_pop <- ObsFishAgeComps_discard_pop
  input_list$data$UseFishAgeComps_discard_pop <- UseFishAgeComps_discard_pop
  input_list$data$ObsFishLenComps_discard_pop <- ObsFishLenComps_discard_pop
  input_list$data$UseFishLenComps_discard_pop <- UseFishLenComps_discard_pop
  input_list$data$FishAgeComps_discard_LikeType <- comp_fishage_discard_like_vals
  input_list$data$FishLenComps_discard_LikeType <- comp_fishlen_discard_like_vals
  input_list$data$FishAgeComps_discard_pop_LikeType <- comp_fishage_discard_pop_like_vals
  input_list$data$FishLenComps_discard_pop_LikeType <- comp_fishlen_discard_pop_like_vals
  input_list$data$FishAgeComps_discard_Type <- FishAgeComps_discard_Type_Mat
  input_list$data$FishLenComps_discard_Type <- FishLenComps_discard_Type_Mat
  input_list$data$FishAgeComps_discard_pop_Type <- FishAgeComps_discard_pop_Type_Mat
  input_list$data$FishLenComps_discard_pop_Type <- FishLenComps_discard_pop_Type_Mat

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the fishery age comps
  if("ln_FishAge_discard_theta" %in% names(starting_values)) input_list$par$ln_FishAge_discard_theta <- starting_values$ln_FishAge_discard_theta
  else input_list$par$ln_FishAge_discard_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery age comps
  if("FishAge_discard_corr_pars" %in% names(starting_values)) input_list$par$FishAge_discard_corr_pars <- starting_values$FishAge_discard_corr_pars
  else input_list$par$FishAge_discard_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishAge_discard_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_discard_theta_agg <- starting_values$ln_FishAge_discard_theta_agg
  else input_list$par$ln_FishAge_discard_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  # aggregated correlation parameters
  if("FishAge_discard_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_discard_corr_pars_agg <- starting_values$FishAge_discard_corr_pars_agg
  else input_list$par$FishAge_discard_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for fishery length comps
  if("ln_FishLen_discard_theta" %in% names(starting_values)) input_list$par$ln_FishLen_discard_theta <- starting_values$ln_FishLen_discard_theta
  else input_list$par$ln_FishLen_discard_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery length comps
  if("FishLen_discard_corr_pars" %in% names(starting_values)) input_list$par$FishLen_discard_corr_pars <- starting_values$FishLen_discard_corr_pars
  else input_list$par$FishLen_discard_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishLen_discard_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_discard_theta_agg <- starting_values$ln_FishLen_discard_theta_agg
  else input_list$par$ln_FishLen_discard_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  if("FishLen_discard_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_discard_corr_pars_agg <- starting_values$FishLen_discard_corr_pars_agg
  else input_list$par$FishLen_discard_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for the population fishery age comps
  if("ln_FishAge_discard_pop_theta" %in% names(starting_values)) input_list$par$ln_FishAge_discard_pop_theta <- starting_values$ln_FishAge_discard_pop_theta
  else input_list$par$ln_FishAge_discard_pop_theta <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery age comps
  if("FishAge_discard_pop_corr_pars" %in% names(starting_values)) input_list$par$FishAge_discard_pop_corr_pars <- starting_values$FishAge_discard_pop_corr_pars
  else input_list$par$FishAge_discard_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishAge_discard_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_discard_pop_theta_agg <- starting_values$ln_FishAge_discard_pop_theta_agg
  else input_list$par$ln_FishAge_discard_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # aggregated population correlation parameters
  if("FishAge_discard_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_discard_pop_corr_pars_agg <- starting_values$FishAge_discard_pop_corr_pars_agg
  else input_list$par$FishAge_discard_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # Dispersion parameters for population fishery length comps
  if("ln_FishLen_discard_pop_theta" %in% names(starting_values)) input_list$par$ln_FishLen_discard_pop_theta <- starting_values$ln_FishLen_discard_pop_theta
  else input_list$par$ln_FishLen_discard_pop_theta <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery length comps
  if("FishLen_discard_pop_corr_pars" %in% names(starting_values)) input_list$par$FishLen_discard_pop_corr_pars <- starting_values$FishLen_discard_pop_corr_pars
  else input_list$par$FishLen_discard_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishLen_discard_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_discard_pop_theta_agg <- starting_values$ln_FishLen_discard_pop_theta_agg
  else input_list$par$ln_FishLen_discard_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  if("FishLen_discard_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_discard_pop_corr_pars_agg <- starting_values$FishLen_discard_pop_corr_pars_agg
  else input_list$par$FishLen_discard_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishAge", discard = TRUE)
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishLen", discard = TRUE)
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishAge", discard = TRUE)
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishLen", discard = TRUE)

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishAge", discard = TRUE, has_pop = TRUE)
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishLen", discard = TRUE, has_pop = TRUE)
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishAge", discard = TRUE, has_pop = TRUE)
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishLen", discard = TRUE, has_pop = TRUE)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Set up discards, fishery index, age composition, and length composition inputs
#'
#' Populates \code{input_list} with observed fishery indices, age compositions,
#' and length compositions (both pooled and population-specific) along with
#' their usage indicators, likelihood types, composition structure types, input
#' sample sizes, and overdispersion and correlation parameter starting values
#' and mappings. Must be called after \code{\link{Setup_Mod_Catch_and_F}}.
#'
#' When \code{ISS_FishAgeComps}, \code{ISS_FishLenComps},
#' \code{ISS_FishAgeComps_pop}, or \code{ISS_FishLenComps_pop} are \code{NULL},
#' input sample sizes are derived automatically by summing the observed
#' composition arrays within each year-fleet-season-region cell, consistent
#' with the specified composition type.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsFishIdx Observed fishery CPUE or biomass index array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#' @param ObsFishIdx_SE Standard errors of \code{ObsFishIdx} on the log scale,
#'   same dimensions as \code{ObsFishIdx}.
#' @param UseFishIdx Binary indicator array \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = include index in the likelihood; \code{0} = exclude.
#' @param fish_idx_ages Per-fleet selection of which ages contribute to the
#'   index total. Either a list with one element per fishery fleet, where each
#'   element is a vector of ages or \code{NULL} for all ages, or an array
#'   \code{[n_ages x n_fish_fleets]} of 0/1 weights. Default \code{NULL} uses
#'   every age for every fleet. The fleet's compositions are unaffected.
#' @param FishIdx_LikeType Character vector \code{[n_fish_fleets]} giving the
#'   error structure of each fishery index. Options are \code{"lognormal"}
#'   (default, the observation standard errors are on the log scale),
#'   \code{"normal"} (arithmetic scale), and \code{"mvn"} (multivariate normal
#'   on the arithmetic scale using a fixed covariance supplied through
#'   \code{FishIdx_Cov}). One-step-ahead residuals are available only for
#'   lognormal fleets. A fleet's population-specific index stream follows the
#'   same choice for \code{"lognormal"} and \code{"normal"}, but stays
#'   lognormal under \code{"mvn"}, whose covariance describes the regional
#'   series only.
#' @param FishIdx_Cov List with one element per fishery fleet holding the fixed
#'   covariance matrix for fleets using \code{"mvn"}, and \code{NULL}
#'   otherwise. Each matrix must be square with one row per observation the
#'   fleet fits, ordered as the observations appear when scanning that fleet's
#'   \code{UseFishIdx} slice in array order.
#' @param fish_idx_type Character vector of length \code{n_fish_fleets} specifying
#'   the index type for each fleet. \code{"biom"} = biomass; \code{"abd"} =
#'   abundance; \code{"none"} = no index for this fleet.
#' @param t_fish Array \code{[n_regions x n_seas x n_fish_fleets]} giving the
#'   fishery index timing: the fraction of the season elapsed when each index
#'   is observed. Numbers at age are decayed by \code{exp(-t_fish * ZAA)}
#'   before the index is formed, the same convention \code{t_srv} uses for
#'   surveys. Defaults to \code{0} (start of season), which is what the model
#'   did before this argument existed; set \code{0.5} for a mid-season index.
#' @param FishLenComps_sel Character vector \code{[n_fish_fleets]}, whether a
#'   length-based selectivity is applied before or after the fish are spread
#'   over lengths. \code{"age"} (default) selects at age and spreads the catch
#'   afterwards, so every fish of an age is equally catchable and the length
#'   composition within an age is just the key's. \code{"length"} spreads the
#'   fish at each age over the key first and selects them length by length, so
#'   the long fish of an age are taken more often. The key is the fleet's own,
#'   at \code{t_fish}. Requires length-based fishery selectivity. Use
#'   \code{"length"} when selectivity is length based and the length
#'   compositions are what inform it. The two give different expected
#'   compositions, not two roundings of the same one.
#' @param fish_waa_selected Integer vector \code{[n_fish_fleets]} (0/1). With
#'   weight at age derived from growth and length-based selectivity, \code{1}
#'   makes the fleet's catch biomass use the mean weight of the fish it takes at
#'   each age, \eqn{\sum_l P(l \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)},
#'   instead of the population mean weight at that age. Use it when the gear
#'   selects strongly within an age. With flat or age-based selectivity the two
#'   are the same.
#' @param ObsFishIdx_pop Observed population-specific fishery index array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#' @param ObsFishIdx_pop_SE Lognormal standard errors for \code{ObsFishIdx_pop},
#'   same dimensions \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#' @param UseFishIdx_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}. \code{1} =
#'   include population-specific index in likelihood; \code{0} = exclude.
#'   Default: all zeros.
#' @param ObsFishAgeComps Observed fishery age composition array
#'   \code{[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   Values may be raw counts or proportions; if proportions, supply
#'   \code{ISS_FishAgeComps} explicitly.
#' @param UseFishAgeComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit age compositions; \code{0} = exclude.
#' @param ISS_FishAgeComps Input sample size array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), computed automatically by summing
#'   \code{ObsFishAgeComps} within each year-fleet-season-region cell
#'   according to \code{FishAgeComps_Type}.
#' @param ObsFishLenComps Observed fishery length composition array
#'   \code{[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'   Only required when \code{input_list$data$fit_lengths == 1}.
#' @param UseFishLenComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit length compositions; \code{0} = exclude.
#' @param ISS_FishLenComps Input sample size array for length compositions
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), derived automatically from \code{ObsFishLenComps}.
#' @param FishAgeComps_LikeType Character vector of length \code{n_fish_fleets}
#'   specifying the likelihood for fishery age compositions. Options:
#'   \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#'   \code{"iid-Logistic-Normal"}, \code{"1d-Logistic-Normal"},
#'   \code{"2d-Logistic-Normal"}, \code{"none"}.
#' @param FishLenComps_LikeType Same as \code{FishAgeComps_LikeType} but for
#'   length compositions.
#' @param FishAgeComps_Type Character vector defining the age composition
#'   structure (aggregation level) for each fleet and time period. Each element
#'   must follow the format \code{"<type>_Year_<start>-<end>_Fleet_<f>"} or
#'   \code{"<type>_Year_<start>-terminal_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"agg"}}{Aggregated across regions and sexes
#'       (incompatible with \code{"2d-Logistic-Normal"}).}
#'     \item{\code{"spltRspltS"}}{Split by region and sex.}
#'     \item{\code{"spltRjntS"}}{Split by region, summed jointly across sexes.}
#'     \item{\code{"none"}}{No composition data for this fleet and period.}
#'   }
#'   Example: \code{c("spltRjntS_Year_1-10_Fleet_1", "agg_Year_11-terminal_Fleet_1")}.
#' @param FishLenComps_Type Same format and options as \code{FishAgeComps_Type}
#'   but applied to length compositions.
#' @param ObsFish_caal Observed conditional age-at-length array
#'   \code{[n_regions x n_years x n_seas x n_lens x n_ages x n_sexes x
#'   n_fish_fleets]}. A CAAL observation is the age composition of the fish aged
#'   from one length bin, so the age margin of each length row is what gets fit.
#'   \code{NULL} (default) for a model with no CAAL data.
#' @param UseFish_caal Use flags \code{[n_regions x n_years x n_seas x n_lens x
#'   n_fish_fleets]}. Length bins with no aged fish carry a zero and are skipped.
#' @param ISS_Fish_caal Input sample sizes \code{[n_regions x n_years x n_seas x
#'   n_lens x n_sexes x n_fish_fleets]}, the number aged within each length bin
#'   rather than the number measured. Summed from \code{ObsFish_caal} when
#'   \code{NULL}.
#' @param Fish_caal_LikeType Character vector of length \code{n_fish_fleets}.
#'   One of \code{"none"}, \code{"Multinomial"} or
#'   \code{"Dirichlet-Multinomial"}. The logistic-normal families are not
#'   available for CAAL, since a single length bin's age sample is small and
#'   mostly zeros, which the additive log-ratio transform cannot handle.
#' @param Fish_caal_Type Composition type specification, using the same
#'   \code{"CompType_Year_x-y_Fleet_z"} convention as the marginal compositions.
#' @param FishAgeComps_bins Which age bins each fishery fleet's age composition
#'   is fitted over. Supply a list with one element per fleet, each a vector of
#'   age indices or \code{NULL} for all ages, or an
#'   \code{[n_ages x n_fish_fleets]} array of 0/1 weights. Both observed and
#'   expected compositions are restricted to the named bins and renormalized
#'   within them, so excluded bins are left out of the likelihood rather than
#'   being forced to be explained; this is how a fleet that only ages part of
#'   its age range is fitted. Indices refer to observed bins, that is after any
#'   ageing error has mapped model ages onto observed ones. Every fleet must
#'   retain at least one bin. Default \code{NULL}, which fits all ages for all
#'   fleets.
#' @param ObsFishAgeComps_pop Observed population-specific fishery age
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   Required when any element of \code{UseFishAgeComps_pop} is \code{1}.
#' @param UseFishAgeComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit population-specific age compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_FishAgeComps_pop Input sample size array for population-specific
#'   age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), computed automatically by summing
#'   \code{ObsFishAgeComps_pop} within each population-year-fleet-season-region
#'   cell according to \code{FishAgeComps_pop_Type}.
#' @param ObsFishLenComps_pop Observed population-specific fishery length
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'   Required when \code{input_list$data$fit_lengths == 1} and any element of
#'   \code{UseFishLenComps_pop} is \code{1}.
#' @param UseFishLenComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit population-specific length compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_FishLenComps_pop Input sample size array for population-specific
#'   length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), derived automatically from
#'   \code{ObsFishLenComps_pop}.
#' @param FishAgeComps_pop_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying the likelihood for population-specific
#'   fishery age compositions. Same options as \code{FishAgeComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param FishLenComps_pop_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying the likelihood for population-specific
#'   fishery length compositions. Same options as \code{FishLenComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param FishAgeComps_pop_Type Character vector defining the composition
#'   structure for population-specific age compositions. Same format and options
#'   as \code{FishAgeComps_Type}. Default: \code{"none"} for all fleets across
#'   all years.
#' @param FishLenComps_pop_Type Character vector defining the composition
#'   structure for population-specific length compositions. Same format and
#'   options as \code{FishLenComps_Type}. Default: \code{"none"} for all fleets
#'   across all years.
#' @param ... Optional starting value overrides for overdispersion and
#'   correlation parameters.
#' @param ObsFishAgeComps_discard Observed fishery age composition from discards
#'   \code{[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   Structure must match \code{ObsFishAgeComps}.
#'
#' @param UseFishAgeComps_discard Binary indicator array for discard age compositions
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = include discard age compositions in likelihood; \code{0} = exclude.
#'
#' @param ISS_FishAgeComps_discard Input sample size array for discard age compositions
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, derived automatically from \code{ObsFishAgeComps_discard}
#'   using \code{FishAgeComps_discard_Type}.
#'
#' @param ObsFishLenComps_discard Observed fishery length composition from discards
#'   \code{[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'   Required if \code{input_list$data$fit_lengths == 1}.
#'
#' @param UseFishLenComps_discard Binary indicator array for discard length compositions
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = include discard length compositions in likelihood; \code{0} = exclude.
#'
#' @param ISS_FishLenComps_discard Input sample size array for discard length compositions
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, derived automatically from \code{ObsFishLenComps_discard}.
#'
#' @param FishAgeComps_discard_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying likelihood type for discard age compositions.
#'   Options:
#'   \describe{
#'     \item{\code{"Multinomial"}}{Standard multinomial likelihood}
#'     \item{\code{"Dirichlet-Multinomial"}}{Overdispersed multinomial}
#'     \item{\code{"iid-Logistic-Normal"}}{Independent logistic-normal}
#'     \item{\code{"1d-Logistic-Normal"}}{1D correlated logistic-normal}
#'     \item{\code{"2d-Logistic-Normal"}}{2D correlated logistic-normal}
#'     \item{\code{"none"}}{No discard age composition likelihood}
#'   }
#'
#' @param FishLenComps_discard_LikeType Same specification as
#'   \code{FishAgeComps_discard_LikeType}, but for discard length compositions.
#'
#' @param FishAgeComps_discard_Type Character vector defining discard age composition
#'   structure by fleet and year block.
#'   Format:
#'   \code{"<type>_Year_<start>-<end>_Fleet_<f>"} or
#'   \code{"<type>_Year_<start>-terminal_Fleet_<f>"}.
#'   Valid types:
#'   \describe{
#'     \item{\code{"agg"}}{Aggregated across regions and sexes}
#'     \item{\code{"spltRspltS"}}{Split by region and sex}
#'     \item{\code{"spltRjntS"}}{Split by region, joint across sexes}
#'     \item{\code{"none"}}{No discard age composition}
#'   }
#'
#' @param FishLenComps_discard_Type Same format and options as
#'   \code{FishAgeComps_discard_Type}, applied to discard length compositions.
#'
#' @param ObsFishAgeComps_discard_pop Observed population-specific discard age
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'
#' @param UseFishAgeComps_discard_pop Binary indicator array for population-specific
#'   discard age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'
#' @param ISS_FishAgeComps_discard_pop Input sample size array for population-specific
#'   discard age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, computed from \code{ObsFishAgeComps_discard_pop}.
#'
#' @param ObsFishLenComps_discard_pop Observed population-specific discard length
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'
#' @param UseFishLenComps_discard_pop Binary indicator array for population-specific
#'   discard length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'
#' @param ISS_FishLenComps_discard_pop Input sample size array for population-specific
#'   discard length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, derived from \code{ObsFishLenComps_discard_pop}.
#'
#' @param FishAgeComps_discard_pop_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying likelihood type for population-specific
#'   discard age compositions. Same options as \code{FishAgeComps_discard_LikeType}.
#'
#' @param FishLenComps_discard_pop_LikeType Same as above but for discard length compositions.
#'
#' @param FishAgeComps_discard_pop_Type Character vector defining structure for
#'   population-specific discard age compositions. Same format as
#'   \code{FishAgeComps_discard_Type}.
#'
#' @param FishLenComps_discard_pop_Type Character vector defining structure for
#'   population-specific discard length compositions. Same format as
#'   \code{FishLenComps_discard_Type}.
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated with all fishery index and composition fields, including
#'   pooled and population-specific observed arrays, computed or supplied ISS
#'   arrays, integer-coded likelihood and composition type matrices,
#'   overdispersion parameters, and their factor maps.
#'
#' @export Setup_Mod_FishIdx_and_Comps
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_FishIdx_and_Comps <- function(input_list,
                                        ObsFishIdx,
                                        ObsFishIdx_SE,
                                        ObsFishIdx_pop = NULL,
                                        ObsFishIdx_pop_SE = NULL,
                                        UseFishIdx_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        fish_idx_type,
                                        t_fish = array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        FishLenComps_sel = rep("age", input_list$data$n_fish_fleets),
                                        fish_waa_selected = rep(0, input_list$data$n_fish_fleets),
                                        UseFishIdx,

                                        # Retained Compositions
                                        ObsFishAgeComps,
                                        UseFishAgeComps,
                                        ISS_FishAgeComps = NULL,
                                        ObsFishLenComps,
                                        UseFishLenComps,
                                        ISS_FishLenComps = NULL,
                                        FishAgeComps_LikeType,
                                        FishLenComps_LikeType,
                                        FishAgeComps_Type,
                                        FishLenComps_Type,
                                        ObsFishAgeComps_pop = NULL,
                                        UseFishAgeComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishAgeComps_pop = NULL,
                                        ObsFishLenComps_pop = NULL,
                                        UseFishLenComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishLenComps_pop = NULL,
                                        FishAgeComps_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishLenComps_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishAgeComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        FishLenComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        fish_idx_ages = NULL,
                                        FishAgeComps_bins = NULL,
                                        FishIdx_LikeType = rep("lognormal", input_list$data$n_fish_fleets),
                                        FishIdx_Cov = NULL,

                                        # Conditional Age-at-Length
                                        ObsFish_caal = NULL,
                                        UseFish_caal = NULL,
                                        ISS_Fish_caal = NULL,
                                        Fish_caal_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        Fish_caal_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),

                                        # Discard Compositions (forwarded to Setup_Mod_Discard_Comps)
                                        ObsFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                        UseFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishAgeComps_discard = NULL,
                                        ObsFishLenComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, obs_len_bins(input_list), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                        UseFishLenComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishLenComps_discard = NULL,
                                        FishAgeComps_discard_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishLenComps_discard_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishAgeComps_discard_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        FishLenComps_discard_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        ObsFishAgeComps_discard_pop = NULL,
                                        UseFishAgeComps_discard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishAgeComps_discard_pop = NULL,
                                        ObsFishLenComps_discard_pop = NULL,
                                        UseFishLenComps_discard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishLenComps_discard_pop = NULL,
                                        FishAgeComps_discard_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishLenComps_discard_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishAgeComps_discard_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        FishLenComps_discard_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        ...
                                        ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_FishIdx_and_Comps <- mget(names(formals()))[-1]

  # Input Validation ---------------------------------------------------------

  # Fishery indices
  check_data_dimensions(ObsFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx')
  check_data_dimensions(ObsFishIdx_SE, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_SE')
  check_data_dimensions(UseFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishIdx')

  if(any(UseFishIdx_pop == 1)) {
    check_data_dimensions(ObsFishIdx_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_pop')
    check_data_dimensions(ObsFishIdx_pop_SE, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_pop_SE')
    check_data_dimensions(UseFishIdx_pop, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishIdx_pop')
  }

  if(!all(fish_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for fish_idx_type. Should be either abd, biom, or none")

  # Fishery compositions
  check_data_dimensions(ObsFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps')
  check_data_dimensions(UseFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps')
  check_data_dimensions(UseFishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsFishLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = obs_len_bins(input_list), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps')
  if(!is.null(ISS_FishAgeComps)) check_data_dimensions(ISS_FishAgeComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps')
  if(!is.null(ISS_FishLenComps)) check_data_dimensions(ISS_FishLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps')
  check_data_dimensions(FishAgeComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_LikeType')
  check_data_dimensions(FishLenComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_LikeType')
  if(!all(FishAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

   # Fishery compositions (population-specific)
  if(any(UseFishAgeComps_pop == 1)) check_data_dimensions(ObsFishAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps_pop')
  check_data_dimensions(UseFishAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps_pop')
  check_data_dimensions(UseFishLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps_pop')
  if(input_list$data$fit_lengths == 1 && any(UseFishLenComps_pop == 1)) check_data_dimensions(ObsFishLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = obs_len_bins(input_list), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps_pop')
  if(!is.null(ISS_FishAgeComps_pop)) check_data_dimensions(ISS_FishAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps_pop')
  if(!is.null(ISS_FishLenComps_pop)) check_data_dimensions(ISS_FishLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps_pop')
  check_data_dimensions(FishAgeComps_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_pop_LikeType')
  check_data_dimensions(FishLenComps_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_pop_LikeType')
  if(!all(FishAgeComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # checking to make sure defaults are not applied
  if(any(UseFishAgeComps_pop == 1)) {
    if(is.null(ObsFishAgeComps_pop)) stop("ObsFishAgeComps_pop is NULL, but UseFishAgeComps_pop contains 1s!")
    if(any(str_detect(FishAgeComps_pop_LikeType, "none"))) warning("FishAgeComps_pop_LikeType has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(FishAgeComps_pop_Type, "none"))) warning("FishAgeComps_pop_Type has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
  }

  if(any(UseFishLenComps_pop == 1)) {
    if(is.null(ObsFishLenComps_pop)) stop("ObsFishLenComps_pop is NULL, but UseFishLenComps_pop contains 1s!")
    if(any(str_detect(FishLenComps_pop_LikeType, "none"))) warning("FishLenComps_pop_LikeType has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(FishLenComps_pop_Type, "none"))) warning("FishLenComps_pop_Type has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
  }


  # Fishery Index Options ---------------------------------------------------

  check_data_dimensions(t_fish, n_regions = input_list$data$n_regions,
                        n_seas = input_list$data$n_seas,
                        n_fish_fleets = input_list$data$n_fish_fleets, what = "t_fish")
  if(any(t_fish < 0 | t_fish > 1)) stop("t_fish must be a fraction of the season in [0, 1]")

  fish_idx_type_vals <- array(NA, dim = c(input_list$data$n_fish_fleets))
  for(f in 1:input_list$data$n_fish_fleets) {
    if(fish_idx_type[f] == 'biom') fish_idx_type_vals[f] <- 1 # biomass
    if(fish_idx_type[f] == 'abd') fish_idx_type_vals[f] <- 0 # abundance
    if(fish_idx_type[f] == 'none') fish_idx_type_vals[f] <- 999 # none
    collect_message(paste("Fishery Index", "for fishery fleet", f, "specified as:" , fish_idx_type[f]))
  } # end f loop


  # Fishery Age Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishage_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_LikeType[f] == 'none') comp_fishage_like_vals <- c(comp_fishage_like_vals, 999)
    if(FishAgeComps_LikeType[f] == "Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 0)
    if(FishAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 1)
    if(FishAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 2)
    if(FishAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 3)
    if(FishAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 4)
    collect_message(paste("Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_like_vals[fleet] == 4) stop("Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_Type_Mat))) stop("FishAgeComps_Type is returning an NA. Did you update the year range of FishAgeComps_Type?")

  # Specifying composition likelihood for population-specific data
  comp_fishage_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_pop_LikeType[f] == 'none') comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 999)
    if(FishAgeComps_pop_LikeType[f] == "Multinomial") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 0)
    if(FishAgeComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 1)
    if(FishAgeComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 2)
    if(FishAgeComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 3)
    if(FishAgeComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 4)
    collect_message(paste("Population Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_pop_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_pop_like_vals[fleet] == 4) stop("Population Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_pop_Type_Mat))) stop("FishAgeComps_pop_Type is returning an NA. Did you update the year range of FishAgeComps_pop_Type?")

  # Fishery Length Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishlen_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_LikeType[f] == 'none') comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 999)
    if(FishLenComps_LikeType[f] == "Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 0)
    if(FishLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 1)
    if(FishLenComps_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 2)
    if(FishLenComps_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 3)
    if(FishLenComps_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 4)
    collect_message(paste("Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # define composition types
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_like_vals[fleet] == 4) stop("Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_Type_Mat))) stop("FishLenComps_Type_Mat is returning an NA. Did you update the year range of FishLenComps_Type_Mat?")


  # Specifying composition likelihood for population-specific data
  comp_fishlen_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_pop_LikeType[f] == 'none') comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 999)
    if(FishLenComps_pop_LikeType[f] == "Multinomial") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 0)
    if(FishLenComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 1)
    if(FishLenComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 2)
    if(FishLenComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 3)
    if(FishLenComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 4)
    collect_message(paste("Population Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_pop_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_pop_like_vals[fleet] == 4) stop("Population Len composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_pop_Type_Mat))) stop("FishLenComps_pop_Type is returning an NA. Did you update the year range of FishLenComps_pop_Type?")

  # whether length selectivity is applied at length or through the size-age key.
  # Whether the selectivity is length based is only known once Setup_Mod_Fishsel_and_Q
  # runs, which checks this against it
  if(length(FishLenComps_sel) != input_list$data$n_fish_fleets || !all(FishLenComps_sel %in% c("age", "length"))) stop("FishLenComps_sel must be one of age or length for each fishery fleet")
  fish_len_comp_sel_vals <- as.numeric(FishLenComps_sel == "length")
  for(f in 1:input_list$data$n_fish_fleets) if(FishLenComps_sel[f] == "length") collect_message("Fishery length compositions for fleet ", f, " apply selectivity at length")

  # Fishery Weight at Age Options ---------------------------------------------

  if(length(fish_waa_selected) != input_list$data$n_fish_fleets || !all(fish_waa_selected %in% c(0, 1))) stop("fish_waa_selected must be 0 or 1 for each fishery fleet")
  if(any(fish_waa_selected == 1) && (is.null(input_list$data$derive_waa) || input_list$data$derive_waa != 1)) stop("fish_waa_selected = 1 needs waa_model = 'wt_len' in Setup_Mod_Biologicals")
  for(f in which(fish_waa_selected == 1)) collect_message("Fishery fleet ", f, " takes its catch biomass on the selection-weighted weight at age")

  # ISS Munging -------------------------------------------------------------

  # Fishery Ages
  if(is.null(ISS_FishAgeComps)) {
    collect_message("No ISS is specified for FishAgeComps. ISS weighting is calculated by summing up values from ObsFishAgeComps each year")
    ISS_FishAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0) or joint across sexes
          if(FishAgeComps_Type_Mat[y,f] == 0) ISS_FishAgeComps[1,y,seas,1,f] <- sum(ObsFishAgeComps[,y,seas,,,f])
          # if split by region and sex
          if(FishAgeComps_Type_Mat[y,f] == 1) ISS_FishAgeComps[,y,seas,,f] <- apply(ObsFishAgeComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishAgeComps_Type_Mat[y,f] == 2) ISS_FishAgeComps[,y,seas,1,f] <- apply(ObsFishAgeComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps)) {
    collect_message("No ISS is specified for FishLenComps. ISS weighting is calculated by summing up values from ObsFishLenComps each year")
    ISS_FishLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(FishLenComps_Type_Mat[y,f] == 0) ISS_FishLenComps[1,y,seas,1,f] <- sum(ObsFishLenComps[,y,seas,,,f])
          # if split by region and sex
          if(FishLenComps_Type_Mat[y,f] == 1) ISS_FishLenComps[,y,seas,,f] <- apply(ObsFishLenComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishLenComps_Type_Mat[y,f] == 2) ISS_FishLenComps[,y,seas,1,f] <- apply(ObsFishLenComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Ages
  if(is.null(ISS_FishAgeComps_pop)) {
    collect_message("No ISS is specified for pop_FishAgeComps. ISS weighting is calculated by summing up values from ObsFishAgeComps_pop each year")
    ISS_FishAgeComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0) or joint across sexes
            if(FishAgeComps_pop_Type_Mat[y,f] == 0) ISS_FishAgeComps_pop[p,1,y,seas,1,f] <- sum(ObsFishAgeComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishAgeComps_pop_Type_Mat[y,f] == 1) ISS_FishAgeComps_pop[p,,y,seas,,f] <- apply(ObsFishAgeComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishAgeComps_pop_Type_Mat[y,f] == 2) ISS_FishAgeComps_pop[p,,y,seas,1,f] <- apply(ObsFishAgeComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps_pop)) {
    collect_message("No ISS is specified for pop_FishLenComps. ISS weighting is calculated by summing up values from ObsFishLenComps_pop each year")
    ISS_FishLenComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0)
            if(FishLenComps_pop_Type_Mat[y,f] == 0) ISS_FishLenComps_pop[p,1,y,seas,1,f] <- sum(ObsFishLenComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishLenComps_pop_Type_Mat[y,f] == 1) ISS_FishLenComps_pop[p,,y,seas,,f] <- apply(ObsFishLenComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishLenComps_pop_Type_Mat[y,f] == 2) ISS_FishLenComps_pop[p,,y,seas,1,f] <- apply(ObsFishLenComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_FishAgeComps <- ISS_FishAgeComps
  input_list$data$ISS_FishLenComps <- ISS_FishLenComps
  input_list$data$ISS_FishAgeComps_pop <- ISS_FishAgeComps_pop
  input_list$data$ISS_FishLenComps_pop <- ISS_FishLenComps_pop
  input_list$data$ObsFishIdx <- ObsFishIdx
  input_list$data$ObsFishIdx_SE <- ObsFishIdx_SE
  input_list$data$UseFishIdx <- UseFishIdx
  input_list$data$ObsFishIdx_pop <- ObsFishIdx_pop
  input_list$data$ObsFishIdx_pop_SE <- ObsFishIdx_pop_SE
  input_list$data$UseFishIdx_pop <- UseFishIdx_pop
  input_list$data$fish_idx_type <- fish_idx_type_vals

  # Index age selection and error structure ---------------------------------
  if(!all(FishIdx_LikeType %in% c("lognormal", "normal", "mvn"))) stop("Invalid specification for FishIdx_LikeType. Should be lognormal, normal, or mvn")
  if(length(FishIdx_LikeType) != input_list$data$n_fish_fleets) stop("FishIdx_LikeType is not length n_fish_fleets")

  fish_idx_like_vals <- convert_to_numeric(FishIdx_LikeType, list(lognormal = 0, normal = 1, mvn = 2))
  fish_idx_ages_arr <- parse_idx_ages(fish_idx_ages, length(input_list$data$ages), input_list$data$n_fish_fleets, "fish_idx_ages")
  fish_idx_cov_parsed <- parse_idx_cov(FishIdx_Cov, fish_idx_like_vals, UseFishIdx, input_list$data$n_fish_fleets, "FishIdx_Cov")

  for(f in 1:input_list$data$n_fish_fleets) {
    collect_message(paste("Fishery Index likelihood for fishery fleet", f, "specified as:", FishIdx_LikeType[f]))
    if(sum(fish_idx_ages_arr[,f]) != length(input_list$data$ages)) {
      collect_message(paste("Fishery Index for fishery fleet", f, "is restricted to ages:", paste(which(fish_idx_ages_arr[,f] == 1), collapse = ", ")))
    }
  } # end f loop

  input_list$data$FishIdx_LikeType <- fish_idx_like_vals
  input_list$data$fish_idx_ages <- fish_idx_ages_arr
  input_list$data$FishAgeComps_bins <- parse_idx_ages(FishAgeComps_bins, length(input_list$data$ages), input_list$data$n_fish_fleets, "FishAgeComps_bins")
  input_list$data$FishIdx_Cov <- fish_idx_cov_parsed
  input_list$data$t_fish <- t_fish
  input_list$data$fish_len_comp_sel <- fish_len_comp_sel_vals
  input_list$data$fish_waa_selected <- fish_waa_selected
  input_list$data$ObsFishAgeComps <- ObsFishAgeComps
  input_list$data$UseFishAgeComps <- UseFishAgeComps
  input_list$data$ObsFishLenComps <- ObsFishLenComps
  input_list$data$UseFishLenComps <- UseFishLenComps
  input_list$data$ObsFishAgeComps_pop <- ObsFishAgeComps_pop
  input_list$data$UseFishAgeComps_pop <- UseFishAgeComps_pop
  input_list$data$ObsFishLenComps_pop <- ObsFishLenComps_pop
  input_list$data$UseFishLenComps_pop <- UseFishLenComps_pop
  input_list$data$FishAgeComps_LikeType <- comp_fishage_like_vals
  input_list$data$FishLenComps_LikeType <- comp_fishlen_like_vals
  input_list$data$FishAgeComps_pop_LikeType <- comp_fishage_pop_like_vals
  input_list$data$FishLenComps_pop_LikeType <- comp_fishlen_pop_like_vals
  input_list$data$FishAgeComps_Type <- FishAgeComps_Type_Mat
  input_list$data$FishLenComps_Type <- FishLenComps_Type_Mat
  input_list$data$FishAgeComps_pop_Type <- FishAgeComps_pop_Type_Mat
  input_list$data$FishLenComps_pop_Type <- FishLenComps_pop_Type_Mat

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the fishery age comps
  if("ln_FishAge_theta" %in% names(starting_values)) input_list$par$ln_FishAge_theta <- starting_values$ln_FishAge_theta
  else input_list$par$ln_FishAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery age comps
  if("FishAge_corr_pars" %in% names(starting_values)) input_list$par$FishAge_corr_pars <- starting_values$FishAge_corr_pars
  else input_list$par$FishAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishAge_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_theta_agg <- starting_values$ln_FishAge_theta_agg
  else input_list$par$ln_FishAge_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  # aggregated correlation parameters
  if("FishAge_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_corr_pars_agg <- starting_values$FishAge_corr_pars_agg
  else input_list$par$FishAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for fishery length comps
  if("ln_FishLen_theta" %in% names(starting_values)) input_list$par$ln_FishLen_theta <- starting_values$ln_FishLen_theta
  else input_list$par$ln_FishLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery length comps
  if("FishLen_corr_pars" %in% names(starting_values)) input_list$par$FishLen_corr_pars <- starting_values$FishLen_corr_pars
  else input_list$par$FishLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishLen_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_theta_agg <- starting_values$ln_FishLen_theta_agg
  else input_list$par$ln_FishLen_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  if("FishLen_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_corr_pars_agg <- starting_values$FishLen_corr_pars_agg
  else input_list$par$FishLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for the population fishery age comps
  if("ln_FishAge_pop_theta" %in% names(starting_values)) input_list$par$ln_FishAge_pop_theta <- starting_values$ln_FishAge_pop_theta
  else input_list$par$ln_FishAge_pop_theta <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery age comps
  if("FishAge_pop_corr_pars" %in% names(starting_values)) input_list$par$FishAge_pop_corr_pars <- starting_values$FishAge_pop_corr_pars
  else input_list$par$FishAge_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishAge_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_pop_theta_agg <- starting_values$ln_FishAge_pop_theta_agg
  else input_list$par$ln_FishAge_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # aggregated population correlation parameters
  if("FishAge_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_pop_corr_pars_agg <- starting_values$FishAge_pop_corr_pars_agg
  else input_list$par$FishAge_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # Dispersion parameters for population fishery length comps
  if("ln_FishLen_pop_theta" %in% names(starting_values)) input_list$par$ln_FishLen_pop_theta <- starting_values$ln_FishLen_pop_theta
  else input_list$par$ln_FishLen_pop_theta <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery length comps
  if("FishLen_pop_corr_pars" %in% names(starting_values)) input_list$par$FishLen_pop_corr_pars <- starting_values$FishLen_pop_corr_pars
  else input_list$par$FishLen_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishLen_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_pop_theta_agg <- starting_values$ln_FishLen_pop_theta_agg
  else input_list$par$ln_FishLen_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  if("FishLen_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_pop_corr_pars_agg <- starting_values$FishLen_pop_corr_pars_agg
  else input_list$par$FishLen_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishAge")
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishLen")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishAge")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishLen")

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishAge", has_pop = TRUE)
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "FishLen", has_pop = TRUE)
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishAge", has_pop = TRUE)
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "FishLen", has_pop = TRUE)

  # Discard Compositions (forwarded to Setup_Mod_Discard_Comps) ---------------
  input_list <- Setup_Mod_Discard_Comps(
    input_list,
    ObsFishAgeComps_discard        = ObsFishAgeComps_discard,
    UseFishAgeComps_discard        = UseFishAgeComps_discard,
    ISS_FishAgeComps_discard       = ISS_FishAgeComps_discard,
    ObsFishLenComps_discard        = ObsFishLenComps_discard,
    UseFishLenComps_discard        = UseFishLenComps_discard,
    ISS_FishLenComps_discard       = ISS_FishLenComps_discard,
    FishAgeComps_discard_LikeType  = FishAgeComps_discard_LikeType,
    FishLenComps_discard_LikeType  = FishLenComps_discard_LikeType,
    FishAgeComps_discard_Type      = FishAgeComps_discard_Type,
    FishLenComps_discard_Type      = FishLenComps_discard_Type,
    ObsFishAgeComps_discard_pop    = ObsFishAgeComps_discard_pop,
    UseFishAgeComps_discard_pop    = UseFishAgeComps_discard_pop,
    ISS_FishAgeComps_discard_pop   = ISS_FishAgeComps_discard_pop,
    ObsFishLenComps_discard_pop    = ObsFishLenComps_discard_pop,
    UseFishLenComps_discard_pop    = UseFishLenComps_discard_pop,
    ISS_FishLenComps_discard_pop   = ISS_FishLenComps_discard_pop,
    FishAgeComps_discard_pop_LikeType = FishAgeComps_discard_pop_LikeType,
    FishLenComps_discard_pop_LikeType = FishLenComps_discard_pop_LikeType,
    FishAgeComps_discard_pop_Type     = FishAgeComps_discard_pop_Type,
    FishLenComps_discard_pop_Type     = FishLenComps_discard_pop_Type,
    ...
  )

  # Conditional Age-at-Length --------------------------------------------------
  # A CAAL observation is an age composition within a length bin, so it belongs
  # with the other compositions for this fleet rather than in a block of its own.
  input_list <- setup_caal_stream(
    input_list,
    ObsCAAL = ObsFish_caal, UseCAAL = UseFish_caal,
    ISS_CAAL = ISS_Fish_caal,
    CAAL_LikeType = Fish_caal_LikeType, CAAL_Type = Fish_caal_Type,
    fleet_type = "Fish"
  )

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
