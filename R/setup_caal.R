# Stage 1 of 3: setup
#
# Conditional age-at-length data streams. These live inside the fishery and survey
# composition blocks rather than in a block of their own, since a CAAL observation
# is an age composition and belongs with the other compositions for that fleet.
# Setup_Mod_FishIdx_and_Comps and Setup_Mod_SrvIdx_and_Comps both call
# setup_caal_stream, which validates the arrays, translates the likelihood and
# composition type specifications, and builds the overdispersion parameters and
# their maps.

#' Translate a CAAL composition type specification into a year by fleet matrix
#'
#' Uses the same \code{"CompType_Year_x-y_Fleet_z"} string convention as the
#' marginal compositions, so a CAAL type is specified exactly the way an age or
#' length composition type is. \code{"terminal"} in place of the upper year runs
#' the block to the last model year.
#'
#' @param CAAL_Type Character vector of composition type specifications.
#' @param n_yrs Number of model years.
#' @param n_fleets Number of fleets.
#' @param what Name used in error messages.
#'
#' @return Integer matrix \code{[year x fleet]} of composition type codes
#'   (0 aggregated, 1 split by region and sex, 2 joint by sex and split by
#'   region, 999 none).
#'
#' @keywords internal
parse_caal_type <- function(CAAL_Type, n_yrs, n_fleets, what) {

  type_mat <- array(999, dim = c(n_yrs, n_fleets))

  for(i in seq_along(CAAL_Type)) {

    tmp <- CAAL_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1]
    fleet <- as.numeric(tmp_vec[5])

    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", "none"))
      stop(paste(what, "not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none"))
    if(!fleet %in% 1:n_fleets)
      stop(paste("Invalid fleet specified for", what, ". This needs to be specified as CompType_Year_x-y_Fleet_x"))

    if(!stringr::str_detect(tmp, "terminal")) {
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2]
    } else {
      years <- as.numeric(unlist(strsplit(tmp_vec[3], "-"))[1]):n_yrs
    }

    if(comps_type_tmp == "agg") comps_type_val <- 0
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    type_mat[years, fleet] <- comps_type_val

  } # end i loop

  return(type_mat)
}


#' Set up a conditional age-at-length data stream
#'
#' Shared by \code{\link{Setup_Mod_FishIdx_and_Comps}} and
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}}. Validates the CAAL arrays, converts
#' the likelihood and composition type specifications, fills in a default input
#' sample size when none is supplied, and creates the Dirichlet-multinomial
#' overdispersion parameters together with their maps. The likelihood weights
#' (\code{Wt_Fish_caal}, \code{Wt_Srv_caal}) belong to
#' \code{\link{Setup_Mod_Weighting}} with every other weight.
#'
#' Only the multinomial and Dirichlet-multinomial families are available. A CAAL
#' row is the age composition of the otoliths taken from one length bin, so it is
#' usually a small sample that is mostly zeros, and the logistic-normal forms need
#' an additive log-ratio transform that such a row cannot support.
#'
#' Supplying any CAAL data switches on \code{do_caal}, since the likelihood reads
#' the joint arrays at length and age that flag builds.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#' @param ObsCAAL Observed CAAL array \code{[n_regions x n_years x n_seas x
#'   n_lens x n_ages x n_sexes x n_fleets]}, or \code{NULL} for none.
#' @param UseCAAL Use flags \code{[n_regions x n_years x n_seas x n_lens x
#'   n_fleets]}. A length bin with no aged fish carries a zero and is skipped.
#' @param ISS_CAAL Input sample sizes \code{[n_regions x n_years x n_seas x
#'   n_lens x n_sexes x n_fleets]}, the number aged within each length bin. When
#'   \code{NULL} it is summed from \code{ObsCAAL}.
#' @param CAAL_LikeType Character vector of length \code{n_fleets}. One of
#'   \code{"none"}, \code{"Multinomial"} or \code{"Dirichlet-Multinomial"}.
#' @param CAAL_Type Character vector of composition type specifications, using
#'   the same \code{"CompType_Year_x-y_Fleet_z"} convention as the marginal
#'   compositions.
#' @param fleet_type Character, either \code{"Fish"} or \code{"Srv"}.
#'
#' @return The updated \code{input_list}.
#'
#' @keywords internal
setup_caal_stream <- function(input_list, ObsCAAL, UseCAAL, ISS_CAAL,
                              CAAL_LikeType, CAAL_Type, fleet_type) {

  n_regions <- input_list$data$n_regions
  n_yrs <- length(input_list$data$years)
  n_seas <- input_list$data$n_seas
  n_lens <- length(input_list$data$lens)
  n_ages <- length(input_list$data$ages)
  n_sexes <- input_list$data$n_sexes
  n_fleets <- if(fleet_type == "Fish") input_list$data$n_fish_fleets else input_list$data$n_srv_fleets

  stub <- paste0(fleet_type, "_caal")
  obs_nm <- paste0("Obs", stub)
  use_nm <- paste0("Use", stub)
  iss_nm <- paste0("ISS_", stub)
  type_nm <- paste0(stub, "_Type")
  like_nm <- paste0(stub, "_LikeType")
  theta_nm <- paste0("ln_", stub, "_theta")
  theta_agg_nm <- paste0("ln_", stub, "_theta_agg")

  obs_dim <- c(n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_fleets)
  iss_dim <- c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fleets)
  use_dim <- c(n_regions, n_yrs, n_seas, n_lens, n_fleets)

  # Defaults for a model that carries no CAAL data at all
  if(is.null(UseCAAL)) UseCAAL <- array(0, dim = use_dim)
  if(is.null(ObsCAAL)) ObsCAAL <- array(0, dim = obs_dim)
  if(is.null(CAAL_LikeType)) CAAL_LikeType <- rep("none", n_fleets)
  if(is.null(CAAL_Type)) CAAL_Type <- paste0("none_Year_1-terminal_Fleet_", 1:n_fleets)

  # Input Validation --------------------------------------------------------
  # Checked here rather than through check_data_dimensions, which dispatches on a
  # fixed list of names and has no CAAL cases
  dim_msg <- function(x, want, what, labels) {
    if(length(dim(x)) != length(want) || !all(dim(x) == want))
      stop(paste0("Dimensions of ", what, " are not correct. Should be ", paste(labels, collapse = ", "),
                  " (", paste(want, collapse = " x "), "), but are (", paste(dim(x), collapse = " x "), ")"))
  }

  if(any(UseCAAL == 1)) {
    if(input_list$data$fit_lengths != 1)
      stop(paste0(use_nm, " has observations, but fit_lengths is 0. Conditional age-at-length needs the size-age transition matrix, so set fit_lengths = 1 in Setup_Mod_Biologicals"))
    # The observed age bins need not be the model ages; the ageing error matrix maps
    # one onto the other, as it does for the marginal age compositions
    n_obs_ages <- dim(ObsCAAL)[5]
    if(is.null(n_obs_ages) || n_obs_ages < 1) stop(paste0(obs_nm, " must carry at least one observed age bin"))
    dim_msg(ObsCAAL, replace(obs_dim, 5, n_obs_ages), obs_nm, c("n_regions", "n_years", "n_seas", "n_lens", "n_obs_ages", "n_sexes", "n_fleets"))
  }
  dim_msg(UseCAAL, use_dim, use_nm, c("n_regions", "n_years", "n_seas", "n_lens", "n_fleets"))

  if(length(CAAL_LikeType) != n_fleets) stop(paste("Dimensions of", like_nm, "are not correct. Should be a vector of length n_fleets"))
  if(!all(CAAL_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial")))
    stop(paste("Invalid specification for", like_nm, ". Should be either none, Multinomial, or Dirichlet-Multinomial. The logistic-normal families are not available for conditional age-at-length"))

  # Likelihood and composition types -----------------------------------------
  like_vals <- rep(999, n_fleets)
  for(f in 1:n_fleets) {
    if(CAAL_LikeType[f] == "Multinomial") like_vals[f] <- 0
    if(CAAL_LikeType[f] == "Dirichlet-Multinomial") like_vals[f] <- 1
    collect_message(paste("Conditional Age-at-Length Likelihoods for", fleet_type, "fleet", f, "specified as:", CAAL_LikeType[f]))
  } # end f loop

  type_mat <- parse_caal_type(CAAL_Type, n_yrs = n_yrs, n_fleets = n_fleets, what = type_nm)

  # comp_const_obs isn't final until Setup_Mod_Weighting runs, so its Dirichlet-multinomial
  # sanity check lives there instead.

  # A fleet with no likelihood contributes nothing, so make sure its use flags are off
  for(f in 1:n_fleets) if(like_vals[f] == 999) UseCAAL[,,,,f] <- 0

  # Input sample size and weighting -------------------------------------------
  if(is.null(ISS_CAAL)) {
    collect_message(paste0("No ISS is specified for ", stub, ". ISS is calculated by summing up values from ", obs_nm, " within each length bin"))
    ISS_CAAL <- array(0, dim = iss_dim)
    for(y in 1:n_yrs) {
      for(f in 1:n_fleets) {
        for(seas in 1:n_seas) {
          for(l in 1:n_lens) {
            # aggregated (0) and joint across sexes (2) draw one sample per cell, so the
            # count belongs in the first slot the likelihood reads
            if(type_mat[y,f] == 0) ISS_CAAL[1,y,seas,l,1,f] <- sum(ObsCAAL[,y,seas,l,,,f])
            if(type_mat[y,f] == 2) for(r in 1:n_regions) ISS_CAAL[r,y,seas,l,1,f] <- sum(ObsCAAL[r,y,seas,l,,,f])
            if(type_mat[y,f] == 1) for(r in 1:n_regions) for(s in 1:n_sexes) ISS_CAAL[r,y,seas,l,s,f] <- sum(ObsCAAL[r,y,seas,l,,s,f])
          } # end l loop
        } # end seas loop
      } # end f loop
    } # end y loop
  } else dim_msg(ISS_CAAL, iss_dim, iss_nm, c("n_regions", "n_years", "n_seas", "n_lens", "n_sexes", "n_fleets"))

  # Reconcile the use flags with any bin restriction, so the fitting likelihood
  # and the residual machinery agree on which length bins carry aged fish. The
  # bins array is set by the caller before this runs.
  UseCAAL <- drop_empty_fitted_blocks(ObsCAAL, UseCAAL, input_list$data[[paste0(stub, "_bins")]], 5, stub)

  # Populate Data List ------------------------------------------------------
  input_list$data[[obs_nm]] <- ObsCAAL
  input_list$data[[use_nm]] <- UseCAAL
  input_list$data[[iss_nm]] <- ISS_CAAL
  input_list$data[[type_nm]] <- type_mat
  input_list$data[[like_nm]] <- like_vals

  # The likelihood reads the joint arrays at length and age, so any CAAL data
  # switches the flag that builds them
  if(any(UseCAAL == 1)) input_list$data$do_caal <- 1

  # Populate Parameter List -------------------------------------------------
  # One theta per region and sex is shared across length bins, since the bins come
  # from one length-stratified sample rather than from independent surveys
  input_list$par[[theta_nm]] <- array(0, dim = c(n_regions, n_sexes, n_fleets))
  input_list$par[[theta_agg_nm]] <- array(0, dim = n_fleets)

  # Mapping Options ---------------------------------------------------------
  # only the Dirichlet-multinomial estimates an overdispersion, so every other
  # theta is held
  map_theta <- input_list$par[[theta_nm]]
  map_theta_agg <- input_list$par[[theta_agg_nm]]
  map_theta[] <- NA
  map_theta_agg[] <- NA
  counter <- 1
  counter_agg <- 1

  for(f in 1:n_fleets) {

    comp_type <- unique(type_mat[,f])
    like_type <- like_vals[f]

    if(any(comp_type == 0) && like_type == 1) {
      map_theta_agg[f] <- counter_agg
      counter_agg <- counter_agg + 1
    }

    for(r in 1:n_regions) {

      # a region with no CAAL for this fleet never contributes, so its theta must stay
      # fixed, otherwise it is a free unidentifiable parameter
      region_has_data <- sum(UseCAAL[r,,,,f]) > 0

      for(s in 1:n_sexes) {
        if(any(comp_type == 1) && like_type == 1 && region_has_data) {
          map_theta[r,s,f] <- counter
          counter <- counter + 1
        }
        # joint across sexes reads only the first sex slot
        if(any(comp_type == 2) && like_type == 1 && region_has_data && s == 1) {
          map_theta[r,s,f] <- counter
          counter <- counter + 1
        }
      } # end s loop
    } # end r loop
  } # end f loop

  input_list$map[[theta_nm]] <- factor(map_theta)
  input_list$map[[theta_agg_nm]] <- factor(map_theta_agg)

  return(input_list)
}
