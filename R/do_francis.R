#' Computes Francis weights, which is used internally by do_francis_reweighting
#'
#' @param n_regions Number of regions
#' @param n_sexes Number of sexes
#' @param n_fleets Number of fleets (fishery or survey)
#' @param Use Array from data list that specifies whether to use data that year
#' @param ISS Input sample size array
#' @param Pred_array Predicted values array dimensioned by n_regions, n_years, n_seas, n_bins, n_sexes, n_fleets
#' @param Obs_array Observed values array dimensioned by n_regions, n_years, n_seas, n_bins, n_sexes, n_fleets
#' @param bins Vector of bins used (age or length)
#' @param weights Array of francis weights (NAs) to apply dimensioned by n_regions, n_years, n_seas, n_sexes, n_fleets
#' @param comp_type Matrix of composition structure types dimensioned by year and fleet
#' @param n_years Number of years
#' @param n_seas Number of seasons
#' @param n_bins Number of bins
#'
#' @returns List of values for calculated francis weight, and a dataframe of observed and expected means
#' @keywords internal
get_francis_weights <- function(n_regions,
                                n_sexes,
                                n_fleets,
                                n_years,
                                n_seas,
                                n_bins,
                                Use,
                                ISS,
                                Pred_array,
                                Obs_array,
                                weights,
                                bins,
                                comp_type
) {

# Inverse variance calculations used in Francis reweighting (makes it so that it doesn't take variance from vector of length 1)
safe_inv_var <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) < 2) return(1)
  v <- stats::var(x)
  if(!is.finite(v) || v == 0) return(1)
  1 / v
}

  mean_francis <- data.frame()

  for(f in 1:n_fleets) {

    data_indices <- matrix(nrow=0, ncol=3)  # storage container for data indices

    for(r in 1:n_regions) {
      for(y in 1:n_years) {
        for(seas in 1:n_seas) {
          if(Use[r, y, seas, f] == 1) {
            data_indices <- rbind(data_indices, c(r, y, seas)) # get data indices by regionn and year
          } # end if
        } # end seas loop
      } # end y loop
    } # end r loop

    if(nrow(data_indices) == 0) next  # skip if no usable data

    # get year ranges with data to loop through
    data_yrs <- sort(unique(data_indices[,2]))

    # Set up reweighting vectors
    exp_bar <- array(NA, dim = c(n_regions, length(data_yrs), n_seas, n_sexes),  dimnames = list(NULL, data_yrs, NULL, NULL))
    obs_bar <- array(NA, dim = c(n_regions, length(data_yrs), n_seas, n_sexes), dimnames = list(NULL, data_yrs, NULL, NULL))
    v_y <- array(NA, dim = c(n_regions, length(data_yrs), n_seas, n_sexes))
    w_denom <- array(NA, dim = c(n_regions, length(data_yrs), n_seas, n_sexes))

    for(y in data_yrs) {

      for(seas in 1:n_seas) {

        # Extract out temporary variables
        tmp_iss_obs <- ISS[,y,seas,,f, drop = FALSE] # temporary ISS
        tmp_exp <- Pred_array[,y,seas,,,f, drop = FALSE] # temporary expected values
        tmp_obs <- Obs_array[,y,seas,,,f, drop = FALSE] # temporary observed values
        yr_alt_idx <- which(data_yrs == y) # get indexing to start from 1
        use_regions <- data_indices[which(data_indices[,2] == y & data_indices[,3] == seas),1] # get regions with data

        # If compositions are aggregated across regions and sexes
        if(comp_type[y,f] == 0) {
          exp_bar[1,yr_alt_idx,seas,1] <- sum(bins * as.vector(tmp_exp[1,1,1,,1,1])) # get mean pred comps
          obs_bar[1,yr_alt_idx,seas,1] <- sum(bins * as.vector(tmp_obs[1,1,1,,1,1])) # get mean obs comps
          v_y[1,yr_alt_idx,seas,1] <- sum(bins^2*tmp_exp[1,1,1,,1,1])-exp_bar[1,yr_alt_idx,seas,1]^2 # get variance
          w_denom[1,yr_alt_idx,seas,1] <- (obs_bar[1,yr_alt_idx,seas,1]-exp_bar[1,yr_alt_idx,seas,1])/sqrt(v_y[1,yr_alt_idx,seas,1]/tmp_iss_obs[1,1,1,1,1]) # get weights
        } # end if aggregated

        # If compositions are split by region and sex
        if(comp_type[y,f] == 1) {
          for(r in use_regions) {
            for(s in 1:n_sexes) {
              exp_bar[r,yr_alt_idx,seas,s] <- sum(bins * as.vector(tmp_exp[r,1,1,,s,1])) # get mean pred comps
              obs_bar[r,yr_alt_idx,seas,s] <- sum(bins * as.vector(tmp_obs[r,1,1,,s,1])) # get mean obs comps
              v_y[r,yr_alt_idx,seas,s] <- sum(bins^2*tmp_exp[r,1,1,,s,1])-exp_bar[r,yr_alt_idx,seas,s]^2 # get variance
              w_denom[r,yr_alt_idx,seas,s] <- (obs_bar[r,yr_alt_idx,seas,s]-exp_bar[r,yr_alt_idx,seas,s])/sqrt(v_y[r,yr_alt_idx,seas,s]/tmp_iss_obs[r,1,1,s,1]) # get weights
            } # end s loop
          } # end r loop
        } # end if split by region and sex

        # If compositions are split by region, but joint by sex
        if(comp_type[y,f] == 2) {
          for(r in use_regions) {
            mat_exp <- matrix(tmp_exp[r,1,1,,,1], nrow = n_bins)
            mat_obs <- matrix(tmp_obs[r,1,1,,,1], nrow = n_bins)
            exp_bar[r,yr_alt_idx,seas,1] <- sum(bins * rowSums(mat_exp)) # input mean pred comps
            obs_bar[r,yr_alt_idx,seas,1] <- sum(bins * rowSums(mat_obs)) # input mean obs comps
            v_y[r,yr_alt_idx,seas,1] <- sum(bins^2*rowSums(mat_exp))-exp_bar[r,yr_alt_idx,seas,1]^2  # variance
            w_denom[r,yr_alt_idx,seas,1] <- (obs_bar[r,yr_alt_idx,seas,1]-exp_bar[r,yr_alt_idx,seas,1])/sqrt(v_y[r,yr_alt_idx,seas,1]/tmp_iss_obs[r,1,1,1,1]) # get weights
          } # end r loop
        } # end if split by region, joint by sex

      } # end seas loop
    } # end y loop

    # get unique composition types
    unique_comp_type <- unique(comp_type[,f])

    for(j in 1:length(unique_comp_type)) {

      # get year pointer index to subset w_denom and calculate weights separately for each composition type
      year_pointer <- which(comp_type[data_yrs,f] == unique_comp_type[j])

      for(seas in 1:n_seas) {
        # if aggregated or joint by region and sex
        if(unique_comp_type[j] == 0) weights[1,data_yrs,seas,1,f] <- safe_inv_var(w_denom[1,year_pointer,seas,1])

        # if split by sex and region
        if(unique_comp_type[j] == 1) for(r in 1:n_regions) for(s in 1:n_sexes) weights[r,data_yrs,seas,s,f] <- safe_inv_var(w_denom[r,year_pointer,seas,s])

        # if split by region, joint by sex
        if(unique_comp_type[j] == 2) for(r in 1:n_regions) weights[r,data_yrs,seas,1,f] <- safe_inv_var(w_denom[r,year_pointer,seas,1])
      } # end seas loop

    } # end j loop

    # Summarize observed and expected means
    tmp_means <- reshape2::melt(obs_bar) %>%
      dplyr::rename(Region = Var1, Comp_Year = Var2, Comp_Seas = Var3, Sex = Var4, obs = value) %>%
      dplyr::left_join(reshape2::melt(exp_bar) %>%
                         dplyr::rename(Region = Var1, Comp_Year = Var2, Comp_Seas = Var3, Sex = Var4, pred = value),
                       by = c("Region", "Comp_Year", "Comp_Seas", "Sex")) %>%
      dplyr::mutate(Fleet = f)

    mean_francis <- rbind(mean_francis, tmp_means)

  } # end f loop

  return(list(weights = weights, mean_francis = mean_francis))
} # end function

#' Get Francis Weights
#'
#' Computes Francis composition reweighting factors for pooled and
#' population-specific fishery and survey age and length compositions.
#' Used inside \code{\link{run_francis}} or directly in a user-defined
#' reweighting loop.
#'
#' @param data List of model data inputs containing observed compositions,
#'   usage flags, input sample sizes, weight arrays, and composition type
#'   matrices for both pooled and population-specific data streams.
#' @param rep Report list from a fitted RTMB model containing predicted
#'   compositions (\code{CAA}, \code{CAL}, \code{SrvIAA}, \code{SrvIAL})
#'   with a leading population dimension.
#' @param age_labels Vector of observed age bin labels.
#' @param len_labels Vector of observed length bin labels.
#' @param year_labels Vector of assessment year labels.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{new_fish_age_wts}, \code{new_fish_len_wts},
#'       \code{new_srv_age_wts}, \code{new_srv_len_wts}}{Updated Francis
#'       weight arrays for pooled compositions, same dimensions as the
#'       corresponding \code{Wt_*} arrays in \code{data}.}
#'     \item{\code{new_fish_age_pop_wts}, \code{new_fish_len_pop_wts},
#'       \code{new_srv_age_pop_wts}, \code{new_srv_len_pop_wts}}{Updated
#'       Francis weight arrays for population-specific compositions, same
#'       dimensions as the corresponding \code{Wt_*_pop} arrays in
#'       \code{data}. Cells remain \code{NA} when the corresponding
#'       \code{Use*_pop} flag contains no ones.}
#'     \item{\code{mean_francis}}{Long-format dataframe of observed and
#'       expected composition means across all data streams and populations,
#'       with columns \code{Region}, \code{Comp_Year}, \code{Comp_Seas},
#'       \code{Sex}, \code{Fleet}, \code{obs}, \code{pred}, \code{Type},
#'       and \code{Pop} (pooled rows have \code{Pop = NA}).}
#'   }
#'
#' @export do_francis_reweighting
#' @family Francis Reweighting
#'
#' @examples
#' \dontrun{
#' for(j in 1:5) {
#'   if(j == 1) {
#'     data$Wt_FishAgeComps[] <- 1; data$Wt_FishLenComps[] <- 1
#'     data$Wt_SrvAgeComps[]  <- 1; data$Wt_SrvLenComps[]  <- 1
#'     data$Wt_FishAgeComps_pop[] <- 0; data$Wt_FishLenComps_pop[] <- 0
#'     data$Wt_SrvAgeComps_pop[]  <- 0; data$Wt_SrvLenComps_pop[]  <- 0
#'   } else {
#'     data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
#'     data$Wt_FishLenComps[] <- wts$new_fish_len_wts
#'     data$Wt_SrvAgeComps[]  <- wts$new_srv_age_wts
#'     data$Wt_SrvLenComps[]  <- wts$new_srv_len_wts
#'     if(any(data$UseFishAgeComps_pop == 1)) data$Wt_FishAgeComps_pop[] <- wts$new_fish_age_pop_wts
#'     if(any(data$UseFishLenComps_pop == 1)) data$Wt_FishLenComps_pop[] <- wts$new_fish_len_pop_wts
#'     if(any(data$UseSrvAgeComps_pop  == 1)) data$Wt_SrvAgeComps_pop[]  <- wts$new_srv_age_pop_wts
#'     if(any(data$UseSrvLenComps_pop  == 1)) data$Wt_SrvLenComps_pop[]  <- wts$new_srv_len_pop_wts
#'   }
#'   obj <- fit_model(data, parameters, mapping, random = NULL,
#'                    newton_loops = 3, silent = TRUE)
#'   rep <- obj$report(obj$env$last.par.best)
#'   wts <- do_francis_reweighting(data = data, rep = rep, age_labels = 2:31,
#'                                 len_labels = seq(41, 99, 2), year_labels = 1960:2024)
#' }
#' }
do_francis_reweighting <- function(data,
                                   rep,
                                   age_labels,
                                   len_labels,
                                   year_labels
) {

  # Get indexing
  n_regions <- data$n_regions
  n_fish_fleets <- data$n_fish_fleets
  n_srv_fleets <- data$n_srv_fleets
  n_sexes <- data$n_sexes
  n_ages <- dim(data$ObsFishAgeComps)[4]
  n_lens <- dim(data$ObsFishLenComps)[4]
  n_years <- length(year_labels)
  n_seas <- data$n_seas

  # Get composition proportions
  comp_prop <- get_comp_prop(data = data,
                             rep = rep,
                             age_labels = age_labels,
                             len_labels = len_labels,
                             year_labels = year_labels)

  # Fishery Ages ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseFishAgeComps
  tmp_ISS <- data$ISS_FishAgeComps
  tmp_Pred <- comp_prop$Pred_FishAge_mat
  tmp_Obs <- comp_prop$Obs_FishAge_mat
  tmp_comp_type <- data$FishAgeComps_Type

  new_fish_age_wts <- data$Wt_FishAgeComps # matrix to store new weights
  new_fish_age_wts[] <- NA

  # Get francis weights here
  fish_age_info <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets,
    n_years = n_years,
    n_seas = n_seas,
    n_bins = n_ages,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_age_wts,
    bins = age_labels,
    comp_type = tmp_comp_type
  )

  new_fish_age_wts[] <- fish_age_info$weights

  # Fishery Ages Discards ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseFishAgeComps_discard
  tmp_ISS <- data$ISS_FishAgeComps_discard
  tmp_Pred <- comp_prop$Pred_FishAge_discard_mat
  tmp_Obs <- comp_prop$Obs_FishAge_discard_mat
  tmp_comp_type <- data$FishAgeComps_discard_Type

  new_fish_age_discard_wts <- data$Wt_FishAgeComps_discard # matrix to store new weights
  new_fish_age_discard_wts[] <- NA

  # Get francis weights here
  fish_age_discard_info <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets,
    n_years = n_years,
    n_seas = n_seas,
    n_bins = n_ages,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_age_discard_wts,
    bins = age_labels,
    comp_type = tmp_comp_type
  )

  new_fish_age_discard_wts[] <- fish_age_discard_info$weights

  # Fishery Lengths ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseFishLenComps
  tmp_ISS <- data$ISS_FishLenComps
  tmp_Pred <- comp_prop$Pred_FishLen_mat
  tmp_Obs <- comp_prop$Obs_FishLen_mat
  tmp_comp_type <- data$FishLenComps_Type

  new_fish_len_wts <- data$Wt_FishLenComps # matrix to store new weights
  new_fish_len_wts[] <- NA

  # Get francis weights here
  fish_len_info <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets,
    n_years = n_years,
    n_seas = n_seas,
    n_bins = n_lens,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_len_wts,
    bins = len_labels,
    comp_type = tmp_comp_type
  )

  new_fish_len_wts[] <- fish_len_info$weights

  # Fishery Lengths Discards ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseFishLenComps_discard
  tmp_ISS <- data$ISS_FishLenComps_discard
  tmp_Pred <- comp_prop$Pred_FishLen_discard_mat
  tmp_Obs <- comp_prop$Obs_FishLen_discard_mat
  tmp_comp_type <- data$FishLenComps_discard_Type

  new_fish_len_discard_wts <- data$Wt_FishLenComps_discard # matrix to store new weights
  new_fish_len_discard_wts[] <- NA

  # Get francis weights here
  fish_len_discard_info <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets,
    n_years = n_years,
    n_seas = n_seas,
    n_bins = n_lens,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_len_discard_wts,
    bins = len_labels,
    comp_type = tmp_comp_type
  )

  new_fish_len_discard_wts[] <- fish_len_discard_info$weights


  # Survey Ages ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseSrvAgeComps
  tmp_ISS <- data$ISS_SrvAgeComps
  tmp_Pred <- comp_prop$Pred_SrvAge_mat
  tmp_Obs <- comp_prop$Obs_SrvAge_mat
  tmp_comp_type <- data$SrvAgeComps_Type

  new_srv_age_wts <- data$Wt_SrvAgeComps # matrix to store new weights
  new_srv_age_wts[] <- NA

  # Get francis weights here
  srv_age_info <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_srv_fleets,
    n_years = n_years,
    n_seas = n_seas,
    n_bins = n_ages,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_srv_age_wts,
    bins = age_labels,
    comp_type = tmp_comp_type
  )

  new_srv_age_wts[] <- srv_age_info$weights

  # Survey Lengths ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseSrvLenComps
  tmp_ISS <- data$ISS_SrvLenComps
  tmp_Pred <- comp_prop$Pred_SrvLen_mat
  tmp_Obs <- comp_prop$Obs_SrvLen_mat
  tmp_comp_type <- data$SrvLenComps_Type

  new_srv_len_wts <- data$Wt_SrvLenComps # matrix to store new weights
  new_srv_len_wts[] <- NA

  # Get francis weights here
  srv_len_info <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_srv_fleets,
    n_years = n_years,
    n_seas = n_seas,
    n_bins = n_lens,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_srv_len_wts,
    bins = len_labels,
    comp_type = tmp_comp_type
  )

  new_srv_len_wts[] <- srv_len_info$weights


  # Population-specific comps ---------------------------------------------
  n_pop <- data$n_pop
  pop_dim_fish <- dim(data$Wt_FishAgeComps_pop)[-1]  # [n_regions, n_years, n_seas, n_sexes, n_fish_fleets]
  pop_dim_srv  <- dim(data$Wt_SrvAgeComps_pop)[-1]   # [n_regions, n_years, n_seas, n_sexes, n_srv_fleets]

  ### Fishery Ages (Pop) ------------------------------------------------------
  new_fish_age_pop_wts <- data$Wt_FishAgeComps_pop
  new_fish_age_pop_wts[] <- NA
  mean_francis_fish_age_pop <- data.frame()
  if(any(data$UseFishAgeComps_pop == 1)) {
    for(p in 1:n_pop) {
      tmp_info <- get_francis_weights(
        n_regions = n_regions, n_sexes = n_sexes, n_fleets = n_fish_fleets,
        n_years = n_years, n_seas = n_seas,
        Use      = array(data$UseFishAgeComps_pop[p,,,,],     dim = dim(data$UseFishAgeComps_pop)[-1]),
        ISS      = array(data$ISS_FishAgeComps_pop[p,,,,,],   dim = dim(data$ISS_FishAgeComps_pop)[-1]),
        Pred_array = array(comp_prop$Pred_FishAge_pop_mat[p,,,,,,], dim = dim(comp_prop$Pred_FishAge_pop_mat)[-1]),
        Obs_array  = array(comp_prop$Obs_FishAge_pop_mat[p,,,,,,],  dim = dim(comp_prop$Obs_FishAge_pop_mat)[-1]),
        weights    = array(new_fish_age_pop_wts[p,,,,,],      dim = pop_dim_fish),
        n_bins       = n_ages,
        bins = age_labels,
        comp_type  = data$FishAgeComps_pop_Type
      )
      new_fish_age_pop_wts[p,,,,,] <- tmp_info$weights
      mean_francis_fish_age_pop <- rbind(mean_francis_fish_age_pop, tmp_info$mean_francis %>% dplyr::mutate(Pop = p))
    }
  }


  ### Fishery Lengths (Pop) ---------------------------------------------------
  new_fish_len_pop_wts <- data$Wt_FishLenComps_pop
  new_fish_len_pop_wts[] <- NA
  mean_francis_fish_len_pop <- data.frame()
  if(any(data$UseFishLenComps_pop == 1)) {
    for(p in 1:n_pop) {
      tmp_info <- get_francis_weights(
        n_regions = n_regions, n_sexes = n_sexes, n_fleets = n_fish_fleets,
        n_years = n_years, n_seas = n_seas,
        Use      = array(data$UseFishLenComps_pop[p,,,,],     dim = dim(data$UseFishLenComps_pop)[-1]),
        ISS      = array(data$ISS_FishLenComps_pop[p,,,,,],   dim = dim(data$ISS_FishLenComps_pop)[-1]),
        Pred_array = array(comp_prop$Pred_FishLen_pop_mat[p,,,,,,], dim = dim(comp_prop$Pred_FishLen_pop_mat)[-1]),
        Obs_array  = array(comp_prop$Obs_FishLen_pop_mat[p,,,,,,],  dim = dim(comp_prop$Obs_FishLen_pop_mat)[-1]),
        weights    = array(new_fish_len_pop_wts[p,,,,,],      dim = pop_dim_fish),
        n_bins       = n_lens,
        bins = len_labels,
        comp_type  = data$FishLenComps_pop_Type
      )
      new_fish_len_pop_wts[p,,,,,] <- tmp_info$weights
      mean_francis_fish_len_pop <- rbind(mean_francis_fish_len_pop, tmp_info$mean_francis %>% dplyr::mutate(Pop = p))
    }
  }

  ### Fishery Age Discards (Pop) ------------------------------------------------------
  new_fish_age_discard_pop_wts <- data$Wt_FishAgeComps_discard_pop
  new_fish_age_discard_pop_wts[] <- NA
  mean_francis_fish_age_discard_pop <- data.frame()
  if(any(data$UseFishAgeComps_discard_pop == 1)) {
    for(p in 1:n_pop) {
      tmp_info <- get_francis_weights(
        n_regions = n_regions, n_sexes = n_sexes, n_fleets = n_fish_fleets,
        n_years = n_years, n_seas = n_seas,
        Use      = array(data$UseFishAgeComps_discard_pop[p,,,,],     dim = dim(data$UseFishAgeComps_discard_pop)[-1]),
        ISS      = array(data$ISS_FishAgeComps_discard_pop[p,,,,,],   dim = dim(data$ISS_FishAgeComps_discard_pop)[-1]),
        Pred_array = array(comp_prop$Pred_FishAge_discard_pop_mat[p,,,,,,], dim = dim(comp_prop$Pred_FishAge_discard_pop_mat)[-1]),
        Obs_array  = array(comp_prop$Obs_FishAge_discard_pop_mat[p,,,,,,],  dim = dim(comp_prop$Obs_FishAge_discard_pop_mat)[-1]),
        weights    = array(new_fish_age_discard_pop_wts[p,,,,,],      dim = pop_dim_fish),
        n_bins       = n_ages,
        bins = age_labels,
        comp_type  = data$FishAgeComps_discard_pop_Type
      )
      new_fish_age_discard_pop_wts[p,,,,,] <- tmp_info$weights
      mean_francis_fish_age_discard_pop <- rbind(mean_francis_fish_age_discard_pop, tmp_info$mean_francis %>% dplyr::mutate(Pop = p))
    }
  }

  ### Fishery Length Discards (Pop) ---------------------------------------------------
  new_fish_len_discard_pop_wts <- data$Wt_FishLenComps_discard_pop
  new_fish_len_discard_pop_wts[] <- NA
  mean_francis_fish_len_discard_pop <- data.frame()
  if(any(data$UseFishLenComps_discard_pop == 1)) {
    for(p in 1:n_pop) {
      tmp_info <- get_francis_weights(
        n_regions = n_regions, n_sexes = n_sexes, n_fleets = n_fish_fleets,
        n_years = n_years, n_seas = n_seas,
        Use      = array(data$UseFishLenComps_discard_pop[p,,,,],     dim = dim(data$UseFishLenComps_discard_pop)[-1]),
        ISS      = array(data$ISS_FishLenComps_discard_pop[p,,,,,],   dim = dim(data$ISS_FishLenComps_discard_pop)[-1]),
        Pred_array = array(comp_prop$Pred_FishLen_discard_pop_mat[p,,,,,,], dim = dim(comp_prop$Pred_FishLen_discard_pop_mat)[-1]),
        Obs_array  = array(comp_prop$Obs_FishLen_discard_pop_mat[p,,,,,,],  dim = dim(comp_prop$Obs_FishLen_discard_pop_mat)[-1]),
        weights    = array(new_fish_len_discard_pop_wts[p,,,,,],      dim = pop_dim_fish),
        n_bins       = n_lens,
        bins = len_labels,
        comp_type  = data$FishLenComps_discard_pop_Type
      )
      new_fish_len_discard_pop_wts[p,,,,,] <- tmp_info$weights
      mean_francis_fish_len_discard_pop <- rbind(mean_francis_fish_len_discard_pop, tmp_info$mean_francis %>% dplyr::mutate(Pop = p))
    }
  }

  ### Survey Ages (Pop) -------------------------------------------------------
  new_srv_age_pop_wts <- data$Wt_SrvAgeComps_pop
  new_srv_age_pop_wts[] <- NA
  mean_francis_srv_age_pop <- data.frame()
  if(any(data$UseSrvAgeComps_pop == 1)) {
    for(p in 1:n_pop) {
      tmp_info <- get_francis_weights(
        n_regions = n_regions, n_sexes = n_sexes, n_fleets = n_srv_fleets,
        n_years = n_years, n_seas = n_seas,
        Use      = array(data$UseSrvAgeComps_pop[p,,,,],     dim = dim(data$UseSrvAgeComps_pop)[-1]),
        ISS      = array(data$ISS_SrvAgeComps_pop[p,,,,,],   dim = dim(data$ISS_SrvAgeComps_pop)[-1]),
        Pred_array = array(comp_prop$Pred_SrvAge_pop_mat[p,,,,,,], dim = dim(comp_prop$Pred_SrvAge_pop_mat)[-1]),
        Obs_array  = array(comp_prop$Obs_SrvAge_pop_mat[p,,,,,,],  dim = dim(comp_prop$Obs_SrvAge_pop_mat)[-1]),
        weights    = array(new_srv_age_pop_wts[p,,,,,],      dim = pop_dim_srv),
        n_bins       = n_ages,
        bins = age_labels,
        comp_type  = data$SrvAgeComps_pop_Type
      )
      new_srv_age_pop_wts[p,,,,,] <- tmp_info$weights
      mean_francis_srv_age_pop <- rbind(mean_francis_srv_age_pop, tmp_info$mean_francis %>% dplyr::mutate(Pop = p))
    }
  }


  ### Survey Lengths (Pop) ----------------------------------------------------
  new_srv_len_pop_wts <- data$Wt_SrvLenComps_pop
  new_srv_len_pop_wts[] <- NA
  mean_francis_srv_len_pop <- data.frame()
  if(any(data$UseSrvLenComps_pop == 1)) {
    for(p in 1:n_pop) {
      tmp_info <- get_francis_weights(
        n_regions = n_regions, n_sexes = n_sexes, n_fleets = n_srv_fleets,
        n_years = n_years, n_seas = n_seas,
        Use      = array(data$UseSrvLenComps_pop[p,,,,],     dim = dim(data$UseSrvLenComps_pop)[-1]),
        ISS      = array(data$ISS_SrvLenComps_pop[p,,,,,],   dim = dim(data$ISS_SrvLenComps_pop)[-1]),
        Pred_array = array(comp_prop$Pred_SrvLen_pop_mat[p,,,,,,], dim = dim(comp_prop$Pred_SrvLen_pop_mat)[-1]),
        Obs_array  = array(comp_prop$Obs_SrvLen_pop_mat[p,,,,,,],  dim = dim(comp_prop$Obs_SrvLen_pop_mat)[-1]),
        weights    = array(new_srv_len_pop_wts[p,,,,,],      dim = pop_dim_srv),
        n_bins       = n_lens,
        bins = len_labels,
        comp_type  = data$SrvLenComps_pop_Type
      )
      new_srv_len_pop_wts[p,,,,,] <- tmp_info$weights
      mean_francis_srv_len_pop <- rbind(mean_francis_srv_len_pop, tmp_info$mean_francis %>% dplyr::mutate(Pop = p))
    }
  }

  # Get francis mean fits
  mean_francis <- rbind(
    fish_age_info$mean_francis %>% dplyr::mutate(Type = "Fishery Ages", Pop = NA),
    fish_age_discard_info$mean_francis %>% dplyr::mutate(Type = "Fishery Discard Ages", Pop = NA),
    srv_age_info$mean_francis  %>% dplyr::mutate(Type = "Survey Ages", Pop = NA),
    fish_len_info$mean_francis %>% dplyr::mutate(Type = "Fishery Lengths", Pop = NA),
    fish_len_discard_info$mean_francis %>% dplyr::mutate(Type = "Fishery Discard Lengths", Pop = NA),
    srv_len_info$mean_francis  %>% dplyr::mutate(Type = "Survey Lengths", Pop = NA),
    mean_francis_fish_age_pop  %>% dplyr::mutate(Type = "Pop Fishery Ages"),
    mean_francis_fish_len_pop  %>% dplyr::mutate(Type = "Pop Fishery Lengths"),
    mean_francis_fish_age_discard_pop  %>% dplyr::mutate(Type = "Pop Fishery Discard Ages"),
    mean_francis_fish_len_discard_pop  %>% dplyr::mutate(Type = "Pop Fishery Discard Lengths"),
    mean_francis_srv_age_pop   %>% dplyr::mutate(Type = "Pop Survey Ages"),
    mean_francis_srv_len_pop   %>% dplyr::mutate(Type = "Pop Survey Lengths")
  )

  return(list(new_fish_age_wts     = new_fish_age_wts,
              new_fish_len_wts     = new_fish_len_wts,
              new_fish_age_discard_wts     = new_fish_age_discard_wts,
              new_fish_len_discard_wts     = new_fish_len_discard_wts,
              new_srv_age_wts      = new_srv_age_wts,
              new_srv_len_wts      = new_srv_len_wts,
              new_fish_age_pop_wts = new_fish_age_pop_wts,
              new_fish_len_pop_wts = new_fish_len_pop_wts,
              new_fish_age_discard_pop_wts = new_fish_age_discard_pop_wts,
              new_fish_len_discard_pop_wts = new_fish_len_discard_pop_wts,
              new_srv_age_pop_wts  = new_srv_age_pop_wts,
              new_srv_len_pop_wts  = new_srv_len_pop_wts,
              mean_francis         = mean_francis
              ))


} # end function

#' Run Iterative Francis Reweighting Procedure
#'
#' Runs an iterative Francis reweighting procedure for pooled and
#' population-specific fishery and survey age and length compositions.
#' Repeatedly fits the model and updates composition weights until the
#' specified number of iterations is reached.
#'
#' @param data A list of model input data containing observed compositions,
#'   usage flags, input sample sizes, and weight arrays for both pooled
#'   (\code{Wt_FishAgeComps}, \code{Wt_FishLenComps}, \code{Wt_SrvAgeComps},
#'   \code{Wt_SrvLenComps}, \code{Wt_FishAgeComps_discard},
#'   \code{Wt_FishLenComps_discard}) and population-specific
#'   (\code{Wt_FishAgeComps_pop}, \code{Wt_FishLenComps_pop},
#'   \code{Wt_SrvAgeComps_pop}, \code{Wt_SrvLenComps_pop},
#'   \code{Wt_FishAgeComps_discard_pop}, \code{Wt_FishLenComps_discard_pop})
#'   data streams.
#' @param parameters A list of model parameters passed to \code{\link{fit_model}}.
#' @param mapping A list or mapping object passed to \code{\link{fit_model}}.
#' @param random Character vector of random effects passed to
#'   \code{\link{fit_model}}. Default \code{NULL}.
#' @param n_francis_iter Integer. Number of Francis reweighting iterations.
#'   Default \code{10}.
#' @param newton_loops Integer. Number of Newton refinement steps passed to
#'   \code{\link{fit_model}}. Default \code{0}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{obj}}{The fitted model object from the final iteration,
#'       augmented with \code{$data}, \code{$parameters}, \code{$mapping},
#'       \code{$random}, and \code{$rep}.}
#'     \item{\code{start_mean_francis}}{Mean composition fits from the first
#'       iteration, as returned by \code{\link{do_francis_reweighting}}.}
#'     \item{\code{end_mean_francis}}{Mean composition fits from the final
#'       iteration.}
#'     \item{\code{recorded_weights}}{Long-format dataframe of Francis weights
#'       from every iteration for all pooled and population-specific data
#'       streams, with columns \code{Region}, \code{Year}, \code{Seas},
#'       \code{Sex}, \code{Fleet}, \code{Weight}, \code{Type}, \code{Pop}
#'       (pooled rows have \code{Pop = NA}), and \code{iter}.}
#'   }
#'
#' @family Francis Reweighting
#' @export run_francis
#'
#' @examples
#' \dontrun{
#'   out <- run_francis(data = data, parameters = parameters,
#'                      mapping = mapping, random = NULL,
#'                      n_francis_iter = 5, newton_loops = 3)
#'   out$obj
#'   out$end_mean_francis
#'   out$recorded_weights
#' }
run_francis <- function(data,
                        parameters,
                        mapping,
                        random = NULL,
                        n_francis_iter = 10,
                        newton_loops = 0
                        ) {

  wts_df <- data.frame() # empty dataframe for recorded weights to bind to

  # run francis
  for(j in 1:n_francis_iter) {

    if(j == 1) {
      data$Wt_FishAgeComps[] <- 1
      data$Wt_FishLenComps[] <- 1
      data$Wt_FishAgeComps_discard[] <- 1
      data$Wt_FishLenComps_discard[] <- 1
      data$Wt_SrvAgeComps[]  <- 1
      data$Wt_SrvLenComps[]  <- 1
      data$Wt_FishAgeComps_pop[] <- 1
      data$Wt_FishLenComps_pop[] <- 1
      data$Wt_FishAgeComps_discard_pop[] <- 1
      data$Wt_FishLenComps_discard_pop[] <- 1
      data$Wt_SrvAgeComps_pop[]  <- 1
      data$Wt_SrvLenComps_pop[]  <- 1
    } else {
      data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
      data$Wt_FishLenComps[] <- wts$new_fish_len_wts
      data$Wt_SrvAgeComps[]  <- wts$new_srv_age_wts
      data$Wt_SrvLenComps[]  <- wts$new_srv_len_wts
      if(any(data$UseFishAgeComps_pop == 1)) data$Wt_FishAgeComps_pop[] <- wts$new_fish_age_pop_wts
      if(any(data$UseFishLenComps_pop == 1)) data$Wt_FishLenComps_pop[] <- wts$new_fish_len_pop_wts
      if(any(data$UseFishAgeComps_discard_pop == 1)) data$Wt_FishAgeComps_discard_pop[] <- wts$new_fish_age_discard_pop_wts
      if(any(data$UseFishLenComps_discard_pop == 1)) data$Wt_FishLenComps_discard_pop[] <- wts$new_fish_len_discard_pop_wts
      if(any(data$UseFishAgeComps_discard == 1)) data$Wt_FishAgeComps_discard[] <- wts$new_fish_age_discard_wts
      if(any(data$UseFishLenComps_discard == 1)) data$Wt_FishLenComps_discard[] <- wts$new_fish_len_discard_wts
      if(any(data$UseSrvAgeComps_pop  == 1)) data$Wt_SrvAgeComps_pop[]  <- wts$new_srv_age_pop_wts
      if(any(data$UseSrvLenComps_pop  == 1)) data$Wt_SrvLenComps_pop[]  <- wts$new_srv_len_pop_wts
    }

    # run model
    obj <- fit_model(data,
                     parameters,
                     mapping,
                     random = random,
                     newton_loops = newton_loops,
                     silent = TRUE
    )

    rep <- obj$report(obj$env$last.par.best) # Get report

    # get francis weights
    wts <- do_francis_reweighting(
      data = data, rep = rep,
      # uses fishery ages to index, because of potential for uneven number of observed and modelled ages
      age_labels = 1:dim(data$ObsFishAgeComps)[4],
      len_labels = data$lens,
      year_labels = data$years
    )

    # save first iteration run
    if(j == 1) wts_1 <- wts

    # record weights
    fish_age_wts_df <- reshape2::melt(wts$new_fish_age_wts); fish_age_wts_df$Type <- "Fishery Ages"
    fish_len_wts_df <- reshape2::melt(wts$new_fish_len_wts); fish_len_wts_df$Type <- "Fishery Lengths"
    fish_age_discard_wts_df <- reshape2::melt(wts$new_fish_age_discard_wts); fish_age_discard_wts_df$Type <- "Fishery Discard Ages"
    fish_len_discard_wts_df <- reshape2::melt(wts$new_fish_len_discard_wts); fish_len_discard_wts_df$Type <- "Fishery Discard Lengths"
    srv_age_wts_df  <- reshape2::melt(wts$new_srv_age_wts);  srv_age_wts_df$Type  <- "Survey Ages"
    srv_len_wts_df  <- reshape2::melt(wts$new_srv_len_wts);  srv_len_wts_df$Type  <- "Survey Lengths"
    colnames(fish_age_wts_df) <- c("Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(fish_len_wts_df) <- c("Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(fish_age_discard_wts_df) <- c("Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(fish_len_discard_wts_df) <- c("Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(srv_age_wts_df)  <- c("Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(srv_len_wts_df)  <- c("Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    fish_age_wts_df$Pop <- NA; fish_len_wts_df$Pop <- NA; fish_age_discard_wts_df$Pop <- NA; fish_len_discard_wts_df$Pop <- NA
    srv_age_wts_df$Pop  <- NA; srv_len_wts_df$Pop  <- NA

    fish_age_pop_wts_df <- reshape2::melt(wts$new_fish_age_pop_wts); fish_age_pop_wts_df$Type <- "Pop Fishery Ages"
    fish_len_pop_wts_df <- reshape2::melt(wts$new_fish_len_pop_wts); fish_len_pop_wts_df$Type <- "Pop Fishery Lengths"
    fish_age_discard_pop_wts_df <- reshape2::melt(wts$new_fish_age_discard_pop_wts); fish_age_discard_pop_wts_df$Type <- "Pop Fishery Discard Ages"
    fish_len_discard_pop_wts_df <- reshape2::melt(wts$new_fish_len_discard_pop_wts); fish_len_discard_pop_wts_df$Type <- "Pop Fishery Discard Lengths"
    srv_age_pop_wts_df  <- reshape2::melt(wts$new_srv_age_pop_wts);  srv_age_pop_wts_df$Type  <- "Pop Survey Ages"
    srv_len_pop_wts_df  <- reshape2::melt(wts$new_srv_len_pop_wts);  srv_len_pop_wts_df$Type  <- "Pop Survey Lengths"
    colnames(fish_age_pop_wts_df) <- c("Pop", "Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(fish_len_pop_wts_df) <- c("Pop", "Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(fish_age_discard_pop_wts_df) <- c("Pop", "Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(fish_len_discard_pop_wts_df) <- c("Pop", "Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(srv_age_pop_wts_df)  <- c("Pop", "Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")
    colnames(srv_len_pop_wts_df)  <- c("Pop", "Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type")

    tmp_wts_df <- rbind(fish_age_wts_df, fish_len_wts_df, fish_age_discard_wts_df, fish_len_discard_wts_df, srv_age_wts_df, srv_len_wts_df,
                        fish_age_pop_wts_df, fish_len_pop_wts_df, fish_age_discard_pop_wts_df, fish_len_discard_pop_wts_df, srv_age_pop_wts_df, srv_len_pop_wts_df)
    tmp_wts_df$iter <- j
    wts_df <- rbind(wts_df, tmp_wts_df)
    print(paste("Francis iteration:", j))

  } # end j loop

  obj$data <- data
  obj$parameters <- parameters
  obj$mapping <- mapping
  obj$random <- random
  obj$rep <- rep

  return(list(obj = obj, start_mean_francis = wts_1$mean_francis,
              end_mean_francis = wts$mean_francis, recorded_weights = wts_df))

}

