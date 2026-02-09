#' Computes Francis weights, which is used internally by do_francis_reweighting
#'
#' @param n_regions Number of regions
#' @param n_sexes Number of sexes
#' @param n_fleets Number of fleets (fishery or survey)
#' @param Use Array from data list that specifies whether to use data that year
#' @param ISS Input sample size array
#' @param Pred_array Predicted values array dimensioned by n_regions, n_years, n_seas, n_ages, n_sexes, n_fleets
#' @param Obs_array Observed values array dimensioned by n_regions, n_years, n_seas, n_ages, n_sexes, n_fleets
#' @param bins Vector of bins used (age or length)
#' @param weights Array of francis weights (NAs) to apply dimensioned by n_regions, n_years, n_seas, n_sexes, n_fleets
#' @param comp_type Matrix of composition structure types dimensioned by year and fleet
#'
#' @returns List of values for calculated francis weight, and a dataframe of observed and expected means
#' @keywords internal
get_francis_weights <- function(n_regions,
                                n_sexes,
                                n_fleets,
                                Use,
                                ISS,
                                Pred_array,
                                Obs_array,
                                weights,
                                bins,
                                comp_type
                                ) {

  mean_francis <- data.frame()
  n_years <- dim(Use)[2] # get n_years
  n_seas <- dim(Use)[3] # get n_seas

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
            exp_bar[r,yr_alt_idx,seas,1] <- sum(bins * rowSums(tmp_exp[r,1,1,,,1])) # input mean pred comps
            obs_bar[r,yr_alt_idx,seas,1] <- sum(bins * rowSums(tmp_obs[r,1,1,,,1])) # input mean obs comps
            v_y[r,yr_alt_idx,seas,1] <- sum(bins^2*rowSums(tmp_exp[r,1,1,,,1]))-exp_bar[r,yr_alt_idx,seas,1]^2  # variance
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
        if(unique_comp_type[j] == 0) weights[1,data_yrs,seas,1,f] <- 1 / var(w_denom[1,year_pointer,seas,1], na.rm = TRUE)

        # if split by sex and region
        if(unique_comp_type[j] == 1) for(r in 1:n_regions) for(s in 1:n_sexes) weights[r,data_yrs,seas,s,f] <- 1 / var(w_denom[r,year_pointer,seas,s], na.rm = TRUE)

        # if split by region, joint by sex
        if(unique_comp_type[j] == 2) for(r in 1:n_regions) weights[r,data_yrs,seas,1,f] <- 1 / var(w_denom[r,year_pointer,seas,1], na.rm = TRUE)
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
#' @param rep Report file list
#' @param age_labels Age labels
#' @param len_labels Length labels
#' @param year_labels Year labels
#' @param data List of data inputs
#'
#' @return A list object of francis weights (note that it will be NAs for some, if using jnt composition approaches - i.e., only uses one dimension), as well as a dataframe of francis mean fits
#' @export do_francis_reweighting
#' @family Francis Reweighting
#' @details
#' Function to get francis weights. Used inside the wrapper function run_francis(), or can be defined by the user as a loop to extract Francis weights (see example).
#'
#'
#' @examples
#' \dontrun{
#' for(j in 1:5) {
#'
#'   if(j == 1) { # reset weights at 1
#'     data$Wt_FishAgeComps[] <- 1
#'     data$Wt_FishLenComps[] <- 1
#'     data$Wt_SrvAgeComps[] <- 1
#'     data$Wt_SrvLenComps[] <- 1
#'   } else {
#'     data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
#'     data$Wt_FishLenComps[] <- wts$new_fish_len_wts
#'     data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
#'     data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
#'   }
#'
#'   sabie_rtmb_model <- fit_model(data,
#'                                 parameters,
#'                                 mapping,
#'                                 random = NULL,
#'                                 newton_loops = 3,
#'                                 silent = TRUE
#'   )
#'
#'   rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report
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
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_age_wts,
    bins = age_labels,
    comp_type = tmp_comp_type
  )

  new_fish_age_wts[] <- fish_age_info$weights

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
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_len_wts,
    bins = len_labels,
    comp_type = tmp_comp_type
  )

  new_fish_len_wts[] <- fish_len_info$weights

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
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_srv_len_wts,
    bins = len_labels,
    comp_type = tmp_comp_type
  )

  new_srv_len_wts[] <- srv_len_info$weights

  # Get francis mean fits
  mean_francis <- rbind(
    fish_age_info$mean_francis %>% dplyr::mutate(Type = "Fishery Ages"),
    srv_age_info$mean_francis %>% dplyr::mutate(Type = "Survey Ages"),
    fish_len_info$mean_francis %>% dplyr::mutate(Type = "Fishery Lengths"),
    srv_len_info$mean_francis %>% dplyr::mutate(Type = "Survey Lengths")
  )

  return(list(new_fish_age_wts = new_fish_age_wts, new_fish_len_wts = new_fish_len_wts,
              new_srv_age_wts = new_srv_age_wts, new_srv_len_wts = new_srv_len_wts, mean_francis = mean_francis))

} # end function

#' Run Iterative Francis Reweighting Procedure
#'
#' Runs an iterative Francis reweighting procedure for composition data
#' (fishery and survey age- and length-compositions). The function
#' reweights input data, repeatedly fits the model, and computes
#' updated Francis weights.
#'
#' @param data A list of model input data, including at least observed
#'   compositions (`ObsFishAgeComps`, `ObsFishLenComps`, `ObsSrvAgeComps`,
#'   `ObsSrvLenComps`) and corresponding weights
#'   (`Wt_FishAgeComps`, `Wt_FishLenComps`, `Wt_SrvAgeComps`,
#'   `Wt_SrvLenComps`).
#' @param parameters A list of model parameters to be passed to
#'   [fit_model()].
#' @param mapping A list or mapping object used to specify fixed or
#'   estimated parameters in [fit_model()].
#' @param n_francis_iter Integer. Number of Francis reweighting
#'   iterations to perform. Default is `10`.
#' @param newton_loops Integer. Number of Newton loops passed to
#'   [fit_model()]. Default is `0`.
#' @param random A character string of random effects passed to [fit_model()].
#'
#' @returns A list with three elements:
#' \describe{
#'   \item{obj}{The fitted model object returned by [fit_model()],
#'   including all elements of a TMB object, data, parameters, mapping, random effects specified, and report.}
#'   \item{end_mean_francis}{A summary of the mean Francis weights from the first iteration.}
#'   \item{end_mean_francis}{A summary of the mean Francis weights from the final iteration.}
#'   \item{recorded_weights}{A summary of recorded francis weights from each iteartion.}
#' }
#' @family Francis Reweighting
#' @export run_francis
#'
#' @examples
#' \dontrun{
#'   out <- run_francis(data = data,
#'                      parameters = parameters,
#'                      mapping = mapping,
#'                      random = NULL,
#'                      n_francis_iter = 5,
#'                      newton_loops = 3)
#'   out$obj
#'   out$mean_francis
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

    if(j == 1) { # reset weights at 1
      data$Wt_FishAgeComps[] <- 1
      data$Wt_FishLenComps[] <- 1
      data$Wt_SrvAgeComps[] <- 1
      data$Wt_SrvLenComps[] <- 1
    } else {
      # iterate francis weights
      data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
      data$Wt_FishLenComps[] <- wts$new_fish_len_wts
      data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
      data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
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
    # reshape matrices and label
    fish_age_wts_df <- reshape2::melt(wts$new_fish_age_wts); fish_age_wts_df$Type <- "Fishery Ages"
    fish_len_wts_df <- reshape2::melt(wts$new_fish_len_wts); fish_len_wts_df$Type <- "Fishery Lengths"
    srv_age_wts_df <- reshape2::melt(wts$new_srv_age_wts); srv_age_wts_df$Type <- "Survey Ages"
    srv_len_wts_df <- reshape2::melt(wts$new_srv_len_wts); srv_len_wts_df$Type <- "Survey Lengths"
    tmp_wts_df <- rbind(fish_age_wts_df, fish_len_wts_df, srv_age_wts_df, srv_len_wts_df) # bind dataframes together
    colnames(tmp_wts_df) <- c("Region", "Year", "Seas", "Sex", "Fleet", "Weight", "Type") # rename columns
    tmp_wts_df$iter <- j # indicate francis iteration
    wts_df <- rbind(wts_df, tmp_wts_df) # update dataframe
    print(paste("Francis iteration:", j)) # print francis iteration

  } # end j loop

  obj$data <- data
  obj$parameters <- parameters
  obj$mapping <- mapping
  obj$random <- random
  obj$rep <- rep

  return(list(obj = obj, start_mean_francis = wts_1$mean_francis,
              end_mean_francis = wts$mean_francis, recorded_weights = wts_df))

}

