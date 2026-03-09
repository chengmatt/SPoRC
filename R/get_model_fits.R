#' Extract Index Fit Results
#'
#' Generates a tidy dataframe of observed and predicted survey and fishery
#' indices from a fitted RTMB model, including standard errors, confidence
#' intervals, residuals, and catchability blocks.
#'
#' @param data List; input data used in the RTMB model. Must contain
#'   \code{ObsSrvIdx}, \code{ObsSrvIdx_SE}, \code{ObsFishIdx}, \code{ObsFishIdx_SE},
#'   \code{Wt_SrvIdx}, \code{Wt_FishIdx}, \code{srv_q_blocks}, and \code{fish_q_blocks}.
#' @param rep List; RTMB report output containing \code{PredSrvIdx} and \code{PredFishIdx}.
#' @param year_labs Vector; year labels to assign to the third dimension of the
#'   predicted indices and the columns of observed indices.
#'
#' @return A dataframe containing combined survey and fishery indices with the
#' following columns:
#' \itemize{
#'   \item \code{Region} – Region label (prefixed with "Region")
#'   \item \code{Year} – Year
#'   \item \code{Seas} – Season
#'   \item \code{Fleet} – Fleet identifier
#'   \item \code{Type} – Either "Survey" or "Fishery"
#'   \item \code{obs} – Observed index
#'   \item \code{value} – Predicted index
#'   \item \code{se} – Standard error of the observed index
#'   \item \code{lci}, \code{uci} – 95% confidence interval for the observed index
#'   \item \code{q_block} – Catchability block value
#'   \item \code{resid} – Log-scale residual (\eqn{\log(obs) - \log(predicted)})
#'   \item \code{Category} – Combined Type, Fleet, Season, and Q-block label
#' }
#'
#' @examples
#' \dontrun{
#' idx_fits <- get_idx_fits(
#'   data = data,
#'   rep = rep,
#'   year_labs = seq(1960, 2024, 1)
#' )
#' }
#'
#' @export get_idx_fits
#' @family Model Diagnostics
#' @import dplyr
#' @importFrom tidyr drop_na
get_idx_fits <- function(data,
                         rep,
                         year_labs
                         ) {

  colnames(data$ObsSrvIdx) <- year_labs
  colnames(data$ObsSrvIdx_SE) <- year_labs
  dimnames(rep$PredSrvIdx)[3] <- list(year_labs)
  colnames(data$ObsFishIdx) <- year_labs
  colnames(data$ObsFishIdx_SE) <- year_labs
  colnames(data$srv_q_blocks) <- year_labs
  colnames(data$fish_q_blocks) <- year_labs
  dimnames(rep$PredFishIdx)[3] <- list(year_labs)

  # Remove data not used(i.e., removing ghost fits)
  data$ObsFishIdx[which(data$UseFishIdx == 0)] <- NA
  data$ObsSrvIdx[which(data$UseSrvIdx == 0)] <- NA
  data$ObsFishIdx_SE[which(data$UseFishIdx == 0)] <- NA
  data$ObsSrvIdx_SE[which(data$UseSrvIdx == 0)] <- NA

  # Observed survey index
  obs_srv <- reshape2::melt(data$ObsSrvIdx) %>% dplyr::rename(obs = value) %>%
    dplyr::left_join(reshape2::melt(data$ObsSrvIdx_SE / sqrt(data$Wt_SrvIdx)) %>%  dplyr::rename(se = value), by = c("Var1", "Var2", "Var3", "Var4")) %>%
    dplyr::mutate(lci = exp(log(obs) - (1.96 * se)), uci = exp(log(obs) + (1.96 * se)), Type = 'Survey') %>%
    tidyr::drop_na() %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4)

  # Predicted survey index
  pred_srv <- reshape2::melt(rep$PredSrvIdx) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
    dplyr::mutate(Type = 'Survey') %>%
    dplyr::filter(Year %in% unique(obs_srv$Year), value != 0) %>%
    dplyr::group_by(Region, Year, Seas, Fleet, Type) %>%
    dplyr::summarise(value = sum(value))

  # Get survey catchability
  srv_q <- reshape2::melt(data$srv_q_blocks) %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3, q_block = value) %>%
    dplyr::mutate(Type = 'Survey')

  # combine survey results
  all_srv <- obs_srv %>%
    dplyr::left_join(pred_srv, by = c("Region", "Year", "Seas", "Fleet", 'Type')) %>%
    dplyr::left_join(srv_q, by = c("Region", "Year", "Fleet", 'Type')) %>%
    dplyr::mutate(resid = log(obs) - log(value))

  # Observed fishery index
  obs_fish <- reshape2::melt(data$ObsFishIdx) %>%
    dplyr::rename(obs = value) %>%
    dplyr::left_join(reshape2::melt(data$ObsFishIdx_SE / sqrt(data$Wt_FishIdx)) %>%
                       dplyr::rename(se = value), by = c("Var1", "Var2", "Var3", "Var4")) %>%
    dplyr::mutate(lci = exp(log(obs) - (1.96 * se)), uci = exp(log(obs) + (1.96 * se)), Type = 'Fishery') %>%
    tidyr::drop_na() %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4)

  # Predicted fishery index
  pred_fish <- reshape2::melt(rep$PredFishIdx) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
    dplyr::mutate(Type = 'Fishery') %>%
    dplyr::filter(Year %in% unique(obs_fish$Year), value != 0) %>%
    dplyr::group_by(Region, Year, Seas, Fleet, Type) %>%
    dplyr::summarise(value = sum(value))

  # Get fishery catchability
  fish_q <- reshape2::melt(data$fish_q_blocks) %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3, q_block = value) %>%
    dplyr::mutate(Type = 'Fishery')

  # combine survey results
  all_fish <- obs_fish %>%
    dplyr::left_join(pred_fish, by = c("Region", "Year", "Seas", "Fleet", 'Type')) %>%
    dplyr::left_join(fish_q, by = c("Region", "Year", "Fleet", 'Type')) %>%
    dplyr::mutate(resid = log(obs) - log(value))

  all_idx <- rbind(all_fish, all_srv) %>%
    dplyr::mutate(Category = paste(Type, Fleet, ", Seas", Seas, ", Q", q_block, sep = ''),
                  Region = paste("Region", Region))

  return(all_idx)
}

#' Restructure Composition Values
#'
#' Restructures observed and expected composition values (catch-at-age or survey
#' index-at-age) for use in Francis reweighting or other composition-based analyses.
#' The function can handle aggregated, split, or joint sex composition parameterizations
#' and can apply ageing error if provided.
#'
#' @param Exp Array; expected composition values indexed by
#'   \code{[region, pop, year, age/length bins, sex, fleet]}.
#' @param Obs Array; observed composition values indexed similarly to \code{Exp}.
#' @param Comp_Type Integer; composition parameterization type:
#'   \itemize{
#'     \item 0 = aggregated across sexes
#'     \item 1 = split by sex (no implicit sex ratio information)
#'     \item 2 = joint across sexes (implicit sex ratio information)
#'   }
#' @param age_or_len Integer; 0 for age compositions, 1 for length compositions.
#' @param AgeingError Matrix; ageing error transition matrix (applied to age compositions only).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{Exp} – array of expected composition values in observed bins
#'   \item \code{Obs} – array of observed composition values in observed bins
#' }
#'
#' @keywords internal
Restrc_Comps <- function(Exp,
                         Obs,
                         Comp_Type,
                         age_or_len,
                         AgeingError) {

  const <- 0

  # Add constant to observed and expected (gets normalized later on)
  Obs <- Obs + const
  Exp <- Exp + const

  # Dimensions
  n_regions <- dim(Exp)[1]
  n_model_bins <- dim(Exp)[4]
  n_obs_bins <- dim(Obs)[4]
  n_sexes <- dim(Exp)[5]

  # Storage (Expected values get converted to n_obs_bins dimensions)
  Exp_mat = array(NA, c(n_regions, n_obs_bins, n_sexes))
  Obs_mat = array(NA, c(n_regions, n_obs_bins, n_sexes))

  # Aggregated comps by sex and region
  if(Comp_Type == 0) {

    tmp_Exp = matrix(rowSums(matrix(Exp, nrow = n_model_bins)) / (n_sexes * n_regions), nrow = 1) # aggregate
    # Expected age bins get collapsed to observed age bins if ageing error is non-square
    if(age_or_len == 0) {
      tmp_Exp = tmp_Exp %*% AgeingError # apply ageing error
      tmp_Exp = tmp_Exp / sum(tmp_Exp) # renormalize
    }
    if(age_or_len == 1) tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize (lengths)

    # Normalize observed
    tmp_Obs =  (Obs[1,1,1,,1,1]) / sum( Obs[1,1,1,,1,1])

    # Input into storage matrix
    Exp_mat[1,,1] = tmp_Exp
    Obs_mat[1,,1] = tmp_Obs

  } # end if aggregated comps across sex

  # 'Split' comps by sex (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      for(r in 1:n_regions) {

        # Expected Values
        if(age_or_len == 0) {
          tmp_Exp = ((Exp[r,1,1,,s,1]) / sum(Exp[r,1,1,,s,1])) %*% AgeingError # Normalize temporary variable (ages), collapses to observed age bins if ageing error is non square
          tmp_Exp = tmp_Exp / sum(tmp_Exp) # renormalize
        }
        if(age_or_len == 1) tmp_Exp = (Exp[r,1,1,,s,1]) / sum(Exp[r,1,1,,s,1]) # normalize lengths

        tmp_Obs = (Obs[r,1,1,,s,1]) / sum(Obs[r,1,1,,s,1]) # Normalize observed temporary variable

        # Input into storage matrix
        Exp_mat[r,,s] = tmp_Exp
        Obs_mat[r,,s] = tmp_Obs
      } # end r loop
    } # end s loop
  } # end if 'Split' comps by sex

  if(Comp_Type == 2) {
    for(r in 1:n_regions) {
      # Expected values
      if(age_or_len == 0) { # if ages
        tmp_Exp = t(as.vector((Exp[r,1,1,,,1])/ sum(Exp[r,1,1,,,1]))) %*% kronecker(diag(n_sexes), AgeingError) # apply ageing error, collapses to observed age bins if ageing error is non square
        tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize to make sure sum to 1
      } # if ages
      if(age_or_len == 1) tmp_Exp = as.vector((Exp[r,1,1,,,1]) / sum((Exp[r,1,1,,,1]))) # Normalize temporary variable (lengths)

      tmp_Obs = (Obs[r,1,1,,,1]) / sum(Obs[r,1,1,,,1]) # Normalize observed temporary variable

      # Input into storage matrix
      Exp_mat[r,,] = array(tmp_Exp, dim = c(n_obs_bins, n_sexes))
      Obs_mat[r,,] = array(tmp_Obs, dim = c(n_obs_bins, n_sexes))
    } # end r loop
  } # end if 'Joint' comps by sex

  return(list(Exp = Exp_mat, Obs = Obs_mat))

} # end function


#' Get Composition Proportions from RTMB Output
#'
#' Extracts and normalizes age and length composition data for fishery and survey
#' fleets according to the assessment specifications from RTMB. This includes
#' both observed and expected compositions and returns them in array and
#' long-dataframe formats.
#'
#' @param data List; data inputs used by RTMB, containing observed compositions,
#'   composition types, aggregation types, ageing errors, and fleet/region/season info.
#' @param rep List; report output from RTMB containing predicted compositions (CAA, CAL, SrvIAA, SrvIAL).
#' @param year_labels Vector; year labels corresponding to the assessment years.
#' @param age_labels Vector; observed age labels used in the assessment.
#' @param len_labels Vector; observed length labels used in the assessment.
#'
#' @return List containing:
#' \itemize{
#'   \item \strong{Fishery_Ages, Fishery_Lens, Survey_Ages, Survey_Lens}: Dataframes with
#'         observed (`obs`) and expected (`pred`) compositions and associated metadata.
#'   \item \strong{Obs_FishAge_mat, Obs_FishLen_mat, Obs_SrvAge_mat, Obs_SrvLen_mat}:
#'         Arrays of observed compositions (dimensioned by region, year, season, bin, sex, fleet).
#'   \item \strong{Pred_FishAge_mat, Pred_FishLen_mat, Pred_SrvAge_mat, Pred_SrvLen_mat}:
#'         Arrays of expected compositions (dimensioned by region, year, season, bin, sex, fleet).
#' }
#'
#' @import dplyr
#' @importFrom tidyr drop_na
#' @export get_comp_prop
#' @family Model Diagnostics
#'
#' @examples
#' \dontrun{
#' comp_props <- get_comp_prop(
#'   data = data, rep = rep,
#'   age_labels = 2:31,
#'   len_labels = seq(41, 99, 2),
#'   year_labels = 1960:2024
#' )
#' }
get_comp_prop <- function(data,
                          rep,
                          age_labels,
                          len_labels,
                          year_labels
                          ) {

  # dimensinoing
  n_regions <- data$n_regions
  n_yrs <- length(data$years)
  n_seas <- data$n_seas
  n_fish_ages <- dim(data$ObsFishAgeComps)[4]
  n_srv_ages <- dim(data$ObsSrvAgeComps)[4]
  n_lens <- length(data$lens)
  n_sexes <- data$n_sexes
  n_fish_fleets <- data$n_fish_fleets
  n_srv_fleets <- data$n_srv_fleets

  # storage containers
  Obs_FishAge <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_fish_ages, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL, age_labels, NULL, NULL)) # Obs fishery ages
  Obs_FishLen <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_lens, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL, len_labels, NULL, NULL)) # Obs fishery lengths
  Obs_SrvAge <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_srv_ages, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, age_labels, NULL, NULL)) # Obs survey ages
  Obs_SrvLen <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_lens, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, len_labels, NULL, NULL)) # Obs survey lengths
  Pred_FishAge <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_fish_ages, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL,  age_labels, NULL, NULL)) # Predicted fishery ages
  Pred_FishLen <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL,  len_labels, NULL, NULL)) # Predicted fishery lengths
  Pred_SrvAge <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_srv_ages, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, age_labels, NULL, NULL)) # Predicted survey ages
  Pred_SrvLen <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, len_labels, NULL, NULL)) # Predicted survey lengths

  # Get quantities
  # setup ageing error if user-supplied is not year specific
  if(length(dim(data$AgeingError)) == 2) {
    AgeingError_t <- array(0, dim = c(length(data$years), dim(data$AgeingError)))
    for(i in 1:length(data$years)) AgeingError_t[i,,] <-  data$AgeingError
  }
  # ageing error if it is year specific (just reassigning)
  if(length(dim(data$AgeingError)) == 3) AgeingError_t <-  data$AgeingError

  AgeingError <- AgeingError_t # ageing errors
  CAA <- apply(rep$CAA, 2:7, sum) # catch at age
  CAL <- apply(rep$CAL, 2:7, sum) # catch at len
  SrvIAA <- apply(rep$SrvIAA, 2:7, sum) # survey at age
  SrvIAL <- apply(rep$SrvIAL, 2:7, sum) # survey at length

  # Observed quantities
  ObsFishAgeComps <- data$ObsFishAgeComps
  ObsFishLenComps <- data$ObsFishLenComps
  ObsSrvAgeComps <- data$ObsSrvAgeComps
  ObsSrvLenComps <- data$ObsSrvLenComps

  # Composition Types
  FishAge_CompType <- data$FishAgeComps_Type
  SrvAge_CompType <- data$SrvAgeComps_Type
  FishLen_CompType <- data$FishLenComps_Type
  SrvLen_CompType <- data$SrvLenComps_Type

  # Aggregation Types
  FishAge_comp_agg_type <- data$FishAge_comp_agg_type
  FishLen_comp_agg_type <- data$FishLen_comp_agg_type
  SrvAge_comp_agg_type <- data$SrvAge_comp_agg_type
  SrvLen_comp_agg_type <- data$SrvLen_comp_agg_type

  # Whether ouse comp data
  UseFishAgeComps <- data$UseFishAgeComps
  UseFishLenComps <- data$UseFishLenComps
  UseSrvAgeComps <- data$UseSrvAgeComps
  UseSrvLenComps <- data$UseSrvLenComps

  # Fishery Ages
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(seas in 1:n_seas) {

        # figure out regions with obs
        use_regions <- which(UseFishAgeComps[,y,seas,f] == 1)

        if(sum(use_regions) > 0) {
          Exp <- CAA[,y,seas,,,f, drop = FALSE] # expected
          Obs <- ObsFishAgeComps[,y,seas,,,f, drop = FALSE] # observed
          Comp_Type <- FishAge_CompType[y,f] # composition type

          # reformat expected compositions
          tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                    age_or_len = 0, AgeingError = AgeingError_t[y,,])
          # Input into storage
          Obs_FishAge[,y,seas,,,f] <- tmp_comps$Obs
          Pred_FishAge[,y,seas,,,f] <- tmp_comps$Exp
        }

      } # end seas loop
    } # end f
  } # end y

  # Fishery Lengths
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(seas in 1:n_seas) {

        use_regions <- which(UseFishLenComps[,y,seas,f] == 1) # figure out regions with obs

        if(sum(use_regions) > 0) {
          Exp <- CAL[,y,seas,,,f, drop = FALSE] # expected
          Obs <- ObsFishLenComps[,y,seas,,,f, drop = FALSE] # observed
          Comp_Type <- FishLen_CompType[y,f] # composition type

          # get compositions
          tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                    age_or_len = 1, AgeingError = NA)
          # Input into storage
          Obs_FishLen[,y,seas,,,f] <- tmp_comps$Obs
          Pred_FishLen[,y,seas,,,f] <- tmp_comps$Exp
        }

      } # end seas
    } # end f
  } # end y

  # Survey Ages
  for(y in 1:n_yrs) {
    for(f in 1:n_srv_fleets) {
      for(seas in 1:n_seas) {

        # figure out regions with obs
        use_regions <- which(UseSrvAgeComps[,y,seas,f] == 1)

        if(sum(use_regions) > 0) {
          Exp <- SrvIAA[,y,seas,,,f, drop = FALSE] # expected
          Obs <- ObsSrvAgeComps[,y,seas,,,f, drop = FALSE] # observed
          Comp_Type <- SrvAge_CompType[y,f] # composition type

          # reformat expected compositions
          tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                    age_or_len = 0, AgeingError = AgeingError_t[y,,])
          # Input into storage
          Obs_SrvAge[,y,seas,,,f] <- tmp_comps$Obs
          Pred_SrvAge[,y,seas,,,f] <- tmp_comps$Exp
        }

      } # end seas loop
    } # end f
  } # end y

  # Survey Lengths
  for(y in 1:n_yrs) {
    for(f in 1:n_srv_fleets) {
      for(seas in 1:n_seas) {

        # figure out regions with obs
        use_regions <- which(UseSrvLenComps[,y,seas,f] == 1)

        if(sum(use_regions) > 0) {
          Exp <- SrvIAL[,y,seas,,,f, drop = FALSE] # expected
          Obs <- ObsSrvLenComps[,y,seas,,,f, drop = FALSE] # observed
          Comp_Type <- SrvLen_CompType[y,f] # composition type

          # reformat expected compositions
          tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                    age_or_len = 1, AgeingError = NA)
          # Input into storage
          Obs_SrvLen[,y,seas,,,f] <- tmp_comps$Obs
          Pred_SrvLen[,y,seas,,,f] <- tmp_comps$Exp
        }

      } # end seas
    } # end f
  } # end y

  # Process outputs
  all_fishages <- reshape2::melt(Obs_FishAge) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishAge) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Age = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Fishery Ages')

  # Fish Lengths
  all_fishlens <- reshape2::melt(Obs_FishLen) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishLen) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Len = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Fishery Lengths')

  # Survey Ages
  all_srvages <- reshape2::melt(Obs_SrvAge) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvAge) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Age = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Survey Ages')

  # Survey Lengths
  all_srvlens <- reshape2::melt(Obs_SrvLen) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvLen) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Len = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Survey Lengths')

  return(list(# data frames of observed and expected comps
              Fishery_Ages = all_fishages,
              Fishery_Lens = all_fishlens,
              Survey_Ages = all_srvages,
              Survey_Lens = all_srvlens,

              # Arrays of observed comps
              Obs_FishAge_mat = Obs_FishAge,
              Obs_FishLen_mat = Obs_FishLen,
              Obs_SrvAge_mat = Obs_SrvAge,
              Obs_SrvLen_mat = Obs_SrvLen,

              # Arrays of expected comps
              Pred_FishAge_mat = Pred_FishAge,
              Pred_FishLen_mat = Pred_FishLen,
              Pred_SrvAge_mat = Pred_SrvAge,
              Pred_SrvLen_mat = Pred_SrvLen))

} # end function

#' Internal function to compute OSA residuals for a single composition slice
#'
#' Backend function used by [get_osa()] to calculate one-step-ahead (OSA) residuals
#' from observed and expected composition data.
#'
#' @param obs Matrix of observed compositions (rows = years, columns = bins)
#' @param exp Matrix of expected compositions (same shape as obs)
#' @param N Sample size for multinomial/Dirichlet-multinomial
#' @param DM_theta Dirichlet-multinomial overdispersion parameter(s)
#' @param LN_Sigma Logistic-normal covariance matrix
#' @param fleet Fleet identifier
#' @param index Vector of composition bins (last entry dropped for logistic-normal)
#' @param years Vector of years corresponding to rows of obs/exp
#' @param index_label Character describing ages or lengths
#' @param comp_like Integer specifying likelihood:
#'   \describe{
#'     \item{0}{multinomial}
#'     \item{1}{Dirichlet-multinomial}
#'     \item{2-4}{logistic-normal variants}
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{res}{Data frame of OSA residuals with columns fleet, index_label, year, index, resid}
#' }
#'
#' @keywords internal
run_osa <- function(obs,
                    exp,
                    N = NULL,
                    DM_theta = NULL,
                    LN_Sigma = NULL,
                    fleet,
                    index,
                    years,
                    index_label,
                    comp_like) {

  if (!requireNamespace("compResidual", quietly = TRUE)) {
    stop("Package 'compResidual' is required for get_osa(). Please install it with remotes::install_github('fishfollower/compResidual/compResidual').")
  } else{
    set.seed(722533474)

    # Multinomial
    if(comp_like == 0) {
      if(is.null(N)) stop("N is NULL. Please provide the appropriate values for the Multinomial!")
      o <- round(N * obs/rowSums(obs), 0) # get observed (needs to be integers)
      p <- exp/rowSums(exp) # get expected
      res <- compResidual::resMulti(t(o), t(p)) # get residuals
      # clean up residual dataframe
      mat <- t(matrix(res, nrow = nrow(res), ncol = ncol(res))) # coerce into matrix
      dimnames(mat) <- list(year = years, index = index[1:(length(index) - 1)]) # name matrix
    }

    # Dirichlet-Multinomial
    if(comp_like == 1) {
      if(is.null(N) || is.null(DM_theta)) stop("N or DM_theta is NULL. Please provide the appropriate values for the Dirichlet-multinomial!")
      o <- round(N * obs/rowSums(obs), 0) # get observed (needs to be integers)
      p <- N * DM_theta * exp/rowSums(exp) # get expected
      res <- compResidual::resDirM(obs = t(o), alpha = t(p)) # get residuals
      # clean up residual dataframe
      mat <- t(matrix(res, nrow = nrow(res), ncol = ncol(res))) # coerce into matrix
      dimnames(mat) <- list(year = years, index = index[1:(length(index) - 1)]) # name matrix
    }

    # Logistic Normal
    if(comp_like %in% c(2:4)) {
      if(is.null(LN_Sigma)) stop("LN_Sigma is NULL. Please provide the appropriate values for the Logistic Normal!")
      # create residual dataframe
      mat <-  matrix(NA, ncol = ncol(t(obs)), nrow = nrow(t(obs)) - 1)
      # Transpose
      obs <- t(obs)
      exp <- t(exp)

      # loop through to normalize compositions and get OSAs
      for(i in 1:length(years)) {

        # normalize compositions
        tmp_obs <- obs[,i] / sum(obs[,i])
        tmp_exp <- exp[,i] / sum(exp[,i])

        # figure out zeros and keep track of original indices
        zeros <- which(tmp_obs == 0)
        original_length <- length(tmp_obs)

        if(length(zeros) > 0) {
          # Keep track of non-zero indices for mapping back
          non_zero_indices <- setdiff(1:original_length, zeros)
          tmp_obs <- tmp_obs[-zeros] / sum(tmp_obs[-zeros]) # renormalize w/o zeros
          tmp_exp <- tmp_exp[-zeros] / sum(tmp_exp[-zeros]) # renormalize w/o zeros
          tmp_Sigma <- LN_Sigma[-zeros, -zeros]
        } else {
          tmp_Sigma <- LN_Sigma
          non_zero_indices <- 1:original_length
        }

        # transform and drop last bin after filtering zeros
        tmp_obs <- compResidual::logistictransf(tmp_obs, FALSE)
        tmp_exp <- compResidual::logistictransf(tmp_exp, FALSE)
        tmp_Sigma <- tmp_Sigma[-nrow(tmp_Sigma), -ncol(tmp_Sigma)] # remove last bins

        # Update non_zero_indices to account for dropped last bin
        non_zero_indices <- non_zero_indices[-length(non_zero_indices)]

        # set up TMB OSA lists
        dat <- list()
        dat$code <- 4
        dat$obs <- tmp_obs
        dat$mu <- tmp_exp
        dat$Sigma <- tmp_Sigma
        param <- list(dummy = 0)

        # get OSAs
        obj <- TMB::MakeADFun(dat, param, DLL = "compResidual", silent = F)
        opt <- nlminb(obj$par, obj$fn, obj$gr)
        tmp <- TMB::oneStepPredict(obj, observation.name = "obs",
                                   data.term.indicator = "keep",
                                   method = "oneStepGaussianOffMode",
                                   trace = FALSE, reverse = T)

        # store OSAs
        mat[non_zero_indices, i] <- tmp$residual
      } # end i loop

      # clean up amtrix
      mat <- t(mat) # transpose to year x bins
      dimnames(mat) <- list(year = years, index = index[1:(length(index) - 1)]) # name matrix

    } # end logistic normal

    res_df <- reshape2::melt(mat, value.name = "resid") %>% # turn into dataframe
      dplyr::mutate(fleet = fleet, index_label = index_label) %>%
      dplyr::relocate(fleet, index_label, .before = year)

    return(list(res = res_df))
  }

}

#' Compute OSA residuals for composition data
#'
#' Formats observed and expected composition data and calculates one-step-ahead
#' (OSA) residuals using multinomial, Dirichlet-multinomial, or logistic-normal
#' likelihoods. This function is the main interface for residual diagnostics,
#' internally calling [run_osa()] to perform the residual calculations.
#'
#' @param obs_mat Array of observed compositions, dimensioned by
#'   \code{[region, year, bin, sex, fleet]}. May contain \code{NA}s, which are
#'   removed when filtering by \code{years}.
#' @param exp_mat Array of expected compositions, dimensioned the same as
#'   \code{obs_mat}. May contain \code{NA}s, which are removed when filtering by
#'   \code{years}.
#' @param N Input (or effective if Multinomial) sample size. Dimensions depend on \code{comp_type}:
#'   \itemize{
#'     \item \code{comp_type = 0} (aggregated): vector of length \code{n_years}.
#'     \item \code{comp_type = 1} (split by region and sex): array
#'       \code{[n_regions, n_years, n_sexes]}.
#'     \item \code{comp_type = 2} (split by region, joint by sex): matrix
#'       \code{[n_regions, n_years]}.
#'   }
#'   For years without data, users can simply input an NA or any abritary number (it gets filtered out within the function).
#' @param years Vector of years to filter to if composition type is aggregated (0). Otherwise, this expects a list where each list element is a vector of years for each region where compositions are available for use (split by region and sex, or split by region, joint by sex).
#' @param fleet Fleet identifier (character or numeric) to filter to.
#' @param bins Vector of age or length bin labels corresponding to the
#'   composition categories.
#' @param comp_type Integer specifying how compositions are structured:
#'   \itemize{
#'     \item 0 = aggregated across regions and sexes
#'     \item 1 = split by region and sex
#'     \item 2 = split by region, joint by sex
#'   }
#' @param bin_label Character label describing whether bins represent ages or
#'   lengths.
#' @param comp_like Integer specifying the likelihood type (defaults to 0):
#'   \itemize{
#'     \item 0 = multinomial
#'     \item 1 = Dirichlet-multinomial
#'     \item 2–4 = logistic-normal variants
#'   }
#' @param DM_theta Dirichlet-multinomial overdispersion parameter(s). Dimensions
#'   must match \code{N}:
#'   \itemize{
#'     \item aggregated: scalar
#'     \item split by sex: matrix \code{[n_regions, n_sexes]}
#'     \item joint by sex: vector of length \code{n_regions}
#'   }
#' @param LN_Sigma Logistic-normal covariance matrix. Dimensions depend on
#'   \code{comp_type}:
#'   \itemize{
#'     \item aggregated: matrix \code{[n_bins, n_bins]}
#'     \item split by region and sex: array \code{[n_regions, n_bins, n_bins, n_sexes]}
#'     \item joint by sex: array \code{[n_regions, n_bins, n_bins]}
#'   }
#'   Use [get_logistN_Sigma()] to help construct this input.
#' @param addtocomp Constant that is added to compositions
#' @param seas Season index
#'
#' @return A list with one element:
#' \describe{
#'   \item{res}{Data frame of OSA residuals. Columns include:
#'     \code{fleet}, \code{index_label}, \code{year}, \code{index},
#'     \code{resid}, \code{region}, \code{seas}, \code{sex}, and \code{comp_type}.}
#' }
#'
#' @family Model Diagnostics
#' @import dplyr
#' @export get_osa
get_osa <- function(obs_mat,
                    exp_mat,
                    N = NULL,
                    DM_theta = NULL,
                    LN_Sigma = NULL,
                    years,
                    seas,
                    fleet,
                    bins,
                    comp_type,
                    bin_label,
                    comp_like = 0,
                    addtocomp = 0
                    ) {

  if (!requireNamespace("compResidual", quietly = TRUE)) {
    stop("Package 'compResidual' is required for get_osa(). Please follow installation instructions from https://github.com/fishfollower/compResidual/compResidual")
  } else{

    # get dimensions
    n_regions <- dim(obs_mat)[1]
    n_sexes <- dim(obs_mat)[5]

    # if comps are aggregated
    if(comp_type == 0) {

      obs <- obs_mat[,years,seas,,,fleet, drop = FALSE] # get filtered observed matrix
      exp <- exp_mat[,years,seas,,,fleet, drop = FALSE] # get filtered expected matrix
      tmp_obs <- obs[1,,1,,1,1] # only get a single sex and single region out since aggregated
      tmp_exp <- exp[1,,1,,1,1] # only get a single sex and single region out since aggregated

      # normalize
      tmp_obs <- (tmp_obs + addtocomp) / rowSums(tmp_obs + addtocomp)
      tmp_exp <- (tmp_exp + addtocomp) / rowSums(tmp_exp + addtocomp)

      # compute OSA
      tmp_osa <- run_osa(obs = tmp_obs, exp = tmp_exp, N = N, DM_theta = DM_theta,
                         years = years, comp_like = comp_like, LN_Sigma = LN_Sigma,
                         index = bins, fleet = as.character(fleet), index_label = bin_label)

      # Doing some naming stuff
      tmp_osa$res$region <- 1 # 1s below b/c aggregated across all dimensions
      tmp_osa$res$sex <- 1
      tmp_osa$res$seas <- seas
      tmp_osa$res$comp_type <- "Aggregated"
      osa_all <- tmp_osa
    }

    # if comps are split by region and sex
    if(comp_type == 1) {

      # empty dataframes to bind to
      res_all <- data.frame()
      agg_all <- data.frame()

      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {

          obs <- obs_mat[,years[[r]],seas,,,fleet, drop = FALSE] # get filtered observed matrix
          exp <- exp_mat[,years[[r]],seas,,,fleet, drop = FALSE] # get filtered expected matrix

          tmp_obs <- obs[r,,1,,s,1] # get observations
          tmp_exp <- exp[r,,1,,s,1] # get expected

          # normalize
          tmp_obs <- (tmp_obs + addtocomp) / rowSums(tmp_obs + addtocomp)
          tmp_exp <- (tmp_exp + addtocomp) / rowSums(tmp_exp + addtocomp)

          # compute OSA
          tmp_osa <- run_osa(obs = tmp_obs, exp = tmp_exp, N = N[r,years[[r]],s], DM_theta = DM_theta[r,s],
                             years = years[[r]], comp_like = comp_like, LN_Sigma = LN_Sigma[r,,,s],
                             index = bins, fleet = as.character(fleet), index_label = bin_label)

          # Doing some naming stuff
          tmp_osa$res$region <- r
          tmp_osa$res$sex <- s
          tmp_osa$res$seas <- seas
          tmp_osa$res$comp_type <- "SpltR_SpltS"

          res_all <- rbind(res_all, tmp_osa$res)

        } # end s loop
      } # end r loop

      osa_all <- list(res = res_all)

    } # end split region and sex

    # if comp types are join by sex, split by region
    if(comp_type == 2) {

      # empty dataframes to bind to
      res_all <- data.frame()
      agg_all <- data.frame()

      for(r in 1:n_regions) {

        obs <- obs_mat[,years[[r]],seas,,,fleet, drop = FALSE] # get filtered observed matrix
        exp <- exp_mat[,years[[r]],seas,,,fleet, drop = FALSE] # get filtered expected matrix

        # initialize to cbind
        tmp_obs <- NULL
        tmp_exp <- NULL

        for(s in 1:n_sexes) {
          tmp_obs <- cbind(tmp_obs, obs[r,,1,,s,1]) # get observations
          tmp_exp <- cbind(tmp_exp, exp[r,,1,,s,1]) # get expected
        } # end s loop

        # normalize
        tmp_obs <- (tmp_obs + addtocomp) / rowSums(tmp_obs + addtocomp)
        tmp_exp <- (tmp_exp + addtocomp) / rowSums(tmp_exp + addtocomp)

        # compute OSA
        tmp_osa <- run_osa(obs = tmp_obs, exp = tmp_exp, N = N[r,years[[r]]],
                           DM_theta = DM_theta[r], years = years[[r]], comp_like = comp_like, LN_Sigma = LN_Sigma[r,,],
                           index = paste(rep(1:n_sexes, each = length(bins)), "_", rep(bins, times = n_sexes), sep = ""),
                           fleet = as.character(fleet), index_label = bin_label)

        # Doing some naming stuff
        tmp_osa$res$region <- r
        tmp_osa$res$seas <- seas
        tmp_osa$res <- tmp_osa$res %>% dplyr::mutate(split_index = str_split(index, "_"),  # Split once and store as list
                                                     sex = sapply(split_index, `[`, 1),
                                                     index = sapply(split_index, `[`, 2)) %>% dplyr::select(-split_index)
        tmp_osa$res$comp_type <- "SpltR_JntS"

        res_all <- rbind(res_all, tmp_osa$res)
      } # end r loop

      osa_all <- list(res = res_all)

    } # end split region, joint by sex

    return(osa_all)
  }


}

#' Plot OSA residuals from outputs of get_osa
#'
#' Generates diagnostic plots for one-step-ahead (OSA) residuals. Includes
#' QQ-plots with SDNR annotations and bubble plots showing residual magnitude and sign.
#'
#' @param osa_results List obtained from get_osa() containing residuals dataframe.
#'
#' @return List of plots: \code{sdnr_plot} (QQ-plot) and \code{bubble_plot} (bubble residuals).
#' @export plot_resids
#' @family Model Diagnostics
#' @import dplyr
#' @import ggplot2
plot_resids <- function(osa_results) {

  # extract results
  res <- osa_results$res %>% dplyr::mutate(sign = ifelse(resid < 0, "Neg", "Pos"))

  # Aggregated Comps
  if(unique(res$comp_type) == "Aggregated") {

    # Get standarized normal residuals
    sdnr <- res %>% dplyr::summarise(sdnr = paste0("SDNR = ", formatC(round(sd(resid, na.rm = TRUE), 3), format = "f", digits = 2)))

    # sdnr plot
    sdnr_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      theme_bw(base_size = 20) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 4)
  }

  # Split Sex and Split Region
  if(unique(res$comp_type) == "SpltR_SpltS") {

    # Get standarized normal residuals
    sdnr <- res %>% dplyr::group_by(region, sex) %>%
      dplyr::summarise(sdnr = paste0("SDNR = ", formatC(round(sd(resid, na.rm = TRUE), 3), format = "f", digits = 2)))

    # sdnr plot
    sdnr_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      facet_grid(region ~ sex, labeller = labeller(
        region = function(x) paste0("Region ", x),
        sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 4)
  }


  # Joint Sex and Split Region
  if(unique(res$comp_type) == "SpltR_JntS") {

    # Get standarized normal residuals
    sdnr <- res %>% dplyr::group_by(region) %>%
      dplyr::summarise(sdnr = paste0("SDNR = ", formatC(round(sd(resid, na.rm = TRUE), 3), format = "f", digits = 2)))

    # sdnr plot
    sdnr_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      facet_grid(~region, labeller = labeller(
        region = function(x) paste0("Region ", x)
      )) +
      theme_bw(base_size = 20) +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 4)
  }

  # bubble plot
  bubble_plot <- ggplot(data = res, aes(x = year, y = as.numeric(index),
                                        color = sign, size = abs(resid), alpha = abs(resid))) +
    geom_point() +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "Year", y = unique(res$index_label), color = "Sign", size = "abs(Resid)", alpha = "abs(Resid)") +
    facet_grid(region ~ sex, labeller = labeller(
      region = function(x) paste0("Region ", x),
      sex = function(x) paste0("Sex ", x)
    )) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'top')

  return(list(sdnr_plot, bubble_plot))

}


