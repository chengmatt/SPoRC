#' Extract Index Fit Results
#'
#' Generates a tidy dataframe of observed and predicted survey and fishery
#' indices from a fitted RTMB model, including standard errors, confidence
#' intervals, a raw log-scale (Pearson-style) residual, and catchability
#' blocks. Both pooled and population-specific indices are returned when the
#' corresponding \code{Use*_pop} flags contain any ones.
#'
#' The \code{resid} column here is the simple log-scale residual
#' \eqn{\log(\text{obs}) - \log(\text{predicted})} -- \emph{not} a
#' one-step-ahead (OSA) residual. For properly decorrelated OSA index
#' residuals (with QQ-plots and SDNR diagnostics via \code{\link{plot_resids}}),
#' use \code{\link{get_osa}} with \code{index_source = }. The observed-vs-
#' predicted \emph{fit} plot built from this function's output is
#' \code{\link{get_idx_fits_plot}}.
#'
#' @param data List. Input data used in the RTMB model. Must contain
#'   \code{ObsSrvIdx}, \code{ObsSrvIdx_SE}, \code{ObsFishIdx},
#'   \code{ObsFishIdx_SE}, \code{Wt_SrvIdx}, \code{Wt_FishIdx},
#'   \code{srv_q_blocks}, \code{fish_q_blocks}, \code{UseFishIdx}, and
#'   \code{UseSrvIdx}. For population-specific indices, also requires
#'   \code{ObsSrvIdx_pop}, \code{ObsSrvIdx_pop_SE}, \code{ObsFishIdx_pop},
#'   \code{ObsFishIdx_pop_SE}, \code{Wt_SrvIdx_pop}, \code{Wt_FishIdx_pop},
#'   \code{UseSrvIdx_pop}, and \code{UseFishIdx_pop}.
#' @param rep List. RTMB report output containing \code{PredSrvIdx} and
#'   \code{PredFishIdx}, both dimensioned
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fleets]}. Pooled indices
#'   are obtained by summing across the population dimension; population-specific
#'   indices use each population slice directly.
#' @param year_labs Vector. Year labels assigned to the year dimension of
#'   predicted and observed index arrays.
#'
#' @return A dataframe containing pooled and, when active, population-specific
#'   survey and fishery indices with the following columns:
#'   \describe{
#'     \item{\code{Region}}{Region label (prefixed with \code{"Region"}).}
#'     \item{\code{Year}}{Year.}
#'     \item{\code{Seas}}{Season.}
#'     \item{\code{Fleet}}{Fleet identifier.}
#'     \item{\code{Type}}{One of \code{"Survey"}, \code{"Fishery"},
#'       \code{"Pop Survey"}, or \code{"Pop Fishery"}.}
#'     \item{\code{obs}}{Observed index value.}
#'     \item{\code{value}}{Predicted index value.}
#'     \item{\code{se}}{Standard error of the observed index (weight-adjusted).}
#'     \item{\code{lci}, \code{uci}}{95\% log-normal confidence interval for
#'       the observed index.}
#'     \item{\code{q_block}}{Catchability block identifier.}
#'     \item{\code{resid}}{Log-scale residual
#'       (\eqn{\log(\text{obs}) - \log(\text{predicted})}).}
#'     \item{\code{Category}}{Combined label of Type, population (for pop
#'       rows), Fleet, Season, and Q-block.}
#'   }
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

  # Pop-specific dimnames — year is dim 3
  if(any(data$UseFishIdx_pop == 1)) {
    dimnames(data$ObsFishIdx_pop)[3]    <- list(year_labs)
    dimnames(data$ObsFishIdx_pop_SE)[3] <- list(year_labs)
  }
  if(any(data$UseSrvIdx_pop == 1)) {
    dimnames(data$ObsSrvIdx_pop)[3]    <- list(year_labs)
    dimnames(data$ObsSrvIdx_pop_SE)[3] <- list(year_labs)
  }

  # Remove data not used
  data$ObsFishIdx[which(data$UseFishIdx == 0)] <- NA
  data$ObsSrvIdx[which(data$UseSrvIdx == 0)] <- NA
  data$ObsFishIdx_SE[which(data$UseFishIdx == 0)] <- NA
  data$ObsSrvIdx_SE[which(data$UseSrvIdx == 0)] <- NA

  if(any(data$UseFishIdx_pop == 1)) {
    data$ObsFishIdx_pop[which(data$UseFishIdx_pop == 0)]    <- NA
    data$ObsFishIdx_pop_SE[which(data$UseFishIdx_pop == 0)] <- NA
  }
  if(any(data$UseSrvIdx_pop == 1)) {
    data$ObsSrvIdx_pop[which(data$UseSrvIdx_pop == 0)]    <- NA
    data$ObsSrvIdx_pop_SE[which(data$UseSrvIdx_pop == 0)] <- NA
  }

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

  # Population-specific Survey Index
  if(any(data$UseSrvIdx_pop == 1)) {

    data$ObsSrvIdx_pop[which(data$UseSrvIdx_pop == 0)] <- NA
    data$ObsSrvIdx_pop_SE[which(data$UseSrvIdx_pop == 0)] <- NA

    obs_srv_pop <- reshape2::melt(data$ObsSrvIdx_pop) %>%
      dplyr::rename(obs = value) %>%
      dplyr::left_join(reshape2::melt(data$ObsSrvIdx_pop_SE / sqrt(data$Wt_SrvIdx_pop)) %>%
                         dplyr::rename(se = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
      dplyr::mutate(lci = exp(log(obs) - (1.96 * se)), uci = exp(log(obs) + (1.96 * se)), Type = 'Pop Survey') %>%
      tidyr::drop_na() %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5)

    pred_srv_pop <- reshape2::melt(rep$PredSrvIdx) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
      dplyr::mutate(Type = 'Pop Survey') %>%
      dplyr::filter(Year %in% unique(obs_srv_pop$Year), value != 0)

    all_srv_pop <- obs_srv_pop %>%
      dplyr::left_join(pred_srv_pop, by = c("Pop", "Region", "Year", "Seas", "Fleet", "Type")) %>%
      dplyr::left_join(srv_q %>% dplyr::mutate(Type = 'Pop Survey'),
                       by = c("Region", "Year", "Fleet", "Type")) %>% # modifying to allow joining
      dplyr::mutate(resid = log(obs) - log(value),
                    Category = paste(Type, " Pop", Pop, ", Fleet", Fleet, ", Seas", Seas, ", Q", q_block, sep = ''),
                    Region = paste("Region", Region))

    all_idx <- rbind(all_idx, all_srv_pop %>% dplyr::select(-Pop))
  }

  # Population-specific Fishery Index
  if(any(data$UseFishIdx_pop == 1)) {

    data$ObsFishIdx_pop[which(data$UseFishIdx_pop == 0)] <- NA
    data$ObsFishIdx_pop_SE[which(data$UseFishIdx_pop == 0)] <- NA

    obs_fish_pop <- reshape2::melt(data$ObsFishIdx_pop) %>%
      dplyr::rename(obs = value) %>%
      dplyr::left_join(reshape2::melt(data$ObsFishIdx_pop_SE / sqrt(data$Wt_FishIdx_pop)) %>%
                         dplyr::rename(se = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
      dplyr::mutate(lci = exp(log(obs) - (1.96 * se)), uci = exp(log(obs) + (1.96 * se)), Type = 'Pop Fishery') %>%
      tidyr::drop_na() %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5)

    pred_fish_pop <- reshape2::melt(rep$PredFishIdx) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
      dplyr::mutate(Type = 'Pop Fishery') %>%
      dplyr::filter(Year %in% unique(obs_fish_pop$Year), value != 0)

    all_fish_pop <- obs_fish_pop %>%
      dplyr::left_join(pred_fish_pop, by = c("Pop", "Region", "Year", "Seas", "Fleet", "Type")) %>%
      dplyr::left_join(fish_q %>% dplyr::mutate(Type = 'Pop Fishery'),
                       by = c("Region", "Year", "Fleet", "Type")) %>% # modifying to allow joining
      dplyr::mutate(resid = log(obs) - log(value),
                    Category = paste(Type, " Pop", Pop, ", Fleet", Fleet, ", Seas", Seas, ", Q", q_block, sep = ''),
                    Region = paste("Region", Region))

    all_idx <- rbind(all_idx, all_fish_pop %>% dplyr::select(-Pop))
  }

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
#' Extracts and standardizes age and length composition data for fishery and
#' survey fleets from RTMB model output. The function processes both observed
#' and expected compositions for pooled (all populations combined) and
#' population-specific data streams, returning results in both long-format
#' data frames and array formats suitable for diagnostics and likelihood
#' evaluation.
#'
#' Compositions include retained catch, discard catch, and survey observations,
#' and are aligned across regions, years, seasons, fleets, sexes, and age/length
#' bins. Observed compositions are optionally filtered by usage flags defined in
#' the input data object.
#'
#' @param data List. RTMB input data containing observed compositions, usage
#'   flags, composition types, ageing error matrices, and fleet/region/season/
#'   population dimensions.
#'
#' @param rep List. RTMB model output containing predicted compositions
#'   (\code{CAA}, \code{CAL}, \code{DAA}, \code{DAL}, \code{SrvIAA}, \code{SrvIAL})
#'   with a leading population dimension that is summed to produce pooled outputs.
#'
#' @param year_labels Vector. Labels for model years corresponding to projection
#'   or estimation periods.
#'
#' @param age_labels Vector. Age bin labels used in composition data.
#'
#' @param len_labels Vector. Length bin labels used in composition data.
#'
#' @return A named list containing:
#'
#' \describe{
#'   \item{Fishery_Ages, Fishery_Lens}{Long-format data frames of observed and
#'   predicted pooled fishery age and length compositions, including metadata
#'   (Region, Year, Season, Age/Len, Sex, Fleet, Type).}
#'
#'   \item{Fishery_Discard_Ages, Fishery_Discard_Lens}{Long-format data frames
#'   for discard compositions.}
#'
#'   \item{Survey_Ages, Survey_Lens}{Long-format data frames for survey-based
#'   compositions.}
#'
#'   \item{Pop_Fishery_Ages, Pop_Fishery_Lens, Pop_Survey_Ages, Pop_Survey_Lens}{Population-specific
#'   long-format composition outputs with an additional Pop column.}
#'
#'   \item{Obs_*_mat}{Arrays of observed compositions indexed by
#'   Region × Year × Season × Bin × Sex × Fleet (and Pop where applicable).}
#'
#'   \item{Pred_*_mat}{Arrays of model-predicted compositions with matching
#'   dimensions to observed arrays.}
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
#'   data = data,
#'   rep = rep,
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
  n_pop <- data$n_pop
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
  Obs_FishAge_discard <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_fish_ages, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL, age_labels, NULL, NULL)) # Obs discard fishery ages
  Obs_FishLen_discard <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_lens, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL, len_labels, NULL, NULL)) # Obs discard fishery lengths
  Obs_SrvAge <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_srv_ages, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, age_labels, NULL, NULL)) # Obs survey ages
  Obs_SrvLen <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_lens, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, len_labels, NULL, NULL)) # Obs survey lengths
  Pred_FishAge <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_fish_ages, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL,  age_labels, NULL, NULL)) # Predicted fishery ages
  Pred_FishLen <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL,  len_labels, NULL, NULL)) # Predicted fishery lengths
  Pred_FishAge_discard <- array(data = NA, dim = c(n_regions, n_yrs,  n_seas, n_fish_ages, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL,  age_labels, NULL, NULL)) # Predicted discard fishery ages
  Pred_FishLen_discard <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, NULL,  len_labels, NULL, NULL)) # Predicted discard fishery lengths
  Pred_SrvAge <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_srv_ages, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, age_labels, NULL, NULL)) # Predicted survey ages
  Pred_SrvLen <- array(data = NA, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, NULL, len_labels, NULL, NULL)) # Predicted survey lengths
  Obs_FishAge_pop  <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_ages, n_sexes, n_fish_fleets))
  Pred_FishAge_pop <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_ages, n_sexes, n_fish_fleets))
  Obs_FishAge_discard_pop  <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_ages, n_sexes, n_fish_fleets))
  Pred_FishAge_discard_pop <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_ages, n_sexes, n_fish_fleets))
  Obs_FishLen_pop  <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens,      n_sexes, n_fish_fleets))
  Pred_FishLen_pop <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens,      n_sexes, n_fish_fleets))
  Obs_FishLen_discard_pop  <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens,      n_sexes, n_fish_fleets))
  Pred_FishLen_discard_pop <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens,      n_sexes, n_fish_fleets))
  Obs_SrvAge_pop   <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_srv_ages,  n_sexes, n_srv_fleets))
  Pred_SrvAge_pop  <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_srv_ages,  n_sexes, n_srv_fleets))
  Obs_SrvLen_pop   <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens,      n_sexes, n_srv_fleets))
  Pred_SrvLen_pop  <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens,      n_sexes, n_srv_fleets))

  # Get quantities
  # setup ageing error if user-supplied is not year specific
  if(length(dim(data$AgeingError)) == 2) {
    AgeingError_t <- array(0, dim = c(length(data$years), dim(data$AgeingError)))
    for(i in 1:length(data$years)) AgeingError_t[i,,] <-  data$AgeingError
  }
  # ageing error if it is year specific (just reassigning)
  if(length(dim(data$AgeingError)) == 3) AgeingError_t <-  data$AgeingError

  AgeingError <- AgeingError_t # ageing errors

  # sum across pops
  CAA <- apply(rep$CAA, 2:7, sum) # retained catch at age
  CAL <- apply(rep$CAL, 2:7, sum) # retained catch at len
  DAA <- apply(rep$DAA, 2:7, sum) # discarded catch at age
  DAL <- apply(rep$DAL, 2:7, sum) # discarded catch at len
  SrvIAA <- apply(rep$SrvIAA, 2:7, sum) # survey at age
  SrvIAL <- apply(rep$SrvIAL, 2:7, sum) # survey at length

  # pop-specific quantities
  CAA_pop <- rep$CAA
  CAL_pop <- rep$CAL
  DAA_pop <- rep$DAA
  DAL_pop <- rep$DAL
  SrvIAA_pop <- rep$SrvIAA
  SrvIAL_pop <- rep$SrvIAL

  # Observed quantities
  ObsFishAgeComps <- data$ObsFishAgeComps
  ObsFishLenComps <- data$ObsFishLenComps
  ObsFishAgeComps_discard <- data$ObsFishAgeComps_discard
  ObsFishLenComps_discard <- data$ObsFishLenComps_discard
  ObsSrvAgeComps <- data$ObsSrvAgeComps
  ObsSrvLenComps <- data$ObsSrvLenComps

  ObsFishAgeComps_pop <- data$ObsFishAgeComps_pop
  ObsFishLenComps_pop <- data$ObsFishLenComps_pop
  ObsFishAgeComps_discard_pop <- data$ObsFishAgeComps_discard_pop
  ObsFishLenComps_discard_pop <- data$ObsFishLenComps_discard_pop
  ObsSrvAgeComps_pop <- data$ObsSrvAgeComps_pop
  ObsSrvLenComps_pop <- data$ObsSrvLenComps_pop

  # Composition Types
  FishAge_CompType <- data$FishAgeComps_Type
  FishAgeComps_discard_Type <- data$FishAgeComps_discard_Type
  SrvAge_CompType <- data$SrvAgeComps_Type
  FishLen_CompType <- data$FishLenComps_Type
  FishLenComps_discard_Type <- data$FishLenComps_discard_Type
  SrvLen_CompType <- data$SrvLenComps_Type

  FishAge_pop_CompType <- data$FishAgeComps_pop_Type
  FishAgeComps_discard_pop_Type <- data$FishAgeComps_discard_pop_Type
  SrvAge_pop_CompType <- data$SrvAgeComps_pop_Type
  FishLen_pop_CompType <- data$FishLenComps_pop_Type
  FishLenComps_discard_pop_Type <- data$FishLenComps_discard_pop_Type
  SrvLen_pop_CompType <- data$SrvLenComps_pop_Type

  # Whether ouse comp data
  UseFishAgeComps <- data$UseFishAgeComps
  UseFishLenComps <- data$UseFishLenComps
  UseFishAgeComps_discard <- data$UseFishAgeComps_discard
  UseFishLenComps_discard <- data$UseFishLenComps_discard
  UseSrvAgeComps <- data$UseSrvAgeComps
  UseSrvLenComps <- data$UseSrvLenComps

  UseFishAgeComps_pop <- data$UseFishAgeComps_pop
  UseFishLenComps_pop <- data$UseFishLenComps_pop
  UseFishAgeComps_discard_pop <- data$UseFishAgeComps_discard_pop
  UseFishLenComps_discard_pop <- data$UseFishLenComps_discard_pop
  UseSrvAgeComps_pop <- data$UseSrvAgeComps_pop
  UseSrvLenComps_pop <- data$UseSrvLenComps_pop

  # Fishery Ages
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(seas in 1:n_seas) {

        # figure out regions with obs
        use_regions <- which(UseFishAgeComps[,y,seas,f] == 1)

        if(length(use_regions) > 0) {
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

        if(length(use_regions) > 0) {
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


  # Fishery Ages Discards
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(seas in 1:n_seas) {

        # figure out regions with obs
        use_regions <- which(UseFishAgeComps_discard[,y,seas,f] == 1)

        if(length(use_regions) > 0) {
          Exp <- DAA[,y,seas,,,f, drop = FALSE] # expected
          Obs <- ObsFishAgeComps_discard[,y,seas,,,f, drop = FALSE] # observed
          Comp_Type <- FishAgeComps_discard_Type[y,f] # composition type

          # reformat expected compositions
          tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                    age_or_len = 0, AgeingError = AgeingError_t[y,,])
          # Input into storage
          Obs_FishAge_discard[,y,seas,,,f] <- tmp_comps$Obs
          Pred_FishAge_discard[,y,seas,,,f] <- tmp_comps$Exp
        }

      } # end seas loop
    } # end f
  } # end y

  # Fishery Lengths Discards
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(seas in 1:n_seas) {

        use_regions <- which(UseFishLenComps_discard[,y,seas,f] == 1) # figure out regions with obs

        if(length(use_regions) > 0) {
          Exp <- DAL[,y,seas,,,f, drop = FALSE] # expected
          Obs <- ObsFishLenComps_discard[,y,seas,,,f, drop = FALSE] # observed
          Comp_Type <- FishLenComps_discard_Type[y,f] # composition type

          # get compositions
          tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                    age_or_len = 1, AgeingError = NA)
          # Input into storage
          Obs_FishLen_discard[,y,seas,,,f] <- tmp_comps$Obs
          Pred_FishLen_discard[,y,seas,,,f] <- tmp_comps$Exp
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

        if(length(use_regions) > 0) {
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

        if(length(use_regions) > 0) {
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

  # Pop Fishery Ages
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        for(seas in 1:n_seas) {
          use_regions <- which(UseFishAgeComps_pop[p,,y,seas,f] == 1)
          if(length(use_regions) > 0) {
            Exp <- array(CAA_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_fish_ages, n_sexes, 1))
            Obs <- array(ObsFishAgeComps_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_fish_ages, n_sexes, 1))
            tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = FishAge_pop_CompType[y,f],
                                      age_or_len = 0, AgeingError = AgeingError_t[y,,])
            Obs_FishAge_pop[p,,y,seas,,,f]  <- tmp_comps$Obs
            Pred_FishAge_pop[p,,y,seas,,,f] <- tmp_comps$Exp
          }
        }
      }
    }
  }

  # Pop Fishery Lengths
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        for(seas in 1:n_seas) {
          use_regions <- which(UseFishLenComps_pop[p,,y,seas,f] == 1)
          if(length(use_regions) > 0) {
            Exp <- array(CAL_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_lens, n_sexes, 1))
            Obs <- array(ObsFishLenComps_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_lens, n_sexes, 1))
            tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = FishLen_pop_CompType[y,f],
                                      age_or_len = 1, AgeingError = NA)
            Obs_FishLen_pop[p,,y,seas,,,f]  <- tmp_comps$Obs
            Pred_FishLen_pop[p,,y,seas,,,f] <- tmp_comps$Exp
          }
        }
      }
    }
  }

  # Pop Fishery Ages Discards
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        for(seas in 1:n_seas) {
          use_regions <- which(UseFishAgeComps_discard_pop[p,,y,seas,f] == 1)
          if(length(use_regions) > 0) {
            Exp <- array(DAA_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_fish_ages, n_sexes, 1))
            Obs <- array(ObsFishAgeComps_discard_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_fish_ages, n_sexes, 1))
            tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = FishAgeComps_discard_pop_Type[y,f],
                                      age_or_len = 0, AgeingError = AgeingError_t[y,,])
            Obs_FishAge_discard_pop[p,,y,seas,,,f]  <- tmp_comps$Obs
            Pred_FishAge_discard_pop[p,,y,seas,,,f] <- tmp_comps$Exp
          }
        }
      }
    }
  }

  # Pop Fishery Lengths Discards
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        for(seas in 1:n_seas) {
          use_regions <- which(UseFishLenComps_discard_pop[p,,y,seas,f] == 1)
          if(length(use_regions) > 0) {
            Exp <- array(DAL_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_lens, n_sexes, 1))
            Obs <- array(ObsFishLenComps_discard_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_lens, n_sexes, 1))
            tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = FishLenComps_discard_pop_Type[y,f],
                                      age_or_len = 1, AgeingError = NA)
            Obs_FishLen_discard_pop[p,,y,seas,,,f]  <- tmp_comps$Obs
            Pred_FishLen_discard_pop[p,,y,seas,,,f] <- tmp_comps$Exp
          }
        }
      }
    }
  }

  # Pop Survey Ages
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_srv_fleets) {
        for(seas in 1:n_seas) {
          use_regions <- which(UseSrvAgeComps_pop[p,,y,seas,f] == 1)
          if(length(use_regions) > 0) {
            Exp <- array(SrvIAA_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_srv_ages, n_sexes, 1))
            Obs <- array(ObsSrvAgeComps_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_srv_ages, n_sexes, 1))
            tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = SrvAge_pop_CompType[y,f],
                                      age_or_len = 0, AgeingError = AgeingError_t[y,,])
            Obs_SrvAge_pop[p,,y,seas,,,f]  <- tmp_comps$Obs
            Pred_SrvAge_pop[p,,y,seas,,,f] <- tmp_comps$Exp
          }
        }
      }
    }
  }

  # Pop Survey Lengths
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_srv_fleets) {
        for(seas in 1:n_seas) {
          use_regions <- which(UseSrvLenComps_pop[p,,y,seas,f] == 1)
          if(length(use_regions) > 0) {
            Exp <- array(SrvIAL_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_lens, n_sexes, 1))
            Obs <- array(ObsSrvLenComps_pop[p,,y,seas,,,f], dim = c(n_regions, 1, 1, n_lens, n_sexes, 1))
            tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = SrvLen_pop_CompType[y,f],
                                      age_or_len = 1, AgeingError = NA)
            Obs_SrvLen_pop[p,,y,seas,,,f]  <- tmp_comps$Obs
            Pred_SrvLen_pop[p,,y,seas,,,f] <- tmp_comps$Exp
          }
        }
      }
    }
  }


  # Process outputs
  all_fishages <- reshape2::melt(Obs_FishAge) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishAge) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Age = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Fishery Ages', Pop = NA)

  # Fish Lengths
  all_fishlens <- reshape2::melt(Obs_FishLen) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishLen) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Len = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Fishery Lengths', Pop = NA)

  # Fish Age Discard
  all_fishages_discard <- reshape2::melt(Obs_FishAge_discard) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishAge_discard) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Age = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Fishery Discard Ages', Pop = NA)

  # Fish Lengths Discard
  all_fishlens_discard <- reshape2::melt(Obs_FishLen_discard) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishLen_discard) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Len = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Fishery Discard Lengths', Pop = NA)

  # Survey Ages
  all_srvages <- reshape2::melt(Obs_SrvAge) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvAge) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Age = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Survey Ages', Pop = NA)

  # Survey Lengths
  all_srvlens <- reshape2::melt(Obs_SrvLen) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvLen) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Len = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(Type = 'Survey Lengths', Pop = NA)

  # Pop Fishery Ages
  all_fishages_pop <- reshape2::melt(Obs_FishAge_pop) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishAge_pop) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7")) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Age = Var5, Sex = Var6, Fleet = Var7) %>%
    dplyr::mutate(Type = 'Pop Fishery Ages')

  # Pop Fishery Lengths
  all_fishlens_pop <- reshape2::melt(Obs_FishLen_pop) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishLen_pop) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7")) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Len = Var5, Sex = Var6, Fleet = Var7) %>%
    dplyr::mutate(Type = 'Pop Fishery Lengths')

  # Pop Fishery Ages Discard
  all_fishages_discard_pop <- reshape2::melt(Obs_FishAge_discard_pop) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishAge_discard_pop) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7")) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Age = Var5, Sex = Var6, Fleet = Var7) %>%
    dplyr::mutate(Type = 'Pop Fishery Discard Ages')

  # Pop Fishery Lengths Discard
  all_fishlens_discard_pop <- reshape2::melt(Obs_FishLen_discard_pop) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishLen_discard_pop) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7")) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Len = Var5, Sex = Var6, Fleet = Var7) %>%
    dplyr::mutate(Type = 'Pop Fishery Discard Lengths')

  # Pop Survey Ages
  all_srvages_pop <- reshape2::melt(Obs_SrvAge_pop) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvAge_pop) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7")) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Age = Var5, Sex = Var6, Fleet = Var7) %>%
    dplyr::mutate(Type = 'Pop Survey Ages')

  # Pop Survey Lengths
  all_srvlens_pop <- reshape2::melt(Obs_SrvLen_pop) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvLen_pop) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7")) %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Len = Var5, Sex = Var6, Fleet = Var7) %>%
    dplyr::mutate(Type = 'Pop Survey Lengths')

  return(list(Fishery_Ages     = all_fishages,
              Fishery_Lens     = all_fishlens,
              Fishery_Discard_Ages     = all_fishages_discard,
              Fishery_Discard_Lens     = all_fishlens_discard,
              Survey_Ages      = all_srvages,
              Survey_Lens      = all_srvlens,
              Pop_Fishery_Ages = all_fishages_pop,
              Pop_Fishery_Lens = all_fishlens_pop,
              Pop_Fishery_Discard_Ages = all_fishages_discard_pop,
              Pop_Fishery_Discard_Lens = all_fishlens_discard_pop,
              Pop_Survey_Ages  = all_srvages_pop,
              Pop_Survey_Lens  = all_srvlens_pop,

              Obs_FishAge_mat      = Obs_FishAge,
              Obs_FishLen_mat      = Obs_FishLen,
              Obs_FishAge_discard_mat      = Obs_FishAge_discard,
              Obs_FishLen_discard_mat      = Obs_FishLen_discard,
              Obs_SrvAge_mat       = Obs_SrvAge,
              Obs_SrvLen_mat       = Obs_SrvLen,
              Obs_FishAge_pop_mat  = Obs_FishAge_pop,
              Obs_FishLen_pop_mat  = Obs_FishLen_pop,
              Obs_FishAge_discard_pop_mat  = Obs_FishAge_discard_pop,
              Obs_FishLen_discard_pop_mat  = Obs_FishLen_discard_pop,
              Obs_SrvAge_pop_mat   = Obs_SrvAge_pop,
              Obs_SrvLen_pop_mat   = Obs_SrvLen_pop,

              Pred_FishAge_mat     = Pred_FishAge,
              Pred_FishLen_mat     = Pred_FishLen,
              Pred_FishAge_discard_mat     = Pred_FishAge_discard,
              Pred_FishLen_discard_mat     = Pred_FishLen_discard,
              Pred_SrvAge_mat      = Pred_SrvAge,
              Pred_SrvLen_mat      = Pred_SrvLen,
              Pred_FishAge_pop_mat = Pred_FishAge_pop,
              Pred_FishLen_pop_mat = Pred_FishLen_pop,
              Pred_FishAge_discard_pop_mat = Pred_FishAge_discard_pop,
              Pred_FishLen_discard_pop_mat = Pred_FishLen_discard_pop,
              Pred_SrvAge_pop_mat  = Pred_SrvAge_pop,
              Pred_SrvLen_pop_mat  = Pred_SrvLen_pop
              ))

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
run_external_comp_osa <- function(obs,
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

#' Map a composition data source to its internal-OSA field names
#'
#' Translates a composition data source identifier into the exact field names
#' used in the model \code{data} list and the RTMB-tracked OSA vector name,
#' following the naming convention used throughout \code{SPoRC_rtmb.R} (e.g.
#' \code{ObsFishAgeComps}, \code{ISS_FishAgeComps}, \code{FishAgeComps_Type},
#' \code{ObsFishAgeComps_osa_discrete}).
#'
#' @param comp_source One of \code{"FishAge"}, \code{"FishLen"}, \code{"SrvAge"}, \code{"SrvLen"}.
#' @param pop Logical; population-specific composition source.
#' @param discard Logical; discard composition source (only valid for Fish* sources).
#'
#' @return A named list of field names: \code{Obs}, \code{ISS}, \code{Wt}, \code{Use},
#'   \code{Type}, \code{LikeType}, \code{n_fleets_field}.
#' @keywords internal
comp_osa_field_map <- function(comp_source, pop = FALSE, discard = FALSE) {
  valid_sources <- c("FishAge", "FishLen", "SrvAge", "SrvLen")
  if(!comp_source %in% valid_sources) {
    stop("`comp_source` must be one of: ", paste(valid_sources, collapse = ", "))
  }
  if(discard && grepl("^Srv", comp_source)) {
    stop("`discard = TRUE` is only valid for Fish* sources (FishAge, FishLen).")
  }
  suffix <- paste0(if(discard) "_discard" else "", if(pop) "_pop" else "")
  list(
    Obs      = paste0("Obs", comp_source, "Comps", suffix),
    ISS      = paste0("ISS_", comp_source, "Comps", suffix),
    Wt       = paste0("Wt_", comp_source, "Comps", suffix),
    Use      = paste0("Use", comp_source, "Comps", suffix),
    Type     = paste0(comp_source, "Comps", suffix, "_Type"),
    LikeType = paste0(comp_source, "Comps", suffix, "_LikeType"),
    n_fleets_field = if(grepl("^Srv", comp_source)) "n_srv_fleets" else "n_fish_fleets"
  )
}

#' Keep-subset for internal OSA residuals
#'
#' Elements flagged \code{last_in_group == TRUE} are the statistically
#' determined/reference cell of their group (excluded from the discrete OSA
#' evaluation).
#'
#' @param last_in_group Logical vector (with possible NAs), as produced by
#'   \code{pack_comp_osa(..., return_labels = TRUE)} or
#'   \code{pack_tag_osa(..., return_labels = TRUE)}.
#' @return Integer vector of positions to keep.
#' @keywords internal
osa_keep_subset <- function(last_in_group) {
  which(is.na(last_in_group) | !last_in_group)
}

#' Validate an internal-OSA \code{method}
#'
#' The internal OSA path deliberately restricts \code{RTMB::oneStepPredict()}'s
#' \code{method} to the generic/Gaussian family. In particular the \code{"cdf"}
#' method is disallowed: it is numerically fragile for the discrete
#' (multinomial / count) likelihoods used here and can silently return
#' mis-calibrated residuals, so only \code{"oneStepGeneric"},
#' \code{"oneStepGaussian"}, and \code{"oneStepGaussianOffMode"} are accepted.
#'
#' @param method Character scalar method name.
#' @return \code{method} invisibly, if valid; otherwise an error is raised.
#' @keywords internal
validate_osa_method <- function(method) {
  allowed <- c("oneStepGeneric", "oneStepGaussianOffMode", "oneStepGaussian")
  if(!method %in% allowed) {
    stop("OSA method '", method, "' is invalid! Must be one of: ",
         paste(allowed, collapse = ", "),
         " (the 'cdf' method is not permitted for internal OSA residuals).")
  }
  invisible(method)
}

#' Run internal (model-based) OSA residuals for a composition data source
#'
#' Internal counterpart to \code{\link{run_external_comp_osa}}'s external (post-hoc,
#' compResidual-based) path, called by \code{\link{get_osa}} when a fitted
#' \code{model} is supplied. Calls \code{RTMB::oneStepPredict()} directly on
#' the model's internally tracked OSA vector (built via
#' \code{do_internal_comp_osa = TRUE} in \code{\link{Setup_Mod_Dim}}), and
#' relabels the resulting residuals using \code{\link{pack_comp_osa}}'s
#' \code{return_labels = TRUE} output so the result matches the same
#' \code{res} schema produced by the external path.
#'
#' @param model A fitted RTMB model object from \code{\link{fit_model}}, built
#'   with \code{do_internal_comp_osa = TRUE}.
#' @param data The model \code{data} list (e.g. \code{input_list$data}) used
#'   to build \code{model}.
#' @param comp_source One of \code{"FishAge"}, \code{"FishLen"}, \code{"SrvAge"}, \code{"SrvLen"}.
#' @param family Character, \code{"discrete"} or \code{"continuous"}.
#' @param pop Logical; population-specific composition source. Default \code{FALSE}.
#' @param discard Logical; discard composition source. Default \code{FALSE}.
#' @param bins Vector of age or length bin labels for display.
#' @param bin_label Character label describing whether bins represent ages or lengths.
#' @param parallel Whether or not to parallelize OSA computation. Defaults to \code{FALSE}.
#' @param osa_method Optional override for \code{RTMB::oneStepPredict}'s \code{method}.
#'   Must be one of \code{"oneStepGeneric"}, \code{"oneStepGaussianOffMode"}, or
#'   \code{"oneStepGaussian"}; the \code{"cdf"} method is not permitted (it is
#'   numerically fragile for the discrete likelihoods used here). Defaults to
#'   \code{"oneStepGeneric"} for discrete data and \code{"oneStepGaussianOffMode"}
#'   for continuous (logistic-normal) data. See
#'   \code{\link[TMB:oneStepPredict]{TMB::oneStepPredict}} for further details.
#'   Note that if data are discrete, the only valid option is \code{"oneStepGeneric"}.
#'
#' @return A list with one element \code{res}, matching \code{\link{get_osa}}'s
#'   external-mode schema (columns \code{fleet}, \code{index_label}, \code{year},
#'   \code{index}, \code{resid}, \code{region}, \code{sex}, \code{seas},
#'   \code{comp_type}) plus a \code{pop} column (population index; always 1 for
#'   \code{pop = FALSE} sources), or \code{NULL} if no data of the requested
#'   family/source is present.
#' @keywords internal
run_internal_comp_osa <- function(model, data, comp_source, family,
                                  pop = FALSE, discard = FALSE, parallel = FALSE,
                                  bins, bin_label, osa_method = NULL) {

  fm <- comp_osa_field_map(comp_source, pop = pop, discard = discard)
  n_pop <- if(pop) data$n_pop else 1

  packed <- pack_comp_osa(
    ObsArr = data[[fm$Obs]], ISSArr = data[[fm$ISS]], WtArr = data[[fm$Wt]],
    UseArr = data[[fm$Use]], TypeMat = data[[fm$Type]], LikeTypeVec = data[[fm$LikeType]],
    n_yrs = length(data$years), n_seas = data$n_seas, n_fleets = data[[fm$n_fleets_field]],
    n_sexes = data$n_sexes, addtocomp = data$addtocomp, family = family,
    pop = pop, n_pop = n_pop, return_labels = TRUE
  )

  if(is.null(packed)) {
    warning("No '", family, "' family data found for comp_source '", comp_source, "'; returning NULL.")
    return(NULL)
  }

  tracked_name <- paste0(fm$Obs, "_osa_", family)
  discrete <- (family == "discrete")
  method <- if(!is.null(osa_method)) osa_method else if(discrete) "oneStepGeneric" else "oneStepGaussianOffMode"
  validate_osa_method(method)
  subset_idx <- osa_keep_subset(packed$labels$last_in_group)

  osa <- RTMB::oneStepPredict(model, observation.name = tracked_name, method = method,
                              discrete = discrete, parallel = parallel,
                              subset = subset_idx, trace = FALSE,
                              discreteSupport = if(discrete) 0:max(model$env$obs[[tracked_name]]) else NULL)

  lab <- packed$labels[subset_idx, ]
  lab$resid <- osa$residual

  res <- lab %>%
    dplyr::transmute(
      fleet = as.character(fleet),
      index_label = bin_label,
      year = year,
      index = bins[bin],
      resid = resid,
      region = region,
      sex = sex,
      seas = season,
      pop = pop,
      comp_type = dplyr::case_when(
        comp_type == 0 ~ "Aggregated",
        comp_type == 1 ~ "SpltR_SpltS",
        TRUE           ~ "SpltR_JntS"
      )
    )

  list(res = res)
}

#' Run internal (model-based) OSA residuals for conventional tagging data
#'
#' Internal counterpart to \code{\link{run_internal_comp_osa}} for
#' conventional tag recapture data packed via \code{\link{pack_tag_osa}}
#' (requires \code{do_internal_conv_tag_osa = TRUE} in
#' \code{\link{Setup_Mod_Dim}}). Called by \code{\link{get_osa}} when
#' \code{tag = TRUE} and a fitted \code{model} is supplied.
#'
#' @param model A fitted RTMB model object from \code{\link{fit_model}}, built
#'   with \code{do_internal_conv_tag_osa = TRUE}.
#' @param data The model \code{data} list (e.g. \code{input_list$data}) used
#'   to build \code{model}.
#' @param osa_method Optional override for \code{RTMB::oneStepPredict}'s \code{method}.
#'   Must be one of \code{"oneStepGeneric"}, \code{"oneStepGaussianOffMode"}, or
#'   \code{"oneStepGaussian"} (the \code{"cdf"} method is not permitted). Defaults
#'   to \code{"oneStepGeneric"} (tag recapture data are always discrete/count-valued).
#' @param parallel Whether or not to parallelize OSA computation. Defaults to \code{FALSE}.
#'
#' @return A list with one element \code{res}: columns \code{fleet}, \code{region},
#'   \code{pop_pool}, \code{age_pool}, \code{sex_pool}, \code{cohort},
#'   \code{release_year}, \code{release_region}, \code{release_season},
#'   \code{recovery_year}, \code{recovery_season}, \code{years_at_liberty},
#'   \code{is_tail}, \code{resid}, \code{family}, \code{comp_type = "Tag"},
#'   or \code{NULL} if no tagging data is present.
#' @keywords internal
run_internal_tag_osa <- function(model, data, osa_method = NULL, parallel = FALSE) {

  family <- tag_fam_of(data$conv_fish_tag_like)

  packed <- pack_tag_osa(
    family = family, like_type = data$conv_fish_tag_like,
    obs_recap = data$obs_conv_tag_fish_recap, pred_recap = data$obs_conv_tag_fish_recap,
    tagged_fish = data$conv_tagged_fish,
    conv_tag_release_indicator = data$conv_tag_release_indicator,
    conv_tag_max_liberty = data$conv_tag_max_liberty,
    n_conv_tag_cohorts = data$n_conv_tag_cohorts,
    n_yrs = length(data$years), n_seas = data$n_seas, n_regions = data$n_regions,
    n_fish_fleets = data$n_fish_fleets,
    n_pop_pool = length(data$conv_tag_pop_pool), n_age_pool = length(data$conv_tag_age_pool),
    n_sex_pool = length(data$conv_tag_sex_pool),
    pop_pool = data$conv_tag_pop_pool, age_pool = data$conv_tag_age_pool, sex_pool = data$conv_tag_sex_pool,
    use_fish_tagging = data$use_conv_fish_tagging, conv_tag_mixing_period = data$conv_tag_mixing_period,
    addtotag = data$addtotag, return_labels = TRUE
  )

  if(is.null(packed$vec)) {
    warning("No conventional tagging data found; returning NULL.")
    return(NULL)
  }

  tracked_name <- paste0("ObsConvTag_osa_", family)
  method <- if(!is.null(osa_method)) osa_method else "oneStepGeneric"
  validate_osa_method(method)
  subset_idx <- osa_keep_subset(packed$labels$last_in_group)

  osa <- RTMB::oneStepPredict(model, observation.name = tracked_name, method = method,
                              discreteSupport = 0:max(model$env$obs[[tracked_name]]),
                              discrete = TRUE, parallel = parallel, subset = subset_idx, trace = FALSE)

  lab <- packed$labels[subset_idx, ]
  lab$resid <- osa$residual

  res <- lab %>%
    dplyr::mutate(
      fleet = as.character(fleet),
      region = region,
      pop_pool = pop_pool,
      age_pool = age_pool,
      sex_pool = sex_pool,
      cohort = tc,
      release_year = data$years[ty],
      release_region = tr,
      release_season = tseas,
      recovery_year = data$years[ty + ry - 1],
      recovery_season = rseas,
      years_at_liberty = ry,
      is_tail = is_tail,
      resid = resid,
      family = family,
      comp_type = "Tag",
      .keep = 'none'
    )

  list(res = res)
}

#' Map an index-type data source to its internal-OSA field names
#'
#' Translates an index-type data source identifier into the exact field names
#' used in the model \code{data} list, following the naming convention used
#' throughout \code{SPoRC_rtmb.R} (e.g. \code{ObsFishIdx}, \code{UseFishIdx}).
#'
#' @param index_source One of \code{"Catch"}, \code{"Discard"}, \code{"FishIdx"}, \code{"SrvIdx"}.
#' @param pop Logical; population-specific index source.
#'
#' @return A named list of field names: \code{Obs}, \code{Use}.
#' @keywords internal
index_osa_field_map <- function(index_source, pop = FALSE) {
  valid_sources <- c("Catch", "Discard", "FishIdx", "SrvIdx")
  if(!index_source %in% valid_sources) {
    stop("`index_source` must be one of: ", paste(valid_sources, collapse = ", "))
  }
  suffix <- if(pop) "_pop" else ""
  list(
    Obs = paste0("Obs", index_source, suffix),
    Use = paste0("Use", index_source, suffix)
  )
}

#' Run internal (model-based) OSA residuals for an index-type data source
#'
#' Internal counterpart to \code{\link{run_internal_comp_osa}} for
#' catch/discard/index data (\code{ObsCatch}, \code{ObsDiscard},
#' \code{ObsFishIdx}, \code{ObsSrvIdx}, and their \code{_pop} variants). These
#' are always continuous observations.
#'
#' @param model A fitted RTMB model object from \code{\link{fit_model}}.
#' @param data The model \code{data} list (e.g. \code{input_list$data}) used
#'   to build \code{model}.
#' @param index_source One of \code{"Catch"}, \code{"Discard"}, \code{"FishIdx"}, \code{"SrvIdx"}.
#' @param pop Logical; population-specific index source. Default \code{FALSE}.
#' @param osa_method Optional override for \code{RTMB::oneStepPredict}'s
#'   \code{method}. Must be one of \code{"oneStepGeneric"},
#'   \code{"oneStepGaussianOffMode"}, or \code{"oneStepGaussian"} (the
#'   \code{"cdf"} method is not permitted). Defaults to \code{"oneStepGeneric"}.
#' @param parallel Whether or not to parallelize OSA computation. Defaults to \code{FALSE}.
#'
#' @return A list with one element \code{res}: columns \code{fleet},
#'   \code{region}, \code{year}, \code{season}, \code{pop}, \code{resid}, and
#'   \code{idx_type} (set to \code{index_source}; named \code{idx_type} rather
#'   than \code{comp_type} because index-type sources are not compositions), or
#'   \code{NULL} if no data of the requested source is present.
#' @keywords internal
run_internal_index_osa <- function(model, data, index_source, pop = FALSE,
                                   osa_method = NULL, parallel = FALSE) {

  fm <- index_osa_field_map(index_source, pop = pop)
  use_arr <- data[[fm$Use]]

  if(is.null(use_arr) || !any(use_arr == 1, na.rm = TRUE)) {
    warning("No data found for index_source '", index_source, "'", if(pop) " (pop)" else "", "; returning NULL.")
    return(NULL)
  }

  valid_idx <- which(use_arr == 1)
  dim_names <- if(pop) c("pop", "region", "year", "season", "fleet") else c("region", "year", "season", "fleet")
  map <- as.data.frame(arrayInd(valid_idx, dim(use_arr)))
  colnames(map) <- dim_names

  tracked_name <- fm$Obs
  method <- if(!is.null(osa_method)) osa_method else "oneStepGeneric"
  validate_osa_method(method)

  osa <- RTMB::oneStepPredict(model, observation.name = tracked_name, method = method,
                              discrete = FALSE, parallel = parallel, trace = FALSE)

  res <- data.frame(
    fleet = as.character(map$fleet),
    region = map$region,
    year = data$years[map$year],
    season = map$season,
    pop = if(pop) map$pop else 1L,
    resid = osa$residual,
    idx_type = index_source
  )

  list(res = res)
}

#' Compute OSA residuals for composition data
#'
#' Formats observed and expected composition data and calculates one-step-ahead
#' (OSA) residuals using multinomial, Dirichlet-multinomial, or logistic-normal
#' likelihoods. This function is the main interface for residual diagnostics,
#' internally calling [run_external_comp_osa()] to perform the residual calculations.
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
#' @param model A fitted RTMB model object from \code{\link{fit_model}}
#'   (built with \code{do_internal_comp_osa = TRUE} or
#'   \code{do_internal_conv_tag_osa = TRUE}). Supplying \code{model} switches
#'   \code{get_osa()} from the default external (post-hoc, compResidual-based)
#'   path to the internal path, which calls \code{RTMB::oneStepPredict()}
#'   directly on the model's internally tracked OSA vector. All \code{obs_mat}/
#'   \code{exp_mat}/\code{N}/\code{DM_theta}/\code{LN_Sigma}/\code{years}/
#'   \code{comp_type}/\code{comp_like} arguments above are ignored in this mode.
#' @param data The model \code{data} list (e.g. \code{input_list$data}) used to
#'   build \code{model}. Required when \code{model} is supplied.
#' @param comp_source One of \code{"FishAge"}, \code{"FishLen"}, \code{"SrvAge"},
#'   \code{"SrvLen"}, identifying which composition data source to pull
#'   internal OSA residuals for. Required when \code{model} is supplied and
#'   \code{index_source} is \code{NULL} and \code{tag = FALSE}.
#' @param index_source One of \code{"Catch"}, \code{"Discard"}, \code{"FishIdx"},
#'   \code{"SrvIdx"}, identifying which continuous (log-normal) index-type
#'   data source to pull internal OSA residuals for. When supplied, takes
#'   precedence over \code{comp_source}/\code{tag}. Only used when
#'   \code{model} is supplied.
#' @param family Character, \code{"discrete"} or \code{"continuous"}; which of
#'   the two internally-tracked OSA vectors to use for \code{comp_source} (a
#'   source can have both, e.g. some fleets multinomial and others
#'   logistic-normal). Only used when \code{model} is supplied, \code{tag =
#'   FALSE}, and \code{index_source} is \code{NULL}.
#' @param pop Logical; population-specific composition or index source. Only
#'   used when \code{model} is supplied and \code{tag = FALSE}. Default \code{FALSE}.
#' @param discard Logical; discard composition source (only valid for
#'   \code{comp_source \%in\% c("FishAge","FishLen")}). Only used when
#'   \code{model} is supplied, \code{tag = FALSE}, and \code{index_source} is
#'   \code{NULL}. Default \code{FALSE}.
#' @param tag Logical; if \code{TRUE} (and \code{model} is supplied, and
#'   \code{index_source} is \code{NULL}), compute internal OSA residuals for
#'   conventional tag recapture data instead of composition data. Default
#'   \code{FALSE}.
#' @param osa_method Optional override for \code{RTMB::oneStepPredict}'s
#'   \code{method}, used only in internal mode. Must be one of
#'   \code{"oneStepGeneric"}, \code{"oneStepGaussianOffMode"}, or
#'   \code{"oneStepGaussian"}; the \code{"cdf"} method is not permitted (it is
#'   numerically fragile for the discrete likelihoods used here). Defaults to
#'   \code{"oneStepGeneric"} for discrete families/tags and
#'   \code{"oneStepGaussianOffMode"} for continuous (logistic-normal)
#'   composition families.
#' @param parallel Whether or not to parallelize OSA computation in internal
#'   mode. Defaults to \code{FALSE}.
#'
#' @details
#' When computing OSA residuals for population-specific composition data,
#' slice the leading population dimension from the \code{obs_mat} and
#' \code{exp_mat} arrays before passing them to this function. For example,
#' to compute residuals for population \code{p}:
#'
#' \preformatted{
#' get_osa(obs_mat = Obs_FishAge_pop_mat[p,,,,,,],
#'         exp_mat = Pred_FishAge_pop_mat[p,,,,,,],
#'         ...)
#' }
#'
#' Population-specific composition arrays returned by
#' \code{\link{get_comp_prop}} are dimensioned
#' \code{[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fleets]}.
#' Slicing on \code{p} yields a 6D array matching the expected input dimensions.
#'
#' For internal (model-based) OSA residuals, fit the model with
#' \code{do_internal_comp_osa = TRUE} and/or \code{do_internal_conv_tag_osa =
#' TRUE} (set in \code{\link{Setup_Mod_Dim}}), then call, e.g.:
#'
#' \preformatted{
#' get_osa(model = fitted_obj, data = input_list$data, comp_source = "FishAge",
#'         family = "discrete", bins = input_list$data$ages, bin_label = "Age")
#' get_osa(model = fitted_obj, data = input_list$data, tag = TRUE)
#' get_osa(model = fitted_obj, data = input_list$data, index_source = "SrvIdx")
#' }
#'
#' @return A list with one element:
#' \describe{
#'   \item{res}{Data frame of OSA residuals. Columns include:
#'     \code{fleet}, \code{index_label}, \code{year}, \code{index},
#'     \code{resid}, \code{region}, \code{seas}, \code{sex}, and \code{comp_type}
#'     (composition sources, external or internal); \code{fleet}, \code{region},
#'     \code{cohort}, \code{release_year}/\code{release_region}/\code{release_season},
#'     \code{recovery_year}/\code{recovery_season}, \code{years_at_liberty}, \code{resid},
#'     and \code{comp_type = "Tag"} (\code{tag = TRUE}); or \code{fleet}, \code{region},
#'     \code{year}, \code{season}, \code{pop}, \code{resid}, and \code{idx_type}
#'     set to \code{index_source} (\code{index_source} supplied.}
#' }
#'
#' @family Model Diagnostics
#' @import dplyr
#' @export get_osa
get_osa <- function(obs_mat = NULL,
                    exp_mat = NULL,
                    N = NULL,
                    DM_theta = NULL,
                    LN_Sigma = NULL,
                    years = NULL,
                    seas = NULL,
                    fleet = NULL,
                    bins = NULL,
                    comp_type = NULL,
                    bin_label = NULL,
                    comp_like = 0,
                    addtocomp = 0,
                    model = NULL,
                    data = NULL,
                    comp_source = NULL,
                    index_source = NULL,
                    family = "discrete",
                    pop = FALSE,
                    discard = FALSE,
                    tag = FALSE,
                    osa_method = NULL,
                    parallel = FALSE
                    ) {

  # Internal (model-based) OSA path, via RTMB::oneStepPredict
  if(!is.null(model)) {
    if(!is.null(index_source)) {
      return(run_internal_index_osa(model = model, data = data, index_source = index_source,
                                    pop = pop, osa_method = osa_method, parallel = parallel))
    } else if(tag) {
      return(run_internal_tag_osa(model = model, data = data, osa_method = osa_method, parallel = parallel))
    } else {
      return(run_internal_comp_osa(model = model, data = data, comp_source = comp_source, family = family,
                                   pop = pop, discard = discard, bins = bins, bin_label = bin_label,
                                   osa_method = osa_method, parallel = parallel))
    }
  }

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
      tmp_osa <- run_external_comp_osa(obs = tmp_obs, exp = tmp_exp, N = N, DM_theta = DM_theta,
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
          tmp_osa <- run_external_comp_osa(obs = tmp_obs, exp = tmp_exp, N = N[r,years[[r]],s], DM_theta = DM_theta[r,s],
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
        tmp_osa <- run_external_comp_osa(obs = tmp_obs, exp = tmp_exp, N = N[r,years[[r]]],
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
#' QQ-plots with SDNR annotations and bubble plots showing residual magnitude
#' and sign.
#'
#' Panels are faceted by every structural dimension the residual data frame
#' actually spans. Composition plots always facet by \code{region}/\code{sex} and additionally facet by \code{fleet}, \code{pop}, and \code{seas}
#' whenever \code{osa_results$res} contains more than one of each (\code{seas} matters because
#' \code{year} + bin alone don't uniquely place a bubble-plot point when compositions are
#' collected in more than one season). Tagging plots only show QQ plots given the number of dimensions in tagging data.
#' Index-type residuals (from \code{get_osa(..., index_source = ...)}, carrying
#' an \code{idx_type} column \code{\%in\% c("Catch","Discard","FishIdx","SrvIdx")}
#' instead of \code{comp_type}) facet by \code{region},
#' \code{season}, \code{fleet}, and \code{pop} whenever those span more than
#' one level, and pair the QQ-plot with a residual-vs-year point plot instead
#' of a bubble plot (there is no bin/age/length dimension to plot against).
#' Note: these are one-step-ahead residuals; for the simpler raw log-scale
#' (Pearson-style) index residual and the observed-vs-predicted index fit, see
#' \code{\link{get_idx_fits}} / \code{\link{get_idx_fits_plot}} instead.
#'
#' @param osa_results List obtained from get_osa() containing residuals dataframe.
#'
#' @return List of plots: \code{sdnr_plot} (QQ-plot) and a second element
#'   that is a \code{bubble_plot} (composition/tag residual magnitude and
#'   sign) or a \code{resid_plot} (index-type residual vs. year).
#' @export plot_resids
#' @family Model Diagnostics
#' @import dplyr
#' @import ggplot2
plot_resids <- function(osa_results) {

  # extract results
  res <- osa_results$res %>% dplyr::mutate(sign = ifelse(resid < 0, "Neg", "Pos"))

  # Single dispatch discriminator. Composition/tag residuals carry a
  # `comp_type` column; index-type residuals (Catch/Discard/FishIdx/SrvIdx)
  # carry an `idx_type` column instead (they are not compositions).
  res_type <- if("idx_type" %in% names(res)) as.character(unique(res$idx_type)) else as.character(unique(res$comp_type))

  # Which optional structural dimensions does this result actually span?
  has_multi <- function(v) v %in% names(res) && dplyr::n_distinct(res[[v]], na.rm = TRUE) > 1
  multi_fleet <- has_multi("fleet")
  multi_pop   <- has_multi("pop")
  multi_comp_seas <- has_multi("seas")

  # Facet labels for every dimension we might facet on (only used when present).
  lab_fns <- list(
    region          = function(x) paste0("Region ", x),
    sex             = function(x) paste0("Sex ", x),
    fleet           = function(x) paste0("Fleet ", x),
    pop             = function(x) paste0("Pop ", x),
    recovery_season = function(x) paste0("Seas ", x),
    season          = function(x) paste0("Seas ", x),
    seas            = function(x) paste0("Seas ", x),
    pop_pool        = function(x) ifelse(x == "tail", "Tail (non-recap)", paste0("Pool ", x))
  )

  # Build a facet_grid() from character vectors of row/column facet variables
  build_facet <- function(row_vars, col_vars) {
    used <- c(row_vars, col_vars)
    if(length(used) == 0) return(NULL)
    lhs <- if(length(row_vars)) paste(row_vars, collapse = " + ") else "."
    rhs <- if(length(col_vars)) paste(col_vars, collapse = " + ") else "."
    facet_grid(stats::as.formula(paste(lhs, "~", rhs)),
               labeller = do.call(labeller, lab_fns[used]))
  }

  # SDNR annotation table
  sdnr_table <- function(res, grp_vars) {
    grouped <- if(length(grp_vars)) dplyr::group_by(res, dplyr::across(dplyr::all_of(grp_vars))) else res
    grouped %>%
      dplyr::summarise(
        df  = n() - 1,
        HCI = sqrt(qchisq(.975, df) / df),
        LCI = sqrt(qchisq(.025, df) / df),
        est = sd(resid), .groups = "drop") %>%
      dplyr::mutate(
        sdnr = paste0('SDNR=', sprintf('%.2f', est)),
        sdnr = paste0(sdnr, '\n(', sprintf('%.2f', LCI), '-', sprintf('%.2f', HCI), ')')
      )
  }

  qq_base <- function(res, sdnr) {
    ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      theme_bw(base_size = 20) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 4)
  }

  # Extra (non-structural) composition facet columns
  comp_extra_cols <- c(if(multi_fleet) "fleet", if(multi_pop) "pop", if(multi_comp_seas) "seas")

  # Conventional Tagging OSA Residuals ----------------------------------------
  if(res_type == "Tag") {

    multi_region   <- has_multi("region")
    multi_seas     <- has_multi("recovery_season")
    multi_pop_pool <- has_multi("pop_pool")

    if(multi_pop_pool) {
      # the release-conditioned "tail" (non-recap) has no pop_pool of its own
      res <- res %>% dplyr::mutate(pop_pool = ifelse(is.na(pop_pool), "tail", as.character(pop_pool)))
    }
    tag_row_vars <- if(multi_region) "region"
    tag_col_vars <- c(if(multi_seas) "recovery_season", if(multi_fleet) "fleet", if(multi_pop_pool) "pop_pool")

    sdnr <- sdnr_table(res, c(tag_row_vars, tag_col_vars))
    sdnr_plot <- qq_base(res, sdnr) + build_facet(tag_row_vars, tag_col_vars)

    return(list(sdnr_plot))
  }

  # Index-type OSA Residuals (Catch/Discard/FishIdx/SrvIdx) -------------------
  if(res_type %in% c("Catch", "Discard", "FishIdx", "SrvIdx")) {

    multi_region <- has_multi("region")
    multi_seas   <- has_multi("season")

    idx_row_vars <- if(multi_region) "region"
    idx_col_vars <- c(if(multi_seas) "season", if(multi_fleet) "fleet", if(multi_pop) "pop")

    sdnr <- sdnr_table(res, c(idx_row_vars, idx_col_vars))
    sdnr_plot <- qq_base(res, sdnr) + build_facet(idx_row_vars, idx_col_vars)

    resid_plot <- ggplot(data = res, aes(x = year, y = resid, color = sign)) +
      geom_point(size = 2.5) +
      geom_hline(yintercept = 0, lty = 2) +
      scale_color_manual(values = c("blue", "red")) +
      labs(x = "Year", y = "OSA Residual", color = "Sign") +
      theme_bw(base_size = 20) +
      theme(legend.position = 'top') +
      build_facet(idx_row_vars, idx_col_vars)

    return(list(sdnr_plot, resid_plot))
  }

  # Aggregated Comps ----------------------------------------------------------
  if(res_type == "Aggregated") {
    sdnr <- sdnr_table(res, comp_extra_cols)
    sdnr_plot <- qq_base(res, sdnr) + build_facet(character(0), comp_extra_cols)
  }

  # Split Sex and Split Region ------------------------------------------------
  if(res_type == "SpltR_SpltS") {
    grp <- c("region", "sex", comp_extra_cols)
    sdnr <- sdnr_table(res, grp)
    sdnr_plot <- qq_base(res, sdnr) + build_facet("region", c("sex", comp_extra_cols))
  }

  # Joint Sex and Split Region ------------------------------------------------
  if(res_type == "SpltR_JntS") {
    grp <- c("region", comp_extra_cols)
    sdnr <- sdnr_table(res, grp)
    sdnr_plot <- qq_base(res, sdnr) + build_facet(character(0), c("region", comp_extra_cols))
  }

  # bubble plot (shared across composition comp types) ------------------------
  bubble_plot <- ggplot(data = res, aes(x = year, y = as.numeric(index),
                                        color = sign, size = abs(resid), alpha = abs(resid))) +
    geom_point() +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "Year", y = unique(res$index_label), color = "Sign", size = "abs(Resid)", alpha = "abs(Resid)") +
    theme_bw(base_size = 20) +
    theme(legend.position = 'top') +
    build_facet("region", c("sex", comp_extra_cols))

  return(list(sdnr_plot, bubble_plot))

}


