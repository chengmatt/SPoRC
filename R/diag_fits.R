# Stage 3 of 3: post fit
#
# Observed against predicted, as tidy data frames: index fits and composition
# proportions. The extraction layer underneath the plots in
# plot_figures_tables.R and the reweighting in diag_francis.R.

#' Extract Index Fit Results
#'
#' Generates a tidy dataframe of observed and predicted survey and fishery
#' indices from a fitted RTMB model, including standard errors, confidence
#' intervals, a raw log-scale (Pearson-style) residual, and catchability
#' blocks. Both pooled and population-specific indices are returned when the
#' corresponding \code{Use*_pop} flags contain any ones.
#'
#' The \code{resid} column here is the simple log-scale residual
#' \eqn{\log(\text{obs}) - \log(\text{predicted})}, \emph{not} a
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

  # Pop-specific dimnames, year is dim 3
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
#'   \item \code{Exp}: array of expected composition values in observed bins
#'   \item \code{Obs}: array of observed composition values in observed bins
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
