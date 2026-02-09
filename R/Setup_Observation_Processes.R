#' Set up SPoRC model weighting
#'
#' @param input_list List containing data, parameter, and map lists.
#' @param Wt_Catch Either a numeric scalar (lambda) applied to the overall catch dataset or an array of lambdas dimensioned by n_regions, n_years, n_seas, n_fish_fleets.
#' @param Wt_FishIdx Either a numeric scalar (lambda) applied to the overall fishery index dataset  or an array of lambdas dimensioned by n_regions, n_years, n_seas, n_fish_fleets.
#' @param Wt_SrvIdx Either a numeric scalar (lambda) applied to the overall survey index dataset or an array of lambdas dimensioned by n_regions, n_years, n_seas, n_srv_fleets.
#' @param Wt_Rec Numeric weight (lambda) applied to the recruitment penalty.
#' @param Wt_F Numeric weight (lambda) applied to the fishing mortality penalty.
#' @param Wt_FishAgeComps Numeric weight (lambda) applied to fishery age composition data.
#' @param Wt_SrvAgeComps Numeric weight (lambda) applied to survey age composition data.
#' @param Wt_FishLenComps Numeric weight (lambda) applied to fishery length composition data.
#' @param Wt_SrvLenComps Numeric weight (lambda) applied to survey length composition data.
#' @param Wt_Tagging Numeric weight (lambda) applied to tagging data.
#'
#' @export Setup_Mod_Weighting
#' @family Model Setup
Setup_Mod_Weighting <- function(input_list,
                                Wt_Catch = 1,
                                Wt_FishIdx = 1,
                                Wt_SrvIdx = 1,
                                Wt_Rec = 1,
                                Wt_F = 1,
                                Wt_Tagging = 1,
                                Wt_FishAgeComps,
                                Wt_SrvAgeComps,
                                Wt_FishLenComps,
                                Wt_SrvLenComps
                                ) {

  messages_list <<- character(0) # string to attach to for printing messages

  input_list$data$Wt_Catch <- Wt_Catch
  input_list$data$Wt_FishIdx <- Wt_FishIdx
  input_list$data$Wt_SrvIdx <- Wt_SrvIdx
  input_list$data$Wt_Rec <- Wt_Rec
  input_list$data$Wt_F <- Wt_F
  input_list$data$Wt_FishAgeComps<- Wt_FishAgeComps
  input_list$data$Wt_SrvAgeComps<- Wt_SrvAgeComps
  input_list$data$Wt_FishLenComps<- Wt_FishLenComps
  input_list$data$Wt_SrvLenComps<- Wt_SrvLenComps
  input_list$data$Wt_Tagging <- Wt_Tagging

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
