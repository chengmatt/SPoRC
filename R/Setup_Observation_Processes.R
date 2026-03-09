#' Set likelihood and penalty weights for the estimation model
#'
#' Assigns lambda (\eqn{\lambda}) multipliers to each likelihood component and
#' penalty term in the TMB/RTMB objective function. Weights scale the relative
#' contribution of each data source during estimation and can be used to
#' down-weight noisy data, implement iterative reweighting schemes (e.g.,
#' Francis or McAllister–Ianelli), or disable a component entirely by setting
#' its weight to \code{0}. Must be called after all data setup functions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param Wt_Catch Weight applied to the catch likelihood. Either a scalar
#'   applied uniformly across all fleets, regions, years, and seasons, or a
#'   numeric array \code{[n_regions × n_years × n_seas × n_fish_fleets]} for
#'   fleet- or time-specific weighting. Default \code{1}.
#' @param Wt_FishIdx Weight applied to the fishery index likelihood. Accepts
#'   the same scalar or array format as \code{Wt_Catch}, dimensioned
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}. Default \code{1}.
#' @param Wt_SrvIdx Weight applied to the survey index likelihood. Accepts
#'   the same scalar or array format, dimensioned
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. Default \code{1}.
#' @param Wt_FishAgeComps Scalar weight applied to the fishery age composition
#'   likelihood across all fleets. No default; must be provided.
#' @param Wt_SrvAgeComps Scalar weight applied to the survey age composition
#'   likelihood across all fleets. No default; must be provided.
#' @param Wt_FishLenComps Scalar weight applied to the fishery length
#'   composition likelihood across all fleets. No default; must be provided.
#' @param Wt_SrvLenComps Scalar weight applied to the survey length composition
#'   likelihood across all fleets. No default; must be provided.
#' @param Wt_Rec Scalar weight applied to the recruitment deviation penalty
#'   (\code{ln_RecDevs}). Default \code{1}.
#' @param Wt_F Scalar weight applied to the fishing mortality deviation penalty
#'   (\code{ln_F_devs}). Default \code{1}.
#' @param Wt_Tagging Scalar weight applied to the tag-recovery likelihood.
#'   Default \code{1}.
#'
#' @return The input \code{input_list} with all weight values stored in
#'   \code{$data} under their respective names (\code{Wt_Catch},
#'   \code{Wt_FishIdx}, \code{Wt_SrvIdx}, \code{Wt_FishAgeComps},
#'   \code{Wt_SrvAgeComps}, \code{Wt_FishLenComps}, \code{Wt_SrvLenComps},
#'   \code{Wt_Rec}, \code{Wt_F}, \code{Wt_Tagging}).
#'
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
