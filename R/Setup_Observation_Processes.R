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
#' @param Wt_Catch_pop Weight applied to the population-specific catch
#'   likelihood. Either a scalar applied uniformly or a numeric array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}. Default \code{0}.
#' @param Wt_FishIdx_pop Weight applied to the population-specific fishery
#'   index likelihood. Same scalar or array format as \code{Wt_Catch_pop},
#'   dimensioned \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   Default \code{0}.
#' @param Wt_SrvIdx_pop Weight applied to the population-specific survey index
#'   likelihood. Same scalar or array format as \code{Wt_Catch_pop}, dimensioned
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}. Default \code{0}.
#' @param Wt_FishAgeComps Weight applied to the fishery age composition
#'   likelihood. Either a scalar or a numeric array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_SrvAgeComps Weight applied to the survey age composition
#'   likelihood. Either a scalar or a numeric array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_FishLenComps Weight applied to the fishery length composition
#'   likelihood. Same format as \code{Wt_FishAgeComps},
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{0}s.
#' @param Wt_SrvLenComps Weight applied to the survey length composition
#'   likelihood. Same format as \code{Wt_SrvAgeComps},
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   Default: array of \code{0}s.
#' @param Wt_FishAgeComps_pop Weight applied to the population-specific fishery
#'   age composition likelihood. Either a scalar or a numeric array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{0}s.
#' @param Wt_SrvAgeComps_pop Weight applied to the population-specific survey
#'   age composition likelihood. Either a scalar or a numeric array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   Default: array of \code{0}s.
#' @param Wt_FishLenComps_pop Weight applied to the population-specific fishery
#'   length composition likelihood. Same format as \code{Wt_FishAgeComps_pop},
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{0}s.
#' @param Wt_SrvLenComps_pop Weight applied to the population-specific survey
#'   length composition likelihood. Same format as \code{Wt_SrvAgeComps_pop},
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   Default: array of \code{0}s.
#' @param Wt_Rec Scalar weight applied to the recruitment deviation penalty
#'   (\code{ln_RecDevs}). Default \code{1}.
#' @param Wt_F Scalar weight applied to the fishing mortality deviation penalty
#'   (\code{ln_F_devs}). Default \code{1}.
#' @param Wt_Tagging Scalar weight applied to the tag-recovery likelihood.
#'   Default \code{1}.
#'
#' @return The input \code{input_list} with all weight values stored in
#'   \code{$data} under their respective names (\code{Wt_Catch},
#'   \code{Wt_FishIdx}, \code{Wt_SrvIdx}, \code{Wt_Catch_pop},
#'   \code{Wt_FishIdx_pop}, \code{Wt_SrvIdx_pop}, \code{Wt_FishAgeComps},
#'   \code{Wt_SrvAgeComps}, \code{Wt_FishLenComps}, \code{Wt_SrvLenComps},
#'   \code{Wt_FishAgeComps_pop}, \code{Wt_SrvAgeComps_pop},
#'   \code{Wt_FishLenComps_pop}, \code{Wt_SrvLenComps_pop},
#'   \code{Wt_Rec}, \code{Wt_F}, \code{Wt_Tagging}).
#'
#' @export Setup_Mod_Weighting
#' @family Model Setup
Setup_Mod_Weighting <- function(input_list,
                                Wt_Catch = 1,
                                Wt_FishIdx = 1,
                                Wt_SrvIdx = 1,
                                Wt_Catch_pop = 0,
                                Wt_FishIdx_pop = 0,
                                Wt_SrvIdx_pop = 0,
                                Wt_Rec = 1,
                                Wt_F = 1,
                                Wt_Tagging = 1,
                                Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                Wt_FishLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                Wt_FishAgeComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvAgeComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                Wt_FishLenComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvLenComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
                                ) {

  messages_list <<- character(0) # string to attach to for printing messages

  input_list$data$Wt_Catch <- Wt_Catch
  input_list$data$Wt_FishIdx <- Wt_FishIdx
  input_list$data$Wt_SrvIdx <- Wt_SrvIdx
  input_list$data$Wt_Catch_pop <- Wt_Catch_pop
  input_list$data$Wt_FishIdx_pop <- Wt_FishIdx_pop
  input_list$data$Wt_SrvIdx_pop <- Wt_SrvIdx_pop
  input_list$data$Wt_Rec <- Wt_Rec
  input_list$data$Wt_F <- Wt_F
  input_list$data$Wt_FishAgeComps<- Wt_FishAgeComps
  input_list$data$Wt_SrvAgeComps<- Wt_SrvAgeComps
  input_list$data$Wt_FishLenComps<- Wt_FishLenComps
  input_list$data$Wt_SrvLenComps<- Wt_SrvLenComps
  input_list$data$Wt_FishAgeComps_pop <- Wt_FishAgeComps_pop
  input_list$data$Wt_SrvAgeComps_pop <- Wt_SrvAgeComps_pop
  input_list$data$Wt_FishLenComps_pop <- Wt_FishLenComps_pop
  input_list$data$Wt_SrvLenComps_pop <- Wt_SrvLenComps_pop
  input_list$data$Wt_Tagging <- Wt_Tagging

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
