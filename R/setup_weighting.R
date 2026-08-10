# Stage 1 of 3: model setup
#
# Relative weighting across data sources, applied as multipliers on each
# likelihood component. Runs last in the setup chain, once every data source it
# weights has been declared.

#' Set likelihood and penalty weights for the estimation model
#'
#' Assigns lambda (\eqn{\lambda}) multipliers to each likelihood component and
#' penalty term in the TMB/RTMB objective function. Weights scale the relative
#' contribution of each data source during estimation and can be used to
#' down-weight noisy data, implement iterative reweighting schemes (e.g.,
#' Francis), or disable a component entirely by setting
#' its weight to \code{1}. Must be called after all data setup functions.
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
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}. Default \code{1}.
#' @param Wt_FishIdx_pop Weight applied to the population-specific fishery
#'   index likelihood. Same scalar or array format as \code{Wt_Catch_pop},
#'   dimensioned \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   Default \code{1}.
#' @param Wt_SrvIdx_pop Weight applied to the population-specific survey index
#'   likelihood. Same scalar or array format as \code{Wt_Catch_pop}, dimensioned
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}. Default \code{1}.
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
#'   Default: array of \code{1}s.
#' @param Wt_SrvLenComps Weight applied to the survey length composition
#'   likelihood. Same format as \code{Wt_SrvAgeComps},
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_FishAgeComps_pop Weight applied to the population-specific fishery
#'   age composition likelihood. Either a scalar or a numeric array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_SrvAgeComps_pop Weight applied to the population-specific survey
#'   age composition likelihood. Either a scalar or a numeric array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_FishLenComps_pop Weight applied to the population-specific fishery
#'   length composition likelihood. Same format as \code{Wt_FishAgeComps_pop},
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_SrvLenComps_pop Weight applied to the population-specific survey
#'   length composition likelihood. Same format as \code{Wt_SrvAgeComps_pop},
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_Rec Weight applied to the recruitment deviation penalty
#'   (\code{ln_RecDevs}). Either a scalar applied uniformly or a numeric array
#'   \code{[n_pop × n_regions × n_est_rec_devs]} for deviation-specific
#'   weighting, where \code{n_est_rec_devs} is the third dimension of
#'   \code{ln_RecDevs} rather than the number of years, since
#'   \code{dont_est_recdev_last} and \code{n_proj_yrs_devs} both move it.
#'   Default \code{1}. A weight of zero on a deviation leaves it estimated but
#'   removes it from the penalty entirely, which is how a stock-recruit
#'   relationship is fit over a window of years while recruitment stays free in
#'   every year. That is distinct from \code{dont_est_recdev_last}, which
#'   removes the deviations themselves so recruitment reverts to the
#'   deterministic prediction in those years.
#' @param Wt_Init_Rec Weight applied to the initial age deviation penalty
#'   (\code{ln_InitDevs}). Either a scalar or a numeric array
#'   \code{[n_pop × n_regions × (n_ages - 1)]}. Defaults to \code{NULL}, which
#'   takes whatever \code{Wt_Rec} is when \code{Wt_Rec} is a scalar; supply it
#'   explicitly when \code{Wt_Rec} is an array, since the two penalties are
#'   dimensioned differently.
#' @param Wt_F Scalar weight applied to the fishing mortality deviation penalty
#'   (\code{ln_F_devs}). Default \code{1}.
#' @param Wt_Tagging Scalar weight applied to the tag-recovery likelihood.
#'   Default \code{1}.
#' @param Wt_Discard Weight applied to the aggregated discard amount or
#'   fraction likelihood. Either a scalar applied uniformly or a numeric
#'   array \code{[n_regions × n_years × n_seas × n_fish_fleets]}. Default
#'   \code{1}.
#' @param Wt_Discard_pop Weight applied to the population-specific discard
#'   likelihood. Either a scalar or a numeric array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}. Default
#'   \code{1}.
#' @param Wt_D Scalar weight applied to the discard mortality rate
#'   deviation penalty (\code{logit_dmr_devs}). Default \code{1}.
#' @param Wt_FishAgeComps_discard Weight applied to the discard fishery age
#'   composition likelihood. Either a scalar or a numeric array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_FishLenComps_discard Weight applied to the discard fishery
#'   length composition likelihood. Same format as
#'   \code{Wt_FishAgeComps_discard},
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_FishAgeComps_discard_pop Weight applied to the
#'   population-specific discard fishery age composition likelihood. Either
#'   a scalar or a numeric array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{1}s.
#' @param Wt_FishLenComps_discard_pop Weight applied to the
#'   population-specific discard fishery length composition likelihood.
#'   Same format as \code{Wt_FishAgeComps_discard_pop},
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   Default: array of \code{1}s.
#' @param fish_sel_pen_wts \code{NULL} (default), or a named numeric vector/list
#'   with independent weights for any subset of six selectivity smoothness
#'   penalty terms (see \code{\link{resolve_sel_pen_wts}} and
#'   \code{\link{Get_Selex_Smoothness_Penalty}}), evaluated directly on the
#'   fleet's realized selectivity-at-bin-at-year surface and so applicable to
#'   any selectivity functional form:
#'   \describe{
#'     \item{\code{"smooth_bin_curve"}}{Second-difference (curvature) penalty across bins.}
#'     \item{\code{"smooth_bin_diff"}}{Unconditional first-difference penalty across bins.}
#'     \item{\code{"smooth_yr_diff"}}{First-difference penalty across years.}
#'     \item{\code{"smooth_yr_curve"}}{Second-difference (curvature) penalty across years.}
#'     \item{\code{"smooth_dome"}}{Dome-shape (non-monotonicity) penalty across bins.}
#'     \item{\code{"smooth_mean_center"}}{Per-year mean-centering regularization.}
#'   }
#'   Any name not supplied defaults to \code{0} (off). Must be called after
#'   \code{Setup_Mod_Fishsel_and_Q}.
#'   Each weight may instead be a vector with one value per model year, so a
#'   penalty can act only in some years or with a different strength in each.
#'   The specification may also carry \code{"bin_range"}, a length-two vector
#'   giving the first and last bin the penalties act over. To give each fleet
#'   its own penalties, pass an unnamed list with one named specification per
#'   fleet instead of a single specification.
#' @param ret_sel_pen_wts Same format as \code{fish_sel_pen_wts}, for the
#'   retained fishery selectivity penalty.
#' @param srv_sel_pen_wts Same format as \code{fish_sel_pen_wts}, for the
#'   survey selectivity penalty. Must be called
#'   after \code{Setup_Mod_Srvsel_and_Q}.
#'
#' @return The input \code{input_list} with all weight values stored in
#'   \code{$data} under their respective names (\code{Wt_Catch},
#'   \code{Wt_FishIdx}, \code{Wt_SrvIdx}, \code{Wt_Catch_pop},
#'   \code{Wt_FishIdx_pop}, \code{Wt_SrvIdx_pop}, \code{Wt_FishAgeComps},
#'   \code{Wt_SrvAgeComps}, \code{Wt_FishLenComps}, \code{Wt_SrvLenComps},
#'   \code{Wt_FishAgeComps_pop}, \code{Wt_SrvAgeComps_pop},
#'   \code{Wt_FishLenComps_pop}, \code{Wt_SrvLenComps_pop},
#'   \code{Wt_Rec}, \code{Wt_F}, \code{Wt_Tagging}, \code{Wt_Discard}, \code{Wt_Discard_pop}, \code{Wt_D},
#'   \code{Wt_FishAgeComps_discard}, \code{Wt_FishLenComps_discard},
#'   \code{Wt_FishAgeComps_discard_pop}, \code{Wt_FishLenComps_discard_pop}).
#'
#' @export Setup_Mod_Weighting
#' @family Model Setup
Setup_Mod_Weighting <- function(input_list,
                                Wt_Catch = 1,
                                Wt_FishIdx = 1,
                                Wt_SrvIdx = 1,
                                Wt_Catch_pop = 1,
                                Wt_FishIdx_pop = 1,
                                Wt_SrvIdx_pop = 1,
                                Wt_Rec = 1,
                                Wt_Init_Rec = NULL,
                                Wt_F = 1,
                                Wt_Tagging = 1,

                                # Retained Catch Stuff
                                Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                Wt_FishAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                Wt_FishLenComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                   input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_SrvLenComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                  input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),

                                # Discard Stuff
                                Wt_Discard = 1,
                                Wt_Discard_pop = 1,
                                Wt_D = 1,
                                Wt_FishAgeComps_discard = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                           input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_FishLenComps_discard = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                           input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_FishAgeComps_discard_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                               input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_FishLenComps_discard_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                                                               input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),

                                # Selectivity penalty weights
                                fish_sel_pen_wts = NULL,
                                ret_sel_pen_wts = NULL,
                                srv_sel_pen_wts = NULL

                                ) {

  messages_list <<- character(0) # string to attach to for printing messages
  if(input_list$store_config) input_list$config$Setup_Mod_Weighting <- mget(names(formals()))[-1]

  input_list$data$Wt_Catch <- Wt_Catch
  input_list$data$Wt_FishIdx <- Wt_FishIdx
  input_list$data$Wt_SrvIdx <- Wt_SrvIdx

  # A multivariate normal index cannot be split across observations, so its whole
  # fleet nLL is returned in a single cell (model_distributions.R). The jnLL then
  # applies the weights elementwise, which means a year-varying weight on an MVN
  # fleet silently applies only the first observed cell's value to the entire
  # fleet, and a zero there deletes the fleet outright. Constant weights are the
  # only ones that behave, so anything else is an error rather than a surprise.
  check_mvn_weight <- function(wt, like_type, use, n_flt, nm) {
    if(is.null(like_type) || length(wt) == 1) return(invisible(NULL))
    for(fl in seq_len(n_flt)) {
      if(is.na(like_type[fl]) || like_type[fl] != 2) next
      w <- wt[,,,fl][use[,,,fl] == 1]
      if(length(w) > 1 && length(unique(w)) > 1)
        stop(nm, " varies across observations for fleet ", fl, ", which uses a multivariate normal likelihood. ",
             "An MVN fleet's likelihood is a single number over the whole observation vector, so only a constant weight is meaningful.")
    }
    invisible(NULL)
  }
  check_mvn_weight(Wt_SrvIdx, input_list$data$SrvIdx_LikeType, input_list$data$UseSrvIdx,
                   input_list$data$n_srv_fleets, "Wt_SrvIdx")
  check_mvn_weight(Wt_FishIdx, input_list$data$FishIdx_LikeType, input_list$data$UseFishIdx,
                   input_list$data$n_fish_fleets, "Wt_FishIdx")
  input_list$data$Wt_Catch_pop <- Wt_Catch_pop
  input_list$data$Wt_FishIdx_pop <- Wt_FishIdx_pop
  input_list$data$Wt_SrvIdx_pop <- Wt_SrvIdx_pop
  # The recruitment and initial-age penalties are dimensioned differently, so an
  # array Wt_Rec cannot also serve the initial ages. A scalar still covers both.
  rec_dev_dim <- dim(input_list$par$ln_RecDevs)
  init_dev_dim <- dim(input_list$par$ln_InitDevs)
  if(length(Wt_Rec) != 1 && !identical(as.integer(dim(Wt_Rec)), as.integer(rec_dev_dim))) {
    stop("Wt_Rec must be a scalar or an array of dimension ", paste(rec_dev_dim, collapse = " x "),
         " (n_pop x n_regions x the third dimension of ln_RecDevs, which dont_est_recdev_last and n_proj_yrs_devs both change).")
  }
  if(is.null(Wt_Init_Rec)) {
    if(length(Wt_Rec) != 1) stop("Wt_Rec is an array, so Wt_Init_Rec must be supplied explicitly; the recruitment and initial age penalties are dimensioned differently.")
    Wt_Init_Rec <- Wt_Rec
  }
  if(length(Wt_Init_Rec) != 1 && !identical(as.integer(dim(Wt_Init_Rec)), as.integer(init_dev_dim))) {
    stop("Wt_Init_Rec must be a scalar or an array of dimension ", paste(init_dev_dim, collapse = " x "), ".")
  }
  if(length(Wt_Rec) > 1 && any(Wt_Rec == 0)) collect_message("Recruitment deviations excluded from the recruitment penalty but still estimated: ", sum(Wt_Rec[1,1,] == 0), " of ", rec_dev_dim[3])

  input_list$data$Wt_Rec <- Wt_Rec
  input_list$data$Wt_Init_Rec <- Wt_Init_Rec
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

  # Retained catch stuff
  input_list$data$Wt_Discard <- Wt_Discard
  input_list$data$Wt_Discard_pop <- Wt_Discard_pop
  input_list$data$Wt_D <- Wt_D
  input_list$data$Wt_FishAgeComps_discard <- Wt_FishAgeComps_discard
  input_list$data$Wt_FishLenComps_discard <- Wt_FishLenComps_discard
  input_list$data$Wt_FishAgeComps_discard_pop <- Wt_FishAgeComps_discard_pop
  input_list$data$Wt_FishLenComps_discard_pop <- Wt_FishLenComps_discard_pop
  input_list$data$fish_sel_pen_wts <- resolve_sel_pen_wts(fish_sel_pen_wts, input_list$data$n_fish_fleets)
  input_list$data$ret_sel_pen_wts <- resolve_sel_pen_wts(ret_sel_pen_wts, input_list$data$n_fish_fleets)
  input_list$data$srv_sel_pen_wts <- resolve_sel_pen_wts(srv_sel_pen_wts, input_list$data$n_srv_fleets)

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
