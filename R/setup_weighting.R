# Stage 1 of 3: model setup
#
# Relative weighting across data sources, applied as multipliers on each likelihood component. Runs
# last in the setup chain, once every data source it weights has been declared.

#' Set likelihood and penalty weights for the estimation model
#'
#' Assigns the \eqn{\lambda} multiplier on each likelihood component and penalty
#' term, which is how a noisy data source is down-weighted, how an iterative
#' reweighting such as Francis is applied, and how a component is switched off
#' entirely, with a weight of \code{0}. Call after every data setup function.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}.
#' @param addtocomp Small constant added to the composition proportions to avoid
#'   \code{log(0)}. Default \code{1e-3}. Ignored by the logistic normal, which
#'   handles zeros itself.
#' @param comp_const_obs Integer switch for where \code{addtocomp} enters the
#'   multinomial, not a constant to tune. \code{1} (default) adds it to the
#'   observed proportions that weight the likelihood as well as inside the
#'   logarithms, so the likelihood is stationary at \code{pred = obs}; \code{0}
#'   weights by the raw observed proportions. With a Dirichlet-multinomial
#'   conditional age-at-length fleet, \code{1} warns, since the constant biases
#'   theta upward when most age bins in a length bin are structurally empty.
#' @param addtofishidx,addtosrvidx Small constants added to the fishery and survey
#'   indices. Default \code{1e-4}.
#' @param addtotag Small constant added to the tag recovery observations. Default
#'   \code{1e-10}.
#' @param Wt_Catch,Wt_FishIdx Weights on the catch and fishery index likelihoods,
#'   a scalar or an array \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   Default \code{1}.
#' @param Wt_SrvIdx Weight on the survey index likelihood, a scalar or an array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. Default \code{1}.
#' @param Wt_Catch_pop,Wt_FishIdx_pop The population-specific catch and fishery
#'   index weights, a scalar or an array \code{[n_pop × n_regions × n_years ×
#'   n_seas × n_fish_fleets]}. Default \code{1}.
#' @param Wt_SrvIdx_pop The population-specific survey index weight, a scalar or an
#'   array \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}. Default
#'   \code{1}.
#' @param Wt_FishAgeComps,Wt_FishLenComps,Wt_FishAgeComps_discard,Wt_FishLenComps_discard
#'   Weights on the fishery and discard composition likelihoods, a scalar or an
#'   array \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}. Default
#'   one everywhere.
#' @param Wt_SrvAgeComps,Wt_SrvLenComps Weights on the survey composition
#'   likelihoods, a scalar or an array \code{[n_regions × n_years × n_seas ×
#'   n_sexes × n_srv_fleets]}. Default one everywhere.
#' @param Wt_FishAgeComps_pop,Wt_FishLenComps_pop,Wt_FishAgeComps_discard_pop,Wt_FishLenComps_discard_pop
#'   The population-specific fishery and discard composition weights, a scalar or
#'   an array \code{[n_pop × n_regions × n_years × n_seas × n_sexes ×
#'   n_fish_fleets]}. Default one everywhere.
#' @param Wt_SrvAgeComps_pop,Wt_SrvLenComps_pop The population-specific survey
#'   composition weights, a scalar or an array \code{[n_pop × n_regions × n_years ×
#'   n_seas × n_sexes × n_srv_fleets]}. Default one everywhere.
#' @param Wt_Rec Weight on the recruitment deviation penalty, a scalar or an array
#'   \code{[n_pop × n_regions × n_est_rec_devs]}, where the third dim is
#'   \code{ln_RecDevs}'s own rather than the number of years, since
#'   \code{dont_est_recdev_last} and \code{n_proj_yrs_devs} both move it. Default
#'   \code{1}. A zero leaves a deviation estimated but takes it out of the penalty,
#'   which is how a stock-recruit relationship is fit over a window of years while
#'   recruitment stays free in every year; \code{dont_est_recdev_last} instead
#'   removes the deviations, so recruitment reverts to the deterministic
#'   prediction.
#' @param Wt_Init_Rec Weight on the initial age deviation penalty, a scalar or an
#'   array \code{[n_pop × n_regions × (n_ages - 1) × n_sexes]}. \code{NULL}
#'   (default) takes \code{Wt_Rec} when that is a scalar; supply it explicitly when
#'   \code{Wt_Rec} is an array, since the two penalties are dimensioned
#'   differently.
#' @param Wt_F Scalar weight on the fishing mortality deviation penalty. Default
#'   \code{1}.
#' @param Wt_Tagging Scalar weight on the tag recovery likelihood. Default \code{1}.
#' @param Wt_Discard Weight on the aggregated discard amount or fraction
#'   likelihood, a scalar or an array \code{[n_regions × n_years × n_seas ×
#'   n_fish_fleets]}. Default \code{1}.
#' @param Wt_Discard_pop The population-specific discard weight, a scalar or an
#'   array with a leading \code{n_pop} dim. Default \code{1}.
#' @param Wt_D Scalar weight on the discard mortality rate deviation penalty.
#'   Default \code{1}.
#' @param Wt_Fish_caal Weight on the fishery conditional age-at-length likelihood,
#'   multiplying each length bin's input sample size. Array \code{[n_regions x
#'   n_years x n_seas x n_lens x n_sexes x n_fish_fleets]}, the shape of
#'   \code{ISS_Fish_caal}. Default one everywhere.
#' @param Wt_Srv_caal The survey counterpart, with \code{n_srv_fleets} last.
#'   Default one everywhere.
#' @param fish_sel_pen_wts \code{NULL} (default), or a named numeric vector or list
#'   weighting any subset of six selectivity smoothness penalties, which are
#'   evaluated on the fleet's realized selectivity by bin and year surface and so
#'   apply to any functional form: \code{"smooth_bin_curve"} and
#'   \code{"smooth_bin_diff"} are the second and first difference across bins,
#'   \code{"smooth_yr_diff"} and \code{"smooth_yr_curve"} the same across years,
#'   \code{"smooth_dome"} penalizes non-monotonicity across bins, and
#'   \code{"smooth_mean_center"} regularizes each year's mean. See
#'   \code{\link{resolve_sel_pen_wts}} and
#'   \code{\link{Get_Selex_Smoothness_Penalty}}. A name left out is \code{0}. Each
#'   weight may instead be a vector with one value per model year, so a penalty can
#'   act in some years only or at a different strength in each. The specification
#'   may also hold \code{"bin_range"}, the first and last bin the penalties act
#'   over. Pass an unnamed list of per-fleet specifications to give each fleet its
#'   own. Call after \code{Setup_Mod_Fishsel_and_Q}.
#' @param ret_sel_pen_wts As \code{fish_sel_pen_wts}, for retained fishery
#'   selectivity.
#' @param srv_sel_pen_wts As \code{fish_sel_pen_wts}, for survey selectivity. Call
#'   after \code{Setup_Mod_Srvsel_and_Q}.
#'
#' @return \code{input_list} with every weight stored in \code{$data} under its own
#'   name.
#'
#' @export Setup_Mod_Weighting
#' @family Model Setup
Setup_Mod_Weighting <- function(input_list,
                                addtocomp = 1e-3,
                                comp_const_obs = 1,
                                addtofishidx = 1e-4,
                                addtosrvidx = 1e-4,
                                addtotag = 1e-10,
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

                                # Conditional age-at-length
                                Wt_Fish_caal = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                                                length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                Wt_Srv_caal = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                                               length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_srv_fleets)),

                                # Selectivity penalty weights
                                fish_sel_pen_wts = NULL,
                                ret_sel_pen_wts = NULL,
                                srv_sel_pen_wts = NULL

                                ) {

  messages_list <<- character(0) # string to attach to for printing messages # nolint: object_usage_linter.
  if(input_list$store_config) input_list$config$Setup_Mod_Weighting <- mget(names(formals()))[-1]

  # A value still passed to the deprecated Setup_Mod_Biologicals arguments wins
  # over this function's own default, but not over a value supplied here
  legacy <- input_list$.legacy_weighting
  resolve_legacy <- function(val, is_missing, legacy_val, arg_name) {
    if(is_missing && !is.null(legacy_val)) {
      collect_message(arg_name, " taken from the deprecated Setup_Mod_Biologicals argument of the same name.")
      return(legacy_val)
    }
    val
  }

  addtocomp      <- resolve_legacy(addtocomp, missing(addtocomp), legacy$addtocomp, "addtocomp")
  comp_const_obs <- resolve_legacy(comp_const_obs, missing(comp_const_obs), legacy$comp_const_obs, "comp_const_obs")
  addtofishidx   <- resolve_legacy(addtofishidx, missing(addtofishidx), legacy$addtofishidx, "addtofishidx")
  addtosrvidx    <- resolve_legacy(addtosrvidx, missing(addtosrvidx),  legacy$addtosrvidx,  "addtosrvidx")
  addtotag       <- resolve_legacy(addtotag, missing(addtotag), legacy$addtotag, "addtotag")
  input_list$.legacy_weighting <- NULL
  if(!comp_const_obs %in% c(0, 1)) stop("comp_const_obs must be 0 or 1. It is a switch, not a constant: 1 weights the multinomial by obs + addtocomp (unbiased), 0 weights by the raw observed proportions (the ADMB convention).")
  input_list$data$addtocomp <- addtocomp
  input_list$data$comp_const_obs <- comp_const_obs
  input_list$data$addtofishidx <- addtofishidx
  input_list$data$addtosrvidx <- addtosrvidx
  input_list$data$addtotag <- addtotag

  # comp_const_obs = 1 inflates Dirichlet-multinomial theta on CAAL's structurally empty
  # bins. Fish_caal_LikeType/Srv_caal_LikeType code 1 for that likelihood.
  dm_caal_fleet_types <- c(if(isTRUE(any(input_list$data$Fish_caal_LikeType == 1))) "fishery",
                            if(isTRUE(any(input_list$data$Srv_caal_LikeType == 1))) "survey")
  if(length(dm_caal_fleet_types) > 0 && comp_const_obs == 1)
    warning(paste0("Conditional age-at-length uses the Dirichlet-Multinomial for at least one ", paste(dm_caal_fleet_types, collapse = " and "),
                   " fleet with comp_const_obs = 1. The constant added to the observed proportions biases theta upward ",
                   "when most age bins in a length bin are structurally empty; set comp_const_obs = 0 for conditional age-at-length."), call. = FALSE)

  input_list$data$Wt_Catch <- Wt_Catch
  input_list$data$Wt_FishIdx <- Wt_FishIdx
  input_list$data$Wt_SrvIdx <- Wt_SrvIdx

  # Input Validation --------------------------------------------------------
  # Checking to see if sigma is identifiable ...
  check_sigma_weight_confound <- function(wt, form, arg_name, spec_name) {
    if(!is.null(form) && form > 0 && any(wt != 1)) {
      warning(arg_name, " is not 1 everywhere while ", spec_name, " estimates the index ",
              "observation error. A likelihood weight and an estimated standard ",
              "deviation are confounded, so the estimate will absorb the weight. ",
              "Set ", arg_name, " to 1 when estimating, or fix the sigma when weighting.")
    }
  }
  check_sigma_weight_confound(Wt_SrvIdx, input_list$data$sigmaSrvIdx_form, "Wt_SrvIdx", "sigmaSrvIdx_spec")
  check_sigma_weight_confound(Wt_FishIdx, input_list$data$sigmaFishIdx_form, "Wt_FishIdx", "sigmaFishIdx_spec")

  # Checking validity of MVN weights
  check_mvn_weight <- function(wt, like_type, use, n_flt, arg_name) {
    if(is.null(like_type) || length(wt) == 1) return(invisible(NULL))
    for(fl in seq_len(n_flt)) {
      if(is.na(like_type[fl]) || like_type[fl] != 2) next
      w <- wt[,,,fl][use[,,,fl] == 1]
      if(length(w) > 1 && length(unique(w)) > 1)
        stop(arg_name, " varies across observations for fleet ", fl, ", which uses a multivariate normal likelihood. ",
             "An MVN fleet's likelihood is a single number over the whole observation vector, so only a constant weight is meaningful.")
    }
    invisible(NULL)
  }

  check_mvn_weight(Wt_SrvIdx, input_list$data$SrvIdx_LikeType, input_list$data$UseSrvIdx, input_list$data$n_srv_fleets, "Wt_SrvIdx")
  check_mvn_weight(Wt_FishIdx, input_list$data$FishIdx_LikeType, input_list$data$UseFishIdx, input_list$data$n_fish_fleets, "Wt_FishIdx")

  input_list$data$Wt_Catch_pop <- Wt_Catch_pop
  input_list$data$Wt_FishIdx_pop <- Wt_FishIdx_pop
  input_list$data$Wt_SrvIdx_pop <- Wt_SrvIdx_pop

  # The recruitment and initial-age penalties are dimensioned differently, so an array Wt_Rec cannot also serve the initial ages. A scalar still covers both.
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

  # Populate Data List ------------------------------------------------------
  # input into list
  input_list$data$Wt_Rec <- Wt_Rec
  input_list$data$Wt_Init_Rec <- Wt_Init_Rec
  input_list$data$Wt_F <- Wt_F
  input_list$data$Wt_FishAgeComps <- Wt_FishAgeComps
  input_list$data$Wt_SrvAgeComps <- Wt_SrvAgeComps
  input_list$data$Wt_FishLenComps <- Wt_FishLenComps
  input_list$data$Wt_SrvLenComps <- Wt_SrvLenComps
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

  # Checking for conditional age-at-length stuff
  caal_dim <- function(n_fleets) c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                   length(input_list$data$lens), input_list$data$n_sexes, n_fleets)
  if(!identical(as.integer(dim(Wt_Fish_caal)), as.integer(caal_dim(input_list$data$n_fish_fleets))))
    stop("Wt_Fish_caal must be an array of dimension n_regions x n_years x n_seas x n_lens x n_sexes x n_fish_fleets (", paste(caal_dim(input_list$data$n_fish_fleets), collapse = " x "), ").")
  if(!identical(as.integer(dim(Wt_Srv_caal)), as.integer(caal_dim(input_list$data$n_srv_fleets))))
    stop("Wt_Srv_caal must be an array of dimension n_regions x n_years x n_seas x n_lens x n_sexes x n_srv_fleets (", paste(caal_dim(input_list$data$n_srv_fleets), collapse = " x "), ").")
  input_list$data$Wt_Fish_caal <- Wt_Fish_caal
  input_list$data$Wt_Srv_caal <- Wt_Srv_caal

  input_list$data$fish_sel_pen_wts <- resolve_sel_pen_wts(fish_sel_pen_wts, input_list$data$n_fish_fleets)
  input_list$data$ret_sel_pen_wts <- resolve_sel_pen_wts(ret_sel_pen_wts, input_list$data$n_fish_fleets)
  input_list$data$srv_sel_pen_wts <- resolve_sel_pen_wts(srv_sel_pen_wts, input_list$data$n_srv_fleets)

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
