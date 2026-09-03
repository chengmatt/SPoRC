# Stage 3 of 3: post fit
#
# One step ahead residuals for compositions, indices and tags. The internal
# routines read the residual bookkeeping that the objective function stored during
# its evaluation; run_external_comp_osa recomputes residuals from scratch for
# compositions the objective did not carry.

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
    Bins     = paste0(comp_source, "Comps", suffix, "_bins"),
    n_fleets_field = if(grepl("^Srv", comp_source)) "n_srv_fleets" else "n_fish_fleets"
  )
}

#' Default bin labels for an internal OSA call
#'
#' \code{bins} and \code{bin_label} label the residual's bin axis, and
#' \code{\link{plot_resids}} draws its second panel against them. Left empty
#' the label columns never reach the residual frame and that panel errors when it
#' is printed, so they are filled from the data here: the observed age bins for
#' an age or conditional age-at-length source and the length bins for a length
#' source. The observed bins are read off the observation array, since ageing
#' error can leave fewer of them than the model carries.
#'
#' @param data Data list from the fitted model.
#' @param comp_source Composition source name, as \code{\link{get_osa}} takes it.
#' @param pop,discard Whether the source is population-specific or the discard
#'   stream, as \code{\link{get_osa}} takes them.
#'
#' @returns List with \code{bins} and \code{bin_label}.
#'
#' @keywords internal
osa_default_bins <- function(data, comp_source, pop = FALSE, discard = FALSE) {

  is_caal <- comp_source %in% c("Fish_caal", "Srv_caal")
  is_len <- grepl("Len$", comp_source)

  if(is_caal) {
    obs_nm <- paste0("Obs", comp_source)
    n_bins <- dim(data[[obs_nm]])[5] # region, year, season, len, AGE, sex, fleet
  } else {
    obs_nm <- comp_osa_field_map(comp_source, pop = pop, discard = discard)$Obs
    n_bins <- dim(data[[obs_nm]])[if(pop) 5 else 4]
  }

  if(is_len) {
    bins <- if(length(data$lens) == n_bins) data$lens else seq_len(n_bins)
    return(list(bins = bins, bin_label = "Length"))
  }

  # ages: the model's own when the observation carries all of them, and the bin
  # index otherwise, since which model ages an observed bin holds is the ageing
  # error's business rather than something to guess at here
  bins <- if(length(data$ages) == n_bins) data$ages else seq_len(n_bins)
  return(list(bins = bins, bin_label = "Age"))
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

#' Call \code{RTMB::oneStepPredict()} with the model's TMB DLL resolved
#'
#' Wrapper used by the internal OSA runners, which works around two quirks of
#' \code{\link[RTMB]{oneStepPredict}}:
#' \itemize{
#'   \item Its \code{parallel} branch calls \code{TMB::openmp()} without a
#'     \code{DLL} argument, so TMB falls back to guessing the DLL and errors
#'     with "Multiple TMB models loaded" whenever a session has more than one
#'     TMB DLL loaded (e.g. RTMB alongside compResidual, which
#'     \code{\link{run_external_comp_osa}} loads). 
#'   \item \code{discreteSupport} is detected with \code{missing()}, so
#'     supplying it as \code{NULL} is not the same as omitting it: a \code{NULL}
#'     sends continuous families down the mixed discrete/continuous branch,
#'     which rejects the Gaussian methods. It is forwarded here only when it is
#'     non-\code{NULL}.
#' }
#'
#' @param model A fitted RTMB model object from \code{\link{fit_model}}.
#' @param ... Further arguments passed to \code{\link[RTMB]{oneStepPredict}}.
#' @param discreteSupport Support of the discrete observations, or \code{NULL}
#'   (the default) to omit the argument entirely.
#' @param parallel Whether or not to parallelize OSA computation. Defaults to
#'   \code{FALSE}.
#'
#' @return The \code{\link[RTMB]{oneStepPredict}} result.
#' @keywords internal
osa_one_step_predict <- function(model, ..., discreteSupport = NULL, parallel = FALSE) {

  run_osa <- function(use_parallel) {
    if(is.null(discreteSupport)) {
      RTMB::oneStepPredict(model, ..., parallel = use_parallel)
    } else {
      RTMB::oneStepPredict(model, ..., discreteSupport = discreteSupport, parallel = use_parallel)
    }
  }

  # only the parallel branch consults TMB's DLL guess, so serial calls need no patching
  if(!isTRUE(parallel)) return(run_osa(parallel))

  # without a resolvable DLL there is nothing to bind, so take oneStepPredict's usual path
  dll <- model$env$DLL
  if(!is.character(dll) || length(dll) != 1 || !nzchar(dll)) return(run_osa(parallel))

  openmp_orig <- get("openmp", envir = asNamespace("TMB"))
  openmp_dll <- openmp_orig
  formals(openmp_dll)$DLL <- dll

  bound <- tryCatch({
    utils::assignInNamespace("openmp", openmp_dll, ns = "TMB")
    TRUE
  }, error = function(e) FALSE)

  if(!bound) {
    warning("Unable to bind TMB DLL '", dll, "' for parallel OSA computation; running serially instead.")
    return(run_osa(FALSE))
  }

  on.exit(utils::assignInNamespace("openmp", openmp_orig, ns = "TMB"), add = TRUE)
  run_osa(TRUE)
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
#' @param bins Vector of age or length bin labels for display. Must span every observed bin of the stream, not just
#'   the ones a \code{*_bins} restriction fits: residuals are labeled by true
#'   observed bin number, so a subset here shifts every label.
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
    pop = pop, n_pop = n_pop, return_labels = TRUE,
    # must match what the objective packed, or the tracked vector is a different length
    BinsArr = bins_or_null(data[[fm$Bins]])
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

  osa <- osa_one_step_predict(model, observation.name = tracked_name, method = method,
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

#' Run internal (model-based) OSA residuals for conditional age-at-length data
#'
#' Internal counterpart to \code{\link{run_internal_comp_osa}} for CAAL data
#' packed via \code{\link{pack_caal_osa}} (requires
#' \code{do_internal_comp_osa = TRUE} in \code{\link{Setup_Mod_Dim}}). Called by
#' \code{\link{get_osa}} when \code{comp_source} is \code{"Fish_caal"} or
#' \code{"Srv_caal"}. CAAL carries only the discrete families, so there is no
#' family argument; the residuals come back with an extra \code{len} column
#' giving the length bin each age composition was conditioned on.
#'
#' @param model Fitted model object.
#' @param data The data list the model was built from.
#' @param comp_source Either \code{"Fish_caal"} or \code{"Srv_caal"}.
#' @param bins Age bins, used to label the residuals. Must span every observed bin of the stream, not just
#'   the ones a \code{*_bins} restriction fits: residuals are labeled by true
#'   observed bin number, so a subset here shifts every label.
#' @param bin_label Label for the bin axis, typically \code{"Age"}.
#' @param osa_method Optional override for the \code{oneStepPredict} method.
#' @param parallel Logical, passed to \code{oneStepPredict}.
#'
#' @return A list with one element \code{res}, matching
#'   \code{\link{run_internal_comp_osa}}'s schema plus a \code{len} column.
#'
#' @keywords internal
run_internal_caal_osa <- function(model, data, comp_source, bins, bin_label,
                                  osa_method = NULL, parallel = FALSE) {

  if(!comp_source %in% c("Fish_caal", "Srv_caal")) stop("`comp_source` for CAAL must be one of: Fish_caal, Srv_caal")

  obs_nm <- paste0("Obs", comp_source)
  n_fleets_field <- if(comp_source == "Srv_caal") "n_srv_fleets" else "n_fish_fleets"

  packed <- pack_caal_osa(
    ObsArr = data[[obs_nm]], ISSArr = data[[paste0("ISS_", comp_source)]],
    WtArr = data[[paste0("Wt_", comp_source)]], UseArr = data[[paste0("Use", comp_source)]],
    TypeMat = data[[paste0(comp_source, "_Type")]], LikeTypeVec = data[[paste0(comp_source, "_LikeType")]],
    n_yrs = length(data$years), n_seas = data$n_seas, n_lens = length(data$lens),
    n_fleets = data[[n_fleets_field]], n_sexes = data$n_sexes, addtocomp = data$addtocomp,
    return_labels = TRUE,
    # must match what the objective packed, or the tracked vector is a different length
    BinsArr = bins_or_null(data[[paste0(comp_source, "_bins")]])
  )

  if(is.null(packed)) {
    warning("No conditional age-at-length data found for comp_source '", comp_source, "'; returning NULL.")
    return(NULL)
  }

  tracked_name <- paste0(obs_nm, "_osa")
  method <- if(!is.null(osa_method)) osa_method else "oneStepGeneric"
  validate_osa_method(method)
  subset_idx <- osa_keep_subset(packed$labels$last_in_group)

  osa <- osa_one_step_predict(model, observation.name = tracked_name, method = method,
                              discrete = TRUE, parallel = parallel,
                              subset = subset_idx, trace = FALSE,
                              discreteSupport = 0:max(model$env$obs[[tracked_name]]))

  lab <- packed$labels[subset_idx, ]
  lab$resid <- osa$residual

  res <- lab %>%
    dplyr::transmute(
      fleet = as.character(fleet),
      index_label = bin_label,
      year = year,
      index = bins[bin],
      len = data$lens[len],
      resid = resid,
      region = region,
      sex = sex,
      seas = season,
      pop = 1,
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

  osa <- osa_one_step_predict(model, observation.name = tracked_name, method = method,
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
#' @param index_source One of \code{"Catch"}, \code{"Discard"}, \code{"FishIdx"},
#'   \code{"SrvIdx"}, or their at-age forms \code{"CatchAA"}, \code{"DiscardAA"},
#'   \code{"SrvIdxAA"}. At-age sources return extra
#'   \code{age} and \code{sex} columns.
#' @param pop Logical; population-specific index source.
#'
#' @return A named list of field names: \code{Obs}, \code{Use}.
#' @keywords internal
index_osa_field_map <- function(index_source, pop = FALSE) {
  valid_sources <- c("Catch", "Discard", "FishIdx", "SrvIdx",
                     "CatchAA", "DiscardAA", "SrvIdxAA")
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
#' @param index_source One of \code{"Catch"}, \code{"Discard"}, \code{"FishIdx"},
#'   \code{"SrvIdx"}, or their at-age forms \code{"CatchAA"}, \code{"DiscardAA"},
#'   \code{"SrvIdxAA"}. At-age sources return extra
#'   \code{age} and \code{sex} columns.
#' @param pop Logical; population-specific index source. Default \code{FALSE}.
#' @param osa_method Optional override for \code{RTMB::oneStepPredict}'s
#'   \code{method}. Must be one of \code{"oneStepGeneric"},
#'   \code{"oneStepGaussianOffMode"}, or \code{"oneStepGaussian"} (the
#'   \code{"cdf"} method is not permitted). Defaults to \code{"oneStepGeneric"}.
#' @param parallel Whether or not to parallelize OSA computation. Defaults to \code{FALSE}.
#'
#' @return A list with one element \code{res}: columns \code{fleet},
#'   \code{region}, \code{year}, \code{season}, \code{pop}, \code{age},
#'   \code{sex}, \code{resid}, and \code{idx_type} (set to
#'   \code{index_source}; named \code{idx_type} rather than \code{comp_type}
#'   because index-type sources are not compositions), or \code{NULL} if no data
#'   of the requested source is present. \code{age} and \code{sex} are
#'   \code{NA} for the aggregated sources.
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
  # at-age sources carry an age and a sex dimension before the fleet
  at_age <- grepl("AA$", index_source)
  dim_names <- c(if(pop) "pop", "region", "year", "season", if(at_age) c("age", "sex"), "fleet")
  map <- as.data.frame(arrayInd(valid_idx, dim(use_arr)))
  colnames(map) <- dim_names

  # a margin the fleet sums over holds its observation in slot one, and that slot
  # is not the first region or the first sex. Labelling it as such would face
  # aggregated residuals alongside split ones once the setting varies by year
  aa_split <- list(region = TRUE, sex = TRUE)
  if(at_age) {
    aa_type <- data[[paste0(index_source, if(pop) "_pop", "_Type")]]
    code <- if(is.null(aa_type)) rep(3, nrow(map)) else
            if(is.null(dim(aa_type))) aa_type[map$fleet] else
            aa_type[cbind(map$year, map$fleet)]
    aa_split <- at_age_split(code)
  }

  tracked_name <- fm$Obs
  method <- if(!is.null(osa_method)) osa_method else "oneStepGeneric"
  validate_osa_method(method)

  osa <- osa_one_step_predict(model, observation.name = tracked_name, method = method,
                              discrete = FALSE, parallel = parallel, trace = FALSE)

  res <- data.frame(
    fleet = as.character(map$fleet),
    region = if(at_age) ifelse(aa_split$region, as.character(map$region), "summed") else map$region,
    year = data$years[map$year],
    season = map$season,
    pop = if(pop) map$pop else 1L,
    age = if(at_age) data$ages[map$age] else NA_integer_,
    sex = if(at_age) ifelse(aa_split$sex, as.character(map$sex), "summed") else NA_integer_,
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
#' @param N Input (or effective if Multinomial) sample size, at the model's full
#'   year dimension in every case: \code{years} (or \code{years_by_region})
#'   selects the years actually used, the same way it selects them from
#'   \code{obs_mat}/\code{exp_mat}, so \code{N} is never pre-filtered by the
#'   caller. Dimensions depend on \code{comp_type}:
#'   \itemize{
#'     \item \code{comp_type = 0} (aggregated): vector of length \code{n_years}
#'       (the model's, not the length of \code{years}).
#'     \item \code{comp_type = 1} (split by region and sex): array
#'       \code{[n_regions, n_years, n_sexes]}.
#'     \item \code{comp_type = 2} (split by region, joint by sex): matrix
#'       \code{[n_regions, n_years]}.
#'   }
#'   For years without data, users can simply input an NA or any abritary number (it gets filtered out within the function). This is the same array \code{ISS_*Comps} already is in the model's data list, so it can usually be passed straight through.
#' @param years Years with composition data. Either a plain vector, used for every region, or a list with one vector of years per region when the regions carry different years. Both forms work for every composition type. A region with no years is skipped.
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
#'     \item 2-4 = logistic-normal variants
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
#'   internal OSA residuals for. Conditional age-at-length uses
#'   \code{"Fish_caal"} or \code{"Srv_caal"}, which return an extra
#'   \code{len} column giving the length bin each age composition was
#'   conditioned on and ignore \code{family}, since CAAL carries only the
#'   discrete likelihoods. Required when \code{model} is supplied and
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

    # plot_resids draws its second panel against the bin labels, so fill them
    # from the data when the caller leaves them out
    if(!is.null(comp_source) && (is.null(bins) || is.null(bin_label))) {
      defaults <- osa_default_bins(data, comp_source, pop = pop, discard = discard)
      if(is.null(bins)) bins <- defaults$bins
      if(is.null(bin_label)) bin_label <- defaults$bin_label
    }

    if(!is.null(index_source)) {
      return(run_internal_index_osa(model = model, data = data, index_source = index_source,
                                    pop = pop, osa_method = osa_method, parallel = parallel))
    } else if(tag) {
      return(run_internal_tag_osa(model = model, data = data, osa_method = osa_method, parallel = parallel))
    } else if(!is.null(comp_source) && comp_source %in% c("Fish_caal", "Srv_caal")) {
      return(run_internal_caal_osa(model = model, data = data, comp_source = comp_source,
                                   bins = bins, bin_label = bin_label,
                                   osa_method = osa_method, parallel = parallel))
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

    # The split composition types read one year vector per region, the aggregated
    # type a single vector for the one composition it has. Accept either form for
    # any type, so a caller does not have to know which branch it will land in.
    years_by_region <- if(is.list(years)) years else replicate(n_regions, years, simplify = FALSE)

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
      tmp_osa <- run_external_comp_osa(obs = tmp_obs, exp = tmp_exp, N = N[years], DM_theta = DM_theta,
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

        # a region a fleet never sampled has nothing to compute a residual from
        if(length(years_by_region[[r]]) == 0) next

        for(s in 1:n_sexes) {

          obs <- obs_mat[,years_by_region[[r]],seas,,,fleet, drop = FALSE] # get filtered observed matrix
          exp <- exp_mat[,years_by_region[[r]],seas,,,fleet, drop = FALSE] # get filtered expected matrix

          tmp_obs <- obs[r,,1,,s,1] # get observations
          tmp_exp <- exp[r,,1,,s,1] # get expected

          # A plain year vector says "these years" without saying which regions
          # sampled them, so a cell with nothing in it is skipped here rather
          # than normalized into a residual with no data behind it
          if(!any(is.finite(tmp_obs)) || sum(tmp_obs, na.rm = TRUE) == 0) next

          # normalize
          tmp_obs <- (tmp_obs + addtocomp) / rowSums(tmp_obs + addtocomp)
          tmp_exp <- (tmp_exp + addtocomp) / rowSums(tmp_exp + addtocomp)

          # compute OSA
          tmp_osa <- run_external_comp_osa(obs = tmp_obs, exp = tmp_exp, N = N[r,years_by_region[[r]],s], DM_theta = DM_theta[r,s],
                             years = years_by_region[[r]], comp_like = comp_like, LN_Sigma = LN_Sigma[r,,,s],
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

        # a region a fleet never sampled has nothing to compute a residual from
        if(length(years_by_region[[r]]) == 0) next

        obs <- obs_mat[,years_by_region[[r]],seas,,,fleet, drop = FALSE] # get filtered observed matrix
        exp <- exp_mat[,years_by_region[[r]],seas,,,fleet, drop = FALSE] # get filtered expected matrix

        # initialize to cbind
        tmp_obs <- NULL
        tmp_exp <- NULL

        for(s in 1:n_sexes) {
          tmp_obs <- cbind(tmp_obs, obs[r,,1,,s,1]) # get observations
          tmp_exp <- cbind(tmp_exp, exp[r,,1,,s,1]) # get expected
        } # end s loop

        # A plain year vector says "these years" without saying which regions
        # sampled them, so a region with nothing in it is skipped here rather
        # than normalized into a residual with no data behind it
        if(!any(is.finite(tmp_obs)) || sum(tmp_obs, na.rm = TRUE) == 0) next

        # normalize
        tmp_obs <- (tmp_obs + addtocomp) / rowSums(tmp_obs + addtocomp)
        tmp_exp <- (tmp_exp + addtocomp) / rowSums(tmp_exp + addtocomp)

        # compute OSA
        tmp_osa <- run_external_comp_osa(obs = tmp_obs, exp = tmp_exp, N = N[r,years_by_region[[r]]],
                           DM_theta = DM_theta[r], years = years_by_region[[r]], comp_like = comp_like, LN_Sigma = LN_Sigma[r,,],
                           index = paste(rep(1:n_sexes, each = length(bins)), "_", rep(bins, times = n_sexes), sep = ""),
                           fleet = as.character(fleet), index_label = bin_label)

        # Doing some naming stuff
        tmp_osa$res$region <- r
        tmp_osa$res$seas <- seas
        tmp_osa$res <- tmp_osa$res %>% dplyr::mutate(split_index = stringr::str_split(index, "_"),  # Split once and store as list
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
#' @return List of two plots, named and also safe to index by position:
#'   \code{sdnr_plot} (QQ-plot) first, then \code{bubble_plot}
#'   (composition/tag residual magnitude and sign) or \code{resid_plot}
#'   (index-type residual vs. year).
#' @export plot_resids
#' @family Model Diagnostics
#' @import dplyr
#' @import ggplot2
plot_resids <- function(osa_results) {

  # extract results
  res <- osa_results$res %>% dplyr::mutate(sign = ifelse(resid < 0, "Neg", "Pos"))

  # comp_type = composition/tag residuals, idx_type = index-type residuals
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
    age             = function(x) paste0("Age ", x),
    len             = function(x) paste0("Len ", x),
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
      dplyr::summarize(
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
  if(res_type %in% c("Catch", "Discard", "FishIdx", "SrvIdx",
                     "CatchAA", "DiscardAA", "SrvIdxAA")) {

    multi_region <- has_multi("region")
    multi_seas   <- has_multi("season")
    multi_age    <- has_multi("age")   # at-age sources only
    multi_sex    <- has_multi("sex")   # and only when the stream splits sexes

    idx_row_vars <- c(if(multi_region) "region", if(multi_age) "age")
    idx_col_vars <- c(if(multi_seas) "season", if(multi_fleet) "fleet", if(multi_pop) "pop",
                      if(multi_sex) "sex")

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

    return(list(sdnr_plot = sdnr_plot, resid_plot = resid_plot))
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
  # Conditional age-at-length residuals carry the length bin they were aged from,
  # a dimension the marginal compositions do not have. Without it in the facets
  # every length bin lands on the same year and age, so a model with any number
  # of bins draws as one solid stack. A model with many bins gives many panels;
  # subset osa_results$res on len first when that is too much to read.
  has_len <- has_multi("len")

  bubble_plot <- ggplot(data = res, aes(x = year, y = as.numeric(index),
                                        color = sign, size = abs(resid), alpha = abs(resid))) +
    geom_point() +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "Year", y = unique(res$index_label), color = "Sign", size = "abs(Resid)", alpha = "abs(Resid)") +
    theme_bw(base_size = 20) +
    theme(legend.position = 'top') +
    build_facet(c("region", if(has_len) "len"), c("sex", comp_extra_cols))

  return(list(sdnr_plot = sdnr_plot, bubble_plot = bubble_plot))

}


