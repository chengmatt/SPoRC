# Shared helpers
#
# Helpers for building and inspecting the input_list during setup: message
# collection, type coercion, safe extraction, year extension, and switching a
# data source off.

#' Resolve a modular selectivity penalty weight vector
#'
#' Converts a user-supplied selectivity penalty weight specification into the
#' complete named weight vector consumed by
#' \code{\link{Get_Selex_Smoothness_Penalty}}.
#'
#' @param pen_wts \code{NULL}, or a named numeric vector/list giving independent
#'   weights for any subset of \code{"smooth_bin_curve"}, \code{"smooth_bin_diff"},
#'   \code{"smooth_yr_diff"}, \code{"smooth_yr_curve"}, \code{"smooth_dome"},
#'   \code{"smooth_mean_center"}. Any name not supplied defaults to \code{0}.
#'   Each weight is either a scalar applied to every year, or a vector with one
#'   value per model year, which lets a penalty act only in some years or act
#'   with a different strength in each. The list may also carry
#'   \code{"bin_range"}, a length-two vector giving the first and last age or
#'   length bin the penalties apply over; the default of every bin is what a
#'   missing \code{bin_range} reproduces.
#'
#' @param n_fleets Integer. Number of fleets the specification covers. A single
#'   named specification is applied to every fleet; an unnamed list of length
#'   \code{n_fleets} gives each fleet its own.
#'
#' @return List of length \code{n_fleets}. Each element is a named list of
#'   length 7: the six penalty terms, each a scalar or a per-year vector, plus
#'   \code{bin_range} (\code{NULL} for all bins).
#'
#' @keywords internal
resolve_sel_pen_wts <- function(pen_wts, n_fleets = 1) {

  term_names <- c("smooth_bin_curve", "smooth_bin_diff", "smooth_yr_diff", "smooth_yr_curve", "smooth_dome", "smooth_mean_center")
  allowed <- c(term_names, "bin_range", "normalize", "yr_diff_ref")

  resolve_one <- function(spec) {

    out <- stats::setNames(vector("list", length(allowed)), allowed)
    for(nm in term_names) out[[nm]] <- 0
    out$bin_range <- NULL
    out$normalize <- TRUE
    out$yr_diff_ref <- NULL

    # NULL or an empty list both mean "this fleet has no penalties", which is how
    # a per-fleet specification names the fleets it does not constrain.
    if(is.null(spec) || length(spec) == 0) return(out)

    if(is.null(names(spec)) || !all(names(spec) %in% allowed))
      stop("pen_wts must be a named numeric vector/list with names in: ", paste(allowed, collapse = ", "))

    spec <- as.list(spec)
    for(nm in names(spec)) {
      if(is.null(spec[[nm]])) next
      if(nm == "bin_range" && !is.list(spec[[nm]])) {
        if(length(spec[[nm]]) != 2) stop("bin_range must be a length-two vector giving the first and last bin the selectivity penalties apply over.")
        if(spec[[nm]][1] > spec[[nm]][2]) stop("bin_range must be increasing.")
      }
      out[[nm]] <- spec[[nm]]
    } # end nm loop

    return(out)
  }

  # A single named specification covers every fleet; an unnamed list gives each
  # fleet its own, which is what surveys with different smoothing needs require.
  per_fleet <- is.list(pen_wts) && is.null(names(pen_wts))

  if(per_fleet) {
    if(length(pen_wts) != n_fleets) stop("pen_wts is an unnamed list of length ", length(pen_wts), " but there are ", n_fleets, " fleets. Supply one specification per fleet, or a single named specification for all of them.")
    return(lapply(pen_wts, resolve_one))
  }

  return(rep(list(resolve_one(pen_wts)), n_fleets))
}

#' Append a message to the global messages list
#'
#' Concatenates its arguments into a single string and appends the result to
#' \code{messages_list} in the calling environment via \code{<<-}. Used
#' internally to accumulate validation and setup messages for deferred display.
#'
#' @param ... Character strings passed to \code{paste(..., sep = "")}.
#'
#' @return \code{NULL} invisibly. Side effect: \code{messages_list} is
#'   updated in the parent environment.
#'
#' @keywords internal
collect_message <- function(...) {
  messages_list <<- c(messages_list, paste(..., sep = ""))
}

#' Safely extract a named element from a list object
#'
#' Returns the named element if it exists and is non-\code{NULL}; returns
#' \code{0} otherwise. Used to guard against missing or \code{NULL} report
#' fields when accumulating likelihood components.
#'
#' @param obj Named list, typically \code{obj$rep} from a fitted RTMB model.
#' @param name Character string. Name of the element to extract.
#'
#' @return The value of \code{obj[[name]]} if present and non-\code{NULL};
#'   \code{0} otherwise.
#'
#' @keywords internal
safe_extract <- function(obj, name) {
  if (name %in% names(obj) && !is.null(obj[[name]])) {
    return(obj[[name]])
  } else {
    return(0)
  }
}

#' Extend an array along its year dimension
#'
#' Appends additional year slices to an array along a specified dimension,
#' using one of several fill strategies. Used in SPoRC to extend biological,
#' selectivity, and sample-size arrays from the conditioning period into
#' projection years before running closed-loop MSE simulations.
#'
#' @param arr Array of any dimensionality to extend.
#' @param n_years Integer. Total number of years in the extended array
#'   (i.e., the new size of dimension \code{yr_dim}). Must be greater than
#'   \code{dim(arr)[yr_dim]}.
#' @param yr_dim Integer. Index of the dimension in \code{arr} corresponding
#'   to years.
#' @param fill Character string or numeric. Fill strategy for the appended
#'   year slices:
#'   \describe{
#'     \item{\code{"zeros"}}{Fill with zeros.}
#'     \item{\code{"last"}}{Repeat the last year slice that contains at least
#'       one non-\code{NA}, non-\code{NaN} value. If no valid slice exists,
#'       fills with \code{NA}.}
#'     \item{\code{"mean"}}{Fill with the per-element mean across years,
#'       excluding zeros, \code{NA}, and \code{NaN} values. Elements with no
#'       valid values are set to zero.}
#'     \item{\code{"F_pattern"}}{Fill with zeros; signals to downstream
#'       functions (e.g., \code{\link{predict_sim_fish_iss_fmort}}) that
#'       sample sizes should be dynamically updated based on projected fishing
#'       mortality during the closed-loop simulation.}
#'     \item{Numeric scalar or array}{Fill all appended slices with the
#'       supplied constant value, recycled via \code{array()}.}
#'   }
#'
#' @return Array with the same dimensions as \code{arr} except that
#'   \code{dim(result)[yr_dim] == dim(arr)[yr_dim] + n_years}, formed by
#'   binding \code{arr} and the fill array along \code{yr_dim} via
#'   \code{abind::abind}.
#'
#' @keywords internal
extend_years <- function(arr, n_years, yr_dim, fill = "zeros") {
  new_dims <- dim(arr); new_dims[yr_dim] <- n_years
  if(fill %in% c("zeros", "F_pattern")) {
    fill_array <- array(0, dim = new_dims)
  } else if(fill == "last") {
    # Get last non-NaN year slice along yr_dim
    # First, find the last year index that contains at least some non-NaN values
    last_valid_idx <- NULL
    for(i in dim(arr)[yr_dim]:1) {
      indices <- rep(list(quote(expr=)), length(dim(arr)))
      indices[[yr_dim]] <- i
      year_slice <- do.call(`[`, c(list(arr), indices, drop = FALSE))
      # check if this slice has any non-NaN values
      if(any(!is.na(year_slice) & !is.nan(year_slice))) {
        last_valid_idx <- i
        break
      }
    }
    # use NA if no valid year found
    if(is.null(last_valid_idx)) {
      fill_array <- array(NA, dim = new_dims)
    } else {
      # Get the last valid year slice
      indices <- rep(list(quote(expr=)), length(dim(arr)))
      indices[[yr_dim]] <- last_valid_idx
      last_year_slice <- do.call(`[`, c(list(arr), indices, drop = FALSE))

      # repeat slice n_years times
      fill_array <- array(0, dim = new_dims)
      for(i in 1:n_years) {
        fill_indices <- rep(list(quote(expr=)), length(dim(arr)))
        fill_indices[[yr_dim]] <- i
        fill_array <- do.call(`[<-`, c(list(fill_array), fill_indices, list(last_year_slice)))
      }
    }
  } else if (fill == "mean") {
    # get mean along the year dimension, excluding zeros and NaNs
    margins <- setdiff(seq_along(dim(arr)), yr_dim)
    # get mean excluding zeros and NaNs
    mean_slice <- apply(arr, margins, function(x) {
      valid_values <- x[!is.na(x) & !is.nan(x) & x != 0]
      if(length(valid_values) == 0) return(0) else mean(valid_values)
    })
    # extend mean_slice along year dimension
    fill_array <- array(mean_slice, dim = new_dims)
  } else if (is.numeric(fill)) {
    fill_array <- array(fill, dim = new_dims)
  }
  return(abind::abind(arr, fill_array, along = yr_dim))
}

#' Set data indicators to unused for specified years
#'
#' Sets \code{Use*} binary indicator arrays to \code{0} for the specified
#' years across one or more data types, and optionally removes conventional
#' tag cohorts released in those years. Used in MSE closed-loop simulations
#' and retrospective analyses to exclude future or withheld data from the
#' estimation model without modifying the underlying observation arrays.
#' Only years present in \code{data$years} are affected; out-of-range values
#' in \code{unused_years} are silently ignored.
#'
#' @param data Named list of model data as constructed by the
#'   \code{Setup_Mod_*} family of functions.
#' @param unused_years Integer vector. Year indices (relative to
#'   \code{data$years}) to mark as unused. Values not present in
#'   \code{1:length(data$years)} are dropped.
#' @param what Character vector. Data types to modify. Any combination of:
#'   \describe{
#'     \item{\code{"Catch"}}{Sets \code{UseCatch[, unused_years, , ] <- 0}.}
#'     \item{\code{"Catch_pop"}}{Sets \code{UseCatch_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"Discard"}}{Sets \code{UseDiscard[, unused_years, , ] <- 0}.}
#'     \item{\code{"Discard_pop"}}{Sets \code{UseDiscard_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"FishIdx"}}{Sets \code{UseFishIdx[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishIdx_pop"}}{Sets \code{UseFishIdx_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"FishAgeComps"}}{Sets \code{UseFishAgeComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishAgeComps_pop"}}{Sets \code{UseFishAgeComps_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"FishLenComps"}}{Sets \code{UseFishLenComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishLenComps_pop"}}{Sets \code{UseFishLenComps_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"FishAgeComps_discard"}}{Sets \code{UseFishAgeComps_discard[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishAgeComps_discard_pop"}}{Sets \code{UseFishAgeComps_discard_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"FishLenComps_discard"}}{Sets \code{UseFishLenComps_discard[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishLenComps_discard_pop"}}{Sets \code{UseFishLenComps_discard_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"SrvIdx"}}{Sets \code{UseSrvIdx[, unused_years, , ] <- 0}.}
#'     \item{\code{"SrvIdx_pop"}}{Sets \code{UseSrvIdx_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"SrvAgeComps"}}{Sets \code{UseSrvAgeComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"SrvAgeComps_pop"}}{Sets \code{UseSrvAgeComps_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"SrvLenComps"}}{Sets \code{UseSrvLenComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"SrvLenComps_pop"}}{Sets \code{UseSrvLenComps_pop[, , unused_years, , ] <- 0}.}
#'     \item{\code{"conv_tagging"}}{Removes tag cohorts whose release year
#'       falls in \code{unused_years} from \code{conv_tagged_fish},
#'       \code{obs_conv_tag_fish_recap}, \code{conv_tag_release_indicator},
#'       and updates \code{n_conv_tag_cohorts}. Only applied when
#'       \code{any(data$use_conv_fish_tagging == 1)}.}
#'   }
#'
#' @return The modified \code{data} list with \code{Use*} indicators set to
#'   \code{0} for \code{unused_years} and, if applicable, tagging cohorts
#'   from those years removed.
#'
#' @export set_data_indicator_unused
#' @family Utility
set_data_indicator_unused <- function(data,
                                      unused_years,
                                      what = c("Catch", "Catch_pop",
                                               "Discard", "Discard_pop",
                                               "FishIdx", "FishIdx_pop",
                                               "FishAgeComps", "FishAgeComps_pop",
                                               "FishLenComps", "FishLenComps_pop",
                                               "FishAgeComps_discard", "FishAgeComps_discard_pop",
                                               "FishLenComps_discard", "FishLenComps_discard_pop",
                                               "SrvIdx", "SrvIdx_pop",
                                               "SrvAgeComps", "SrvAgeComps_pop",
                                               "SrvLenComps", "SrvLenComps_pop",
                                               "conv_tagging")) {

  # figure out year dimensions
  data_years <- 1:length(data$years)
  unused_years <- unused_years[which(unused_years %in% data_years)]

  if(length(unused_years) > 0) {
    # set to not use
    if("Catch" %in% what) data$UseCatch[,unused_years,,] <- 0
    if("Catch_pop" %in% what) data$UseCatch_pop[,,unused_years,,] <- 0
    if("Discard" %in% what) data$UseDiscard[,unused_years,,] <- 0
    if("Discard_pop" %in% what) data$UseDiscard_pop[,,unused_years,,] <- 0
    if("FishIdx" %in% what) data$UseFishIdx[,unused_years,,] <- 0
    if("FishIdx_pop" %in% what) data$UseFishIdx_pop[,,unused_years,,] <- 0
    if("FishAgeComps" %in% what) data$UseFishAgeComps[,unused_years,,] <- 0
    if("FishAgeComps_pop" %in% what) data$UseFishAgeComps_pop[,,unused_years,,] <- 0
    if("FishLenComps" %in% what) data$UseFishLenComps[,unused_years,,] <- 0
    if("FishLenComps_pop" %in% what) data$UseFishLenComps_pop[,,unused_years,,] <- 0
    if("FishAgeComps_discard" %in% what) data$UseFishAgeComps_discard[,unused_years,,] <- 0
    if("FishAgeComps_discard_pop" %in% what) data$UseFishAgeComps_discard_pop[,,unused_years,,] <- 0
    if("FishLenComps_discard" %in% what) data$UseFishLenComps_discard[,unused_years,,] <- 0
    if("FishLenComps_discard_pop" %in% what) data$UseFishLenComps_discard_pop[,,unused_years,,] <- 0
    if("SrvIdx" %in% what) data$UseSrvIdx[,unused_years,,] <- 0
    if("SrvIdx_pop" %in% what) data$UseSrvIdx_pop[,,unused_years,,] <- 0
    if("SrvAgeComps" %in% what) data$UseSrvAgeComps[,unused_years,,] <- 0
    if("SrvAgeComps_pop" %in% what) data$UseSrvAgeComps_pop[,,unused_years,,] <- 0
    if("SrvLenComps" %in% what) data$UseSrvLenComps[,unused_years,,] <- 0
    if("SrvLenComps_pop" %in% what) data$UseSrvLenComps_pop[,,unused_years,,] <- 0

  }

  # modify tagging stuff
  if(any(data$use_conv_fish_tagging == 1) && "conv_tagging" %in% what) {
    tags_to_remove <- which(data$conv_tag_release_indicator[,2] %in% unused_years)
    if(length(tags_to_remove) > 0) {
      data$conv_tagged_fish <- data$conv_tagged_fish[-tags_to_remove,,,,drop=FALSE]
      data$obs_conv_tag_fish_recap <- data$obs_conv_tag_fish_recap[,,-tags_to_remove,,,,,,drop=FALSE]
      data$conv_tag_release_indicator <- data$conv_tag_release_indicator[-tags_to_remove,,drop=FALSE]
      data$n_conv_tag_cohorts <- nrow(data$conv_tag_release_indicator)
    }
  }

  return(data)
}

#' Convert character or numeric input to numeric codes
#'
#' Maps a character vector to integer codes via a named lookup, or passes
#' numeric input through unchanged. Character arrays and matrices are converted
#' element-wise and keep their original dimensions. Unrecognised character
#' values raise an informative error listing both the invalid inputs and the
#' valid options.
#'
#' @param x Character vector, numeric vector, or array to convert.
#' @param lookup Named list (or named atomic vector) mapping valid character
#'   strings to numeric codes (e.g., \code{list("none" = 999, "multinomial" = 0)}).
#'
#' @return Numeric vector or array of the same shape as \code{x}, without names.
#'
#' @keywords internal
convert_to_numeric <- function(x, lookup) {

  # Return numeric if already numeric
  if(is.numeric(x)) {
    return(x)
  }

  # if character, validate against the lookup names and convert. Matching on
  # names is required because subsetting a list by a missing name yields a
  # one-element list holding NULL, which is.na() reports as FALSE.
  if(is.character(x)) {
    idx <- match(x, names(lookup))
    if(any(is.na(idx))) {
      invalid <- unique(x[is.na(idx)])
      stop("Invalid character input: ", paste(invalid, collapse = ", "),
           "\nValid options: ", paste(names(lookup), collapse = ", "))
    }
    result <- unlist(lookup[idx], use.names = FALSE)
    # arrays/matrices keep their shape
    if(is.array(x)) result <- array(result, dim = dim(x), dimnames = dimnames(x))
    return(result)
  }

  stop("Input must be numeric or character")
}

#' Build an index age-selection array
#'
#' Turns a per-fleet specification of which ages contribute to an index into
#' the \code{[age, fleet]} array of 0/1 weights the objective function uses.
#' Accepts a list with one element per fleet, where each element is a vector of
#' ages or \code{NULL} for all ages, or an array already in \code{[age, fleet]}
#' form.
#'
#' @param idx_ages List, array, or \code{NULL}. Per-fleet age selection.
#' @param n_ages Integer. Number of model ages.
#' @param n_fleets Integer. Number of fleets.
#' @param what Character. Name used in error messages.
#'
#' @return Array \code{[n_ages x n_fleets]} of 0/1 weights.
#'
#' @keywords internal
parse_idx_ages <- function(idx_ages, n_ages, n_fleets, what) {

  if(is.null(idx_ages)) return(array(1, dim = c(n_ages, n_fleets)))

  if(is.list(idx_ages)) {
    if(length(idx_ages) != n_fleets) stop(what, " is a list of length ", length(idx_ages), " but there are ", n_fleets, " fleets.")
    out <- array(1, dim = c(n_ages, n_fleets))
    for(f in seq_len(n_fleets)) {
      if(is.null(idx_ages[[f]])) next
      ages_f <- idx_ages[[f]]
      if(!all(ages_f %in% seq_len(n_ages))) stop(what, " for fleet ", f, " refers to ages outside 1:", n_ages, ".")
      out[,f] <- 0
      out[ages_f,f] <- 1
    } # end f loop
    if(any(colSums(out) == 0)) stop(what, " leaves at least one fleet with no ages contributing to its index.")
    return(out)
  }

  if(length(idx_ages) != n_ages * n_fleets) stop(what, " must have ", n_ages, " x ", n_fleets, " elements when supplied as an array.")
  out <- array(as.numeric(idx_ages), dim = c(n_ages, n_fleets))
  if(any(!out %in% c(0, 1))) stop(what, " must contain only 0 and 1.")
  if(any(colSums(out) == 0)) stop(what, " leaves at least one fleet with no ages contributing to its index.")
  return(out)
}

#' Validate fixed index covariance matrices for a multivariate normal likelihood
#'
#' Checks each supplied covariance once at setup and returns it ready for
#' \code{\link[RTMB]{dmvnorm}}. Each matrix must be square with one row per
#' observation the fleet actually fits, ordered the way the observations appear
#' when scanning the fleet's use flags in array order (region varies fastest,
#' then year, then season).
#'
#' The checks are not decorative. \code{RTMB::dmvnorm} reads only the lower
#' triangle without verifying symmetry, and returns \code{NaN} silently when the
#' covariance is not positive definite, so either mistake would otherwise
#' surface as an unexplained \code{NaN} objective rather than a setup error.
#'
#' @param cov_list List with one element per fleet, each either a covariance
#'   matrix or \code{NULL}, or \code{NULL} for no matrices at all.
#' @param like_type_vals Integer vector of index likelihood codes; only fleets
#'   coded \code{2} require a matrix.
#' @param use_arr Array \code{[region, year, season, fleet]} of use flags.
#' @param n_fleets Integer. Number of fleets.
#' @param what Character. Name used in error messages.
#'
#' @return List with one element per fleet, holding the validated covariance
#'   matrix for fleets using the multivariate normal and \code{NULL} otherwise.
#'
#' @keywords internal
parse_idx_cov <- function(cov_list, like_type_vals, use_arr, n_fleets, what) {

  out <- vector("list", n_fleets)

  for(f in seq_len(n_fleets)) {
    if(like_type_vals[f] != 2) next

    if(is.null(cov_list) || is.null(cov_list[[f]])) stop(what, " must supply a covariance matrix for fleet ", f, ", which uses a multivariate normal index likelihood.")

    cov_f <- as.matrix(cov_list[[f]])
    n_obs <- sum(array(use_arr[,,,f], dim = dim(use_arr)[1:3]) == 1)

    if(nrow(cov_f) != ncol(cov_f)) stop(what, " for fleet ", f, " is not square.")
    if(nrow(cov_f) != n_obs) stop(what, " for fleet ", f, " is ", nrow(cov_f), " x ", ncol(cov_f), " but that fleet fits ", n_obs, " observations.")
    if(!isTRUE(all.equal(cov_f, t(cov_f), check.attributes = FALSE))) stop(what, " for fleet ", f, " is not symmetric. Only its lower triangle would be used.")

    tryCatch(chol(cov_f), error = function(e) stop(what, " for fleet ", f, " is not positive definite."))
    dimnames(cov_f) <- NULL
    out[[f]] <- cov_f

  } # end f loop

  return(out)
}

#' Validate a selectivity fixed-effect centering penalty table
#'
#' Checks the table consumed by \code{\link{get_selex_fixed_penalty}} and
#' normalises its \code{par} column to a list of integer vectors, so a row may
#' name either one parameter or a whole set.
#'
#' @param selex_penalty Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par}, and \code{wt}, or \code{NULL}.
#' @param use_flag Integer (0/1). When \code{0} the table is returned unchanged
#'   and never validated, matching how the prior tables are guarded.
#' @param what Character. Name used in error messages.
#'
#' @return The validated table with \code{par} as a list column, or the input
#'   unchanged when \code{use_flag} is \code{0}.
#'
#' @keywords internal
validate_selex_penalty <- function(selex_penalty, use_flag, what) {

  if(is.null(use_flag) || use_flag != 1) return(selex_penalty)
  if(is.null(selex_penalty)) stop(what, " is NULL but its Use_ flag is 1. Please provide a penalty table.")

  required_cols <- c("region", "fleet", "block", "sex", "par", "wt")
  missing_cols <- setdiff(required_cols, names(selex_penalty))
  if(length(missing_cols) > 0) stop(what, " is missing required columns: ", paste(missing_cols, collapse = ", "))

  if(!is.list(selex_penalty$par)) selex_penalty$par <- as.list(selex_penalty$par)
  if(any(!is.finite(selex_penalty$wt)) || any(selex_penalty$wt < 0)) stop(what, "$wt must be finite and non-negative.")
  if(any(lengths(selex_penalty$par) == 0)) stop(what, "$par names an empty set of parameters in at least one row.")

  return(selex_penalty)
}

#' Validate the type column of a selectivity prior table
#'
#' Checks the optional \code{type} column consumed by
#' \code{\link{get_selex_prior}}. A table without the column is all
#' \code{"par"} rows (the original behaviour) and passes untouched.
#' \code{"value"} rows are range-checked here because their \code{par} column
#' indexes the selectivity grid rather than the parameter vector, and their
#' \code{block} must exist in the fleet's block map to resolve to a year.
#'
#' @param selex_prior Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par}, \code{mu}, \code{sd}, and optionally
#'   \code{type}, or \code{NULL}.
#' @param use_flag Integer (0/1). When \code{0} the table is returned unchanged
#'   and never validated, matching how the prior tables are guarded.
#' @param what Character. Name used in error messages.
#' @param sel_blocks Integer array \code{[region, year, fleet]} mapping model
#'   years to selectivity blocks.
#' @param n_bins Integer. Size of the grid the stream's selectivity is
#'   parameterized on (ages or lengths per its selectivity type).
#'
#' @return The validated table, or the input unchanged when \code{use_flag} is
#'   \code{0}.
#'
#' @keywords internal
validate_selex_prior_types <- function(selex_prior, use_flag, what, sel_blocks, n_bins) {

  if(is.null(use_flag) || use_flag != 1) return(selex_prior)
  if(is.null(selex_prior$type)) return(selex_prior)

  if(any(!selex_prior$type %in% c("par", "value"))) stop(what, "$type must be 'par' or 'value'.")

  val <- selex_prior$type == "value"
  if(any(!selex_prior$par[val] %in% 1:n_bins)) stop(what, "$par must index the selectivity grid (1:", n_bins, " for this stream's selectivity type) on rows with type = 'value'.")
  if(any(!is.finite(selex_prior$sd[val])) || any(selex_prior$sd[val] <= 0)) stop(what, "$sd must be finite and positive on rows with type = 'value'.")
  for(i in which(val)) {
    if(!selex_prior$block[i] %in% sel_blocks[selex_prior$region[i], , selex_prior$fleet[i]])
      stop(what, "$block does not exist in the block map for region ", selex_prior$region[i], ", fleet ", selex_prior$fleet[i], " on a row with type = 'value'.")
  } # end i loop

  return(selex_prior)
}

#' Set up bin-override selectivity deviations for one selectivity stream
#'
#' Creates the bin-override deviation parameter array, its factor map, and its
#' process-error hyperparameters, and records which bins each fleet overrides.
#' Bins named here take a free annual selectivity value instead of whatever the
#' fleet's functional form produces, which lets an otherwise parametric curve
#' carry a handful of freely estimated bins.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}.
#' @param bin_dev_bins List with one element per fleet, each a vector of bins to
#'   override or \code{NULL} for none, or \code{NULL} for no overrides anywhere.
#' @param pe_model Character vector \code{[n_fleets]} giving the process-error
#'   structure of the override deviations: \code{"none"}, \code{"iid"}, or
#'   \code{"rw"}.
#' @param prefix One of \code{"fish"}, \code{"ret"}, \code{"srv"}.
#' @param n_fleets Integer. Number of fleets in this stream.
#' @param bins Integer. Number of age or length bins.
#' @param starting_values Named list of user-supplied starting values.
#'
#' @return \code{input_list} with the parameter, map, and data entries added.
#'
#' @keywords internal
setup_sel_bin_devs <- function(input_list, bin_dev_bins, pe_model, prefix, n_fleets, bins, starting_values = list()) {

  dev_nm <- paste0("ln_", prefix, "sel_bin_devs")
  pe_nm <- paste0(prefix, "sel_bin_devs_pe_pars")
  bins_nm <- paste0(prefix, "_sel_bin_dev_bins")
  pe_data_nm <- paste0("cont_tv_", prefix, "sel_bin_devs")

  n_yrs_dev <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs

  if(!all(pe_model %in% c("none", "iid", "rw"))) stop(pe_data_nm, " must be one of none, iid, or rw")
  if(length(pe_model) != n_fleets) stop(pe_data_nm, " is not length ", n_fleets)
  pe_vals <- convert_to_numeric(pe_model, list(none = 0, iid = 1, rw = 2))

  # 0/1 array of which bins each fleet overrides; no overrides is all zeros
  bins_arr <- array(0, dim = c(bins, n_fleets))
  if(!is.null(bin_dev_bins)) {
    if(!is.list(bin_dev_bins) || length(bin_dev_bins) != n_fleets) stop(bins_nm, " must be a list with one element per fleet (use NULL for a fleet with no overrides).")
    for(f in seq_len(n_fleets)) {
      if(is.null(bin_dev_bins[[f]])) next
      if(!all(bin_dev_bins[[f]] %in% seq_len(bins))) stop(bins_nm, " for fleet ", f, " refers to bins outside 1:", bins)
      bins_arr[bin_dev_bins[[f]], f] <- 1
    } # end f loop
  }

  if(dev_nm %in% names(starting_values)) input_list$par[[dev_nm]] <- starting_values[[dev_nm]]
  else input_list$par[[dev_nm]] <- array(0, dim = c(input_list$data$n_regions, n_yrs_dev, bins, input_list$data$n_sexes, n_fleets))

  if(pe_nm %in% names(starting_values)) input_list$par[[pe_nm]] <- starting_values[[pe_nm]]
  else input_list$par[[pe_nm]] <- array(0, dim = c(input_list$data$n_regions, bins, input_list$data$n_sexes, n_fleets))

  # Only the named bins are estimated; everything else is fixed at zero and
  # never reaches the selectivity curve, so it costs nothing.
  map_dev <- array(NA_real_, dim = dim(input_list$par[[dev_nm]]))
  counter <- 1
  for(f in seq_len(n_fleets)) {
    for(r in seq_len(input_list$data$n_regions)) {
      for(s in seq_len(input_list$data$n_sexes)) {
        for(b in which(bins_arr[,f] == 1)) {
          map_dev[r, 1:n_yrs_dev, b, s, f] <- counter:(counter + n_yrs_dev - 1)
          counter <- counter + n_yrs_dev
        } # end b loop
      } # end s loop
    } # end r loop
  } # end f loop
  input_list$map[[dev_nm]] <- factor(map_dev)

  # A process-error hyperparameter only exists where the deviations are both
  # estimated and given a structure to be penalized against.
  map_pe <- array(NA_real_, dim = dim(input_list$par[[pe_nm]]))
  counter <- 1
  for(f in seq_len(n_fleets)) {
    if(pe_vals[f] == 0) next
    for(r in seq_len(input_list$data$n_regions)) {
      for(s in seq_len(input_list$data$n_sexes)) {
        for(b in which(bins_arr[,f] == 1)) {
          map_pe[r, b, s, f] <- counter
          counter <- counter + 1
        } # end b loop
      } # end s loop
    } # end r loop
  } # end f loop
  input_list$map[[pe_nm]] <- factor(map_pe)

  input_list$data[[bins_nm]] <- bins_arr
  input_list$data[[pe_data_nm]] <- pe_vals
  input_list$data[[paste0("map_", dev_nm)]] <- map_dev

  return(input_list)
}

#' Assign a value to every region x year cell of one fleet belonging to a selectivity block
#'
#' Selectivity block arrays are \code{[region, year, fleet]}, so a single fleet's slice is a
#' \code{region x year} MATRIX. \code{which(slice == block)} on that matrix returns LINEAR
#' positions running down the columns, in \code{1:(n_regions * n_years)} -- not year indices.
#' Using them as a year subscript (\code{arr[, which(...), fleet] <- value}) is therefore wrong
#' whenever \code{n_regions > 1}: it either errors with a subscript out of bounds, or, when the
#' block is early enough that the linear positions stay within \code{n_years}, SILENTLY writes
#' the wrong years. With three regions and 35 years, a block covering years 1-5 produces linear
#' positions 1-15 and quietly overwrites years 1-15. At \code{n_regions == 1} the linear
#' position equals the column index, which is why this only shows up in spatial models.
#'
#' Indexing with the logical matrix directly is correct in both cases, and stays correct if
#' blocks are ever allowed to differ between regions.
#'
#' @param arr Array \code{[region, year, fleet]} to write into.
#' @param blocks_arr Block array \code{[region, year, fleet]}, same first three dims as \code{arr}.
#' @param fleet Fleet index.
#' @param block Block value to match.
#' @param value Scalar to assign to the matching cells.
#' @return \code{arr} with the matching cells of \code{fleet} set to \code{value}.
#' @keywords internal
assign_sel_block <- function(arr, blocks_arr, fleet, block, value) {
  slice <- arr[, , fleet]
  slice[blocks_arr[, , fleet] == block] <- value
  arr[, , fleet] <- slice
  arr
}
