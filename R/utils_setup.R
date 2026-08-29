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
  # fleet its own, as surveys with different smoothing needs require.
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
  # The setup entry points open a fresh messages_list, but the helpers they call
  # are also reachable on their own, so start one here rather than failing on a
  # binding that has not been created yet.
  if(!exists("messages_list", inherits = TRUE)) messages_list <<- character(0)
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

#' Check every parameter block against the map that indexes it
#'
#' RTMB pairs a parameter with its map by position, so the two must be the same
#' length. They come apart when a starting value is supplied at the wrong shape,
#' or when a map is built from dimensions the parameter does not have. Neither is
#' caught where it happens: the objective reads the shorter of the two past its
#' end, and RTMB reports an invalid advector from somewhere unrelated.
#'
#' @param parameters Parameter list.
#' @param mapping Map list. Entries naming a parameter that is absent are
#'   reported too, since a map with no parameter is silently ignored.
#'
#' @return \code{invisible(NULL)}. Called for its error.
#'
#' @keywords internal
check_par_map_lengths <- function(parameters, mapping) {
  if(is.null(mapping) || length(mapping) == 0) return(invisible(NULL))

  problems <- character()
  for(nm in names(mapping)) {
    m <- mapping[[nm]]
    if(is.null(m) || length(m) == 0) next
    if(!nm %in% names(parameters)) {
      problems <- c(problems, sprintf("map '%s' has no parameter of that name", nm))
      next
    }
    n_par <- length(parameters[[nm]])
    if(length(m) != n_par)
      problems <- c(problems, sprintf("'%s' is length %d against a map of length %d",
                                      nm, n_par, length(m)))
  }

  if(length(problems) > 0)
    stop("parameters and their maps disagree on length. This is usually a starting ",
         "value supplied at the wrong shape.\n  ", paste(problems, collapse = "\n  "))

  invisible(NULL)
}

#' Substitute a user-supplied starting value, checking its shape first
#'
#' Every \code{Setup_Mod_*} stage builds a default for each parameter and then
#' lets \code{starting_values} replace it. The replacement used to be taken as
#' given. A value of the wrong shape is not rejected by that: it is carried into
#' the objective, read position by position, and indexes past its own end
#' somewhere far from the argument that caused it. What comes back is RTMB's
#' \code{'*this' is not a valid 'advector'}, which names nothing useful.
#'
#' The default already carries the shape the model expects, so it is the
#' reference the supplied value is measured against.
#'
#' @param default The parameter as the stage built it.
#' @param starting_values The user's list of starting values.
#' @param nm Name of the parameter.
#'
#' @return \code{starting_values[[nm]]} when it is supplied and correctly
#'   shaped, otherwise \code{default}.
#'
#' @keywords internal
use_starting_value <- function(default, starting_values, nm) {
  if(!nm %in% names(starting_values)) return(default)
  supplied <- starting_values[[nm]]

  # a default the stage has not built yet carries no shape to check against
  if(is.null(default) || length(default) == 0) return(supplied)

  shape <- function(x) if(is.null(dim(x))) paste0("length ", length(x)) else
    paste(dim(x), collapse = " by ")

  # A single value where the model wants many is the common slip and has an
  # obvious repair, so the message carries the call rather than only the
  # measurement.
  hint <- if(length(supplied) == 1 && length(default) > 1)
    paste0(" Recycle it with rep(", nm, ", ", length(default), ") if every element takes the same value.")
  else ""

  wrong <- if(!is.null(dim(default)) && !is.null(dim(supplied)))
    !identical(as.integer(dim(supplied)), as.integer(dim(default)))
  else length(supplied) != length(default)

  if(wrong)
    stop("starting value for ", nm, " is ", shape(supplied), " where the model expects ",
         shape(default), ".", hint)

  supplied
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

  # return "array" if just NULL and don't do anything with it
  if(is.null(arr) || length(arr) == 0) return(arr)

  # fill array based on specified option
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
#' element-wise and keep their original dimensions. Unrecognized character
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

#' Build a per-fleet bin-selection array
#'
#' Turns a per-fleet specification of which bins are used into the
#' \code{[bin, fleet]} array of 0/1 weights the objective function reads. Shared
#' by the ages that contribute to an index (\code{fish_idx_ages},
#' \code{srv_idx_ages}) and by the observed bins a composition is fitted over
#' (the \code{*_bins} arguments), so both spellings behave identically. Accepts a
#' list with one element per fleet, where each element is a vector of bin indices
#' or \code{NULL} for all bins, or an array already in \code{[bin, fleet]} form.
#'
#' @param idx_bins List, array, or \code{NULL}. Per-fleet bin selection.
#' @param n_bins Integer. Number of bins the selection indexes into: model ages
#'   for the index arguments, observed composition bins for the \code{*_bins}
#'   arguments.
#' @param n_fleets Integer. Number of fleets.
#' @param what Character. Name used in error messages.
#'
#' @return Array \code{[n_bins x n_fleets]} of 0/1 weights.
#'
#' @keywords internal
parse_bin_subset <- function(idx_bins, n_bins, n_fleets, what) {

  if(length(n_bins) != 1 || is.na(n_bins) || n_bins < 1) {
    stop(what, ": could not work out how many bins to index into. This is an internal error, not something the argument caused.")
  }
  if(is.null(idx_bins)) return(array(1, dim = c(n_bins, n_fleets)))

  if(is.list(idx_bins)) {
    if(length(idx_bins) != n_fleets) stop(what, " is a list of length ", length(idx_bins), " but there are ", n_fleets, " fleets.")
    out <- array(1, dim = c(n_bins, n_fleets))
    for(f in seq_len(n_fleets)) {
      if(is.null(idx_bins[[f]])) next
      bins_f <- idx_bins[[f]]
      if(!all(bins_f %in% seq_len(n_bins))) stop(what, " for fleet ", f, " refers to bins outside 1:", n_bins, ".")
      out[,f] <- 0
      out[bins_f,f] <- 1
    } # end f loop
    if(any(colSums(out) == 0)) stop(what, " leaves at least one fleet with no bins at all.")
    return(out)
  }

  if(!is.null(dim(idx_bins)) && length(dim(idx_bins)) == 2 && !all(dim(idx_bins) == c(n_bins, n_fleets))) {
    stop(what, " is ", paste(dim(idx_bins), collapse = " x "), " but must be ", n_bins, " x ", n_fleets,
         " (bins down the rows, fleets across the columns).")
  }
  if(length(idx_bins) != n_bins * n_fleets) stop(what, " must have ", n_bins, " x ", n_fleets, " elements when supplied as an array.")
  out <- array(as.numeric(idx_bins), dim = c(n_bins, n_fleets))
  if(any(!out %in% c(0, 1))) stop(what, " must contain only 0 and 1.")
  if(any(colSums(out) == 0)) stop(what, " leaves at least one fleet with no bins at all.")
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
#' normalizes its \code{par} column to a list of integer vectors, so a row may
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
#' \code{"par"} rows (the original behavior) and passes untouched.
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

#' Build the selectivity standardization window
#'
#' Records which bins the mean-one standardization averages over for each fleet.
#' Every bin is the default; a fleet whose catchability is defined against only
#' part of the bin range standardizes over that part instead.
#'
#' @param input_list Input list to append to
#' @param sel_norm_bins List with one element per fleet, or NULL for every bin
#' @param prefix One of "fish", "ret" or "srv"
#' @param n_fleets Number of fleets
#' @param bins Number of bins
#' @keywords internal
setup_sel_norm_bins <- function(input_list, sel_norm_bins, prefix, n_fleets, bins) {

  nm <- paste0(prefix, "_sel_norm_bins")

  # 0/1 array of which bins each fleet standardizes over; all ones is every bin
  arr <- array(1, dim = c(bins, n_fleets))
  if(!is.null(sel_norm_bins)) {
    if(!is.list(sel_norm_bins) || length(sel_norm_bins) != n_fleets) stop(nm, " must be a list with one element per fleet (use NULL for a fleet standardizing over every bin).")
    for(f in seq_len(n_fleets)) {
      if(is.null(sel_norm_bins[[f]])) next
      if(!all(sel_norm_bins[[f]] %in% seq_len(bins))) stop(nm, " for fleet ", f, " refers to bins outside 1:", bins)
      if(length(sel_norm_bins[[f]]) == 0) stop(nm, " for fleet ", f, " is empty; a standardization needs at least one bin.")
      arr[, f] <- 0
      arr[sel_norm_bins[[f]], f] <- 1
    } # end f loop
  }

  input_list$data[[nm]] <- arr
  return(input_list)
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

#' Sensible starting values for the double normal's peak
#'
#' The double normal carries the bin at which its plateau begins on the bin
#' scale, so a parameter left at zero puts the peak at bin zero. Where the first
#' bin is itself zero the ascending limb has no extent, its rescaling divides by
#' zero, and the curve comes back as \code{NaN}. A starting value in the middle
#' of the bin range is a curve, which is what a default should be.
#'
#' Only a sex whose slot holds a peak is seeded. Under a par sex offset the
#' slots of every sex beyond the first hold offsets on the first sex's
#' parameters, where zero is the right starting value and the middle of the bin
#' range would shift that sex's curve by half the bin range.
#'
#' @param pars Array \code{[region, par, block, sex, fleet]} of fixed-effect
#'   selectivity parameters.
#' @param sel_model_arr Integer array \code{[region, year, fleet]} of functional
#'   forms.
#' @param bin_vec Numeric vector of the bins selectivity is evaluated over.
#' @param sex_offset Character vector \code{[n_fleets]} of sex-offset
#'   specifications, used to tell a peak from an offset.
#'
#' @return \code{pars} with the peak slot of every double normal fleet set to
#'   the middle of the bin range, for the sexes that carry a peak there.
#'
#' @keywords internal
seed_dbnrml_peak <- function(pars, sel_model_arr, bin_vec, sex_offset = NULL) {
  if(!any(sel_model_arr == 4, na.rm = TRUE)) return(pars)
  mid <- min(bin_vec) + 0.5 * (max(bin_vec) - min(bin_vec))
  n_fleets <- dim(pars)[5]
  if(is.null(sex_offset)) sex_offset <- rep("none", n_fleets)
  for(f in seq_len(n_fleets)) {
    if(!any(sel_model_arr[,,f] == 4, na.rm = TRUE)) next
    # under a par offset only the first sex's slot is a peak; the rest are
    # offsets on it and belong at zero
    sexes <- if(sex_offset[f] %in% c("par", "par_scale", "par_apical")) 1 else seq_len(dim(pars)[4])
    pars[,1,,sexes,f] <- mid
  } # end f loop
  pars
} # end seed_dbnrml_peak

#' Validate the raw double normal limb flags
#'
#' @param x \code{NULL} or a 0/1 matrix \code{[n_fleets x 2]}.
#' @param n_fleets Number of fleets.
#' @param what Argument name for messages.
#'
#' @return An integer matrix \code{[n_fleets x 2]}, all zero for \code{NULL}.
#'
#' @keywords internal
setup_dbnrml_raw <- function(x, n_fleets, what) {
  if(is.null(x)) return(array(0, dim = c(n_fleets, 2)))
  x <- array(as.numeric(x), dim = c(n_fleets, 2))
  if(!all(x %in% c(0, 1))) stop(what, " must hold 0/1 values")
  x
}

#' Set up sex offsets on selectivity for one selectivity stream
#'
#' Parses the per-fleet sex-offset specification, stores the model flags, and
#' creates the curve scale-offset parameters with their factor map. Under a
#' \code{"par"} offset the sexes beyond the first store additive offsets on the
#' first sex's transformed-scale fixed-effect parameters, which needs no new
#' parameters, only the flag. Under a \code{"scale"} offset each sex beyond the
#' first carries a constant log-scale offset on its whole realized curve,
#' estimated per region and block.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}.
#' @param sex_offset Character vector \code{[n_fleets]}: \code{"none"},
#'   \code{"par"}, \code{"scale"}, \code{"par_scale"}, \code{"apical"}, or
#'   \code{"par_apical"}.
#' @param prefix One of \code{"fish"}, \code{"ret"}, or \code{"srv"}.
#' @param n_fleets Integer. Number of fleets in this stream.
#' @param fleet_label Character used in messages.
#' @param sel_model_arr Integer array \code{[region, year, fleet]} of functional
#'   forms, used to refuse a scale offset on forms whose standardization would
#'   cancel it (non-parametric forms 5 and 9).
#' @param cont_tv_mat Integer matrix \code{[region, fleet]} of continuous
#'   time-varying structures, used to refuse a scale offset under the
#'   semi-parametric structures (3-5) for the same reason.
#' @param max_blks Integer. Maximum number of selectivity blocks, sizing the
#'   scale parameter array.
#' @param sel_blocks Integer array \code{[region, year, fleet]} of selectivity
#'   block indices, so a scale offset is only estimated for the blocks a fleet
#'   actually has (the array is padded to \code{max_blks}).
#' @param fixed_spec Character vector \code{[n_fleets]} of fixed-parameter
#'   sharing specifications, or \code{NULL}. A \code{"par"} offset reads the
#'   later sexes' slots as offsets on the first sex's, so a specification that
#'   shares those slots across sexes (\code{"est_shared_s"},
#'   \code{"est_shared_r_s"}) would silently double the first sex's
#'   parameters and is refused.
#' @param starting_values Named list of user-supplied starting values.
#'
#' @return \code{input_list} with the flags, scale parameters, and map added.
#'
#' @keywords internal
setup_sel_sex_offset <- function(input_list, sex_offset, prefix, n_fleets, fleet_label,
                                 sel_model_arr, cont_tv_mat, max_blks, sel_blocks = NULL,
                                 fixed_spec = NULL, starting_values = list()) {

  scale_nm <- paste0("ln_", prefix, "sel_sex_scale")
  par_flag_nm <- paste0(prefix, "sel_sex_par_offset")
  scale_flag_nm <- paste0(prefix, "sel_sex_scale_offset")
  apical_flag_nm <- paste0(prefix, "sel_sex_apical_offset")

  if(length(sex_offset) != n_fleets) stop(prefix, "_sel_sex_offset is not length ", n_fleets)
  valid <- c("none", "par", "scale", "par_scale", "apical", "par_apical")
  if(!all(sex_offset %in% valid)) stop(prefix, "_sel_sex_offset must be one of ", paste(valid, collapse = ", "))
  if(input_list$data$n_sexes == 1 && any(sex_offset != "none")) stop(prefix, "_sel_sex_offset links sexes, so it requires n_sexes > 1")

  par_flag <- as.numeric(sex_offset %in% c("par", "par_scale", "par_apical"))
  scale_flag <- as.numeric(sex_offset %in% c("scale", "par_scale"))
  apical_flag <- as.numeric(sex_offset %in% c("apical", "par_apical"))

  # the apical offset is the height the double normal builds its limbs up to,
  # so it has no meaning for any other functional form
  for(f in seq_len(n_fleets)) {
    if(apical_flag[f] == 0) next
    if(any(sel_model_arr[,,f] != 4)) stop(prefix, "_sel_sex_offset for ", fleet_label, " ", f, " requests an apical offset, which is the height the double normal builds its limbs up to. That fleet is not on the double normal. Use the scale offset, which multiplies whatever curve the form returns.")
  } # end f loop

  # Under a par offset the later sexes' stored slots ARE the offsets, so a
  # specification sharing those slots across sexes would make them the first
  # sex's parameters and double them.
  if(!is.null(fixed_spec)) for(f in seq_len(n_fleets)) {
    if(par_flag[f] == 1 && fixed_spec[f] %in% c("est_shared_s", "est_shared_r_s"))
      stop(prefix, "_sel_sex_offset for ", fleet_label, " ", f, " requests a par offset, but its fixed-parameter specification shares the sex slots, which would read the first sex's parameters as their own offsets. Use est_all or est_shared_r with a par offset.")
  } # end f loop

  # A constant multiplier on the curve is canceled by the mean standardization
  # the non-parametric forms and semi-parametric structures apply afterwards.
  for(f in seq_len(n_fleets)) {
    if(scale_flag[f] == 0) next
    if(any(sel_model_arr[,,f] %in% c(5, 9))) stop(prefix, "_sel_sex_offset for ", fleet_label, " ", f, " requests a scale offset, but its non-parametric form is mean-standardized, which cancels a constant multiplier. Use the par offset instead.")
    if(any(cont_tv_mat[,f] %in% 3:5)) stop(prefix, "_sel_sex_offset for ", fleet_label, " ", f, " requests a scale offset, but its semi-parametric time variation is mean-standardized, which cancels a constant multiplier. Use the par offset instead.")
  } # end f loop

  if(scale_nm %in% names(starting_values)) input_list$par[[scale_nm]] <- starting_values[[scale_nm]]
  else input_list$par[[scale_nm]] <- array(0, dim = c(input_list$data$n_regions, max_blks, input_list$data$n_sexes, n_fleets))

  # The first sex is the reference and never carries a scale or apical offset, and a
  # fleet only carries one for the blocks it actually has
  map_scale <- array(NA_real_, dim = dim(input_list$par[[scale_nm]]))
  counter <- 1
  for(f in seq_len(n_fleets)) {
    if(scale_flag[f] == 0 && apical_flag[f] == 0) next
    for(r in seq_len(input_list$data$n_regions)) {
      n_blks_rf <- if(is.null(sel_blocks)) max_blks else length(unique(sel_blocks[r,,f]))
      for(b in seq_len(n_blks_rf)) {
        for(s in seq_len(input_list$data$n_sexes)[-1]) {
          map_scale[r, b, s, f] <- counter
          counter <- counter + 1
        } # end s loop
      } # end b loop
    } # end r loop
  } # end f loop
  input_list$map[[scale_nm]] <- factor(map_scale)

  input_list$data[[par_flag_nm]] <- par_flag
  input_list$data[[scale_flag_nm]] <- scale_flag
  input_list$data[[apical_flag_nm]] <- apical_flag
  for(f in seq_len(n_fleets)) if(sex_offset[f] != "none") collect_message("Selectivity sex offset for ", fleet_label, " ", f, " is: ", sex_offset[f])

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


#' Expand a fleet-specific ageing error specification to its full array
#'
#' Ageing error is a property of the sampling program, so a fishery that ages
#' its catch from otoliths and a survey that reads scales do not misclassify the
#' same way. \code{AgeingError_fish} and \code{AgeingError_srv} let each fleet
#' carry its own matrix, defaulting to the shared \code{AgeingError} so a model
#' written before they existed behaves exactly as it did.
#'
#' Every fleet must land on the same observed age bins, because the observed
#' composition arrays carry a single age dimension shared across fleets.
#'
#' @param x \code{NULL}, a 3D array \code{[n_ages x n_obs_ages x n_fleets]} for a
#'   time-invariant fleet-specific matrix, or a 4D array
#'   \code{[n_years x n_ages x n_obs_ages x n_fleets]} for a time-varying one.
#' @param shared The shared \code{[n_years x n_ages x n_obs_ages]} array to fall
#'   back on, already expanded over years.
#' @param n_fleets Integer. Number of fleets.
#' @param what Character. Argument name, used in messages and errors.
#'
#' @return Array \code{[n_years x n_ages x n_obs_ages x n_fleets]}.
#' @keywords internal
expand_fleet_ageing_error <- function(x, shared, n_fleets, what) {
  n_years <- dim(shared)[1]
  n_ages <- dim(shared)[2]
  n_obs_ages <- dim(shared)[3]
  out <- array(0, dim = c(n_years, n_ages, n_obs_ages, n_fleets))

  if(is.null(x)) {   # no fleet-specific matrices, so every fleet reads the shared one
    for(f in seq_len(n_fleets)) out[,,,f] <- shared
    return(out)
  }

  d <- dim(x)
  if(length(d) == 3) {          # time-invariant, fleet-specific
    if(d[3] != n_fleets) stop(what, " has ", d[3], " fleets but the model has ", n_fleets, ".")
    if(d[1] != n_ages || d[2] != n_obs_ages) {
      stop(what, " slices must be ", n_ages, " model ages by ", n_obs_ages, " observed ages, matching AgeingError, but are ", d[1], " by ", d[2], ".")
    }
    for(f in seq_len(n_fleets)) {
      check_bin_map(x[,,f], n_ages, paste0(what, " fleet ", f), strict = FALSE, tol = 0.05)
      for(i in seq_len(n_years)) out[i,,,f] <- x[,,f]
    } # end f loop
    collect_message(what, " is fleet specific and time invariant")
    return(out)
  }

  if(length(d) == 4) {          # time-varying, fleet-specific
    if(d[4] != n_fleets) stop(what, " has ", d[4], " fleets but the model has ", n_fleets, ".")
    if(d[1] != n_years) stop(what, " has ", d[1], " years but the model has ", n_years, ".")
    if(d[2] != n_ages || d[3] != n_obs_ages) {
      stop(what, " slices must be ", n_ages, " model ages by ", n_obs_ages, " observed ages, matching AgeingError, but are ", d[2], " by ", d[3], ".")
    }
    for(f in seq_len(n_fleets)) for(i in seq_len(n_years)) check_bin_map(x[i,,,f], n_ages, paste0(what, " fleet ", f, " year ", i), strict = FALSE, tol = 0.05)
    collect_message(what, " is fleet specific and time varying")
    return(array(as.numeric(x), dim = c(n_years, n_ages, n_obs_ages, n_fleets)))
  }

  stop(what, " must be NULL, a 3D array [n_ages x n_obs_ages x n_fleets], or a 4D array [n_years x n_ages x n_obs_ages x n_fleets].")
}

#' Read a fitted model's fleet-specific ageing error, falling back on the shared one
#'
#' The post-fit diagnostics have to reproduce the expected compositions the
#' likelihood built, which means reading the same ageing error the objective
#' read. Models fitted before \code{AgeingError_fish} and \code{AgeingError_srv}
#' existed carry only the shared matrix, so it is replicated across the fleets
#' and the diagnostics come out exactly as they did.
#'
#' @param data Data list from the fitted model.
#' @param shared The shared \code{[n_years x n_ages x n_obs_ages]} array.
#' @param which Either \code{"fish"} or \code{"srv"}.
#'
#' @return Array \code{[n_years x n_ages x n_obs_ages x n_fleets]}.
#' @keywords internal
fleet_ageing_error <- function(data, shared, which) {
  nm <- paste0("AgeingError_", which)
  n_fleets <- if(which == "srv") data$n_srv_fleets else data$n_fish_fleets
  if(!is.null(data[[nm]])) return(data[[nm]])
  out <- array(0, dim = c(dim(shared), n_fleets))
  for(f in seq_len(n_fleets)) out[,,,f] <- shared
  return(out)
}

#' A bin selection array, or NULL when it restricts nothing
#'
#' The composition machinery treats \code{NULL} as "fit every bin", which lets
#' the likelihood and the OSA packers skip the restriction entirely and lets a
#' backwards-compatible all-ones array of the wrong length never be indexed
#' into. Both the objective and \code{\link{get_osa}} decide that here, so they
#' cannot disagree about which bins were fitted.
#'
#' @param x A \code{[n_obs_bins x n_fleets]} 0/1 array, or \code{NULL}.
#'
#' @return \code{x}, or \code{NULL} if it is absent or selects every bin.
#' @keywords internal
bins_or_null <- function(x) {
  if(is.null(x) || all(x == 1)) return(NULL)
  return(x)
}

#' Validate a model-bin to observed-bin map
#'
#' \code{AgeingError} and \code{LenBinMap} are the same operation on different
#' axes: an \code{[n_model_bins x n_obs_bins]} matrix that the expected
#' composition is multiplied through so it lands on the bins the observations
#' were recorded on. The likelihood does not distinguish them, and neither does
#' this check, so a mistake in either one is reported the same way.
#'
#' A row is one model bin's share across the observed bins, so it sums to one. A
#' row of zeros is allowed and drops that model bin from the observations
#' entirely, which is how observed bins that start above the first model bin are
#' expressed (a shifted identity such as \code{diag(1, 10)[, 2:10]}).
#'
#' The row-sum tolerance is a caller's choice. Published ageing error matrices
#' are rounded at source, and real ones come in with rows summing to 0.997 or
#' 1.002; the likelihood renormalizes the expectation after the multiply, so a
#' row off by that much reweights nothing, and \code{AgeingError} passes
#' \code{tol = 0.05}. A length bin map is written by hand rather than read from a
#' rounded table, so \code{LenBinMap} keeps the \code{1e-8} it has always been
#' held to. Only a row off by more than \code{tol} is reported, since that means
#' the matrix is not the map its author thought it was.
#'
#' \code{strict} decides whether that is fatal. \code{LenBinMap} has always
#' rejected such a matrix outright and keeps doing so. \code{AgeingError} has
#' not been checked before, so a bad row is reported through the setup messages
#' rather than stopping a model that ran yesterday.
#'
#' A column of zeros is an observed bin nothing maps into, whose expected
#' proportion is a structural zero the composition likelihood cannot fit. It
#' follows \code{strict} for the same reason the row sums do. A negative entry is
#' fatal either way, since nothing downstream can interpret one.
#'
#' @param x The matrix to check.
#' @param n_model_bins Integer. Number of model bins, the required row count.
#' @param what Character. Argument name, used in messages.
#' @param strict Logical. \code{TRUE} (default) makes a bad row sum an error,
#'   \code{FALSE} reports it through \code{\link{collect_message}}.
#' @param tol Numeric. How far a row sum may sit from one before it is reported.
#'
#' @return \code{x} invisibly, as a matrix.
#' @keywords internal
check_bin_map <- function(x, n_model_bins, what, strict = TRUE, tol = 1e-8) {
  x <- as.matrix(x)
  if(nrow(x) != n_model_bins) stop(what, " must have one row per model bin (", n_model_bins, "), but has ", nrow(x), ".")
  if(any(x < 0)) stop(what, " must not contain negative values.")

  rs <- rowSums(x)
  bad <- which(abs(rs - 1) > tol & abs(rs) > 1e-8)
  if(length(bad) > 0) {
    msg <- paste0(what, " rows must each sum to one, spreading a model bin over the observed bins, or to zero to drop that model bin from the observations. Rows ",
                  paste(utils::head(bad, 10), collapse = ", "), " sum to neither (worst is ",
                  signif(rs[bad][which.max(abs(rs[bad] - 1))], 6), ").")
    if(strict) stop(msg)
    collect_message(msg, " Left as supplied, since the likelihood renormalizes the expectation after the multiply.")
  }

  empty <- which(colSums(x) <= 1e-8)
  if(length(empty) > 0) {
    msg <- paste0(what, " leaves observed bins ", paste(utils::head(empty, 10), collapse = ", "),
                  " with nothing mapped into them. Their expected proportion is a structural zero the composition likelihood cannot fit. Drop those bins from the observations instead, using the matching *_bins argument.")
    if(strict) stop(msg)
    collect_message(msg)
  }
  invisible(x)
}

#' Reject a bin restriction that leaves a stream nothing to fit
#'
#' A composition fitted over a single bin carries no information: the normalized
#' proportion in that bin is identically one whatever the model says. Every
#' family degenerates, and the machinery around them degenerates further. The
#' logistic-normal families spend one bin as the additive log-ratio reference and
#' so have no free element left, which gives a zero-length packed block, a
#' zero-row label frame, and a zero-length slice request in
#' \code{\link{eval_comp_osa}}. The discrete families mark their one bin as the
#' determined cell of the multinomial, which leaves \code{get_osa} with nothing
#' to keep and it fails inside \code{RTMB::oneStepPredict}.
#'
#' Two bins is therefore the minimum, and it is checked at setup where the
#' message can name the argument and the fleet. Fleets whose likelihood is
#' \code{"none"} are skipped, since their bins are never read.
#'
#' @param bins_arr \code{[n_obs_bins x n_fleets]} 0/1 array from
#'   \code{\link{parse_comp_bins}}.
#' @param like_vals Integer vector of likelihood codes, one per fleet.
#'   \code{999} marks a fleet that is not fitted. \code{NULL} checks every fleet.
#' @param what Character. Name of the bins argument, used in the error.
#'
#' @return \code{bins_arr} invisibly.
#' @keywords internal
check_comp_bins_min <- function(bins_arr, like_vals, what) {
  n_fleets <- ncol(bins_arr)
  for(f in seq_len(n_fleets)) {
    if(!is.null(like_vals) && f <= length(like_vals) && like_vals[f] == 999) next
    n_fit <- sum(bins_arr[,f])
    if(n_fit >= 2) next
    stop(what, " leaves fleet ", f, " with ", n_fit, " fitted bin",
         if(n_fit == 1) "" else "s",
         ". A composition needs at least two bins to say anything, since the proportion in a lone bin is one whatever the model predicts. Name two or more bins for that fleet.")
  } # end f loop
  invisible(bins_arr)
}

#' Drop blocks a bin restriction has emptied
#'
#' A restriction can leave a region, year and season with no observations at all
#' in the bins being fitted, even though the full composition had plenty. The
#' fitting likelihood already skips such a block, since normalizing it would
#' divide by zero, but the one-step-ahead packer has no way to know: it sees the
#' use flag say "there is data here", normalizes \code{(0 + addtocomp)} into a
#' flat composition and fits that. The packer and the evaluator cannot agree on
#' an emptiness test between themselves, because the evaluator never sees the
#' observations, so the two are reconciled here instead by clearing the use flag.
#'
#' Clearing it changes nothing about the fit: the likelihood was already
#' contributing zero for those blocks. It only stops the residual machinery
#' inventing an observation that was never there. Anything cleared is reported,
#' so a restriction that guts a stream is visible rather than silent.
#'
#' A block counts as empty when it holds no finite values at all, exactly as when
#' it sums to zero, since that is the test the likelihood's own guard applies.
#'
#' @param obs Observation array for the stream.
#' @param use Use-flag array for the stream, whose margins are \code{obs} without
#'   its bin and sex dimensions.
#' @param bins_arr \code{[n_obs_bins x n_fleets]} 0/1 array.
#' @param bin_dim Integer. Which dimension of \code{obs} holds the bins. The sex
#'   dimension is taken to be the next one, and fleets the last.
#' @param what Character. Stream name, used in the message.
#'
#' @return \code{use}, with emptied blocks cleared.
#' @keywords internal
drop_empty_fitted_blocks <- function(obs, use, bins_arr, bin_dim, what) {
  if(is.null(obs) || is.null(use) || is.null(bins_arr)) return(use)
  if(all(bins_arr == 1)) return(use)                 # nothing restricted, nothing to check
  d <- dim(obs)
  if(is.null(d) || length(d) < bin_dim + 1) return(use)
  if(d[bin_dim] != nrow(bins_arr)) return(use)       # sized for a different stream, left to the packers to reject

  n_dims <- length(d)
  margins <- setdiff(seq_len(n_dims), c(bin_dim, bin_dim + 1))   # everything but bins and sexes
  fleet_pos <- length(margins)                                    # fleets are last in obs, so last in margins
  cleared <- 0

  for(f in seq_len(ncol(bins_arr))) {
    fit <- which(bins_arr[,f] == 1)
    if(length(fit) == dim(obs)[bin_dim]) next                     # this fleet fits everything
    idx <- lapply(seq_len(n_dims), function(i) seq_len(d[i]))
    idx[[bin_dim]] <- fit
    idx[[n_dims]] <- f
    sub <- do.call(`[`, c(list(obs), idx, list(drop = FALSE)))
    # Same predicate the likelihood's own guard uses: a block with no finite values
    # counts as empty just as an all-zero one does, so the two cannot disagree
    tot <- apply(sub, margins, function(v) if(!any(is.finite(v))) 0 else sum(v, na.rm = TRUE))
    dim(tot) <- dim(tot)                                          # keep it an array for the assignment below
    uidx <- lapply(seq_along(dim(use)), function(i) seq_len(dim(use)[i]))
    uidx[[length(dim(use))]] <- f
    cur <- do.call(`[`, c(list(use), uidx, list(drop = FALSE)))
    hit <- which(as.vector(cur) == 1 & as.vector(tot) == 0)
    if(length(hit) == 0) next
    cur[hit] <- 0
    cleared <- cleared + length(hit)
    use <- do.call(`[<-`, c(list(use), uidx, list(value = cur)))
  } # end f loop

  if(cleared > 0) {
    collect_message(what, ": ", cleared, " region/year/season block",
                    if(cleared == 1) "" else "s",
                    " hold no observations inside the fitted bins, so their use flag was cleared. The likelihood already skipped them; this keeps the one-step-ahead residuals from inventing a flat composition in their place.")
  }
  return(use)
}

#' Re-reconcile use flags against freshly simulated observations
#'
#' \code{\link{drop_empty_fitted_blocks}} runs once at setup, against the
#' observations the model was built with. A self test or closed loop replaces
#' those observations replicate by replicate while carrying the setup's use flags
#' forward, so under a bin restriction a simulated replicate can hold nothing in
#' the fitted bins of a block the flags still call used. This walks the marginal
#' composition streams of a data list and reconciles them again.
#'
#' A no-op when no stream is restricted, which is the usual case.
#'
#' @param data A data list whose \code{Obs*} arrays have just been replaced.
#'
#' @return \code{data}, with its use flags reconciled.
#' @keywords internal
resync_fitted_blocks <- function(data) {
  streams <- list(
    list(obs = "ObsFishAgeComps", use = "UseFishAgeComps", bins = "FishAgeComps_bins", d = 4),
    list(obs = "ObsFishLenComps", use = "UseFishLenComps", bins = "FishLenComps_bins", d = 4),
    list(obs = "ObsSrvAgeComps", use = "UseSrvAgeComps", bins = "SrvAgeComps_bins", d = 4),
    list(obs = "ObsSrvLenComps", use = "UseSrvLenComps", bins = "SrvLenComps_bins", d = 4),
    list(obs = "ObsFishAgeComps_pop", use = "UseFishAgeComps_pop", bins = "FishAgeComps_pop_bins", d = 5),
    list(obs = "ObsFishLenComps_pop", use = "UseFishLenComps_pop", bins = "FishLenComps_pop_bins", d = 5),
    list(obs = "ObsSrvAgeComps_pop", use = "UseSrvAgeComps_pop", bins = "SrvAgeComps_pop_bins", d = 5),
    list(obs = "ObsSrvLenComps_pop", use = "UseSrvLenComps_pop", bins = "SrvLenComps_pop_bins", d = 5),
    list(obs = "ObsFish_caal", use = "UseFish_caal", bins = "Fish_caal_bins", d = 5),
    list(obs = "ObsSrv_caal", use = "UseSrv_caal", bins = "Srv_caal_bins", d = 5)
  )
  for(st in streams) {
    if(is.null(data[[st$bins]]) || all(data[[st$bins]] == 1)) next
    data[[st$use]] <- drop_empty_fitted_blocks(data[[st$obs]], data[[st$use]], data[[st$bins]], st$d, st$obs)
  } # end st loop
  return(data)
}

#' Number of observed bins a composition stream is recorded on
#'
#' The \code{*_bins} arguments index into observed bins, so they need the bin
#' count of the array they will be applied to. That is normally read straight
#' off the supplied observation array, but a model carrying no data for a stream
#' can hand in an array with no dimensions at all, so the model's own observed
#' bin count stands in: the ageing error's observed-age dimension for ages, and
#' \code{\link{obs_len_bins}} for lengths.
#'
#' @param input_list Input list, used for the fallback.
#' @param obs The observation array for the stream, possibly dimensionless.
#' @param dim_i Integer. Which dimension of \code{obs} holds the bins.
#' @param axis Either \code{"age"} or \code{"len"}.
#'
#' @return A single positive integer.
#' @keywords internal
obs_bin_count <- function(input_list, obs, dim_i, axis) {
  d <- dim(obs)
  if(!is.null(d) && length(d) >= dim_i && !is.na(d[dim_i]) && d[dim_i] >= 1) return(d[dim_i])
  if(axis == "len") return(obs_len_bins(input_list))
  ae <- input_list$data$AgeingError
  if(is.null(ae) || is.null(dim(ae))) return(length(input_list$data$ages))
  return(dim(ae)[length(dim(ae))])   # last dimension is the observed ages, 2D or 3D alike
}

#' Parse and report the observed bins a composition stream is fitted over
#'
#' Wraps \code{\link{parse_bin_subset}} and records which fleets ended up
#' restricted, so the setup messages name the bins a stream is fitted over rather
#' than leaving it implicit. Used for every \code{*_bins} argument, age and
#' length, retained and discarded, marginal and conditional, so a restriction
#' reads the same way whichever stream it was set on.
#'
#' @param bins List, array, or \code{NULL}. Per-fleet bin selection.
#' @param n_bins Integer. Number of observed bins the stream is recorded on,
#'   that is after any ageing error or length bin map.
#' @param n_fleets Integer. Number of fleets.
#' @param what Character. Argument name, used in messages and errors.
#'
#' @return Array \code{[n_bins x n_fleets]} of 0/1 weights.
#' @keywords internal
parse_comp_bins <- function(bins, n_bins, n_fleets, what) {
  out <- parse_bin_subset(bins, n_bins, n_fleets, what)
  for(f in seq_len(n_fleets)) {
    if(sum(out[,f]) != n_bins) collect_message(what, " for fleet ", f, " is fitted over observed bins: ", paste(which(out[,f] == 1), collapse = ", "))
  } # end f loop
  return(out)
}

#' Number of length bins the observed compositions are recorded on
#'
#' The model's own bins unless a length bin map was given, in which case the
#' observed compositions sit on the map's columns.
#'
#' @param input_list Input list after \code{\link{Setup_Mod_Biologicals}}.
#' @keywords internal
obs_len_bins <- function(input_list) {
  if(is.null(input_list$data$LenBinMap)) length(input_list$data$lens) else ncol(input_list$data$LenBinMap)
}


#' Validate the double normal start bin per fleet
#'
#' @param x \code{NULL} or an integer vector of start bins, one per fleet.
#' @param n_fleets Number of fleets.
#' @param n_bins Number of selectivity bins.
#' @param what Argument name for messages.
#' @return Integer vector \code{[n_fleets]}, ones when \code{x} is \code{NULL}.
#' @keywords internal
setup_dbnrml_startbin <- function(x, n_fleets, n_bins, what) {
  if(is.null(x)) return(rep(1, n_fleets))
  if(length(x) != n_fleets) stop(what, " must have one entry per fleet (", n_fleets, ")")
  if(any(x < 1 | x > n_bins | x != round(x))) stop(what, " entries must be bin indices between 1 and ", n_bins)
  as.integer(x)
}
