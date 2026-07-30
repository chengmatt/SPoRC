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
#'
#' @return Named numeric vector of length 6: \code{c(smooth_bin_curve,
#'   smooth_bin_diff, smooth_yr_diff, smooth_yr_curve, smooth_dome,
#'   smooth_mean_center)}.
#'
#' @keywords internal
resolve_sel_pen_wts <- function(pen_wts) {

  term_names <- c("smooth_bin_curve", "smooth_bin_diff", "smooth_yr_diff", "smooth_yr_curve", "smooth_dome", "smooth_mean_center")

  if(is.null(pen_wts)) {
    out <- stats::setNames(rep(0, length(term_names)), term_names)
    return(out)
  }

  if(is.null(names(pen_wts)) || !all(names(pen_wts) %in% term_names))
    stop("pen_wts must be a named numeric vector/list with names in: ", paste(term_names, collapse = ", "))

  out <- stats::setNames(rep(0, length(term_names)), term_names)
  out[names(pen_wts)] <- unlist(pen_wts)
  return(out)
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
#' Maps a character vector to integer codes via a named lookup list, or
#' passes numeric input through unchanged. Arrays and matrices are flattened,
#' converted element-wise, and restored to their original dimensions.
#' Unrecognised character values raise an informative error listing both the
#' invalid inputs and the valid options.
#'
#' @param x Character vector, numeric vector, or array to convert.
#' @param lookup Named list mapping valid character strings to numeric codes
#'   (e.g., \code{list("none" = 999, "multinomial" = 0)}).
#'
#' @return Numeric vector or array of the same shape as \code{x}.
#'
#' @keywords internal
convert_to_numeric <- function(x, lookup) {

  # Return numberic if already numeric
  if (is.numeric(x)) {
    return(x)
  }

  # if character, return numeric and convert
  if (is.character(x)) {
    result <- lookup[x]
    if (any(is.na(result))) {
      invalid <- x[is.na(lookup[x])]
      stop("Invalid character input: ", paste(invalid, collapse = ", "),
           "\nValid options: ", paste(names(lookup), collapse = ", "))
    }
    return(unlist(result))
  }

  # Handle arrays/matrices
  if (is.array(x)) {
    dims <- dim(x)
    result <- convert_to_numeric(as.vector(x), lookup)
    return(array(result, dim = dims))
  }

  stop("Input must be numeric or character")
}
