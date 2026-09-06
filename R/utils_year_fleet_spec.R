# Stage 1 of 3: model setup
#
# The year-by-fleet specification grammar, which twelve composition data sources each parsed with their own
# copy of. One parser means a source cannot gain a setting in a form the others do not accept.

#' Parse a year-by-fleet specification into a matrix
#'
#' Settings that change part way through a series are given as strings of the
#' form \code{Value_Year_x-y_Fleet_f}, where \code{x} and \code{y} are positions
#' in the year vector rather than calendar years and \code{y} may be
#' \code{"terminal"}. A vector of them describes a whole model, one entry per
#' block, and every year of every fleet must be covered by some entry.
#'
#' Later entries overwrite earlier ones where they overlap, so a general setting
#' can be given first and a period carved out of it afterwards.
#'
#' @param spec Character vector of specifications.
#' @param arg_name Name of the argument being parsed, used in error messages.
#' @param n_fleets Number of fleets the matrix must cover.
#' @param n_yrs Number of years the matrix must cover.
#' @param codes Named numeric vector mapping each accepted value to the code the
#'   model reads. The names are the accepted vocabulary, so a data source with its own
#'   set of values passes its own mapping.
#' @param check Optional function called as \code{check(value, fleet)} for each
#'   entry, returning a character string to stop with or \code{NULL} to accept.
#'   Used for constraints that depend on other settings, such as a likelihood
#'   that cannot take an aggregated observation.
#'
#' @return Numeric matrix \code{[n_yrs, n_fleets]} of codes.
#'
#' @keywords internal
parse_year_fleet_spec <- function(spec, arg_name, n_fleets, n_yrs, codes, check = NULL) {

  valid <- names(codes)
  grammar <- paste0("Give it as Value_Year_x-y_Fleet_f, where Value is one of ",
                    paste(valid, collapse = ", "), ", the years are positions in ",
                    "the year vector and the last year may be 'terminal', for ",
                    "example ", valid[1], "_Year_1-terminal_Fleet_1.")

  out <- array(NA_real_, dim = c(n_yrs, n_fleets))

  for(i in seq_along(spec)) {

    entry <- spec[i]
    parts <- unlist(strsplit(entry, "_"))

    if(length(parts) < 5) {
      stop(arg_name, " is '", entry, "', which is not a year and fleet ",
           "specification. ", grammar)
    }

    value <- parts[1]
    if(!value %in% valid) {
      stop(arg_name, " is '", entry, "'. ", grammar)
    }

    fleet <- suppressWarnings(as.numeric(parts[5]))
    if(is.na(fleet) || !fleet %in% seq_len(n_fleets)) {
      stop(arg_name, " is '", entry, "', naming fleet ", parts[5], " of ",
           n_fleets, " fleets. ", grammar)
    }

    # the year block, given as positions rather than calendar years
    if(grepl("terminal", entry, fixed = TRUE)) {
      first <- suppressWarnings(as.numeric(unlist(strsplit(parts[3], "-"))[1]))
      last <- n_yrs
    } else {
      bounds <- suppressWarnings(as.numeric(unlist(strsplit(parts[3], "-"))))
      first <- bounds[1]
      last <- if(length(bounds) > 1) bounds[2] else NA_real_
    }

    if(is.na(first) || is.na(last)) {
      stop(arg_name, " is '", entry, "', whose year range '", parts[3],
           "' could not be read. ", grammar)
    }
    if(first < 1 || last > n_yrs || first > last) {
      stop(arg_name, " is '", entry, "', naming years ", first, " to ", last,
           " of a ", n_yrs, " year model. The years are positions in the year ",
           "vector, so they run from 1 to ", n_yrs, ".")
    }

    if(!is.null(check)) {
      msg <- check(value, fleet)
      if(!is.null(msg)) stop(arg_name, " is '", entry, "'. ", msg)
    }

    out[first:last, fleet] <- codes[[value]]

  } # end i loop

  if(anyNA(out)) {
    gaps <- which(is.na(out), arr.ind = TRUE)
    stop(arg_name, " leaves ", nrow(gaps), " year and fleet combinations unset, ",
         "starting at year ", gaps[1, 1], " fleet ", gaps[1, 2],
         ". Every year of every fleet needs an entry covering it.")
  }

  out
}

#' Expand an at-age aggregation setting to a year by fleet matrix
#'
#' The at-age data sources accept their aggregation either as a bare value, standing
#' for the whole series, or as year and fleet specifications in the same form the
#' composition data sources take. Both arrive here and leave as a matrix, so
#' everything downstream reads one shape.
#'
#' @param type Character. Bare values, one for all fleets or one per fleet, or
#'   \code{Value_Year_x-y_Fleet_f} specifications.
#' @param n_fleets Number of fleets.
#' @param n_yrs Number of years.
#' @param arg_name Name of the argument being set, used in error messages.
#'
#' @return Numeric matrix \code{[n_yrs, n_fleets]} of codes.
#'
#' @keywords internal
at_age_type_matrix <- function(type, n_fleets, n_yrs, arg_name = "at-age Type") {

  codes_map <- c(agg = 0, spltRaggS = 1, aggRspltS = 2, spltRspltS = 3)
  valid <- names(codes_map)

  # a model already set up hands its settings back as codes rather than labels, which is what
  # conditioning a closed loop does, so those pass through rather than being read as labels
  if(is.numeric(type)) {
    if(!all(type %in% unname(codes_map))) {
      stop(arg_name, " has codes outside ", paste(unname(codes_map), collapse = ", "),
           ". A setting given as numbers is one this package produced, so a value ",
           "outside that set means it has been altered since.")
    }
    if(!is.null(dim(type)) && ncol(type) == n_fleets)
      return(base::matrix(as.numeric(type), nrow = nrow(type), ncol = n_fleets))
    return(base::matrix(rep(rep_len(as.numeric(type), n_fleets), each = n_yrs), nrow = n_yrs))
  }

  # the year grammar announces itself, so a bare value keeps meaning what it
  # always did: this setting, for the whole series
  if(any(grepl("_Year_", type, fixed = TRUE)))
    return(parse_year_fleet_spec(type, arg_name, n_fleets, n_yrs, codes_map))

  if(length(type) == 1) type <- rep(type, n_fleets)
  if(length(type) != n_fleets) {
    stop(arg_name, " has ", length(type), " entries for ", n_fleets, " fleets. Supply ",
         "one setting per fleet, one setting for all of them, or year and fleet ",
         "specifications such as '", valid[2], "_Year_1-terminal_Fleet_1'.")
  }
  if(!all(type %in% valid)) {
    stop(arg_name, " is '", paste(unique(type[!type %in% valid]), collapse = "', '"),
         "'. Valid options: ", paste(valid, collapse = ", "))
  }

  base::matrix(rep(unname(codes_map[type]), each = n_yrs), nrow = n_yrs)
}
