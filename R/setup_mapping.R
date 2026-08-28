# Stage 1 of 3: model setup
#
# Parameter map builders shared by more than one Setup_* file. build_pe_map and
# build_shared_spec_map turn a dimension sharing spec string into an RTMB factor
# map; the do_*_mapping functions here each serve several parameter blocks at
# once, selected by a prefix argument. A mapping helper belongs in this file only
# if more than one Setup_* file calls it; single caller helpers stay with their
# caller.

#' Build a generic sharing/process-error parameter map
#'
#' Assigns a unique estimation ID to every combination of the "key"
#' dimensions (those \emph{not} listed in \code{share_over}), so cells that
#' agree on the key dimensions get the same ID regardless of their value
#' along the \code{share_over} dimensions. This is the general form behind
#' every \code{"est_shared_<dims>"} spec used throughout SPoRC's
#' \code{Setup_*} mapping functions (e.g. \code{est_shared_r},
#' \code{est_shared_r_seas}, ...): rather than hand-enumerating one branch
#' per combination of dimensions, callers just say which dimensions to
#' collapse.
#'
#' @param dims Named integer vector of array dimensions, e.g.
#'   \code{c(region = 3, season = 4, fleet = 2)}. Names must be unique and
#'   non-empty.
#' @param share_over Character vector, subset of \code{names(dims)}, giving
#'   the dimensions across which a single parameter is shared. Use
#'   \code{character(0)} (default) to estimate a unique parameter per cell
#'   (equivalent to \code{"est_all"}), or \code{names(dims)} to share one
#'   value across everything.
#'
#' @return Integer array with \code{dim = dims} and \code{names(dim(.))
#'   = names(dims)}, containing sequential estimation IDs from 1 to the
#'   number of unique groups. No \code{NA} handling is done here; apply
#'   fixing (e.g. \code{"fix"} or use-flags) to the result afterwards.
#'
#' @keywords internal
build_pe_map <- function(dims, share_over = character(0)) {

  if(is.null(names(dims)) || any(names(dims) == "")) stop("build_pe_map: 'dims' must be a fully named vector.")
  if(!all(share_over %in% names(dims))) stop("build_pe_map: 'share_over' must be a subset of names(dims).")

  key_dims <- setdiff(names(dims), share_over)
  grid <- expand.grid(lapply(dims, seq_len), KEEP.OUT.ATTRS = FALSE)

  key <- if(length(key_dims) == 0) rep(1L, nrow(grid)) else do.call(paste, c(grid[key_dims], sep = "_"))
  id <- match(key, unique(key))

  arr <- array(id, dim = dims)
  names(dim(arr)) <- names(dims)
  arr
}

#' Build a factor map from an "est_all"/"fix"/"est_shared_..." spec string
#'
#' Convenience wrapper around \code{\link{build_pe_map}} for the common case
#' of a single spec string (as used by e.g. \code{sigmaC_spec},
#' \code{sigmaF_spec}) governing a fixed-effect array with no additional
#' use/fix masking. Validates \code{spec} against every dimension
#' combination implied by \code{dim_abbrev} before building the map, so
#' invalid specs still fail with the same style of error message as the
#' hand-written mapping functions.
#'
#' @param dims Named integer vector of array dimensions (see
#'   \code{\link{build_pe_map}}).
#' @param spec Character scalar: \code{"est_all"}, \code{"fix"}, or
#'   \code{"est_shared_<abbrev>[_<abbrev>...]"}.
#' @param dim_abbrev Named character vector mapping abbreviation to
#'   dimension name, given in canonical order, e.g.
#'   \code{c(r = "region", y = "year", seas = "season", f = "fleet")}.
#'   Values must match \code{names(dims)} exactly.
#' @param use Optional array of the same dimensions as \code{dims}, non-zero
#'   where an observation informs that cell (e.g. \code{UseCatch},
#'   \code{UseSrvIdx}). When supplied, the resulting map is checked for
#'   observation-error parameters that no data can identify. See
#'   \code{\link{check_spec_map_identifiable}}. Defaults to \code{NULL}, which
#'   skips the check and leaves behavior unchanged.
#' @param what Character label for the parameter, used in the check's messages.
#' @param min_obs Integer. A group informed by fewer than this many
#'   observations raises an error. Default \code{2}.
#' @param warn_obs Integer. A group informed by fewer than this many
#'   observations raises a warning. Default \code{5}.
#'
#' @return Factor vector of length \code{prod(dims)}, suitable for direct
#'   assignment to \code{input_list$map$<par>}.
#'
#' @keywords internal
build_shared_spec_map <- function(dims, spec, dim_abbrev, use = NULL,
                                  what = "parameter", min_obs = 2, warn_obs = 5) {

  shared_specs <- unlist(lapply(seq_along(dim_abbrev), function(k) {
    combs <- utils::combn(names(dim_abbrev), k)
    apply(combs, 2, function(x) paste0("est_shared_", paste(x, collapse = "_")))
  }))
  valid_specs <- c("fix", "est_all", shared_specs)

  if(!spec %in% valid_specs) stop("spec '", spec, "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  if(spec == "fix") return(factor(rep(NA, prod(dims))))

  share_over <- if(spec == "est_all") character(0) else {
    parts <- strsplit(sub("^est_shared_", "", spec), "_")[[1]]
    unname(dim_abbrev[parts])
  }

  map <- factor(as.vector(build_pe_map(dims, share_over = share_over)))

  if(!is.null(use)) check_spec_map_identifiable(map, use, spec, dims, dim_abbrev, what, min_obs, warn_obs)

  map
}

#' Check that an observation-error spec leaves every parameter identified
#'
#' An observation-error standard deviation needs more than one observation to
#' be estimable. Given one, the likelihood is unbounded: the standard deviation
#' collapses towards zero on whatever residual the model can fit exactly and the
#' \code{log(sigma)} term runs to negative infinity. The optimizer reports
#' convergence, so nothing about the fit announces the problem.
#'
#' This is not confined to \code{"est_all"}. Any spec that leaves a dimension
#' free when the corresponding observation array is one cell deep in the other
#' dimensions has the same failure, so which specs are safe depends on the
#' model's dimensions rather than on the spec string alone. The check therefore
#' counts the observations actually informing each estimation group instead of
#' rejecting particular spec names.
#'
#' @param map Factor map returned by \code{\link{build_shared_spec_map}}.
#' @param use Array of the same dimensions as the parameter, non-zero where an
#'   observation informs that cell.
#' @param spec The spec string, used in messages.
#' @param dims Named integer vector of array dimensions.
#' @param dim_abbrev Named character vector of dimension abbreviations.
#' @param what Character label for the parameter, used in messages.
#' @param min_obs Integer below which a group raises an error.
#' @param warn_obs Integer below which a group raises a warning.
#'
#' @return \code{invisible(NULL)}. Called for its error and warning side effects.
#'
#' @keywords internal
check_spec_map_identifiable <- function(map, use, spec, dims, dim_abbrev,
                                        what = "parameter", min_obs = 2, warn_obs = 5) {

  if(length(use) != length(map)) {
    warning("Could not check identifiability of ", what, ": the 'use' array has ",
            length(use), " cells but the parameter has ", length(map),
            ". Skipping the check.")
    return(invisible(NULL))
  }

  n_obs <- tapply(as.numeric(use) != 0, map, sum)
  n_obs <- n_obs[!is.na(n_obs)]
  if(length(n_obs) == 0) return(invisible(NULL))

  # A spec that shares over every dimension always collapses to one parameter,
  # so it is the safe fallback to point people at.
  all_shared <- paste0("est_shared_", paste(names(dim_abbrev), collapse = "_"))

  if(any(n_obs < min_obs)) {
    stop(what, " is specified as '", spec, "', which creates ", length(n_obs),
         " parameters for ", sum(as.numeric(use) != 0), " observations. ",
         sum(n_obs < min_obs), " of them are informed by fewer than ", min_obs,
         " observations, so they cannot be estimated: an observation error ",
         "standard deviation with a single observation drives the likelihood to ",
         "negative infinity rather than failing outright. Share this parameter ",
         "over more dimensions (for example '", all_shared, "'), or fix it with 'fix'.")
  }

  if(any(n_obs < warn_obs)) {
    warning(what, " is specified as '", spec, "', which leaves ",
            sum(n_obs < warn_obs), " of ", length(n_obs), " parameters informed by ",
            "fewer than ", warn_obs, " observations. These are estimable but poorly ",
            "determined, and the smallest is informed by ", min(n_obs), ". Consider ",
            "sharing over more dimensions or fixing it.")
  }

  invisible(NULL)
}

#' Dimensions of one at-age observation array
#'
#' Every at-age stream is stored region by year by season by age by sex by fleet,
#' with a leading population dimension for the population-specific form. This is
#' the layout of the prediction arrays the likelihood reads, so the two line up
#' margin for margin.
#'
#' @param input_list Named list with \code{$data}.
#' @param fleet_field \code{"n_fish_fleets"} or \code{"n_srv_fleets"}.
#' @param pop Logical. \code{TRUE} for the population-specific stream.
#'
#' @return An integer vector of dimensions.
#'
#' @keywords internal
at_age_dims <- function(input_list, fleet_field, pop = FALSE) {
  d <- c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
         length(input_list$data$ages), at_age_n_sexes(input_list), input_list$data[[fleet_field]])
  n_pop <- input_list$data$n_pop
  return(as.integer(if(pop) c(if(is.null(n_pop)) 1L else n_pop, d) else d))
}

#' How many sexes an input list carries
#'
#' A list assembled without \code{\link{Setup_Mod_Dim}}, as the mapping unit
#' tests do, has no sex dimension to read. One sex is what such a list means, and
#' is what the at-age arrays were shaped as before the margin existed.
#'
#' @param input_list Named list with \code{$data}.
#'
#' @return An integer.
#'
#' @keywords internal
at_age_n_sexes <- function(input_list) {
  n_sexes <- input_list$data$n_sexes
  return(if(is.null(n_sexes)) 1L else n_sexes)
}

#' Refuse an at-age array or parameter that is missing a dimension
#'
#' The at-age streams carry a sex margin, and nothing promotes an array into it:
#' an array one dimension short would otherwise be indexed by position and read
#' the wrong age or sex. This reports what was supplied against what is wanted.
#'
#' @param x The array to check, or \code{NULL} to skip.
#' @param want Integer vector of the dimensions expected.
#' @param what Name used in the message.
#'
#' @return \code{invisible(NULL)}. Called for its error.
#'
#' @keywords internal
check_at_age_shape <- function(x, want, what) {

  if(is.null(x)) return(invisible(NULL))
  got <- dim(as.array(x))

  if(is.null(got) || !identical(as.integer(got), as.integer(want))) {
    stop(what, " is ", if(is.null(got)) paste0("a length ", length(x), " vector")
                       else paste(got, collapse = " by "),
         " where ", paste(want, collapse = " by "), " is expected. The at-age streams carry ",
         "a sex margin, so supply the full array rather than one summed over sexes.")
  }

  return(invisible(NULL))
}

#' Store one at-age stream's observations, use flags and standard errors
#'
#' Shared by the four streams and their population-specific forms. An absent
#' stream is given a zeroed array so the objective can index it unconditionally.
#'
#' @param input_list Named list with \code{$data}.
#' @param obs,use,se Observation, use and reported standard error arrays, any of
#'   them \code{NULL}.
#' @param stream Stream tag naming the data elements, e.g. \code{"CatchAA"}.
#' @param fleet_field \code{"n_fish_fleets"} or \code{"n_srv_fleets"}.
#' @param pop Logical. \code{TRUE} for the population-specific stream.
#'
#' @return \code{input_list} with the stream's three data elements set.
#'
#' @keywords internal
do_at_age_data_setup <- function(input_list, obs, use, se, stream, fleet_field, pop = FALSE) {

  n_sexes <- at_age_n_sexes(input_list)
  tag <- if(pop) paste0(stream, "_pop") else stream
  dims <- at_age_dims(input_list, fleet_field, pop)

  supplied <- list(obs = !is.null(obs), use = !is.null(use), se = !is.null(se))

  if(is.null(obs)) {
    obs <- array(0, dim = dims)
    use <- array(0, dim = dims)
  } else if(is.null(use)) {
    stop("Obs", tag, " was supplied without Use", tag, ".")
  }
  if(is.null(se)) se <- array(0, dim = dims)

  # only what the caller supplied is checked: the zeroed placeholders built
  # above are correct by construction, and a model carrying no at-age data need
  # not have an age or sex dimension for them to be checked against
  fleet_args <- list(input_list$data[[fleet_field]])
  names(fleet_args) <- fleet_field
  check_one <- function(x, what) {
    do.call(check_data_dimensions, c(list(x, n_regions = input_list$data$n_regions,
                                          n_years = length(input_list$data$years),
                                          n_seas = input_list$data$n_seas,
                                          n_ages = length(input_list$data$ages),
                                          n_sexes = n_sexes,
                                          n_pop = input_list$data$n_pop,
                                          what = what), fleet_args))
  }
  if(supplied$obs) {
    check_one(obs, paste0("Obs", tag))
    check_one(use, paste0("Use", tag))
  }
  if(supplied$se) check_one(se, paste0("Obs", tag, "_SE"))

  input_list$data[[paste0("Obs", tag)]] <- obs
  input_list$data[[paste0("Use", tag)]] <- use
  input_list$data[[paste0("Obs", tag, "_SE")]] <- se

  return(input_list)
}

#' Should an at-age stream's observation error parameter be estimated?
#'
#' A stream nobody fits has nothing to inform its standard deviation, and a
#' stream taking its error from reported standard errors alone has no parameter
#' to read, so both are held fixed whatever the spec says.
#'
#' @param spec \code{"est"} or \code{"fix"} as supplied by the caller.
#' @param form The stream's \code{sigma_form}.
#' @param any_used \code{TRUE} when any fleet fits this stream.
#'
#' @return \code{"est"} or \code{"fix"}.
#'
#' @keywords internal
at_age_sigma_spec <- function(spec, form, any_used) {
  forms <- convert_to_numeric(form, list(none = 0, data = 1, est_additive = 2, est_quadrature = 3))
  if(all(forms == 1)) return("fix")
  if(!any_used) return("fix")
  return(spec)
}

#' Name the split margins of one at-age stream
#'
#' An at-age observation is stored over regions and sexes whether or not it is
#' reported that way, and the Type code names which of those margins the fleet
#' reports separately. The vocabulary is the composition one, so a model stating
#' both kinds of data states them the same way.
#'
#' A margin the fleet sums over carries its observation in slot one, and a use
#' flag anywhere else on that margin is refused rather than quietly ignored.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#' @param type Character, one of \code{"agg"}, \code{"spltRaggS"},
#'   \code{"aggRspltS"} or \code{"spltRspltS"}, either one setting for every
#'   fleet or one per fleet.
#' @param stream Stream tag naming the data element, e.g. \code{"CatchAA"}.
#' @param fleet_field \code{"n_fish_fleets"} or \code{"n_srv_fleets"}.
#' @param use_field Name of the use array for this stream.
#' @param pop Logical. \code{TRUE} for the population-specific stream.
#'
#' @return \code{input_list} with \code{$data$<stream>_Type} set.
#'
#' @keywords internal
do_at_age_type_setup <- function(input_list, type, stream, fleet_field, use_field, pop = FALSE) {

  valid <- c("agg", "spltRaggS", "aggRspltS", "spltRspltS")
  n_fleets <- input_list$data[[fleet_field]]
  tag <- if(pop) paste0(stream, "_pop") else stream

  if(length(type) == 1) type <- rep(type, n_fleets)
  if(length(type) != n_fleets) {
    stop(tag, "_Type has ", length(type), " entries for ", n_fleets, " fleets. Supply ",
         "one setting per fleet, or one setting for all of them.")
  }
  if(!all(type %in% valid)) {
    stop(tag, "_Type is '", paste(unique(type[!type %in% valid]), collapse = "', '"),
         "'. Valid options: ", paste(valid, collapse = ", "))
  }

  codes <- convert_to_numeric(type, list(agg = 0, spltRaggS = 1, aggRspltS = 2, spltRspltS = 3))
  input_list$data[[paste0(tag, "_Type")]] <- codes

  # a summed margin holds the observation in slot one, so nothing else on that
  # margin may be flagged: the prediction has already added those cells in
  use_arr <- input_list$data[[use_field]]
  nd <- length(dim(use_arr))
  i_r <- if(pop) 2 else 1
  i_s <- nd - 1
  by_region <- apply(use_arr, c(i_r, nd), function(x) sum(x != 0))
  by_sex <- apply(use_arr, c(i_s, nd), function(x) sum(x != 0))

  for(f in seq_len(n_fleets)) {
    split <- at_age_split(codes[f])
    if(!split$region && dim(use_arr)[i_r] > 1 && any(by_region[-1,f] != 0)) {
      stop(tag, " fleet ", f, " is '", type[f], "', which sums over regions, but it flags ",
           "observations in more than one region. An observation summed over regions is ",
           "one number, so it belongs in region 1. Use a type splitting regions if the ",
           "fleet reports them separately.")
    }
    if(!split$sex && dim(use_arr)[i_s] > 1 && any(by_sex[-1,f] != 0)) {
      stop(tag, " fleet ", f, " is '", type[f], "', which sums over sexes, but it flags ",
           "observations in more than one sex. An observation summed over sexes is one ",
           "number, so it belongs in sex 1. Use a type splitting sexes if the fleet ",
           "reports them separately.")
    }
  } # end f loop

  return(input_list)
}

#' Set the error scale and likelihood of one at-age stream
#'
#' An at-age observation may be lognormal or normal, and its standard deviation
#' may come from an estimated parameter, from reported standard errors, or from
#' both. This is the parity the aggregated index streams already have, stated per
#' fleet.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#' @param like_type \code{"lognormal"} or \code{"normal"}, one setting for every
#'   fleet or one per fleet.
#' @param sigma_form \code{"none"} for the parameter alone, \code{"data"} for the
#'   reported standard errors alone, \code{"est_additive"} or
#'   \code{"est_quadrature"} for both.
#' @param stream,fleet_field,pop See \code{\link{do_at_age_type_setup}}.
#'
#' @return \code{input_list} with \code{$data$<stream>_LikeType} and
#'   \code{$data$<stream>_sigma_form} set.
#'
#' @keywords internal
do_at_age_like_setup <- function(input_list, like_type, sigma_form, stream, fleet_field, pop = FALSE) {

  n_fleets <- input_list$data[[fleet_field]]
  tag <- if(pop) paste0(stream, "_pop") else stream

  expand <- function(x, nm, valid) {
    if(length(x) == 1) x <- rep(x, n_fleets)
    if(length(x) != n_fleets) {
      stop(tag, nm, " has ", length(x), " entries for ", n_fleets, " fleets. Supply one ",
           "setting per fleet, or one setting for all of them.")
    }
    if(!all(x %in% valid)) {
      stop(tag, nm, " is '", paste(unique(x[!x %in% valid]), collapse = "', '"),
           "'. Valid options: ", paste(valid, collapse = ", "))
    }
    return(x)
  }

  like_type <- expand(like_type, "_LikeType", c("lognormal", "normal"))
  sigma_form <- expand(sigma_form, "_sigma_form", c("none", "data", "est_additive", "est_quadrature"))

  input_list$data[[paste0(tag, "_LikeType")]] <- convert_to_numeric(like_type, list(lognormal = 0, normal = 1))
  input_list$data[[paste0(tag, "_sigma_form")]] <- convert_to_numeric(
    sigma_form, list(none = 0, data = 1, est_additive = 2, est_quadrature = 3))

  return(input_list)
}

#' Set the correlation structure for one at-age stream
#'
#' Each stream is configured where its data are configured, so the catch and
#' discard streams are set in \code{\link{Setup_Mod_Catch_and_F}} and the index
#' streams in their own setup functions. The population-specific form carries its
#' own setting rather than borrowing the aggregated one.
#'
#' Four structures are available, per fleet. \code{"iid"} treats ages as
#' independent. \code{"1dar1"} correlates them as an AR(1) in age distance, so a
#' fleet skipping ages is spaced correctly rather than treated as consecutive.
#' \code{"us"} estimates an unstructured correlation across ages, the third
#' structure ICES age-structured assessments offer. \code{"2dar1"} correlates
#' over ages and years jointly through a separable AR(1), which is defined on a
#' complete grid and so requires the fleet's observed ages and years to form one.
#'
#' How the correlations are shared follows the package's spec strings rather than
#' a structure of its own. They sit over region, sex and fleet, with a leading
#' population margin for the population-specific streams, so
#' \code{"est_shared_r_s"} (the default, \code{"est_shared_p_r_s"} for the
#' population form) gives one per fleet, \code{"est_shared_r_s_f"} a single
#' value, \code{"est_all"} a free one per cell, and \code{"fix"} holds them all.
#' One spec governs the stream's across-age correlation, its across-year
#' correlation and its unstructured matrix together, so two fleets sharing a
#' correlation share a whole matrix under \code{"us"}. A region, sex or
#' population a fleet never observes carries no parameter whatever the spec says,
#' which is what holds the unused slots of a summed margin out.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#' @param corr \code{"iid"}, \code{"1dar1"}, \code{"us"} or \code{"2dar1"},
#'   either one setting for every fleet or one per fleet.
#' @param stream Stream tag: \code{"catch"}, \code{"discard"},
#'   \code{"fish_idx"} or \code{"srv_idx"}.
#' @param fleet_field \code{"n_fish_fleets"} or \code{"n_srv_fleets"}.
#' @param use_field Name of the use array for this stream.
#' @param starting_values Named list from the caller's \code{...}.
#' @param rho_spec Character string controlling how the correlation parameters
#'   are shared: \code{"est_all"}, \code{"fix"}, or \code{"est_shared_"}
#'   followed by any combination of \code{r}, \code{s} and \code{f}, gaining
#'   \code{p} for the population-specific streams. \code{NULL} (the default)
#'   takes \code{"est_shared_r_s"}, or \code{"est_shared_p_r_s"} when
#'   \code{pop}, both of which give one correlation per fleet.
#' @param pop Logical. \code{TRUE} for the population-specific stream.
#'
#' @return \code{input_list} with the stream's correlation flag and its
#'   correlation parameters set.
#'
#' @keywords internal
do_age_corr_setup <- function(input_list, corr, stream, fleet_field, use_field,
                              starting_values = list(), rho_spec = NULL, pop = FALSE) {

  valid <- c("iid", "1dar1", "us", "2dar1")
  n_fleets <- input_list$data[[fleet_field]]
  n_ages <- length(input_list$data$ages)
  n_sexes <- at_age_n_sexes(input_list)

  tag <- if(pop) paste0(stream, "_pop") else stream
  rho_name <- paste0("trans_rho_", tag)
  yr_name <- paste0("trans_rho_", tag, "_year")
  us_name <- paste0("trans_rho_", tag, "_us")

  if(length(corr) == 1) corr <- rep(corr, n_fleets)
  if(length(corr) != n_fleets) {
    stop("AgeObsCorr_", tag, " has ", length(corr), " entries for ", n_fleets, " fleets. ",
         "Supply one setting per fleet, or one setting for all of them.")
  }
  if(!all(corr %in% valid)) {
    stop("AgeObsCorr_", tag, " is '", paste(unique(corr[!corr %in% valid]), collapse = "', '"),
         "'. Valid options: ", paste(valid, collapse = ", "))
  }

  codes <- convert_to_numeric(corr, list(iid = 0, `1dar1` = 1, us = 2, `2dar1` = 3))
  input_list$data[[paste0("AgeObsCorr_", tag)]] <- codes

  use_arr <- input_list$data[[use_field]]
  nd <- length(dim(use_arr))
  i_a <- nd - 2
  i_y <- nd - 4
  # the correlations sit over the margins the observations are split by, so a
  # slot a fleet never observes drops out on its own
  n_regions <- input_list$data$n_regions
  n_pop <- input_list$data$n_pop
  obs_margins <- if(pop) c(1, 2, nd - 1, nd) else c(1, nd - 1, nd)
  obs_by <- apply(use_arr, obs_margins, function(x) sum(x != 0)) # (pop,) region, sex, fleet

  spec_dims <- if(pop) c(pop = n_pop, region = n_regions, sex = n_sexes, fleet = n_fleets)
               else c(region = n_regions, sex = n_sexes, fleet = n_fleets)
  spec_abbrev <- if(pop) c(p = "pop", r = "region", s = "sex", f = "fleet")
                 else c(r = "region", s = "sex", f = "fleet")
  if(is.null(rho_spec)) rho_spec <- if(pop) "est_shared_p_r_s" else "est_shared_r_s"

  n_pairs <- max(1, n_ages * (n_ages - 1) / 2)
  rho_dims <- as.integer(spec_dims)
  us_dims <- as.integer(c(n_pairs, spec_dims))

  for(nm in c(rho_name, yr_name)) {
    input_list$par[[nm]] <- if(nm %in% names(starting_values)) array(starting_values[[nm]], dim = rho_dims)
                            else array(0, dim = rho_dims)
  } # end nm loop
  input_list$par[[us_name]] <- if(us_name %in% names(starting_values)) array(starting_values[[us_name]], dim = us_dims)
                               else array(0, dim = us_dims)

  # sharing follows the same spec strings the rest of the package uses, so a
  # correlation is not a new idea to learn, only a new place to apply one
  base_map <- array(as.integer(build_shared_spec_map(
    dims = spec_dims, spec = rho_spec, dim_abbrev = spec_abbrev)), dim = rho_dims)
  base_map[obs_by == 0] <- NA_integer_   # a slot this fleet never observes

  renumber <- function(m) {
    if(any(!is.na(m))) m[!is.na(m)] <- as.integer(factor(m[!is.na(m)]))
    return(m)
  }
  hold_unless <- function(m, keep) {
    for(f in seq_len(n_fleets)) {
      if(!codes[f] %in% keep) {
        if(pop) m[,,,f] <- NA_integer_ else m[,,f] <- NA_integer_
      }
    } # end f loop
    return(renumber(m))
  }

  input_list$map[[rho_name]] <- factor(hold_unless(base_map, c(1, 3)))  # across ages
  input_list$map[[yr_name]] <- factor(hold_unless(base_map, 3))        # across years

  # an unstructured correlation shares whole matrices: cells the spec groups
  # together get one matrix between them, each pair its own parameter within it
  us_map <- array(NA_integer_, dim = us_dims)
  grp <- hold_unless(base_map, 2)
  if(any(!is.na(grp))) {
    n_grp <- max(grp, na.rm = TRUE)
    # the pair margin runs fastest, so a group's matrix is n_pairs strided ids
    us_map[] <- rep((seq_len(n_pairs) - 1) * n_grp, times = length(grp)) +
                rep(as.vector(grp), each = n_pairs)
  }
  input_list$map[[us_name]] <- factor(us_map)

  # guard rails: an unstructured correlation grows with the square of the ages,
  # and a separable one needs a complete grid
  if(any(codes %in% c(2, 3))) {
    cell_margins <- setdiff(seq_len(nd), i_a)
    cells <- apply(use_arr == 1, cell_margins, any)
    n_cells <- apply(cells, length(dim(cells)), sum)

    for(f in which(codes == 2)) {
      if(n_ages < 2) stop("AgeObsCorr_", tag, " is 'us' for fleet ", f, ", which needs at ",
                          "least two ages to describe a correlation.")
      if(n_cells[f] <= n_pairs) {
        stop("AgeObsCorr_", tag, " is 'us' for fleet ", f, ", which estimates ", n_pairs,
             " correlations from ", n_cells[f], " observed cells. An unstructured ",
             "correlation across ", n_ages, " ages needs more cells than it has parameters. ",
             "Use '1dar1', which spends one.")
      }
      if(n_cells[f] < 3 * n_pairs) {
        warning("AgeObsCorr_", tag, " is 'us' for fleet ", f, ", estimating ", n_pairs,
                " correlations from ", n_cells[f], " observed cells.")
      }
    } # end f loop

    for(f in which(codes == 3)) {
      if(!at_age_block_complete(use_arr, f, nd, i_y, i_a)) {
        stop("AgeObsCorr_", tag, " is '2dar1' for fleet ", f, ", which correlates over ages ",
             "and years jointly and is defined on a complete grid. That fleet's observed ",
             "ages and years do not form one. Use '1dar1' or 'us', which are defined over ",
             "whatever ages a cell observes.")
      }
    } # end f loop
  }

  return(input_list)
}

#' Do a fleet's at-age observations fill a complete age by year grid?
#'
#' A separable correlation over ages and years is a statement about a rectangle.
#' This reports whether every combination of the years and ages a fleet observes
#' is present, within each population, region, season and sex.
#'
#' @param use_arr Use array for the stream.
#' @param f Fleet index.
#' @param nd Number of dimensions of \code{use_arr}.
#' @param i_y,i_a Positions of the year and age dimensions.
#'
#' @return \code{TRUE} when every block is complete.
#'
#' @keywords internal
at_age_block_complete <- function(use_arr, f, nd, i_y, i_a) {

  d <- dim(use_arr)
  other <- setdiff(seq_len(nd), c(i_y, i_a, nd))
  grid <- if(length(other)) expand.grid(lapply(d[other], seq_len)) else data.frame(x = 1)

  for(g in seq_len(nrow(grid))) {
    idx <- vector("list", nd)
    for(k in seq_along(other)) idx[[other[k]]] <- grid[g,k]
    idx[[i_y]] <- seq_len(d[i_y]); idx[[i_a]] <- seq_len(d[i_a]); idx[[nd]] <- f
    blk <- base::matrix(do.call("[", c(list(use_arr), idx)), nrow = d[i_y])
    obs_yrs <- which(rowSums(blk == 1) > 0)
    obs_ages <- which(colSums(blk == 1) > 0)
    if(length(obs_yrs) == 0) next
    if(!all(blk[obs_yrs, obs_ages] == 1)) return(FALSE)
  } # end g loop

  return(TRUE)
}

#' Is a fleet's selectivity informed by data in some year?
#'
#' Reads both the aggregated use array and its at-age counterpart. A fleet
#' fitting catch at age or an index at age carries no aggregated observations,
#' so keying only off the aggregated array maps its selectivity off and silently
#' holds it at the starting value.
#'
#' @param data Model data list.
#' @param use_field Stub, \code{"Catch"}, \code{"Discard"}, \code{"FishIdx"} or
#'   \code{"SrvIdx"}.
#' @param r,f Region and fleet.
#'
#' @return \code{TRUE} when any stream for that region and fleet is fit.
#'
#' @keywords internal
sel_has_data <- function(data, use_field, r, f) {

  agg <- paste0("Use", use_field)              # aggregated stream
  agg_pop <- paste0(agg, "_pop")
  at_age <- paste0("Use", use_field, "AA")     # at-age stream
  at_age_pop <- paste0(at_age, "_pop")

  any_of <- function(nm, pop) {
    arr <- data[[nm]]
    if(is.null(arr)) return(FALSE)
    if(pop) sum(arr[, r, , , , , f]) > 0 else sum(arr[r, , , , , f]) > 0
  }

  if(!is.null(data[[agg]]) && sum(data[[agg]][r,,,f]) > 0) return(TRUE)
  if(!is.null(data[[agg_pop]]) && sum(data[[agg_pop]][,r,,,f]) > 0) return(TRUE)
  if(any_of(at_age, FALSE)) return(TRUE)
  if(any_of(at_age_pop, TRUE)) return(TRUE)

  return(FALSE)
}

#' Map an at-age observation error or catchability from a key matrix
#'
#' The key is an integer matrix \code{[n_ages, n_fleets]} in which equal entries
#' share a parameter and \code{NA} excludes one. This is the key matrix
#' convention ICES age-structured assessments use for coupling. One structure
#' covers every sharing pattern that would otherwise need its own spec string:
#' one parameter per age, one per age group, or one for the whole fleet.
#'
#' Shared by every at-age stream: the catch, discard and index observation
#' errors, the age-specific catchabilities, and their population-specific
#' counterparts.
#'
#' Coupled parameters are checked against the observations informing them. A
#' standard deviation with a single observation is not merely poorly determined:
#' the likelihood is unbounded, since it collapses onto whatever residual the
#' model can fit exactly and the \code{log(sigma)} term runs to negative
#' infinity. The optimizer reports convergence either way.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#' @param key Integer array \code{[n_ages, n_sexes, n_fleets]}, or \code{NULL}
#'   for the default given by \code{default_shared}. Gains a leading population
#'   dimension when \code{pop} is \code{TRUE}. The sex margin is required: a key
#'   coupling the sexes says so by repeating its entries across them.
#' @param spec \code{"est"} or \code{"fix"}.
#' @param par_name Name of the parameter to map, e.g. \code{"ln_sigmaCAA"}.
#' @param fleet_field \code{"n_fish_fleets"} or \code{"n_srv_fleets"}.
#' @param use_field Name of the use array informing this parameter, e.g.
#'   \code{"UseCatchAA"}. An age a fleet never observes is excluded whatever the
#'   key says.
#' @param starting_values Named list from the caller's \code{...}, read for
#'   \code{par_name}.
#' @param pop Logical. \code{TRUE} for the population-specific stream.
#' @param default_shared Logical. When \code{key} is \code{NULL}, \code{TRUE}
#'   gives one parameter per fleet shared across ages and \code{FALSE} gives one
#'   per age and fleet. Catchability defaults to the latter.
#'
#' @return \code{input_list} with \code{$par$<par_name>} and
#'   \code{$map$<par_name>} set.
#'
#' @keywords internal

do_key_mapping <- function(input_list, key, spec, par_name, fleet_field, use_field,
                           starting_values = list(), pop = FALSE, default_shared = TRUE) {

  n_fleets <- input_list$data[[fleet_field]]
  n_ages <- length(input_list$data$ages)
  n_sexes <- at_age_n_sexes(input_list)
  n_pop <- input_list$data$n_pop
  dims <- if(pop) as.integer(c(n_pop, n_ages, n_sexes, n_fleets)) else as.integer(c(n_ages, n_sexes, n_fleets))
  shape_msg <- paste0(if(pop) "n_pop, " else "", "n_ages, n_sexes, ", fleet_field)

  if(!spec %in% c("est", "fix")) stop(par_name, " spec is '", spec, "', which is not recognized. Valid options: est, fix")

  if(par_name %in% names(starting_values)) input_list$par[[par_name]] <- starting_values[[par_name]]
  else input_list$par[[par_name]] <- array(log(0.5), dim = dims)

  if(!identical(dim(input_list$par[[par_name]]), dims)) {
    stop(par_name, " is not the correct dimension. Should be ", shape_msg)
  }

  if(spec == "fix") {
    input_list$map[[par_name]] <- factor(array(NA_integer_, dim = dims))
    collect_message(par_name, " is specified as: fix")
    return(input_list)
  }

  if(is.null(key)) {
    key <- if(default_shared) array(rep(seq_len(dims[length(dims)]), each = prod(dims[-length(dims)])), dim = dims)
           else array(seq_len(prod(dims)), dim = dims)
  }
  if(!identical(dim(as.array(key)), dims)) {
    stop(par_name, " key is not the correct dimension. Should be ", shape_msg)
  }
  key <- array(as.integer(key), dim = dims)

  # an age and sex a fleet never observes carries no parameter, whatever the key
  # says. This is what holds the unused sexes of a sex-aggregated stream out
  use_arr <- input_list$data[[use_field]]
  nd <- length(dim(use_arr))
  margins <- if(pop) c(1, nd - 2, nd - 1, nd) else c(nd - 2, nd - 1, nd)
  obs_by <- apply(use_arr, margins, function(x) sum(x != 0))
  key[obs_by == 0] <- NA_integer_

  keep <- !is.na(key)
  if(any(keep)) {
    n_obs <- tapply(obs_by[keep], key[keep], sum)
    if(any(n_obs < 2)) {
      stop(par_name, " couples ", length(n_obs), " parameters, and ", sum(n_obs < 2),
           " of them are informed by fewer than 2 observations, so they cannot be ",
           "estimated. Couple more ages or fleets to the same entry, or set that entry ",
           "to NA to hold it fixed.")
    }
    if(any(n_obs < 5)) {
      warning(par_name, " leaves ", sum(n_obs < 5), " of ", length(n_obs),
              " parameters informed by fewer than 5 observations. The smallest has ",
              min(n_obs), ".")
    }
    key[keep] <- as.integer(factor(key[keep]))
  }

  input_list$map[[par_name]] <- factor(key)
  collect_message(par_name, " is specified as: est, with ",
                  length(unique(stats::na.omit(as.vector(key)))), " parameters")
  return(input_list)
}

#' Map an estimated index observation error standard deviation
#'
#' Shared by \code{\link{Setup_Mod_FishIdx_and_Comps}} and
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}}, which carry the same parameter under
#' different names and over a different fleet dimension.
#'
#' The parameter is one value per fleet. Index standard errors are reported per
#' observation, so a year-resolved estimated standard deviation would place one
#' parameter on one observation and drive the likelihood to negative infinity;
#' see \code{\link{check_spec_map_identifiable}} for the same problem in the
#' catch and fishing mortality sigmas, which are dimensioned that way for other
#' reasons.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param spec Character scalar. One of \code{"fix"}, \code{"est_additive"},
#'   \code{"est_quadrature"} or \code{"est_replace"}.
#' @param fleet_field Character. Name of the \code{$data} field giving the
#'   number of fleets, \code{"n_fish_fleets"} or \code{"n_srv_fleets"}.
#' @param par_name Character. Name of the parameter to map, e.g.
#'   \code{"ln_sigmaSrvIdx"}.
#' @param fleet_map Optional integer vector of length \code{n_fleets}. Fleets
#'   sharing a value share a parameter, and \code{NA} fixes that fleet at its
#'   starting value. Defaults to one free parameter per fleet. Useful when a
#'   reference assessment estimated some fleets and pinned others at a bound.
#'
#' @return The input \code{input_list} with \code{$map$<par_name>} set.
#'
#' @keywords internal
do_sigmaIdx_mapping <- function(input_list, spec, fleet_field, par_name, fleet_map = NULL) {

  valid_specs <- c("fix", "est_additive", "est_quadrature", "est_replace")
  if(!spec %in% valid_specs) {
    stop(par_name, " spec '", spec, "' not recognized. Valid options: ",
         paste(valid_specs, collapse = ", "))
  }

  n_fleets <- input_list$data[[fleet_field]]

  if(spec == "fix") {
    input_list$map[[par_name]] <- factor(rep(NA, n_fleets))
    collect_message(par_name, " is specified as: ", spec)
    return(input_list)
  }

  if(is.null(fleet_map)) fleet_map <- seq_len(n_fleets)
  if(length(fleet_map) != n_fleets) {
    stop("The fleet map supplied for ", par_name, " has length ", length(fleet_map),
         " but there are ", n_fleets, " fleets.")
  }

  input_list$map[[par_name]] <- factor(fleet_map)
  collect_message(par_name, " is specified as: ", spec, ", estimating ",
                  length(unique(stats::na.omit(fleet_map))), " of ", n_fleets, " fleets")

  return(input_list)
}

#' Map fishery composition overdispersion (theta) parameters
#'
#' Constructs factor maps for a composition overdispersion parameter (e.g.
#' \code{ln_FishAge_theta}) and its aggregated counterpart (e.g.
#' \code{ln_FishAge_theta_agg}) based on the composition type and likelihood
#' specified in the corresponding \code{$data} fields. Parameters are mapped
#' to \code{NA} for fleets using multinomial likelihoods (\code{LikeType ==
#' 0}), with no observed compositions, or (per-region/per-pop-region) with no
#' active comps for that cell.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param comp_prefix Character, either \code{"FishAge"} or \code{"FishLen"}.
#' @param discard Logical. If \code{TRUE}, maps the discard variant
#'   (\code{..._discard_...} fields). Default \code{FALSE}.
#' @param has_pop Logical. If \code{TRUE}, maps the population-specific
#'   variant (\code{..._pop_...} fields, with an added population dimension
#'   and loop). Default \code{FALSE}.
#' @param fleet_field Character. Name of the \code{$data} field giving the
#'   number of fleets to loop over. Default \code{"n_fish_fleets"}; pass
#'   \code{"n_srv_fleets"} for survey composition types (\code{"SrvAge"},
#'   \code{"SrvLen"}), which have no \code{discard} variant.
#'
#' @return The input \code{input_list} with the corresponding
#'   \code{$map$ln_<...>_theta} and \code{$map$ln_<...>_theta_agg} set to
#'   factor vectors. Active parameters receive sequential integer indices;
#'   inactive parameters are \code{NA}.
#'
#' @keywords internal
do_comp_theta_mapping <- function(input_list, comp_prefix, discard = FALSE, has_pop = FALSE, fleet_field = "n_fish_fleets") {

  suffix <- paste0(if(discard) "_discard" else "", if(has_pop) "_pop" else "")
  data_stub <- paste0(comp_prefix, "Comps", suffix) # e.g. "FishAgeComps_discard_pop"
  par_stub <- paste0(comp_prefix, suffix) # e.g. "FishAge_discard_pop"

  Type_nm <- paste0(data_stub, "_Type")
  LikeType_nm <- paste0(data_stub, "_LikeType")
  Use_nm <- paste0("Use", data_stub)
  theta_nm <- paste0("ln_", par_stub, "_theta")
  theta_agg_nm <- paste0("ln_", par_stub, "_theta_agg")
  n_fleets <- input_list$data[[fleet_field]]

  # setup counters
  counter_agg <- 1
  counter <- 1

  # initialize array to set up mapping
  map_theta <- input_list$par[[theta_nm]]
  map_theta_agg <- input_list$par[[theta_agg_nm]]
  map_theta[] <- NA
  map_theta_agg[] <- NA

  pop_range <- if(has_pop) 1:input_list$data$n_pop else 1

  for(p in pop_range) {
    for(f in 1:n_fleets) {

      # get unique fishery comp types (not itself pop-specific, even for the has_pop variants)
      comp_type <- unique(input_list$data[[Type_nm]][,f])
      like_type <- input_list$data[[LikeType_nm]][f]

      # If aggregated
      if(any(comp_type == 0) && like_type != 0) {
        if(has_pop) map_theta_agg[p,f] <- counter_agg else map_theta_agg[f] <- counter_agg
        counter_agg <- counter_agg + 1 # aggregated
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {

        # a (pop,) region with no active comps for this fleet never contributes to the
        # likelihood, so its theta cell must stay NA (fixed), otherwise it's a free,
        # unidentifiable parameter
        region_has_data <- if(has_pop) {
          sum(input_list$data[[Use_nm]][p,r,,,f]) > 0
        } else {
          sum(input_list$data[[Use_nm]][r,,,f]) > 0
        }

        for(s in 1:input_list$data$n_sexes) {

          # if split by sex and region
          if(any(comp_type == 1) && like_type != 0 && region_has_data) {
            if(has_pop) map_theta[p,r,s,f] <- counter else map_theta[r,s,f] <- counter
            counter <- counter + 1 # split by sex and region
          }

          # joint by sex, split by region
          if(any(comp_type == 2) && like_type != 0 && s == 1 && region_has_data) {
            if(has_pop) map_theta[p,r,1,f] <- counter else map_theta[r,1,f] <- counter
            counter <- counter + 1 # joint by sex, split by region
          }

        } # end s loop
      } # end r loop

      # If we are using a multinomial or there aren't any comps for a given fleet
      use_sum <- if(has_pop) sum(input_list$data[[Use_nm]][p,,,,f]) else sum(input_list$data[[Use_nm]][,,,f])
      if(like_type == 0 || use_sum == 0) {
        if(has_pop) {
          map_theta[p,,,f] <- NA
          map_theta_agg[p,f] <- NA
        } else {
          map_theta[,,f] <- NA
          map_theta_agg[f] <- NA
        }
      }

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map[[theta_nm]] <- factor(map_theta)
  input_list$map[[theta_agg_nm]] <- factor(map_theta_agg)

  return(input_list)
}

#' Map fishery composition correlation parameters
#'
#' Constructs factor maps for a composition correlation-parameter array (e.g.
#' \code{FishAge_corr_pars}, region- and sex-specific AR1/sex correlation) and
#' its aggregated counterpart (e.g. \code{FishAge_corr_pars_agg}) for 1D and
#' 2D logistic-normal composition likelihoods. Parameters are activated only
#' when the corresponding \code{LikeType} is in \code{c(3, 4)} (1D / 2D
#' logistic-normal); all other likelihoods, or fleets with no observed
#' compositions, map correlation parameters to \code{NA}. For the 2D
#' logistic-normal (\code{LikeType == 4}), both trailing elements of the
#' \code{[...,2]} slice are activated: element 1 for the AR1 coefficient and
#' element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' One function serves every composition block. \code{comp_prefix} selects age
#' versus length and fishery versus survey, \code{discard} selects the retained
#' or discarded stream, \code{has_pop} selects the aggregated or population
#' specific stream, and \code{fleet_field} names the fleet count to size the map
#' by. Same parameterization as \code{\link{do_comp_theta_mapping}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param comp_prefix Character, either \code{"FishAge"} or \code{"FishLen"}.
#' @param discard Logical. If \code{TRUE}, maps the discard variant. Default
#'   \code{FALSE}.
#' @param has_pop Logical. If \code{TRUE}, maps the population-specific
#'   variant. Default \code{FALSE}.
#' @param fleet_field Character. Name of the \code{$data} field giving the
#'   number of fleets to loop over. Default \code{"n_fish_fleets"}; pass
#'   \code{"n_srv_fleets"} for survey composition types (\code{"SrvAge"},
#'   \code{"SrvLen"}), which have no \code{discard} variant.
#'
#' @return The input \code{input_list} with the corresponding
#'   \code{$map$<...>_corr_pars} and \code{$map$<...>_corr_pars_agg} set to
#'   factor vectors.
#'
#' @keywords internal
do_comp_corr_pars_mapping <- function(input_list, comp_prefix, discard = FALSE, has_pop = FALSE, fleet_field = "n_fish_fleets") {

  suffix <- paste0(if(discard) "_discard" else "", if(has_pop) "_pop" else "")
  data_stub <- paste0(comp_prefix, "Comps", suffix) # e.g. "FishAgeComps_discard_pop"
  par_stub <- paste0(comp_prefix, suffix) # e.g. "FishAge_discard_pop"

  Type_nm <- paste0(data_stub, "_Type")
  LikeType_nm <- paste0(data_stub, "_LikeType")
  Use_nm <- paste0("Use", data_stub)
  corr_nm <- paste0(par_stub, "_corr_pars")
  corr_agg_nm <- paste0(par_stub, "_corr_pars_agg")
  n_fleets <- input_list$data[[fleet_field]]

  # setup counters
  counter_corr <- 1
  counter_corr_agg <- 1

  # initialize array to set up mapping
  map_corr <- input_list$par[[corr_nm]]
  map_corr_agg <- input_list$par[[corr_agg_nm]]
  map_corr[] <- NA
  map_corr_agg[] <- NA

  pop_range <- if(has_pop) 1:input_list$data$n_pop else 1

  for(p in pop_range) {
    for(f in 1:n_fleets) {

      like_type <- input_list$data[[LikeType_nm]][f]
      use_sum <- if(has_pop) sum(input_list$data[[Use_nm]][p,,,,f]) else sum(input_list$data[[Use_nm]][,,,f])

      # No overdispersion parameters estimated
      if(like_type == 0 || use_sum == 0) {
        if(has_pop) {
          map_corr[p,,,f,] <- NA
          map_corr_agg[p,f] <- NA
        } else {
          map_corr[,,f,] <- NA
          map_corr_agg[f] <- NA
        }
        next # skip if none
      }

      # get unique fishery comp types
      comp_type <- unique(input_list$data[[Type_nm]][,f])

      # Aggregated Correlation Parameters
      if(any(comp_type == 0) && like_type != 0) {
        if(like_type == 3) {
          if(has_pop) map_corr_agg[p,f] <- counter_corr_agg else map_corr_agg[f] <- counter_corr_agg
          counter_corr_agg <- counter_corr_agg + 1 # aggregated
        }
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # Split by region and sex
          if(any(comp_type == 1) && like_type != 0) {
            if(like_type == 3) {
              if(has_pop) map_corr[p,r,s,f,1] <- counter_corr else map_corr[r,s,f,1] <- counter_corr
              counter_corr <- counter_corr + 1
            }
          }

          # Joint by sex, split by region
          if(any(comp_type == 2) && like_type != 0 && s == 1) {

            # 1dar1 correlation
            if(like_type == 3) {
              if(has_pop) map_corr[p,r,1,f,1] <- counter_corr else map_corr[r,1,f,1] <- counter_corr
              counter_corr <- counter_corr + 1
            }

            # 2dar1 correlation
            if(like_type == 4) {
              for(i in 1:2) {
                if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
                if(has_pop) map_corr[p,r,1,f,i] <- counter_corr else map_corr[r,1,f,i] <- counter_corr
                counter_corr <- counter_corr + 1
              } # end i
            } # end if

          }
        } # end s loop
      } # end r loop

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map[[corr_agg_nm]] <- factor(map_corr_agg)
  input_list$map[[corr_nm]] <- factor(map_corr)

  return(input_list)
}

#' Map fishery or survey catchability parameters
#'
#' Constructs the factor map for \code{ln_fish_q} or \code{ln_srv_q},
#' controlling whether catchability parameters are estimated independently
#' per region and time block or shared across regions. Cells with no index
#' observations (aggregated or population-specific) are automatically mapped
#' to \code{NA}. Catchability scales the aggregated index alone, so an at-age
#' index stream does not switch one on: its age multiplier lives in the
#' selectivity. Serves both the fishery and the survey; \code{prefix} is
#' \code{"fish"} or \code{"srv"} and picks which parameter and data names to
#' read.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param q_spec Character vector of length \code{n_fish_fleets}/\code{n_srv_fleets}.
#'   Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate catchability per region × block × fleet.}
#'     \item{\code{"est_shared_r"}}{Single catchability shared across regions,
#'       unique per block × fleet.}
#'     \item{\code{"fix"}}{All catchability parameters fixed (mapped to \code{NA}).}
#'   }
#' @param prefix Character, either \code{"fish"} or \code{"srv"}. Derives
#'   the data/parameter field names: \code{ln_<prefix>_q}, \code{<prefix>_q_blocks},
#'   \code{Use<Fish/Srv>Idx}, \code{Use<Fish/Srv>Idx_pop}.
#' @param fleet_field Character. Name of the \code{$data} field giving the
#'   number of fleets. \code{"n_fish_fleets"} or \code{"n_srv_fleets"}.
#' @param fleet_label Character. Used only in the collected setup message,
#'   e.g. \code{"fishery fleet"} or \code{"survey fleet"}.
#'
#' @return The input \code{input_list} with \code{$map$ln_<prefix>_q} set to a
#'   factor vector.
#'
#' @keywords internal
do_q_mapping <- function(input_list, q_spec, prefix, fleet_field, fleet_label) {

  cap_prefix <- paste0(toupper(substring(prefix, 1, 1)), substring(prefix, 2))
  q_nm <- paste0("ln_", prefix, "_q")
  blocks_nm <- paste0(prefix, "_q_blocks")
  Use_nm <- paste0("Use", cap_prefix, "Idx")
  Use_pop_nm <- paste0(Use_nm, "_pop")
  n_fleets <- input_list$data[[fleet_field]]

  # Initialize counter and mapping array for catchability
  q_counter <- 1
  map_q <- input_list$par[[q_nm]]
  map_q[] <- NA

  for(f in 1:n_fleets) {

    # Validate options
    if(!is.null(q_spec)) {
      if(!q_spec[f] %in% c("est_all", "est_shared_r", "fix"))
        stop(prefix, "_q_spec not correctly specfied. Should be one of these: est_all, est_shared_r, fix")
    }

    for(r in 1:input_list$data$n_regions) {

      if(sum(input_list$data[[Use_nm]][r,,,f]) == 0 && sum(input_list$data[[Use_pop_nm]][,r,,,f]) == 0) {
        map_q[r,,f] <- NA # fix parameters if we are not using indices for these fleets and regions
      } else {

        # Extract number of catchability blocks
        q_blocks_tmp <- unique(as.vector(input_list$data[[blocks_nm]][r,,f]))

        for(b in 1:length(q_blocks_tmp)) {

          # Estimate for all regions
          if(q_spec[f] == 'est_all') {
            map_q[r,b,f] <- q_counter
            q_counter <- q_counter + 1
          } # end if

          # Estimate but share q across regions
          if(q_spec[f] == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              if(q_blocks_tmp[b] %in% input_list$data[[blocks_nm]][rr,,f]) {
                map_q[rr, b, f] <- q_counter
              } # end if
            } # end rr loop
            q_counter <- q_counter + 1
          } # end if

        } # end b loop
      } # end else loop
    } # end r loop

    # fix all parameters
    if(q_spec[f] == 'fix') map_q[,,f] <- NA
    collect_message(prefix, "_q_spec is specified as: ", q_spec[f], " for ", fleet_label, " ", f)
  } # end f loop

  # input into mapping list
  input_list$map[[q_nm]] <- factor(map_q)

  return(input_list)
}

#' Map fixed selectivity parameters (fishery, retention, or survey)
#'
#' Constructs the factor map for fixed-effects selectivity parameters
#' (\code{fish_fixed_sel_pars}, \code{ret_fixed_sel_pars}, or
#' \code{srv_fixed_sel_pars}) across region, bin, block, sex, and fleet,
#' handling every parametric selectivity form (logistic, gamma, exponential,
#' double normal, asymptotic logistic, bicubic spline) as well as
#' non-parametric selectivity with user-defined bin groupings. Fleet sharing
#' (\code{"est_shared_f_x"}) is handled in a second pass, copying the full
#' mapping from a reference fleet.
#'
#' Serves fishery, retention, and survey selectivity. \code{prefix} is
#' \code{"fish"}, \code{"ret"}, or \code{"srv"}; retention takes its own
#' \code{use_field} because it is switched on by the fishery side rather than
#' by a retention specific flag.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param sel_pars_spec Character vector of length \code{n_<fleet_field>}.
#'   Options: \code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"fix"}, \code{"fix_<prefix>_sel_input"}, or
#'   \code{"est_shared_f_x"} (copy from fleet \code{x}).
#' @param bins Number of selectivity bins.
#' @param sel_nonpar_est_bins Optional list \code{[[fleet]][[block]]} of bin
#'   index groupings for non-parametric selectivity (model 5).
#' @param prefix Character, one of \code{"fish"}, \code{"ret"}, or \code{"srv"}.
#'   Drives the domain-specific field names: \code{use_fixed_<prefix>_sel},
#'   \code{<prefix>_sel_blocks}, \code{<prefix>_sel_model},
#'   \code{<prefix>_sel_bicubic_binnodes}/\code{_yrnodes},
#'   \code{<prefix>_fixed_sel_pars} (par/map name).
#' @param fleet_field Character. Name of the \code{$data} field giving the
#'   number of fleets to loop over (\code{"n_fish_fleets"} for both
#'   \code{"fish"} and \code{"ret"}, since retention shares fishery fleets;
#'   \code{"n_srv_fleets"} for \code{"srv"}).
#' @param use_field Character. Stub for the usage-indicator fields read by
#'   \code{\link{sel_has_data}}: \code{Use<use_field>} and its at-age
#'   counterpart \code{Use<use_field>AA}, each with a \code{_pop} variant.
#'   \code{"Catch"} for \code{"fish"}/\code{"ret"}; \code{"SrvIdx"} for
#'   \code{"srv"}.
#' @param fleet_label Character. Used only in the collected setup message,
#'   e.g. \code{"fishery fleet"} or \code{"survey fleet"}.
#'
#' @return The input \code{input_list} with \code{$map$<prefix>_fixed_sel_pars}
#'   set to a factor vector.
#'
#' @keywords internal
#' @importFrom stringr str_detect str_extract_all
do_fixed_sel_pars_mapping <- function(input_list, sel_pars_spec, bins, sel_nonpar_est_bins,
                                       prefix, fleet_field, use_field, fleet_label) {

  par_nm <- paste0(prefix, "_fixed_sel_pars")
  use_fixed_nm <- paste0("use_fixed_", prefix, "_sel")
  blocks_nm <- paste0(prefix, "_sel_blocks")
  model_nm <- paste0(prefix, "_sel_model")
  bicubic_bin_nm <- paste0(prefix, "_sel_bicubic_binnodes")
  bicubic_yr_nm <- paste0(prefix, "_sel_bicubic_yrnodes")
  fix_input_valid <- paste0("fix_", prefix, "_sel_input")
  fix_input_check <- paste0("fixed_", prefix, "_sel_input")
  Use_nm <- paste0("Use", use_field)
  Use_pop_nm <- paste0(Use_nm, "_pop")
  n_fleets <- input_list$data[[fleet_field]]

  # Initialize counter and mapping array for fixed effects selectivity
  sel_pars_counter <- 1
  map_sel_pars <- input_list$par[[par_nm]]
  map_sel_pars[] <- NA

  for(f in 1:n_fleets) {

    # Validate Options
    if(!sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix", "est_shared_s", fix_input_valid) &&
       !stringr::str_detect(sel_pars_spec[f], "est_shared_f_\\d+"))
      stop(prefix, "_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    # checking fixed selex options
    if(input_list$data[[use_fixed_nm]][f] == 1 && stringr::str_detect(sel_pars_spec[f], 'est'))
      stop(use_fixed_nm, " has 1s for a given fleet, but ", prefix, "_fixed_sel_pars_spec is specified at an est variant.")
    if(input_list$data[[use_fixed_nm]][f] == 0 && sel_pars_spec[f] == fix_input_valid)
      stop(use_fixed_nm, " has 0s for a given fleet, but ", prefix, "_fixed_sel_pars_spec is specified at ", fix_input_valid)

    # Skip fleet sharing specs in first pass
    if(stringr::str_detect(sel_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # Only add a counter if catches/index data are avaliable in some years for a given region and fleet combination
      if(sel_has_data(input_list$data, use_field, r, f)) {

        # Extract number of selectivity blocks
        sel_blocks_tmp <- unique(as.vector(input_list$data[[blocks_nm]][r,,f]))

        for(s in 1:input_list$data$n_sexes) {
          for(b in 1:length(sel_blocks_tmp)) {

            block_years <- which(input_list$data[[blocks_nm]][r,,f] == sel_blocks_tmp[b]) # figure out block years
            sel_model_this_block <- unique(input_list$data[[model_nm]][r, block_years, f]) # get selectivity form for a given block
            if(length(sel_model_this_block) > 1) stop("Block ", sel_blocks_tmp[b], " for fleet ", f, " region ", r, " has multiple selectivity models assigned to it")

            # determine maximum selectivity parameters
            if(sel_model_this_block == 2) max_sel_pars <- 1 # exponential
            if(sel_model_this_block %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(sel_model_this_block == 4) max_sel_pars <- 6 # double normal
            if(sel_model_this_block %in% c(6,7)) max_sel_pars <- 3 # logistic w/ asymptotic selectivity
            if(sel_model_this_block == 8) { # bicubic spline: flattened bin-node x year-node grid (group_bins below reduces to a plain 1:max_sel_pars mapping, same as other parametric forms)
              n_bin_nodes_this <- unique(input_list$data[[bicubic_bin_nm]][r, block_years, f])
              n_yr_nodes_this <- unique(input_list$data[[bicubic_yr_nm]][r, block_years, f])
              max_sel_pars <- n_bin_nodes_this * n_yr_nodes_this
            }

            # non-parametric selectivity
            if(sel_model_this_block %in% c(5,9,10)) {

              if(is.null(sel_nonpar_est_bins)) stop("Non-parametric selectivtiy specified, but ", prefix, "_sel_nonpar_est_bins is NULL. Please specify bins!")
              bin_groups <- sel_nonpar_est_bins[[f]][[b]]
              max_sel_pars <- length(bin_groups)  # number of groups = number of estimated pars

              # validate
              all_bins <- unlist(bin_groups)
              if(any(all_bins < 1) || any(all_bins > bins))
                stop(prefix, "_sel_nonpar_est_bins[[", f, "]][[", b, "]] contains indices outside 1:", bins)
              if(length(all_bins) != length(unique(all_bins)))
                stop(prefix, "_sel_nonpar_est_bins[[", f, "]][[", b, "]] has duplicate bin indices")
            }

            for(i in 1:max_sel_pars) {

              # get non-parametric selectivity bins
              group_bins <- if(sel_model_this_block %in% c(5,9,10)) bin_groups[[i]] else i

              # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
              if(sel_pars_spec[f] == "est_all") {
                for(bi in group_bins) map_sel_pars[r,bi,b,s,f] <- sel_pars_counter
                sel_pars_counter <- sel_pars_counter + 1
              } # end if

              # Estimating parameters shared across regions (but unique for each sex, fleet, parameter)
              if(sel_pars_spec[f] == 'est_shared_r' && r == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  if(sel_blocks_tmp[b] %in% input_list$data[[blocks_nm]][rr,,f]) {
                    for(bi in group_bins) map_sel_pars[rr, bi, b, s, f] <- sel_pars_counter
                  } # end if
                } # end rr loop
                sel_pars_counter <- sel_pars_counter + 1
              } # end if

              # Estimating parameters shared across sexes (but unique for each region, fleet, parameter)
              if(sel_pars_spec[f] == 'est_shared_s' && s == 1) {
                for(ss in 1:input_list$data$n_sexes) {
                  for(bi in group_bins) map_sel_pars[r, bi, b, ss, f] <- sel_pars_counter
                } # end ss loop
                sel_pars_counter <- sel_pars_counter + 1
              } # end if

              # Estimating parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  for(ss in 1:input_list$data$n_sexes) {
                    if(sel_blocks_tmp[b] %in% input_list$data[[blocks_nm]][rr,,f]) {
                      for(bi in group_bins) map_sel_pars[rr, bi, b, ss, f] <- sel_pars_counter
                    } # end if
                  } # end ss loop
                } #end rr loop
                sel_pars_counter <- sel_pars_counter + 1
              } # end if

            } # end i loop
          } # end b loop
        } # end s loop
      } # end if statement
    } # end r loop

    # fix all parameters
    if(sel_pars_spec[f] %in% c("fix", fix_input_check)) map_sel_pars[,,,,f] <- NA
    collect_message(prefix, "_fixed_sel_pars_spec is specified as: ", sel_pars_spec[f], " for ", fleet_label, " ", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:n_fleets) {
    if(stringr::str_detect(sel_pars_spec[f], "est_shared_f")) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(sel_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > n_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(sel_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_sel_pars[,,,,f] <- map_sel_pars[,,,,flt_shared]
      collect_message(prefix, "_fixed_sel_pars_spec is specified as: ", sel_pars_spec[f], " for ", fleet_label, " ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map[[par_nm]] <- factor(map_sel_pars)
  return(input_list)
}

#' Map selectivity process error hyperparameters (fishery, retention, or survey)
#'
#' Constructs the factor map for the variance/correlation hyperparameters
#' governing continuous time-varying selectivity (\code{fishsel_pe_pars},
#' \code{retsel_pe_pars}, or \code{srvsel_pe_pars}). The set of active
#' parameters depends on the time-variation type: iid/random-walk forms use
#' up to 2 parameters (log-sigma); 3D GMRF forms use up to 4 (partial
#' correlations for age, year, cohort dimensions plus log-sigma); the 2D AR1
#' form uses 3 (bin AR1, year AR1, log-sigma). Correlation components can be
#' selectively suppressed via \code{corr_opt_semipar}. Fleet sharing
#' (\code{"est_shared_f_x"}) is handled in a second pass.
#'
#' Serves fishery, retention, and survey selectivity, selected by
#' \code{prefix} exactly as in \code{\link{do_fixed_sel_pars_mapping}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param pe_pars_spec Character vector of length \code{n_<fleet_field>}.
#'   Options: \code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"fix"}/\code{"none"}, or
#'   \code{"est_shared_f_x"}.
#' @param corr_opt_semipar Character vector of length \code{n_<fleet_field>}
#'   specifying which correlation components to suppress for semi-parametric
#'   models (\code{NA}, \code{"corr_zero_y"}, \code{"corr_zero_b"},
#'   \code{"corr_zero_y_b"}, \code{"corr_zero_c"}, \code{"corr_zero_y_c"},
#'   \code{"corr_zero_b_c"}, \code{"corr_zero_y_b_c"}). Cohort options are
#'   only valid for 3D GMRF forms.
#' @param bins Number of selectivity bins.
#' @param prefix Character, one of \code{"fish"}, \code{"ret"}, or \code{"srv"}.
#'   Drives the domain-specific field names: \code{cont_tv_<prefix>_sel},
#'   \code{<prefix>_sel_model}, \code{<prefix>_selex_type},
#'   \code{<prefix>sel_pe_pars} (par/map name, no underscore before "sel").
#' @param fleet_field Character. Name of the \code{$data} field giving the
#'   number of fleets (\code{"n_fish_fleets"} for \code{"fish"}/\code{"ret"};
#'   \code{"n_srv_fleets"} for \code{"srv"}).
#' @param use_field Character. Stub for the usage-indicator fields read by
#'   \code{\link{sel_has_data}}: \code{Use<use_field>} and its at-age
#'   counterpart \code{Use<use_field>AA}, each with a \code{_pop} variant.
#'   \code{"Catch"} for \code{"fish"}/\code{"ret"}; \code{"SrvIdx"} for
#'   \code{"srv"}.
#' @param fleet_label Character. Used only in the collected setup message.
#'
#' @return The input \code{input_list} with \code{$map$<prefix>sel_pe_pars}
#'   set to a factor vector.
#'
#' @keywords internal
#' @importFrom stringr str_detect str_extract_all
do_sel_pe_pars_mapping <- function(input_list, pe_pars_spec, corr_opt_semipar, bins,
                                    prefix, fleet_field, use_field, fleet_label) {

  par_nm <- paste0(prefix, "sel_pe_pars")
  cont_tv_nm <- paste0("cont_tv_", prefix, "_sel")
  model_nm <- paste0(prefix, "_sel_model")
  selex_type_nm <- paste0(prefix, "_selex_type")
  n_fleets <- input_list$data[[fleet_field]]

  # Initialize counter and mapping array for process errors
  pe_pars_counter <- 1 # initalize counter
  map_pe_pars <- input_list$par[[par_nm]] # initalize array
  map_pe_pars[] <- NA

  for(f in 1:n_fleets) {

    # Validate options
    if(!is.null(pe_pars_spec)) {
      if(!pe_pars_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s") &&
         !stringr::str_detect(pe_pars_spec[f], "est_shared_f_\\d+"))
        stop(prefix, "sel_pe_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    }

    # Skip fleet sharing specs in first pass
    if(!is.null(pe_pars_spec)) if(stringr::str_detect(pe_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data[[cont_tv_nm]][r,f] == 0 || !sel_has_data(input_list$data, use_field, r, f)) {
        map_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation

        # Figure out max number of selectivity parameters for a given region and fleet
        if(unique(input_list$data[[model_nm]][r,,f]) %in% 2) max_sel_pars <- 1 # exponential
        if(unique(input_list$data[[model_nm]][r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
        if(unique(input_list$data[[model_nm]][r,,f]) == 4) max_sel_pars <- 6 # double normal
        if(unique(input_list$data[[model_nm]][r,,f]) %in% c(5,9,10)) max_sel_pars <- bins # non-parametric selectivity
        if(unique(input_list$data[[model_nm]][r,,f]) %in% c(6,7)) max_sel_pars <- 3 # logistic selectivity w/ asmyptote

        for(s in 1:input_list$data$n_sexes) {

          # If iid time-variation or random walk for this fleet
          if(input_list$data[[cont_tv_nm]][r,f] %in% c(1,2)) {

            for(i in 1:max_sel_pars) {

              # either fixing parameters or not used for a given fleet
              if(pe_pars_spec[f] %in% c("none", "fix")) map_pe_pars[r,i,s,f] <- NA

              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(pe_pars_spec[f] == "est_all") {
                map_pe_pars[r,i,s,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_pe_pars[,i,s,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_pe_pars[r,i,,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_pe_pars[,i,,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              }

            } # end i loop
          } # end iid or random walk variation

          # If 3d gmrf or 2dar1
          if(input_list$data[[cont_tv_nm]][r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data[[cont_tv_nm]][r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data[[cont_tv_nm]][r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_bin, 2 = pcorr_year, 4 = log_sigma)
            if(input_list$data[[cont_tv_nm]][r,f] %in% c(3,4) && input_list$data[[selex_type_nm]] == 1) stop("Cohort-based selectivity deviations are specified, but selectivity is specified as length-based. Please choose another deviation form!")

            for(i in idx) {

              # either fixing parameters or not used for a given fleet
              if(pe_pars_spec[f] %in% c("none", "fix")) map_pe_pars[r,i,s,f] <- NA

              # Estimating all process error parameters
              if(pe_pars_spec[f] == "est_all") {
                map_pe_pars[r,i,s,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_pe_pars[,i,s,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_pe_pars[r,i,,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_pe_pars[,i,,f] <- pe_pars_counter
                pe_pars_counter <- pe_pars_counter + 1
              }

            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(corr_opt_semipar)) {

              opt <- input_list$data[[cont_tv_nm]][r,f] # get random effects options

              # Validate options
              if(!corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_b", "corr_zero_y_b", "corr_zero_c", "corr_zero_y_c", "corr_zero_b_c", "corr_zero_y_b_c"))
                stop("corr_opt_semipar not correctly specfied. Should be one of these: corr_zero_y, corr_zero_b, corr_zero_y_b, corr_zero_c, corr_zero_y_c, corr_zero_b_c, corr_zero_y_b_c, NA")
              if(opt == 5 && corr_opt_semipar[f] %in% c("corr_zero_c","corr_zero_y_c","corr_zero_b_c","corr_zero_y_b_c"))
                stop("Invalid corr_opt_semipar for 2dar1 (opt=5): cohort correlations are not allowed.")

              if (opt %in% c(3,4,5)) {
                # 2d and 3d options
                if (corr_opt_semipar[f] == "corr_zero_y")    map_pe_pars[,2,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_b")    map_pe_pars[,1,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b")  map_pe_pars[,1:2,,f]   <- NA
              }

              if(opt %in% c(3,4)) {
                # 3d gmrf options only (adds the cohort dimension)
                if (corr_opt_semipar[f] == "corr_zero_c")      map_pe_pars[,3,,f]   <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_c")    map_pe_pars[,2:3,,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_b_c")    map_pe_pars[,c(1,3),,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b_c")  map_pe_pars[,1:3,,f] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_pe_pars))
              map_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("corr_opt_semipar is specified as: ", corr_opt_semipar[f], "for ", fleet_label, f)

            }
          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(pe_pars_spec[f] == "fix") map_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(pe_pars_spec)) collect_message(prefix, "sel_pe_pars_spec is specified as: ", pe_pars_spec[f], "for ", fleet_label, f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:n_fleets) {
    if(stringr::str_detect(pe_pars_spec[f], "est_shared_f") && !is.null(pe_pars_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(pe_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > n_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(pe_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_pe_pars[,,,f] <- map_pe_pars[,,,flt_shared]
      collect_message(prefix, "sel_pe_pars_spec is specified as: ", pe_pars_spec[f], " for ", fleet_label, " ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map[[par_nm]] <- factor(map_pe_pars)

  return(input_list)
}

#' Map selectivity deviation parameters (fishery, retention, or survey)
#'
#' Constructs the factor map for continuous time-varying selectivity
#' deviations (\code{ln_fishsel_devs}, \code{ln_retsel_devs}, or
#' \code{ln_srvsel_devs}) across region, year, bin, sex, and fleet. For
#' iid/random-walk forms, active bins are governed by the fitted selectivity
#' model's parameter count; for 3D GMRF/2D AR1 forms, every age bin is active,
#' optionally shared via \code{sel_devs_shared_bins} groupings
#' (\code{"est_shared_b"} and its combinations). Fleet sharing
#' (\code{"est_shared_f_x"}) is handled in a second pass.
#'
#' Serves fishery, retention, and survey selectivity, selected by
#' \code{prefix} exactly as in \code{\link{do_fixed_sel_pars_mapping}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param sel_devs_spec Character vector of length \code{n_<fleet_field>}.
#'   Options: \code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"est_shared_b"}, \code{"est_shared_r_b"},
#'   \code{"est_shared_b_s"}, \code{"est_shared_r_b_s"},
#'   \code{"fix"}/\code{"none"}, or \code{"est_shared_f_x"}.
#' @param sel_devs_shared_bins List of integer vectors, each defining a group
#'   of bins that share a single estimated deviation. Required when
#'   \code{sel_devs_spec} includes \code{"est_shared_b"} or its variants.
#' @param bins Number of selectivity bins.
#' @param prefix Character, one of \code{"fish"}, \code{"ret"}, or \code{"srv"}.
#'   Drives the domain-specific field names: \code{cont_tv_<prefix>_sel},
#'   \code{<prefix>_sel_model}, \code{ln_<prefix>sel_devs} (par/map name, no
#'   underscore before "sel").
#' @param fleet_field Character. Name of the \code{$data} field giving the
#'   number of fleets (\code{"n_fish_fleets"} for \code{"fish"}/\code{"ret"};
#'   \code{"n_srv_fleets"} for \code{"srv"}).
#' @param use_field Character. Stub for the usage-indicator fields read by
#'   \code{\link{sel_has_data}}: \code{Use<use_field>} and its at-age
#'   counterpart \code{Use<use_field>AA}, each with a \code{_pop} variant.
#'   \code{"Catch"} for \code{"fish"}/\code{"ret"}; \code{"SrvIdx"} for
#'   \code{"srv"}.
#' @param fleet_label Character. Used only in the collected setup message.
#'
#' @return The input \code{input_list} with \code{$map$ln_<prefix>sel_devs}
#'   set to a factor vector, and \code{$data$map_ln_<prefix>sel_devs} set to
#'   the equivalent integer array.
#'
#' @keywords internal
#' @importFrom stringr str_detect str_extract_all
do_sel_devs_mapping <- function(input_list, sel_devs_spec, sel_devs_shared_bins, bins,
                                 prefix, fleet_field, use_field, fleet_label) {

  par_nm <- paste0("ln_", prefix, "sel_devs")
  cont_tv_nm <- paste0("cont_tv_", prefix, "_sel")
  model_nm <- paste0(prefix, "_sel_model")
  n_fleets <- input_list$data[[fleet_field]]

  # Initialize counter and mapping array for selectivity deviations
  sel_devs_counter <- 1
  map_sel_devs <- input_list$par[[par_nm]]
  map_sel_devs[] <- NA

  # How many deviation slots a form reads under iid/rw time variation: its
  # parameter count. Slots beyond it are never read by the model and stay
  # unmapped, so a parametric fleet carries no dead deviation parameters.
  sel_dev_slot_count <- function(code) {
    if(code == 2) return(1)                 # exponential
    if(code %in% c(0, 1, 3)) return(2)      # logistic or gamma
    if(code == 4) return(6)                 # double normal
    if(code %in% c(5, 9, 10)) return(bins)      # non-parametric
    if(code %in% c(6, 7)) return(3)         # logistic w/ asymptote
    stop("Continuous time-varying (iid/rw) deviations are not supported for selectivity model code ", code, ".")
  }

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:n_fleets) {

      # Validate options
      if(!is.null(sel_devs_spec)) {
        if(!sel_devs_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "est_shared_b", "est_shared_r_b", "est_shared_r_b_s", "est_shared_b_s") &&
           !stringr::str_detect(sel_devs_spec[f], "est_shared_f_\\d+"))
          stop(prefix, "_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, est_shared_b, est_shared_r_b, est_shared_r_b_s, est_shared_r_s, fix, or est_shared_f_# (where # is fleet number)")
        # Sharing across bins needs the deviations to be indexed by bin. That is
        # true for the GMRF and AR1 time-varying forms whatever the functional
        # form, and also for the non-parametric forms under iid or a random walk,
        # where each deviation already belongs to one bin. It is not true for a
        # parametric form under iid or a random walk, whose deviation slots are
        # its parameters rather than its bins.
        nonpar_fleet <- all(input_list$data[[model_nm]][r,,f] %in% c(5,9,10))
        if(sel_devs_spec[f] %in% c("est_shared_b", "est_shared_r_b", "est_shared_r_b_s", "est_shared_b_s") &&
           !input_list$data[[cont_tv_nm]][r,f] %in% c(3,4,5) && !nonpar_fleet)
          stop("Sharing bin deviations is only supported when the deviations are indexed by bin: either a GMRF or AR1 time-varying form, or a non-parametric selectivity form. A parametric form under iid or a random walk indexes its deviations by parameter instead.")
       }

      # Skip fleet sharing specs in first pass
      if(!is.null(sel_devs_spec)) if(stringr::str_detect(sel_devs_spec[f], "est_shared_f")) next

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:(length(input_list$data$years) + input_list$data$n_proj_yrs_devs)) {

          # Which regions actually have data for this fleet
          reg_has_dat <- sapply(1:input_list$data$n_regions, function(rr)
            sel_has_data(input_list$data, use_field, rr, f))
          r_anchor <- if(any(reg_has_dat)) min(which(reg_has_dat)) else 1L
          shares_r <- !is.null(sel_devs_spec) &&
            sel_devs_spec[f] %in% c('est_shared_r', 'est_shared_r_s', 'est_shared_r_b', 'est_shared_r_b_s')
          dat_ok <- if(shares_r) any(reg_has_dat) else reg_has_dat[r]

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data[[cont_tv_nm]][r,f] == 0 || !dat_ok) {
            map_sel_devs[r,y,,s,f] <- NA
          } else {

            # Figure out max number of selectivity parameters for a given region and
            # fleet. A fleet changing forms across blocks takes the most any of its
            # forms reads, since the deviation slots are shared across years.
            max_sel_pars <- max(sapply(unique(input_list$data[[model_nm]][r,,f]), sel_dev_slot_count))

            # If iid or random walk time-variation for this fleet
            if(input_list$data[[cont_tv_nm]][r,f] %in% c(1,2)) {

              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(sel_devs_spec[f] == 'est_all') {
                  map_sel_devs[r,y,i,s,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(sel_devs_spec[f] == 'est_shared_r' && r == r_anchor) {
                  map_sel_devs[,y,i,s,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_sel_devs[r,y,i,,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(sel_devs_spec[f] == 'est_shared_r_s' && r == r_anchor && s == 1) {
                  map_sel_devs[,y,i,,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data[[cont_tv_nm]][r,f] %in% c(3,4,5)) {

              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and bins
                if(sel_devs_spec[f] == 'est_all') {
                  map_sel_devs[r,y,i,s,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across regions
                if(sel_devs_spec[f] == 'est_shared_r' && r == r_anchor) {
                  map_sel_devs[,y,i,s,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes
                if(sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_sel_devs[r,y,i,,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes and regions
                if(sel_devs_spec[f] == 'est_shared_r_s' && s == 1 && r == r_anchor) {
                  map_sel_devs[,y,i,,f] <- sel_devs_counter
                  sel_devs_counter <- sel_devs_counter + 1
                }

                if(sel_devs_spec[f] == 'est_shared_b') {
                  for(k in 1:length(sel_devs_shared_bins)) {
                    map_sel_devs[r,y,sel_devs_shared_bins[[k]],s,f] <- sel_devs_counter
                    sel_devs_counter <- sel_devs_counter + 1
                  } # end k loop
                }

                if(sel_devs_spec[f] == 'est_shared_r_b' && r == r_anchor) {
                  for(k in 1:length(sel_devs_shared_bins)) {
                    map_sel_devs[,y,sel_devs_shared_bins[[k]],s,f] <- sel_devs_counter
                    sel_devs_counter <- sel_devs_counter + 1
                  } # end k loop
                }

                if(sel_devs_spec[f] == 'est_shared_b_s' && s == 1) {
                  for(k in 1:length(sel_devs_shared_bins)) {
                    map_sel_devs[r,y,sel_devs_shared_bins[[k]],,f] <- sel_devs_counter
                    sel_devs_counter <- sel_devs_counter + 1
                  } # end k loop
                }

                if(sel_devs_spec[f] == 'est_shared_r_b_s' && s == 1 && r == r_anchor) {
                  for(k in 1:length(sel_devs_shared_bins)) {
                    map_sel_devs[,y,sel_devs_shared_bins[[k]],,f] <- sel_devs_counter
                    sel_devs_counter <- sel_devs_counter + 1
                  } # end k loop
                }

              } # end i loop
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

      if(!is.null(sel_devs_spec)) collect_message(prefix, "_sel_devs_spec is specified as: ", sel_devs_spec[f], "for ", fleet_label, f, "and region ", r)

    } # end f loop
  } # end r loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:n_fleets) {
    if(stringr::str_detect(sel_devs_spec[f], "est_shared_f") && !is.null(sel_devs_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(sel_devs_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > n_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(sel_devs_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_sel_devs[,,,,f] <- map_sel_devs[,,,,flt_shared]
      collect_message(prefix, "_sel_devs_spec is specified as: ", sel_devs_spec[f], " for ", fleet_label, " ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map[[par_nm]] <- factor(map_sel_devs)
  input_list$data[[paste0("map_", par_nm)]] <- array(as.numeric(input_list$map[[par_nm]]), dim = dim(input_list$par[[par_nm]]))

  return(input_list)
}

#' Refresh the map mirrors held in the data list
#'
#' Several deviation penalties key on a copy of their parameter's factor map
#' carried in the data list under \code{map_<par>}, since the map itself is
#' applied by \code{RTMB::MakeADFun} and is invisible inside the objective.
#' Those copies are written at setup, so a map edited by hand afterwards would
#' otherwise leave the penalty evaluating deviations that are no longer
#' estimated. Rebuilding the mirrors from the map immediately before the model
#' is constructed keeps the two in step, with the map treated as authoritative.
#'
#' A mirror whose parameter has no entry in \code{mapping}, or whose length no
#' longer matches (as when a caller has truncated one but not the other), is
#' left untouched.
#'
#' @param data Named list of model data, as passed to \code{RTMB::MakeADFun}.
#' @param mapping Named list of factor maps, as passed to
#'   \code{RTMB::MakeADFun}.
#'
#' @return \code{data} with every \code{map_<par>} element refreshed from
#'   \code{mapping[[par]]}.
#'
#' @keywords internal
sync_dev_map_data <- function(data, mapping) {

  for(nm in grep("^map_", names(data), value = TRUE)) {
    par_nm <- sub("^map_", "", nm)
    if(is.null(mapping[[par_nm]])) next
    if(length(mapping[[par_nm]]) != length(data[[nm]])) next
    data[[nm]] <- array(as.numeric(mapping[[par_nm]]), dim = dim(data[[nm]]))
  }

  return(data)
}
