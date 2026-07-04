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
#'
#' @return Factor vector of length \code{prod(dims)}, suitable for direct
#'   assignment to \code{input_list$map$<par>}.
#'
#' @keywords internal
build_shared_spec_map <- function(dims, spec, dim_abbrev) {

  shared_specs <- unlist(lapply(seq_along(dim_abbrev), function(k) {
    combs <- combn(names(dim_abbrev), k)
    apply(combs, 2, function(x) paste0("est_shared_", paste(x, collapse = "_")))
  }))
  valid_specs <- c("fix", "est_all", shared_specs)

  if(!spec %in% valid_specs) stop("spec '", spec, "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  if(spec == "fix") return(factor(rep(NA, prod(dims))))

  share_over <- if(spec == "est_all") character(0) else {
    parts <- strsplit(sub("^est_shared_", "", spec), "_")[[1]]
    unname(dim_abbrev[parts])
  }

  factor(as.vector(build_pe_map(dims, share_over = share_over)))
}
