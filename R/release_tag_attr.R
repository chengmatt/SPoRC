#' Distribute Tagged Fish Releases to Full Population Dimensions
#'
#' When tag release data are not recorded at full population resolution
#' (i.e., when one or more of population, age, or sex dimensions are not
#' attended in \code{tag_attr}), this function distributes the known tag
#' totals out to full \code{[n_pop, n_ages, n_sexes]} resolution using
#' apportionment weights derived from the release platform (population
#' abundance, catch, or survey index). If all three dimensions are attended
#' (\code{tag_attr = "p_a_s"}), the function returns \code{tagged_fish}
#' unchanged.
#'
#' Apportionment is performed conditionally on the attended dimensions. For
#' each cell \code{[p, a, s]}, the normalized weight is the raw weight
#' divided by the sum of raw weights over all cells that share the same
#' attended-dimension indices. For example, if only age is attended, the
#' denominator is the total weight at that age across all populations and
#' sexes, so that the age totals in \code{tagged_fish} are preserved while
#' tags are distributed across the unattended population and sex dimensions.
#'
#' @param tagged_fish Numeric vector or array of released tagged fish for a
#'   single tag cohort, with values placed into index 1 of any unattended
#'   dimensions. Will be reshaped internally to
#'   \code{[n_pop, n_ages, n_sexes]}.
#' @param tag_attr Character string specifying which population dimensions
#'   are attended. Built from any combination of \code{"p"} (population),
#'   \code{"a"} (age), and \code{"s"} (sex), joined by underscores (e.g.,
#'   \code{"p_a_s"}, \code{"a"}, \code{"p_a"}). See
#'   \code{\link{Setup_Mod_Tagging}} for full details.
#' @param tag_release_platform Character vector of length 2 giving the
#'   release platform (\code{"population"}, \code{"fishery"}, or
#'   \code{"survey"}) and the fleet index (or \code{NA} for
#'   \code{"population"}) for this tag cohort.
#' @param SrvIAA Numeric array of survey index-at-age used as apportionment
#'   weights when \code{platform = "survey"}.
#' @param CAA Numeric array of catch-at-age used as apportionment weights
#'   when \code{platform = "fishery"}.
#' @param NAA Numeric array of numbers-at-age used as apportionment weights
#'   when \code{platform = "population"}.
#' @param ty Integer. Current model year index.
#' @param tseas Integer. Current season index.
#' @param tr Integer. Release region index.
#' @param n_pop Integer. Number of populations.
#' @param n_ages Integer. Number of age classes.
#' @param n_sexes Integer. Number of sexes.
#'
#' @return Numeric array of dimensions \code{[n_pop, n_ages, n_sexes]} with
#'   tagged fish distributed to full population resolution according to the
#'   apportionment weights.
#'
#' @keywords internal
release_conv_tag_attr <- function(tagged_fish,
                                  tag_attr,
                                  tag_release_platform,
                                  SrvIAA,
                                  CAA,
                                  NAA,
                                  ty,
                                  tseas,
                                  tr,
                                  n_pop,
                                  n_ages,
                                  n_sexes
                                  ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # reshape into array
  tagged_fish = array(tagged_fish, dim = c(n_pop, n_ages, n_sexes))

  # get tagging attributes
  attr_parts = strsplit(tag_attr, "_")[[1]]
  attended_p = "p" %in% attr_parts
  attended_a = "a" %in% attr_parts
  attended_s = "s" %in% attr_parts

  # return tagged_fish if p_a_s
  if(attended_p && attended_a && attended_s) return(tagged_fish)

  # get platform and fleet release
  platform = tag_release_platform[1]
  fleet = as.integer(tag_release_platform[2])

  # Get raw tag apportionment weights
  weights = if(platform == "population") {
    array(NAA[, tr, ty, tseas, , ], dim = c(n_pop, n_ages, n_sexes))
  } else if(platform == "fishery") {
    array(CAA[, tr, ty, tseas, , , fleet], dim = c(n_pop, n_ages, n_sexes))
  } else if(platform == "survey") {
    array(SrvIAA[, tr, ty, tseas, , , fleet], dim = c(n_pop, n_ages, n_sexes))
  }

  # Normalize weights within attended dimensions.
  # The denominator sums over unattended dims only, preserving totals
  # within each attended dim combination.
  norm_weights = array(0, dim = c(n_pop, n_ages, n_sexes))
  for(p in seq_len(n_pop)) {
    for(a in seq_len(n_ages)) {
      for(s in seq_len(n_sexes)) {
        # Denominator: sum weights over the unattended dims that share the
        # same attended-dim indices as this [p, a, s] cell
        denom = 0
        for(pp in seq_len(n_pop)) {
          for(aa in seq_len(n_ages)) {
            for(ss in seq_len(n_sexes)) {
              # Include this cell in the denominator only if it matches on
              # all attended dimensions
              same_p = !attended_p || pp == p
              same_a = !attended_a || aa == a
              same_s = !attended_s || ss == s
              if(same_p && same_a && same_s) denom = denom + weights[pp, aa, ss]
            } # end ss loop
          } # end aa loop
        } # end pp loop
        norm_weights[p, a, s] = weights[p, a, s] / denom
      } # end s loop
    } # end a loop
  } # end p loop

  out = array(0, dim = c(n_pop, n_ages, n_sexes))
  for(p in seq_len(n_pop)) {
    for(a in seq_len(n_ages)) {
      for(s in seq_len(n_sexes)) {
        pi = if(attended_p) p else 1
        ai = if(attended_a) a else 1
        si = if(attended_s) s else 1
        out[p, a, s] = tagged_fish[pi, ai, si] * norm_weights[p, a, s]
      } # end s loop
    } # end a loop
  } # end p loop

  return(out)

} # end function
