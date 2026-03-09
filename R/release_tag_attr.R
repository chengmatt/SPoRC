#' Distribute Tagged Fish Releases to Full Population Dimensions
#'
#' When tag release data are not recorded at full population resolution —
#' i.e. when one or more of the population, age, or sex dimensions are
#' unattended in \code{tag_attr} — this function distributes the known tag
#' totals to full \code{[n_pop, n_ages, n_sexes]} resolution using
#' apportionment weights derived from the release platform (population
#' abundance, fishery catch-at-age, or survey index-at-age). If all three
#' dimensions are attended (\code{tag_attr = "p_a_s"}), \code{tagged_fish}
#' is returned unchanged with no computation performed.
#'
#' Apportionment weights are constructed from numbers-at-age
#' (\code{platform = "population"}), numbers-at-age multiplied by fishery
#' selectivity (\code{platform = "fishery"}), or numbers-at-age multiplied
#' by survey selectivity (\code{platform = "survey"}), all evaluated at the
#' release region, year, and season. Weights are then normalised
#' conditionally on the attended dimensions: the denominator for cell
#' \code{[p, a, s]} is the sum of raw weights across all cells that share
#' the same indices in the attended dimensions. This ensures that the
#' marginal totals of \code{tagged_fish} are preserved exactly along every
#' attended dimension. For example, if only age is attended
#' (\code{tag_attr = "a"}), age-specific totals in \code{tagged_fish} are
#' preserved while tags are distributed across population and sex in
#' proportion to the platform weights.
#'
#' @param tagged_fish Numeric vector or array of released tagged fish for a
#'   single tag cohort. Unattended dimensions are expected to be collapsed
#'   to index 1. Reshaped internally to \code{[n_pop, n_ages, n_sexes]}.
#' @param tag_attr Character string specifying which population dimensions
#'   are attended in \code{tagged_fish}. Constructed from any combination of
#'   \code{"p"} (population), \code{"a"} (age), and \code{"s"} (sex),
#'   joined by underscores (e.g. \code{"p_a_s"}, \code{"a"}, \code{"p_a"}).
#'   See \code{\link{Setup_Mod_Tagging}} for the full set of valid strings.
#' @param tag_release_platform Character vector of length 2. Element 1 is
#'   the release platform: one of \code{"population"}, \code{"fishery"}, or
#'   \code{"survey"}. Element 2 is the fleet index as a character string
#'   (coerced to integer internally), or \code{NA} when
#'   \code{platform = "population"}.
#' @param srv_sel Numeric array of survey selectivity
#'   \code{[n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets]}. Used as the
#'   age-sex apportionment weight when \code{platform = "survey"}.
#' @param fish_sel Numeric array of fishery selectivity
#'   \code{[n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets]}. Used as the
#'   age-sex apportionment weight when \code{platform = "fishery"}.
#' @param NAA Numeric array of numbers-at-age prior to movement
#'   \code{[n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes]}. Used
#'   directly as weights when \code{platform = "population"}, and multiplied
#'   by selectivity for fishery and survey platforms.
#' @param ty Integer. Model year index of the tag release cohort.
#' @param tseas Integer. Season index of the tag release cohort.
#' @param tr Integer. Region index of the tag release cohort.
#' @param n_pop Integer. Number of populations.
#' @param n_ages Integer. Number of age classes.
#' @param n_sexes Integer. Number of sexes (1 or 2).
#'
#'
#' @keywords internal
release_conv_tag_attr <- function(tagged_fish,
                                  tag_attr,
                                  tag_release_platform,
                                  srv_sel,
                                  fish_sel,
                                  NAA,
                                  ty,
                                  tseas,
                                  tr,
                                  n_pop,
                                  n_ages,
                                  n_sexes
                                  ) {

  "c" = RTMB::ADoverload("c")
  "[<-" = RTMB::ADoverload("[<-")

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
    NAA_slice  =NAA[, tr, ty, tseas, , , drop = FALSE]
    dim(NAA_slice) =c(n_pop, n_ages, n_sexes)
    sel_slice  =fish_sel[tr, ty, , , fleet, drop = FALSE]
    dim(sel_slice) =c(n_ages, n_sexes)
    sweep(NAA_slice, 2:3, sel_slice, "*")
  } else if(platform == "survey") {
    NAA_slice  =NAA[, tr, ty, tseas, , , drop = FALSE]
    dim(NAA_slice) =c(n_pop, n_ages, n_sexes)
    sel_slice  =srv_sel[tr, ty, , , fleet, drop = FALSE]
    dim(sel_slice) =c(n_ages, n_sexes)
    sweep(NAA_slice, 2:3, sel_slice, "*")
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
