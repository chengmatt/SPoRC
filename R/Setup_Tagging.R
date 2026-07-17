#' Set up conventional tagging dynamics for the operating model simulation
#'
#' Populates \code{sim_list} with all conventional tagging inputs needed by
#' the operating model: tag release cohort definitions, release platform,
#' tagging timing, tag-induced mortality, chronic shedding, recapture
#' likelihood settings, fishery reporting rates, and zero-filled containers
#' for predicted and observed recaptures. Must be called after
#' \code{\link{Setup_Sim_Dim}}.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}}.
#' @param n_tags Numeric scalar. Constant number of tags released per release
#'   event. Cannot be specified simultaneously with \code{n_tags_rel_input}.
#'   Default \code{NULL}.
#' @param n_tags_rel_input Array specifying cohort-specific tag release
#'   numbers (e.g., varying by length, age, or other structure). Overrides
#'   \code{n_tags} when supplied. Cannot be specified simultaneously with
#'   \code{n_tags}. Default \code{NULL}.
#' @param use_conv_fish_tagging Integer (0/1). Whether conventional tag
#'   recaptures from fisheries are simulated. Default \code{0}.
#' @param conv_tag_max_liberty Integer. Maximum number of years-at-liberty
#'   tracked per tag cohort. Controls the size of recapture containers.
#'   Default: \code{n_ages / 2}.
#' @param conv_tag_release_indicator Data frame defining tag release cohorts.
#'   Required columns: \code{regions} (release region), \code{tag_years}
#'   (release year), \code{tag_seas} (release season). Default: all
#'   combinations of regions, years, and seasons from \code{sim_list}.
#' @param conv_tag_release_platform Character matrix
#'   \code{[n_conv_tag_cohorts × 2]} specifying the release platform and
#'   associated fleet index for each cohort. Rows must align with
#'   \code{conv_tag_release_indicator}. Column \code{platform} must be one
#'   of:
#'   \describe{
#'     \item{\code{"population"}}{Tags released directly into the population,
#'       independent of any fleet or survey sampling process.
#'       \code{fleet} column should be \code{NA}.}
#'     \item{\code{"fishery"}}{Tags released through a fishery fleet.
#'       \code{fleet} column gives the fleet index.}
#'     \item{\code{"survey"}}{Tags released through a survey fleet.
#'       \code{fleet} column gives the fleet index.}
#'   }
#'   Default: \code{"survey"} platform with fleet \code{1} for every cohort.
#' @param conv_tag_t_tagging Numeric scalar or vector of length
#'   \code{n_tag_rel_events} (one value per row of
#'   \code{conv_tag_release_indicator}), each in \eqn{[0, 1]}. Fraction of the
#'   season remaining at the time of tag release for that release event.
#'   \code{1} = start of season; \code{0.5} = mid-season; \code{0} = end of
#'   season. A scalar is recycled to all release events. Default \code{1}.
#' @param ln_init_conv_tag_mort Log-scale initial tag-induced mortality
#'   applied at the moment of release. Numeric scalar or vector of length
#'   \code{n_tag_rel_events}, one value per release event. A scalar is
#'   recycled to all release events. Default \code{-1000} (approximately
#'   zero mortality).
#' @param ln_conv_tag_shed Log-scale annual chronic tag shedding rate.
#'   Numeric scalar or vector of length \code{n_tag_rel_events}, one value per
#'   release event. A scalar is recycled to all release events. Default
#'   \code{-1000} (approximately no shedding).
#' @param conv_fish_tag_like Integer or character scalar specifying the tag
#'   recapture likelihood. Default \code{0} (Poisson). Options:
#'   \describe{
#'     \item{\code{0}/\code{"Poisson"}}{Poisson.}
#'     \item{\code{1}/\code{"NegBin"}}{Negative binomial.}
#'     \item{\code{2}/\code{"Multinomial_Release"}}{Multinomial conditioned
#'       on releases.}
#'     \item{\code{3}/\code{"Multinomial_Recapture"}}{Multinomial conditioned
#'       on recaptures.}
#'     \item{\code{4}/\code{"Dirichlet-Multinomial_Release"}}{Dirichlet-multinomial
#'       conditioned on releases.}
#'     \item{\code{5}/\code{"Dirichlet-Multinomial_Recapture"}}{Dirichlet-multinomial
#'       conditioned on recaptures.}
#'   }
#' @param ln_conv_fish_tag_theta Log-scale overdispersion parameter for
#'   negative-binomial (\code{1}) or Dirichlet-multinomial (\code{4},
#'   \code{5}) likelihoods. Ignored for Poisson and multinomial likelihoods.
#'   Default \code{log(1)}.
#' @param conv_fish_tag_attr Character scalar or character vector of length
#'   \code{n_tag_rel_events} specifying which biological
#'   dimensions are resolved at release (and hence retained at recapture) for
#'   each tag release event. A scalar is recycled to every event; a vector
#'   lets different events resolve different dimensions (e.g. event 1 known at
#'   \code{"p_a_s"}, event 2 only at \code{"p_a"}). Each element is built from
#'   any combination of \code{"p"} (population), \code{"a"} (age), and
#'   \code{"s"} (sex), joined by underscores. Region and fleet are always
#'   retained. Valid values: \code{"p_a_s"}, \code{"a_s"}, \code{"p_a"},
#'   \code{"p_s"}, \code{"a"}, \code{"s"}, \code{"p"}, \code{"none"}. Default
#'   \code{"p_a_s"}.
#' @param conv_tag_fish_reporting_input Fishery tag reporting rate array
#'   \code{[n_regions × n_yrs × n_fish_fleets × n_sims]}. Values represent
#'   the probability that a recaptured tag is reported and must be in
#'   \eqn{[0, 1]}. Default: \code{0.5} for all cells.
#'
#' @return The input \code{sim_list} with tagging-related fields appended:
#'   \code{$n_tags} or \code{$n_tags_rel_input} (depending on which is
#'   provided), \code{$conv_tag_max_liberty}, \code{$conv_tag_t_tagging}
#'   (length \code{n_tag_rel_events}), \code{$ln_init_conv_tag_mort} (length
#'   \code{n_tag_rel_events}), \code{$ln_conv_tag_shed} (length
#'   \code{n_tag_rel_events}),
#'   \code{$conv_tag_release_indicator}, \code{$conv_tag_release_platform},
#'   \code{$n_tag_rel_events}, \code{$use_conv_fish_tagging},
#'   \code{$conv_fish_tag_like}, \code{$conv_fish_tag_attr},
#'   \code{$ln_conv_fish_tag_theta}, \code{$conv_tag_fish_reporting}, and
#'   zero-filled containers \code{$conv_tagged_fish},
#'   \code{$conv_tagged_fish_attr}, \code{$conv_tag_fish_avail},
#'   \code{$obs_conv_tag_fish_recap}, \code{$pred_conv_tag_fish_recap}.
#'
#'
#' @export Setup_Sim_Tagging
#' @family Simulation Setup
Setup_Sim_Tagging <- function(
                              sim_list,
                              n_tags = NULL,
                              n_tags_rel_input = NULL,
                              use_conv_fish_tagging = 0,
                              conv_tag_max_liberty = sim_list$n_ages / 2,
                              conv_tag_release_indicator = expand.grid(regions = 1:sim_list$n_regions, tag_years = 1:sim_list$n_yrs, tag_seas = 1:sim_list$n_seas),
                              conv_tag_release_platform = matrix(c("survey", "1"), nrow = nrow(conv_tag_release_indicator),  ncol = 2, byrow = TRUE,dimnames = list(NULL, c("platform", "fleet"))),
                              conv_tag_t_tagging = 1,
                              ln_init_conv_tag_mort = -1000,
                              ln_conv_tag_shed = -1000,
                              conv_fish_tag_attr = 'p_a_s',
                              conv_tag_fish_reporting_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              conv_fish_tag_like = 0,
                              ln_conv_fish_tag_theta = log(1)
                              ) {

  if(any(use_conv_fish_tagging == 1)) {
    check_sim_dimensions(conv_tag_fish_reporting_input, n_regions = sim_list$n_regions,
                         n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets,
                         n_sims = sim_list$n_sims, what = "conv_tag_fish_reporting_input")
  }

  # Convert codes to numeric
  conv_fish_tag_like <- convert_to_numeric(conv_fish_tag_like, list(Poisson = 0, NegBin = 1,
                                                Multinomial_Release = 2,
                                                Multinomial_Recapture = 3,
                                                `Dirichlet-Multinomial_Release` = 4,
                                                `Dirichlet-Multinomial_Recapture` = 5))

  # Output variables into list
  if(!is.null(n_tags) && !is.null(n_tags_rel_input)) stop("n_tags and n_tags_rel_input cannot be specified simultaneously. n_tags is a scalar, while n_tags_rel_input specifies cohort-specific tags!")
  if(!is.null(n_tags)) sim_list$n_tags <- n_tags
  if(!is.null(n_tags_rel_input)) sim_list$n_tags_rel_input <- n_tags_rel_input
  sim_list$conv_tag_max_liberty <- conv_tag_max_liberty
  sim_list$conv_tag_release_indicator <- conv_tag_release_indicator # tag release indicator (by tag years and regions = a tag cohort)
  sim_list$conv_tag_release_platform <- conv_tag_release_platform # how tags are released
  # number of tag release events - tag years x tag region (tag cohorts).
  # When tagging is inactive, conv_tag_release_indicator may arrive as NA/NULL
  # (no rows); fall back to a length-1 placeholder so scalar recycling is valid.
  sim_list$n_tag_rel_events <- if(is.null(nrow(conv_tag_release_indicator))) 1L else nrow(conv_tag_release_indicator)

  # Per-release-event timing / mortality / shedding (scalars are recycled to all events)
  sim_list$conv_tag_t_tagging <- recycle_tag_event_par(conv_tag_t_tagging, sim_list$n_tag_rel_events, "conv_tag_t_tagging") # time of tagging
  sim_list$ln_init_conv_tag_mort <- recycle_tag_event_par(ln_init_conv_tag_mort, sim_list$n_tag_rel_events, "ln_init_conv_tag_mort") # tag induced mortality
  sim_list$ln_conv_tag_shed <- recycle_tag_event_par(ln_conv_tag_shed, sim_list$n_tag_rel_events, "ln_conv_tag_shed") # tag shedding

  # Containers
  sim_list$conv_tagged_fish <- sim_list$conv_tagged_fish_attr <-
    array(0, dim = c(sim_list$n_tag_rel_events, sim_list$n_pop,
                     sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # number of tagged fish

  sim_list$conv_tag_fish_avail <-
    array(0, dim = c(sim_list$conv_tag_max_liberty + 1, sim_list$n_seas, sim_list$n_tag_rel_events,
                     sim_list$n_pop, sim_list$n_regions,
                     sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # tags availiable for recapture every year

  sim_list$obs_conv_tag_fish_recap <- sim_list$pred_conv_tag_fish_recap <-
    array(0, dim = c(sim_list$conv_tag_max_liberty, sim_list$n_seas, sim_list$n_tag_rel_events,
                     sim_list$n_pop, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes,
                     sim_list$n_fish_fleets, sim_list$n_sims)) # predicted tag recaptures

  sim_list$conv_tag_fish_reporting <- conv_tag_fish_reporting_input # output this into list
  sim_list$use_conv_fish_tagging <- use_conv_fish_tagging # output into list
  sim_list$conv_fish_tag_like <- conv_fish_tag_like # tag likelihood
  sim_list$conv_fish_tag_attr <- recycle_tag_event_par(conv_fish_tag_attr, sim_list$n_tag_rel_events, "conv_fish_tag_attr") # tag release/recapture attributes for fishery conventional tags (one per release event; scalar recycled)
  sim_list$ln_conv_fish_tag_theta <- ln_conv_fish_tag_theta # tag likelihood overdispersion parameter

  return(sim_list)
}

#' Recycle a per-tag-release-event scalar or vector to the full event count
#'
#' Internal helper shared by \code{\link{Setup_Sim_Tagging}} and
#' \code{\link{Setup_Mod_Tagging}}. Accepts either a scalar (recycled to all
#' release events) or a vector already of length \code{n_events}, and errors
#' on any other length.
#'
#' @param x Numeric scalar or vector.
#' @param n_events Integer. Number of tag release events (cohorts).
#' @param what Character string used in the error message to identify the
#'   offending argument.
#'
#' @return Numeric vector of length \code{n_events}.
#'
#' @keywords internal
recycle_tag_event_par <- function(x, n_events, what) {
  if(length(x) == 1) return(rep(x, n_events))
  if(length(x) != n_events) stop(what, " must be length 1 (recycled to all release events) or length n_tag_rel_events (", n_events, "). Got length ", length(x), ".")
  x
}

#' Map initial tag-induced mortality parameter
#'
#' Internal helper called by \code{\link{Setup_Mod_Tagging}} to construct the
#' TMB/RTMB factor map for \code{ln_init_conv_tag_mort}, the log-scale
#' mortality applied to fish at the moment of tag release. This parameter is
#' a vector of length \code{n_conv_tag_cohorts} (one element per tag release
#' event / row of \code{conv_tag_release_indicator}); when tagging is
#' inactive it is a length-1 placeholder.
#'
#' When no fishery fleet uses conventional tagging
#' (\code{all(use_conv_fish_tagging == 0)}), the parameter is automatically
#' mapped to \code{NA} regardless of \code{init_conv_tag_mort_spec}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and
#'   \code{$map} sublists. Requires \code{$data$use_conv_fish_tagging} and
#'   \code{$par$ln_init_conv_tag_mort}.
#' @param init_conv_tag_mort_spec Character string. One of:
#'   \describe{
#'     \item{\code{"fix"}}{Fix \code{ln_init_conv_tag_mort} at its starting
#'       values for every release event (mapped to \code{NA}).}
#'     \item{\code{"est_shared"}}{Estimate a single value of
#'       \code{ln_init_conv_tag_mort} shared across all release events
#'       (mapped to factor level \code{1} for every event).}
#'     \item{\code{"est_all"}}{Estimate an independent value of
#'       \code{ln_init_conv_tag_mort} for every release event (mapped to
#'       distinct factor levels \code{1:n_conv_tag_cohorts}). Issues a
#'       warning, as per-event initial tag mortality is frequently
#'       non-identifiable.}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_init_conv_tag_mort}
#'   set to a factor vector the same length as
#'   \code{$par$ln_init_conv_tag_mort}. Fixed: all \code{NA}; shared: all
#'   factor level \code{1}; independent: distinct levels per event.
#'
#'
#' @keywords internal
do_conv_init_tag_mort_mapping <- function(input_list, init_conv_tag_mort_spec) {
  n_events <- length(input_list$par$ln_init_conv_tag_mort)
  if(all(input_list$data$use_conv_fish_tagging == 0)) input_list$map$ln_init_conv_tag_mort <- factor(rep(NA, n_events)) # initial tag mortality
  if(any(input_list$data$use_conv_fish_tagging == 1)) {
    # Validate input
    if(!init_conv_tag_mort_spec %in% c("fix", "est_shared", "est_all")) stop("init_conv_tag_mort_spec is incorrectly specified. Should be one of these: fix, est_shared, est_all")
    # Initial tag mortality
    if(init_conv_tag_mort_spec == "fix") input_list$map$ln_init_conv_tag_mort <- factor(rep(NA, n_events))
    if(init_conv_tag_mort_spec == "est_shared") input_list$map$ln_init_conv_tag_mort <- factor(rep(1, n_events))
    if(init_conv_tag_mort_spec == "est_all") {
      input_list$map$ln_init_conv_tag_mort <- factor(1:n_events)
      warning("init_conv_tag_mort_spec = 'est_all' estimates an independent initial tag-induced mortality for each of the ", n_events, " tag release events. This is frequently non-identifiable (initial tag mortality is confounded with natural and fishing mortality, and per-event recaptures are often sparse). Consider 'est_shared' or fixing the parameter unless you have strong prior/data support for event-specific values.")
    }
    collect_message("Conventional Initial Tag Mortality is specified as: ", init_conv_tag_mort_spec)
  }
  return(input_list)
}

#' Map chronic tag shedding parameter
#'
#' Internal helper called by \code{\link{Setup_Mod_Tagging}} to construct the
#' TMB/RTMB factor map for \code{ln_conv_tag_shed}, the log-scale annual
#' chronic tag shedding rate. This parameter is a vector of length
#' \code{n_conv_tag_cohorts} (one element per tag release event / row of
#' \code{conv_tag_release_indicator}); when tagging is inactive it is a
#' length-1 placeholder.
#'
#' When no fishery fleet uses conventional tagging
#' (\code{all(use_conv_fish_tagging == 0)}), the parameter is automatically
#' mapped to \code{NA} regardless of \code{conv_tag_shed_spec}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and
#'   \code{$map} sublists. Requires \code{$data$use_conv_fish_tagging} and
#'   \code{$par$ln_conv_tag_shed}.
#' @param conv_tag_shed_spec Character string. One of:
#'   \describe{
#'     \item{\code{"fix"}}{Fix \code{ln_conv_tag_shed} at its starting values
#'       for every release event (mapped to \code{NA}).}
#'     \item{\code{"est_shared"}}{Estimate a single value of
#'       \code{ln_conv_tag_shed} shared across all release events (mapped to
#'       factor level \code{1} for every event).}
#'     \item{\code{"est_all"}}{Estimate an independent value of
#'       \code{ln_conv_tag_shed} for every release event (mapped to distinct
#'       factor levels \code{1:n_conv_tag_cohorts}). Issues a warning, as
#'       per-event chronic shedding is frequently non-identifiable.}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_conv_tag_shed} set
#'   to a factor vector the same length as \code{$par$ln_conv_tag_shed}.
#'   Fixed: all \code{NA}; shared: all factor level \code{1}; independent:
#'   distinct levels per event.
#'
#'
#' @keywords internal
do_conv_tag_shed_mapping <- function(input_list, conv_tag_shed_spec) {
  n_events <- length(input_list$par$ln_conv_tag_shed)
  if(all(input_list$data$use_conv_fish_tagging == 0)) input_list$map$ln_conv_tag_shed <- factor(rep(NA, n_events)) # chronic tag shedding
  if(any(input_list$data$use_conv_fish_tagging == 1)) {
    # Validate input
    if(!conv_tag_shed_spec %in% c("fix", "est_shared", "est_all")) stop("conv_tag_shed_spec is incorrectly specified. Should be one of these: fix, est_shared, est_all")
    # Tag Shedding
    if(conv_tag_shed_spec == "fix") input_list$map$ln_conv_tag_shed <- factor(rep(NA, n_events))
    if(conv_tag_shed_spec == "est_shared") input_list$map$ln_conv_tag_shed <- factor(rep(1, n_events))
    if(conv_tag_shed_spec == "est_all") {
      input_list$map$ln_conv_tag_shed <- factor(1:n_events)
      warning("conv_tag_shed_spec = 'est_all' estimates an independent chronic tag shedding rate for each of the ", n_events, " tag release events. This is frequently non-identifiable (shedding is confounded with natural and fishing mortality, and per-event recaptures are often sparse). Consider 'est_shared' or fixing the parameter unless you have strong prior/data support for event-specific values.")
    }
    collect_message("Conventional Chronic Tag Shedding is specified as: ", conv_tag_shed_spec)
  }
  return(input_list)
}

#' Map tag recapture overdispersion parameter
#'
#' Internal helper called by \code{\link{Setup_Mod_Tagging}} to construct the
#' TMB/RTMB factor map for \code{ln_conv_fish_tag_theta}, the log-scale
#' overdispersion parameter for the tag recapture likelihood. This is a scalar
#' parameter active only for negative-binomial and Dirichlet-multinomial
#' likelihoods.
#'
#' When no fishery fleet uses conventional tagging, the parameter is mapped to
#' \code{NA}. When tagging is active, activation depends on
#' \code{conv_fish_tag_like}: \code{NA} for Poisson (\code{0}) and multinomial
#' (\code{2}, \code{3}); factor level \code{1} for negative binomial
#' (\code{1}) and Dirichlet-multinomial (\code{4}, \code{5}).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and
#'   \code{$map} sublists. Requires \code{$data$use_conv_fish_tagging} and
#'   \code{$data$conv_fish_tag_like}.
#'
#' @return The input \code{input_list} with \code{$map$ln_conv_fish_tag_theta}
#'   set to a length-1 factor. Active: factor level \code{1}; fixed: \code{NA}.
#'
#'
#' @keywords internal
do_conv_tag_theta_mapping <- function(input_list) {
  if(all(input_list$data$use_conv_fish_tagging == 0)) input_list$map$ln_conv_fish_tag_theta <- factor(NA) # tag overdispersion
  if(any(input_list$data$use_conv_fish_tagging == 1)) {
    # Tag Overdispersion
    if(input_list$data$conv_fish_tag_like %in% c(0,2,3)) input_list$map$ln_conv_fish_tag_theta <- factor(NA)
    if(input_list$data$conv_fish_tag_like %in% c(1,4,5)) input_list$map$ln_conv_fish_tag_theta <- factor(1)
  }
  return(input_list)
}

#' Map fishery tag reporting rate parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Tagging}} to construct the
#' TMB/RTMB factor map for \code{conv_tag_fish_reporting_pars}
#' \code{[n_regions × max_tagrep_blocks × n_fish_fleets]}, the logit-scale
#' tag reporting rate parameters. Parameters are only active for fleets with
#' \code{use_conv_fish_tagging = 1} and likelihoods that use reporting rates
#' (Poisson, negative binomial, or Dirichlet-multinomial conditioned on
#' releases: codes \code{0}, \code{1}, \code{2}, \code{4}).
#'
#' Sharing options that impose region- or fleet-wide pooling
#' (\code{"est_shared_r"}, \code{"est_shared_f"}, \code{"est_shared_r_f"})
#' require that all pooled cells share an identical block structure; a
#' mismatch raises an error. When all fleets have tagging disabled, all
#' parameters are mapped to \code{NA}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and
#'   \code{$map} sublists. Requires \code{$data$use_conv_fish_tagging},
#'   \code{$data$conv_fish_tag_like}, \code{$data$n_regions},
#'   \code{$data$n_fish_fleets}, and \code{$data$conv_tag_fish_reporting_blocks}.
#' @param conv_tagrep_spec Character string specifying the sharing structure
#'   for reporting rate parameters. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate parameter per region, fleet, and
#'       block.}
#'     \item{\code{"est_shared_r"}}{Shared across regions within each fleet
#'       and block. Requires identical block structure across regions within
#'       each fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across active fleets within each
#'       region and block. Requires identical block structure across fleets
#'       within each region.}
#'     \item{\code{"est_shared_r_f"}}{Shared across all regions and active
#'       fleets within each block. Requires identical block structure across
#'       all regions and fleets.}
#'     \item{\code{"fix"}}{All reporting rate parameters fixed at starting
#'       values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with
#'   \code{$map$conv_tag_fish_reporting_pars} set to a factor vector of length
#'   \code{prod(dim(par$conv_tag_fish_reporting_pars))}. Active parameters
#'   receive sequential integer indices; inactive fleet slots and fixed
#'   parameters are \code{NA}.
#'
#'
#' @keywords internal
do_conv_tag_fish_reporting_pars_mapping <- function(input_list, conv_tagrep_spec) {

  if(all(input_list$data$use_conv_fish_tagging == 0)) {
    input_list$map$conv_tag_fish_reporting_pars <- factor(rep(NA, length(input_list$par$conv_tag_fish_reporting_pars)))
  }

  if(any(input_list$data$use_conv_fish_tagging == 1)) {

    first_active_fleet <- min(which(input_list$data$use_conv_fish_tagging == 1))
    map_TagRep <- input_list$par$conv_tag_fish_reporting_pars
    map_TagRep[] <- NA
    tagrep_counter <- 1

    if(input_list$data$conv_fish_tag_like %in% c(0,1,2,4)) {

      # Validate inputs
      if(!conv_tagrep_spec %in% c("est_all", "est_shared_r", "est_shared_f", "est_shared_r_f", "fix"))
        stop("Tag Reporting Specification is not correctly specified. Needs to be fix, est_all, est_shared_r, est_shared_f, or est_shared_r_f")

      # Block consistency checks for sharing specs
      if(conv_tagrep_spec == "est_shared_r") {
        # blocks must match across regions within each fleet
        for(f in 1:input_list$data$n_fish_fleets) {
          ref_blocks <- input_list$data$conv_tag_fish_reporting_blocks[1,,f]
          for(rr in 1:input_list$data$n_regions) {
            if(!identical(as.vector(input_list$data$conv_tag_fish_reporting_blocks[rr,,f]), as.vector(ref_blocks)))
              stop("est_shared_r requires consistent tag reporting block structure across all regions within each fleet.")
          }
        }
      }

      if(conv_tagrep_spec == "est_shared_f") {
        # blocks must match across fleets within each region
        for(r in 1:input_list$data$n_regions) {
          ref_blocks <- input_list$data$conv_tag_fish_reporting_blocks[r,,1]
          for(ff in 1:input_list$data$n_fish_fleets) {
            if(!identical(as.vector(input_list$data$conv_tag_fish_reporting_blocks[r,,ff]), as.vector(ref_blocks)))
              stop("est_shared_f requires consistent tag reporting block structure across all fleets within each region.")
          }
        }
      }

      if(conv_tagrep_spec == "est_shared_r_f") {
        # blocks must match across all regions and fleets
        ref_blocks <- input_list$data$conv_tag_fish_reporting_blocks[1,,1]
        for(r in 1:input_list$data$n_regions) {
          for(f in 1:input_list$data$n_fish_fleets) {
            if(!identical(as.vector(input_list$data$conv_tag_fish_reporting_blocks[r,,f]), as.vector(ref_blocks)))
              stop("est_shared_r_f requires consistent tag reporting block structure across all regions and fleets.")
          }
        }
      }

      if(conv_tagrep_spec == 'fix') map_TagRep[] <- NA

      for(r in 1:input_list$data$n_regions) {
        for(f in 1:input_list$data$n_fish_fleets) {
          tagrep_blocks_tmp <- unique(as.vector(input_list$data$conv_tag_fish_reporting_blocks[r,,f]))

          if(input_list$data$use_conv_fish_tagging[f] == 1) {
            for(b in 1:length(tagrep_blocks_tmp)) {

              if(conv_tagrep_spec == 'est_all') {
                map_TagRep[r,b,f] <- tagrep_counter
                tagrep_counter <- tagrep_counter + 1
              }

              if(conv_tagrep_spec == 'est_shared_r' && r == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  if(tagrep_blocks_tmp[b] %in% input_list$data$conv_tag_fish_reporting_blocks[rr,,f])
                    map_TagRep[rr,b,f] <- tagrep_counter
                }
                tagrep_counter <- tagrep_counter + 1
              }

              if(conv_tagrep_spec == 'est_shared_f' && f == first_active_fleet) {
                for(ff in 1:input_list$data$n_fish_fleets) {
                  if(input_list$data$use_conv_fish_tagging[ff] == 1 &&
                     tagrep_blocks_tmp[b] %in% input_list$data$conv_tag_fish_reporting_blocks[r,,ff])
                    map_TagRep[r,b,ff] <- tagrep_counter
                }
                tagrep_counter <- tagrep_counter + 1
              }

              if(conv_tagrep_spec == 'est_shared_r_f' && r == 1 && f == first_active_fleet) {
                for(rr in 1:input_list$data$n_regions) {
                  for(ff in 1:input_list$data$n_fish_fleets) {
                    if(input_list$data$use_conv_fish_tagging[ff] == 1 &&
                       tagrep_blocks_tmp[b] %in% input_list$data$conv_tag_fish_reporting_blocks[rr,,ff])
                      map_TagRep[rr,b,ff] <- tagrep_counter
                  }
                }
                tagrep_counter <- tagrep_counter + 1
              }

            } # end b loop

          } # end if
        } # end f loop
      } # end r loop

      collect_message("Conventional Tag Reporting is specified as: ", conv_tagrep_spec)
    }

    input_list$map$conv_tag_fish_reporting_pars <- factor(map_TagRep)
  }

  return(input_list)
}


#' Set up the conventional tagging module for model fitting
#'
#' Configures all conventional tagging components of the estimation model:
#' tag release cohort definitions, recapture data, tag recapture likelihood,
#' mixing period, release platform, dimension-attendance and pooling
#' structure, reporting rate time blocks and sharing, and optional reporting
#' rate priors. Delegates parameter mapping to four internal helpers
#' (\code{\link{do_conv_init_tag_mort_mapping}},
#' \code{\link{do_conv_tag_shed_mapping}},
#' \code{\link{do_conv_tag_theta_mapping}},
#' \code{\link{do_conv_tag_fish_reporting_pars_mapping}}). Must be called
#' after \code{\link{Setup_Mod_Dim}} and before model compilation.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param use_conv_fish_tagging Integer vector \code{[n_fish_fleets]} (0/1).
#'   Whether conventional tagging data are included in the likelihood for
#'   each fishery fleet. Default: \code{0} for all fleets.
#' @param conv_tag_release_indicator Integer matrix
#'   \code{[n_conv_tag_cohorts × 3]} giving the release region, year, and
#'   season for each tag cohort. Required when any
#'   \code{use_conv_fish_tagging = 1}. Default \code{NULL}.
#' @param conv_tag_max_liberty Integer. Maximum years-at-liberty included in
#'   the likelihood; recaptures beyond this horizon are ignored. Must be
#'   \eqn{> 0} when tagging is active. Default \code{0}.
#' @param conv_tagged_fish Array \code{[n_conv_tag_cohorts × n_pop × n_ages ×
#'   n_sexes]} of tagged fish released per cohort. Required when any
#'   \code{use_conv_fish_tagging = 1}. Dimensions not attended in
#'   \code{conv_fish_tag_attr} should have all fish placed into index 1 with
#'   remaining indices set to zero. Default \code{NA}.
#' @param obs_conv_tag_fish_recap Array of observed recaptures
#'   \code{[conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop ×
#'   n_regions × n_ages × n_sexes × n_fish_fleets]}. Required when any
#'   \code{use_conv_fish_tagging = 1}. Default \code{NA}.
#' @param conv_fish_tag_like Character string specifying the tag recapture
#'   likelihood. One of \code{"Poisson"}, \code{"NegBin"},
#'   \code{"Multinomial_Release"}, \code{"Multinomial_Recapture"},
#'   \code{"Dirichlet-Multinomial_Release"},
#'   \code{"Dirichlet-Multinomial_Recapture"}. Converted to integer codes
#'   (\code{0}–\code{5}) before storage. Default \code{NA}.
#' @param conv_tag_mixing_period Integer. Minimum number of years (or seasons
#'   in seasonal models) post-release before recaptures contribute to the
#'   likelihood. Allows time for tags to mix within the population before
#'   informing movement estimation. Default \code{1}.
#' @param conv_tag_t_tagging Numeric scalar or vector of length
#'   \code{n_conv_tag_cohorts} (one value per row of
#'   \code{conv_tag_release_indicator}), each in \eqn{[0, 1]}. Fraction of the
#'   season remaining at tag release for that release event. \code{1} = start
#'   of season; \code{0.5} = mid-season; \code{0} = end of season. A scalar is
#'   recycled to all release events. Default \code{1}.
#' @param conv_fish_tag_attr Character scalar or character vector of length
#'   \code{n_conv_tag_cohorts} specifying which biological
#'   dimensions are attended (resolved) at release for each tag release event.
#'   A scalar is recycled to every event; a vector lets different events
#'   resolve different dimensions (e.g. event 1 known at \code{"p_a_s"},
#'   event 2 only at \code{"p_a"}). Each element is built
#'   from any combination of \code{"p"} (population), \code{"a"} (age), and
#'   \code{"s"} (sex), joined by underscores. Region and fleet are always
#'   retained. When a dimension is not attended for an event, all released
#'   fish for that event are placed into index 1 of that dimension and
#'   apportioned to full resolution via the release platform. A dimension may
#'   only be split into more than one pooling group if it is attended in
#'   \emph{every} release event; otherwise the corresponding pooling argument
#'   is overridden to a single group (with a warning). Valid
#'   values: \code{"p_a_s"}, \code{"a_s"}, \code{"p_a"}, \code{"p_s"},
#'   \code{"a"}, \code{"s"}, \code{"p"}, \code{"none"}. Default
#'   \code{"p_a_s"}.
#' @param conv_tag_release_platform Character matrix
#'   \code{[n_conv_tag_cohorts × 2]} specifying the release platform and
#'   fleet index per cohort. Same format as in \code{\link{Setup_Sim_Tagging}}.
#'   Default \code{NULL}.
#' @param conv_tag_pop_pool List of integer vectors defining population pooling
#'   groups for the tagging likelihood. When \code{"p"} is not attended in
#'   \code{conv_fish_tag_attr}, use \code{list(1:n_pop)}. If the pooling
#'   structure is inconsistent with \code{conv_fish_tag_attr}, a warning is
#'   issued and the structure is automatically overridden to a single group.
#'   Default: \code{as.list(1:n_pop)} (population-specific).
#' @param conv_tag_age_pool List of integer vectors defining age pooling
#'   groups. When \code{"a"} is not attended, use \code{list(1:n_ages)}.
#'   Custom groupings (e.g., \code{list(1:5, 6:10)}) are supported when
#'   \code{"a"} is attended. Default: \code{as.list(1:n_ages)}.
#' @param conv_tag_sex_pool List of integer vectors defining sex pooling
#'   groups. When \code{"s"} is not attended, use \code{list(1:n_sexes)}.
#'   Default: \code{as.list(1:n_sexes)}.
#' @param init_conv_tag_mort_spec Character string (\code{"fix"},
#'   \code{"est_shared"}, or \code{"est_all"}). Whether initial tag-induced
#'   mortality is fixed at its starting values, estimated as a single value
#'   shared across all release events, or estimated independently for every
#'   release event. See \code{\link{do_conv_init_tag_mort_mapping}}. Default
#'   \code{NULL}.
#' @param conv_tag_shed_spec Character string (\code{"fix"},
#'   \code{"est_shared"}, or \code{"est_all"}). Whether chronic tag shedding
#'   is fixed at its starting values, estimated as a single value shared
#'   across all release events, or estimated independently for every release
#'   event. See \code{\link{do_conv_tag_shed_mapping}}. Default \code{NULL}.
#' @param conv_tag_fish_reporting_blocks Character vector defining time blocks
#'   for fishery tag reporting rates. Each element follows one of:
#'   \describe{
#'     \item{\code{"none_Region_r_Fleet_f"}}{Constant reporting rate for
#'       region \code{r} and fleet \code{f}.}
#'     \item{\code{"Block_b_Year_y1-y2_Region_r_Fleet_f"}}{Block \code{b}
#'       applies to years \code{y1}–\code{y2}. Use \code{"terminal"} for the
#'       end year to extend to the final model year.}
#'   }
#'   Parsed into an array \code{[n_regions × n_years × n_fish_fleets]}.
#'   If \code{NULL}, a single constant block is used for all region–fleet
#'   combinations. Default \code{NULL}.
#' @param conv_tagrep_spec Character string. Sharing structure for reporting
#'   rate parameters \code{conv_tag_fish_reporting_pars}
#'   \code{[n_regions × max_tagrep_blocks × n_fish_fleets]}. See
#'   \code{\link{do_conv_tag_fish_reporting_pars_mapping}} for full option
#'   descriptions. Default \code{"fix"} (a warning is issued if this was
#'   unintentional).
#' @param use_conv_tag_fishrep_prior Integer (0/1). Whether priors are applied
#'   to reporting rate parameters. Default \code{0}.
#' @param conv_tag_fishrep_prior Data frame of prior specifications for
#'   reporting rates. Required columns: \code{region}, \code{block},
#'   \code{fleet}, \code{mu}, \code{sd}, \code{type}. Ignored when
#'   \code{use_conv_tag_fishrep_prior = 0}. Default \code{NULL}.
#' @param ... Optional named starting values for parameters. Supported names
#'   and defaults:
#'   \code{ln_init_conv_tag_mort} (scalar or length \code{n_conv_tag_cohorts}
#'     vector; a scalar is recycled to all release events; default
#'     \code{-1000}),
#'   \code{ln_conv_tag_shed} (scalar or length \code{n_conv_tag_cohorts}
#'     vector; a scalar is recycled to all release events; default
#'     \code{-1000}),
#'   \code{ln_conv_fish_tag_theta} (scalar, default \code{0}),
#'   \code{conv_tag_fish_reporting_pars}
#'     \code{[n_regions × max_tagrep_blocks × n_fish_fleets]},
#'   default \code{0} (logit scale ≈ 0.5 reporting probability; inactive
#'   fleet slots overwritten to \code{-1000}).
#'
#' @return The input \code{input_list} with tagging configuration stored in
#'   \code{$data} (\code{use_conv_fish_tagging}, \code{conv_tag_release_indicator},
#'   \code{n_conv_tag_cohorts}, \code{conv_tag_max_liberty},
#'   \code{conv_tagged_fish}, \code{obs_conv_tag_fish_recap},
#'   \code{conv_fish_tag_like}, \code{conv_tag_mixing_period},
#'   \code{conv_tag_t_tagging}, \code{use_conv_tag_fishrep_prior},
#'   \code{conv_tag_fishrep_prior}, \code{conv_tag_pop_pool},
#'   \code{conv_tag_age_pool}, \code{conv_tag_sex_pool},
#'   \code{conv_tag_fish_reporting_blocks}, \code{conv_fish_tag_attr},
#'   \code{conv_tag_release_platform}); starting values in \code{$par} for
#'   \code{ln_init_conv_tag_mort} and \code{ln_conv_tag_shed} (each length
#'   \code{n_conv_tag_cohorts}, or length 1 when tagging is inactive),
#'   \code{ln_conv_fish_tag_theta}, and \code{conv_tag_fish_reporting_pars};
#'   and factor maps in \code{$map} for all four parameter arrays.
#'
#'
#' @export Setup_Mod_Tagging
#' @family Model Setup
Setup_Mod_Tagging <- function(input_list,
                              use_conv_fish_tagging = rep(0, input_list$data$n_fish_fleets),
                              conv_tag_release_indicator = NULL,
                              conv_tag_max_liberty = 0,
                              conv_tagged_fish = NA,
                              obs_conv_tag_fish_recap = NA,
                              conv_fish_tag_like = NA,
                              conv_tag_mixing_period = 1,
                              conv_tag_t_tagging = 1,
                              use_conv_tag_fishrep_prior = 0,
                              conv_tag_fishrep_prior = NULL,
                              conv_tag_pop_pool = as.list(1:input_list$data$n_pop),
                              conv_tag_age_pool = as.list(1:length(input_list$data$ages)),
                              conv_tag_sex_pool = as.list(1:input_list$data$n_sexes),
                              init_conv_tag_mort_spec = NULL,
                              conv_tag_shed_spec = NULL,
                              conv_tagrep_spec = 'fix',
                              conv_tag_fish_reporting_blocks = NULL,
                              conv_fish_tag_attr = 'p_a_s',
                              conv_tag_release_platform = NULL,
                              ...
                              ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Tagging <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Checking data and other specifications
  if(any(use_conv_fish_tagging == 1)) {

    # check specifications
    if(is.na(sum(conv_tagged_fish))) stop("No data is provided for conv_tagged_fish")
    if(is.na(sum(obs_conv_tag_fish_recap))) stop("No data is provided for obs_conv_tag_fish_recap")
    if(is.na(conv_fish_tag_like)) stop("No likelihood is provided for conv_fish_tag_like")
    if(conv_tag_max_liberty == 0) stop("conv_tag_max_liberty must be greater than 0")
    if(conv_tagrep_spec == 'fix') warning("Note that tag reporting rates is fixed. Specify est_all or est_shared_r if this was not the intention.")

    # Check data
    check_data_dimensions(conv_tagged_fish, n_conv_tag_cohorts = nrow(conv_tag_release_indicator), n_ages = length(input_list$data$ages),
                          n_sexes = input_list$data$n_sexes, n_pop = input_list$data$n_pop, what = 'conv_tagged_fish')

    check_data_dimensions(obs_conv_tag_fish_recap, conv_tag_max_liberty = conv_tag_max_liberty,  n_pop = input_list$data$n_pop,
                          n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_fish_fleets = input_list$data$n_fish_fleets,
                          n_conv_tag_cohorts = nrow(conv_tag_release_indicator), n_ages = length(input_list$data$ages),
                          n_sexes = input_list$data$n_sexes, what = 'obs_conv_tag_fish_recap')
  }

  # Checking tagging priors
  if(use_conv_tag_fishrep_prior == 1) {
    required_cols <- c("region", "block", "fleet", "mu", "sd", 'type')
    missing_cols <- setdiff(required_cols, names(conv_tag_fishrep_prior))
    if(length(missing_cols) > 0) {
      stop("conv_tag_fishrep_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    collect_message("Tagging priors are used")
  }

  # Checking tag likelihoods
  tag_like_map <- data.frame(type = c("Poisson", "NegBin", "Multinomial_Release", "Multinomial_Recapture",
                                      "Dirichlet-Multinomial_Release", "Dirichlet-Multinomial_Recapture"), num = c(0,1,2,3,4,5))

  if(is.na(conv_fish_tag_like)) conv_fish_tag_like_vals <- 999
  else {
    if(!conv_fish_tag_like %in% c(tag_like_map$type)) stop("Tag Likelihood not correctly specified. Should be one of these: Poisson, NegBin, Multinomial_Release, Multinomial_Recapture, Dirichlet-Multinomial_Release, Dirichlet-Multinomial_Recapture")
    conv_fish_tag_like_vals <- tag_like_map$num[tag_like_map$type == conv_fish_tag_like]
    collect_message("Conventional Tag Likelihood specified as: ", conv_fish_tag_like)
  }


  # Tag Pooling Options -----------------------------------------------------
  # If movement is pooled either across pops, sexes or ages
  if(is.character(conv_tag_pop_pool)){
    if(conv_tag_pop_pool == "all") move_pop_tag_pool_vals = list(1:input_list$data$n_pop)
  } else move_pop_tag_pool_vals = conv_tag_pop_pool

  if(is.character(conv_tag_age_pool)){
    if(conv_tag_age_pool == "all") move_age_tag_pool_vals = list(input_list$data$ages)
  } else move_age_tag_pool_vals = conv_tag_age_pool

  if(is.character(conv_tag_sex_pool)){
    if(conv_tag_sex_pool == "all") move_sex_tag_pool_vals = list(1:input_list$data$n_sexes)
  } else move_sex_tag_pool_vals = conv_tag_sex_pool

  # Recycle conv_fish_tag_attr to one attended-dimension string per release event
  # (scalar recycled to all events). Different release events may resolve
  # different dimensions (e.g. event 1 known at p_a_s, event 2 only at p_a).
  # n_conv_tag_cohorts is populated further below, so derive the event count
  # locally here (matching that later logic) since the pooling check needs it.
  n_events_attr <- if(all(use_conv_fish_tagging == 0)) 1L else nrow(conv_tag_release_indicator)
  conv_fish_tag_attr <- recycle_tag_event_par(conv_fish_tag_attr, max(n_events_attr, 1), "conv_fish_tag_attr")

  # Enforce consistency between pooling structure and conv_fish_tag_attr.
  # A dimension may only be split into >1 pooling group if it is attended in
  # EVERY release event; if any event does not resolve that dimension, the
  # recapture likelihood cannot distinguish groups along it, so pool fully.
  # Warn the user and correct automatically if not.
  attr_parts_list <- strsplit(conv_fish_tag_attr, "_")
  attended_all <- function(dim) all(vapply(attr_parts_list, function(x) dim %in% x, logical(1)))

  if(!attended_all("p") && length(move_pop_tag_pool_vals) > 1) {
    warning("conv_tag_pop_pool has more than one group but 'p' is not attended in every release event's conv_fish_tag_attr. ",
            "Overriding to list(1:n_pop) to pool all populations.")
    move_pop_tag_pool_vals <- list(1:input_list$data$n_pop)
  }

  if(!attended_all("a") && length(move_age_tag_pool_vals) > 1) {
    warning("conv_tag_age_pool has more than one group but 'a' is not attended in every release event's conv_fish_tag_attr. ",
            "Overriding to list(1:n_ages) to pool all ages.")
    move_age_tag_pool_vals <- list(seq_along(input_list$data$ages))
  }

  if(!attended_all("s") && length(move_sex_tag_pool_vals) > 1) {
    warning("conv_tag_sex_pool has more than one group but 's' is not attended in every release event's conv_fish_tag_attr. ",
            "Overriding to list(1:n_sexes) to pool all sexes.")
    move_sex_tag_pool_vals <- list(1:input_list$data$n_sexes)
  }

  collect_message("Conventional Tagging data are fit to ", length(move_pop_tag_pool_vals), " population groups")
  collect_message("Conventional Tagging data are fit to ", length(move_age_tag_pool_vals), " age groups")
  collect_message("Conventional Tagging data are fit to ", length(move_sex_tag_pool_vals), " sex groups")

  # Tag Reporting Rates Options ---------------------------------------------
  conv_tag_fish_reporting_blocks_mat <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))

  if(!is.null(conv_tag_fish_reporting_blocks)) {
    for(i in 1:length(conv_tag_fish_reporting_blocks)) {

      # Extract out components from list
      tmp <- conv_tag_fish_reporting_blocks[i]
      tmp_vec <- unlist(strsplit(tmp, "_"))

      if(!tmp_vec[1] %in% c("none", "Block")) stop("Tag Reporting Blocks not correctly specified. This should be either none_Region_x_Fleet_x or Block_x_Year_x-y_Region_x_Fleet_x")

      # extract out fleets if constant
      if(tmp_vec[1] == "none") {
        region <- as.numeric(tmp_vec[3]) # get region index
        fleet <- as.numeric(tmp_vec[5]) # get fleet index
        conv_tag_fish_reporting_blocks_mat[region,,fleet] <- 1 # input tag reporting time block
      }

      if(tmp_vec[1] == "Block") {
        block_val <- as.numeric(tmp_vec[2]) # get block value
        region <- as.numeric(tmp_vec[6]) # get region value
        fleet <- as.numeric(tmp_vec[8]) # get fleet value

        # get year ranges
        if(!str_detect(tmp, "terminal")) { # if not terminal year
          year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
          years <- year_range[1]:year_range[2] # get sequence of years
        } else { # if terminal year
          year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
          years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
        }

        conv_tag_fish_reporting_blocks_mat[region,years,fleet] <- block_val # input tag reporting time block
      }

    } # end i loop
  } else conv_tag_fish_reporting_blocks_mat[] <- 1

   for(f in 1:input_list$data$n_fish_fleets) for(r in 1:input_list$data$n_regions)
     collect_message("Conventional Tag Reporting estimated with ", length(unique(conv_tag_fish_reporting_blocks_mat[r,,f])), " blocks for region ", r, " and fleet ", f)


  # Populate Data List ------------------------------------------------------

  input_list$data$use_conv_fish_tagging <- use_conv_fish_tagging
  input_list$data$conv_tag_release_indicator <- conv_tag_release_indicator
  if(all(use_conv_fish_tagging == 0)) input_list$data$n_conv_tag_cohorts <- 0
  if(any(use_conv_fish_tagging == 1)) input_list$data$n_conv_tag_cohorts <- nrow(conv_tag_release_indicator)
  input_list$data$conv_tag_max_liberty <- conv_tag_max_liberty
  input_list$data$conv_tagged_fish <- conv_tagged_fish
  input_list$data$obs_conv_tag_fish_recap <- obs_conv_tag_fish_recap
  input_list$data$conv_fish_tag_like <- conv_fish_tag_like_vals
  input_list$data$conv_tag_mixing_period <- conv_tag_mixing_period
  input_list$data$conv_tag_t_tagging <- recycle_tag_event_par(conv_tag_t_tagging, max(input_list$data$n_conv_tag_cohorts, 1), "conv_tag_t_tagging")
  input_list$data$use_conv_tag_fishrep_prior <- use_conv_tag_fishrep_prior
  input_list$data$conv_tag_fishrep_prior <- conv_tag_fishrep_prior
  input_list$data$conv_tag_pop_pool <- move_pop_tag_pool_vals
  input_list$data$conv_tag_age_pool <- move_age_tag_pool_vals
  input_list$data$conv_tag_sex_pool <- move_sex_tag_pool_vals
  input_list$data$conv_tag_fish_reporting_blocks <- conv_tag_fish_reporting_blocks_mat
  input_list$data$conv_fish_tag_attr <- conv_fish_tag_attr
  input_list$data$conv_tag_release_platform <- conv_tag_release_platform

  # Populate Parameter List ------------------------------------------------------

  n_tag_par_events <- max(input_list$data$n_conv_tag_cohorts, 1) # one value per release event; length-1 placeholder when tagging is inactive

  # Initial tag induced mortality
  if("ln_init_conv_tag_mort" %in% names(starting_values)) input_list$par$ln_init_conv_tag_mort <- recycle_tag_event_par(starting_values$ln_init_conv_tag_mort, n_tag_par_events, "ln_init_conv_tag_mort")
  else input_list$par$ln_init_conv_tag_mort <- rep(-1000, n_tag_par_events)

  # Chronic tag shedding
  if("ln_conv_tag_shed" %in% names(starting_values)) input_list$par$ln_conv_tag_shed <- recycle_tag_event_par(starting_values$ln_conv_tag_shed, n_tag_par_events, "ln_conv_tag_shed")
  else input_list$par$ln_conv_tag_shed <- rep(-1000, n_tag_par_events)

  # tag overdispersion parameter
  if("ln_conv_fish_tag_theta" %in% names(starting_values)) input_list$par$ln_conv_fish_tag_theta <- starting_values$ln_conv_fish_tag_theta
  else input_list$par$ln_conv_fish_tag_theta <- 0

  # tag reporting rate parameters
  max_tagrep_blks <- max(apply(input_list$data$conv_tag_fish_reporting_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of tag reporting rate blocks for each region and fleet
  if("conv_tag_fish_reporting_pars" %in% names(starting_values)) input_list$par$conv_tag_fish_reporting_pars <- starting_values$conv_tag_fish_reporting_pars
  else input_list$par$conv_tag_fish_reporting_pars <- array(0, dim = c(input_list$data$n_regions, max_tagrep_blks, input_list$data$n_fish_fleets)) # specified at 0.5 in inverse logit space

  # For fleets with no tagging data, set reporting pars to -1000 on the logit scale
  # (effectively zero reporting probability; these will be fixed via map)
  if(any(input_list$data$use_conv_fish_tagging == 0)) input_list$par$conv_tag_fish_reporting_pars[,, which(input_list$data$use_conv_fish_tagging == 0)] <- -1000

  # Mapping Options ---------------------------------------------------------
  input_list <- do_conv_init_tag_mort_mapping(input_list, init_conv_tag_mort_spec)
  input_list <- do_conv_tag_shed_mapping(input_list, conv_tag_shed_spec)
  input_list <- do_conv_tag_theta_mapping(input_list)
  input_list <- do_conv_tag_fish_reporting_pars_mapping(input_list, conv_tagrep_spec)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose && any(use_conv_fish_tagging == 1)) for(msg in messages_list) message(msg)

  return(input_list)
}
