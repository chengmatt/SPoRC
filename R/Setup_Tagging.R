#' Set Up Simulated Tagging Dynamics
#'
#' Initializes conventional tagging processes for a simulation, including tag
#' releases, tag-induced mortality, tag shedding, reporting rates, and
#' containers for predicted and observed recaptures.
#'
#' This function defines tag release cohorts (by region, year, and season),
#' sets assumptions about tag mortality and shedding, and allocates arrays
#' required to track tagged fish and recapture observations through time.
#'
#' @param n_tags Numeric scalar. Total number of tags released per release
#'   event. Cannot be specified simultaneously with `n_tags_rel_input`.
#'
#' @param n_tags_rel_input Optional array specifying the number of tags
#'   released by cohort (e.g., by length, age, or other structure defined in
#'   the operating model). Overrides `n_tags` if supplied. Cannot be specified
#'   simultaneously with `n_tags`.
#'
#' @param use_conv_fish_tagging Integer (0/1). Indicator for whether
#'   conventional tag recaptures from fisheries are simulated.
#'
#' @param max_liberty Maximum number of years-at-liberty tracked for each tag
#'   cohort. Default is `sim_list$n_ages / 2`.
#'
#' @param tag_release_indicator Data frame defining tag release events
#'   (i.e., tag cohorts). Default expands across all regions, years, and
#'   seasons in `sim_list`. Must contain columns:
#'   \describe{
#'     \item{regions}{Region of release}
#'     \item{tag_years}{Year of release}
#'     \item{tag_seas}{Season of release}
#'   }
#'
#' @param tag_release_platform Character matrix specifying the release
#'   platform and fleet number associated with each tag release event.
#'   Rows must align with `tag_release_indicator`. The first column
#'   (`platform`) must be one of:
#'   \describe{
#'     \item{"population"}{Tags are released directly into the population
#'     (i.e., independent of a fleet or survey sampling process).}
#'     \item{"fishery"}{Tags are released through a fishery fleet. The
#'     corresponding fleet number must be provided in the second column.}
#'     \item{"survey"}{Tags are released through a survey fleet. The
#'     corresponding fleet number must be provided in the second column.}
#'   }
#'   The second column (`fleet`) identifies the fleet index when
#'   `platform` is `"fishery"` or `"survey"`, and may be set to `NA`
#'   when `platform` is `"population"`.
#'
#' @param t_tagging Numeric scalar in [0,1]. Fraction of the season remaining
#'   when tags are released:
#'   \itemize{
#'     \item 1 = start of season
#'     \item 0.5 = mid-season
#'     \item 0 = end of season
#'   }
#'
#' @param ln_init_conv_tag_mort Log initial tag-induced mortality applied at
#'   release. Default of -1000 approximates zero mortality.
#'
#' @param ln_conv_tag_shed Log annual chronic tag shedding rate. Default of -1000
#'   approximates no shedding.
#'
#' @param conv_fish_tag_like Integer specifying likelihood used for tag
#'   recapture observations:
#'   \describe{
#'     \item Poisson
#'     \item Negative binomial
#'     \item Multinomial (release)
#'     \item Multinomial (recapture)
#'     \item Dirichlet-multinomial (release)
#'     \item Dirichlet-multinomial (recapture)
#'   }
#'
#' @param ln_conv_fish_tag_theta Log overdispersion parameter for
#'   negative-binomial or Dirichlet-multinomial likelihoods.
#'
#' @param sim_list A list passed on from previous Setup functions.
#'
#' @param conv_fish_tag_attr Character string specifying which biological attributes
#'   are recorded at recapture. Constructed from any combination of \code{"p"} (population),
#'   \code{"a"} (age), and \code{"s"} (sex), joined by underscores. Region and fleet are
#'   always retained. Supported values: \code{"p_a_s"}, \code{"a_s"}, \code{"p_a"},
#'   \code{"p_s"}, \code{"a"}, \code{"s"}, \code{"p"}, \code{"none"}.
#'
#' @param conv_tag_fish_reporting_input Numeric array of tag reporting rates for
#'   fishery fleets with dimensions \code{[n_regions, n_yrs, n_fish_fleets, n_sims]}.
#'   Values represent the probability that a recaptured tag is reported by the
#'   fishery and should be betweem \code{[0, 1]}. Reporting rates can vary by region,
#'   year, fleet, and simulation replicate.
#'
#' @export Setup_Sim_Tagging
#' @family Simulation Setup
Setup_Sim_Tagging <- function(n_tags = NULL,
                              n_tags_rel_input = NULL,
                              use_conv_fish_tagging = 0,
                              max_liberty = sim_list$n_ages / 2,
                              tag_release_indicator = expand.grid(regions = 1:sim_list$n_regions, tag_years = 1:sim_list$n_yrs, tag_seas = 1:sim_list$n_seas),
                              tag_release_platform = matrix(c("survey", "1"), nrow = nrow(tag_release_indicator),  ncol = 2, byrow = TRUE,dimnames = list(NULL, c("platform", "fleet"))),
                              t_tagging = 1,
                              ln_init_conv_tag_mort = -1000,
                              ln_conv_tag_shed = -1000,
                              conv_fish_tag_attr = 'p_a_s',
                              conv_tag_fish_reporting_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              conv_fish_tag_like = 0,
                              ln_conv_fish_tag_theta = log(1),
                              sim_list
                              ) {

  check_sim_dimensions(conv_tag_fish_reporting_input, n_regions = sim_list$n_regions,
                       n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = "conv_tag_fish_reporting_input")

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
  sim_list$max_liberty <- max_liberty
  sim_list$t_tagging <- t_tagging # time of tagging
  sim_list$ln_init_conv_tag_mort <- ln_init_conv_tag_mort # tag induced mortality
  sim_list$ln_conv_tag_shed <- ln_conv_tag_shed # tag shedding
  sim_list$tag_release_indicator <- tag_release_indicator # tag release indicator (by tag years and regions = a tag cohort)
  sim_list$tag_release_platform <- tag_release_platform # how tags are released
  sim_list$n_tag_rel_events <- nrow(tag_release_indicator) # number of tag release events - tag years x tag region (tag cohorts)

  # Containers
  sim_list$conv_conv_tagged_fish <-
    array(0, dim = c(sim_list$n_tag_rel_events, sim_list$n_pop,
                     sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # number of tagged fish

  sim_list$conv_tag_fish_avail <-
    array(0, dim = c(sim_list$max_liberty + 1, sim_list$n_seas, sim_list$n_tag_rel_events,
                     sim_list$n_pop, sim_list$n_regions,
                     sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # tags availiable for recapture every year

  sim_list$obs_conv_tag_fish_recap <- sim_list$pred_conv_tag_fish_recap <-
    array(0, dim = c(sim_list$max_liberty, sim_list$n_seas, sim_list$n_tag_rel_events,
                     sim_list$n_pop, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes,
                     sim_list$n_fish_fleets, sim_list$n_sims)) # predicted tag recaptures

  sim_list$conv_tag_fish_reporting <- conv_tag_fish_reporting_input # output this into list
  sim_list$use_conv_fish_tagging <- use_conv_fish_tagging # output into list
  sim_list$conv_fish_tag_like <- conv_fish_tag_like # tag likelihood
  sim_list$conv_fish_tag_attr <- conv_fish_tag_attr # tag recpature attributes for fishery conventional tags
  sim_list$ln_conv_fish_tag_theta <- ln_conv_fish_tag_theta # tag likelihood overdispersion parameter

  return(sim_list)
}

#' Helper function to setup initial tag mortality parameters
#'
#' @param input_list Input list
#' @param init_conv_tag_mort_spec Character vector specifying initial tag mortality parameterization
#' @keywords internal
do_conv_init_tag_mort_mapping <- function(input_list, init_conv_tag_mort_spec) {
  if(input_list$data$use_conv_fish_tagging == 0) input_list$map$ln_init_conv_tag_mort <- factor(NA) # initial tag mortality
  if(input_list$data$use_conv_fish_tagging == 1) {
    # Validate input
    if(!init_conv_tag_mort_spec %in% c("fix", "est")) stop("init_conv_tag_mort_spec is incorrectly specified. Should be one of these: fix, est")
    # Initial tag mortality
    if(init_conv_tag_mort_spec == "fix") input_list$map$ln_init_conv_tag_mort <- factor(NA)
    if(init_conv_tag_mort_spec == "est") input_list$map$ln_init_conv_tag_mort <- factor(1)
    collect_message("Conventional Initial Tag Mortality is specified as: ", init_conv_tag_mort_spec)
  }
  return(input_list)
}

#' Helper function to set up tag shedding parameters
#'
#' @param input_list Input list
#' @param conv_tag_shed_spec Character specifying tag shedding parameterization
#' @keywords internal
do_conv_tag_shed_mapping <- function(input_list, conv_tag_shed_spec) {
  if(input_list$data$use_conv_fish_tagging == 0) input_list$map$ln_conv_tag_shed <- factor(NA) # chronic tag shedding
  if(input_list$data$use_conv_fish_tagging == 1) {
    # Validate input
    if(!conv_tag_shed_spec %in% c("fix", "est")) stop("conv_tag_shed_spec is incorrectly specified. Should be one of these: fix, est")
    # Tag Shedding
    if(conv_tag_shed_spec == "fix" || use_conv_fish_tagging == 0) input_list$map$ln_conv_tag_shed <- factor(NA)
    if(conv_tag_shed_spec == "est") input_list$map$ln_conv_tag_shed <- factor(1)
    collect_message("Conventional Chronic Tag Shedding is specified as: ", conv_tag_shed_spec)
  }
  return(input_list)
}

#' Helper function to set up tag overdispersion parameters
#'
#' @param input_list Input list
#' @keywords internal
do_conv_tag_theta_mapping <- function(input_list) {
  if(input_list$data$use_conv_fish_tagging == 0) input_list$map$ln_conv_fish_tag_theta <- factor(NA) # tag overdispersion
  if(input_list$data$use_conv_fish_tagging == 1) {
    # Tag Overdispersion
    if(input_list$data$conv_fish_tag_like %in% c(0,2,3)) input_list$map$ln_conv_fish_tag_theta <- factor(NA)
    if(input_list$data$conv_fish_tag_like %in% c(1,4,5)) input_list$map$ln_conv_fish_tag_theta <- factor(1)
  }
  return(input_list)
}

#' Helper function to set up tag reporting rate parameters
#'
#' @param input_list Input list
#' @param conv_tagrep_spec Charcacter specifying tag reporting parameterization
#' @keywords internal
do_conv_tag_fish_reporting_pars_mapping <- function(input_list, conv_tagrep_spec) {

  if(input_list$data$use_conv_fish_tagging == 0) {
    input_list$map$conv_tag_fish_reporting_pars <- factor(rep(NA, length(input_list$par$conv_tag_fish_reporting_pars)))
  }

  if(input_list$data$use_conv_fish_tagging == 1) {

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

            if(conv_tagrep_spec == 'est_shared_f' && f == 1) {
              for(ff in 1:input_list$data$n_fish_fleets) {
                if(tagrep_blocks_tmp[b] %in% input_list$data$conv_tag_fish_reporting_blocks[r,,ff])
                  map_TagRep[r,b,ff] <- tagrep_counter
              }
              tagrep_counter <- tagrep_counter + 1
            }

            if(conv_tagrep_spec == 'est_shared_r_f' && r == 1 && f == 1) {
              for(rr in 1:input_list$data$n_regions) {
                for(ff in 1:input_list$data$n_fish_fleets) {
                  if(tagrep_blocks_tmp[b] %in% input_list$data$conv_tag_fish_reporting_blocks[rr,,ff])
                    map_TagRep[rr,b,ff] <- tagrep_counter
                }
              }
              tagrep_counter <- tagrep_counter + 1
            }

          } # end b loop
        } # end f loop
      } # end r loop

      collect_message("Conventional Tag Reporting is specified as: ", conv_tagrep_spec)
    }

    input_list$map$conv_tag_fish_reporting_pars <- factor(map_TagRep)
  }

  return(input_list)
}

#' Set Up Model-Based Conventional Tagging
#'
#' Configures conventional tagging data, likelihood structure, parameter
#' mappings, and optional priors for model fitting.
#'
#' This function:
#' \itemize{
#'   \item Validates tagging data inputs,
#'   \item Defines likelihood type,
#'   \item Specifies parameterizations for tag mortality and shedding,
#'   \item Configures reporting-rate blocks and sharing structure,
#'   \item Applies optional priors to reporting parameters.
#' }
#'
#' @param input_list List containing \code{$data}, \code{$par}, and \code{$map}.
#'
#' @param use_conv_fish_tagging Integer (0/1) indicating whether tagging data
#'   are included in model fitting.
#'
#' @param tag_release_indicator Matrix [n_tag_cohorts x 3] giving release
#'   region, year, and season.
#'
#' @param max_tag_liberty Maximum number of years-at-liberty included in fitting.
#'
#' @param conv_tagged_fish Array describing tagged fish releases with
#'   dimensions:
#'   \code{[n_tag_cohorts, n_pop, n_ages, n_sexes]}.
#'
#' @param obs_conv_tag_fish_recap Array of observed recaptures with dimensions:
#'   \code{[max_tag_liberty, n_seas, n_tag_cohorts, n_pop,
#'           n_regions, n_ages, n_sexes]}.
#'
#' @param conv_fish_tag_like Character specifying likelihood type.
#'   One of:
#'   \code{"Poisson"}, \code{"NegBin"},
#'   \code{"Multinomial_Release"}, \code{"Multinomial_Recapture"},
#'   \code{"Dirichlet-Multinomial_Release"},
#'   \code{"Dirichlet-Multinomial_Recapture"}.
#'
#' @param mixing_period Minimum years (or seasons if seasonal model)
#'   post-release included in fitting.
#'
#' @param t_tagging Numeric scalar in [0,1]. Fraction of the season remaining
#'   when tags are released:
#'   \itemize{
#'     \item 1 = start of season
#'     \item 0.5 = mid-season
#'     \item 0 = end of season
#'   }
#'
#' @param use_conv_tag_fishrep_prior Integer (0/1) indicating whether priors
#'   are applied to reporting parameters.
#'
#' @param conv_tag_fishrep_prior Data frame specifying priors on reporting
#'   parameters. Must include columns:
#'   \code{region}, \code{block}, \code{fleet},
#'   \code{mu}, \code{sd}, and \code{type}.
#'
#' @param move_age_tag_pool List defining age pooling structure for tagging
#'   likelihood. Examples:
#'   \itemize{
#'     \item \code{list(1:5, 6:10)}
#'     \item \code{"all"}
#'     \item \code{as.list(1:n_ages)}
#'   }
#'
#' @param move_sex_tag_pool List defining sex pooling structure.
#'
#' @param init_conv_tag_mort_spec Character string \code{"fix"} or \code{"est"}
#'   specifying whether initial tag mortality is fixed or estimated.
#'
#' @param conv_tag_shed_spec Character string \code{"fix"} or \code{"est"}
#'   specifying whether chronic tag shedding is fixed or estimated.
#'
#' @param conv_tag_fish_reporting_blocks Character vector defining reporting
#'   rate blocks. Examples:
#'   \itemize{
#'     \item \code{"none_Region_1_Fleet_1"}
#'     \item \code{"Block_2_Year_1-20_Region_1_Fleet_1"}
#'     \item \code{"Block_3_Year_21-terminal_Region_2_Fleet_2"}
#'   }
#'
#' @param conv_tagrep_spec Character specifying reporting-rate sharing scheme:
#'   \describe{
#'     \item{est_all}{Estimate independently for all regions, fleets, and blocks}
#'     \item{est_shared_r}{Share across regions within fleet}
#'     \item{est_shared_f}{Share across fleets within region}
#'     \item{est_shared_r_f}{Share across all regions and fleets}
#'     \item{fix}{Fix reporting rates}
#'   }
#'
#' @param ... Optional starting values for tagging parameters
#'   (\code{ln_init_conv_tag_mort}, \code{ln_conv_tag_shed},
#'   \code{ln_conv_fish_tag_theta}, \code{conv_tag_fish_reporting_pars}).
#'
#'
#' @export Setup_Mod_Tagging
#' @family Model Setup
Setup_Mod_Tagging <- function(input_list,
                              use_conv_fish_tagging = 0,
                              tag_release_indicator = NULL,
                              max_tag_liberty = 0,
                              conv_tagged_fish = NA,
                              obs_conv_tag_fish_recap = NA,
                              conv_fish_tag_like = NA,
                              mixing_period = 1,
                              t_tagging = 1,
                              use_conv_tag_fishrep_prior = 0,
                              conv_tag_fishrep_prior = NULL,
                              move_age_tag_pool = as.list(1:length(input_list$data$ages)),
                              move_sex_tag_pool = as.list(1:input_list$data$n_sexes),
                              init_conv_tag_mort_spec = NULL,
                              conv_tag_shed_spec = NULL,
                              conv_tagrep_spec = 'fix',
                              conv_tag_fish_reporting_blocks = NULL,
                              ...
                              ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Checking data and other specifications
  if(use_conv_fish_tagging == 1) {

    # check specifications
    if(is.na(sum(conv_tagged_fish))) stop("No data is provided for conv_tagged_fish")
    if(is.na(sum(obs_conv_tag_fish_recap))) stop("No data is provided for obs_conv_tag_fish_recap")
    if(is.na(conv_fish_tag_like)) stop("No likelihood is provided for conv_fish_tag_like")
    if(max_tag_liberty == 0) stop("max_tag_liberty must be greater than 0")
    if(conv_tagrep_spec == 'fix') warning("Note that tag reporting rates is fixed. Specify est_all or est_shared_r if this was not the intention.")

    # Check data
    check_data_dimensions(conv_tagged_fish, n_tag_cohorts = nrow(tag_release_indicator), n_ages = length(input_list$data$ages),
                          n_sexes = input_list$data$n_sexes, n_pop = input_list$data$n_pop, what = 'conv_tagged_fish')

    check_data_dimensions(obs_conv_tag_fish_recap, max_tag_liberty = max_tag_liberty,  n_pop = input_list$data$n_pop,
                          n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions,
                          n_tag_cohorts = nrow(tag_release_indicator), n_ages = length(input_list$data$ages),
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
  # If movement is pooled either across sexes or ages
  if(is.character(move_age_tag_pool)){
    if(move_age_tag_pool == "all") move_age_tag_pool_vals = list(input_list$data$ages)
  } else move_age_tag_pool_vals = move_age_tag_pool

  if(is.character(move_sex_tag_pool)){
    if(move_sex_tag_pool == "all") move_sex_tag_pool_vals = list(1:input_list$data$n_sexes)
  } else move_sex_tag_pool_vals = move_sex_tag_pool

  collect_message("Conventional Tagging data are fit to ", length(move_age_tag_pool_vals), " age groups")
  collect_message("Conventional Tagging data are fit to ", length(move_sex_tag_pool_vals), " sex groups")

  # Tag Reporting Rates Options ---------------------------------------------
  conv_tag_fish_reporting_blocks_mat <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))

  if(!is.null(conv_tag_fish_reporting_blocks)) {
    for(i in 1:length(conv_tag_fish_reporting_blocks)) {

      # Extract out components from list
      # conv_tag_fish_reporting_blocks = c('none_Region_1_Fleet_1', 'Block_2_Year_1-35_Region2_Fleet_2')
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
  input_list$data$tag_release_indicator <- tag_release_indicator
  if(use_conv_fish_tagging == 0) input_list$data$n_tag_cohorts <- 0
  if(use_conv_fish_tagging == 1) input_list$data$n_tag_cohorts <- nrow(tag_release_indicator)
  input_list$data$max_tag_liberty <- max_tag_liberty
  input_list$data$conv_tagged_fish <- conv_tagged_fish
  input_list$data$obs_conv_tag_fish_recap <- obs_conv_tag_fish_recap
  input_list$data$conv_fish_tag_like <- conv_fish_tag_like_vals
  input_list$data$mixing_period <- mixing_period
  input_list$data$t_tagging <- t_tagging
  input_list$data$use_conv_tag_fishrep_prior <- use_conv_tag_fishrep_prior
  input_list$data$conv_tag_fishrep_prior <- conv_tag_fishrep_prior
  input_list$data$move_age_tag_pool <- move_age_tag_pool_vals
  input_list$data$move_sex_tag_pool <- move_sex_tag_pool_vals
  input_list$data$conv_tag_fish_reporting_blocks <- conv_tag_fish_reporting_blocks_mat

  # Populate Parameter List ------------------------------------------------------

  # Initial tag induced mortality
  if("ln_init_conv_tag_mort" %in% names(starting_values)) input_list$par$ln_init_conv_tag_mort <- starting_values$ln_init_conv_tag_mort
  else input_list$par$ln_init_conv_tag_mort <- -1000

  # Chronic tag shedding
  if("ln_conv_tag_shed" %in% names(starting_values)) input_list$par$ln_conv_tag_shed <- starting_values$ln_conv_tag_shed
  else input_list$par$ln_conv_tag_shed <- -1000

  # tag overdispersion parameter
  if("ln_conv_fish_tag_theta" %in% names(starting_values)) input_list$par$ln_conv_fish_tag_theta <- starting_values$ln_conv_fish_tag_theta
  else input_list$par$ln_conv_fish_tag_theta <- 0

  # tag reporting rate parameters
  max_tagrep_blks <- max(apply(input_list$data$conv_tag_fish_reporting_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of tag reporting rate blocks for each region and fleet
  if("conv_tag_fish_reporting_pars" %in% names(starting_values)) input_list$par$conv_tag_fish_reporting_pars <- starting_values$conv_tag_fish_reporting_pars
  else input_list$par$conv_tag_fish_reporting_pars <- array(0, dim = c(input_list$data$n_regions, max_tagrep_blks, input_list$data$n_fish_fleets)) # specified at 0.5 in inverse logit space


  # Mapping Options ---------------------------------------------------------
  input_list <- do_conv_init_tag_mort_mapping(input_list, init_conv_tag_mort_spec)
  input_list <- do_conv_tag_shed_mapping(input_list, conv_tag_shed_spec)
  input_list <- do_conv_tag_theta_mapping(input_list)
  input_list <- do_conv_tag_fish_reporting_pars_mapping(input_list, conv_tagrep_spec)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose && use_conv_fish_tagging == 1) for(msg in messages_list) message(msg)

  return(input_list)
}
