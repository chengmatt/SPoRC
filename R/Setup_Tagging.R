#' Set up simulated tagging dynamics
#'
#' @param n_tags Number of tags to release in a given year (scalar, default = NULL)
#' @param max_liberty Maximum liberty (years) to track cohorts (default = sim_list$n_ages / 2)
#' @param t_tagging Fraction of season remaining when tags are released (e.g., start of season == 1, mid season == 0.5, end of season == 0; default = 1)
#' @param ln_Init_Tag_Mort Log initial tag-induced mortality (default = -1000)
#' @param ln_Tag_Shed Log chronic tag shedding rate (default = -1000) annual rate
#' @param sim_list Simulation list (required)
#' @param UseTagging Boolean to use tagging (default = 0):
#'   \itemize{
#'     \item 0: Do not simulate tagging
#'     \item 1: Simulate tagging
#'   }
#' @param Tag_Reporting_input Tag reporting input [n_regions × n_yrs × n_sims]
#'   (default = 0.5)
#' @param n_tags_rel_input Number of tag releases by tag cohort length (default = NULL)
#' @param tag_selex Tag selectivity type (default = 5):
#'   \itemize{
#'     \item \code{0} or \code{"Uniform_DomFleet"}: Uniform by age/sex, dominant fleet
#'     \item \code{1} or \code{"SexAgg_DomFleet"}: Sex-aggregated, dominant fleet
#'     \item \code{2} or \code{"SexSp_DomFleet"}: Sex-specific, dominant fleet
#'     \item \code{3} or \code{"Uniform_AllFleet"}: Uniform by age/sex, all fleets
#'     \item \code{4} or \code{"SexAgg_AllFleet"}: Sex-aggregated, all fleets
#'     \item \code{5} or \code{"SexSp_AllFleet"}: Sex-specific, all fleets
#'   }
#' @param tag_natmort Tag natural mortality type (default = 3):
#'   \itemize{
#'     \item \code{0} or \code{"AgeAgg_SexAgg"}: Age-aggregated, sex-aggregated
#'     \item \code{1} or \code{"AgeSp_SexAgg"}: Age-specific, sex-aggregated
#'     \item \code{2} or \code{"AgeAgg_SexSp"}: Age-aggregated, sex-specific
#'     \item \code{3} or \code{"AgeSp_SexSp"}: Age-specific, sex-specific
#'   }
#' @param tag_like Tag likelihood type (default = 0):
#'   \itemize{
#'     \item \code{0} or \code{"Poisson"}: Poisson
#'     \item \code{1} or \code{"NegBin"}: Negative Binomial
#'     \item \code{2} or \code{"Multinomial_Release"}: Multinomial by release cohort
#'     \item \code{3} or \code{"Multinomial_Recapture"}: Multinomial by recapture event
#'     \item \code{4} or \code{"Dirichlet-Multinomial_Release"}: Dirichlet-Multinomial by release cohort
#'     \item \code{5} or \code{"Dirichlet-Multinomial_Recapture"}: Dirichlet-Multinomial by recapture event
#'   }
#' @param ln_tag_theta Scalar in log space describing tag likelihood overdispersion (default = log(1))
#'
#' @export Setup_Sim_Tagging
#' @family Simulation Setup
Setup_Sim_Tagging <- function(n_tags = NULL,
                              n_tags_rel_input = NULL,
                              UseTagging = 0,
                              max_liberty = sim_list$n_ages / 2,
                              tag_release_indicator = expand.grid(regions = 1:sim_list$n_regions, tag_years = 1:sim_list$n_yrs, tag_seas = 1:sim_list$n_seas),
                              t_tagging = 1,
                              ln_Init_Tag_Mort = -1000,
                              ln_Tag_Shed = -1000,
                              Tag_Reporting_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
                              tag_selex = 5,
                              tag_natmort = 3,
                              tag_like = 0,
                              ln_tag_theta = log(1),
                              sim_list
                              ) {

  # Convert codes to numeric
  tag_selex <- convert_to_numeric(tag_selex, list(Uniform_DomFleet = 0, SexAgg_DomFleet = 1, SexSp_DomFleet = 2, Uniform_AllFleet = 3, SexAgg_AllFleet = 4, SexSp_AllFleet = 5))
  tag_natmort <- convert_to_numeric(tag_natmort, list(AgeAgg_SexAgg = 0, AgeSp_SexAgg = 1, AgeAgg_SexSp = 2, AgeSp_SexSp = 3))
  tag_like <- convert_to_numeric(tag_like, list(Poisson = 0, NegBin = 1, Multinomial_Release = 2,
                                                Multinomial_Recapture = 3, `Dirichlet-Multinomial_Release` = 4, `Dirichlet-Multinomial_Recapture` = 5))

  if(!is.null(Tag_Reporting_input)) check_sim_dimensions(Tag_Reporting_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_sims = sim_list$n_sims, what = "Tag_Reporting_input")

  # Output variables into list
  if(!is.null(n_tags) && !is.null(n_tags_rel_input)) stop("n_tags and n_tags_rel_input cannot be specified simultaneously. n_tags is a scalar, while n_tags_rel_input specifies cohort-specific tags!")
  if(!is.null(n_tags)) sim_list$n_tags <- n_tags
  if(!is.null(n_tags_rel_input)) sim_list$n_tags_rel_input <- n_tags_rel_input
  sim_list$max_liberty <- max_liberty
  sim_list$t_tagging <- t_tagging # time of tagging
  sim_list$ln_Init_Tag_Mort <- ln_Init_Tag_Mort # tag induced mortality
  sim_list$ln_Tag_Shed <- ln_Tag_Shed # tag shedding
  sim_list$tag_release_indicator <- tag_release_indicator # tag release indicator (by tag years and regions = a tag cohort)
  sim_list$n_tag_rel_events <- nrow(tag_release_indicator) # number of tag release events - tag years x tag region (tag cohorts)

  # Containers
  sim_list$Tagged_Fish <- array(0, dim = c(sim_list$n_tag_rel_events, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # number of tagged fish
  sim_list$Tags_Avail <- array(0, dim = c(sim_list$max_liberty + 1, sim_list$n_seas, sim_list$n_tag_rel_events, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # tags availiable for recapture every year
  sim_list$Pred_Tag_Recap <- array(0, dim = c(sim_list$max_liberty, sim_list$n_seas, sim_list$n_tag_rel_events, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # predicted tag recaptures
  sim_list$Obs_Tag_Recap <- array(0, dim = c(sim_list$max_liberty, sim_list$n_seas, sim_list$n_tag_rel_events, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # observed tag recaptures

  sim_list$Tag_Reporting <- Tag_Reporting_input # output this into list
  sim_list$UseTagging <- UseTagging # output into list
  sim_list$tag_selex <- tag_selex # output into list
  sim_list$tag_natmort <- tag_natmort # otuput into list
  sim_list$tag_like <- tag_like # tag likelihood
  sim_list$ln_tag_theta <- ln_tag_theta # tag likelihood overdispersion parameter

  return(sim_list)
}

#' Helper function to setup initial tag mortality parameters
#'
#' @param input_list Input list
#' @param Init_Tag_Mort_spec Character vector specifying initial tag mortality parameterization
#' @keywords internal
do_Init_Tag_Mort_mapping <- function(input_list, Init_Tag_Mort_spec) {
  if(input_list$data$UseTagging == 0) input_list$map$ln_Init_Tag_Mort <- factor(NA) # initial tag mortality
  if(input_list$data$UseTagging == 1) {
    # Validate input
    if(!Init_Tag_Mort_spec %in% c("fix", "est")) stop("Init_Tag_Mort_spec is incorrectly specified. Should be one of these: fix, est")
    # Initial tag mortality
    if(Init_Tag_Mort_spec == "fix") input_list$map$ln_Init_Tag_Mort <- factor(NA)
    if(Init_Tag_Mort_spec == "est") input_list$map$ln_Init_Tag_Mort <- factor(1)
    collect_message("Initial Tag Mortality is specified as: ", Init_Tag_Mort_spec)
  }
  return(input_list)
}

#' Helper function to set up tag shedding parameters
#'
#' @param input_list Input list
#' @param Tag_Shed_spec Character specifying tag shedding parameterization
#' @keywords internal
do_Tag_Shed_mapping <- function(input_list, Tag_Shed_spec) {
  if(input_list$data$UseTagging == 0) input_list$map$ln_Tag_Shed <- factor(NA) # chronic tag shedding
  if(input_list$data$UseTagging == 1) {
    # Validate input
    if(!Tag_Shed_spec %in% c("fix", "est")) stop("Tag_Shed_spec is incorrectly specified. Should be one of these: fix, est")
    # Tag Shedding
    if(Tag_Shed_spec == "fix" || UseTagging == 0) input_list$map$ln_Tag_Shed <- factor(NA)
    if(Tag_Shed_spec == "est") input_list$map$ln_Tag_Shed <- factor(1)
    collect_message("Chronic Tag Shedding is specified as: ", Tag_Shed_spec)
  }
  return(input_list)
}

#' Helper function to set up tag overdispersion parameters
#'
#' @param input_list Input list
#' @keywords internal
do_tag_theta_mapping <- function(input_list) {
  if(input_list$data$UseTagging == 0) input_list$map$ln_tag_theta <- factor(NA) # tag overdispersion
  if(input_list$data$UseTagging == 1) {
    # Tag Overdispersion
    if(input_list$data$Tag_LikeType %in% c(0,2,3)) input_list$map$ln_tag_theta <- factor(NA)
    if(input_list$data$Tag_LikeType %in% c(1,4,5)) input_list$map$ln_tag_theta <- factor(1)
  }
  return(input_list)
}

#' Helper function to set up tag reporting rate parameters
#'
#' @param input_list Input list
#' @param TagRep_spec Charcacter specifying tag reporting parameterization
#' @keywords internal
do_Tag_Reporting_Pars_mapping <- function(input_list, TagRep_spec) {

  # If not using tagging data
  if(input_list$data$UseTagging == 0) {
    input_list$map$Tag_Reporting_Pars <- factor(rep(NA, length(input_list$par$Tag_Reporting_Pars))) # tag reporting rates
  }

  if(input_list$data$UseTagging == 1) {

    # Initialize arrays and counters
    map_TagRep <- input_list$par$Tag_Reporting_Pars
    map_TagRep[] <- NA
    tagrep_counter <- 1

    if(input_list$data$Tag_LikeType %in% c(0,1,2,4)) { # If this is a poisson, negative binomial, multinomial or dirichlet-multinomial release conditioned

      # Validate inputs here
      if(!TagRep_spec %in% c("est_all", "est_shared_r", "fix")) stop("Tag Reporting Specificaiton is not correctly specified. Needs to be fix, est_all, or est_shared_r")

      # if we want to fix
      if(TagRep_spec == 'fix') map_TagRep[] <- NA

      for(r in 1:input_list$data$n_regions) {

        # Get number of tag reporting rate blocks
        tagrep_blocks_tmp <- unique(as.vector(input_list$data$Tag_Reporting_blocks[r,]))

        for(b in 1:length(tagrep_blocks_tmp)) {

          # Estimate for all regions
          if(TagRep_spec == 'est_all') {
            map_TagRep[r,b] <- tagrep_counter
            tagrep_counter <- tagrep_counter + 1
          }

          # Estimate but share tag reporting across regions
          if(TagRep_spec == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              # only assign if this value exists for this region
              if(tagrep_blocks_tmp[b] %in% input_list$data$Tag_Reporting_blocks[rr,]) {
                map_TagRep[rr, b] <- tagrep_counter
              } # end if
            } # end rr loop
            tagrep_counter <- tagrep_counter + 1
          }

        } # end b loop
      } # end r loop

      collect_message("Tag Reporting is specified as: ", TagRep_spec)

    } # end if

    # input tag reporting rates into mapping list
    input_list$map$Tag_Reporting_Pars <- factor(map_TagRep) # tag reporting rates

  }

  return(input_list)
}

#' Setup tagging processes and parameters
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param UseTagging Numeric (0 or 1) indicating whether to use tagging data (1) or not (0)
#' @param tag_release_indicator Matrix [n_tag_cohorts x 3], where columns are release region, release year, and release season
#' @param max_tag_liberty Maximum number of years to track a tagged cohort
#' @param Tagged_Fish Array [n_tag_cohorts x n_ages x n_sexes] describing tagged fish releases
#' @param Obs_Tag_Recap Array [max_tag_liberty x n_seas x n_tag_cohorts x n_regions x n_ages x n_sexes] observed tag recaptures
#' @param Tag_LikeType Character string specifying tag likelihood type. One of:
#'   \itemize{
#'     \item \code{"Poisson"}
#'     \item \code{"NegBin"}
#'     \item \code{"Multinomial_Release"}
#'     \item \code{"Multinomial_Recapture"}
#'     \item \code{"Dirichlet-Multinomial_Release"}
#'     \item \code{"Dirichlet-Multinomial_Recapture"}
#'   }
#'   Example: \code{Tag_LikeType = "NegBin"}
#' @param mixing_period Numeric indicating minimum years post-release to include in fitting (or minimum seasons post-release if model is seasonal)
#' @param t_tagging Fraction of season remaining when tags are released (e.g., start of season == 1, mid season == 0.5, end of season == 0; default = 1)
#' @param tag_selex Character string specifying tag recovery selectivity. One of:
#'   \itemize{
#'     \item \code{"Uniform_DomFleet"}
#'     \item \code{"SexAgg_DomFleet"}
#'     \item \code{"SexSp_DomFleet"}
#'     \item \code{"Uniform_AllFleet"}
#'     \item \code{"SexAgg_AllFleet"}
#'     \item \code{"SexSp_AllFleet"}
#'   }
#'   Example: \code{tag_selex = "SexSp_AllFleet"}
#' @param tag_natmort Character string specifying tag natural mortality parameterization. One of:
#'   \itemize{
#'     \item \code{"AgeAgg_SexAgg"}
#'     \item \code{"AgeSp_SexAgg"}
#'     \item \code{"AgeAgg_SexSp"}
#'     \item \code{"AgeSp_SexSp"}
#'   }
#'   Example: \code{tag_natmort = "AgeSp_SexSp"}
#' @param Use_TagRep_Prior Numeric (0 or 1) whether to use tag reporting rate prior
#' @param move_age_tag_pool List or character specifying pooling of tagging data by age groups. Default does not pool ages. Examples:
#'   \itemize{
#'     \item \code{list(1:5, 6:11, 12:20)} pools these age groups together
#'     \item \code{"all"} pools all ages together (internally converted to \code{list(1:n_ages)})
#'     \item \code{as.list(1:n_ages)} fits each sex separately
#'   }
#' @param move_sex_tag_pool List or character specifying pooling of tagging data by sex groups. Default do not pool sexes. Examples:
#'   \itemize{
#'     \item \code{list(1:2)} pools sexes together
#'     \item \code{"all"} pools all sexes together (internally converted to \code{list(1:n_sexes)})
#'     \item \code{list(1, 2)} fits each sex separately
#'   }
#' @param Init_Tag_Mort_spec Character string \code{"fix"} or \code{"est"} specifying if initial tag mortality is fixed or estimated
#' @param Tag_Shed_spec Character string \code{"fix"} or \code{"est"} specifying if chronic tag shedding is fixed or estimated
#' @param Tag_Reporting_blocks Character vector specifying blocks of years and regions for tag reporting rates. Default is a single block for all regions. Format examples:
#'   \itemize{
#'     \item \code{"Block_1_Year_1-15_Region_1"}
#'     \item \code{"Block_2_Year_16-terminal_Region_2"}
#'     \item \code{"none_Region_3"} (means no block, constant for that region)
#'   }
#' @param TagRep_spec Character string specifying tag reporting rate estimation scheme:
#'   \itemize{
#'     \item \code{"est_all"} estimates rates for all blocks and regions independently
#'     \item \code{"est_shared_r"} estimates rates shared across regions but varying by block
#'     \item \code{"fix"} fixes all reporting rates (no estimation)
#'   }
#' @param ... Additional starting values for tagging parameters such as \code{ln_Init_Tag_Mort}, \code{ln_Tag_Shed}, \code{ln_tag_theta}, \code{Tag_Reporting_Pars}
#' @param TagRep_Prior Data frame containing prior specifications for tag reporting parameters.
#'   Must include columns: \code{region} (region index), \code{block} (time block index),
#'   \code{mu} (Numeric mean for tag reporting prior (normal space); \code{NA} if symmetric beta is used),
#'   \code{sd} (Numeric standard deviation for tag reporting prior (normal space)), and \code{type} (0 == symmetric beta, 1 == regular beta).
#'   Each row specifies a beta prior for one tag reporting parameter.
#'   Only parameters with rows in this data frame will have priors applied.
#'
#' @export Setup_Mod_Tagging
#' @family Model Setup
Setup_Mod_Tagging <- function(input_list,
                              UseTagging = 0,
                              tag_release_indicator = NULL,
                              max_tag_liberty = 0,
                              Tagged_Fish = NA,
                              Obs_Tag_Recap = NA,
                              Tag_LikeType = NA,
                              mixing_period = 1,
                              t_tagging = 1,
                              tag_selex = NA,
                              tag_natmort = NA,
                              Use_TagRep_Prior = 0,
                              TagRep_Prior = NULL,
                              move_age_tag_pool = as.list(1:length(input_list$data$ages)),
                              move_sex_tag_pool = as.list(1:input_list$data$n_sexes),
                              Init_Tag_Mort_spec = NULL,
                              Tag_Shed_spec = NULL,
                              TagRep_spec = 'fix',
                              Tag_Reporting_blocks = NULL,
                              ...
                              ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Checking data and other specifications
  if(UseTagging == 1) {

    # check specifications
    if(is.na(sum(Tagged_Fish))) stop("No data is provided for Tagged_Fish")
    if(is.na(sum(Obs_Tag_Recap))) stop("No data is provided for Obs_Tag_Recap")
    if(is.na(Tag_LikeType)) stop("No likelihood is provided for Tag_LikeType")
    if(is.na(tag_selex)) stop("No specfication is provided for tag_selex")
    if(is.na(tag_natmort)) stop("No specfication is provided for tag_natmort")
    if(max_tag_liberty == 0) stop("max_tag_liberty must be greater than 0")
    if(TagRep_spec == 'fix') warning("Note that tag reporting rates is fixed. Specify est_all or est_shared_r if this was not the intention.")

    # Check data
    check_data_dimensions(Tagged_Fish, n_tag_cohorts = nrow(tag_release_indicator), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'Tagged_Fish')
    check_data_dimensions(Obs_Tag_Recap, max_tag_liberty = max_tag_liberty, n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_tag_cohorts = nrow(tag_release_indicator), n_ages = length(input_list$data$ages),
                          n_sexes = input_list$data$n_sexes, what = 'Obs_Tag_Recap')
  }

  # Checking tagging priors
  if(Use_TagRep_Prior == 1) {
    required_cols <- c("region", "block", "mu", "sd", 'type')
    missing_cols <- setdiff(required_cols, names(TagRep_Prior))
    if(length(missing_cols) > 0) {
      stop("TagRep_Prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    collect_message("Tagging priors are used")
  }

  # Checking tag likelihoods
  tag_like_map <- data.frame(type = c("Poisson", "NegBin", "Multinomial_Release", "Multinomial_Recapture",
                                      "Dirichlet-Multinomial_Release", "Dirichlet-Multinomial_Recapture"), num = c(0,1,2,3,4,5))

  if(is.na(Tag_LikeType)) Tag_LikeType_vals <- 999
  else {
    if(!Tag_LikeType %in% c(tag_like_map$type)) stop("Tag Likelihood not correctly specified. Should be one of these: Poisson, NegBin, Multinomial_Release, Multinomial_Recapture, Dirichlet-Multinomial_Release, Dirichlet-Multinomial_Recapture")
    Tag_LikeType_vals <- tag_like_map$num[tag_like_map$type == Tag_LikeType]
    collect_message("Tag Likelihood specified as: ", Tag_LikeType)
  }

  # Setup tagging selectivity
  tag_selex_map <- data.frame(type = c("Uniform_DomFleet", "SexAgg_DomFleet", "SexSp_DomFleet", "Uniform_AllFleet", "SexAgg_AllFleet", "SexSp_AllFleet"), num = c(0,1,2,3,4,5))

  if(is.na(Tag_LikeType)) tag_selex_vals <- 999
  else {
    if(!tag_selex %in% c(tag_selex_map$type)) stop("Tag Selectivity not correctly specified. Should be one of these: Uniform_DomFleet, SexAgg_DomFleet, SexSp_DomFleet, Uniform_AllFleet, SexAgg_AllFleet, SexSp_AllFleet")
    tag_selex_vals <- tag_selex_map$num[tag_selex_map$type == tag_selex]
    collect_message("Tag Selectivity specified as: ", tag_selex)
  }

  # Checking tagging natural moratlity
  tag_natmort_map <- data.frame(type = c("AgeAgg_SexAgg", "AgeSp_SexAgg", "AgeAgg_SexSp", "AgeSp_SexSp"), num = c(0,1,2,3))

  if(is.na(Tag_LikeType)) tag_natmort_vals <- 999
  else {
    if(!tag_natmort %in% c(tag_natmort_map$type)) stop("Tag Natural Mortality not correctly specified. Should be one of these: AgeAgg_SexAgg, AgeSp_SexAgg, AgeAgg_SexSp, AgeSp_SexSp")
    tag_natmort_vals <- tag_natmort_map$num[tag_natmort_map$type == tag_natmort]
     collect_message("Tag Natural Mortality specified as: ", tag_natmort)
  }


  # Tag Pooling Options -----------------------------------------------------
  # If movement is pooled either across sexes or ages
  if(is.character(move_age_tag_pool)){
    if(move_age_tag_pool == "all") move_age_tag_pool_vals = list(input_list$data$ages)
  } else move_age_tag_pool_vals = move_age_tag_pool

  if(is.character(move_sex_tag_pool)){
    if(move_sex_tag_pool == "all") move_sex_tag_pool_vals = list(1:input_list$data$n_sexes)
  } else move_sex_tag_pool_vals = move_sex_tag_pool

  collect_message("Tagging data are fit to ", length(move_age_tag_pool_vals), " age groups")
  collect_message("Tagging data are fit to ", length(move_sex_tag_pool_vals), " sex groups")

  # Tag Reporting Rates Options ---------------------------------------------
  Tag_Reporting_blocks_mat <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years)))

  if(!is.null(Tag_Reporting_blocks)) {
    for(i in 1:length(Tag_Reporting_blocks)) {

      # Extract out components from list
      tmp <- Tag_Reporting_blocks[i]
      tmp_vec <- unlist(strsplit(tmp, "_"))

      if(!tmp_vec[1] %in% c("none", "Block")) stop("Tag Reporting Blocks not correctly specified. This should be either none_Region_x or Block_x_Year_x-y_Region_x")

      # extract out fleets if constant
      if(tmp_vec[1] == "none") {
        region <- as.numeric(tmp_vec[3]) # get region index
        Tag_Reporting_blocks_mat[region,] <- 1 # input tag reporting time block
      }

      if(tmp_vec[1] == "Block") {
        block_val <- as.numeric(tmp_vec[2]) # get block value
        region <- as.numeric(tmp_vec[6]) # get region value

        # get year ranges
        if(!str_detect(tmp, "terminal")) { # if not terminal year
          year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
          years <- year_range[1]:year_range[2] # get sequence of years
        } else { # if terminal year
          year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
          years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
        }

        Tag_Reporting_blocks_mat[region,years] <- block_val # input tag reporting time block
      }

    } # end i loop
  } else Tag_Reporting_blocks_mat[] <- 1

   for(r in 1:input_list$data$n_regions) collect_message("Tag Reporting estimated with ", length(unique(Tag_Reporting_blocks_mat[r,])), " block for region ", r)


  # Populate Data List ------------------------------------------------------

  input_list$data$UseTagging <- UseTagging
  input_list$data$tag_release_indicator <- tag_release_indicator
  if(UseTagging == 0) input_list$data$n_tag_cohorts <- 0
  if(UseTagging == 1) input_list$data$n_tag_cohorts <- nrow(tag_release_indicator)
  input_list$data$max_tag_liberty <- max_tag_liberty
  input_list$data$Tagged_Fish <- Tagged_Fish
  input_list$data$Obs_Tag_Recap <- Obs_Tag_Recap
  input_list$data$Tag_LikeType <- Tag_LikeType_vals
  input_list$data$mixing_period <- mixing_period
  input_list$data$t_tagging <- t_tagging
  input_list$data$tag_selex <- tag_selex_vals
  input_list$data$tag_natmort <- tag_natmort_vals
  input_list$data$Use_TagRep_Prior <- Use_TagRep_Prior
  input_list$data$TagRep_Prior <- TagRep_Prior
  input_list$data$move_age_tag_pool <- move_age_tag_pool_vals
  input_list$data$move_sex_tag_pool <- move_sex_tag_pool_vals
  input_list$data$Tag_Reporting_blocks <- Tag_Reporting_blocks_mat

  # Populate Parameter List ------------------------------------------------------

  # Initial tag induced mortality
  if("ln_Init_Tag_Mort" %in% names(starting_values)) input_list$par$ln_Init_Tag_Mort <- starting_values$ln_Init_Tag_Mort
  else input_list$par$ln_Init_Tag_Mort <- -1000

  # Chronic tag shedding
  if("ln_Tag_Shed" %in% names(starting_values)) input_list$par$ln_Tag_Shed <- starting_values$ln_Tag_Shed
  else input_list$par$ln_Tag_Shed <- -1000

  # tag overdispersion parameter
  if("ln_tag_theta" %in% names(starting_values)) input_list$par$ln_tag_theta <- starting_values$ln_tag_theta
  else input_list$par$ln_tag_theta <- 0

  # tag reporting rate parameters
  max_tagrep_blks <- max(apply(input_list$data$Tag_Reporting_blocks, 1, FUN = function(x) length(unique(x)))) # figure out maximum number of tag reporting rate blocks for each region
  if("Tag_Reporting_Pars" %in% names(starting_values)) input_list$par$Tag_Reporting_Pars <- starting_values$Tag_Reporting_Pars
  else input_list$par$Tag_Reporting_Pars <- array(0, dim = c(input_list$data$n_regions, max_tagrep_blks)) # specified at 0.5 in inverse logit space


  # Mapping Options ---------------------------------------------------------
  input_list <- do_Init_Tag_Mort_mapping(input_list, Init_Tag_Mort_spec)
  input_list <- do_Tag_Shed_mapping(input_list, Tag_Shed_spec)
  input_list <- do_tag_theta_mapping(input_list)
  input_list <- do_Tag_Reporting_Pars_mapping(input_list, TagRep_spec)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose && UseTagging == 1) for(msg in messages_list) message(msg)

  return(input_list)
}
