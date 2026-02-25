#' Construct Parameter Mapping for Movement Fixed Effects
#'
#' Internal helper used by \code{Setup_Mod_Movement()} to construct the
#' TMB parameter mapping for movement parameters under the unstructured
#' Markov formulation (\code{move_type = 0}).
#'
#' This function defines how movement parameters are shared across
#' population, age, year, season, and sex blocks. When blocks are specified,
#' parameters within each block share a common estimated value via a common
#' factor level in the TMB \code{map} list.
#'
#' @param input_list A model input list containing \code{$data}, \code{$par},
#'   and \code{$map} components. This object is modified and returned.
#' @param Movement_popblk_spec Either \code{"constant"} (population-invariant
#'   movement) or a \code{list} of integer vectors defining population blocks
#'   that share parameters.
#' @param Movement_ageblk_spec Either \code{"constant"} (age-invariant
#'   movement) or a \code{list} of integer vectors defining age blocks.
#' @param Movement_yearblk_spec Either \code{"constant"} (time-invariant
#'   movement) or a \code{list} of integer vectors defining year blocks.
#' @param Movement_sexblk_spec Either \code{"constant"} (sex-invariant
#'   movement) or a \code{list} of integer vectors defining sex blocks.
#' @param Movement_seasblk_spec Either \code{"constant"} (season-invariant
#'   movement) or a \code{list} of integer vectors defining seasonal blocks.
#' @param use_fixed_movement Integer flag (0 or 1). If \code{1}, movement is
#'   treated as fixed and parameters are not estimated.
#'
#' @details
#' For \code{move_type = 0}, movement is parameterized as an unstructured
#' Markov transition among regions. Parameters are defined for each
#' population × region_from × region_to × year × season × age × sex
#' combination (excluding self-transitions).
#'
#' Block specifications determine which indices share a common parameter.
#' Internally, unique integer identifiers are assigned to each block
#' combination and converted to factors for TMB.
#'
#' For \code{move_type = 1} (CTMC), all unstructured movement parameters are
#' mapped to \code{NA}, and only CTMC-specific parameters are estimated.
#'
#' @return The updated \code{input_list} with movement-related entries in
#'   \code{$map} populated.
#'
#' @keywords internal
do_move_pars_mapping <- function(input_list, Movement_popblk_spec,
                                 Movement_ageblk_spec, Movement_yearblk_spec,
                                 Movement_sexblk_spec, Movement_seasblk_spec, use_fixed_movement) {

  # Setup mapping list
  map_Movement_Pars <- input_list$par$move_pars # initialize array with same dimensions as parameters
  map_log_move_diffusion_pars <- input_list$par$log_move_diffusion_pars # initialize array with same dimensions as parameters
  map_move_preference_pars <- input_list$par$move_preference_pars # initialize array with same dimensions as parameters

  if(input_list$data$move_type == 0) {

    # Setup dimensions
    n_regions_from <- dim(map_Movement_Pars)[2]
    n_regions_to <- dim(map_Movement_Pars)[3]

    # If movement is constant for populations
    if(is.character(Movement_popblk_spec)){
      if(Movement_popblk_spec == "constant") Movement_popblk_spec_vals <- list(1:input_list$data$n_pop)
    } else Movement_popblk_spec_vals <- Movement_popblk_spec


    # If movement is constant for ages
    if(is.character(Movement_ageblk_spec)){
      if(Movement_ageblk_spec == "constant") Movement_ageblk_spec_vals <- list(input_list$data$ages)
    } else Movement_ageblk_spec_vals <- Movement_ageblk_spec

    # If movement is constant across years
    if(is.character(Movement_yearblk_spec)){
      if(Movement_yearblk_spec == "constant") Movement_yearblk_spec_vals = list(input_list$data$years)
    } else Movement_yearblk_spec_vals = Movement_yearblk_spec

    # If movement is constant across sexes
    if(is.character(Movement_sexblk_spec)){
      if(Movement_sexblk_spec == "constant") Movement_sexblk_spec_vals <- list(1:input_list$data$n_sexes)
    } else Movement_sexblk_spec_vals <- Movement_sexblk_spec

    # If movement is constant across seasons
    if(is.character(Movement_seasblk_spec)){
      if(Movement_seasblk_spec == "constant") Movement_seasblk_spec_vals <- list(1:input_list$data$n_seas)
    } else Movement_seasblk_spec_vals <- Movement_seasblk_spec

    # If spatial model
    if(input_list$data$n_regions > 1 &&
       input_list$data$use_fixed_movement == 0 # if not using fixed movement matrix
       ){

      # Initialize counter
      counter <- 1

      for(popblk in 1:length(Movement_popblk_spec_vals)) {
        # get populations to block and map off
        map_p <- Movement_popblk_spec_vals[[popblk]]

        for(ageblk in 1:length(Movement_ageblk_spec_vals)) {
          # get ages to block and map off
          map_a <- Movement_ageblk_spec_vals[[ageblk]]

          for(yearblk in 1:length(Movement_yearblk_spec_vals)) {
            # get years to block and map off
            map_y <- Movement_yearblk_spec_vals[[yearblk]]

            for(seasblk in 1:length(Movement_seasblk_spec_vals)) {
              # get seasons to block and map off
              map_seas <- Movement_seasblk_spec_vals[[seasblk]]

              for(sexblk in 1:length(Movement_sexblk_spec_vals)) {
                # get sexes to block and map off
                map_s <- Movement_sexblk_spec_vals[[sexblk]]

                # Now, loop through each combination and increment get unique indices
                map_idx <- array(0, dim = c(n_regions_from, n_regions_to))

                # Each region from and to has a new counter variable
                for(i in 1:n_regions_from) {
                  for(j in 1:n_regions_to) {
                    map_idx[i,j] <- counter
                    counter <- counter + 1 # increment counter
                  } # end j loop
                } # end i loop

                # Input unique counters into unique pop, age, year, season, and sex blocks
                for(p in map_p) for(a in map_a) for(y in map_y) for(seas in map_seas) for(s in map_s) map_Movement_Pars[p,,,y,seas,a,s] <- map_idx

              } # end sex block
            } # end season block
          } # end year block
        } # end age block
      } # end pop block

    } else map_Movement_Pars <- factor(rep(NA, length(input_list$par$move_pars))) # don't estimate movement

    # Turn off parameters for CTMC
    map_log_move_diffusion_pars <- factor(rep(NA, length(map_log_move_diffusion_pars)))
    map_move_preference_pars <- factor(rep(NA, length(map_move_preference_pars)))
  }

  # CTMC movement
  if(input_list$data$move_type == 1) {
    # turn off parameters for unstructured markov
    map_Movement_Pars <- factor(rep(NA, length(input_list$par$move_pars))) # don't estimate movement
    # estimate parameters for CTMC
    if(length(map_log_move_diffusion_pars) != 0 ) map_log_move_diffusion_pars <- factor(1:length(map_log_move_diffusion_pars))
    if(length(map_move_preference_pars) != 0 ) map_move_preference_pars <- factor(1:length(map_move_preference_pars))
  }

  # Input into mapping list
  input_list$map$move_pars <- factor(map_Movement_Pars)
  input_list$map$log_move_diffusion_pars <- factor(map_log_move_diffusion_pars)
  input_list$map$move_preference_pars <- factor(map_move_preference_pars)

  return(input_list)
}

#' Construct Mapping for Continuous Movement Deviations and Process Error
#'
#' Internal helper used by \code{Setup_Mod_Movement()} to define the TMB
#' parameter mapping for continuous time-varying movement deviations and
#' associated process-error parameters.
#'
#' @param input_list Model input list containing \code{$data}, \code{$par},
#'   and \code{$map}. Modified in-place and returned.
#' @param cont_vary_movement Character string specifying the structure of
#'   independent deviations (e.g., \code{"iid_y"}, \code{"iid_y_a_s"},
#'   \code{"iid_p_y_seas_a_s"}).
#' @param Movement_cont_pe_pars_spec Character string specifying how movement
#'   process-error parameters are treated. Options:
#'   \itemize{
#'     \item \code{"none"}: no process-error parameters.
#'     \item \code{"fix"}: parameters fixed (not estimated).
#'     \item \code{"est_shared"}: single shared parameter estimated.
#'     \item \code{"est_all"}: all process-error parameters estimated independently.
#'   }
#'
#' @details
#' When continuous movement variation is enabled
#' (\code{cont_vary_movement != "none"}) and movement is estimated,
#' this function assigns unique factor levels to movement deviation
#' parameters according to the specified structure.
#'
#' Supported structures include:
#' \itemize{
#'   \item Year-only (\code{"iid_y"})
#'   \item Age-only (\code{"iid_a"})
#'   \item Year × Age (\code{"iid_y_a"})
#'   \item Year × Age × Sex (\code{"iid_y_a_s"})
#'   \item Year × Season × Age × Sex (\code{"iid_y_seas_a_s"})
#'   \item Population-interacted analogues (e.g., \code{"iid_p_y_a_s"})
#' }
#'
#' For CTMC movement, deviations are assigned only to region pairs that are
#' connected in the adjacency matrix.
#'
#' @return The updated \code{input_list} with \code{$map$move_devs} and
#'   \code{$map$move_pe_pars} defined.
#'
#' @keywords internal
do_cont_vary_move_mapping <- function(input_list, cont_vary_movement, Movement_cont_pe_pars_spec) {

  # Setup mapping list
  n_regions_to <- dim(input_list$par$move_devs)[3] # get movement to
  map_move_devs <- array(NA, dim = dim(input_list$par$move_devs))
  map_move_pe_pars <- array(NA, dim = dim(input_list$par$move_pe_pars))

  # Whether or not recruits move
  age_start <- ifelse(input_list$data$do_recruits_move == 0 && length(input_list$data$ages) >= 2, 2, 1)

  # Movement Deviations -----------------------------------------------
  if(input_list$data$n_regions > 1 && # if spatial model
     input_list$data$cont_vary_movement > 0 && # if continuous varying movement
     input_list$data$use_fixed_movement == 0 # if not using fixed movement matrix
  ) {

    # Dimensions
    n_pop <- input_list$data$n_pop
    n_yrs_devs <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs
    n_seas <- input_list$data$n_seas
    n_ages <- length(input_list$data$ages)
    n_sexes <- input_list$data$n_sexes

    counter <- 1 # setup counter

    for(r in 1:input_list$data$n_regions) {
      for(rr in 1:n_regions_to) {

        # if regions are not adjacent or residency (CTMC), input NA and restart loop
        if(input_list$data$move_type == 1 & input_list$data$adjacency_collapsed[r,rr] == 0) {
          map_move_devs[,r,rr,,,,] <- NA
          next
        }

        # --- Year only ---
        if(cont_vary_movement %in% c('iid_y')) {
          for(y in 1:n_yrs_devs) {
            map_move_devs[,r,rr,y,,,] <- counter
            counter <- counter + 1
          } # end y loop
        } # end if iid_y

        # --- Age only ---
        if(cont_vary_movement %in% c('iid_a')) {
          for(a in age_start:n_ages) {
            map_move_devs[,r,rr,,,a,] <- counter
            counter <- counter + 1
          } # end a loop
        } # end if iid_a

        # --- Year x Age ---
        if(cont_vary_movement %in% c('iid_y_a')) {
          for(y in 1:n_yrs_devs) {
            for(a in age_start:n_ages) {
              map_move_devs[,r,rr,y,,a,] <- counter
              counter <- counter + 1
            } # end a loop
          } # end y loop
        } # end if iid_y_a

        # --- Year x Age x Sex ---
        if(cont_vary_movement %in% c('iid_y_a_s')) {
          for(y in 1:n_yrs_devs) {
            for(a in age_start:n_ages) {
              for(s in 1:n_sexes) {
                map_move_devs[,r,rr,y,,a,s] <- counter
                counter <- counter + 1
              } # end s loop
            } # end a loop
          } # end y loop
        } # end if iid_y_a_s

        # --- Year x Season x Age x Sex ---
        if(cont_vary_movement %in% c('iid_y_seas_a_s')) {
          for(y in 1:n_yrs_devs) {
            for(seas in 1:n_seas) {
              for(a in age_start:n_ages) {
                for(s in 1:n_sexes) {
                  map_move_devs[,r,rr,y,seas,a,s] <- counter
                  counter <- counter + 1
                } # end s loop
              } # end a loop
            } # end seas loop
          } # end y loop
        } # end if iid_y_seas_a_s

        # --- pop x Year only ---
        if(cont_vary_movement %in% c('iid_p_y')) {
          for(p in 1:n_pop) {
            for(y in 1:n_yrs_devs) {
              map_move_devs[p,r,rr,y,,,] <- counter
              counter <- counter + 1
            } # end y loop
          } # end p loop
        } # end if iid_p_y

        # --- pop x Age only ---
        if(cont_vary_movement %in% c('iid_p_a')) {
          for(p in 1:n_pop) {
            for(a in age_start:n_ages) {
              map_move_devs[p,r,rr,,,a,] <- counter
              counter <- counter + 1
            } # end a loop
          } # end p loop
        } # end iid_p_a

        # --- pop x Year x Age ---
        if(cont_vary_movement %in% c('iid_p_y_a')) {
          for(p in 1:n_pop) {
            for(y in 1:n_yrs_devs) {
              for(a in age_start:n_ages) {
                map_move_devs[p,r,rr,y,,a,] <- counter
                counter <- counter + 1
              } # end a loop
            } # end y loop
          } # end p loop
        } # end iid_p_y_a

        # --- popp x Year x Age x Sex ---
        if(cont_vary_movement %in% c('iid_p_y_a_s')) {
          for(p in 1:n_pop) {
            for(y in 1:n_yrs_devs) {
              for(a in age_start:n_ages) {
                for(s in 1:n_sexes) {
                  map_move_devs[p,r,rr,y,,a,s] <- counter
                  counter <- counter + 1
                } # end s loop
              } # end a loop
            } # end y loop
          } # end p loop
        } # end iid_p_y_a_s

        # --- pop x Year x Season x Age x Sex ---
        if(cont_vary_movement %in% c('iid_p_y_seas_a_s')) {
          for(p in 1:n_pop) {
            for(y in 1:n_yrs_devs) {
              for(seas in 1:n_seas) {
                for(a in age_start:n_ages) {
                  for(s in 1:n_sexes) {
                    map_move_devs[p,r,rr,y,seas,a,s] <- counter
                    counter <- counter + 1
                  } # end s loop
                } # end a loop
              } # end seas loop
            } # end y loop
          } # end p loop
        } # end iid_p_y_seas_a_s

      } # end rr
    } # end r loop
  }

    # Movement Process Error Parameters ---------------------------------------

    # Mapping for movement process error deviations
    if(Movement_cont_pe_pars_spec %in% c("fix", "none")) map_move_pe_pars <- map_move_pe_pars
    if(Movement_cont_pe_pars_spec == 'est_all') map_move_pe_pars[] <- 1:length(map_move_pe_pars)
    if(Movement_cont_pe_pars_spec == 'est_shared') map_move_pe_pars[] <- 1

    # return to input list
    input_list$map$move_devs <- factor(map_move_devs)
    input_list$data$map_move_devs <- array(as.numeric(input_list$map$move_devs), dim = dim(input_list$par$move_devs))
    input_list$map$move_pe_pars <- factor(map_move_pe_pars)
    return(input_list)

}

#' Configure Movement Processes for SPoRC
#'
#' Define and initialize movement model components for SPoRC, including
#' unstructured Markov movement (\code{move_type = 0}) and Continuous Time
#' Markov Chain (CTMC) movement (\code{move_type = 1}). This function
#' validates inputs, initializes parameter arrays, constructs TMB mapping
#' structures, and updates the model input list accordingly.
#'
#' @param input_list A list containing \code{$data}, \code{$par}, and
#'   \code{$map} components. This object is updated and returned.
#'
#' @param move_type Integer specifying the movement model:
#'   \itemize{
#'     \item \code{0}: Unstructured Markov transition model (multinomial logit parameterization).
#'     \item \code{1}: Continuous Time Markov Chain (CTMC).
#'   }
#'
#' @param do_recruits_move Integer (0/1) indicating whether recruits
#'   participate in movement.
#'
#' @param use_fixed_movement Integer (0/1). If \code{1}, movement rates are
#'   supplied via \code{Fixed_Movement} and are not estimated.
#'
#' @param Fixed_Movement Numeric array of dimension
#'   \code{[n_pop, n_regions, n_regions, n_years, n_seas, n_ages, n_sexes]}.
#'   Required when \code{use_fixed_movement = 1}. If \code{NA}, a neutral
#'   array of ones is created.
#'
#' @param Use_Movement_Prior Integer (0/1) indicating whether Dirichlet
#'   priors are applied to movement probabilities.
#'
#' @param Movement_prior Optional data.frame defining Dirichlet prior
#'   concentration parameters. Must include columns:
#'   \code{pop}, \code{region_from}, \code{year}, \code{seas},
#'   \code{age}, \code{sex}, and \code{alpha}, where \code{alpha}
#'   is a list-column of length-\code{n_regions} numeric vectors.
#'
#' @param Movement_popblk_spec,Movement_ageblk_spec,Movement_yearblk_spec,
#'   Movement_seasblk_spec,Movement_sexblk_spec Block specifications for
#'   parameter sharing under \code{move_type = 0}. Each may be either
#'   \code{"constant"} or a \code{list} of integer vectors defining blocks.
#'
#' @param cont_vary_movement Character string specifying the structure of
#'   independent movement deviations. Options include:
#'   \code{"none"}, \code{"iid_y"}, \code{"iid_a"},
#'   \code{"iid_y_a"}, \code{"iid_y_a_s"},
#'   \code{"iid_y_seas_a_s"}, and population-interacted variants (e.g., \code{"iid_p_y_a}).
#'
#' @param Movement_cont_pe_pars_spec Character string controlling estimation
#'   of process-error parameters for continuous deviations:
#'   \code{"none"}, \code{"fix"}, \code{"est_shared"}, or \code{"est_all"}.
#'
#' @param ctmc_move_dat Data.frame required when \code{move_type = 1}.
#'   Must include columns:
#'   \code{pop}, \code{regions}, \code{years}, \code{seas},
#'   \code{ages}, \code{sexes}, plus any covariates referenced in
#'   \code{diffusion_formula} or \code{preference_formula}.
#'
#' @param adjacency_mat Square \code{n_regions × n_regions} matrix defining
#'   allowable transitions for CTMC movement.
#'
#' @param area_r Numeric vector of region areas (length \code{n_regions}),
#'   required for CTMC movement.
#'
#' @param diffusion_formula R formula defining the linear predictor for
#'   CTMC diffusion parameters.
#'
#' @param preference_formula R formula defining the linear predictor for
#'   CTMC preference parameters.
#'
#' @param ctmc_diffusion_bounds Integer (0/1) indicating whether bounds are
#'   applied to ensure a Metzler generator matrix.
#'
#' @param ... Optional named starting values for parameters (e.g.,
#'   \code{move_pars}, \code{move_devs}, \code{move_pe_pars},
#'   \code{log_move_diffusion_pars}, \code{move_preference_pars}).
#'
#' @details
#' For \strong{unstructured Markov movement}, region-to-region transition
#' probabilities are parameterized independently (subject to multinomial
#' constraints). Block specifications determine parameter sharing across
#' dimensions.
#'
#' For \strong{CTMC movement}, diffusion and preference components are
#' constructed from user-specified design matrices derived from the supplied
#' formulas. The adjacency matrix restricts allowable transitions.
#'
#' Continuous movement deviations (if enabled) introduce independent
#' random effects structured according to \code{cont_vary_movement}.
#'
#' @export Setup_Mod_Movement
#' @family Model Setup
Setup_Mod_Movement <- function(input_list,
                               move_type = 0,
                               do_recruits_move = 0,
                               use_fixed_movement = 0,
                               Fixed_Movement = NA,
                               Use_Movement_Prior = 0,
                               Movement_prior = NULL,
                               Movement_popblk_spec = 'constant',
                               Movement_ageblk_spec = 'constant',
                               Movement_yearblk_spec = 'constant',
                               Movement_seasblk_spec = 'constant',
                               Movement_sexblk_spec = 'constant',
                               cont_vary_movement = 'none',
                               Movement_cont_pe_pars_spec = 'none',
                               ctmc_move_dat = NULL,
                               adjacency_mat = NULL,
                               area_r = rep(1, input_list$data$n_regions),
                               diffusion_formula = NULL,
                               preference_formula = NULL,
                               ctmc_diffusion_bounds = 0,
                               ...
) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...) # get starting values if there are any

  # Input Validation --------------------------------------------------------

  # If no movement matrix is provided
  if(is.na(sum(Fixed_Movement))) {
    Fixed_Movement <- array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_regions,
                                       length(input_list$data$years), input_list$data$n_seas,
                                       length(input_list$data$ages), input_list$data$n_sexes))
  }

  # Check fixed movement matrix
  if(!use_fixed_movement %in% c(0,1)) stop('Options for fixing movement are not correctly specified. The options are use_fixed_movement == 0 (dont use and estiamte movement parameters), or == 1 (use)')
  else collect_message("Movement is: ", ifelse(use_fixed_movement == 0, "Estimated", "Fixed"))
  if(use_fixed_movement == 1) check_data_dimensions(Fixed_Movement, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_seas = input_list$data$n_seas, what = 'Fixed_Movement')

  # Check for movement priors
  if(!Use_Movement_Prior %in% c(0,1)) stop('Options for movement priors not correctly specified. The options are Use_Movement_Prior == 0 (dont use), or == 1 (use)')
  else collect_message("Movement priors are: ", ifelse(Use_Movement_Prior == 0, "Not Used", "Used"))

  # Check for recruits moving
  if(!do_recruits_move %in% c(0,1)) stop('Movement for recruits is not correctly specified. The options are do_recruits_move == 0 (they dont move), or == 1 (they move)')
  else collect_message("Recruits are: ", ifelse(do_recruits_move == 0, "Not Moving", "Moving"))

  # Check movement continuous varying parameterization
  if(!cont_vary_movement %in% c("none", "iid_y", "iid_a", "iid_y_a", "iid_y_a_s", "iid_y_seas_a_s",
                                "iid_p_y", "iid_p_a", "iid_p_y_a", "iid_p_y_a_s", "iid_p_y_seas_a_s" ))
    stop('Options for continuous movement is not correctly specified. The options are none,
         iid_y, iid_a, iid_y_a, iid_y_a_s, iid_y_seas_a_s, iid_p_y, iid_p_a, iid_p_y_a, iid_p_y_a_s, iid_p_y_seas_a_s')
  else collect_message("Continuous movement specification is: ", cont_vary_movement)

  # Check movement process error estimation (no change needed here)
  if(!Movement_cont_pe_pars_spec %in% c('none', 'fix', 'est_all', 'est_shared'))
    stop('Options for continuous movement process error is not correctly specified.')
  else collect_message("Continuous movement process error specification is: ", Movement_cont_pe_pars_spec)

  if(!move_type %in% c(0, 1)) stop('move_type must be 0 (unstructured) or 1 (Continuous Time Markov Chain)')
  collect_message("Movement type is: ", ifelse(move_type == 0, "Unstructured Markov", "Continuous Time Markov Chain"))

  # Check movement blocks (for unstructured markov)
  if(move_type == 0) {
    if(!is.null(Movement_popblk_spec)) if(!typeof(Movement_popblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects population blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 years and wanted 2 population blocks, this would be list(1,2).")
    if(!is.null(Movement_ageblk_spec)) if(!typeof(Movement_ageblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects age blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 ages and wanted 2 age blocks, this would be list(c(1:5), c(6:10)) such that ages 1 - 5 are a block, and ages 6 - 10 are a block.")
    if(!is.null(Movement_yearblk_spec)) if(!typeof(Movement_yearblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects year blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 years and wanted 2 year blocks, this would be list(c(1:5), c(6:10)) such that years 1 - 5 are a block, and years 6 - 10 are a block.")
    if(!is.null(Movement_sexblk_spec)) if(!typeof(Movement_sexblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects sex blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 2 sexes and wanted sex-specific movement, this would be list(1, 2).")
    if(!is.null(Movement_seasblk_spec)) if(!typeof(Movement_seasblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects season blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 4 seasons and wanted 2 season blocks, this would be list(c(1:2), c(3:4)) such that seasons 1 - 2 are a block, and seasons 3 - 4 are a block.")
    if(is.list(Movement_popblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_popblk_spec), " population blocks") else collect_message("Movement fixed effect blocks are population-invariant")
    if(is.list(Movement_seasblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_seasblk_spec), " seas blocks") else collect_message("Movement fixed effect blocks are season-invariant")
    if(is.list(Movement_sexblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_sexblk_spec), " sex blocks") else collect_message("Movement fixed effect blocks are sex-invariant")
    if(is.list(Movement_yearblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_yearblk_spec), " year blocks") else collect_message("Movement fixed effect blocks are time-invariant")
    if(is.list(Movement_ageblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_ageblk_spec), " age blocks") else collect_message("Movement fixed effect blocks are age-invariant")
    # create fully connected adjacency matrix
    adjacency_mat <- base::matrix(1, nrow = input_list$data$n_regions, ncol = input_list$data$n_regions)
    diag(adjacency_mat) <- 0
  }

  # Check CTMC movement
  if(move_type == 1) {

    # Make sure blocks are not specified
    if ((Movement_popblk_spec != "constant") ||
        (Movement_ageblk_spec != "constant") ||
        (Movement_yearblk_spec != "constant") ||
        (Movement_sexblk_spec != "constant") ||
        (Movement_seasblk_spec != "constant")) {
      stop("Movement blocks (pop, age, year, seas, or sex) must be NULL or 'constant' when CTMC movement is used.")
    }

    # check adjacency matrix
    if(is.null(adjacency_mat)) stop("adjacency_mat is required for CTMC movement (move_type = 1)")
    if(nrow(adjacency_mat) != input_list$data$n_regions ||
       ncol(adjacency_mat) != input_list$data$n_regions) {
      stop("adjacency_mat must be a square matrix with dimensions n_regions x n_regions")
    }

    # check area sizes
    if(is.null(area_r)) stop("area_r is required for CTMC movement (move_type = 1)")
    if(length(area_r) != input_list$data$n_regions) stop("area_r must have length n_regions")

    # check ctmc data frame
    required_cols <- c("pop", "regions", "years", "seas", "ages", "sexes")
    if(!all(required_cols %in% names(ctmc_move_dat))) {
      missing <- setdiff(required_cols, names(ctmc_move_dat))
      stop("ctmc_move_dat must have columns: ", paste(required_cols, collapse = ", "),
           "\n  Missing: ", paste(missing, collapse = ", "))
    }

    # Extract variables from formulas and check they exist in ctmc_move_dat
    diffusion_vars <- all.vars(diffusion_formula)
    preference_vars <- all.vars(preference_formula)

    # Check formulas
    if(is.null(diffusion_formula)) stop("diffusion_formula is required for CTMC movement")
    if(is.null(preference_formula)) stop("preference_formula is required for CTMC movement")

    # Check diffusion formula variables
    missing_diff <- setdiff(diffusion_vars, names(ctmc_move_dat))
    if(length(missing_diff) > 0) {
      stop("Variables in diffusion_formula not found in ctmc_move_dat:\n",
           "  Missing: ", paste(missing_diff, collapse = ", "), "\n",
           "  Available: ", paste(names(ctmc_move_dat), collapse = ", "))
    }

    # Check preference formula variables
    missing_pref <- setdiff(preference_vars, names(ctmc_move_dat))
    if(length(missing_pref) > 0) {
      stop("Variables in preference_formula not found in ctmc_move_dat:\n",
           "  Missing: ", paste(missing_pref, collapse = ", "), "\n",
           "  Available: ", paste(names(ctmc_move_dat), collapse = ", "))
    }

    # Coerce CTMC years > n_yrs to equal n_yrs if that is the case; makes sure we are not extrapolating splines
    # while allowing for covariate projections
    proj_year_idx <- which(ctmc_move_dat$years > length(input_list$data$years))
    if(length(proj_year_idx) > 0) {
      ctmc_move_dat$years[proj_year_idx] <- length(input_list$data$years)
      collect_message("ctmc_move_dat has years > n_yrs for projections. These years are capped at n_yrs to prevent spline extrapolation, while allowing for covariate projections.")
    }
  }

  # check movement prior
  if(!is.null(Movement_prior)) {
    required_cols <- c("pop", "region_from", 'year', 'seas', "age", "sex", "alpha")
    missing_cols <- setdiff(required_cols, names(Movement_prior))
    if(length(missing_cols) > 0) stop("Movement_prior is missing required columns: ", paste(missing_cols, collapse = ", "))

    # check dimensions for alpha
    for(i in 1:nrow(Movement_prior)) {
      alpha_vec <- Movement_prior$alpha[[i]]
      if(length(alpha_vec) != input_list$data$n_regions) stop("Row ", i, ": alpha vector has length ", length(alpha_vec), " but should have length ", input_list$data$n_regions)
    } # end i loop
  }

  # make collapsed adjacency matrix
  adjacency_collapsed = base::matrix(NA, nrow = input_list$data$n_regions, ncol = input_list$data$n_regions - 1) # get collapsed adjacency matrix
  # create collapsed adjacency matrix for indexing devs that should be penalized
  for(r in 1:input_list$data$n_regions) {
    counter_col <- 1
    for(rr_full in 1:input_list$data$n_regions) {
      if(r != rr_full) {  # Skip diagonal
        adjacency_collapsed[r, counter_col] = as.numeric(adjacency_mat[r, rr_full])
        counter_col = counter_col + 1
      } # end if
    } # end rr_full
  } # end r loop

  # Populate Data List ------------------------------------------------------
  input_list$data$move_type <- move_type
  input_list$data$do_recruits_move <- do_recruits_move
  input_list$data$use_fixed_movement <- use_fixed_movement
  input_list$data$Fixed_Movement <- Fixed_Movement
  input_list$data$Use_Movement_Prior <- Use_Movement_Prior
  input_list$data$Movement_prior <- Movement_prior

  # define things for CTMC movement
  input_list$data$adjacency_mat <- adjacency_mat
  input_list$data$adjacency_collapsed <- adjacency_collapsed
  input_list$data$area_r <- area_r
  input_list$data$ctmc_move_dat <- ctmc_move_dat
  input_list$data$diffusion_formula <- diffusion_formula
  input_list$data$preference_formula <- preference_formula
  input_list$data$ctmc_diffusion_bounds <- ctmc_diffusion_bounds

  # define for continuous varying movement
  cont_move_map <- data.frame(
    type = c("none", "iid_y", "iid_a", "iid_y_a", "iid_y_a_s", "iid_y_seas_a_s", "iid_p_y", "iid_p_a", "iid_p_y_a", "iid_p_y_a_s", "iid_p_y_seas_a_s"),
    num = 0:10
  )
  cont_vary_movement_val <- cont_move_map$num[cont_move_map$type == cont_vary_movement] # look for number corresponding to specified option
  input_list$data$cont_vary_movement <- cont_vary_movement_val

  # Populate Parameter List -------------------------------------------------

  # Movement Parameters (for unstructured markov; move_type == 0)
  if("move_pars" %in% names(starting_values)) input_list$par$move_pars <- starting_values$move_pars
  else input_list$par$move_pars <- array(0, dim = c(input_list$data$n_pop,
                                                    input_list$data$n_regions, input_list$data$n_regions - 1,
                                                    length(input_list$data$years), input_list$data$n_seas,
                                                    length(input_list$data$ages), input_list$data$n_sexes))

  # Movement Parameters (for CTMTC; move_type == 1)
  # get design matrix to figure out number of parameters needed
  if(move_type == 0) n_gamma <- n_theta <- 1 # if unstructered markov, then use 1 as place holder
  if(move_type == 1) {
    if(do_recruits_move == 0) {
      recruit_idx <- which(input_list$data$ctmc_move_dat$ages == min(input_list$data$ctmc_move_dat$ages))
      if(length(recruit_idx) == 1) input_list$data$ctmc_move_dat <- input_list$data$ctmc_move_dat[-recruit_idx,] # remove recruits
    }
    designs = get_movement_dp_design_matrix(data = input_list$data$ctmc_move_dat,
                                            preference_formula = input_list$data$preference_formula,
                                            diffusion_formula = input_list$data$diffusion_formula)
    n_theta <- designs$n_theta # extract out number of pars
    n_gamma <- designs$n_gamma # extract out number of pars
  }

  # diffusion parameters
  if("log_move_diffusion_pars" %in% names(starting_values)) input_list$par$log_move_diffusion_pars <- starting_values$log_move_diffusion_pars
  else input_list$par$log_move_diffusion_pars <- rep(log(0.1), n_theta)

  # preference parameters
  if("move_preference_pars" %in% names(starting_values)) input_list$par$move_preference_pars <- starting_values$move_preference_pars
  else input_list$par$move_preference_pars <- rep(0, n_gamma)

  # Movement deviations
  if("move_devs" %in% names(starting_values)) input_list$par$move_devs <- starting_values$move_devs
  else {
    input_list$par$move_devs <- array(0, c(input_list$data$n_pop,
                                           input_list$data$n_regions, input_list$data$n_regions - 1,
                                           length(input_list$data$years) + input_list$data$n_proj_yrs_devs,
                                           input_list$data$n_seas,
                                           length(input_list$data$ages),
                                           input_list$data$n_sexes))
  }

  # Movement process error parameters
  if("move_pe_pars" %in% names(starting_values)) input_list$par$move_pe_pars <- starting_values$move_pe_pars
  else input_list$par$move_pe_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                       input_list$data$n_seas, length(input_list$data$ages),
                                                       input_list$data$n_sexes)) # max 4 parameters or the ages


  # Mapping Options ---------------------------------------------------------
  input_list <- do_move_pars_mapping(input_list, Movement_popblk_spec, Movement_ageblk_spec, Movement_yearblk_spec, Movement_sexblk_spec, Movement_seasblk_spec, use_fixed_movement)
  input_list <- do_cont_vary_move_mapping(input_list, cont_vary_movement, Movement_cont_pe_pars_spec)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
