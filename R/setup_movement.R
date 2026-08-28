# Stage 1 of 3: model setup
#
# Movement inputs: how fish redistribute among regions, and when in the season
# that happens. Chooses between unstructured transition matrices (move_type 0)
# and a continuous time generator (move_type 1), sets move_timing, and builds the
# movement parameter map. Borrows get_movement_dp_design_matrix from
# model_movement.R so the setup and the objective agree on the design matrix.

#' Map unstructured Markov movement parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Movement}} to construct the
#' TMB/RTMB factor maps for \code{move_pars} (unstructured Markov transitions),
#' \code{log_move_diffusion_pars} (CTMC diffusion), and
#' \code{move_preference_pars} (CTMC taxis). Under the unstructured Markov
#' formulation (\code{move_type = 0}), parameters within each block combination
#' share a common estimation index; for CTMC (\code{move_type = 1}),
#' \code{move_pars} is mapped entirely to \code{NA} and only the CTMC-specific
#' arrays are activated.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param Movement_popblk_spec \code{"constant"} for population-invariant
#'   movement, or a list of integer vectors partitioning populations into blocks
#'   that share parameters.
#' @param Movement_ageblk_spec \code{"constant"} for age-invariant movement, or
#'   a list of integer vectors defining age blocks.
#' @param Movement_yearblk_spec \code{"constant"} for time-invariant movement,
#'   or a list of integer vectors defining year blocks.
#' @param Movement_sexblk_spec \code{"constant"} for sex-invariant movement, or
#'   a list of integer vectors defining sex blocks.
#' @param Movement_seasblk_spec \code{"constant"} for season-invariant movement,
#'   or a list of integer vectors defining season blocks.
#' @param use_fixed_movement Integer flag. \code{1} = movement rates are
#'   externally supplied; all \code{move_pars} are mapped to \code{NA} and not
#'   estimated. \code{0} = movement is estimated.
#'
#' @return The input \code{input_list} with three \code{$map} entries updated:
#'   \describe{
#'     \item{\code{$map$move_pars}}{Factor vector for unstructured Markov
#'       transition parameters. Under \code{move_type = 0} with estimated
#'       movement, cells within the same block receive the same integer index;
#'       cells outside a spatial model or with fixed movement are \code{NA}.
#'       Entirely \code{NA} under \code{move_type = 1}.}
#'     \item{\code{$map$log_move_diffusion_pars}}{Factor vector for CTMC
#'       diffusion parameters. Active (sequential integers) under
#'       \code{move_type = 1}; entirely \code{NA} under \code{move_type = 0}.}
#'     \item{\code{$map$move_preference_pars}}{Factor vector for CTMC
#'       preference (taxis) parameters. Active under \code{move_type = 1};
#'       entirely \code{NA} under \code{move_type = 0}.}
#'   }
#'
#' @keywords internal
do_move_pars_mapping <- function(input_list, Movement_popblk_spec,
                                 Movement_ageblk_spec, Movement_yearblk_spec,
                                 Movement_sexblk_spec, Movement_seasblk_spec, use_fixed_movement) {

  # Setup mapping list
  map_Movement_Pars <- input_list$par$move_pars # initialize array with same dimensions as parameters
  map_Movement_Pars[] <- NA # any cell no block covers stays fixed rather than sharing a level
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
      if(Movement_ageblk_spec == "constant") Movement_ageblk_spec_vals <- list(seq_along(input_list$data$ages))
    } else Movement_ageblk_spec_vals <- Movement_ageblk_spec

    # If movement is constant across years
    if(is.character(Movement_yearblk_spec)){
      if(Movement_yearblk_spec == "constant") Movement_yearblk_spec_vals = list(seq_along(input_list$data$years))
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

#' Map continuous movement deviation and process-error parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Movement}} to construct the
#' TMB/RTMB factor maps for \code{move_devs} (iid deviations on the movement
#' logit or log-rate surface) and \code{move_pe_pars} (process-error variance
#' parameters). Deviations are only activated when the model is spatial
#' (\code{n_regions > 1}), continuous variation is requested
#' (\code{cont_vary_movement > 0}), and movement is estimated
#' (\code{use_fixed_movement == 0}). For CTMC movement, deviations are only
#' assigned to region pairs that are connected in the adjacency matrix. The
#' resulting integer map is also stored as \code{$data$map_move_devs} for use
#' in the C++ template.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param cont_vary_movement Character string specifying the deviation structure.
#'   One of \code{"none"}, \code{"iid_y"}, \code{"iid_a"}, \code{"iid_y_a"},
#'   \code{"iid_y_a_s"}, \code{"iid_y_seas_a_s"}, or the population-specific
#'   analogues \code{"iid_p_y"}, \code{"iid_p_a"}, \code{"iid_p_y_a"},
#'   \code{"iid_p_y_a_s"}, \code{"iid_p_y_seas_a_s"}. Dimensions present in
#'   the string receive unique estimation indices; absent dimensions share a
#'   single index. \code{"none"} maps all deviations to \code{NA}.
#' @param Movement_cont_pe_pars_spec Character string controlling estimation of
#'   the process-error variance for movement deviations. One of:
#'   \describe{
#'     \item{\code{"none"} or \code{"fix"}}{All \code{move_pe_pars} held fixed
#'       (mapped to \code{NA} or at starting values).}
#'     \item{\code{"est_shared"}}{Single variance parameter shared across all
#'       dimensions (all elements mapped to index 1).}
#'     \item{\code{"est_all"}}{All \code{move_pe_pars} estimated independently
#'       with dimensions
#'       \code{[n_pop × n_regions × n_seas × n_ages × n_sexes]}.}
#'   }
#'
#' @return The input \code{input_list} with three entries updated:
#'   \describe{
#'     \item{\code{$map$move_devs}}{Factor vector for movement deviations.
#'       Active cells receive sequential integer indices; non-adjacent CTMC
#'       pairs and inactive configurations are \code{NA}.}
#'     \item{\code{$data$map_move_devs}}{Integer array (same dimensions as
#'       \code{$par$move_devs}) storing the numeric version of the factor map
#'       for use in the C++ objective function.}
#'     \item{\code{$map$move_pe_pars}}{Factor vector for process-error
#'       variance parameters, following \code{Movement_cont_pe_pars_spec}.}
#'   }
#'
#' @keywords internal
do_cont_vary_move_mapping <- function(input_list, cont_vary_movement, Movement_cont_pe_pars_spec) {

  # Setup mapping list
  n_regions_to <- dim(input_list$par$move_devs)[3] # get movement to
  map_move_devs <- array(NA, dim = dim(input_list$par$move_devs))
  map_move_pe_pars <- array(NA, dim = dim(input_list$par$move_pe_pars))

  # Movement Deviations -----------------------------------------------
  if(input_list$data$n_regions > 1 && # if spatial model
     input_list$data$cont_vary_movement > 0 && # if continuous varying movement
     input_list$data$use_fixed_movement == 0 # if not using fixed movement matrix
  ) {

    # Dimensions (region_from/region_to are always "key": every region pair
    # gets its own set of estimated deviations, never shared with another pair)
    dims <- c(pop = input_list$data$n_pop,
              region_from = input_list$data$n_regions,
              region_to = n_regions_to,
              year = length(input_list$data$years) + input_list$data$n_proj_yrs_devs,
              season = input_list$data$n_seas,
              age = length(input_list$data$ages),
              sex = input_list$data$n_sexes)

    # dims named in the spec (e.g. "iid_y_a_s" -> year, age, sex) are "key"
    # dims (unique value per combination); dims not named are shared/broadcast
    dim_abbrev <- c(p = "pop", y = "year", seas = "season", a = "age", s = "sex")
    key_extra <- unname(dim_abbrev[strsplit(sub("^iid_", "", cont_vary_movement), "_")[[1]]])
    share_over <- setdiff(dim_abbrev, key_extra)

    map_move_devs <- build_pe_map(dims, share_over = share_over)

    # Whether or not recruits (age 1) move: only relevant when age is itself
    # a key dim of the spec, when age is a broadcast dim, age 1 shares the
    # same tied deviation as every other age instead of being masked out.
    if("age" %in% key_extra && input_list$data$do_recruits_move == 0 && dims["age"] >= 2) {
      map_move_devs[,,,,,1,] <- NA
    }

    # if regions are not adjacent or residency (CTMC), no deviation is estimated
    if(input_list$data$move_type == 1) {
      for(r in 1:input_list$data$n_regions) {
        for(rr in 1:n_regions_to) {
          if(input_list$data$adjacency_collapsed[r,rr] == 0) map_move_devs[,r,rr,,,,] <- NA
        } # end rr
      } # end r
    }
  }

    # Movement Process Error Parameters ---------------------------------------

    # Mapping for movement process error deviations
    if(Movement_cont_pe_pars_spec %in% c("fix", "none")) map_move_pe_pars <- map_move_pe_pars
    if(Movement_cont_pe_pars_spec == 'est_all') map_move_pe_pars[] <- 1:length(map_move_pe_pars)
    if(Movement_cont_pe_pars_spec == 'est_shared') map_move_pe_pars[] <- 1

    # return to input list
    input_list$map$move_devs <- factor(as.vector(map_move_devs))
    input_list$data$map_move_devs <- array(as.numeric(input_list$map$move_devs), dim = dim(input_list$par$move_devs))
    input_list$map$move_pe_pars <- factor(map_move_pe_pars)
    return(input_list)

}

#' Set up movement model inputs and parameter structures
#'
#' Configures all aspects of spatial movement for the estimation model,
#' supporting both unstructured Markov transition (\code{move_type = 0}) and
#' Continuous Time Markov Chain (\code{move_type = 1}) formulations, with
#' optional continuous iid deviations on the movement surface. Validates all
#' inputs, initialises parameter arrays, constructs TMB/RTMB factor maps, and
#' populates \code{input_list$data} accordingly. Must be called after
#' \code{\link{Setup_Mod_Biologicals}}.
#'
#' @section Unstructured Markov movement (\code{move_type = 0}):
#' Transition probabilities from region \eqn{r} to all other regions are
#' parameterised via a multinomial logit with a reference-cell constraint.
#' The parameter array \code{move_pars} has dimensions
#' \code{[n_pop × n_regions × (n_regions - 1) × n_years × n_seas × n_ages × n_sexes]}.
#' Block specifications (\code{Movement_*blk_spec}) control sharing: indices
#' within the same block receive the same TMB factor level and are estimated
#' as a single free parameter. A fully connected adjacency matrix is
#' constructed automatically. Blocked and continuous time-varying movement
#' can be combined: use \code{Movement_yearblk_spec} for discrete structural
#' breaks and \code{cont_vary_movement} for residual year-to-year variation.
#'
#' @section CTMC movement (\code{move_type = 1}):
#' The instantaneous rate matrix \eqn{Q} is decomposed into diffusion
#' (\eqn{\theta}) and preference (\eqn{\gamma}) components following
#' Thorson et al. Design matrices for both components are derived from
#' \code{diffusion_formula} and \code{preference_formula} evaluated on
#' \code{ctmc_move_dat}. The discrete-time movement matrix for each time step
#' is \eqn{\exp(Q \Delta t)}. Parameter blocking is not supported for CTMC;
#' all \code{Movement_*blk_spec} arguments must remain \code{"constant"}.
#' Structural variation across populations, ages, sexes, or seasons should
#' instead be introduced through formula covariates in \code{ctmc_move_dat}.
#'
#' @section Continuous movement deviations:
#' IID deviations (\code{move_devs}) are added to the movement logit surface
#' (unstructured Markov) or log-rate surface (CTMC) before computing
#' probabilities. Deviations are penalised as normal random effects; the
#' variance is optionally estimated via \code{Movement_cont_pe_pars_spec}.
#' If \code{do_recruits_move = 0}, age-1 deviations are fixed at zero.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param move_type Integer. Movement model formulation: \code{0} =
#'   unstructured Markov; \code{1} = CTMC. Default \code{0}.
#' @param do_recruits_move Integer flag. \code{0} = age-1 fish do not move
#'   (default); movement deviations and CTMC rows for the minimum age are
#'   fixed at zero. \code{1} = recruits participate in movement.
#' @param use_fixed_movement Integer flag. \code{0} = estimate movement
#'   (default); \code{1} = fix movement rates to \code{Fixed_Movement} and
#'   map all movement parameters to \code{NA}.
#' @param Fixed_Movement Numeric array of externally supplied movement
#'   probability matrices, dimensioned
#'   \code{[n_pop × n_regions × n_regions × n_years × n_seas × n_ages × n_sexes]}.
#'   Each \code{[n_regions × n_regions]} slice must be row-stochastic (rows
#'   sum to 1). Required when \code{use_fixed_movement = 1}. If \code{NA}
#'   (default), an identity matrix (no movement) is constructed internally.
#' @param Use_Movement_Prior Integer flag. \code{1} = apply Dirichlet priors
#'   to movement row probabilities; \code{0} = no priors (default). Requires
#'   \code{Movement_prior}.
#' @param Movement_prior Data frame of Dirichlet prior concentration
#'   parameters. Required columns: \code{pop}, \code{region_from},
#'   \code{year}, \code{seas}, \code{age}, \code{sex}, and \code{alpha},
#'   where \code{alpha} is a list-column with each element a numeric vector
#'   of length \code{n_regions} giving the Dirichlet concentration for
#'   transitions out of \code{region_from}. Values near 1 are uninformative;
#'   larger values concentrate the prior toward equal movement. Only used
#'   when \code{Use_Movement_Prior = 1}.
#' @param Movement_popblk_spec \code{"constant"} (default, shared across all
#'   populations) or a list of integer vectors partitioning populations into
#'   blocks. Example: \code{list(c(1, 2), 3)} shares parameters for
#'   populations 1 and 2 and estimates a separate parameter for population 3.
#'   Ignored when \code{move_type = 1}.
#' @param Movement_ageblk_spec \code{"constant"} (default) or a list of
#'   integer vectors defining age blocks. Example: \code{list(1:4, 5:10)}
#'   creates a juvenile block (ages 1-4) and an adult block (ages 5-10).
#'   Ignored when \code{move_type = 1}.
#' @param Movement_yearblk_spec \code{"constant"} (default) or a list of
#'   integer vectors defining year blocks for discrete structural breaks in
#'   movement. For residual annual variation, use \code{cont_vary_movement}
#'   instead. Ignored when \code{move_type = 1}.
#' @param Movement_seasblk_spec \code{"constant"} (default) or a list of
#'   integer vectors defining season blocks. Example: \code{list(c(1, 2),
#'   c(3, 4))} groups winter/spring and summer/fall. Ignored when
#'   \code{move_type = 1}.
#' @param Movement_sexblk_spec \code{"constant"} (default, sex-invariant) or
#'   a list of integer vectors defining sex blocks. Example: \code{list(1, 2)}
#'   estimates sex-specific movement independently. Ignored when
#'   \code{move_type = 1}.
#' @param cont_vary_movement Character string specifying the structure of
#'   continuous iid movement deviations added on top of the fixed-effect
#'   movement surface. Default \code{"none"}. Options:
#'   \describe{
#'     \item{\code{"none"}}{No deviations.}
#'     \item{\code{"iid_y"}}{Year-varying; shared across pop, age, sex, season.}
#'     \item{\code{"iid_a"}}{Age-varying; shared across pop, year, sex, season.}
#'     \item{\code{"iid_y_a"}}{Year \eqn{\times} age.}
#'     \item{\code{"iid_y_a_s"}}{Year \eqn{\times} age \eqn{\times} sex.}
#'     \item{\code{"iid_y_seas_a_s"}}{Year \eqn{\times} season \eqn{\times} age \eqn{\times} sex.}
#'     \item{\code{"iid_p_y"}, \code{"iid_p_a"}, \code{"iid_p_y_a"}, \code{"iid_p_y_a_s"}, \code{"iid_p_y_seas_a_s"}}{Population-specific analogues of the above.}
#'   }
#' @param Movement_cont_pe_pars_spec Character string specifying estimation of
#'   process-error variance for \code{cont_vary_movement} deviations. One of:
#'   \describe{
#'     \item{\code{"none"}}{No process-error parameters; use with \code{cont_vary_movement = "none"}.}
#'     \item{\code{"fix"}}{Parameters initialised but not estimated; fixes
#'       deviation variance at its starting value.}
#'     \item{\code{"est_shared"}}{Single variance estimated, shared across all
#'       dimensions.}
#'     \item{\code{"est_all"}}{All variance parameters estimated independently,
#'       dimensioned \code{[n_pop × n_regions × n_seas × n_ages × n_sexes]}.}
#'   }
#' @param ctmc_move_dat Data frame required when \code{move_type = 1}. Each
#'   row corresponds to a unique pop-region-year-season-age-sex combination.
#'   Required columns: \code{pop}, \code{regions}, \code{years}, \code{seas},
#'   \code{ages}, \code{sexes}, plus any covariate columns referenced in
#'   \code{diffusion_formula} or \code{preference_formula}. Projection years
#'   exceeding \code{n_years} are automatically capped to the final estimation
#'   year to prevent spline extrapolation.
#' @param adjacency_mat Square numeric matrix \code{[n_regions × n_regions]}
#'   with 1 indicating an allowed transition and 0 indicating no direct
#'   connection; diagonal should be 0. Required for \code{move_type = 1}.
#'   For \code{move_type = 0} a fully connected matrix is constructed
#'   automatically.
#' @param area_r Numeric vector of length \code{n_regions} giving the area of
#'   each region, used to scale CTMC diffusion rates. Required for
#'   \code{move_type = 1}. Default: \code{rep(1, n_regions)}.
#' @param diffusion_formula R \code{formula} defining the linear predictor for
#'   the CTMC diffusion (\eqn{\theta}) component (e.g.,
#'   \code{~ bs(depth, df = 4)}). All right-hand-side variables must be
#'   present in \code{ctmc_move_dat}. Required for \code{move_type = 1}.
#' @param preference_formula R \code{formula} defining the linear predictor for
#'   the CTMC habitat-preference (taxis, \eqn{\gamma}) component. All
#'   variables must be present in \code{ctmc_move_dat}. Required for
#'   \code{move_type = 1}.
#' @param ctmc_diffusion_bounds Integer flag. \code{1} = apply bounds to
#'   diffusion parameters to ensure the CTMC generator matrix is a valid
#'   Metzler matrix (non-negative off-diagonal entries). \code{0} = no bounds
#'   (default).
#' @param move_timing Integer flag setting how movement and mortality are
#'   sequenced within a season. \code{0} = movement then mortality (default,
#'   historical SPoRC behaviour); \code{1} = mortality then movement;
#'   \code{2} = continuous, with movement and mortality acting simultaneously via
#'   the matrix exponential of \eqn{Q\Delta - \mathrm{diag}(Z)}. \code{move_timing = 2}
#'   requires an estimated CTMC generator, i.e. \code{move_type = 1} and
#'   \code{use_fixed_movement = 0}.
#' @param ctmc_scale_by_seasdur Integer flag controlling the time units of the CTMC
#'   generator. \code{1} (default) treats \eqn{Q} as an annual rate, exponentiating
#'   \eqn{Q \cdot \mathrm{seasdur}[s]} in each season so that movement and mortality
#'   share time units. \code{0} exponentiates
#'   \eqn{Q} once per season regardless of duration. Only has an effect when
#'   \code{move_type = 1} and \code{n_seas > 1}; forced to \code{1} when
#'   \code{move_timing = 2}.
#' @param move_expm_nsub Integer controlling how matrix exponentials of the CTMC
#'   generator are evaluated, both when converting \eqn{Q} to movement fractions and
#'   inside the \code{move_timing = 2} seasonal operators. \code{0} (default)
#'   uses \code{Matrix::expm}. A power of two \eqn{n \ge 1} instead uses \eqn{n} implicit
#'   (backward Euler) substeps, \eqn{(I - A/n)^{-n}}, evaluated as one linear solve plus
#'   \eqn{\log_2 n} squarings, which is why \eqn{n} must be a power of two. The implicit form has a much cheaper reverse-mode adjoint
#'   than a matrix exponential, so the gradient is several times faster, but it is a
#'   first-order approximation: \eqn{n = 1} is plain \code{solve(I - A)} and is an approximation.
#' @param ... Optional starting value overrides, passed by name. Recognised
#'   arguments:
#'   \describe{
#'     \item{\code{move_pars}}{Array \code{[n_pop × n_regions × (n_regions-1) × n_years × n_seas × n_ages × n_sexes]}. Default: \code{0} (equal movement on logit scale).}
#'     \item{\code{log_move_diffusion_pars}}{Vector of length \code{n_theta}. Default: \code{log(0.1)}.}
#'     \item{\code{move_preference_pars}}{Vector of length \code{n_gamma}. Default: \code{0}.}
#'     \item{\code{move_devs}}{Array \code{[n_pop × n_regions × (n_regions-1) × (n_years + n_proj_yrs_devs) × n_seas × n_ages × n_sexes]}. Default: \code{0}.}
#'     \item{\code{move_pe_pars}}{Array \code{[n_pop × n_regions × n_seas × n_ages × n_sexes]}. Default: \code{0}.}
#'   }
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated. Key additions to \code{$data} include
#'   \code{move_type}, \code{use_fixed_movement}, \code{Fixed_Movement},
#'   \code{adjacency_mat}, \code{adjacency_collapsed}, \code{area_r},
#'   \code{ctmc_move_dat}, \code{diffusion_formula}, \code{preference_formula},
#'   and \code{cont_vary_movement} (stored as an integer code). Parameter
#'   arrays \code{move_pars}, \code{log_move_diffusion_pars},
#'   \code{move_preference_pars}, \code{move_devs}, and \code{move_pe_pars}
#'   are added to \code{$par}, with corresponding factor maps in \code{$map}.
#'
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
                               move_timing = 0,
                               ctmc_scale_by_seasdur = 1,
                               move_expm_nsub = 0,
                               ...
) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...) # get starting values if there are any
  if(input_list$store_config) input_list$config$Setup_Mod_Movement <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # If no movement matrix is provided
  if(is.na(sum(Fixed_Movement))) {
    Fixed_Movement <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_regions,
                                       length(input_list$data$years), input_list$data$n_seas,
                                       length(input_list$data$ages), input_list$data$n_sexes))
    for(p in 1:input_list$data$n_pop) Fixed_Movement[p,,,,,,] <- diag(1, input_list$data$n_regions)
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

  # Check movement / mortality sequencing
  if(!move_timing %in% c(0, 1, 2))
    stop('move_timing is not correctly specified. The options are move_timing == 0 (movement then mortality),
         == 1 (mortality then movement), or == 2 (continuous, simultaneous movement and mortality)')

  collect_message("Movement timing is: ", c("Movement then mortality",
                                            "Mortality then movement",
                                            "Continuous (simultaneous)")[move_timing + 1])

  # Continuous movement needs an instantaneous rate matrix, which only exists for an
  # estimated CTMC. A discrete multinomial-logit matrix has no guaranteed real generator
  # (the Markov embedding problem), so we refuse rather than attempt a matrix logarithm.
  if(move_timing == 2) {
    if(move_type != 1)
      stop("move_timing == 2 (continuous movement) requires move_type == 1 (CTMC). ",
           "Unstructured multinomial-logit movement has no instantaneous rate matrix, and one ",
           "cannot in general be recovered from the movement fractions.")
    if(use_fixed_movement == 1)
      stop("move_timing == 2 (continuous movement) requires use_fixed_movement == 0. ",
           "A fixed movement matrix supplies transition fractions, not the instantaneous rates ",
           "that continuous movement needs.")
  }

  if(!ctmc_scale_by_seasdur %in% c(0, 1))
    stop('ctmc_scale_by_seasdur is not correctly specified. The options are 0 (unscaled, one exponentiation per season) or 1 (scale by season duration)')

  # Matrix exponential evaluation. Backward Euler trades accuracy for a quicker eval
  if(!is.numeric(move_expm_nsub) || length(move_expm_nsub) != 1 || is.na(move_expm_nsub) ||
     move_expm_nsub != as.integer(move_expm_nsub) || move_expm_nsub < 0)
    stop('move_expm_nsub is not correctly specified. It must be a single non-negative integer: 0 (exact, Matrix::expm) or the number of implicit backward Euler substeps.')

  move_expm_nsub <- as.integer(move_expm_nsub)

  # Substeps are applied by repeated squaring, which reaches powers of two exactly and
  # nothing else. The scheme is first order regardless, so a finer ladder would buy nothing.
  if(move_expm_nsub > 0 && bitwAnd(move_expm_nsub, move_expm_nsub - 1L) != 0)
    stop('move_expm_nsub must be a power of two (1, 2, 4, 8, ... ), since substeps are applied by repeated squaring. Got ', move_expm_nsub, '.')

  if(move_type == 1 && use_fixed_movement == 0) {
    if(move_expm_nsub == 0) collect_message("Matrix exponential is: exact (Matrix::expm)")
    else collect_message("Matrix exponential is: implicit backward Euler with ", move_expm_nsub,
                         " substep(s), (I - A/n)^-n. First order in 1/n; the exponential is approximated, not reproduced.")
  } else if(move_expm_nsub != 0) {
    collect_message("move_expm_nsub ignored: matrix exponentials are only taken for an estimated CTMC generator (move_type = 1 with use_fixed_movement = 0).")
    move_expm_nsub <- 0L
  }

  # Mixing an unscaled generator with seasdur-scaled mortality is dimensionally inconsistent,
  # so continuous movement forces the scaling on.
  if(move_timing == 2 && ctmc_scale_by_seasdur == 0) {
    collect_message("ctmc_scale_by_seasdur forced to 1: continuous movement requires the generator and mortality to share time units.")
    ctmc_scale_by_seasdur <- 1
  }

  if(move_type == 1 && input_list$data$n_seas > 1) {
    collect_message("CTMC generator time units: ", ifelse(ctmc_scale_by_seasdur == 1,
                                                          "annual rate, scaled by seasdur each season",
                                                          "per-season rate, unscaled"))
  }

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
  input_list$data$move_timing <- move_timing
  input_list$data$ctmc_scale_by_seasdur <- ctmc_scale_by_seasdur
  input_list$data$move_expm_nsub <- move_expm_nsub

  # define for continuous varying movement
  cont_move_map <- data.frame(
    type = c("none", "iid_y", "iid_a", "iid_y_a", "iid_y_a_s", "iid_y_seas_a_s",
             "iid_p_y", "iid_p_a", "iid_p_y_a", "iid_p_y_a_s", "iid_p_y_seas_a_s"),
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
  # A preference formula with no terms (e.g. ~ 0) means pure diffusion with no taxis.
  # RTMB does not accept a zero-length parameter, so keep a length-1 placeholder; it is
  # mapped off below and never enters the generator.
  else input_list$par$move_preference_pars <- rep(0, max(n_gamma, 1))

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

  # Pure diffusion (preference formula with no terms): fix the placeholder preference
  # parameter so it is never estimated
  if(move_type == 1 && n_gamma == 0) {
    input_list$map$move_preference_pars <- factor(rep(NA, length(input_list$par$move_preference_pars)))
    collect_message("Preference formula has no terms: movement is pure diffusion (no taxis).")
  }

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
