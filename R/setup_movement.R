# Stage 1 of 3: model setup
#
# Movement inputs: how fish redistribute among regions and when in the season. Chooses unstructured matrices
# (move_type 0) or a continuous time generator (move_type 1), sets move_timing, and builds the map.

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
    if(is.character(Movement_popblk_spec)) {
      if(!identical(Movement_popblk_spec, "constant")) stop("Movement_popblk_spec must be \"constant\" or a list of population blocks, but was: ", Movement_popblk_spec)
      Movement_popblk_spec_vals <- list(1:input_list$data$n_pop)
    } else Movement_popblk_spec_vals <- Movement_popblk_spec


    # If movement is constant for ages
    if(is.character(Movement_ageblk_spec)) {
      if(!identical(Movement_ageblk_spec, "constant")) stop("Movement_ageblk_spec must be \"constant\" or a list of age blocks, but was: ", Movement_ageblk_spec)
      Movement_ageblk_spec_vals <- list(seq_along(input_list$data$ages))
    } else Movement_ageblk_spec_vals <- Movement_ageblk_spec

    # If movement is constant across years
    if(is.character(Movement_yearblk_spec)) {
      if(!identical(Movement_yearblk_spec, "constant")) stop("Movement_yearblk_spec must be \"constant\" or a list of year blocks, but was: ", Movement_yearblk_spec)
      Movement_yearblk_spec_vals <- list(seq_along(input_list$data$years))
    } else Movement_yearblk_spec_vals <- Movement_yearblk_spec

    # If movement is constant across sexes
    if(is.character(Movement_sexblk_spec)) {
      if(!identical(Movement_sexblk_spec, "constant")) stop("Movement_sexblk_spec must be \"constant\" or a list of sex blocks, but was: ", Movement_sexblk_spec)
      Movement_sexblk_spec_vals <- list(1:input_list$data$n_sexes)
    } else Movement_sexblk_spec_vals <- Movement_sexblk_spec

    # If movement is constant across seasons
    if(is.character(Movement_seasblk_spec)) {
      if(!identical(Movement_seasblk_spec, "constant")) stop("Movement_seasblk_spec must be \"constant\" or a list of season blocks, but was: ", Movement_seasblk_spec)
      Movement_seasblk_spec_vals <- list(1:input_list$data$n_seas)
    } else Movement_seasblk_spec_vals <- Movement_seasblk_spec

    # If spatial model
    if(input_list$data$n_regions > 1 &&
       input_list$data$use_fixed_movement == 0 # if not using fixed movement matrix
       ) {

      # Initialize counter
      counter <- 1

      for(popblk in seq_along(Movement_popblk_spec_vals)) {
        # get populations to block and map off
        map_p <- Movement_popblk_spec_vals[[popblk]]

        for(ageblk in seq_along(Movement_ageblk_spec_vals)) {
          # get ages to block and map off
          map_a <- Movement_ageblk_spec_vals[[ageblk]]

          for(yearblk in seq_along(Movement_yearblk_spec_vals)) {
            # get years to block and map off
            map_y <- Movement_yearblk_spec_vals[[yearblk]]

            for(seasblk in seq_along(Movement_seasblk_spec_vals)) {
              # get seasons to block and map off
              map_seas <- Movement_seasblk_spec_vals[[seasblk]]

              for(sexblk in seq_along(Movement_sexblk_spec_vals)) {
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
    if(length(map_log_move_diffusion_pars) != 0) map_log_move_diffusion_pars <- factor(seq_along(map_log_move_diffusion_pars))
    if(length(map_move_preference_pars) != 0) map_move_preference_pars <- factor(seq_along(map_move_preference_pars))
  }

  # Input into mapping list
  input_list$map$move_pars <- factor(map_Movement_Pars)
  input_list$map$log_move_diffusion_pars <- factor(map_log_move_diffusion_pars)
  input_list$map$move_preference_pars <- factor(map_move_preference_pars)

  return(input_list)
}

#' Preference terms that vary over years
#'
#' Internal helper called by \code{\link{Setup_Mod_DSEM}}. A dsem on the movement
#' deviations describes the year to year part of habitat preference, so a
#' preference covariate that also varies over years writes that part a second
#' time, and the formula's version of it is not penalized.
#'
#' @param input_list Named list with \code{$data}, built through
#'   \code{\link{Setup_Mod_Movement}}.
#'
#' @return Character vector naming the preference design columns whose values
#'   change across years within a population, region, season, age and sex. Empty
#'   when movement is not CTMC, when the preference formula has no terms, or when
#'   none of its terms vary over years.
#'
#' @keywords internal
get_yr_varying_pref_terms <- function(input_list) {

  dat <- input_list$data$ctmc_move_dat
  if(!isTRUE(input_list$data$move_type == 1) || is.null(dat)) return(character(0))

  design <- get_movement_dp_design_matrix(dat, input_list$data$preference_formula, input_list$data$diffusion_formula)
  if(design$n_gamma == 0) return(character(0)) # pure diffusion, so there is no preference term to write twice
  W_zk <- design$W_zk

  # rows of one stratum differ only in their year, so comparing each row with the first row of its
  # stratum finds every column whose value the formula lets change from year to year
  key <- do.call(paste, c(lapply(c("pop", "regions", "seas", "ages", "sexes"), function(v) dat[,v]), sep = "_"))
  first <- match(key, key)
  varies <- apply(W_zk != W_zk[first,,drop = FALSE], 2, any)

  return(colnames(W_zk)[varies])

}

#' Map continuous movement deviation and process-error parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Movement}} to construct the
#' TMB/RTMB factor maps for \code{move_devs} (iid deviations on the movement
#' logit or log-rate surface) and \code{move_pe_pars} (process-error variance
#' parameters). Deviations are only activated when the model is spatial
#' (\code{n_regions > 1}), continuous variation is requested
#' (\code{cont_vary_movement} is not \code{"none"}), and movement is estimated
#' (\code{use_fixed_movement == 0}). For CTMC movement, deviations sit on each
#' region's preference, so the only ones left unestimated belong to a region no
#' edge of the adjacency matrix touches. The resulting integer map is also
#' stored as \code{$data$map_move_devs} for use in the C++ template.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param cont_vary_movement Character string specifying the deviation structure:
#'   \code{"none"}, or \code{"iid_"} followed by the dims the deviations vary
#'   over, any of p, y, seas, a, s in any order (\code{"iid_y"}, \code{"iid_y_a_s"},
#'   \code{"iid_p_y_seas_a_s"}, any combination). Dimensions present in the
#'   string receive unique estimation indices; absent dimensions share a single
#'   index. \code{"none"} maps all deviations to \code{NA}.
#' @param Movement_cont_pe_pars_spec Character string controlling estimation of
#'   the process-error variance for movement deviations. One of:
#'   \describe{
#'     \item{\code{"none"} or \code{"fix"}}{All \code{move_pe_pars} kept fixed
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
#'       Active cells receive sequential integer indices; CTMC regions no edge
#'       touches and inactive configurations are \code{NA}.}
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
     input_list$data$cont_vary_movement != "none" && # if continuous varying movement
     input_list$data$use_fixed_movement == 0 # if not using fixed movement matrix
  ) {

    # Dimensions (every region pair, or every region under the CTMC, gets its own set of estimated deviations, never shared with another)
    dims <- c(pop = input_list$data$n_pop,
              region_from = input_list$data$n_regions,
              region_to = n_regions_to,
              year = length(input_list$data$years) + input_list$data$n_proj_yrs_devs,
              season = input_list$data$n_seas,
              age = length(input_list$data$ages),
              sex = input_list$data$n_sexes)

    # dims named in the spec (e.g. "iid_y_a_s" -> year, age, sex) get unique value per combination and those dims not named are shared/broadcast
    dim_abbrev <- c(p = "pop", y = "year", seas = "season", a = "age", s = "sex")
    key_extra <- unname(dim_abbrev[strsplit(sub("^iid_", "", cont_vary_movement), "_")[[1]]])
    share_over <- setdiff(dim_abbrev, key_extra)

    map_move_devs <- build_pe_map(dims, share_over = share_over)

    # whether recruits (age 1) move
    if("age" %in% key_extra && input_list$data$do_recruits_move == 0 && dims["age"] >= 2) {
      map_move_devs[,,,,,1,] <- NA
    }

    # ctmc deviations
    if(input_list$data$move_type == 1) {
      for(r in 1:input_list$data$n_regions) {
        isolated <- all(input_list$data$adjacency_mat[r,] == 0) && all(input_list$data$adjacency_mat[,r] == 0)
        if(isolated) map_move_devs[,r,,,,,] <- NA
      } # end r
    }
  }

    # Movement Process Error Parameters ---------------------------------------

    # Mapping for movement process error deviations
    if(Movement_cont_pe_pars_spec %in% c("fix", "none")) map_move_pe_pars <- map_move_pe_pars
    if(Movement_cont_pe_pars_spec == 'est_all') map_move_pe_pars[] <- seq_along(map_move_pe_pars)
    if(Movement_cont_pe_pars_spec == 'est_shared') map_move_pe_pars[] <- 1

    # return to input list
    input_list$map$move_devs <- factor(as.vector(map_move_devs))
    input_list$data$map_move_devs <- array(as.numeric(input_list$map$move_devs), dim = dim(input_list$par$move_devs))
    input_list$map$move_pe_pars <- factor(map_move_pe_pars)
    return(input_list)

}

#' Set up movement model inputs and parameter structures
#'
#' Sets up unstructured Markov transition movement (\code{move_type = 0}) or a
#' continuous time Markov chain (\code{move_type = 1}), with optional iid
#' deviations on the movement surface, and builds the parameter arrays and factor
#' maps. Call after \code{\link{Setup_Mod_Biologicals}}.
#'
#' @section Unstructured Markov movement (\code{move_type = 0}):
#' Transitions out of region \eqn{r} are a multinomial logit with a reference cell,
#' so \code{move_pars} is \code{[n_pop × n_regions × (n_regions - 1) × n_years ×
#' n_seas × n_ages × n_sexes]}. The \code{Movement_*blk_spec} arguments share
#' parameters: indices in one block take the same factor level. A fully connected
#' adjacency matrix is built automatically. Blocks and continuous time variation
#' combine: use \code{Movement_yearblk_spec} for structural breaks and
#' \code{cont_vary_movement} for residual year-to-year variation.
#'
#' @section CTMC movement (\code{move_type = 1}):
#' The rate matrix \eqn{Q} is decomposed into diffusion (\eqn{\theta}) and
#' preference (\eqn{\gamma}), with design matrices from \code{diffusion_formula}
#' and \code{preference_formula} evaluated on \code{ctmc_move_dat}, and each time
#' step's movement matrix is \eqn{\exp(Q \Delta t)}. Blocking is not supported, so
#' every \code{Movement_*blk_spec} must stay \code{"constant"}; put structure
#' across populations, ages, sexes or seasons into formula covariates instead.
#'
#' @section Continuous movement deviations:
#' Deviations are added to the movement logit surface under unstructured movement,
#' or to each region's preference under CTMC, before probabilities are computed,
#' and are penalized as normal random effects whose variance
#' \code{Movement_cont_pe_pars_spec} can estimate. Age-1 deviations are fixed at
#' zero when \code{do_recruits_move = 0}.
#'
#' The two types size the deviations differently. Unstructured movement holds one
#' per origin and destination pair, \code{[n_regions x (n_regions - 1)]}, while the
#' CTMC holds one per region, \code{[n_regions x 1]}, since preference is a surface
#' over regions rather than a rate along an edge. A CTMC deviation raises the rates
#' into its region and lowers those out of it on every edge at once. Only
#' differences in preference reach the generator, so a constant added to every
#' region's deviation leaves movement unchanged and nothing but the process error
#' penalty holds the level of the field down.
#'
#' Because the deviations are preference they enter the generator additively, and
#' one larger than an edge's diffusion rate drives that rate negative. Under
#' \code{ctmc_diffusion_bounds = "none"} the movement fractions then leave
#' \code{[0, 1]} while still summing to one, so use \code{"upwind"} or
#' \code{"softplus"} whenever the deviations are estimated.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}.
#' @param move_type Integer. \code{0} (default) unstructured Markov, \code{1} CTMC.
#' @param do_recruits_move Integer flag. \code{0} (default) fixes the movement
#'   deviations and CTMC rows at the minimum age to zero, \code{1} moves recruits.
#' @param use_fixed_movement Integer flag. \code{0} (default) estimates movement,
#'   \code{1} fixes it at \code{Fixed_Movement} and maps every movement parameter
#'   to \code{NA}.
#' @param Fixed_Movement Movement probability array \code{[n_pop × n_regions ×
#'   n_regions × n_years × n_seas × n_ages × n_sexes]}, each \code{[n_regions ×
#'   n_regions]} slice row-stochastic. Required when
#'   \code{use_fixed_movement = 1}. \code{NA} (default) builds an identity matrix.
#' @param Use_Movement_Prior Integer flag, \code{1} for Dirichlet priors on the
#'   movement rows. Default \code{0}.
#' @param Movement_prior Data frame with columns \code{pop}, \code{region_from},
#'   \code{year}, \code{seas}, \code{age}, \code{sex} and \code{alpha}, the last a
#'   list-column of length-\code{n_regions} concentrations for transitions out of
#'   \code{region_from}. Values near 1 are uninformative, larger ones concentrate
#'   toward equal movement. Read when \code{Use_Movement_Prior = 1}.
#' @param Movement_popblk_spec,Movement_ageblk_spec,Movement_yearblk_spec,Movement_seasblk_spec,Movement_sexblk_spec
#'   Blocking across populations, ages, years, seasons and sexes: \code{"constant"}
#'   (default) or a list of integer vectors, e.g. \code{list(c(1, 2), 3)} for
#'   populations, \code{list(1:4, 5:10)} for a juvenile and an adult block, or
#'   \code{list(1, 2)} for sex-specific movement. Use
#'   \code{Movement_yearblk_spec} for structural breaks and
#'   \code{cont_vary_movement} for residual annual variation. All are ignored when
#'   \code{move_type = 1}.
#' @param cont_vary_movement Structure of the continuous deviations on the
#'   fixed-effect movement surface. \code{"none"} (default), or \code{"iid_"}
#'   followed by the dims they vary over, any of p (population), y (year), seas
#'   (season), a (age) and s (sex) in any order: \code{"iid_y"} is one deviation per
#'   year and region pair, or per year and region under CTMC movement, shared
#'   across everything else, and \code{"iid_p_y_seas_a_s"} varies by every dim. A
#'   dim left out shares one deviation across it. They are random effects with
#'   \code{Movement_cont_pe_pars_spec} estimating the sd and
#'   \code{random = "move_devs"} in \code{\link{fit_model}}. \code{"dsem"} instead
#'   hands their density to the arrows given to \code{\link{Setup_Mod_DSEM}}, one
#'   series per origin and destination (per region under CTMC, whose deviations
#'   hold no destination) and per level of every other dim with more than one,
#'   which it names itself since a deviation shared across a dim cannot be linked;
#'   \code{move_pe_pars} are then read by nothing.
#' @param Movement_cont_pe_pars_spec Estimation of the process error variance for
#'   the \code{cont_vary_movement} deviations. \code{"none"} creates no parameters
#'   and pairs with \code{cont_vary_movement = "none"}, \code{"fix"} holds the
#'   variance at its starting value, \code{"est_shared"} estimates one shared
#'   value, and \code{"est_all"} estimates \code{[n_pop × n_regions × n_seas ×
#'   n_ages × n_sexes]} independently.
#' @param ctmc_move_dat Data frame required when \code{move_type = 1}, one row per
#'   population, region, year, season, age and sex, with columns \code{pop},
#'   \code{regions}, \code{years}, \code{seas}, \code{ages}, \code{sexes} and any
#'   covariates the formulas name. Projection years beyond \code{n_years} are
#'   capped at the final estimation year to prevent spline extrapolation.
#' @param adjacency_mat Square \code{[n_regions × n_regions]} matrix, 1 for an
#'   allowed transition and 0 for none. The diagonal must be 0: residency falls out
#'   of the generator, and a non-zero diagonal leaves the generator columns summing
#'   to something other than zero, so the movement matrix loses abundance rather
#'   than redistributing it. A fully connected matrix is \code{1 - diag(n_regions)}
#'   (\code{diag(1, n_regions)} is the identity, not an adjacency matrix). Required
#'   under \code{move_type = 1}, where it is validated for dimension, 0/1 entries,
#'   a zero diagonal and at least one connection; built automatically under
#'   \code{move_type = 0}.
#' @param area_r Numeric vector \code{[n_regions]} of region areas, used to scale
#'   the CTMC diffusion rates. Required under \code{move_type = 1}. Default
#'   \code{rep(1, n_regions)}.
#' @param diffusion_formula Formula for the CTMC diffusion (\eqn{\theta}) linear
#'   predictor, e.g. \code{~ bs(depth, df = 4)}. Every right-hand-side variable
#'   must be in \code{ctmc_move_dat}. Required under \code{move_type = 1}.
#' @param preference_formula Formula for the CTMC preference (taxis, \eqn{\gamma})
#'   linear predictor, on the same terms. Required under \code{move_type = 1}.
#' @param ctmc_diffusion_bounds How the CTMC generator is kept a valid Metzler
#'   matrix when taxis outweighs diffusion. \code{"softplus"} takes a softplus of
#'   \eqn{\theta_j + d} of width \code{ctmc_diffusion_eps}; \code{"upwind"} takes
#'   the finite volume flux \eqn{\theta_j + \max(d, 0)}, which keeps diffusion
#'   whole and adds only the down-gradient taxis, so positivity never depends on
#'   the two cancelling.
#' @param ctmc_diffusion_eps Positive width of the softplus under
#'   \code{ctmc_diffusion_bounds = "softplus"}. Default \code{0.1}. An edge where
#'   taxis exactly cancels diffusion has \code{eps * log(2)}, so this is a floor on
#'   exchange as well as a smoothing constant.
#' @param move_timing How movement and mortality are sequenced within a season.
#'   \code{0} (default) moves then kills, \code{1} kills then moves, and \code{2}
#'   runs the two together through the matrix exponential of
#'   \eqn{Q\Delta - \mathrm{diag}(Z)}. \code{2} needs an estimated CTMC generator,
#'   so \code{move_type = 1} and \code{use_fixed_movement = 0}.
#' @param ctmc_scale_by_seasdur Integer flag for the time units of the CTMC
#'   generator. \code{1} (default) treats \eqn{Q} as an annual rate and
#'   exponentiates \eqn{Q \cdot \mathrm{seasdur}[s]} each season, so movement and
#'   mortality share time units; \code{0} exponentiates \eqn{Q} once per season
#'   whatever its duration. Only matters under \code{move_type = 1} with
#'   \code{n_seas > 1}, and is forced to \code{1} under \code{move_timing = 2}.
#' @param move_expm_nsub How matrix exponentials of the generator are evaluated,
#'   both converting \eqn{Q} to movement fractions and inside the
#'   \code{move_timing = 2} seasonal operators. \code{0} (default) uses
#'   \code{Matrix::expm}. A power of two \eqn{n \ge 1} uses \eqn{n} implicit
#'   backward Euler substeps, \eqn{(I - A/n)^{-n}}, as one linear solve plus
#'   \eqn{\log_2 n} squarings, which is why \eqn{n} must be a power of two. Its
#'   reverse-mode derivative is much cheaper, so the gradient is several times
#'   faster, but it is a first-order approximation and \eqn{n = 1} is plain
#'   \code{solve(I - A)}.
#' @param ... Optional starting values by name: \code{move_pars} \code{[n_pop ×
#'   n_regions × (n_regions-1) × n_years × n_seas × n_ages × n_sexes]}, default
#'   \code{0}; \code{log_move_diffusion_pars} of length \code{n_theta}, default
#'   \code{log(0.1)}; \code{move_preference_pars} of length \code{n_gamma}, default
#'   \code{0}; \code{move_devs}, shaped as \code{move_pars} with
#'   \code{n_years + n_proj_yrs_devs} years and the third dim \code{1} under
#'   \code{move_type = 1}, default \code{0}; and \code{move_pe_pars} \code{[n_pop ×
#'   n_regions × n_seas × n_ages × n_sexes]}, default \code{0}.
#'
#' @return \code{input_list} with \code{$data}, \code{$par} and \code{$map}
#'   updated. \code{$data} gains \code{move_type}, \code{use_fixed_movement},
#'   \code{Fixed_Movement}, \code{adjacency_mat}, \code{adjacency_collapsed},
#'   \code{area_r}, \code{ctmc_move_dat}, \code{diffusion_formula},
#'   \code{preference_formula} and \code{cont_vary_movement} as its form string.
#'   \code{move_pars}, \code{log_move_diffusion_pars}, \code{move_preference_pars},
#'   \code{move_devs} and \code{move_pe_pars} go into \code{$par} with their factor
#'   maps in \code{$map}.
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
                               ctmc_diffusion_eps = 0.1,
                               move_timing = 0,
                               ctmc_scale_by_seasdur = 1,
                               move_expm_nsub = 0,
                               ...
) {

  move_pe_spec_given <- !missing(Movement_cont_pe_pars_spec) # read before anything assigns it
  messages_list <<- character(0) # string to attach to for printing messages # nolint: object_usage_linter.
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
  if(use_fixed_movement == 1) check_data_dimensions(
    Fixed_Movement,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_ages = length(input_list$data$ages),
    n_sexes = input_list$data$n_sexes,
    n_seas = input_list$data$n_seas,
    what = 'Fixed_Movement'
  )

  # Check for movement priors
  if(!Use_Movement_Prior %in% c(0,1)) stop('Options for movement priors not correctly specified. The options are Use_Movement_Prior == 0 (dont use), or == 1 (use)')
  else collect_message("Movement priors are: ", ifelse(Use_Movement_Prior == 0, "Not Used", "Used"))

  # Check for recruits moving
  if(!do_recruits_move %in% c(0,1)) stop('Movement for recruits is not correctly specified. The options are do_recruits_move == 0 (they dont move), or == 1 (they move)')
  else collect_message("Recruits are: ", ifelse(do_recruits_move == 0, "Not Moving", "Moving"))

  # Check movement continuous varying parameterization. "dsem" reuses the iid form
  # the dsem sets the linked cells to NA in map_move_devs, so the penalty doesn't use them and the dsem supplies their density
  if(identical(cont_vary_movement, "dsem")) cont_vary_movement <- paste(c("dsem", if(input_list$data$n_pop > 1) "p", "y", if(input_list$data$n_seas > 1) "seas", if(length(input_list$data$ages) > 1) "a", if(input_list$data$n_sexes > 1) "s"), collapse = "_") # names every dim with more than one level, since a deviation shared across a dim cannot be linked
  move_dsem <- grepl("^dsem_", cont_vary_movement)
  cont_vary_movement <- sub("^dsem_", "iid_", cont_vary_movement)
  dim_order <- c("p", "y", "seas", "a", "s") # the dims may be written in any order and are read in this one
  named <- strsplit(sub("^iid_", "", cont_vary_movement), "_")[[1]]
  form_ok <- identical(cont_vary_movement, "none") ||
    (grepl("^iid_", cont_vary_movement) && length(named) > 0 && !any(duplicated(named)) && all(named %in% dim_order))
  if(!form_ok)
    stop("cont_vary_movement should be 'none', 'iid_' followed by the dims the deviations vary over, any of p, y, seas, a, s in any order (iid_y, iid_y_a_s, iid_p_y_seas_a_s, ...), or 'dsem'.")
  if(grepl("^iid_", cont_vary_movement)) cont_vary_movement <- paste(c("iid", dim_order[dim_order %in% named]), collapse = "_")
  collect_message("Continuous movement specification is: ", if(move_dsem) sub("^iid_", "dsem_", cont_vary_movement) else cont_vary_movement)

  # Check movement process error estimation (no change needed here)
  if(!Movement_cont_pe_pars_spec %in% c('none', 'fix', 'est_all', 'est_shared'))
    stop('Options for continuous movement process error is not correctly specified.')
  else collect_message("Continuous movement process error specification is: ", Movement_cont_pe_pars_spec)

  # under the handover the process error sd is read by nothing, so it is fixed here
  if(move_dsem) {
    if(move_pe_spec_given && Movement_cont_pe_pars_spec %in% c("est_all", "est_shared")) stop("cont_vary_movement = 'dsem_...' takes the movement deviations' density from the dsem arrows, so move_pe_pars are read by nothing and cannot be estimated. Leave Movement_cont_pe_pars_spec out or set it to 'fix'.")
    # a dim the form leaves out shares one deviation across it, and a shared deviation cannot sit under the dsem
    named <- strsplit(sub("^iid_", "", cont_vary_movement), "_")[[1]]
    needed <- c(if(input_list$data$n_pop > 1) "p", "y", if(input_list$data$n_seas > 1) "seas", if(length(input_list$data$ages) > 1) "a", if(input_list$data$n_sexes > 1) "s")
    if(!all(needed %in% named)) stop(paste0("cont_vary_movement = 'dsem_", paste(named, collapse = "_"), "' shares a deviation across ", paste(setdiff(needed, named), collapse = ", "), ", and a shared deviation cannot be linked. Write cont_vary_movement = 'dsem', which names every dim itself, or name every dim the deviations vary over: dsem_", paste(needed, collapse = "_"), "."))
    Movement_cont_pe_pars_spec <- "fix"
    input_list$data$dsem_declared <- union(input_list$data$dsem_declared, "move")
    collect_message("cont_vary_movement = 'dsem_...': the movement deviations' density comes from Setup_Mod_DSEM, and move_pe_pars stay at their start.")
  }
  input_list$data$move_dsem <- as.numeric(move_dsem)

  if(!move_type %in% c(0, 1)) stop('move_type must be 0 (unstructured) or 1 (Continuous Time Markov Chain)')
  collect_message("Movement type is: ", ifelse(move_type == 0, "Unstructured Markov", "Continuous Time Markov Chain"))

  # Check movement / mortality sequencing
  if(!move_timing %in% c(0, 1, 2))
    stop('move_timing is not correctly specified. The options are move_timing == 0 (movement then mortality),
         == 1 (mortality then movement), or == 2 (continuous, simultaneous movement and mortality)')

  collect_message("Movement timing is: ", c("Movement then mortality",
                                            "Mortality then movement",
                                            "Continuous (simultaneous)")[move_timing + 1])

  # continuous movement needs an instantaneous rate matrix, which only exists for an estimated CTMC.
  # a discrete multinomial-logit matrix has no guaranteed real generator, so refuse rather than try
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
    if(anyNA(adjacency_mat)) stop("adjacency_mat cannot contain NA values")
    if(!all(adjacency_mat %in% c(0, 1))) stop("adjacency_mat entries must be 0 (not connected) or 1 (connected)")

    # Guard for adjacency matrix
    if(any(diag(adjacency_mat) != 0)) {
      stop("adjacency_mat must have a zero diagonal, since residency is implied by the generator rather than specified. ",
           "Non-zero diagonal entries in region(s): ", paste(which(diag(adjacency_mat) != 0), collapse = ", "), ". ",
           "Note that diag(1, n_regions) is the identity matrix, not an adjacency matrix; ",
           "a fully connected matrix is 1 - diag(n_regions).")
    }

    # with no off-diagonal connections nothing can move and diffusion is unidentifiable
    if(use_fixed_movement == 0 && all(adjacency_mat == 0)) {
      stop("adjacency_mat has no off-diagonal connections, so no transitions are possible and the CTMC diffusion ",
           "parameters are unidentifiable. Set use_fixed_movement = 1 to run without movement.")
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
    for(i in seq_len(nrow(Movement_prior))) {
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

  # taxis makes Q = D + Z, so without a bound an off diagonal of Q can go negative
  bound_form <- get_ctmc_bound_form(ctmc_diffusion_bounds)
  if(is.na(bound_form)) stop('ctmc_diffusion_bounds must be 0, 1, 2, or the matching name "none", "softplus", "upwind"')

  if(move_type == 1 && !is.null(preference_formula) && bound_form == "none") {
    pref_terms <- length(attr(stats::terms(preference_formula), "term.labels")) + attr(stats::terms(preference_formula), "intercept")
    if(pref_terms > 0) warning('preference_formula has terms but ctmc_diffusion_bounds is "none"; the generator can go invalid where taxis outweighs diffusion. A bounded form is recommended.')
    # the deviations are preference too, so they make taxis out of a model whose formula has none
    if(pref_terms == 0 && cont_vary_movement != "none") collect_message('cont_vary_movement puts deviations on preference, but ctmc_diffusion_bounds is "none"; a deviation larger than the diffusion rate makes the generator invalid. A bounded form is recommended.')
  }

  input_list$data$ctmc_diffusion_bounds <- ctmc_diffusion_bounds

  # softplus width for the bounds
  if(!is.numeric(ctmc_diffusion_eps) || length(ctmc_diffusion_eps) != 1 || !is.finite(ctmc_diffusion_eps) || ctmc_diffusion_eps <= 0) stop('ctmc_diffusion_eps must be a single positive number: the softplus width used when ctmc_diffusion_bounds = "softplus" (default 0.1)')
  if(move_type == 1) {
    msg <- switch(
      bound_form,
      none = "unbounded generator",
      softplus = paste0("softplus of width ", ctmc_diffusion_eps, " on the adjacency edges, so a cancelled edge has a floor of ", signif(ctmc_diffusion_eps * log(2), 3)),
      upwind = "discontinuous Galerkin upwind flux, diffusion kept whole"
    )
    collect_message("CTMC generator bounds: ", msg)
    if(bound_form == "upwind" && ctmc_diffusion_eps != 0.1) collect_message('ctmc_diffusion_eps is not read by the "upwind" form.')
  }
  input_list$data$ctmc_diffusion_eps <- ctmc_diffusion_eps
  input_list$data$move_timing <- move_timing
  input_list$data$ctmc_scale_by_seasdur <- ctmc_scale_by_seasdur
  input_list$data$move_expm_nsub <- move_expm_nsub

  input_list$data$cont_vary_movement <- cont_vary_movement

  # Populate Parameter List -------------------------------------------------

  # Movement Parameters (for unstructured markov; move_type == 0)
  input_list$par$move_pars <- array(0, dim = c(input_list$data$n_pop,
                                                    input_list$data$n_regions, input_list$data$n_regions - 1,
                                                    length(input_list$data$years), input_list$data$n_seas,
                                                    length(input_list$data$ages), input_list$data$n_sexes))
  input_list$par$move_pars <- use_starting_value(input_list$par$move_pars, starting_values, "move_pars")

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
  input_list$par$log_move_diffusion_pars <- rep(log(0.1), n_theta)
  input_list$par$log_move_diffusion_pars <- use_starting_value(input_list$par$log_move_diffusion_pars, starting_values, "log_move_diffusion_pars")

  # preference parameters
  input_list$par$move_preference_pars <- rep(0, max(n_gamma, 1))
  input_list$par$move_preference_pars <- use_starting_value(input_list$par$move_preference_pars, starting_values, "move_preference_pars")

  # Movement deviations. The unstructured model holds one per origin-destination pair, while the
  # CTMC holds one per region, on that region's preference, so its destination axis has length one
  {
    n_dev_to <- if(move_type == 1) 1 else input_list$data$n_regions - 1
    input_list$par$move_devs <- array(0, c(input_list$data$n_pop,
                                           input_list$data$n_regions, n_dev_to,
                                           length(input_list$data$years) + input_list$data$n_proj_yrs_devs,
                                           input_list$data$n_seas,
                                           length(input_list$data$ages),
                                           input_list$data$n_sexes))
  }
  input_list$par$move_devs <- use_starting_value(input_list$par$move_devs, starting_values, "move_devs")

  # Movement process error parameters
  input_list$par$move_pe_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                       input_list$data$n_seas, length(input_list$data$ages),
                                                       input_list$data$n_sexes)) # max 4 parameters or the ages
  input_list$par$move_pe_pars <- use_starting_value(input_list$par$move_pe_pars, starting_values, "move_pe_pars")


  # Mapping Options ---------------------------------------------------------
  input_list <- do_move_pars_mapping(input_list, Movement_popblk_spec, Movement_ageblk_spec, Movement_yearblk_spec, Movement_sexblk_spec, Movement_seasblk_spec, use_fixed_movement)
  input_list <- do_cont_vary_move_mapping(input_list, cont_vary_movement, Movement_cont_pe_pars_spec)

  # Pure diffusion (preference formula with no terms)
  if(move_type == 1 && n_gamma == 0) {
    input_list$map$move_preference_pars <- factor(rep(NA, length(input_list$par$move_preference_pars)))
    collect_message("Preference formula has no terms: movement is pure diffusion (no taxis).")
  }

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

