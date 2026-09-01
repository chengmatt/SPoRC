# Stage 1 of 3: model setup
#
# Growth, weight at age, maturity and natural mortality inputs, for both the
# estimation model (Setup_Mod_Biologicals) and the operating model
# (Setup_Sim_Biologicals).

#' Set up biological parameter inputs for closed-loop simulation
#'
#' Populates a simulation list (created by \code{\link{Setup_Sim_Dim}}) with
#' biological arrays needed to run the operating model: natural mortality,
#' weight-at-age (spawning, fishery, and survey), maturity-at-age, ageing
#' error, and an optional size-age transition matrix for length compositions.
#' All arrays must conform to the dimension structure stored in \code{sim_list}.
#'
#' @param sim_list Simulation list object returned by \code{\link{Setup_Sim_Dim}},
#'   which defines the dimension sizes used to validate all input arrays.
#' @param natmort_input Natural mortality array with dimensions
#'   \code{[n_pop × n_regions × n_yrs × n_ages × n_sexes × n_sims]}.
#'   Note: natural mortality is not season-specific and therefore lacks an
#'   \code{n_seas} dimension.
#' @param WAA_input Spawning weight-at-age array with dimensions
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]}.
#'   Used to compute spawning stock biomass.
#' @param WAA_fish_input Fishery weight-at-age array with dimensions
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims]}.
#'   Used to compute fishery biomass and catch in weight.
#' @param WAA_srv_input Survey weight-at-age array with dimensions
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]}.
#'   Used to compute survey biomass indices.
#' @param MatAA_input Maturity-at-age array with dimensions
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]}.
#'   Values should be proportions in \eqn{[0, 1]}. When \code{rec_lag = 0}
#'   (age-0 recruitment, set via \code{\link{Setup_Sim_Rec}}), maturity at
#'   the recruit age (the first age class) must be exactly \code{0} for all
#'   populations, regions, years, seasons, and sexes. An error is raised
#'   otherwise.
#' @param AgeingError_fish_input Optional fleet-specific ageing error for the
#'   simulated fishery fleets, dimensioned
#'   \code{[n_yrs × n_ages × n_obs_ages × n_fish_fleets × n_sims]}, or
#'   \code{NULL} (default) to give every fishery fleet \code{AgeingError_input}.
#' @param AgeingError_srv_input Optional fleet-specific ageing error for the
#'   simulated survey fleets, dimensioned
#'   \code{[n_yrs × n_ages × n_obs_ages × n_srv_fleets × n_sims]}, or
#'   \code{NULL} (default) to give every survey fleet \code{AgeingError_input}.
#' @param AgeingError_input Ageing error (age-length transition) array with
#'   dimensions \code{[n_yrs × n_model_ages × n_obs_ages × n_sims]}, where each
#'   \code{[n_model_ages × n_obs_ages]} slice is a row-stochastic matrix mapping
#'   true modeled ages to observed age bins. If \code{NULL} (default), an
#'   identity matrix is constructed for each year and simulation, which assumes
#'   that modeled and observed age bins are identical in number and alignment.
#'   \strong{If observed age bins are a subset of modeled ages} (e.g., observed
#'   ages 2-10 vs. modeled ages 1-10), the default identity matrix will cause a
#'   dimensional mismatch. In that case, supply a shifted identity matrix such as
#'   \code{diag(1, n_model_ages)[, obs_age_index]} to correctly drop or collapse
#'   model ages into observed bins.
#' @param SizeAgeTrans_fish_input,SizeAgeTrans_srv_input Optional size-age
#'   transition arrays per fleet, \code{[n_pop x n_regions x n_yrs x n_seas x
#'   n_lens x n_ages x n_sexes x n_fleets x n_sims]}, each read at that fleet's
#'   own timing. When supplied they are used for that fleet type in place of
#'   \code{SizeAgeTrans_input}; the self-test passes the fitted model's own keys
#'   here when growth was estimated.
#' @param SizeAgeTrans_input Size-age transition matrix array with dimensions
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_lens × n_ages × n_sexes × n_sims]}.
#'   Each slice maps age classes to length bins and should be column-stochastic
#'   (columns sum to 1). Only required when fitting length compositions;
#'   defaults to \code{NULL}.
#'
#' @return The input \code{sim_list} with the following fields added or updated:
#'   \describe{
#'     \item{\code{$natmort}}{Natural mortality array.}
#'     \item{\code{$WAA}}{Spawning weight-at-age array.}
#'     \item{\code{$WAA_fish}}{Fishery weight-at-age array.}
#'     \item{\code{$WAA_srv}}{Survey weight-at-age array.}
#'     \item{\code{$MatAA}}{Maturity-at-age array.}
#'     \item{\code{$AgeingError}}{Ageing error array (identity matrix if not supplied).}
#'     \item{\code{$SizeAgeTrans}}{Size-age transition array (only added if supplied).}
#'   }
#'
#'
#' @export Setup_Sim_Biologicals
#' @family Simulation Setup
Setup_Sim_Biologicals <- function(
                                  sim_list,
                                  natmort_input,
                                  WAA_input,
                                  WAA_fish_input,
                                  WAA_srv_input,
                                  MatAA_input,
                                  AgeingError_input = NULL,
                                  AgeingError_fish_input = NULL,
                                  AgeingError_srv_input = NULL,
                                  SizeAgeTrans_input = NULL,
                                  SizeAgeTrans_fish_input = NULL,
                                  SizeAgeTrans_srv_input = NULL
                                  ) {

  check_sim_dimensions(natmort_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_pop = sim_list$n_pop,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_sims = sim_list$n_sims, what = 'natmort_input')
  check_sim_dimensions(WAA_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,  n_pop = sim_list$n_pop,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_sims = sim_list$n_sims, what = 'WAA_input')
  check_sim_dimensions(WAA_fish_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_pop = sim_list$n_pop,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = 'WAA_fish_input')
  check_sim_dimensions(WAA_srv_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_pop = sim_list$n_pop,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets,
                       n_sims = sim_list$n_sims, what = 'WAA_srv_input')
  check_sim_dimensions(MatAA_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_pop = sim_list$n_pop,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_sims = sim_list$n_sims, what = 'MatAA_input')

  # Age-0 (rec_lag = 0) recruitment requires the recruit age class (the first
  # age) to be immature everywhere, since age-0 fish can't spawn the year they're
  # born. Requires Setup_Sim_Rec() to have run first so $rec_lag is already set.
  if(!is.null(sim_list$rec_lag) && sim_list$rec_lag == 0 && any(MatAA_input[,,,,1,,] != 0)) {
    stop("rec_lag = 0 (age-0 recruitment) requires MatAA_input to be zero at the recruit age (the first age class) for all populations, regions, years, seasons, and sexes, since age-0 fish cannot be mature.")
  }
  if(!is.null(SizeAgeTrans_input)) check_sim_dimensions(SizeAgeTrans_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_lens = sim_list$n_lens, n_seas = sim_list$n_seas, n_pop = sim_list$n_pop,
                                                        n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_sims = sim_list$n_sims, what = 'SizeAgeTrans_input')

  # output into list
  sim_list$natmort <- natmort_input
  sim_list$WAA <- WAA_input
  sim_list$WAA_fish <- WAA_fish_input
  sim_list$WAA_srv <- WAA_srv_input
  sim_list$MatAA <- MatAA_input
  if(!is.null(SizeAgeTrans_input)) sim_list$SizeAgeTrans <- SizeAgeTrans_input
  # keys per fleet, read at each fleet's own timing, take precedence over the shared one
  if(!is.null(SizeAgeTrans_fish_input)) {
    if(length(dim(SizeAgeTrans_fish_input)) != 9) stop("SizeAgeTrans_fish_input must be dimensioned n_pop x n_regions x n_yrs x n_seas x n_lens x n_ages x n_sexes x n_fish_fleets x n_sims")
    sim_list$SizeAgeTrans_fish <- SizeAgeTrans_fish_input
  }
  if(!is.null(SizeAgeTrans_srv_input)) {
    if(length(dim(SizeAgeTrans_srv_input)) != 9) stop("SizeAgeTrans_srv_input must be dimensioned n_pop x n_regions x n_yrs x n_seas x n_lens x n_ages x n_sexes x n_srv_fleets x n_sims")
    sim_list$SizeAgeTrans_srv <- SizeAgeTrans_srv_input
  }
  if(!is.null(AgeingError_input)) sim_list$AgeingError <- AgeingError_input
  else {
    # if null, create an identity matrix
    identity_AgeingError <- array(0, dim = c(sim_list$n_yrs, sim_list$n_ages, sim_list$n_ages, sim_list$n_sims))
    for(i in 1:sim_list$n_yrs) for(sim in 1:sim_list$n_sims) diag(identity_AgeingError[i,,,sim]) <- 1 # create identity matrix for each year
    sim_list$AgeingError <- identity_AgeingError
    warning("No ageing error matrix was provided. A default identity matrix was used, which assumes that the number and structure of modeled age bins exactly match the observed age bins. If the observed age composition data includes fewer age bins than the model (e.g., observed ages 2-10 while modeled ages are 1-10), this default assumption will cause a dimensional mismatch and potentially misalign the modeled and observed compositions. To avoid this, please provide an ageing error matrix of dimension n_model_ages x n_obs_ages that correctly maps modeled ages to observed age bins. For example, if observed ages are 2-10, supply a matrix that drops the first model age by using a shifted identity matrix: diag(1, 10)[, 2:10]. This will ensure the age bins are correctly aligned for likelihood calculations.")
  }

  # expand fleet-specific ageing error for the OM for backwards compatbility
  expand_sim_ae <- function(x, n_fleets, what) {
    shared <- sim_list$AgeingError
    d <- dim(shared)   # [n_yrs, n_ages, n_obs_ages, n_sims]
    if(is.null(x)) {
      out <- array(0, dim = c(d[1], d[2], d[3], n_fleets, d[4]))
      for(f in seq_len(n_fleets)) out[,,,f,] <- shared
      return(out)
    }
    if(length(dim(x)) != 5) stop(what, " must be dimensioned n_yrs x n_ages x n_obs_ages x n_fleets x n_sims")
    if(dim(x)[4] != n_fleets) stop(what, " has ", dim(x)[4], " fleets but the operating model has ", n_fleets, ".")
    if(!all(dim(x)[c(1,2,3,5)] == d)) stop(what, " must match AgeingError on years, ages, observed ages and sims.")
    return(x)
  }
  sim_list$AgeingError_fish <- expand_sim_ae(AgeingError_fish_input, sim_list$n_fish_fleets, "AgeingError_fish_input")
  sim_list$AgeingError_srv <- expand_sim_ae(AgeingError_srv_input, sim_list$n_srv_fleets, "AgeingError_srv_input")

  return(sim_list)

}

#' Map natural mortality parameters to a block structure
#'
#' Constructs the \code{M_blocks} index array and the \code{ln_M} factor map used
#' by the TMB/RTMB objective function to share or fix natural mortality parameters
#' across population, region, year, age, and sex dimensions. Each unique combination
#' of blocks is assigned a sequential integer ID; all cells within a block share the
#' same \code{ln_M} parameter.
#'
#' @param input_list Named list containing \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param M_spec Character string controlling whether \code{ln_M} is estimated or
#'   fixed. One of:
#'   \describe{
#'     \item{\code{"est_ln_M"}}{Freely estimate \code{ln_M} across all defined blocks.}
#'     \item{\code{"fix"}}{Fix all \code{ln_M} parameters by mapping them to \code{NA}.}
#'   }
#' @param M_popblk_spec_vals List of integer vectors assigning population indices to
#'   blocks, e.g., \code{list(1, 2)} for two population-specific blocks or
#'   \code{list(1:2)} for a single shared block.
#' @param M_regionblk_spec_vals List of integer vectors assigning region indices to
#'   blocks, e.g., \code{list(1:3, 4:5)} for two region blocks.
#' @param M_yearblk_spec_vals List of integer vectors assigning year indices to
#'   blocks, e.g., \code{list(1:10, 11:30)} for two time periods.
#' @param M_ageblk_spec_vals List of integer vectors assigning age indices to
#'   blocks, e.g., \code{list(1:5, 6:10)} for two age groups.
#' @param M_sexblk_spec_vals List of integer vectors assigning sex indices to
#'   blocks. Use \code{list(1:2)} for a sex-invariant block or \code{list(1, 2)}
#'   for sex-specific mortality.
#'
#' @return The input \code{input_list} with two fields updated:
#'   \describe{
#'     \item{\code{$map$ln_M}}{Factor vector of length equal to \code{prod(dim(par$ln_M))}.
#'       Each element is an integer estimation index when \code{M_spec = "est_ln_M"},
#'       or \code{NA} when \code{M_spec = "fix"}.}
#'     \item{\code{$data$M_blocks}}{Integer array of dimensions
#'       \code{[n_pop × n_regions × n_years × n_ages × n_sexes]} mapping each
#'       population-region-year-age-sex cell to its corresponding \code{ln_M}
#'       parameter index.}
#'   }
#'
#' @keywords internal
do_natmort_mapping <- function(input_list,
                         M_spec,
                         M_popblk_spec_vals,
                         M_regionblk_spec_vals,
                         M_yearblk_spec_vals,
                         M_ageblk_spec_vals,
                         M_sexblk_spec_vals) {

  # Validate options
  if(!M_spec %in% c('est_ln_M', 'fix')) stop("M_spec needs to be specified as either est_ln_M or fix")

  # set up whether fixing M or estimating
  if(M_spec == 'est_ln_M') input_list$map$ln_M <- factor(1:length(input_list$par$ln_M))
  if(M_spec == 'fix') input_list$map$ln_M <- factor(rep(NA, length(input_list$par$ln_M)))

  # create array for blocks
  M_blocks <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))

  # loop through to get counters for blocking structure for indexing
  counter <- 1
  for(popblk in 1:length(M_popblk_spec_vals)) {
    map_p <- M_popblk_spec_vals[[popblk]]

    for (regionblk in 1:length(M_regionblk_spec_vals)) {
      map_r <- M_regionblk_spec_vals[[regionblk]]

      for (yearblk in 1:length(M_yearblk_spec_vals)) {
        map_y <- M_yearblk_spec_vals[[yearblk]]

        for (ageblk in 1:length(M_ageblk_spec_vals)) {
          map_a <- M_ageblk_spec_vals[[ageblk]]

          for (sexblk in 1:length(M_sexblk_spec_vals)) {
            map_s <- M_sexblk_spec_vals[[sexblk]]

            # Assign the current counter to this block
            M_blocks[map_p, map_r, map_y, map_a, map_s] <- counter
            counter <- counter + 1

          } # end sexblk
        } # end ageblk
      } # end yearblk
    } # end regionblk
  }

  collect_message("Natural Mortality specified as: ", M_spec)
  collect_message("Natural Mortality Population Blocks is specified as: ", length(M_popblk_spec_vals))
  collect_message("Natural Mortality Region Blocks is specified as: ", length(M_regionblk_spec_vals))
  collect_message("Natural Mortality Year Blocks is specified as: ", length(M_yearblk_spec_vals))
  collect_message("Natural Mortality Age Blocks is specified as: ", length(M_ageblk_spec_vals))
  collect_message("Natural Mortality Sex Blocks is specified as: ", length(M_sexblk_spec_vals))

  input_list$data$M_blocks <- M_blocks

  return(input_list)
}

#' Set up mapping for growth
#'
#' Builds the estimation maps for the growth parameters, for the deviations of
#' any parameter that varies over time, and for the semi-parametric surface on
#' mean length at age. Called by \code{\link{Setup_Mod_Biologicals}} once the
#' parameter list is populated, since every map is dimensioned off its parameter.
#'
#' @param input_list List containing data, parameter, and map lists, with the
#'   growth parameters already populated.
#' @param growth_spec Character. How the growth parameters are estimated, one of
#'   \code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"} or \code{"fix"}.
#' @param growth_fix Logical vector, one entry per growth parameter, naming any
#'   held at its starting value while the others are estimated.
#' @param tv_vals Integer vector, one entry per growth parameter, of the time
#'   variation each carries (0 none, 1 iid, 2 random walk).
#' @param tv_active Matrix \code{[n_years x n_gpars]} of ones in the years each
#'   parameter's deviations are estimated in.
#' @param growth_tv_spec Character. How the deviation series are shared across
#'   regions and sexes.
#' @param growth_tv_sigma_spec Character. \code{"fix"} holds the process error
#'   standard deviations, \code{"est"} estimates them.
#' @param semipar_val Integer code of the semi-parametric form (0 none, 1 iid,
#'   2 random walk, 3 \code{3dmarg}, 4 \code{3dcond}, 5 \code{2dar1}).
#' @param growth_semipar_spec Character. Whether the surface's process error
#'   parameters are estimated.
#' @param semipar_age_idx,semipar_yr_idx Integer vectors of the age and year
#'   indices the surface is estimated over.
#'
#' @return The input \code{input_list} with \code{$map$ln_growth_pars},
#'   \code{$map$ln_growth_devs}, \code{$map$growth_pe_pars} and
#'   \code{$map$ln_growth_semipar_devs} set, along with the
#'   \code{$data$map_ln_growth_devs} and
#'   \code{$data$map_ln_growth_semipar_devs} mirrors the deviation penalties
#'   read.
#'
#' @keywords internal
do_growth_mapping <- function(input_list,
                              growth_spec,
                              growth_fix,
                              tv_vals,
                              tv_active,
                              growth_tv_spec,
                              growth_tv_sigma_spec,
                              semipar_val,
                              growth_semipar_spec,
                              semipar_age_idx,
                              semipar_yr_idx) {

  n_pop <- input_list$data$n_pop
  n_regions <- input_list$data$n_regions
  n_sexes <- input_list$data$n_sexes
  n_gpars <- dim(input_list$par$ln_growth_pars)[4]

  # Growth parameters ---------------------------------------------------------
  # one index per estimated cell, shared across regions and sexes as growth_spec
  # says, NA wherever the parameter is held
  map_growth <- array(NA_integer_, dim = dim(input_list$par$ln_growth_pars))
  counter <- 1

  for(k in 1:n_gpars) {

    if(growth_spec == "fix" || isTRUE(growth_fix[k])) next

    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {

          share_r <- growth_spec %in% c("est_shared_r", "est_shared_r_s") && r > 1
          share_s <- growth_spec %in% c("est_shared_s", "est_shared_r_s") && s > 1

          if(share_r) {
            map_growth[p, r, s, k] <- map_growth[p, 1, s, k]
          } else if(share_s) {
            map_growth[p, r, s, k] <- map_growth[p, r, 1, k]
          } else {
            map_growth[p, r, s, k] <- counter
            counter <- counter + 1
          }

        } # end s loop
      } # end r loop
    } # end p loop
  } # end k loop

  input_list$map$ln_growth_pars <- factor(map_growth)

  # Time-varying growth parameters --------------------------------------------
  # one deviation per year a parameter is active in, and one process error
  # standard deviation per varying parameter in the first stream of the shared
  # process error array, both shared as growth_tv_spec says
  map_devs <- array(NA_integer_, dim = dim(input_list$par$ln_growth_devs))
  map_pe <- array(NA_integer_, dim = dim(input_list$par$growth_pe_pars))
  counter <- 1
  counter_pe <- 1

  for(k in which(tv_vals > 0)) {
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {

          share_r <- growth_tv_spec %in% c("est_shared_r", "est_shared_r_s") && r > 1
          share_s <- growth_tv_spec %in% c("est_shared_s", "est_shared_r_s") && s > 1

          if(share_r) {

            map_devs[p, r, , k, s] <- map_devs[p, 1, , k, s]
            map_pe[p, r, k, s, 1] <- map_pe[p, 1, k, s, 1]

          } else if(share_s) {

            map_devs[p, r, , k, s] <- map_devs[p, r, , k, 1]
            map_pe[p, r, k, s, 1] <- map_pe[p, r, k, 1, 1]

          } else {

            for(y in which(tv_active[, k] == 1)) {
              map_devs[p, r, y, k, s] <- counter
              counter <- counter + 1
            } # end y loop

            if(growth_tv_sigma_spec == "est") {
              map_pe[p, r, k, s, 1] <- counter_pe
              counter_pe <- counter_pe + 1
            }

          }

        } # end s loop
      } # end r loop
    } # end p loop
  } # end k loop

  input_list$map$ln_growth_devs <- factor(map_devs)
  input_list$data$map_ln_growth_devs <- array(as.numeric(input_list$map$ln_growth_devs), dim = dim(map_devs))

  # Semi-parametric growth surface --------------------------------------------
  # a deviation in every age and year the surface is estimated over, and process
  # error parameters in the second stream's slots the form reads
  map_semi <- array(NA_integer_, dim = dim(input_list$par$ln_growth_semipar_devs))

  if(semipar_val > 0) {

    counter <- 1

    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          for(y in semipar_yr_idx) {
            for(a in semipar_age_idx) {

              map_semi[p, r, y, a, s] <- counter
              counter <- counter + 1

            } # end a loop
          } # end y loop
        } # end s loop
      } # end r loop
    } # end p loop

    if(growth_semipar_spec == "est") {

      # the correlated forms read three correlations and a scale, iid and the
      # random walk one sigma per age
      slots <- if(semipar_val %in% 1:2) semipar_age_idx else if(semipar_val %in% 3:4) 1:4 else c(1, 2, 4)

      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            for(k in slots) {

              map_pe[p, r, k, s, 2] <- counter_pe
              counter_pe <- counter_pe + 1

            } # end k loop
          } # end s loop
        } # end r loop
      } # end p loop
    }
  }

  input_list$map$ln_growth_semipar_devs <- factor(map_semi)
  input_list$map$growth_pe_pars <- factor(map_pe)
  input_list$data$map_ln_growth_semipar_devs <- array(as.numeric(input_list$map$ln_growth_semipar_devs), dim = dim(map_semi))

  return(input_list)
}

#' Set up biological inputs for the estimation model
#'
#' Populates \code{input_list} with biological arrays and parameter structures
#' needed by the TMB/RTMB objective function: weight-at-age (spawning, fishery,
#' and survey), maturity-at-age, ageing error, the size-age transition matrix
#' (optional), small constants for numerical stability, and the natural mortality
#' block structure and mapping. Called after \code{\link{Setup_Mod_Dim}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map}, and
#'   \code{$verbose} sublists, as returned by \code{\link{Setup_Mod_Dim}}.
#' @param WAA Numeric array of spawning weight-at-age with dimensions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]}.
#'   Used to compute spawning stock biomass. Also serves as the fallback for
#'   \code{WAA_fish} and \code{WAA_srv} when those are \code{NULL}.
#' @param WAA_fish Numeric array of fishery weight-at-age with dimensions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), \code{WAA} is broadcast across all fishery fleets.
#' @param WAA_srv Numeric array of survey weight-at-age with dimensions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]}.
#'   If \code{NULL} (default), \code{WAA} is broadcast across all survey fleets.
#' @param MatAA Numeric array of maturity-at-age proportions (\eqn{\in [0,1]}) with
#'   dimensions \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]}.
#'   When \code{rec_lag = 0} (age-0 recruitment, set via \code{\link{Setup_Mod_Rec}}),
#'   maturity at the recruit age (the first age class) must be exactly \code{0}
#'   for all populations, regions, years, seasons, and sexes, an error is
#'   raised otherwise. Requires \code{\link{Setup_Mod_Rec}} to have been called
#'   first so \code{rec_lag} is already set.
#' @param addtocomp \strong{Deprecated here}, pass it to \code{\link{Setup_Mod_Weighting}}
#'   instead, which now owns this constant along with every other likelihood weight.
#'   Still accepted for backward compatibility: if supplied, it is forwarded to
#'   \code{Setup_Mod_Weighting} with a message rather than applied here directly.
#'   (Small constant added to composition proportions before likelihood evaluation
#'   to avoid \code{log(0)}; default \code{1e-3} in \code{Setup_Mod_Weighting}.
#'   Ignored when a logistic-normal likelihood is specified, as that family handles
#'   zeros internally.)
#' @param comp_const_obs \strong{Deprecated here}, pass it to
#'   \code{\link{Setup_Mod_Weighting}} instead. Still accepted for backward
#'   compatibility (forwarded with a message). Integer switch (\code{0} or
#'   \code{1}) controlling where \code{addtocomp} is applied in the multinomial
#'   likelihood, not a constant to be tuned. \code{1} (default in
#'   \code{Setup_Mod_Weighting}) adds it to the observed proportions that weight
#'   the multinomial as well as inside the logarithms, so the likelihood is
#'   stationary exactly at \code{pred = obs}. \code{0} weights by the raw
#'   observed proportions. The Dirichlet-multinomial sanity check that used to
#'   read it here (inside \code{Setup_Mod_FishIdx_and_Comps}/
#'   \code{Setup_Mod_SrvIdx_and_Comps}) now runs inside \code{Setup_Mod_Weighting}
#'   once the final value is known.
#' @param addtofishidx \strong{Deprecated here}, pass it to
#'   \code{\link{Setup_Mod_Weighting}} instead. Still accepted for backward
#'   compatibility (forwarded with a message). Small constant added to fishery
#'   indices; default \code{1e-4} in \code{Setup_Mod_Weighting}.
#' @param addtosrvidx \strong{Deprecated here}, pass it to
#'   \code{\link{Setup_Mod_Weighting}} instead. Still accepted for backward
#'   compatibility (forwarded with a message). Small constant added to survey
#'   indices; default \code{1e-4} in \code{Setup_Mod_Weighting}.
#' @param addtotag \strong{Deprecated here}, pass it to
#'   \code{\link{Setup_Mod_Weighting}} instead. Still accepted for backward
#'   compatibility (forwarded with a message). Small constant added to tag
#'   recovery observations; default \code{1e-10} in \code{Setup_Mod_Weighting}.
#' @param AgeingError Ageing error (age-age transition) array mapping true modeled ages
#'   to observed age bins. Each row is one model age's share across the observed
#'   bins and sums to one, or to zero to drop that model age from the
#'   observations. This is the age-axis twin of \code{LenBinMap}: the likelihood
#'   applies the two identically and validates them identically, so read either
#'   one for the other. It changes which bins the compositions are recorded on;
#'   to leave observed bins out of the likelihood without changing the bins
#'   themselves, use the \code{*_bins} arguments instead. Accepted forms:
#'   \describe{
#'     \item{2D matrix \code{[n_model_ages × n_obs_ages]}}{Time-invariant ageing error;
#'       replicated internally across all years.}
#'     \item{3D array \code{[n_years × n_model_ages × n_obs_ages]}}{Time-varying ageing error.}
#'     \item{\code{NULL} (default)}{An identity matrix is constructed, assuming modeled
#'       and observed age bins are identical. If observed bins are a subset of modeled
#'       ages (e.g., observed ages 2-10 vs. modeled ages 1-10), supply a shifted
#'       identity matrix such as \code{diag(1, n_model_ages)[, obs_age_index]} to
#'       avoid a dimensional mismatch.}
#'   }
#' @param AgeingError_fish Optional fleet-specific ageing error for the fishery
#'   fleets, for when the fleets do not read ages the same way. Accepted forms:
#'   a 3D array \code{[n_model_ages × n_obs_ages × n_fish_fleets]} for a
#'   time-invariant matrix per fleet, a 4D array
#'   \code{[n_years × n_model_ages × n_obs_ages × n_fish_fleets]} for a
#'   time-varying one, or \code{NULL} (default), which gives every fishery fleet
#'   the shared \code{AgeingError}. Each fleet's slice is validated the same way
#'   \code{AgeingError} is, and every fleet must land on the same observed age
#'   bins, since the observed composition arrays carry one age dimension shared
#'   across fleets.
#' @param AgeingError_srv Optional fleet-specific ageing error for the survey
#'   fleets, in the same forms as \code{AgeingError_fish}, with
#'   \code{n_srv_fleets} in place of \code{n_fish_fleets}. \code{NULL}
#'   (default) gives every survey fleet the shared \code{AgeingError}.
#' @param Use_M_prior Integer flag to apply a lognormal prior on natural mortality.
#'   \code{0} = no prior (default); \code{1} = apply prior.
#' @param M_prior Data frame of prior hyperparameters for natural mortality, with one
#'   row per unique block combination. Required columns:
#'   \describe{
#'     \item{\code{popblk}, \code{regionblk}, \code{yearblk}, \code{ageblk}, \code{sexblk}}{Block indices identifying which parameter the prior applies to.}
#'     \item{\code{mu}}{Prior mean in natural (untransformed) space.}
#'     \item{\code{sd}}{Prior standard deviation.}
#'   }
#'   Example for a single shared prior:
#'   \preformatted{M_prior <- data.frame(
#'     popblk = 1, regionblk = 1, yearblk = 1,
#'     ageblk = 1, sexblk = 1,
#'     mu = 0.085, sd = 0.05
#'   )}
#'   Only used when \code{Use_M_prior = 1}.
#' @param fit_lengths Integer flag for fitting length compositions. \code{0} = no
#'   (default); \code{1} = yes. Requires a valid \code{SizeAgeTrans} array.
#' @param SizeAgeTrans Numeric array of size-at-age transition probabilities
#'   (column-stochastic; each age column sums to 1) with dimensions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_ages × n_sexes]}.
#'   Required when \code{fit_lengths = 1}; ignored otherwise. The shared key
#'   every fleet reads unless \code{SizeAgeTrans_fish}/\code{SizeAgeTrans_srv}
#'   override it for that fleet type.
#' @param SizeAgeTrans_fish,SizeAgeTrans_srv Optional per-fleet size-at-age
#'   transition arrays, dimensioned like \code{SizeAgeTrans} with an added
#'   trailing fleet dimension (\code{n_fish_fleets}/\code{n_srv_fleets}).
#'   \code{NULL} (default) reads every fleet's key from the shared
#'   \code{SizeAgeTrans}. Only meaningful with \code{growth_model = "none"};
#'   a growth model already derives one key per fleet, at that fleet's own
#'   timing, and rejects these to avoid mixing two sources for the same key.
#'   This is the fixed-data counterpart of \code{Setup_Sim_Biologicals}'s
#'   \code{SizeAgeTrans_fish_input}/\code{SizeAgeTrans_srv_input}, and of
#'   \code{WAA_fish}/\code{WAA_srv} overriding the shared \code{WAA}.
#' @param do_caal Integer flag for building the joint arrays at length and age.
#'   \code{0} = no (default); \code{1} = yes. Requires \code{fit_lengths = 1}.
#'   Turning this on adds \code{Fish_caal}, \code{Fish_caal_discard} and \code{Srv_caal} to the
#'   report, holding predicted retained catch, discards and survey index jointly
#'   by length and age.
#' @param growth_model Character. \code{"none"} (default) keeps \code{SizeAgeTrans}
#'   and the weight-at-age arrays as data. \code{"vb_schnute"} builds the size-age
#'   transition from estimable von Bertalanffy parameters in Schnute's form:
#'   length \code{L1} at reference age \code{growth_A1}, length \code{L2} at
#'   \code{growth_A2}, rate \code{K}, and CVs of length at age \code{CV1} and
#'   \code{CV2} at the two reference ages. Growth below \code{growth_A1} is linear
#'   from \code{growth_L0} at age zero, the CV interpolates between the two
#'   references, and the plus group carries an adjustment for fish older than
#'   the accumulator age. \code{"richards"} is the same curve with a sixth
#'   parameter, the Richards coefficient \code{rho}, applied to the lengths
#'   raised to that power (\code{rho = 1} recovers the von Bertalanffy form).
#'   Requires \code{fit_lengths = 1}; \code{SizeAgeTrans} is then ignored and
#'   may be \code{NA}.
#' @param growth_spec Character. How the growth parameters are estimated:
#'   \code{"est_all"} (default, one set per population, region and sex),
#'   \code{"est_shared_r"} (shared across regions), \code{"est_shared_s"}
#'   (shared across sexes), \code{"est_shared_r_s"} (one set per population),
#'   or \code{"fix"}.
#' @param growth_fix Logical vector, one entry per growth parameter, naming
#'   which of L1, L2, K, CV1, CV2 (and rho) stay at their starting values
#'   whatever \code{growth_spec} says.
#' @param growth_tv_model Time variation of the growth parameters. \code{NULL}
#'   (default) holds every parameter constant. Otherwise a character vector
#'   naming a structure per parameter, either of length \code{n_gpars} in the
#'   parameter order or named by parameter (\code{L1}, \code{L2}, \code{K},
#'   \code{CV1}, \code{CV2}, \code{rho}) with the rest constant, each one of
#'   \code{"none"}, \code{"iid"} (independent annual deviations) or
#'   \code{"rw"} (a random walk). A varying parameter gets a deviation series
#'   \code{ln_growth_devs} and a log sigma in the first stream of
#'   \code{growth_pe_pars}.
#' @param growth_tv_years Years the deviations are active in, calendar years.
#'   \code{NULL} (default) for every model year, a vector applied to every
#'   varying parameter, or a list named by parameter. Deviations outside the
#'   range are held at zero.
#' @param growth_tv_link Character, the scale a deviation enters on.
#'   \code{"log"} (default) multiplies the parameter by \eqn{e^{\delta}};
#'   \code{"logit"} keeps it inside \code{growth_par_bounds},
#'   \eqn{P_y = lo + (hi - lo)\,\mathrm{logit}^{-1}(\mathrm{logit}((P - lo)/(hi - lo)) + \delta_y)},
#'   so the parameter approaches a bound however large the deviation instead of
#'   crossing it.
#' @param growth_par_bounds Matrix \code{[n_gpars x 2]} of lower and upper
#'   bounds, natural scale, required under the logit link.
#' @param growth_tv_sigma_spec Character, \code{"fix"} (default) holds the
#'   process error standard deviations of the deviations at their starting
#'   values, \code{"est"} estimates them. Both read the first stream of
#'   \code{growth_pe_pars}, one slot per growth parameter.
#' @param growth_tv_spec Character, how the deviations are shared across
#'   strata, with the same vocabulary as \code{growth_spec}: \code{"est_all"}
#'   (default), \code{"est_shared_r"}, \code{"est_shared_s"} or
#'   \code{"est_shared_r_s"}.
#' @param growth_tv_type Character. \code{"curve"} (default) reads every
#'   year's size at age off that year's curve. \code{"cohort"} carries size at
#'   age forward cohort by cohort: each year every cohort grows by the increment
#'   the current year's parameters imply from the size it reached, ages still in
#'   the linear phase keep the length at \code{growth_A1} their birth year's
#'   parameters gave them, the first age past \code{growth_A1} is placed on the
#'   current year's curve, and the plus group's size blends the cohort entering it
#'   with the fish already there by their numbers at age. The CV at age is
#'   then held at the first year's sizes. The propagation starts in the first
#'   year any deviation is active; every earlier year sits on the first year's
#'   curve.
#' @param growth_rw_init_sigma Standard deviation given to the first year of a
#'   random walk on a growth parameter, as \code{srvsel_rw_init_sigma} for
#'   selectivity. Default \code{5}.
#' @param growth_semipar Character. Semi-parametric growth: a year-by-age
#'   surface of deviations on mean length at age, multiplying the parametric
#'   curve, so the curve stays the parametric part and the deviations hold departures
#'   from it. \code{"none"} (default) keeps growth purely parametric; otherwise
#'   one of \code{"iid"}, \code{"rw"} (a random walk over years within an age),
#'   \code{"3dmarg"} or \code{"3dcond"} (a three-dimensional Gaussian Markov
#'   random field over age, year and cohort, on the marginal or conditional
#'   variance), or \code{"2dar1"} (a separable first-order autoregression over
#'   ages and years). The same process error forms the selectivity deviations
#'   use, so a growth surface and a selectivity surface are scored the same way.
#'   The spread at age follows the deviated mean, which leaves the coefficient of
#'   variation at age to the parametric part.
#' @param growth_semipar_spec Character, whether the second stream of
#'   \code{growth_pe_pars} is estimated. Whether the process error
#'   hyperparameters are estimated (\code{"est"}) or held at their starting
#'   values (\code{"fix"}, the default). The deviations themselves are always
#'   estimated.
#' @param growth_semipar_ages Ages the deviations are estimated over, as ages
#'   (not indices). \code{NULL} (default) uses every age. Ages outside the set
#'   are held at zero, which is how a surface is restricted to the ages the
#'   length data actually inform.
#' @param growth_semipar_years Years the deviations are estimated over, calendar
#'   years. \code{NULL} (default) uses every year.
#' @param LenBinMap Optional matrix \code{[n_lens x n_obs_lens]} mapping the
#'   model's length bins onto the bins the length compositions are recorded on,
#'   for compositions on coarser bins than the model carries (a population of
#'   1 cm bins fit to 5 cm compositions, say). Observed length compositions are
#'   then dimensioned by \code{n_obs_lens} and the expected compositions are
#'   mapped through it inside the likelihood. This is the length-axis twin of
#'   \code{AgeingError}: the likelihood applies the two identically and
#'   validates them identically, so read either one for the other. Each row is
#'   one model bin's share across the observed bins and sums to one, or to zero
#'   to drop that model bin from the observations. It changes which bins the
#'   compositions are recorded on; to leave observed bins out of the likelihood
#'   without changing the bins themselves, use the \code{*LenComps_bins}
#'   arguments instead. \code{NULL} (default) fits the compositions on the model
#'   bins.
#' @param growth_A1,growth_A2 Reference ages for \code{L1} and \code{L2}.
#'   \code{growth_A2 = "Linf"} instead makes \code{L2} the asymptotic length
#'   itself, with no second reference age to solve it from.
#' @param growth_len_lower Numeric vector of the lower edges of the length bins.
#'   \code{lens} in \code{Setup_Mod_Dim} are bin midpoints; the key is built on
#'   the edges.
#' @param growth_L0 Length at age zero anchoring the linear phase. Defaults to
#'   \code{growth_len_lower[1]}.
#' @param growth_cv_type Character, \code{"len"} (default) interpolates the CV
#'   on mean length between \code{L1} and \code{L2}, \code{"age"} on age.
#' @param growth_sd_type Character, \code{"cv"} (default) scales the mean by the
#'   CV parameters, \code{"sd"} reads them as standard deviations.
#' @param growth_dist Character, \code{"normal"} (default) or \code{"lognormal"}
#'   distribution of length at age.
#' @param growth_plus_group Character. \code{"mixture"} (default) takes the plus
#'   group's mean length as the survivorship-weighted mixture of the ages it
#'   holds, their numbers declining at an assumed 0.2 per year and their length
#'   rising from the curve at the accumulator age to the asymptote; \code{"curve"}
#'   reads the curve at the accumulator age.
#' @param waa_model Character. Where weight at age comes from.
#'   \code{"data"} (default) reads \code{WAA}, \code{WAA_fish} and
#'   \code{WAA_srv} from the arguments of the same name. \code{"wt_len"} builds
#'   them from the size-age key and the weight-length relationship
#'   \eqn{W = a L^b} applied at the bin midpoints, so weight at age carries the
#'   spread of length at age rather than being the weight of the mean length;
#'   the spawning weight uses the key at spawning time and each fleet's weight
#'   the key at that fleet's timing, \code{t_fish} or \code{t_srv}. Under
#'   \code{"wt_len"}, \code{WAA} may be \code{NULL}, and reference point and
#'   projection code still read \code{data$WAA}, so copy the reported arrays
#'   into the data list before calling them.
#' @param wt_len_pars Weight-length parameters \eqn{a, b} in \eqn{W = a L^b},
#'   a vector of two or an array \code{[n_pop x n_regions x n_sexes x 2]}.
#'   Required when \code{waa_model = "wt_len"}.
#' @param M_spec Character string controlling natural mortality estimation. One of:
#'   \describe{
#'     \item{\code{"est_ln_M"} (default)}{Estimate \code{ln_M} across the defined blocks.}
#'     \item{\code{"fix"}}{Fix mortality to \code{Fixed_natmort}; \code{ln_M} parameters
#'       are mapped to \code{NA} and not passed to the optimizer.}
#'   }
#' @param NAA_re Character. State-space numbers at age: the log numbers themselves
#'   become parameters for ages two and older, including the plus group, and the
#'   deterministic mortality and ageing step becomes a prediction they are scored
#'   against. \code{"none"} (default) leaves the numbers deterministic. Otherwise
#'   one of \code{"iid"},
#'   \code{"1dar1_a"} (an autoregression over ages, years independent),
#'   \code{"1dar1_y"} (an autoregression over years, ages independent),
#'   \code{"2dar1"} (a separable
#'   autoregression over ages and years), or \code{"3dcond"} and \code{"3dmarg"}
#'   (a three-dimensional Gaussian Markov random field over age, year and cohort,
#'   on the conditional or the marginal variance). These are the same process
#'   error forms the selectivity and growth surfaces use, so a state on the
#'   numbers and a surface on selectivity are scored the same way. Age one belongs
#'   to \code{ln_RecDevs} and year one at ages two and older to
#'   \code{ln_InitDevs}, so the three partition the numbers at age rather than
#'   overlapping. The state is the log numbers, not a deviation, so \code{ln_NAA}
#'   starts from an equilibrium decay from \eqn{R_0} at the mean \eqn{M}, split
#'   over regions and sexes, with the plus group accumulated. Pass \code{ln_NAA}
#'   through \code{starting_values} to start somewhere else.
#' @param NAA_re_ages Ages the state is estimated over, as ages rather than
#'   indices. \code{NULL} (default) uses every age from the second onward. Must be
#'   a contiguous run.
#' @param NAA_re_years Calendar years the state is estimated over. \code{NULL}
#'   (default) uses every year from the second onward. Must be a contiguous run.
#' @param NAA_re_region Character. Correlation across regions, composed with
#'   whatever \code{NAA_re} gives over the age and year grid. \code{"iid"} (the
#'   default) leaves regions independent; \code{"us"} estimates an unstructured
#'   correlation, \eqn{n_r(n_r-1)/2} parameters, placing no shape on how regions
#'   covary. Independence is the default deliberately: a flexible correlation
#'   manufactures structure from independent data far more readily than it misses
#'   real structure.
#' @param NAA_re_region_spec Character controlling how the region correlations are
#'   shared, following the package's spec strings: \code{"est_all"} (the default)
#'   gives a free correlation matrix per population and sex, \code{"est_shared_p"}
#'   and \code{"est_shared_s"} share it over one of those margins,
#'   \code{"est_shared_p_s"} gives a single matrix for the whole model, and
#'   \code{"fix"} holds them all.
#' @param NAA_re_pop,NAA_re_sex Character. Correlation across populations and
#'   across sexes, composed with the region, age and year structures the same way.
#'   \code{"iid"} (the default) leaves them independent; \code{"us"} estimates an
#'   unstructured correlation. Both are global to the model rather than varying
#'   over the other margins, so a two-sex model spends exactly one parameter on
#'   \code{NAA_re_sex = "us"}. Note that a survival correlation between sexes is
#'   weakly identified in most configurations: spawning biomass reads the first
#'   sex while abundance indices sum over sexes, so an antisymmetric shock moves
#'   spawning biomass at almost no cost in the indices, and only sex-structured
#'   compositions push back on it.
#' @param NAA_sigma_spec Character, whether the process error standard deviations
#'   are estimated (\code{"est"}, the default) or held at their starting values
#'   (\code{"fix"}). The states themselves are always estimated.
#' @param NAA_sigma_popblk_spec,NAA_sigma_regionblk_spec,NAA_sigma_yearblk_spec,NAA_sigma_ageblk_spec,NAA_sigma_sexblk_spec
#'   Blocking for the process error standard deviation, each either
#'   \code{"constant"} (the default) or a list of integer vectors assigning
#'   indices to blocks, exactly as the \code{M_*blk_spec} arguments do. Blocking
#'   shares a standard deviation; it never removes a cell from the state. Only
#'   \code{NAA_re = "iid"} admits a standard deviation that varies over years or
#'   ages: every other form is separable or Markov in a margin, so it carries one
#'   standard deviation per population, region and sex.
#' @param Fixed_natmort Numeric array of fixed natural mortality rates with dimensions
#'   \code{[n_pop × n_regions × n_years × n_ages × n_sexes]}. Note the absence of an
#'   \code{n_seas} dimension. Required when \code{M_spec = "fix"}; ignored otherwise.
#' @param M_popblk_spec Blocking structure for \code{ln_M} across populations. Either
#'   \code{"constant"} (default; single shared value) or a list of integer index
#'   vectors defining population groups, e.g., \code{list(1, 2)} for
#'   population-specific M.
#' @param M_regionblk_spec Blocking structure across regions. Either \code{"constant"}
#'   (default) or a list of integer index vectors, e.g., \code{list(1:3, 4:5)}.
#' @param M_yearblk_spec Blocking structure across years. Either \code{"constant"}
#'   (default) or a list of integer index vectors, e.g., \code{list(1:10, 11:30)}.
#' @param M_ageblk_spec Blocking structure across ages. Either \code{"constant"}
#'   (default) or a list of integer index vectors, e.g., \code{list(1:5, 6:10)}.
#' @param M_sexblk_spec Blocking structure across sexes. Either \code{"constant"}
#'   (default; shared across sexes) or a list of integer index vectors, e.g.,
#'   \code{list(1, 2)} for sex-specific M.
#' @param ... Optional starting value overrides passed by name. Currently recognized:
#'   \describe{
#'     \item{\code{ln_M}}{Array of log-scale starting values for natural mortality,
#'       dimensioned \code{[n_popblks × n_regionblks × n_yearblks × n_ageblks × n_sexblks]}.
#'       Defaults to \code{log(0.5)} for all blocks if not supplied.}
#'     \item{\code{ln_growth_pars}}{Array of log-scale starting values for the
#'       growth parameters, dimensioned \code{[n_pop × n_regions × n_sexes × n_gpars]}
#'       in the order \code{L1, L2, K, CV1, CV2} and, under the Richards form,
#'       \code{rho}. Defaults to the ends of the length bins with a rate of
#'       \code{0.15} and CVs of \code{0.1}, so supply your own for any real
#'       model.}
#'     \item{\code{growth_pe_pars}}{Array of process error starting values for
#'       both growth deviation streams, dimensioned
#'       \code{[n_pop × n_regions × max(4, n_ages, n_gpars) × n_sexes × 2]}.
#'       The first stream holds one log sigma per growth parameter for the
#'       time-varying deviations; the second holds the semi-parametric surface's
#'       correlations by age, year and cohort in slots one to three and a log
#'       scale in slot four for the correlated forms, or one log sigma per age
#'       for \code{"iid"} and \code{"rw"}. Defaults to \code{log(0.1)} for the
#'       first stream and \code{log(0.05)} with correlations of \code{0.3} for
#'       the second. Slots a form does not read are mapped off.}
#'   }
#'   All \code{...} arguments are silently ignored when \code{M_spec = "fix"}.
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and \code{$map}
#'   sublists updated. Key additions include \code{$data$WAA}, \code{$data$WAA_fish},
#'   \code{$data$WAA_srv}, \code{$data$MatAA}, \code{$data$AgeingError},
#'   \code{$data$M_blocks}, \code{$par$ln_M}, and \code{$map$ln_M}.
#'
#' @export Setup_Mod_Biologicals
#' @family Model Setup
Setup_Mod_Biologicals <- function(input_list,
                                  WAA,
                                  WAA_fish = NULL,
                                  WAA_srv = NULL,
                                  MatAA,
                                  addtocomp = NULL,
                                  comp_const_obs = NULL,
                                  addtofishidx = NULL,
                                  addtosrvidx = NULL,
                                  addtotag = NULL,
                                  AgeingError = NULL,
                                  AgeingError_fish = NULL,
                                  AgeingError_srv = NULL,
                                  Use_M_prior = 0,
                                  M_prior = NA,
                                  fit_lengths = 0,
                                  SizeAgeTrans = NA,
                                  SizeAgeTrans_fish = NULL,
                                  SizeAgeTrans_srv = NULL,
                                  do_caal = 0,
                                  growth_model = "none",
                                  growth_spec = "est_all",
                                  growth_fix = NULL,
                                  growth_tv_model = NULL,
                                  growth_tv_years = NULL,
                                  growth_tv_link = "log",
                                  growth_par_bounds = NULL,
                                  growth_tv_sigma_spec = "fix",
                                  growth_tv_spec = "est_all",
                                  growth_tv_type = "curve",
                                  growth_rw_init_sigma = 5,
                                  growth_semipar = "none",
                                  growth_semipar_spec = "fix",
                                  growth_semipar_ages = NULL,
                                  growth_semipar_years = NULL,
                                  LenBinMap = NULL,
                                  growth_A1 = NULL,
                                  growth_A2 = NULL,
                                  growth_len_lower = NULL,
                                  growth_L0 = NULL,
                                  growth_cv_type = "len",
                                  growth_sd_type = "cv",
                                  growth_dist = "normal",
                                  growth_plus_group = "mixture",
                                  waa_model = "data",
                                  wt_len_pars = NULL,
                                  M_spec = "est_ln_M",
                                  M_popblk_spec = 'constant',
                                  M_ageblk_spec = 'constant',
                                  M_regionblk_spec = 'constant',
                                  M_yearblk_spec = 'constant',
                                  M_sexblk_spec = 'constant',
                                  Fixed_natmort = NULL,
                                  NAA_re = "none",
                                  NAA_re_ages = NULL,
                                  NAA_re_years = NULL,
                                  NAA_sigma_spec = "est",
                                  NAA_re_region = "iid",
                                  NAA_re_region_spec = "est_all",
                                  NAA_re_pop = "iid",
                                  NAA_re_sex = "iid",
                                  NAA_sigma_popblk_spec = 'constant',
                                  NAA_sigma_regionblk_spec = 'constant',
                                  NAA_sigma_yearblk_spec = 'constant',
                                  NAA_sigma_ageblk_spec = 'constant',
                                  NAA_sigma_sexblk_spec = 'constant',
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Biologicals <- mget(names(formals()))[-1]

  # Growth Options ---------------------------------------------------------
  n_pop <- input_list$data$n_pop; n_regions <- input_list$data$n_regions; n_sexes <- input_list$data$n_sexes
  n_yrs <- length(input_list$data$years); n_ages <- length(input_list$data$ages); n_seas <- input_list$data$n_seas
  if(!growth_model %in% c("none", "vb_schnute", "richards")) stop("growth_model must be one of: none, vb_schnute, richards")
  growth_model_val <- c(none = 0, vb_schnute = 1, richards = 2)[[growth_model]]
  gpar_names <- c("L1", "L2", "K", "CV1", "CV2", "rho")
  n_gpars <- if(growth_model_val == 2) 6 else 5
  if(growth_model_val != 0) {
    if(fit_lengths != 1) stop("growth_model = '", growth_model, "' builds the size-age transition inside the model, so fit_lengths must be 1")
    if(is.null(growth_len_lower)) stop("growth_len_lower (lower edges of the length bins) is required when growth is estimated")
    if(length(growth_len_lower) != length(input_list$data$lens)) stop("growth_len_lower must have one entry per length bin in Setup_Mod_Dim")
    if(is.null(growth_A1) || is.null(growth_A2)) stop("growth_A1 and growth_A2 (reference ages for L1 and L2) are required when growth is estimated")
    # L2 is either the length at a second reference age or the asymptote itself,
    # in which case the CV interpolation and the plus group read the accumulator age
    growth_L2_asymptote <- as.numeric(identical(growth_A2, "Linf"))
    if(growth_L2_asymptote == 1) growth_A2 <- max(input_list$data$ages)
    if(!is.numeric(growth_A1) || length(growth_A1) != 1) stop("growth_A1 must be a single reference age")
    if(!is.numeric(growth_A2) || length(growth_A2) != 1) stop("growth_A2 must be a single reference age, or \"Linf\" when L2 is the asymptote")
    if(is.null(growth_L0)) growth_L0 <- growth_len_lower[1]
    if(!growth_cv_type %in% c("len", "age")) stop("growth_cv_type must be len or age")
    if(!growth_sd_type %in% c("cv", "sd")) stop("growth_sd_type must be cv or sd")
    if(!growth_dist %in% c("normal", "lognormal")) stop("growth_dist must be normal or lognormal")
    if(!growth_plus_group %in% c("mixture", "curve")) stop("growth_plus_group must be mixture or curve")
    if(!growth_spec %in% c("est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "fix")) stop("growth_spec must be one of: est_all, est_shared_r, est_shared_s, est_shared_r_s, fix")
    if(is.null(growth_fix)) growth_fix <- rep(FALSE, n_gpars)
    if(length(growth_fix) != n_gpars) stop("growth_fix must be a logical vector of length ", n_gpars, " (", paste(gpar_names[1:n_gpars], collapse = ", "), ")")
    # starting values come through starting_values, as every other parameter's
    # do. the default puts the reference lengths at the ends of the length bins
    # with a middling rate. kept on the natural scale here for the bound checks
    gp_arr <- array(NA_real_, dim = c(n_pop, n_regions, n_sexes, n_gpars))
    gp_default <- c(L1 = min(input_list$data$lens), L2 = max(input_list$data$lens), K = 0.15, CV1 = 0.1, CV2 = 0.1, rho = 1)

    if("ln_growth_pars" %in% names(starting_values)) {
      sv_growth <- starting_values$ln_growth_pars
      if(is.null(dim(sv_growth)) || !all(dim(sv_growth) == dim(gp_arr))) stop("starting_values$ln_growth_pars must be an array [n_pop, n_regions, n_sexes, ", n_gpars, "]")
      gp_arr[] <- exp(sv_growth)
      if(any(!is.finite(gp_arr))) stop("starting_values$ln_growth_pars must be finite; the growth parameters are estimated on the log scale")
    } else {
      for(k in 1:n_gpars) gp_arr[, , , k] <- gp_default[[gpar_names[k]]]
      collect_message("No starting_values$ln_growth_pars supplied; starting from the length bins with K = 0.15")
    }
    # the size-age transition is built inside the model, so a placeholder stands in for the data checks
    SizeAgeTrans <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, length(input_list$data$lens), n_ages, n_sexes))
    collect_message("Growth is estimated (", if(growth_model_val == 1) "von Bertalanffy, Schnute form" else "Richards", "); SizeAgeTrans is built inside the model")

    # Time variation of the growth parameters ---------------------------------
    tv_vals <- rep(0, n_gpars); names(tv_vals) <- gpar_names[1:n_gpars]
    if(!is.null(growth_tv_model)) {
      tv_codes <- c(none = 0, iid = 1, rw = 2)
      if(!all(growth_tv_model %in% names(tv_codes))) stop("growth_tv_model entries must be one of: none, iid, rw")
      if(!is.null(names(growth_tv_model)) && all(names(growth_tv_model) != "")) {
        bad <- setdiff(names(growth_tv_model), gpar_names[1:n_gpars])
        if(length(bad) > 0) stop("growth_tv_model names not growth parameters: ", paste(bad, collapse = ", "), ". Use ", paste(gpar_names[1:n_gpars], collapse = ", "))
        for(name in names(growth_tv_model)) tv_vals[name] <- tv_codes[[growth_tv_model[[name]]]]
      } else {
        if(length(growth_tv_model) != n_gpars) stop("an unnamed growth_tv_model must have one entry per growth parameter (", n_gpars, "), or be named by parameter")
        for(k in 1:n_gpars) tv_vals[k] <- tv_codes[[growth_tv_model[k]]]
      }
    }
    if(!growth_tv_link %in% c("log", "logit")) stop("growth_tv_link must be log or logit")
    growth_tv_link_val <- c(log = 0, logit = 1)[[growth_tv_link]]
    if(growth_tv_link_val == 1) {
      if(is.null(growth_par_bounds)) stop("growth_par_bounds ([n_gpars x 2], natural scale) is required under the logit link")
      growth_par_bounds <- matrix(growth_par_bounds, ncol = 2)
      if(nrow(growth_par_bounds) != n_gpars) stop("growth_par_bounds must have one row per growth parameter (", n_gpars, ")")
      for(k in which(tv_vals > 0)) if(any(gp_arr[,,,k] <= growth_par_bounds[k, 1] | gp_arr[,,,k] >= growth_par_bounds[k, 2])) stop("growth_pars for ", gpar_names[k], " must lie strictly inside growth_par_bounds under the logit link")
    } else growth_par_bounds <- matrix(0, n_gpars, 2)
    if(!growth_tv_type %in% c("curve", "cohort")) stop("growth_tv_type must be curve or cohort")
    growth_tv_type_val <- c(curve = 0, cohort = 1)[[growth_tv_type]]
    if(!growth_tv_spec %in% c("est_all", "est_shared_r", "est_shared_s", "est_shared_r_s")) stop("growth_tv_spec must be one of: est_all, est_shared_r, est_shared_s, est_shared_r_s")
    if(!growth_tv_sigma_spec %in% c("fix", "est")) stop("growth_tv_sigma_spec must be fix or est")
    # active years per parameter, calendar years into indices
    tv_active <- matrix(0, n_yrs, n_gpars)
    for(k in which(tv_vals > 0)) {
      yrs_k <- if(is.null(growth_tv_years)) input_list$data$years else if(is.list(growth_tv_years)) growth_tv_years[[gpar_names[k]]] else growth_tv_years
      if(is.null(yrs_k)) yrs_k <- input_list$data$years
      if(!all(yrs_k %in% input_list$data$years)) stop("growth_tv_years for ", gpar_names[k], " has years outside the model years")
      tv_active[match(yrs_k, input_list$data$years), k] <- 1
    }
    growth_cohort_styr <- if(any(tv_vals > 0)) min(which(rowSums(tv_active) > 0)) else 1
    if(any(tv_vals > 0)) collect_message("Growth parameters varying over time: ", paste(paste0(gpar_names[tv_vals > 0], " (", c("none", "iid", "rw")[tv_vals[tv_vals > 0] + 1], ")"), collapse = ", "),
                                         "; link ", growth_tv_link, "; size at age read from ", if(growth_tv_type_val == 1) paste0("cohort propagation from ", input_list$data$years[growth_cohort_styr]) else "each year's curve")
    if(!waa_model %in% c("data", "wt_len")) stop("waa_model must be data or wt_len")
    if(waa_model == "wt_len") {
      if(is.null(wt_len_pars)) stop("wt_len_pars (a, b in W = a L^b) are required when waa_model = 'wt_len'")
      wl_arr <- array(NA_real_, dim = c(n_pop, n_regions, n_sexes, 2))
      if(is.null(dim(wt_len_pars))) {
        if(length(wt_len_pars) != 2) stop("wt_len_pars must be two values (a, b) or an array [n_pop, n_regions, n_sexes, 2]")
        for(k in 1:2) wl_arr[,,,k] <- wt_len_pars[k]
      } else {
        if(!all(dim(wt_len_pars) == c(n_pop, n_regions, n_sexes, 2))) stop("wt_len_pars array must be [n_pop, n_regions, n_sexes, 2]")
        wl_arr[] <- wt_len_pars
      }
      wt_len_pars <- wl_arr
      # weight at age is derived, so the inputs become placeholders too
      if(is.null(WAA)) WAA <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
      collect_message("Weight at age (spawning, fishery, survey) is derived from growth and the weight-length relationship")
    }

    # Semi-parametric growth: a deviation surface over years and ages ---------
    semipar_codes <- c(none = 0, iid = 1, rw = 2, `3dmarg` = 3, `3dcond` = 4, `2dar1` = 5)
    if(length(growth_semipar) != 1 || !growth_semipar %in% names(semipar_codes)) stop("growth_semipar must be one of: ", paste(names(semipar_codes), collapse = ", "))
    semipar_val <- semipar_codes[[growth_semipar]]
    if(!growth_semipar_spec %in% c("fix", "est")) stop("growth_semipar_spec must be fix or est")
    # the unconstrained scale a correlation is read on, 2/(1+exp(-2x))-1 inverted
    rho_untrans <- function(x) 0.5 * log((1 + x) / (1 - x))
    semipar_age_idx <- seq_len(n_ages)
    semipar_yr_idx <- seq_len(n_yrs)
    if(semipar_val > 0) {
      ages_use <- if(is.null(growth_semipar_ages)) input_list$data$ages else growth_semipar_ages
      yrs_use <- if(is.null(growth_semipar_years)) input_list$data$years else growth_semipar_years
      if(!all(ages_use %in% input_list$data$ages)) stop("growth_semipar_ages has ages outside the model ages")
      if(!all(yrs_use %in% input_list$data$years)) stop("growth_semipar_years has years outside the model years")
      semipar_age_idx <- match(ages_use, input_list$data$ages)
      semipar_yr_idx <- match(yrs_use, input_list$data$years)
      collect_message("Semi-parametric growth: ", growth_semipar, " deviations on mean length at age over ",
                      length(yrs_use), " years and ", length(ages_use), " ages, process error ", growth_semipar_spec)
    }
  } else {
    if(waa_model == "wt_len") stop("waa_model = 'wt_len' needs growth_model = 'vb_schnute' or 'richards'")
  }

  # Length Bin Map Options --------------------------------------------------
  if(!is.null(LenBinMap)) {
    LenBinMap <- as.matrix(LenBinMap)
    check_bin_map(LenBinMap, length(input_list$data$lens), "LenBinMap")
    collect_message("Length compositions are recorded on ", ncol(LenBinMap), " bins, mapped from the model's ", nrow(LenBinMap), " bins inside the likelihood")
  }

  # Input Validation --------------------------------------------------------

  # Weight at age checking

  check_data_dimensions(WAA, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_seas = input_list$data$n_seas, what = 'WAA')
  if(!is.null(WAA_fish)) check_data_dimensions(WAA_fish, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'WAA_fish')
  if(!is.null(WAA_srv)) check_data_dimensions(WAA_srv, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'WAA_srv')

  # Maturity at age checking
  check_data_dimensions(MatAA, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas,  n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'MatAA')

  # Age-0 (rec_lag = 0) recruitment requires the recruit age class (the first
  # age) to be immature everywhere, since age-0 fish can't spawn the year they're
  # born. This is relied on (rather than special-cased) when excluding age-0
  # from spawning biomass per recruit. Requires Setup_Mod_Rec() to have run
  # first so $data$rec_lag is already set.
  if(!is.null(input_list$data$rec_lag) && input_list$data$rec_lag == 0 && any(MatAA[,,,,1,] != 0)) {
    stop("rec_lag = 0 (age-0 recruitment) requires MatAA to be zero at the recruit age (the first age class) for all populations, regions, years, seasons, and sexes, since age-0 fish cannot be mature.")
  }

  # Length checking
  if(!fit_lengths %in% c(0,1)) stop("Values for fit_lengths are not valid. They are == 0 (not used), or == 1 (used)")
  collect_message("Length Composition data are: ", ifelse(fit_lengths == 0, "Not Used", "Used"))

  # Size Age Transition checking
  if(fit_lengths == 1) check_data_dimensions(SizeAgeTrans, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'SizeAgeTrans')
  if(fit_lengths == 1 & is.na(sum(SizeAgeTrans))) stop("Length composition are fit to, but the size-age transition matrix is NA")

  # Per-fleet fixed keys: only meaningful without a growth module, which already
  # derives one key per fleet and would leave two sources for the same quantity
  if(!is.null(SizeAgeTrans_fish) || !is.null(SizeAgeTrans_srv)) {
    if(growth_model_val != 0) stop("SizeAgeTrans_fish/SizeAgeTrans_srv are for growth_model = 'none'; a growth model already derives one key per fleet")
    if(!is.null(SizeAgeTrans_fish)) check_data_dimensions(SizeAgeTrans_fish, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'SizeAgeTrans_fish')
    if(!is.null(SizeAgeTrans_srv)) check_data_dimensions(SizeAgeTrans_srv, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SizeAgeTrans_srv')
    if(!is.null(SizeAgeTrans_fish)) collect_message("Fishery keys read per fleet from SizeAgeTrans_fish rather than the shared SizeAgeTrans")
    if(!is.null(SizeAgeTrans_srv)) collect_message("Survey keys read per fleet from SizeAgeTrans_srv rather than the shared SizeAgeTrans")
  }

  # Joint length and age array checking
  if(!do_caal %in% c(0,1)) stop("Values for do_caal are not valid. They are == 0 (not used), or == 1 (used)")
  if(do_caal == 1 && fit_lengths == 0) stop("do_caal == 1 requires fit_lengths == 1, since the joint arrays are built from the size-age transition matrix")
  if(do_caal == 1) collect_message("Joint arrays at length and age (Fish_caal, Fish_caal_discard, Srv_caal) are: Reported")

  # Natural Mortality checking
  if(!is.null(M_spec)) {
    if(M_spec == 'fix') {
      if(is.null(Fixed_natmort)) stop("Please provide a fixed natural mortality array dimensioned by n_pop, n_regions, n_years, n_ages, and n_sexes!")
      check_data_dimensions(Fixed_natmort, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'Fixed_natmort')
    }
  }

  # Check M blocks
  if(!is.null(M_ageblk_spec)) if(!typeof(M_ageblk_spec) %in% c("list", "character", NULL)) stop("M fixed effects age blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 ages and wanted 2 age blocks, this would be list(c(1:5), c(6:10)) such that ages 1 - 5 are a block, and ages 6 - 10 are a block.")
  if(!is.null(M_yearblk_spec)) if(!typeof(M_yearblk_spec) %in% c("list", "character", NULL)) stop("M fixed effects year blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 years and wanted 2 year blocks, this would be list(c(1:5), c(6:10)) such that years 1 - 5 are a block, and years 6 - 10 are a block.")
  if(!is.null(M_sexblk_spec)) if(!typeof(M_sexblk_spec) %in% c("list", "character", NULL)) stop("M fixed effects sex blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 2 sexes and wanted sex-specific M, this would be list(1, 2).")
  if(!is.null(M_regionblk_spec)) if(!typeof(M_regionblk_spec) %in% c("list", "character", NULL)) stop("M fixed effects region blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 2 regions and wanted region-specific M, this would be list(1, 2).")
  if(!is.null(M_popblk_spec)) if(!typeof(M_popblk_spec) %in% c("list", "character", NULL)) stop("M fixed effects population blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 2 populations and wanted population-specific M, this would be list(1, 2).")

  # Natural Mortality prior checking
  if(!Use_M_prior %in% c(0,1)) stop("Values for Use_M_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  collect_message("Natural Mortality priors are: ", ifelse(Use_M_prior == 0, "Not Used", "Used"))

  if(Use_M_prior == 1) {
    required_cols <- c("popblk", "regionblk", "yearblk", "ageblk", "sexblk", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(M_prior))
    if(length(missing_cols) > 0) {
      stop("M_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }

  # Checking ageing error dimensions. AgeingError and LenBinMap are the same
  # model-bin to observed-bin map on different axes, so both go through
  # check_bin_map and a mistake in either reads the same way.
  if(!is.null(AgeingError)) {
    if(length(dim(AgeingError)) == 2) { # user supplied ageing error is not time-varying
      check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), what = 'AgeingError')
      check_bin_map(AgeingError, length(input_list$data$ages), "AgeingError", strict = FALSE, tol = 0.05)
    }
    if(length(dim(AgeingError)) == 3) { # user supplied ageing error is time-varying
      check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), n_years = length(input_list$data$years), what = 'AgeingError_t')
      for(i in 1:dim(AgeingError)[1]) check_bin_map(AgeingError[i,,], length(input_list$data$ages), paste0("AgeingError year ", i), strict = FALSE, tol = 0.05)
    } # end i loop
  }

  # Weight at Age Options ---------------------------------------------------

  # setup fishery and survey specific weight at age (if not specified - just uses the WAA (spawning) already supplied)
  if(is.null(WAA_fish)) { # if no fishery WAA provided, use spawning WAA supplied
    WAA_fish <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), n_seas = input_list$data$n_seas, length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(f in 1:input_list$data$n_fish_fleets) WAA_fish[,,,,,,f] <- WAA
    collect_message("WAA_fish was specified at NULL. Using the spawning WAA for WAA_fish")
  }

  # if no survey WAA provided, use spawning WAA supplied
  if(is.null(WAA_srv)) {
    WAA_srv <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), n_seas = input_list$data$n_seas, length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(f in 1:input_list$data$n_srv_fleets) WAA_srv[,,,,,,f] <- WAA
    collect_message("WAA_srv was specified at NULL. Using the spawning WAA for WAA_srv")
  }

  # Ageing Error Options ----------------------------------------------------

  # setup ageing error if not provided
  if(is.null(AgeingError)) {
    AgeingError <- diag(1, length(input_list$data$ages)) # if no inputs for ageing error, then create identity matrix
    AgeingError_t <- array(0, dim = c(length(input_list$data$years), dim(AgeingError)))
    for(i in 1:length(input_list$data$years)) AgeingError_t[i,,] <- AgeingError
    warning("No ageing error matrix was provided. A default identity matrix was used, which assumes that the number and structure of modeled age bins exactly match the observed age bins. If the observed age composition data includes fewer age bins than the model (e.g., observed ages 2-10 while modeled ages are 1-10), this default assumption will cause a dimensional mismatch and potentially misalign the modeled and observed compositions. To avoid this, please provide an ageing error matrix of dimension n_model_ages x n_obs_ages that correctly maps modeled ages to observed age bins. For example, if observed ages are 2-10, supply a matrix that drops the first model age by using a shifted identity matrix: diag(1, 10)[, 2:10]. This will ensure the age bins are correctly aligned for likelihood calculations.")
  } else if(length(dim(AgeingError)) == 2) {   # setup ageing error if user-supplied is not year specific
    AgeingError_t <- array(0, dim = c(length(input_list$data$years), dim(AgeingError)))
    for(i in 1:length(input_list$data$years)) AgeingError_t[i,,] <- AgeingError
    collect_message("Ageing Error is specified to be time-invariant")
  } else if(length(dim(AgeingError)) == 3) {   # ageing error if it is year specific (just reassigning)
    AgeingError_t <- AgeingError
    collect_message("Ageing Error is specified to be time-varying")
  }

  # expand fleet-specific ageing error
  AgeingError_fish_t <- expand_fleet_ageing_error(AgeingError_fish, AgeingError_t, input_list$data$n_fish_fleets, "AgeingError_fish")
  AgeingError_srv_t <- expand_fleet_ageing_error(AgeingError_srv, AgeingError_t, input_list$data$n_srv_fleets, "AgeingError_srv")

  # Natural Mortality Options -----------------------------------------------
  # Input indicator for estimating or not estimating M
  if(is.null(M_spec) || M_spec == "est_ln_M") input_list$data$use_fixed_natmort <- 0
  else if(M_spec == "fix") input_list$data$use_fixed_natmort <- 1

  # Populate Data List ------------------------------------------------------

  input_list$data$WAA <- WAA
  input_list$data$WAA_fish <- WAA_fish
  input_list$data$WAA_srv <- WAA_srv
  input_list$data$MatAA <- MatAA
  input_list$data$AgeingError <- AgeingError_t
  input_list$data$AgeingError_fish <- AgeingError_fish_t
  input_list$data$AgeingError_srv <- AgeingError_srv_t
  input_list$data$fit_lengths <- fit_lengths
  input_list$data$SizeAgeTrans <- SizeAgeTrans
  input_list$data$SizeAgeTrans_fish <- SizeAgeTrans_fish
  input_list$data$SizeAgeTrans_srv <- SizeAgeTrans_srv
  input_list$data$do_caal <- do_caal

  # Growth module storage
  input_list$data$growth_model <- growth_model_val
  if(growth_model_val != 0) {
    input_list$data$growth_A1 <- growth_A1
    input_list$data$growth_A2 <- growth_A2
    input_list$data$growth_L2_asymptote <- growth_L2_asymptote
    input_list$data$growth_L0 <- growth_L0
    input_list$data$growth_len_lower <- growth_len_lower
    input_list$data$growth_cv_type <- c(len = 0, age = 1)[[growth_cv_type]]
    input_list$data$growth_sd_type <- c(cv = 0, sd = 1)[[growth_sd_type]]
    input_list$data$growth_dist <- c(normal = 0, lognormal = 1)[[growth_dist]]
    input_list$data$growth_plus_group <- c(mixture = 1, curve = 0)[[growth_plus_group]]
    input_list$data$derive_waa <- c(data = 0, wt_len = 1)[[waa_model]]
    input_list$data$wt_len_pars <- if(waa_model == "wt_len") wt_len_pars else array(0, dim = c(n_pop, n_regions, n_sexes, 2))
    # time variation of the growth parameters
    input_list$data$growth_tv_model <- as.numeric(tv_vals)
    input_list$data$growth_tv_link <- growth_tv_link_val
    input_list$data$growth_par_bounds <- growth_par_bounds
    input_list$data$growth_tv_type <- growth_tv_type_val
    input_list$data$growth_cohort_styr <- growth_cohort_styr
    input_list$data$growth_rw_init_sigma <- growth_rw_init_sigma
    # semi-parametric surface
    input_list$data$growth_semipar <- semipar_val
    input_list$data$growth_semipar_bins <- semipar_age_idx
  }
  input_list$data$LenBinMap <- LenBinMap
  input_list$data$Use_M_prior <- Use_M_prior
  input_list$data$M_prior <- M_prior
  input_list$data$Fixed_natmort <- Fixed_natmort
  # addtocomp/comp_const_obs/addtofishidx/addtosrvidx/addtotag now belong to
  # Setup_Mod_Weighting; a value still passed here is stashed and picked up
  # there, so old scripts keep working
  legacy_weighting <- list(addtocomp = addtocomp, comp_const_obs = comp_const_obs, addtofishidx = addtofishidx, addtosrvidx = addtosrvidx, addtotag = addtotag)
  if(any(!vapply(legacy_weighting, is.null, logical(1))))
    collect_message("addtocomp/comp_const_obs/addtofishidx/addtosrvidx/addtotag passed to Setup_Mod_Biologicals are deprecated; pass them to Setup_Mod_Weighting instead.")
  input_list$.legacy_weighting <- legacy_weighting

  # Populate Parameter List -------------------------------------------------

  # If M is constant for ages
  if(is.character(M_ageblk_spec)) {
    if(!identical(M_ageblk_spec, "constant")) stop("M_ageblk_spec must be \"constant\" or a list of age blocks, but was: ", M_ageblk_spec)
    M_ageblk_spec_vals <- list(1:length(input_list$data$ages))
  } else M_ageblk_spec_vals <- M_ageblk_spec

  # If M is constant across years
  if(is.character(M_yearblk_spec)) {
    if(!identical(M_yearblk_spec, "constant")) stop("M_yearblk_spec must be \"constant\" or a list of year blocks, but was: ", M_yearblk_spec)
    M_yearblk_spec_vals <- list(1:length(input_list$data$years))
  } else M_yearblk_spec_vals <- M_yearblk_spec

  # If M is constant across sexes
  if(is.character(M_sexblk_spec)) {
    if(!identical(M_sexblk_spec, "constant")) stop("M_sexblk_spec must be \"constant\" or a list of sex blocks, but was: ", M_sexblk_spec)
    M_sexblk_spec_vals <- list(1:input_list$data$n_sexes)
  } else M_sexblk_spec_vals <- M_sexblk_spec

  # If M is constant across regions
  if(is.character(M_regionblk_spec)) {
    if(!identical(M_regionblk_spec, "constant")) stop("M_regionblk_spec must be \"constant\" or a list of region blocks, but was: ", M_regionblk_spec)
    M_regionblk_spec_vals <- list(1:input_list$data$n_regions)
  } else M_regionblk_spec_vals <- M_regionblk_spec

  # If M is constant across populations
  if(is.character(M_popblk_spec)) {
    if(!identical(M_popblk_spec, "constant")) stop("M_popblk_spec must be \"constant\" or a list of population blocks, but was: ", M_popblk_spec)
    M_popblk_spec_vals <- list(1:input_list$data$n_pop)
  } else M_popblk_spec_vals <- M_popblk_spec

  input_list$par$ln_M <- array(log(0.5), dim = c(length(M_popblk_spec_vals),
                                                      length(M_regionblk_spec_vals),
                                                      length(M_yearblk_spec_vals),
                                                      length(M_ageblk_spec_vals),
                                                      length(M_sexblk_spec_vals)))
  input_list$par$ln_M <- use_starting_value(input_list$par$ln_M, starting_values, "ln_M")

  # Growth parameters, the deviations of any that vary over time, and the
  # semi-parametric surface on mean length at age
  if(growth_model_val != 0) {

    input_list$par$ln_growth_pars <- log(gp_arr)
    input_list$par$ln_growth_pars <- use_starting_value(input_list$par$ln_growth_pars, starting_values, "ln_growth_pars")

    input_list$par$ln_growth_devs <- array(0, dim = c(n_pop, n_regions, n_yrs, n_gpars, n_sexes))
    input_list$par$ln_growth_semipar_devs <- array(0, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes))

    # One process error array carries both growth deviation streams, in the same
    # slots the selectivity forms use. Stream one holds a log sigma per growth
    # parameter for the time-varying deviations. Stream two holds the surface's
    # correlations by age, year and cohort in slots one to three and a log scale
    # in slot four for the correlated forms, or one log sigma per age for iid
    # and the random walk. The array is as wide as whichever stream needs more,
    # and every slot a form does not read is mapped off and never reaches the
    # objective.
    if("growth_pe_pars" %in% names(starting_values)) {

      input_list$par$growth_pe_pars <- starting_values$growth_pe_pars

    } else {

      pe <- array(0, dim = c(n_pop, n_regions, max(4, n_ages, n_gpars), n_sexes, 2))
      pe[, , , , 1] <- log(0.1)  # time-varying growth parameters
      pe[, , , , 2] <- log(0.05) # the surface's scale, and its per-age sigmas

      # the correlated forms read correlations in the first three slots
      if(semipar_val %in% 3:5) {
        pe[, , 1, , 2] <- rho_untrans(0.3)
        pe[, , 2, , 2] <- rho_untrans(0.3)
        if(semipar_val %in% 3:4) pe[, , 3, , 2] <- rho_untrans(0.3)
      }

      input_list$par$growth_pe_pars <- pe

    }
  }

  # Mapping Options ---------------------------------------------------------
  input_list <- do_natmort_mapping(input_list, M_spec, M_popblk_spec_vals, M_regionblk_spec_vals,
                             M_yearblk_spec_vals, M_ageblk_spec_vals, M_sexblk_spec_vals) # natural mortality mapping

  # State-space numbers at age
  NAA_sigma_blk_vals <- lapply(
    list(pop = NAA_sigma_popblk_spec, region = NAA_sigma_regionblk_spec, year = NAA_sigma_yearblk_spec,
         age = NAA_sigma_ageblk_spec, sex = NAA_sigma_sexblk_spec),
    function(x) x)
  n_by_dim <- list(pop = input_list$data$n_pop, region = input_list$data$n_regions,
                   year = length(input_list$data$years), age = length(input_list$data$ages),
                   sex = input_list$data$n_sexes)
  for(name in names(NAA_sigma_blk_vals)) {
    if(is.character(NAA_sigma_blk_vals[[name]])) {
      if(!identical(NAA_sigma_blk_vals[[name]], "constant"))
        stop("NAA_sigma_", name, "blk_spec must be \"constant\" or a list of blocks, but was: ",
             NAA_sigma_blk_vals[[name]])
      NAA_sigma_blk_vals[[name]] <- list(1:n_by_dim[[name]])
    }
  } # end name loop

  input_list <- do_NAAstate_mapping(input_list, NAA_re, NAA_re_ages, NAA_re_years, NAA_sigma_spec,
                                    NAA_re_region, NAA_re_region_spec, NAA_re_pop, NAA_re_sex,
                                    starting_values,
                                    NAA_sigma_blk_vals$pop, NAA_sigma_blk_vals$region,
                                    NAA_sigma_blk_vals$year, NAA_sigma_blk_vals$age,
                                    NAA_sigma_blk_vals$sex)
  if(growth_model_val != 0) input_list <- do_growth_mapping(input_list, growth_spec, growth_fix, tv_vals, tv_active, growth_tv_spec,
                                                            growth_tv_sigma_spec, semipar_val, growth_semipar_spec,
                                                            semipar_age_idx, semipar_yr_idx) # growth mapping

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}


#' Set up the state-space numbers at age
#'
#' Builds the \code{ln_NAA} parameter array and its map, the process error
#' standard deviations and their blocking index, and the data fields the
#' dynamics and the penalty read.
#'
#' The state covers ages 2 and older, including the plus group, over years 2 and
#' later. Year 1 at those ages belongs to \code{ln_InitDevs} and age 1 belongs to
#' \code{ln_RecDevs} in every year, so the three parameterizations partition the
#' numbers at age rather than overlapping. The plus group is inside the state and
#' not optional: it is the only cell whose influence never decays, so leaving it
#' deterministic gives up the conditional independence the state is worth having
#' for.
#'
#' Ages and years must each be a contiguous run. The penalty scores the active
#' cells as one rectangular slice, and a gap would either leave scored cells the
#' dynamics never wrote or force a per-cell branch onto the tape.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#' @param NAA_re Character. \code{"none"} (default) leaves the numbers at age
#'   deterministic. \code{"iid"} gives every active cell an independent Gaussian
#'   innovation.
#' @param NAA_re_ages Ages the state is active over, as ages rather than indices.
#'   \code{NULL} (default) uses every age from the second onward.
#' @param NAA_re_years Calendar years the state is active over. \code{NULL}
#'   (default) uses every year from the second onward.
#' @param NAA_sigma_spec Character, \code{"est"} or \code{"fix"}, whether the
#'   process error standard deviations are estimated.
#' @param NAA_sigma_popblk_spec_vals,NAA_sigma_regionblk_spec_vals,NAA_sigma_yearblk_spec_vals,NAA_sigma_ageblk_spec_vals,NAA_sigma_sexblk_spec_vals
#'   Lists of integer vectors assigning indices to blocks, exactly as the
#'   \code{M_*blk_spec_vals} arguments do. Blocking shares the standard
#'   deviation; it never removes a cell from the state.
#'
#' @return \code{input_list} with \code{$par$ln_NAA}, \code{$par$ln_sigmaNAA},
#'   their maps, and the data fields \code{NAA_re}, \code{n_est_naa_re},
#'   \code{naa_re_ages}, \code{naa_re_yrs} and \code{naa_sigma_blocks}.
#'
#' @keywords internal
do_NAAstate_mapping <- function(input_list,
                                NAA_re,
                                NAA_re_ages,
                                NAA_re_years,
                                NAA_sigma_spec,
                                NAA_re_region = "iid",
                                NAA_re_region_spec = "est_all",
                                NAA_re_pop = "iid",
                                NAA_re_sex = "iid",
                                starting_values = list(),
                                NAA_sigma_popblk_spec_vals,
                                NAA_sigma_regionblk_spec_vals,
                                NAA_sigma_yearblk_spec_vals,
                                NAA_sigma_ageblk_spec_vals,
                                NAA_sigma_sexblk_spec_vals) {

  n_pop <- input_list$data$n_pop
  n_regions <- input_list$data$n_regions
  n_sexes <- input_list$data$n_sexes
  ages <- input_list$data$ages
  years <- input_list$data$years
  n_ages <- length(ages)
  n_yrs <- length(years)

  # 1dar1 without a suffix means across ages, matching do_age_corr_setup, where age is the only
  # margin those observations have. The state has two, so both are also spellable explicitly.
  naa_codes <- c(none = 0, iid = 1, `1dar1_a` = 2, `1dar1_y` = 3, `2dar1` = 4, `3dcond` = 5, `3dmarg` = 6)
  if(length(NAA_re) != 1 || !NAA_re %in% names(naa_codes))
    stop("NAA_re is '", NAA_re, "'. Valid options: ", paste(names(naa_codes), collapse = ", "))
  naa_val <- naa_codes[[NAA_re]]

  # the array always exists so the objective can index it unconditionally; whether any of it
  # is estimated is carried by n_est_naa_re, never inferred from the array's dimensions.
  # Unseeded, the cold start is an equilibrium decay from R0 at the mean M, split over regions and
  # sexes, with the plus group accumulated. A constant would put every age at the same abundance,
  # which is nowhere near any plausible solution; the decay at least has the right shape. A seed
  # from a deterministic pass, passed as starting_values$ln_NAA, is always better and overrides it.
  M_bar <- mean(exp(as.vector(input_list$par$ln_M)))
  ln_R0 <- if(is.null(input_list$par$ln_global_R0)) 5 else as.vector(input_list$par$ln_global_R0)[1]
  eq_naa <- ln_R0 - log(n_regions) - log(n_sexes) - M_bar * (seq_len(n_ages) - 1)
  eq_naa[n_ages] <- eq_naa[n_ages] - log(1 - exp(-M_bar)) # plus group at equilibrium
  if(!all(is.finite(eq_naa))) eq_naa <- rep(5, n_ages)
  input_list$par$ln_NAA <- use_starting_value(
    array(rep(eq_naa, each = n_pop * n_regions * n_yrs),
          dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes)),
    starting_values, "ln_NAA")

  # map all of this stuff off if no process error
  if(NAA_re == "none") {
    input_list$map$ln_NAA <- factor(rep(NA, length(input_list$par$ln_NAA)))
    input_list$par$ln_sigmaNAA <- array(log(0.3), dim = c(1, 1, 1, 1, 1))
    input_list$map$ln_sigmaNAA <- factor(NA)
    input_list$data$NAA_re <- 0
    input_list$data$n_est_naa_re <- 0
    input_list$data$naa_re_ages <- integer(0)
    input_list$data$naa_re_yrs <- integer(0)
    input_list$data$naa_sigma_blocks <- array(1, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes))
    input_list$par$NAA_pe_pars <- array(0, dim = c(n_pop, n_regions, 3, n_sexes))
    input_list$map$NAA_pe_pars <- factor(rep(NA, length(input_list$par$NAA_pe_pars)))
    input_list$data$NAA_re_region <- 0
    input_list$par$NAA_region_corr_pars <- array(0, dim = c(n_pop, max(1, n_regions * (n_regions - 1) / 2), n_sexes))
    input_list$map$NAA_region_corr_pars <- factor(rep(NA, length(input_list$par$NAA_region_corr_pars)))
    input_list$data$NAA_re_pop <- 0
    input_list$data$NAA_re_sex <- 0
    input_list$par$NAA_pop_corr_pars <- rep(0, max(1, n_pop * (n_pop - 1) / 2))
    input_list$map$NAA_pop_corr_pars <- factor(rep(NA, length(input_list$par$NAA_pop_corr_pars)))
    input_list$par$NAA_sex_corr_pars <- rep(0, max(1, n_sexes * (n_sexes - 1) / 2))
    input_list$map$NAA_sex_corr_pars <- factor(rep(NA, length(input_list$par$NAA_sex_corr_pars)))
    return(input_list)
  }

  # A time-blocked M that is freely estimated is the same unpenalized log-survival surface the
  # state already carries, so the two trade off one for one and neither is pinned down.
  m_map <- input_list$map$ln_M
  m_estimated <- is.null(m_map) || any(!is.na(m_map))
  if(m_estimated && length(dim(input_list$par$ln_M)) >= 3 && dim(input_list$par$ln_M)[3] > 1)
    stop("A freely estimated time-blocked natural mortality and the numbers-at-age state are not ",
         "separately identified: both are an unpenalized log-survival surface over years and ages. ",
         "Fix M with M_spec = 'fix', or drop the year blocking, or turn the state off.")

  # Ages and years the state runs over
  age_use <- if(is.null(NAA_re_ages)) ages[-1] else NAA_re_ages
  yr_use <- if(is.null(NAA_re_years)) years[-1] else NAA_re_years
  if(!all(age_use %in% ages)) stop("NAA_re_ages has ages outside the model ages.")
  if(!all(yr_use %in% years)) stop("NAA_re_years has years outside the model years.")

  age_idx <- sort(match(age_use, ages))
  yr_idx <- sort(match(yr_use, years))

  if(any(age_idx < 2))
    stop("NAA_re_ages includes the first age. Age one is recruitment and belongs to ln_RecDevs; ",
         "the state covers ages two and older.")
  if(any(yr_idx < 2))
    stop("NAA_re_years includes the first year. Year one at ages two and older belongs to ",
         "ln_InitDevs; the state covers years two and later.")
  if(length(age_idx) > 1 && !all(diff(age_idx) == 1))
    stop("NAA_re_ages must be a contiguous run of ages. The state is scored as one rectangular ",
         "slice, so a gap would leave scored cells the dynamics never wrote.")
  if(length(yr_idx) > 1 && !all(diff(yr_idx) == 1))
    stop("NAA_re_years must be a contiguous run of years. The state is scored as one rectangular ",
         "slice, so a gap would leave scored cells the dynamics never wrote.")

  # Map: estimate the active rectangle, hold everything else
  map_naa <- array(NA_integer_, dim = dim(input_list$par$ln_NAA))
  n_active <- n_pop * n_regions * length(yr_idx) * length(age_idx) * n_sexes
  map_naa[,,yr_idx,age_idx,] <- seq_len(n_active)
  input_list$map$ln_NAA <- factor(map_naa)
  input_list$data$map_ln_NAA <- array(as.numeric(input_list$map$ln_NAA), dim = dim(map_naa))

  # Process error standard deviations, blocked the way natural mortality is
  sigma_blocks <- array(0, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes))
  counter <- 1
  for(popblk in 1:length(NAA_sigma_popblk_spec_vals)) {
    map_p <- NAA_sigma_popblk_spec_vals[[popblk]]
    for(regionblk in 1:length(NAA_sigma_regionblk_spec_vals)) {
      map_r <- NAA_sigma_regionblk_spec_vals[[regionblk]]
      for(yearblk in 1:length(NAA_sigma_yearblk_spec_vals)) {
        map_y <- NAA_sigma_yearblk_spec_vals[[yearblk]]
        for(ageblk in 1:length(NAA_sigma_ageblk_spec_vals)) {
          map_a <- NAA_sigma_ageblk_spec_vals[[ageblk]]
          for(sexblk in 1:length(NAA_sigma_sexblk_spec_vals)) {
            map_s <- NAA_sigma_sexblk_spec_vals[[sexblk]]
            sigma_blocks[map_p, map_r, map_y, map_a, map_s] <- counter
            counter <- counter + 1
          } # end sexblk loop
        } # end ageblk loop
      } # end yearblk loop
    } # end regionblk loop
  } # end popblk loop

  if(any(sigma_blocks == 0))
    stop("The NAA_sigma block specifications do not cover every cell. Each of the population, ",
         "region, year, age and sex block lists must partition its dimension.")

  input_list$par$ln_sigmaNAA <- array(log(0.3), dim = c(length(NAA_sigma_popblk_spec_vals),
                                                        length(NAA_sigma_regionblk_spec_vals),
                                                        length(NAA_sigma_yearblk_spec_vals),
                                                        length(NAA_sigma_ageblk_spec_vals),
                                                        length(NAA_sigma_sexblk_spec_vals)))

  if(!NAA_sigma_spec %in% c("est", "fix"))
    stop("NAA_sigma_spec is '", NAA_sigma_spec, "'. Valid options: est, fix")
  if(NAA_sigma_spec == "est") input_list$map$ln_sigmaNAA <- factor(1:length(input_list$par$ln_sigmaNAA))
  if(NAA_sigma_spec == "fix") input_list$map$ln_sigmaNAA <- factor(rep(NA, length(input_list$par$ln_sigmaNAA)))

  # Only the independent form admits a standard deviation that varies cell by cell. Every other
  # structure is separable or Markov in a margin, so one standard deviation has to serve the whole
  # age-year grid within a population, region and sex; a per-cell variance would break separability.
  if(naa_val > 1 && (length(NAA_sigma_yearblk_spec_vals) > 1 || length(NAA_sigma_ageblk_spec_vals) > 1)) {
    stop("NAA_re = \"", NAA_re, "\" is a correlated structure over the age and year grid, but the ",
         "process error standard deviation is blocked over ", length(NAA_sigma_yearblk_spec_vals),
         " year and ", length(NAA_sigma_ageblk_spec_vals), " age blocks. A separable or Markov ",
         "correlation carries one standard deviation per population, region and sex. Hold ",
         "NAA_sigma_yearblk_spec and NAA_sigma_ageblk_spec at \"constant\", or use NAA_re = \"iid\", ",
         "which is the only form a cell-varying standard deviation is defined for.")
  }

  # Correlation parameters, read as age, year and cohort
  input_list$par$NAA_pe_pars <- array(0, dim = c(n_pop, n_regions, 3, n_sexes))
  input_list$par$NAA_pe_pars <- use_starting_value(input_list$par$NAA_pe_pars, starting_values, "NAA_pe_pars")
  # which of the three slots (age, year, cohort) a form reads; the rest stay at zero and mapped off
  pe_slots <- list(`1` = integer(0), `2` = 1L, `3` = 2L, `4` = c(1L, 2L),
                   `5` = 1:3, `6` = 1:3)[[as.character(naa_val)]]
  map_pe <- array(NA_integer_, dim = dim(input_list$par$NAA_pe_pars))
  if(length(pe_slots)) map_pe[,,pe_slots,] <- seq_len(n_pop * n_regions * length(pe_slots) * n_sexes)
  input_list$map$NAA_pe_pars <- factor(map_pe)

  # Correlation across regions
  region_codes <- c(iid = 0, us = 1)
  if(length(NAA_re_region) != 1 || !NAA_re_region %in% names(region_codes))
    stop("NAA_re_region is '", NAA_re_region, "'. Valid options: ", paste(names(region_codes), collapse = ", "))
  region_val <- region_codes[[NAA_re_region]]

  if(region_val > 0 && n_regions == 1)
    stop("NAA_re_region = \"", NAA_re_region, "\" needs more than one region, but the model has one. ",
         "There is nothing for a cross-region correlation to describe; leave it at \"iid\".")

  n_pairs <- max(1, n_regions * (n_regions - 1) / 2)
  input_list$par$NAA_region_corr_pars <- array(0, dim = c(n_pop, n_pairs, n_sexes))
  input_list$par$NAA_region_corr_pars <- use_starting_value(input_list$par$NAA_region_corr_pars,
                                                            starting_values, "NAA_region_corr_pars")

  # The correlation sits over population and sex, so "est_shared_p_s" gives one correlation matrix for the whole
  # model and "est_all" a free one per population and sex.
  valid_spec <- c("est_all", "est_shared_p", "est_shared_s", "est_shared_p_s", "fix")
  if(length(NAA_re_region_spec) != 1 || !NAA_re_region_spec %in% valid_spec)
    stop("NAA_re_region_spec is '", NAA_re_region_spec, "'. Valid options: ",
         paste(valid_spec, collapse = ", "))

  map_rc <- array(NA_integer_, dim = c(n_pop, n_pairs, n_sexes))
  if(region_val > 0 && NAA_re_region_spec != "fix") {
    share_p <- NAA_re_region_spec %in% c("est_shared_p", "est_shared_p_s")
    share_s <- NAA_re_region_spec %in% c("est_shared_s", "est_shared_p_s")
    for(p1 in 1:n_pop) {
      for(k in 1:n_pairs) {
        for(s1 in 1:n_sexes) {
          pi_ <- if(share_p) 1 else p1
          si_ <- if(share_s) 1 else s1
          map_rc[p1,k,s1] <- (pi_ - 1) * n_pairs * n_sexes + (k - 1) * n_sexes + si_
        } # end s1 loop
      } # end k loop
    } # end p1 loop
    map_rc[] <- as.integer(factor(map_rc)) # renumber to a dense sequence
  }
  input_list$map$NAA_region_corr_pars <- factor(map_rc)
  input_list$data$NAA_re_region <- region_val

  # Population and sex correlations
  for(mg in c("pop", "sex")) {
    val_chr <- if(mg == "pop") NAA_re_pop else NAA_re_sex
    n_lvl <- if(mg == "pop") n_pop else n_sexes
    if(length(val_chr) != 1 || !val_chr %in% names(region_codes))
      stop("NAA_re_", mg, " is '", val_chr, "'. Valid options: ", paste(names(region_codes), collapse = ", "))
    val <- region_codes[[val_chr]]
    if(val > 0 && n_lvl == 1)
      stop("NAA_re_", mg, " = \"", val_chr, "\" needs more than one ", mg, ", but the model has one. ",
           "There is nothing for the correlation to describe; leave it at \"iid\".")
    name <- paste0("NAA_", mg, "_corr_pars")
    input_list$par[[name]] <- rep(0, max(1, n_lvl * (n_lvl - 1) / 2))
    input_list$par[[name]] <- use_starting_value(input_list$par[[name]], starting_values, name)
    input_list$map[[name]] <- factor(if(val > 0) seq_along(input_list$par[[name]])
                                   else rep(NA, length(input_list$par[[name]])))
    input_list$data[[paste0("NAA_re_", mg)]] <- val
  } # end mg loop

  # A correlated survival shock between sexes
  if(input_list$data$NAA_re_sex > 0)
    collect_message("NAA_re_sex is estimating a survival correlation between sexes. Spawning biomass ",
                    "reads sex one while indices sum over sexes, so this is weakly identified unless ",
                    "sex-structured compositions carry the contrast. Profile it before trusting it.")

  input_list$data$NAA_re <- naa_val
  input_list$data$n_est_naa_re <- n_active
  input_list$data$naa_re_ages <- age_idx
  input_list$data$naa_re_yrs <- yr_idx
  input_list$data$naa_sigma_blocks <- sigma_blocks

  collect_message("State-space numbers at age: ", NAA_re, " innovations over ",
                  length(yr_idx), " years and ", length(age_idx), " ages (",
                  n_active, " states), process error ", NAA_sigma_spec)

  input_list
}
