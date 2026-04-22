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
#'   Values should be proportions in \eqn{[0, 1]}.
#' @param AgeingError_input Ageing error (age-length transition) array with
#'   dimensions \code{[n_yrs × n_model_ages × n_obs_ages × n_sims]}, where each
#'   \code{[n_model_ages × n_obs_ages]} slice is a row-stochastic matrix mapping
#'   true modelled ages to observed age bins. If \code{NULL} (default), an
#'   identity matrix is constructed for each year and simulation, which assumes
#'   that modelled and observed age bins are identical in number and alignment.
#'   \strong{If observed age bins are a subset of modelled ages} (e.g., observed
#'   ages 2–10 vs. modelled ages 1–10), the default identity matrix will cause a
#'   dimensional mismatch. In that case, supply a shifted identity matrix such as
#'   \code{diag(1, n_model_ages)[, obs_age_index]} to correctly drop or collapse
#'   model ages into observed bins.
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
                                  SizeAgeTrans_input = NULL
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
  if(!is.null(SizeAgeTrans_input)) check_sim_dimensions(SizeAgeTrans_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_lens = sim_list$n_lens, n_seas = sim_list$n_seas, n_pop = sim_list$n_pop,
                                                        n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_sims = sim_list$n_sims, what = 'SizeAgeTrans_input')

  # output into list
  sim_list$natmort <- natmort_input
  sim_list$WAA <- WAA_input
  sim_list$WAA_fish <- WAA_fish_input
  sim_list$WAA_srv <- WAA_srv_input
  sim_list$MatAA <- MatAA_input
  if(!is.null(SizeAgeTrans_input)) sim_list$SizeAgeTrans <- SizeAgeTrans_input
  if(!is.null(AgeingError_input)) sim_list$AgeingError <- AgeingError_input
  else {
    # if null, create an identity matrix
    identity_AgeingError <- array(0, dim = c(sim_list$n_yrs, sim_list$n_ages, sim_list$n_ages, sim_list$n_sims))
    for(i in 1:sim_list$n_yrs) for(sim in 1:sim_list$n_sims) diag(identity_AgeingError[i,,,sim]) <- 1 # create identity matrix for each year
    sim_list$AgeingError <- identity_AgeingError
    warning("No ageing error matrix was provided. A default identity matrix was used, which assumes that the number and structure of modelled age bins exactly match the observed age bins. If the observed age composition data includes fewer age bins than the model (e.g., observed ages 2-10 while modelled ages are 1-10), this default assumption will cause a dimensional mismatch and potentially misalign the modelled and observed compositions. To avoid this, please provide an ageing error matrix of dimension n_model_ages x n_obs_ages that correctly maps modelled ages to observed age bins. For example, if observed ages are 2-10, supply a matrix that drops the first model age by using a shifted identity matrix: diag(1, 10)[, 2:10]. This will ensure the age bins are correctly aligned for likelihood calculations.")
  }

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
#'       population–region–year–age–sex cell to its corresponding \code{ln_M}
#'       parameter index.}
#'   }
#'
#' @keywords internal
do_M_mapping <- function(input_list,
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
#' @param addtocomp Small constant added to composition proportions before likelihood
#'   evaluation to avoid \code{log(0)}. Default \code{1e-3}. Ignored when a
#'   logistic-normal likelihood is specified, as that family handles zeros internally.
#' @param addtofishidx Small constant added to fishery indices. Default \code{1e-4}.
#' @param addtosrvidx Small constant added to survey indices. Default \code{1e-4}.
#' @param addtotag Small constant added to tag recovery observations. Default \code{1e-10}.
#' @param AgeingError Ageing error (age–age transition) array mapping true modelled ages
#'   to observed age bins. Accepted forms:
#'   \describe{
#'     \item{2D matrix \code{[n_model_ages × n_obs_ages]}}{Time-invariant ageing error;
#'       replicated internally across all years.}
#'     \item{3D array \code{[n_years × n_model_ages × n_obs_ages]}}{Time-varying ageing error.}
#'     \item{\code{NULL} (default)}{An identity matrix is constructed, assuming modelled
#'       and observed age bins are identical. If observed bins are a subset of modelled
#'       ages (e.g., observed ages 2–10 vs. modelled ages 1–10), supply a shifted
#'       identity matrix such as \code{diag(1, n_model_ages)[, obs_age_index]} to
#'       avoid a dimensional mismatch.}
#'   }
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
#'   Required when \code{fit_lengths = 1}; ignored otherwise.
#' @param Selex_Type Character string specifying the selectivity basis. One of:
#'   \describe{
#'     \item{\code{"age"} (default)}{Age-based selectivity.}
#'     \item{\code{"length"}}{Length-based selectivity. Requires \code{fit_lengths = 1}.}
#'   }
#' @param M_spec Character string controlling natural mortality estimation. One of:
#'   \describe{
#'     \item{\code{"est_ln_M"} (default)}{Estimate \code{ln_M} across the defined blocks.}
#'     \item{\code{"fix"}}{Fix mortality to \code{Fixed_natmort}; \code{ln_M} parameters
#'       are mapped to \code{NA} and not passed to the optimizer.}
#'   }
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
#' @param ... Optional starting value overrides passed by name. Currently recognised:
#'   \describe{
#'     \item{\code{ln_M}}{Array of log-scale starting values for natural mortality,
#'       dimensioned \code{[n_popblks × n_regionblks × n_yearblks × n_ageblks × n_sexblks]}.
#'       Defaults to \code{log(0.5)} for all blocks if not supplied.}
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
                                  addtocomp = 1e-3,
                                  addtofishidx = 1e-4,
                                  addtosrvidx = 1e-4,
                                  addtotag = 1e-10,
                                  AgeingError = NULL,
                                  Use_M_prior = 0,
                                  M_prior = NA,
                                  fit_lengths = 0,
                                  SizeAgeTrans = NA,
                                  Selex_Type = 'age',
                                  M_spec = "est_ln_M",
                                  M_popblk_spec = 'constant',
                                  M_ageblk_spec = 'constant',
                                  M_regionblk_spec = 'constant',
                                  M_yearblk_spec = 'constant',
                                  M_sexblk_spec = 'constant',
                                  Fixed_natmort = NULL,
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Biologicals <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Weight at age checking
  check_data_dimensions(WAA, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_seas = input_list$data$n_seas, what = 'WAA')
  if(!is.null(WAA_fish)) check_data_dimensions(WAA_fish, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'WAA_fish')
  if(!is.null(WAA_srv)) check_data_dimensions(WAA_srv, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'WAA_srv')

  # Maturity at age checking
  check_data_dimensions(MatAA, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas,  n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'MatAA')

  # Length checking
  if(!fit_lengths %in% c(0,1)) stop("Values for fit_lengths are not valid. They are == 0 (not used), or == 1 (used)")
  collect_message("Length Composition data are: ", ifelse(fit_lengths == 0, "Not Used", "Used"))

  # Size Age Transition checking
  if(fit_lengths == 1) check_data_dimensions(SizeAgeTrans, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'SizeAgeTrans')
  if(fit_lengths == 1 & is.na(sum(SizeAgeTrans))) stop("Length composition are fit to, but the size-age transition matrix is NA")

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

  # Checking ageing error dimensions
  if(!is.null(AgeingError)) {
    if(length(dim(AgeingError)) == 2) check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), what = 'AgeingError') # user supplied ageing error is not time-varying
    if(length(dim(AgeingError)) == 3) check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), n_years = length(input_list$data$years), what = 'AgeingError_t') # user supplied ageing error is time-varying
  }

  # Selectivity Options -----------------------------------------------------

  # Age based selectivity
  if(Selex_Type == 'age') {
    Selex_Type <- 0
    collect_message("Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(Selex_Type == 'length') {
    if(fit_lengths == 0) stop("Length composition data are not fit, but selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    Selex_Type <- 1
    collect_message("Selectivity is length-based")
  } # if length based


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
    warning("No ageing error matrix was provided. A default identity matrix was used, which assumes that the number and structure of modelled age bins exactly match the observed age bins. If the observed age composition data includes fewer age bins than the model (e.g., observed ages 2-10 while modelled ages are 1-10), this default assumption will cause a dimensional mismatch and potentially misalign the modelled and observed compositions. To avoid this, please provide an ageing error matrix of dimension n_model_ages x n_obs_ages that correctly maps modelled ages to observed age bins. For example, if observed ages are 2-10, supply a matrix that drops the first model age by using a shifted identity matrix: diag(1, 10)[, 2:10]. This will ensure the age bins are correctly aligned for likelihood calculations.")
  } else if(length(dim(AgeingError)) == 2) {   # setup ageing error if user-supplied is not year specific
    AgeingError_t <- array(0, dim = c(length(input_list$data$years), dim(AgeingError)))
    for(i in 1:length(input_list$data$years)) AgeingError_t[i,,] <- AgeingError
    collect_message("Ageing Error is specified to be time-invariant")
  } else if(length(dim(AgeingError)) == 3) {   # ageing error if it is year specific (just reassigning)
    AgeingError_t <- AgeingError
    collect_message("Ageing Error is specified to be time-varying")
  }


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
  input_list$data$fit_lengths <- fit_lengths
  input_list$data$SizeAgeTrans <- SizeAgeTrans
  input_list$data$Use_M_prior <- Use_M_prior
  input_list$data$M_prior <- M_prior
  input_list$data$Fixed_natmort <- Fixed_natmort
  input_list$data$Selex_Type <- Selex_Type
  input_list$data$addtocomp <- addtocomp
  input_list$data$addtofishidx <- addtofishidx
  input_list$data$addtosrvidx <- addtosrvidx
  input_list$data$addtotag <- addtotag

  # Populate Parameter List -------------------------------------------------

  # If M is constant for ages
  if(is.character(M_ageblk_spec)){
    if(M_ageblk_spec == "constant") M_ageblk_spec_vals <- list(1:length(input_list$data$ages))
  } else M_ageblk_spec_vals <- M_ageblk_spec

  # If M is constant across years
  if(is.character(M_yearblk_spec)){
    if(M_yearblk_spec == "constant") M_yearblk_spec_vals = list(1:length(input_list$data$years))
  } else M_yearblk_spec_vals = M_yearblk_spec

  # If M is constant across sexes
  if(is.character(M_sexblk_spec)){
    if(M_sexblk_spec == "constant") M_sexblk_spec_vals <- list(1:input_list$data$n_sexes)
  } else M_sexblk_spec_vals <- M_sexblk_spec

  # If M is constant across regions
  if(is.character(M_regionblk_spec)){
    if(M_regionblk_spec == "constant") M_regionblk_spec_vals <- list(1:input_list$data$n_regions)
  } else M_regionblk_spec_vals <- M_regionblk_spec

  # If M is constant across populations
  if(is.character(M_popblk_spec)){
    if(M_popblk_spec == "constant") M_popblk_spec_vals <- list(1:input_list$data$n_pop)
  } else M_popblk_spec_vals <- M_popblk_spec

  if("ln_M" %in% names(starting_values)) input_list$par$ln_M <- starting_values$ln_M
  else input_list$par$ln_M <- array(log(0.5), dim = c(length(M_popblk_spec_vals),
                                                      length(M_regionblk_spec_vals),
                                                      length(M_yearblk_spec_vals),
                                                      length(M_ageblk_spec_vals),
                                                      length(M_sexblk_spec_vals)))

  # Mapping Options ---------------------------------------------------------
  input_list <- do_M_mapping(input_list, M_spec, M_popblk_spec_vals, M_regionblk_spec_vals,
                             M_yearblk_spec_vals, M_ageblk_spec_vals, M_sexblk_spec_vals) # natural mortality mapping

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

