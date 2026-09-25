# Stage 1 of 3: model setup
#
# Fishery selectivity, retention and catchability inputs: Setup_Mod_Fishsel_and_Q for the gear's
# selectivity on the stock, Setup_Mod_Retsel for what fraction of encountered fish are kept.

#' Set up retained fishery selectivity
#'
#' Sets the retention selectivity forms, time blocks, continuous time variation,
#' process error, annual deviations and fixed or estimated parameters. Time
#' variation through \code{cont_tv_ret_sel} and blocks through
#' \code{ret_sel_blocks} are mutually exclusive within a fleet. Call after
#' \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}.
#' @param cont_tv_ret_sel Character vector of continuous time variation per fleet,
#'   each \code{"<type>_Fleet_<f>"}: \code{"none"}, \code{"iid"}, \code{"rw"},
#'   \code{"3dmarg"}, \code{"3dcond"} or \code{"2dar1"}. Stored as integer codes.
#' @param ret_sel_blocks Character vector of discrete blocks, either
#'   \code{"none_Fleet_<f>"} or \code{"Block_<b>_Year_<start>-<end>_Fleet_<f>"} with
#'   \code{"terminal"} allowed as the end year. Expanded into an \code{[n_regions ×
#'   n_years × n_fish_fleets]} array.
#' @param ret_sel_model Character vector of the retention form per fleet, and
#'   optionally per block: \code{"<type>_Fleet_<f>"} or
#'   \code{"<type>_Fleet_<f>_Block_<b>"}. The forms and their syntax are those of
#'   \code{fish_sel_model} in \code{\link{Setup_Mod_Fishsel_and_Q}}.
#' @param retsel_pe_pars_spec Process error parameters for time-varying retention,
#'   one entry per fishery fleet.
#' @param ret_fixed_sel_pars_spec Which fixed retention parameters are estimated.
#' @param ret_sel_devs_spec Structure of the annual retention deviations per fleet.
#' @param ret_sel_corr_opt_semipar Optional correlation structure for
#'   semi-parametric retention deviations, one entry per fishery fleet.
#' @param Use_ret_selex_prior Integer flag (0/1) for retention selectivity priors.
#' @param ret_selex_prior Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par}, \code{mu}, \code{sd} and an optional
#'   \code{type} (\code{"par"} or \code{"value"}), as \code{fish_selex_prior} in
#'   \code{\link{Setup_Mod_Fishsel_and_Q}}.
#' @param retsel_devs_shared_bins Bins sharing one retention deviation series.
#' @param ret_sel_bin_dev_bins List with one element per fleet naming the bins that
#'   fleet overrides, or \code{NULL} for none. An overridden bin takes a freely
#'   estimated annual value in place of what the functional form produced, applied
#'   after every other transformation including standardization. Default
#'   \code{NULL}.
#' @param cont_tv_retsel_bin_devs Character vector \code{[n_fish_fleets]} of the
#'   process error on the bin-override deviations: \code{"none"} (default),
#'   \code{"iid"} or \code{"rw"}.
#' @param retsel_pe_wt Numeric vector \code{[n_fish_fleets]} multiplying the
#'   retention process error likelihood. Default \code{1}. \code{0} skips that
#'   fleet's process error, so the deviations stay estimated but enter the objective
#'   only through the data and any smoothness or centering penalties. Values other
#'   than 0 or 1 make an estimated process error sigma reinterpretable. Applies to
#'   \code{ln_retsel_devs} only.
#' @param retsel_rw_init_sigma Numeric vector \code{[n_fish_fleets]} giving the
#'   standard deviation of the first year of an \code{"rw"} deviation series.
#'   Default \code{5}, which leaves that year effectively free. \code{NA} instead
#'   starts the walk at zero under the walk's own estimated sigma.
#' @param ret_selex_type Character scalar, \code{"age"} or \code{"length"}.
#' @param use_fixed_ret_sel Binary array \code{[n_pop × n_regions × n_years ×
#'   n_seas × n_fish_fleets]}, \code{1} to use \code{ret_sel_input}.
#' @param ret_sel_input Fixed retention array \code{[n_pop × n_regions × n_years ×
#'   n_seas × n_bins × n_sexes × n_fish_fleets]}.
#' @param ret_sel_nonpar_est_bins Estimated bins for non-parametric retention.
#' @param ret_sel_sex_offset Character vector \code{[n_fish_fleets]} linking the
#'   sexes of a fleet's retention curve: \code{"none"} (default), \code{"par"},
#'   \code{"scale"} or \code{"par_scale"}, as \code{fish_sel_sex_offset} in
#'   \code{\link{Setup_Mod_Fishsel_and_Q}}. Retention is a fraction, so a scale
#'   offset only makes sense where the scaled curve stays at or below one, and
#'   nothing enforces that.
#' @param ... Optional starting values for the selectivity parameters and
#'   deviations.
#'
#' @return \code{input_list} with the parsed retention structure arrays in
#'   \code{$data}, the starting values in \code{$par} and the factor maps in
#'   \code{$map}.
#'
#' @keywords internal
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Retsel <- function(
  input_list,
  cont_tv_ret_sel,
  ret_sel_blocks,
  ret_sel_model,
  retsel_pe_pars_spec,
  ret_fixed_sel_pars_spec,
  ret_sel_devs_spec,
  ret_sel_corr_opt_semipar,
  Use_ret_selex_prior,
  ret_selex_prior,
  retsel_devs_shared_bins,
  ret_selex_type,
  use_fixed_ret_sel,
  ret_sel_input,
  ret_sel_bin_dev_bins = NULL,
  cont_tv_retsel_bin_devs = rep("none", input_list$data$n_fish_fleets),
  retsel_pe_wt = rep(1, input_list$data$n_fish_fleets),
  retsel_rw_init_sigma = rep(5, input_list$data$n_fish_fleets),
  retsel_dont_est_dev_first = rep(0, input_list$data$n_fish_fleets),
  ret_sel_nonpar_est_bins,
  ret_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
  ...
) {

  messages_list <<- character(0) # string to attach to for printing messages # nolint: object_usage_linter.
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Fishsel_and_Q <- c(input_list$config$Setup_Mod_Fishsel_and_Q, mget(names(formals()))[-1])


  # Input Validation --------------------------------------------------------
  # Continuous Selectivity Deviations
  check_fleet_spec_length(retsel_pe_pars_spec, input_list$data$n_fish_fleets, "retsel_pe_pars_spec", allow_null = TRUE)
  check_fleet_spec_length(ret_sel_devs_spec, input_list$data$n_fish_fleets, "ret_sel_devs_spec", allow_null = TRUE)
  check_fleet_spec_length(ret_sel_corr_opt_semipar, input_list$data$n_fish_fleets, "ret_sel_corr_opt_semipar", allow_null = TRUE)

  # A short vector here is read per fleet in the objective, so a length mismatch
  # silently becomes NA rather than being recycled.
  check_fleet_spec_length(retsel_pe_wt, input_list$data$n_fish_fleets, "retsel_pe_wt")
  check_fleet_spec_length(retsel_rw_init_sigma, input_list$data$n_fish_fleets, "retsel_rw_init_sigma")
  check_fleet_spec_length(retsel_dont_est_dev_first, input_list$data$n_fish_fleets, "retsel_dont_est_dev_first")
  if(!all(retsel_dont_est_dev_first %in% c(0, 1))) stop("retsel_dont_est_dev_first must be 0 or 1 for every fleet")

  # Selectivity Priors
  if(!Use_ret_selex_prior %in% c(0,1)) stop("Values for Use_ret_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_ret_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(ret_selex_prior))
    if(length(missing_cols) > 0) {
      stop("ret_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Retained Fishery Selectivity priors are: ", ifelse(Use_ret_selex_prior == 0, "Not Used", "Used"))

  if(any(use_fixed_ret_sel == 1) && is.null(ret_sel_input)) stop("ret_sel_input is NULL, please provide an input array.")
  if(any(use_fixed_ret_sel == 1) && ret_selex_type == 'age') check_data_dimensions(
    ret_sel_input,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_ages = length(input_list$data$ages),
    n_sexes = input_list$data$n_sexes,
    n_fish_fleets = input_list$data$n_fish_fleets,
    what = 'ret_sel_input_age'
  )
  if(any(use_fixed_ret_sel == 1) && ret_selex_type == 'length') check_data_dimensions(
    ret_sel_input,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_lens = length(input_list$data$lens),
    n_sexes = input_list$data$n_sexes,
    n_fish_fleets = input_list$data$n_fish_fleets,
    what = 'ret_sel_input_len'
  )

  # Selectivity Options -----------------------------------------------------
  # The bin vector is kept as well as its length. Starting values stated on the
  # bin scale are seeded further down, by which point ret_selex_type holds the
  # numeric code rather than the name it arrived as.
  if(ret_selex_type == 'age') {
    ret_selex_type <- 0
    ret_sel_bin_vec <- input_list$data$ages
    collect_message("Retained Fishery Selectivity is aged-based.")
  } else if(ret_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but retained selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    ret_selex_type <- 1
    ret_sel_bin_vec <- input_list$data$lens
    collect_message("Retained Fishery Selectivity is length-based")
  } else stop("ret_selex_type must be 'age' or 'length', but was: ", ret_selex_type)

  bins <- length(ret_sel_bin_vec)


  # Continuous Retained Time-Varying Selectivity Options -----------------------------
  cont_tv_ret_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in seq_along(cont_tv_ret_sel)) {
    # Extract out components from list
    tmp <- cont_tv_ret_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for cont_tv_ret_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_ret_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_ret_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous retained fishery time-varying selectivity specified as: ", cont_tv_type, " for fishery fleet ", fleet)
  }

  if(any(cont_tv_ret_sel_mat > 0) && (is.null(retsel_pe_pars_spec) || is.null(ret_sel_devs_spec))) stop("Continuous time-varying selectivity specified, but retsel_pe_pars_spec and/or ret_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Retained Time-Varying Selectivity Options --------------------------------
  ret_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in seq_along(ret_sel_blocks)) {

    # Extract out components from list
    tmp <- ret_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Retained Fishery Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      ret_sel_blocks_arr[,,fleet] <- 1 # input only 1 fishery time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # extract fleet index

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      ret_sel_blocks_arr[,years,fleet] <- block_val
    }

  }

  if(any(is.na(ret_sel_blocks_arr))) stop("Retained Fishery Selectivtiy Blocks are returning an NA. Did you forget to specify the year range of ret_sel_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Retained Fishery Selectivity Time Blocks for fishery", f, "is specified at:", length(unique(ret_sel_blocks_arr[,,f]))))

  # Retained Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml", 'nonpar', 'asymplogist1', "asymplogist2", "bicubic", "nonparlog", "nonparfree"), num = c(0,1,2,3,4,5,6,7,8,9,10)) # set up values we can map to
  ret_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  ret_sel_bicubic_binnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bin nodes, only set where ret_sel_model == 8 (bicubic)
  ret_sel_bicubic_yrnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of year nodes, only set where ret_sel_model == 8 (bicubic)
  ret_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year); years before this are edge-kept, matching fish_sel_model's SelStyr
  ret_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bins the bicubic surface is actually fit over (0 = all bins); bins beyond this are edge-kept, matching fish_sel_model's NSelBins

  for(i in seq_along(ret_sel_model)) {

    # Extract out retained fishery selectivity components from vector
    tmp_sel_form <- ret_sel_model[i]
    tmp_sel_form_vec <- unlist(strsplit(tmp_sel_form, "_")) # split string
    sel_form <- tmp_sel_form_vec[1] # get selectivity type

    if(sel_form == "bicubic") {
      # bicubic spline: bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f>[_Block_<b>][_SelStyr_<year>][_NSelBins_<n>]
      bin_pos <- which(tmp_sel_form_vec == "Bin")
      yr_pos <- which(tmp_sel_form_vec == "Yr")
      fleet_pos <- which(tmp_sel_form_vec == "Fleet")
      block_pos <- which(tmp_sel_form_vec == "Block")
      selstyr_pos <- which(tmp_sel_form_vec == "SelStyr")
      nselbins_pos <- which(tmp_sel_form_vec == "NSelBins")
      if(length(bin_pos) != 1 || length(yr_pos) != 1 || length(fleet_pos) != 1)
        stop("ret_sel_model 'bicubic' entries must be specified as bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f> or bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Block_<b>_Fleet_<f>, optionally with _SelStyr_<year> and/or _NSelBins_<n>")
      tmp_n_bin_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[bin_pos + 1]))
      tmp_n_yr_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[yr_pos + 1]))
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_selstyr <- if(length(selstyr_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[selstyr_pos + 1])) else 0
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(is.na(tmp_n_bin_nodes) || tmp_n_bin_nodes < 2) stop("bicubic ret_sel_model requires at least 2 bin nodes (n_bin_nodes >= 2)")
      if(is.na(tmp_n_yr_nodes) || tmp_n_yr_nodes < 1) stop("bicubic ret_sel_model requires at least 1 year node (n_yr_nodes >= 1). Use n_yr_nodes == 1 for a time-invariant bin-only spline.")
      if(length(selstyr_pos) == 1 && (is.na(tmp_selstyr) || !tmp_selstyr %in% input_list$data$years)) stop("bicubic ret_sel_model SelStyr must be a calendar year within the modeled years")
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("bicubic ret_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
    } else {
      # optional _NSelBins_<n> plateau suffix is joined with  _Block_<b> in either order for any parametric form
      fleet_pos <- which(tmp_sel_form_vec == "Fleet")
      block_pos <- which(tmp_sel_form_vec == "Block")
      nselbins_pos <- which(tmp_sel_form_vec == "NSelBins")
      if(length(fleet_pos) != 1) stop("ret_sel_model entries must name their fleet exactly once, as _Fleet_<f>")
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("ret_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
      if(length(nselbins_pos) == 1 && sel_form %in% c("nonpar", "nonparlog", "nonparfree")) stop("ret_sel_model NSelBins is for the parametric forms; a non-parametric form would keep estimating the bins it then overwrites. Group those bins through ret_sel_nonpar_est_bins instead.")
    }

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("ret_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml, nonpar, nonparlog, nonparfree, asymplogist1, asymplogist2, bicubic (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for ret_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) {
      ret_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
      if(sel_form == "bicubic") {
        ret_sel_bicubic_binnodes_arr[,,tmp_fleet] <- tmp_n_bin_nodes
        ret_sel_bicubic_yrnodes_arr[,,tmp_fleet] <- tmp_n_yr_nodes
        ret_sel_bicubic_selstyr_arr[,,tmp_fleet] <- tmp_selstyr
      }
      ret_sel_bicubic_nselbins_arr[,,tmp_fleet] <- tmp_nselbins # plateau bin; read by every form (0 = none)
    } else {
      ret_sel_model_arr <- assign_sel_block(ret_sel_model_arr, ret_sel_blocks_arr, tmp_fleet, tmp_block, sel_map$num[which(sel_map$sel == sel_form)])
      if(sel_form == "bicubic") {
        ret_sel_bicubic_binnodes_arr <- assign_sel_block(ret_sel_bicubic_binnodes_arr, ret_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_bin_nodes)
        ret_sel_bicubic_yrnodes_arr <- assign_sel_block(ret_sel_bicubic_yrnodes_arr, ret_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_yr_nodes)
        ret_sel_bicubic_selstyr_arr <- assign_sel_block(ret_sel_bicubic_selstyr_arr, ret_sel_blocks_arr, tmp_fleet, tmp_block, tmp_selstyr)
      }
      ret_sel_bicubic_nselbins_arr <- assign_sel_block(ret_sel_bicubic_nselbins_arr, ret_sel_blocks_arr, tmp_fleet, tmp_block, tmp_nselbins)
    }
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Retained Fishery selectivity functional form specified as:", sel_form, " for fishery fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_fish_fleets) {
    has_blocks <- length(unique(ret_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_ret_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$cont_tv_ret_sel <- cont_tv_ret_sel_mat
  input_list$data$ret_sel_blocks <- ret_sel_blocks_arr
  input_list$data$ret_sel_model <- ret_sel_model_arr
  input_list$data$ret_sel_bicubic_binnodes <- ret_sel_bicubic_binnodes_arr
  input_list$data$ret_sel_bicubic_yrnodes <- ret_sel_bicubic_yrnodes_arr
  input_list$data$ret_sel_bicubic_selstyr <- ret_sel_bicubic_selstyr_arr
  input_list$data$ret_sel_bicubic_nselbins <- ret_sel_bicubic_nselbins_arr
  input_list$data$Use_ret_selex_prior <- Use_ret_selex_prior
  input_list$data$ret_selex_prior <- validate_selex_prior_types(ret_selex_prior, Use_ret_selex_prior, "ret_selex_prior",
                                                                sel_blocks = ret_sel_blocks_arr, n_bins = bins)
  input_list$data$ret_selex_type <- ret_selex_type
  input_list$data$use_fixed_ret_sel <- use_fixed_ret_sel
  input_list$data$ret_sel_input <- ret_sel_input
  input_list$data$retsel_pe_wt <- retsel_pe_wt
  input_list$data$retsel_rw_init_sigma <- retsel_rw_init_sigma
  input_list <- setup_sel_bin_devs(
    input_list,
    ret_sel_bin_dev_bins,
    cont_tv_retsel_bin_devs,
    prefix = "ret",
    n_fleets = input_list$data$n_fish_fleets,
    bins = bins,
    starting_values = starting_values
  )
  input_list$data$retsel_devs_min_shared_bins <- if(!is.null(retsel_devs_shared_bins)) unlist(lapply(retsel_devs_shared_bins, min)) else seq_along(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------

  # Figure out number of selectivity parameters for a given functional form
  unique_retsel_vals <- unique(as.vector(input_list$data$ret_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in seq_along(unique_retsel_vals)) {
    if(unique_retsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_retsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_retsel_vals[i] %in% c(4)) sel_pars_vec[i] <- 6 # double normal
    if(unique_retsel_vals[i] %in% c(5,9,10)) sel_pars_vec[i] <- bins # non-parametric
    if(unique_retsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic w/ asymptote parameter
    if(unique_retsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$ret_sel_bicubic_binnodes * input_list$data$ret_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  # figure out maximum number of retained fishery selectivity blocks for a given reigon and fleet
  max_retsel_blks <- max(apply(input_list$data$ret_sel_blocks, c(1,3), FUN = function(x) length(unique(x))))
  # maximum number of selectivity parameters across all forms
  max_retsel_pars <- max(sel_pars_vec)
  input_list$par$ret_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_retsel_pars, max_retsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))
  input_list$par$ret_fixed_sel_pars <- use_starting_value(input_list$par$ret_fixed_sel_pars, starting_values, "ret_fixed_sel_pars")

  # a double normal has its peak on the bin scale, so a default of zero
  # would put it at bin zero, where the ascending limb has no extent
  if(!"ret_fixed_sel_pars" %in% names(starting_values)) {
    input_list$par$ret_fixed_sel_pars <- seed_dbnrml_peak(input_list$par$ret_fixed_sel_pars, ret_sel_model_arr,
                                          ret_sel_bin_vec,
                                          ret_sel_sex_offset)
  }

  # Bicubic spline interpolation weight matrices (bin node x year node grid) for retention selectivity,
  # mirroring fish_sel_bicubic_Wbin/Wyr above (see Get_Selex documentation).
  has_bicubic_ret_sel <- any(input_list$data$ret_sel_model == 8)
  max_bin_nodes_bicubic_ret <- if(has_bicubic_ret_sel) max(input_list$data$ret_sel_bicubic_binnodes) else 1
  max_yr_nodes_bicubic_ret <- if(has_bicubic_ret_sel) max(input_list$data$ret_sel_bicubic_yrnodes) else 1
  n_yrs_total_bicubic_ret <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs

  ret_sel_bicubic_Wbin <- array(0, dim = c(input_list$data$n_regions, bins, max_bin_nodes_bicubic_ret, max_retsel_blks, input_list$data$n_fish_fleets))
  ret_sel_bicubic_Wyr <- array(0, dim = c(input_list$data$n_regions, n_yrs_total_bicubic_ret, max_yr_nodes_bicubic_ret, max_retsel_blks, input_list$data$n_fish_fleets))

  if(has_bicubic_ret_sel) {
    for(f in 1:input_list$data$n_fish_fleets) {
      for(r in 1:input_list$data$n_regions) {

        retsel_blocks_tmp <- unique(as.vector(input_list$data$ret_sel_blocks[r,,f]))

        for(b in seq_along(retsel_blocks_tmp)) {

          block_years <- which(input_list$data$ret_sel_blocks[r,,f] == retsel_blocks_tmp[b])
          if(unique(input_list$data$ret_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$ret_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$ret_sel_bicubic_yrnodes[r, block_years, f])

          # Bin dimension: see fish_sel_bicubic_Wbin construction above for NSelBins truncation/plateau details
          nselbins_this <- unique(input_list$data$ret_sel_bicubic_nselbins[r, block_years, f])
          n_fit_bins <- if(nselbins_this == 0) bins else nselbins_this

          bin_nodes_scaled <- seq(0, 1, length.out = n_bin_nodes_this)
          fit_bin_scaled <- seq(0, 1, length.out = n_fit_bins)
          Wbin_fit <- Get_Natural_Cubic_Spline_Weights(bin_nodes_scaled, fit_bin_scaled)

          Wbin_this <- matrix(0, nrow = bins, ncol = n_bin_nodes_this)
          Wbin_this[1:n_fit_bins, ] <- Wbin_fit
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(
            Wbin_fit[nrow(Wbin_fit), ],
            nrow = bins - n_fit_bins,
            ncol = n_bin_nodes_this,
            byrow = TRUE
          )

          ret_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # Year dimension: see fish_sel_bicubic_Wyr construction above for SelStyr truncation/edge-hold details
          selstyr_this <- unique(input_list$data$ret_sel_bicubic_selstyr[r, block_years, f])
          selstyr_idx <- if(selstyr_this == 0) min(block_years) else which(input_list$data$years == selstyr_this)
          fit_years <- block_years[block_years >= selstyr_idx]
          pre_fit_years <- block_years[block_years < selstyr_idx]

          yr_nodes_scaled <- seq(0, 1, length.out = n_yr_nodes_this)
          fit_yr_scaled <- seq(0, 1, length.out = length(fit_years))
          Wyr_block <- Get_Natural_Cubic_Spline_Weights(yr_nodes_scaled, fit_yr_scaled)

          Wyr_this <- matrix(0, nrow = n_yrs_total_bicubic_ret, ncol = n_yr_nodes_this)
          Wyr_this[fit_years, ] <- Wyr_block
          if(length(pre_fit_years) > 0) Wyr_this[pre_fit_years, ] <- matrix(Wyr_block[1, ], nrow = length(pre_fit_years), ncol = n_yr_nodes_this, byrow = TRUE)
          if(min(block_years) > 1) Wyr_this[1:(min(block_years) - 1), ] <- matrix(Wyr_block[1, ], nrow = min(block_years) - 1, ncol = n_yr_nodes_this, byrow = TRUE)
          if(max(block_years) < n_yrs_total_bicubic_ret) Wyr_this[(max(block_years) + 1):n_yrs_total_bicubic_ret, ] <- matrix(
            Wyr_block[nrow(Wyr_block), ],
            nrow = n_yrs_total_bicubic_ret - max(block_years),
            ncol = n_yr_nodes_this,
            byrow = TRUE
          )

          ret_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this
        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_ret_sel

  input_list$data$ret_sel_bicubic_Wbin <- ret_sel_bicubic_Wbin
  input_list$data$ret_sel_bicubic_Wyr <- ret_sel_bicubic_Wyr

  # Retained Fishery selectivity process error parameters
  input_list$par$retsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_retsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using
  input_list$par$retsel_pe_pars <- use_starting_value(input_list$par$retsel_pe_pars, starting_values, "retsel_pe_pars")

  # Retained Fishery selectivity deviations
  if(input_list$data$ret_selex_type == 0) bins <- length(input_list$data$ages) # age based deviations
  if(input_list$data$ret_selex_type == 1) bins <- length(input_list$data$lens) # length based deviations
  input_list$par$ln_retsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_fish_fleets))
  input_list$par$ln_retsel_devs <- use_starting_value(input_list$par$ln_retsel_devs, starting_values, "ln_retsel_devs")

  # Sex offsets on retention selectivity (parameter offsets and/or a curve scale offset)
  input_list <- setup_sel_sex_offset(
    input_list,
    ret_sel_sex_offset,
    prefix = "ret",
    n_fleets = input_list$data$n_fish_fleets,
    fleet_label = "fishery fleet (retention)",
    sel_model_arr = input_list$data$ret_sel_model,
    cont_tv_mat = cont_tv_ret_sel_mat,
    max_blks = max_retsel_blks,
    sel_blocks = input_list$data$ret_sel_blocks,
    fixed_spec = ret_fixed_sel_pars_spec,
    starting_values = starting_values
  )

  # Mapping Options ---------------------------------------------------------
  input_list <- do_fixed_sel_pars_mapping(
    input_list,
    ret_fixed_sel_pars_spec,
    bins,
    ret_sel_nonpar_est_bins,
    prefix = "ret",
    fleet_field = "n_fish_fleets",
    use_field = "Catch",
    fleet_label = "fishery fleet"
  )
  input_list <- do_sel_pe_pars_mapping(
    input_list,
    retsel_pe_pars_spec,
    ret_sel_corr_opt_semipar,
    bins,
    sel_devs_spec = ret_sel_devs_spec,
    sel_devs_shared_bins = retsel_devs_shared_bins,
    prefix = "ret",
    fleet_field = "n_fish_fleets",
    use_field = "Catch",
    fleet_label = "fishery fleet"
  )
  input_list <- do_sel_devs_mapping(
    input_list,
    ret_sel_devs_spec,
    retsel_devs_shared_bins,
    bins,
    dont_est_dev_first = retsel_dont_est_dev_first,
    prefix = "ret",
    fleet_field = "n_fish_fleets",
    use_field = "Catch",
    fleet_label = "fishery fleet"
  )

   return(input_list)
}

#' Set up total and retained fishery selectivity and catchability specifications
#'
#' Sets the selectivity functional forms, time blocks, continuous time variation,
#' process error hyperparameters, annual deviations, and the catchability blocks
#' and estimation structure, for total and retained selectivity alike. Time
#' variation and blocked selectivity are mutually exclusive within a fleet. Call
#' after \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}.
#' @param fish_sel_model Character vector of the selectivity form per fleet, and
#'   optionally per block: \code{"<model>_Fleet_<f>"} or
#'   \code{"<model>_Fleet_<f>_Block_<b>"}. The forms are \code{"logist1"}
#'   (\eqn{a_{50}} and slope), \code{"logist2"} (\eqn{a_{50}} and \eqn{a_{95}}),
#'   \code{"gamma"} (dome, \eqn{a_{max}} and \eqn{\delta}), \code{"exponential"}
#'   (one power), \code{"dbnrml"} (double normal, six parameters),
#'   \code{"asymplogist1"} and \code{"asymplogist2"} (the two logistics with an
#'   asymptote), the three non-parametric forms and \code{"bicubic"}.
#'
#'   \code{"nonpar"} is on the logit scale, mean-standardized jointly over years
#'   and bins so the grand mean of the surface is one. \code{"nonparlog"} is on the
#'   log scale, standardized so each year averages to one over
#'   \code{*_sel_norm_bins}, leaving only within-year contrasts identified.
#'   \code{"nonparfree"} is on the log scale with no standardization,
#'   \eqn{\exp(\theta)}, so the values hold the height of the curve as well as its
#'   shape; this is the form for a data source fit age by age, where a free
#'   catchability per age and a selectivity estimated at age are one quantity
#'   written twice, so no catchability is set. Pin one bin, by leaving it out of
#'   the estimated bins, whenever the mean it multiplies is also free.
#'
#'   \code{"bicubic"} is a spline over a bin-node by year-node grid, written
#'   \code{"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"} with an optional
#'   \code{_Block_k}. One form covers a smooth bin by year surface
#'   (\code{n_yr_nodes > 1}), a time-invariant bin-only spline
#'   (\code{n_yr_nodes == 1}), and a bin-only spline re-fit per block. An optional
#'   \code{_SelStyr_<year>} restricts the fit to \code{SelStyr}:block-end, holding
#'   earlier years of the block at the \code{SelStyr} curve, and an optional
#'   \code{_NSelBins_<n>} restricts it to the first \code{n} bins, holding the rest
#'   at the last fitted bin. See \code{\link{Get_Selex}} and the model equations
#'   vignette.
#' @param cont_tv_fish_sel Character vector \code{[n_fish_fleets]} of continuous
#'   time variation per fleet, each \code{"<type>_Fleet_<f>"}: \code{"none"}
#'   (default), \code{"iid"} or \code{"rw"} on the selectivity parameters, or
#'   \code{"3dmarg"}, \code{"3dcond"} (3D GMRF on the marginal or conditional
#'   variance) and \code{"2dar1"} (separable over bin and year). Any fleet other
#'   than \code{"none"} also needs \code{fishsel_pe_pars_spec} and
#'   \code{fish_sel_devs_spec}.
#' @param fish_sel_blocks Character vector of discrete selectivity time blocks per
#'   fleet, each \code{"Block_<b>_Year_<s>-<e>_Fleet_<f>"} with \code{"terminal"}
#'   allowed as the end year, or \code{"none_Fleet_<f>"} (default) for one constant
#'   block. Blocks must not overlap and together must span every model year for
#'   that fleet. Mutually exclusive with \code{cont_tv_fish_sel != "none"}.
#' @param fish_q_blocks Catchability time blocks per fleet, in the same format as
#'   \code{fish_sel_blocks}. Default one constant block.
#' @param fish_fixed_sel_pars_spec Character vector \code{[n_fish_fleets]} of how
#'   the fixed-effect selectivity parameters are estimated: \code{"est_all"},
#'   \code{"est_shared_r"}, \code{"est_shared_s"}, \code{"est_shared_r_s"},
#'   \code{"est_shared_f_x"} or \code{"fix"}. See
#'   \code{\link{do_fixed_sel_pars_mapping}}.
#' @param fish_q_spec Character vector \code{[n_fish_fleets]} of the catchability
#'   estimation structure: \code{"est_all"}, \code{"est_shared_r"} or \code{"fix"}.
#'   See \code{\link{do_q_mapping}}.
#' @param fishsel_pe_pars_spec Character vector \code{[n_fish_fleets]} of the
#'   estimation structure for the selectivity process error hyperparameters,
#'   required when any fleet varies continuously. See
#'   \code{\link{do_sel_pe_pars_mapping}}.
#' @param fish_sel_devs_spec Character vector \code{[n_fish_fleets]} of the
#'   estimation structure for the annual selectivity deviations, required when any
#'   fleet varies continuously. See \code{\link{do_sel_devs_mapping}}.
#' @param Use_fish_selex_penalty Integer (0/1). Whether a centering penalty is
#'   applied to sets of fishery selectivity fixed-effect parameters. Default
#'   \code{0}.
#' @param fish_selex_penalty Data frame of centering penalties with columns
#'   \code{region}, \code{fleet}, \code{block}, \code{sex}, \code{par} and
#'   \code{wt}, required when \code{Use_fish_selex_penalty = 1}. Each row penalizes
#'   \code{wt * (log(mean(exp(pars))))^2} over the parameters named in \code{par},
#'   a single index or a list column of integer vectors. This pins the scalar of a
#'   non-parametric curve that catchability or fishing mortality would otherwise
#'   absorb, and is softer than fixing a bin. Meant for parameter sets on the log
#'   scale. Default \code{NULL}.
#' @param fishsel_devs_shared_bins List of integer vectors grouping the bins that
#'   share one deviation series, e.g. \code{list(1:5, 6:10, 11:30)}. Only read when
#'   \code{fish_sel_devs_spec} names an \code{"est_shared_b"} variant.
#' @param corr_opt_semipar Character vector \code{[n_fish_fleets]} of which
#'   correlation components to suppress under 3D GMRF or 2D AR1 time variation.
#'   \code{NA} (default) suppresses none, and the cohort options are invalid for
#'   \code{"2dar1"}. See \code{\link{do_sel_pe_pars_mapping}}.
#' @param fish_q_type Character vector \code{[n_fish_fleets]} of how catchability
#'   is obtained. \code{"est"} (default) estimates \code{ln_fish_q},
#'   \code{"arith"} concentrates it out as the ratio of mean observed to mean
#'   predicted index, and \code{"geo"} does the same on the log scale as
#'   \code{exp(mean(log(obs) - log(pred)))}. Both analytic forms use the years with
#'   observations only and fix that fleet's \code{ln_fish_q} whatever
#'   \code{fish_q_spec} says. The solve runs within each \code{fish_q_blocks}
#'   block, so a blocked catchability gets one solved value per block.
#' @param fish_q_model Character vector \code{[n_fish_fleets]} of the process error
#'   on annual catchability deviations: \code{"none"} (default), \code{"iid"},
#'   \code{"rw"}, \code{"ar1"} or \code{"dsem"}, which hands the series to
#'   \code{\link{Setup_Mod_DSEM}}. Catchability is then \eqn{\exp(\ln q_{r,b,f} +
#'   \epsilon_{r,y,f})}. A fleet with deviations cannot also have
#'   \code{fish_q_blocks} or an analytically solved \code{fish_q_type}.
#' @param sigma_fish_q_spec Sharing string for the deviation standard deviation
#'   over region and fleet: \code{"est_all"} (default), \code{"est_shared_r"},
#'   \code{"est_shared_f"}, \code{"est_shared_r_f"} or \code{"fix"}.
#' @param fish_q_rho_spec Sharing string for the AR1 correlation, with the same
#'   options as \code{sigma_fish_q_spec}. Default \code{"est_all"}. Only read under
#'   \code{fish_q_model = "ar1"}.
#' @param fish_q_rw_init_sigma Standard deviation of the first estimated year of a
#'   random walk. \code{NA} (default) starts the walk at zero under its own sigma,
#'   which keeps \code{ln_fish_q} as the level of the series.
#' @param fishsel_pe_wt Numeric vector \code{[n_fish_fleets]} multiplying the
#'   fishery selectivity process error likelihood. Default \code{1}. \code{0} skips
#'   that fleet's process error, so the deviations stay estimated but enter the
#'   objective only through the data and any smoothness or centering penalties.
#'   Values other than 0 or 1 make an estimated process error sigma
#'   reinterpretable. Applies to \code{ln_fishsel_devs} only; the bin-override
#'   deviations have their own process error.
#' @param fishsel_rw_init_sigma Numeric vector \code{[n_fish_fleets]} giving the
#'   standard deviation of the first year of an \code{"rw"} deviation series.
#'   Default \code{5}, which leaves that year effectively free. \code{NA} instead
#'   starts the walk at zero under the walk's own estimated sigma, which suits a
#'   base curve that already describes the first year well.
#' @param fish_sel_norm_bins List with one element per fleet naming the bins the
#'   mean-one standardization averages over, or \code{NULL} for fleets
#'   standardizing over every bin. Read under \code{"nonparlog"} only. A gear whose
#'   catchability is defined against part of the bin range standardizes over that
#'   part, and catchability absorbs the difference in scale. Default \code{NULL}.
#' @param fish_sel_bin_dev_bins List with one element per fleet naming the bins
#'   that fleet overrides, or \code{NULL} for none, e.g. \code{list(1, NULL)}. An
#'   overridden bin takes a freely estimated annual value
#'   \eqn{\exp(\epsilon_{y,b})} in place of what the functional form produced,
#'   applied after every other transformation including standardization, while the
#'   rest of the curve keeps its parametric shape. Default \code{NULL}.
#' @param cont_tv_fishsel_bin_devs Character vector \code{[n_fish_fleets]} of the
#'   process error on the bin-override deviations: \code{"none"} (default),
#'   \code{"iid"} or \code{"rw"}. A walk has its own estimated sigma per bin.
#' @param fishsel_dont_est_dev_first Integer vector \code{[n_fish_fleets]} of 0/1,
#'   default \code{0}. Where \code{1}, that fleet's deviations start in year two and
#'   the fixed parameters hold year one. A non-parametric form has one free base
#'   parameter per bin, so year one's deviation is that same value written twice
#'   with only \code{fishsel_rw_init_sigma} between them, a prior on a level usually
#'   meant to be free. Refused for the GMRF and 2D AR1 forms, whose deviations are a
#'   field over years and bins rather than a walk anchored at year one.
#' @param Use_fish_q_prior Integer flag, \code{1} for lognormal priors on
#'   catchability. Default \code{0}.
#' @param fish_q_prior Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{mu} on the natural scale and \code{sd} on the log scale,
#'   one row per \eqn{\text{Normal}(\log(\mu), \sigma)} prior. Read when
#'   \code{Use_fish_q_prior = 1}.
#' @param Use_fish_selex_prior Integer flag, \code{1} for priors on the selectivity
#'   parameters. Default \code{0}.
#' @param fish_selex_prior Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par}, \code{mu}, \code{sd} and an optional
#'   \code{type}. \code{"par"} (the default) is a lognormal prior on one fixed
#'   selectivity parameter, with \code{mu} on the natural scale and \code{sd} on the
#'   log scale. \code{"value"} is a normal prior on the realized selectivity at one
#'   bin, both on the natural scale, where \code{par} names the bin and the value is
#'   read at the first model year of \code{block}; that is the ADMB convention of
#'   pinning selectivity at a reference age near one, which no set of independent
#'   parameter priors can express. Read when \code{Use_fish_selex_prior = 1}.
#' @param cont_tv_ret_sel Continuous time variation on retention, with the options
#'   and requirements of \code{cont_tv_fish_sel}. Any fleet other than
#'   \code{"none"} also needs \code{retsel_pe_pars_spec} and
#'   \code{ret_sel_devs_spec}.
#' @param ret_sel_blocks Discrete retention time blocks per fleet, in the format of
#'   \code{fish_sel_blocks} and mutually exclusive with
#'   \code{cont_tv_ret_sel != "none"}.
#' @param ret_sel_model Retention selectivity form per fleet and block, with the
#'   same syntax and the same set of forms as \code{fish_sel_model}.
#' @param retsel_pe_pars_spec Estimation structure for the retention process error
#'   hyperparameters, as \code{fishsel_pe_pars_spec}.
#' @param ret_fixed_sel_pars_spec How the retention fixed-effect parameters are
#'   estimated, with the options of \code{fish_fixed_sel_pars_spec}.
#' @param ret_sel_devs_spec Estimation structure for the annual retention
#'   deviations, as \code{fish_sel_devs_spec}.
#' @param retsel_devs_shared_bins Bins sharing one retention deviation series, as
#'   \code{fishsel_devs_shared_bins}.
#' @param retsel_pe_wt Per-fleet multiplier on the retention process error
#'   likelihood, as \code{fishsel_pe_wt}. Default \code{1}.
#' @param retsel_rw_init_sigma Standard deviation of the first year of an
#'   \code{"rw"} retention deviation series, as \code{fishsel_rw_init_sigma}.
#'   Default \code{5}.
#' @param retsel_dont_est_dev_first Whether each fleet's retention deviations start
#'   in year two, as \code{fishsel_dont_est_dev_first}. Default \code{0}.
#' @param ret_sel_corr_opt_semipar Which correlation components to suppress under
#'   semi-parametric retention time variation, as \code{corr_opt_semipar}.
#' @param Use_ret_selex_prior Integer flag, \code{1} for priors on the retention
#'   selectivity parameters. Default \code{0}.
#' @param ret_selex_prior Data frame with the columns and the optional \code{type}
#'   of \code{fish_selex_prior}.
#' @param use_fixed_ret_sel Integer vector \code{[n_fish_fleets]}, \code{1} to fix
#'   retention selectivity and \code{0} to estimate it.
#' @param ret_sel_input Array of fixed retention values \code{[n_pop × n_regions ×
#'   n_years × n_seas × n_bins × n_sexes × n_fish_fleets]}.
#' @param ret_sel_nonpar_est_bins Optional bin groupings for non-parametric
#'   retention, structured \code{[[fleet]][[block]]}, each element a list of bin
#'   index vectors defining grouped parameters.
#' @param ret_sel_sex_offset Character vector \code{[n_fish_fleets]} linking the
#'   sexes of a fleet's retention curve, with the options of
#'   \code{fish_sel_sex_offset}. Default \code{"none"}. Retention is a fraction, so
#'   a scale offset only makes sense where the scaled curve stays at or below one.
#' @param fish_selex_type Character scalar, \code{"age"} or \code{"length"}, the
#'   bin dim every fishery selectivity function is defined over.
#' @param use_fixed_fish_sel Integer vector \code{[n_fish_fleets]}, \code{1} to fix
#'   fishery selectivity and \code{0} to estimate it.
#' @param fish_sel_input Array of fixed fishery selectivity values \code{[n_pop ×
#'   n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]}. Required
#'   when any \code{use_fixed_fish_sel == 1}.
#' @param fish_sel_nonpar_est_bins Optional bin groupings for non-parametric
#'   fishery selectivity, structured \code{[[fleet]][[block]]}, each element a list
#'   of integer vectors naming the bins that share one estimated parameter. Indices
#'   are on the bin dim \code{fish_selex_type} names.
#' @param fish_sel_dbnrml_startbin \code{NULL} (default) or an integer vector
#'   \code{[n_fish_fleets]}, the bin each fleet's double normal anchors its
#'   ascending limb at. Bins below it take the squared ratio of their bin to it
#'   times the selectivity there, which is Stock Synthesis's convention when the
#'   compositions start above the population's first length bin.
#' @param fish_sel_dbnrml_raw \code{NULL} (default) or a 0/1 matrix
#'   \code{[n_fish_fleets x 2]} for fleets on the double normal: column one leaves
#'   the ascending limb a raw Gaussian instead of anchoring it to \code{p5} at the
#'   first bin, column two does the same for the descending limb and \code{p6}.
#' @param fish_sel_sex_offset Character vector \code{[n_fish_fleets]} linking the
#'   sexes of a fleet's selectivity when \code{n_sexes > 1}. \code{"none"} (default)
#'   keeps each sex's stored parameters its own. \code{"par"} makes every sex beyond
#'   the first hold additive offsets on the first sex's stored parameters, so a
#'   log-scale parameter's natural value is the first sex's times \eqn{e^{\delta}};
#'   offsets fixed at zero reproduce sex-shared parameters. \code{"scale"} keeps
#'   each sex's own parameters and adds a constant log-scale offset on the whole
#'   realized curve, \code{exp(ln_fishsel_sex_scale)}, per region, block and sex,
#'   which may exceed one and is refused for the non-parametric forms and
#'   semi-parametric time variation, whose standardization would cancel it.
#'   \code{"apical"} has the double normal build its limbs up to
#'   \code{exp(ln_*sel_sex_scale)} rather than one, so the offset moves the middle
#'   of the curve and leaves its ends where that sex's own parameters put them.
#'   \code{"par_apical"} and \code{"par_scale"} combine a par offset with each.
#' @param ret_selex_type Character scalar, \code{"age"} or \code{"length"}, the bin
#'   dim every retention selectivity function is defined over.
#' @param ... Optional starting values for the selectivity parameters.
#'
#' @return \code{input_list} with \code{$data}, \code{$par} and \code{$map}
#'   updated: the parsed integer arrays for \code{cont_tv_fish_sel},
#'   \code{fish_sel_blocks}, \code{fish_sel_model} and \code{fish_q_blocks}, the
#'   starting values for all four parameter groups, and their factor maps.
#'
#' @export Setup_Mod_Fishsel_and_Q
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Fishsel_and_Q <- function(input_list,

                                    # Total Selectivity
                                    cont_tv_fish_sel = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fish_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fish_sel_model,
                                    Use_fish_q_prior = 0,
                                    fish_q_prior = NA,
                                    fish_q_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fish_q_type = rep("est", input_list$data$n_fish_fleets),
                                    fish_q_model = NULL,
                                    sigma_fish_q_spec = "est_all",
                                    fish_q_rho_spec = "est_all",
                                    fish_q_rw_init_sigma = NA,
                                    fishsel_pe_pars_spec = NULL,
                                    fish_fixed_sel_pars_spec = NULL,
                                    fish_q_spec = NULL,
                                    fish_sel_devs_spec = NULL,
                                    corr_opt_semipar = NULL,
                                    Use_fish_selex_prior = 0,
                                    fish_selex_prior = NULL,
                                    Use_fish_selex_penalty = 0,
                                    fish_sel_norm_bins = NULL,
                                    fish_sel_bin_dev_bins = NULL,
                                    fishsel_pe_wt = rep(1, input_list$data$n_fish_fleets),
                                    fishsel_rw_init_sigma = rep(5, input_list$data$n_fish_fleets),
                                    fishsel_dont_est_dev_first = rep(0, input_list$data$n_fish_fleets),
                                    cont_tv_fishsel_bin_devs = rep("none", input_list$data$n_fish_fleets),
                                    fish_selex_penalty = NULL,
                                    fishsel_devs_shared_bins = NULL,
                                    fish_selex_type = 'age',
                                    use_fixed_fish_sel = rep(0, input_list$data$n_fish_fleets),
                                    fish_sel_input = NULL,
                                    fish_sel_nonpar_est_bins = NULL,
                                    fish_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
                                    fish_sel_dbnrml_raw = NULL,
                                    fish_sel_dbnrml_startbin = NULL,

                                    # Retained Selectivity
                                    cont_tv_ret_sel = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    ret_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    ret_sel_model = paste("logist1_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    retsel_pe_pars_spec = NULL,
                                    ret_fixed_sel_pars_spec = rep("fix_ret_sel_input", input_list$data$n_fish_fleets),
                                    ret_sel_devs_spec = NULL,
                                    ret_sel_corr_opt_semipar = NULL,
                                    Use_ret_selex_prior = 0,
                                    ret_selex_prior = NULL,
                                    retsel_devs_shared_bins = NULL,
                                    retsel_pe_wt = rep(1, input_list$data$n_fish_fleets),
                                    retsel_rw_init_sigma = rep(5, input_list$data$n_fish_fleets),
                                    retsel_dont_est_dev_first = rep(0, input_list$data$n_fish_fleets),
                                    ret_selex_type = 'age',
                                    use_fixed_ret_sel = rep(1, input_list$data$n_fish_fleets),
                                    ret_sel_input = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                    ret_sel_nonpar_est_bins = NULL,
                                    ret_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
                                    ...
                                    ) {

  messages_list <<- character(0) # string to attach to for printing messages # nolint: object_usage_linter.
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Fishsel_and_Q <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Continuous Selectivity Deviations
  check_fleet_spec_length(fishsel_pe_pars_spec, input_list$data$n_fish_fleets, "fishsel_pe_pars_spec", allow_null = TRUE)
  check_fleet_spec_length(fish_sel_devs_spec, input_list$data$n_fish_fleets, "fish_sel_devs_spec", allow_null = TRUE)
  check_fleet_spec_length(corr_opt_semipar, input_list$data$n_fish_fleets, "corr_opt_semipar", allow_null = TRUE)

  # A short vector here is read per fleet in the objective, so a length mismatch
  # silently becomes NA rather than being recycled.
  check_fleet_spec_length(fishsel_pe_wt, input_list$data$n_fish_fleets, "fishsel_pe_wt")
  check_fleet_spec_length(fishsel_rw_init_sigma, input_list$data$n_fish_fleets, "fishsel_rw_init_sigma")
  check_fleet_spec_length(fishsel_dont_est_dev_first, input_list$data$n_fish_fleets, "fishsel_dont_est_dev_first")
  if(!all(fishsel_dont_est_dev_first %in% c(0, 1))) stop("fishsel_dont_est_dev_first must be 0 or 1 for every fleet")

  # Catchability Priors
  if(!Use_fish_q_prior %in% c(0,1)) stop("Values for Use_fish_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking catchability priors
  if(Use_fish_q_prior == 1) {
    required_cols <- c("region", "fleet", "block", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(fish_q_prior))
    if(length(missing_cols) > 0) {
      stop("fish_q_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Fishery Catchability priors are: ", ifelse(Use_fish_q_prior == 0, "Not Used", "Used"))

  # fishery catchability type, mirroring srv_q_type. "est" estimates ln_fish_q; "arith" and "geo"
  # concentrate it out, and a concentrated fleet's ln_fish_q is never read, so map it off as well
  if(!all(fish_q_type %in% c("est", "arith", "geo"))) stop("Invalid specification for fish_q_type. Should be one of est, arith, or geo.")
  check_fleet_spec_length(fish_q_type, input_list$data$n_fish_fleets, "fish_q_type")
  fish_q_type_val <- match(fish_q_type, c("est", "arith", "geo")) - 1
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste0("Fishery Catchability for fishery ", f, " is: ", fish_q_type[f]))

  # Selectivity Priors
  if(!Use_fish_selex_prior %in% c(0,1)) stop("Values for Use_fish_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_fish_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(fish_selex_prior))
    if(length(missing_cols) > 0) {
      stop("fish_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Fishery Selectivity priors are: ", ifelse(Use_fish_selex_prior == 0, "Not Used", "Used"))

  # Fixed selectivity options
  if(any(use_fixed_fish_sel == 1) && is.null(fish_sel_input)) stop("fish_sel_input is NULL, please provide an input array.")
  if(any(use_fixed_fish_sel == 1) && fish_selex_type == 'age') check_data_dimensions(
    fish_sel_input,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_ages = length(input_list$data$ages),
    n_sexes = input_list$data$n_sexes,
    n_fish_fleets = input_list$data$n_fish_fleets,
    what = 'fish_sel_input_age'
  )
  if(any(use_fixed_fish_sel == 1) && fish_selex_type == 'length') check_data_dimensions(
    fish_sel_input,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_lens = length(input_list$data$lens),
    n_sexes = input_list$data$n_sexes,
    n_fish_fleets = input_list$data$n_fish_fleets,
    what = 'fish_sel_input_len'
  )

  # Selectivity Options -----------------------------------------------------
  # The bin vector is kept as well as its length. Starting values stated on the
  # bin scale are seeded further down, by which point fish_selex_type holds the
  # numeric code rather than the name it arrived as.
  if(fish_selex_type == 'age') {
    fish_selex_type <- 0
    fish_sel_bin_vec <- input_list$data$ages
    collect_message("Total Fishery Selectivity is aged-based.")
  } else if(fish_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but total selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    fish_selex_type <- 1
    fish_sel_bin_vec <- input_list$data$lens
    collect_message("Total Fishery Selectivity is length-based")
  } else stop("fish_selex_type must be 'age' or 'length', but was: ", fish_selex_type)

  bins <- length(fish_sel_bin_vec)

  # Continuous Time-Varying Selectivity Options -----------------------------
  cont_tv_fish_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in seq_along(cont_tv_fish_sel)) {
    # Extract out components from list
    tmp <- cont_tv_fish_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for cont_tv_fish_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_fish_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_fish_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous fishery time-varying selectivity specified as: ", cont_tv_type, " for fishery fleet ", fleet)
  }

  if(any(cont_tv_fish_sel_mat > 0) && (is.null(fishsel_pe_pars_spec) || is.null(fish_sel_devs_spec))) stop("Continuous time-varying selectivity specified, but fishsel_pe_pars_spec and/or fish_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  fish_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in seq_along(fish_sel_blocks)) {

    # Extract out components from list
    tmp <- fish_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_sel_blocks_arr[,,fleet] <- 1 # input only 1 fishery time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # extract fleet index

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fish_sel_blocks_arr[,years,fleet] <- block_val
    }

  }

  if(any(is.na(fish_sel_blocks_arr))) stop("Fishery Selectivtiy Blocks are returning an NA. Did you forget to specify the year range of fish_sel_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Selectivity Time Blocks for fishery", f, "is specified at:", length(unique(fish_sel_blocks_arr[,,f]))))

  # Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml", 'nonpar', 'asymplogist1', "asymplogist2", "bicubic", "nonparlog", "nonparfree"), num = c(0,1,2,3,4,5,6,7,8,9,10)) # set up values we can map to
  fish_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  fish_sel_bicubic_binnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bin nodes, only set where fish_sel_model == 8 (bicubic)
  fish_sel_bicubic_yrnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of year nodes, only set where fish_sel_model == 8 (bicubic)
  fish_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year, i.e. no offset); years within the block before this are edge-kept at this year's fitted curve
  fish_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bins (starting from the first) the bicubic surface is actually fit over (0 = all bins, i.e. no truncation); bins beyond this are kept flat at the last fitted bin's value

  for(i in seq_along(fish_sel_model)) {

    # Extract out fishery selectivity components from vector
    tmp_sel_form <- fish_sel_model[i]
    tmp_sel_form_vec <- unlist(strsplit(tmp_sel_form, "_")) # split string
    sel_form <- tmp_sel_form_vec[1] # get selectivity type

    if(sel_form == "bicubic") {
      # bicubic spline: bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f>[_Block_<b>][_SelStyr_<year>][_NSelBins_<n>]
      bin_pos <- which(tmp_sel_form_vec == "Bin")
      yr_pos <- which(tmp_sel_form_vec == "Yr")
      fleet_pos <- which(tmp_sel_form_vec == "Fleet")
      block_pos <- which(tmp_sel_form_vec == "Block")
      selstyr_pos <- which(tmp_sel_form_vec == "SelStyr")
      nselbins_pos <- which(tmp_sel_form_vec == "NSelBins")
      if(length(bin_pos) != 1 || length(yr_pos) != 1 || length(fleet_pos) != 1)
        stop("fish_sel_model 'bicubic' entries must be specified as bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f> or bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Block_<b>_Fleet_<f>, optionally with _SelStyr_<year> and/or _NSelBins_<n>")
      tmp_n_bin_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[bin_pos + 1]))
      tmp_n_yr_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[yr_pos + 1]))
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_selstyr <- if(length(selstyr_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[selstyr_pos + 1])) else 0
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(is.na(tmp_n_bin_nodes) || tmp_n_bin_nodes < 2) stop("bicubic fish_sel_model requires at least 2 bin nodes (n_bin_nodes >= 2)")
      if(is.na(tmp_n_yr_nodes) || tmp_n_yr_nodes < 1) stop("bicubic fish_sel_model requires at least 1 year node (n_yr_nodes >= 1). Use n_yr_nodes == 1 for a time-invariant bin-only spline.")
      if(length(selstyr_pos) == 1 && (is.na(tmp_selstyr) || !tmp_selstyr %in% input_list$data$years)) stop("bicubic fish_sel_model SelStyr must be a calendar year within the modeled years")
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("bicubic fish_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
    } else {
      # optional _NSelBins_<n> plateau suffix is joined with  _Block_<b> in either order for any parametric form
      fleet_pos <- which(tmp_sel_form_vec == "Fleet")
      block_pos <- which(tmp_sel_form_vec == "Block")
      nselbins_pos <- which(tmp_sel_form_vec == "NSelBins")
      if(length(fleet_pos) != 1) stop("fish_sel_model entries must name their fleet exactly once, as _Fleet_<f>")
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("fish_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
      if(length(nselbins_pos) == 1 && sel_form %in% c("nonpar", "nonparlog", "nonparfree")) stop("fish_sel_model NSelBins is for the parametric forms; a non-parametric form would keep estimating the bins it then overwrites. Group those bins through fish_sel_nonpar_est_bins instead.")
    }

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("fish_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml, nonpar, nonparlog, nonparfree, asymplogist1, asymplogist2, bicubic (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for fish_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) {
      fish_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
      if(sel_form == "bicubic") {
        fish_sel_bicubic_binnodes_arr[,,tmp_fleet] <- tmp_n_bin_nodes
        fish_sel_bicubic_yrnodes_arr[,,tmp_fleet] <- tmp_n_yr_nodes
        fish_sel_bicubic_selstyr_arr[,,tmp_fleet] <- tmp_selstyr
      }
      fish_sel_bicubic_nselbins_arr[,,tmp_fleet] <- tmp_nselbins # plateau bin; read by every form (0 = none)
    } else {
      fish_sel_model_arr <- assign_sel_block(fish_sel_model_arr, fish_sel_blocks_arr, tmp_fleet, tmp_block, sel_map$num[which(sel_map$sel == sel_form)])
      if(sel_form == "bicubic") {
        fish_sel_bicubic_binnodes_arr <- assign_sel_block(fish_sel_bicubic_binnodes_arr, fish_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_bin_nodes)
        fish_sel_bicubic_yrnodes_arr <- assign_sel_block(fish_sel_bicubic_yrnodes_arr, fish_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_yr_nodes)
        fish_sel_bicubic_selstyr_arr <- assign_sel_block(fish_sel_bicubic_selstyr_arr, fish_sel_blocks_arr, tmp_fleet, tmp_block, tmp_selstyr)
      }
      fish_sel_bicubic_nselbins_arr <- assign_sel_block(fish_sel_bicubic_nselbins_arr, fish_sel_blocks_arr, tmp_fleet, tmp_block, tmp_nselbins)
    }
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Fishery selectivity functional form specified as:", sel_form, " for fishery fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_fish_fleets) {
    has_blocks <- length(unique(fish_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_fish_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Blocked Catchability Options --------------------------------------------
  fish_q_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in seq_along(fish_q_blocks)) {
    # Extract out components from list
    tmp <- fish_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_q_blocks_arr[,,fleet] <- 1 # input only 1 fishery catchability time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # get fleet number

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fish_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }

  if(any(is.na(fish_q_blocks_arr))) stop("Fishery Catchability Blocks are returning an NA. Did you forget to specify the year range of fish_q_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Catchability Time Blocks for fishery", f, "is specified at:", length(unique(fish_q_blocks_arr[,,f]))))

  # Populate Data List ------------------------------------------------------

  input_list$data$cont_tv_fish_sel <- cont_tv_fish_sel_mat
  input_list$data$fish_sel_blocks <- fish_sel_blocks_arr
  input_list$data$fish_sel_model <- fish_sel_model_arr
  input_list$data$fish_sel_bicubic_binnodes <- fish_sel_bicubic_binnodes_arr
  input_list$data$fish_sel_bicubic_yrnodes <- fish_sel_bicubic_yrnodes_arr
  input_list$data$fish_sel_bicubic_selstyr <- fish_sel_bicubic_selstyr_arr
  input_list$data$fish_sel_bicubic_nselbins <- fish_sel_bicubic_nselbins_arr
  input_list$data$fish_q_type <- fish_q_type_val
  input_list$data$fish_q_blocks <- fish_q_blocks_arr
  input_list$data$fish_q_prior <- fish_q_prior
  input_list$data$Use_fish_q_prior <- Use_fish_q_prior
  input_list$data$Use_fish_selex_prior <- Use_fish_selex_prior
  input_list$data$fish_selex_prior <- validate_selex_prior_types(fish_selex_prior, Use_fish_selex_prior, "fish_selex_prior",
                                                                 sel_blocks = fish_sel_blocks_arr, n_bins = bins)
  input_list$data$Use_fish_selex_penalty <- Use_fish_selex_penalty
  input_list$data$fishsel_pe_wt <- fishsel_pe_wt
  input_list$data$fishsel_rw_init_sigma <- fishsel_rw_init_sigma
  input_list <- setup_sel_bin_devs(
    input_list,
    fish_sel_bin_dev_bins,
    cont_tv_fishsel_bin_devs,
    prefix = "fish",
    n_fleets = input_list$data$n_fish_fleets,
    bins = bins,
    starting_values = starting_values
  )
  input_list <- setup_sel_norm_bins(
    input_list,
    fish_sel_norm_bins,
    prefix = "fish",
    n_fleets = input_list$data$n_fish_fleets,
    bins = bins
  )
  input_list$data$fish_selex_penalty <- validate_selex_penalty(fish_selex_penalty, Use_fish_selex_penalty, "fish_selex_penalty")
  input_list$data$fish_selex_type <- fish_selex_type
  if(!is.null(input_list$data$fish_len_comp_sel) && any(input_list$data$fish_len_comp_sel == 1) && fish_selex_type != 1) stop("FishLenComps_sel = 'length' in Setup_Mod_FishIdx_and_Comps applies the length selectivity at length, so fish_selex_type must be length")
  if(!is.null(input_list$data$fish_waa_selected) && any(input_list$data$fish_waa_selected == 1) && fish_selex_type != 1) stop("fish_waa_selected = 1 in Setup_Mod_FishIdx_and_Comps weights the fishery weight at age by length selectivity, so fish_selex_type must be length")
  input_list$data$use_fixed_fish_sel <- use_fixed_fish_sel
  input_list$data$fish_sel_input <- fish_sel_input
  input_list$data$fishsel_devs_min_shared_bins <- if(!is.null(fishsel_devs_shared_bins)) unlist(lapply(fishsel_devs_shared_bins, min)) else seq_along(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------

  # Figure out number of selectivity parameters for a given functional form
  unique_fishsel_vals <- unique(as.vector(input_list$data$fish_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in seq_along(unique_fishsel_vals)) {
    if(unique_fishsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_fishsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_fishsel_vals[i] %in% c(4)) sel_pars_vec[i] <- 6 # double normal
    if(unique_fishsel_vals[i] %in% c(5,9,10)) sel_pars_vec[i] <- bins # non-parametric selex
    if(unique_fishsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic selex w/ asymptote
    if(unique_fishsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$fish_sel_bicubic_binnodes * input_list$data$fish_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  # figure out maximum number of fishery selectivity blocks for a given reigon and fleet
  max_fishsel_blks <- max(apply(input_list$data$fish_sel_blocks, c(1,3), FUN = function(x) length(unique(x))))

  # bicubic spline interpolation weight matrices (bin node by year node grid), passed through the
  # model with the flattened node parameters. zero-weight rows and columns never contribute
  has_bicubic_fish_sel <- any(input_list$data$fish_sel_model == 8)
  max_bin_nodes_bicubic <- if(has_bicubic_fish_sel) max(input_list$data$fish_sel_bicubic_binnodes) else 1
  max_yr_nodes_bicubic <- if(has_bicubic_fish_sel) max(input_list$data$fish_sel_bicubic_yrnodes) else 1
  n_yrs_total_bicubic <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs

  fish_sel_bicubic_Wbin <- array(0, dim = c(input_list$data$n_regions, bins, max_bin_nodes_bicubic, max_fishsel_blks, input_list$data$n_fish_fleets))
  fish_sel_bicubic_Wyr <- array(0, dim = c(input_list$data$n_regions, n_yrs_total_bicubic, max_yr_nodes_bicubic, max_fishsel_blks, input_list$data$n_fish_fleets))

  if(has_bicubic_fish_sel) {
    for(f in 1:input_list$data$n_fish_fleets) {
      for(r in 1:input_list$data$n_regions) {

        fishsel_blocks_tmp <- unique(as.vector(input_list$data$fish_sel_blocks[r,,f]))

        for(b in seq_along(fishsel_blocks_tmp)) {

          block_years <- which(input_list$data$fish_sel_blocks[r,,f] == fishsel_blocks_tmp[b])
          if(unique(input_list$data$fish_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$fish_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$fish_sel_bicubic_yrnodes[r, block_years, f])

          # bin nodes evenly spaced over [0,1]. the spline is fit over all bins unless NSelBins is
          # set, beyond which bins are edge-kept at the last fitted bin's weights
          nselbins_this <- unique(input_list$data$fish_sel_bicubic_nselbins[r, block_years, f])
          n_fit_bins <- if(nselbins_this == 0) bins else nselbins_this

          bin_nodes_scaled <- seq(0, 1, length.out = n_bin_nodes_this)
          fit_bin_scaled <- seq(0, 1, length.out = n_fit_bins)
          Wbin_fit <- Get_Natural_Cubic_Spline_Weights(bin_nodes_scaled, fit_bin_scaled)

          Wbin_this <- matrix(0, nrow = bins, ncol = n_bin_nodes_this)
          Wbin_this[1:n_fit_bins, ] <- Wbin_fit
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(
            Wbin_fit[nrow(Wbin_fit), ],
            nrow = bins - n_fit_bins,
            ncol = n_bin_nodes_this,
            byrow = TRUE
          )

          fish_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # year nodes are evenly spaced over the block's fit range, which is the whole block unless
          # SelStyr is set. years outside that range hold the boundary node weights constant
          selstyr_this <- unique(input_list$data$fish_sel_bicubic_selstyr[r, block_years, f])
          selstyr_idx <- if(selstyr_this == 0) min(block_years) else which(input_list$data$years == selstyr_this)
          fit_years <- block_years[block_years >= selstyr_idx]
          pre_fit_years <- block_years[block_years < selstyr_idx]

          yr_nodes_scaled <- seq(0, 1, length.out = n_yr_nodes_this)
          fit_yr_scaled <- seq(0, 1, length.out = length(fit_years))
          Wyr_block <- Get_Natural_Cubic_Spline_Weights(yr_nodes_scaled, fit_yr_scaled)

          Wyr_this <- matrix(0, nrow = n_yrs_total_bicubic, ncol = n_yr_nodes_this)
          Wyr_this[fit_years, ] <- Wyr_block
          if(length(pre_fit_years) > 0) Wyr_this[pre_fit_years, ] <- matrix(Wyr_block[1, ], nrow = length(pre_fit_years), ncol = n_yr_nodes_this, byrow = TRUE)
          if(min(block_years) > 1) Wyr_this[1:(min(block_years) - 1), ] <- matrix(Wyr_block[1, ], nrow = min(block_years) - 1, ncol = n_yr_nodes_this, byrow = TRUE)
          if(max(block_years) < n_yrs_total_bicubic) Wyr_this[(max(block_years) + 1):n_yrs_total_bicubic, ] <- matrix(
            Wyr_block[nrow(Wyr_block), ],
            nrow = n_yrs_total_bicubic - max(block_years),
            ncol = n_yr_nodes_this,
            byrow = TRUE
          )

          fish_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this

        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_fish_sel

  input_list$data$fish_sel_bicubic_Wbin <- fish_sel_bicubic_Wbin
  input_list$data$fish_sel_bicubic_Wyr <- fish_sel_bicubic_Wyr

  # maximum number of selectivity parameters across all forms
  max_fishsel_pars <- max(sel_pars_vec)
  input_list$par$fish_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_fishsel_pars, max_fishsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))
  input_list$par$fish_fixed_sel_pars <- use_starting_value(input_list$par$fish_fixed_sel_pars, starting_values, "fish_fixed_sel_pars")

  # a double normal has its peak on the bin scale, so a default of zero
  # would put it at bin zero, where the ascending limb has no extent
  if(!"fish_fixed_sel_pars" %in% names(starting_values)) {
    input_list$par$fish_fixed_sel_pars <- seed_dbnrml_peak(input_list$par$fish_fixed_sel_pars, fish_sel_model_arr,
                                          fish_sel_bin_vec,
                                          fish_sel_sex_offset)
  }

  # Fishery catchability
  max_fishq_blks <- max(apply(input_list$data$fish_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of fishery catchability blocks for a given reigon and fleet
  input_list$par$ln_fish_q <- array(0, dim = c(input_list$data$n_regions, max_fishq_blks, input_list$data$n_fish_fleets))
  input_list$par$ln_fish_q <- use_starting_value(input_list$par$ln_fish_q, starting_values, "ln_fish_q")

  # Fishery selectivity process error parameters
  input_list$par$fishsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_fishsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using
  input_list$par$fishsel_pe_pars <- use_starting_value(input_list$par$fishsel_pe_pars, starting_values, "fishsel_pe_pars")

  # Fishery selectivity deviations
  input_list$par$ln_fishsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_fish_fleets))
  input_list$par$ln_fishsel_devs <- use_starting_value(input_list$par$ln_fishsel_devs, starting_values, "ln_fishsel_devs")

  # Sex offsets on selectivity (parameter offsets and/or a curve scale offset)
  input_list <- setup_sel_sex_offset(
    input_list,
    fish_sel_sex_offset,
    prefix = "fish",
    n_fleets = input_list$data$n_fish_fleets,
    fleet_label = "fishery fleet",
    sel_model_arr = input_list$data$fish_sel_model,
    cont_tv_mat = cont_tv_fish_sel_mat,
    max_blks = max_fishsel_blks,
    sel_blocks = input_list$data$fish_sel_blocks,
    fixed_spec = fish_fixed_sel_pars_spec,
    starting_values = starting_values
  )

  # Raw (unanchored) double normal limbs; retention keeps anchored limbs
  input_list$data$fish_dbnrml_raw <- setup_dbnrml_raw(fish_sel_dbnrml_raw, input_list$data$n_fish_fleets, "fish_sel_dbnrml_raw")
  input_list$data$ret_dbnrml_raw <- array(0, dim = c(input_list$data$n_fish_fleets, 2))
  input_list$data$fish_dbnrml_startbin <- setup_dbnrml_startbin(fish_sel_dbnrml_startbin, input_list$data$n_fish_fleets, bins, "fish_sel_dbnrml_startbin")
  input_list$data$ret_dbnrml_startbin <- rep(1, input_list$data$n_fish_fleets)

  # Mapping Options ---------------------------------------------------------
  input_list <- do_fixed_sel_pars_mapping(
    input_list,
    fish_fixed_sel_pars_spec,
    bins,
    fish_sel_nonpar_est_bins,
    prefix = "fish",
    fleet_field = "n_fish_fleets",
    use_field = "Catch",
    fleet_label = "fishery fleet"
  )
  input_list <- do_q_mapping(
    input_list,
    fish_q_spec,
    prefix = "fish",
    fleet_field = "n_fish_fleets",
    fleet_label = "fishery fleet"
  )
  input_list <- setup_q_devs(
    input_list,
    q_model = fish_q_model,
    sigma_q_spec = sigma_fish_q_spec,
    q_rho_spec = fish_q_rho_spec,
    q_rw_init_sigma = fish_q_rw_init_sigma,
    q_type = fish_q_type,
    prefix = "fish",
    fleet_field = "n_fish_fleets",
    use_field = "UseFishIdx",
    fleet_label = "fishery fleet",
    starting_values = starting_values
  )
  input_list <- do_sel_pe_pars_mapping(
    input_list,
    fishsel_pe_pars_spec,
    corr_opt_semipar,
    bins,
    sel_devs_spec = fish_sel_devs_spec,
    sel_devs_shared_bins = fishsel_devs_shared_bins,
    prefix = "fish",
    fleet_field = "n_fish_fleets",
    use_field = "Catch",
    fleet_label = "fishery fleet"
  )
  input_list <- do_sel_devs_mapping(
    input_list,
    fish_sel_devs_spec,
    fishsel_devs_shared_bins,
    bins,
    dont_est_dev_first = fishsel_dont_est_dev_first,
    prefix = "fish",
    fleet_field = "n_fish_fleets",
    use_field = "Catch",
    fleet_label = "fishery fleet"
  )


  # Retained Selectivity ---------------------------------------------
  input_list <- Setup_Mod_Retsel(input_list,
                                 cont_tv_ret_sel = cont_tv_ret_sel,
                                 ret_sel_blocks = ret_sel_blocks,
                                 ret_sel_model = ret_sel_model,
                                 retsel_pe_pars_spec = retsel_pe_pars_spec,
                                 ret_fixed_sel_pars_spec = ret_fixed_sel_pars_spec,
                                 ret_sel_devs_spec = ret_sel_devs_spec,
                                 ret_sel_corr_opt_semipar = ret_sel_corr_opt_semipar,
                                 Use_ret_selex_prior = Use_ret_selex_prior,
                                 ret_selex_prior = ret_selex_prior,
                                 retsel_devs_shared_bins = retsel_devs_shared_bins,
                                 retsel_pe_wt = retsel_pe_wt,
                                 retsel_rw_init_sigma = retsel_rw_init_sigma,
                                 retsel_dont_est_dev_first = retsel_dont_est_dev_first,
                                 ret_selex_type = ret_selex_type,
                                 use_fixed_ret_sel = use_fixed_ret_sel,
                                 ret_sel_input = ret_sel_input,
                                 ret_sel_nonpar_est_bins = ret_sel_nonpar_est_bins,
                                 ret_sel_sex_offset = ret_sel_sex_offset,
                                 ...)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
