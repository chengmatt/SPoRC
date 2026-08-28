# Stage 1 of 3: model setup
#
# Fishery selectivity, retention and catchability inputs.
# Setup_Mod_Fishsel_and_Q covers the selectivity of the gear on the stock;
# Setup_Mod_Retsel covers what fraction of the encountered fish are kept.
# Both delegate their parameter maps to the shared builders in setup_mapping.R.

#' Set up retained fishery selectivity
#'
#' Configures all aspects of retained fishery selectivity and catchability for the
#' estimation model, including functional forms, time blocks, continuous time-varying
#' selectivity, process error structure, annual deviations, and fixed or estimated
#' selectivity parameters. Must be called after
#' \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' Selectivity time-variation via \code{cont_tv_ret_sel} and blocked selectivity via
#' \code{ret_sel_blocks} are mutually exclusive within a fleet. Specifying both for
#' the same fleet will result in an error.
#'
#' @param input_list Named list containing \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists as created by upstream setup functions.
#'
#' @param cont_tv_ret_sel Character vector defining continuous time-varying selectivity
#'   structure per fleet. Format:
#'   \code{"<type>_Fleet_<f>"} where \code{type} is one of:
#'   \describe{
#'     \item{\code{"none"}}{No time variation}
#'     \item{\code{"iid"}}{Independent year effects}
#'     \item{\code{"rw"}}{Random walk process}
#'     \item{\code{"3dmarg"}}{3D marginal structure}
#'     \item{\code{"3dcond"}}{3D conditional structure}
#'     \item{\code{"2dar1"}}{2D AR1 structure}
#'   }
#'   Output is mapped to integer codes internally.
#'
#' @param ret_sel_blocks Character vector defining discrete selectivity blocks.
#'   Format:
#'   \code{"none_Fleet_<f>"} or
#'   \code{"Block_<b>_Year_<start>-<end>_Fleet_<f>"} (or terminal year version).
#'   Values are expanded into a 3D array:
#'   \code{[n_regions × n_years × n_fish_fleets]}.
#'
#' @param ret_sel_model Character vector defining selectivity functional forms.
#'   Format:
#'   \code{"<type>_Fleet_<f>"} or
#'   \code{"<type>_Fleet_<f>_Block_<b>"}.
#'   Supported types:
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'     \item{\code{"nonpar"}}{Non-parametric over discrete age or length bins, on
#'       the logit scale, then mean-standardized jointly over years and bins so the
#'       grand mean of the surface is one. Bins may be grouped through the
#'       non-parametric bin mapping. No fixed functional form is imposed.}
#'     \item{\code{"nonparlog"}}{Non-parametric on the log scale, standardized so
#'       each year's selectivity averages to one over \code{*_sel_norm_bins}. Only
#'       within-year contrasts are identified; the level is absorbed by
#'       catchability or fishing mortality.}
#'     \item{\code{"nonparfree"}}{Non-parametric on the log scale with no
#'       standardization, \eqn{\exp(\theta)}, so the values carry the height of the
#'       curve as well as its shape. This is the form for a stream fit age by age:
#'       a free catchability per age and a selectivity estimated at age are one
#'       quantity written two ways, so the whole age multiplier lives here and no
#'       catchability is set. Pin one bin, by leaving it out of the estimated bins,
#'       whenever the mean it multiplies is also free.}
#'     \item{\code{"asymplogist1"}}{Logistic selectivity with \eqn{a_{50}} and slope \eqn{k} and asymptotic control (3 parameters).}
#'     \item{\code{"asymplogist2"}}{Logistic selectivity with with \eqn{a_{50}} and \eqn{a_{95}} and asymptotic control (3 parameters).}
#'     \item{\code{"bicubic"}}{Bicubic spline over a bin-node x year-node grid; see \code{ret_sel_model} in \code{\link{Setup_Mod_Fishsel_and_Q}} for full syntax.}
#'   }
#'
#' @param retsel_pe_pars_spec Specification of process error parameters for
#'   time-varying selectivity. Length must equal \code{n_fish_fleets}.
#'
#' @param ret_fixed_sel_pars_spec Specification controlling which fixed
#'   selectivity parameters are estimated vs fixed.
#'
#' @param ret_sel_devs_spec Specification of annual selectivity deviations
#'   structure per fleet.
#'
#' @param ret_sel_corr_opt_semipar Optional correlation structure for semi-parametric
#'   selectivity deviations. Length must equal \code{n_fish_fleets}.
#'
#' @param Use_ret_selex_prior Integer flag (\code{0/1}) indicating whether
#'   selectivity priors are used.
#'
#' @param ret_selex_prior Data frame of priors for selectivity parameters.
#'   Must include columns: \code{region}, \code{fleet}, \code{block}, \code{sex},
#'   \code{par}, \code{mu}, \code{sd}, plus an optional \code{type}
#'   (\code{"par"}/\code{"value"}; see \code{fish_selex_prior} in
#'   \code{Setup_Mod_Fishsel_and_Q}).
#'
#' @param retsel_devs_shared_bins Vector defining shared bins for selectivity
#'   deviation estimation (e.g., age or length grouping structure).
#'
#' @param ret_sel_bin_dev_bins List with one element per fishery fleet naming the
#'   bins that fleet overrides, or \code{NULL} for fleets with no overrides. An
#'   overridden bin takes a freely estimated annual value in place of whatever
#'   the functional form produced, applied after every other transformation
#'   including standardization. Default \code{NULL}.
#'
#' @param cont_tv_retsel_bin_devs Character vector of length
#'   \code{n_fish_fleets} giving the process error on the bin-override
#'   deviations for each fleet: \code{"none"} (default), \code{"iid"}, or
#'   \code{"rw"}.
#'
#' @param retsel_pe_wt Numeric vector of length \code{n_fish_fleets}. Per-fleet
#'   multiplier on the retention selectivity process error likelihood. Default
#'   \code{1} for every fleet. \code{0} skips that fleet's process error
#'   likelihood altogether, so the deviations stay estimated but enter the
#'   objective only through the data and any explicit smoothness or centering
#'   penalties. Values other than 0 or 1 make an estimated process error sigma
#'   reinterpretable, so prefer 0 or 1 unless deliberately down-weighting.
#'   Applies only to \code{ln_retsel_devs}; the bin-override deviations carry
#'   their own process error and are not affected.
#'
#' @param retsel_rw_init_sigma Numeric vector of length \code{n_fish_fleets}.
#'   Standard deviation given to the first year of an \code{"rw"} deviation
#'   series. Default \code{5}, which leaves that year effectively free.
#'   \code{NA} instead starts the walk at zero under the walk's own estimated
#'   sigma, making the first year as smooth as every later step. Appropriate
#'   when the base parametric curve already describes the first year well.
#'
#' @param ret_selex_type Character string indicating selectivity domain:
#'   \code{"age"} or \code{"length"}.
#'
#' @param use_fixed_ret_sel Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = use fixed selectivity input; \code{0} = estimate.
#'
#' @param ret_sel_input Fixed selectivity input array:
#'   \code{[n_pop × n_regions × n_years × n_seas × (ages or lengths) × n_sexes × n_fish_fleets]}.
#'
#' @param ret_sel_nonpar_est_bins Vector defining estimated bins for
#'   non-parametric selectivity.
#'
#' @param ret_sel_sex_offset Character vector of length \code{n_fish_fleets}
#'   linking the sexes of a fleet's retention curve: \code{"none"} (default),
#'   \code{"par"}, \code{"scale"}, or \code{"par_scale"}, with the same
#'   meaning as \code{fish_sel_sex_offset} in
#'   \code{\link{Setup_Mod_Fishsel_and_Q}}. Retention is a fraction, so a
#'   scale offset is only sensible where the scaled curve stays at or below
#'   one (a negative offset, the less-retained sex); nothing enforces that.
#'
#' @param ... Optional starting values for selectivity parameters and deviations.
#'
#' @return The updated \code{input_list} with:
#' \itemize{
#'   \item Parsed selectivity structure arrays in \code{$data}
#'   \item Starting parameter arrays in \code{$par}
#'   \item Factor mapping objects in \code{$map}
#' }
#'
#' @keywords internal
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Retsel <- function(input_list,
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
                             ret_sel_nonpar_est_bins,
                             ret_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
                             ...
) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Fishsel_and_Q <- c(input_list$config$Setup_Mod_Fishsel_and_Q, mget(names(formals()))[-1])


  # Input Validation --------------------------------------------------------
  # Continuous Selectivity Deviations
  if(!is.null(retsel_pe_pars_spec)) if(length(retsel_pe_pars_spec) != input_list$data$n_fish_fleets) stop("retsel_pe_pars_spec is not length n_fish_fleets")
  if(!is.null(ret_sel_devs_spec)) if(length(ret_sel_devs_spec) != input_list$data$n_fish_fleets) stop("ret_sel_devs_spec is not length n_fish_fleets")
  if(!is.null(ret_sel_corr_opt_semipar)) if(length(ret_sel_corr_opt_semipar) != input_list$data$n_fish_fleets) stop("ret_sel_corr_opt_semipar is not length n_fish_fleets")

  # A short vector here is read per fleet in the objective, so a length mismatch
  # silently becomes NA rather than being recycled.
  if(length(retsel_pe_wt) != input_list$data$n_fish_fleets) stop("retsel_pe_wt is not length n_fish_fleets")
  if(length(retsel_rw_init_sigma) != input_list$data$n_fish_fleets) stop("retsel_rw_init_sigma is not length n_fish_fleets")

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
  if(any(use_fixed_ret_sel == 1) && ret_selex_type == 'age') check_data_dimensions(ret_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ret_sel_input_age')
  if(any(use_fixed_ret_sel == 1) && ret_selex_type == 'length') check_data_dimensions(ret_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ret_sel_input_len')

  # Selectivity Options -----------------------------------------------------
  # Age based selectivity
  if(ret_selex_type == 'age') {
    ret_selex_type <- 0
    bins <- length(input_list$data$ages)
    collect_message("Retained Fishery Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(ret_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but retained selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    ret_selex_type <- 1
    bins <- length(input_list$data$lens)
    collect_message("Retained Fishery Selectivity is length-based")
  } # if length based


  # Continuous Retained Time-Varying Selectivity Options -----------------------------
  cont_tv_ret_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_ret_sel)) {
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

  if(any(cont_tv_ret_sel_mat > 0) && is.null(retsel_pe_pars_spec) && is.null(ret_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but retsel_pe_pars_spec and/or ret_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Retained Time-Varying Selectivity Options --------------------------------
  ret_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(ret_sel_blocks)) {

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
  ret_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year); years before this are edge-held, matching fish_sel_model's SelStyr
  ret_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bins the bicubic surface is actually fit over (0 = all bins); bins beyond this are edge-held, matching fish_sel_model's NSelBins

  for(i in 1:length(ret_sel_model)) {

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
  input_list <- setup_sel_bin_devs(input_list, ret_sel_bin_dev_bins, cont_tv_retsel_bin_devs,
                                   prefix = "ret", n_fleets = input_list$data$n_fish_fleets,
                                   bins = bins, starting_values = starting_values)
  input_list$data$retsel_devs_min_shared_bins <- if(!is.null(retsel_devs_shared_bins)) unlist(lapply(retsel_devs_shared_bins, min)) else 1:length(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------

  # Figure out number of selectivity parameters for a given functional form
  unique_retsel_vals <- unique(as.vector(input_list$data$ret_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_retsel_vals)) {
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
  if("ret_fixed_sel_pars" %in% names(starting_values)) input_list$par$ret_fixed_sel_pars <- starting_values$ret_fixed_sel_pars
  else input_list$par$ret_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_retsel_pars, max_retsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # a double normal carries its peak on the bin scale, so a default of zero
  # would put it at bin zero, where the ascending limb has no extent
  if(!"ret_fixed_sel_pars" %in% names(starting_values)) {
    input_list$par$ret_fixed_sel_pars <- seed_dbnrml_peak(input_list$par$ret_fixed_sel_pars, ret_sel_model_arr,
                                          if(ret_selex_type == 'length') input_list$data$lens else input_list$data$ages,
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

        for(b in 1:length(retsel_blocks_tmp)) {

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
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(Wbin_fit[nrow(Wbin_fit), ], nrow = bins - n_fit_bins, ncol = n_bin_nodes_this, byrow = TRUE)

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
          if(max(block_years) < n_yrs_total_bicubic_ret) Wyr_this[(max(block_years) + 1):n_yrs_total_bicubic_ret, ] <- matrix(Wyr_block[nrow(Wyr_block), ], nrow = n_yrs_total_bicubic_ret - max(block_years), ncol = n_yr_nodes_this, byrow = TRUE)

          ret_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this
        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_ret_sel

  input_list$data$ret_sel_bicubic_Wbin <- ret_sel_bicubic_Wbin
  input_list$data$ret_sel_bicubic_Wyr <- ret_sel_bicubic_Wyr

  # Retained Fishery selectivity process error parameters
  if("retsel_pe_pars" %in% names(starting_values)) input_list$par$retsel_pe_pars <- starting_values$retsel_pe_pars
  else input_list$par$retsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_retsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Retained Fishery selectivity deviations
  if(input_list$data$ret_selex_type == 0) bins <- length(input_list$data$ages) # age based deviations
  if(input_list$data$ret_selex_type == 1) bins <- length(input_list$data$lens) # length based deviations
  if("ln_retsel_devs" %in% names(starting_values)) input_list$par$ln_retsel_devs <- starting_values$ln_retsel_devs
  else input_list$par$ln_retsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Sex offsets on retention selectivity (parameter offsets and/or a curve scale offset)
  input_list <- setup_sel_sex_offset(input_list, ret_sel_sex_offset, prefix = "ret",
                                     n_fleets = input_list$data$n_fish_fleets, fleet_label = "fishery fleet (retention)",
                                     sel_model_arr = input_list$data$ret_sel_model, cont_tv_mat = cont_tv_ret_sel_mat,
                                     max_blks = max_retsel_blks, sel_blocks = input_list$data$ret_sel_blocks,
                                     fixed_spec = ret_fixed_sel_pars_spec, starting_values = starting_values)

  # Mapping Options ---------------------------------------------------------
  input_list <- do_fixed_sel_pars_mapping(input_list, ret_fixed_sel_pars_spec, bins, ret_sel_nonpar_est_bins,
                                          prefix = "ret", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet")
  input_list <- do_sel_pe_pars_mapping(input_list, retsel_pe_pars_spec, ret_sel_corr_opt_semipar, bins,
                                       prefix = "ret", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet")
  input_list <- do_sel_devs_mapping(input_list, ret_sel_devs_spec, retsel_devs_shared_bins, bins,
                                    prefix = "ret", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet")

  return(input_list)
}

#' Set up total and retained fishery selectivity and catchability specifications
#'
#' Configures all aspects of fishery selectivity and catchability for the
#' estimation model: functional forms, time blocks, continuous time-varying
#' structures, process error hyperparameters, annual deviations, and
#' catchability blocks and estimation structure. Must be called after
#' \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' Selectivity time-variation and blocked selectivity are mutually exclusive
#' within a fleet. Specifying both for the same fleet will raise an error.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists.
#' @param fish_sel_model Character vector specifying the selectivity functional
#'   form for each fleet (and optionally each time block). Each element must
#'   follow one of:
#'   \itemize{
#'     \item \code{"<model>_Fleet_<f>"}: single form for all years of fleet \code{f}.
#'     \item \code{"<model>_Fleet_<f>_Block_<b>"}: form specific to block \code{b}
#'       of fleet \code{f}, as defined in \code{fish_sel_blocks}.
#'   }
#'   Available models:
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'     \item{\code{"nonpar"}}{Non-parametric over discrete age or length bins, on
#'       the logit scale, then mean-standardized jointly over years and bins so the
#'       grand mean of the surface is one. Bins may be grouped through the
#'       non-parametric bin mapping. No fixed functional form is imposed.}
#'     \item{\code{"nonparlog"}}{Non-parametric on the log scale, standardized so
#'       each year's selectivity averages to one over \code{*_sel_norm_bins}. Only
#'       within-year contrasts are identified; the level is absorbed by
#'       catchability or fishing mortality.}
#'     \item{\code{"nonparfree"}}{Non-parametric on the log scale with no
#'       standardization, \eqn{\exp(\theta)}, so the values carry the height of the
#'       curve as well as its shape. This is the form for a stream fit age by age:
#'       a free catchability per age and a selectivity estimated at age are one
#'       quantity written two ways, so the whole age multiplier lives here and no
#'       catchability is set. Pin one bin, by leaving it out of the estimated bins,
#'       whenever the mean it multiplies is also free.}
#'     \item{\code{"asymplogist1"}}{Logistic selectivity with \eqn{a_{50}} and slope \eqn{k} and asymptotic control (3 parameters).}
#'     \item{\code{"asymplogist2"}}{Logistic selectivity with with \eqn{a_{50}} and \eqn{a_{95}} and asymptotic control (3 parameters).}
#'     \item{\code{"bicubic"}}{Bicubic spline over a bin-node x year-node grid
#'       (see \code{\link{Get_Selex}}, \code{Selex_Model == 8}). Specified as
#'       \code{"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"} (optionally
#'       with \code{_Block_k}). One generalized form covers a smooth bin x year
#'       surface (\code{n_yr_nodes > 1}), a time-invariant bin-only spline
#'       (\code{n_yr_nodes == 1}), or a bin-only spline re-fit independently
#'       per year-block (\code{n_yr_nodes == 1} within each of several blocks
#'       defined via \code{fish_sel_blocks}). An optional \code{_SelStyr_<year>}
#'       suffix (a calendar year within the block) restricts the actual spline
#'       fit to \code{SelStyr}:block-end; years within the block before
#'       \code{SelStyr} are held constant at the \code{SelStyr} year's fitted
#'       curve, rather than fitting the surface over the whole block. An
#'       optional \code{_NSelBins_<n>} suffix restricts the actual spline fit
#'       to the first \code{n} bins (ages or lengths, per \code{fish_selex_type});
#'       bins beyond \code{n} are held constant at the last fitted bin's
#'       curve.}
#'   }
#'   See the model equations vignette for mathematical definitions.
#' @param cont_tv_fish_sel Character vector of length \code{n_fish_fleets}
#'   specifying continuous time-varying selectivity per fleet. Each element
#'   must be \code{"<type>_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"none"}}{No continuous time-variation (default).}
#'     \item{\code{"iid"}}{IID annual deviations on selectivity parameters.}
#'     \item{\code{"rw"}}{Random walk in selectivity parameters over time.}
#'     \item{\code{"3dmarg"}}{3D GMRF with marginal variance parameterization.}
#'     \item{\code{"3dcond"}}{3D GMRF with conditional variance parameterization.}
#'     \item{\code{"2dar1"}}{2D separable AR1 in bin and year dimensions.}
#'   }
#'   If any fleet has \code{cont_tv_fish_sel != "none"}, both
#'   \code{fishsel_pe_pars_spec} and \code{fish_sel_devs_spec} must also be
#'   provided.
#' @param fish_sel_blocks Character vector defining discrete selectivity time
#'   blocks per fleet. Each element follows \code{"Block_<b>_Year_<s>-<e>_Fleet_<f>"}
#'   or \code{"Block_<b>_Year_<s>-terminal_Fleet_<f>"}. Use
#'   \code{"none_Fleet_<f>"} (default) for a single constant block. Blocks must
#'   be non-overlapping and together span all model years for the specified fleet.
#'   Mutually exclusive with \code{cont_tv_fish_sel != "none"} for the same fleet.
#' @param fish_q_blocks Character vector defining catchability time blocks per
#'   fleet, using the same format as \code{fish_sel_blocks}. Default
#'   \code{"none_Fleet_<f>"} gives a single constant block.
#' @param fish_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying how fixed-effect selectivity parameters are estimated. See
#'   \code{\link{do_fixed_sel_pars_mapping}} for all options
#'   (\code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"est_shared_f_x"}, \code{"fix"}).
#' @param fish_q_spec Character vector of length \code{n_fish_fleets} specifying
#'   catchability estimation structure. See \code{\link{do_q_mapping}} for
#'   options (\code{"est_all"}, \code{"est_shared_r"}, \code{"fix"}).
#' @param fishsel_pe_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for selectivity process error
#'   hyperparameters. Required when any fleet has continuous time-variation.
#'   See \code{\link{do_sel_pe_pars_mapping}} for all options.
#' @param fish_sel_devs_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for annual selectivity deviations.
#'   Required when any fleet has continuous time-variation. See
#'   \code{\link{do_sel_devs_mapping}} for all options including age-sharing
#'   options for semi-parametric forms.
#' @param Use_fish_selex_penalty Integer (0/1). Whether a centering penalty is
#'   applied to sets of fishery selectivity fixed-effect parameters. Default
#'   \code{0}.
#' @param fish_selex_penalty Data frame of centering penalty specifications,
#'   required when \code{Use_fish_selex_penalty = 1}. Required columns:
#'   \code{region}, \code{fleet}, \code{block}, \code{sex}, \code{par}, and
#'   \code{wt}. Each row penalizes \code{wt * (log(mean(exp(pars))))^2} over the
#'   set of parameters named in \code{par}, which may be a single index or a
#'   list column of integer vectors naming a whole set. This pins the scalar of
#'   a non-parametric curve that catchability or fishing mortality would
#'   otherwise absorb, and is softer than fixing a bin outright. Intended for
#'   parameter sets held on the log scale. Default \code{NULL}.
#' @param fishsel_devs_shared_bins List of integer vectors grouping age or length
#'   bins that share a single deviation series. Only used when
#'   \code{fish_sel_devs_spec} contains one of the \code{"est_shared_b"} variants.
#'   Example: \code{list(1:5, 6:10, 11:30)}.
#' @param corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   controlling which correlation components to suppress in semi-parametric
#'   (3D GMRF or 2D AR1) time-varying selectivity. Set to \code{NA} (default)
#'   for no suppression. See \code{\link{do_sel_pe_pars_mapping}} for valid
#'   suppression codes. Cohort-correlation options are invalid for \code{"2dar1"}.
#' @param fish_q_type Character vector of length \code{n_fish_fleets} controlling
#'   how catchability is obtained. \code{"est"} (default) estimates
#'   \code{ln_fish_q}. \code{"arith"} concentrates it out of the likelihood as
#'   the ratio of mean observed to mean predicted index, and \code{"geo"} does
#'   the same on the log scale as \code{exp(mean(log(obs) - log(pred)))}. Both
#'   analytic forms solve one catchability per region and fleet using only the
#'   years with observations, ignore any block structure, and fix that fleet's
#'   \code{ln_fish_q} regardless of \code{fish_q_spec}.
#' @param fish_q_formula Named list of one-sided formulas, one element per fleet
#'   requiring catchability covariates, referencing series in
#'   \code{fish_q_cov_dat}. \code{NULL} (default) excludes covariate effects.
#' @param fish_q_cov_dat Named list of numeric vectors (length =
#'   \code{n_years}) containing the covariate time series referenced in
#'   \code{fish_q_formula}. All vectors must be the same length and contain no
#'   missing values; set values to \code{0} for years when the fishery index is
#'   not active. Default \code{NULL}.
#' @param fishsel_pe_wt Numeric vector of length \code{n_fish_fleets}. Per-fleet
#'   multiplier on the fishery selectivity process error likelihood. Default
#'   \code{1} for every fleet. \code{0} skips that fleet's process error
#'   likelihood altogether, so the deviations stay estimated but enter the
#'   objective only through the data and any explicit smoothness or centering
#'   penalties, which is how several existing assessments constrain them. Values
#'   other than 0 or 1 make an estimated process error sigma reinterpretable, so
#'   prefer 0 or 1 unless deliberately down-weighting. Applies only to
#'   \code{ln_fishsel_devs}; the bin-override deviations carry their own process
#'   error and are not affected.
#' @param fishsel_rw_init_sigma Numeric vector of length \code{n_fish_fleets}.
#'   Standard deviation given to the first year of an \code{"rw"} deviation
#'   series. Default \code{5}, which leaves that year effectively free.
#'   \code{NA} instead starts the walk at zero under the walk's own estimated
#'   sigma, making the first year as smooth as every later step. Appropriate when
#'   the base parametric curve already describes the first year well.
#' @param fish_sel_norm_bins List with one element per fishery fleet naming the
#'   bins the mean-one standardization averages over, or \code{NULL} for fleets
#'   standardizing over every bin (\code{Selex_Model = 9} only). A gear whose
#'   catchability is defined against part of the bin range standardizes over
#'   that part, and catchability absorbs the difference in scale. Default
#'   \code{NULL}.
#' @param fish_sel_bin_dev_bins List with one element per fishery fleet naming
#'   the bins that fleet overrides, or \code{NULL} for fleets with no overrides
#'   (e.g. \code{list(1, NULL)} frees bin 1 of fleet 1 only). An overridden bin
#'   takes a freely estimated annual value \eqn{\exp(\epsilon_{y,b})} in place of
#'   whatever the functional form produced, applied after every other
#'   transformation including standardization. The rest of the curve keeps its
#'   parametric shape. Default \code{NULL}.
#' @param cont_tv_fishsel_bin_devs Character vector of length
#'   \code{n_fish_fleets} giving the process error on the bin-override
#'   deviations for each fleet: \code{"none"} (default), \code{"iid"}, or
#'   \code{"rw"}. A random walk carries its own estimated sigma per bin, with
#'   \code{fishsel_bin_devs_rw_init_sigma} governing its first year.
#' @param Use_fish_q_prior Integer flag. \code{1} = apply lognormal priors to
#'   catchability; \code{0} = no priors (default). Requires \code{fish_q_prior}.
#' @param fish_q_prior Data frame of catchability prior hyperparameters. Required
#'   columns: \code{region}, \code{fleet}, \code{block} (block index), \code{mu}
#'   (prior mean on natural scale), \code{sd} (prior SD on log scale). Each row
#'   specifies a \eqn{\text{Normal}(\log(\mu), \sigma)} prior for one catchability
#'   parameter. Only used when \code{Use_fish_q_prior = 1}.
#' @param Use_fish_selex_prior Integer flag. \code{1} = apply lognormal priors to
#'   selectivity parameters; \code{0} = no priors (default). Requires
#'   \code{fish_selex_prior}.
#' @param fish_selex_prior Data frame of selectivity prior hyperparameters, one
#'   row per prior. Required columns: \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par} (parameter index within the
#'   functional form), \code{mu}, \code{sd}, plus an optional \code{type}
#'   giving each row's target: \code{"par"} (the default when the column is
#'   absent) is a lognormal prior on one fixed selectivity parameter, with
#'   \code{mu} on the natural scale and \code{sd} on the log scale;
#'   \code{"value"} is a normal prior on the realized selectivity value at one
#'   bin, with both on the natural scale, where \code{par} instead names the
#'   bin (on ages or lengths per \code{fish_selex_type}) and the value is read
#'   at the first model year of \code{block}. A \code{"value"} row constrains
#'   the derived selectivity value rather than the parameters, matching the
#'   ADMB convention of pinning selectivity at a reference age near one, which
#'   no set of independent parameter priors can express. Only used when
#'   \code{Use_fish_selex_prior = 1}.
#' @param ... Optional starting value overrides for selectivity parameters.
#' @param cont_tv_ret_sel Character vector of length \code{n_fish_fleets}
#'   specifying continuous time-varying selectivity per fleet. Each element
#'   must be \code{"<type>_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"none"}}{No continuous time-variation (default).}
#'     \item{\code{"iid"}}{IID annual deviations on selectivity parameters.}
#'     \item{\code{"rw"}}{Random walk in selectivity parameters over time.}
#'     \item{\code{"3dmarg"}}{3D GMRF with marginal variance parameterization.}
#'     \item{\code{"3dcond"}}{3D GMRF with conditional variance parameterization.}
#'     \item{\code{"2dar1"}}{2D separable AR1 in bin and year dimensions.}
#'   }
#'   If any fleet has \code{cont_tv_ret_sel != "none"}, both
#'   \code{retsel_pe_pars_spec} and \code{ret_sel_devs_spec} must also be
#'   provided.
#'
#' @param ret_sel_blocks Character vector defining discrete selectivity time
#'   blocks per fleet. Each element follows \code{"Block_<b>_Year_<s>-<e>_Fleet_<f>"}
#'   or \code{"Block_<b>_Year_<s>-terminal_Fleet_<f>"}. Use
#'   \code{"none_Fleet_<f>"} (default) for a single constant block. Blocks must
#'   be non-overlapping and together span all model years for the specified fleet.
#'   Mutually exclusive with \code{cont_tv_ret_sel != "none"} for the same fleet.
#'
#' @param ret_sel_model Character vector specifying the selectivity functional
#'   form for each fleet (and optionally each time block). Each element must
#'   follow one of:
#'   \itemize{
#'     \item \code{"<model>_Fleet_<f>"}: single form for all years of fleet \code{f}.
#'     \item \code{"<model>_Fleet_<f>_Block_<b>"}: form specific to block \code{b}
#'       of fleet \code{f}, as defined in \code{ret_sel_blocks}.
#'   }
#'   Available models:
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'     \item{\code{"nonpar"}}{Non-parametric over discrete age or length bins, on
#'       the logit scale, then mean-standardized jointly over years and bins so the
#'       grand mean of the surface is one. Bins may be grouped through the
#'       non-parametric bin mapping. No fixed functional form is imposed.}
#'     \item{\code{"nonparlog"}}{Non-parametric on the log scale, standardized so
#'       each year's selectivity averages to one over \code{*_sel_norm_bins}. Only
#'       within-year contrasts are identified; the level is absorbed by
#'       catchability or fishing mortality.}
#'     \item{\code{"nonparfree"}}{Non-parametric on the log scale with no
#'       standardization, \eqn{\exp(\theta)}, so the values carry the height of the
#'       curve as well as its shape. This is the form for a stream fit age by age:
#'       a free catchability per age and a selectivity estimated at age are one
#'       quantity written two ways, so the whole age multiplier lives here and no
#'       catchability is set. Pin one bin, by leaving it out of the estimated bins,
#'       whenever the mean it multiplies is also free.}
#'     \item{\code{"asymplogist1"}}{Logistic selectivity with \eqn{a_{50}} and slope \eqn{k} and asymptotic control (3 parameters).}
#'     \item{\code{"asymplogist2"}}{Logistic selectivity with with \eqn{a_{50}} and \eqn{a_{95}} and asymptotic control (3 parameters).}
#'     \item{\code{"bicubic"}}{Bicubic spline over a bin-node x year-node grid, specified as
#'       \code{"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"} (optionally with \code{_Block_k},
#'       \code{_SelStyr_<year>}, and/or \code{_NSelBins_<n>}); see \code{fish_sel_model} above for the
#'       full syntax and \code{\link{Get_Selex}} (\code{Selex_Model == 8}) for the underlying math.}
#'   }
#'   See the model equations vignette for mathematical definitions.
#'
#' @param retsel_pe_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for selectivity process error
#'   hyperparameters. Required when any fleet has continuous time-variation.
#'   See \code{\link{do_sel_pe_pars_mapping}} for all options.
#'
#' @param ret_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying how fixed-effect selectivity parameters are estimated. See
#'   \code{\link{do_fixed_sel_pars_mapping}} for all options
#'   (\code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"est_shared_f_x"}, \code{"fix"}).
#'
#' @param ret_sel_devs_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for annual selectivity deviations.
#'   Required when any fleet has continuous time-variation. See
#'   \code{\link{do_sel_devs_mapping}} for all options including age-sharing
#'   options for semi-parametric forms.
#'
#' @param retsel_devs_shared_bins List of integer vectors grouping age or length
#'   bins that share a single deviation series. Only used when
#'   \code{ret_sel_devs_spec} contains one of the \code{"est_shared_b"} variants.
#'   Example: \code{list(1:5, 6:10, 11:30)}.
#'
#' @param retsel_pe_wt Numeric vector of length \code{n_fish_fleets}. Per-fleet
#'   multiplier on the retention selectivity process error likelihood, the
#'   retention counterpart of \code{fishsel_pe_wt}. Default \code{1} for every
#'   fleet, and \code{0} skips that fleet's process error likelihood so its
#'   deviations stay estimated but are constrained only by the data and any
#'   explicit smoothness or centering penalties.
#'
#' @param retsel_rw_init_sigma Numeric vector of length \code{n_fish_fleets}.
#'   Standard deviation given to the first year of an \code{"rw"} retention
#'   deviation series, the retention counterpart of
#'   \code{fishsel_rw_init_sigma}. Default \code{5}; \code{NA} instead starts
#'   the walk at zero under the walk's own estimated sigma.
#'
#' @param ret_sel_corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   controlling which correlation components to suppress in semi-parametric
#'   (3D GMRF or 2D AR1) time-varying selectivity. Set to \code{NA} (default)
#'   for no suppression. See \code{\link{do_sel_pe_pars_mapping}} for valid
#'   suppression codes. Cohort-correlation options are invalid for \code{"2dar1"}.
#'
#' @param Use_ret_selex_prior Integer flag. \code{1} = apply lognormal priors to
#'   selectivity parameters; \code{0} = no priors (default). Requires
#'   \code{ret_selex_prior}.
#'
#' @param ret_selex_prior Data frame of selectivity prior hyperparameters.
#'   Required columns: \code{region}, \code{fleet}, \code{block}, \code{sex},
#'   \code{par}, \code{mu}, \code{sd}, plus an optional \code{type}
#'   (\code{"par"}/\code{"value"}; see \code{fish_selex_prior}).
#'
#' @param use_fixed_ret_sel Integer vector of length \code{n_fish_fleets}
#'   indicating whether to fix selectivity (\code{1}) or estimate it (\code{0}).
#'
#' @param ret_sel_input Array of fixed selectivity values with dimensions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]}.
#'
#' @param ret_sel_nonpar_est_bins Optional list specifying bin groupings for
#' non-parametric retained selectivity. Structure is \code{[[fleet]][[block]]}, where each
#' element is a list of bin index vectors defining grouped parameters.
#' @param ret_sel_sex_offset Character vector of length \code{n_fish_fleets}
#'   linking the sexes of a fleet's retention curve, with the options and
#'   meaning of \code{fish_sel_sex_offset}. Default \code{"none"}. Retention
#'   is a fraction, so a scale offset is only sensible where the scaled curve
#'   stays at or below one.
#' @param fish_selex_type Character scalar specifying whether selectivity is
#'   age- or length-based. Options:
#'   \describe{
#'     \item{\code{"age"}}{Selectivity is defined over age bins.}
#'     \item{\code{"length"}}{Selectivity is defined over length bins.}
#'   }
#'   Determines the bin dimension used for all fishery selectivity functions,
#'   including parametric, time-varying, and non-parametric forms.

#' @param use_fixed_fish_sel Integer vector of length \code{n_fish_fleets}
#'   indicating whether fishery selectivity is fixed (\code{1}) or estimated
#'   (\code{0}) for each fleet.

#' @param fish_sel_input Array of fixed fishery selectivity values with
#'   dimensions:
#'   \code{[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]}.
#'   Required when any element of \code{use_fixed_fish_sel == 1}.

#' @param fish_sel_nonpar_est_bins Optional list defining bin groupings for
#' non-parametric fishery selectivity. Structure is \code{[[fleet]][[block]]}, where each
#' element is a list of integer vectors. Each vector defines a group of bins that
#' share a single estimated selectivity parameter. Indices must correspond to the
#' bin dimension defined by \code{fish_selex_type}.

#' @param fish_sel_dbnrml_startbin \code{NULL} (default) or an integer vector
#'   \code{[n_fish_fleets]}, the bin each fleet's double normal anchors its
#'   ascending limb at (\code{1} is the first bin). Bins below it take the
#'   squared ratio of their bin to it times the selectivity there, Stock
#'   Synthesis's convention when the compositions start above the population's
#'   first length bin.
#' @param fish_sel_dbnrml_raw \code{NULL} (default) or a 0/1 matrix
#'   \code{[n_fish_fleets x 2]} for fleets on the double normal: column one
#'   leaves the ascending limb as a raw Gaussian instead of anchoring it to
#'   \code{p5} at the first bin, column two does the same for the descending
#'   limb and \code{p6}.
#' @param fish_sel_sex_offset Character vector of length \code{n_fish_fleets}
#'   linking the sexes of a fleet's selectivity, for models with
#'   \code{n_sexes > 1}. Options per fleet:
#'   \describe{
#'     \item{\code{"none"} (default)}{Each sex's stored parameters are its own,
#'       exactly as before this option existed.}
#'     \item{\code{"par"}}{The stored fixed-effect parameter slots of every sex
#'       beyond the first hold additive offsets on the first sex's stored
#'       (transformed-scale) parameters, so for log-scale parameters the
#'       sex-\eqn{s} natural value is the first sex's times
#'       \eqn{e^{\delta}}. Offsets fixed at zero reproduce sex-shared
#'       parameters; estimating them links the sexes through the offset the way
#'       several existing assessments parameterize male selectivity.}
#'     \item{\code{"scale"}}{Each sex keeps its own parameters, and every sex
#'       beyond the first additionally carries a constant log-scale offset on
#'       the whole realized curve, \code{exp(ln_fishsel_sex_scale)}, estimated
#'       per region, block, and sex. The scaled curve may exceed one. Refused
#'       for non-parametric forms and semi-parametric time variation, whose
#'       post-hoc standardization would cancel a constant multiplier.}
#'     \item{\code{"apical"}}{Each sex keeps its own parameters, and for every
#'       sex beyond the first the double normal builds its limbs up to
#'       \code{exp(ln_*sel_sex_scale)} rather than to one. Selectivity at the
#'       first and last bins stays where that sex's own parameters put it, so
#'       the offset moves the middle of the curve and leaves its ends anchored.
#'       Requires the double normal.}
#'     \item{\code{"par_apical"}}{Both a par offset and an apical offset.}
#'     \item{\code{"par_scale"}}{Both a par offset and a scale offset.}
#'   }

#' @param ret_selex_type Character scalar specifying whether retained selectivity
#'   is age- or length-based. Options:
#'   \describe{
#'     \item{\code{"age"}}{Selectivity is defined over age bins.}
#'     \item{\code{"length"}}{Selectivity is defined over length bins.}
#'   }
#'   Determines the bin dimension used for all retained selectivity functions,
#'   including parametric, time-varying, and non-parametric forms.
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated. Key additions include the parsed integer arrays for
#'   \code{cont_tv_fish_sel}, \code{fish_sel_blocks}, \code{fish_sel_model}, and
#'   \code{fish_q_blocks}; starting value arrays for all four parameter groups;
#'   and their corresponding factor maps.
#'
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
                                    fish_q_cov_dat = NULL,
                                    fish_q_formula = NULL,
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
                                    ret_selex_type = 'age',
                                    use_fixed_ret_sel = rep(1, input_list$data$n_fish_fleets),
                                    ret_sel_input = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets )),
                                    ret_sel_nonpar_est_bins = NULL,
                                    ret_sel_sex_offset = rep("none", input_list$data$n_fish_fleets),
                                    ...
                                    ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Fishsel_and_Q <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Continuous Selectivity Deviations
  if(!is.null(fishsel_pe_pars_spec)) if(length(fishsel_pe_pars_spec) != input_list$data$n_fish_fleets) stop("fishsel_pe_pars_spec is not length n_fish_fleets")
  if(!is.null(fish_sel_devs_spec)) if(length(fish_sel_devs_spec) != input_list$data$n_fish_fleets) stop("fish_sel_devs_spec is not length n_fish_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_fish_fleets) stop("corr_opt_semipar is not length n_fish_fleets")

  # A short vector here is read per fleet in the objective, so a length mismatch
  # silently becomes NA rather than being recycled.
  if(length(fishsel_pe_wt) != input_list$data$n_fish_fleets) stop("fishsel_pe_wt is not length n_fish_fleets")
  if(length(fishsel_rw_init_sigma) != input_list$data$n_fish_fleets) stop("fishsel_rw_init_sigma is not length n_fish_fleets")

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

  # Fishery catchability type, mirroring srv_q_type. "est" estimates ln_fish_q;
  # "arith" and "geo" concentrate it out of the likelihood, solving it from the
  # observed and predicted index over the years with observations. A concentrated
  # fleet's ln_fish_q is never read, so it should also be mapped off.
  if(!all(fish_q_type %in% c("est", "arith", "geo"))) stop("Invalid specification for fish_q_type. Should be one of est, arith, or geo.")
  if(length(fish_q_type) != input_list$data$n_fish_fleets) stop("fish_q_type is not length n_fish_fleets")
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
  if(any(use_fixed_fish_sel == 1) && fish_selex_type == 'age') check_data_dimensions(fish_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'fish_sel_input_age')
  if(any(use_fixed_fish_sel == 1) && fish_selex_type == 'length') check_data_dimensions(fish_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'fish_sel_input_len')

  # Selectivity Options -----------------------------------------------------
  # Age based selectivity
  if(fish_selex_type == 'age') {
    fish_selex_type <- 0
    bins <- length(input_list$data$ages)
    collect_message("Total Fishery Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(fish_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but total selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    fish_selex_type <- 1
    bins <- length(input_list$data$lens)
    collect_message("Total Fishery Selectivity is length-based")
  } # if length based

  # Continuous Time-Varying Selectivity Options -----------------------------
  cont_tv_fish_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_fish_sel)) {
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

  if(any(cont_tv_fish_sel_mat > 0) && is.null(fishsel_pe_pars_spec) && is.null(fish_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but fishsel_pe_pars_spec and/or fish_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  fish_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_sel_blocks)) {

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
  fish_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year, i.e. no offset); years within the block before this are edge-held at this year's fitted curve
  fish_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bins (starting from the first) the bicubic surface is actually fit over (0 = all bins, i.e. no truncation); bins beyond this are held flat at the last fitted bin's value

  for(i in 1:length(fish_sel_model)) {

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
  for(i in 1:length(fish_q_blocks)) {
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

  if(any(is.na(fish_q_blocks))) stop("Fishery Catchability Blocks are returning an NA. Did you forget to specify the year range of fish_q_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Catchability Time Blocks for fishery", f, "is specified at:", length(unique(fish_q_blocks_arr[,,f]))))

  # Populate Data List ------------------------------------------------------

  input_list$data$cont_tv_fish_sel <- cont_tv_fish_sel_mat
  input_list$data$fish_sel_blocks <- fish_sel_blocks_arr
  input_list$data$fish_sel_model <- fish_sel_model_arr
  input_list$data$fish_sel_bicubic_binnodes <- fish_sel_bicubic_binnodes_arr
  input_list$data$fish_sel_bicubic_yrnodes <- fish_sel_bicubic_yrnodes_arr
  input_list$data$fish_sel_bicubic_selstyr <- fish_sel_bicubic_selstyr_arr
  input_list$data$fish_sel_bicubic_nselbins <- fish_sel_bicubic_nselbins_arr
  # Catchability covariate containers, mirroring the survey
  if(!is.null(fish_q_cov_dat) && !is.null(fish_q_formula)) {
    n_fish_q_cov <- max(sapply(names(fish_q_formula), function(key) {
      tmp_formula <- fish_q_formula[[key]]
      var_names <- all.vars(tmp_formula)
      tmp_dat <- data.frame(fish_q_cov_dat[var_names])
      ncol(stats::model.matrix(tmp_formula, data = tmp_dat))
    }))
  } else {
    do_fish_q_cov <- 0
    n_fish_q_cov <- 1 # dummy to initialize the array
  }

  fish_q_cov <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets, n_fish_q_cov)) # environmental time series
  fish_q_coeff <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets, n_fish_q_cov)) # coefficients to be estimated
  map_fish_q_coeff <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets, n_fish_q_cov)) # coefficients to be mapped off

  if(!is.null(fish_q_cov_dat) && !is.null(fish_q_formula)) {

    cov_lengths <- lengths(fish_q_cov_dat)
    if (length(unique(cov_lengths)) != 1) stop("All covariates in 'fish_q_cov_dat' must have the same length. If some years are missing data, either impute some value, or set at 0 (if it is not used in the calculation).")
    if (unique(cov_lengths) != length(input_list$data$years)) stop(paste0("Covariate length mismatch: expected ",  length(input_list$data$years),  " years but got ", unique(cov_lengths),  "."))

    do_fish_q_cov <- 1
    coeff_counter <- 0

    for(r in 1:input_list$data$n_regions) {
      for(f in 1:input_list$data$n_fish_fleets) {

        key <- paste(paste("Region", r, sep = "_"), "_Fleet_", f, sep = "")
        tmp_formula <- fish_q_formula[[key]]
        var_names <- all.vars(tmp_formula)
        if(length(var_names) == 0) next
        tmp_dat <- data.frame(fish_q_cov_dat[var_names])
        tmp_design_mat <- stats::model.matrix(tmp_formula, data = tmp_dat)
        fish_q_cov[r,,f,1:ncol(tmp_design_mat)] <- tmp_design_mat

        for(i in 1:ncol(tmp_design_mat)) {
          coeff_counter <- coeff_counter + 1
          map_fish_q_coeff[r,f,i] <- coeff_counter
        } # end i loop

      } # end f loop
    } # end r loop
  } # if using covariates

  input_list$data$fish_q_type <- fish_q_type_val
  input_list$data$do_fish_q_cov <- do_fish_q_cov
  input_list$data$fish_q_cov <- fish_q_cov
  input_list$par$fish_q_coeff <- fish_q_coeff
  input_list$map$fish_q_coeff <- factor(map_fish_q_coeff)
  input_list$data$fish_q_blocks <- fish_q_blocks_arr
  input_list$data$fish_q_prior <- fish_q_prior
  input_list$data$Use_fish_q_prior <- Use_fish_q_prior
  input_list$data$Use_fish_selex_prior <- Use_fish_selex_prior
  input_list$data$fish_selex_prior <- validate_selex_prior_types(fish_selex_prior, Use_fish_selex_prior, "fish_selex_prior",
                                                                 sel_blocks = fish_sel_blocks_arr, n_bins = bins)
  input_list$data$Use_fish_selex_penalty <- Use_fish_selex_penalty
  input_list$data$fishsel_pe_wt <- fishsel_pe_wt
  input_list$data$fishsel_rw_init_sigma <- fishsel_rw_init_sigma
  input_list <- setup_sel_bin_devs(input_list, fish_sel_bin_dev_bins, cont_tv_fishsel_bin_devs,
                                   prefix = "fish", n_fleets = input_list$data$n_fish_fleets,
                                   bins = bins, starting_values = starting_values)
  input_list <- setup_sel_norm_bins(input_list, fish_sel_norm_bins, prefix = "fish", n_fleets = input_list$data$n_fish_fleets, bins = bins)
  input_list$data$fish_selex_penalty <- validate_selex_penalty(fish_selex_penalty, Use_fish_selex_penalty, "fish_selex_penalty")
  input_list$data$fish_selex_type <- fish_selex_type
  if(!is.null(input_list$data$fish_len_comp_sel) && any(input_list$data$fish_len_comp_sel == 1) && fish_selex_type != 1) stop("FishLenComps_sel = 'length' in Setup_Mod_FishIdx_and_Comps applies the length selectivity at length, so fish_selex_type must be length")
  if(!is.null(input_list$data$fish_waa_selected) && any(input_list$data$fish_waa_selected == 1) && fish_selex_type != 1) stop("fish_waa_selected = 1 in Setup_Mod_FishIdx_and_Comps weights the fishery weight at age by length selectivity, so fish_selex_type must be length")
  input_list$data$use_fixed_fish_sel <- use_fixed_fish_sel
  input_list$data$fish_sel_input <- fish_sel_input
  input_list$data$fishsel_devs_min_shared_bins <- if(!is.null(fishsel_devs_shared_bins)) unlist(lapply(fishsel_devs_shared_bins, min)) else 1:length(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------

  # Figure out number of selectivity parameters for a given functional form
  unique_fishsel_vals <- unique(as.vector(input_list$data$fish_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_fishsel_vals)) {
    if(unique_fishsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_fishsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_fishsel_vals[i] %in% c(4)) sel_pars_vec[i] <- 6 # double normal
    if(unique_fishsel_vals[i] %in% c(5,9,10)) sel_pars_vec[i] <- bins # non-parametric selex
    if(unique_fishsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic selex w/ asymptote
    if(unique_fishsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$fish_sel_bicubic_binnodes * input_list$data$fish_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  # figure out maximum number of fishery selectivity blocks for a given reigon and fleet
  max_fishsel_blks <- max(apply(input_list$data$fish_sel_blocks, c(1,3), FUN = function(x) length(unique(x))))

  # Bicubic spline interpolation weight matrices (bin node x year node grid), built here so they can be
  # passed through the model  alongside the flattened node parameters (see Get_Selex, Selex_Model == 8).
  # Padded with zeros to a common width across regions/blocks/fleets; padding doesn't really do anything because unused
  # (zero-weight) columns/rows never contribute to the resulting selectivity (see Get_Selex documentation).
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

        for(b in 1:length(fishsel_blocks_tmp)) {

          block_years <- which(input_list$data$fish_sel_blocks[r,,f] == fishsel_blocks_tmp[b])
          if(unique(input_list$data$fish_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$fish_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$fish_sel_bicubic_yrnodes[r, block_years, f])

          # Age/length dimension: nodes evenly spaced over [0,1]. By default (NSelBins unset, i.e. 0)
          # the spline is evaluated over all bins, as before. When NSelBins is set , the spline surface is only actually fit over the first NSelBins bins;
          # bins beyond that are edge-held at the last fitted bin's weights ("plateau")
          nselbins_this <- unique(input_list$data$fish_sel_bicubic_nselbins[r, block_years, f])
          n_fit_bins <- if(nselbins_this == 0) bins else nselbins_this

          bin_nodes_scaled <- seq(0, 1, length.out = n_bin_nodes_this)
          fit_bin_scaled <- seq(0, 1, length.out = n_fit_bins)
          Wbin_fit <- Get_Natural_Cubic_Spline_Weights(bin_nodes_scaled, fit_bin_scaled)

          Wbin_this <- matrix(0, nrow = bins, ncol = n_bin_nodes_this)
          Wbin_this[1:n_fit_bins, ] <- Wbin_fit
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(Wbin_fit[nrow(Wbin_fit), ], nrow = bins - n_fit_bins, ncol = n_bin_nodes_this, byrow = TRUE)

          fish_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # Year dimension: nodes evenly spaced over the block's own contiguous fit range. By default
          # (SelStyr unset, i.e. 0) the fit range is the whole block, as before. When SelStyr is set, 
          # only years from SelStyr through the block's end are
          # actually spline-fit; years within the block before SelStyr are edge-held at the SelStyr
          # row's weights ("previous years are filled"). Years outside the block entirely (before it,
          # after it, and any projection years, since projections reuse the terminal modeled year's
          # block) hold the boundary node weights constant, which for a spline evaluated exactly at
          # its first/last node reduces to full weight on that node.
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
          if(max(block_years) < n_yrs_total_bicubic) Wyr_this[(max(block_years) + 1):n_yrs_total_bicubic, ] <- matrix(Wyr_block[nrow(Wyr_block), ], nrow = n_yrs_total_bicubic - max(block_years), ncol = n_yr_nodes_this, byrow = TRUE)

          fish_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this

        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_fish_sel

  input_list$data$fish_sel_bicubic_Wbin <- fish_sel_bicubic_Wbin
  input_list$data$fish_sel_bicubic_Wyr <- fish_sel_bicubic_Wyr

  # maximum number of selectivity parameters across all forms
  max_fishsel_pars <- max(sel_pars_vec)
  if("fish_fixed_sel_pars" %in% names(starting_values)) input_list$par$fish_fixed_sel_pars <- starting_values$fish_fixed_sel_pars
  else input_list$par$fish_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_fishsel_pars, max_fishsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # a double normal carries its peak on the bin scale, so a default of zero
  # would put it at bin zero, where the ascending limb has no extent
  if(!"fish_fixed_sel_pars" %in% names(starting_values)) {
    input_list$par$fish_fixed_sel_pars <- seed_dbnrml_peak(input_list$par$fish_fixed_sel_pars, fish_sel_model_arr,
                                          if(fish_selex_type == 'length') input_list$data$lens else input_list$data$ages,
                                          fish_sel_sex_offset)
  }

  # Fishery catchability
  max_fishq_blks <- max(apply(input_list$data$fish_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of fishery catchability blocks for a given reigon and fleet
  if("ln_fish_q" %in% names(starting_values)) input_list$par$ln_fish_q <- starting_values$ln_fish_q
  else input_list$par$ln_fish_q <- array(0, dim = c(input_list$data$n_regions, max_fishq_blks, input_list$data$n_fish_fleets))

  # Fishery selectivity process error parameters
  if("fishsel_pe_pars" %in% names(starting_values)) input_list$par$fishsel_pe_pars <- starting_values$fishsel_pe_pars
  else input_list$par$fishsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_fishsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Fishery selectivity deviations
  if("ln_fishsel_devs" %in% names(starting_values)) input_list$par$ln_fishsel_devs <- starting_values$ln_fishsel_devs
  else input_list$par$ln_fishsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Sex offsets on selectivity (parameter offsets and/or a curve scale offset)
  input_list <- setup_sel_sex_offset(input_list, fish_sel_sex_offset, prefix = "fish",
                                     n_fleets = input_list$data$n_fish_fleets, fleet_label = "fishery fleet",
                                     sel_model_arr = input_list$data$fish_sel_model, cont_tv_mat = cont_tv_fish_sel_mat,
                                     max_blks = max_fishsel_blks, sel_blocks = input_list$data$fish_sel_blocks,
                                     fixed_spec = fish_fixed_sel_pars_spec, starting_values = starting_values)

  # Raw (unanchored) double normal limbs; retention keeps anchored limbs
  input_list$data$fish_dbnrml_raw <- setup_dbnrml_raw(fish_sel_dbnrml_raw, input_list$data$n_fish_fleets, "fish_sel_dbnrml_raw")
  input_list$data$ret_dbnrml_raw <- array(0, dim = c(input_list$data$n_fish_fleets, 2))
  input_list$data$fish_dbnrml_startbin <- setup_dbnrml_startbin(fish_sel_dbnrml_startbin, input_list$data$n_fish_fleets, bins, "fish_sel_dbnrml_startbin")
  input_list$data$ret_dbnrml_startbin <- rep(1, input_list$data$n_fish_fleets)

  # Mapping Options ---------------------------------------------------------
  input_list <- do_fixed_sel_pars_mapping(input_list, fish_fixed_sel_pars_spec, bins, fish_sel_nonpar_est_bins,
                                          prefix = "fish", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet")
  input_list <- do_q_mapping(input_list, fish_q_spec, prefix = "fish", fleet_field = "n_fish_fleets", fleet_label = "fishery fleet")
  input_list <- do_sel_pe_pars_mapping(input_list, fishsel_pe_pars_spec, corr_opt_semipar, bins,
                                       prefix = "fish", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet")
  input_list <- do_sel_devs_mapping(input_list, fish_sel_devs_spec, fishsel_devs_shared_bins, bins,
                                    prefix = "fish", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet")


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
