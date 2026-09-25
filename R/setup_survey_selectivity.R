# Stage 1 of 3: model setup
#
# Survey selectivity and catchability inputs. Mirrors setup_fishery_selectivity.R but has no
# retention, and delegates its parameter maps to the shared builders in setup_mapping.R.

#' Set up survey selectivity and catchability specifications
#'
#' Sets the survey selectivity forms, time blocks, continuous time variation,
#' process error and deviations, the catchability blocks and estimation structure,
#' and the selectivity and catchability priors. Call after
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}. \code{$data$srv_selex_type} must already be set by
#'   \code{\link{Setup_Mod_Biologicals}}.
#' @param srv_sel_model Character vector of the selectivity form per fleet, and
#'   optionally per block: \code{"<model>_Fleet_x"} or
#'   \code{"<model>_Fleet_x_Block_k"}, the latter required when a fleet has several
#'   blocks. The forms and their syntax, including the \code{"bicubic"} suffixes,
#'   are those of \code{fish_sel_model} in \code{\link{Setup_Mod_Fishsel_and_Q}}.
#'   No default.
#' @param srv_fixed_sel_pars_spec Character vector \code{[n_srv_fleets]} of the
#'   sharing structure for the fixed-effect selectivity parameters. No default. See
#'   \code{\link{do_fixed_sel_pars_mapping}}.
#' @param srv_sel_blocks Character vector of discrete selectivity time blocks, each
#'   \code{"Block_k_Year_a-b_Fleet_x"} with \code{"terminal"} allowed as the end
#'   year, or \code{"none_Fleet_x"} (default) for constant selectivity. Parsed into
#'   an \code{[n_regions × n_years × n_srv_fleets]} array. Mutually exclusive with
#'   continuous time variation for the same fleet.
#' @param cont_tv_srv_sel Character vector of continuous time variation per fleet,
#'   each \code{"<type>_Fleet_x"}: \code{"none"} (default), \code{"iid"},
#'   \code{"rw"}, \code{"3dmarg"}, \code{"3dcond"} (3D GMRF over age, year and
#'   cohort) or \code{"2dar1"} (separable over bin and year). Any fleet other than
#'   \code{"none"} also needs \code{srvsel_pe_pars_spec} and
#'   \code{srv_sel_devs_spec}.
#' @param srvsel_pe_pars_spec Character vector \code{[n_srv_fleets]} or \code{NULL}
#'   (default) of the sharing structure for the process error hyperparameters. See
#'   \code{\link{do_sel_pe_pars_mapping}}.
#' @param srv_sel_devs_spec Character vector \code{[n_srv_fleets]} or \code{NULL}
#'   (default) of the sharing structure for the deviation series. See
#'   \code{\link{do_sel_devs_mapping}}.
#' @param Use_srv_selex_penalty Integer (0/1). Whether a centering penalty is
#'   applied to sets of survey selectivity fixed-effect parameters. Default
#'   \code{0}.
#' @param srv_selex_penalty Data frame of centering penalties with columns
#'   \code{region}, \code{fleet}, \code{block}, \code{sex}, \code{par} and
#'   \code{wt}, required when \code{Use_srv_selex_penalty = 1}. Each row penalizes
#'   \code{wt * (log(mean(exp(pars))))^2} over the parameters named in \code{par},
#'   a single index or a list column of integer vectors. This pins the scalar of a
#'   non-parametric curve that catchability would otherwise absorb, and is softer
#'   than fixing a bin. Meant for parameter sets on the log scale. Default
#'   \code{NULL}.
#' @param srvsel_devs_shared_bins List of integer vectors grouping the bins that
#'   share one deviation series, e.g. \code{list(1:5, 6:10, 11:30)}. Required when
#'   \code{srv_sel_devs_spec} names an \code{"est_shared_b"} variant. Default
#'   \code{NULL}.
#' @param corr_opt_semipar Character vector \code{[n_srv_fleets]} or \code{NULL}
#'   (default) of which correlation components to suppress under the 3D GMRF or 2D
#'   AR1 forms. See \code{\link{do_sel_pe_pars_mapping}}.
#' @param srvsel_pe_wt Numeric vector \code{[n_srv_fleets]} multiplying the survey
#'   selectivity process error likelihood. Default \code{1}. \code{0} skips that
#'   fleet's process error, so the deviations stay estimated but enter the objective
#'   only through the data and any smoothness or centering penalties. Values other
#'   than 0 or 1 make an estimated process error sigma reinterpretable. Applies to
#'   \code{ln_srvsel_devs} only; the bin-override deviations have their own process
#'   error.
#' @param srvsel_rw_init_sigma Numeric vector \code{[n_srv_fleets]} giving the
#'   standard deviation of the first year of an \code{"rw"} deviation series.
#'   Default \code{5}, which leaves that year effectively free. \code{NA} instead
#'   starts the walk at zero under the walk's own estimated sigma, which suits a
#'   base curve that already describes the first year well.
#' @param srv_sel_norm_bins List with one element per fleet naming the bins the
#'   mean-one standardization averages over, or \code{NULL} for fleets
#'   standardizing over every bin. Read under \code{"nonparlog"} only. A gear whose
#'   catchability is defined against part of the bin range standardizes over that
#'   part, and catchability absorbs the difference in scale. Default \code{NULL}.
#' @param srv_sel_bin_dev_bins List with one element per fleet naming the bins that
#'   fleet overrides, or \code{NULL} for none, e.g. \code{list(1, NULL)}. An
#'   overridden bin takes a freely estimated annual value
#'   \eqn{\exp(\epsilon_{y,b})} in place of what the functional form produced,
#'   applied after every other transformation including standardization, while the
#'   rest of the curve keeps its parametric shape. Default \code{NULL}.
#' @param cont_tv_srvsel_bin_devs Character vector \code{[n_srv_fleets]} of the
#'   process error on the bin-override deviations: \code{"none"} (default),
#'   \code{"iid"} or \code{"rw"}. A walk has its own estimated sigma per bin.
#' @param srv_q_blocks Character vector of discrete catchability time blocks, in
#'   the format of \code{srv_sel_blocks}. Parsed into an \code{[n_regions × n_years
#'   × n_srv_fleets]} array. Default \code{"none_Fleet_x"}.
#' @param srv_q_spec Character vector \code{[n_srv_fleets]} or \code{NULL}
#'   (default) of the sharing structure for catchability. See
#'   \code{\link{do_q_mapping}}.
#' @param srv_q_type Character vector \code{[n_srv_fleets]} of how catchability is
#'   obtained. \code{"est"} (default) estimates \code{ln_srv_q}, \code{"arith"}
#'   concentrates it out as the ratio of mean observed to mean predicted index, and
#'   \code{"geo"} does the same on the log scale as
#'   \code{exp(mean(log(obs) - log(pred)))}. Both analytic forms use the years with
#'   observations only and fix that fleet's \code{ln_srv_q} whatever
#'   \code{srv_q_spec} says. The solve runs within each \code{srv_q_blocks} block.
#' @param Use_srv_q_prior Integer (0/1) for lognormal priors on survey
#'   catchability. Default \code{0}.
#' @param srv_q_prior Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{mu} on the natural scale and \code{sd} on the log scale,
#'   one row per \eqn{\log\text{N}(\log(\mu), \text{sd})} prior. Default \code{NA}.
#' @param srv_q_model Character vector \code{[n_srv_fleets]} of the process error
#'   on annual catchability deviations: \code{"none"} (default), \code{"iid"},
#'   \code{"rw"}, \code{"ar1"} or \code{"dsem"}, which hands the series to
#'   \code{\link{Setup_Mod_DSEM}}. Catchability is then \eqn{\exp(\ln q_{r,b,f} +
#'   \epsilon_{r,y,f})}. A fleet with deviations cannot also have
#'   \code{srv_q_blocks} or an analytically solved \code{srv_q_type}.
#' @param sigma_srv_q_spec Sharing string for the deviation standard deviation over
#'   region and fleet: \code{"est_all"} (default), \code{"est_shared_r"},
#'   \code{"est_shared_f"}, \code{"est_shared_r_f"} or \code{"fix"}.
#' @param srv_q_rho_spec Sharing string for the AR1 correlation, with the same
#'   options as \code{sigma_srv_q_spec}. Default \code{"est_all"}. Only read under
#'   \code{srv_q_model = "ar1"}.
#' @param srv_q_rw_init_sigma Standard deviation of the first estimated year of a
#'   random walk. \code{NA} (default) starts the walk at zero under its own sigma,
#'   which keeps \code{ln_srv_q} as the level of the series; a wide value leaves
#'   the level free and confounds it with \code{ln_srv_q}.
#' @param t_srv Survey timing as a fraction of the year (annual models) or the
#'   season (seasonal models), array \code{[n_regions × n_seas × n_srv_fleets]}.
#'   Default \code{1}, the end of the period.
#' @param Use_srv_selex_prior Integer (0/1) for priors on the survey selectivity
#'   parameters. Default \code{0}.
#' @param srv_selex_prior Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par}, \code{mu}, \code{sd} and an optional
#'   \code{type}. \code{"par"} (the default) is a lognormal prior on one fixed
#'   selectivity parameter, with \code{mu} on the natural scale and \code{sd} on
#'   the log scale. \code{"value"} is a normal prior on the realized selectivity at
#'   one bin, both on the natural scale, where \code{par} names the bin and the
#'   value is read at the first model year of \code{block}; that is the ADMB
#'   convention of pinning survey selectivity at a reference age near one, which no
#'   set of independent parameter priors can express. Default \code{NULL}.
#' @param srv_selex_type Character scalar, \code{"age"} (default) or
#'   \code{"length"}.
#' @param use_fixed_srv_sel Integer vector \code{[n_srv_fleets]}, \code{1} to fix
#'   survey selectivity and \code{0} to estimate it.
#' @param srv_sel_input Array of fixed survey selectivity values \code{[n_pop ×
#'   n_regions × n_years × n_seas × n_bins × n_sexes × n_srv_fleets]}, required
#'   whenever any survey has fixed selectivity.
#' @param srv_sel_nonpar_est_bins Optional bin groupings for non-parametric survey
#'   selectivity, structured \code{[[survey]][[block]]}, each element a list of
#'   integer vectors naming the bins that share one estimated parameter. Indices
#'   are on the bin dim the survey selectivity type names.
#' @param srvsel_dont_est_dev_first Integer vector \code{[n_srv_fleets]} of 0/1,
#'   default \code{0}. Where \code{1}, that fleet's deviations start in year two and
#'   the fixed parameters hold year one. A non-parametric form has one free base
#'   parameter per bin, so year one's deviation is that same value written twice
#'   with only \code{srvsel_rw_init_sigma} between them, a prior on a level usually
#'   meant to be free. Refused for the GMRF and 2D AR1 forms, whose deviations are a
#'   field over years and bins rather than a walk anchored at year one.
#' @param srv_sel_dbnrml_startbin \code{NULL} (default) or an integer vector
#'   \code{[n_srv_fleets]}, the bin each survey's double normal anchors its
#'   ascending limb at. See \code{fish_sel_dbnrml_startbin} in
#'   \code{\link{Setup_Mod_Fishsel_and_Q}}.
#' @param srv_sel_dbnrml_raw \code{NULL} (default) or a 0/1 matrix
#'   \code{[n_srv_fleets x 2]} for fleets on the double normal: column one leaves
#'   the ascending limb a raw Gaussian instead of anchoring it to \code{p5} at the
#'   first bin, column two does the same for the descending limb and \code{p6}.
#' @param srv_sel_sex_offset Character vector \code{[n_srv_fleets]} linking the
#'   sexes of a fleet's selectivity when \code{n_sexes > 1}: \code{"none"}
#'   (default), \code{"par"}, \code{"scale"}, \code{"par_scale"}, \code{"apical"}
#'   or \code{"par_apical"}. See \code{fish_sel_sex_offset} in
#'   \code{\link{Setup_Mod_Fishsel_and_Q}} for what each one does.
#' @param ... Optional starting values for the selectivity and catchability
#'   parameters.
#'
#' @return \code{input_list} with the selectivity and catchability configuration in
#'   \code{$data} (\code{cont_tv_srv_sel}, \code{srv_sel_blocks},
#'   \code{srv_sel_model}, \code{srv_q_blocks}, \code{srv_q_prior},
#'   \code{Use_srv_q_prior}, \code{srv_q_model}, \code{Use_srv_selex_prior},
#'   \code{srv_selex_prior}, \code{t_srv}), the starting values in \code{$par} for
#'   \code{srv_fixed_sel_pars}, \code{ln_srv_q}, \code{srvsel_pe_pars},
#'   \code{ln_srvsel_devs}, \code{ln_srv_q_devs}, \code{ln_sigma_srv_q} and
#'   \code{srv_q_rho}, and their factor maps in \code{$map}.
#'
#' @export Setup_Mod_Srvsel_and_Q
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Srvsel_and_Q <- function(
  input_list,
  cont_tv_srv_sel = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
  srv_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
  srv_sel_model,
  Use_srv_q_prior = 0,
  srv_q_prior = NA,
  srv_q_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
  srvsel_pe_pars_spec = NULL,
  srv_fixed_sel_pars_spec,
  srv_q_spec = NULL,
  srv_q_type = rep("est", input_list$data$n_srv_fleets),
  srv_sel_devs_spec = NULL,
  corr_opt_semipar = NULL,
  srv_q_model = NULL,
  sigma_srv_q_spec = "est_all",
  srv_q_rho_spec = "est_all",
  srv_q_rw_init_sigma = NA,
  Use_srv_selex_prior = 0,
  srv_selex_prior = NULL,
  Use_srv_selex_penalty = 0,
  srv_sel_norm_bins = NULL,
  srv_sel_bin_dev_bins = NULL,
  srvsel_pe_wt = rep(1, input_list$data$n_srv_fleets),
  srvsel_rw_init_sigma = rep(5, input_list$data$n_srv_fleets),
  cont_tv_srvsel_bin_devs = rep("none", input_list$data$n_srv_fleets),
  srv_selex_penalty = NULL,
  t_srv = array(1, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_srv_fleets)),
  srvsel_devs_shared_bins = NULL,
  srv_selex_type = 'age',
  use_fixed_srv_sel = rep(0, input_list$data$n_srv_fleets),
  srv_sel_input = NULL,
  srv_sel_nonpar_est_bins = NULL,
  srvsel_dont_est_dev_first = rep(0, input_list$data$n_srv_fleets),
  srv_sel_sex_offset = rep("none", input_list$data$n_srv_fleets),
  srv_sel_dbnrml_raw = NULL,
  srv_sel_dbnrml_startbin = NULL,
  ...
) {

  messages_list <<- character(0) # string to attach to for printing messages # nolint: object_usage_linter.
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Srvsel_and_Q <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Continuous Selectivity Deviations
  check_fleet_spec_length(srvsel_pe_pars_spec, input_list$data$n_srv_fleets, "srvsel_pe_pars_spec", allow_null = TRUE)
  check_fleet_spec_length(srv_sel_devs_spec, input_list$data$n_srv_fleets, "srv_sel_devs_spec", allow_null = TRUE)
  check_fleet_spec_length(corr_opt_semipar, input_list$data$n_srv_fleets, "corr_opt_semipar", allow_null = TRUE)

  # A short vector here is read per fleet in the objective, so a length mismatch
  # silently becomes NA rather than being recycled.
  check_fleet_spec_length(srvsel_pe_wt, input_list$data$n_srv_fleets, "srvsel_pe_wt")
  check_fleet_spec_length(srvsel_rw_init_sigma, input_list$data$n_srv_fleets, "srvsel_rw_init_sigma")
  check_fleet_spec_length(srvsel_dont_est_dev_first, input_list$data$n_srv_fleets, "srvsel_dont_est_dev_first")
  if(!all(srvsel_dont_est_dev_first %in% c(0, 1))) stop("srvsel_dont_est_dev_first must be 0 or 1 for every fleet")

  # Catchability Priors
  if(!Use_srv_q_prior %in% c(0,1)) stop("Values for Use_srv_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking catchability priors
  if(Use_srv_q_prior == 1) {
    required_cols <- c("region", "fleet", "block", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(srv_q_prior))
    if(length(missing_cols) > 0) {
      stop("srv_q_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Survey Catchability priors are: ", ifelse(Use_srv_q_prior == 0, "Not Used", "Used"))

  # Selectivity Priors
  if(!Use_srv_selex_prior %in% c(0,1)) stop("Values for Use_srv_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_srv_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(srv_selex_prior))
    if(length(missing_cols) > 0) {
      stop("srv_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Survey Selectivity priors are: ", ifelse(Use_srv_selex_prior == 0, "Not Used", "Used"))

  if(any(use_fixed_srv_sel == 1) && is.null(srv_sel_input)) stop("srv_sel_input is NULL, please provide an input array.")
  if(any(use_fixed_srv_sel == 1) && srv_selex_type == 'age') check_data_dimensions(
    srv_sel_input,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_ages = length(input_list$data$ages),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'srv_sel_input_age'
  )
  if(any(use_fixed_srv_sel == 1) && srv_selex_type == 'length') check_data_dimensions(
    srv_sel_input,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_lens = length(input_list$data$lens),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'srv_sel_input_len'
  )

  # Selectivity Options -----------------------------------------------------
  # The bin vector is kept as well as its length. Starting values stated on the
  # bin scale are seeded further down, by which point srv_selex_type holds the
  # numeric code rather than the name it arrived as.
  if(srv_selex_type == 'age') {
    srv_selex_type <- 0
    srv_sel_bin_vec <- input_list$data$ages
    collect_message("Survey Selectivity is aged-based.")
  } else if(srv_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but survey selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    srv_selex_type <- 1
    srv_sel_bin_vec <- input_list$data$lens
    collect_message("Survey Selectivity is length-based")
  } else stop("srv_selex_type must be 'age' or 'length', but was: ", srv_selex_type)

  bins <- length(srv_sel_bin_vec)

  # Continuous Time-Varying Selectivity Options -----------------------------
  # define for continuous time-varying selectivity
  cont_tv_srv_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in seq_along(cont_tv_srv_sel)) {
    # Extract out components from list
    tmp <- cont_tv_srv_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for cont_tv_srv_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_srv_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_srv_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous survey time-varying selectivity specified as: ", cont_tv_type, " for survey fleet ", fleet)
  }

  if(any(cont_tv_srv_sel_mat > 0) && (is.null(srvsel_pe_pars_spec) || is.null(srv_sel_devs_spec))) stop("Continuous time-varying selectivity specified, but srvsel_pe_pars_spec and/or srv_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  srv_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in seq_along(srv_sel_blocks)) {

    # Extract out components from list
    tmp <- srv_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_sel_blocks_arr[,,fleet] <- 1 # input only 1 survey time block
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

      srv_sel_blocks_arr[,years,fleet] <- block_val
    }
  }

  if(any(is.na(srv_sel_blocks_arr))) stop("Survey Selectivity Blocks are returning an NA. Did you update the year range of srv_sel_blocks?")
  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Selectivity Time Blocks for survey", f, "is specified at:", length(unique(srv_sel_blocks_arr[,,f]))))

  # Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml", 'nonpar', 'asymplogist1', "asymplogist2", "bicubic", "nonparlog", "nonparfree"), num = c(0,1,2,3,4,5,6,7,8,9,10)) # set up values we can map to
  srv_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  srv_sel_bicubic_binnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of bin nodes, only set where srv_sel_model == 8 (bicubic)
  srv_sel_bicubic_yrnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of year nodes, only set where srv_sel_model == 8 (bicubic)
  srv_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year, i.e. no offset); years within the block before this are edge-kept at this year's fitted curve
  srv_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of bins (starting from the first) the bicubic surface is actually fit over (0 = all bins, i.e. no truncation); bins beyond this are kept flat at the last fitted bin's value

  for(i in seq_along(srv_sel_model)) {

    # Extract out survey selectivity components from vector
    tmp_sel_form <- srv_sel_model[i]
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
        stop("srv_sel_model 'bicubic' entries must be specified as bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f> or bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Block_<b>_Fleet_<f>, optionally with _SelStyr_<year> and/or _NSelBins_<n>")
      tmp_n_bin_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[bin_pos + 1]))
      tmp_n_yr_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[yr_pos + 1]))
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_selstyr <- if(length(selstyr_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[selstyr_pos + 1])) else 0
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(is.na(tmp_n_bin_nodes) || tmp_n_bin_nodes < 2) stop("bicubic srv_sel_model requires at least 2 bin nodes (n_bin_nodes >= 2)")
      if(is.na(tmp_n_yr_nodes) || tmp_n_yr_nodes < 1) stop("bicubic srv_sel_model requires at least 1 year node (n_yr_nodes >= 1). Use n_yr_nodes == 1 for a time-invariant bin-only spline.")
      if(length(selstyr_pos) == 1 && (is.na(tmp_selstyr) || !tmp_selstyr %in% input_list$data$years)) stop("bicubic srv_sel_model SelStyr must be a calendar year within the modeled years")
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("bicubic srv_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
    } else {
      # optional _NSelBins_<n> plateau suffix is joined with  _Block_<b> in either order for any parametric form
      fleet_pos <- which(tmp_sel_form_vec == "Fleet")
      block_pos <- which(tmp_sel_form_vec == "Block")
      nselbins_pos <- which(tmp_sel_form_vec == "NSelBins")
      if(length(fleet_pos) != 1) stop("srv_sel_model entries must name their fleet exactly once, as _Fleet_<f>")
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("srv_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
      if(length(nselbins_pos) == 1 && sel_form %in% c("nonpar", "nonparlog", "nonparfree")) stop("srv_sel_model NSelBins is for the parametric forms; a non-parametric form would keep estimating the bins it then overwrites. Group those bins through srv_sel_nonpar_est_bins instead.")
    }

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("srv_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml, nonpar, nonparlog, nonparfree, asymplogist1, asymplogist2, bicubic (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for srv_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) {
      srv_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
      if(sel_form == "bicubic") {
        srv_sel_bicubic_binnodes_arr[,,tmp_fleet] <- tmp_n_bin_nodes
        srv_sel_bicubic_yrnodes_arr[,,tmp_fleet] <- tmp_n_yr_nodes
        srv_sel_bicubic_selstyr_arr[,,tmp_fleet] <- tmp_selstyr
      }
      srv_sel_bicubic_nselbins_arr[,,tmp_fleet] <- tmp_nselbins # plateau bin; read by every form (0 = none)
    } else {
      srv_sel_model_arr <- assign_sel_block(srv_sel_model_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, sel_map$num[which(sel_map$sel == sel_form)])
      if(sel_form == "bicubic") {
        srv_sel_bicubic_binnodes_arr <- assign_sel_block(srv_sel_bicubic_binnodes_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_bin_nodes)
        srv_sel_bicubic_yrnodes_arr <- assign_sel_block(srv_sel_bicubic_yrnodes_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_yr_nodes)
        srv_sel_bicubic_selstyr_arr <- assign_sel_block(srv_sel_bicubic_selstyr_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_selstyr)
      }
      srv_sel_bicubic_nselbins_arr <- assign_sel_block(srv_sel_bicubic_nselbins_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_nselbins)
    }
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Survey selectivity functional form specified as:", sel_form, " for survey fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_srv_fleets) {
    has_blocks <- length(unique(srv_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_srv_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Blocked Catchability Options --------------------------------------------
  srv_q_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in seq_along(srv_q_blocks)) {

    # Extract out components from list
    tmp <- srv_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Vakudate option
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_q_blocks_arr[,,fleet] <- 1 # input only 1 survey catchability time block
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
      srv_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }
  if(any(is.na(srv_q_blocks_arr))) stop("Survey Catchability Blocks are returning an NA. Did you update the year range of srv_q_blocks?")
  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Catchability Time Blocks for survey", f, "is specified at:", length(unique(srv_q_blocks_arr[,,f]))))

  # Populate Data List ------------------------------------------------------
  input_list$data$cont_tv_srv_sel <- cont_tv_srv_sel_mat
  input_list$data$srv_sel_blocks <- srv_sel_blocks_arr
  input_list$data$srv_sel_model <- srv_sel_model_arr
  input_list$data$srv_sel_bicubic_binnodes <- srv_sel_bicubic_binnodes_arr
  input_list$data$srv_sel_bicubic_yrnodes <- srv_sel_bicubic_yrnodes_arr
  input_list$data$srv_sel_bicubic_selstyr <- srv_sel_bicubic_selstyr_arr
  input_list$data$srv_sel_bicubic_nselbins <- srv_sel_bicubic_nselbins_arr
  input_list$data$srv_q_blocks <- srv_q_blocks_arr
  input_list$data$srv_q_prior <- srv_q_prior
  input_list$data$Use_srv_q_prior <- Use_srv_q_prior
  input_list$data$Use_srv_selex_prior <- Use_srv_selex_prior
  input_list$data$srv_selex_prior <- validate_selex_prior_types(srv_selex_prior, Use_srv_selex_prior, "srv_selex_prior",
                                                                sel_blocks = srv_sel_blocks_arr, n_bins = bins)
  input_list$data$Use_srv_selex_penalty <- Use_srv_selex_penalty
  input_list$data$srvsel_pe_wt <- srvsel_pe_wt
  input_list$data$srvsel_rw_init_sigma <- srvsel_rw_init_sigma
  input_list <- setup_sel_bin_devs(
    input_list,
    srv_sel_bin_dev_bins,
    cont_tv_srvsel_bin_devs,
    prefix = "srv",
    n_fleets = input_list$data$n_srv_fleets,
    bins = bins,
    starting_values = starting_values
  )
  input_list <- setup_sel_norm_bins(
    input_list,
    srv_sel_norm_bins,
    prefix = "srv",
    n_fleets = input_list$data$n_srv_fleets,
    bins = bins
  )
  input_list$data$srv_selex_penalty <- validate_selex_penalty(srv_selex_penalty, Use_srv_selex_penalty, "srv_selex_penalty")
  input_list$data$t_srv <- t_srv
  if(!is.null(input_list$data$srv_len_comp_sel) && any(input_list$data$srv_len_comp_sel == 1) && srv_selex_type != 1) stop("SrvLenComps_sel = 'length' in Setup_Mod_SrvIdx_and_Comps applies the length selectivity at length, so srv_selex_type must be length")
  input_list$data$srv_selex_type <- srv_selex_type
  input_list$data$use_fixed_srv_sel <- use_fixed_srv_sel
  input_list$data$srv_sel_input <- srv_sel_input
  input_list$data$srvsel_devs_min_shared_bins <- if(!is.null(srvsel_devs_shared_bins)) unlist(lapply(srvsel_devs_shared_bins, min)) else seq_along(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------
  # Figure out number of selectivity parameters for a given functional form
  unique_srvsel_vals <- unique(as.vector(input_list$data$srv_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in seq_along(unique_srvsel_vals)) {
    if(unique_srvsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_srvsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_srvsel_vals[i] == 4) sel_pars_vec[i] <- 6 # double normal
    if(unique_srvsel_vals[i] %in% c(5,9,10)) sel_pars_vec[i] <- bins # non-parametric selex
    if(unique_srvsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic selex w/ asymptote
    if(unique_srvsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$srv_sel_bicubic_binnodes * input_list$data$srv_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  max_srvsel_blks <- max(apply(input_list$data$srv_sel_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey selectivity blocks for a given reigon and fleet

  # bicubic spline interpolation weight matrices (bin node by year node grid), passed through the
  # model with the flattened node parameters. zero-weight rows and columns never contribute
  has_bicubic_srv_sel <- any(input_list$data$srv_sel_model == 8)
  max_bin_nodes_bicubic <- if(has_bicubic_srv_sel) max(input_list$data$srv_sel_bicubic_binnodes) else 1
  max_yr_nodes_bicubic <- if(has_bicubic_srv_sel) max(input_list$data$srv_sel_bicubic_yrnodes) else 1
  n_yrs_total_bicubic <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs

  srv_sel_bicubic_Wbin <- array(0, dim = c(input_list$data$n_regions, bins, max_bin_nodes_bicubic, max_srvsel_blks, input_list$data$n_srv_fleets))
  srv_sel_bicubic_Wyr <- array(0, dim = c(input_list$data$n_regions, n_yrs_total_bicubic, max_yr_nodes_bicubic, max_srvsel_blks, input_list$data$n_srv_fleets))

  if(has_bicubic_srv_sel) {
    for(f in 1:input_list$data$n_srv_fleets) {
      for(r in 1:input_list$data$n_regions) {

        srvsel_blocks_tmp <- unique(as.vector(input_list$data$srv_sel_blocks[r,,f]))

        for(b in seq_along(srvsel_blocks_tmp)) {

          block_years <- which(input_list$data$srv_sel_blocks[r,,f] == srvsel_blocks_tmp[b])
          if(unique(input_list$data$srv_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$srv_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$srv_sel_bicubic_yrnodes[r, block_years, f])

          # bin nodes evenly spaced over [0,1]. the spline is fit over all bins unless NSelBins is
          # set, beyond which bins are edge-kept at the last fitted bin's weights
          nselbins_this <- unique(input_list$data$srv_sel_bicubic_nselbins[r, block_years, f])
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

          srv_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # year nodes are evenly spaced over the block's fit range, which is the whole block unless
          # SelStyr is set. years outside that range hold the boundary node weights constant

          selstyr_this <- unique(input_list$data$srv_sel_bicubic_selstyr[r, block_years, f])
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

          srv_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this

        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_srv_sel

  input_list$data$srv_sel_bicubic_Wbin <- srv_sel_bicubic_Wbin
  input_list$data$srv_sel_bicubic_Wyr <- srv_sel_bicubic_Wyr

  max_srvsel_pars <- max(sel_pars_vec) # maximum number of selectivity parameters across all forms
  input_list$par$srv_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_srvsel_pars, max_srvsel_blks, input_list$data$n_sexes, input_list$data$n_srv_fleets))
  input_list$par$srv_fixed_sel_pars <- use_starting_value(input_list$par$srv_fixed_sel_pars, starting_values, "srv_fixed_sel_pars")

  # a double normal has its peak on the bin scale, so a default of zero
  # would put it at bin zero, where the ascending limb has no extent
  if(!"srv_fixed_sel_pars" %in% names(starting_values)) {
    input_list$par$srv_fixed_sel_pars <- seed_dbnrml_peak(input_list$par$srv_fixed_sel_pars, srv_sel_model_arr,
                                          srv_sel_bin_vec,
                                          srv_sel_sex_offset)
  }

  # Survey catchability
  max_srvq_blks <- max(apply(input_list$data$srv_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey catchability blocks for a given reigon and fleet
  input_list$par$ln_srv_q <- array(0, dim = c(input_list$data$n_regions, max_srvq_blks, input_list$data$n_srv_fleets))
  input_list$par$ln_srv_q <- use_starting_value(input_list$par$ln_srv_q, starting_values, "ln_srv_q")

  # Survey selectivity process error parameters
  input_list$par$srvsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_srvsel_pars, 4), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using
  input_list$par$srvsel_pe_pars <- use_starting_value(input_list$par$srvsel_pe_pars, starting_values, "srvsel_pe_pars")

  # Survey selectivity deviations
  input_list$par$ln_srvsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_srv_fleets))
  input_list$par$ln_srvsel_devs <- use_starting_value(input_list$par$ln_srvsel_devs, starting_values, "ln_srvsel_devs")

  # Sex offsets on selectivity (parameter offsets and/or a curve scale offset)
  input_list <- setup_sel_sex_offset(
    input_list,
    srv_sel_sex_offset,
    prefix = "srv",
    n_fleets = input_list$data$n_srv_fleets,
    fleet_label = "survey fleet",
    sel_model_arr = input_list$data$srv_sel_model,
    cont_tv_mat = cont_tv_srv_sel_mat,
    max_blks = max_srvsel_blks,
    sel_blocks = input_list$data$srv_sel_blocks,
    fixed_spec = srv_fixed_sel_pars_spec,
    starting_values = starting_values
  )

  if(!is.null(input_list$data$srv_waa_selected) && any(input_list$data$srv_waa_selected == 1) && srv_selex_type != 1) stop("srv_waa_selected = 1 in Setup_Mod_SrvIdx_and_Comps weights the survey weight at age by length selectivity, so srv_selex_type must be length")

  # Raw (unanchored) double normal limbs
  input_list$data$srv_dbnrml_raw <- setup_dbnrml_raw(srv_sel_dbnrml_raw, input_list$data$n_srv_fleets, "srv_sel_dbnrml_raw")
  input_list$data$srv_dbnrml_startbin <- setup_dbnrml_startbin(srv_sel_dbnrml_startbin, input_list$data$n_srv_fleets, bins, "srv_sel_dbnrml_startbin")

  ## Parameter Maps ---------------------------------------------------------
  ## Catchability solving ---------------------------------------------------
  # A fleet whose catchability is concentrated out of the likelihood has no
  # free ln_srv_q, so its mapping is fixed regardless of what srv_q_spec specifies.
  if(!all(srv_q_type %in% c("est", "arith", "geo"))) stop("Invalid specification for srv_q_type. Should be est, arith, or geo")
  check_fleet_spec_length(srv_q_type, input_list$data$n_srv_fleets, "srv_q_type")

  srv_q_type_vals <- convert_to_numeric(srv_q_type, list(est = 0, arith = 1, geo = 2))
  input_list$data$srv_q_type <- srv_q_type_vals

  for(f in 1:input_list$data$n_srv_fleets) {
    collect_message(paste("Survey Catchability for survey fleet", f, "is:",
                          switch(srv_q_type[f],
                                 est = "estimated",
                                 arith = "solved analytically as the ratio of mean observed to mean predicted",
                                 geo = "solved analytically on the log scale")))
  } # end f loop

  if(any(srv_q_type_vals != 0)) {
    if(is.null(srv_q_spec)) srv_q_spec <- rep("est_all", input_list$data$n_srv_fleets)
    srv_q_spec[srv_q_type_vals != 0] <- "fix"
  }

  # Mapping Options ---------------------------------------------------------
  input_list <- do_fixed_sel_pars_mapping(
    input_list,
    srv_fixed_sel_pars_spec,
    bins,
    srv_sel_nonpar_est_bins,
    prefix = "srv",
    fleet_field = "n_srv_fleets",
    use_field = "SrvIdx",
    fleet_label = "survey fleet"
  )
  input_list <- do_q_mapping(
    input_list,
    srv_q_spec,
    prefix = "srv",
    fleet_field = "n_srv_fleets",
    fleet_label = "survey fleet"
  )
  input_list <- setup_q_devs(
    input_list,
    q_model = srv_q_model,
    sigma_q_spec = sigma_srv_q_spec,
    q_rho_spec = srv_q_rho_spec,
    q_rw_init_sigma = srv_q_rw_init_sigma,
    q_type = srv_q_type,
    prefix = "srv",
    fleet_field = "n_srv_fleets",
    use_field = "UseSrvIdx",
    fleet_label = "survey fleet",
    starting_values = starting_values
  )
  input_list <- do_sel_pe_pars_mapping(
    input_list,
    srvsel_pe_pars_spec,
    corr_opt_semipar,
    bins,
    sel_devs_spec = srv_sel_devs_spec,
    sel_devs_shared_bins = srvsel_devs_shared_bins,
    prefix = "srv",
    fleet_field = "n_srv_fleets",
    use_field = "SrvIdx",
    fleet_label = "survey fleet"
  )
  input_list <- do_sel_devs_mapping(
    input_list,
    srv_sel_devs_spec,
    srvsel_devs_shared_bins,
    bins,
    dont_est_dev_first = srvsel_dont_est_dev_first,
    prefix = "srv",
    fleet_field = "n_srv_fleets",
    use_field = "SrvIdx",
    fleet_label = "survey fleet"
  )


  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
