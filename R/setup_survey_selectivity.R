# Stage 1 of 3: model setup
#
# Survey selectivity and catchability inputs. Setup_Mod_Srvsel_and_Q covers the
# selectivity of the gear on the stock and the catchability blocks, and
# delegates its parameter maps to the shared builders in setup_mapping.R.
# Mirrors setup_fishery_selectivity.R but has no retention stream.

#' Set up survey selectivity and catchability specifications
#'
#' Configures all survey selectivity and catchability components of the
#' estimation model: continuous and blocked time-varying selectivity,
#' selectivity functional forms, catchability blocks and optional
#' environmental covariate effects, process error and deviation mapping, and
#' selectivity/catchability priors. Delegates parameter mapping to four
#' internal helpers (\code{\link{do_fixed_sel_pars_mapping}},
#' \code{\link{do_q_mapping}}, \code{\link{do_sel_pe_pars_mapping}},
#' \code{\link{do_sel_devs_mapping}}). Must be called after
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}} and before model compilation.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#'   \code{$data$srv_selex_type} must already be set by
#'   \code{\link{Setup_Mod_Biologicals}}.
#'
#' @param srv_sel_model Character vector specifying the selectivity functional
#'   form per fleet, and optionally per time block. Each element follows one
#'   of:
#'   \describe{
#'     \item{\code{"<model>_Fleet_x"}}{Single form applied across all years for
#'       fleet \code{x}.}
#'     \item{\code{"<model>_Fleet_x_Block_k"}}{Form applied only to block
#'       \code{k} for fleet \code{x}, as defined in \code{srv_sel_blocks}.
#'       Required when multiple blocks are defined for a fleet.}
#'   }
#'   Available models (see the model equations vignette for parameterisations):
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'     \item{\code{"nonpar"}}{Non-parametric selectivity defined over discrete age or length bins, where selectivity is estimated as independent parameters (or grouped bins if specified via nonparametric bin mapping). No fixed functional form is imposed.}
#'     \item{\code{"asymplogist1"}}{Logistic selectivity with \eqn{a_{50}} and slope \eqn{k} and asymptotic control (3 parameters).}
#'     \item{\code{"asymplogist2"}}{Logistic selectivity with with \eqn{a_{50}} and \eqn{a_{95}} and asymptotic control (3 parameters).}
#'     \item{\code{"bicubic"}}{Bicubic spline over a bin-node x year-node grid
#'       (see \code{\link{Get_Selex}}, \code{Selex_Model == 8}). Specified as
#'       \code{"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"} (optionally
#'       with \code{_Block_k}). One generalized form covers a smooth bin x year
#'       surface (\code{n_yr_nodes > 1}), a time-invariant bin-only spline
#'       (\code{n_yr_nodes == 1}), or a bin-only spline re-fit independently
#'       per year-block (\code{n_yr_nodes == 1} within each of several blocks
#'       defined via \code{srv_sel_blocks}). An optional \code{_SelStyr_<year>}
#'       suffix (a calendar year within the block) restricts the actual spline
#'       fit to \code{SelStyr}:block-end; years within the block before
#'       \code{SelStyr} are held constant at the \code{SelStyr} year's fitted
#'       curve, rather than fitting the surface over the whole block. An
#'       optional \code{_NSelBins_<n>} suffix restricts the actual spline fit
#'       to the first \code{n} bins (ages or lengths, per \code{srv_selex_type});
#'       bins beyond \code{n} are held constant at the last fitted bin's
#'       curve.}
#'   }
#'   No default; must be provided.
#' @param srv_fixed_sel_pars_spec Character vector \code{[n_srv_fleets]}.
#'   Sharing structure for fixed-effect selectivity parameters. See
#'   \code{\link{do_fixed_sel_pars_mapping}} for full option descriptions.
#'   No default; must be provided.
#' @param srv_sel_blocks Character vector defining discrete time blocks for
#'   survey selectivity. Each element follows \code{"Block_k_Year_a-b_Fleet_x"}
#'   or \code{"none_Fleet_x"} (constant selectivity). Use \code{"terminal"} in
#'   place of the end year to extend to the final model year. Parsed into an
#'   internal array \code{[n_regions × n_years × n_srv_fleets]}. Blocked and
#'   continuous time-varying selectivity are mutually exclusive for a given
#'   fleet. Default: \code{"none_Fleet_x"} for each fleet.
#'
#' @param cont_tv_srv_sel Character vector defining the continuous
#'   time-variation form per fleet. Each element follows
#'   \code{"<type>_Fleet_x"}. Options:
#'   \describe{
#'     \item{\code{"none"}}{No continuous time variation (default).}
#'     \item{\code{"iid"}}{IID deviations across years.}
#'     \item{\code{"rw"}}{Random walk in time.}
#'     \item{\code{"3dmarg"}}{3D marginal GMRF (age × year × cohort).}
#'     \item{\code{"3dcond"}}{3D conditional GMRF.}
#'     \item{\code{"2dar1"}}{2D AR1 (bin × year).}
#'   }
#'   When any fleet uses a non-\code{"none"} type, both
#'   \code{srvsel_pe_pars_spec} and \code{srv_sel_devs_spec} must be
#'   specified. Default: \code{"none_Fleet_x"} for each fleet.
#' @param srvsel_pe_pars_spec Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Sharing structure for process error hyperparameters. See
#'   \code{\link{do_sel_pe_pars_mapping}} for full option descriptions.
#'   Default \code{NULL}.
#' @param srv_sel_devs_spec Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Sharing structure for selectivity deviation time series.
#'   See \code{\link{do_sel_devs_mapping}} for full option descriptions.
#'   Default \code{NULL}.
#' @param Use_srv_selex_penalty Integer (0/1). Whether a centering penalty is
#'   applied to sets of survey selectivity fixed-effect parameters. Default
#'   \code{0}.
#' @param srv_selex_penalty Data frame of centering penalty specifications,
#'   required when \code{Use_srv_selex_penalty = 1}. Required columns:
#'   \code{region}, \code{fleet}, \code{block}, \code{sex}, \code{par}, and
#'   \code{wt}. Each row penalizes \code{wt * (log(mean(exp(pars))))^2} over the
#'   set of parameters named in \code{par}, which may be a single index or a
#'   list column of integer vectors naming a whole set. This pins the scalar of
#'   a non-parametric curve that catchability or fishing mortality would
#'   otherwise absorb, and is softer than fixing a bin outright. Intended for
#'   parameter sets held on the log scale. Default \code{NULL}.
#' @param srvsel_devs_shared_bins List of integer vectors defining bin groups
#'   for age/length-sharing of deviations under semi-parametric forms (e.g.,
#'   \code{list(1:5, 6:10, 11:30)}). Required when \code{srv_sel_devs_spec}
#'   includes any \code{"est_shared_b"} variant. Default \code{NULL}.
#' @param corr_opt_semipar Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Specifies correlation components to suppress for 3D GMRF or
#'   2D AR1 forms. See \code{\link{do_sel_pe_pars_mapping}} for valid
#'   values. Default \code{NULL}.
#' @param srvsel_pe_wt Numeric vector \code{[n_srv_fleets]}. Per-fleet
#'   multiplier on the survey selectivity process error likelihood. Default
#'   \code{1} for every fleet. \code{0} skips that fleet's process error
#'   likelihood altogether, so the deviations stay estimated but enter the
#'   objective only through the data and any explicit smoothness or centering
#'   penalties, which is how several existing assessments constrain them. Values
#'   other than 0 or 1 make an estimated process error sigma reinterpretable, so
#'   prefer 0 or 1 unless deliberately down-weighting. Applies only to
#'   \code{ln_srvsel_devs}; the bin-override deviations carry their own process
#'   error and are not affected.
#' @param srvsel_rw_init_sigma Numeric vector \code{[n_srv_fleets]}. Standard
#'   deviation given to the first year of an \code{"rw"} deviation series.
#'   Default \code{5}, which leaves that year effectively free. \code{NA}
#'   instead starts the walk at zero under the walk's own estimated sigma,
#'   making the first year as smooth as every later step. Appropriate when the
#'   base parametric curve already describes the first year well.
#' @param srv_sel_bin_dev_bins List with one element per survey fleet naming the
#'   bins that fleet overrides, or \code{NULL} for fleets with no overrides
#'   (e.g. \code{list(1, NULL)} frees bin 1 of fleet 1 only). An overridden bin
#'   takes a freely estimated annual value \eqn{\exp(\epsilon_{y,b})} in place of
#'   whatever the functional form produced, applied after every other
#'   transformation including standardization. The rest of the curve keeps its
#'   parametric shape. Default \code{NULL}.
#' @param cont_tv_srvsel_bin_devs Character vector \code{[n_srv_fleets]} giving
#'   the process error on the bin-override deviations for each fleet:
#'   \code{"none"} (default), \code{"iid"}, or \code{"rw"}. A random walk
#'   carries its own estimated sigma per bin, with
#'   \code{srvsel_bin_devs_rw_init_sigma} governing its first year.
#'
#' @param srv_q_blocks Character vector defining discrete time blocks for
#'   survey catchability. Same format as \code{srv_sel_blocks}:
#'   \code{"Block_k_Year_a-b_Fleet_x"} or \code{"none_Fleet_x"}. Parsed into
#'   an array \code{[n_regions × n_years × n_srv_fleets]}. Default:
#'   \code{"none_Fleet_x"} for each fleet.
#' @param srv_q_spec Character vector \code{[n_srv_fleets]} or \code{NULL}.
#'   Sharing structure for catchability. See \code{\link{do_q_mapping}}
#'   for full option descriptions. Default \code{NULL}.
#' @param srv_q_type Character vector \code{[n_srv_fleets]} controlling how
#'   catchability is obtained. \code{"est"} (default) estimates
#'   \code{ln_srv_q}. \code{"arith"} concentrates it out of the likelihood as
#'   the ratio of mean observed to mean predicted index, and \code{"geo"} does
#'   the same on the log scale as \code{exp(mean(log(obs) - log(pred)))}. Both
#'   analytic forms solve one catchability per region and fleet using only the
#'   years with observations, ignore any block structure, and fix that fleet's
#'   \code{ln_srv_q} regardless of \code{srv_q_spec}.
#' @param Use_srv_q_prior Integer (0/1). Whether log-normal priors are applied
#'   to survey catchability parameters. Default \code{0}.
#' @param srv_q_prior Data frame of catchability prior specifications. Required
#'   columns: \code{region}, \code{fleet}, \code{block}, \code{mu} (prior mean
#'   on natural scale), \code{sd} (prior SD on log scale). Each row specifies
#'   a \eqn{\log\text{N}(\log(\mu), \text{sd})} prior. Ignored when
#'   \code{Use_srv_q_prior = 0}. Default \code{NA}.
#' @param srv_q_formula Named list of R formulas specifying environmental
#'   covariate relationships for catchability per region-fleet combination.
#'   Names follow the convention \code{"Region_r_Fleet_f"}. Covariates must
#'   be present in \code{srv_q_cov_dat}. If \code{NULL}, no covariate effects
#'   are included. Default \code{NULL}.
#' @param srv_q_cov_dat Named list of numeric vectors (length = \code{n_years})
#'   containing covariate time series referenced in \code{srv_q_formula}.
#'   All vectors must be the same length and contain no missing values; set
#'   values to \code{0} for years when the survey is not active. If
#'   \code{NULL}, covariate effects are excluded. Default \code{NULL}.
#'
#' @param t_srv Survey timing fraction within a given year (annual models) or
#'   season (seasonal models), array
#'   \code{[n_regions × n_seas × n_srv_fleets]}. Default: \code{1}
#'   (end of period).
#' @param Use_srv_selex_prior Integer (0/1). Whether log-normal priors are
#'   applied to survey selectivity parameters. Default \code{0}.
#' @param srv_selex_prior Data frame of selectivity prior specifications, one
#'   row per prior. Required columns: \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par}, \code{mu}, \code{sd}, plus an
#'   optional \code{type} giving each row's target: \code{"par"} (the default
#'   when the column is absent) is a lognormal prior on one fixed selectivity
#'   parameter, with \code{mu} on the natural scale and \code{sd} on the log
#'   scale; \code{"value"} is a normal prior on the realized selectivity value
#'   at one bin, with both on the natural scale, where \code{par} instead names
#'   the bin (on ages or lengths per \code{srv_selex_type}) and the value is
#'   read at the first model year of \code{block}. A \code{"value"} row
#'   constrains the derived selectivity value rather than the parameters,
#'   matching the ADMB convention of pinning survey selectivity at a reference
#'   age near one, which no set of independent parameter priors can express.
#'   Ignored when \code{Use_srv_selex_prior = 0}. Default \code{NULL}.
#'
#' @param ... Optional named starting values for selectivity and catchability
#'   parameters.
#' @param srv_selex_type Character. Whether survey selectivity type is 'age' or 'length' based. Default: \code{age}.
#' @param use_fixed_srv_sel Integer vector of length \code{n_srv_fleets}
#'   indicating whether survey selectivity is fixed (\code{1}) or estimated
#'   (\code{0}) for each survey index.
#'
#' @param srv_sel_input Array of fixed survey selectivity values used when
#'   \code{use_fixed_srv_sel == 1}. Dimensions:
#'   \code{[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_srv_fleets]}.
#'   Required whenever any survey has fixed selectivity specified.
#'
#' @param srv_sel_nonpar_est_bins Optional list defining bin groupings for
#'   non-parametric survey selectivity. Structure is
#'   \code{[[survey]][[block]]}, where each element is a list of integer vectors.
#'   Each vector defines a group of bins that share a single estimated
#'   selectivity parameter. Indices must correspond to the bin dimension
#'   defined by the survey selectivity type (age or length).
#'
#' @return The input \code{input_list} with selectivity and catchability
#'   configuration stored in \code{$data} (\code{cont_tv_srv_sel}, \code{srv_sel_blocks},
#'   \code{srv_sel_model}, \code{srv_q_blocks}, \code{srv_q_prior},
#'   \code{Use_srv_q_prior}, \code{do_srv_q_cov}, \code{srv_q_cov},
#'   \code{Use_srv_selex_prior}, \code{srv_selex_prior}, \code{t_srv});
#'   starting values in \code{$par} for \code{srv_fixed_sel_pars},
#'   \code{ln_srv_q}, \code{srvsel_pe_pars}, \code{ln_srvsel_devs}, and
#'   \code{srv_q_coeff}; and factor maps in \code{$map} for all five
#'   parameter arrays plus \code{srv_q_coeff}.
#'
#' @export Setup_Mod_Srvsel_and_Q
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Srvsel_and_Q <- function(input_list,
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
                                   srv_q_formula = NULL,
                                   srv_q_cov_dat = NULL,
                                   Use_srv_selex_prior = 0,
                                   srv_selex_prior = NULL,
                                   Use_srv_selex_penalty = 0,
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
                                   ...
                                   ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Srvsel_and_Q <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Continuous Selectivity Deviations
  if(!is.null(srvsel_pe_pars_spec)) if(length(srvsel_pe_pars_spec) != input_list$data$n_srv_fleets) stop("srvsel_pe_pars_spec is not length n_srv_fleets")
  if(!is.null(srv_sel_devs_spec)) if(length(srv_sel_devs_spec) != input_list$data$n_srv_fleets) stop("srv_sel_devs_spec is not length n_srv_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_srv_fleets) stop("corr_opt_semipar is not length n_srv_fleets")

  # A short vector here is read per fleet in the objective, so a length mismatch
  # silently becomes NA rather than being recycled.
  if(length(srvsel_pe_wt) != input_list$data$n_srv_fleets) stop("srvsel_pe_wt is not length n_srv_fleets")
  if(length(srvsel_rw_init_sigma) != input_list$data$n_srv_fleets) stop("srvsel_rw_init_sigma is not length n_srv_fleets")

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
  if(any(use_fixed_srv_sel == 1) && srv_selex_type == 'age') check_data_dimensions(srv_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'srv_sel_input_age')
  if(any(use_fixed_srv_sel == 1) && srv_selex_type == 'length') check_data_dimensions(srv_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'srv_sel_input_len')

  # Selectivity Options -----------------------------------------------------
  # Age based selectivity
  if(srv_selex_type == 'age') {
    srv_selex_type <- 0
    bins <- length(input_list$data$ages)
    collect_message("Survey Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(srv_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but survey selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    srv_selex_type <- 1
    bins <- length(input_list$data$lens)
    collect_message("Survey Selectivity is length-based")
  } # if length based

  # Continuous Time-Varying Selectivity Options -----------------------------
  # define for continuous time-varying selectivity
  cont_tv_srv_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_srv_sel)) {
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

  if(any(cont_tv_srv_sel_mat > 0) && is.null(srvsel_pe_pars_spec) && is.null(srv_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but srvsel_pe_pars_spec and/or srv_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  srv_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_sel_blocks)) {

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
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml", 'nonpar', 'asymplogist1', "asymplogist2", "bicubic", "nonparlog"), num = c(0,1,2,3,4,5,6,7,8,9)) # set up values we can map to
  srv_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  srv_sel_bicubic_binnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of bin nodes, only set where srv_sel_model == 8 (bicubic)
  srv_sel_bicubic_yrnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of year nodes, only set where srv_sel_model == 8 (bicubic)
  srv_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year, i.e. no offset); years within the block before this are edge-held at this year's fitted curve
  srv_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of bins (starting from the first) the bicubic surface is actually fit over (0 = all bins, i.e. no truncation); bins beyond this are held flat at the last fitted bin's value

  for(i in 1:length(srv_sel_model)) {

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
      # get fleet index
      tmp_fleet <- if(length(tmp_sel_form_vec) == 3) as.numeric(tmp_sel_form_vec[3]) else as.numeric(tmp_sel_form_vec[5]) # fleet index changes if block is included in character vector
      # get block index
      tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[3]) else NULL
    }

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("srv_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml, nonpar, asymplogist1, asymplogist2, bicubic (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for srv_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) {
      srv_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
      if(sel_form == "bicubic") {
        srv_sel_bicubic_binnodes_arr[,,tmp_fleet] <- tmp_n_bin_nodes
        srv_sel_bicubic_yrnodes_arr[,,tmp_fleet] <- tmp_n_yr_nodes
        srv_sel_bicubic_selstyr_arr[,,tmp_fleet] <- tmp_selstyr
        srv_sel_bicubic_nselbins_arr[,,tmp_fleet] <- tmp_nselbins
      }
    } else {
      srv_sel_model_arr <- assign_sel_block(srv_sel_model_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, sel_map$num[which(sel_map$sel == sel_form)])
      if(sel_form == "bicubic") {
        srv_sel_bicubic_binnodes_arr <- assign_sel_block(srv_sel_bicubic_binnodes_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_bin_nodes)
        srv_sel_bicubic_yrnodes_arr <- assign_sel_block(srv_sel_bicubic_yrnodes_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_n_yr_nodes)
        srv_sel_bicubic_selstyr_arr <- assign_sel_block(srv_sel_bicubic_selstyr_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_selstyr)
        srv_sel_bicubic_nselbins_arr <- assign_sel_block(srv_sel_bicubic_nselbins_arr, srv_sel_blocks_arr, tmp_fleet, tmp_block, tmp_nselbins)
      }
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
  for(i in 1:length(srv_q_blocks)) {

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

  # Covariate Catchability Options ------------------------------------------
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) collect_message("Using covariates to predict survey catchability")

  # Figure out the total number of regression coefficients that could be estimated
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) {
    n_srv_q_cov <- max(sapply(names(srv_q_formula), function(key) { # sapply to extract names from formula
      tmp_formula <- srv_q_formula[[key]] # get formula
      var_names <- all.vars(tmp_formula) # get var names
      tmp_dat <- data.frame(srv_q_cov_dat[var_names]) # make dataframe
      ncol(stats::model.matrix(tmp_formula, data = tmp_dat)) # figure out number of columns for formula (max number of coefficients to estimate)
    }))
  } else {
    do_srv_q_cov <- 0 # Indicator for whether covariates are included into survey catchability
    n_srv_q_cov <- 1 # dummy to initialize the array
  }

  # Catchability covariate containers
  srv_q_cov <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets, n_srv_q_cov)) # environmental time series
  srv_q_coeff <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets, n_srv_q_cov)) # coefficients to be estimated
  map_srv_q_coeff <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets, n_srv_q_cov)) # coefficients to be mapped off

  # Loop through to map stuff off and populate containers
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) {

    # Validate covariate length
    cov_lengths <- lengths(srv_q_cov_dat)

    # Validate options
    # Check all covariates are the same length
    if (length(unique(cov_lengths)) != 1) stop("All covariates in 'srv_q_cov_dat' must have the same length. If some years are missing data, either impute some value, or set at 0 (if it is not used in the calculation).")
    # Check that length matches the model year structure
    if (unique(cov_lengths) != length(input_list$data$years)) stop(paste0("Covariate length mismatch: expected ",  length(input_list$data$years),  " years but got ", unique(cov_lengths),  "."))

    do_srv_q_cov <- 1 # Indicator for whether covariates are included into survey catchability
    coeff_counter <- 0 # setup counter for mapping

    for(r in 1:input_list$data$n_regions) {
      for(f in 1:input_list$data$n_srv_fleets) {

        # Get key to index
        key <- paste(paste("Region", r, sep = "_"), "_Fleet_", f, sep = "")
        # get temporary formula
        tmp_formula <- srv_q_formula[[key]]
        # extract variable names
        var_names <- all.vars(tmp_formula)
        if(length(var_names) == 0) next # skip if no variables
        # get environmental covariates from environmental data list, based on model formula
        tmp_dat <- data.frame(srv_q_cov_dat[var_names])
        # Generate design matrix
        tmp_design_mat <- stats::model.matrix(tmp_formula, data = tmp_dat)
        # store covariate effects into container
        srv_q_cov[r,,f,1:ncol(tmp_design_mat)] <- tmp_design_mat

        # setup mapping - assign unique counter values for each coefficient
        for(i in 1:ncol(tmp_design_mat)) {
          coeff_counter <- coeff_counter + 1
          map_srv_q_coeff[r,f,i] <- coeff_counter
        } # end i loop

      } # end sf loop
    } # end r loop
  } # if using covariates

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
  input_list$data$do_srv_q_cov <- do_srv_q_cov
  input_list$data$srv_q_cov <- srv_q_cov
  input_list$data$Use_srv_selex_prior <- Use_srv_selex_prior
  input_list$data$srv_selex_prior <- validate_selex_prior_types(srv_selex_prior, Use_srv_selex_prior, "srv_selex_prior",
                                                                sel_blocks = srv_sel_blocks_arr, n_bins = bins)
  input_list$data$Use_srv_selex_penalty <- Use_srv_selex_penalty
  input_list$data$srvsel_pe_wt <- srvsel_pe_wt
  input_list$data$srvsel_rw_init_sigma <- srvsel_rw_init_sigma
  input_list <- setup_sel_bin_devs(input_list, srv_sel_bin_dev_bins, cont_tv_srvsel_bin_devs,
                                   prefix = "srv", n_fleets = input_list$data$n_srv_fleets,
                                   bins = bins, starting_values = starting_values)
  input_list$data$srv_selex_penalty <- validate_selex_penalty(srv_selex_penalty, Use_srv_selex_penalty, "srv_selex_penalty")
  input_list$data$t_srv <- t_srv
  input_list$data$srv_selex_type <- srv_selex_type
  input_list$data$use_fixed_srv_sel <- use_fixed_srv_sel
  input_list$data$srv_sel_input <- srv_sel_input
  input_list$data$srvsel_devs_min_shared_bins <- if(!is.null(srvsel_devs_shared_bins)) unlist(lapply(srvsel_devs_shared_bins, min)) else 1:length(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------
  # Figure out number of selectivity parameters for a given functional form
  unique_srvsel_vals <- unique(as.vector(input_list$data$srv_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_srvsel_vals)) {
    if(unique_srvsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_srvsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_srvsel_vals[i] == 4) sel_pars_vec[i] <- 6 # double normal
    if(unique_srvsel_vals[i] %in% c(5,9)) sel_pars_vec[i] <- bins # non-parametric selex
    if(unique_srvsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic selex w/ asymptote
    if(unique_srvsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$srv_sel_bicubic_binnodes * input_list$data$srv_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  max_srvsel_blks <- max(apply(input_list$data$srv_sel_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey selectivity blocks for a given reigon and fleet

  # Bicubic spline interpolation weight matrices (bin node x year node grid), built here so they can be
  # threaded through SPoRC_rtmb.R alongside the flattened node parameters (see Get_Selex, Selex_Model == 8).
  # Padded with zeros to a common width across regions/blocks/fleets; padding is harmless because unused
  # (zero-weight) columns/rows never contribute to the resulting selectivity (see Get_Selex documentation).
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

        for(b in 1:length(srvsel_blocks_tmp)) {

          block_years <- which(input_list$data$srv_sel_blocks[r,,f] == srvsel_blocks_tmp[b])
          if(unique(input_list$data$srv_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$srv_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$srv_sel_bicubic_yrnodes[r, block_years, f])

          # Bin dimension: nodes evenly spaced over [0,1]. By default (NSelBins unset, i.e. 0) the
          # spline is evaluated over all bins, as before. When NSelBins is set, the spline surface is only actually fit over the first NSelBins bins;
          # bins beyond that are edge-held at the last fitted bin's weights ("plateau").
          nselbins_this <- unique(input_list$data$srv_sel_bicubic_nselbins[r, block_years, f])
          n_fit_bins <- if(nselbins_this == 0) bins else nselbins_this

          bin_nodes_scaled <- seq(0, 1, length.out = n_bin_nodes_this)
          fit_bin_scaled <- seq(0, 1, length.out = n_fit_bins)
          Wbin_fit <- Get_Natural_Cubic_Spline_Weights(bin_nodes_scaled, fit_bin_scaled)

          Wbin_this <- matrix(0, nrow = bins, ncol = n_bin_nodes_this)
          Wbin_this[1:n_fit_bins, ] <- Wbin_fit
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(Wbin_fit[nrow(Wbin_fit), ], nrow = bins - n_fit_bins, ncol = n_bin_nodes_this, byrow = TRUE)

          srv_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # Year dimension: nodes evenly spaced over the block's own contiguous fit range. By default
          # (SelStyr unset, i.e. 0) the fit range is the whole block, as before. When SelStyr is set , only years from SelStyr through the block's end are
          # actually spline-fit; years within the block before SelStyr are edge-held at the SelStyr
          # row's weights ("previous years are filled"). Years outside the block entirely (before it,
          # after it, and any projection years, since projections reuse the terminal modeled year's
          # block) hold the boundary node weights constant, which for a spline evaluated exactly at
          # its first/last node reduces to full weight on that node.
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
          if(max(block_years) < n_yrs_total_bicubic) Wyr_this[(max(block_years) + 1):n_yrs_total_bicubic, ] <- matrix(Wyr_block[nrow(Wyr_block), ], nrow = n_yrs_total_bicubic - max(block_years), ncol = n_yr_nodes_this, byrow = TRUE)

          srv_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this

        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_srv_sel

  input_list$data$srv_sel_bicubic_Wbin <- srv_sel_bicubic_Wbin
  input_list$data$srv_sel_bicubic_Wyr <- srv_sel_bicubic_Wyr

  max_srvsel_pars <- max(sel_pars_vec) # maximum number of selectivity parameters across all forms
  if("srv_fixed_sel_pars" %in% names(starting_values)) input_list$par$srv_fixed_sel_pars <- starting_values$srv_fixed_sel_pars
  else input_list$par$srv_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_srvsel_pars, max_srvsel_blks, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability
  max_srvq_blks <- max(apply(input_list$data$srv_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey catchability blocks for a given reigon and fleet
  if("ln_srv_q" %in% names(starting_values)) input_list$par$ln_srv_q <- starting_values$ln_srv_q
  else input_list$par$ln_srv_q <- array(0, dim = c(input_list$data$n_regions, max_srvq_blks, input_list$data$n_srv_fleets))

  # Survey selectivity process error parameters
  if("srvsel_pe_pars" %in% names(starting_values)) input_list$par$srvsel_pe_pars <- starting_values$srvsel_pe_pars
  else input_list$par$srvsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_srvsel_pars, 4), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Survey selectivity deviations
  if("ln_srvsel_devs" %in% names(starting_values)) input_list$par$ln_srvsel_devs <- starting_values$ln_srvsel_devs
  else input_list$par$ln_srvsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability covariate effects
  if("srv_q_coeff" %in% names(starting_values)) input_list$par$srv_q_coeff <- starting_values$srv_q_coeff
  else input_list$par$srv_q_coeff <- srv_q_coeff # input parameter array

  # Catchability solving ----------------------------------------------------
  # A fleet whose catchability is concentrated out of the likelihood carries no
  # free ln_srv_q, so its mapping is fixed regardless of what srv_q_spec asks for.
  if(!all(srv_q_type %in% c("est", "arith", "geo"))) stop("Invalid specification for srv_q_type. Should be est, arith, or geo")
  if(length(srv_q_type) != input_list$data$n_srv_fleets) stop("srv_q_type is not length n_srv_fleets")

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
  input_list$map$srv_q_coeff <- factor(map_srv_q_coeff) # set up mapping for catchability covariate
  input_list <- do_fixed_sel_pars_mapping(input_list, srv_fixed_sel_pars_spec, bins, srv_sel_nonpar_est_bins,
                                          prefix = "srv", fleet_field = "n_srv_fleets", use_field = "SrvIdx", fleet_label = "survey fleet")
  input_list <- do_q_mapping(input_list, srv_q_spec, prefix = "srv", fleet_field = "n_srv_fleets", fleet_label = "survey fleet")
  input_list <- do_sel_pe_pars_mapping(input_list, srvsel_pe_pars_spec, corr_opt_semipar, bins,
                                       prefix = "srv", fleet_field = "n_srv_fleets", use_field = "SrvIdx", fleet_label = "survey fleet")
  input_list <- do_sel_devs_mapping(input_list, srv_sel_devs_spec, srvsel_devs_shared_bins, bins,
                                    prefix = "srv", fleet_field = "n_srv_fleets", use_field = "SrvIdx", fleet_label = "survey fleet")


  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
