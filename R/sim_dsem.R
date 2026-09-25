# Operating model
#
# dsem series in the operating model: covariates and linked deviations drawn jointly for every
# replicate in Setup_sim_env, then written into the arrays the population loop reads.

#' Set up dsem series in a simulation list
#'
#' Two ways in. From a fit: pass \code{data} and \code{pars} and the operating model draws
#' from the fitted arrows at their estimates, either every year or with the fitted years kept
#' and only later years drawn. From scratch: pass \code{dsem_arrows} and \code{dsem_values}
#' and it draws from arrows you wrote at values you chose, which is how a self test states its
#' truth.
#'
#' @section What a drawn series feeds:
#' Each linked series is written into the array the operating model reads, at the cells the
#' fit linked, so those cells are never copied from the report:
#' \itemize{
#'   \item recruitment into \code{ln_RecDevs}, which \code{generate_recruitment} reads.
#'   \item catchability into \code{ln_fish_q_devs} or \code{ln_srv_q_devs}, and
#'     \code{\link{draw_sim_q_devs}} scales that fleet's catchability by them whatever the
#'     operating model's own \code{fish_q_model} or \code{srv_q_model} says. The catchability
#'     the deviations scale is the block mean, so the series drives the conditioning years and
#'     every year past them alike.
#'   \item the numbers at age into \code{naa_eta_all}, as the log state less its deterministic
#'     prediction, replacing those cells of the replicate's own innovations.
#'   \item growth and movement into arrays built here, since the operating model has none of
#'     its own: fitted years at their estimates and later years at zero. The fit's own
#'     \code{Get_Growth} or \code{Get_Movement} is then rerun per replicate to rebuild weight
#'     at age, the size-age keys and the movement matrix, reading every argument from
#'     \code{data} and \code{pars} by name. Under selectivity at length a fleet's selectivity
#'     at age is its curve read through the rebuilt key, as in the fit.
#' }
#' Growth propagated cohort by cohort (\code{growth_tv_type = 1}) runs the way the fit runs
#' it: the years before \code{growth_cohort_styr} are built up front, and from there
#' \code{run_annual_cycle} advances each replicate one year at a time from its own start of
#' year numbers. From scratch only recruitment and the numbers at age can be linked, since
#' those are the arrays the operating model has itself.
#'
#' @param sim_list Simulation list, for example from
#'   \code{condition_closed_loop_simulations} or the \code{Setup_Sim_*} calls.
#' @param data Data list of the fitted model, set up by \code{Setup_Mod_DSEM}.
#'   \code{NULL} from scratch.
#' @param pars Parameter list at the fitted values, as \code{get_optim_param_list}
#'   or \code{obj$env$parList()} returns it. \code{NULL} from scratch.
#' @param rep Report of the fit. Needed for a numbers at age series, whose \code{NAA_pred}
#'   turns the fitted log state into innovations, and for a growth series linked under
#'   selectivity at length, which the rebuilt keys are read through.
#' @param dsem_cov_use Optional matrix \code{[year, covariate]} of 0/1 flags for the years
#'   each covariate is observed in the simulation. Defaults to the fitted pattern for fitted
#'   years and every year after, or every year from scratch.
#' @param condition_on_fit \code{TRUE} keeps every fitted year at its estimate and draws
#'   later years given them. \code{FALSE} (default) draws every year.
#' @param dsem_arrows From scratch: the arrow lines, one per element or one per line of a
#'   string, naming covariates and deviation series.
#' @param dsem_values From scratch: named vector giving every arrow parameter its value, sds
#'   on the natural scale, for example
#'   \code{c(b_env = 0.5, rho_env = 0.6, sd_env = 1, sd_rec = 0.8)}.
#' @param dsem_processes From scratch: processes the arrows may name, \code{"rec"} (default)
#'   or \code{"NAA"}, or both.
#' @param dsem_cov_mu From scratch: named vector of covariate means, one per covariate the
#'   arrows name. Default zero for each.
#' @param dsem_cov_obs_sd From scratch: named vector, one per covariate, of the spread of its
#'   observations on its family's own scale: the sd for \code{"normal"}, the CV for
#'   \code{"gamma"}, the sd of the log for \code{"lognormal"}, the dispersion for
#'   \code{"tweedie"}; \code{NA} for a covariate observed without error. Default \code{NA}.
#' @param dsem_cov_family From scratch: named character vector, one per covariate, of the
#'   family the replicates' observations are drawn from, as in \code{\link{Setup_Mod_DSEM}}.
#'   Default \code{"normal"} where \code{dsem_cov_obs_sd} is given and \code{"fixed"}
#'   otherwise. A link-scale family reads \code{dsem_cov_mu} and the arrows on that scale.
#' @param dsem_cov_tweedie_p From scratch: named vector of tweedie powers, strictly between
#'   1 and 2. Default 1.5.
#' @param dsem_cov_link From scratch: named character vector of links, as in
#'   \code{\link{Setup_Mod_DSEM}}. Default each family's own. A \code{"gaussian_fixed_sd"}
#'   covariate from scratch takes its one sd from \code{dsem_cov_obs_sd} for every year.
#' @param mod_var_logscale From scratch: passed to \code{read_dsem_arrows}.
#' @param dsem_variance From scratch: passed to \code{read_dsem_arrows} as \code{variance},
#'   \code{"conditional"} (default) or \code{"diagonal"}. From a fit the data list's own
#'   setting is used.
#'
#' @return \code{sim_list} with the dsem arrows, values and observation pattern added, plus
#'   the deviation arrays and the growth or movement arguments the operating model needs to
#'   consume a drawn series.
#'
#' @seealso \code{\link{Setup_Mod_DSEM}}
#'
#' @export
Setup_Sim_DSEM <- function(sim_list,
                           data = NULL,
                           pars = NULL,
                           rep = NULL,
                           dsem_cov_use = NULL,
                           condition_on_fit = FALSE,
                           dsem_arrows = NULL,
                           dsem_values = NULL,
                           dsem_processes = "rec",
                           dsem_cov_mu = NULL,
                           dsem_cov_obs_sd = NULL,
                           dsem_cov_family = NULL,
                           dsem_cov_tweedie_p = NULL,
                           dsem_cov_link = NULL,
                           mod_var_logscale = FALSE,
                           dsem_variance = "conditional") {

  # From Scratch ------------------------------------------------------------

  # arrows and values stand in for a fit: a data list and parameter list are built from them, with the
  # link table read off the operating model's own arrays, and the rest of the function runs as from a fit
  if(is.null(data)) {

    if(is.null(dsem_arrows) || is.null(dsem_values)) {
      stop("Give data and pars from a fit, or dsem_arrows and dsem_values to write the dsem from scratch.")
    }

    if(condition_on_fit) stop("There is no fit to condition on when the dsem is written from scratch.")

    scratch <- scratch_dsem_fit(sim_list = sim_list,
                                dsem_arrows = dsem_arrows,
                                dsem_values = dsem_values,
                                dsem_processes = dsem_processes,
                                dsem_cov_mu = dsem_cov_mu,
                                dsem_cov_obs_sd = dsem_cov_obs_sd,
                                mod_var_logscale = mod_var_logscale,
                                dsem_variance = dsem_variance,
                                dsem_cov_family = dsem_cov_family,
                                dsem_cov_tweedie_p = dsem_cov_tweedie_p,
                                dsem_cov_link = dsem_cov_link)

    data <- scratch$data
    pars <- scratch$pars

  } # end if written from scratch

  # Options -----------------------------------------------------------------

  n_fit_yrs <- length(data$years) # years the fit covers
  n_sim_yrs <- sim_list$n_yrs # years the operating model runs
  n_cov <- length(data$dsem_cov_var_idx)

  process_table <- dsem_process_table()
  table_par <- vapply(process_table, function(entry) entry$par, "")
  link_entry <- process_table[match(data$dsem_link_par, table_par)] # the process each linked series belongs to

  growth_linked <- any(data$dsem_link_par %in% c("ln_growth_devs", "ln_growth_semipar_devs"))
  move_linked <- any(data$dsem_link_par == "move_devs")
  naa_linked <- any(data$dsem_link_par == "ln_NAA")
  sel_at_length <- c(fish = isTRUE(data$fish_selex_type == 1), # absent from older data lists, so isTRUE
                     ret = isTRUE(data$ret_selex_type == 1),
                     srv = isTRUE(data$srv_selex_type == 1))

  # Input Validation --------------------------------------------------------

  if(is.null(data$dsem_model)) {
    stop("obj was not set up with Setup_Mod_DSEM, so there is no dsem to simulate from.")
  }

  if(sim_list$n_pop != data$n_pop || sim_list$n_regions != data$n_regions) {
    stop("sim_list and obj have different populations or regions.")
  }

  if(condition_on_fit && n_sim_yrs < n_fit_yrs) {
    stop("Conditioning on the fit needs the operating model to run at least the fitted years.")
  }

  if(!is.null(dsem_cov_use) && !identical(dim(dsem_cov_use), c(as.integer(n_sim_yrs), as.integer(n_cov)))) {
    stop("dsem_cov_use should be [", n_sim_yrs, ", ", n_cov,
         "], one row per simulated year and one column per covariate.")
  }

  if(naa_linked && !isTRUE(sim_list$NAA_re > 0)) {
    stop("The dsem links the numbers at age state, but the operating model runs it deterministic ",
         "(NAA_re is off). Condition the simulation with the state on so the drawn cells have ",
         "somewhere to go.")
  }

  if(naa_linked && is.null(rep$NAA_pred) && !isTRUE(data$dsem_from_scratch)) {
    stop("The dsem links the numbers at age state, so Setup_Sim_DSEM needs rep, the fit's report, ",
         "to measure the fitted state against its prediction (ln_NAA - log(NAA_pred)).")
  }

  if(growth_linked && any(sel_at_length) && is.null(rep)) {
    stop("The dsem links a growth series and the fit has selectivity at length, so Setup_Sim_DSEM ",
         "needs rep, the fit's report, for the selectivity at length the rebuilt keys are read through.")
  }

  # Defaults for Unsupplied Inputs ------------------------------------------

  if(is.null(dsem_cov_use)) {

    dsem_cov_use <- matrix(1, n_sim_yrs, n_cov) # observed every year after the fit, and every year from scratch
    fit_rows <- 1:min(n_fit_yrs, n_sim_yrs)

    # the fitted years instead take the pattern the fit was handed
    if(!isTRUE(data$dsem_from_scratch)) {
      dsem_cov_use[fit_rows,] <- as.numeric(!is.na(data$dsem_cov_obs[fit_rows,,drop = FALSE]))
    }

  } # end if no pattern was given

  # Populate Sim List ------------------------------------------------

  sim_list$dsem_model <- data$dsem_model
  sim_list$dsem_beta <- pars$dsem_beta
  sim_list$ln_dsem_sd <- pars$ln_dsem_sd
  sim_list$dsem_mu <- pars$dsem_mu
  sim_list$ln_dsem_obs_sd <- pars$ln_dsem_obs_sd
  sim_list$logit_dsem_tweedie_p <- pars$logit_dsem_tweedie_p # absent from a fit before the tweedie family existed

  # RecDevs_model = 'dsem' reads sigmaR off the arrows' recruitment sd, so the operating model's
  # initial age deviations read the same value
  if("rec" %in% data$dsem_declared && !is.null(sim_list$ln_sigmaR)) {

    arrow_value <- as.numeric(get_dsem_arrow_values(pars$dsem_beta, pars$ln_dsem_sd, data$dsem_model))

    for(s in which(data$dsem_link_par == "ln_RecDevs")) {
      if(!isTRUE(data$dsem_link_sd_arrow[s] > 0)) next # a moderated sd (arrow 0) leaves the list's own sigmaR in place
      pop <- data$dsem_link_idx[[s]][1]
      region <- data$dsem_link_idx[[s]][2]
      sim_list$ln_sigmaR[,pop,region] <- log(arrow_value[data$dsem_link_sd_arrow[s]])
    } # end s loop

  } # end if recruitment is declared

  sim_list$dsem_cov_family <- data$dsem_cov_family
  sim_list$dsem_cov_link <- if(is.null(data$dsem_cov_link)) dsem_default_link(data$dsem_cov_family) else data$dsem_cov_link
  sim_list$dsem_cov_var_idx <- data$dsem_cov_var_idx

  # a known sd per year for the fixed-sd normal, reused past the fit at its last given value
  sim_list$dsem_cov_fixed_sd <- matrix(NA_real_, n_sim_yrs, n_cov)

  if(!is.null(data$dsem_cov_fixed_sd)) {
    for(k in seq_len(n_cov)) {

      given <- data$dsem_cov_fixed_sd[,k] # the sds the fit was handed for this covariate
      if(all(is.na(given))) next # not a fixed-sd covariate

      fit_rows <- 1:min(nrow(data$dsem_cov_fixed_sd), n_sim_yrs)
      sim_list$dsem_cov_fixed_sd[fit_rows,k] <- given[fit_rows]
      sim_list$dsem_cov_fixed_sd[is.na(sim_list$dsem_cov_fixed_sd[,k]),k] <- given[max(which(!is.na(given)))]

    } # end k loop
  } # end if the fit has fixed sds

  # which cells of which arrays the dsem writes into
  for(field_name in c("dsem_link_par", "dsem_link_col", "dsem_link_row", "dsem_link_cell", "dsem_link_idx", "dsem_link_yr_dim")) {
    sim_list[[field_name]] <- data[[field_name]]
  } # end field_name loop

  sim_list$dsem_link_sim_par <- vapply(link_entry, function(entry) entry$sim_par, "") # the array each series is written into
  sim_list$dsem_cov_use <- dsem_cov_use
  sim_list$dsem_n_cond_yrs <- if(condition_on_fit) n_fit_yrs else 0 # leading years kept at the fit

  # What a Conditioned Draw Is Given ----------------------------------------

  # a linked recruitment cell's correction conditions on the covariate years the estimation model is
  # handed, and on the rows before a linked series starts
  x_known <- matrix(FALSE, n_sim_yrs, length(data$dsem_model$variables))
  for(k in seq_len(n_cov)) x_known[dsem_cov_use[,k] == 1, data$dsem_cov_var_idx[k]] <- TRUE
  for(s in seq_along(data$dsem_link_par)) x_known[seq_len(min(data$dsem_link_row[[s]]) - 1), data$dsem_link_col[s]] <- TRUE
  sim_list$dsem_x_known <- x_known

  # the fit says whether its own penalty takes a correction at all
  fit_ramp <- 1
  if(!isTRUE(data$dsem_from_scratch) && !is.null(data$do_rec_bias_ramp)) {
    fit_ramp <- get_rec_bias_ramp(data$do_rec_bias_ramp, data$bias_year, n_fit_yrs, data$max_bias_ramp_fct)
  }
  sim_list$dsem_rec_corr_on <- any(fit_ramp != 0) && !isTRUE(data$RecDevs_model != 1) && !isTRUE(data$RecDevs_pen_center == 1)

  # fitted values of every series, the cells a conditioned draw is given
  x_fit <- matrix(0, n_fit_yrs, length(data$dsem_model$variables))
  for(k in seq_len(n_cov)) x_fit[,data$dsem_cov_var_idx[k]] <- pars$dsem_x[1:n_fit_yrs,data$dsem_cov_var_idx[k]]

  for(s in seq_along(data$dsem_link_par)) {

    fit_rows <- data$dsem_link_row[[s]][data$dsem_link_row[[s]] <= n_fit_yrs]
    cells <- data$dsem_link_cell[[s]][seq_along(fit_rows)]
    value <- pars[[data$dsem_link_par[s]]][cells]

    # the arrows describe the state less its prediction, not the state
    if(data$dsem_link_par[s] == "ln_NAA" && !is.null(rep$NAA_pred)) value <- value - log(rep$NAA_pred[cells])

    x_fit[fit_rows,data$dsem_link_col[s]] <- value

  } # end s loop

  sim_list$dsem_x_fit <- x_fit

  # Deviation Arrays the Operating Model Lacks ------------------------------

  # growth reads both of its arrays, so both come across when either is linked. fitted years keep
  # the fit's values, years past what the fit holds sit at zero, and every replicate starts the same
  build_pars <- c(if(growth_linked) c("ln_growth_devs", "ln_growth_semipar_devs"), if(move_linked) "move_devs")

  for(par_name in build_pars) {

    if(!is.null(sim_list[[par_name]]) || is.null(pars[[par_name]])) next # already there, or the fit has none
    dev_array <- pars[[par_name]]
    yr_dim <- which(process_table[[match(par_name, table_par)]]$dim_names == "Yr")
    n_fit_dev_yrs <- dim(dev_array)[yr_dim]

    # cut the years the operating model does not run
    if(n_fit_dev_yrs > n_sim_yrs) {
      keep_levels <- lapply(dim(dev_array), seq_len)
      keep_levels[[yr_dim]] <- seq_len(n_sim_yrs)
      dev_array <- do.call("[", c(list(dev_array), keep_levels, list(drop = FALSE)))
    }

    if(n_fit_dev_yrs < n_sim_yrs) dev_array <- extend_years(dev_array, n_sim_yrs - n_fit_dev_yrs, yr_dim, "zeros")
    sim_list[[par_name]] <- array(dev_array, dim = c(dim(dev_array), sim_list$n_sims)) # replicate dim last

  } # end par_name loop

  # the fit's own growth function, with every argument read from data and pars by name
  if(growth_linked) sim_list$dsem_growth_args <- match_model_args(Get_Growth, data, pars, n_yrs = n_sim_yrs)

  # what the fit does after growth needs the selectivity at length the report holds, run out to the
  # operating model's years the way the closed loop extends its other inputs
  if(growth_linked) {

    sel_len <- list()

    for(sel_name in names(sel_at_length)[sel_at_length]) {

      sel_array <- rep[[paste0(sel_name, "_sel_l")]] # [region, year, len, sex, fleet]
      n_rep_yrs <- dim(sel_array)[2]
      if(n_rep_yrs > n_sim_yrs) sel_array <- sel_array[,seq_len(n_sim_yrs),,,,drop = FALSE]
      if(n_rep_yrs < n_sim_yrs) sel_array <- extend_years(sel_array, n_sim_yrs - n_rep_yrs, 2, "last")
      sel_len[[sel_name]] <- sel_array

    } # end sel_name loop

    # and which fleets weigh their catch by what they select
    sim_list$dsem_length_sel <- list(fish_selex_type = as.integer(sel_at_length[["fish"]]),
                                     ret_selex_type = as.integer(sel_at_length[["ret"]]),
                                     srv_selex_type = as.integer(sel_at_length[["srv"]]),
                                     fish_waa_selected = if(is.null(data$fish_waa_selected)) 0 else data$fish_waa_selected,
                                     srv_waa_selected = if(is.null(data$srv_waa_selected)) 0 else data$srv_waa_selected,
                                     fish_sel_l = sel_len$fish,
                                     ret_sel_l = sel_len$ret,
                                     srv_sel_l = sel_len$srv)

  } # end if a growth series is linked

  # movement keeps the fit's year count so covariates stop at the data, as they do in a projection
  if(move_linked) {
    sim_list$dsem_move_args <- match_model_args(Get_Movement, data, pars,
                                               n_yrs = min(n_fit_yrs, n_sim_yrs),
                                               n_proj_yrs_devs = max(0, n_sim_yrs - n_fit_yrs),
                                               n_ages = length(data$ages),
                                               expm_nsub = if(is.null(data$move_expm_nsub)) 0 else data$move_expm_nsub)
  }

  # whether or not q dsem can be drawn ... (since can have 0 sd in q series)
  check_q_dsem_drawable(sim_list)

  return(sim_list)

} # end function

#' A data list and parameter list for a dsem written from scratch
#'
#' Stands in for a fit when an operating model states its own dsem: the arrows are read
#' against the operating model's arrays (every cell estimated, no map), and the values given
#' become the parameters. Only recruitment and the numbers at age can be named, since those
#' are the arrays the operating model has.
#'
#' @param sim_list Simulation list holding the dims and the arrays.
#' @param dsem_arrows,dsem_values,dsem_processes,dsem_cov_mu,dsem_cov_obs_sd,mod_var_logscale,dsem_variance,dsem_cov_family,dsem_cov_tweedie_p,dsem_cov_link
#'   As in \code{\link{Setup_Sim_DSEM}}.
#'
#' @return List with \code{data} and \code{pars} holding the dsem fields
#'   \code{Setup_Sim_DSEM} reads from a fit.
#'
#' @keywords internal
scratch_dsem_fit <- function(sim_list,
                             dsem_arrows,
                             dsem_values,
                             dsem_processes,
                             dsem_cov_mu,
                             dsem_cov_obs_sd,
                             mod_var_logscale,
                             dsem_variance = "conditional",
                             dsem_cov_family = NULL,
                             dsem_cov_tweedie_p = NULL,
                             dsem_cov_link = NULL) {

  # The Operating Model's Arrays, Read as a Fit's ---------------------------

  if(!all(dsem_processes %in% c("rec", "NAA"))) {
    stop("From scratch the arrows can name recruitment and the numbers at age only; ",
         "growth and movement need a fit's data and pars.")
  }

  n_yrs <- sim_list$n_yrs

  # the arrays without their replicate dim, laid out as a fit's input_list so get_dsem_link can read them
  fit_input <- list(data = list(n_pop = sim_list$n_pop,
                                n_regions = sim_list$n_regions,
                                years = seq_len(n_yrs),
                                n_proj_yrs_devs = 0),
                    par = list(),
                    map = list()) # no map, so every cell is estimated

  fit_input$par$ln_RecDevs <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, n_yrs))

  if("NAA" %in% dsem_processes) {

    if(!isTRUE(sim_list$NAA_re > 0)) {
      stop("The arrows name the numbers at age state, but the operating model runs it deterministic. ",
           "Call Setup_Sim_NAA_state first.")
    }

    fit_input$par$ln_NAA <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, n_yrs,
                                             sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))
    fit_input$data$naa_re_yrs <- sim_list$naa_re_yrs # the state's own years, seasons and ages
    fit_input$data$naa_re_seas <- sim_list$naa_re_seas
    fit_input$data$naa_re_ages <- sim_list$naa_re_ages

  } # end if the state is named

  dsem_arrows <- paste(dsem_arrows, collapse = "\n") # one string with a line per arrow, or one arrow per element
  arrow_text <- paste(sub("#.*$", "", strsplit(dsem_arrows, "\n")[[1]]), collapse = " ") # comments off each line
  link <- get_dsem_link(fit_input, dsem_processes, arrow_text, n_yrs)

  if(length(link$name) == 0) {
    stop("No deviation series appears in dsem_arrows. The series available are: ",
         paste(utils::head(link$offered, 20), collapse = ", "), ".")
  }

  # Covariates the Arrows Name ----------------------------------------------

  # every name in the arrows that is not a series is a covariate, given a mean and, when asked, a measurement sd
  tokens <- unique(unlist(regmatches(arrow_text, gregexpr("[A-Za-z][A-Za-z0-9_.]*", arrow_text))))
  cov_names <- if(is.null(dsem_cov_mu)) setdiff(tokens, c(link$offered, "NA", names(dsem_values))) else names(dsem_cov_mu)

  if(is.null(dsem_cov_mu)) dsem_cov_mu <- stats::setNames(rep(0, length(cov_names)), cov_names)
  if(is.null(dsem_cov_obs_sd)) dsem_cov_obs_sd <- stats::setNames(rep(NA_real_, length(cov_names)), cov_names)

  if(!all(cov_names %in% names(dsem_cov_obs_sd))) {
    stop("dsem_cov_obs_sd needs one entry per covariate: ", paste(cov_names, collapse = ", "))
  }

  variables <- c(cov_names, link$name) # the grid's columns, covariates first
  dsem_model <- read_dsem_arrows(dsem_arrows, variables, mod_var_logscale = mod_var_logscale, variance = dsem_variance)

  # Values for Every Arrow --------------------------------------------------

  par_names <- c(dsem_model$beta_names, dsem_model$ln_sd_names)

  if(!all(par_names %in% names(dsem_values))) {
    stop("dsem_values needs every arrow parameter: ", paste(setdiff(par_names, names(dsem_values)), collapse = ", "))
  }

  if(any(dsem_values[dsem_model$ln_sd_names] <= 0)) {
    stop("An sd in dsem_values must be positive; sds are given on the natural scale.")
  }

  # the family each covariate is observed through: without one, an sd means normal error and NA means known
  if(is.null(dsem_cov_family)) {
    dsem_cov_family <- stats::setNames(ifelse(is.na(dsem_cov_obs_sd[cov_names]), "fixed", "normal"), cov_names)
  }

  family_code <- unname(dsem_family_codes()[dsem_cov_family[cov_names]])

  if(any(is.na(family_code))) {
    stop("dsem_cov_family should name every covariate with one of '",
         paste(names(dsem_family_codes()), collapse = "', '"), "'.")
  }

  needs_sd <- family_code %in% c(1, 4, 5, 6, 7) & is.na(dsem_cov_obs_sd[cov_names]) # every family but fixed, binomial and poisson

  if(any(needs_sd)) {
    stop("dsem_cov_obs_sd is needed for a normal (sd), gamma (CV), gaussian_fixed_sd (sd), ",
         "lognormal (sd of the log) or tweedie (dispersion) covariate: ",
         paste(cov_names[needs_sd], collapse = ", "))
  }

  link_code <- if(is.null(dsem_cov_link)) dsem_default_link(family_code) else unname(dsem_link_codes()[dsem_cov_link[cov_names]])

  if(any(is.na(link_code))) {
    stop("dsem_cov_link should name every covariate with one of '",
         paste(names(dsem_link_codes()), collapse = "', '"), "'.")
  }

  if(is.null(dsem_cov_tweedie_p)) dsem_cov_tweedie_p <- stats::setNames(rep(1.5, length(cov_names)), cov_names)
  tweedie_p <- ifelse(is.na(dsem_cov_tweedie_p[cov_names]), 1.5, dsem_cov_tweedie_p[cov_names])

  if(any(family_code == 7 & (tweedie_p <= 1 | tweedie_p >= 2))) {
    stop("dsem_cov_tweedie_p needs a power strictly between 1 and 2 for every tweedie covariate.")
  }

  # Populate the Data List --------------------------------------------------

  link_col <- match(link$name, variables) # the grid column each linked series sits in

  # the fixed-sd normal takes one sd in every year, and no other family has one
  cov_fixed_sd <- matrix(NA_real_, n_yrs, length(cov_names))
  for(k in seq_along(cov_names)) if(family_code[k] == 5) cov_fixed_sd[,k] <- dsem_cov_obs_sd[cov_names[k]]

  data <- list(years = seq_len(n_yrs),
               n_pop = sim_list$n_pop,
               n_regions = sim_list$n_regions,
               dsem_from_scratch = TRUE, # read wherever a fitted quantity would otherwise be looked for
               dsem_model = dsem_model,
               dsem_n_grid_yrs = n_yrs,
               dsem_var_names = variables,
               dsem_cov_obs = matrix(NA_real_, n_yrs, length(cov_names), dimnames = list(NULL, cov_names)),
               dsem_cov_var_idx = match(cov_names, variables),
               dsem_cov_family = family_code,
               dsem_cov_link = link_code,
               dsem_cov_fixed_sd = cov_fixed_sd,
               dsem_link_par = link$par,
               dsem_link_col = link_col,
               dsem_link_row = link$grid_row,
               dsem_link_cell = link$cell,
               dsem_link_idx = link$idx,
               dsem_link_yr_dim = link$yr_dim,
               dsem_link_sd_arrow = get_dsem_link_sd_arrow(dsem_model, link_col),
               dsem_declared = unique(link$label))

  # Populate the Parameter List ---------------------------------------------

  mu <- rep(0, length(variables)) # a deviation series has mean zero, a covariate the mean given
  mu[match(cov_names, variables)] <- dsem_cov_mu[cov_names]

  pars <- list(dsem_beta = unname(dsem_values[dsem_model$beta_names]),
               ln_dsem_sd = unname(log(dsem_values[dsem_model$ln_sd_names])), # given on the natural scale
               dsem_mu = mu,
               dsem_x = matrix(rep(mu, each = n_yrs), n_yrs, length(variables)), # the grid at its means
               ln_dsem_obs_sd = unname(ifelse(is.na(dsem_cov_obs_sd[cov_names]), 0, log(dsem_cov_obs_sd[cov_names]))),
               logit_dsem_tweedie_p = unname(stats::qlogis(tweedie_p - 1)))

  pars$ln_RecDevs <- fit_input$par$ln_RecDevs
  if(!is.null(fit_input$par$ln_NAA)) pars$ln_NAA <- fit_input$par$ln_NAA

  return(list(data = data, pars = pars))

} # end function

#' Draw the dsem grid year by year
#'
#' A moderated arrow puts the grid's own values inside \eqn{B}, so there is no one covariance
#' to draw from. Each year is drawn instead in the order \code{series_order} gives, every cell
#' from what points into it and an innovation, which is the model the density evaluates. Years
#' already known are held as is.
#'
#' @param dsem_model Output of \code{read_dsem_arrows}.
#' @param arrow_value Values of the arrows, from \code{get_dsem_arrow_values}.
#' @param mu_grid Matrix \code{[year, series]} of series means.
#' @param n_sims Number of replicates.
#' @param n_cond Leading years kept at \code{x_known} rather than drawn.
#' @param x_known Matrix \code{[year, series]} read for those years.
#'
#' @return Array \code{[year, series, sim]}.
#'
#' @keywords internal
draw_dsem_recursive <- function(dsem_model,
                                arrow_value,
                                mu_grid,
                                n_sims,
                                n_cond = 0,
                                x_known = NULL) {

  arrows <- dsem_model$arrows
  n_yrs <- nrow(mu_grid)
  n_vars <- length(dsem_model$variables)

  if(is.null(dsem_model$series_order)) stop("Same-year paths form a loop, so there is no order to draw the series in.")
  if(any(arrows$lag < 0)) stop("An arrow reads a later year, so the grid cannot be drawn forward in time. Draw from the precision instead.")

  x <- array(0, dim = c(n_yrs, n_vars, n_sims))
  if(n_cond > 0) for(sim in 1:n_sims) x[1:n_cond,,sim] <- x_known[1:n_cond,]

  for(yr in seq_len(n_yrs)) {

    if(yr <= n_cond) next # kept at the values given

    for(j in dsem_model$series_order) {

      value <- rep(mu_grid[yr,j], n_sims)

      # each from series' deviation from its own mean, times this arrow's value that year
      for(i in which(arrows$type == "path" & arrows$to_idx == j)) {
        if(yr - arrows$lag[i] < 1) next # no earlier value to read
        coef <- if(arrows$mod_idx[i] > 0) x[yr,arrows$mod_idx[i],] else arrow_value[i] # a moderated path reads the grid
        value <- value + coef * (x[yr - arrows$lag[i],arrows$from_idx[i],] - mu_grid[yr - arrows$lag[i],arrows$from_idx[i]])
      } # end i loop

      i_sd <- which(arrows$type == "sd" & arrows$to_idx == j) # this series' innovation sd

      # a moderated sd line reads its sd off the grid too, on the log scale when asked
      if(arrows$mod_idx[i_sd] > 0) {
        sd_yr <- x[yr,arrows$mod_idx[i_sd],]
        if(isTRUE(dsem_model$mod_var_logscale)) sd_yr <- exp(sd_yr)
      } else sd_yr <- arrow_value[i_sd]

      x[yr,j,] <- value + stats::rnorm(n_sims, 0, abs(sd_yr))

    } # end j loop
  } # end yr loop

  return(x)

} # end function

#' Draw dsem series for every replicate
#'
#' Draws the year by series grid for all replicates up front, which is the same in
#' distribution as drawing year by year because the series do not depend on the harvest, and
#' keeps the same draws across management procedures. A linked recruitment cell is drawn about
#' minus half its variance given the known cells (\code{\link{get_dsem_margvar}}, with a
#' moderating series at its mean) when \code{rec_bias_correct} is on.
#'
#' @section What it fills:
#' \code{dsem_x_sim} \code{[year, series, sim]}, every linked cell of every linked array from
#' it with \code{dsem_drawn} marking those cells, and \code{dsem_cov_obs_sim} \code{[year,
#' covariate, sim]} (\code{NA} where not observed). Growth and movement are then rebuilt from
#' the arrays a drawn series went into.
#'
#' @param sim_env Simulation environment from \code{Setup_sim_env}, built from a list set up
#'   by \code{Setup_Sim_DSEM}.
#'
#' @return \code{invisible(NULL)}; \code{sim_env} is modified in place.
#'
#' @keywords internal
draw_dsem_sim <- function(sim_env) {

  n_yrs <- sim_env$n_yrs
  n_sims <- sim_env$n_sims
  n_vars <- length(sim_env$dsem_model$variables)
  n_cov <- length(sim_env$dsem_cov_var_idx)
  n_cond <- sim_env$dsem_n_cond_yrs

  # Means Over the Grid -----------------------------------------------------

  # zero for a linked series, since it is a deviation (for the state series, its innovation)
  mu_grid <- matrix(0, n_yrs, n_vars)
  for(k in seq_len(n_cov)) mu_grid[,sim_env$dsem_cov_var_idx[k]] <- sim_env$dsem_mu[sim_env$dsem_cov_var_idx[k]]

  # a linked recruitment cell is drawn about minus half its variance given the known cells, the full
  # correction the density gives a random effect, when the fit's penalty takes a correction at all
  rec_col <- sim_env$dsem_link_col[sim_env$dsem_link_par == "ln_RecDevs"]

  if(length(rec_col) > 0 && isTRUE(sim_env$rec_bias_correct == 1) && isTRUE(sim_env$dsem_rec_corr_on)) {

    margvar <- get_dsem_margvar(sim_env$dsem_beta,
                                sim_env$ln_dsem_sd,
                                mu_grid,
                                sim_env$dsem_model,
                                get_dsem_cells(sim_env$dsem_model, n_yrs),
                                as.vector(sim_env$dsem_x_known))

    mu_grid[,rec_col] <- mu_grid[,rec_col] - 0.5 * as.matrix(margvar[,rec_col,drop = FALSE])

  } # end if the correction is on

  # Draw Every Replicate ----------------------------------------------------

  x_known <- matrix(0, n_yrs, n_vars) # the fitted years, when conditioning on them
  if(n_cond > 0) x_known[1:n_cond,] <- sim_env$dsem_x_fit[1:n_cond,]

  if(any(sim_env$dsem_model$arrows$mod_idx > 0)) {

    # a moderated arrow holds the grid's own values, so the field has no one covariance to draw from
    arrow_value <- get_dsem_arrow_values(sim_env$dsem_beta, sim_env$ln_dsem_sd, sim_env$dsem_model)
    x_sim <- draw_dsem_recursive(sim_env$dsem_model, arrow_value, mu_grid, n_sims, n_cond, x_known)

  } else {

    Q <- get_dsem_precision(sim_env$dsem_beta, sim_env$ln_dsem_sd, sim_env$dsem_model, get_dsem_cells(sim_env$dsem_model, n_yrs))

    known <- matrix(FALSE, n_yrs, n_vars)
    if(n_cond > 0) known[1:n_cond,] <- TRUE
    known_cell <- which(as.vector(known)) # as.vector stacks series, matching the cell numbering

    dsem_cond <- get_dsem_conditional(Q, as.vector(mu_grid), known_cell, as.vector(x_known)[known_cell])

    x_sim <- array(as.vector(x_known), dim = c(n_yrs * n_vars, n_sims))
    x_sim[dsem_cond$unknown_cell,] <- draw_dsem_conditional(dsem_cond, n_sims)
    x_sim <- array(x_sim, dim = c(n_yrs, n_vars, n_sims))

  } # end if any moderated arrow

  sim_env$dsem_x_sim <- x_sim

  # Write Each Series Into Its Array ----------------------------------------

  # straight into the array the operating model reads, the way the objective gathers them. that array
  # has the replicate dim last, so each cell is rebuilt at its dims, and dsem_drawn marks it
  for(s in seq_along(sim_env$dsem_link_par)) {

    par_name <- sim_env$dsem_link_sim_par[s]

    if(is.null(sim_env[[par_name]])) {
      stop("The dsem links ", sim_env$dsem_link_par[s], " but the operating model holds no ",
           par_name, " array to write it into.")
    }

    dims <- dim(sim_env[[par_name]])
    fit_dims <- dims[-length(dims)] # everything but the replicate
    stride <- cumprod(c(1, fit_dims[-length(fit_dims)])) # step in the flattened array per dim
    yr_dim <- sim_env$dsem_link_yr_dim[s]
    idx <- sim_env$dsem_link_idx[[s]]

    grid_row <- sim_env$dsem_link_row[[s]]
    grid_row <- c(grid_row, if(max(grid_row) < n_yrs) (max(grid_row) + 1):n_yrs) # and on through the simulated years
    grid_row <- grid_row[grid_row <= fit_dims[yr_dim]]

    for(y in grid_row) {
      idx[yr_dim] <- y
      cell <- 1L + sum((idx - 1L) * stride) # this year's cell in the first replicate
      sim_env[[par_name]][cell + (seq_len(n_sims) - 1L) * prod(fit_dims)] <- x_sim[y,sim_env$dsem_link_col[s],]
      sim_env$dsem_drawn[[par_name]][cell] <- TRUE
    } # end y loop

  } # end s loop

  # Covariate Observations --------------------------------------------------

  # the state itself when known, otherwise a draw from the family about the link-scale state
  cov_obs <- array(NA_real_, dim = c(n_yrs, n_cov, n_sims))

  for(k in seq_len(n_cov)) {

    state <- x_sim[,sim_env$dsem_cov_var_idx[k],,drop = FALSE] # [year, 1, sim]
    tweedie_p <- if(is.null(sim_env$logit_dsem_tweedie_p)) 1.5 else 1 + stats::plogis(sim_env$logit_dsem_tweedie_p[k])
    link <- if(is.null(sim_env$dsem_cov_link)) dsem_default_link(sim_env$dsem_cov_family[k]) else sim_env$dsem_cov_link[k]
    fixed_sd <- if(is.null(sim_env$dsem_cov_fixed_sd)) NULL else sim_env$dsem_cov_fixed_sd[,k]

    observed <- draw_dsem_cov_obs(state, sim_env$dsem_cov_family[k], link, exp(sim_env$ln_dsem_obs_sd[k]), tweedie_p, fixed_sd)
    observed[sim_env$dsem_cov_use[,k] == 0,,] <- NA # years this covariate is not observed
    cov_obs[,k,] <- observed[,1,]

  } # end k loop

  sim_env$dsem_cov_obs_sim <- cov_obs

  # the population loop reads weight at age and movement, not their deviations, so a drawn growth
  # or movement series is turned into those here by the fit's own functions
  if(!is.null(sim_env$dsem_growth_args)) derive_sim_growth(sim_env)
  if(!is.null(sim_env$dsem_move_args)) derive_sim_movement(sim_env)

  return(invisible(NULL))

} # end function

#' Conditional distribution of unknown cells in a Gaussian Markov random field
#'
#' For \eqn{x \sim N(\mu, Q^{-1})} split into unknown cells \eqn{u} and known cells \eqn{k},
#' \eqn{x_u \mid x_k \sim N(\mu_u - Q_{uu}^{-1} Q_{uk} (x_k - \mu_k), Q_{uu}^{-1})}. Only the
#' unknown block of the precision is factored.
#'
#' @param Q Sparse precision matrix.
#' @param mu_cell Mean of every cell.
#' @param known_cell Indices of the known cells.
#' @param x_known Values of the known cells, in the order of \code{known_cell}.
#'
#' @return List with \code{unknown_cell}, \code{cond_mean} and \code{chol_uu}, the sparse
#'   Cholesky factor of \eqn{Q_{uu}}.
#'
#' @keywords internal
get_dsem_conditional <- function(Q,
                                 mu_cell,
                                 known_cell,
                                 x_known) {

  unknown_cell <- setdiff(seq_len(nrow(Q)), known_cell)
  Q_uu <- Matrix::forceSymmetric(Q[unknown_cell, unknown_cell, drop = FALSE]) # precision of the unknown cells given the known
  chol_uu <- Matrix::Cholesky(Q_uu, LDL = FALSE, perm = TRUE) # P Q_uu t(P) = L t(L)

  cond_mean <- mu_cell[unknown_cell] # with nothing known the draw is the unconditional field about its mean

  if(length(known_cell) > 0) {
    Q_uk <- Q[unknown_cell, known_cell, drop = FALSE] # how the known cells pull on the unknown
    shift <- Matrix::solve(chol_uu, Q_uk %*% (x_known - mu_cell[known_cell]), system = "A") # Q_uu^-1 Q_uk (x_k - mu_k)
    cond_mean <- cond_mean - as.vector(shift)
  } # end if any cell is known

  return(list(unknown_cell = unknown_cell, cond_mean = cond_mean, chol_uu = chol_uu))

} # end function

#' Draw unknown cells from their conditional distribution
#'
#' @param dsem_cond Output of \code{get_dsem_conditional}.
#' @param n_sims Number of draws.
#'
#' @return Matrix \code{[unknown cell, sim]}.
#'
#' @keywords internal
draw_dsem_conditional <- function(dsem_cond,
                                  n_sims) {

  n_unknown <- length(dsem_cond$cond_mean)
  z <- matrix(stats::rnorm(n_unknown * n_sims), n_unknown, n_sims) # independent standard normals
  w <- Matrix::solve(dsem_cond$chol_uu, z, system = "Lt") # t(L)^-1 z
  dev <- Matrix::solve(dsem_cond$chol_uu, w, system = "Pt") # t(P) t(L)^-1 z, covariance Q_uu^-1

  return(dsem_cond$cond_mean + as.matrix(dev))

} # end function

#' Draw covariate observations from their family and link
#'
#' One draw per cell: the state through the link is the mean, and the family draws about it, as
#' \code{\link{get_dsem_obs_nLL}} evaluates it. Fixed returns the state itself.
#'
#' @param state Array \code{[year, 1, sim]} of link-scale cells.
#' @param family Family code as in \code{\link{Setup_Mod_DSEM}}.
#' @param link Link code.
#' @param obs_sd Spread parameter on the family's own scale.
#' @param tweedie_p Tweedie power in (1, 2).
#' @param fixed_sd Known sd per year for gaussian_fixed_sd, recycled over replicates;
#'   NULL otherwise.
#'
#' @return Array shaped like \code{state}.
#'
#' @keywords internal
draw_dsem_cov_obs <- function(state,
                              family,
                              link,
                              obs_sd,
                              tweedie_p,
                              fixed_sd = NULL) {

  n <- length(state)

  # the mean on the covariate's own scale: identity unless a link was asked for
  mu <- switch(as.character(link),
               "1" = exp(state), # log
               "2" = stats::plogis(state), # logit
               "3" = 1 - exp(-exp(state)), # complementary log log
               state)

  drawn <- switch(as.character(family),
                  "0" = state, # fixed, so the state is the observation
                  "1" = mu + stats::rnorm(n, 0, obs_sd), # normal
                  "2" = stats::rbinom(n, 1, mu), # bernoulli
                  "3" = stats::rpois(n, mu), # poisson
                  "4" = stats::rgamma(n, shape = 1 / obs_sd^2, scale = mu * obs_sd^2), # gamma, obs_sd the CV
                  "5" = mu + stats::rnorm(n, 0, rep(fixed_sd, length.out = n)), # normal with a known sd per year
                  "6" = exp(log(mu) + stats::rnorm(n, 0, obs_sd)), # lognormal, obs_sd the sd of the log
                  "7" = draw_tweedie(mu, obs_sd, tweedie_p), # tweedie, obs_sd the dispersion
                  stop("Unknown dsem covariate family code ", family))

  return(array(as.numeric(drawn), dim = dim(state)))

} # end function

#' Draw from a tweedie as a poisson sum of gammas
#'
#' A tweedie with power between 1 and 2 is a poisson number of gamma jumps: the count has rate
#' \eqn{\mu^{2-p} / (\phi (2-p))}, each jump shape \eqn{(2-p)/(p-1)} and scale
#' \eqn{\phi (p-1) \mu^{p-1}}, and no jumps is a zero. Mean \eqn{\mu}, variance
#' \eqn{\phi \mu^p}.
#'
#' @param mu Means, one per draw.
#' @param phi Dispersion.
#' @param p Power in (1, 2).
#'
#' @return Numeric vector, one draw per mean.
#'
#' @keywords internal
draw_tweedie <- function(mu,
                         phi,
                         p) {

  mu <- as.numeric(mu)
  n_jump <- stats::rpois(length(mu), mu^(2 - p) / (phi * (2 - p))) # how many jumps this cell gets
  drawn <- stats::rgamma(length(mu), shape = n_jump * (2 - p) / (p - 1), scale = phi * (p - 1) * mu^(p - 1))
  drawn[n_jump == 0] <- 0 # no jumps is an exact zero, which is the mass a tweedie puts there

  return(drawn)

} # end function
