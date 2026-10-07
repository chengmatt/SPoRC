# Operating model
#
# Self testing: simulate from known parameters, fit back, and compare. simulation_data_to_SPoRC
# turns operating model output into the data and par lists fit_model expects.

#' Warn when the estimation and operating models start from different R0
#'
#' The operating model takes \code{R0_input} by year, and year one of it does
#' two jobs: it is the equilibrium \code{\link{generate_initial_age_structure}}
#' solves from, and it is the R0 that generates year one's recruitment. So the
#' operating model necessarily starts from the block in force at year one.
#'
#' The estimation model instead builds its initial age structure from
#' \code{R0[R0_ref_block]}. The two agree whenever \code{R0_ref_block} is the
#' block covering year one, which is the default and the case that makes
#' physical sense, since the equilibrium the series starts from is the one the
#' first block describes. Under any other \code{R0_ref_block} the estimation
#' model cannot reproduce the operating model's initial numbers no matter how
#' well it fits, so a self-test would be measuring that gap rather than the
#' feature under test.
#'
#' Only matters under \code{use_rinit = 0}; with \code{use_rinit = 1} both
#' sides initialize from \code{rinit} and the blocks never enter.
#'
#' @param data The model data list.
#' @param where Name of the calling routine, for the message.
#'
#' @return \code{NULL}, invisibly. Called for the warning.
#'
#' @keywords internal
warn_R0_ref_block_om <- function(data, where) {
  if(is.null(data$R0_blocks) || isTRUE(data$use_rinit == 1)) return(invisible(NULL))
  ref <- if(is.null(data$R0_ref_block)) 1 else data$R0_ref_block
  yr1 <- unique(as.vector(data$R0_blocks[, 1, , drop = FALSE]))
  if(all(yr1 == ref)) return(invisible(NULL))
  warning(where, ": R0_ref_block is ", ref, " but year one sits in block ",
          paste(yr1, collapse = "/"), ". The operating model starts from the block in force at year one, ",
          "since that R0 both solves its equilibrium and generates its first year's recruitment, while ",
          "the estimation model initializes from the reference block. The two therefore start from ",
          "different numbers at age and the fit cannot recover the operating model. Set R0_ref_block to ",
          "the block covering year one, which is the default, or use_rinit = 1.")
  invisible(NULL)
}

#' Stack one array per replicate
#'
#' The operating model reads every input with the replicates on the last dim.
#'
#' @param parts List of arrays, one per replicate, all the same shape.
#'
#' @return One array, replicates on the last dim.
#'
#' @keywords internal
bind_sims <- function(parts) {
  d <- dim(parts[[1]])
  if(is.null(d)) d <- length(parts[[1]])
  array(unlist(parts), dim = c(d, length(parts)))
}

#' Make the data effectively exact
#'
#' Every observation gets a standard deviation of \code{obs_sd} and every composition
#' a sample size of \code{iss}. Process error is left alone, being part of the truth.
#' That covers the population-specific at-age sources, the estimated part of an
#' index sd, and a multivariate normal index, whose covariance keeps its
#' correlations at a marginal sd of \code{obs_sd}. Tag recaptures are counts and
#' keep their error.
#'
#' @param sim_list The simulation list, once the setup routines have filled it.
#' @param obs_sd Observation standard deviation to impose. Default \code{1e-3}.
#' @param iss Input sample size to impose on every composition. Default \code{1e6}.
#'
#' @return \code{sim_list} with its observation error replaced.
#'
#' @keywords internal
make_data_perfect <- function(sim_list, obs_sd = 1e-3, iss = 1e6) {

  ln_sigma_names <- c("ln_sigmaC", "ln_sigmaC_pop", "ln_sigmaD", "ln_sigmaD_pop",
                      "ln_sigmaCAA", "ln_sigmaDAA", "ln_sigmaSrvIdxAA",
                      "ln_sigmaCAA_pop", "ln_sigmaDAA_pop", "ln_sigmaSrvIdxAA_pop",
                      "ln_sigmaFishIdx", "ln_sigmaFishIdx_pop", "ln_sigmaSrvIdx", "ln_sigmaSrvIdx_pop") # the last four are an index sd's estimated part
  se_names <- c("ObsFishIdx_SE", "ObsFishIdx_pop_SE", "ObsSrvIdx_SE", "ObsSrvIdx_pop_SE",
                "ObsCatchAA_SE", "ObsDiscardAA_SE", "ObsSrvIdxAA_SE",
                "ObsCatchAA_pop_SE", "ObsDiscardAA_pop_SE", "ObsSrvIdxAA_pop_SE", "dsem_cov_obs_sd")
  iss_names <- grep("^ISS_", names(sim_list), value = TRUE)

  # catch, discards and the at-age sources take their sd on the log scale
  for(ln_sigma_name in ln_sigma_names) {
    if(is.null(sim_list[[ln_sigma_name]])) next
    sim_list[[ln_sigma_name]][] <- log(obs_sd)
  } # end ln_sigma_name loop

  # the indices take theirs on the natural scale, and a year with no survey stays NA
  for(se_name in se_names) {
    if(is.null(sim_list[[se_name]])) next
    survey_se <- sim_list[[se_name]]
    survey_se[!is.na(survey_se)] <- obs_sd
    sim_list[[se_name]] <- survey_se
  } # end se_name loop

  # a multivariate normal index keeps its correlations, each cell at a marginal sd of obs_sd
  for(mvn_name in c("fish_idx_mvn", "srv_idx_mvn")) {
    for(f in seq_along(sim_list[[mvn_name]])) {
      mvn <- sim_list[[mvn_name]][[f]]
      if(is.null(mvn)) next
      if(!is.null(mvn$chol_lower)) mvn$chol_lower <- obs_sd * mvn$chol_lower / mvn$d # rows of the factor over each cell's sd
      mvn$d[] <- obs_sd
      mvn$d_mean <- obs_sd
      sim_list[[mvn_name]][[f]] <- mvn
    } # end f loop
  } # end mvn_name loop

  # a composition with no fish aged stays at zero, so only the ones with a sample change
  for(iss_name in iss_names) {
    if(is.null(sim_list[[iss_name]])) next
    sample_size <- sim_list[[iss_name]]
    sample_size[sample_size > 0] <- iss
    sim_list[[iss_name]] <- sample_size
  } # end iss_name loop

  sim_list

}

#' Give failed replicates the shape of the ones that refit
#'
#' A replicate whose refit fails is left empty or a single \code{NA}, and
#' \code{simplify2array} returns a list rather than an array as soon as one
#' replicate does. Each failed replicate becomes \code{NA} in the shape of a
#' replicate that refit, so the results keep their replicate dim and arithmetic
#' against the truth runs, \code{NA} where a refit failed.
#'
#' @param x List of one result per replicate.
#'
#' @return \code{x}, failed replicates filled.
#'
#' @keywords internal
fill_failed_replicates <- function(x) {

  failed <- vapply(x, function(res) length(res) == 0 || (length(res) == 1 && is.na(res)), logical(1))
  if(!any(failed) || all(failed)) return(x)
  template <- x[[which(!failed)[1]]]
  for(i in which(failed)) x[[i]] <- array(NA, dim = if(is.null(dim(template))) length(template) else dim(template))
  x

} # end fill_failed_replicates

#' The priors a self test's refits read
#'
#' A prior's mean is data the refit reads, but a self test does not redraw it with the
#' observations, so where the data pulled the fit away from a prior every replicate is pulled
#' back toward it, a bias of the design rather than of the estimator. \code{"assessment"} keeps
#' the priors as the assessment gives them; \code{"truth"} moves the mean of every catchability,
#' natural mortality, steepness, R0 and selectivity prior to the value the replicate ran on,
#' keeping its sd, so the prior keeps its information and loses its pull; \code{"off"} turns every
#' prior off, which leaves a parameter only a prior identifies unidentified. A Dirichlet (on
#' movement and on the recruitment apportionment) or a mean and sd beta (on the stray rate and the
#' tag reporting rate) is moved the same way, its mode put on the value the replicate ran on and
#' its concentration kept, so it has no gradient there (\code{\link{dirichlet_mode_at}},
#' \code{\link{beta_mode_at}}); one with no interior mode is left as given, as is the symmetric
#' beta on the tag reporting rate, which sits at one half by construction.
#'
#' @param data Data list of the fit.
#' @param prior_means \code{"assessment"}, \code{"truth"} or \code{"off"}.
#' @param pars,rep The parameter list and report the replicate ran on.
#'
#' @return \code{data} with its priors set.
#'
#' @keywords internal
self_test_priors <- function(data, prior_means, pars, rep) {

  if(prior_means == "assessment") return(data)

  if(prior_means == "off") {
    prior_flags <- c("Use_fish_q_prior", "Use_srv_q_prior", "Use_M_prior", "Use_h_prior", "use_r0_prior",
                     "Use_fish_selex_prior", "Use_srv_selex_prior", "Use_ret_selex_prior", "Use_Movement_Prior",
                     "use_conv_tag_fishrep_prior", "use_rec_region_prop_prior", "use_rec_seas_prop_prior", "use_stray_rate_prior")
    for(flag in intersect(prior_flags, names(data))) data[[flag]][] <- 0
    return(data)
  }

  # catchability, normal on log q with its mean on the natural scale
  for(prefix in c("fish", "srv")) {
    if(!isTRUE(any(data[[paste0("Use_", prefix, "_q_prior")]] == 1))) next
    q_prior <- data[[paste0(prefix, "_q_prior")]]
    ln_q <- pars[[paste0("ln_", prefix, "_q")]]
    for(i in seq_len(nrow(q_prior))) q_prior$mu[i] <- exp(ln_q[q_prior$region[i], q_prior$block[i], q_prior$fleet[i]])
    data[[paste0(prefix, "_q_prior")]] <- q_prior
  } # end prefix loop

  # natural mortality, read through its blocks
  if(isTRUE(any(data$Use_M_prior == 1))) {
    for(i in seq_len(nrow(data$M_prior))) {
      seas <- if(is.null(data$M_prior$seasblk)) 1 else data$M_prior$seasblk[i]
      M_idx <- data$M_blocks[data$M_prior$popblk[i], data$M_prior$regionblk[i], data$M_prior$yearblk[i], seas, data$M_prior$ageblk[i], data$M_prior$sexblk[i]]
      data$M_prior$mu[i] <- exp(pars$ln_M[M_idx])
    } # end i loop
  }

  # steepness and R0, both with their means on the natural scale
  if(isTRUE(any(data$Use_h_prior == 1))) {
    for(i in seq_len(nrow(data$h_prior))) data$h_prior$mu[i] <- rep$h_trans[data$h_prior$pop[i], data$h_prior$region[i]]
  }
  if(isTRUE(any(data$use_r0_prior == 1))) {
    for(i in seq_len(nrow(data$r0_prior))) data$r0_prior$mu[i] <- rep$R0[data$r0_prior$pop[i]]
  }

  # selectivity, a parameter on the log scale or the selectivity in the first year of its block
  for(prefix in c("fish", "srv", "ret")) {
    if(!isTRUE(any(data[[paste0("Use_", prefix, "_selex_prior")]] == 1))) next
    selex_prior <- data[[paste0(prefix, "_selex_prior")]]
    row_type <- if(is.null(selex_prior$type)) rep("par", nrow(selex_prior)) else selex_prior$type
    for(i in seq_len(nrow(selex_prior))) {
      r <- selex_prior$region[i]
      p <- selex_prior$par[i]
      b <- selex_prior$block[i]
      s <- selex_prior$sex[i]
      f <- selex_prior$fleet[i]
      if(row_type[i] == "value") {
        y <- min(which(data[[paste0(prefix, "_sel_blocks")]][r,,f] == b))
        selex_prior$mu[i] <- if(data[[paste0(prefix, "_selex_type")]] == 0) rep[[paste0(prefix, "_sel")]][1,r,y,1,p,s,f] else rep[[paste0(prefix, "_sel_l")]][r,y,p,s,f]
      } else selex_prior$mu[i] <- exp(pars[[paste0(prefix, "_fixed_sel_pars")]][r,p,b,s,f])
    } # end i loop
    data[[paste0(prefix, "_selex_prior")]] <- selex_prior
  } # end prefix loop

  # movement, a Dirichlet on the annual fractions out of a region
  if(isTRUE(data$Use_Movement_Prior == 1) && !is.null(data$Movement_prior)) {
    move_prior <- data$Movement_prior
    for(i in seq_len(nrow(move_prior))) {
      p <- move_prior$pop[i]
      r_from <- move_prior$region_from[i]
      y <- move_prior$year[i]
      seas <- move_prior$seas[i]
      a <- move_prior$age[i]
      s <- move_prior$sex[i]
      frac <- if(is.null(rep$Mrate)) rep$Movement[p,r_from,,y,seas,a,s] else as.matrix(Matrix::expm(methods::as(rep$Mrate[p,,,y,seas,a,s], "sparseMatrix")))[r_from,]
      move_prior$alpha[[i]] <- dirichlet_mode_at(move_prior$alpha[[i]], frac)
    } # end i loop
    data$Movement_prior <- move_prior
  }

  # recruitment apportionment over regions, and over the seasons the penalty reads
  if(isTRUE(data$use_rec_region_prop_prior == 1) && !is.null(data$rec_region_prop_prior)) {
    for(i in seq_len(nrow(data$rec_region_prop_prior))) {
      p <- data$rec_region_prop_prior$pop[i]
      data$rec_region_prop_prior$alpha[[i]] <- dirichlet_mode_at(data$rec_region_prop_prior$alpha[[i]], rep$rec_region_prop[p,])
    } # end i loop
  }
  if(isTRUE(data$use_rec_seas_prop_prior == 1) && isTRUE(data$use_fixed_rec_seas_prop == 0) && !is.null(data$rec_seas_prop_prior)) {
    seas_read <- if(isTRUE(data$rec_lag == 0) && isTRUE(data$spawn_seas > 1)) data$spawn_seas:data$n_seas else seq_len(data$n_seas)
    for(i in seq_len(nrow(data$rec_seas_prop_prior))) {
      p <- data$rec_seas_prop_prior$pop[i]
      data$rec_seas_prop_prior$alpha[[i]] <- dirichlet_mode_at(data$rec_seas_prop_prior$alpha[[i]], rep$rec_seas_prop[p,seas_read])
    } # end i loop
  }

  # the stray rate and the tag reporting rate, each a mean and sd beta on the rate the penalty reads
  if(isTRUE(data$use_stray_rate_prior == 1) && !is.null(data$stray_rate_prior)) {
    for(i in seq_len(nrow(data$stray_rate_prior))) {
      stray <- 1e-4 + (1 - 2 * 1e-4) * stats::plogis(pars$stray_rate_pars[data$stray_rate_prior$pop[i], data$stray_rate_prior$block[i]])
      moved <- beta_mode_at(data$stray_rate_prior$mu[i], data$stray_rate_prior$sd[i], stray)
      data$stray_rate_prior$mu[i] <- moved[1]
      data$stray_rate_prior$sd[i] <- moved[2]
    } # end i loop
  }
  if(isTRUE(data$use_conv_tag_fishrep_prior == 1) && !is.null(data$conv_tag_fishrep_prior)) {
    tag_prior <- data$conv_tag_fishrep_prior
    for(i in which(tag_prior$type == 1)) {
      reporting <- stats::plogis(pars$conv_tag_fish_reporting_pars[tag_prior$region[i], tag_prior$block[i], tag_prior$fleet[i]])
      moved <- beta_mode_at(tag_prior$mu[i], tag_prior$sd[i], reporting)
      tag_prior$mu[i] <- moved[1]
      tag_prior$sd[i] <- moved[2]
    } # end i loop
    data$conv_tag_fishrep_prior <- tag_prior
  }

  data

} # end self_test_priors

#' A Dirichlet moved so its mode sits at given fractions
#'
#' The concentration (the sum of \code{alpha}) is kept, so the prior keeps its
#' strength and loses its pull at \code{frac}. A Dirichlet whose concentration is
#' at or below its number of cells has no interior mode and is returned as given.
#'
#' @param alpha Concentration vector.
#' @param frac Fractions summing to one, the same length.
#'
#' @return The moved \code{alpha}.
#'
#' @keywords internal
dirichlet_mode_at <- function(alpha, frac) {
  if(sum(alpha) <= length(alpha)) return(alpha)
  1 + (sum(alpha) - length(alpha)) * frac
}

#' A mean and sd beta moved so its mode sits at a given value
#'
#' The concentration \code{mu * (1 - mu) / sd^2 - 1} is kept and the mean and sd
#' are read back from the moved shape parameters, as
#' \code{\link{get_tagrep_prior}} forms them. A beta whose concentration is at
#' or below two has no interior mode and is returned as given.
#'
#' @param mu,sd The prior's mean and sd on the natural scale.
#' @param x Value in (0, 1) the mode is put at.
#'
#' @return The moved \code{c(mu, sd)}.
#'
#' @keywords internal
beta_mode_at <- function(mu, sd, x) {
  conc <- mu * (1 - mu) / sd^2 - 1
  if(conc <= 2) return(c(mu, sd))
  mu_new <- (1 + (conc - 2) * x) / conc
  c(mu_new, sqrt(mu_new * (1 - mu_new) / (conc + 1)))
}

#' One replicate's true value of a reported quantity or parameter
#'
#' The operating model's own value where it holds the quantity under the same
#' name with a replicate dim last, so each replicate is compared with what it
#' ran on rather than with the fit. Where the fit's
#' array runs longer along one dim, a projected year say, the operating model's
#' cells fill the leading part and the fit's the rest; where the operating
#' model's runs longer, as recruitment deviations that stop before the last year
#' do, the fit's cells take its leading part. Anything else is the fit's.
#'
#' @param om_arr The operating model's array, or \code{NULL}.
#' @param fit_arr The fit's value, or under a joint self test the replicate's draw.
#' @param sim Replicate.
#' @param n_sims Number of replicates, which the operating model's last dim must be.
#'
#' @return \code{fit_arr} with the operating model's values in it.
#'
#' @keywords internal
om_truth <- function(om_arr, fit_arr, sim, n_sims) {

  om_dim <- dim(om_arr)
  n_om <- length(om_dim)
  fit_dim <- if(is.null(dim(fit_arr))) length(fit_arr) else dim(fit_arr)
  if(is.null(om_dim) || is.null(fit_arr) || om_dim[n_om] != n_sims || length(fit_dim) != n_om - 1) return(fit_arr)

  rep_dim <- om_dim[-n_om] # one replicate's shape
  replicate_cells <- array(om_arr[slice.index(om_arr, n_om) == sim], dim = rep_dim)
  if(all(fit_dim == rep_dim)) return(array(replicate_cells, dim = fit_dim))

  # one array runs longer along one dim: the operating model's leading part when it is the longer one
  longer <- which(fit_dim != rep_dim)
  if(length(longer) != 1) return(fit_arr)
  if(fit_dim[longer] < rep_dim[longer]) {
    replicate_arr <- array(replicate_cells, dim = rep_dim)
    return(array(replicate_arr[slice.index(replicate_arr, longer) <= fit_dim[longer]], dim = fit_dim))
  }

  # otherwise the fit's array runs longer, and the operating model fills its leading part
  truth <- array(fit_arr, dim = fit_dim)
  truth[slice.index(truth, longer) <= rep_dim[longer]] <- replicate_cells
  truth

} # end om_truth

#' Seasons a fit holds its year totals in
#'
#' For each data source a fleet reports once a year, the season its \code{Use}
#' array turns on in each year, which is where the operating model draws that
#' year's total (\code{seas_agg_slot} in \code{\link{Setup_Sim_Fishing}} and
#' \code{\link{Setup_Sim_Survey}}). A year with no observation takes season one,
#' and years past the fit repeat the last fitted year.
#'
#' @param data Data list of the fitted model.
#' @param n_yrs Number of years the operating model runs.
#' @param platform \code{"fish"} or \code{"srv"}.
#'
#' @return Named list of integer matrices \code{[n_yrs, n_fleets]}, one for each
#'   data source with a fleet reporting once a year.
#'
#' @keywords internal
seas_agg_slot_list <- function(data, n_yrs, platform) {

  data_names <- if(platform == "fish") c("Catch", "Catch_pop", "Discard", "Discard_pop", "FishIdx", "FishIdx_pop",
                                         "FishAgeComps", "FishLenComps", "FishAgeComps_pop", "FishLenComps_pop",
                                         "FishAgeComps_discard", "FishLenComps_discard",
                                         "FishAgeComps_discard_pop", "FishLenComps_discard_pop")
                else c("SrvIdx", "SrvIdx_pop", "SrvAgeComps", "SrvLenComps", "SrvAgeComps_pop", "SrvLenComps_pop")
  n_fleets <- if(platform == "fish") data$n_fish_fleets else data$n_srv_fleets
  slots <- list()

  for(data_name in data_names) {

    seas_agg <- data[[paste0(data_name, "_seas_Type")]]
    use <- data[[paste0("Use", data_name)]]
    if(is.null(seas_agg) || is.null(use) || !any(seas_agg == 1)) next
    pop_source <- grepl("_pop$", data_name) # population leads the use array
    n_fit_yrs <- dim(use)[if(pop_source) 3 else 2]
    slot <- matrix(1, nrow = n_yrs, ncol = n_fleets) # season one unless the fit says otherwise

    for(f in which(seas_agg == 1)) {
      for(y in seq_len(n_yrs)) {

        y_fit <- min(y, n_fit_yrs) # past the fit, the last fitted year
        season_on <- if(pop_source) apply(use[,,y_fit,,f,drop = FALSE] == 1, 4, any)
                     else apply(use[,y_fit,,f,drop = FALSE] == 1, 3, any)
        seasons <- which(season_on)
        if(length(seasons) > 1)
          stop("Use", data_name, " puts fleet ", f, "'s year total for year ", y_fit, " in seasons ",
               paste(seasons, collapse = " and "), " across its regions or populations. The operating model ",
               "draws a fleet's year total into one season, so put every region's and population's in the same one.")
        if(length(seasons) == 1) slot[y,f] <- seasons

      } # end y loop
    } # end f loop

    slots[[data_name]] <- slot

  } # end data_name loop

  slots

} # end seas_agg_slot_list

#' Parameters each replicate's operating model runs on
#'
#' Conditional gives every replicate the fitted parameters, so one truth covers
#' them all and only the data change. Joint gives each replicate its own draw, so
#' each has its own truth.
#'
#' The draw is taken at the fitted parameter vector, the order the joint precision
#' is in, and \code{parList} puts the map back. A fit with no random effects has no
#' joint precision, so the fixed effect covariance is inverted in its place.
#'
#' @param sim_type Either \code{"conditional"} or \code{"joint"}.
#' @param n_sims Number of replicates.
#' @param fit_rep Report list from the fitted model.
#' @param parameters Parameter list the model was built with.
#' @param mapping Factor maps the model was built with.
#' @param sd_rep \code{sdreport} object from the fitted model.
#' @param random Character vector of random effect names.
#' @param obj Fitted object, needed only under \code{"joint"}.
#'
#' @return \code{reps} and \code{pars}, one per replicate, plus \code{fit_pars}
#'   at the fit.
#'
#' @keywords internal
sim_draw_views <- function(sim_type, n_sims, fit_rep, parameters, mapping, sd_rep, random, obj) {

  # get fitted pars
  fit_pars <- get_optim_param_list(parameters, mapping, sd_rep, random)

  # conditional sims
  if(sim_type == "conditional")  return(list(reps = rep(list(fit_rep), n_sims), pars = rep(list(fit_pars), n_sims), fit_pars = fit_pars))
  if(is.null(obj)) stop("simulation_self_test: sim_type = 'joint' draws at the fitted parameter vector, so pass obj = the fitted object.")

  # joint sims
  prec <- sd_rep$jointPrecision
  if(is.null(prec)) {
    # sdreport only assembles a joint precision when the fit has random effects. without them the
    # joint distribution is the fixed effects on their own, whose precision is the inverse covariance
    if(length(obj$env$random) > 0) stop("simulation_self_test: sim_type = 'joint' needs sd_rep$jointPrecision, so call RTMB::sdreport(obj, getJointPrecision = TRUE).")
    prec <- methods::as(Matrix::forceSymmetric(solve(sd_rep$cov.fixed)), "CsparseMatrix")
  }

  # make draws
  draws <- rmvnorm_prec(obj$env$last.par.best, prec, n_sims = n_sims)

  # output
  list(reps = lapply(seq_len(n_sims), function(i) obj$report(draws[, i])),
       pars = lapply(seq_len(n_sims), function(i) obj$env$parList(par = draws[, i])),
       fit_pars = fit_pars)
}

#' Run a simulation self-test of a fitted RTMB estimation model
#'
#' Validates model performance by: (1) generating \code{n_sims} new datasets
#' from the fitted model parameters using \code{\link{Simulate_Pop_Static}},
#' (2) re-fitting the estimation model to each simulated dataset, and (3)
#' storing user-specified report quantities for comparison against true values.
#' Supports sequential or parallel execution via
#' \code{future}/\code{future.apply}. Likelihood weights from the original fit
#' are propagated into the simulation (e.g., ISS scaled by
#' \code{Wt_FishAgeComps}; \code{ObsSrvIdx_SE} divided by
#' \code{sqrt(Wt_SrvIdx)}), and the data weights are reset to 1 when refitting.
#' A penalty weight (\code{Wt_Rec}, \code{Wt_Init_Rec}, \code{Wt_F},
#' \code{Wt_D}, \code{*_pe_wt}) is kept, the operating model drawing that
#' process at its sd over the root of the weight, so each refit is the fit's
#' own estimator. A weighted penalty whose sd is estimated comes back at the
#' sd over the root of the weight, a weighted penalty being a density only up
#' to its normalizing constant.
#' Failed replicates are silently stored as \code{NA}. The operating model takes
#' the fit's \code{bias_correct_pe}, \code{bias_correct_oe} and
#' \code{sigmaR_switch}, so both sides center process and observation error alike.
#'
#' @param data Named list of model data from a fitted RTMB object
#'   (\code{$data}).
#' @param parameters Named list of fitted parameter values (\code{$par} or
#'   equivalent).
#' @param mapping Named list of parameter factor maps (\code{$map}).
#' @param random Character vector of random effect names passed to RTMB.
#' @param rep Named list of report values from the fitted model
#'   (\code{obj$rep}).
#' @param sd_rep \code{sdreport} object from the fitted model, used to
#'   extract optimized parameter values in list format via
#'   \code{get_optim_param_list}.
#' @param n_sims Integer. Number of simulation replicates.
#' @param newton_loops Integer. Number of Newton refinement steps applied
#'   during re-fitting. Default \code{3}.
#' @param do_sdrep Logical. Whether to compute \code{sdreport} for each
#'   fitted replicate. Results stored as \code{$sd_rep} in the output list;
#'   failed \code{sdreport} calls stored as \code{NA}. Default \code{FALSE}.
#' @param do_par Logical. Whether to run replicates in parallel via
#'   \code{future::multisession}. Default \code{FALSE}.
#' @param n_cores Integer. Number of parallel workers. If \code{NULL}
#'   (default), \code{parallel::detectCores() - 1} is used.
#' @param output_path Character string. Path to save the simulated dataset
#'   RDS file. Passed to \code{\link{Simulate_Pop_Static}}. Default
#'   \code{NULL}.
#' @param obj Fitted object the self test is run from. Needed under
#'   \code{sim_type = "joint"}, which draws at its parameter vector and reports
#'   through it. Default \code{NULL}.
#' @param what Character vector. Names of report elements (keys of
#'   \code{rep}) to extract and store from each replicate. An error is raised
#'   if any name is not found in \code{rep}. Default \code{c("SSB", "Rec")}.
#' @param what_par Character vector. Names of parameters (keys of
#'   \code{parameters}) to extract and store from each replicate, read off the
#'   refit's own parameter list so that mapped elements come back at the values
#'   the map gave them. An error is raised if any name is not found in
#'   \code{parameters}. Default \code{NULL}, which stores none.
#' @param perfect_data Logical. Whether to shrink the observation error before
#'   simulating, sds to 0.001 and sample sizes to 1e6, leaving process error alone.
#'   A correct model whose data identify the population then returns the
#'   operating model to several decimals. A state-space model need not: exact
#'   catch and survey indices at age still leave each age's scale to cohort
#'   continuity under process error, so recovery there is to a few percent. A
#'   process error sd conditioned on the fit comes back low, since the fitted
#'   value includes the posterior variance of its own deviations and data this
#'   clean remove it; on NEA cod the numbers at age sd returns 0.12 against 0.19.
#'   The truth for an estimated observation sd is the 0.001 the data were drawn at.
#'   Default \code{FALSE}.
#' @param prior_means Where each refit's priors are centered (see
#'   \code{\link{self_test_priors}}). \code{"assessment"} (default) keeps them as
#'   the assessment gives them, so every replicate is pulled toward a prior the fit
#'   sits away from; \code{"truth"} centers each catchability, natural mortality,
#'   steepness, R0 and selectivity prior on the value the replicate ran on,
#'   keeping its sd, and puts the mode of each Dirichlet (movement, recruitment
#'   apportionment) and mean and sd beta (stray rate, tag reporting) there,
#'   keeping its concentration; \code{"off"} turns every prior off.
#' @param sim_type Character. Which self test to run. \code{"conditional"}
#'   (default) runs every replicate at the fitted parameters and every process
#'   deviation at the fit's estimate, so one population covers every replicate
#'   and only the observations are new. It measures how well the data determine
#'   this history; a process sd comes back low, since the fit's deviations are
#'   shrunk estimates that vary less than the sd describes. \code{"joint"} draws
#'   each replicate's parameters from the fit's uncertainty
#'   (\code{sd_rep$jointPrecision}) and then every process fresh at them:
#'   recruitment and the initial ages from their penalty
#'   (\code{\link{rec_devs_past_fit}}, \code{\link{init_devs_past_fit}}), the
#'   numbers at age, selectivity, catchability, movement, growth and a linked
#'   dsem from their processes, and F and discard mortality deviations when the
#'   fit integrates them out (fixed effect F deviations come from the parameter
#'   draw). Each replicate is a new population with its own truth, so the self
#'   test measures bias in every estimate, process sds included, and with
#'   \code{do_sdrep = TRUE} whether the standard errors cover the truth. A fit
#'   with no random effects has no joint precision, and the fixed effect
#'   covariance is inverted in its place.
#'
#'   Joint moves F, both selectivities, catchability, natural mortality, weight and
#'   size at age, movement, steepness, sex ratio, recruitment, the initial
#'   deviations, the numbers at age and a linked dsem. Observation error and the
#'   composition parameters stay at the fit, having no replicate dim.
#'
#' @return Named list with one element per entry in \code{what} and then one per
#'   entry in \code{what_par}, each an array with the last dimension indexing
#'   simulation replicates (via \code{simplify2array}). If \code{do_sdrep = TRUE},
#'   an additional element \code{"sd_rep"} contains a list of \code{sdreport}
#'   objects (or \code{NA} for failed replicates). A final element \code{"truth"}
#'   holds each replicate's true values for the same names, the operating model's
#'   own wherever it holds the quantity (\code{\link{om_truth}}), so replicates
#'   drawn under \code{sim_type = "joint"} are each compared with what they ran on.
#'   Observation error, composition, at-age correlation and tag loss parameters
#'   take the fit's values, which every replicate runs at.
#'
#'
#' @export
#' @family Simulation Setup
#'
#' @examples
#' \dontrun{
#' res <- simulation_self_test(
#'   data = obj$data, parameters = obj$par, mapping = obj$map,
#'   random = obj$random, rep = obj$rep, sd_rep = obj$sd_rep,
#'   n_sims = 100, what = c("SSB", "Rec", "Fmort")
#' )
#' str(res$SSB)
#'
#' # parameters drawn from the fit's uncertainty, each replicate compared with its own truth
#' sd_rep <- RTMB::sdreport(fit, getJointPrecision = TRUE)
#' res <- simulation_self_test(
#'   data = fit$data, parameters = par, mapping = map, random = NULL,
#'   rep = fit$rep, sd_rep = sd_rep, obj = fit, n_sims = 100,
#'   what = "SSB", what_par = "ln_global_R0", sim_type = "joint"
#' )
#' rel_err <- (res$SSB - res$truth$SSB) / res$truth$SSB
#' }
simulation_self_test <- function(
  data,
  parameters,
  mapping,
  random,
  rep,
  sd_rep,
  n_sims,
  newton_loops = 3,
  do_sdrep = FALSE,
  do_par = FALSE,
  obj = NULL,
  n_cores = NULL,
  output_path = NULL,
  what = c('SSB', 'Rec'),
  what_par = NULL,
  perfect_data = FALSE,
  sim_type = c("conditional", "joint"),
  prior_means = c("assessment", "truth", "off")
) {

  sim_type <- match.arg(sim_type)
  prior_means <- match.arg(prior_means)
  n_cond_yrs <- if(sim_type == "joint") 0 else length(data$years) # conditional keeps every process deviation at the fit; joint draws every one fresh
  warn_R0_ref_block_om(data, "simulation_self_test")

  missing_names <- setdiff(what, names(rep))
  if(length(missing_names) > 0)  stop(paste("The following elements in 'what' are not found in rep:",  paste(missing_names, collapse = ", ")))
  missing_pars <- setdiff(what_par, names(parameters))
  if(length(missing_pars) > 0) stop(paste("The following elements in 'what_par' are not found in parameters:", paste(missing_pars, collapse = ", ")))

  # make draws if joint sim, if conditional, just return the report values
  views <- sim_draw_views(sim_type, n_sims, rep, parameters, mapping, sd_rep, random, obj)
  optim_parameters_list <- views$fit_pars

  # weights become simulation sds as sd / sqrt(wt). a weight of zero means excluded, not infinite
  # error, so excluded cells keep their nominal sd rather than giving Inf and then NaN
  deweight <- function(sd, wt) {
    w <- wt
    w[!is.finite(w) | w <= 0] <- 1
    sd / sqrt(w)
  }

  # an index sd with an estimated part (sigma*Idx_form) is drawn unweighted: a weight scales the likelihood, so it
  # cancels from that sd's estimate, which already matches the fit's residual spread
  idx_draw_se <- function(se, wt, form) if(is.null(form) || form == 0) deweight(se, wt) else se

  # Modify any data weights that are NA to 0
  if(any(is.na(data$Wt_Catch))) data$Wt_Catch[is.na(data$Wt_Catch)] <- 0
  if(any(is.na(data$Wt_Catch_pop))) data$Wt_Catch_pop[is.na(data$Wt_Catch_pop)] <- 0
  if(any(is.na(data$Wt_FishAgeComps))) data$Wt_FishAgeComps[is.na(data$Wt_FishAgeComps)] <- 0
  if(any(is.na(data$Wt_FishLenComps))) data$Wt_FishLenComps[is.na(data$Wt_FishLenComps)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_discard))) data$Wt_FishAgeComps_discard[is.na(data$Wt_FishAgeComps_discard)] <- 0
  if(any(is.na(data$Wt_FishLenComps_discard))) data$Wt_FishLenComps_discard[is.na(data$Wt_FishLenComps_discard)] <- 0
  if(any(is.na(data$Wt_FishIdx))) data$Wt_FishIdx[is.na(data$Wt_FishIdx)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_pop))) data$Wt_FishAgeComps_pop[is.na(data$Wt_FishAgeComps_pop)] <- 0
  if(any(is.na(data$Wt_FishLenComps_pop))) data$Wt_FishLenComps_pop[is.na(data$Wt_FishLenComps_pop)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_discard_pop))) data$Wt_FishAgeComps_discard_pop[is.na(data$Wt_FishAgeComps_discard_pop)] <- 0
  if(any(is.na(data$Wt_FishLenComps_discard_pop))) data$Wt_FishLenComps_discard_pop[is.na(data$Wt_FishLenComps_discard_pop)] <- 0
  if(any(is.na(data$Wt_FishIdx_pop))) data$Wt_FishIdx_pop[is.na(data$Wt_FishIdx_pop)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps))) data$Wt_SrvAgeComps[is.na(data$Wt_SrvAgeComps)] <- 0
  if(any(is.na(data$Wt_SrvLenComps))) data$Wt_SrvLenComps[is.na(data$Wt_SrvLenComps)] <- 0
  if(any(is.na(data$Wt_SrvIdx))) data$Wt_SrvIdx[is.na(data$Wt_SrvIdx)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps_pop))) data$Wt_SrvAgeComps_pop[is.na(data$Wt_SrvAgeComps_pop)] <- 0
  if(any(is.na(data$Wt_SrvLenComps_pop))) data$Wt_SrvLenComps_pop[is.na(data$Wt_SrvLenComps_pop)] <- 0
  if(any(is.na(data$Wt_SrvIdx_pop))) data$Wt_SrvIdx_pop[is.na(data$Wt_SrvIdx_pop)] <- 0
  if(any(is.na(data$Wt_Tagging))) data$Wt_Tagging[is.na(data$Wt_Tagging)] <- 0

  # Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims, # number of simulations
                            n_yrs = length(data$years), # number of years
                            n_regions = data$n_regions,  # number of regions
                            n_ages = length(data$ages), # number of ages
                            # Use fishery or survey observed ages depending on what is available
                            n_obs_ages = if(any(data$UseFishAgeComps == 1)) {
                              dim(data$ObsFishAgeComps)[4]
                            } else if(any(data$UseFishAgeComps_pop == 1)) {
                              dim(data$ObsFishAgeComps_pop)[5]
                            } else if(any(data$UseSrvAgeComps == 1)) {
                              dim(data$ObsSrvAgeComps)[4]
                            } else if(any(data$UseSrvAgeComps_pop == 1)) {
                              dim(data$ObsSrvAgeComps_pop)[5]
                            } else if(!is.null(data$UseFish_caal) && any(data$UseFish_caal == 1)) {
                              dim(data$ObsFish_caal)[5] # conditional age-at-length holds the observed ages
                            } else if(!is.null(data$UseSrv_caal) && any(data$UseSrv_caal == 1)) {
                              dim(data$ObsSrv_caal)[5]
                            } else {
                              dim(data$AgeingError)[3] # otherwise the ageing error's observed ages, which the at-age data sit on
                            },
                            n_lens = length(data$lens), # number of lengths
                            n_obs_lens = if(is.null(data$LenBinMap)) length(data$lens) else ncol(data$LenBinMap), # length bins the comps are recorded on
                            n_caal_lens = length(caal_row_lens(data)), # age-at-length rows, each covering the bins CAAL_LenBinMap gives it
                            n_sexes = data$n_sexes, # number of sexes
                            n_fish_fleets = data$n_fish_fleets, # number of fishery fleets
                            n_srv_fleets = data$n_srv_fleets, # number of survey fleets
                            # Seasonal stuff
                            n_seas = data$n_seas,
                            seasdur = data$seasdur,
                            # Population stuff
                            n_pop = data$n_pop,
                            natal_region = data$natal_region,
                            # the fit's centering of process and observation error, so both sides draw what the other evaluates
                            bias_correct_pe = if(is.null(data$bias_correct_pe)) "rec" else data$bias_correct_pe,
                            bias_correct_oe = if(is.null(data$bias_correct_oe)) 0 else data$bias_correct_oe
  )

  # Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list)

  # Catchability: the reported value is the mean times the fit's deviation, and the OM wants those split
  fit_yrs <- seq_along(data$years)
  fish_q_fit <- lapply(seq_len(n_sims), function(i)
    split_reported_q(views$reps[[i]]$fish_q[,fit_yrs,,drop = FALSE], views$pars[[i]]$ln_fish_q,
                     data$fish_q_blocks[,fit_yrs,,drop = FALSE], data$fish_q_type))
  srv_q_fit <- lapply(seq_len(n_sims), function(i)
    split_reported_q(views$reps[[i]]$srv_q[,fit_yrs,,drop = FALSE], views$pars[[i]]$ln_srv_q,
                     data$srv_q_blocks[,fit_yrs,,drop = FALSE], data$srv_q_type))

  # Setup Fishing Processes -------------------------------------------------

  # Region-specific sigmaC
  ln_sigmaC <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch)) ln_sigmaC[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC[r,,,f]), data$Wt_Catch[r,,,f]))
    else ln_sigmaC[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC[r,,,f]), data$Wt_Catch))
  }

  # Population-specific sigmaC
  ln_sigmaC_pop <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch_pop)) ln_sigmaC_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]), data$Wt_Catch_pop[p,r,,,f]))
    else ln_sigmaC_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]), data$Wt_Catch_pop))
  }

  # Region-specific sigmaD
  ln_sigmaD <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaD container
  # Loop through to populate ln_sigmaD with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Discard)) ln_sigmaD[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD[r,,,f]), data$Wt_Discard[r,,,f]))
    else ln_sigmaD[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD[r,,,f]), data$Wt_Discard))
  }

  # Population-specific sigmaD
  ln_sigmaD_pop <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaD container
  # Loop through to populate ln_sigmaD with associated weights
  for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Discard_pop)) ln_sigmaD_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]), data$Wt_Discard_pop[p,r,,,f]))
    else ln_sigmaD_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]), data$Wt_Discard_pop))
  }

  n_obs_om <- sim_list$n_obs_ages # observed ages the operating model draws on
  catch_aa_used <- any(data$UseCatchAA == 1) # whether the fit observes catch at age
  discard_aa_used <- any(data$UseDiscardAA == 1) # discards at age
  srv_idx_aa_used <- any(data$UseSrvIdxAA == 1) # survey index at age
  catch_aa_pop_used <- any(data$UseCatchAA_pop == 1) # population-specific catch at age
  discard_aa_pop_used <- any(data$UseDiscardAA_pop == 1) # population-specific discards at age
  srv_idx_aa_pop_used <- any(data$UseSrvIdxAA_pop == 1) # population-specific survey index at age

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(sim_list = sim_list,
                                ln_sigmaC = ln_sigmaC,
                                ln_sigmaC_pop = ln_sigmaC_pop,
                                ln_sigmaD = ln_sigmaD,
                                ln_sigmaD_pop = ln_sigmaD_pop,
                                ln_sigmaCAA = unused_at_age_on_obs_ages(optim_parameters_list$ln_sigmaCAA, catch_aa_used, 1, n_obs_om, log(0.5)),
                                ln_sigmaDAA = unused_at_age_on_obs_ages(optim_parameters_list$ln_sigmaDAA, discard_aa_used, 1, n_obs_om, log(0.5)),
                                UseCatchAA = unused_at_age_on_obs_ages(data$UseCatchAA, catch_aa_used, 4, n_obs_om),
                                UseDiscardAA = unused_at_age_on_obs_ages(data$UseDiscardAA, discard_aa_used, 4, n_obs_om),
                                use_catch_aa = data$use_catch_aa,
                                use_discard_aa = data$use_discard_aa,
                                # the fit's correlation of each fleet's at-age residuals, under the estimation model's names
                                AgeObsCorr_catch = data$AgeObsCorr_catch,
                                AgeObsCorr_discard = data$AgeObsCorr_discard,
                                trans_rho_catch = if(catch_aa_used) optim_parameters_list$trans_rho_catch,
                                trans_rho_catch_year = if(catch_aa_used) optim_parameters_list$trans_rho_catch_year,
                                trans_rho_catch_us = if(catch_aa_used) optim_parameters_list$trans_rho_catch_us,
                                trans_rho_discard = if(discard_aa_used) optim_parameters_list$trans_rho_discard,
                                trans_rho_discard_year = if(discard_aa_used) optim_parameters_list$trans_rho_discard_year,
                                trans_rho_discard_us = if(discard_aa_used) optim_parameters_list$trans_rho_discard_us,
                                # the population-specific at-age data sources and the year totals, as the fit reads them
                                CatchAA_seas_Type = data$CatchAA_seas_Type,
                                DiscardAA_seas_Type = data$DiscardAA_seas_Type,
                                UseCatchAA_pop = if(catch_aa_pop_used) data$UseCatchAA_pop,
                                UseDiscardAA_pop = if(discard_aa_pop_used) data$UseDiscardAA_pop,
                                ln_sigmaCAA_pop = if(catch_aa_pop_used) optim_parameters_list$ln_sigmaCAA_pop,
                                ln_sigmaDAA_pop = if(discard_aa_pop_used) optim_parameters_list$ln_sigmaDAA_pop,
                                ObsCatchAA_pop_SE = if(catch_aa_pop_used) data$ObsCatchAA_pop_SE,
                                ObsDiscardAA_pop_SE = if(discard_aa_pop_used) data$ObsDiscardAA_pop_SE,
                                CatchAA_pop_Type = data$CatchAA_pop_Type,
                                DiscardAA_pop_Type = data$DiscardAA_pop_Type,
                                CatchAA_pop_LikeType = data$CatchAA_pop_LikeType,
                                DiscardAA_pop_LikeType = data$DiscardAA_pop_LikeType,
                                CatchAA_pop_sigma_form = data$CatchAA_pop_sigma_form,
                                DiscardAA_pop_sigma_form = data$DiscardAA_pop_sigma_form,
                                CatchAA_pop_seas_Type = data$CatchAA_pop_seas_Type,
                                DiscardAA_pop_seas_Type = data$DiscardAA_pop_seas_Type,
                                AgeObsCorr_catch_pop = data$AgeObsCorr_catch_pop,
                                AgeObsCorr_discard_pop = data$AgeObsCorr_discard_pop,
                                trans_rho_catch_pop = if(catch_aa_pop_used) optim_parameters_list$trans_rho_catch_pop,
                                trans_rho_catch_pop_year = if(catch_aa_pop_used) optim_parameters_list$trans_rho_catch_pop_year,
                                trans_rho_catch_pop_us = if(catch_aa_pop_used) optim_parameters_list$trans_rho_catch_pop_us,
                                trans_rho_discard_pop = if(discard_aa_pop_used) optim_parameters_list$trans_rho_discard_pop,
                                trans_rho_discard_pop_year = if(discard_aa_pop_used) optim_parameters_list$trans_rho_discard_pop_year,
                                trans_rho_discard_pop_us = if(discard_aa_pop_used) optim_parameters_list$trans_rho_discard_pop_us,
                                ObsCatchAA_SE = unused_at_age_on_obs_ages(data$ObsCatchAA_SE, catch_aa_used, 4, n_obs_om),
                                ObsDiscardAA_SE = unused_at_age_on_obs_ages(data$ObsDiscardAA_SE, discard_aa_used, 4, n_obs_om),
                                CatchAA_Type = data$CatchAA_Type,
                                DiscardAA_Type = data$DiscardAA_Type,
                                CatchAA_LikeType = data$CatchAA_LikeType,
                                DiscardAA_LikeType = data$DiscardAA_LikeType,
                                CatchAA_sigma_form = data$CatchAA_sigma_form,
                                DiscardAA_sigma_form = data$DiscardAA_sigma_form,
                                catch_units = data$catch_units,
                                discard_units = data$discard_units,
                                # data sources the fitted model reports once a year stay annual
                                Catch_seas_Type = data$Catch_seas_Type,
                                Catch_pop_seas_Type = data$Catch_pop_seas_Type,
                                FishIdx_seas_Type = data$FishIdx_seas_Type,
                                FishIdx_pop_seas_Type = data$FishIdx_pop_seas_Type,
                                FishAgeComps_seas_Type = data$FishAgeComps_seas_Type,
                                Discard_seas_Type = data$Discard_seas_Type,
                                Discard_pop_seas_Type = data$Discard_pop_seas_Type,
                                FishLenComps_seas_Type = data$FishLenComps_seas_Type,
                                FishAgeComps_pop_seas_Type = data$FishAgeComps_pop_seas_Type,
                                FishLenComps_pop_seas_Type = data$FishLenComps_pop_seas_Type,
                                FishAgeComps_discard_seas_Type = data$FishAgeComps_discard_seas_Type,
                                FishLenComps_discard_seas_Type = data$FishLenComps_discard_seas_Type,
                                FishAgeComps_discard_pop_seas_Type = data$FishAgeComps_discard_pop_seas_Type,
                                FishLenComps_discard_pop_seas_Type = data$FishLenComps_discard_pop_seas_Type,
                                seas_agg_slot = seas_agg_slot_list(data, length(data$years), "fish"), # the season each year total sits in
                                Fmort_input = bind_sims(lapply(views$reps, function(rp) rp$Fmort[,seq_along(data$years),,,drop = FALSE])),
                                dmr_input = bind_sims(lapply(views$reps, function(rp) rp$dmr[,seq_along(data$years),,,drop = FALSE])),
                                fish_sel_input = bind_sims(lapply(views$reps, function(rp) rp$fish_sel[,,seq_along(data$years),,,,,drop = FALSE])),
                                ret_sel_input = bind_sims(lapply(views$reps, function(rp) rp$ret_sel[,,seq_along(data$years),,,,,drop = FALSE])),
                                # length comps selected at length read the fit's selectivity at length
                                FishLenComps_sel = if(is.null(data$fish_len_comp_sel)) rep("age", data$n_fish_fleets) else ifelse(data$fish_len_comp_sel == 1, "length", "age"),
                                fish_sel_l_input = if(any(data$fish_len_comp_sel == 1)) bind_sims(lapply(views$reps, function(rp) rp$fish_sel_l[,seq_along(data$years),,,,drop = FALSE])) else NULL,
                                ret_sel_l_input = if(any(data$fish_len_comp_sel == 1) && isTRUE(data$ret_selex_type == 1)) bind_sims(lapply(views$reps, function(rp) rp$ret_sel_l[,seq_along(data$years),,,,drop = FALSE])) else NULL,
                                fish_q_input = bind_sims(lapply(fish_q_fit, function(q) q$q_mean)),
                                # the reported errors, with any estimated part of the sd drawn on top as the fit forms it
                                ObsFishIdx_SE = idx_draw_se(data$ObsFishIdx_SE, data$Wt_FishIdx, data$sigmaFishIdx_form),
                                sigmaFishIdx_form = if(is.null(data$sigmaFishIdx_form)) 0 else data$sigmaFishIdx_form,
                                ln_sigmaFishIdx = optim_parameters_list$ln_sigmaFishIdx,
                                sigmaFishIdx_pop_form = if(is.null(data$sigmaFishIdx_pop_form)) 0 else data$sigmaFishIdx_pop_form,
                                ln_sigmaFishIdx_pop = optim_parameters_list$ln_sigmaFishIdx_pop,
                                ObsFishIdx_pop_SE = if(any(data$UseFishIdx_pop == 1)) {
                                  idx_draw_se(data$ObsFishIdx_pop_SE, data$Wt_FishIdx_pop, data$sigmaFishIdx_pop_form)
                                } else {
                                  array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
                                },
                                fish_idx_type = data$fish_idx_type,
                                fish_idx_ages = data$fish_idx_ages, # ages in each index total
                                t_fish = if(is.null(data$t_fish)) array(0, dim = c(data$n_regions, data$n_seas, data$n_fish_fleets)) else data$t_fish, # fishery index timing
                                FishIdx_LikeType = if(is.null(data$FishIdx_LikeType)) rep(0, data$n_fish_fleets) else data$FishIdx_LikeType,
                                FishIdx_Cov = data$FishIdx_Cov,
                                UseFishIdx = data$UseFishIdx,
                                init_F_val = rep$init_F,

                                # fishery age composition specifications
                                comp_fishage_like = data$FishAgeComps_LikeType,
                                FishAgeComps_Type = data$FishAgeComps_Type,
                                ISS_FishAgeComps = replicate(sim_list$n_sims, data$ISS_FishAgeComps[,,,,,drop = FALSE] * data$Wt_FishAgeComps),
                                ln_FishAge_theta = optim_parameters_list$ln_FishAge_theta[,,,drop = FALSE],
                                ln_FishAge_theta_agg = optim_parameters_list$ln_FishAge_theta_agg,
                                FishAge_corr_pars_agg = optim_parameters_list$FishAge_corr_pars_agg,
                                FishAge_corr_pars = optim_parameters_list$FishAge_corr_pars[,,,,drop = FALSE],

                                # fishery length composition specifications
                                comp_fishlen_like = data$FishLenComps_LikeType,
                                FishLenComps_Type = data$FishLenComps_Type,
                                ISS_FishLenComps = replicate(sim_list$n_sims, data$ISS_FishLenComps[,,,,,drop = FALSE] * data$Wt_FishLenComps),
                                ln_FishLen_theta = optim_parameters_list$ln_FishLen_theta[,,,drop = FALSE],
                                ln_FishLen_theta_agg = optim_parameters_list$ln_FishLen_theta_agg,
                                FishLen_corr_pars_agg = optim_parameters_list$FishLen_corr_pars_agg,
                                FishLen_corr_pars = optim_parameters_list$FishLen_corr_pars[,,,,drop = FALSE],

                                # population-specific age composition specifications
                                comp_fishage_pop_like = data$FishAgeComps_pop_LikeType,
                                FishAgeComps_pop_Type = data$FishAgeComps_pop_Type,
                                ISS_FishAgeComps_pop = if(any(data$UseFishAgeComps_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishAgeComps_pop[,,,,,,drop = FALSE] * data$Wt_FishAgeComps_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishAge_pop_theta = optim_parameters_list$ln_FishAge_pop_theta[,,,,drop = FALSE],
                                ln_FishAge_pop_theta_agg = optim_parameters_list$ln_FishAge_pop_theta_agg,
                                FishAge_pop_corr_pars_agg = optim_parameters_list$FishAge_pop_corr_pars_agg,
                                FishAge_pop_corr_pars = optim_parameters_list$FishAge_pop_corr_pars[,,,,,drop = FALSE],

                                # population-specific length composition specifications
                                comp_fishlen_pop_like = data$FishLenComps_pop_LikeType,
                                FishLenComps_pop_Type = data$FishLenComps_pop_Type,
                                ISS_FishLenComps_pop = if(any(data$UseFishLenComps_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishLenComps_pop[,,,,,,drop = FALSE] * data$Wt_FishLenComps_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishLen_pop_theta = optim_parameters_list$ln_FishLen_pop_theta[,,,,drop = FALSE],
                                ln_FishLen_pop_theta_agg = optim_parameters_list$ln_FishLen_pop_theta_agg,
                                FishLen_pop_corr_pars_agg = optim_parameters_list$FishLen_pop_corr_pars_agg,
                                FishLen_pop_corr_pars = optim_parameters_list$FishLen_pop_corr_pars[,,,,,drop = FALSE],

                                # discarded fishery age composition specifications
                                comp_fishage_discard_like = data$FishAgeComps_discard_LikeType,
                                FishAgeComps_discard_Type = data$FishAgeComps_discard_Type,
                                ISS_FishAgeComps_discard = replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard[,,,,,drop = FALSE] * data$Wt_FishAgeComps_discard),
                                ln_FishAge_discard_theta = optim_parameters_list$ln_FishAge_discard_theta[,,,drop = FALSE],
                                ln_FishAge_discard_theta_agg = optim_parameters_list$ln_FishAge_discard_theta_agg,
                                FishAge_discard_corr_pars_agg = optim_parameters_list$FishAge_discard_corr_pars_agg,
                                FishAge_discard_corr_pars = optim_parameters_list$FishAge_discard_corr_pars[,,,,drop = FALSE],

                                # discarded fishery length composition specifications
                                comp_fishlen_discard_like = data$FishLenComps_discard_LikeType,
                                FishLenComps_discard_Type = data$FishLenComps_discard_Type,
                                ISS_FishLenComps_discard = replicate(sim_list$n_sims, data$ISS_FishLenComps_discard[,,,,,drop = FALSE] * data$Wt_FishLenComps_discard),
                                ln_FishLen_discard_theta = optim_parameters_list$ln_FishLen_discard_theta[,,,drop = FALSE],
                                ln_FishLen_discard_theta_agg = optim_parameters_list$ln_FishLen_discard_theta_agg,
                                FishLen_discard_corr_pars_agg = optim_parameters_list$FishLen_discard_corr_pars_agg,
                                FishLen_discard_corr_pars = optim_parameters_list$FishLen_discard_corr_pars[,,,,drop = FALSE],

                                # discarded population-specific age composition specifications
                                comp_fishage_discard_pop_like = data$FishAgeComps_discard_pop_LikeType,
                                FishAgeComps_discard_pop_Type = data$FishAgeComps_discard_pop_Type,
                                ISS_FishAgeComps_discard_pop = if(any(data$UseFishAgeComps_discard_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard_pop[,,,,,,drop = FALSE] * data$Wt_FishAgeComps_discard_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishAge_discard_pop_theta = optim_parameters_list$ln_FishAge_discard_pop_theta[,,,,drop = FALSE],
                                ln_FishAge_discard_pop_theta_agg = optim_parameters_list$ln_FishAge_discard_pop_theta_agg,
                                FishAge_discard_pop_corr_pars_agg = optim_parameters_list$FishAge_discard_pop_corr_pars_agg,
                                FishAge_discard_pop_corr_pars = optim_parameters_list$FishAge_discard_pop_corr_pars[,,,,,drop = FALSE],

                                # discarded population-specific length composition specifications
                                comp_fishlen_discard_pop_like = data$FishLenComps_discard_pop_LikeType,
                                FishLenComps_discard_pop_Type = data$FishLenComps_discard_pop_Type,
                                ISS_FishLenComps_discard_pop = if(any(data$UseFishLenComps_discard_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishLenComps_discard_pop[,,,,,,drop = FALSE] * data$Wt_FishLenComps_discard_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishLen_discard_pop_theta = optim_parameters_list$ln_FishLen_discard_pop_theta[,,,,drop = FALSE],
                                ln_FishLen_discard_pop_theta_agg = optim_parameters_list$ln_FishLen_discard_pop_theta_agg,
                                FishLen_discard_pop_corr_pars_agg = optim_parameters_list$FishLen_discard_pop_corr_pars_agg,
                                FishLen_discard_pop_corr_pars = optim_parameters_list$FishLen_discard_pop_corr_pars[,,,,,drop = FALSE],

                                # conditional age-at-length specifications; absent on models built
                                # before the data source existed, which the defaults leave off
                                comp_fish_caal_like = if(is.null(data$Fish_caal_LikeType)) rep(999, sim_list$n_fish_fleets) else data$Fish_caal_LikeType,
                                Fish_caal_Type = if(is.null(data$Fish_caal_Type)) array(999, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)) else data$Fish_caal_Type,
                                ISS_Fish_caal = if(!is.null(data$UseFish_caal) && any(data$UseFish_caal == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_Fish_caal[,,,,,,drop = FALSE] * data$Wt_Fish_caal)
                                } else NULL,
                                ln_Fish_caal_theta = optim_parameters_list$ln_Fish_caal_theta,
                                ln_Fish_caal_theta_agg = optim_parameters_list$ln_Fish_caal_theta_agg
  )

  # Setup Survey Processes --------------------------------------------------
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = bind_sims(lapply(views$reps, function(rp) rp$srv_sel[,,seq_along(data$years),,,,,drop = FALSE])),
    # length comps selected at length read the fit's selectivity at length
    SrvLenComps_sel = if(is.null(data$srv_len_comp_sel)) rep("age", data$n_srv_fleets) else ifelse(data$srv_len_comp_sel == 1, "length", "age"),
    srv_sel_l_input = if(any(data$srv_len_comp_sel == 1)) bind_sims(lapply(views$reps, function(rp) rp$srv_sel_l[,seq_along(data$years),,,,drop = FALSE])) else NULL,
    srv_q_input = bind_sims(lapply(srv_q_fit, function(q) q$q_mean)),
    # the reported errors, with any estimated part of the sd drawn on top as the fit forms it
    ObsSrvIdx_SE = idx_draw_se(data$ObsSrvIdx_SE, data$Wt_SrvIdx, data$sigmaSrvIdx_form),
    sigmaSrvIdx_form = if(is.null(data$sigmaSrvIdx_form)) 0 else data$sigmaSrvIdx_form,
    ln_sigmaSrvIdx = optim_parameters_list$ln_sigmaSrvIdx,
    sigmaSrvIdx_pop_form = if(is.null(data$sigmaSrvIdx_pop_form)) 0 else data$sigmaSrvIdx_pop_form,
    ln_sigmaSrvIdx_pop = optim_parameters_list$ln_sigmaSrvIdx_pop,
    # the index at age has its own error by age and fleet, so no weight is applied to it
    ln_sigmaSrvIdxAA = unused_at_age_on_obs_ages(optim_parameters_list$ln_sigmaSrvIdxAA, srv_idx_aa_used, 1, n_obs_om, log(0.5)),
    UseSrvIdxAA = unused_at_age_on_obs_ages(data$UseSrvIdxAA, srv_idx_aa_used, 4, n_obs_om),
    use_srv_idx_aa = data$use_srv_idx_aa,
    AgeObsCorr_srv_idx = data$AgeObsCorr_srv_idx,
    trans_rho_srv_idx = if(srv_idx_aa_used) optim_parameters_list$trans_rho_srv_idx,
    trans_rho_srv_idx_year = if(srv_idx_aa_used) optim_parameters_list$trans_rho_srv_idx_year,
    trans_rho_srv_idx_us = if(srv_idx_aa_used) optim_parameters_list$trans_rho_srv_idx_us,
    SrvIdxAA_seas_Type = data$SrvIdxAA_seas_Type,
    UseSrvIdxAA_pop = if(srv_idx_aa_pop_used) data$UseSrvIdxAA_pop,
    ln_sigmaSrvIdxAA_pop = if(srv_idx_aa_pop_used) optim_parameters_list$ln_sigmaSrvIdxAA_pop,
    ObsSrvIdxAA_pop_SE = if(srv_idx_aa_pop_used) data$ObsSrvIdxAA_pop_SE,
    SrvIdxAA_pop_Type = data$SrvIdxAA_pop_Type,
    SrvIdxAA_pop_LikeType = data$SrvIdxAA_pop_LikeType,
    SrvIdxAA_pop_sigma_form = data$SrvIdxAA_pop_sigma_form,
    SrvIdxAA_pop_seas_Type = data$SrvIdxAA_pop_seas_Type,
    AgeObsCorr_srv_idx_pop = data$AgeObsCorr_srv_idx_pop,
    trans_rho_srv_idx_pop = if(srv_idx_aa_pop_used) optim_parameters_list$trans_rho_srv_idx_pop,
    trans_rho_srv_idx_pop_year = if(srv_idx_aa_pop_used) optim_parameters_list$trans_rho_srv_idx_pop_year,
    trans_rho_srv_idx_pop_us = if(srv_idx_aa_pop_used) optim_parameters_list$trans_rho_srv_idx_pop_us,
    ObsSrvIdxAA_SE = unused_at_age_on_obs_ages(data$ObsSrvIdxAA_SE, srv_idx_aa_used, 4, n_obs_om),
    SrvIdxAA_Type = data$SrvIdxAA_Type,
    SrvIdxAA_LikeType = data$SrvIdxAA_LikeType,
    SrvIdxAA_sigma_form = data$SrvIdxAA_sigma_form,
    ObsSrvIdx_pop_SE = if(any(data$UseSrvIdx_pop == 1)) {
      idx_draw_se(data$ObsSrvIdx_pop_SE, data$Wt_SrvIdx_pop, data$sigmaSrvIdx_pop_form)
    } else {
      array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))
    },
    srv_idx_type = data$srv_idx_type,
    srv_idx_ages = data$srv_idx_ages, # ages in each index total
    SrvIdx_LikeType = if(is.null(data$SrvIdx_LikeType)) rep(0, data$n_srv_fleets) else data$SrvIdx_LikeType,
    SrvIdx_Cov = data$SrvIdx_Cov,
    UseSrvIdx = data$UseSrvIdx,
    # data sources the fitted model reports once a year stay annual
    SrvIdx_seas_Type = data$SrvIdx_seas_Type,
    SrvIdx_pop_seas_Type = data$SrvIdx_pop_seas_Type,
    SrvAgeComps_seas_Type = data$SrvAgeComps_seas_Type,
    SrvLenComps_seas_Type = data$SrvLenComps_seas_Type,
    SrvAgeComps_pop_seas_Type = data$SrvAgeComps_pop_seas_Type,
    SrvLenComps_pop_seas_Type = data$SrvLenComps_pop_seas_Type,
    seas_agg_slot = seas_agg_slot_list(data, length(data$years), "srv"), # the season each year total sits in
    t_srv = data$t_srv,

    # survey age composition specifications
    comp_srvage_like = data$SrvAgeComps_LikeType,
    SrvAgeComps_Type = data$SrvAgeComps_Type,
    ISS_SrvAgeComps = replicate(sim_list$n_sims, data$ISS_SrvAgeComps[,,,,,drop = FALSE] * data$Wt_SrvAgeComps),
    ln_SrvAge_theta = optim_parameters_list$ln_SrvAge_theta[,,,drop = FALSE],
    ln_SrvAge_theta_agg = optim_parameters_list$ln_SrvAge_theta_agg,
    SrvAge_corr_pars_agg = optim_parameters_list$SrvAge_corr_pars_agg,
    SrvAge_corr_pars = optim_parameters_list$SrvAge_corr_pars[,,,,drop = FALSE],

    # survey length composition specifications
    comp_srvlen_like = data$SrvLenComps_LikeType,
    SrvLenComps_Type = data$SrvLenComps_Type,
    ISS_SrvLenComps = replicate(sim_list$n_sims, data$ISS_SrvLenComps[,,,,,drop = FALSE] * data$Wt_SrvLenComps),
    ln_SrvLen_theta = optim_parameters_list$ln_SrvLen_theta[,,,drop = FALSE],
    ln_SrvLen_theta_agg = optim_parameters_list$ln_SrvLen_theta_agg,
    SrvLen_corr_pars_agg = optim_parameters_list$SrvLen_corr_pars_agg,
    SrvLen_corr_pars = optim_parameters_list$SrvLen_corr_pars[,,,,drop = FALSE],

    # population-specific age composition specifications
    comp_srvage_pop_like = data$SrvAgeComps_pop_LikeType,
    SrvAgeComps_pop_Type = data$SrvAgeComps_pop_Type,
    ISS_SrvAgeComps_pop = if(any(data$UseSrvAgeComps_pop == 1)) {
      replicate(sim_list$n_sims, data$ISS_SrvAgeComps_pop[,,,,,,drop = FALSE] * data$Wt_SrvAgeComps_pop)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    },
    ln_SrvAge_pop_theta = optim_parameters_list$ln_SrvAge_pop_theta[,,,,drop = FALSE],
    ln_SrvAge_pop_theta_agg = optim_parameters_list$ln_SrvAge_pop_theta_agg,
    SrvAge_pop_corr_pars_agg = optim_parameters_list$SrvAge_pop_corr_pars_agg,
    SrvAge_pop_corr_pars = optim_parameters_list$SrvAge_pop_corr_pars[,,,,,drop = FALSE],

    # population-specific length composition specifications
    comp_srvlen_pop_like = data$SrvLenComps_pop_LikeType,
    SrvLenComps_pop_Type = data$SrvLenComps_pop_Type,
    ISS_SrvLenComps_pop = if(any(data$UseSrvLenComps_pop == 1)) {
      replicate(sim_list$n_sims, data$ISS_SrvLenComps_pop[,,,,,,drop = FALSE] * data$Wt_SrvLenComps_pop)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    },
    ln_SrvLen_pop_theta = optim_parameters_list$ln_SrvLen_pop_theta[,,,,drop = FALSE],
    ln_SrvLen_pop_theta_agg = optim_parameters_list$ln_SrvLen_pop_theta_agg,
    SrvLen_pop_corr_pars_agg = optim_parameters_list$SrvLen_pop_corr_pars_agg,
    SrvLen_pop_corr_pars = optim_parameters_list$SrvLen_pop_corr_pars[,,,,,drop = FALSE],

    # conditional age-at-length specifications
    comp_srv_caal_like = if(is.null(data$Srv_caal_LikeType)) rep(999, sim_list$n_srv_fleets) else data$Srv_caal_LikeType,
    Srv_caal_Type = if(is.null(data$Srv_caal_Type)) array(999, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)) else data$Srv_caal_Type,
    ISS_Srv_caal = if(!is.null(data$UseSrv_caal) && any(data$UseSrv_caal == 1)) {
      replicate(sim_list$n_sims, data$ISS_Srv_caal[,,,,,,drop = FALSE] * data$Wt_Srv_caal)
    } else NULL,
    ln_Srv_caal_theta = optim_parameters_list$ln_Srv_caal_theta,
    ln_Srv_caal_theta_agg = optim_parameters_list$ln_Srv_caal_theta_agg
  )

  # Setup Biological Dynamics -----------------------------------------------
  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list, # simualtion list
    natmort_input = bind_sims(lapply(views$reps, function(rp) truncate_years(expand_natmort_seasons(rp$natmort, data$n_seas), length(data$years)))), # natural mortality
    # derived by the growth module when present, otherwise the data the model was given
    WAA_input = bind_sims(lapply(views$reps, function(rp) (if(is.null(rp$WAA)) data$WAA else rp$WAA)[,,seq_along(data$years),,,,drop = FALSE])), # weight at age
    WAA_fish_input = bind_sims(lapply(views$reps, function(rp) (if(is.null(rp$WAA_fish)) data$WAA_fish else rp$WAA_fish)[,,seq_along(data$years),,,,,drop = FALSE])), # fishery weight at age
    WAA_srv_input = bind_sims(lapply(views$reps, function(rp) (if(is.null(rp$WAA_srv)) data$WAA_srv else rp$WAA_srv)[,,seq_along(data$years),,,,,drop = FALSE])), # survey weight at age
    MatAA_input = replicate(n = sim_list$n_sims, data$MatAA[,,seq_along(data$years),,,,drop = FALSE]), # maturity at age
    AgeingError_input = replicate(n = sim_list$n_sims, data$AgeingError[seq_along(data$years),,,drop = FALSE]), # ageing error
    # fleet-specific ageing error, absent from data lists written before it existed, in which case the operating model falls back on the shared matrix
    AgeingError_fish_input = if(is.null(data$AgeingError_fish)) NULL else replicate(n = sim_list$n_sims, data$AgeingError_fish[seq_along(data$years),,,,drop = FALSE]),
    AgeingError_srv_input = if(is.null(data$AgeingError_srv)) NULL else replicate(n = sim_list$n_sims, data$AgeingError_srv[seq_along(data$years),,,,drop = FALSE]),
    LenBinMap_input = data$LenBinMap, # model length bins onto the recorded ones, as the fit maps them
    CAAL_LenBinMap_input = data$CAAL_LenBinMap, # model length bins each age-at-length row covers
    SizeAgeTrans_input = if(data$fit_lengths == 0 || is.null(data$SizeAgeTrans) || all(is.na(data$SizeAgeTrans))) NULL else replicate(n = sim_list$n_sims, data$SizeAgeTrans[,,seq_along(data$years),,,,,drop = FALSE]),
    # keys per fleet from the growth module, each at its fleet's own timing
    SizeAgeTrans_fish_input = if(is.null(rep$SizeAgeTrans_fish)) NULL else bind_sims(lapply(views$reps, function(rp) rp$SizeAgeTrans_fish[,,seq_along(data$years),,,,,,drop = FALSE])),
    SizeAgeTrans_srv_input = if(is.null(rep$SizeAgeTrans_srv)) NULL else bind_sims(lapply(views$reps, function(rp) rp$SizeAgeTrans_srv[,,seq_along(data$years),,,,,,drop = FALSE])) # size age transition matrix, derived by the growth module when present
  )

  # growth deviations, the fit's under conditional, drawn fresh at each replicate's own parameters under joint
  sim_list <- Setup_Sim_Growth_RE(sim_list, data, optim_parameters_list, rep = rep,
                                  pars_by_sim = if(sim_type == "joint") views$pars,
                                  rep_by_sim = if(sim_type == "joint") views$reps)

  # Movement
  sim_list$Movement <- bind_sims(lapply(views$reps, function(rp) rp$Movement[,,,seq_along(data$years),,,,drop = FALSE]))
  sim_list$sgl_seas_spawning_movement <- bind_sims(lapply(views$reps, function(rp) rp$sgl_seas_spawning_movement[,,,seq_along(data$years),,,drop = FALSE]))
  # Movement / mortality sequencing; absent for models built before this option existed
  sim_list$move_timing <- if(is.null(data$move_timing)) 0 else data$move_timing
  # How the matrix exponential is evaluated
  sim_list$expm_nsub <- if(is.null(data$move_expm_nsub)) 0 else data$move_expm_nsub
  # The instantaneous rate matrix only exists for an estimated CTMC, and is only needed for continuous movement
  sim_list$Mrate <- if(sim_list$move_timing == 2) bind_sims(lapply(views$reps, function(rp) rp$Mrate[,,,seq_along(data$years),,,,drop = FALSE])) else NULL
  # movement deviations, the fit's under conditional, drawn fresh at each replicate's own parameters under joint
  sim_list <- Setup_Sim_Movement(sim_list, data, optim_parameters_list,
                                 pars_by_sim = if(sim_type == "joint") views$pars)

  # selectivity deviations, and F and discard mortality deviations the fit integrates out, the fit's under conditional and
  # drawn fresh at each replicate's own parameters under joint
  sim_list <- Setup_Sim_Fleet_Devs(sim_list, data, optim_parameters_list, random = random,
                                   pars_by_sim = if(sim_type == "joint") views$pars)

  # Setup Recruitment Processes ---------------------------------------------

  # under joint recruitment and the initial ages are drawn fresh from each replicate's penalty, unless a dsem links them
  draw_rec <- sim_type == "joint" && !("rec" %in% data$dsem_declared)

  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    spawn_seas = data$spawn_seas, # spawning season
    do_recruits_move = data$do_recruits_move, # whether recruits move
    t_spawn = data$t_spawn, # spawn timing
    init_age_strc = data$init_age_strc, # initilaizing age structure
    h_input = bind_sims(lapply(views$reps, function(rp) array(rp$h_trans, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs)))), # steepness
    R0_input = {
      # R0 can have time blocks, so the operating model takes the year-by-year value rather than
      # R0's single reference-block value. identical in an unblocked model
      tmp = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
      for(i in seq_len(n_sims)) {
        rp = views$reps[[i]]
        R0_by_yr = if(is.null(rp$R0_yr)) matrix(rp$R0, sim_list$n_pop, sim_list$n_yrs)
                   else matrix(rp$R0_yr, sim_list$n_pop, ncol(rp$R0_yr))[, pmin(1:sim_list$n_yrs, ncol(rp$R0_yr)), drop = FALSE]
        for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) tmp[p,r,,i] = R0_by_yr[p,] * rp$rec_region_prop[p,r]
      }
      tmp
    },
    rinit_input = {
      tmp = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims))
      for(i in seq_len(n_sims)) for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) {
        tmp[p,r,i] = views$reps[[i]]$rinit[p] * views$reps[[i]]$rec_region_prop[p,r]
      }
      tmp
    },
    use_rinit = data$use_rinit,
    sexratio_input = bind_sims(lapply(views$reps, function(rp) rp$sexratio[,,seq_along(data$years),,drop = FALSE])), # sex ratio
    # the fit's own; recruitment and the initial ages are supplied below, drawn under Wt_Rec and Wt_Init_Rec where they are drawn
    ln_sigmaR = optim_parameters_list$ln_sigmaR,
    sigmaR_switch = if(is.null(data$sigmaR_switch)) 1 else data$sigmaR_switch, # first year on the late sigmaR, as the fit reads it
    # under conditional recruitment goes in year by year as the fit has it; under joint drawn deviations go on the stock-recruit curve
    Rec_input = if(sim_type == "conditional") bind_sims(lapply(views$reps, function(rp) rp$Rec[,,seq_along(data$years),drop = FALSE])),
    ln_RecDevs_input = if(draw_rec) {
      bind_sims(lapply(seq_len(n_sims), function(i) rec_devs_past_fit(data, views$pars[[i]], views$reps[[i]], 0, length(data$years), mapping$ln_RecDevs)))
    },
    ln_InitDevs_input = if(draw_rec) {
      bind_sims(lapply(seq_len(n_sims), function(i) init_devs_past_fit(data, views$pars[[i]], views$reps[[i]], mapping$ln_InitDevs)))
    } else bind_sims(lapply(views$pars, function(pr) pr$ln_InitDevs)),
    stray_rate_input = replicate(sim_list$n_sims, data$stray_rate[,seq_along(data$years), drop = FALSE]),
    rec_seas_prop_input = array(
      unlist(lapply(views$reps, function(rp) rp$rec_seas_prop)),
      dim = c(data$n_pop, data$n_seas, sim_list$n_sims)), # seasonal recruitment apportionment

    # Not needed; already specified in Rec_input and ln_InitDevs_input
    recruitment_opt = data$rec_model,
    rec_dd = data$rec_dd,
    init_dd = data$rec_dd,
    rec_lag = data$rec_lag,
    # the per-recruit reference year has to match the fit, or the operating model and the
    # estimation model build S0 from different biology and the self test measures that gap
    SR_ref_yr = if(is.null(data$SR_ref_yr)) 1 else data$SR_ref_yr,
    # the operating model draws recruitment under the process the fit was estimated with
    RecDevs_model = c("iid", "rw", "ar1")[if(is.null(data$RecDevs_model)) 1 else data$RecDevs_model],
    RecDevs_rho = if(is.null(optim_parameters_list$RecDevs_rho)) array(0, dim = c(sim_list$n_pop, sim_list$n_regions))
                  else rho_trans(array(optim_parameters_list$RecDevs_rho, dim = c(sim_list$n_pop, sim_list$n_regions)))
  )

  # Setup Numbers At Age State ----------------------------------------------
  if(isTRUE(data$NAA_re > 0)) {

    # figure out eta from NAA PE
    naa_codes <- c("none", "iid", "1dar1_a", "1dar1_y", "2dar1", "3dcond", "3dmarg")
    ny_state <- dim(optim_parameters_list$ln_NAA)[3]
    draw_eta <- function(rp, pr) {
      pred <- rp$NAA_pred[,,seq_len(ny_state),,,,drop = FALSE]
      out <- array(0, dim = dim(pr$ln_NAA))
      out[pred > 0] <- pr$ln_NAA[pred > 0] - log(pred[pred > 0])
      out
    }
    eta <- bind_sims(lapply(seq_len(n_sims), function(i) draw_eta(views$reps[[i]], views$pars[[i]])))

    sim_list <- Setup_Sim_NAA_state(
      sim_list = sim_list,
      NAA_re = naa_codes[data$NAA_re + 1],
      NAA_re_ages = data$naa_re_ages,
      NAA_re_years = data$naa_re_yrs,
      NAA_re_seasons = data$naa_re_seas,
      naa_eta_input = if(sim_type == "conditional") eta
    )

    # the fit's process, its sd and correlations, and under joint each replicate's own, which its fresh draw reads
    naa_process <- naa_process_from_fit(data, optim_parameters_list)
    sim_list[names(naa_process)] <- naa_process
    if(sim_type == "joint") sim_list$naa_process_by_sim <- lapply(views$pars, function(pr) naa_process_from_fit(data, pr))
  }

  # Catchability Stuff -------------------------------------------------
  sim_list$n_cond_yrs <- n_cond_yrs # the leading years whose catchability, movement and growth deviations are the fit's
  sim_list$ln_fish_q_devs <- bind_sims(lapply(fish_q_fit, function(q) q$devs))
  sim_list$ln_srv_q_devs <- bind_sims(lapply(srv_q_fit, function(q) q$devs))

  # under joint the process the fit penalizes: its forms, each replicate's own sigma and correlation, and only the cells it estimates
  if(sim_type == "joint" && any(c(data$fish_q_model, data$srv_q_model) %in% c(2, 3, 4))) {
    q_forms <- c("none", "iid", "rw", "ar1", "dsem") # the fit's codes, in order
    sim_list <- Setup_Sim_q_devs(
      sim_list = sim_list,
      fish_q_model = q_forms[data$fish_q_model],
      srv_q_model = q_forms[data$srv_q_model],
      sigma_fish_q = bind_sims(lapply(views$pars, function(pr) exp(pr$ln_sigma_fish_q))),
      sigma_srv_q = bind_sims(lapply(views$pars, function(pr) exp(pr$ln_sigma_srv_q))),
      fish_q_rho = bind_sims(lapply(views$pars, function(pr) rho_trans(pr$fish_q_rho))),
      srv_q_rho = bind_sims(lapply(views$pars, function(pr) rho_trans(pr$srv_q_rho))),
      fish_q_rw_init_sigma = if(is.null(data$fish_q_rw_init_sigma)) NA else data$fish_q_rw_init_sigma,
      srv_q_rw_init_sigma = if(is.null(data$srv_q_rw_init_sigma)) NA else data$srv_q_rw_init_sigma
    )
    est_cells <- function(map) if(is.null(map)) NULL else !is.na(map[,fit_yrs,,drop = FALSE]) # a year outside q_re_years or a region with no index keeps the fit's value
    sim_list$fish_q_devs_est <- est_cells(data$map_ln_fish_q_devs)
    sim_list$srv_q_devs_est <- est_cells(data$map_ln_srv_q_devs)
  }

  # Setup DSEM --------------------------------------------------------------
  # conditional keeps the dsem's states at the fit; joint draws them fresh from each replicate's own arrows
  if(!is.null(data$dsem_model)) {
    sim_list <- Setup_Sim_DSEM(sim_list, data, optim_parameters_list, rep = rep, condition_on_fit = sim_type == "conditional",
                               pars_by_sim = if(sim_type == "joint") views$pars else NULL,
                               rep_by_sim = if(sim_type == "joint") views$reps else NULL)
  }

  # Setup Tagging -----------------------------------------------------------
  if(!is.na(sum(data$conv_tagged_fish))) n_tags_rel_input <- apply(data$conv_tagged_fish, 1, sum) else n_tags_rel_input <- NA
  if(exists("conv_tag_release_indicator", data)) conv_tag_release_indicator <- data$conv_tag_release_indicator  else conv_tag_release_indicator <- NA
  conv_tag_fish_reporting_input <- if(!is.null(rep$conv_tag_fish_reporting)) bind_sims(lapply(views$reps, function(rp) rp$conv_tag_fish_reporting)) else NULL

  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list, # simulation list
    conv_tag_max_liberty = data$conv_tag_max_liberty, # maximum tag liberty
    conv_tag_t_tagging = data$conv_tag_t_tagging, # time of tagging
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,  # number of tags to release per event
    conv_tag_release_indicator = conv_tag_release_indicator,  # tag release indicator
    ln_init_conv_tag_mort = optim_parameters_list$ln_init_conv_tag_mort,  # inital tagging mortality
    ln_conv_tag_shed = optim_parameters_list$ln_conv_tag_shed, # chronic tag shedding
    conv_tag_fish_reporting_input = conv_tag_fish_reporting_input, # tag reporting rates
    use_conv_fish_tagging = data$use_conv_fish_tagging, # whether or not tagging is used / simulated
    conv_fish_tag_like = data$conv_fish_tag_like, # tag likelihood
    ln_conv_fish_tag_theta = optim_parameters_list$ln_conv_fish_tag_theta, # tag overdispersion, the fit's
    conv_tag_pop_pool = data$conv_tag_pop_pool, # recaptures the fit pools into one count
    conv_tag_age_pool = data$conv_tag_age_pool,
    conv_tag_sex_pool = data$conv_tag_sex_pool,
    conv_tagged_fish_input = if(!anyNA(data$conv_tagged_fish)) data$conv_tagged_fish * data$Wt_Tagging, # the fit's releases by age, released as they were
    conv_tag_release_platform = if(is.null(data$conv_tag_release_platform)) default_tag_release_platform(conv_tag_release_indicator) else data$conv_tag_release_platform, # a data list built before the setup stored one
    conv_fish_tag_attr = data$conv_fish_tag_attr # tag attributes
  )

  # the refits read these back out of sim_list, so both sides use the same error
  if(isTRUE(perfect_data)) sim_list <- make_data_perfect(sim_list)


  # Run Simulation ----------------------------------------------------------

  # storage for the report values, then the parameters, then the sdreport
  n_store <- length(what) + length(what_par) + 1
  store_res_list <- vector("list", n_store) # get list
  names(store_res_list) <- c(what, what_par, "sd_rep") # name list
  for(j in seq_along(c(what, what_par))) store_res_list[[j]] <- vector("list", n_sims) # stick in n_sims lists into storage
  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = output_path) # get simulated datasets

  # simulated comps on other bins than the data would be poured into the wrong cells below, so refuse them
  comp_sources <- c("FishAgeComps", "FishAgeComps_pop", "FishAgeComps_discard", "FishAgeComps_discard_pop", "SrvAgeComps", "SrvAgeComps_pop")
  if(data$fit_lengths != 0) comp_sources <- c(comp_sources, "FishLenComps", "FishLenComps_pop", "FishLenComps_discard", "FishLenComps_discard_pop", "SrvLenComps", "SrvLenComps_pop", "Fish_caal", "Srv_caal")
  for(comp_source in comp_sources) {
    if(!any(data[[paste0("Use", comp_source)]] == 1)) next # not fit, so not written back
    sim_dims <- dim(sim_obj[[paste0("Obs", comp_source)]])
    data_dims <- dim(data[[paste0("Obs", comp_source)]])
    if(!identical(as.numeric(sim_dims[-length(sim_dims)]), as.numeric(data_dims))) {
      stop("simulation_self_test: the operating model drew Obs", comp_source, " as ", paste(sim_dims[-length(sim_dims)], collapse = " x "),
           " but the data are ", paste(data_dims, collapse = " x "), ". The observed age or length bins of the two do not match.")
    }
  } # end comp_source loop

  if(do_par == FALSE) {

    for(i in 1:n_sims) {

      tryCatch({

        # set up data stuff
        tmp_data <- self_test_priors(data, prior_means, views$pars[[i]], views$reps[[i]])
        tmp_pars <- parameters
        tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx))
        tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx))
        tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch))
        tmp_data$ObsCatch[is.na(data$ObsCatch)] <- NA # a catch the fit had missing stays missing, which keeps its fleet fishing
        tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps))
        tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps))
        if(tmp_data$fit_lengths != 0) {
          tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps))
          tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps))
        }

        # population-specific observations
        if(any(tmp_data$UseFishIdx_pop == 1)) {
          tmp_data$ObsFishIdx_pop <- array(sim_obj$ObsFishIdx_pop[,,,,,i], dim = dim(tmp_data$ObsFishIdx_pop))
        }
        if(any(tmp_data$UseSrvIdx_pop == 1)) {
          tmp_data$ObsSrvIdx_pop <- array(sim_obj$ObsSrvIdx_pop[,,,,,i], dim = dim(tmp_data$ObsSrvIdx_pop))
        }
        if(any(tmp_data$UseFishAgeComps_pop == 1)) {
          tmp_data$ObsFishAgeComps_pop <- array(sim_obj$ObsFishAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_pop))
        }
        if(any(tmp_data$UseSrvAgeComps_pop == 1)) {
          tmp_data$ObsSrvAgeComps_pop <- array(sim_obj$ObsSrvAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps_pop))
        }
        if(tmp_data$fit_lengths != 0) {
          if(any(tmp_data$UseFishLenComps_pop == 1)) {
            tmp_data$ObsFishLenComps_pop <- array(sim_obj$ObsFishLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_pop))
          }
          if(any(tmp_data$UseSrvLenComps_pop == 1)) {
            tmp_data$ObsSrvLenComps_pop <- array(sim_obj$ObsSrvLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps_pop))
          }
        }
        if(any(tmp_data$UseCatch_pop == 1)) {
          tmp_data$ObsCatch_pop <- array(sim_obj$ObsCatch_pop[,,,,,i], dim = dim(tmp_data$ObsCatch_pop))
        }

        # set up discarding stuff
        if(any(tmp_data$UseDiscard == 1)) tmp_data$ObsDiscard <- array(sim_obj$ObsDiscard[,,,,i], dim = dim(tmp_data$ObsDiscard))
        if(any(tmp_data$UseDiscard_pop == 1)) tmp_data$ObsDiscard_pop <- array(sim_obj$ObsDiscard_pop[,,,,,i], dim = dim(tmp_data$ObsDiscard_pop))
        if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ObsFishAgeComps_discard <- array(sim_obj$ObsFishAgeComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard))
        if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ObsFishLenComps_discard <- array(sim_obj$ObsFishLenComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard))
        if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ObsFishAgeComps_discard_pop <- array(sim_obj$ObsFishAgeComps_discard_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard_pop))
        if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ObsFishLenComps_discard_pop <- array(sim_obj$ObsFishLenComps_discard_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard_pop))

        # setup tagging data stuff if tagging is done
        if(any(tmp_data$use_conv_fish_tagging == 1)) {
          tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,,i], dim = dim(tmp_data$conv_tagged_fish))
          tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap))
          tmp_data$conv_tag_release_indicator <- sim_obj$conv_tag_release_indicator
        }

        # reset weights
        # the data weights go to one, the operating model having drawn at the error they imply; the penalty weights
        # (Wt_Rec, Wt_Init_Rec, Wt_F, Wt_D, *_pe_wt) stay, the operating model having drawn from the weighted penalty's density
        tmp_data$Wt_Tagging <- 1
        tmp_data$Wt_Catch[] <- 1
        tmp_data$Wt_Discard[] <- 1
        tmp_data$Wt_FishAgeComps[] <- 1
        tmp_data$Wt_FishAgeComps_discard[] <- 1
        tmp_data$Wt_FishIdx[] <- 1
        tmp_data$Wt_FishLenComps[] <- 1
        tmp_data$Wt_FishLenComps_discard[] <- 1
        tmp_data$Wt_SrvAgeComps[] <- 1
        tmp_data$Wt_SrvIdx[] <- 1
        tmp_data$Wt_SrvLenComps[] <- 1
        tmp_data$Wt_Catch_pop[] <- 1
        tmp_data$Wt_Discard_pop[] <- 1
        tmp_data$Wt_FishIdx_pop[] <- 1
        tmp_data$Wt_SrvIdx_pop[] <- 1
        tmp_data$Wt_FishAgeComps_pop[] <- 1
        tmp_data$Wt_FishAgeComps_discard_pop[] <- 1
        tmp_data$Wt_SrvAgeComps_pop[] <- 1
        tmp_data$Wt_FishLenComps_pop[] <- 1
        tmp_data$Wt_FishLenComps_discard_pop[] <- 1
        tmp_data$Wt_SrvLenComps_pop[] <- 1

        # input simulated uncertainty
        tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
        tmp_pars$ln_sigmaC_pop[] <- sim_list$ln_sigmaC_pop
        tmp_pars$ln_sigmaD[] <- sim_list$ln_sigmaD
        tmp_pars$ln_sigmaD_pop[] <- sim_list$ln_sigmaD_pop
        tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
        tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
        # only a data source the model fits takes the draws, so an unused one keeps its own shape
        if(!is.null(tmp_data$ObsCatchAA) && !is.null(sim_obj$ObsCatchAA) && any(tmp_data$UseCatchAA == 1))
          tmp_data$ObsCatchAA[] <- sim_obj$ObsCatchAA[,,,,,,i]
        if(!is.null(tmp_data$ObsDiscardAA) && !is.null(sim_obj$ObsDiscardAA) && any(tmp_data$UseDiscardAA == 1))
          tmp_data$ObsDiscardAA[] <- sim_obj$ObsDiscardAA[,,,,,,i]
        if(!is.null(tmp_data$ObsSrvIdxAA) && !is.null(sim_obj$ObsSrvIdxAA) && any(tmp_data$UseSrvIdxAA == 1))
          tmp_data$ObsSrvIdxAA[] <- sim_obj$ObsSrvIdxAA[,,,,,,i]
        # population-specific at-age data sources, wherever the fit observes them
        for(name in c("CatchAA_pop", "DiscardAA_pop", "SrvIdxAA_pop")) {
          if(is.null(tmp_data[[paste0("Obs", name)]]) || !any(tmp_data[[paste0("Use", name)]] == 1)) next
          tmp_data[[paste0("Obs", name)]][] <- sim_obj[[paste0("Obs", name)]][,,,,,,,i]
        } # end name loop
        if(!is.null(tmp_pars$ln_sigmaCAA)) tmp_pars$ln_sigmaCAA[] <- parameters$ln_sigmaCAA
        if(!is.null(tmp_pars$ln_sigmaSrvIdxAA)) tmp_pars$ln_sigmaSrvIdxAA[] <- parameters$ln_sigmaSrvIdxAA
        if(!is.null(tmp_pars$ln_sigmaFishIdx)) tmp_pars$ln_sigmaFishIdx[] <- parameters$ln_sigmaFishIdx
        if(!is.null(tmp_pars$ln_sigmaSrvIdx)) tmp_pars$ln_sigmaSrvIdx[] <- parameters$ln_sigmaSrvIdx
        tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,,i]
        tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,,i]
        if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ISS_FishAgeComps_discard[] <- sim_list$ISS_FishAgeComps_discard[,,,,,i]
        if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ISS_FishLenComps_discard[] <- sim_list$ISS_FishLenComps_discard[,,,,,i]
        tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,,i]
        tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,,i]
        if(any(tmp_data$UseFishIdx_pop == 1)) tmp_data$ObsFishIdx_pop_SE[] <- sim_list$ObsFishIdx_pop_SE
        if(any(tmp_data$UseSrvIdx_pop == 1)) tmp_data$ObsSrvIdx_pop_SE[] <- sim_list$ObsSrvIdx_pop_SE
        if(any(tmp_data$UseFishAgeComps_pop == 1)) tmp_data$ISS_FishAgeComps_pop[] <- sim_list$ISS_FishAgeComps_pop[,,,,,,i]
        if(any(tmp_data$UseFishLenComps_pop == 1)) tmp_data$ISS_FishLenComps_pop[] <- sim_list$ISS_FishLenComps_pop[,,,,,,i]
        if(any(tmp_data$UseSrvAgeComps_pop == 1)) tmp_data$ISS_SrvAgeComps_pop[] <- sim_list$ISS_SrvAgeComps_pop[,,,,,,i]
        if(any(tmp_data$UseSrvLenComps_pop == 1)) tmp_data$ISS_SrvLenComps_pop[] <- sim_list$ISS_SrvLenComps_pop[,,,,,,i]
        if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ISS_FishAgeComps_discard_pop[] <- sim_list$ISS_FishAgeComps_discard_pop[,,,,,,i]
        if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ISS_FishLenComps_discard_pop[] <- sim_list$ISS_FishLenComps_discard_pop[,,,,,,i]

        # conditional age-at-length observations
        if(!is.null(tmp_data$UseFish_caal) && any(tmp_data$UseFish_caal == 1)) {
          tmp_data$ObsFish_caal <- array(sim_obj$ObsFish_caal[,,,,,,,i], dim = dim(tmp_data$ObsFish_caal))
          tmp_data$ISS_Fish_caal[] <- sim_list$ISS_Fish_caal[,,,,,,i]
          tmp_data$Wt_Fish_caal[] <- 1
        }
        if(!is.null(tmp_data$UseSrv_caal) && any(tmp_data$UseSrv_caal == 1)) {
          tmp_data$ObsSrv_caal <- array(sim_obj$ObsSrv_caal[,,,,,,,i], dim = dim(tmp_data$ObsSrv_caal))
          tmp_data$ISS_Srv_caal[] <- sim_list$ISS_Srv_caal[,,,,,,i]
          tmp_data$Wt_Srv_caal[] <- 1
        }

        # do dsem covariate stuff here
        if(!is.null(tmp_data$dsem_model)) {
          tmp_data$dsem_cov_obs <- array(sim_obj$dsem_cov_obs_sim[,,i], dim = dim(tmp_data$dsem_cov_obs))
          for(k in seq_along(tmp_data$dsem_cov_var_idx)) {
            if(tmp_data$dsem_cov_family[k] != 0) next
            seen <- !is.na(tmp_data$dsem_cov_obs[,k])
            tmp_pars$dsem_x[seen,tmp_data$dsem_cov_var_idx[k]] <- tmp_data$dsem_cov_obs[seen,k]
          } # end k loop
        }

        # update setup stuff if needed
        tmp_data <- resync_fitted_blocks(tmp_data)

        # Fit model
        fit_i <- fit_model(
          data = tmp_data,
          parameters = tmp_pars,
          mapping = mapping,
          random = random,
          newton_loops = newton_loops,
          silent = TRUE
        )

        # Populate results into store list
        for(j in seq_along(what)) store_res_list[[j]][[i]] <- fit_i$rep[[what[j]]]
        if(length(what_par) > 0) {
          est_pars <- fit_i$env$parList(par = fit_i$env$last.par.best) # name the argument, the bare first one is fixed effects only
          for(j in seq_along(what_par)) store_res_list[[length(what) + j]][[i]] <- est_pars[[what_par[j]]]
        }

        if(do_sdrep == TRUE) {
          tryCatch({
            fit_i$sd_rep <- RTMB::sdreport(fit_i)
            store_res_list[[n_store]][[i]] <- fit_i$sd_rep # input sd report
          }, error = function(e) {
            store_res_list[[n_store]][[i]] <- NA
          })
        }

      }, error = function(e) {
        # Skip failed simulations, saying why
        warning(sprintf("simulation %d failed: %s", i, conditionMessage(e)), call. = FALSE)
        for(j in seq_along(c(what, what_par))) store_res_list[[j]][[i]] <- NA
        if(do_sdrep == TRUE) store_res_list[[n_store]][[i]] <- NA
      })

    } # end i loop

    # Convert result lists to array
    for(j in seq_along(c(what, what_par))) store_res_list[[j]] <- simplify2array(fill_failed_replicates(store_res_list[[j]]))

  } # not doing parallelization

  if(do_par == TRUE) {

    # initialize cores
    if(is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    options(future.globals.maxSize = 5e3 * 1024^2)  # increase parrlalel size
    future::plan(future::multisession, workers = n_cores) # set up cores
    progressr::with_progress({
      p <- progressr::progressor(along = 1:n_sims) # progress bar

      sim_results <- future.apply::future_lapply(1:n_sims, function(i) {

        tryCatch({

          # set up data stuff
          tmp_data <- self_test_priors(data, prior_means, views$pars[[i]], views$reps[[i]])
          tmp_pars <- parameters
          tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx))
          tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx))
          tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch))
          tmp_data$ObsCatch[is.na(data$ObsCatch)] <- NA # a catch the fit had missing stays missing, which keeps its fleet fishing
          tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps))
          tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps))
          if(tmp_data$fit_lengths != 0) {
            tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps))
            tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps))
          }

          # population-specific observations
          if(any(tmp_data$UseFishIdx_pop == 1)) {
            tmp_data$ObsFishIdx_pop <- array(sim_obj$ObsFishIdx_pop[,,,,,i], dim = dim(tmp_data$ObsFishIdx_pop))
          }
          if(any(tmp_data$UseSrvIdx_pop == 1)) {
            tmp_data$ObsSrvIdx_pop <- array(sim_obj$ObsSrvIdx_pop[,,,,,i], dim = dim(tmp_data$ObsSrvIdx_pop))
          }
          if(any(tmp_data$UseFishAgeComps_pop == 1)) {
            tmp_data$ObsFishAgeComps_pop <- array(sim_obj$ObsFishAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_pop))
          }
          if(any(tmp_data$UseSrvAgeComps_pop == 1)) {
            tmp_data$ObsSrvAgeComps_pop <- array(sim_obj$ObsSrvAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps_pop))
          }
          if(tmp_data$fit_lengths != 0) {
            if(any(tmp_data$UseFishLenComps_pop == 1)) {
              tmp_data$ObsFishLenComps_pop <- array(sim_obj$ObsFishLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_pop))
            }
            if(any(tmp_data$UseSrvLenComps_pop == 1)) {
              tmp_data$ObsSrvLenComps_pop <- array(sim_obj$ObsSrvLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps_pop))
            }
          }
          if(any(tmp_data$UseCatch_pop == 1)) {
            tmp_data$ObsCatch_pop <- array(sim_obj$ObsCatch_pop[,,,,,i], dim = dim(tmp_data$ObsCatch_pop))
          }

          # set up discarding stuff
          if(any(tmp_data$UseDiscard == 1)) tmp_data$ObsDiscard <- array(sim_obj$ObsDiscard[,,,,i], dim = dim(tmp_data$ObsDiscard))
          if(any(tmp_data$UseDiscard_pop == 1)) tmp_data$ObsDiscard_pop <- array(sim_obj$ObsDiscard_pop[,,,,,i], dim = dim(tmp_data$ObsDiscard_pop))
          if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ObsFishAgeComps_discard <- array(sim_obj$ObsFishAgeComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard))
          if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ObsFishLenComps_discard <- array(sim_obj$ObsFishLenComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard))
          if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ObsFishAgeComps_discard_pop <- array(sim_obj$ObsFishAgeComps_discard_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard_pop))
          if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ObsFishLenComps_discard_pop <- array(sim_obj$ObsFishLenComps_discard_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard_pop))

          # setup tagging data stuff if tagging is done
          if(any(tmp_data$use_conv_fish_tagging == 1)) {
            tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,,i], dim = dim(tmp_data$conv_tagged_fish))
            tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap))
            tmp_data$conv_tag_release_indicator <- sim_obj$conv_tag_release_indicator
          }

          # reset weights
          # the data weights go to one, the operating model having drawn at the error they imply; the penalty weights
          # (Wt_Rec, Wt_Init_Rec, Wt_F, Wt_D, *_pe_wt) stay, the operating model having drawn from the weighted penalty's density
          tmp_data$Wt_Tagging <- 1
          tmp_data$Wt_Catch[] <- 1
          tmp_data$Wt_Discard[] <- 1
          tmp_data$Wt_FishAgeComps[] <- 1
          tmp_data$Wt_FishAgeComps_discard[] <- 1
          tmp_data$Wt_FishIdx[] <- 1
          tmp_data$Wt_FishLenComps[] <- 1
          tmp_data$Wt_FishLenComps_discard[] <- 1
          tmp_data$Wt_SrvAgeComps[] <- 1
          tmp_data$Wt_SrvIdx[] <- 1
          tmp_data$Wt_SrvLenComps[] <- 1
          tmp_data$Wt_Catch_pop[] <- 1
          tmp_data$Wt_Discard_pop[] <- 1
          tmp_data$Wt_FishIdx_pop[] <- 1
          tmp_data$Wt_SrvIdx_pop[] <- 1
          tmp_data$Wt_FishAgeComps_pop[] <- 1
          tmp_data$Wt_FishAgeComps_discard_pop[] <- 1
          tmp_data$Wt_SrvAgeComps_pop[] <- 1
          tmp_data$Wt_FishLenComps_pop[] <- 1
          tmp_data$Wt_FishLenComps_discard_pop[] <- 1
          tmp_data$Wt_SrvLenComps_pop[] <- 1

          # input simulated uncertainty
          tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
          tmp_pars$ln_sigmaC_pop[] <- sim_list$ln_sigmaC_pop
          tmp_pars$ln_sigmaD[] <- sim_list$ln_sigmaD
          tmp_pars$ln_sigmaD_pop[] <- sim_list$ln_sigmaD_pop
          tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
          tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
          # only a data source the model fits takes the draws, so an unused one keeps its own shape
          if(!is.null(tmp_data$ObsCatchAA) && !is.null(sim_obj$ObsCatchAA) && any(tmp_data$UseCatchAA == 1))
            tmp_data$ObsCatchAA[] <- sim_obj$ObsCatchAA[,,,,,,i]
          if(!is.null(tmp_data$ObsDiscardAA) && !is.null(sim_obj$ObsDiscardAA) && any(tmp_data$UseDiscardAA == 1))
            tmp_data$ObsDiscardAA[] <- sim_obj$ObsDiscardAA[,,,,,,i]
          if(!is.null(tmp_data$ObsSrvIdxAA) && !is.null(sim_obj$ObsSrvIdxAA) && any(tmp_data$UseSrvIdxAA == 1))
            tmp_data$ObsSrvIdxAA[] <- sim_obj$ObsSrvIdxAA[,,,,,,i]
          # population-specific at-age data sources, wherever the fit observes them
          for(name in c("CatchAA_pop", "DiscardAA_pop", "SrvIdxAA_pop")) {
            if(is.null(tmp_data[[paste0("Obs", name)]]) || !any(tmp_data[[paste0("Use", name)]] == 1)) next
            tmp_data[[paste0("Obs", name)]][] <- sim_obj[[paste0("Obs", name)]][,,,,,,,i]
          } # end name loop
          if(!is.null(tmp_pars$ln_sigmaCAA)) tmp_pars$ln_sigmaCAA[] <- parameters$ln_sigmaCAA
          if(!is.null(tmp_pars$ln_sigmaSrvIdxAA)) tmp_pars$ln_sigmaSrvIdxAA[] <- parameters$ln_sigmaSrvIdxAA
          if(!is.null(tmp_pars$ln_sigmaFishIdx)) tmp_pars$ln_sigmaFishIdx[] <- parameters$ln_sigmaFishIdx
          if(!is.null(tmp_pars$ln_sigmaSrvIdx)) tmp_pars$ln_sigmaSrvIdx[] <- parameters$ln_sigmaSrvIdx
          tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,,i]
          tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,,i]
          if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ISS_FishAgeComps_discard[] <- sim_list$ISS_FishAgeComps_discard[,,,,,i]
          if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ISS_FishLenComps_discard[] <- sim_list$ISS_FishLenComps_discard[,,,,,i]
          tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,,i]
          tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,,i]
          if(any(tmp_data$UseFishIdx_pop == 1)) tmp_data$ObsFishIdx_pop_SE[] <- sim_list$ObsFishIdx_pop_SE
          if(any(tmp_data$UseSrvIdx_pop == 1)) tmp_data$ObsSrvIdx_pop_SE[] <- sim_list$ObsSrvIdx_pop_SE
          if(any(tmp_data$UseFishAgeComps_pop == 1)) tmp_data$ISS_FishAgeComps_pop[] <- sim_list$ISS_FishAgeComps_pop[,,,,,,i]
          if(any(tmp_data$UseFishLenComps_pop == 1)) tmp_data$ISS_FishLenComps_pop[] <- sim_list$ISS_FishLenComps_pop[,,,,,,i]
          if(any(tmp_data$UseSrvAgeComps_pop == 1)) tmp_data$ISS_SrvAgeComps_pop[] <- sim_list$ISS_SrvAgeComps_pop[,,,,,,i]
          if(any(tmp_data$UseSrvLenComps_pop == 1)) tmp_data$ISS_SrvLenComps_pop[] <- sim_list$ISS_SrvLenComps_pop[,,,,,,i]
          if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ISS_FishAgeComps_discard_pop[] <- sim_list$ISS_FishAgeComps_discard_pop[,,,,,,i]
          if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ISS_FishLenComps_discard_pop[] <- sim_list$ISS_FishLenComps_discard_pop[,,,,,,i]

          # conditional age-at-length observations
          if(!is.null(tmp_data$UseFish_caal) && any(tmp_data$UseFish_caal == 1)) {
            tmp_data$ObsFish_caal <- array(sim_obj$ObsFish_caal[,,,,,,,i], dim = dim(tmp_data$ObsFish_caal))
            tmp_data$ISS_Fish_caal[] <- sim_list$ISS_Fish_caal[,,,,,,i]
            tmp_data$Wt_Fish_caal[] <- 1
          }
          if(!is.null(tmp_data$UseSrv_caal) && any(tmp_data$UseSrv_caal == 1)) {
            tmp_data$ObsSrv_caal <- array(sim_obj$ObsSrv_caal[,,,,,,,i], dim = dim(tmp_data$ObsSrv_caal))
            tmp_data$ISS_Srv_caal[] <- sim_list$ISS_Srv_caal[,,,,,,i]
            tmp_data$Wt_Srv_caal[] <- 1
          }

          # dsem covariates, as in the single-model path
          if(!is.null(tmp_data$dsem_model)) {
            tmp_data$dsem_cov_obs <- array(sim_obj$dsem_cov_obs_sim[,,i], dim = dim(tmp_data$dsem_cov_obs))
            for(k in seq_along(tmp_data$dsem_cov_var_idx)) {
              if(tmp_data$dsem_cov_family[k] != 0) next
              seen <- !is.na(tmp_data$dsem_cov_obs[,k])
              tmp_pars$dsem_x[seen,tmp_data$dsem_cov_var_idx[k]] <- tmp_data$dsem_cov_obs[seen,k]
            } # end k loop
          }

          # see the note at the single-model path above
          tmp_data <- resync_fitted_blocks(tmp_data)

          # Fit model
          fit_i <- fit_model(
            data = tmp_data,
            parameters = tmp_pars,
            mapping = mapping,
            random = random,
            newton_loops = newton_loops,
            silent = TRUE
          )

          # Extract what we need and return
          result <- list()
          for(j in seq_along(what)) result[[what[j]]] <- fit_i$rep[[what[j]]]
          if(length(what_par) > 0) {
            est_pars <- fit_i$env$parList(par = fit_i$env$last.par.best) # name the argument, the bare first one is fixed effects only
            for(j in seq_along(what_par)) result[[what_par[j]]] <- est_pars[[what_par[j]]]
          }

          if(do_sdrep == TRUE) {
            tryCatch({
              fit_i$sd_rep <- RTMB::sdreport(fit_i) # get sdreport
              result[[n_store]] <- fit_i$sd_rep # input sd report
            }, error = function(e) {
              result[[n_store]] <- NA
            })
          }

          p() # update progress
          return(result)

        }, error = function(e) {
          # Skip failed simulations, saying why
          warning(sprintf("simulation %d failed: %s", i, conditionMessage(e)), call. = FALSE)
          result <- list()
          for(j in seq_along(what)) result[[what[j]]] <- NA
          for(j in seq_along(what_par)) result[[what_par[j]]] <- NA
          if(do_sdrep == TRUE) result[[n_store]] <- NA

          p() # update progress
          return(result)
        })

      }, future.seed = TRUE)

      future::plan(future::sequential)  # Reset
    })

    # Populate results from parallel run
    for(i in 1:n_sims) for(j in seq_along(c(what, what_par))) store_res_list[[j]][[i]] <- sim_results[[i]][[c(what, what_par)[j]]]
    if(do_sdrep == TRUE) for(i in 1:n_sims) store_res_list[[n_store]][[i]] <- sim_results[[i]][[n_store]]
    for(j in seq_along(c(what, what_par))) store_res_list[[j]] <- simplify2array(fill_failed_replicates(store_res_list[[j]]))  # Convert lists to array
  }

  # the values the OM ran on, its own wherever it holds them; anything it does not hold is the fit's, or under
  # joint the replicate's draw
  # observation error, composition, at-age correlation and tag loss parameters have no replicate dim in the OM,
  # so every replicate ran at the fit's
  # a catch or discard sd the fit weighted was drawn, and is refit, at the sd over the root of the weight
  fit_valued <- grepl("^ln_sigma(C|D|CAA|DAA|SrvIdxAA|FishIdx|SrvIdx)(_pop)?$|_theta(_agg)?$|_corr_pars(_agg)?$|^trans_rho_|^ln_conv_tag_shed$|^ln_init_conv_tag_mort$", what_par)
  deweighted <- what_par %in% c("ln_sigmaC", "ln_sigmaC_pop", "ln_sigmaD", "ln_sigmaD_pop")
  fit_value <- function(w, j) {
    if(deweighted[j] && !is.null(sim_list[[w]]) && identical(dim(sim_list[[w]]), dim(views$fit_pars[[w]]))) return(sim_list[[w]])
    views$fit_pars[[w]]
  }
  store_res_list$truth <- c(
    stats::setNames(lapply(what, function(w) simplify2array(lapply(seq_len(n_sims), function(i) om_truth(sim_obj[[w]], views$reps[[i]][[w]], i, n_sims)))), what),
    stats::setNames(lapply(seq_along(what_par), function(j) {
      w <- what_par[j]
      simplify2array(lapply(seq_len(n_sims), function(i) om_truth(sim_obj[[w]], if(fit_valued[j]) fit_value(w, j) else views$pars[[i]][[w]], i, n_sims)))
    }), what_par)
  )

  # perfect data were drawn at make_data_perfect's sd rather than the fit's, so that is the true observation sd
  if(isTRUE(perfect_data)) {
    perfect_sd_names <- c("ln_sigmaC", "ln_sigmaC_pop", "ln_sigmaD", "ln_sigmaD_pop",
                          "ln_sigmaCAA", "ln_sigmaDAA", "ln_sigmaSrvIdxAA",
                          "ln_sigmaCAA_pop", "ln_sigmaDAA_pop", "ln_sigmaSrvIdxAA_pop",
                          "ln_sigmaFishIdx", "ln_sigmaFishIdx_pop", "ln_sigmaSrvIdx", "ln_sigmaSrvIdx_pop")
    for(w in intersect(what_par, perfect_sd_names)) store_res_list$truth[[w]][] <- log(1e-3)
  }

  return(store_res_list)

}

#' Extract simulation outputs into SPoRC estimation model format
#'
#' Subsets and reshapes the biological, tagging, fishery and survey arrays of a
#' simulation environment to years \code{1:y} and replicate \code{sim}, giving a
#' list ready for \code{\link{Setup_Mod_Biologicals}},
#' \code{\link{Setup_Mod_Catch_and_F}}, \code{\link{Setup_Mod_SrvIdx_and_Comps}}
#' and \code{\link{Setup_Mod_Tagging}}. The \code{Use*} flags are derived from the
#' extracted observations, 1 where a value is present and positive.
#'
#' Population-specific arrays are extracted when \code{sim_env} holds them, with
#' their own flags derived the same way. The length composition outputs and
#' \code{SizeAgeTrans} are \code{NULL} when the environment holds no size-age
#' transition matrix, and the tagging outputs are \code{NULL} under
#' \code{use_conv_fish_tagging = 0}; otherwise only cohorts released in
#' \code{1:y} are kept.
#'
#' @param sim_env Simulation environment or list, from
#'   \code{\link{Simulate_Pop_Static}} or \code{\link{Setup_sim_env}}, holding the
#'   operating model arrays.
#' @param y Integer. Last year to include; years \code{1:y} are kept.
#' @param sim Integer. Simulation replicate to extract.
#'
#' @return Named list, every array with \code{y} in its year dim. The
#'   biological elements are \code{WAA} \code{[n_pop x n_regions x y x n_seas x
#'   n_ages x n_sexes]}, \code{WAA_fish} and \code{WAA_srv} with a trailing fleet
#'   dim, \code{MatAA}, \code{SizeAgeTrans}, \code{AgeingError} \code{[y x n_ages x
#'   n_obs_ages]}, whose columns are the observed ages the model ages are read
#'   onto, and \code{AgeingError_fish} and \code{AgeingError_srv}, \code{NULL} when
#'   the fleets share one matrix. \code{LenBinMap} is the map onto the length bins
#'   the length compositions are recorded on, \code{NULL} when those are the
#'   model's bins, and goes to \code{\link{Setup_Mod_Biologicals}} with the rest,
#'   as does \code{CAAL_LenBinMap}, the model bins each age-at-length row covers.
#'   The tagging elements are
#'   \code{use_conv_fish_tagging}, \code{conv_tag_release_indicator},
#'   \code{obs_conv_tag_fish_recap}, \code{conv_tagged_fish},
#'   \code{conv_tagged_fish_attr} and \code{n_tag_cohorts}.
#'
#'   Each fishery and survey data source contributes its observations, its
#'   observation error or input sample sizes, and its use flag: \code{ObsCatch}
#'   with \code{ln_sigmaC} and \code{UseCatch}, \code{ObsDiscard} with
#'   \code{ln_sigmaD} and \code{UseDiscard}, \code{ObsFishIdx} and
#'   \code{ObsSrvIdx} with their \code{_SE} and use arrays, and the four
#'   composition streams (fishery and survey age and length) with their
#'   \code{ISS_*} and \code{Use*} arrays. Each has a \code{_pop} counterpart, and
#'   the fishery compositions also have \code{_discard} and \code{_discard_pop}
#'   ones. The length composition elements and their input sample sizes are
#'   \code{NULL} when no size-age transition matrix is present.
#'
#'   The at-age data sources, \code{ObsCatchAA}, \code{ObsDiscardAA} and
#'   \code{ObsSrvIdxAA} and their \code{_pop} counterparts, come with their
#'   \code{_SE} arrays and the operating model's own use flags, since a normal
#'   likelihood can draw a value at or below zero. Each is \code{NULL} when the
#'   operating model draws none of it.
#'
#' @seealso \code{\link{Setup_Mod_Biologicals}},
#'   \code{\link{Setup_Mod_Catch_and_F}},
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}, \code{\link{Setup_Mod_Tagging}},
#'   \code{\link{Simulate_Pop_Static}}, \code{\link{Setup_sim_env}}
#'
#' @export simulation_data_to_SPoRC
#' @family Simulation Utilities
simulation_data_to_SPoRC <- function(sim_env,
                                     y,
                                     sim) {

  # Biologicals
  WAA <- array(sim_env$WAA[,,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
  WAA_fish <- array(sim_env$WAA_fish[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets))
  WAA_srv <- array(sim_env$WAA_srv[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_srv_fleets))
  MatAA <- array(sim_env$MatAA[,,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
  SizeAgeTrans <- if(!is.null(sim_env$SizeAgeTrans)) {
    array(sim_env$SizeAgeTrans[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_ages, sim_env$n_sexes))
  } else NULL
  AgeingError <- array(sim_env$AgeingError[1:y,,,sim, drop = FALSE],
                       dim = c(length(1:y), dim(sim_env$AgeingError)[2], dim(sim_env$AgeingError)[3]))
  n_obs_lens <- if(is.null(sim_env$n_obs_lens)) sim_env$n_lens else sim_env$n_obs_lens # length bins the comps are recorded on
  n_caal_lens <- if(is.null(sim_env$n_caal_lens)) sim_env$n_lens else sim_env$n_caal_lens # age-at-length rows
  # the fleet-specific matrices the estimation model reads, peeled to the same years
  AgeingError_fish <- if(is.null(sim_env$AgeingError_fish)) NULL else {
    array(sim_env$AgeingError_fish[1:y,,,,sim, drop = FALSE], dim = c(length(1:y), dim(sim_env$AgeingError_fish)[2:4]))
  }
  AgeingError_srv <- if(is.null(sim_env$AgeingError_srv)) NULL else {
    array(sim_env$AgeingError_srv[1:y,,,,sim, drop = FALSE], dim = c(length(1:y), dim(sim_env$AgeingError_srv)[2:4]))
  }

  # Tagging
  if(sim_env$use_conv_fish_tagging == 1) {
    keep_tag_cohorts <- which(sim_env$conv_tag_release_indicator[,2] %in% 1:y)
    conv_tag_release_indicator <- sim_env$conv_tag_release_indicator[keep_tag_cohorts,,drop = FALSE]
    obs_conv_tag_fish_recap <- array(sim_env$obs_conv_tag_fish_recap[,,keep_tag_cohorts,,,,,,sim], dim = dim(sim_env$obs_conv_tag_fish_recap)[-length(dim(sim_env$obs_conv_tag_fish_recap))])
    conv_tagged_fish <- array(sim_env$conv_tagged_fish[keep_tag_cohorts,,,,sim], dim = c(dim(sim_env$conv_tagged_fish)[-length(dim(sim_env$conv_tagged_fish))]))
    conv_tagged_fish_attr <- array(sim_env$conv_tagged_fish_attr[keep_tag_cohorts,,,,sim], dim = c(dim(sim_env$conv_tagged_fish_attr)[-length(dim(sim_env$conv_tagged_fish_attr))]))
    n_tag_cohorts <- nrow(conv_tag_release_indicator)
  } else {
    conv_tag_release_indicator = obs_conv_tag_fish_recap = conv_tagged_fish = conv_tagged_fish_attr = n_tag_cohorts = NULL
  }

  # Fishery Landed Catches
  ObsCatch <- array(sim_env$ObsCatch[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaC <- array(sim_env$ln_sigmaC[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseCatch <- array(0, dim = dim(ObsCatch))
  UseCatch[!is.na(ObsCatch) & ObsCatch > 0] <- 1

  # Population-specific catches
  ObsCatch_pop <- array(sim_env$ObsCatch_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseCatch_pop <- array(0, dim = dim(ObsCatch_pop))
  UseCatch_pop[!is.na(ObsCatch_pop) & ObsCatch_pop > 0] <- 1
  ln_sigmaC_pop <- array(sim_env$ln_sigmaC_pop[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Discards
  ObsDiscard <- array(sim_env$ObsDiscard[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaD <- array(sim_env$ln_sigmaD[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseDiscard <- array(0, dim = dim(ObsDiscard))
  UseDiscard[!is.na(ObsDiscard) & ObsDiscard > 0] <- 1

  # Population-specific discards
  ObsDiscard_pop <- array(sim_env$ObsDiscard_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaD_pop <- array(sim_env$ln_sigmaD_pop[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseDiscard_pop <- array(0, dim = dim(ObsDiscard_pop))
  UseDiscard_pop[!is.na(ObsDiscard_pop) & ObsDiscard_pop > 0] <- 1

  # Fishery Indices
  ObsFishIdx <- array(sim_env$ObsFishIdx[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ObsFishIdx_SE <- array(sim_env$ObsFishIdx_SE[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseFishIdx <- array(0, dim = dim(ObsFishIdx))
  UseFishIdx[!is.na(ObsFishIdx) & ObsFishIdx > 0] <- 1

  # Population-specific fishery indices
  ObsFishIdx_pop <- array(sim_env$ObsFishIdx_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseFishIdx_pop <- array(0, dim = dim(ObsFishIdx_pop))
  UseFishIdx_pop[!is.na(ObsFishIdx_pop) & ObsFishIdx_pop > 0] <- 1
  ObsFishIdx_pop_SE <- array(sim_env$ObsFishIdx_pop_SE[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Retained Fishery Compositions
  ObsFishAgeComps <- array(sim_env$ObsFishAgeComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps <- array(sim_env$ISS_FishAgeComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps <- apply(ObsFishAgeComps, c(1,2,3,6), sum)
  UseFishAgeComps[!is.na(UseFishAgeComps) & UseFishAgeComps > 0] <- 1

  ObsFishLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps <- apply(ObsFishLenComps, c(1,2,3,6), sum)
    UseFishLenComps[!is.na(UseFishLenComps) & UseFishLenComps > 0] <- 1
  } else UseFishLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Population-specific retained fishery compositions
  ObsFishAgeComps_pop <- array(sim_env$ObsFishAgeComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_pop <- array(sim_env$ISS_FishAgeComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_pop <- apply(ObsFishAgeComps_pop, c(1,2,3,4,7), sum)
  UseFishAgeComps_pop[!is.na(UseFishAgeComps_pop) & UseFishAgeComps_pop > 0] <- 1

  ObsFishLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_pop <- apply(ObsFishLenComps_pop, c(1,2,3,4,7), sum)
    UseFishLenComps_pop[!is.na(UseFishLenComps_pop) & UseFishLenComps_pop > 0] <- 1
  } else UseFishLenComps_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Discarded Fishery Compositions
  ObsFishAgeComps_discard <- array(sim_env$ObsFishAgeComps_discard[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_discard <- array(sim_env$ISS_FishAgeComps_discard[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_discard <- apply(ObsFishAgeComps_discard, c(1,2,3,6), sum)
  UseFishAgeComps_discard[!is.na(UseFishAgeComps_discard) & UseFishAgeComps_discard > 0] <- 1

  ObsFishLenComps_discard <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_discard[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_discard <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_discard[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_discard <- apply(ObsFishLenComps_discard, c(1,2,3,6), sum)
    UseFishLenComps_discard[!is.na(UseFishLenComps_discard) & UseFishLenComps_discard > 0] <- 1
  } else UseFishLenComps_discard <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Population-specific discarded fishery compositions
  ObsFishAgeComps_discard_pop <- array(sim_env$ObsFishAgeComps_discard_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_discard_pop <- array(sim_env$ISS_FishAgeComps_discard_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_discard_pop <- apply(ObsFishAgeComps_discard_pop, c(1,2,3,4,7), sum)
  UseFishAgeComps_discard_pop[!is.na(UseFishAgeComps_discard_pop) & UseFishAgeComps_discard_pop > 0] <- 1

  ObsFishLenComps_discard_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_discard_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_discard_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_discard_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_discard_pop <- apply(ObsFishLenComps_discard_pop, c(1,2,3,4,7), sum)
    UseFishLenComps_discard_pop[!is.na(UseFishLenComps_discard_pop) & UseFishLenComps_discard_pop > 0] <- 1
  } else UseFishLenComps_discard_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Survey Indices
  ObsSrvIdx <- array(sim_env$ObsSrvIdx[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  ObsSrvIdx_SE <- array(sim_env$ObsSrvIdx_SE[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  UseSrvIdx <- array(0, dim = dim(ObsSrvIdx))
  UseSrvIdx[!is.na(ObsSrvIdx) & ObsSrvIdx > 0] <- 1

  # Population-specific survey indices
  ObsSrvIdx_pop <- array(sim_env$ObsSrvIdx_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  UseSrvIdx_pop <- array(0, dim = dim(ObsSrvIdx_pop))
  UseSrvIdx_pop[!is.na(ObsSrvIdx_pop) & ObsSrvIdx_pop > 0] <- 1
  ObsSrvIdx_pop_SE <- array(sim_env$ObsSrvIdx_pop_SE[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Survey Compositions
  ObsSrvAgeComps <- array(sim_env$ObsSrvAgeComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
  ISS_SrvAgeComps <- array(sim_env$ISS_SrvAgeComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  UseSrvAgeComps <- apply(ObsSrvAgeComps, c(1,2,3,6), sum)
  UseSrvAgeComps[!is.na(UseSrvAgeComps) & UseSrvAgeComps > 0] <- 1

  ObsSrvLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsSrvLenComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  ISS_SrvLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_SrvLenComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseSrvLenComps <- apply(ObsSrvLenComps, c(1,2,3,6), sum)
    UseSrvLenComps[!is.na(UseSrvLenComps) & UseSrvLenComps > 0] <- 1
  } else UseSrvLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Population-specific survey compositions
  ObsSrvAgeComps_pop <- array(sim_env$ObsSrvAgeComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
  ISS_SrvAgeComps_pop <- array(sim_env$ISS_SrvAgeComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  UseSrvAgeComps_pop <- apply(ObsSrvAgeComps_pop, c(1,2,3,4,7), sum)
  UseSrvAgeComps_pop[!is.na(UseSrvAgeComps_pop) & UseSrvAgeComps_pop > 0] <- 1

  ObsSrvLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsSrvLenComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  ISS_SrvLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_SrvLenComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseSrvLenComps_pop <- apply(ObsSrvLenComps_pop, c(1,2,3,4,7), sum)
    UseSrvLenComps_pop[!is.na(UseSrvLenComps_pop) & UseSrvLenComps_pop > 0] <- 1
  } else UseSrvLenComps_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Conditional age-at-length, present only when the simulation drew it. The use
  # flag marks length bins that received at least one aged fish.
  n_obs_ages <- dim(sim_env$AgeingError)[3]
  if(isTRUE(sim_env$do_fish_caal)) {
    ObsFish_caal <- array(sim_env$ObsFish_caal[,1:y,,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_caal_lens, n_obs_ages, sim_env$n_sexes, sim_env$n_fish_fleets))
    ISS_Fish_caal <- array(sim_env$ISS_Fish_caal[,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_caal_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
    UseFish_caal <- apply(ObsFish_caal, c(1,2,3,4,7), sum)
    UseFish_caal[] <- as.numeric(!is.na(UseFish_caal) & UseFish_caal > 0)
  } else ObsFish_caal <- ISS_Fish_caal <- UseFish_caal <- NULL
  if(isTRUE(sim_env$do_srv_caal)) {
    ObsSrv_caal <- array(sim_env$ObsSrv_caal[,1:y,,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_caal_lens, n_obs_ages, sim_env$n_sexes, sim_env$n_srv_fleets))
    ISS_Srv_caal <- array(sim_env$ISS_Srv_caal[,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_caal_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
    UseSrv_caal <- apply(ObsSrv_caal, c(1,2,3,4,7), sum)
    UseSrv_caal[] <- as.numeric(!is.na(UseSrv_caal) & UseSrv_caal > 0)
  } else ObsSrv_caal <- ISS_Srv_caal <- UseSrv_caal <- NULL

  # At-age data keep the operating model's use flags, since a normal likelihood can draw a value at or
  # below zero. A data source the operating model never draws comes back NULL.
  fish_aa_dim <- c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_ages, sim_env$n_sexes, sim_env$n_fish_fleets)
  srv_aa_dim <- c(sim_env$n_regions, length(1:y), sim_env$n_seas, n_obs_ages, sim_env$n_sexes, sim_env$n_srv_fleets)
  fish_aa_pop_dim <- c(sim_env$n_pop, fish_aa_dim)
  srv_aa_pop_dim <- c(sim_env$n_pop, srv_aa_dim)

  # Catch at age
  UseCatchAA <- ObsCatchAA <- ObsCatchAA_SE <- NULL
  if(isTRUE(any(sim_env$UseCatchAA == 1))) {
    UseCatchAA <- array(sim_env$UseCatchAA[,1:y,,,,, drop = FALSE], dim = fish_aa_dim)
    ObsCatchAA <- array(sim_env$ObsCatchAA[,1:y,,,,,sim, drop = FALSE], dim = fish_aa_dim)
    ObsCatchAA_SE <- array(sim_env$ObsCatchAA_SE[,1:y,,,,, drop = FALSE], dim = fish_aa_dim)
  }

  # Discards at age
  UseDiscardAA <- ObsDiscardAA <- ObsDiscardAA_SE <- NULL
  if(isTRUE(any(sim_env$UseDiscardAA == 1))) {
    UseDiscardAA <- array(sim_env$UseDiscardAA[,1:y,,,,, drop = FALSE], dim = fish_aa_dim)
    ObsDiscardAA <- array(sim_env$ObsDiscardAA[,1:y,,,,,sim, drop = FALSE], dim = fish_aa_dim)
    ObsDiscardAA_SE <- array(sim_env$ObsDiscardAA_SE[,1:y,,,,, drop = FALSE], dim = fish_aa_dim)
  }

  # Survey index at age
  UseSrvIdxAA <- ObsSrvIdxAA <- ObsSrvIdxAA_SE <- NULL
  if(isTRUE(any(sim_env$UseSrvIdxAA == 1))) {
    UseSrvIdxAA <- array(sim_env$UseSrvIdxAA[,1:y,,,,, drop = FALSE], dim = srv_aa_dim)
    ObsSrvIdxAA <- array(sim_env$ObsSrvIdxAA[,1:y,,,,,sim, drop = FALSE], dim = srv_aa_dim)
    ObsSrvIdxAA_SE <- array(sim_env$ObsSrvIdxAA_SE[,1:y,,,,, drop = FALSE], dim = srv_aa_dim)
  }

  # Population-specific catch at age
  UseCatchAA_pop <- ObsCatchAA_pop <- ObsCatchAA_pop_SE <- NULL
  if(isTRUE(any(sim_env$UseCatchAA_pop == 1))) {
    UseCatchAA_pop <- array(sim_env$UseCatchAA_pop[,,1:y,,,,, drop = FALSE], dim = fish_aa_pop_dim)
    ObsCatchAA_pop <- array(sim_env$ObsCatchAA_pop[,,1:y,,,,,sim, drop = FALSE], dim = fish_aa_pop_dim)
    ObsCatchAA_pop_SE <- array(sim_env$ObsCatchAA_pop_SE[,,1:y,,,,, drop = FALSE], dim = fish_aa_pop_dim)
  }

  # Population-specific discards at age
  UseDiscardAA_pop <- ObsDiscardAA_pop <- ObsDiscardAA_pop_SE <- NULL
  if(isTRUE(any(sim_env$UseDiscardAA_pop == 1))) {
    UseDiscardAA_pop <- array(sim_env$UseDiscardAA_pop[,,1:y,,,,, drop = FALSE], dim = fish_aa_pop_dim)
    ObsDiscardAA_pop <- array(sim_env$ObsDiscardAA_pop[,,1:y,,,,,sim, drop = FALSE], dim = fish_aa_pop_dim)
    ObsDiscardAA_pop_SE <- array(sim_env$ObsDiscardAA_pop_SE[,,1:y,,,,, drop = FALSE], dim = fish_aa_pop_dim)
  }

  # Population-specific survey index at age
  UseSrvIdxAA_pop <- ObsSrvIdxAA_pop <- ObsSrvIdxAA_pop_SE <- NULL
  if(isTRUE(any(sim_env$UseSrvIdxAA_pop == 1))) {
    UseSrvIdxAA_pop <- array(sim_env$UseSrvIdxAA_pop[,,1:y,,,,, drop = FALSE], dim = srv_aa_pop_dim)
    ObsSrvIdxAA_pop <- array(sim_env$ObsSrvIdxAA_pop[,,1:y,,,,,sim, drop = FALSE], dim = srv_aa_pop_dim)
    ObsSrvIdxAA_pop_SE <- array(sim_env$ObsSrvIdxAA_pop_SE[,,1:y,,,,, drop = FALSE], dim = srv_aa_pop_dim)
  }

  # Return
  return(list(
    # Biologicals
    WAA = WAA,
    WAA_fish = WAA_fish,
    WAA_srv = WAA_srv,
    MatAA = MatAA,
    SizeAgeTrans = SizeAgeTrans,
    AgeingError = AgeingError,
    AgeingError_fish = AgeingError_fish,
    AgeingError_srv = AgeingError_srv,
    LenBinMap = sim_env$LenBinMap,
    CAAL_LenBinMap = sim_env$CAAL_LenBinMap,

    # Tagging
    use_conv_fish_tagging = sim_env$use_conv_fish_tagging,
    conv_tag_release_indicator = conv_tag_release_indicator,
    obs_conv_tag_fish_recap = obs_conv_tag_fish_recap,
    conv_tagged_fish = conv_tagged_fish,
    conv_tagged_fish_attr = conv_tagged_fish_attr,
    n_tag_cohorts = n_tag_cohorts,

    # Aggregated catches
    ObsCatch = ObsCatch,
    ln_sigmaC = ln_sigmaC,
    UseCatch = UseCatch,
    ObsDiscard = ObsDiscard,
    ln_sigmaD = ln_sigmaD,
    UseDiscard = UseDiscard,

    # Population-specific catches
    ObsCatch_pop = ObsCatch_pop,
    ln_sigmaC_pop = ln_sigmaC_pop,
    UseCatch_pop = UseCatch_pop,
    ObsDiscard_pop = ObsDiscard_pop,
    ln_sigmaD_pop = ln_sigmaD_pop,
    UseDiscard_pop = UseDiscard_pop,

    # Aggregated fishery indices
    ObsFishIdx = ObsFishIdx,
    ObsFishIdx_SE = ObsFishIdx_SE,
    UseFishIdx = UseFishIdx,

    # Population-specific fishery indices
    ObsFishIdx_pop = ObsFishIdx_pop,
    ObsFishIdx_pop_SE = ObsFishIdx_pop_SE,
    UseFishIdx_pop = UseFishIdx_pop,

    # Aggregated retained fishery compositions
    ObsFishAgeComps = ObsFishAgeComps,
    ISS_FishAgeComps = ISS_FishAgeComps,
    UseFishAgeComps = UseFishAgeComps,
    ObsFishLenComps = ObsFishLenComps,
    ISS_FishLenComps = ISS_FishLenComps,
    UseFishLenComps = UseFishLenComps,

    # Population-specific retained fishery compositions
    ObsFishAgeComps_pop = ObsFishAgeComps_pop,
    ISS_FishAgeComps_pop = ISS_FishAgeComps_pop,
    UseFishAgeComps_pop = UseFishAgeComps_pop,
    ObsFishLenComps_pop = ObsFishLenComps_pop,
    ISS_FishLenComps_pop = ISS_FishLenComps_pop,
    UseFishLenComps_pop = UseFishLenComps_pop,

    # Aggregated discarded fishery compositions
    ObsFishAgeComps_discard = ObsFishAgeComps_discard,
    ISS_FishAgeComps_discard = ISS_FishAgeComps_discard,
    UseFishAgeComps_discard = UseFishAgeComps_discard,
    ObsFishLenComps_discard = ObsFishLenComps_discard,
    ISS_FishLenComps_discard = ISS_FishLenComps_discard,
    UseFishLenComps_discard = UseFishLenComps_discard,

    # Population-specific discarded fishery compositions
    ObsFishAgeComps_discard_pop = ObsFishAgeComps_discard_pop,
    ISS_FishAgeComps_discard_pop = ISS_FishAgeComps_discard_pop,
    UseFishAgeComps_discard_pop = UseFishAgeComps_discard_pop,
    ObsFishLenComps_discard_pop = ObsFishLenComps_discard_pop,
    ISS_FishLenComps_discard_pop = ISS_FishLenComps_discard_pop,
    UseFishLenComps_discard_pop = UseFishLenComps_discard_pop,

    # Aggregated survey indices
    ObsSrvIdx = ObsSrvIdx,
    ObsSrvIdx_SE = ObsSrvIdx_SE,
    UseSrvIdx = UseSrvIdx,

    # Population-specific survey indices
    ObsSrvIdx_pop = ObsSrvIdx_pop,
    ObsSrvIdx_pop_SE = ObsSrvIdx_pop_SE,
    UseSrvIdx_pop = UseSrvIdx_pop,

    # Aggregated survey compositions
    ObsSrvAgeComps = ObsSrvAgeComps,
    ISS_SrvAgeComps = ISS_SrvAgeComps,
    UseSrvAgeComps = UseSrvAgeComps,
    ObsSrvLenComps = ObsSrvLenComps,
    ISS_SrvLenComps = ISS_SrvLenComps,
    UseSrvLenComps = UseSrvLenComps,

    # Population-specific survey compositions
    ObsSrvAgeComps_pop = ObsSrvAgeComps_pop,
    ISS_SrvAgeComps_pop = ISS_SrvAgeComps_pop,
    UseSrvAgeComps_pop = UseSrvAgeComps_pop,
    ObsSrvLenComps_pop = ObsSrvLenComps_pop,
    ISS_SrvLenComps_pop = ISS_SrvLenComps_pop,
    UseSrvLenComps_pop = UseSrvLenComps_pop,

    # Conditional age-at-length
    ObsFish_caal = ObsFish_caal,
    ISS_Fish_caal = ISS_Fish_caal,
    UseFish_caal = UseFish_caal,
    ObsSrv_caal = ObsSrv_caal,
    ISS_Srv_caal = ISS_Srv_caal,
    UseSrv_caal = UseSrv_caal,

    # At-age data
    ObsCatchAA = ObsCatchAA,
    ObsCatchAA_SE = ObsCatchAA_SE,
    UseCatchAA = UseCatchAA,
    ObsDiscardAA = ObsDiscardAA,
    ObsDiscardAA_SE = ObsDiscardAA_SE,
    UseDiscardAA = UseDiscardAA,
    ObsSrvIdxAA = ObsSrvIdxAA,
    ObsSrvIdxAA_SE = ObsSrvIdxAA_SE,
    UseSrvIdxAA = UseSrvIdxAA,
    ObsCatchAA_pop = ObsCatchAA_pop,
    ObsCatchAA_pop_SE = ObsCatchAA_pop_SE,
    UseCatchAA_pop = UseCatchAA_pop,
    ObsDiscardAA_pop = ObsDiscardAA_pop,
    ObsDiscardAA_pop_SE = ObsDiscardAA_pop_SE,
    UseDiscardAA_pop = UseDiscardAA_pop,
    ObsSrvIdxAA_pop = ObsSrvIdxAA_pop,
    ObsSrvIdxAA_pop_SE = ObsSrvIdxAA_pop_SE,
    UseSrvIdxAA_pop = UseSrvIdxAA_pop
  ))

}
