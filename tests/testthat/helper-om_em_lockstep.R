# The operating model a self test builds, run at a model's current parameters and set beside the estimation
# model's report there. Conditioned on the model's own processes, every true state and observation is the prediction.

#' The self test's operating model at a model's current parameters
#'
#' @param il Input list, built and seeded however the caller likes; no fit is needed.
#' @param random Random effects the model integrates out.
#' @param perfect_data Passed to the self test, so the draws sit on the truth.
#' @param draw_rec Whether the recruitment and initial age deviations are drawn from the model's penalty, as a
#'   joint self test draws them, rather than kept at the model's own.
#' @param n_sims Replicates the operating model runs.
#'
#' @return List of the operating model's output (\code{om}), the estimation model's report at the same
#'   parameters (\code{rep}), its data list and the simulation list the operating model ran on.
#'
#' @keywords internal
lockstep_om <- function(il, random = NULL, perfect_data = TRUE, draw_rec = FALSE, n_sims = 1) {

  obj <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, random = random, do_optim = FALSE, silent = TRUE)))
  full <- obj$env$last.par.best # random effects at their inner mode, so the report and the conditioning agree
  rep <- obj$report(full)
  is_random <- seq_along(full) %in% obj$env$random
  sd_rep <- list(par.fixed = full[!is_random], par.random = full[is_random])

  # the self test builds the operating model; stop it before it simulates and refits
  captured <- new.env()
  attempt <- testthat::with_mocked_bindings(
    Simulate_Pop_Static = function(sim_list, ...) {
      captured$sim_list <- sim_list
      stop("captured")
    },
    try(suppressWarnings(suppressMessages(
      simulation_self_test(data = obj$data, parameters = il$par, mapping = il$map, random = random, rep = rep,
                           sd_rep = sd_rep, n_sims = n_sims, newton_loops = 0, what = "SSB", perfect_data = perfect_data))),
      silent = TRUE),
    .package = "SPoRC")
  if(is.null(captured$sim_list)) stop("the self test stopped before its operating model was built: ", attempt)

  # recruitment and the initial ages drawn from the penalty at the model's parameters, as joint draws them
  if(draw_rec) {
    pars <- obj$env$parList(par = full)
    n_yrs <- length(obj$data$years)
    captured$sim_list$Rec_input <- NULL
    captured$sim_list$ln_RecDevs_input <- simplify2array(lapply(seq_len(n_sims), function(i) rec_devs_past_fit(obj$data, pars, rep, 0, n_yrs, il$map$ln_RecDevs)))
    captured$sim_list$ln_InitDevs_input <- simplify2array(lapply(seq_len(n_sims), function(i) init_devs_past_fit(obj$data, pars, rep, il$map$ln_InitDevs)))
  }

  om <- suppressWarnings(suppressMessages(Simulate_Pop_Static(captured$sim_list)))
  list(om = om, rep = rep, data = obj$data, sim_list = captured$sim_list, il = il, random = random)
}

#' Largest relative difference between the operating model's truth and the estimation model's prediction
#'
#' @param truth,pred Numeric vectors, matched cell by cell.
#'
#' @return The largest |truth / pred - 1| over cells the prediction is not zero, or the largest |truth| where
#'   it is; \code{NA} with no cells.
#'
#' @keywords internal
lockstep_rel <- function(truth, pred) {
  truth <- as.vector(truth)
  pred <- as.vector(pred)
  if(!length(pred)) return(NA)
  nonzero <- pred != 0
  max(c(abs(truth[nonzero] / pred[nonzero] - 1), abs(truth[!nonzero])))
}

#' Every true state and fitted observation of the operating model against the estimation model's prediction
#'
#' The population state over the fitted years, then each data source over the cells its use flags fit: the
#' totals through the estimation model's own season and population sums (\code{get_seas_pred},
#' \code{get_seas_pred_pop}, \code{get_discard_pred}), the at-age sources cell by cell, and the tag
#' recaptures the operating model predicts.
#'
#' @param run Output of \code{lockstep_om}.
#'
#' @return Named numeric vector of largest relative differences, \code{NA} for anything the model does not have.
#'
#' @keywords internal
lockstep_diffs <- function(run) {

  om <- run$om
  rep <- run$rep
  data <- run$data
  n_yrs <- length(data$years)
  yrs <- seq_len(n_yrs)
  out <- c()

  # population state and the numbers behind every composition
  out["NAA"] <- lockstep_rel(om$NAA[,,yrs,,,,1], rep$NAA[,,yrs,,,])
  out["SSB"] <- lockstep_rel(om$SSB[,,yrs,1], rep$SSB[,,yrs])
  out["Rec"] <- lockstep_rel(om$Rec[,,yrs,1], rep$Rec[,,yrs])
  out["Fmort"] <- lockstep_rel(om$Fmort[,yrs,,,1], rep$Fmort[,yrs,,])
  out["tot_FAA"] <- lockstep_rel(om$tot_FAA[,,yrs,,,,,1], rep$tot_FAA[,,yrs,,,,])
  out["ZAA"] <- lockstep_rel(om$ZAA[,,yrs,,,,1], rep$ZAA[,,yrs,,,])
  out["CAA"] <- lockstep_rel(om$CAA[,,yrs,,,,,1], rep$CAA[,,yrs,,,,])
  out["DAA"] <- lockstep_rel(om$DAA[,,yrs,,,,,1], rep$DAA[,,yrs,,,,])
  out["SrvIAA"] <- lockstep_rel(om$SrvIAA[,,yrs,,,,,1], rep$SrvIAA[,,yrs,,,,])
  if(!is.null(om$CAL) && !is.null(rep$CAL) && isTRUE(data$fit_lengths == 1)) {
    out["CAL"] <- lockstep_rel(om$CAL[,,yrs,,,,,1], rep$CAL[,,yrs,,,,])
    out["SrvIAL"] <- lockstep_rel(om$SrvIAL[,,yrs,,,,,1], rep$SrvIAL[,,yrs,,,,])
  }

  # totals, each flagged cell against the estimation model's prediction for it
  n_seas <- data$n_seas
  pop_seq <- seq_len(data$n_pop)
  total_sources <- list(Catch = "PredCatch", FishIdx = "PredFishIdx", SrvIdx = "PredSrvIdx")
  for(data_name in names(total_sources)) {
    for(pop_source in c(FALSE, TRUE)) {
      name <- if(pop_source) paste0(data_name, "_pop") else data_name
      use <- data[[paste0("Use", name)]]
      truth_arr <- om[[paste0("True", name)]]
      if(is.null(use) || !any(use == 1) || is.null(truth_arr)) next
      seas_agg <- data[[paste0(name, "_seas_Type")]]
      if(is.null(seas_agg)) seas_agg <- rep(0, dim(use)[length(dim(use))])
      cells <- which(use == 1, arr.ind = TRUE)
      truth <- pred <- numeric(nrow(cells))
      for(i in seq_len(nrow(cells))) {
        cell <- cells[i,]
        if(pop_source) {
          truth[i] <- truth_arr[cell[1], cell[2], cell[3], cell[4], cell[5], 1]
          pred[i] <- get_seas_pred_pop(rep[[total_sources[[data_name]]]], cell[1], cell[2], cell[3], cell[4], cell[5], seas_agg[cell[5]])
        } else {
          truth[i] <- truth_arr[cell[1], cell[2], cell[3], cell[4], 1]
          pred[i] <- get_seas_pred(rep[[total_sources[[data_name]]]], cell[1], cell[2], cell[3], cell[4], seas_agg[cell[4]])
        }
      } # end i loop
      out[name] <- lockstep_rel(truth, pred)
    } # end pop_source loop
  } # end data_name loop

  # discards, through the estimation model's fraction when the units are one
  for(pop_source in c(FALSE, TRUE)) {
    name <- if(pop_source) "Discard_pop" else "Discard"
    use <- data[[paste0("Use", name)]]
    if(is.null(use) || !any(use == 1)) next
    seas_agg <- data[[paste0(name, "_seas_Type")]]
    if(is.null(seas_agg)) seas_agg <- rep(0, data$n_fish_fleets)
    cells <- which(use == 1, arr.ind = TRUE)
    truth <- pred <- numeric(nrow(cells))
    for(i in seq_len(nrow(cells))) {
      cell <- cells[i,]
      f <- cell[length(cell)]
      seasons <- if(seas_agg[f] == 1) seq_len(n_seas) else cell[length(cell) - 1]
      if(pop_source) {
        truth[i] <- om$TrueDiscard_pop[cell[1], cell[2], cell[3], cell[4], f, 1]
        pred[i] <- get_discard_pred(rep$PredDiscard, rep$CAA, rep$DAA, rep$dmr, data$WAA_fish, data$discard_units[f], cell[1], cell[2], cell[3], seasons, f)
      } else {
        truth[i] <- om$TrueDiscard[cell[1], cell[2], cell[3], f, 1]
        pred[i] <- get_discard_pred(rep$PredDiscard, rep$CAA, rep$DAA, rep$dmr, data$WAA_fish, data$discard_units[f], pop_seq, cell[1], cell[2], seasons, f)
      }
    } # end i loop
    out[name] <- lockstep_rel(truth, pred)
  } # end pop_source loop

  # at-age sources, cell by cell over the flags
  for(name in c("CatchAA", "DiscardAA", "SrvIdxAA", "CatchAA_pop", "DiscardAA_pop", "SrvIdxAA_pop")) {
    use <- data[[paste0("Use", name)]]
    if(is.null(use) || !any(use == 1)) next
    truth_arr <- om[[paste0("True", name)]]
    n_dims <- length(dim(truth_arr))
    truth_obs <- array(truth_arr, dim = dim(truth_arr))[slice.index(truth_arr, n_dims) == 1] # replicate one
    out[name] <- lockstep_rel(truth_obs[use == 1], rep[[paste0("Pred", name)]][use == 1])
  } # end name loop

  # tag recaptures, wherever the estimation model predicts any
  if(any(data$use_conv_fish_tagging == 1) && !is.null(rep$pred_conv_tag_fish_recap)) {
    om_recap <- om$pred_conv_tag_fish_recap
    n_dims <- length(dim(om_recap))
    om_recap <- om_recap[slice.index(om_recap, n_dims) == 1]
    out["TagRecap"] <- lockstep_rel(om_recap, rep$pred_conv_tag_fish_recap)
  }

  out

}

#' Each multinomial composition's likelihood at the truth, given the operating model's perfect-data draws
#'
#' Writes replicate one's compositions and sample sizes into the estimation model's data, every composition
#' weight at one, and evaluates the model at its own parameters. A composition drawn from the expectation the
#' model builds has an nLL of about (bins - 1) / 2 at a sample size of 1e6, one drawn from any other thousands.
#'
#' @param run Output of \code{lockstep_om} with \code{perfect_data = TRUE}.
#'
#' @return Data frame of each multinomial source fit: its likelihood per composition and (bins - 1) / 2.
#'
#' @keywords internal
lockstep_comp_nll <- function(run) {

  data <- run$data
  sources <- c("FishAgeComps", "FishLenComps", "SrvAgeComps", "SrvLenComps", "FishAgeComps_discard", "FishLenComps_discard",
               "FishAgeComps_pop", "FishLenComps_pop", "SrvAgeComps_pop", "SrvLenComps_pop",
               "FishAgeComps_discard_pop", "FishLenComps_discard_pop")
  first_rep <- function(x) { d <- dim(x); array(x[slice.index(x, length(d)) == 1], dim = d[-length(d)]) }
  fit_sources <- c()
  for(name in sources) {
    use <- data[[paste0("Use", name)]]
    if(is.null(use) || !any(use == 1) || is.null(run$om[[paste0("Obs", name)]])) next
    data[[paste0("Obs", name)]][] <- first_rep(run$om[[paste0("Obs", name)]])
    data[[paste0("ISS_", name)]][] <- first_rep(run$om[[paste0("ISS_", name)]])
    data[[paste0("Wt_", name)]][] <- 1
    fit_sources <- c(fit_sources, name)
  } # end name loop

  rep <- suppressWarnings(suppressMessages(fit_model(data, run$il$par, run$il$map, random = run$random, do_optim = FALSE, silent = TRUE)))$rep
  out <- NULL
  for(name in fit_sources) {
    if(!all(data[[paste0(name, "_LikeType")]] == 0)) next # the multinomial alone has this reference value
    obs <- data[[paste0("Obs", name)]]
    n_bins <- dim(obs)[if(grepl("_pop$", name)) 5 else 4]
    n_comps <- sum(data[[paste0("Use", name)]] == 1)
    out <- rbind(out, data.frame(source = name, nll_per_comp = sum(rep[[paste0(name, "_nLL")]]) / n_comps, reference = (n_bins - 1) / 2))
  } # end name loop
  out

}
