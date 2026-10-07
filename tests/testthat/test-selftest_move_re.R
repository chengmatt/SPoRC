# Movement random effects end to end: the operating model draws AR1 deviations around the fit's mean
# movement and rebuilds it, the estimation model at those deviations agrees, and a refit recovers the sd.

library(SPoRC)
library(testthat)

# a two region, one fleet, one sex sweep model, fitted once without deviations as the mean movement
move_re_setup <- local({
  cached <- NULL
  function() {
    if(!is.null(cached)) return(cached)
    small <- list(n_regions = 2, n_fish_fleets = 1, n_srv_fleets = 1, n_sexes = 1)
    base <- suppressMessages(sweep_input(dims = small, move = list(use_fixed_movement = 0, Fixed_Movement = NA)))
    base_fit <- suppressWarnings(fit_model(base$data, base$par, base$map, random = NULL, silent = TRUE, newton_loops = 1))
    em <- suppressMessages(sweep_input(dims = small, move = list(use_fixed_movement = 0, Fixed_Movement = NA, move_year_re = "ar1")))

    # the base fit's parameters, deviations at zero, and a process error of our choosing
    truth <- em$par
    base_pl <- base_fit$env$parList()
    for(nm in intersect(names(truth), names(base_pl))) if(length(truth[[nm]]) == length(base_pl[[nm]])) truth[[nm]] <- base_pl[[nm]]
    truth$move_devs[] <- 0
    truth$move_pe_pars[] <- 0
    truth$move_pe_pars[,,1] <- log(0.3) # conditional sd
    truth$move_pe_pars[,,3] <- 0.5 # year correlation of rho_trans(0.5)

    cached <<- list(em = em, truth = truth, base_fit = base_fit, n_yrs = length(em$data$years))
    cached
  }
})

test_that("the operating model draws the process the estimation model penalizes", {

  s <- move_re_setup()
  n_sims <- 400
  rho <- SPoRC:::rho_trans(0.5)
  marginal <- 0.3 / sqrt(1 - rho^2)

  sim_list <- list(n_pop = 1, n_regions = 2, n_yrs = s$n_yrs, n_seas = 1, n_ages = 7, n_sexes = 1, n_sims = n_sims,
                   Movement = replicate(n_sims, s$base_fit$rep$Movement), move_timing = 0, expm_nsub = 0)
  sim_list <- Setup_Sim_Movement(sim_list, s$em$data, s$truth)
  expect_equal(sim_list$move_year_re, 2)
  set.seed(5)
  se <- Setup_sim_env(sim_list)
  d <- se$move_devs
  expect_equal(dim(d), c(dim(s$truth$move_devs), n_sims))

  # the marginal sd and the lag one correlation over years, from one pair's year series across replicates
  x <- d[1,1,1,,1,3,1,]
  expect_equal(sd(as.vector(x)), marginal, tolerance = 0.1)
  expect_equal(cor(as.vector(x[-s$n_yrs,]), as.vector(x[-1,])), rho, tolerance = 0.15)
  expect_true(all(d[1,,,,1,1,1,] == 0)) # recruits do not move
  expect_equal(d[1,1,1,,1,2,1,], d[1,1,1,,1,7,1,]) # ages share the year deviation
  expect_false(isTRUE(all.equal(d[1,1,1,,1,3,1,1], d[1,2,1,,1,3,1,1]))) # pairs draw their own

  # movement rebuilt from the draw is Get_Movement at the draw, and the estimation model reports the same
  dev1 <- array(d[,,,,,,,1], dim = dim(s$truth$move_devs))
  args <- sim_list$move_args
  args$move_devs <- dev1
  by_hand <- do.call(Get_Movement, args)$Movement
  expect_equal(as.numeric(se$Movement[,,,,,,,1]), as.numeric(by_hand), tolerance = 1e-12)
  pars <- s$truth
  pars$move_devs <- dev1
  obj <- fit_model(s$em$data, pars, s$em$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_equal(as.numeric(obj$rep$Movement), as.numeric(by_hand), tolerance = 1e-12)
  for(y in 1:s$n_yrs) expect_equal(unname(rowSums(se$Movement[1,,,y,1,3,1,1])), rep(1, 2), tolerance = 1e-12) # still rows of probabilities
  expect_gt(max(abs(se$Movement[1,1,,,1,3,1,1] - s$base_fit$rep$Movement[1,1,,,1,3,1])), 1e-3) # and the draw moved them

  # the penalty at the draw is the ar1 density by hand, summed over the two pairs
  mvn <- function(v, S) { L <- chol(S); z <- backsolve(L, v, transpose = TRUE); -0.5 * (length(v) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2)) }
  S <- marginal^2 * rho^abs(outer(1:s$n_yrs, 1:s$n_yrs, "-"))
  by_hand_ll <- mvn(dev1[1,1,1,,1,2,1], S) + mvn(dev1[1,2,1,,1,2,1], S)
  expect_equal(-sum(obj$rep$Movement_nLL), by_hand_ll, tolerance = 1e-8)

})

test_that("a refit on simulated data recovers the movement process error", {

  # the three region CTMC operating model of helper-spatial_ctmc_om.R, with every other parameter pinned at
  # the truth, so the diffusion rate, the deviations and their process error are what the refit estimates
  om <- build_om(move_timing = 0)
  em <- build_em(om, 0, move_extra = list(move_year_re = "ar1", ctmc_diffusion_bounds = "upwind"))
  pinned <- pin_at_truth(em, TRUE_LOG_THETA, keep_free = c("move_devs", "move_pe_pars"))
  truth <- pinned$par
  truth$move_pe_pars[] <- 0
  truth$move_pe_pars[,,1] <- log(0.15) # conditional sd of a region's preference deviations
  truth$move_pe_pars[,,3] <- 0.5 # year correlation of rho_trans(0.5)
  at_truth <- fit_model(pinned$data, truth, pinned$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_equal(nrow(pinned$data$move_pairs), N_REGIONS) # one surface per region under the CTMC

  set.seed(31)
  res <- suppressWarnings(simulation_self_test(data = pinned$data, parameters = truth, mapping = pinned$map, random = "move_devs",
                                               rep = at_truth$rep, sd_rep = exact_pars_sd_rep(at_truth, "move_devs"), obj = at_truth,
                                               n_sims = 3, newton_loops = 1, sim_type = "joint",
                                               what = "SSB", what_par = c("move_pe_pars", "log_move_diffusion_pars")))
  expect_equal(sum(is.na(res$SSB)), 0) # every replicate refit

  # the diffusion rate, the conditional sd and the year correlation come back
  theta <- exp(res$log_move_diffusion_pars)
  expect_lt(max(abs(theta / exp(TRUE_LOG_THETA) - 1)), 0.25)
  pe <- res$move_pe_pars # [region, 1, 3, sim]
  est_sd <- exp(pe[,1,1,])
  est_rho <- SPoRC:::rho_trans(pe[,1,3,])
  expect_true(all(is.finite(est_sd)))
  expect_lt(abs(log(median(est_sd)) - log(0.15)), log(2)) # within a factor of two over nine region by replicate estimates
  expect_lt(abs(median(est_rho) - SPoRC:::rho_trans(0.5)), 0.4)

})

test_that("a conditioned closed loop keeps the fit's deviations and continues them", {

  # the small sweep model with ar1 year deviations at values of our own, so the fit's years are recognizable
  s <- move_re_setup()
  truth <- s$truth
  set.seed(8)
  fitted_series <- matrix(rnorm(2 * s$n_yrs, 0, 0.3), 2, s$n_yrs) # one series per pair
  for(a in 2:7) truth$move_devs[1,,1,,1,a,1] <- fitted_series # shared across the ages that move
  obj <- fit_model(s$em$data, truth, s$em$map, random = NULL, do_optim = FALSE, silent = TRUE)
  rho <- SPoRC:::rho_trans(0.5)
  n_proj <- 6
  n_sims <- 1000
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = n_proj, n_sims = n_sims, data = obj$data, parameters = obj$parameters, mapping = obj$mapping,
                                                sd_rep = list(par.fixed = obj$par, par.random = NULL), rep = obj$rep, random = NULL)
  expect_equal(sim_list$n_cond_yrs, s$n_yrs)
  set.seed(9)
  sim_env <- Setup_sim_env(sim_list)
  d <- sim_env$move_devs # [pop, from, to, year, season, age, sex, sim]

  # the fit's years are the fit's deviations in every replicate, and its movement
  for(sim in c(1, n_sims)) expect_equal(as.numeric(d[,,,1:s$n_yrs,,,,sim]), as.numeric(truth$move_devs), tolerance = 1e-12)
  expect_equal(as.numeric(sim_env$Movement[,,,1:s$n_yrs,,,,1]), as.numeric(obj$rep$Movement), tolerance = 1e-10)

  # the projection years are drawn, differ across replicates, and continue the ar1 from the last fitted year
  proj <- d[1,1,1,s$n_yrs + 1:n_proj,1,2,1,]
  expect_gt(sd(proj[1,]), 0.1)
  last_fit <- truth$move_devs[1,1,1,s$n_yrs,1,2,1]
  expect_lt(abs(mean(proj[1,]) - rho * last_fit), 0.04) # the first drawn year is centered on rho times the last fitted one
  expect_lt(abs(sd(proj[1,]) - 0.3), 0.03) # with the conditional sd as its spread
  expect_lt(abs(cor(proj[2,], proj[1,]) - rho), 0.1)
  expect_equal(d[1,1,1,s$n_yrs + 2,1,7,1,], proj[2,]) # ages still share the year's deviation

})

test_that("replicates on their own parameter draws rebuild movement from their own", {

  # two draws of the deviations, one value per map level so sharing across ages holds
  s <- move_re_setup()
  move_level <- as.integer(s$em$map$move_devs)
  n_levels <- max(move_level, na.rm = TRUE)
  set.seed(9)
  pars_by_sim <- vector("list", 2)
  for(sim in 1:2) {
    pars_sim <- s$truth
    level_draw <- stats::rnorm(n_levels, 0, 0.3)
    pars_sim$move_devs[!is.na(move_level)] <- level_draw[move_level[!is.na(move_level)]]
    pars_by_sim[[sim]] <- pars_sim
  } # end sim loop

  # every year conditioned, so each replicate keeps its own deviations throughout
  sim_list <- list(n_pop = 1, n_regions = 2, n_yrs = s$n_yrs, n_seas = 1, n_ages = 7, n_sexes = 1, n_sims = 2, n_cond_yrs = s$n_yrs,
                   Movement = replicate(2, s$base_fit$rep$Movement), move_timing = 0, expm_nsub = 0)
  sim_list <- Setup_Sim_Movement(sim_list, s$em$data, s$truth, pars_by_sim = pars_by_sim)
  se <- Setup_sim_env(sim_list)

  for(sim in 1:2) {
    expect_equal(as.numeric(se$move_devs[,,,,,,,sim]), as.numeric(pars_by_sim[[sim]]$move_devs))
    at_draw <- fit_model(s$em$data, pars_by_sim[[sim]], s$em$map, random = NULL, do_optim = FALSE, silent = TRUE)
    expect_equal(as.numeric(se$Movement[,,,,,,,sim]), as.numeric(at_draw$rep$Movement), tolerance = 1e-12)
  } # end sim loop
  expect_gt(max(abs(se$Movement[,,,,,,,1] - se$Movement[,,,,,,,2])), 1e-3)

})
