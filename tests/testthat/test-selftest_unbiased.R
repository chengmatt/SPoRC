# A joint self test, every process drawn fresh, returns the population without bias: over replicates, the mean log
# error of spawning biomass, recruitment, R0 and catchability sits within sampling error of zero.

library(SPoRC)
library(testthat)

test_that("a joint self test returns unbiased spawning biomass, recruitment, R0 and catchability", {

  # one region, sex and fleet, recruitment deviations integrated out under an estimated sigmaR
  il <- suppressWarnings(suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1),
                                                      rec = list(sigmaR_spec = "est_all"))))
  il$par$ln_sigmaR[] <- log(0.6)
  obj <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, random = "ln_RecDevs", do_optim = FALSE, silent = TRUE)))
  full <- obj$env$last.par.best

  set.seed(3)
  res <- suppressWarnings(suppressMessages(simulation_self_test(
    data = obj$data,
    parameters = il$par,
    mapping = il$map,
    random = "ln_RecDevs",
    rep = obj$report(full),
    sd_rep = exact_pars_sd_rep(obj, "ln_RecDevs"),
    obj = obj,
    n_sims = 80,
    newton_loops = 1,
    sim_type = "joint",
    what = c("SSB", "Rec"),
    what_par = c("ln_global_R0", "ln_fish_q", "ln_srv_q", "ln_sigmaR")
  )))
  expect_equal(sum(is.na(res$SSB)), 0) # every replicate refit

  # each replicate's mean log error over years, then their mean against its standard error
  log_err <- list(
    SSB = apply(log(res$SSB[1,1,,] / res$truth$SSB[1,1,,]), 2, mean),
    Rec = apply(log(res$Rec[1,1,,] / res$truth$Rec[1,1,,]), 2, mean),
    ln_global_R0 = as.vector(res$ln_global_R0) - as.vector(res$truth$ln_global_R0),
    ln_fish_q = as.vector(res$ln_fish_q) - as.vector(res$truth$ln_fish_q),
    ln_srv_q = as.vector(res$ln_srv_q) - as.vector(res$truth$ln_srv_q)
  )
  for(name in names(log_err)) {
    err <- log_err[[name]]
    expect_lt(abs(mean(err)) / (stats::sd(err) / sqrt(length(err))), 3.5, label = name)
  } # end name loop

  # a variance from 13 deviations comes back low by maximum likelihood's own bias, half the expected log of a
  # chi-square on 12 degrees of freedom over 13, so the sd is checked between that and zero
  sd_err <- res$ln_sigmaR[2,1,1,] - array(res$truth$ln_sigmaR, dim = dim(res$ln_sigmaR))[2,1,1,]
  sd_se <- stats::sd(sd_err) / sqrt(length(sd_err))
  ml_bias <- 0.5 * (digamma(12 / 2) + log(2) - log(13))
  expect_gt(mean(sd_err), ml_bias - 3 * sd_se)
  expect_lt(mean(sd_err), 3 * sd_se)

})

test_that("a replicate that fails to refit is NA in the results, which stay arrays", {

  il <- suppressWarnings(suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1))))
  obj <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))

  # the second refit fails, the others run as usual
  real_fit_model <- SPoRC::fit_model
  calls <- new.env()
  calls$n <- 0
  res <- testthat::with_mocked_bindings(
    fit_model = function(...) {
      calls$n <- calls$n + 1
      if(calls$n == 2) stop("refit failed")
      real_fit_model(...)
    },
    suppressWarnings(suppressMessages(simulation_self_test(data = obj$data, parameters = il$par, mapping = il$map, random = NULL,
                                                           rep = obj$rep, sd_rep = list(par.fixed = obj$par, par.random = NULL),
                                                           n_sims = 3, newton_loops = 0, what = "SSB", what_par = "ln_global_R0"))),
    .package = "SPoRC"
  )

  expect_true(is.array(res$SSB))
  expect_equal(dim(res$SSB), dim(res$truth$SSB))
  expect_equal(which(apply(is.na(res$SSB), 4, all)), 2)
  expect_true(is.na(res$ln_global_R0[2]) && all(!is.na(res$ln_global_R0[-2])))
  expect_true(all(is.finite((res$SSB / res$truth$SSB)[,,,-2]))) # arithmetic against the truth runs

})

test_that("a prior the truth sits away from pulls every refit unless its mean is moved to the truth or it is off", {

  # an informative survey catchability prior centered half again above the value the replicates run on
  il <- suppressWarnings(suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1))))
  il$data$Use_srv_q_prior <- 1
  il$data$srv_q_prior <- data.frame(region = 1, block = 1, fleet = 1, mu = 1.5 * exp(il$par$ln_srv_q[1,1,1]), sd = 0.1)
  obj <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))
  full <- obj$env$last.par.best

  # each setting's mean error in log q across replicates, against its standard error
  q_t <- function(prior_means) {
    set.seed(5)
    res <- suppressWarnings(suppressMessages(simulation_self_test(
      data = obj$data, parameters = il$par, mapping = il$map, random = NULL, rep = obj$report(full),
      sd_rep = exact_pars_sd_rep(obj), obj = obj, n_sims = 30, newton_loops = 0, sim_type = "joint",
      what = "SSB", what_par = "ln_srv_q", prior_means = prior_means
    )))
    err <- matrix(res$ln_srv_q, ncol = 30)[1,] - matrix(res$truth$ln_srv_q, ncol = 30)[1,] # the one survey's q, by replicate
    mean(err) / (stats::sd(err) / sqrt(length(err)))
  }

  expect_gt(q_t("assessment"), 3.5) # pulled toward the prior's larger mean
  expect_lt(abs(q_t("truth")), 3.5)
  expect_lt(abs(q_t("off")), 3.5)

})

test_that("a joint self test compares an observation sd with the fit's value, which every replicate runs at", {

  il <- suppressWarnings(suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1),
                                                      catch = list(sigmaC_spec = "est_shared_r_y_seas_f"))))
  obj <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))
  full <- obj$env$last.par.best
  expect_equal(sum(names(full) == "ln_sigmaC"), 1)

  # the catch sd drawn at an sd of one, every other parameter kept at the fit
  prec <- Matrix::.sparseDiagonal(length(full), ifelse(names(full) == "ln_sigmaC", 1, 1e12), shape = "s")
  set.seed(4)
  res <- suppressWarnings(suppressMessages(simulation_self_test(
    data = obj$data, parameters = il$par, mapping = il$map, random = NULL, rep = obj$report(full),
    sd_rep = list(par.fixed = full, par.random = NULL, jointPrecision = prec), obj = obj, n_sims = 3,
    newton_loops = 0, sim_type = "joint", what = "SSB", what_par = "ln_sigmaC"
  )))
  expect_true(all(res$truth$ln_sigmaC == full[["ln_sigmaC"]])) # every year and replicate, one shared level

})

test_that("a refit keeps the fit's penalty weights and resets its data weights", {

  # recruitment penalized at half weight, F at twice, the catch at three times
  il <- suppressWarnings(suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1),
                                                      wt = list(Wt_Rec = 0.5, Wt_F = 2, Wt_Catch = 3))))
  obj <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))
  expect_true(all(obj$data$Wt_Rec == 0.5) && all(obj$data$Wt_Catch == 3))

  # the data each refit is handed
  seen <- new.env()
  testthat::with_mocked_bindings(
    fit_model = function(data, parameters, ...) { assign("data", data, envir = seen); assign("pars", parameters, envir = seen); stop("intercepted") },
    try(suppressWarnings(suppressMessages(simulation_self_test(data = obj$data, parameters = il$par, mapping = il$map, random = NULL,
                                                               rep = obj$rep, sd_rep = list(par.fixed = obj$par, par.random = NULL),
                                                               n_sims = 1, newton_loops = 0, what = "SSB"))), silent = TRUE),
    .package = "SPoRC"
  )
  refit_data <- get("data", envir = seen)
  refit_pars <- get("pars", envir = seen)
  expect_true(all(refit_data$Wt_Rec == 0.5)) # the penalty weights stay, so the refit is the fit's own estimator
  expect_equal(refit_data$Wt_F, 2)
  expect_true(all(refit_data$Wt_Catch == 1)) # the data weights go to one, the catch having been drawn at its sd over root three
  expect_equal(as.vector(exp(refit_pars$ln_sigmaC)), as.vector(exp(il$par$ln_sigmaC) / sqrt(3)), tolerance = 1e-12)

})

test_that("the movement Dirichlet prior is moved so its mode sits on the fractions the replicate ran on", {

  # two regions with estimated movement and a Dirichlet prior on each region's annual fractions, concentration 6
  prior <- expand.grid(pop = 1, region_from = 1:2, year = 1, age = 1, seas = 1, sex = 1, alpha = I(list(c(4, 2))))
  il <- suppressWarnings(suppressMessages(sweep_input(dims = list(n_regions = 2, n_sexes = 1, n_fish_fleets = 1),
                                                      move = list(use_fixed_movement = 0, Fixed_Movement = NA, Use_Movement_Prior = 1, Movement_prior = prior))))
  obj <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))
  pars <- obj$env$parList()
  rep <- obj$rep
  moved <- SPoRC:::self_test_priors(obj$data, "truth", pars, rep)$Movement_prior

  # the annual fractions the prior is evaluated on, as get_movement_dirichlet_prior reads them
  frac_true <- function(i) {
    if(is.null(rep$Mrate)) rep$Movement[1, prior$region_from[i], , 1, 1, 1, 1]
    else as.matrix(Matrix::expm(methods::as(rep$Mrate[1, , , 1, 1, 1, 1], "sparseMatrix")))[prior$region_from[i], ]
  }
  for(i in 1:2) {
    alpha <- moved$alpha[[i]]
    expect_equal(sum(alpha), 6) # the concentration is kept
    expect_equal((alpha - 1) / (sum(alpha) - 2), frac_true(i), tolerance = 1e-12) # the mode is the truth
  } # end i loop

  # so the prior has no gradient at the truth along the simplex, where the assessment's did
  nll_at <- function(prior_df, shift) {
    movement <- rep$Movement
    movement[1, 1, , 1, 1, 1, 1] <- movement[1, 1, , 1, 1, 1, 1] + c(shift, -shift)
    mrate <- rep$Mrate
    if(!is.null(mrate)) return(NA) # the CTMC branch is covered by the mode identity above
    SPoRC:::get_movement_dirichlet_prior(prior_df, movement, NULL)
  }
  if(is.null(rep$Mrate)) {
    h <- 1e-4
    expect_lt(abs((nll_at(moved, h) - nll_at(moved, -h)) / (2 * h)), 1e-6)
    expect_gt(abs((nll_at(obj$data$Movement_prior, h) - nll_at(obj$data$Movement_prior, -h)) / (2 * h)), 1e-2)
  }

  # a flat Dirichlet has no interior mode to move and is left as given
  flat <- obj$data
  flat$Movement_prior$alpha <- I(list(c(1, 1), c(1, 1)))
  expect_equal(SPoRC:::self_test_priors(flat, "truth", pars, rep)$Movement_prior$alpha[[1]], c(1, 1))

})

test_that("a Dirichlet or beta prior moved to a value has its mode there and keeps its concentration", {

  # Dirichlet: the mode (alpha - 1) / (sum(alpha) - K) is the fractions, the sum is kept, a flat one is left alone
  alpha <- SPoRC:::dirichlet_mode_at(c(5, 3, 2), c(0.6, 0.3, 0.1))
  expect_equal(sum(alpha), 10)
  expect_equal((alpha - 1) / (sum(alpha) - 3), c(0.6, 0.3, 0.1))
  expect_equal(SPoRC:::dirichlet_mode_at(c(1, 1, 1), c(0.6, 0.3, 0.1)), c(1, 1, 1))

  # beta: the shapes get_tagrep_prior forms from the moved mean and sd have their mode at x and the same concentration
  moved <- SPoRC:::beta_mode_at(0.3, 0.1, 0.7)
  conc <- function(mu, sd) mu * (1 - mu) / sd^2 - 1
  expect_equal(conc(moved[1], moved[2]), conc(0.3, 0.1))
  shape1 <- moved[1] * conc(moved[1], moved[2])
  shape2 <- (1 - moved[1]) * conc(moved[1], moved[2])
  expect_equal((shape1 - 1) / (shape1 + shape2 - 2), 0.7)
  expect_equal(stats::dbeta(0.7, shape1, shape2), max(stats::dbeta(seq(0.01, 0.99, 0.001), shape1, shape2)), tolerance = 1e-4)
  expect_equal(SPoRC:::beta_mode_at(0.5, 0.5, 0.7), c(0.5, 0.5)) # a concentration below two has no interior mode

})
