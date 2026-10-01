# Growth deviations end to end: the operating model draws each varying parameter's series and the
# semi-parametric surface from the process the estimation model penalizes, rebuilds growth from them,
# the estimation model at those deviations agrees, and a refit recovers the sd.

library(SPoRC)
library(testthat)

# the small one region model of helper-selftest_growth_semipar.R, its biologicals overridden per test
growth_re_obs <- local({
  cached <- NULL
  function() {
    if(is.null(cached)) cached <<- semipar_simulate(seed = 11)$obs
    cached
  }
})

growth_re_input <- function(biol) {
  with_overrides <- function(stage, overrides) function(...) { given <- list(...); given[names(overrides)] <- NULL; do.call(stage, c(given, overrides)) }
  build <- semipar_input
  body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Biologicals = with_overrides(Setup_Mod_Biologicals, biol))))
  suppressMessages(build("none", growth_re_obs()))
}

# the smallest operating model list the draw needs: the dims and the report's growth for every replicate
growth_re_sim_list <- function(rep, n_sims) {
  list(n_pop = 1, n_regions = 1, n_yrs = spcfg$n_yrs, n_seas = 1, n_ages = spcfg$n_ages, n_sexes = 1, n_sims = n_sims,
       WAA = replicate(n_sims, rep$WAA), WAA_fish = replicate(n_sims, rep$WAA_fish), WAA_srv = replicate(n_sims, rep$WAA_srv),
       SizeAgeTrans_fish = replicate(n_sims, rep$SizeAgeTrans_fish), SizeAgeTrans_srv = replicate(n_sims, rep$SizeAgeTrans_srv))
}

rho_untrans <- function(x) 0.5 * log((1 + x) / (1 - x))
mvn_logdens <- function(v, cov_by_hand) { L <- chol(cov_by_hand); z <- backsolve(L, v, transpose = TRUE); -0.5 * (length(v) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2)) }
ar1_corr <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

test_that("the operating model draws iid and random walk series on the growth parameters", {

  input_list <- growth_re_input(list(growth_tv_model = c(L1 = "iid", K = "rw"), growth_tv_sigma_spec = "fix", growth_semipar = "none"))
  truth <- input_list$par
  truth$growth_pe_pars[1,1,1,1,1] <- log(0.1) # L1, iid
  truth$growth_pe_pars[1,1,3,1,1] <- log(0.05) # K, random walk
  obj <- fit_model(input_list$data, truth, input_list$map, random = NULL, do_optim = FALSE, silent = TRUE)
  n_yrs <- spcfg$n_yrs
  n_sims <- 200

  sim_list <- Setup_Sim_Growth_RE(growth_re_sim_list(obj$rep, n_sims), input_list$data, truth)
  expect_equal(sim_list$growth_tv_model, c(1, 0, 2, 0, 0))
  set.seed(5)
  sim_env <- Setup_sim_env(sim_list)
  devs <- sim_env$ln_growth_devs
  expect_equal(dim(devs), c(dim(truth$ln_growth_devs), n_sims))

  # the iid sd, the walk's first year at its own sd and its steps at that sd, across replicates
  expect_equal(sd(as.vector(devs[1,1,,1,1,])), 0.1, tolerance = 0.1)
  expect_equal(sd(devs[1,1,1,3,1,]), 0.05, tolerance = 0.25)
  expect_equal(sd(as.vector(apply(devs[1,1,,3,1,], 2, diff))), 0.05, tolerance = 0.1)
  expect_true(all(devs[1,1,,c(2, 4, 5),1,] == 0)) # the parameters that do not vary
  expect_gt(max(abs(devs[1,1,,1,1,1] - devs[1,1,,1,1,2])), 0.05) # replicates draw their own

  # weight at age rebuilt from the draw is Get_Growth at the draw, and the estimation model reports the same
  devs_rep1 <- array(devs[,,,,,1], dim = dim(truth$ln_growth_devs))
  growth_args <- sim_list$growth_args
  growth_args$ln_growth_devs <- devs_rep1
  growth_args$ln_growth_semipar_devs <- array(sim_env$ln_growth_semipar_devs[,,,,,1], dim = dim(truth$ln_growth_semipar_devs))
  by_hand <- do.call(Get_Growth, growth_args)
  expect_equal(as.numeric(sim_env$WAA[,,,,,,1]), as.numeric(by_hand$WAA), tolerance = 1e-12)
  expect_equal(as.numeric(sim_env$SizeAgeTrans_srv[,,,,,,,,1]), as.numeric(by_hand$SizeAgeTrans_srv), tolerance = 1e-12)
  pars_at_draw <- truth
  pars_at_draw$ln_growth_devs <- devs_rep1
  at_draw <- fit_model(input_list$data, pars_at_draw, input_list$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_equal(as.numeric(at_draw$rep$WAA), as.numeric(by_hand$WAA), tolerance = 1e-12)
  expect_equal(as.numeric(at_draw$rep$SizeAgeTrans_fish), as.numeric(by_hand$SizeAgeTrans_fish), tolerance = 1e-12)
  expect_gt(max(abs(sim_env$WAA[1,1,,1,,1,1] - obj$rep$WAA[1,1,,1,,1])), 1e-3) # and the draw moved growth

  # the penalty at the draw is the density by hand; the estimation model gives the walk's first year
  # the diffuse sd growth_rw_init_sigma, the one convention the draw does not copy
  by_hand_ll <- sum(dnorm(devs_rep1[1,1,,1,1], 0, 0.1, TRUE)) + dnorm(devs_rep1[1,1,1,3,1], 0, input_list$data$growth_rw_init_sigma, TRUE) +
    sum(dnorm(diff(devs_rep1[1,1,,3,1]), 0, 0.05, TRUE))
  expect_equal(-at_draw$rep$growth_tv_nLL, by_hand_ll, tolerance = 1e-8)

})

test_that("the operating model draws the semi-parametric surface under the separable AR1", {

  input_list <- growth_re_input(list(growth_semipar = "2dar1", growth_semipar_spec = "fix"))
  truth <- input_list$par
  truth$growth_pe_pars[1,1,1,1,2] <- rho_untrans(0.5) # across ages
  truth$growth_pe_pars[1,1,2,1,2] <- rho_untrans(0.3) # across years
  truth$growth_pe_pars[1,1,4,1,2] <- log(0.08) # conditional sd
  obj <- fit_model(input_list$data, truth, input_list$map, random = NULL, do_optim = FALSE, silent = TRUE)
  n_yrs <- spcfg$n_yrs
  n_ages <- spcfg$n_ages
  n_sims <- 100
  marginal <- 0.08 / sqrt(1 - 0.3^2) / sqrt(1 - 0.5^2)

  sim_list <- Setup_Sim_Growth_RE(growth_re_sim_list(obj$rep, n_sims), input_list$data, truth)
  set.seed(6)
  sim_env <- Setup_sim_env(sim_list)
  devs <- sim_env$ln_growth_semipar_devs
  expect_equal(dim(devs), c(dim(truth$ln_growth_semipar_devs), n_sims))
  expect_true(all(sim_env$ln_growth_devs == 0)) # no parameter varies

  # the marginal sd and the lag one correlations over years and ages, across replicates
  surface <- devs[1,1,,,1,]
  expect_equal(sd(as.vector(surface)), marginal, tolerance = 0.1)
  expect_equal(cor(as.vector(surface[-n_yrs,,]), as.vector(surface[-1,,])), 0.3, tolerance = 0.15)
  expect_equal(cor(as.vector(surface[,-n_ages,]), as.vector(surface[,-1,])), 0.5, tolerance = 0.1)

  # growth rebuilt from the draw is what the estimation model reports, and the penalty is the Kronecker density by hand
  devs_rep1 <- array(devs[,,,,,1], dim = dim(truth$ln_growth_semipar_devs))
  pars_at_draw <- truth
  pars_at_draw$ln_growth_semipar_devs <- devs_rep1
  at_draw <- fit_model(input_list$data, pars_at_draw, input_list$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_equal(as.numeric(sim_env$WAA[,,,,,,1]), as.numeric(at_draw$rep$WAA), tolerance = 1e-12)
  expect_equal(as.numeric(sim_env$SizeAgeTrans_fish[,,,,,,,,1]), as.numeric(at_draw$rep$SizeAgeTrans_fish), tolerance = 1e-12)
  cov_by_hand <- marginal^2 * kronecker(ar1_corr(n_ages, 0.5), ar1_corr(n_yrs, 0.3)) # year fastest, as the surface is stored
  expect_equal(-at_draw$rep$growth_semipar_nLL, mvn_logdens(as.vector(devs_rep1[1,1,,,1]), cov_by_hand), tolerance = 1e-8)

})

test_that("the three dimensional field is drawn conditional on the cells the map fixes", {

  for(form in c("3dmarg", "3dcond")) {

    # ages and years left out of the surface stay at zero and sit inside the density at zero
    input_list <- growth_re_input(list(growth_semipar = form, growth_semipar_spec = "fix", growth_semipar_ages = 2:11, growth_semipar_years = 3:30))
    truth <- input_list$par
    truth$growth_pe_pars[1,1,1:3,1,2] <- c(0.3, 0.4, 0.2) # partial correlations over age, year and cohort
    truth$growth_pe_pars[1,1,4,1,2] <- log(0.1^2) # log variance
    obj <- fit_model(input_list$data, truth, input_list$map, random = NULL, do_optim = FALSE, silent = TRUE)
    n_yrs <- spcfg$n_yrs

    sim_list <- Setup_Sim_Growth_RE(growth_re_sim_list(obj$rep, 3), input_list$data, truth)
    set.seed(7)
    sim_env <- Setup_sim_env(sim_list)
    devs <- sim_env$ln_growth_semipar_devs
    expect_true(all(devs[1,1,1:2,,1,] == 0))
    expect_true(all(devs[1,1,,c(1, 12),1,] == 0))
    expect_true(all(devs[1,1,3:n_yrs,2:11,1,] != 0))

    devs_rep1 <- array(devs[,,,,,1], dim = dim(truth$ln_growth_semipar_devs))
    pars_at_draw <- truth
    pars_at_draw$ln_growth_semipar_devs <- devs_rep1
    at_draw <- fit_model(input_list$data, pars_at_draw, input_list$map, random = NULL, do_optim = FALSE, silent = TRUE)
    expect_equal(as.numeric(sim_env$WAA[,,,,,,1]), as.numeric(at_draw$rep$WAA), tolerance = 1e-12)
    precision <- as.matrix(SPoRC:::Get_3d_precision(10, n_yrs, 0.3, 0.4, 0.2, log(0.1^2), Var_Type = if(form == "3dmarg") 0 else 1))
    nodes <- as.vector(t(devs_rep1[1,1,,2:11,1])) # age fastest, as the precision numbers its nodes
    by_hand_ll <- 0.5 * as.numeric(determinant(precision, logarithm = TRUE)$modulus) - 0.5 * as.numeric(t(nodes) %*% precision %*% nodes) - 0.5 * length(nodes) * log(2 * pi)
    expect_equal(-at_draw$rep$growth_semipar_nLL, by_hand_ll, tolerance = 1e-8)

  } # end form loop

})

test_that("cohort growth advances from drawn deviations through the closed loop", {

  # pcod varies L1 and K over part of its years, propagates cohort by cohort and selects at length
  data("sgl_rg_ebs_pcod_data", envir = environment())
  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(sgl_rg_ebs_pcod_data))), sgl_rg_ebs_pcod_data)
  n_yrs <- length(input_list$data$years)
  cohort_styr <- input_list$data$growth_cohort_styr
  obj <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL, do_optim = FALSE, silent = TRUE)
  fit_pars <- obj$env$parList()
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = 3, n_sims = 2, data = obj$data, parameters = obj$parameters, mapping = obj$mapping,
                                                sd_rep = list(par.fixed = obj$par, par.random = NULL), rep = obj$rep, random = NULL)
  expect_false(is.null(sim_list$growth_args))
  expect_equal(sim_list$growth_length_sel$fish_selex_type, 1)
  set.seed(8)
  sim_env <- Setup_sim_env(sim_list) # built here, since the annual cycle finds the environment by this name in the frame that built it
  devs <- sim_env$ln_growth_devs
  devs_map <- input_list$data$map_ln_growth_devs
  expect_equal(dim(devs), c(1, 1, n_yrs + 3, 6, 1, 2))

  # the fit's devs_map decides the fitted years, the closed loop years are active, and the draws are each replicate's own
  expect_true(all(devs[,,1:n_yrs,,,][is.na(devs_map)] == 0))
  expect_true(all(devs[,,1:n_yrs,,,1][!is.na(devs_map)] != 0))
  expect_true(all(devs[1,1,n_yrs + 1:3,c(1, 3),1,] != 0))
  expect_true(all(devs[1,1,,c(2, 4:6),1,] == 0))
  expect_gt(max(abs(devs[,,,,,1] - devs[,,,,,2])), 0.01)
  expect_length(sim_env$growth_state, 2)

  for(y in 1:(n_yrs + 3)) run_annual_cycle(y, 1, sim_env)

  # the same chain by hand, each year's plus group blended by the numbers this replicate had at the start of it
  growth_args <- sim_list$growth_args
  growth_args$ln_growth_devs <- array(devs[,,,,,1], dim = dim(devs)[-6])
  growth_args$ln_growth_semipar_devs <- array(sim_env$ln_growth_semipar_devs[,,,,,1], dim = dim(sim_env$ln_growth_semipar_devs)[-6])
  growth <- do.call(Get_Growth, growth_args)
  year_args <- growth_args[names(growth_args) %in% names(formals(Get_Growth_Year))]
  for(y in cohort_styr:(n_yrs + 3)) {
    growth <- do.call(Get_Growth_Year, c(year_args, list(growth = growth, y = y, NAA_y = array(sim_env$NAA[,,y,1,,,1], dim = c(1, 1, length(input_list$data$ages), 1)))))
  } # end y loop
  expect_equal(as.numeric(sim_env$WAA[1,1,,,,,1]), as.numeric(growth$WAA[1,1,,,,]), tolerance = 1e-10)
  expect_equal(as.numeric(sim_env$SizeAgeTrans_fish[1,1,,,,,,,1]), as.numeric(growth$SizeAgeTrans_fish[1,1,,,,,,]), tolerance = 1e-10)
  expect_true(all(is.finite(sim_env$SSB[,,,1])))
  expect_gt(max(abs(sim_env$WAA[1,1,cohort_styr:n_yrs,1,,1,1] - obj$rep$WAA[1,1,cohort_styr:n_yrs,1,,1])), 1e-4) # the propagated years are not the report's

})

test_that("a refit on simulated data recovers the growth process error", {

  # iid deviations on K with their sd estimated; every other parameter but the recruitment deviations is
  # pinned at the truth, and both deviation sets are integrated
  input_list <- growth_re_input(list(growth_tv_model = c(K = "iid"), growth_tv_sigma_spec = "est", growth_semipar = "none"))
  truth <- input_list$par
  truth$ln_growth_devs[] <- 0
  truth$growth_pe_pars[1,1,3,1,1] <- log(0.15)
  free_pars <- c("ln_growth_devs", "growth_pe_pars", "ln_RecDevs")
  pinned_map <- input_list$map
  for(nm in setdiff(names(truth), free_pars)) pinned_map[[nm]] <- factor(rep(NA, length(truth[[nm]])))
  expect_equal(sum(!is.na(pinned_map$growth_pe_pars)), 1)
  at_truth <- fit_model(input_list$data, truth, pinned_map, random = NULL, do_optim = FALSE, silent = TRUE)

  set.seed(31)
  self_test <- suppressWarnings(simulation_self_test(data = input_list$data, parameters = truth, mapping = pinned_map, random = c("ln_growth_devs", "ln_RecDevs"),
                                               rep = at_truth$rep, sd_rep = NULL, n_sims = 3, newton_loops = 1,
                                               what = "SSB", what_par = "growth_pe_pars"))
  expect_equal(sum(is.na(self_test$SSB)), 0) # every replicate refit
  est_sd <- exp(self_test$growth_pe_pars[1,1,3,1,1,])
  expect_true(all(is.finite(est_sd)))
  expect_lt(max(abs(log(est_sd / 0.15))), log(2)) # within a factor of two in every replicate

})
