# The operating model draws length compositions on the bins the data are recorded on, through the fit's
# LenBinMap, and keeps the fit's growth deviations over the conditioned years. The EBS Pacific cod
# self test on perfect data recovers the operating model to within the 0.001 observation error.

library(SPoRC)
library(testthat)

# six model length bins summed in pairs onto three recorded bins
len_bin_map <- kronecker(diag(3), matrix(1, 2, 1))

test_that("length compositions are drawn on the recorded bins through the map the fit uses", {

  set.seed(3)
  n_sexes <- 2
  exp_len <- array(c(5, 9, 14, 11, 6, 2, 3, 8, 12, 12, 7, 4), dim = c(1, 1, 1, 1, 6, n_sexes, 1, 1)) # [pop, region, year, season, len, sex, fleet, sim]
  iss <- array(1e8, dim = c(1, 1, 1, n_sexes, 1, 1))
  draw <- function(comp_type) {
    simulate_comps(r = 1, y = 1, f = 1, seas = 1, sim = 1, Exp = exp_len, ISS = iss, AgeingError = len_bin_map,
                   comp_like = 0, comp_type = matrix(comp_type, 1, 1), n_sexes = n_sexes, n_regions = 1, n_cat = 6,
                   Obs = array(0, dim = c(1, 1, 1, 3, n_sexes, 1, 1)), age_or_len = 1)
  }
  by_sex <- matrix(exp_len, 6, n_sexes) # [len, sex]

  # split by sex: each sex's composition normalized, then mapped, as Get_Comp_Likelihoods does
  split_obs <- draw(1)
  for(s in 1:n_sexes) {
    expected <- as.vector(by_sex[,s] %*% len_bin_map)
    expect_equal(split_obs[1,1,1,,s,1,1] / sum(split_obs[1,1,1,,s,1,1]), expected / sum(expected), tolerance = 1e-3)
  }

  # joint across sexes: one composition over the length by sex stack, each sex mapped on its own
  joint_obs <- draw(2)
  expected <- as.vector(as.vector(by_sex) %*% kronecker(diag(n_sexes), len_bin_map))
  expect_equal(as.vector(joint_obs[1,1,1,,,1,1]) / sum(joint_obs), expected / sum(expected), tolerance = 1e-3)

  # aggregated over sexes and regions
  agg_obs <- draw(0)
  expected <- as.vector(rowSums(by_sex) %*% len_bin_map)
  expect_equal(agg_obs[1,1,1,,1,1,1] / sum(agg_obs[1,1,1,,1,1,1]), expected / sum(expected), tolerance = 1e-3)

})

test_that("the operating model sizes length comps on the recorded bins and refuses a map that does not match", {

  n_yrs <- 2
  n_ages <- 3
  at_age <- array(0.5, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1)) # [pop, region, year, season, age, sex, sim]
  at_age_fleet <- array(0.5, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1, 1))
  ageing_error <- array(0, dim = c(n_yrs, n_ages, n_ages, 1)) # [year, model age, observed age, sim]
  for(y in 1:n_yrs) ageing_error[y,,,1] <- diag(n_ages)
  biologicals <- function(sim_list, map) {
    Setup_Sim_Biologicals(sim_list = sim_list, natmort_input = at_age, WAA_input = at_age, WAA_fish_input = at_age_fleet,
                          WAA_srv_input = at_age_fleet, MatAA_input = at_age * 0, AgeingError_input = ageing_error,
                          LenBinMap_input = map)
  }
  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = 6, n_obs_lens = 3,
                            n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1)
  sim_list <- Setup_Sim_Containers(sim_list)

  # observed comps on the three recorded bins, the true index at length on the six model bins
  expect_equal(dim(sim_list$ObsSrvLenComps)[4], 3)
  expect_equal(dim(sim_list$ObsFishLenComps_pop)[5], 3)
  expect_equal(dim(sim_list$SrvIAL)[5], 6)
  expect_equal(biologicals(sim_list, len_bin_map)$LenBinMap, len_bin_map)

  # no map for bins that differ, a map onto the wrong number of bins, and rows that sum to neither
  expect_error(biologicals(sim_list, NULL), "supply LenBinMap_input")
  expect_error(biologicals(sim_list, kronecker(diag(2), matrix(1, 3, 1))), "maps onto 2 length bins")
  expect_error(biologicals(sim_list, len_bin_map * 0.5), "sum to neither")

  # simulation lists that never set n_obs_lens record lengths on the model's bins
  old_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = 6, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1)
  old_list$n_obs_lens <- NULL
  expect_equal(dim(Setup_Sim_Containers(old_list)$ObsSrvLenComps)[4], 6)

})

test_that("growth deviations keep the fit's values over the conditioned years and are drawn given them after", {

  n_years <- 6
  n_cond <- 3
  map <- array(1:n_years, dim = c(1, 1, n_years, 1, 1))
  fit_devs <- array(c(0.2, -0.1, 0.5, 0, 0, 0), dim = c(1, 1, n_years, 1, 1))

  # separable AR1 over years: past the conditioned years the mean decays from the fit's last value
  rho_year <- 0.8
  cond_sd <- 0.3
  ar1_pars <- array(0, dim = c(1, 1, 4, 1))
  ar1_pars[1,1,2,1] <- 0.5 * log((1 + rho_year) / (1 - rho_year)) # year correlation on the scale rho_trans reads
  ar1_pars[1,1,4,1] <- log(cond_sd)
  set.seed(5)
  draws <- replicate(4000, draw_growth_pe_surface(5, map, ar1_pars, 1, fit_devs = fit_devs, n_cond = n_cond)[1,1,,1,1])
  expect_true(all(draws[1:n_cond,] == fit_devs[1,1,1:n_cond,1,1]))
  expect_equal(mean(draws[4,]), rho_year * 0.5, tolerance = 0.02)
  expect_equal(sd(draws[4,]), cond_sd, tolerance = 0.02)

  # a random walk steps on from the fit's last value
  walk_pars <- array(log(1e-6), dim = c(1, 1, 1, 1))
  walk <- draw_growth_pe_surface(2, map, walk_pars, 1, fit_devs = fit_devs, n_cond = n_cond)[1,1,,1,1]
  expect_equal(walk[4:6], rep(0.5, 3), tolerance = 1e-4)

  # every year conditioned gives back the fit under each form
  for(PE_model in c(1, 2, 5)) {
    pe_pars <- if(PE_model == 5) ar1_pars else walk_pars
    expect_identical(draw_growth_pe_surface(PE_model, map, pe_pars, 1, fit_devs = fit_devs, n_cond = n_years), fit_devs)
  }

})

test_that("the EBS Pacific cod self test on perfect data recovers the operating model", {

  skip_on_cran()
  data("sgl_rg_ebs_pcod_data")
  dat <- sgl_rg_ebs_pcod_data
  n_yrs <- length(dat$years)
  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = TRUE, newton_loops = 2, silent = TRUE)

  set.seed(17)
  sim_path <- tempfile(fileext = ".rds")
  st <- simulation_self_test(data = fit$data, parameters = fit$parameters, mapping = fit$mapping, random = NULL,
                             rep = fit$rep, sd_rep = list(par.fixed = fit$env$last.par.best), n_sims = 1,
                             newton_loops = 2, what = c("SSB", "Rec", "srv_q"), what_par = "ln_growth_pars",
                             perfect_data = TRUE, output_path = sim_path)
  sim <- readRDS(sim_path)

  # the operating model runs the fit's population, growth deviations included
  expect_equal(sim$SSB[1,1,,1], fit$rep$SSB[1,1,1:n_yrs], tolerance = 1e-10)

  # its length comps sit on the 24 recorded bins and match the fit's expected comps through the map
  norm <- function(v) v / sum(v)
  worst_len <- 0
  for(y in which(fit$data$UseSrvLenComps[1,,1,1] == 1)) {
    expected <- norm(as.vector(fit$rep$SrvIAL[1,1,y,1,,1,1] %*% dat$LenBinMap))
    worst_len <- max(worst_len, abs(norm(sim$ObsSrvLenComps[1,y,1,,1,1,1]) - expected))
  }
  for(y in which(fit$data$UseFishLenComps[1,,1,1] == 1)) {
    expected <- norm(as.vector(fit$rep$CAL[1,1,y,1,,1,1] %*% dat$LenBinMap))
    worst_len <- max(worst_len, abs(norm(sim$ObsFishLenComps[1,y,1,,1,1,1]) - expected))
  }
  expect_equal(dim(sim$ObsSrvLenComps)[4], ncol(dat$LenBinMap))
  expect_lt(worst_len, 0.005)

  # and the refit recovers it. before the bin map was applied SSB was off by a factor of about 1800
  rel_err <- function(est, truth) max(abs(est / truth - 1))
  expect_lt(rel_err(st$SSB, st$truth$SSB), 0.01)
  expect_lt(rel_err(st$srv_q, st$truth$srv_q), 0.005)
  # the last two year classes are seen in very few observations, so they get a looser bound
  rec_err <- abs(as.vector(st$Rec) / as.vector(st$truth$Rec) - 1)
  expect_lt(max(head(rec_err, -2)), 0.01)
  expect_lt(max(tail(rec_err, 2)), 0.05)
  expect_lt(max(abs(st$ln_growth_pars - st$truth$ln_growth_pars)), 0.005)

})
