# Length compositions selected at length (FishLenComps_sel / SrvLenComps_sel = "length") in the operating
# model: conditioned on the EBS Pacific cod fit, catch and index at length match the fit exactly, age-at-length
# rows follow the length-selected joint, and retention at length splits catch and discards as the fit does.

library(SPoRC)
library(testthat)

test_that("the operating model refuses selectivity at length without the selectivity at length", {

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = 2, n_regions = 1, n_ages = 3, n_lens = 4, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  sel_at_age <- array(1, dim = c(1, 1, 2, 1, 3, 1, 1, 1)) # [pop, region, year, season, age, sex, fleet, sim]
  wrong_dims <- array(1, dim = c(1, 2, 3, 1, 1, 1)) # three length bins where the model has four

  expect_error(Setup_Sim_Fishing(sim_list, fish_sel_input = sel_at_age, FishLenComps_sel = "length"), "fish_sel_l_input must be")
  expect_error(Setup_Sim_Fishing(sim_list, fish_sel_input = sel_at_age, FishLenComps_sel = "length", fish_sel_l_input = wrong_dims), "fish_sel_l_input must be")
  expect_error(Setup_Sim_Fishing(sim_list, fish_sel_input = sel_at_age, FishLenComps_sel = "lengths"), "one of age or length")
  expect_error(Setup_Sim_Survey(sim_list, srv_sel_input = sel_at_age, SrvLenComps_sel = "length"), "srv_sel_l_input must be")

})

test_that("selectivity at length forms catch, index and age-at-length as the fit does", {

  # pcod selects both its fishery and its survey at length
  data("sgl_rg_ebs_pcod_data", envir = environment())
  dat <- sgl_rg_ebs_pcod_data
  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)
  obj <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  sim_list <- suppressMessages(condition_closed_loop_simulations(closed_loop_yrs = 1, n_sims = 2, data = obj$data, parameters = obj$parameters,
                                                                 mapping = obj$mapping, sd_rep = list(par.fixed = obj$par, par.random = NULL),
                                                                 rep = obj$rep, random = NULL))
  expect_equal(c(sim_list$fish_len_comp_sel, sim_list$srv_len_comp_sel), c(1, 1))
  n_sim_yrs <- sim_list$n_yrs
  n_lens <- sim_list$n_lens
  caal_yr <- 30

  # age-at-length draws in one year, ten million fish per length bin, for the fishery and the survey
  sim_list$do_fish_caal <- sim_list$do_srv_caal <- TRUE
  sim_list$comp_fish_caal_like <- sim_list$comp_srv_caal_like <- 0 # multinomial
  sim_list$Fish_caal_Type <- array(999, dim = c(n_sim_yrs, 1))
  sim_list$Fish_caal_Type[caal_yr, 1] <- 1
  sim_list$Srv_caal_Type <- sim_list$Fish_caal_Type
  sim_list$ISS_Fish_caal <- sim_list$ISS_Srv_caal <- array(1e7, dim = c(1, n_sim_yrs, 1, n_lens, 1, 1, 2))
  sim_list$ln_Fish_caal_theta <- sim_list$ln_Srv_caal_theta <- array(0, dim = c(1, 1, 1))
  sim_list$ln_Fish_caal_theta_agg <- sim_list$ln_Srv_caal_theta_agg <- 0

  # the second replicate keeps fish by length and kills half of what it discards
  ret_at_len <- plogis((dat$lens - 55) / 4)
  sim_list$ret_sel_l <- array(1, dim = c(1, n_sim_yrs, n_lens, 1, 1, 2)) # [region, year, len, sex, fleet, sim]
  for(y in 1:n_sim_yrs) sim_list$ret_sel_l[1,y,,1,1,2] <- ret_at_len
  sim_list$dmr[,,,,2] <- 0.5

  set.seed(2)
  sim_env <- Setup_sim_env(sim_list)
  for(sim in 1:2) for(y in 1:caal_yr) run_annual_cycle(y, sim, sim_env)

  # conditioned on the fit, catch and index at length are the fit's in every year
  yrs <- 1:caal_yr
  expect_equal(sim_env$CAL[1,1,yrs,1,,1,1,1], obj$rep$CAL[1,1,yrs,1,,1,1], tolerance = 1e-10)
  expect_equal(sim_env$SrvIAL[1,1,yrs,1,,1,1,1], obj$rep$SrvIAL[1,1,yrs,1,,1,1], tolerance = 1e-10)

  # an age-at-length row is the fish available at each age in that length bin, read through the ageing error.
  # selectivity at age would weight the ages by how selected they are, which moves the rows a long way
  norm_rows <- function(m) m / rowSums(m)
  y <- caal_yr
  n_at_age <- sim_env$NAA[1,1,y,1,,1,1]
  z_at_age <- sim_env$ZAA[1,1,y,1,,1,1]
  fish_key <- sim_env$SizeAgeTrans_fish[1,1,y,1,,,1,1,1]
  fish_rows <- norm_rows(sim_env$ObsFish_caal[1,y,1,,,1,1,1])
  fish_by_length <- norm_rows((fish_key * rep(n_at_age * (1 - exp(-z_at_age)) / z_at_age, each = n_lens)) %*% sim_env$AgeingError_fish[y,,,1,1])
  fish_by_age <- norm_rows((fish_key * rep(sim_env$CAA[1,1,y,1,,1,1,1], each = n_lens)) %*% sim_env$AgeingError_fish[y,,,1,1])
  expect_lt(max(abs(fish_rows - fish_by_length)), 0.002)
  expect_gt(max(abs(fish_rows - fish_by_age)), 0.1)

  srv_key <- sim_env$SizeAgeTrans_srv[1,1,y,1,,,1,1,1]
  srv_rows <- norm_rows(sim_env$ObsSrv_caal[1,y,1,,,1,1,1])
  srv_by_length <- norm_rows((srv_key * rep(n_at_age * exp(-sim_env$t_srv[1,1,1] * z_at_age), each = n_lens)) %*% sim_env$AgeingError_srv[y,,,1,1])
  srv_by_age <- norm_rows((srv_key * rep(sim_env$SrvIAA[1,1,y,1,,1,1,1], each = n_lens)) %*% sim_env$AgeingError_srv[y,,,1,1])
  expect_lt(max(abs(srv_rows - srv_by_length)), 0.002)
  expect_gt(max(abs(srv_rows - srv_by_age)), 0.1)

  # retention at length: what is selected at each length is kept or discarded by length, and half the discards die
  n_at_age <- sim_env$NAA[1,1,y,1,,1,2]
  z_at_age <- sim_env$ZAA[1,1,y,1,,1,2]
  selected <- as.vector(sim_env$SizeAgeTrans_fish[1,1,y,1,,,1,1,2] %*% (n_at_age * (1 - exp(-z_at_age)) / z_at_age)) *
    sim_env$Fmort[1,y,1,1,2] * sim_env$fish_sel_l[1,y,,1,1,2]
  expect_equal(sim_env$CAL[1,1,y,1,,1,1,2], selected * ret_at_len, tolerance = 1e-10)
  expect_equal(sim_env$DAL[1,1,y,1,,1,1,2], selected * (1 - ret_at_len) * 0.5, tolerance = 1e-10)
  expect_gt(sum(sim_env$DAL[1,1,y,1,,1,1,2]), 0)

})
