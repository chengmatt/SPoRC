# Index totals over a subset of ages (srv_idx_ages / fish_idx_ages) and fishery index timing in the operating
# model: conditioned on the EBS Pacific cod fit with both indices counting age 1 only, as the EBS pollock acoustic
# index does, the operating model's indices equal the fit's predictions in every year.

library(SPoRC)
library(testthat)

test_that("the operating model counts only the index ages, at the fit's index timing", {

  data("sgl_rg_ebs_pcod_data", envir = environment())
  dat <- sgl_rg_ebs_pcod_data
  n_yrs <- length(dat$years)
  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)

  # both indices count age 1 only, and the fishery index is abundance at mid season so the fit predicts it
  input_list$data$srv_idx_ages[] <- 0
  input_list$data$srv_idx_ages[2, ] <- 1
  input_list$data$fish_idx_ages[] <- 0
  input_list$data$fish_idx_ages[2, ] <- 1
  input_list$data$fish_idx_type[] <- 0
  expect_equal(as.numeric(input_list$data$t_fish), 0.5)

  obj <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  sim_list <- suppressMessages(condition_closed_loop_simulations(closed_loop_yrs = 1, n_sims = 1, data = obj$data, parameters = obj$parameters,
                                                                 mapping = obj$mapping, sd_rep = list(par.fixed = obj$par, par.random = NULL),
                                                                 rep = obj$rep, random = NULL))
  expect_equal(sim_list$srv_idx_ages, obj$data$srv_idx_ages)
  sim_env <- Setup_sim_env(sim_list) # built here, since the annual cycle finds the environment by this name in the frame that built it
  for(y in 1:n_yrs) run_annual_cycle(y, 1, sim_env)

  expect_equal(sim_env$TrueSrvIdx[1,1:n_yrs,1,1,1], obj$rep$PredSrvIdx[1,1,,1,1], tolerance = 1e-10)
  expect_equal(sim_env$TrueFishIdx[1,1:n_yrs,1,1,1], obj$rep$PredFishIdx[1,1,,1,1], tolerance = 1e-10)

  # counting every age instead would have made the survey index several times larger
  all_ages <- sim_env$srv_q[1,1:n_yrs,1,1] * rowSums(sim_env$SrvIAA[1,1,1:n_yrs,1,,1,1,1])
  expect_gt(min(all_ages / obj$rep$PredSrvIdx[1,1,,1,1]), 2)

  # setup refuses an index age array of the wrong shape
  expect_error(Setup_Sim_Survey(sim_list, srv_sel_input = sim_list$srv_sel, srv_idx_ages = matrix(1, 3, 1)), "srv_idx_ages must be")

})
