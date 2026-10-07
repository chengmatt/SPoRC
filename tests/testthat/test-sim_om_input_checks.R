# Operating model inputs checked and drawn the way the estimation model reads them: the population composition
# likelihoods are length checked, and an index whose sd is partly estimated is drawn at the fit's sd.

library(SPoRC)
library(testthat)

test_that("the population composition likelihoods are length checked", {

  expect_error(check_sim_dimensions(rep(0, 2), n_fish_fleets = 3, what = "comp_fishage_pop_like"), "n_fish_fleets")
  expect_error(check_sim_dimensions(rep(0, 2), n_fish_fleets = 3, what = "comp_fishlen_discard_pop_like"), "n_fish_fleets")
  expect_error(check_sim_dimensions(rep(0, 2), n_srv_fleets = 1, what = "comp_srvage_pop_like"), "n_srv_fleets")
  expect_null(check_sim_dimensions(rep(0, 3), n_fish_fleets = 3, what = "comp_fishage_pop_like"))
})

# a self test of one survey with reported errors of 0.2 and a likelihood weight of 4, keeping the operating model
idx_weight_self_test <- function(sigma_spec, n_sims) {

  il <- suppressWarnings(suppressMessages(sweep_input(dims = list(n_regions = 1, n_yrs = 15),
                                                      srvidx = list(sigmaSrvIdx_spec = sigma_spec, ln_sigmaSrvIdx = log(0.3)))))
  il$data$Wt_SrvIdx[] <- 4
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  sim_file <- tempfile(fileext = ".RDS")
  set.seed(9)
  suppressWarnings(suppressMessages(simulation_self_test(data = obj$data, parameters = obj$parameters, mapping = obj$mapping,
                                                         random = NULL, rep = obj$rep,
                                                         sd_rep = list(par.fixed = obj$par, par.random = NULL),
                                                         n_sims = n_sims, newton_loops = 0, output_path = sim_file, what = "SrvIdx_nLL")))
  list(om = readRDS(sim_file), use = il$data$UseSrvIdx == 1)
}

test_that("an index whose sd is partly estimated is drawn at the fit's sd whatever its weight", {

  # the weight scales the likelihood, so it cancels from the estimated sd: drawn at 0.2 + 0.3, not 0.2 / 2 + 0.3
  run <- idx_weight_self_test("est_additive", n_sims = 40)
  dev <- log(run$om$ObsSrvIdx[,,,1,] / run$om$TrueSrvIdx[,,,1,])
  expect_equal(stats::sd(as.vector(dev)), 0.5, tolerance = 0.06)
  expect_true(all(run$om$ObsSrvIdx_SE[run$use] == 0.2)) # the reported errors the refit reads

  # a fixed sd still takes its weight, as before
  run <- idx_weight_self_test("fix", n_sims = 1)
  expect_true(all(run$om$ObsSrvIdx_SE[run$use] == 0.1))
})
