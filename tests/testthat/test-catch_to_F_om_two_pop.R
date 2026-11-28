# catch_to_F_om on an operating model with two populations in two regions: the F it solves takes the target catch
# by region once the year is run, with both populations contributing. Relative catch within 1e-6.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

# a closed loop conditioned on the small two-population model at its starting values, run through the fitted years
two_pop_sim_env <- function() {

  model <- suppressWarnings(suppressMessages(fitted_small_model(nr = 2, nx = 1, np = 2)))

  # the report at the starting values, where both populations are the same size (the fit empties one of them)
  start <- suppressWarnings(suppressMessages(fit_model(model$il$data, model$il$par, model$il$map, do_optim = FALSE, silent = TRUE)))
  closed_loop <- suppressWarnings(suppressMessages(condition_closed_loop_simulations(
    closed_loop_yrs = 3,
    n_sims = 1,
    model$il$data,
    model$il$par,
    model$il$map,
    sd_rep = model$sdrep,
    rep = start$rep,
    random = NULL
  )))
  sim_env <- suppressWarnings(suppressMessages(Setup_sim_env(closed_loop)))

  n_fit_yrs <- length(model$il$data$years)
  set.seed(31)
  for(y in 1:n_fit_yrs) run_annual_cycle(y, 1, sim_env)

  list(sim_env = sim_env, n_fit_yrs = n_fit_yrs)
}

# retained catch in biomass by region in year y, for one population or summed over both
two_pop_catch <- function(sim_env, y, pops) {
  catch <- rep(0, sim_env$n_regions)
  for(r in 1:sim_env$n_regions) {
    for(p in pops) catch[r] <- catch[r] + sum(sim_env$CAA[p,r,y,1,,,1,1] * sim_env$WAA_fish[p,r,y,1,,,1,1])
  } # end r loop
  catch
}

# Tests ----------------------------------------------------------------------

test_that("the solved F takes the target catch by region with two populations in the water", {

  om <- two_pop_sim_env()
  sim_env <- om$sim_env
  y <- om$n_fit_yrs + 1 # the first year past the fit

  # both populations are caught in both regions, so neither can be left out of the solve
  for(p in 1:2) expect_gt(min(two_pop_catch(sim_env, om$n_fit_yrs, pops = p)), 0.2 * max(two_pop_catch(sim_env, om$n_fit_yrs, pops = 1:2)))

  # 60 pct of last year's catch in region 1 and 140 pct in region 2
  last_catch <- two_pop_catch(sim_env, om$n_fit_yrs, pops = 1:2)
  target <- array(last_catch * c(0.6, 1.4), dim = c(sim_env$n_regions, 1, 1))

  solved <- catch_to_F_om(target = target, y = y, sim = 1, sim_env = sim_env)
  sim_env$Fmort[,y,,,1] <- solved$Fmort
  run_annual_cycle(y, 1, sim_env)

  expect_equal(two_pop_catch(sim_env, y, pops = 1:2), as.numeric(target), tolerance = 1e-6)
  expect_lt(max(abs(solved$resid)), 1e-6)

  # the two regions were asked for different changes, so their F moved apart
  expect_gt(abs(solved$Fmort[2,1,1] / solved$Fmort[1,1,1] - 1), 0.5)

})
