# Under rec_lag 0 the operating model draws next year's recruitment deviations at the end of each year:
# the drawn deviation is the one used, random numbers used in between leave it alone, and spawning
# biomass from trial numbers matches the stored one. Exact equality throughout.

library(SPoRC)
library(testthat)

# one region, two seasons, Beverton-Holt recruits from the same year's spawning biomass in season 2
rec_lag0_sim_list <- function(n_yrs = 8, n_ages = 6) {

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = NULL,
                            n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 2, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  age_curve <- function(v, dims) array(rep(v, each = prod(dims[1:4])), dim = dims)
  sel <- age_curve(1 / (1 + exp(-1.5 * ((1:n_ages) - 2))), c(1, 1, n_yrs, 2, n_ages, 1, 1))
  waa <- 5 / (1 + exp(-1 * ((1:n_ages) - 3)))

  sim_list <- Setup_Sim_Fishing(sim_list = sim_list,
                                fish_sel_input = replicate(1, sel),
                                Fmort_input = array(0.1, dim = c(1, n_yrs, 2, 1, 1)))
  sim_list <- Setup_Sim_Survey(sim_list = sim_list, srv_sel_input = replicate(1, sel))
  sim_list <- suppressWarnings(Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.25, dim = c(1, 1, n_yrs, 2, n_ages, 1))),
    WAA_input = replicate(1, age_curve(waa, c(1, 1, n_yrs, 2, n_ages, 1))),
    WAA_fish_input = replicate(1, age_curve(waa, c(1, 1, n_yrs, 2, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, age_curve(waa, c(1, 1, n_yrs, 2, n_ages, 1, 1))),
    MatAA_input = replicate(1, age_curve(c(0, 1 / (1 + exp(-2 * ((2:n_ages) - 3)))), c(1, 1, n_yrs, 2, n_ages, 1)))
  ))
  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 2, n_ages, 1, 1))

  rec_seas_prop <- array(0, dim = c(1, 2, 1))
  rec_seas_prop[, 2, ] <- 1 # nothing recruits before the spawning season
  Setup_Sim_Rec(sim_list = sim_list,
                R0_input = array(100, dim = c(1, 1, n_yrs, 1)),
                h_input = array(0.7, dim = c(1, 1, n_yrs, 1)),
                ln_sigmaR = array(log(0.5), dim = c(2, 1, 1)),
                rec_seas_prop_input = rec_seas_prop,
                recruitment_opt = "bh_rec",
                spawn_seas = 2,
                t_spawn = 0.5,
                rec_lag = 0,
                init_age_strc = 1)
}

test_that("under rec_lag 0 next year's deviation is drawn at the end of this year and is the one used", {

  sim_list <- rec_lag0_sim_list()
  set.seed(11)
  sim_env <- suppressMessages(Setup_sim_env(sim_list))
  for(y in 1:4) run_annual_cycle(y, 1, sim_env)

  dev_5 <- sim_env$ln_RecDevs[1,1,5,1]
  expect_true(dev_5 != 0) # drawn before year 5 runs
  expect_equal(sim_env$Rec[1,1,5,1], 0) # recruitment itself still waits for year 5's spawning biomass

  run_annual_cycle(5, 1, sim_env)
  expect_identical(sim_env$ln_RecDevs[1,1,5,1], dev_5)
  # age 1 is immature, so the spawning biomass stored after recruits arrive is the one they came from
  det_rec <- SPoRC:::sim_det_recruitment(5, 1, sim_env, SSB_vals = array(sim_env$SSB[,,,1], dim = c(1, 1, sim_env$n_yrs)))
  expect_equal(sim_env$Rec[1,1,5,1], det_rec[1,1] * exp(dev_5), tolerance = 1e-12)

})

test_that("random numbers used between years leave rec_lag 0 recruitment unchanged", {

  sim_list <- rec_lag0_sim_list()

  set.seed(12)
  sim_env <- suppressMessages(Setup_sim_env(sim_list))
  for(y in 1:5) run_annual_cycle(y, 1, sim_env)
  rec_plain <- sim_env$Rec[1,1,5,1]

  # a management procedure drawing random numbers after year 4 is assessed
  set.seed(12)
  sim_env <- suppressMessages(Setup_sim_env(sim_list))
  for(y in 1:4) run_annual_cycle(y, 1, sim_env)
  invisible(stats::runif(25))
  run_annual_cycle(5, 1, sim_env)
  expect_identical(sim_env$Rec[1,1,5,1], rec_plain)

})

test_that("spawning biomass from trial numbers and mortality matches the stored one", {

  sim_list <- rec_lag0_sim_list()
  set.seed(13)
  sim_env <- suppressMessages(Setup_sim_env(sim_list))
  for(y in 1:5) run_annual_cycle(y, 1, sim_env)

  season_dims <- c(1, 1, 1, 1, sim_env$n_ages, 1)
  stored <- SPoRC:::compute_biom_y_sim(5, 2, 1, sim_env)
  trial <- SPoRC:::compute_biom_y_sim(5, 2, 1, sim_env,
                                      NAA_s = array(sim_env$NAA[,,5,2,,,1], dim = season_dims),
                                      ZAA_s = array(sim_env$ZAA[,,5,2,,,1], dim = season_dims))
  expect_identical(trial, stored)

  # lower total mortality before spawning leaves more spawners
  lower_Z <- SPoRC:::compute_biom_y_sim(5, 2, 1, sim_env, ZAA_s = array(sim_env$ZAA[,,5,2,,,1] * 0.5, dim = season_dims))
  expect_gt(lower_Z$SSB_y[1,1], stored$SSB_y[1,1])

})
