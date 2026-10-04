# catch_to_F_om: the F it solves takes the target catch from the operating model once the year is run,
# by region, season and fleet, under rec_lag 0, in numbers, and capped where no F reaches the target.
# Relative catch within 1e-5.

library(SPoRC)
library(testthat)

# two regions, two seasons, two fleets. fleet 2 takes age 1, discards small fish and records numbers
cf_sim_list <- function(move_timing = 0, rec_lag = 1, n_regions = 2, n_seas = 2, n_fish_fleets = 2, n_yrs = 6, n_ages = 6) {

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = n_regions, n_ages = n_ages, n_lens = NULL,
                            n_sexes = 1, n_fish_fleets = n_fish_fleets, n_srv_fleets = 1, n_seas = n_seas,
                            seasdur = if(n_seas == 1) 1 else c(0.4, 0.6), n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  waa <- 5 / (1 + exp(-0.5 * ((1:n_ages) - 3)))

  fish_sel <- ret_sel <- WAA_fish <- array(1, dim = c(1, n_regions, n_yrs, n_seas, n_ages, 1, n_fish_fleets))
  dmr <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets, 1))
  for(f in 1:n_fish_fleets) {
    for(a in 1:n_ages) {
      fish_sel[,,,,a,,f] <- 1 / (1 + exp(-1.2 * (a - (if(f == 1) 4 else 1.5))))
      WAA_fish[,,,,a,,f] <- waa[a] * (if(f == 2) 1.1 else 1)
      if(f == 2) ret_sel[,,,,a,,f] <- 1 / (1 + exp(-1.5 * (a - 3)))
    } # end a loop
  } # end f loop
  if(n_fish_fleets == 2) dmr[,,,2,] <- 0.5

  sim_list <- Setup_Sim_Fishing(sim_list = sim_list,
                                fish_sel_input = replicate(1, fish_sel),
                                ret_sel_input = replicate(1, ret_sel),
                                dmr_input = dmr,
                                catch_units = if(n_fish_fleets == 2) c(1, 0) else 1,
                                Fmort_input = array(0.1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets, 1)))
  sim_list <- Setup_Sim_Survey(sim_list = sim_list, srv_sel_input = replicate(1, fish_sel[,,,,,,1, drop = FALSE]))

  biol <- array(0, dim = c(1, n_regions, n_yrs, n_seas, n_ages, 1))
  WAA <- MatAA <- biol
  for(a in 1:n_ages) {
    WAA[,,,,a,] <- waa[a]
    MatAA[,,,,a,] <- 1 / (1 + exp(-1 * (a - 3)))
  } # end a loop
  sim_list <- suppressWarnings(Setup_Sim_Biologicals(sim_list = sim_list,
                                                     natmort_input = replicate(1, biol + 0.2),
                                                     WAA_input = replicate(1, WAA),
                                                     WAA_fish_input = replicate(1, WAA_fish),
                                                     WAA_srv_input = replicate(1, WAA_fish[,,,,,,1, drop = FALSE]),
                                                     MatAA_input = replicate(1, MatAA)))
  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)

  rec_seas_prop <- array(0, dim = c(1, n_seas, 1))
  rec_seas_prop[,,1] <- if(n_seas == 1) 1 else c(0.6, 0.4)
  R0 <- array(0, dim = c(1, n_regions, n_yrs, 1))
  for(r in 1:n_regions) R0[,r,,1] <- c(10, 6)[r]
  sim_list <- Setup_Sim_Rec(sim_list = sim_list,
                            R0_input = R0,
                            ln_sigmaR = array(log(0.5), dim = c(2, 1, n_regions)),
                            rec_seas_prop_input = rec_seas_prop,
                            recruitment_opt = "bh_rec", # recruits depend on spawning biomass under rec_lag 0
                            t_spawn = if(rec_lag == 0) 0.5 else 0,
                            rec_lag = rec_lag)

  # recruits stay put; Movement for timings 0 and 1, rates for timing 2
  Movement <- Mrate <- array(0, dim = c(1, n_regions, n_regions, n_yrs, n_seas, n_ages, 1, 1))
  for(a in 1:n_ages) {
    for(y in 1:n_yrs) {
      for(seas in 1:n_seas) {
        if(n_regions == 1 || a == 1) Movement[1,,,y,seas,a,1,1] <- diag(n_regions)
        else {
          Movement[1,,,y,seas,a,1,1] <- matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
          Mrate[1,,,y,seas,a,1,1] <- matrix(c(-0.25, 0.35, 0.25, -0.35), 2, 2)
        }
      } # end seas loop
    } # end y loop
  } # end a loop
  sim_list$Movement <- Movement
  sim_list$Mrate <- Mrate
  sim_list$move_timing <- move_timing
  sim_list$expm_nsub <- 0

  return(sim_list)
}

# retained catch by region, season and fleet in year y, in biomass or numbers
cf_realized <- function(sim_env, y, biomass = TRUE) {
  out <- array(0, dim = c(sim_env$n_regions, sim_env$n_seas, sim_env$n_fish_fleets))
  for(r in 1:sim_env$n_regions) {
    for(seas in 1:sim_env$n_seas) {
      for(f in 1:sim_env$n_fish_fleets) {
        wt <- if(biomass) sim_env$WAA_fish[,r,y,seas,,,f,1] else 1
        out[r,seas,f] <- sum(sim_env$CAA[,r,y,seas,,,f,1] * wt)
      } # end f loop
    } # end seas loop
  } # end r loop
  return(out)
}

test_that("the solved F takes biomass targets by region, season and fleet", {

  for(move_timing in c(0, 2)) {

    set.seed(21)
    sim_env <- suppressMessages(Setup_sim_env(cf_sim_list(move_timing = move_timing)))
    for(y in 1:4) run_annual_cycle(y, 1, sim_env)
    target <- cf_realized(sim_env, 4) * array(seq(0.5, 1.8, length.out = 8), dim = c(2, 2, 2))
    target[1,2,1] <- 0 # closed to fishing

    sol <- catch_to_F_om(target = target, y = 5, sim = 1, sim_env = sim_env)
    sim_env$Fmort[,5,,,1] <- sol$Fmort
    run_annual_cycle(5, 1, sim_env)
    realized <- cf_realized(sim_env, 5)

    pos <- target > 0
    expect_lt(max(abs(realized[pos] / target[pos] - 1)), 1e-5)
    expect_equal(sol$Fmort[1,2,1], 0)
    expect_equal(realized[1,2,1], 0)
    expect_lt(max(abs(sol$resid)), 1e-5)

  } # end move_timing loop

})

test_that("under rec_lag 0 the solve counts the recruits the year's spawning biomass produces", {

  set.seed(22)
  sim_env <- suppressMessages(Setup_sim_env(cf_sim_list(rec_lag = 0)))
  for(y in 1:4) run_annual_cycle(y, 1, sim_env)
  target <- cf_realized(sim_env, 4) * 1.2

  sol <- catch_to_F_om(target = target, y = 5, sim = 1, sim_env = sim_env)
  sim_env$Fmort[,5,,,1] <- sol$Fmort
  run_annual_cycle(5, 1, sim_env)

  expect_gt(sum(sim_env$CAA[,,5,,1,,2,1]), 0) # fleet 2 does catch this year's recruits
  expect_lt(max(abs(cf_realized(sim_env, 5) / target - 1)), 1e-5)

})

test_that("targets in numbers are taken in numbers", {

  set.seed(23)
  sim_env <- suppressMessages(Setup_sim_env(cf_sim_list()))
  for(y in 1:4) run_annual_cycle(y, 1, sim_env)
  target <- array(0, dim = c(2, 2, 2))
  target[,,1] <- cf_realized(sim_env, 4)[,,1] # fleet 1 in biomass
  target[,,2] <- cf_realized(sim_env, 4, biomass = FALSE)[,,2] * 0.8 # fleet 2 in numbers

  sol <- catch_to_F_om(target = target, y = 5, sim = 1, sim_env = sim_env, target_units = c("biom", "abd"))
  sim_env$Fmort[,5,,,1] <- sol$Fmort
  run_annual_cycle(5, 1, sim_env)

  expect_lt(max(abs(cf_realized(sim_env, 5)[,,1] / target[,,1] - 1)), 1e-5)
  expect_lt(max(abs(sim_env$TrueCatch[,5,,2,1] / target[,,2] - 1)), 1e-5) # fleet 2 records numbers

})

test_that("a target no F reaches is capped at catch_f_max with a warning", {

  set.seed(24)
  sim_env <- suppressMessages(Setup_sim_env(cf_sim_list()))
  for(y in 1:4) run_annual_cycle(y, 1, sim_env)
  target <- cf_realized(sim_env, 4)
  target[2,1,2] <- 1e3 * target[2,1,2]

  expect_warning(sol <- catch_to_F_om(target = target, y = 5, sim = 1, sim_env = sim_env, catch_f_max = 3), "not reachable")
  expect_equal(sol$Fmort[2,1,2], 3)
  expect_lt(sol$resid[2,1,2], -0.5) # undershot
  expect_lt(max(abs(sol$resid[-6])), 1e-5) # every other cell still hits its target

})

test_that("one region, season and fleet matches the Baranov equation solved by hand", {

  set.seed(25)
  sim_env <- suppressMessages(Setup_sim_env(cf_sim_list(n_regions = 1, n_seas = 1, n_fish_fleets = 1)))
  for(y in 1:4) run_annual_cycle(y, 1, sim_env)
  target <- 0.8 * cf_realized(sim_env, 4)[1,1,1]
  sol <- catch_to_F_om(target = target, y = 5, sim = 1, sim_env = sim_env)

  N <- sim_env$NAA[1,1,5,1,,1,1]
  sel <- sim_env$fish_sel[1,1,5,1,,1,1,1]
  M <- sim_env$natmort[1,1,5,1,,1,1]
  W <- sim_env$WAA_fish[1,1,5,1,,1,1,1]
  baranov <- function(F) sum(F * sel / (F * sel + M) * N * (1 - exp(-(F * sel + M))) * W) - target
  F_hand <- stats::uniroot(baranov, c(0, 5), tol = 1e-12)$root
  expect_equal(sol$Fmort[1,1,1], F_hand, tolerance = 1e-5)

})
