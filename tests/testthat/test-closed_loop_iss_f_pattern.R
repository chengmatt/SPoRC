# F_pattern input sample sizes in the closed loop: each projection year's sample size is set between
# the smallest and largest earlier one by that year's F over the largest earlier F, capped at one

data("dusky_rtmb_model")

test_that("population-specific fishery sample sizes follow their own F_pattern flag", {
  # the population-specific branch tested the pooled flag, so F_pattern set on the population sample
  # sizes alone stopped the run with "object 'ISS_FishAgeComps_fill' not found"
  m <- dusky_rtmb_model
  closed_loop_yrs <- 4
  n_hist <- length(m$data$years)
  hist <- seq_len(n_hist) # conditioning years
  proj <- n_hist + seq_len(closed_loop_yrs) # projection years
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = 2,
    m$data,
    m$parameters,
    m$mapping,
    sd_rep = m$sdrep,
    rep = m$rep,
    random = NULL,
    ISS_FishAgeComps_pop_fill = "F_pattern",
    ISS_FishLenComps_pop_fill = "F_pattern"
  )
  iss_hist <- 20 + 5 * (hist %% 7) # conditioning sample sizes running from 20 to 50
  cl$ISS_FishAgeComps_pop[1, 1, hist, 1, 1, 1, ] <- iss_hist
  cl$ISS_FishLenComps_pop[1, 1, hist, 1, 1, 1, ] <- iss_hist
  pooled_before <- cl$ISS_FishAgeComps # the pooled fill stays "mean", so the run should leave it alone

  set.seed(3)
  sim_env <- Setup_sim_env(cl)
  sim_env$Fmort[1, proj, 1, 1, ] <- c(0.01, 0.03, 0.2, 0.06) # set in place of a management procedure; 0.2 is above every earlier F
  for(s in seq_len(sim_env$n_sims)) {
    for(y in seq_len(sim_env$n_yrs)) run_annual_cycle(y, s, sim_env)
  } # end s loop

  fmort <- sim_env$Fmort[1, , 1, 1, 1]
  expected <- rep(NA_real_, closed_loop_yrs)
  for(i in seq_len(closed_loop_yrs)) {
    y <- proj[i]
    expected[i] <- 20 + min(fmort[y] / max(fmort[1:(y - 1)]), 1) * (50 - 20)
  } # end i loop
  expect_gt(length(unique(round(expected, 6))), 2) # varies by year, so a constant fill cannot match

  for(s in seq_len(sim_env$n_sims)) {
    expect_equal(sim_env$ISS_FishAgeComps_pop[1, 1, proj, 1, 1, 1, s], expected, tolerance = 1e-10,
                 label = sprintf("projected population age ISS, simulation %d", s))
    expect_equal(sim_env$ISS_FishLenComps_pop[1, 1, proj, 1, 1, 1, s], expected, tolerance = 1e-10,
                 label = sprintf("projected population length ISS, simulation %d", s))
  } # end s loop
  expect_identical(sim_env$ISS_FishAgeComps, pooled_before)
})


test_that("F_pattern sample sizes are predicted in every season", {
  # the prediction copied its history into season seas of an array with a single season, so any
  # season after the first stopped with "subscript out of bounds"
  n_hist <- 15
  closed_loop_yrs <- 3
  il <- sweep_input(dims = list(
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_yrs = n_hist,
    n_ages = 6,
    n_seas = 2
  ))
  # evaluated at its starting values rather than fitted: conditioning reads only the parameter
  # values, and what is checked here is the season dim, not the estimates
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = 1,
    il$data,
    il$par,
    il$map,
    sd_rep = list(par.fixed = obj$par),
    rep = obj$rep,
    random = NULL,
    ISS_FishAgeComps_fill = "F_pattern"
  )
  hist <- seq_len(n_hist) # conditioning years
  proj <- n_hist + seq_len(closed_loop_yrs) # projection years
  iss_lo <- c(20, 100) # smallest conditioning sample size in each season
  iss_hi <- c(50, 140) # largest, with no overlap, so a season reading the other's history cannot match
  cl$ISS_FishAgeComps[1, hist, 1, 1, 1, 1] <- 20 + 5 * (hist %% 7)
  cl$ISS_FishAgeComps[1, hist, 2, 1, 1, 1] <- 100 + 10 * (hist %% 5)

  set.seed(3)
  sim_env <- Setup_sim_env(cl)
  sim_env$Fmort[1, proj, 1, 1, 1] <- c(0.01, 0.5, 0.02) # season one
  sim_env$Fmort[1, proj, 2, 1, 1] <- c(0.03, 0.01, 0.9) # season two, capped in a different year
  for(y in seq_len(sim_env$n_yrs)) run_annual_cycle(y, 1, sim_env)

  for(k in 1:2) {
    fmort <- sim_env$Fmort[1, , k, 1, 1] # this season's F only
    expected <- rep(NA_real_, closed_loop_yrs)
    for(i in seq_len(closed_loop_yrs)) {
      y <- proj[i]
      expected[i] <- iss_lo[k] + min(fmort[y] / max(fmort[1:(y - 1)]), 1) * (iss_hi[k] - iss_lo[k])
    } # end i loop
    expect_equal(sim_env$ISS_FishAgeComps[1, proj, k, 1, 1, 1], expected, tolerance = 1e-10,
                 label = sprintf("projected age ISS, season %d", k))
  } # end k loop
})
