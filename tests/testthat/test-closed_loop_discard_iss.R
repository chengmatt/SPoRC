# discard input sample sizes in the closed loop, from condition_closed_loop_simulations() through to
# the operating model's composition draws

data("dusky_rtmb_model")

test_that("supplied discard input sample sizes reach the returned sim_list", {
  # an array given for a discard ISS is the complete input, and the closed loop used to read it back
  # from the *_fill name, which is never set, so Setup_Sim_Fishing stopped on a NULL array
  m <- dusky_rtmb_model
  closed_loop_yrs <- 2
  n_sims <- 2
  iss_dim <- c(m$data$n_regions, length(m$data$years) + closed_loop_yrs, m$data$n_seas,
               m$data$n_sexes, m$data$n_fish_fleets, n_sims) # region, year, season, sex, fleet, sim
  age_iss <- array(10 + seq_len(prod(iss_dim)), dim = iss_dim) # a different value in every cell, so a reshaped or generated array cannot match
  len_iss <- age_iss + 1000 # differs from the age array, so swapping the two fails

  # passed through the *_fill arguments, since ISS_FishAgeComps_discard = in the call partially
  # matches both the _discard_fill and _discard_pop_fill formals and R stops before the body runs
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = n_sims,
    m$data,
    m$parameters,
    m$mapping,
    sd_rep = m$sdrep,
    rep = m$rep,
    random = NULL,
    ISS_FishAgeComps_discard_fill = age_iss,
    ISS_FishLenComps_discard_fill = len_iss
  )

  expect_identical(cl$ISS_FishAgeComps_discard, age_iss)
  expect_identical(cl$ISS_FishLenComps_discard, len_iss)
})


test_that("population-specific discard sample sizes follow the F_pattern fill", {
  # the closed loop stores this flag as *_discard_pop_fill, and the operating model looked it up as
  # *_pop_discard_fill, so these sample sizes stayed at zero over the projection years
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
    ISS_FishAgeComps_discard_pop_fill = "F_pattern",
    ISS_FishLenComps_discard_pop_fill = "F_pattern"
  )

  # the dusky fit retains every fish, so discarding is switched on here, since the sample sizes are
  # only updated in years with discards
  cl$ret_sel[] <- 0.7 # proportion retained at every age
  cl$dmr[] <- 0.5 # discard mortality rate
  iss_hist <- 20 + 5 * (hist %% 7) # conditioning sample sizes running from 20 to 50
  cl$ISS_FishAgeComps_discard_pop[1, 1, hist, 1, 1, 1, ] <- iss_hist
  cl$ISS_FishLenComps_discard_pop[1, 1, hist, 1, 1, 1, ] <- iss_hist

  set.seed(3)
  sim_env <- Setup_sim_env(cl)
  sim_env$Fmort[1, proj, 1, 1, ] <- c(0.01, 0.03, 0.2, 0.06) # set in place of a management procedure; 0.2 is above every earlier F
  for(s in seq_len(sim_env$n_sims)) {
    for(y in seq_len(sim_env$n_yrs)) run_annual_cycle(y, s, sim_env)
  } # end s loop

  # each projection year sits between the smallest and largest earlier sample size, at F over the
  # largest earlier F, capped at one
  fmort <- sim_env$Fmort[1, , 1, 1, 1]
  expected <- rep(NA_real_, closed_loop_yrs)
  for(i in seq_len(closed_loop_yrs)) {
    y <- proj[i]
    expected[i] <- 20 + min(fmort[y] / max(fmort[1:(y - 1)]), 1) * (50 - 20)
  } # end i loop
  expect_gt(length(unique(round(expected, 6))), 2) # varies by year, so a constant fill cannot match

  for(s in seq_len(sim_env$n_sims)) {
    expect_equal(sim_env$ISS_FishAgeComps_discard_pop[1, 1, proj, 1, 1, 1, s], expected, tolerance = 1e-10,
                 label = sprintf("projected discard age ISS, simulation %d", s))
    expect_equal(sim_env$ISS_FishLenComps_discard_pop[1, 1, proj, 1, 1, 1, s], expected, tolerance = 1e-10,
                 label = sprintf("projected discard length ISS, simulation %d", s))
  } # end s loop
})
