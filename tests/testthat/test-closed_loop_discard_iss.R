# discard input sample sizes supplied to condition_closed_loop_simulations()
#
# an array given for a discard ISS is the complete input, and the closed loop used to read it back
# from the *_fill name, which is never set, so Setup_Sim_Fishing stopped on a NULL array

data("dusky_rtmb_model")

test_that("supplied discard input sample sizes reach the returned sim_list", {
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
