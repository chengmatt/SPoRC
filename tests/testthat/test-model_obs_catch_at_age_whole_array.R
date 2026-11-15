# Catch and discards at age are one product each over every slice, formed ahead of the observation loops. Checked slice
# by slice against the Baranov catch equation on unequal dims, with numbers at age carrying the model's extra year, and
# under move_timing = 2 against the season-integrated abundance, passed in or built by the function itself.

library(SPoRC)
library(testthat)

n_pop <- 2
n_regions <- 3
n_yrs <- 5
n_seas <- 2
n_ages <- 7
n_sexes <- 2
n_fish_fleets <- 3
dim_age <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes) # at age
dim_fish_age <- c(dim_age, n_fish_fleets) # at age by fleet
dim_fish <- c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets) # by fleet

set.seed(11)
ret_FAA <- array(runif(prod(dim_fish_age), 0.01, 0.3), dim = dim_fish_age)
disc_FAA <- array(runif(prod(dim_fish_age), 0, 0.05), dim = dim_fish_age)
ZAA <- array(0.2, dim = dim_age) + apply(ret_FAA + disc_FAA, 1:6, sum) # natural mortality plus fishing
NAA <- array(runif(n_pop * n_regions * (n_yrs + 1) * n_seas * n_ages * n_sexes, 10, 1000),
             dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # one extra year, as the model has it

# movement generators, rows summing to zero
Mrate <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes))
for(p in 1:n_pop) {
  for(y in 1:n_yrs) {
    for(seas in 1:n_seas) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          Q <- matrix(runif(n_regions^2, 0, 0.3), n_regions, n_regions)
          diag(Q) <- 0
          diag(Q) <- -rowSums(Q)
          Mrate[p,,,y,seas,a,s] <- Q
        } # end s loop
      } # end a loop
    } # end seas loop
  } # end y loop
} # end p loop

obs_input <- function(move_timing = 0, NAA_int = NULL) {
  list(n_pop = n_pop,
       n_regions = n_regions,
       n_yrs = n_yrs,
       n_seas = n_seas,
       n_fish_fleets = n_fish_fleets,
       n_sexes = n_sexes,
       fish_q_blocks = array(1L, dim = c(n_regions, n_yrs, n_fish_fleets)),
       ln_fish_q = array(log(0.001), dim = c(n_regions, 1, n_fish_fleets)),
       fish_q = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)),
       ret_FAA = ret_FAA,
       disc_FAA = disc_FAA,
       ZAA = ZAA,
       NAA = NAA,
       CAA = array(0, dim = dim_fish_age),
       DAA = array(0, dim = dim_fish_age),
       CAL = array(0, dim = dim_fish_age),
       DAL = array(0, dim = dim_fish_age),
       PredCatch = array(0, dim = dim_fish),
       PredDiscard = array(0, dim = dim_fish),
       PredFishIdx = array(0, dim = dim_fish),
       fit_lengths = 0,
       SizeAgeTrans = NULL,
       catch_units = rep(0, n_fish_fleets),
       discard_units = rep(0, n_fish_fleets),
       WAA_fish = array(1, dim = dim_fish_age),
       dmr = array(1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)),
       fish_idx_type = rep(0, n_fish_fleets),
       fish_sel = array(1, dim = dim_fish_age),
       ret_sel = array(1, dim = dim_fish_age),
       Mrate = Mrate,
       move_timing = move_timing,
       NAA_int = NAA_int)
}

# catch at age a slice at a time from an abundance available over the season
slice_catch <- function(FAA, Avail) {
  expected <- array(0, dim = dim_fish_age)
  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(seas in 1:n_seas) {
          for(f in 1:n_fish_fleets) expected[p,r,y,seas,,,f] <- FAA[p,r,y,seas,,,f] * Avail[p,r,y,seas,,]
        } # end seas loop
      } # end y loop
    } # end r loop
  } # end p loop
  expected
}

test_that("catch and discards at age match the Baranov equation slice by slice", {

  out <- do.call(SPoRC:::get_fishery_observation_model, obs_input(move_timing = 0))
  baranov_avail <- NAA[,,1:n_yrs,,,, drop = FALSE] * (1 - exp(-ZAA)) / ZAA

  expect_equal(dim(out$CAA), dim_fish_age)
  expect_equal(out$CAA, slice_catch(ret_FAA, baranov_avail), tolerance = 1e-12)
  expect_equal(out$DAA, slice_catch(disc_FAA, baranov_avail), tolerance = 1e-12)

  # the catch totals read the whole arrays fleet by fleet
  expect_equal(out$PredCatch, apply(out$CAA, c(1:4, 7), sum), tolerance = 1e-12)
  expect_equal(out$PredDiscard, apply(out$DAA, c(1:4, 7), sum), tolerance = 1e-12)

})

test_that("under move_timing = 2 the catch at age takes the season-integrated abundance, passed in or built", {

  # the season-integrated abundance over the model years, as the dynamics pass it
  NAA_int <- array(0, dim = dim_age)
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(seas in 1:n_seas) {
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            NAA_int[p,,y,seas,a,s] <- SPoRC:::integrate_seas_abundance(NAA[p,,y,seas,a,s], ZAA[p,,y,seas,a,s], Mrate[p,,,y,seas,a,s])
          } # end s loop
        } # end a loop
      } # end seas loop
    } # end y loop
  } # end p loop

  out_passed <- do.call(SPoRC:::get_fishery_observation_model, obs_input(move_timing = 2, NAA_int = NAA_int))
  out_built <- do.call(SPoRC:::get_fishery_observation_model, obs_input(move_timing = 2))

  expect_equal(dim(out_passed$CAA), dim_fish_age)
  expect_equal(out_passed$CAA, slice_catch(ret_FAA, NAA_int), tolerance = 1e-12)
  expect_equal(out_passed$DAA, slice_catch(disc_FAA, NAA_int), tolerance = 1e-12)
  expect_equal(out_built$CAA, out_passed$CAA, tolerance = 1e-12)
  expect_equal(out_built$DAA, out_passed$DAA, tolerance = 1e-12)
  expect_equal(out_built$PredFishIdx, out_passed$PredFishIdx, tolerance = 1e-12)

})
