library(SPoRC)
library(testthat)

# A minimal 1-pop, 1-region, 2-season, 3-age fixture. n_regions = 1 means the
# movement step is a structural no-op (gated by n_regions > 1), so this
# isolates the mortality/ageing/recruitment mechanics without also depending
# on Get_Movement()'s correctness. rec_model = 0 (mean recruitment) is used
# because it's a fully deterministic, trivial formula
# (R0 * rec_region_prop), so the recruitment insertion step can be checked
# exactly too, without depending on the Beverton-Holt/SSB machinery.
# spawn_seas = 1 and rec_seas_prop = c(1, 0) mean season 2 gets none of the
# year's recruitment, so the season-2 "insert seasonal recruits" step is a
# guaranteed no-op and doesn't need separate verification.
make_pop_proj_input <- function() {

  n_pop <- 1; n_regions <- 1; n_seas <- 2; n_ages <- 3; n_sexes <- 1; n_yrs <- 1; n_fish_fleets <- 1

  list(
    n_pop = n_pop, n_regions = n_regions, n_seas = n_seas, n_ages = n_ages, n_sexes = n_sexes,
    n_yrs = n_yrs, n_fish_fleets = n_fish_fleets, n_est_rec_devs = 0,
    rec_lag = 0, rec_model = 0, rec_dd = 0,
    R0 = c(2000),
    rec_region_prop = matrix(1, 1, 1),
    rec_seas_prop = matrix(c(1, 0), nrow = 1, ncol = 2),
    h_trans = matrix(0.7, 1, 1), # unused by rec_model = 0
    natal_region = c(1), t_spawn = 0.5, spawn_seas = 1, seasdur = c(0.5, 0.5),
    init_F = array(0, dim = 1), # unused by rec_model = 0
    ln_RecDevs = array(0, dim = c(n_pop, n_regions, 1)), # never indexed (n_est_rec_devs = 0)
    sexratio = array(1, dim = c(n_pop, n_regions, n_yrs, n_sexes)),
    WAA = array(2, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)),
    MatAA = array(0.5, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)),
    natmort = array(0.2, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes)),
    Movement = array(1, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes)),
    stray_rate = array(0, dim = c(n_pop, n_yrs)),
    sgl_seas_spawning_movement = array(1, dim = c(n_pop, n_regions, n_regions, n_yrs, n_ages, n_sexes)),
    do_recruits_move = 0,
    fish_sel = array(1, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
    ret_sel = array(1, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
    dmr = array(1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)),
    # ZAA differs by season so the two mortality steps are independently checkable
    ZAA = local({
      arr <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
      arr[1,1,1,1,,1] <- c(0.5, 0.3, 0.4)   # season 1
      arr[1,1,1,2,,1] <- c(0.6, 0.35, 0.45) # season 2
      arr
    }),
    NAA = local({
      arr <- array(0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes))
      arr[1,1,1,1,,1] <- c(999, 40, 20) # age 1 placeholder (overwritten by recruitment); ages 2-3 are the known check values
      arr
    }),
    NAA0 = array(0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)),
    NAA_bef = array(0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)),
    NAA_aft = array(0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)),
    Rec = array(0, dim = c(n_pop, n_regions, n_yrs)),
    SSB = array(0, dim = c(n_pop, n_regions, n_yrs)),
    Total_Biom = array(0, dim = c(n_pop, n_regions, n_yrs)),
    Dynamic_SSB0 = array(0, dim = c(n_pop, n_regions, n_yrs)),
    eff_SSB = array(0, dim = c(n_pop, n_yrs))
  )
}

test_that("get_population_projection inserts mean recruitment (rec_model = 0) at age 1", {

  il <- make_pop_proj_input()
  out <- do.call(SPoRC:::get_population_projection, il)

  # rec_model = 0: rec = R0 * rec_region_prop; spawn_seas = 1 gets all of rec_seas_prop
  expected_rec <- 2000 * 1
  expect_equal(out$NAA[1,1,1,1,1,1], expected_rec, tolerance = 1e-8)
  expect_equal(out$Rec[1,1,1], expected_rec, tolerance = 1e-8)
})

test_that("get_population_projection: within-year seasonal mortality matches exp(-Z) decay", {

  il <- make_pop_proj_input()
  out <- do.call(SPoRC:::get_population_projection, il)

  # Ages 2-3 at season 1 are untouched by recruitment (which only ever writes
  # to age 1), so their season 1 -> 2 transition is pure exponential decay.
  expect_equal(out$NAA[1,1,1,2,2,1], 40 * exp(-0.3), tolerance = 1e-8)
  expect_equal(out$NAA[1,1,1,2,3,1], 20 * exp(-0.4), tolerance = 1e-8)

  # Age 1 is also checkable here since rec_model = 0 makes recruitment exact.
  expect_equal(out$NAA[1,1,1,2,1,1], 2000 * exp(-0.5), tolerance = 1e-8)
})

test_that("get_population_projection: cross-year age advancement and plus-group accumulation", {

  il <- make_pop_proj_input()
  out <- do.call(SPoRC:::get_population_projection, il)

  naa_s2_age2 <- 40 * exp(-0.3) # verified independently above
  naa_s2_age3 <- 20 * exp(-0.4)

  # age 2 -> age 3 advance, plus age 3 -> age 3 plus-group accumulation, both
  # using season 2's Z (0.35 for age 2, 0.45 for age 3)
  expected_plus_group <- naa_s2_age2 * exp(-0.35) + naa_s2_age3 * exp(-0.45)
  expect_equal(out$NAA[1,1,2,1,3,1], expected_plus_group, tolerance = 1e-8)
})

test_that("get_population_projection returns finite, non-negative derived quantities", {

  il <- make_pop_proj_input()
  out <- do.call(SPoRC:::get_population_projection, il)

  expect_true(all(is.finite(out$SSB)))
  expect_true(all(out$SSB >= 0))
  expect_true(all(is.finite(out$Total_Biom)))
  expect_true(all(out$Total_Biom >= 0))
  expect_true(all(is.finite(out$Aggregated_SSB)))
  expect_true(all(is.finite(out$Dynamic_Aggregated_SSB0)))
})
