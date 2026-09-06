library(SPoRC)
library(testthat)

# Do_Population_Projection is ordinary R, so a projection can be taped with RTMB
# and handed to an optimizer with an exact gradient. That only holds while every
# function in the call chain keeps AD values alive: the projection writes into
# preallocated numeric arrays, and without RTMB's replacement operator those
# arrays quietly coerce the AD values back to doubles. The failure is silent, so
# these tests check the gradient itself rather than that the call runs.
#
# Two options cannot be differentiated through and are refused rather than
# allowed to return a gradient that is wrong: inverse gaussian recruitment draws
# random numbers, and a catch target is inverted by a numerical solve.

# Beverton-Holt settings for the packaged sablefish model, which put projected
# spawning biomass back into recruitment and so put the whole feedback on the tape.
sable_srr_opt <- function() {
  d <- sgl_rg_sable_data
  rp <- sgl_rg_sable_rep
  ny <- length(d$years)
  n_pop <- d$n_pop; n_regions <- d$n_regions; n_seas <- d$n_seas
  n_ages <- length(d$ages); n_fish_fleets <- d$n_fish_fleets

  list(
    rec_dd = 1,
    rec_lag = 2,
    do_recruits_move = 0,
    R0 = rp$R0,
    h = array(rp$h_trans, dim = c(n_pop, n_regions)),
    rec_region_prop = rp$rec_region_prop,
    WAA = array(d$WAA[,,ny,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    MatAA = array(d$MatAA[,,ny,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    SSB = rp$SSB,
    Movement = array(1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
    sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)),
    sgl_seas_spawning_movement = NULL,
    stray_rate = array(0, dim = n_pop),
    # the packaged report predates seasonal M, so it has no season margin to slice
    natmort = array(rp$natmort[,,ny,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    fish_sel = array(rp$fish_sel[,,ny,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    ret_sel = array(rp$ret_sel[,,ny,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  )
}


test_that("a mean recruitment projection tapes and its gradient matches finite differences", {

  yield <- function(ln_F) sum(project_at_F(f = exp(ln_F), n_proj_yrs = 8)$proj_Catch)

  ln_F0 <- log(0.06)
  tape <- RTMB::MakeTape(yield, ln_F0)

  # a tape that has lost the AD class returns the value fine and the gradient wrong,
  # so both are checked
  expect_equal(tape(ln_F0), yield(ln_F0), tolerance = 1e-10)
  expect_equal(as.vector(tape$jacfun()(ln_F0)),
               numDeriv::grad(yield, ln_F0), tolerance = 1e-5)
  expect_gt(abs(as.vector(tape$jacfun()(ln_F0))), 1) # not a detached constant
})


test_that("a Beverton-Holt projection carries the stock-recruit feedback onto the tape", {

  srr_opt <- sable_srr_opt()
  yield <- function(ln_F) sum(project_at_F(f = exp(ln_F), n_proj_yrs = 12,
                                           recruitment_opt = "bh_rec",
                                           srr_opt = srr_opt)$proj_Catch)

  ln_F0 <- log(0.06)
  tape <- RTMB::MakeTape(yield, ln_F0)

  expect_equal(tape(ln_F0), yield(ln_F0), tolerance = 1e-10)
  expect_equal(as.vector(tape$jacfun()(ln_F0)),
               numDeriv::grad(yield, ln_F0), tolerance = 1e-5)
})


test_that("bind_proj_SSB keeps the AD class that abind drops", {

  hist <- array(1, dim = c(1, 1, 3))

  seen <- new.env()

  RTMB::MakeTape(function(x) {
    proj <- array(0, dim = c(1, 1, 2))
    proj[1,1,] <- c(x, x)
    seen$kept <- inherits(SPoRC:::bind_proj_SSB(hist, proj), "advector")
    seen$dropped <- inherits(abind::abind(hist, proj, along = 3), "advector")
    sum(x)
  }, 1)

  expect_true(seen$kept)
  expect_false(seen$dropped)
})


test_that("the projection refuses the options that cannot be differentiated through", {

  expect_error(
    RTMB::MakeTape(function(ln_F) sum(project_at_F(f = exp(ln_F), n_proj_yrs = 4,
                                                   recruitment_opt = "inv_gauss")$proj_Catch),
                   log(0.06)),
    "inv_gauss"
  )
})


test_that("adding the AD operators leaves an ordinary projection unchanged", {

  out <- project_at_F(f = 0.06, n_proj_yrs = 8)

  expect_false(inherits(out$proj_Catch, "advector"))
  expect_type(out$proj_Catch, "double")
  expect_true(all(is.finite(out$proj_SSB)))
})
