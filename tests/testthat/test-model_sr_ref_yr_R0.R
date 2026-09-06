library(SPoRC)
library(testthat)

# SR_ref_yr fixes the BIOLOGY that goes into spawning biomass per recruit. The R0 that
# turns that per-recruit quantity into S0 is the year's own value, and it is the same R0
# the curve's numerator uses. The two cannot be drawn from different years: steepness is
# defined as recruitment at 0.2 * S0 as a fraction of R0, so splitting them leaves the
# unfished state off the curve.

bh_rec <- function(SSB, h, R0, S0) 4 * h * R0 * SSB / ((1 - h) * S0 + (5 * h - 1) * SSB)

test_that("unfished is a fixed point only when S0 and the numerator share an R0", {

  phi0 <- 2.5; h <- 0.7
  for(R0 in c(5, 12)) {
    S0 <- R0 * phi0
    expect_equal(bh_rec(S0, h, R0, S0), R0, tolerance = 1e-10)   # matched
  } # end R0 loop

  # numerator on block 2, S0 still built from block 1: the unfished biomass block 2
  # implies no longer returns block 2's R0
  R0_1 <- 5; R0_2 <- 12
  mismatched <- bh_rec(R0_2 * phi0, h, R0_2, R0_1 * phi0)
  expect_false(isTRUE(all.equal(mismatched, R0_2)))
  expect_equal(mismatched, 12.8, tolerance = 1e-8)
})

test_that("steepness keeps its meaning under the matched pairing and loses it otherwise", {

  phi0 <- 2.5; h <- 0.7; R0 <- 12
  S0 <- R0 * phi0
  # h is recruitment at 0.2 * S0 as a fraction of R0
  expect_equal(bh_rec(0.2 * S0, h, R0, S0) / R0, h, tolerance = 1e-10)

  # the same curve with S0 taken from a different block no longer honors that definition
  S0_wrong <- 5 * phi0
  expect_false(isTRUE(all.equal(bh_rec(0.2 * S0, h, R0, S0_wrong) / R0, h, tolerance = 1e-6)))
})

test_that("SR_ref_yr moves the per-recruit biology without touching R0", {

  skip_on_cran()
  # both settings are accepted and recorded; the R0 side is untouched by either
  d <- Setup_Mod_Dim(
    years = 1:30,
    ages = 1:10,
    lens = 1,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    verbose = FALSE
  )
  mk <- function(yr) suppressWarnings(suppressMessages(
    Setup_Mod_Rec(
      d,
      do_rec_bias_ramp = 0,
      sigmaR_switch = 1,
      ln_sigmaR = array(log(1), c(2, 1, 1)),
      rec_model = "mean_rec",
      sigmaR_spec = "fix",
      init_age_strc = 1,
      equil_init_age_strc = 2,
      SR_ref_yr = yr
    )))
  expect_equal(mk(1)$data$SR_ref_yr, 1)
  expect_equal(mk(30)$data$SR_ref_yr, 30)
  expect_error(suppressMessages(mk(31)), "SR_ref_yr must be a single year index")
  # R0 is sized by its own blocks, not by SR_ref_yr
  expect_equal(dim(mk(30)$par$ln_global_R0), dim(mk(1)$par$ln_global_R0))
})
