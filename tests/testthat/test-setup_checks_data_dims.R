library(SPoRC)
library(testthat)

# check_data_dimensions() validates every estimation input in the setup chain.
# The happy path is exercised throughout the suite, so these tests drive the
# error branches instead.
#
# Dimension sizes are all distinct. check_data_dimensions() compares dim(x)
# elementwise against an expected vector, so equal sizes would let a
# transposed array satisfy a check it should fail.

dims <- list(
  n_pop = 2,
  n_regions = 3,
  n_years = 4,
  n_seas = 5,
  n_ages = 6,
  n_lens = 7,
  n_sexes = 8,
  n_fish_fleets = 9,
  n_srv_fleets = 10,
  conv_tag_max_liberty = 12,
  n_conv_tag_cohorts = 13
)

# Observed composition ages are a free dimension the checks deliberately skip,
# so it must differ from n_ages for the skip to be meaningful.
n_obs_ages <- 11

arr <- function(d) array(0L, dim = d)

check_data <- function(x, what) {
  do.call(SPoRC:::check_data_dimensions, c(list(x = x, what = what), dims))
}

# Each case pairs a set of `what` labels with the dimensions they require.
array_cases <- list(
  # Biologicals
  list(what = c("WAA", "MatAA"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes))),
  list(what = "Fixed_natmort",
       d = with(dims, c(n_pop, n_regions, n_years, n_ages, n_sexes))),
  list(what = "WAA_fish",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets))),
  list(what = "WAA_srv",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets))),
  # The ageing-error checks validate only their leading dimensions; the trailing
  # observed-age axis is free, so tests that perturb the tail must skip these.
  list(what = "AgeingError", tail_checked = FALSE,
       d = with(dims, c(n_ages, n_obs_ages))),
  list(
    what = "AgeingError_t",
    tail_checked = FALSE,
    d = with(dims, c(n_years, n_ages, n_obs_ages))
  ),
  list(what = "SizeAgeTrans",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_lens, n_ages, n_sexes))),
  list(what = "Fixed_Movement",
       d = with(dims, c(n_pop, n_regions, n_regions, n_years, n_seas, n_ages, n_sexes))),
  list(what = "sgl_seas_spawning_movement",
       d = with(dims, c(n_pop, n_regions, n_regions, n_years, n_ages, n_sexes))),
  list(what = "stray_rate",
       d = with(dims, c(n_pop, n_years))),

  # Fishery
  list(what = c("ObsCatch", "UseCatch", "ObsDiscard", "UseDiscard",
                "ObsFishIdx", "ObsFishIdx_SE", "UseFishIdx",
                "UseFishAgeComps", "UseFishLenComps",
                "UseFishAgeComps_discard", "UseFishLenComps_discard"),
       d = with(dims, c(n_regions, n_years, n_seas, n_fish_fleets))),
  list(what = c("ObsCatch_pop", "UseCatch_pop", "ObsDiscard_pop", "UseDiscard_pop",
                "ObsFishIdx_pop", "ObsFishIdx_pop_SE", "UseFishIdx_pop",
                "UseFishAgeComps_pop", "UseFishLenComps_pop",
                "UseFishAgeComps_discard_pop", "UseFishLenComps_discard_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_fish_fleets))),
  list(what = c("ObsFishAgeComps", "ObsFishAgeComps_discard"),
       d = with(dims, c(n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_fish_fleets))),
  list(what = c("ObsFishLenComps", "ObsFishLenComps_discard"),
       d = with(dims, c(n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets))),
  list(what = c("ObsFishAgeComps_pop", "ObsFishAgeComps_discard_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_fish_fleets))),
  list(what = c("fish_sel_input_age", "ret_sel_input_age"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets))),
  list(what = c("ObsFishLenComps_pop", "ObsFishLenComps_discard_pop",
                "fish_sel_input_len", "ret_sel_input_len"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets))),
  list(what = c("ISS_FishLenComps", "ISS_FishAgeComps",
                "ISS_FishLenComps_discard", "ISS_FishAgeComps_discard"),
       d = with(dims, c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets))),
  list(what = c("ISS_FishLenComps_pop", "ISS_FishAgeComps_pop",
                "ISS_FishLenComps_discard_pop", "ISS_FishAgeComps_discard_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets))),

  # Survey
  list(what = c("ObsSrvIdx", "ObsSrvIdx_SE", "UseSrvIdx", "UseSrvAgeComps", "UseSrvLenComps"),
       d = with(dims, c(n_regions, n_years, n_seas, n_srv_fleets))),
  list(what = c("ObsSrvIdx_pop", "ObsSrvIdx_pop_SE", "UseSrvIdx_pop",
                "UseSrvAgeComps_pop", "UseSrvLenComps_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_srv_fleets))),
  list(what = "ObsSrvAgeComps",
       d = with(dims, c(n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_srv_fleets))),
  list(what = "ObsSrvLenComps",
       d = with(dims, c(n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets))),
  list(what = c("ISS_SrvLenComps", "ISS_SrvAgeComps"),
       d = with(dims, c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets))),
  list(what = c("ISS_SrvLenComps_pop", "ISS_SrvAgeComps_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets))),
  list(what = "ObsSrvAgeComps_pop",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_srv_fleets))),
  list(what = "srv_sel_input_age",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets))),
  list(what = c("ObsSrvLenComps_pop", "srv_sel_input_len"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets))),

  # Tagging
  list(what = "conv_tagged_fish",
       d = with(dims, c(n_conv_tag_cohorts, n_pop, n_ages, n_sexes))),
  list(what = "obs_conv_tag_fish_recap",
       d = with(dims, c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_pop,
                        n_regions, n_ages, n_sexes, n_fish_fleets)))
)

# `what` labels validated by length() rather than dim().
length_cases <- list(
  list(what = c("FishAgeComps_LikeType", "FishLenComps_LikeType",
                "pop_FishAgeComps_LikeType", "pop_FishLenComps_LikeType",
                "FishAgeComps_discard_LikeType", "FishLenComps_discard_LikeType",
                "pop_FishAgeComps_discard_LikeType", "pop_FishLenComps_discard_LikeType"),
       n = dims$n_fish_fleets),
  list(what = c("SrvAgeComps_LikeType", "SrvLenComps_LikeType",
                "pop_SrvAgeComps_LikeType", "pop_SrvLenComps_LikeType"),
       n = dims$n_srv_fleets)
)

test_that("correctly dimensioned inputs pass validation", {
  for(case in array_cases) {
    for(w in case$what) {
      expect_silent(check_data(arr(case$d), w))
    }
  } # end case loop
})

test_that("a wrong leading dimension is rejected", {
  for(case in array_cases) {
    bad <- case$d
    bad[1] <- bad[1] + 1
    for(w in case$what) {
      expect_error(check_data(arr(bad), w))
    }
  } # end case loop
})

test_that("a wrong trailing dimension is rejected", {
  # The trailing dimension is the fleet or sex axis for most inputs, and is the
  # one most easily got wrong when a fleet is added to a model.
  for(case in array_cases) {
    if(isFALSE(case$tail_checked)) next
    bad <- case$d
    n <- length(bad)
    bad[n] <- bad[n] + 1
    for(w in case$what) {
      expect_error(check_data(arr(bad), w))
    }
  } # end case loop
})

test_that("an array with too few dimensions is rejected", {
  for(case in array_cases) {
    if(length(case$d) < 3 || isFALSE(case$tail_checked)) next
    bad <- case$d[1:2]
    for(w in case$what) {
      expect_error(suppressWarnings(check_data(arr(bad), w)))
    }
  } # end case loop
})

test_that("a dimensionless vector is rejected where an array is required", {
  # dim() is NULL for a plain vector, so the comparison collapses to logical(0)
  # and the sum can never reach the required count.
  for(case in array_cases) {
    for(w in case$what) {
      expect_error(check_data(rep(0L, prod(case$d[1:2])), w))
    }
  } # end case loop
})

test_that("length-validated inputs accept the fleet count and reject others", {
  for(case in length_cases) {
    for(w in case$what) {
      expect_silent(check_data(rep(0, case$n), w))
      expect_error(check_data(rep(0, case$n + 1), w))
      expect_error(check_data(rep(0, case$n - 1), w))
    }
  } # end case loop
})

test_that("composition inputs ignore the observed-age dimension", {
  # ObsFishAgeComps and friends skip the observed-age axis so that comps can be
  # binned differently from the model age range.
  d_fish <- with(dims, c(n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_fish_fleets))
  for(n_obs in c(1, 3, 20)) {
    d_fish[4] <- n_obs
    expect_silent(check_data(arr(d_fish), "ObsFishAgeComps"))
  }

  d_pop <- with(dims, c(n_pop, n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_fish_fleets))
  for(n_obs in c(1, 3, 20)) {
    d_pop[5] <- n_obs
    expect_silent(check_data(arr(d_pop), "ObsFishAgeComps_pop"))
  }

  d_srv <- with(dims, c(n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_srv_fleets))
  for(n_obs in c(1, 3, 20)) {
    d_srv[4] <- n_obs
    expect_silent(check_data(arr(d_srv), "ObsSrvAgeComps"))
  }
})

test_that("an unrecognized `what` label validates nothing", {
  # Every branch is gated on a `what` match, so a typo silently skips
  # validation rather than erroring. Pinning this documents the behavior.
  expect_silent(check_data(arr(c(1, 1)), "ObsCatchh"))
  expect_silent(check_data(arr(c(1, 1)), "not_a_real_input"))
})

test_that("error messages name the offending input", {
  expect_error(check_data(arr(c(1, 1, 1, 1)), "ObsCatch"), "ObsCatch")
  expect_error(check_data(arr(c(1, 1, 1, 1)), "UseSrvIdx"), "UseSrvIdx")
  expect_error(check_data(arr(c(1, 1, 1, 1, 1, 1)), "WAA"), "WAA")
  expect_error(check_data(arr(c(1, 1, 1, 1)), "conv_tagged_fish"), "conv_tagged_fish")
})
