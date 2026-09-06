library(SPoRC)
library(testthat)

# check_sim_dimensions() is the operating-model counterpart to
# check_data_dimensions(). Simulation inputs have a trailing n_sims axis, so
# these tests also confirm that axis is actually validated rather than ignored.
#
# Sizes are all distinct so a transposed array cannot satisfy a check.

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
  n_sims = 11
)

arr <- function(d) array(0L, dim = d)

check_sim <- function(x, what) {
  do.call(SPoRC:::check_sim_dimensions, c(list(x = x, what = what), dims))
}

array_cases <- list(
  # Biologicals
  list(what = c("WAA_input", "MatAA_input"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_sims))),
  list(what = "natmort_input",
       d = with(dims, c(n_pop, n_regions, n_years, n_ages, n_sexes, n_sims))),
  list(what = "WAA_fish_input",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets, n_sims))),
  list(what = "WAA_srv_input",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets, n_sims))),
  list(what = "SizeAgeTrans_input",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_lens, n_ages, n_sexes, n_sims))),

  # Fishing
  list(what = c("Fmort_input", "dmr_input"),
       d = with(dims, c(n_regions, n_years, n_seas, n_fish_fleets, n_sims))),
  list(what = "fish_q_input",
       d = with(dims, c(n_regions, n_years, n_fish_fleets, n_sims))),
  list(what = c("fish_sel_input", "ret_sel_input"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets, n_sims))),
  list(what = c("ln_sigmaC", "ObsFishIdx_SE", "ln_sigmaD"),
       d = with(dims, c(n_regions, n_years, n_seas, n_fish_fleets))),
  list(what = c("ln_sigmaC_pop", "ObsFishIdx_pop_SE", "ln_sigmaD_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_fish_fleets))),
  list(what = c("FishAge_pop_corr_pars_agg", "FishLen_pop_corr_pars_agg",
                "ln_FishAge_pop_theta_agg", "ln_FishLen_pop_theta_agg",
                "FishAge_discard_pop_corr_pars_agg", "FishLen_discard_pop_corr_pars_agg",
                "ln_FishAge_discard_pop_theta_agg", "ln_FishLen_discard_pop_theta_agg"),
       d = with(dims, c(n_pop, n_fish_fleets))),
  list(what = c("ISS_FishAgeComps", "ISS_FishLenComps",
                "ISS_FishAgeComps_discard", "ISS_FishLenComps_discard"),
       d = with(dims, c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims))),
  list(what = c("ISS_FishAgeComps_pop", "ISS_FishLenComps_pop",
                "ISS_FishAgeComps_discard_pop", "ISS_FishLenComps_discard_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets, n_sims))),
  list(what = c("ln_FishAge_theta", "ln_FishLen_theta",
                "ln_FishAge_discard_theta", "ln_FishLen_discard_theta"),
       d = with(dims, c(n_regions, n_sexes, n_fish_fleets))),
  list(what = c("ln_FishAge_pop_theta", "ln_FishLen_pop_theta",
                "ln_FishAge_discard_pop_theta", "ln_FishLen_discard_pop_theta"),
       d = with(dims, c(n_pop, n_regions, n_sexes, n_fish_fleets))),
  list(what = c("FishAge_corr_pars", "FishLen_corr_pars",
                "FishAge_discard_corr_pars", "FishLen_discard_corr_pars"),
       d = with(dims, c(n_regions, n_sexes, n_fish_fleets, 2))),
  list(what = c("FishAge_pop_corr_pars", "FishLen_pop_corr_pars",
                "FishAge_discard_pop_corr_pars", "FishLen_discard_pop_corr_pars"),
       d = with(dims, c(n_pop, n_regions, n_sexes, n_fish_fleets, 2))),
  list(what = c("FishAgeComps_Type", "FishLenComps_Type",
                "pop_FishAgeComps_Type", "pop_FishLenComps_Type",
                "FishAgeComps_discard_Type", "FishLenComps_discard_Type",
                "pop_FishAgeComps_discard_Type", "pop_FishLenComps_discard_Type"),
       d = with(dims, c(n_years, n_fish_fleets))),

  # Survey
  list(what = "srv_q_input",
       d = with(dims, c(n_regions, n_years, n_srv_fleets, n_sims))),
  list(what = "srv_sel_input",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets, n_sims))),
  list(what = "ObsSrvIdx_SE",
       d = with(dims, c(n_regions, n_years, n_seas, n_srv_fleets))),
  list(what = "ObsSrvIdx_pop_SE",
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_srv_fleets))),
  list(what = c("ISS_SrvAgeComps", "ISS_SrvLenComps"),
       d = with(dims, c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims))),
  list(what = c("ISS_SrvAgeComps_pop", "ISS_SrvLenComps_pop"),
       d = with(dims, c(n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets, n_sims))),
  list(what = c("ln_SrvAge_theta", "ln_SrvLen_theta"),
       d = with(dims, c(n_regions, n_sexes, n_srv_fleets))),
  list(what = c("ln_SrvAge_pop_theta", "ln_SrvLen_pop_theta"),
       d = with(dims, c(n_pop, n_regions, n_sexes, n_srv_fleets))),
  list(what = c("SrvAge_corr_pars", "SrvLen_corr_pars"),
       d = with(dims, c(n_regions, n_sexes, n_srv_fleets, 2))),
  list(what = c("SrvAge_pop_corr_pars_agg", "SrvLen_pop_corr_pars_agg",
                "ln_SrvAge_pop_theta_agg", "ln_SrvLen_pop_theta_agg"),
       d = with(dims, c(n_pop, n_srv_fleets))),
  list(what = c("SrvAgeComps_Type", "SrvLenComps_Type",
                "pop_SrvAgeComps_Type", "pop_SrvLenComps_Type"),
       d = with(dims, c(n_years, n_srv_fleets))),
  list(what = "t_srv",
       d = with(dims, c(n_regions, n_seas, n_srv_fleets))),
  list(what = c("SrvAge_pop_corr_pars", "SrvLen_pop_corr_pars"),
       d = with(dims, c(n_pop, n_regions, n_sexes, n_srv_fleets, 2))),

  # Recruitment
  list(what = "sexratio_input",
       d = with(dims, c(n_pop, n_regions, n_years, n_sexes, n_sims))),
  list(what = "stray_rate_input",
       d = with(dims, c(n_pop, n_years, n_sims))),
  list(what = "rec_seas_prop",
       d = with(dims, c(n_pop, n_seas, n_sims))),
  list(what = c("R0_input", "h_input"),
       d = with(dims, c(n_pop, n_regions, n_years, n_sims))),
  list(what = "rinit_input",
       d = with(dims, c(n_pop, n_regions, n_sims))),
  list(what = "ln_InitDevs_input",
       d = with(dims, c(n_pop, n_regions, n_ages - 1, n_sims))),

  # Tagging
  list(what = "conv_tag_fish_reporting_input",
       d = with(dims, c(n_regions, n_years, n_fish_fleets, n_sims)))
)

length_cases <- list(
  list(what = c("comp_fishage_like", "ln_FishAge_theta_agg", "FishAge_corr_pars_agg",
                "comp_fishlen_like", "ln_FishLen_theta_agg", "FishLen_corr_pars_agg",
                "pop_comp_fishage_like", "pop_comp_fishlen_like",
                "comp_fishage_discard_like", "ln_FishAge_discard_theta_agg",
                "FishAge_discard_corr_pars_agg", "comp_fishlen_discard_like",
                "ln_FishLen_discard_theta_agg", "FishLen_discard_corr_pars_agg",
                "pop_comp_fishage_discard_like", "pop_comp_fishlen_discard_like"),
       n = dims$n_fish_fleets),
  list(what = c("comp_srvage_like", "ln_SrvAge_theta_agg", "SrvAge_corr_pars_agg",
                "comp_srvlen_like", "ln_SrvLen_theta_agg", "SrvLen_corr_pars_agg",
                "pop_comp_srvage_like", "pop_comp_srvlen_like"),
       n = dims$n_srv_fleets)
)

test_that("correctly dimensioned simulation inputs pass validation", {
  for(case in array_cases) {
    for(w in case$what) {
      expect_silent(check_sim(arr(case$d), w))
    }
  } # end case loop
})

test_that("a wrong leading dimension is rejected", {
  for(case in array_cases) {
    bad <- case$d
    bad[1] <- bad[1] + 1
    for(w in case$what) {
      expect_error(check_sim(arr(bad), w))
    }
  } # end case loop
})

test_that("a wrong trailing dimension is rejected", {
  # For most simulation inputs this is the n_sims axis, so this also confirms
  # the simulation dimension is checked.
  for(case in array_cases) {
    bad <- case$d
    n <- length(bad)
    bad[n] <- bad[n] + 1
    for(w in case$what) {
      expect_error(check_sim(arr(bad), w))
    }
  } # end case loop
})

test_that("an array with too few dimensions is rejected", {
  for(case in array_cases) {
    if(length(case$d) < 3) next
    bad <- case$d[1:2]
    for(w in case$what) {
      expect_error(suppressWarnings(check_sim(arr(bad), w)))
    }
  } # end case loop
})

test_that("a dimensionless vector is rejected where an array is required", {
  for(case in array_cases) {
    for(w in case$what) {
      expect_error(check_sim(rep(0L, prod(case$d[1:2])), w))
    }
  } # end case loop
})

test_that("length-validated simulation inputs accept only the fleet count", {
  for(case in length_cases) {
    for(w in case$what) {
      expect_silent(check_sim(rep(0, case$n), w))
      expect_error(check_sim(rep(0, case$n + 1), w))
      expect_error(check_sim(rep(0, case$n - 1), w))
    }
  } # end case loop
})

test_that("correlation-parameter inputs require a trailing axis of exactly 2", {
  # The final axis holds the two AR1 correlation parameters and is hard-coded,
  # so it must not track any of the model dimensions.
  d <- with(dims, c(n_regions, n_sexes, n_fish_fleets, 2))
  expect_silent(check_sim(arr(d), "FishAge_corr_pars"))
  d[4] <- 3
  expect_error(check_sim(arr(d), "FishAge_corr_pars"))

  d <- with(dims, c(n_pop, n_regions, n_sexes, n_srv_fleets, 2))
  expect_silent(check_sim(arr(d), "SrvAge_pop_corr_pars"))
  d[5] <- 1
  expect_error(check_sim(arr(d), "SrvAge_pop_corr_pars"))
})

test_that("an unrecognized `what` label validates nothing", {
  expect_silent(check_sim(arr(c(1, 1)), "Fmort_inputt"))
  expect_silent(check_sim(arr(c(1, 1)), "not_a_real_input"))
})

test_that("error messages name the offending input", {
  expect_error(check_sim(arr(c(1, 1, 1, 1, 1)), "Fmort_input"), "Fmort_input")
  expect_error(check_sim(arr(c(1, 1, 1, 1)), "srv_q_input"), "srv_q_input")
  expect_error(check_sim(arr(c(1, 1, 1)), "rec_seas_prop"), "rec_seas_prop")
})
