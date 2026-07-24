library(SPoRC)
library(testthat)

# do_*_mapping helpers call collect_message(); this function doesn't, but
# release_conv_tag_attr and friends live in the same package convention, so
# keep the same guard as the other test files for consistency/safety.
assign("messages_list", character(0), envir = .GlobalEnv)

# A minimal 1-pop, 1-region, 1-fleet, 1-cohort, 2-age, 2-year-of-liberty
# fixture. conv_fish_tag_attr = "p_a_s" (all dimensions attended) makes
# release_conv_tag_attr() a pure pass-through (see release_tag_attr.R: with
# every dim attended it just reshapes and returns the input unchanged), which
# isolates this function's own tag decay/movement/recapture mechanics from
# that separate, pre-existing apportionment logic. do_recruits_move = 0 and
# n_regions = 1 make movement a structural no-op (only ages >= 2 are moved,
# and with one region movement is multiplication by a 1x1 identity anyway).
make_tagging_input <- function() {

  n_pop <- 1; n_regions <- 1; n_fish_fleets <- 1; n_conv_tag_cohorts <- 1
  n_yrs <- 2; n_seas <- 1; n_ages <- 2; n_sexes <- 1; conv_tag_max_liberty <- 2

  Fmort <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  Fmort[1,1,1,1] <- 0.3; Fmort[1,2,1,1] <- 0.3

  natmort <- array(0.2, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes))

  list(
    n_fish_fleets = n_fish_fleets, n_regions = n_regions, n_conv_tag_cohorts = n_conv_tag_cohorts,
    n_yrs = n_yrs, n_seas = n_seas, n_pop = n_pop, n_ages = n_ages, n_sexes = n_sexes,
    conv_tag_fish_reporting_blocks = array(1L, dim = c(n_regions, n_yrs, n_fish_fleets)),
    conv_tag_fish_reporting_pars = array(0, dim = c(n_regions, 1, n_fish_fleets)), # logit(0.5) = 0
    conv_tag_fish_reporting = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)),
    conv_tag_release_indicator = matrix(c(1, 1, 1), nrow = 1, ncol = 3), # region 1, year 1, season 1
    conv_tag_max_liberty = conv_tag_max_liberty,
    use_conv_fish_tagging = c(1),
    Fmort = Fmort,
    fish_sel = array(1, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
    ret_sel = array(1, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)), # fully retained -> no discards
    dmr = array(1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)),
    natmort = natmort,
    seasdur = c(1),
    ln_conv_tag_shed = c(-10), # negligible shedding
    conv_tag_t_tagging = c(1), # tagging at the very start of the release season -> no Z discount
    conv_tagged_fish = array(c(500, 500), dim = c(n_conv_tag_cohorts, n_pop, n_ages, n_sexes)),
    conv_fish_tag_attr = c("p_a_s"),
    conv_tag_release_platform = matrix(0, nrow = n_conv_tag_cohorts, ncol = 1), # unused (p_a_s short-circuits release_conv_tag_attr)
    srv_sel = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, 1)), # unused, same reason
    NAA_bef = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)),   # unused, same reason
    ln_init_conv_tag_mort = c(-10), # negligible tag-induced mortality
    do_recruits_move = 0,
    Movement = array(1, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes)), # n_regions = 1 -> no-op
    conv_tag_fish_avail = array(0, dim = c(conv_tag_max_liberty + 1, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes)),
    pred_conv_tag_fish_recap = array(0, dim = c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets))
  )
}

test_that("get_tagging_observation_model converts reporting-rate parameters via inverse logit", {

  il <- make_tagging_input()
  out <- do.call(SPoRC:::get_tagging_observation_model, il)

  expect_equal(out$conv_tag_fish_reporting[1,1,1], plogis(0), tolerance = 1e-8) # 0.5
  expect_equal(out$conv_tag_fish_reporting[1,2,1], plogis(0), tolerance = 1e-8)
})

test_that("get_tagging_observation_model: release-year tag availability and predicted recapture (Baranov form)", {

  il <- make_tagging_input()
  out <- do.call(SPoRC:::get_tagging_observation_model, il)

  # tag mortality is negligible (ln_init_conv_tag_mort = -10), so available
  # tags at release are ~ the raw release numbers
  released <- 500 * exp(-exp(-10))
  expect_equal(out$conv_tag_fish_avail[1,1,1,1,1,1,1], released, tolerance = 1e-6)
  expect_equal(out$conv_tag_fish_avail[1,1,1,1,1,2,1], released, tolerance = 1e-6)

  # tmp_ZAA = natmort + F + shed (all ~uniform across ages in this fixture)
  Z <- 0.2 + 0.3 + exp(-10)
  SAA <- exp(-Z)
  reporting <- 0.5
  expected_recap <- reporting * (0.3 / Z) * released * (1 - SAA) # Baranov-style predicted recapture

  expect_equal(out$pred_conv_tag_fish_recap[1,1,1,1,1,1,1,1], expected_recap, tolerance = 1e-6)
  expect_equal(out$pred_conv_tag_fish_recap[1,1,1,1,1,2,1,1], expected_recap, tolerance = 1e-6)
})

test_that("get_tagging_observation_model: cross-year tag decay and plus-group accumulation", {

  il <- make_tagging_input()
  out <- do.call(SPoRC:::get_tagging_observation_model, il)

  released <- 500 * exp(-exp(-10))
  Z <- 0.2 + 0.3 + exp(-10)
  SAA <- exp(-Z)

  # n_ages = 2, so age 2 is the plus group: liberty-year 2's age-2 slot
  # accumulates liberty-year 1's age-1 (advanced) and age-2 (plus-group) survivors
  expected_ry2_age2 <- released * SAA + released * SAA
  expect_equal(out$conv_tag_fish_avail[2,1,1,1,1,2,1], expected_ry2_age2, tolerance = 1e-6)
})
