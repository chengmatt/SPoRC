library(SPoRC)
library(testthat)

# The composition likelihood message used to list six of the nine accepted values and never said
# which fleet was wrong. It now names the fleet, its value and the nearest accepted one.

msg <- function(...) tryCatch({ SPoRC:::check_comp_like_type(...); NA_character_ },
                              error = function(e) conditionMessage(e))


test_that("the message names the fleet and every accepted value", {

  m <- msg(c("Multinomial", "multinomial", "none"), "FishAgeComps_LikeType")

  expect_match(m, "FishAgeComps_LikeType has 1 entry")
  expect_match(m, 'fleet 2 was given "multinomial"')
  # the three -miss0 forms are what the old message left out
  for(v in SPoRC:::comp_like_type_options()) expect_match(m, v, fixed = TRUE)
})


test_that("a near miss gets a suggestion and a distant one does not", {

  expect_match(msg("Multinomal", "FishAgeComps_LikeType"), 'did you mean "Multinomial"')
  expect_match(msg("2d-Logistic-Norm", "FishAgeComps_LikeType"), 'did you mean "2d-Logistic-Normal"')
  expect_match(msg("Dirichlet Multinomial", "FishAgeComps_LikeType"), 'did you mean "Dirichlet-Multinomial"')

  # an abbreviation is not close enough to any one option to guess at
  expect_false(grepl("did you mean", msg("DM", "FishAgeComps_LikeType")))
  expect_false(grepl("did you mean", msg(NA_character_, "FishAgeComps_LikeType")))
})


test_that("conditional age-at-length says why its list is shorter", {

  m <- msg("iid-Logistic-Normal", "CAAL_LikeType",
           allowed = c("none", "Multinomial", "Dirichlet-Multinomial"),
           note = "The logistic-normal families are not available for conditional age-at-length.")

  expect_match(m, "not available for conditional age-at-length", fixed = TRUE)
  expect_false(grepl("miss0", m))
})


test_that("a valid setting passes through untouched", {
  x <- c("none", "Multinomial", "2d-Logistic-Normal-miss0")
  expect_identical(SPoRC:::check_comp_like_type(x, "FishAgeComps_LikeType"), x)
  expect_true(is.na(msg(x, "FishAgeComps_LikeType")))
})


test_that("the exported setup functions raise it", {

  n <- sweep_dims$n_fish_fleets
  m <- tryCatch(sweep_input(fishidx = list(FishAgeComps_LikeType = c("Multinomal", rep("none", n - 1))),
                            stop_after = "fishidx"),
                error = function(e) conditionMessage(e))

  expect_match(m, "FishAgeComps_LikeType has 1 entry")
  expect_match(m, 'did you mean "Multinomial"')
})
