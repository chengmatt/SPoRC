library(SPoRC)
library(testthat)

# do_likelihood_profile() refuses to profile catchability when the model solves it
# analytically. Under fish_q_type/srv_q_type of "arith" or "geo" the model recomputes q from
# the index and never reads ln_fish_q/ln_srv_q, so the grid value never reaches the
# likelihood and the profile would come back flat.
#
# ln_fish_q and ln_srv_q are dimensioned [region, block, fleet] and idx has linear
# indices, so the fleet a profile touches only falls out of the third dim once the
# region and block strides are accounted for.

q_pars <- list(
  # 2 regions x 3 blocks x 2 fleets, so fleet 1 is linear index 1:6 and fleet 2 is 7:12
  ln_fish_q = array(0, dim = c(2, 3, 2)),
  ln_srv_q = array(0, dim = c(1, 1, 2))
)

guard <- function(data, what, idx) {
  tryCatch({
    SPoRC:::check_analytic_q(data, q_pars, what, idx)
    NA_character_
  }, error = function(e) conditionMessage(e))
}

test_that("profiling an analytic fishery catchability is refused", {
  msg <- guard(list(fish_q_type = c(0, 1)), "ln_fish_q", 7:12)
  expect_match(msg, "Cannot profile `ln_fish_q` for fleet(s) 2", fixed = TRUE)
  expect_match(msg, "`fish_q_type` is set to arith", fixed = TRUE)
})

test_that("the refusal names the analytic form and the survey parameter", {
  msg <- guard(list(srv_q_type = c(0, 2)), "ln_srv_q", 2)
  expect_match(msg, "Cannot profile `ln_srv_q` for fleet(s) 2", fixed = TRUE)
  expect_match(msg, "`srv_q_type` is set to geo", fixed = TRUE)
})

test_that("an estimated fleet still profiles when another fleet is analytic", {
  # Fleet 1 is "est" and occupies linear indices 1:6. Reading the fleet off the wrong
  # dim would pick up fleet 2 here and refuse a profile that is perfectly valid.
  expect_true(is.na(guard(list(fish_q_type = c(0, 1)), "ln_fish_q", 1:6)))
})

test_that("a NULL idx is treated as targeting every fleet", {
  expect_match(guard(list(fish_q_type = c(0, 1)), "ln_fish_q", NULL),
               "Cannot profile `ln_fish_q` for fleet(s) 2", fixed = TRUE)
})

test_that("estimated catchability everywhere passes the guard", {
  expect_true(is.na(guard(list(fish_q_type = c(0, 0)), "ln_fish_q", 1:12)))
})

test_that("a data list without q types passes the guard", {
  # Older data lists predate fish_q_type, and the model defaults them to estimated.
  expect_true(is.na(guard(list(), "ln_fish_q", 1:12)))
})

test_that("parameters other than catchability are left alone", {
  expect_true(is.na(guard(list(fish_q_type = c(1, 1)), "ln_global_R0", NULL)))
})
