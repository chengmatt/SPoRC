library(SPoRC)
library(testthat)

# Retention is the only selectivity stream whose process error weight and random walk
# initial sigma were read by the objective without being reachable from setup. These
# pin that they arrive through Setup_Mod_Fishsel_and_Q, which is where a user sets
# retention, and that a zero weight actually removes the retention process error.

build <- function(...) suppressWarnings(suppressMessages(objective_fixture_input(...)))

# A time-varying retention curve, so the process error block the weight scales is
# reached at all. Without this the weight is inert whatever it is set to.
tv_retsel <- list(
  cont_tv_ret_sel = "iid_Fleet_1",
  retsel_pe_pars_spec = "fix",
  ret_sel_devs_spec = "est_shared_r"
)


test_that("retention process error controls default and reach the data list", {
  input <- build()
  n_fish_fleets <- input$data$n_fish_fleets

  expect_equal(input$data$retsel_pe_wt, rep(1, n_fish_fleets))
  expect_equal(input$data$retsel_rw_init_sigma, rep(5, n_fish_fleets))

  # set through Setup_Mod_Fishsel_and_Q, which forwards to Setup_Mod_Retsel
  user <- build(fishsel = list(retsel_pe_wt = 0, retsel_rw_init_sigma = NA))
  expect_equal(user$data$retsel_pe_wt, 0)
  expect_equal(user$data$retsel_rw_init_sigma, NA)
})


test_that("a wrong length process error weight errors rather than becoming NA per fleet", {
  # The message names the length rather than the value: a per-fleet setting given
  # at the wrong length used to be reported as an unrecognized value.
  expect_error(build(fishsel = list(retsel_pe_wt = c(1, 1))),
               "retsel_pe_wt has 2 entries for 1 fleet")
  expect_error(build(fishsel = list(retsel_rw_init_sigma = c(5, 5))),
               "retsel_rw_init_sigma has 2 entries for 1 fleet")

  # the same check now guards the fishery and survey streams it was missing from
  expect_error(build(fishsel = list(fishsel_pe_wt = c(1, 1))),
               "fishsel_pe_wt has 2 entries for 1 fleet")
  expect_error(build(srvsel = list(srvsel_pe_wt = c(1, 1))),
               "srvsel_pe_wt has 2 entries for 1 fleet")
})


test_that("a zero weight removes the retention process error from the objective", {
  weighted <- evaluate_input(build(fishsel = tv_retsel))
  unweighted <- evaluate_input(build(fishsel = modifyList(tv_retsel, list(retsel_pe_wt = 0))))

  # same parameters either way, so any jnLL difference is the dropped penalty alone
  expect_equal(weighted$obj$par, unweighted$obj$par)

  expect_equal(dim(weighted$rep$sel_nLL), dim(unweighted$rep$sel_nLL))
  expect_false(isTRUE(all.equal(sum(weighted$rep$sel_nLL), sum(unweighted$rep$sel_nLL))))
  expect_equal(weighted$rep$jnLL - unweighted$rep$jnLL,
               sum(weighted$rep$sel_nLL) - sum(unweighted$rep$sel_nLL),
               tolerance = 1e-10)

  # both still decompose into their reported terms
  expect_jnLL_decomposes(weighted)
  expect_jnLL_decomposes(unweighted)
})
