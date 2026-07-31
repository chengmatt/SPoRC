library(SPoRC)
library(testthat)

# ── dmr_dev_spec / Use_dmr_pen consistency checks ────────────────────────────

# Minimal input_list accepted by Setup_Mod_Catch_and_F: only the dimensions and
# the catch arrays are needed to reach the validation block and the mapping calls.
make_catch_input_list <- function(n_pop = 1, n_regions = 1, n_yrs = 3, n_seas = 1, n_fish_fleets = 1) {
  list(
    data = list(n_pop = n_pop, n_regions = n_regions, years = 1:n_yrs,
                n_seas = n_seas, n_fish_fleets = n_fish_fleets),
    par = list(),
    map = list(),
    verbose = FALSE,
    store_config = FALSE
  )
}

# Setup_Mod_Catch_and_F emits several unrelated "specified as fix, but no
# starting values" warnings, so collect every warning rather than relying on
# the ordering expect_warning() sees.
collect_warnings <- function(expr) {
  warns <- character(0)
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  warns
}

run_setup <- function(dmr_dev_spec, Use_dmr_pen, n_yrs = 3) {
  il <- make_catch_input_list(n_yrs = n_yrs)
  dims <- c(il$data$n_regions, n_yrs, il$data$n_seas, il$data$n_fish_fleets)
  Setup_Mod_Catch_and_F(
    input_list = il,
    ObsCatch = array(100, dim = dims),
    UseCatch = array(1, dim = dims),
    Use_dmr_pen = Use_dmr_pen,
    dmr_dev_spec = dmr_dev_spec
  )
}

dmr_pen_msg <- "dmr_dev_spec is 'est_all' but Use_dmr_pen is 0"

test_that("estimating dmr deviations with no penalty warns", {
  warns <- collect_warnings(run_setup(dmr_dev_spec = "est_all", Use_dmr_pen = 0))
  expect_true(any(grepl(dmr_pen_msg, warns, fixed = TRUE)))
})

test_that("estimating dmr deviations with a penalty does not warn", {
  warns <- collect_warnings(run_setup(dmr_dev_spec = "est_all", Use_dmr_pen = 1))
  expect_false(any(grepl(dmr_pen_msg, warns, fixed = TRUE)))
})

test_that("fixed dmr deviations with no penalty does not warn", {
  warns <- collect_warnings(run_setup(dmr_dev_spec = "fix", Use_dmr_pen = 0))
  expect_false(any(grepl(dmr_pen_msg, warns, fixed = TRUE)))
})

test_that("fixed dmr deviations with a penalty is an error", {
  expect_error(suppressWarnings(run_setup(dmr_dev_spec = "fix", Use_dmr_pen = 1)),
               "Cannot apply a penalty on deviations that are not estimated")
})

test_that("'est' is not an accepted dmr_dev_spec", {
  # The estimating spec is named "est_all"; guarding on "est" would be dead code.
  expect_error(suppressWarnings(run_setup(dmr_dev_spec = "est", Use_dmr_pen = 0)),
               "dmr_dev_spec 'est' not recognized")
})
