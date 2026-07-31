library(SPoRC)
library(testthat)

# ── mapping a deviation off by hand removes its penalty ─────────────────────
#
# The map is applied by RTMB::MakeADFun and is invisible inside the objective,
# so the F and dmr deviation penalties read a mirror of it carried in the data
# list. fit_model() refreshes those mirrors from the map it is about to hand to
# MakeADFun, which is what makes a hand edit to the map take effect on the
# penalty as well as on what is estimated.

build <- function(...) suppressWarnings(suppressMessages(objective_fixture_input(...)))

fit_no_optim <- function(input) {
  fit_model(data = input$data, parameters = input$par, mapping = input$map,
            random = NULL, do_optim = FALSE, silent = TRUE)
}

map_off_years <- function(input, par_nm, yrs) {
  mp <- array(input$map[[par_nm]], dim = dim(input$par[[par_nm]]))
  mp[, yrs, , ] <- NA
  input$map[[par_nm]] <- factor(mp)
  input
}

test_that("dmr deviations mapped off by hand are not penalized", {
  input <- build(catch_f = list(Use_dmr_pen = 1, dmr_dev_spec = "est_all"))
  n_yrs <- dim(input$par$logit_dmr_devs)[2]
  off_yrs <- 10:20

  before <- fit_no_optim(input)
  expect_true(all(before$rep$dmr_nLL != 0)) # every year penalized to begin with

  after <- fit_no_optim(map_off_years(input, "logit_dmr_devs", off_yrs))

  expect_true(all(after$rep$dmr_nLL[, off_yrs, , ] == 0))
  expect_true(all(after$rep$dmr_nLL[, setdiff(1:n_yrs, off_yrs), , ] != 0))

  # the penalized set is exactly the estimated set
  expect_equal(as.vector(after$rep$dmr_nLL != 0),
               !is.na(after$mapping$logit_dmr_devs))
})

test_that("F deviations mapped off by hand are not penalized", {
  input <- build()
  n_yrs <- dim(input$par$ln_F_devs)[2]
  off_yrs <- 4:7

  before <- fit_no_optim(input)
  expect_true(all(before$rep$Fmort_nLL != 0))

  after <- fit_no_optim(map_off_years(input, "ln_F_devs", off_yrs))

  expect_true(all(after$rep$Fmort_nLL[, off_yrs, , ] == 0))
  expect_true(all(after$rep$Fmort_nLL[, setdiff(1:n_yrs, off_yrs), , ] != 0))
})

test_that("a dmr deviation mapped off falls back on the mean at the default starting value", {
  input <- build(catch_f = list(Use_dmr_pen = 1, dmr_dev_spec = "est_all"))
  expect_true(all(input$par$logit_dmr_devs == 0)) # the default Setup_Mod_Catch_and_F start

  fit <- fit_no_optim(map_off_years(input, "logit_dmr_devs", 10:20))
  logit_mean <- as.numeric(fit$env$last.par.best[names(fit$env$last.par.best) == "logit_dmr_mean"])[1]

  expect_equal(fit$rep$dmr[1, 12, 1, 1], stats::plogis(logit_mean), tolerance = 1e-10)
})

test_that("a dmr deviation mapped off stays at a non-zero starting value rather than reverting", {
  # mapping a parameter off pins it at whatever is in $par, so a user who supplies
  # starting deviations and then maps some off freezes them at those values. This
  # is deliberate -- supplied starting values are treated as intentional -- but it
  # is not the same thing as reverting to logit_dmr_mean
  input <- build(catch_f = list(Use_dmr_pen = 1, dmr_dev_spec = "est_all"))
  input$par$logit_dmr_devs[] <- 0.4

  fit <- fit_no_optim(map_off_years(input, "logit_dmr_devs", 10:20))
  logit_mean <- as.numeric(fit$env$last.par.best[names(fit$env$last.par.best) == "logit_dmr_mean"])[1]

  expect_equal(fit$rep$dmr[1, 12, 1, 1], stats::plogis(logit_mean + 0.4), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fit$rep$dmr[1, 12, 1, 1], stats::plogis(logit_mean))))
})

test_that("an input list carrying no map mirrors falls back on the fished set", {
  # maintain_backwards_compatibility() rebuilds the set the penalties used to
  # compute for themselves, so input lists saved by older SPoRC versions still
  # evaluate to the same penalties
  input <- build(catch_f = list(Use_dmr_pen = 1, dmr_dev_spec = "est_all"))
  current <- fit_no_optim(input)

  legacy <- input
  legacy$data$map_ln_F_devs <- NULL
  legacy$data$map_logit_dmr_devs <- NULL
  expect_false(any(c("map_ln_F_devs", "map_logit_dmr_devs") %in% names(legacy$data)))

  fit <- fit_no_optim(legacy)

  expect_equal(fit$rep$dmr_nLL, current$rep$dmr_nLL, tolerance = 1e-10)
  expect_equal(fit$rep$Fmort_nLL, current$rep$Fmort_nLL, tolerance = 1e-10)
  expect_equal(fit$rep$jnLL, current$rep$jnLL, tolerance = 1e-10)
})

# ── sync_dev_map_data ───────────────────────────────────────────────────────

test_that("sync_dev_map_data refreshes mirrors from the map and leaves others alone", {
  input <- build(catch_f = list(Use_dmr_pen = 1, dmr_dev_spec = "est_all"))

  mp <- array(input$map$logit_dmr_devs, dim = dim(input$par$logit_dmr_devs))
  mp[, 2, , ] <- NA
  input$map$logit_dmr_devs <- factor(mp)

  synced <- SPoRC:::sync_dev_map_data(input$data, input$map)
  expect_true(all(is.na(synced$map_logit_dmr_devs[, 2, , ])))
  expect_false(any(is.na(synced$map_logit_dmr_devs[, 1, , ])))

  # a mirror whose map has been truncated elsewhere is left untouched rather
  # than written at a mismatched length
  truncated <- input$map
  truncated$logit_dmr_devs <- factor(mp[, 1:5, , , drop = FALSE])
  untouched <- SPoRC:::sync_dev_map_data(input$data, truncated)
  expect_equal(dim(untouched$map_logit_dmr_devs), dim(input$data$map_logit_dmr_devs))

  # a parameter absent from the map is skipped
  expect_equal(SPoRC:::sync_dev_map_data(input$data, list())$map_logit_dmr_devs,
               input$data$map_logit_dmr_devs)
})
