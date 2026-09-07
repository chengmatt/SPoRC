# Dropping year one of a selectivity deviation series.
#
# A non-parametric form has one free base parameter per bin, so year one's deviation is that value
# written twice. Checks the map, the penalty it removes, and the refusal for the field forms.

library(SPoRC)
library(testthat)

sel_input <- function(...) {
  suppressMessages(suppressWarnings(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
    fishsel = c(list(
      fish_sel_model = "nonparfree_Fleet_1",
      cont_tv_fish_sel = "rw_Fleet_1",
      fishsel_pe_pars_spec = "est_all",
      fish_sel_devs_spec = "est_all",
      fish_sel_nonpar_est_bins = list(list(as.list(1:6)))
    ), list(...))
  )))
}

test_that("fishsel_dont_est_dev_first drops year one and leaves every later year alone", {

  keep <- sel_input()
  drop <- sel_input(fishsel_dont_est_dev_first = 1)

  map_keep <- array(as.numeric(keep$map$ln_fishsel_devs), dim = dim(keep$par$ln_fishsel_devs))
  map_drop <- array(as.numeric(drop$map$ln_fishsel_devs), dim = dim(drop$par$ln_fishsel_devs))

  expect_true(all(!is.na(map_keep[1, 1, , 1, 1])))
  expect_true(all(is.na(map_drop[1, 1, , 1, 1])))

  # one deviation per bin per year goes, and the rest of the series is untouched
  n_bins <- 6
  expect_equal(length(unique(map_keep[!is.na(map_keep)])) - length(unique(map_drop[!is.na(map_drop)])), n_bins)
  expect_true(all(!is.na(map_drop[1, 2:10, , 1, 1])))

  # the data mirror the penalty keys on follows the map
  expect_true(all(is.na(drop$data$map_ln_fishsel_devs[1, 1, , 1, 1])))
})

test_that("dropping year one removes exactly the walk's first-year term", {

  rep_of <- function(il) {
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    obj$report(obj$par)
  }

  keep <- sel_input()
  drop <- sel_input(fishsel_dont_est_dev_first = 1)

  n_bins <- 6
  init_sd <- 5 # fishsel_rw_init_sigma default

  # at zero deviations the first year contributes -log dnorm(0, 0, init_sd) per bin and nothing else
  gap <- sum(rep_of(keep)$sel_nLL) - sum(rep_of(drop)$sel_nLL)
  expect_equal(gap, n_bins * -dnorm(0, 0, init_sd, log = TRUE), tolerance = 1e-10)
})

test_that("the field forms refuse it", {
  expect_error(
    suppressMessages(suppressWarnings(sweep_input(
      dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
      fishsel = list(
        fish_sel_model = "nonparfree_Fleet_1",
        cont_tv_fish_sel = "2dar1_Fleet_1",
        fishsel_pe_pars_spec = "est_all",
        fish_sel_devs_spec = "est_all",
        fish_sel_nonpar_est_bins = list(list(as.list(1:6))),
        fishsel_dont_est_dev_first = 1
      )
    ))),
    "walk anchored at year one"
  )
})

test_that("the survey stage takes it per fleet", {

  il <- suppressMessages(suppressWarnings(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 2, n_yrs = 10, n_ages = 6),
    srvsel = list(
      srv_sel_model = c("nonparfree_Fleet_1", "nonparfree_Fleet_2"),
      cont_tv_srv_sel = c("rw_Fleet_1", "rw_Fleet_2"),
      srvsel_pe_pars_spec = rep("est_all", 2),
      srv_sel_devs_spec = rep("est_all", 2),
      srv_sel_nonpar_est_bins = list(list(as.list(1:6)), list(as.list(1:6))),
      srvsel_dont_est_dev_first = c(1, 0)
    )
  )))

  map_srv <- array(as.numeric(il$map$ln_srvsel_devs), dim = dim(il$par$ln_srvsel_devs))
  expect_true(all(is.na(map_srv[1, 1, , 1, 1])))  # fleet one drops year one
  expect_true(all(!is.na(map_srv[1, 1, , 1, 2]))) # fleet two keeps it
})

test_that("a bad value is refused", {
  expect_error(sel_input(fishsel_dont_est_dev_first = 2), "must be 0 or 1")
})
