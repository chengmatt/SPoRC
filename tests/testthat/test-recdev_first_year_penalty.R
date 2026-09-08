library(SPoRC)
library(testthat)

# The first year's recruitment is the first year's age one abundance, which under an
# equilibrium initialization belongs to the initial condition. dont_pen_recdev_first
# leaves those years estimated while taking their penalty away, which mapping them off
# would not do: that fixes them instead.

rec_input <- function(dont_pen_recdev_first = 0) {
  sim <- seasonal_M_sim()
  input_list <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))
  SPoRC:::do_RecDevs_mapping(input_list, RecDevs_spec = NULL, rec_dd = 0,
                             dont_pen_recdev_first = dont_pen_recdev_first)
}


test_that("the dropped years leave the map alone and only the penalty mirror", {

  base <- rec_input(0)
  dropped <- rec_input(2)

  # still estimated, so the parameter map is untouched
  expect_identical(base$map$ln_RecDevs, dropped$map$ln_RecDevs)

  # the mirror the penalty reads loses the first two years and keeps the rest
  expect_true(all(is.na(dropped$data$map_ln_RecDevs[, , 1:2])))
  expect_equal(dropped$data$map_ln_RecDevs[, , -(1:2)], base$data$map_ln_RecDevs[, , -(1:2)])
  expect_false(any(is.na(base$data$map_ln_RecDevs)))
})


test_that("the penalty drops those years and nothing else", {

  n_yrs <- seasonal_M_cfg$n_yrs
  devs <- array(stats::rnorm(n_yrs, 0, 0.4), dim = c(1, 1, n_yrs))

  penalty <- function(map_mirror) {
    SPoRC:::get_recruitment_penalty(
      n_pop = 1, n_regions = 1, n_ages = seasonal_M_cfg$n_ages, n_est_rec_devs = n_yrs,
      rec_region_prop_spec = 0,
      rec_region_prop = matrix(1, 1, 1), equil_init_age_strc = 0,
      ln_InitDevs = array(0, dim = c(1, 1, seasonal_M_cfg$n_ages - 1, 1)),
      init_age_devs_shared = NA, ln_sigmaR = array(log(0.5), dim = c(2, 1, 1)),
      bias_ramp = rep(0, n_yrs), sigmaR_switch = 1, ln_RecDevs = devs,
      sigmaR2_early = matrix(0.25, 1, 1), sigmaR2_late = matrix(0.25, 1, 1),
      do_rec_bias_ramp = 0, map_ln_RecDevs = map_mirror)$Rec_nLL
  }

  full <- array(seq_len(n_yrs), dim = c(1, 1, n_yrs)) # one level per year, as setup builds it
  cut <- full
  cut[, , 1] <- NA

  pen_full <- penalty(full)
  pen_cut <- penalty(cut)

  expect_equal(as.numeric(pen_cut[1, 1, 1]), 0)
  expect_true(all(is.finite(pen_cut)))
  expect_equal(pen_cut[, , -1], pen_full[, , -1], tolerance = 1e-12)

  # the total falls by exactly the first year's own contribution
  expect_equal(sum(pen_full) - sum(pen_cut), as.numeric(pen_full[1, 1, 1]), tolerance = 1e-12)
})


test_that("the count has to leave at least one year in the penalty", {

  expect_error(rec_input(seasonal_M_cfg$n_yrs), "at least one year")
  expect_error(rec_input(-1), "whole number")
  expect_error(rec_input(1.5), "whole number")
  expect_error(rec_input(c(1, 2)), "whole number")
})


test_that("the setting survives fit_model, which refreshes the mirror from the map", {

  sim <- seasonal_M_sim()
  base <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))
  base$par$ln_RecDevs[] <- stats::rnorm(length(base$par$ln_RecDevs), 0, 0.3)

  cut <- SPoRC:::do_RecDevs_mapping(base, RecDevs_spec = NULL, rec_dd = 0,
                                    dont_pen_recdev_first = 1)

  fit_base <- fit_model(data = base$data, parameters = base$par, mapping = base$map,
                        random = NULL, do_optim = FALSE, silent = TRUE)
  fit_cut <- fit_model(data = cut$data, parameters = cut$par, mapping = cut$map,
                       random = NULL, do_optim = FALSE, silent = TRUE)

  expect_equal(as.numeric(fit_cut$rep$Rec_nLL[1, 1, 1]), 0)
  expect_equal(sum(fit_base$rep$Rec_nLL) - sum(fit_cut$rep$Rec_nLL),
               as.numeric(fit_base$rep$Rec_nLL[1, 1, 1]), tolerance = 1e-10)

  # the numbers at age are untouched, since only a penalty was removed
  expect_equal(fit_cut$rep$NAA, fit_base$rep$NAA, tolerance = 1e-12)
})


test_that("the setting travels in the data list so the mirror refresh cannot undo it", {

  sim <- seasonal_M_sim()
  base <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))
  cut <- SPoRC:::do_RecDevs_mapping(base, RecDevs_spec = NULL, rec_dd = 0, dont_pen_recdev_first = 1)

  expect_equal(cut$data$dont_pen_recdev_first, 1L)

  # fit_model refreshes every map mirror from the map, so the setting has to be re-applied there
  refreshed <- SPoRC:::sync_dev_map_data(cut$data, cut$map)
  expect_true(is.na(refreshed$map_ln_RecDevs[1, 1, 1]))
  expect_false(any(is.na(refreshed$map_ln_RecDevs[1, 1, -1])))

  # and a model that never asked for it is refreshed the ordinary way
  plain <- SPoRC:::sync_dev_map_data(base$data, base$map)
  expect_false(any(is.na(plain$map_ln_RecDevs)))
})
