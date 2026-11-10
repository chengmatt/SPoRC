# Bins a fleet does not select: their selectivity is exactly zero in every year, the rest of the curve
# keeps its shape, and the operating model builds the same curve.

zero_bins_input <- function(zero_bins, fishsel = list()) {
  suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 2),
    fishsel = fishsel,
    srvsel = list(srv_sel_model = c("desclogist1_Fleet_1", "logist1_Fleet_2"), srv_sel_zero_bins = zero_bins)))
}

zero_bins_fit <- function(il) suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))

test_that("the named bins are zero and every other bin is unchanged", {
  zeroed <- zero_bins_input(list(1:2, NULL))
  plain <- zero_bins_input(NULL)
  expect_null(plain$data$srv_sel_zero_bins)
  expect_equal(which(zeroed$data$srv_sel_zero_bins[,1] == 1), 1:2)
  expect_equal(sum(zeroed$data$srv_sel_zero_bins[,2]), 0)

  sel_zeroed <- zero_bins_fit(zeroed)$rep$srv_sel
  sel_plain <- zero_bins_fit(plain)$rep$srv_sel
  n_ages <- dim(sel_plain)[5]
  expect_true(all(sel_zeroed[,,,,1:2,,1] == 0)) # every year
  expect_identical(sel_zeroed[,,,,3:n_ages,,1], sel_plain[,,,,3:n_ages,,1])
  expect_identical(sel_zeroed[,,,,,,2], sel_plain[,,,,,,2]) # the other fleet
  expect_true(all(sel_plain[,,,,1:2,,1] > 0))
})

test_that("the index leaves those ages out and the gradient stays finite", {
  zeroed <- zero_bins_fit(zero_bins_input(list(1:2, NULL)))
  plain <- zero_bins_fit(zero_bins_input(NULL))
  expect_true(all(zeroed$rep$PredSrvIdx[1,1,,1,1] < plain$rep$PredSrvIdx[1,1,,1,1]))
  expect_equal(zeroed$rep$PredSrvIdx[1,1,,1,2], plain$rep$PredSrvIdx[1,1,,1,2], tolerance = 1e-12)
  expect_true(all(is.finite(zeroed$gr(zeroed$par))))
})

test_that("the operating model builds the same selectivity from the same flags", {
  il <- zero_bins_input(list(1:2, NULL))
  sim_args <- SPoRC:::sim_selex_args("srv", il$data, il$par)
  expect_identical(sim_args$zero_bins, il$data$srv_sel_zero_bins)
  expect_equal(do.call(Get_Selex_Array, sim_args)$sel, zero_bins_fit(il)$rep$srv_sel, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("a fishery fleet takes the option too, and a fixed selectivity is left as given", {
  il <- zero_bins_input(NULL, fishsel = list(fish_sel_zero_bins = list(1)))
  expect_true(all(zero_bins_fit(il)$rep$fish_sel[,,,,1,,1] == 0))

  n_yrs <- length(il$data$years)
  n_ages <- length(il$data$ages)
  fixed <- suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 2),
    fishsel = list(fish_sel_zero_bins = list(1), use_fixed_fish_sel = 1, fish_fixed_sel_pars_spec = "fix",
                   fish_sel_input = array(0.5, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1)))))
  expect_true(all(zero_bins_fit(fixed)$rep$fish_sel[,,1:n_yrs,,,,1] == 0.5))
})

test_that("a list of the wrong shape, bins outside the range, or every bin is refused", {
  expect_error(zero_bins_input(1:2), "srv_sel_zero_bins must be a list")
  expect_error(zero_bins_input(list(1:2)), "srv_sel_zero_bins must be a list")
  expect_error(zero_bins_input(list(99, NULL)), "outside 1:")
  expect_error(zero_bins_input(list(1:7, NULL)), "names every bin")
})
