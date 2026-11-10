# Selectivity deviations added on the natural scale: under "additive" a parameter takes par + dev where the
# default is par * exp(dev), in the curve, in the penalty and in the operating model's selectivity alike.

devs_one_year <- function(values) {
  devs <- array(0, dim = c(1, 2, 12, 1, 1)) # region, year, slot, sex, fleet
  devs[1,2,seq_along(values),1,1] <- values
  devs
}

curve_add <- function(model, pars, devs, additive = NULL) {
  as.numeric(Get_Selex(Selex_Model = model, TimeVary_Model = 2, pars = log(pars), ln_seldevs = devs,
                       Region = 1, Year = 2, Bin = 1:12, Sex = 1, devs_additive = additive))
}

# fleet one a double logistic walking on its ascending limb, fleet two a time-invariant logistic
additive_input <- function(devs_form) {
  suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 2, n_srv_fleets = 1),
    fishsel = list(fish_sel_model = c("dbllogist1_Fleet_1", "logist1_Fleet_2"),
                   cont_tv_fish_sel = c("rw_Fleet_1", "none_Fleet_2"),
                   fishsel_pe_pars_spec = c("fix", "none"),
                   fish_sel_devs_spec = c("est_all", "none"),
                   fish_sel_devs_pars = list(1:2, NULL),
                   fish_sel_devs_form = devs_form)))
}

test_that("an additive deviation moves the parameter by the deviation itself", {
  bins <- 1:12
  devs <- devs_one_year(c(0.4, 0.1))

  # logistic and descending logistic, age at 50% added and slope multiplied
  expect_equal(curve_add(0, c(4, 1.2), devs, additive = c(TRUE, FALSE)), 1 / (1 + exp(-1.2 * exp(0.1) * (bins - 4.4))), tolerance = 1e-12)
  expect_equal(curve_add(11, c(4, 1.2), devs, additive = c(TRUE, FALSE)), 1 / (1 + exp(1.2 * exp(0.1) * (bins - 4.4))), tolerance = 1e-12)

  # double logistic, both ascending parameters added
  by_hand <- (1 / (1 + exp(-(1.2 + 0.1) * (bins - 3.4)))) / (1 + exp(0.9 * (bins - 8)))
  expect_equal(curve_add(12, c(3, 1.2, 8, 0.9), devs, additive = c(TRUE, TRUE, FALSE, FALSE)), by_hand, tolerance = 1e-12)

  # no parameter named, every deviation multiplies
  expect_equal(curve_add(0, c(4, 1.2), devs), 1 / (1 + exp(-1.2 * exp(0.1) * (bins - 4 * exp(0.4)))), tolerance = 1e-12)
  expect_identical(curve_add(0, c(4, 1.2), devs), curve_add(0, c(4, 1.2), devs, additive = c(FALSE, FALSE)))
})

test_that("setup flags the additive parameters and the model's curve follows them", {
  ascending_age <- list(c("additive", "multiplicative", "multiplicative", "multiplicative"), "multiplicative")
  il <- additive_input(ascending_age)
  expect_equal(dim(il$data$fish_sel_devs_additive), c(length(il$data$ages), 2))
  expect_equal(which(il$data$fish_sel_devs_additive[,1] == 1), 1)
  expect_equal(sum(il$data$fish_sel_devs_additive[,2]), 0)

  # one string covers every parameter of the fleet, and the default adds none
  expect_true(all(additive_input(c("additive", "multiplicative"))$data$fish_sel_devs_additive[,1] == 1))
  expect_equal(sum(additive_input(c("multiplicative", "multiplicative"))$data$fish_sel_devs_additive), 0)

  il$par$fish_fixed_sel_pars[1,,1,1,1] <- log(c(3, 1.2, 6, 0.9))
  il$par$ln_fishsel_devs[1,2,1:2,1,1] <- c(0.4, 0.1)
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  ages <- il$data$ages
  by_hand <- (1 / (1 + exp(-1.2 * exp(0.1) * (ages - 3.4)))) / (1 + exp(0.9 * (ages - 6)))
  expect_equal(as.numeric(obj$rep$fish_sel[1,1,2,1,,1,1]), by_hand, tolerance = 1e-12)

  # the operating model builds the same curve from the same flags
  sim_args <- SPoRC:::sim_selex_args("fish", il$data, il$par)
  expect_identical(sim_args$devs_additive, il$data$fish_sel_devs_additive)
  expect_equal(do.call(Get_Selex_Array, sim_args)$sel[1,1,2,1,,1,1], as.numeric(obj$rep$fish_sel[1,1,2,1,,1,1]), tolerance = 1e-12)
})

test_that("the walk penalty is on the deviations as given, so its sd is on the parameter's own scale", {
  il <- additive_input(list(c("additive", "multiplicative", "multiplicative", "multiplicative"), "multiplicative"))
  n_yrs <- dim(il$par$ln_fishsel_devs)[2]
  step <- 0.07
  il$par$ln_fishsel_devs[1,,1,1,1] <- step * (1:n_yrs) # the age at 50% rises 0.07 a year
  il$par$fishsel_pe_pars[1,1,1,1] <- log(0.2) # 0.2 years of age
  il$par$fishsel_pe_pars[1,2,1,1] <- log(0.05)
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  first <- il$data$fishsel_rw_init_sigma[1]
  by_hand <- -dnorm(step, 0, first, log = TRUE) - (n_yrs - 1) * dnorm(step, 0, 0.2, log = TRUE) - # age at 50%
    dnorm(0, 0, first, log = TRUE) - (n_yrs - 1) * dnorm(0, 0, 0.05, log = TRUE) # slope, at zero
  expect_equal(obj$rep$sel_nLL, by_hand, tolerance = 1e-10)
})

test_that("a wrong length, an unknown form, or a fleet whose deviations are not on curve parameters is refused", {
  expect_error(additive_input("additive"), "must have one element per fleet")
  expect_error(additive_input(c("additive", "added")), "must be 'multiplicative' or 'additive'")
  expect_error(additive_input(c("additive", "additive")), "not iid or a random walk")
  expect_error(additive_input(list(rep("additive", 99), "multiplicative")), "more entries than")
})

test_that("a survey whose logistic curve walks takes additive deviations the same way", {
  il <- suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1),
    srvsel = list(srv_sel_model = "logist1_Fleet_1",
                  cont_tv_srv_sel = "rw_Fleet_1",
                  srvsel_pe_pars_spec = "fix",
                  srv_sel_devs_spec = "est_all",
                  srv_sel_devs_form = list(c("additive", "multiplicative")))))
  expect_equal(dim(il$data$srv_sel_devs_additive), c(length(il$data$ages), 1))
  expect_equal(which(il$data$srv_sel_devs_additive[,1] == 1), 1)

  # age at 50% moves by the deviation itself, the slope in proportion
  il$par$srv_fixed_sel_pars[1,,1,1,1] <- log(c(3, 1.2))
  il$par$ln_srvsel_devs[1,2,1:2,1,1] <- c(0.4, 0.1)
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  ages <- il$data$ages
  expect_equal(as.numeric(obj$rep$srv_sel[1,1,2,1,,1,1]), 1 / (1 + exp(-1.2 * exp(0.1) * (ages - 3.4))), tolerance = 1e-12)

  # the operating model builds the same curve from the same flags
  sim_args <- SPoRC:::sim_selex_args("srv", il$data, il$par)
  expect_identical(sim_args$devs_additive, il$data$srv_sel_devs_additive)
  expect_equal(do.call(Get_Selex_Array, sim_args)$sel[1,1,2,1,,1,1], as.numeric(obj$rep$srv_sel[1,1,2,1,,1,1]), tolerance = 1e-12)
})
