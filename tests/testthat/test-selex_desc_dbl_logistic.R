# The descending and double logistic forms: the curves against their formulas, the parameter and
# deviation slots setup gives them, and *_sel_devs_pars keeping deviations to the parameters named.

dev_value_dl <- 0.3 # every estimated deviation set here, so the penalty has a closed form

flat_devs_dl <- array(0, dim = c(1, 2, 12, 1, 1)) # region, year, slot, sex, fleet

curve_dl <- function(model, pars, devs = flat_devs_dl, tv = 0, year = 1) {
  as.numeric(Get_Selex(Selex_Model = model, TimeVary_Model = tv, pars = log(pars), ln_seldevs = devs,
                       Region = 1, Year = year, Bin = 1:12, Sex = 1))
}

# fleet one a double logistic on a random walk, fleet two a descending logistic with iid deviations
logistic_forms_input <- function(devs_pars = NULL, tv = c("rw_Fleet_1", "iid_Fleet_2")) {
  suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 2, n_srv_fleets = 1),
    fishsel = list(fish_sel_model = c("dbllogist1_Fleet_1", "desclogist1_Fleet_2"),
                   cont_tv_fish_sel = tv,
                   fishsel_pe_pars_spec = c("est_all", "est_all"),
                   fish_sel_devs_spec = c("est_all", "est_all"),
                   fish_sel_devs_pars = devs_pars)))
}

# which deviation slots of a fleet are estimated in any year
active_slots <- function(il, f) {
  dev_map <- array(as.numeric(il$map$ln_fishsel_devs), dim = dim(il$par$ln_fishsel_devs)) # region, year, slot, sex, fleet
  which(apply(!is.na(dev_map[1,,,1,f]), 2, any))
}

# which sigma slots of a fleet are estimated
sigma_slots <- function(il, f) {
  pe_map <- array(as.numeric(il$map$fishsel_pe_pars), dim = dim(il$par$fishsel_pe_pars)) # region, slot, sex, fleet
  which(!is.na(pe_map[1,,1,f]))
}

test_that("the descending logistic is one minus the ascending logistic", {
  bins <- 1:12
  for(pars in list(c(4, 1.2), c(8, 0.5))) {
    desc <- curve_dl(11, pars)
    expect_equal(desc, 1 - 1 / (1 + exp(-pars[2] * (bins - pars[1]))), tolerance = 1e-12)
    expect_equal(desc, 1 - curve_dl(0, pars), tolerance = 1e-12)
    expect_true(all(diff(desc) < 0)) # falls with every bin
  }
})

test_that("the double logistic is the ascending logistic times the descending one", {
  bins <- 1:12
  pars <- c(3, 1.2, 8, 0.9) # ascending b50 and slope, descending b50 and slope
  dbl <- curve_dl(12, pars)
  by_hand <- (1 / (1 + exp(-1.2 * (bins - 3)))) * (1 - 1 / (1 + exp(-0.9 * (bins - 8))))
  expect_equal(dbl, by_hand, tolerance = 1e-12)
  expect_equal(dbl, curve_dl(0, pars[1:2]) * curve_dl(11, pars[3:4]), tolerance = 1e-12)
  expect_lt(max(dbl), 1) # not rescaled to peak at one

  # a descending limb far past the last bin leaves the ascending logistic
  expect_equal(curve_dl(12, c(3, 1.2, 200, 0.9)), curve_dl(0, c(3, 1.2)), tolerance = 1e-10)
})

test_that("deviations multiply each parameter of both forms", {
  devs <- flat_devs_dl
  devs[1,2,1:4,1,1] <- c(0.1, -0.2, 0.05, 0.3) # year two
  pars <- c(3, 1.2, 8, 0.9)
  for(tv in 1:2) {
    expect_equal(curve_dl(12, pars, devs, tv = tv, year = 2), curve_dl(12, pars * exp(devs[1,2,1:4,1,1])), tolerance = 1e-12)
    expect_equal(curve_dl(11, pars[1:2], devs, tv = tv, year = 2), curve_dl(11, pars[1:2] * exp(devs[1,2,1:2,1,1])), tolerance = 1e-12)
    expect_equal(curve_dl(12, pars, devs, tv = tv, year = 1), curve_dl(12, pars), tolerance = 1e-12) # year one has none
  }
})

test_that("setup gives the two forms their parameters", {
  il <- logistic_forms_input(tv = c("none_Fleet_1", "none_Fleet_2"))
  fixed_map <- array(as.numeric(il$map$fish_fixed_sel_pars), dim = dim(il$par$fish_fixed_sel_pars)) # region, par, block, sex, fleet
  expect_equal(unique(il$data$fish_sel_model[1,,1]), 12)
  expect_equal(unique(il$data$fish_sel_model[1,,2]), 11)
  expect_equal(dim(il$par$fish_fixed_sel_pars)[2], 4)
  expect_equal(sum(!is.na(fixed_map[1,,1,1,1])), 4) # double logistic
  expect_equal(sum(!is.na(fixed_map[1,,1,1,2])), 2) # descending logistic
})

test_that("deviations and sigmas follow the parameters a fleet names", {
  every_par <- logistic_forms_input()
  expect_equal(active_slots(every_par, 1), 1:4)
  expect_equal(active_slots(every_par, 2), 1:2)
  expect_equal(sigma_slots(every_par, 1), 1:4)
  expect_equal(sigma_slots(every_par, 2), 1:2)

  ascending <- logistic_forms_input(devs_pars = list(1:2, NULL)) # fleet one's ascending limb only
  expect_equal(active_slots(ascending, 1), 1:2)
  expect_equal(active_slots(ascending, 2), 1:2)
  expect_equal(sigma_slots(ascending, 1), 1:2)
  expect_equal(sigma_slots(ascending, 2), 1:2)

  slope_only <- logistic_forms_input(devs_pars = list(c(2, 4), 2))
  expect_equal(active_slots(slope_only, 1), c(2, 4))
  expect_equal(active_slots(slope_only, 2), 2)
  expect_equal(sigma_slots(slope_only, 1), c(2, 4))
  expect_equal(sigma_slots(slope_only, 2), 2)
  expect_identical(is.na(slope_only$data$map_ln_fishsel_devs),
                   is.na(array(as.numeric(slope_only$map$ln_fishsel_devs), dim = dim(slope_only$par$ln_fishsel_devs))))
})

test_that("the penalty covers the named parameters only, and every sigma enters it", {
  il <- logistic_forms_input(devs_pars = list(1:2, NULL))
  n_yrs <- dim(il$par$ln_fishsel_devs)[2]
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  pars <- obj$par
  pars[names(pars) == "ln_fishsel_devs"] <- dev_value_dl
  pars[names(pars) == "fishsel_pe_pars"] <- 0 # sigma of one

  # fleet one: two walks, year one under the wide starting sd and no step after it
  walk <- -dnorm(dev_value_dl, 0, il$data$fishsel_rw_init_sigma[1], log = TRUE) - (n_yrs - 1) * dnorm(0, 0, 1, log = TRUE)
  # fleet two: two iid series
  iid <- -n_yrs * dnorm(dev_value_dl, 0, 1, log = TRUE)
  expect_equal(obj$report(pars)$sel_nLL, 2 * walk + 2 * iid, tolerance = 1e-10)

  set.seed(11)
  pars[names(pars) == "ln_fishsel_devs"] <- rnorm(sum(names(pars) == "ln_fishsel_devs"), 0, 0.2)
  grad <- obj$gr(pars)
  expect_equal(sum(names(pars) == "fishsel_pe_pars"), 4) # two for each fleet
  expect_true(all(grad[names(pars) == "fishsel_pe_pars"] != 0))
  expect_true(all(grad[names(pars) == "ln_fishsel_devs"] != 0))
})

test_that("the operating model draws deviations for the named parameters only", {
  il <- logistic_forms_input(devs_pars = list(1:2, NULL))
  pe_pars <- il$par$fishsel_pe_pars
  pe_pars[] <- log(0.2)
  set.seed(3)
  drawn <- draw_sel_dev_surface(PE_codes = il$data$cont_tv_fish_sel,
                                map = il$data$map_ln_fishsel_devs,
                                pe_pars = pe_pars,
                                fit_devs = il$par$ln_fishsel_devs,
                                bins = il$data$fishsel_devs_min_shared_bins,
                                n_cond = 0,
                                pe_wt = il$data$fishsel_pe_wt,
                                rw_init_sigma = il$data$fishsel_rw_init_sigma)
  expect_true(all(drawn[1,-1,1:2,1,1] != 0)) # ascending limb drawn, the walk starting at zero in year one
  expect_true(all(drawn[1,,3:4,1,1] == 0)) # descending limb left alone
  expect_true(all(drawn[1,,1:2,1,2] != 0)) # fleet two's iid series, both parameters
})

test_that("a list of the wrong shape, or parameters a form does not have, is refused", {
  expect_error(logistic_forms_input(devs_pars = 1:2), "fish_sel_devs_pars must be a list")
  expect_error(logistic_forms_input(devs_pars = list(1:2)), "fish_sel_devs_pars must be a list")
  expect_error(logistic_forms_input(devs_pars = list(1:5, NULL)), "between 1 and 4")
  expect_error(logistic_forms_input(devs_pars = list(c(1, 1), NULL)), "distinct parameters")
  expect_error(logistic_forms_input(devs_pars = list(NULL, 3)), "between 1 and 2")
  expect_error(logistic_forms_input(devs_pars = list(1:2, NULL), tv = c("3dmarg_Fleet_1", "iid_Fleet_2")),
               "not iid or a random walk")
})
