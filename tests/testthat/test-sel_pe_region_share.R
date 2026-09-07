# The selectivity process error penalty is evaluated one region at a time, so a deviation series shared
# across regions must be split over the regions holding it rather than penalized once in each.

library(SPoRC)
library(testthat)

n_yrs <- 10
n_ages <- 6
dev_value <- 0.3 # every deviation set here, so the penalty has a closed form

# one non-parametric fishery fleet under continuous time-varying selectivity
sel_pe_input <- function(n_regions, n_sexes, tv, devs_spec) {

  sweep_input(
    dims = list(n_regions = n_regions, n_sexes = n_sexes, n_fish_fleets = 1,
                n_srv_fleets = 1, n_yrs = n_yrs, n_ages = n_ages),
    fishsel = list(
      fish_sel_model = "nonparfree_Fleet_1",
      cont_tv_fish_sel = paste0(tv, "_Fleet_1"),
      fish_sel_nonpar_est_bins = list(list(as.list(1:n_ages))),
      fish_fixed_sel_pars_spec = "est_all",
      fishsel_pe_pars_spec = "est_all",
      fish_sel_devs_spec = devs_spec
    ))
}

# the selectivity penalty with every deviation held at the same value
sel_penalty <- function(il) {
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  pars <- obj$par
  pars[names(obj$par) == "ln_fishsel_devs"] <- dev_value
  obj$fn(pars)
  obj$report(pars)$sel_nLL
}

# number of distinct deviations actually estimated
n_devs <- function(il) length(unique(stats::na.omit(as.numeric(il$map$ln_fishsel_devs))))


test_that("a region-shared deviation series is penalized once, not once per region", {

  # region sharing leaves the deviation count flat, so the penalty has to stay flat with it
  for(tv in c("iid", "rw", "3dmarg", "2dar1")) {
    for(spec in c("est_shared_r", "est_shared_r_s")) {

      inputs <- lapply(1:3, function(n_regions) sel_pe_input(n_regions, 2, tv, spec))
      label <- paste(tv, spec)

      expect_equal(sapply(inputs, n_devs), rep(n_devs(inputs[[1]]), 3), info = label)
      expect_equal(sapply(inputs, sel_penalty), rep(sel_penalty(inputs[[1]]), 3), info = label)

    } # end spec loop
  } # end tv loop
})


test_that("deviations that are not shared still contribute one penalty each", {

  # est_all gives every region its own series, so here the penalty should scale with region count
  for(tv in c("iid", "rw")) {

    inputs <- lapply(1:3, function(n_regions) sel_pe_input(n_regions, 2, tv, "est_all"))
    one_region <- sel_penalty(inputs[[1]])

    expect_equal(sapply(inputs, n_devs), n_devs(inputs[[1]]) * 1:3, info = tv)
    expect_equal(sapply(inputs, sel_penalty), one_region * 1:3, info = tv)

  } # end tv loop
})


test_that("the iid penalty on a shared series matches the density written out by hand", {

  # one series of n_yrs x n_ages deviations at 0.3 under sigma = exp(0), whatever the region count
  by_hand <- -sum(dnorm(rep(dev_value, n_yrs * n_ages), 0, 1, log = TRUE))

  for(n_regions in 1:3) {
    expect_equal(sel_penalty(sel_pe_input(n_regions, 1, "iid", "est_shared_r")), by_hand,
                 info = paste(n_regions, "regions"))
  } # end n_regions loop
})


test_that("summing the per-unit calls gives one penalty however the levels are shared", {

  # Called directly, since the split is the same contract for selectivity, which loops regions, and
  # for growth, which loops populations by region. n_units is whatever the first dim holds.
  n_units <- 3
  n_bins <- 2
  set.seed(42)
  devs <- array(rnorm(n_yrs * n_bins), dim = c(1, n_yrs, n_bins, 1, 1))
  pe_pars <- array(0, dim = c(1, n_bins, 1, 1)) # sigma = exp(0)

  # every unit holds the same levels, which is what sharing over that dim produces
  one_unit <- array(seq_len(n_yrs * n_bins), dim = c(1, n_yrs, n_bins, 1))
  shared <- array(rep(as.vector(one_unit), n_units), dim = c(n_units, n_yrs, n_bins, 1))

  # each unit holds its own levels, so nothing is shared and nothing is split
  unshared <- array(seq_len(n_units * n_yrs * n_bins), dim = c(n_units, n_yrs, n_bins, 1))

  call_unit <- function(map_full, unit) {
    SPoRC:::Get_PE_loglik(
      PE_model = 1,
      PE_pars = pe_pars,
      ln_devs = devs,
      map_sel_devs = array(map_full[unit,,,], dim = c(1, n_yrs, n_bins, 1)),
      map_sel_devs_full = map_full,
      min_sel_devs_shared_bins = 1:n_bins
    )
  }

  # the shared series is one set of deviations, so the units together owe one penalty
  by_hand <- sum(dnorm(as.vector(devs), 0, 1, log = TRUE))
  expect_equal(sum(sapply(1:n_units, function(u) call_unit(shared, u))), by_hand)

  # holding its own levels, each unit owes a full penalty of its own
  expect_equal(sum(sapply(1:n_units, function(u) call_unit(unshared, u))), n_units * by_hand)
})
