library(SPoRC)
library(testthat)

# The est_shared_b family estimates one deviation series per bin group per year. Under iid or a
# random walk it is legal for non-parametric fleets, whose deviation slots are bins rather than parameters.

n_yrs <- 12
n_ages <- 6
dev_groups <- list(1, 2, 3, 4, 5:6) # five groups over six bins, the last holding two
n_groups <- length(dev_groups)

# build the sweep model with one non-parametric fleet under continuous time-varying selectivity
sel_devs_input <- function(prefix, spec, tv, n_regions, n_sexes) {

  dims <- list(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_yrs = n_yrs,
    n_ages = n_ages
  )

  # total fishery selectivity
  if(prefix == "fish") {
    il <- sweep_input(dims = dims, fishsel = list(
      fish_sel_model = "nonparfree_Fleet_1",
      cont_tv_fish_sel = paste0(tv, "_Fleet_1"),
      fish_sel_nonpar_est_bins = list(list(as.list(1:n_ages))),
      fish_fixed_sel_pars_spec = "est_all",
      fishsel_pe_pars_spec = "est_all",
      fish_sel_devs_spec = spec,
      fishsel_devs_shared_bins = dev_groups
    ))
  }

  # retained fishery selectivity
  if(prefix == "ret") {
    il <- sweep_input(dims = dims, fishsel = list(
      use_fixed_ret_sel = 0,
      ret_sel_model = "nonparfree_Fleet_1",
      cont_tv_ret_sel = paste0(tv, "_Fleet_1"),
      ret_sel_nonpar_est_bins = list(list(as.list(1:n_ages))),
      ret_fixed_sel_pars_spec = "est_all",
      retsel_pe_pars_spec = "est_all",
      ret_sel_devs_spec = spec,
      retsel_devs_shared_bins = dev_groups
    ))
  }

  # survey selectivity
  if(prefix == "srv") {
    il <- sweep_input(dims = dims, srvsel = list(
      srv_sel_model = "nonparfree_Fleet_1",
      cont_tv_srv_sel = paste0(tv, "_Fleet_1"),
      srv_sel_nonpar_est_bins = list(list(as.list(1:n_ages))),
      srv_fixed_sel_pars_spec = "est_all",
      srvsel_pe_pars_spec = "est_all",
      srv_sel_devs_spec = spec,
      srvsel_devs_shared_bins = dev_groups
    ))
  }

  il
}

# the deviation map alone, which most of these checks read
sel_devs_map <- function(prefix, spec, tv, n_regions, n_sexes) {
  sel_devs_input(prefix, spec, tv, n_regions, n_sexes)$map[[paste0("ln_", prefix, "sel_devs")]]
}

# number of distinct deviations actually estimated
n_estimated <- function(map) length(unique(stats::na.omit(as.numeric(map))))

test_that("est_shared_b variants estimate one deviation per bin group per year under iid and rw", {

  n_regions <- 2
  n_sexes <- 2

  # each variant drops the dims it shares over from the deviation count
  expected <- c(
    est_shared_b     = n_yrs * n_groups * n_regions * n_sexes,
    est_shared_r_b   = n_yrs * n_groups * n_sexes,
    est_shared_b_s   = n_yrs * n_groups * n_regions,
    est_shared_r_b_s = n_yrs * n_groups
  )

  for(prefix in c("fish", "ret", "srv")) {
    for(tv in c("iid", "rw")) {
      for(spec in names(expected)) {

        map <- sel_devs_map(prefix, spec, tv, n_regions, n_sexes)
        label <- paste(prefix, tv, spec)

        expect_false(all(is.na(map)), info = label)
        expect_equal(n_estimated(map), unname(expected[[spec]]), info = label)

      } # end spec loop
    } # end tv loop
  } # end prefix loop
})

test_that("iid and rw deviation counts match the 2dar1 form, which indexes deviations by bin too", {

  for(spec in c("est_shared_b", "est_shared_r_b", "est_shared_b_s", "est_shared_r_b_s")) {

    gmrf <- n_estimated(sel_devs_map("fish", spec, "2dar1", n_regions = 2, n_sexes = 2))

    for(tv in c("iid", "rw")) {
      expect_equal(n_estimated(sel_devs_map("fish", spec, tv, n_regions = 2, n_sexes = 2)), gmrf,
                   info = paste(tv, spec))
    } # end tv loop

  } # end spec loop
})

test_that("a single region and sex gives one deviation per bin group per year for every variant", {

  for(spec in c("est_shared_b", "est_shared_r_b", "est_shared_b_s", "est_shared_r_b_s")) {
    for(tv in c("iid", "rw")) {

      map <- sel_devs_map("fish", spec, tv, n_regions = 1, n_sexes = 1)
      expect_equal(n_estimated(map), n_yrs * n_groups, info = paste(tv, spec))

    } # end tv loop
  } # end spec loop
})

test_that("est_shared_b estimates one process error sigma per bin group, each with a non-zero gradient", {

  # under iid or a walk the sigmas are indexed by bin, and a bin group leaves one deviation series, so
  # a sigma at any bin above the group's lowest would be estimated without ever reaching the likelihood
  for(prefix in c("fish", "ret", "srv")) {
    for(tv in c("iid", "rw")) {

      il <- sel_devs_input(prefix, "est_shared_b", tv, n_regions = 1, n_sexes = 1)
      par_nm <- paste0(prefix, "sel_pe_pars")
      label <- paste(prefix, tv)

      # one sigma per group, not one per bin
      expect_equal(n_estimated(il$map[[par_nm]]), n_groups, info = label)

      obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
      grad <- as.numeric(obj$gr(obj$par))[names(obj$par) == par_nm]

      # every sigma the map estimates moves the objective
      expect_length(grad, n_groups)
      expect_true(all(grad != 0), info = label)

    } # end tv loop
  } # end prefix loop
})

test_that("a fleet borrowing another fleet's bin-grouped deviations gets one sigma per group too", {

  # est_shared_f copies the reference fleet's deviations wholesale, bin groups included, so the
  # borrowing fleet's own sigmas answer to the reference fleet's spec rather than to its own
  il <- sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 2, n_srv_fleets = 1,
                n_yrs = n_yrs, n_ages = n_ages),
    fishsel = list(
      fish_sel_model = c("nonparfree_Fleet_1", "nonparfree_Fleet_2"),
      cont_tv_fish_sel = c("rw_Fleet_1", "rw_Fleet_2"),
      fish_sel_nonpar_est_bins = list(list(as.list(1:n_ages)), list(as.list(1:n_ages))),
      fish_fixed_sel_pars_spec = c("est_all", "est_all"),
      fishsel_pe_pars_spec = c("est_all", "est_all"),
      fish_sel_devs_spec = c("est_shared_b", "est_shared_f_1"),
      fishsel_devs_shared_bins = dev_groups
    ))

  # one sigma per group for each of the two fleets
  expect_equal(n_estimated(il$map$fishsel_pe_pars), 2 * n_groups)

  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  grad <- as.numeric(obj$gr(obj$par))[names(obj$par) == "fishsel_pe_pars"]

  expect_true(all(grad != 0))
})

test_that("sharing deviations across sexes leaves one sigma per sex-shared series", {

  # the likelihood evaluates a sex-shared deviation series at the first sex only, so a sigma at any
  # sex above it would be estimated without ever reaching the likelihood
  n_regions <- 2
  n_sexes <- 2

  # sharing over sexes halves the sigma count, and over bin groups takes it from bins to groups
  expected <- c(
    est_shared_s     = n_regions * n_ages,
    est_shared_r_s   = n_regions * n_ages,
    est_shared_b_s   = n_regions * n_groups,
    est_shared_r_b_s = n_regions * n_groups
  )

  for(prefix in c("fish", "ret", "srv")) {
    for(tv in c("iid", "rw")) {
      for(spec in names(expected)) {

        il <- sel_devs_input(prefix, spec, tv, n_regions, n_sexes)
        par_name <- paste0(prefix, "sel_pe_pars")
        label <- paste(prefix, tv, spec)

        expect_equal(n_estimated(il$map[[par_name]]), unname(expected[[spec]]), info = label)

        obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
        grad <- as.numeric(obj$gr(obj$par))[names(obj$par) == par_name]

        # every sigma the map estimates moves the objective
        expect_length(grad, unname(expected[[spec]]))
        expect_true(all(grad != 0), info = label)

      } # end spec loop
    } # end tv loop
  } # end prefix loop
})
