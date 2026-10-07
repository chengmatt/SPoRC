# A selectivity deviation series shared over fleets (est_shared_f_x) or regions is one series with one variance:
# penalized once, and its sigmas shared the same way, so no sigma is left that no penalty reads.

dev_value_fs <- 0.3 # every estimated deviation set here, so the penalty has a closed form

fleet_share_input <- function(n_regions, devs_spec, tv = "iid") {
  suppressMessages(sweep_input(dims = list(n_regions = n_regions, n_sexes = 1, n_fish_fleets = 2, n_srv_fleets = 1),
    fishsel = list(cont_tv_fish_sel = paste0(tv, "_Fleet_", 1:2), fishsel_pe_pars_spec = c("est_all", "est_all"),
                   fish_sel_devs_spec = devs_spec)))
}

# the selectivity penalty at every deviation set to dev_value_fs and every sigma at one
fleet_share_penalty <- function(il) {
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  pars <- obj$par
  pars[names(pars) == "ln_fishsel_devs"] <- dev_value_fs
  pars[names(pars) == "fishsel_pe_pars"] <- 0
  obj$report(pars)$sel_nLL
}

n_dev_levels <- function(il) length(unique(stats::na.omit(as.numeric(il$map$ln_fishsel_devs))))

test_that("a deviation series shared over fleets is penalized once", {
  shared <- fleet_share_input(1, c("est_all", "est_shared_f_1"))
  own <- fleet_share_input(1, c("est_all", "est_all"))
  per_dev <- -dnorm(dev_value_fs, 0, 1, log = TRUE)
  expect_equal(n_dev_levels(own), 2 * n_dev_levels(shared))
  expect_equal(fleet_share_penalty(shared), n_dev_levels(shared) * per_dev)
  expect_equal(fleet_share_penalty(own), n_dev_levels(own) * per_dev)
})

test_that("sigmas are shared wherever the deviations are, over regions and over fleets", {
  il <- fleet_share_input(2, c("est_shared_r", "est_shared_f_1"))
  pe_map <- array(il$map$fishsel_pe_pars, dim = dim(il$par$fishsel_pe_pars)) # region, slot, sex, fleet
  expect_identical(pe_map[1,,,1], pe_map[2,,,1]) # both regions read fleet one's
  expect_identical(pe_map[,,,1], pe_map[,,,2]) # fleet two reads fleet one's
  expect_equal(length(unique(stats::na.omit(as.vector(pe_map)))), 2) # one per logistic parameter
})

test_that("every estimated sigma enters the penalty", {
  for(devs_spec in list(c("est_shared_r", "est_shared_f_1"), c("est_shared_r_s", "est_all"))) {
    il <- fleet_share_input(2, devs_spec)
    obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
    pars <- obj$par
    set.seed(5)
    is_dev <- names(pars) == "ln_fishsel_devs"
    pars[is_dev] <- stats::rnorm(sum(is_dev), 0, 0.3) # off zero, where a sigma's gradient can be zero by symmetry
    grad <- obj$gr(pars)[names(pars) == "fishsel_pe_pars"]
    expect_true(all(grad != 0), info = paste(devs_spec, collapse = " "))
  } # end devs_spec loop
})
