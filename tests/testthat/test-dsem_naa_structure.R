# Checks a numbers at age process error written in arrows with no covariate: the iid arrows give the native iid
# penalty's joint density (1e-8) once the year-one constant is removed, and the 2dar1 arrows parse and evaluate.

data("sgl_rg_dusky_data")

naa_series <- function(a) sprintf("NAA_Pop_1_Region_1_Seas_1_Age_%d_Sex_1", a)
naa_arrows <- function(structure, ages) {
  sd_lines <- paste0(naa_series(ages), " <-> ", naa_series(ages), ", 0, sd_naa")
  year_lines <- paste0(naa_series(ages), " -> ", naa_series(ages), ", 1, rho_yr")
  age_lines <- paste0(naa_series(ages[-1] - 1), " -> ", naa_series(ages[-1]), ", 0, rho_age")
  cohort_lines <- paste0(naa_series(ages[-1] - 1), " -> ", naa_series(ages[-1]), ", 1, rho_cohort")
  paste(switch(structure, iid = sd_lines, "2dar1" = c(year_lines, age_lines, cohort_lines, sd_lines)))
}

test_that("iid arrows with no covariate reproduce the native iid penalty", {

  build_naa <- build_goa_dusky_input
  body(build_naa) <- do.call(substitute, list(body(build_naa), list(Setup_Mod_Biologicals = quote(function(...) Setup_Mod_Biologicals(..., NAA_re = "iid", NAA_re_ages = 5:9, NAA_sigma_spec = "fix")))))
  il <- suppressMessages(build_naa(sgl_rg_dusky_data))
  n_yrs <- length(il$data$years)
  ages <- 2:6 # age indices of the state
  base_rep <- fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)$rep
  set.seed(2)
  il$par$ln_NAA[] <- log(base_rep$NAA[,,1:n_yrs,,,]) + rnorm(length(il$par$ln_NAA), 0, 0.2) # a state off its prediction
  native <- fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)

  # the whole process error model is the arrows: no data frame, no covariate column in the grid
  d <- suppressMessages(Setup_Mod_DSEM(il, hold_arrows(naa_arrows("iid", ages), c(sd_naa = 0.3)), dsem_data = NULL, dsem_processes = "NAA"))
  expect_equal(ncol(d$data$dsem_cov_obs), 0L)
  expect_length(d$par$ln_dsem_obs_sd, 0)
  expect_equal(d$data$dsem_var_names, naa_series(ages))
  expect_equal(sum(!is.na(d$map$dsem_x)), 0L) # every grid cell is a linked deviation or a fixed leading row

  # the year-one rows are not state cells, so they sit fixed at the mean and add a known constant
  via_arrows <- fit_model(d$data, d$par, d$map, random = NULL, do_optim = FALSE, silent = TRUE)
  year_one <- length(ages) * (log(0.3) + 0.5 * log(2 * pi))
  expect_equal(via_arrows$fn(via_arrows$par) - year_one, native$fn(native$par), tolerance = 1e-8)
  expect_equal(via_arrows$rep$NAA_state_nLL, 0, tolerance = 1e-12) # the native penalty is off for every state cell
  expect_equal(via_arrows$rep$dsem_nLL - year_one, native$rep$NAA_state_nLL, tolerance = 1e-8)
  expect_equal(sum(is.na(via_arrows$data$map_ln_NAA[,,2:n_yrs,1,ages,1])), (n_yrs - 1) * length(ages))

  # a data frame with no covariate column is still refused, and a covariate can be added on top
  expect_error(Setup_Mod_DSEM(il, naa_arrows("iid", ages), data.frame(year = il$data$years), dsem_processes = "NAA"), "no covariate columns")

})

test_that("the 2dar1 arrows parse with three correlations and evaluate", {

  build_naa <- build_goa_dusky_input
  body(build_naa) <- do.call(substitute, list(body(build_naa), list(Setup_Mod_Biologicals = quote(function(...) Setup_Mod_Biologicals(..., NAA_re = "iid", NAA_re_ages = 5:9, NAA_sigma_spec = "fix")))))
  il <- suppressMessages(build_naa(sgl_rg_dusky_data))
  ages <- 2:6

  expect_equal(read_dsem_arrows(naa_arrows("2dar1", ages), naa_series(ages))$beta_names, c("rho_yr", "rho_age", "rho_cohort"))
  d <- suppressMessages(Setup_Mod_DSEM(il, hold_arrows(naa_arrows("2dar1", ages), c(rho_yr = 0.6, rho_age = 0.4, rho_cohort = -0.24)),
                                       dsem_data = NULL, dsem_processes = "NAA"))
  expect_length(d$data$dsem_model$beta_names, 0) # held through the arrows, so no coefficient parameter is left
  expect_equal(d$data$dsem_model$ln_sd_names, "sd_naa")
  expect_equal(nrow(d$data$dsem_model$arrows), 3 * length(ages) - 2 + length(ages))
  obj <- fit_model(d$data, d$par, d$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(is.finite(obj$rep$dsem_nLL))

})
