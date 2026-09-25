# Checks SPoRC_rtmb against SPoRC_rtmb plus a hand-built grid density (1e-8): linked series lose
# SPoRC's recruitment penalty, unlinked ones keep it, and each series reads the right deviations.

data("sgl_rg_dusky_data")

objective_value <- function(model, input_list) {
  data <- sync_dev_map_data(input_list$data, input_list$map)
  obj <- RTMB::MakeADFun(cmb(model, data), input_list$par, map = input_list$map, silent = TRUE)
  obj$fn(obj$par)
}

test_that("an sd-only dsem at sigmaR reproduces SPoRC's penalty on dusky", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  set.seed(20)
  input_list$par$ln_RecDevs[] <- rnorm(length(input_list$par$ln_RecDevs), 0, 0.6)
  env <- data.frame(year = input_list$data$years, env = rnorm(length(input_list$data$years)))

  sigmaR <- exp(input_list$par$ln_sigmaR[2,1,1])
  dsem_list <- Setup_Mod_DSEM(input_list,
                              dsem_arrows = c(sprintf("rec <-> rec, 0, NA, %.17g", sigmaR), "env <-> env, 0, NA, 1"),
                              dsem_data = env,
                              dsem_mu_spec = "fix")

  env_nLL <- -sum(dnorm(env$env, mean(env$env), 1, log = TRUE)) # the unlinked covariate's own density
  expect_equal(objective_value(SPoRC_rtmb, dsem_list) - env_nLL, objective_value(SPoRC_rtmb, input_list), tolerance = 1e-8)

})

test_that("linked regions read their own deviations and unlinked regions keep SPoRC's penalty", {

  input_list <- sweep_input() # three regions, distinct extents for every dim
  set.seed(21)
  input_list$par$ln_RecDevs[] <- rnorm(length(input_list$par$ln_RecDevs), 0, 0.6)
  n_yrs <- length(input_list$data$years)
  cov_data <- data.frame(year = input_list$data$years, sst = rnorm(n_yrs), prey = rnorm(n_yrs))
  cov_data$prey[c(2, 5, 9)] <- NA # gaps are latent states

  arrows <- "
    sst -> rec_Pop_1_Region_3, 0, b_sst
    prey -> rec_Pop_1_Region_1, 1, b_prey
    sst -> prey, 0, b_sst_prey
    sst -> sst, 1, rho_sst
    sst <-> sst, 0, sd_sst
    prey <-> prey, 0, sd_prey
    rec_Pop_1_Region_1 <-> rec_Pop_1_Region_1, 0, sd_r1
    rec_Pop_1_Region_3 <-> rec_Pop_1_Region_3, 0, sd_r3
  "
  dsem_list <- Setup_Mod_DSEM(input_list, dsem_arrows = arrows, dsem_data = cov_data)
  dsem_list$par$dsem_beta <- c(0.4, -0.3, 0.2, 0.5)
  dsem_list$par$ln_dsem_sd <- log(c(0.9, 0.7, 0.6, 0.8))
  dsem_list$par$dsem_x[is.na(cov_data$prey), 2] <- c(0.1, -0.2, 0.3)
  dsem_list$par$dsem_mu <- c(0.05, -0.1, 0, 0) # one per series now, linked ones fixed at zero

  # SPoRC with its recruitment penalty off for regions 1 and 3 only
  sporc_list <- input_list
  sporc_list$data$Wt_Rec <- array(1, dim = dim(input_list$par$ln_RecDevs))
  sporc_list$data$Wt_Rec[,c(1, 3),] <- 0

  # grid by hand: dsem_x already spans every series, so only the linked columns need filling
  x_grid <- dsem_list$par$dsem_x
  x_grid[,3] <- input_list$par$ln_RecDevs[1,1,]
  x_grid[,4] <- input_list$par$ln_RecDevs[1,3,]
  mu_grid <- matrix(dsem_list$par$dsem_mu, n_yrs, 4, byrow = TRUE)
  expect_equal(dsem_list$data$dsem_model$variables, c("sst", "prey", "rec_Pop_1_Region_1", "rec_Pop_1_Region_3"))

  # the sweep model runs the full correction, so each linked recruitment cell is centered at minus half its variance
  # given the known cells (the observed covariate years), here by the dense Schur complement rather than the precision
  margvar <- dense_dsem_margvar(dsem_list$par$dsem_beta, dsem_list$par$ln_dsem_sd, x_grid, dsem_list$data$dsem_model, as.vector(dsem_list$data$dsem_x_known))
  mu_grid[,3:4] <- mu_grid[,3:4] - 0.5 * margvar[,3:4]
  by_hand <- dense_dsem_nLL(dsem_list$par$dsem_beta, dsem_list$par$ln_dsem_sd, x_grid, mu_grid, dsem_list$data$dsem_model)

  expect_equal(objective_value(SPoRC_rtmb, dsem_list), objective_value(SPoRC_rtmb, sporc_list) + by_hand, tolerance = 1e-8)

})

test_that("covariates measured with error add their normal likelihood", {

  input_list <- sweep_input()
  set.seed(22)
  n_yrs <- length(input_list$data$years)
  cov_data <- data.frame(year = input_list$data$years, sst = rnorm(n_yrs))
  arrows <- c("sst -> rec_Pop_1_Region_2, 0, b", "sst <-> sst, 0, sd_sst", "rec_Pop_1_Region_2 <-> rec_Pop_1_Region_2, 0, sd_r")
  dsem_list <- Setup_Mod_DSEM(input_list, dsem_arrows = arrows, dsem_data = cov_data, dsem_family = c(sst = "normal"))

  map_x <- matrix(as.integer(dsem_list$map$dsem_x), nrow = n_yrs)
  expect_true(all(!is.na(map_x[,1]))) # every covariate state estimated under measurement error
  expect_true(all(is.na(map_x[,2])))  # the linked column comes from ln_RecDevs instead
  dsem_list$par$dsem_x[,1] <- cov_data$sst + 0.1
  dsem_list$par$ln_dsem_obs_sd <- log(0.3)

  data <- sync_dev_map_data(dsem_list$data, dsem_list$map)
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), dsem_list$par, map = dsem_list$map, silent = TRUE)
  obj$fn(obj$par)
  expect_equal(obj$report()$dsem_obs_nLL, -sum(dnorm(cov_data$sst, cov_data$sst + 0.1, 0.3, log = TRUE)), tolerance = 1e-10)

})

test_that("setup refuses linked deviations it cannot handle", {

  input_list <- sweep_input()
  n_yrs <- length(input_list$data$years)
  cov_data <- data.frame(year = input_list$data$years, sst = rnorm(n_yrs))
  arrows <- c("sst -> rec_Pop_1_Region_1, 0, b", "sst <-> sst, 0, s1", "rec_Pop_1_Region_1 <-> rec_Pop_1_Region_1, 0, s2")

  expect_error(Setup_Mod_DSEM(input_list, c("sst -> rec, 0, b", "sst <-> sst, 0, s1", "rec <-> rec, 0, s2"), cov_data), "No deviation series")
  expect_error(Setup_Mod_DSEM(input_list, arrows, data.frame(year = 1800, sst = 1)), "years should fall")

  shared <- input_list
  shared$map$ln_RecDevs <- factor(rep(1, length(input_list$par$ln_RecDevs))) # one deviation shared everywhere
  expect_error(Setup_Mod_DSEM(shared, arrows, cov_data), "share map levels")

})

test_that("a linked recruitment cell takes the correction its penalty would, and enters an index as it is", {

  # the sweep model runs the full correction. an sd-only link at sigmaR is then the same density as the
  # iid penalty, center included, so the objective does not move. a recruitment deviation index adds the
  # penalty's center back to a cell under it, and reads a linked cell as it is
  plain <- suppressMessages(sweep_input(dims = list(n_regions = 1)))
  n_yrs <- length(plain$data$years)
  sigmaR <- exp(plain$par$ln_sigmaR[2,1,1])
  devs <- seq(-0.3, 0.3, length.out = n_yrs)
  value_of <- function(il, srv_idx_type = NULL) {
    if(!is.null(srv_idx_type)) il$data$srv_idx_type[1] <- srv_idx_type # the anomaly is only computed for a recruitment index
    il$par$ln_RecDevs[] <- devs
    obj <- suppressWarnings(fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)) # a negative anomaly makes the index likelihood NaN, not what is read here
    suppressWarnings(list(fn = obj$fn(obj$par), anom = as.numeric(obj$report()$RecDev_anom), margvar = obj$report()$dsem_margvar_grid))
  }
  sd_only <- suppressMessages(Setup_Mod_DSEM(plain, sprintf("rec <-> rec, 0, NA, %.17g", sigmaR), NULL, dsem_processes = "rec"))
  expect_equal(value_of(sd_only)$fn, value_of(plain)$fn, tolerance = 1e-8)

  # with a covariate observed every year, a linked cell's variance given it is the sd line squared
  env <- data.frame(year = plain$data$years, env = rnorm(n_yrs))
  linked <- suppressMessages(Setup_Mod_DSEM(plain, c("env -> rec, 0, b", "env <-> env, 0, s", "rec <-> rec, 0, sr"), env, dsem_processes = "rec", dsem_mu_spec = "fix"))
  sr <- exp(linked$par$ln_dsem_sd[match("sr", linked$data$dsem_model$ln_sd_names)])
  out <- value_of(linked, srv_idx_type = 2)
  expect_equal(as.numeric(out$margvar[,2]), rep(sr^2, n_yrs), tolerance = 1e-10)
  expect_equal(out$anom, devs + sigmaR^2 / 2) # without the declaration sigmaR is the module's own
  expect_equal(value_of(plain, srv_idx_type = 2)$anom, devs + sigmaR^2 / 2)

  # declared, sigmaR is the arrows' sd line, so the anomaly adds that back
  declared <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"), dims = list(n_regions = 1)))
  declared <- suppressMessages(Setup_Mod_DSEM(declared, "rec <-> rec, 0, NA, 0.7", NULL))
  expect_equal(value_of(declared, srv_idx_type = 2)$anom, devs + 0.7^2 / 2)

  # a ramp at zero is SPoRC's way of taking no correction, and the linked cells then take none either
  plain_none <- suppressMessages(sweep_input(rec = list(do_rec_bias_ramp = 1, bias_year = rep(999, 4)), dims = list(n_regions = 1)))
  sd_none <- suppressMessages(Setup_Mod_DSEM(plain_none, sprintf("rec <-> rec, 0, NA, %.17g", sigmaR), NULL, dsem_processes = "rec"))
  expect_equal(value_of(sd_none)$fn, value_of(plain_none)$fn, tolerance = 1e-8)
  expect_equal(value_of(sd_none, srv_idx_type = 2)$anom, devs)
  expect_true(all(value_of(sd_none)$margvar == 0)) # not worked out when nothing takes it

})
