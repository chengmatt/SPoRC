# Initial age deviations under RecDevs_model = 'dsem' read the recruitment series' settled marginal
# variance, not the sd line. Checked against the dense covariance oracle in helper-dsem.R at 1e-8.

data("dusky_rtmb_model")

# dusky with recruitment under the arrows, one self path and one sd line at known values
init_margvar_input <- function(variance = "conditional", rec_dsem = TRUE, rho = 0.5, sd_line = 0.4) {

  base <- list(data = dusky_rtmb_model$data, par = dusky_rtmb_model$parameters,
               map = dusky_rtmb_model$mapping, verbose = FALSE, store_config = FALSE)

  rec_in <- suppressMessages(suppressWarnings(Setup_Mod_Rec(
    base,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(sd_line), dim = c(2, base$data$n_pop, base$data$n_regions)),
    rec_model = "mean_rec",
    init_age_strc = 1,
    ln_global_R0 = log(2.7),
    t_spawn = base$data$t_spawn,
    RecDevs_model = if(rec_dsem) "dsem" else "iid",
    sigmaR_spec = "fix" # once every rec cell is linked the penalty stops reading it, so it cannot stay estimated
  )))

  arrows <- c(sprintf("rec -> rec, 1, rho, %.17g", rho), sprintf("rec <-> rec, 0, sd_rec, %.17g", sd_line))
  suppressMessages(Setup_Mod_DSEM(rec_in, dsem_arrows = arrows, dsem_data = NULL, dsem_variance = variance))

}

# the report at the starting values, and the pieces the init deviation penalty is read from
init_margvar_rep <- function(il) {

  dat <- sync_dev_map_data(il$data, il$map)
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, dat), il$par, map = il$map, silent = TRUE)
  rep <- obj$report(obj$par)

  n_grid <- nrow(rep$dsem_margvar_grid)
  last_row <- max(dat$dsem_link_row[[1]])
  cell <- (dat$dsem_link_col[1] - 1) * n_grid + last_row
  idx <- 1:(length(dat$ages) - 2) # init_age_strc = 1 leaves the plus group out

  list(rep = rep, n_grid = n_grid, cell = cell,
       settled = as.numeric(rep$dsem_margvar_grid)[cell],
       dev = as.numeric(il$par$ln_InitDevs)[idx],
       nLL = as.numeric(rep$Init_Rec_nLL)[idx])

}

init_nLL_at <- function(dev, sigma) -stats::dnorm(dev, -sigma^2 / 2, sigma, log = TRUE)

test_that("the conditional form gives the init devs the settled variance, not the sd line", {

  rho <- 0.5
  sd_line <- 0.4
  il <- init_margvar_input("conditional", rho = rho, sd_line = sd_line)
  out <- init_margvar_rep(il)

  # with only a self path the settled variance is the stationary one, and the dense oracle agrees
  oracle <- dense_dsem_margvar(il$par$dsem_beta, il$par$ln_dsem_sd, matrix(0, out$n_grid, 1),
                               il$data$dsem_model, rep(FALSE, out$n_grid))
  expect_equal(out$settled, sd_line^2 / (1 - rho^2), tolerance = 1e-6)
  expect_equal(out$settled, as.numeric(oracle)[out$cell], tolerance = 1e-8)

  # the penalty sits at that spread and its own bias-corrected center
  expect_equal(out$nLL, init_nLL_at(out$dev, sqrt(out$settled)), tolerance = 1e-10)
  expect_gt(max(abs(out$nLL - init_nLL_at(out$dev, sd_line))), 1e-3)

})

test_that("the marginal form leaves the init devs at the sd line", {

  sd_line <- 0.4
  il <- init_margvar_input("marginal", sd_line = sd_line)
  out <- init_margvar_rep(il)

  expect_equal(out$settled, sd_line^2, tolerance = 1e-8) # every cell is held at the sd line squared
  expect_equal(out$nLL, init_nLL_at(out$dev, sd_line), tolerance = 1e-10)

})

test_that("a rec link without the declaration leaves the init devs on ln_sigmaR", {

  sd_line <- 0.4
  em_sigma <- 0.25 # the model's own sigmaR, which the arrows do not stand in for here
  il <- init_margvar_input("conditional", rec_dsem = FALSE, sd_line = sd_line)
  il$par$ln_sigmaR[] <- log(em_sigma)
  out <- init_margvar_rep(il)

  expect_length(il$data$dsem_declared, 0)
  expect_equal(out$nLL, init_nLL_at(out$dev, em_sigma), tolerance = 1e-10)

})

test_that("the operating model reads the same settled sd as the fit", {

  rho <- 0.5
  sd_line <- 0.4
  il <- init_margvar_input("conditional", rho = rho, sd_line = sd_line)
  out <- init_margvar_rep(il)

  obj <- suppressMessages(suppressWarnings(fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)))
  sim_list <- suppressMessages(suppressWarnings(condition_closed_loop_simulations(
    closed_loop_yrs = 2, n_sims = 1, data = obj$data, parameters = obj$parameters, mapping = obj$mapping,
    sd_rep = list(par.fixed = obj$par, par.random = NULL), rep = obj$rep, random = NULL)))
  sim_list <- suppressMessages(Setup_Sim_DSEM(sim_list, obj$data, obj$env$parList(), rep = obj$rep))

  expect_equal(as.numeric(exp(sim_list$ln_sigmaR[1,1,1])), sqrt(out$settled), tolerance = 1e-6)
  expect_equal(as.numeric(exp(sim_list$ln_sigmaR[2,1,1])), sqrt(out$settled), tolerance = 1e-6)

})
