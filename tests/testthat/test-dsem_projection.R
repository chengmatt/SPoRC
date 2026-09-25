# Checks that n_proj_yrs_devs is the projection: the fitted model's latent projected cells have the
# conditional mean (1e-6) and sd (1e-3) of the dense conditional distribution given the fitted years.

data("sgl_rg_dusky_data")

test_that("the fit's projected deviation cells are the dsem's conditional forecast", {

  build_ext <- build_goa_dusky_input
  body(build_ext) <- do.call(substitute, list(body(build_ext), list(Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 5)))))
  input_list <- build_ext(sgl_rg_dusky_data)
  n_yrs <- length(input_list$data$years); n_proj <- 5
  expect_equal(dim(input_list$par$ln_RecDevs)[3], n_yrs + n_proj)

  set.seed(5)
  env <- data.frame(year = input_list$data$years, env = rnorm(n_yrs)) # observed in model years only
  d <- Setup_Mod_DSEM(input_list, c("env -> rec, 0, NA, 0.4", "env -> env, 1, NA, 0.6", "env <-> env, 0, NA, 1", "rec <-> rec, 0, NA, 0.9"), env,
                      dsem_mu_spec = "fix")
  expect_equal(d$data$dsem_n_grid_yrs, n_yrs + n_proj)

  fit <- fit_model(d$data, d$par, d$map, random = c("ln_RecDevs", "dsem_x"), newton_loops = 1, silent = TRUE)
  sdr <- RTMB::sdreport(fit)
  pars <- fit$env$parList()
  proj <- (n_yrs + 1):(n_yrs + n_proj)

  # the dense conditional given every fitted cell, the same thing the operating model draws from
  n_grid <- d$data$dsem_n_grid_yrs
  x_grid <- cbind(pars$dsem_x[,1], pars$ln_RecDevs[1,1,])
  mu_grid <- cbind(rep(pars$dsem_mu[1], n_grid), 0)
  known <- rep(c(TRUE, TRUE), each = n_grid); known[c(proj, n_grid + proj)] <- FALSE
  Q <- get_dsem_precision(pars$dsem_beta, pars$ln_dsem_sd, d$data$dsem_model, get_dsem_cells(d$data$dsem_model, n_grid))
  cond <- get_dsem_conditional(Q, as.vector(mu_grid), which(known), as.vector(x_grid)[known])
  cond_sd <- sqrt(Matrix::diag(Matrix::solve(cond$chol_uu, Matrix::Diagonal(length(cond$unknown_cell)), system = "A")))
  rec_cells <- which(cond$unknown_cell > n_grid)

  # the fit's own answer for those cells: their modes and Laplace standard errors
  sr <- summary(sdr, "random")
  rec_rows <- which(rownames(sr) == "ln_RecDevs")
  expect_equal(unname(pars$ln_RecDevs[1,1,proj]), unname(cond$cond_mean[rec_cells]), tolerance = 1e-6)
  expect_equal(unname(sr[rec_rows[proj], "Std. Error"]), unname(cond_sd[rec_cells]), tolerance = 1e-3)

  # and the forecast is the arrow: rec = b_env * (env - mu_env) in years with nothing else to read
  expect_equal(unname(pars$ln_RecDevs[1,1,proj]), 0.4 * (pars$dsem_x[proj,1] - pars$dsem_mu[1]), tolerance = 1e-6)

})
