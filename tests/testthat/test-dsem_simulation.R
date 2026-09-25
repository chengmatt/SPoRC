# Checks the dsem operating model: draws against the dense covariance (5 SE), fitted years kept under
# conditioning, recruitment as the deterministic value times exp(deviation) (1e-10), and swapped inputs.

data("sgl_rg_dusky_data")

# dusky at its fixed effects estimates, with a dsem at chosen values standing in for a fit
sim_years <- sgl_rg_dusky_data$years
sim_base_input <- build_goa_dusky_input(sgl_rg_dusky_data)
sim_base_fit <- fit_model(sim_base_input$data, sim_base_input$par, sim_base_input$map, random = NULL, newton_loops = 0, silent = TRUE)
set.seed(40)
sim_cov <- rnorm(length(sim_years))
sim_cov[sim_years < 1990] <- NA

sim_dsem_list <- Setup_Mod_DSEM(sim_base_input,
                                dsem_arrows = c("env -> rec, 0, b_env", "env -> env, 1, rho_env", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec"),
                                dsem_data = data.frame(year = sim_years, env = sim_cov),
                                dsem_family = c(env = "normal"),
                                dsem_mu_spec = "fix")
sim_base_pars <- sim_base_fit$env$parList(par = sim_base_fit$env$last.par.best)
for(par_name in names(sim_base_pars)) sim_dsem_list$par[[par_name]] <- sim_base_pars[[par_name]]
sim_dsem_list$par$dsem_beta <- c(0.5, 0.4)
sim_dsem_list$par$ln_dsem_sd <- log(c(0.9, 1.1))
sim_dsem_list$par$ln_dsem_obs_sd <- log(0.3)

sim_obj <- fit_model(sim_dsem_list$data, sim_dsem_list$par, sim_dsem_list$map, random = c("ln_RecDevs", "dsem_x"),
                     model = SPoRC_rtmb, do_optim = FALSE, silent = TRUE)
sim_best <- sim_obj$env$last.par.best
sim_sdrep <- list(par.fixed = sim_best[-sim_obj$env$random], par.random = sim_best[sim_obj$env$random])

build_sim_env <- function(recruitment_opt = 0, condition_on_fit = FALSE, seed = 41) {
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = 5, n_sims = 3, data = sim_obj$data, parameters = sim_obj$parameters,
                                                mapping = sim_obj$mapping, sd_rep = sim_sdrep, rep = sim_obj$rep, random = sim_obj$random,
                                                Rec_input = NULL, recruitment_opt = recruitment_opt)
  sim_list <- Setup_Sim_DSEM(sim_list, sim_obj$data, sim_obj$env$parList(x = sim_best[-sim_obj$env$random], par = sim_best), condition_on_fit = condition_on_fit)
  sim_env <- Setup_sim_env(sim_list)
  set.seed(seed)
  draw_dsem_sim(sim_env)
  sim_env
}

run_sim <- function(sim_env, sim, seed = 42) {
  set.seed(seed)
  for(y in 1:sim_env$n_yrs) run_annual_cycle(y, sim, sim_env)
  invisible(sim_env)
}

test_that("unconditional draws have the field's mean and covariance", {

  om <- build_sim_env()
  mc_env <- list2env(mget(c("n_pop", "n_regions", "dsem_model", "dsem_beta", "ln_dsem_sd", "dsem_mu", "ln_dsem_obs_sd",
                            "dsem_cov_family", "dsem_cov_var_idx", "dsem_cov_use", "dsem_n_cond_yrs", "dsem_x_fit",
                            "dsem_link_par", "dsem_link_sim_par", "dsem_link_col", "dsem_link_row", "dsem_link_cell", "dsem_link_idx", "dsem_link_yr_dim", "ln_RecDevs"), envir = om))
  mc_env$dsem_drawn <- list(ln_RecDevs = array(FALSE, dim = dim(om$ln_RecDevs)[-4]))
  mc_env$n_yrs <- om$n_yrs
  mc_env$n_sims <- 20000
  set.seed(43)
  draw_dsem_sim(mc_env)

  Q <- get_dsem_precision(mc_env$dsem_beta, mc_env$ln_dsem_sd, mc_env$dsem_model, get_dsem_cells(mc_env$dsem_model, mc_env$n_yrs))
  Sigma <- solve(as.matrix(Q))
  draws <- matrix(mc_env$dsem_x_sim, nrow = mc_env$n_yrs * 2) # [cell, sim]
  mu_cell <- rep(mc_env$dsem_mu, each = mc_env$n_yrs) # dsem_mu spans every series now, linked ones at zero

  mean_z <- (rowMeans(draws) - mu_cell) / sqrt(diag(Sigma) / mc_env$n_sims)
  var_z <- (apply(draws, 1, stats::var) - diag(Sigma)) / (sqrt(2 / mc_env$n_sims) * diag(Sigma))
  expect_lt(max(abs(mean_z)), 5)
  expect_lt(max(abs(var_z)), 5)

  # measurement error on the observed covariate, NA where unobserved
  err <- mc_env$dsem_cov_obs_sim[,1,] - mc_env$dsem_x_sim[,1,]
  expect_identical(which(rowSums(is.na(err)) > 0), which(mc_env$dsem_cov_use[,1] == 0))
  expect_equal(stats::sd(err[!is.na(err)]), 0.3, tolerance = 0.01)

})

test_that("conditioned draws keep every fitted year at the fit", {

  om <- build_sim_env(condition_on_fit = TRUE)
  n_fit <- length(sim_years)
  for(sim in 1:om$n_sims) expect_equal(om$dsem_x_sim[1:n_fit,,sim], om$dsem_x_fit, tolerance = 1e-12)
  expect_gt(stats::sd(om$dsem_x_sim[n_fit + 1,2,]), 0) # later years are drawn

})

test_that("recruitment is R0 times exp(deviation), the draw already about minus half its variance", {

  om <- run_sim(build_sim_env(), sim = 1)
  R0 <- as.numeric(sim_obj$rep$R0)
  expect_equal(om$Rec[1,1,,1], R0 * exp(om$dsem_x_sim[,2,1]), tolerance = 1e-10) # the dsem's rec column, drawn about its correction
  expect_equal(om$ln_RecDevs[1,1,,1], om$dsem_x_sim[,2,1], tolerance = 1e-12) # written straight into the array

})

test_that("under Beverton-Holt a deviation moves its own year's recruitment by its exponential", {

  om_a <- run_sim(build_sim_env(recruitment_opt = 1), sim = 1)
  om_b <- build_sim_env(recruitment_opt = 1)
  om_b$ln_RecDevs[1,1,30,1] <- om_b$ln_RecDevs[1,1,30,1] + 0.5
  run_sim(om_b, sim = 1)

  expect_equal(om_b$Rec[1,1,1:29,1], om_a$Rec[1,1,1:29,1], tolerance = 1e-12)
  expect_equal(om_b$Rec[1,1,30,1], om_a$Rec[1,1,30,1] * exp(0.5), tolerance = 1e-10)
  det_rec <- om_a$Rec[1,1,,1] / exp(om_a$dsem_x_sim[,2,1]) # the stock recruit prediction
  expect_gt(stats::sd(det_rec) / mean(det_rec), 1e-3) # moves with spawning biomass, so the curve is in use

})

test_that("a closed loop operating model keeps fitted recruitment and continues from it", {

  # historical recruitment comes from the fit (the default Rec_input), later years from conditioned draws
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = 5, n_sims = 2, data = sim_obj$data, parameters = sim_obj$parameters,
                                                mapping = sim_obj$mapping, sd_rep = sim_sdrep, rep = sim_obj$rep, random = sim_obj$random)
  sim_list <- Setup_Sim_DSEM(sim_list, sim_obj$data, sim_obj$env$parList(x = sim_best[-sim_obj$env$random], par = sim_best), condition_on_fit = TRUE)
  sim_env <- Setup_sim_env(sim_list) # the annual cycle finds the environment by this name in the calling frame
  set.seed(44)
  draw_dsem_sim(sim_env)
  run_sim(sim_env, sim = 1)

  n_fit <- length(sim_years)
  R0 <- as.numeric(sim_obj$rep$R0)
  expect_equal(sim_env$Rec[1,1,1:n_fit,1], sim_obj$rep$Rec[1,1,], tolerance = 1e-10)
  expect_equal(sim_env$Rec[1,1,n_fit + 1:5,1], R0 * exp(sim_env$dsem_x_sim[n_fit + 1:5,2,1]), tolerance = 1e-10)

})

test_that("an operating model with no dsem draws its own deviations and the fields sit at zero", {

  # the branch generate_recruitment reads has to exist without a dsem, otherwise every plain
  # simulation would fail on it. Setup_sim_env puts it there, all zero
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = 3, n_sims = 2, data = sim_base_input$data,
                                                parameters = sim_base_input$par, mapping = sim_base_input$map,
                                                sd_rep = list(par.fixed = sim_base_fit$env$last.par.best, par.random = NULL),
                                                rep = sim_base_fit$report(sim_base_fit$env$last.par.best), random = NULL,
                                                Rec_input = NULL, recruitment_opt = 0)
  sim_env <- Setup_sim_env(sim_list) # the annual cycle's helpers find the environment by this name
  expect_true(all(!sim_env$dsem_drawn$ln_RecDevs))
  expect_equal(dim(sim_env$dsem_drawn$ln_RecDevs), c(sim_env$n_pop, sim_env$n_regions, sim_env$n_yrs))
  expect_false("ln_NAA" %in% names(sim_env$dsem_drawn)) # the operating model holds innovations, not the log state, so the mask sits on naa_eta_all
  expect_true(all(names(sim_env$dsem_drawn) %in% vapply(dsem_process_table(), function(e) e$sim_par, "")))

  set.seed(5)
  for(y in 1:sim_env$n_yrs) run_annual_cycle(y, 1, sim_env)
  expect_true(any(sim_env$ln_RecDevs[,,,1] != 0)) # drew its own, since nothing is linked

})
