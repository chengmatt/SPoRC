# Checks the package diagnostics run on a dsem model: profile (which needed a fix for plain vector
# parameters), jitter, Francis reweighting, OSA residuals, and the DAG plot.

data("sgl_rg_dusky_data")

dsem_diag_setup <- function() {
  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  set.seed(11)
  env <- data.frame(year = input_list$data$years, env = rnorm(length(input_list$data$years)))
  Setup_Mod_DSEM(input_list, c("env -> rec, 0, b_env", "env -> env, 1, rho", "env <-> env, 0, sd_env", "rec <-> rec, 0, NA, 0.9"), env,
                 dsem_mu_spec = "fix")
}
dsem_random <- c("ln_RecDevs", "dsem_x")

test_that("a likelihood profile over a dsem coefficient runs and bowls", {

  d <- dsem_diag_setup()
  prof <- do_likelihood_profile(d$data, d$par, d$map, random = dsem_random, what = "dsem_beta", idx = 1,
                                min_val = -0.2, max_val = 0.2, inc = 0.2, do_par = FALSE)
  jn <- prof$jnLL_df$value
  expect_length(jn, 3)
  expect_true(all(is.finite(jn)))
  expect_lt(jn[2], jn[1]) # the fit's b_env is near 0.1, so the middle point beats -0.2
  expect_lt(jn[2], jn[3] + 1) # and is close to 0.2

})

test_that("Francis reweighting runs on a dsem model", {

  d <- dsem_diag_setup()
  fr <- run_francis(d$data, d$par, d$map, random = dsem_random, n_francis_iter = 1, newton_loops = 0)
  expect_true(!is.null(fr))

})

test_that("jitter runs on a dsem model, and a run's returned start, passed back with sd = 0, repeats that run", {

  d <- dsem_diag_setup()
  set.seed(5)
  jit <- do_jitter(d$data, d$par, d$map, random = dsem_random, sd = 0.1, n_jitter = 2, n_newton_loops = 0,
                   do_par = FALSE, jitter_random = TRUE, return_start_values = TRUE)
  start <- jit$start_values
  expect_equal(nrow(start), 2)
  expect_true(any(colnames(start) %in% dsem_random)) # random effects sit in their own columns
  expect_false(isTRUE(all.equal(start[1, ], start[2, ]))) # each run drew its own start

  rerun <- do_jitter(d$data, d$par, d$map, random = dsem_random, sd = 0, n_jitter = 1, n_newton_loops = 0,
                     do_par = FALSE, par_vec = start[2, ], return_start_values = TRUE)
  expect_equal(rerun$start_values[1, ], start[2, ])
  expect_equal(unique(rerun$jitter_df$jnLL), unique(jit$jitter_df$jnLL[jit$jitter_df$jitter == 2]))

})

test_that("OSA residuals and the DAG plot run on a fitted dsem model", {

  d <- dsem_diag_setup()
  fit <- fit_model(d$data, d$par, d$map, random = dsem_random, newton_loops = 0, silent = TRUE)
  osa <- get_osa(model = fit, data = fit$data, index_source = "SrvIdx")
  expect_true(!is.null(osa))

  skip_if_not_installed("igraph")
  grDevices::pdf(NULL)
  g <- plot_dsem_dag(fit)
  grDevices::dev.off()
  e <- igraph::as_data_frame(g)
  expect_equal(igraph::vcount(g), 3L) # env, rec and lag(env,1), the lagged source as its own node like dsem
  expect_equal(e$lag[e$from == "lag(env,1)" & e$to == "env"], 1L) # the ar1 on env
  expect_true(is.finite(e$estimate[e$from == "env" & e$to == "rec"])) # the env -> rec path, labelled with its estimate
  expect_equal(e$estimate[e$from == "rec" & e$to == "rec"], 0.9) # the fixed sd shows its value
  expect_true(all(nzchar(e$label)))

})
