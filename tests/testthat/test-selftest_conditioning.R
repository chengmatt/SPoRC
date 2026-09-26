# A conditioning period reproduces the fit: the self test over all its years, the closed loop over the
# fitted ones with the projection still drawn. What conditioning does to a drawn dsem series is checked
# in test-dsem_sim_derived.R.

selftest_dsem_fit <- function(n_yrs = 40) {
  om <- q_devs_sim("none", n_yrs = n_yrs, seed = 20)
  il <- blocks_em(om)
  set.seed(5)
  env <- data.frame(year = 1:n_yrs, env = as.numeric(scale(stats::rnorm(n_yrs))))
  arrows <- c("env -> rec, 0, b_env", "env <-> env, 0, NA, 1", "rec <-> rec, 0, sd_rec")
  d <- suppressMessages(Setup_Mod_DSEM(il, arrows, env, dsem_processes = "rec",
                                       dsem_family = c(env = "fixed"), dsem_mu_spec = c(env = 0)))
  suppressWarnings(suppressMessages(fit_model(d$data, d$par, d$map, random = "ln_RecDevs",
                                              do_optim = TRUE, newton_loops = 1, silent = TRUE)))
}

test_that("the self test asks its dsem to condition on the fit", {

  fit <- selftest_dsem_fit()
  sd_rep <- list(par.fixed = fit$optim$par, par.random = fit$env$last.par.best[fit$env$random])
  seen <- new.env()

  # local_mocked_bindings reaches the namespace binding the self test calls; trace() does not, because
  # it leaves the copy inside the package untouched and the call goes to the original
  local_mocked_bindings(
    Setup_Sim_DSEM = function(sim_list, data = NULL, pars = NULL, rep = NULL, condition_on_fit = FALSE, ...) {
      assign("cond", condition_on_fit, envir = seen)
      stop("probe") # stop here, so none of the refits run
    },
    .package = "SPoRC"
  )

  invisible(try(suppressWarnings(suppressMessages(
    simulation_self_test(data = fit$data, parameters = fit$env$parList(), mapping = fit$mapping,
                         random = "ln_RecDevs", rep = fit$rep, sd_rep = sd_rep, n_sims = 1,
                         newton_loops = 1, what = "SSB"))), silent = TRUE))

  expect_true(exists("cond", envir = seen))          # the branch was reached at all
  expect_true(get("cond", envir = seen))             # and it asked for the fitted years to be kept
})

test_that("a state space fit is reproduced by the self test's operating model", {

  # data generated with a state, so the fit has innovations worth reproducing
  om <- suppressWarnings(suppressMessages(naaom_make_om()))
  il <- suppressWarnings(suppressMessages(naaom_build_em(naaom_om_data(om), NAA_re = "2dar1")))
  fit <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, random = "ln_NAA",
                                                     do_optim = TRUE, newton_loops = 1, silent = TRUE)))
  # parList() reads last.par, the last point evaluated; the report below is at last.par.best
  pars <- fit$env$parList(par = fit$env$last.par.best)
  n_yrs <- length(fit$data$years)

  eta <- as.numeric(pars$ln_NAA - log(fit$rep$NAA_pred[,,seq_len(dim(pars$ln_NAA)[3]),,,,drop = FALSE]))
  expect_gt(stats::sd(eta[is.finite(eta)]), 0.05) # a state worth reproducing, not one collapsed to zero

  # the sim list the self test builds, caught before the refits it would then run. the mock is scoped to
  # this block so it has unwound by the time the operating model is run below for real
  seen <- new.env()
  local({
    local_mocked_bindings(
      Simulate_Pop_Static = function(sim_list, ...) { assign("sl", sim_list, envir = seen); stop("probe") },
      .package = "SPoRC"
    )
    invisible(try(suppressWarnings(suppressMessages(
      simulation_self_test(data = fit$data, parameters = pars, mapping = fit$mapping, random = "ln_NAA",
                           rep = fit$rep, sd_rep = list(par.fixed = fit$optim$par, par.random = NULL),
                           n_sims = 1, newton_loops = 1, what = "SSB"))), silent = TRUE))
  })
  sl <- get("sl", envir = seen)

  expect_true(isTRUE(sl$NAA_re > 0))          # the state is carried at all
  expect_false(is.null(sl$naa_eta_input))     # and at the fit's own innovations

  om_run <- suppressWarnings(suppressMessages(Simulate_Pop_Static(sim_list = sl, output_path = NULL)))
  ssb_fit <- as.numeric(fit$rep$SSB[1, 1, seq_len(n_yrs)])
  ssb_om <- as.numeric(om_run$SSB[1, 1, seq_len(n_yrs), 1])

  # without the state the operating model's biomass ran 33% away from the fit on average
  # 1e-6 rather than machine precision, since the states come out of the inner Laplace solve
  expect_equal(ssb_om, ssb_fit, tolerance = 1e-6)
})

test_that("a state space self test runs end to end and recovers", {

  om <- suppressWarnings(suppressMessages(naaom_make_om()))
  il <- suppressWarnings(suppressMessages(naaom_build_em(naaom_om_data(om), NAA_re = "2dar1")))
  fit <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, random = "ln_NAA",
                                                     do_optim = TRUE, newton_loops = 1, silent = TRUE)))
  n_yrs <- length(fit$data$years)

  res <- suppressWarnings(suppressMessages(
    simulation_self_test(data = fit$data, parameters = fit$env$parList(), mapping = fit$mapping,
                         random = "ln_NAA", rep = fit$rep,
                         sd_rep = list(par.fixed = fit$optim$par, par.random = NULL),
                         n_sims = 3, newton_loops = 1, what = "SSB")))

  ssb_true <- as.numeric(fit$rep$SSB[1, 1, seq_len(n_yrs)])
  est <- res[[1]]
  expect_equal(dim(est)[3], n_yrs)

  for(i in 1:3) {
    fitted_ssb <- as.numeric(est[1, 1, seq_len(n_yrs), i])
    expect_true(all(is.finite(fitted_ssb)))                                    # the replicate refit at all
    # what is left is observation error, now that the truth is the fitted population rather than
    # the deterministic one, which used to sit 33% away from it
    expect_lt(stats::median(abs(fitted_ssb - ssb_true) / ssb_true), 0.2)
  } # end i loop
})

test_that("a closed loop reproduces a state space fit over its conditioning years and draws past them", {

  om <- suppressWarnings(suppressMessages(naaom_make_om()))
  il <- suppressWarnings(suppressMessages(naaom_build_em(naaom_om_data(om), NAA_re = "2dar1")))
  fit <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, random = "ln_NAA",
                                                     do_optim = TRUE, newton_loops = 1, silent = TRUE)))
  n_yrs <- length(fit$data$years)

  sl <- suppressWarnings(suppressMessages(
    condition_closed_loop_simulations(closed_loop_yrs = 5, n_sims = 2, data = fit$data,
                                      parameters = fit$env$parList(par = fit$env$last.par.best),
                                      mapping = fit$mapping,
                                      sd_rep = list(par.fixed = fit$optim$par,
                                                    par.random = fit$env$last.par.best[fit$env$random]),
                                      rep = fit$rep, random = "ln_NAA")))

  expect_equal(sl$n_cond_yrs, n_yrs)
  expect_false(is.null(sl$naa_eta_input))
  expect_equal(dim(sl$naa_eta_input)[3], n_yrs) # the fitted years only, so the projection still draws

  res <- suppressWarnings(suppressMessages(Simulate_Pop_Static(sim_list = sl, output_path = NULL)))
  ssb_fit <- as.numeric(fit$rep$SSB[1, 1, seq_len(n_yrs)])

  for(i in 1:2) expect_equal(as.numeric(res$SSB[1, 1, seq_len(n_yrs), i]), ssb_fit, tolerance = 1e-6)

  # and the years past the fit are the operating model's own, not the fit's
  proj <- (n_yrs + 1):dim(res$SSB)[3]
  expect_gt(max(abs(res$SSB[1, 1, proj, 1] - res$SSB[1, 1, proj, 2]) / res$SSB[1, 1, proj, 1]), 0.01)
})

test_that("conditioning a state space fit without its random effects warns rather than using zeros", {

  om <- suppressWarnings(suppressMessages(naaom_make_om()))
  il <- suppressWarnings(suppressMessages(naaom_build_em(naaom_om_data(om), NAA_re = "2dar1")))
  fit <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, random = "ln_NAA",
                                                     do_optim = TRUE, newton_loops = 1, silent = TRUE)))

  # par.random left out, so the innovations cannot be read; conditioning on the starting zeros would
  # collapse the population while every replicate agreed, which is the worst way to be wrong
  expect_warning(
    suppressMessages(condition_closed_loop_simulations(closed_loop_yrs = 3, n_sims = 1, data = fit$data,
                                                       parameters = fit$env$parList(), mapping = fit$mapping,
                                                       sd_rep = list(par.fixed = fit$optim$par, par.random = NULL),
                                                       rep = fit$rep, random = "ln_NAA")),
    "innovations cannot be read")
})
