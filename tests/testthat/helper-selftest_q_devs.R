# A random walk on survey catchability fitted once and reused, and the simulation list the self test passes
# the operating model, captured before any replicate runs.

q_rw_cache <- new.env(parent = emptyenv())

q_rw_fit <- function() {
  if(is.null(q_rw_cache$fit)) {
    q_rw_cache$fit <- q_devs_fit(q_devs_em(q_devs_sim("rw", sigma_q = 0.25, seed = 20), "rw"))$fit
    q_rw_cache$sd_rep <- RTMB::sdreport(q_rw_cache$fit, getJointPrecision = TRUE)
  }
  q_rw_cache
}

q_selftest_capture <- function(fit, sd_rep, ...) {
  out <- NULL
  testthat::with_mocked_bindings(
    Simulate_Pop_Static = function(sim_list, ...) { out <<- sim_list; stop("captured") },
    try(suppressWarnings(suppressMessages(
      simulation_self_test(data = fit$data, parameters = fit$env$parList(), mapping = fit$mapping,
                           random = "ln_srv_q_devs", rep = fit$rep, sd_rep = sd_rep, obj = fit,
                           newton_loops = 0, what = "SSB", ...))), silent = TRUE),
    .package = "SPoRC")
  out
}
