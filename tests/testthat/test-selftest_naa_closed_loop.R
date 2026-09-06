# A closed loop run end to end with the state-space numbers at age turned on: an assessment is fitted,
# a closed loop conditioned from it, and a fresh estimation model fitted at each assessment year.
#
# This guards routing rather than statistics: the closed loop never calls Simulate_Pop_Static, and a
# conditioning step that drops the state projects a deterministic future rather than failing.

naacl <- local({
  cached <- NULL
  function() {
    if(!is.null(cached)) return(cached)
    proj_yrs <- 8; assess_every <- 4

    om  <- naaom_make_om(NAA_re = "2dar1", sigmaNAA = 0.30, rho_age = 0.5, rho_year = 0.4, seed = 808)
    il  <- naaom_build_em(naaom_om_data(om), NAA_re = "2dar1")
    fit <- suppressWarnings(fit_model(il$data, il$par, il$map, random = "ln_NAA", silent = TRUE))
    sdr <- suppressWarnings(RTMB::sdreport(fit))

    sim_list <- condition_closed_loop_simulations(
      closed_loop_yrs = proj_yrs,
      n_sims = 1,
      data = il$data,
      parameters = il$par,
      mapping = il$map,
      sd_rep = sdr,
      rep = fit$rep,
      random = "ln_NAA"
    )

    set.seed(77)
    sim_env <- Setup_sim_env(sim_list)
    start <- sim_env$feedback_start_yr
    assess_yrs <- seq(start, sim_env$n_yrs, assess_every)
    f_hist <- sim_env$Fmort[1, start, 1, 1, 1]
    rows <- list()

    for(sim in 1:sim_env$n_sims) {
      for(y in 1:sim_env$n_yrs) {
        run_annual_cycle(y = y, sim = sim, sim_env = sim_env)

        if(y >= start && y %in% assess_yrs) {
          il_y <- naaom_build_em(simulation_data_to_SPoRC(sim_env, y, sim), NAA_re = "2dar1")
          obj <- suppressWarnings(fit_model(il_y$data, il_y$par, il_y$map,
                                            random = "ln_NAA", silent = TRUE))
          rp <- get_closed_loop_reference_points(
            use_true_values = FALSE,
            sim_env = sim_env,
            asmt_data = il_y$data,
            asmt_rep = obj$rep,
            y = y,
            sim = sim,
            n_proj_yrs = assess_every + 1,
            reference_points_opt = list(
              n_avg_yrs = 1,
              SPR_x = 0.4,
              calc_rec_st_yr = 1,
              rec_age = 1,
              type = "single_region",
              what = "SPR"
            )
          )

          # a threshold rule on the assessment's own terminal spawning biomass
          ssb_hat <- sum(obj$rep$SSB[, , y])
          ratio <- ssb_hat / rp$b_ref_pt
          f_next <- as.numeric(rp$f_ref_pt)[1] * max(0, min(1, (ratio - 0.05) / 0.95))

          rows[[length(rows) + 1]] <- data.frame(
            y = y,
            gradient = max(abs(obj$gr(naaom_fixed(obj)))),
            n_states = length(obj$env$random),
            sigmaNAA = exp(naaom_fixed(obj, "ln_sigmaNAA")),
            ssb_hat = ssb_hat,
            ssb_true = sum(sim_env$SSB[, , y, sim]),
            f_set = f_next
          )

          if(y < sim_env$n_yrs)
            sim_env$Fmort[, (y + 1):min(y + assess_every, sim_env$n_yrs), , , sim] <- f_next
        } # end assessment year
      } # end y loop
    } # end sim loop

    cached <<- list(
      sim_env = sim_env,
      log = do.call(rbind, rows),
      start = start,
      f_hist = f_hist,
      conditioned_sd = sim_list$sigmaNAA[1, 1, start, 1, 1, 1],
      assess_yrs = assess_yrs
    )
    cached
  }
})

test_that("the closed loop runs to completion with the state on", {
  d <- naacl()
  expect_equal(nrow(d$log), length(d$assess_yrs))
  # the operating model advanced through every projection year rather than stopping at the data
  expect_true(all(d$sim_env$SSB[, , 1:d$sim_env$n_yrs, 1] > 0))
  expect_true(all(is.finite(d$sim_env$NAA[, , 1:d$sim_env$n_yrs, , , , 1])))
})

test_that("each assessment inside the loop estimates the state and converges", {
  d <- naacl()
  expect_true(all(d$log$gradient < 1e-2))
  # the states are integrated out at every assessment, over a window that grows with the loop
  expect_true(all(d$log$n_states > 0))
  expect_true(all(diff(d$log$n_states) > 0))
  # and the process error is recovered near the value the loop was conditioned on
  expect_true(all(abs(d$log$sigmaNAA - d$conditioned_sd) / d$conditioned_sd < 0.5))
})

test_that("the operating model has process error through the projection years", {
  d <- naacl()
  eta <- d$sim_env$naa_eta_all[1, 1, , 1, , 1, 1]
  ages <- d$sim_env$naa_re_ages
  proj <- eta[(d$start + 1):d$sim_env$n_yrs, ages]
  expect_gt(stats::sd(proj), 0)
  expect_equal(stats::sd(proj), stats::sd(eta[2:d$start, ages]), tolerance = 0.5)
  # the recruitment age and the first year stay deterministic here as they do in the past
  expect_equal(max(abs(eta[, 1])), 0)
})

test_that("the assessment advice feeds back into the operating model", {
  d <- naacl()
  # a loop that never closed would leave the projection years at the terminal historical rate
  f_proj <- d$sim_env$Fmort[1, (d$start + 1):d$sim_env$n_yrs, 1, 1, 1]
  expect_true(all(f_proj > 0))
  expect_false(isTRUE(all.equal(as.vector(f_proj), rep(d$f_hist, length(f_proj)))))
  expect_true(any(abs(diff(f_proj)) > 0)) # the rule moved with the stock, not once at the start
})

test_that("the assessment tracks the operating model it is assessing", {
  d <- naacl()
  # loose: one realization of a state-space fit to a stock under active management. What would
  # fail here is a routing break, where the assessment and the operating model are different stocks
  rel <- (d$log$ssb_hat - d$log$ssb_true) / d$log$ssb_true
  expect_true(all(abs(rel) < 0.75))
  expect_gt(stats::cor(d$log$ssb_hat, d$log$ssb_true), 0.9)
})
