# The projection inside the model with projected recruitment deviations and state innovations set to known, nonzero
# values that differ by year: each recruitment and state option reproduces Do_Population_Projection given the same
# deviations (1e-10), and projected recruits are mean recruitment times the deviation each option names.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

proj_devs_n_proj_yrs <- 6
proj_devs_dims <- list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_seas = 1, n_yrs = 22, n_ages = 8)

# a state-space numbers at age model with five projected deviation years, evaluated at values set by hand
proj_devs_input <- local({
  cache <- NULL
  function() {
    if(!is.null(cache)) return(cache)

    n_yrs <- proj_devs_dims$n_yrs

    # the sweep model with five projected deviation years
    build <- sweep_input
    body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 5)))))

    # a deterministic pass seeds the state over the fitted years
    il0 <- suppressWarnings(suppressMessages(build(dims = proj_devs_dims)))
    obj0 <- suppressWarnings(suppressMessages(fit_model(il0$data, il0$par, il0$map, do_optim = FALSE, silent = TRUE)))
    seed <- log(array(obj0$report(obj0$par)$NAA[,,1:n_yrs,1,,], dim = c(1, 1, n_yrs, 1, proj_devs_dims$n_ages, 1)))
    il <- suppressWarnings(suppressMessages(build(dims = proj_devs_dims, biol = list(NAA_re = "iid", NAA_sigma_spec = "fix", M_spec = "fix", ln_NAA = seed))))

    # recruitment deviations: the terminal year, then the five projected years, all different
    il$par$ln_RecDevs[1,1,n_yrs] <- 0.35
    il$par$ln_RecDevs[1,1,n_yrs + 1:5] <- c(0.4, -0.3, 0.2, -0.5, 0.1)

    # state innovations in the projected years, growing by year and alternating in sign by age
    state_ages <- il$data$naa_re_ages
    for(k in 1:5) {
      for(i in seq_along(state_ages)) il$par$ln_NAA[1,1,n_yrs + k,1,state_ages[i],1] <- 0.05 * k * (-1)^i
    } # end k loop

    cache <<- il
    cache
  }
})

# the model's report with the projection set up as asked
proj_devs_internal <- function(...) {
  il <- proj_devs_input()
  ip <- suppressWarnings(suppressMessages(Setup_Mod_Projection(il, n_proj_yrs = proj_devs_n_proj_yrs, f_opt = "terminal", ...)))
  suppressWarnings(suppressMessages(fit_model(ip$data, ip$par, ip$map, do_optim = FALSE, silent = TRUE)))
}

# Do_Population_Projection from the same report, given multipliers on recruitment and on numbers at age
proj_devs_stand_alone <- function(fit, rec_devs, naa_devs = NULL) {
  n_yrs <- proj_devs_dims$n_yrs
  terminal_F <- sum(fit$rep$Fmort[1,n_yrs,,])
  inputs <- build_projection_inputs(data = fit$data, rep = fit$rep, refpts = list(f_ref_pt = terminal_F, b_ref_pt = array(1, dim = c(1, 1))),
                                    n_proj_yrs = proj_devs_n_proj_yrs, t_spawn = 0)
  do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "Input", rec_devs = rec_devs, naa_devs = naa_devs)))
}

# the multipliers each option should apply, built from the parameters by hand
proj_devs_multipliers <- function() {

  il <- proj_devs_input()
  n_yrs <- proj_devs_dims$n_yrs
  n_proj_yrs <- proj_devs_n_proj_yrs
  state_ages <- il$data$naa_re_ages

  rec_mean <- rec_propagate <- rec_terminal <- array(1, dim = c(1, 1, n_proj_yrs)) # year 1 gets no new recruits
  naa_propagate <- array(1, dim = c(1, 1, n_proj_yrs, 1, proj_devs_dims$n_ages, 1))
  for(y in 2:n_proj_yrs) {
    rec_propagate[1,1,y] <- exp(il$par$ln_RecDevs[1,1,n_yrs + y - 1]) # that year's own projected deviation
    rec_terminal[1,1,y] <- exp(il$par$ln_RecDevs[1,1,n_yrs]) # the terminal year's deviation kept
    naa_propagate[1,1,y,1,state_ages,1] <- exp(il$par$ln_NAA[1,1,n_yrs + y - 1,1,state_ages,1]) # that year's innovations
  } # end y loop

  list(rec_mean = rec_mean, rec_propagate = rec_propagate, rec_terminal = rec_terminal, naa_propagate = naa_propagate)
}

# Tests ----------------------------------------------------------------------

test_that("each recruitment and state option reproduces the stand-alone projection given the same deviations", {

  mult <- proj_devs_multipliers()

  # projected recruitment deviations and state innovations
  fit <- proj_devs_internal()
  expect_equal(fit$rep$proj_NAA, proj_devs_stand_alone(fit, mult$rec_propagate, mult$naa_propagate)$proj_NAA, tolerance = 1e-10)

  # the terminal recruitment deviation kept, with state innovations
  fit <- proj_devs_internal(rec_opt = "terminal")
  expect_equal(fit$rep$proj_NAA, proj_devs_stand_alone(fit, mult$rec_terminal, mult$naa_propagate)$proj_NAA, tolerance = 1e-10)

  # mean recruitment, with state innovations
  fit <- proj_devs_internal(rec_opt = "mean")
  expect_equal(fit$rep$proj_NAA, proj_devs_stand_alone(fit, mult$rec_mean, mult$naa_propagate)$proj_NAA, tolerance = 1e-10)

  # projected recruitment deviations, with the terminal state advanced and no innovations
  fit <- proj_devs_internal(naa_opt = "terminal")
  expect_equal(fit$rep$proj_NAA, proj_devs_stand_alone(fit, mult$rec_propagate)$proj_NAA, tolerance = 1e-10)

})

test_that("projected recruits are mean recruitment times the deviation each option names, in every year", {

  mult <- proj_devs_multipliers()
  n_proj_yrs <- proj_devs_n_proj_yrs

  fit <- proj_devs_internal()
  mean_rec <- mean(fit$rep$Rec[1,1,fit$data$proj_calc_rec_yrs])
  expect_equal(fit$rep$proj_NAA[1,1,2:n_proj_yrs,1,1,1] / mean_rec, mult$rec_propagate[1,1,2:n_proj_yrs], tolerance = 1e-10)

  fit <- proj_devs_internal(rec_opt = "terminal")
  expect_equal(fit$rep$proj_NAA[1,1,2:n_proj_yrs,1,1,1] / mean_rec, mult$rec_terminal[1,1,2:n_proj_yrs], tolerance = 1e-10)

  fit <- proj_devs_internal(rec_opt = "mean")
  expect_equal(fit$rep$proj_NAA[1,1,2:n_proj_yrs,1,1,1] / mean_rec, rep(1, n_proj_yrs - 1), tolerance = 1e-10)

})

test_that("the deviations and innovations set here change the projection", {

  # if they did not, the comparisons above could not tell one year's deviation from another's
  mult <- proj_devs_multipliers()
  with_state <- proj_devs_internal()$rep$proj_NAA
  without_state <- proj_devs_internal(naa_opt = "terminal")$rep$proj_NAA

  expect_gt(max(abs(with_state / without_state - 1), na.rm = TRUE), 0.1)
  expect_gt(max(abs(mult$rec_propagate / mult$rec_terminal - 1)), 0.3)
  expect_equal(length(unique(as.numeric(mult$rec_propagate[1,1,2:proj_devs_n_proj_yrs]))), proj_devs_n_proj_yrs - 1)

})
