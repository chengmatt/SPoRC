# The numbers at age state run over the n_proj_yrs_devs projection years. Those cells hold innovations with
# no data behind them, so they have to leave the fit unchanged, sit at the process forecast given the fitted
# years, and, passed to Do_Population_Projection as naa_devs, land where the fit applies its own state.

library(SPoRC)
library(testthat)

naaproj_n_proj <- 5

# the OM/EM builder with projection years switched on
naaproj_build_ext <- naaom_build_em
body(naaproj_build_ext) <- do.call(substitute, list(body(naaproj_build_ext), list(
  Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 5)))))

naaproj_fits <- local({
  cached <- list()
  function(form) {
    if(!is.null(cached[[form]])) return(cached[[form]])
    om <- naaom_make_om(NAA_re = form, sigmaNAA = 0.25, rho_age = 0,
                        rho_year = if(form == "iid") 0 else 0.6, seed = 808)
    dat <- naaom_om_data(om)
    il0 <- naaom_build_em(dat, NAA_re = form)
    il1 <- naaproj_build_ext(dat, NAA_re = form)
    f0 <- suppressWarnings(fit_model(il0$data, il0$par, il0$map, random = "ln_NAA", silent = TRUE))
    f1 <- suppressWarnings(fit_model(il1$data, il1$par, il1$map, random = "ln_NAA", silent = TRUE))
    cached[[form]] <<- list(il0 = il0, il1 = il1, f0 = f0, f1 = f1)
    cached[[form]]
  }
})

# fitted innovation in the terminal year, ages two and older
naaproj_eta_T <- function(fit, n_yrs) {
  ln_NAA <- fit$env$parList()$ln_NAA
  ln_NAA[1,1,n_yrs,1,2:8,1] - log(fit$rep$NAA_pred[1,1,n_yrs,1,2:8,1])
}

test_that("the state runs over the projection years only when it reaches the terminal year", {
  sweep_ext <- sweep_input
  body(sweep_ext) <- do.call(substitute, list(body(sweep_ext), list(
    Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 3)))))
  dims <- list(n_seas = 3, n_yrs = 11, n_ages = 6)
  build <- function(...) suppressMessages(sweep_ext(dims = dims, stop_after = "biol",
                                                    biol = list(M_spec = "fix", ...)))

  il <- build(NAA_re = "iid", NAA_re_seasons = "all")
  expect_equal(dim(il$par$ln_NAA), c(1, 3, 14, 3, 6, 2))
  expect_equal(il$data$naa_re_yrs, 2:14)
  expect_equal(dim(il$data$naa_sigma_blocks)[3], 14)
  expect_equal(dim(il$data$map_ln_NAA)[3], 14)

  # projected cells start at zero, are estimated where the fitted state is, and age one stays out
  map <- array(as.integer(as.character(il$map$ln_NAA)), dim = dim(il$par$ln_NAA))
  expect_true(all(il$par$ln_NAA[,,12:14,,,] == 0))
  expect_true(all(!is.na(map[,,12:14,,2:6,])))
  expect_true(all(is.na(map[,,12:14,,1,])))
  expect_equal(il$data$n_est_naa_re, 3 * 13 * 3 * 5 * 2)

  # the projection years take the terminal year's sigma blocks
  expect_equal(il$data$naa_sigma_blocks[,,12:14,,,], il$data$naa_sigma_blocks[,,rep(11, 3),,,])

  # a seed over the fitted years alone leaves the projected cells at zero
  il_seed <- build(NAA_re = "iid", ln_NAA = array(1, dim = c(1, 3, 11, 3, 6, 2)))
  expect_true(all(il_seed$par$ln_NAA[,,1:11,,,] == 1))
  expect_true(all(il_seed$par$ln_NAA[,,12:14,,,] == 0))

  # a state that stops early, and no state at all, keep the fitted years
  il_early <- build(NAA_re = "iid", NAA_re_years = 2:8)
  expect_equal(dim(il_early$par$ln_NAA)[3], 11)
  expect_equal(il_early$data$naa_re_yrs, 2:8)
  expect_equal(dim(build(NAA_re = "none")$par$ln_NAA)[3], 11)
})

test_that("the projected cells leave the dynamics alone under every form", {
  sweep_ext <- sweep_input
  body(sweep_ext) <- do.call(substitute, list(body(sweep_ext), list(
    Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 3)))))
  dims <- list(n_seas = 3, n_yrs = 11, n_ages = 6)

  for(form in c("iid", "1dar1_a", "1dar1_y", "2dar1", "3dcond", "3dmarg")) {
    biol <- list(NAA_re = form, M_spec = "fix", NAA_re_seasons = "all",
                 NAA_re_region = if(form == "iid") "iid" else "us")
    b0 <- suppressMessages(sweep_input(dims = dims, biol = biol))
    b1 <- suppressMessages(sweep_ext(dims = dims, biol = biol))
    o0 <- fit_model(b0$data, b0$par, b0$map, do_optim = FALSE, silent = TRUE)
    o1 <- fit_model(b1$data, b1$par, b1$map, do_optim = FALSE, silent = TRUE)
    r0 <- o0$report(o0$par)
    r1 <- o1$report(o1$par)

    expect_equal(r1$SSB, r0$SSB, info = form)
    expect_true(all(is.finite(o1$gr(o1$par))), info = form)

    # every correlation starts at zero, so each projected cell at zero adds the independent density there
    n_new <- b1$data$n_est_naa_re - b0$data$n_est_naa_re
    sig <- exp(b1$par$ln_sigmaNAA[1])
    expect_equal(r1$NAA_state_nLL - r0$NAA_state_nLL, -n_new * dnorm(0, 0, sig, log = TRUE),
                 tolerance = 1e-8, info = form)
  } # end form loop
})

test_that("integrating out the projected cells leaves the marginal likelihood unchanged", {
  fits <- naaproj_fits("iid")
  n_yrs <- length(fits$il1$data$years)

  # the fixed effects other than recruitment deviations, which gain projected cells of their own
  fixed_named <- function(fit) {
    p <- naaom_fixed(fit)
    p[names(p) != "ln_RecDevs"]
  }
  expect_equal(fixed_named(fits$f1), fixed_named(fits$f0), tolerance = 1e-6)
  expect_equal(fits$f1$rep$SSB, fits$f0$rep$SSB, tolerance = 1e-8)

  # iid, no bias correction: each projected cell is at zero with the process sd as its standard error
  ln_NAA <- fits$f1$env$parList()$ln_NAA
  expect_equal(max(abs(ln_NAA[1,1,n_yrs + seq_len(naaproj_n_proj),1,,1])), 0, tolerance = 1e-8)
  sr <- summary(RTMB::sdreport(fits$f1), "random")
  se <- array(sr[rownames(sr) == "ln_NAA", "Std. Error"], dim = c(n_yrs - 1 + naaproj_n_proj, 7))
  expect_equal(as.vector(se[n_yrs - 1 + seq_len(naaproj_n_proj), ]),
               rep(exp(naaom_fixed(fits$f1, "ln_sigmaNAA")), naaproj_n_proj * 7), tolerance = 1e-6)
})

test_that("a year autoregression forecasts the terminal innovation forward", {
  fits <- naaproj_fits("1dar1_y")
  n_yrs <- length(fits$il1$data$years)
  pars <- fits$f1$env$parList()

  rho <- rho_trans(pars$NAA_pe_pars[1,1,2,1])
  eta_T <- naaproj_eta_T(fits$f1, n_yrs)
  forecast <- outer(rho^seq_len(naaproj_n_proj), eta_T)
  expect_equal(unname(pars$ln_NAA[1,1,n_yrs + seq_len(naaproj_n_proj),1,2:8,1]), forecast, tolerance = 1e-8)
  expect_equal(fits$f1$rep$SSB, fits$f0$rep$SSB, tolerance = 1e-8)
})

test_that("naa_devs land where the fit applies its state", {
  fits <- naaproj_fits("1dar1_y")
  il <- fits$il1
  rep <- fits$f1$rep
  n_yrs <- length(il$data$years)

  # replay year T - 1 at year T's fishing with year T's innovation, which has to give the fit's year T
  rp_back <- rep
  rp_back$NAA[,,n_yrs,,,] <- rep$NAA[,,n_yrs - 1,,,]
  rp_back$NAA0[,,n_yrs,,,] <- rep$NAA0[,,n_yrs - 1,,,]
  rp_back$Fmort[,n_yrs,,] <- rep$Fmort[,n_yrs - 1,,]
  naa_devs <- array(1, dim = c(1, 1, 3, 1, 8, 1))
  naa_devs[1,1,2,1,2:8,1] <- exp(naaproj_eta_T(fits$f1, n_yrs))
  run <- function(nd) project_at_F(sum(rep$Fmort[1,n_yrs,,]), n_proj_yrs = 3, data = il$data, rep = rp_back,
                                   rec_yrs = n_yrs, naa_devs = nd)

  with_devs <- run(naa_devs)
  without <- run(NULL)
  expect_equal(with_devs$proj_NAA[1,1,2,1,2:8,1], rep$NAA[1,1,n_yrs,1,2:8,1], tolerance = 1e-12)
  expect_equal(with_devs$proj_NAA0[1,1,2,1,2:8,1], rep$NAA0[1,1,n_yrs,1,2:8,1], tolerance = 1e-12)
  expect_equal(without$proj_NAA[1,1,2,1,2:8,1], rep$NAA_pred[1,1,n_yrs,1,2:8,1], tolerance = 1e-12)

  # ones are no innovation at all, and a mis-sized array is refused
  expect_identical(run(array(1, dim = c(1, 1, 3, 1, 8, 1)))$proj_NAA, without$proj_NAA)
  expect_error(run(array(1, dim = c(1, 1, 2, 1, 8, 1))), "naa_devs should be dimensioned")
})

test_that("a retrospective peel carries the projected cells past the new terminal year", {
  fits <- naaproj_fits("1dar1_y")
  il <- fits$il1
  pars <- fits$f1$env$parList()
  n_yrs <- length(il$data$years)
  peel <- 3
  n_keep <- n_yrs - peel
  proj <- seq_len(naaproj_n_proj)

  cut <- SPoRC:::truncate_yr(peel, il$data, pars, il$map)
  expect_equal(dim(cut$retro_parameters$ln_NAA)[3], n_keep + naaproj_n_proj)
  expect_equal(cut$retro_data$naa_re_yrs, 2:(n_keep + naaproj_n_proj))
  expect_equal(dim(cut$retro_data$naa_sigma_blocks)[3], n_keep + naaproj_n_proj)
  expect_equal(length(cut$retro_mapping$ln_NAA), length(cut$retro_parameters$ln_NAA))

  # the projected cells keep the innovations as starts, rather than the peeled years' log numbers
  expect_equal(cut$retro_parameters$ln_NAA[,,n_keep + proj,,,], pars$ln_NAA[,,n_yrs + proj,,,])

  # and the peeled fit forecasts from its own terminal year
  peeled <- suppressWarnings(fit_model(cut$retro_data, cut$retro_parameters, cut$retro_mapping,
                                       random = "ln_NAA", silent = TRUE))
  pp <- peeled$env$parList()
  rho <- rho_trans(pp$NAA_pe_pars[1,1,2,1])
  expect_equal(unname(pp$ln_NAA[1,1,n_keep + proj,1,2:8,1]),
               outer(rho^proj, naaproj_eta_T(peeled, n_keep)), tolerance = 1e-8)
})

test_that("the closed loop conditions on the fitted years alone", {
  fits <- naaproj_fits("1dar1_y")
  il <- fits$il1
  sdr <- suppressWarnings(RTMB::sdreport(fits$f1))
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = 4, n_sims = 1, data = il$data,
                                                parameters = il$par, mapping = il$map, sd_rep = sdr,
                                                rep = fits$f1$rep, random = "ln_NAA")
  expect_equal(dim(sim_list$naa_eta_input)[3], length(il$data$years))
})
