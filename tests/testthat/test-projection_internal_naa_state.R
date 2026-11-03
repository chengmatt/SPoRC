# The projection inside the objective on a state-space numbers at age fit. With n_proj_yrs_devs the state
# runs over the projection years as innovations: their standard error is the process sd (1e-3), the
# projected numbers' variance grows by exactly that much (5%), year one starts from the terminal state, and
# the internal projection matches the post hoc one with naa_devs (1e-8), or without them under naa_opt terminal.

library(SPoRC)
library(testthat)

n_proj_yrs <- 6
naa_dims <- list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_seas = 1, n_yrs = 22, n_ages = 8)
n_yrs <- naa_dims$n_yrs
sigma_naa <- 0.1 # process sd kept fixed below. at the setup's 0.3 this toy has an unidentified direction

# the sweep model with five projected deviation years
build <- sweep_input
body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 5)))))

# a deterministic pass seeds the state over the fitted years
il0 <- build(dims = naa_dims)
obj0 <- fit_model(il0$data, il0$par, il0$map, do_optim = FALSE, silent = TRUE)
seed <- log(array(obj0$report(obj0$par)$NAA[,,1:n_yrs,1,,], dim = c(1, 1, n_yrs, 1, naa_dims$n_ages, 1)))
il <- build(dims = naa_dims, biol = list(NAA_re = "iid", NAA_sigma_spec = "fix", M_spec = "fix", ln_NAA = seed))
il$par$ln_sigmaNAA[] <- log(sigma_naa) # the setup does not read a starting value for it

post_hoc_inputs <- function(fit) {
  f_term <- sum(fit$rep$Fmort[1,n_yrs,,])
  build_projection_inputs(data = fit$data, rep = fit$rep, refpts = list(f_ref_pt = f_term, b_ref_pt = array(1, dim = c(1, 1))), n_proj_yrs = n_proj_yrs, t_spawn = 0)
}

test_that("the state is allocated and mapped over the projection years, seeded at zero there", {

  expect_equal(dim(il$par$ln_NAA)[3], n_yrs + 5)
  expect_equal(dim(il$data$naa_sigma_blocks)[3], n_yrs + 5)
  expect_equal(il$data$naa_re_yrs, 2:(n_yrs + 5))
  map_naa <- array(il$map$ln_NAA, dim = dim(il$par$ln_NAA))
  expect_true(all(!is.na(map_naa[1,1,n_yrs + 1:5,1,il$data$naa_re_ages,1])))
  expect_true(all(il$par$ln_NAA[1,1,n_yrs + 1:5,,,] == 0))
  expect_equal(il$par$ln_NAA[1,1,1:n_yrs,,,], seed[1,1,,,,])

})

test_that("the projected innovations have the process sd, and the projection takes them with their uncertainty", {

  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "terminal")
  fit <- fit_model(ip$data, ip$par, ip$map, random = "ln_NAA", newton_loops = 1, silent = TRUE)
  sdr <- RTMB::sdreport(fit)
  se_of <- function(name) as.numeric(sdr$sd[names(sdr$value) == name])
  pars <- fit$env$parList(par = fit$env$last.par.best)
  state_ages <- ip$data$naa_re_ages

  # with no data and no bias correction the projected innovations sit at zero with the process sd
  proj_cells <- pars$ln_NAA[1,1,n_yrs + 1:5,1,state_ages,1]
  expect_true(all(abs(proj_cells) < 1e-6))
  sd_random <- summary(sdr, "random")
  ln_naa_rows <- which(rownames(sd_random) == "ln_NAA")
  est_cells <- which(!is.na(array(ip$map$ln_NAA, dim = dim(ip$par$ln_NAA))))
  proj_rows <- ln_naa_rows[arrayInd(est_cells, dim(ip$par$ln_NAA))[,3] > n_yrs] # the estimated cells past the fitted years, in array order
  expect_equal(length(proj_rows), 5 * length(state_ages))
  expect_equal(unname(sd_random[proj_rows, "Std. Error"]), rep(sigma_naa, length(proj_rows)), tolerance = 1e-3)

  # year one replays the terminal year, so its standard error is the terminal one
  se_ssb <- se_of("log_proj_SSB")
  expect_equal(se_ssb[1], sdr$sd[names(sdr$value) == "log_SSB"][n_yrs], tolerance = 1e-8)
  expect_true(all(diff(se_ssb) > 0))

  # the numbers at age the projection starts from are the state
  state_terminal <- exp(pars$ln_NAA[1,1,n_yrs,1,,1])
  expect_equal(as.numeric(fit$rep$proj_NAA[1,1,1,1,state_ages,1]), as.numeric(state_terminal[state_ages]), tolerance = 1e-10)

  # the same fit with the projected state cells kept fixed: year two's numbers at a state age then differ in
  # variance by exactly the process variance, since the innovation is independent of everything else
  il_fixed <- il
  map_fixed <- array(il$map$ln_NAA, dim = dim(il$par$ln_NAA))
  map_fixed[,,n_yrs + 1:5,,,] <- NA
  il_fixed$map$ln_NAA <- factor(map_fixed)
  il_fixed$data$map_ln_NAA <- array(as.numeric(il_fixed$map$ln_NAA), dim = dim(map_fixed))
  il_fixed$data$naa_re_yrs <- il$data$naa_re_yrs[il$data$naa_re_yrs <= n_yrs]
  il_fixed$data$n_est_naa_re <- length(il_fixed$data$naa_re_yrs) * length(state_ages)
  ip_fixed <- Setup_Mod_Projection(il_fixed, n_proj_yrs = n_proj_yrs, f_opt = "terminal")
  fit_fixed <- fit_model(ip_fixed$data, ip_fixed$par, ip_fixed$map, random = "ln_NAA", newton_loops = 1, silent = TRUE)
  sdr_fixed <- RTMB::sdreport(fit_fixed)
  expect_equal(fit$optim$objective, fit_fixed$optim$objective, tolerance = 1e-6)

  naa_se <- array(se_of("log_proj_NAA"), dim = c(1, 1, n_proj_yrs, 1, naa_dims$n_ages, 1))
  naa_se_fixed <- array(as.numeric(sdr_fixed$sd[names(sdr_fixed$value) == "log_proj_NAA"]), dim = c(1, 1, n_proj_yrs, 1, naa_dims$n_ages, 1))
  a <- state_ages[2]
  expect_equal(naa_se[1,1,2,1,a,1]^2 - naa_se_fixed[1,1,2,1,a,1]^2, sigma_naa^2, tolerance = 0.05)
  expect_equal(as.numeric(fit$rep$proj_SSB[1,1,]), as.numeric(fit_fixed$rep$proj_SSB[1,1,]), tolerance = 1e-4) # two optimizations, so not to the last digit

  # the post hoc projection with the fit's innovations passed as naa_devs
  rec_devs <- array(1, dim = c(1, 1, n_proj_yrs))
  rec_devs[1,1,2:n_proj_yrs] <- exp(pars$ln_RecDevs[1,1,n_yrs + 1:(n_proj_yrs - 1)])
  naa_devs <- array(1, dim = c(1, 1, n_proj_yrs, 1, naa_dims$n_ages, 1))
  for(k in 1:(n_proj_yrs - 1)) naa_devs[1,1,k + 1,1,state_ages,1] <- exp(pars$ln_NAA[1,1,n_yrs + k,1,state_ages,1])
  ext <- do.call(Do_Population_Projection, c(post_hoc_inputs(fit), list(fmort_opt = "Input", rec_devs = rec_devs, naa_devs = naa_devs)))
  expect_equal(as.numeric(fit$rep$proj_SSB[1,1,]), as.numeric(ext$proj_SSB[1,1,]), tolerance = 1e-8)
  expect_equal(as.numeric(fit$rep$proj_NAA), as.numeric(ext$proj_NAA), tolerance = 1e-8)

  # the terminal state advanced with no innovations is the post hoc projection without naa_devs
  ip_term <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "terminal", naa_opt = "terminal")
  fit_term <- fit_model(ip_term$data, fit$env$parList(par = fit$env$last.par.best), ip_term$map, random = "ln_NAA", do_optim = FALSE, silent = TRUE)
  ext_term <- do.call(Do_Population_Projection, c(post_hoc_inputs(fit_term), list(fmort_opt = "Input", rec_devs = rec_devs)))
  expect_equal(as.numeric(fit_term$rep$proj_NAA), as.numeric(ext_term$proj_NAA), tolerance = 1e-8) # equal to the propagated run here, since iid innovations with no data sit at zero

})

test_that("naa_devs move the post hoc projection where the state applies and nowhere else", {

  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "terminal")
  fit <- fit_model(ip$data, ip$par, ip$map, random = "ln_NAA", newton_loops = 1, silent = TRUE)
  inputs <- post_hoc_inputs(fit)
  plain <- do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "Input")))

  # a doubling of age three at the start of year two
  naa_devs <- array(1, dim = c(1, 1, n_proj_yrs, 1, naa_dims$n_ages, 1))
  naa_devs[1,1,2,1,3,1] <- 2
  bumped <- do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "Input", naa_devs = naa_devs)))
  expect_equal(bumped$proj_NAA[1,1,2,1,3,1], 2 * plain$proj_NAA[1,1,2,1,3,1], tolerance = 1e-10)
  expect_equal(as.numeric(bumped$proj_NAA[1,1,1,1,,1]), as.numeric(plain$proj_NAA[1,1,1,1,,1]), tolerance = 1e-10)
  expect_equal(as.numeric(bumped$proj_NAA[1,1,2,1,-3,1]), as.numeric(plain$proj_NAA[1,1,2,1,-3,1]), tolerance = 1e-10)
  expect_gt(bumped$proj_SSB[1,1,2], plain$proj_SSB[1,1,2])

  expect_error(do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "Input", naa_devs = naa_devs[,,1:2,,,,drop = FALSE]))), "naa_devs should be dimensioned")

})
