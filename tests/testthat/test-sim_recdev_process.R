# Checks the recruitment deviations the operating model draws under each process error, against
# the density get_recdev_pe_nLL penalizes them at. Year one of an ar1 is the case that differed.

recdev_sim_list <- function(sigmaR, RecDevs_model, rho = 0) {
  n_yrs <- naaom_cfg$n_yrs
  n_ages <- naaom_cfg$n_ages
  curve <- function(slope, infl, scale = 1)
    array(rep(scale / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
          dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))

  sl <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = NULL,
                      n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sl <- Setup_Sim_Containers(sl)
  sl <- Setup_Sim_Fishing(sl,
                          fish_sel_input = replicate(1, curve(3, 2)),
                          ret_sel_input = replicate(1, curve(3, 2)),
                          dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
                          Fmort_input = array(naaom_cfg$f_ramp, dim = c(1, n_yrs, 1, 1, 1)),
                          ISS_FishAgeComps = array(naaom_cfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)))
  sl <- Setup_Sim_Survey(sl,
                         srv_sel_input = replicate(1, curve(1, 3)),
                         ObsSrvIdx_SE = array(naaom_cfg$idx_se, dim = c(1, n_yrs, 1, 1)),
                         ISS_SrvAgeComps = array(naaom_cfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)))
  biol <- function(v) array(rep(v, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sl <- Setup_Sim_Biologicals(sl,
    natmort_input = replicate(1, array(naaom_cfg$M, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, biol(naaom_cfg$waa)),
    WAA_fish_input = replicate(1, array(rep(naaom_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(rep(naaom_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol(naaom_cfg$mat))))
  sl <- Setup_Sim_Tagging(sl, use_conv_fish_tagging = 0)
  sl$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))

  Setup_Sim_Rec(sl,
                R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
                ln_sigmaR = array(log(sigmaR), dim = c(2, 1, 1)),
                recruitment_opt = "mean_rec",
                init_age_strc = 1,
                RecDevs_model = RecDevs_model,
                RecDevs_rho = array(rho, dim = c(1, 1)))
}

recdev_draws <- function(sl, nrep) {
  sapply(seq_len(nrep), function(i) Simulate_Pop_Static(sim_list = sl, output_path = NULL)$ln_RecDevs[1, 1, , 1])
}

test_that("year one of an ar1 is drawn from its stationary marginal", {
  set.seed(404)
  sigmaR <- 0.5
  rho <- 0.7
  d <- recdev_draws(recdev_sim_list(sigmaR, "ar1", rho), 150)

  stationary <- sigmaR / sqrt(1 - rho^2)

  # drawing year one at sigmaR instead left the series short of its own variance for a few years
  expect_equal(sd(d[1, ]), stationary, tolerance = 0.2)
  expect_gt(sd(d[1, ]), sigmaR * 1.1)

  # and the years after it hold that same spread
  expect_equal(mean(apply(d[5:naaom_cfg$n_yrs, ], 1, sd)), stationary, tolerance = 0.15)

  # the reversion itself
  lag1 <- mean(sapply(21:naaom_cfg$n_yrs, function(y) cor(d[y, ], d[y - 1, ])))
  expect_equal(lag1, rho, tolerance = 0.12)
})

test_that("independent deviations are drawn at their own sigma in every year", {
  set.seed(405)
  sigmaR <- 0.4
  d <- recdev_draws(recdev_sim_list(sigmaR, "iid"), 150)

  expect_equal(sd(d[1, ]), sigmaR, tolerance = 0.2)
  expect_equal(mean(apply(d, 1, sd)), sigmaR, tolerance = 0.15)

  # nothing links one year to the next
  lag1 <- mean(sapply(21:naaom_cfg$n_yrs, function(y) cor(d[y, ], d[y - 1, ])))
  expect_lt(abs(lag1), 0.15)
})

test_that("an ar1 correlation on the unit circle is refused", {
  expect_error(recdev_sim_list(0.5, "ar1", 1), "inside \\(-1, 1\\)")
})
