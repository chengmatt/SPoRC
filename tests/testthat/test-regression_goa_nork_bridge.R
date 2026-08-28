# Self-validating bridge test. The expectations are not stored SPoRC output: they are
# the 2024 GOA northern rockfish assessment's own reported quantities, shipped in
# sgl_rg_goa_nork_data$admb. Every parameter is set to the assessment's maximum
# likelihood estimate and the model is evaluated there without optimizing, so a failure
# means the population dynamics, a likelihood, or a selectivity form no longer
# reproduces the assessment at a point where it is known to. Do not loosen the
# tolerances to make this pass. See tests/README.md and the GOA northern rockfish case
# study vignette.

library(SPoRC)
library(testthat)
data("sgl_rg_goa_nork_data")

test_that("GOA northern rockfish reproduces the 2024 ADMB assessment at its own MLE", {

  dat <- sgl_rg_goa_nork_data
  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)

  input_list <- seed_goa_nork_mle(build_goa_nork_input(dat), dat)

  obj <- fit_model(input_list$data, input_list$par, input_list$map,
                   do_optim = FALSE, silent = TRUE)
  r <- obj$rep

  # Selectivity is built by SPoRC's own logistic form, so matching the assessment's
  # curves is a check on the form rather than on the data.
  expect_equal(as.vector(r$fish_sel[1, 1, 1, 1, , 1, 1]), dat$admb$sel[, 1],
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(as.vector(r$srv_sel[1, 1, 1, 1, , 1, 1]), dat$admb$sel[, 2],
               tolerance = 1e-6, ignore_attr = TRUE)

  # The population at the assessment's MLE, over every year and model age.
  naa <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
  expect_equal(naa, t(dat$admb$NAA), tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(as.vector(r$SSB)[1:n_yrs], dat$admb$SSB,
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(as.vector(r$Rec)[1:n_yrs], dat$admb$Rec,
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(as.vector(r$Fmort), dat$admb$Fmort,
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(as.vector(r$PredCatch), dat$admb$pred_catch,
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(as.vector(r$PredSrvIdx)[dat$UseSrvIdx[1, , 1, 1] == 1],
               dat$admb$pred_srv, tolerance = 1e-6, ignore_attr = TRUE)

  # Likelihood components. SPoRC writes each component as a proper density while the
  # assessment drops normalizing constants, so each comparison subtracts exactly the
  # constants the assessment omits. What is left is a like for like comparison. The
  # composition tolerances absorb the two templates' different robustifying
  # constants.
  expect_equal(sum(r$FishAgeComps_nLL), dat$admb$like_fish_age, tolerance = 1e-3)
  expect_equal(sum(r$SrvAgeComps_nLL), dat$admb$like_srv_age, tolerance = 1e-3)
  expect_equal(sum(r$FishLenComps_nLL), dat$admb$like_fish_size, tolerance = 1e-3)

  # The assessment's catch statement is a weighted sum of squares; subtracting
  # SPoRC's per observation constants leaves exactly that.
  ln_sigmaC <- input_list$par$ln_sigmaC
  catch_ssq <- sum(as.vector(r$Catch_nLL) - (0.5 * log(2 * pi) + as.vector(ln_sigmaC)))
  expect_equal(catch_ssq, dat$admb$ssqcatch, tolerance = 1e-6)

  # The assessment's survey statement keeps the log sigma term but drops the
  # sqrt(2 pi) constant, and carries a weight of 0.25.
  n_srv_obs <- sum(dat$UseSrvIdx)
  srv_like <- dat$srv_wt * (sum(r$SrvIdx_nLL) - n_srv_obs * 0.5 * log(2 * pi))
  expect_equal(srv_like, dat$admb$like_srv, tolerance = 1e-6)

  # The F penalty is a weighted sum of squares on the F deviations.
  n_catch_obs <- sum(dat$UseCatch)
  f_reg <- dat$fmort_wt *
    (sum(r$Fmort_nLL) - n_catch_obs * (0.5 * log(2 * pi) + as.vector(input_list$par$ln_sigmaF)))
  expect_equal(as.vector(f_reg), dat$admb$f_regularity, tolerance = 1e-6)

  expect_equal(r$M_nLL - log(sqrt(2 * pi) * dat$cv_M), dat$admb$nll_M, tolerance = 1e-6)
  expect_equal(r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q), dat$admb$nll_q, tolerance = 1e-6)

  # The assessment's recruitment penalty is the lognormal bias corrected sum of
  # squares over all recruitment deviations plus the initial age deviations
  # excluding the plus group, with no normalizing constants. SPoRC's Rec_nLL and
  # Init_Rec_nLL carry the same penalties plus the constants subtracted here; the
  # plus group deviation sits at zero in this parameterization so it contributes
  # nothing on either side.
  s <- dat$sigmaR
  n_recdev <- length(dat$mle$log_Rt)
  n_initdev_pen <- n_ages - 2
  rec_like <- sum(r$Rec_nLL) - n_recdev * (0.5 * log(2 * pi) + log(s)) +
    sum(r$Init_Rec_nLL) - n_initdev_pen * (0.5 * log(2 * pi) + log(s))
  expect_equal(rec_like, dat$admb$like_rec, tolerance = 1e-6)

  expect_jnLL_decomposes(obj)
})
