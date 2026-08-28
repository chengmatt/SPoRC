# Self-validating bridge test. The expectations are not stored SPoRC output: they are
# the 2023 BSAI northern rockfish assessment's own reported quantities, shipped in
# sgl_rg_bsai_nork_data$admb. Every parameter is set to the assessment's maximum
# likelihood estimate and the model is evaluated there without optimizing, so a failure
# means the population dynamics, a likelihood, or a selectivity form no longer
# reproduces the assessment at a point where it is known to. Do not loosen the
# tolerances to make this pass. See tests/README.md and the BSAI northern rockfish case
# study vignette.

library(SPoRC)
library(testthat)
data("sgl_rg_bsai_nork_data")

test_that("BSAI northern rockfish reproduces the 2023 ADMB assessment at its own MLE", {

  dat <- sgl_rg_bsai_nork_data
  n_yrs <- length(dat$years)
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

  # The survey curve is supplied as a fixed input with the assessment's age 30 edge
  # hold applied, which SPoRC's logist1 form cannot express. That isolates the
  # likelihoods from the selectivity form; the companion pinned test refits the
  # uncapped logistic instead.
  input_list <- seed_bsai_nork_mle(build_bsai_nork_input(dat), dat)
  input_list$data <- cap_bsai_nork_srv_sel(input_list$data, dat)

  obj <- fit_model(input_list$data, input_list$par, input_list$map,
                   do_optim = FALSE, silent = TRUE)
  r <- obj$rep

  # Selectivity is built by SPoRC's own logistic form over the assessment's observed
  # age range, so matching the curves is a check on the form rather than on the data.
  expect_equal(as.vector(r$fish_sel[1, 1, 1, 1, seq_len(n_obs_ages), 1, 1]),
               dat$admb$sel_fsh, tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$srv_sel[1, 1, 1, 1, seq_len(n_obs_ages), 1, 1]),
               dat$admb$sel_srv, tolerance = 1e-5, ignore_attr = TRUE)

  # The population at the assessment's MLE. The assessment's last reported age is a
  # plus group while SPoRC carries ages past it, so those columns are summed before
  # comparing.
  naa <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
  naa_obs <- cbind(naa[, 1:(n_obs_ages - 1)], rowSums(naa[, n_obs_ages:n_ages]))
  expect_equal(naa_obs, dat$admb$NAA, tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$SSB)[1:n_yrs], dat$admb$SSB,
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$Rec)[1:n_yrs], dat$admb$Rec,
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$Fmort), dat$admb$Fmort,
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$PredSrvIdx), dat$admb$pred_srv,
               tolerance = 1e-3, ignore_attr = TRUE)

  # Likelihood components. SPoRC writes each component as a proper density while the
  # assessment drops normalizing constants, so each comparison subtracts exactly the
  # constants the assessment omits. What is left is a like for like comparison.
  c2pi <- 0.5 * log(2 * pi)

  expect_equal(sum(r$FishAgeComps_nLL), dat$admb$datalikecomp[["fish.unbiased.ac"]],
               tolerance = 1e-5)
  expect_equal(sum(r$FishLenComps_nLL), dat$admb$datalikecomp[["fish.lc"]],
               tolerance = 1e-5)
  expect_equal(sum(r$SrvAgeComps_nLL), dat$admb$datalikecomp[["aisrv.ac"]],
               tolerance = 1e-5)

  # The assessment's survey statement keeps the log sigma term but drops the
  # sqrt(2 pi) constant. Its standard errors vary by year, so the constant is summed
  # over the observations rather than counted.
  srv_like <- sum(r$SrvIdx_nLL) -
    sum(c2pi + log(dat$ObsSrvIdx_SE[dat$UseSrvIdx == 1]))
  expect_equal(srv_like, dat$admb$datalikecomp[["aisurvlike"]], tolerance = 1e-3)

  # The F penalty is a weighted sum of squares on the F deviations, with the weight
  # carried in sigmaF rather than applied outside the sum.
  n_catch_obs <- sum(dat$UseCatch)
  f_pen <- sum(r$Fmort_nLL) -
    n_catch_obs * (c2pi + as.vector(input_list$par$ln_sigmaF))
  expect_equal(as.vector(f_pen), dat$admb$pen_likecomp[["Fmortpen"]], tolerance = 1e-5)

  # The assessment's recruitment penalty keeps its log sigmaR terms and drops only the
  # sqrt(2 pi) constants, over the recruitment deviations and the initial age
  # deviations together. It is negative because sigmaR is 0.75.
  n_recdev <- length(dat$mle$rec_dev)
  n_fydev <- length(dat$mle$fydev)
  rec_like <- sum(r$Rec_nLL) + sum(r$Init_Rec_nLL) - (n_recdev + n_fydev) * c2pi
  expect_equal(rec_like, dat$admb$pen_likecomp[["reclike"]], tolerance = 1e-5)

  # The survey selectivity constraint is a normal prior on the realized selectivity
  # value at age 30, so its constant is the one that statement omits.
  sel_pri <- sum(r$sel_nLL) - log(sqrt(2 * pi) * 0.003)
  expect_equal(sel_pri, dat$admb$pen_likecomp[["prior_sel"]], tolerance = 1e-5)

  expect_equal(r$M_nLL - log(sqrt(2 * pi) * dat$cv_M),
               dat$admb$pen_likecomp[["prior_m"]], tolerance = 1e-5)
  expect_equal(r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q),
               dat$admb$pen_likecomp[["prior_q"]], tolerance = 1e-5)

  # The catch statement is the one component that does not land on the assessment's
  # value. Both are numerically zero against an objective of 555; the assessment fits
  # catch to 2e-05 and SPoRC to 4e-03, which is a difference in how near-exact catch
  # is driven rather than in the statement. Asserted in absolute terms because a
  # relative tolerance on two numbers this close to zero means nothing.
  catch_ssq <- sum(as.vector(r$Catch_nLL) -
                     (c2pi + as.vector(input_list$par$ln_sigmaC)))
  expect_lt(abs(catch_ssq - dat$admb$datalikecomp[["catch.like"]]), 1e-2)

  # The whole objective, like for like. The assessment estimates maturity internally
  # and SPoRC fixes it, so mat_like is removed from the assessment's side. Its
  # obj_fun is reported to six significant figures, which is the tolerance floor here.
  like_for_like <- sum(r$FishAgeComps_nLL) + sum(r$FishLenComps_nLL) +
    sum(r$SrvAgeComps_nLL) + srv_like + catch_ssq + as.vector(f_pen) + rec_like +
    sel_pri + (r$M_nLL - log(sqrt(2 * pi) * dat$cv_M)) +
    (r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q))
  admb_total <- dat$admb$datalikecomp[["obj_fun"]] - dat$admb$datalikecomp[["mat_like"]]
  expect_equal(like_for_like, admb_total, tolerance = 1e-4)

  expect_jnLL_decomposes(obj)
})
