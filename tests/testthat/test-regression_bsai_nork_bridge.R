# Self-validating bridge test: the expectations are the 2023 BSAI northern rockfish assessment's own
# reported quantities in sgl_rg_bsai_nork_data$admb, evaluated at its estimate without optimizing.

library(SPoRC)
library(testthat)
data("sgl_rg_bsai_nork_data")

test_that("BSAI northern rockfish reproduces the 2023 ADMB assessment at its own MLE", {

  dat <- sgl_rg_bsai_nork_data
  n_yrs <- length(dat$years)
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

  # the survey curve is supplied as a fixed input held flat past age 30, which logist1 cannot
  # express, so this compares the likelihoods rather than the selectivity form
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

  # the assessment's last reported age is a plus group while SPoRC has ages past it, so those
  # columns are summed before comparing
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

  # SPoRC writes each component as a proper density while the assessment drops normalizing
  # constants, so each comparison subtracts exactly the constants it omits
  c2pi <- 0.5 * log(2 * pi)

  expect_equal(sum(r$FishAgeComps_nLL), dat$admb$datalikecomp[["fish.unbiased.ac"]],
               tolerance = 1e-5)
  expect_equal(sum(r$FishLenComps_nLL), dat$admb$datalikecomp[["fish.lc"]],
               tolerance = 1e-5)
  expect_equal(sum(r$SrvAgeComps_nLL), dat$admb$datalikecomp[["aisrv.ac"]],
               tolerance = 1e-5)

  # the assessment's survey term keeps the log sigma but drops the sqrt(2 pi), and its standard
  # errors vary by year, so the constant is summed over the observations
  srv_like <- sum(r$SrvIdx_nLL) -
    sum(c2pi + log(dat$ObsSrvIdx_SE[dat$UseSrvIdx == 1]))
  expect_equal(srv_like, dat$admb$datalikecomp[["aisurvlike"]], tolerance = 1e-3)

  # The F penalty is a weighted sum of squares on the F deviations, with the weight
  # kept in sigmaF rather than applied outside the sum.
  n_catch_obs <- sum(dat$UseCatch)
  f_pen <- sum(r$Fmort_nLL) -
    n_catch_obs * (c2pi + as.vector(input_list$par$ln_sigmaF))
  expect_equal(as.vector(f_pen), dat$admb$pen_likecomp[["Fmortpen"]], tolerance = 1e-5)

  # the assessment's recruitment penalty keeps its log sigmaR terms and drops the sqrt(2 pi),
  # over the recruitment and initial age deviations together. negative, sigmaR being 0.75
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

  # the catch term is the one component that does not land on the assessment's value, both
  # being numerically zero against an objective of 555: 2e-05 there against 4e-03 here.
  #
  # that is a difference in how near-exact the catch is driven, so it is asserted in
  # absolute terms, a relative tolerance on two numbers this small meaning little
  catch_ssq <- sum(as.vector(r$Catch_nLL) -
                     (c2pi + as.vector(input_list$par$ln_sigmaC)))
  expect_lt(abs(catch_ssq - dat$admb$datalikecomp[["catch.like"]]), 1e-2)

  # the whole objective, like for like, with the assessment's maturity likelihood removed since
  # SPoRC fixes maturity. its objective is reported to six figures, the tolerance floor
  like_for_like <- sum(r$FishAgeComps_nLL) + sum(r$FishLenComps_nLL) +
    sum(r$SrvAgeComps_nLL) + srv_like + catch_ssq + as.vector(f_pen) + rec_like +
    sel_pri + (r$M_nLL - log(sqrt(2 * pi) * dat$cv_M)) +
    (r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q))
  admb_total <- dat$admb$datalikecomp[["obj_fun"]] - dat$admb$datalikecomp[["mat_like"]]
  expect_equal(like_for_like, admb_total, tolerance = 1e-4)

  expect_jnLL_decomposes(obj)
})
