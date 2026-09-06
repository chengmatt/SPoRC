# Self-validating bridge test: the expectations are the 2024 BSAI blackspotted and rougheye rockfish
# assessment's own quantities in sgl_rg_rebs_data$admb, evaluated at its estimate without optimizing.

library(SPoRC)
library(testthat)
data("sgl_rg_rebs_data")

test_that("BSAI rougheye reproduces the 2024 ADMB assessment at its own MLE", {

  dat <- sgl_rg_rebs_data
  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

  input_list <- seed_rebs_mle(build_rebs_input(dat), dat)

  obj <- fit_model(input_list$data, input_list$par, input_list$map,
                   do_optim = FALSE, silent = TRUE)
  r <- obj$rep

  # Selectivity is built by SPoRC's own logistic form, so matching the assessment's
  # curves is a check on the form rather than on the data.
  expect_equal(as.vector(r$fish_sel[1, 1, 1, 1, 1:n_obs_ages, 1, 1]),
               as.vector(dat$admb$sel_fsh[1, ]), tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$srv_sel[1, 1, 1, 1, 1:n_obs_ages, 1, 1]),
               as.vector(dat$admb$sel_srv), tolerance = 1e-5, ignore_attr = TRUE)

  # The population at the assessment's MLE. The assessment reports numbers at age
  # over the 43 observed ages with the model's ages 45 to 54 pooled into the last
  # column, so SPoRC's 52 ages are pooled the same way before comparing. The
  # tolerance is set by the assessment's own report, which has six significant
  # figures.
  #
  # Age 1 is recruitment and is checked separately below. The three terminal
  # recruits differ by a documented convention, and they age into the age 2 and 3
  # cells of the last two years, so the comparison covers every year at ages 4 and
  # older and every age over the years the assessment estimates a deviation for.
  # Those are the only cells the convention can reach.
  n_est <- length(dat$mle$rec_dev)
  naa <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
  naa_pooled <- cbind(naa[, 1:(n_obs_ages - 1)], rowSums(naa[, n_obs_ages:n_ages]))
  expect_equal(naa_pooled[, -(1:3)], dat$admb$NAA[, -(1:3)],
               tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(naa_pooled[1:n_est, -1], dat$admb$NAA[1:n_est, -1],
               tolerance = 1e-4, ignore_attr = TRUE)

  expect_equal(as.vector(r$SSB)[1:n_yrs], as.vector(dat$admb$SSB),
               tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(as.vector(r$Fmort)[1:n_yrs], as.vector(dat$admb$Fmort),
               tolerance = 1e-4, ignore_attr = TRUE)

  i_srv <- which(yrs %in% dat$yrs_srv)
  expect_equal(as.vector(r$PredSrvIdx)[i_srv], dat$admb$pred_srv[i_srv],
               tolerance = 1e-3, ignore_attr = TRUE)

  # Recruitment splits into the years the assessment estimates a deviation for and
  # the three terminal years it does not. Over the estimated years the two agree
  # outright. Over the terminal three the assessment multiplies mean recruitment by
  # exp(sigmaR^2 / 2) while leaving the estimated recruitments uncorrected, a legacy
  # ADMB convention SPoRC does not reproduce, so SPoRC sits exactly that factor low.
  # This is a documented convention difference, not an error: do not "fix" it by
  # adding a bias correction switch. It moves terminal spawning biomass by 8e-5,
  # because maturity at age 3 is 0.003.
  expect_equal(as.vector(r$Rec)[1:n_est], as.vector(dat$admb$Rec)[1:n_est],
               tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(as.vector(r$Rec)[(n_est + 1):n_yrs],
               as.vector(dat$admb$Rec)[(n_est + 1):n_yrs] * exp(-dat$sigmaR^2 / 2),
               tolerance = 1e-4, ignore_attr = TRUE)

  # Likelihood components. SPoRC writes each component as a proper density while the
  # assessment drops normalizing constants, so each comparison subtracts exactly the
  # constants the assessment omits. What is left is a like for like comparison.
  d <- input_list$data
  n_recdev <- length(dat$mle$rec_dev)
  n_fydev <- length(dat$mle$fydev)

  expect_equal(sum(r$FishAgeComps_nLL), dat$admb$datalikecomp[["fish.ac"]], tolerance = 1e-3)
  expect_equal(sum(r$FishLenComps_nLL), dat$admb$datalikecomp[["fish.lc"]], tolerance = 1e-3)
  expect_equal(sum(r$SrvAgeComps_nLL), dat$admb$datalikecomp[["AI_survey.sac"]], tolerance = 1e-3)
  expect_equal(sum(r$SrvLenComps_nLL), dat$admb$datalikecomp[["AI_survey.slc"]], tolerance = 1e-2)

  srv_idx <- sum(r$SrvIdx_nLL[r$SrvIdx_nLL != 0]) -
    sum(0.5 * log(2 * pi) + log(d$ObsSrvIdx_SE[r$SrvIdx_nLL != 0]))
  expect_equal(srv_idx, dat$admb$datalikecomp[["AI_survey.biom"]], tolerance = 1e-3)

  catch_ssq <- sum(r$Catch_nLL) - n_yrs * (0.5 * log(2 * pi) + log(0.1))
  expect_equal(catch_ssq,
               50 * sum((log(as.vector(dat$ObsCatch) + 1e-4) -
                           log(as.vector(r$PredCatch) + 1e-5))^2),
               tolerance = 1e-3)

  f_pen <- sum(r$Fmort_nLL) - n_yrs * (0.5 * log(2 * pi) + as.vector(input_list$par$ln_sigmaF))
  expect_equal(as.vector(f_pen), dat$admb$pen_likecomp[["Fmortpen"]], tolerance = 1e-5)

  expect_equal(r$M_nLL - log(sqrt(2 * pi) * d$M_prior$sd),
               dat$admb$pen_likecomp[["prior_m"]], tolerance = 1e-5)
  expect_equal(r$srv_q_nLL - log(sqrt(2 * pi) * d$srv_q_prior$sd),
               dat$admb$pen_likecomp[["AI_survey_prior_q"]], tolerance = 1e-5)

  rec_like <- sum(r$Rec_nLL) - n_recdev * 0.5 * log(2 * pi) +
    sum(r$Init_Rec_nLL) - n_fydev * 0.5 * log(2 * pi)
  expect_equal(rec_like, dat$admb$pen_likecomp[["reclike"]], tolerance = 1e-3)

  expect_jnLL_decomposes(obj)
})
