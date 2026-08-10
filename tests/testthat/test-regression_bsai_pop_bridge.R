# Self-validating bridge test. The expectations are not stored SPoRC output: they are
# the 2024 BSAI Pacific ocean perch assessment's own reported quantities, shipped in
# sgl_rg_bsai_pop_data$admb. Every parameter is set to the assessment's maximum
# likelihood estimate and the model is evaluated there without optimizing, so a failure
# means the population dynamics, a likelihood, or the bicubic selectivity form no longer
# reproduces the assessment at a point where it is known to. Do not loosen the
# tolerances to make this pass. See tests/README.md and the BSAI Pacific ocean perch
# case study vignette.

library(SPoRC)
library(testthat)
data("sgl_rg_bsai_pop_data")

test_that("BSAI Pacific ocean perch reproduces the 2024 ADMB assessment at its own MLE", {

  dat <- sgl_rg_bsai_pop_data
  n_yrs <- length(dat$years)
  n_obs_ages <- length(dat$obs_ages)
  nsel <- dat$nselages
  yr_ind <- match(dat$years, dat$admb$yrs)

  input_list <- seed_bsai_pop_mle(build_bsai_pop_input(dat), dat)
  obj <- fit_model(input_list$data, input_list$par, input_list$map,
                   do_optim = FALSE, silent = TRUE)
  r <- obj$rep

  # The assessment normalises selectivity to a maximum of one within each year for
  # reporting and multiplies F by the same factor, leaving their product invariant.
  # Its internal selectivity is the raw exponentiated bicubic surface, so both sides
  # are put on the reported convention before comparing. This is the sharpest single
  # check in the bridge: it covers the node orientation, the edge hold before 1964
  # and the edge hold beyond age 40 across all 65 by 38 cells at once.
  sel <- r$fish_sel[1, 1, 1:n_yrs, 1, 1:nsel, 1, 1]
  expect_equal(norm_bsai_pop_sel_by_year(sel), dat$admb$sel_fsh,
               tolerance = 1e-5, ignore_attr = TRUE)

  sel_max_yr <- apply(sel, 1, max)
  expect_equal(as.vector(r$Fmort[1, , 1, 1]), dat$admb$Fmort[yr_ind] / sel_max_yr,
               tolerance = 1e-5, ignore_attr = TRUE)

  # The population at the assessment's MLE. Its numbers at age are reported over the
  # 38 observed ages, so the final column is an age 40 plus aggregation rather than
  # the model's age 46 plus group; only the ages below it are compared.
  cmp <- 1:(n_obs_ages - 1)
  expect_equal(as.vector(r$NAA[1, 1, 1, 1, cmp, 1]), dat$admb$NAA[1, cmp],
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$SSB)[1:n_yrs], dat$admb$SSB[yr_ind],
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(as.vector(r$Rec)[1:n_yrs], dat$admb$Rec[yr_ind],
               tolerance = 1e-5, ignore_attr = TRUE)

  # SPoRC's Total_Biom is spawning time biomass and carries the exp(-Z t_spawn)
  # discount, while the assessment reports January 1 biomass. The like for like
  # quantity is built from numbers at age, which is a 6 percent difference if the
  # two are compared directly instead.
  naa_all <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
  expect_equal(as.vector(naa_all %*% dat$pop_waa), dat$admb$TotBiom[yr_ind],
               tolerance = 1e-5, ignore_attr = TRUE)

  # Likelihood components. SPoRC writes each component as a proper density while the
  # assessment drops normalising constants, so each comparison subtracts exactly the
  # constants the assessment omits. What is left is a like for like comparison.
  c2pi <- 0.5 * log(2 * pi)
  dlc <- unlist(dat$admb$datalikecomp)
  plc <- unlist(dat$admb$pen_likecomp)

  expect_equal(sum(r$FishAgeComps_nLL), dlc[["fish.ac"]], tolerance = 1e-5)
  expect_equal(sum(r$FishLenComps_nLL), dlc[["fish.lc"]], tolerance = 1e-5)
  expect_equal(sum(r$SrvAgeComps_nLL[, , , , 1]), dlc[["AI_survey.sac"]], tolerance = 1e-4)
  expect_equal(sum(r$SrvAgeComps_nLL[, , , , 2]), dlc[["EBS_survey.sac"]], tolerance = 1e-4)

  # The Aleutian Islands survey contributes a single year of length compositions, and
  # it is the one data component that does not land on the assessment's value. The gap
  # is 0.15 on a term of 7.5, which is a 0.07 percent difference in the predicted
  # proportions carried through an applied sample size of 224. It is consistent with
  # the multinomial offset convention rather than with the size at age transition, and
  # it is 0.017 percent of the objective. Held at its measured size so a real
  # regression in the length composition machinery still trips this.
  expect_equal(sum(r$SrvLenComps_nLL[, , , , 1]), dlc[["AI_survey.slc"]],
               tolerance = 2.5e-2)

  # The assessment's survey statements keep the log sigma term but drop the sqrt(2 pi)
  # constant. Their standard errors vary by year, so the constant is summed over the
  # observations rather than counted.
  srv_like <- function(fl) {
    keep <- dat$UseSrvIdx[, , , fl] == 1
    sum(r$SrvIdx_nLL[, , , fl][keep]) - sum(c2pi + log(dat$ObsSrvIdx_SE[, , , fl][keep]))
  }
  expect_equal(srv_like(1), dlc[["AI_survey.biom"]], tolerance = 2e-3)
  expect_equal(srv_like(2), dlc[["EBS_survey.biom"]], tolerance = 1e-3)

  # The F penalty is a weighted sum of squares on the F deviations, with the weight
  # carried in sigmaF rather than applied outside the sum.
  f_pen <- sum(r$Fmort_nLL) - n_yrs * (c2pi + as.vector(input_list$par$ln_sigmaF))
  expect_equal(as.vector(f_pen), plc[["Fmortpen"]], tolerance = 1e-5)

  # The assessment's recruitment penalty keeps its log sigmaR terms and drops only the
  # sqrt(2 pi) constants. The seeded deviations are shifted down by sigmaR^2 / 2 and
  # the bias ramp centres the penalty on the same shift, so the quadratic term is the
  # assessment's own statement on its raw deviations. Two constants come off. The
  # sqrt(2 pi) the assessment omits, and the (1 - bias_ramp / 2) log sigmaR credit
  # SPoRC applies wherever the ramp is active, which the assessment has no analogue
  # for; the ramp is one over the whole series here. The first year is an equilibrium,
  # so there are no initial age deviations to add.
  n_recdev <- dim(input_list$par$ln_RecDevs)[3]
  rec_like <- sum(r$Rec_nLL) + sum(r$Init_Rec_nLL) -
    n_recdev * c2pi + n_recdev * 0.5 * log(dat$sigmaR)
  expect_equal(rec_like, plc[["reclike"]], tolerance = 1e-5)

  expect_equal(r$M_nLL - log(sqrt(2 * pi) * dat$cv_M),
               plc[["prior_M_bin_1"]], tolerance = 1e-5)
  expect_equal(r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q[1]),
               plc[["AI_survey_prior_q_bin_1"]], tolerance = 1e-5)

  # The selectivity penalty is the sum of the assessment's lambdas 3 to 6 on the
  # bicubic surface plus the mean centering term its template hardcodes. SPoRC reports
  # only the total, so the five terms are checked together.
  expect_equal(sum(r$sel_nLL), dlc[["selpen"]], tolerance = 1e-5)

  # The catch statement is numerically zero on both sides against an objective of 890;
  # the assessment fits catch to 1e-04 and SPoRC to 6e-02, which is a difference in
  # how near-exact catch is driven rather than in the statement. Asserted in absolute
  # terms because a relative tolerance on two numbers this close to zero means nothing.
  catch_ssq <- sum(as.vector(r$Catch_nLL) - (c2pi + as.vector(input_list$par$ln_sigmaC)))
  expect_lt(abs(catch_ssq - dlc[["catch.like"]]), 1e-1)

  # The whole objective, like for like. The assessment estimates maturity internally
  # and SPoRC fixes it, so mat_like is removed from the assessment's side.
  like_for_like <- sum(r$FishAgeComps_nLL) + sum(r$FishLenComps_nLL) +
    sum(r$SrvAgeComps_nLL) + sum(r$SrvLenComps_nLL) +
    srv_like(1) + srv_like(2) + catch_ssq + as.vector(f_pen) + rec_like +
    (r$M_nLL - log(sqrt(2 * pi) * dat$cv_M)) +
    (r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q[1])) + sum(r$sel_nLL)
  expect_equal(like_for_like, dlc[["obj_fun"]] - dlc[["mat_like"]], tolerance = 3e-4)

  expect_jnLL_decomposes(obj)
})
