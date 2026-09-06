# Bridge test for the sandeel case study (vignette("ah_north_sea_sandeel_case_study")). Every parameter
# sits at smsR's own MLE and the model is evaluated there without optimizing.
#
# This stock is why seasonal M exists in SPoRC: the season-1 share of annual M runs 0.38 to 0.55, so no
# single annual rate gets it. The seeded comparison is an identity, not a fit, so it holds to floating point.

library(SPoRC)
library(testthat)

sandeel <- build_sandeel_bridge()

test_that("the smsR seasonal mortality is converted to a rate per year", {

  skip_if(is.null(sandeel), "the dev tree with the sandeel bundle is not present")

  M_sms <- sandeel$M_sms
  seasdur <- sandeel$seasdur
  natmort <- sandeel$natmort

  expect_equal(dim(natmort), c(1, 1, sandeel$n_yrs, 2, 5, 1))

  # smsR stores mortality accumulated in a season, SPoRC stores the rate, so the
  # rate times the duration has to give back smsR's own number
  for(seas in 1:2) {
    applied <- natmort[1, 1, , seas, , 1] * seasdur[seas]
    expect_equal(applied, t(M_sms[, , seas]), tolerance = 1e-12,
                 label = paste("applied mortality, season", seas))
  } # end seas loop

  # the two seasons really do differ, so this isn't a constant rate with a season
  # dim bolted on
  share <- M_sms[2:5, , 1] / (M_sms[2:5, , 1] + M_sms[2:5, , 2])
  expect_gt(max(share) - min(share), 0.15)
  expect_false(isTRUE(all.equal(natmort[1, 1, , 1, , 1], natmort[1, 1, , 2, , 1])))

  # age 0 doesn't exist until season 2 and has no season 1 mortality, so nothing
  # needs correcting on the way in
  expect_equal(unname(natmort[1, 1, , 1, 1, 1]), rep(0, sandeel$n_yrs))
  expect_true(all(natmort[1, 1, , 2, 1, 1] > 0))
})


test_that("North Sea sandeel bridges to smsR at its own estimate", {

  skip_if(is.null(sandeel), "the dev tree with the sandeel bundle is not present")

  rep <- sandeel$rep
  n_yrs <- sandeel$n_yrs

  # fishing mortality at age, every year, season and age cell
  expect_lt(max(sandeel$rel_F), 1e-12)

  # SSB, recruitment and numbers at age. Seeded rather than fitted, so these are
  # the population dynamics rather than the likelihood
  pct <- function(sp, sms) max(abs(100 * (sp - sms) / sms))

  expect_lt(pct(as.numeric(rep$SSB)[1:n_yrs], sandeel$ref$SSB[1:n_yrs]), 1e-8)
  expect_lt(pct(as.numeric(rep$Rec)[1:n_yrs], sandeel$ref$Rsave[1:n_yrs]), 1e-10)

  sp_N <- rep$NAA[1, 1, 1:n_yrs, 1, , 1]
  for(a in 2:5) expect_lt(pct(sp_N[, a], sandeel$ref_N[, a]), 1e-8)

  sp_F <- apply(rep$tot_FAA[1, 1, 1:n_yrs, , , 1, ], c(1, 3), sum)
  expect_lt(pct(rowMeans(sp_F[, 2:3]), rowMeans(sandeel$ref_F[, 2:3])), 1e-10)
})


test_that("total mortality within a season is the seasonal rate times its duration", {

  skip_if(is.null(sandeel), "the dev tree with the sandeel bundle is not present")

  # the thing the whole feature rests on, read off the fitted object rather than
  # off the inputs
  rep <- sandeel$rep
  n_yrs <- sandeel$n_yrs
  FAA <- apply(rep$ret_FAA, 1:6, sum) + apply(rep$disc_FAA, 1:6, sum)

  for(seas in 1:2) {
    M_applied <- rep$ZAA[1, 1, 1:n_yrs, seas, , 1] - FAA[1, 1, 1:n_yrs, seas, , 1]
    expect_equal(M_applied, t(sandeel$M_sms[, , seas]), tolerance = 1e-10,
                 label = paste("mortality folded into Z, season", seas))
  } # end seas loop
})
