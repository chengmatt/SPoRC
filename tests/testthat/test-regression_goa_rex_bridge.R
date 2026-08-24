# Self-validating bridge test. The expectations are the 2025 Gulf of Alaska rex sole
# assessment's own reported quantities (Stock Synthesis Model 25.1), shipped in
# mlt_rg_goa_rex_data$ss3. Every parameter is set to the assessment's maximum
# likelihood estimate and the model is evaluated there without optimizing, so a
# failure means the growth module, the conditional age-at-length likelihood, the
# two-area dynamics, a selectivity form or the ageing error no longer reproduces the
# assessment at a point where it is known to. Do not loosen the tolerances to make
# this pass. The assessment's report and parameter values carry six significant
# digits, which is the floor under every comparison here.

library(SPoRC)
library(testthat)
data("mlt_rg_goa_rex_data")

test_that("GOA rex sole bridges to the 2025 Stock Synthesis assessment at its own estimate", {

  dat <- mlt_rg_goa_rex_data
  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_reg <- dat$n_regions
  n_sex <- dat$n_sexes

  input_list <- seed_goa_rex_mle(suppressWarnings(suppressMessages(build_goa_rex_input(dat))), dat)
  obj <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  r <- obj$rep

  pct <- function(a, b) max(abs(100 * (a - b) / b), na.rm = TRUE)

  # Growth ----
  # Mean length, its spread and weight at age, at the start of the year and at mid
  # season, per area and sex. Linear growth below A1, the plus-group adjustment and
  # the CV interpolation all have to be right for these to match.
  for(a in 1:n_reg) for(s in 1:n_sex) {
    g <- dat$ss3$growth[[a]][[s]]
    expect_lt(pct(r$mean_LAA_spawn[1, a, 1, 1, , s], g$Len_Beg), 1e-2)
    expect_lt(pct(r$mean_LAA_srv[1, a, 1, 1, , s, 1], g$Len_Mid), 1e-2)
    expect_lt(pct(r$sd_LAA_srv[1, a, 1, 1, , s, 1], g$SD_Mid), 1e-2)
    expect_lt(pct(r$WAA[1, a, 1, 1, , s], g$Wt_Beg), 1e-2)
    expect_lt(pct(r$WAA_fish[1, a, 1, 1, , s, 1], g$Wt_Mid), 1e-2)
  } # end a, s loops

  # The mid-season age-length key, which SS3 prints with its lengths largest first
  alk_ss3 <- dat$ss3$ALK[rev(seq_len(dim(dat$ss3$ALK)[1])), , "Seas: 1 Sub_Seas: 2 Morph: 1"]
  expect_lt(max(abs(r$SizeAgeTrans[1, 1, 1, 1, , , 1] - alk_ss3)), 1e-5)

  # Population dynamics ----
  for(a in 1:n_reg) for(s in 1:n_sex) expect_lt(pct(r$NAA[1, a, 1:n_yrs, 1, , s], dat$ss3$NAA[a, , , s]), 1e-2)
  for(a in 1:n_reg) {
    expect_lt(pct(r$SSB[1, a, 1:n_yrs], dat$ss3$SSB[, a]), 1e-2)
    expect_lt(pct(r$Rec[1, a, 1:n_yrs], dat$ss3$Rec[, a]), 1e-2)
    totb <- sapply(1:n_yrs, function(y) sum(r$NAA[1, a, y, 1, , ] * r$WAA[1, a, y, 1, , ]))
    expect_lt(pct(totb, dat$ss3$Bio_all[, a]), 1e-2)
  } # end a loop

  # The bias ramp, built from the four breakpoints in deviation index space, for the
  # model years and for the years before the first one that the initial ages were
  # born in. A wrong ramp is silent in every other quantity, so it is checked outright.
  main_yrs <- as.integer(names(dat$mle$main_recdev))
  expect_equal(as.vector(r$bias_ramp)[match(main_yrs, yrs)], as.vector(dat$mle$biasadj[as.character(main_yrs)]), tolerance = 1e-6)
  early_yrs <- as.integer(names(dat$mle$early_recdev))
  expect_equal(as.vector(r$init_bias_ramp)[yrs[1] - early_yrs], as.vector(dat$mle$biasadj[as.character(early_yrs)]), tolerance = 1e-6)

  # Predicted observations ----
  expect_lt(pct(r$PredCatch[1, 1, , 1, 1], dat$ss3$dead_B), 1e-2)
  for(sf in 1:2) {
    ci <- dat$ss3$cpue[dat$ss3$cpue$Fleet == sf + 1, ]
    expect_lt(pct(r$PredSrvIdx[1, sf, match(ci$Yr, yrs), 1, sf], ci$Exp), 1e-2)
  } # end sf loop

  # Selectivity: the assessment's double normal parameters go straight in, with the
  # raw ascending limb its -999 convention means and the male offsets on the fishery
  for(s in 1:n_sex) expect_lt(max(abs(r$fish_sel[1, 1, 1, 1, , s, 1] - dat$ss3$sel[[1]][[s]])), 1e-5)
  for(sf in 1:2) for(s in 1:n_sex) expect_lt(max(abs(r$srv_sel[1, sf, 1, 1, , s, sf] - dat$ss3$sel[[sf + 1]][[s]])), 1e-5)

  # Conditional age-at-length: one length bin of one survey year, formed the way the
  # likelihood forms it (the joint row conditioned on its length, ageing error, then
  # the composition constant) against the assessment's own table
  cb <- dat$ss3$condbase
  cr <- cb[cb$Fleet == 2 & cb$Yr == 1993 & cb$Lbin_lo == 13 & cb$Sex == 1, ]
  cr <- cr[order(cr$Bin), ]
  l <- match(13, dat$lens_lower)
  e <- r$Srv_caal[1, 1, match(1993, yrs), 1, l, , 1, 1]
  e <- as.vector((e / sum(e)) %*% dat$AgeingError)
  e <- e / sum(e)
  e <- (e + dat$comp$addtocomp_age) / sum(e + dat$comp$addtocomp_age)
  expect_lt(max(abs(e - cr$Exp)), 1e-4)

  # Likelihoods ----
  # The composition likelihoods agree up to SS3 renormalizing after adding its
  # constant, which scales each by (1 + n_bins * addtocomp); the index likelihood
  # up to the 0.5 log(2 pi) SPoRC's normal density carries per observation.
  c_age <- dat$comp$addtocomp_age
  expect_equal(sum(r$FishLenComps_nLL) + sum(r$SrvLenComps_nLL),
               dat$ss3$likelihoods["Length_comp", "values"] * (1 + 2 * length(dat$lens) * c_age), tolerance = 2e-4)
  expect_equal(sum(r$FishAgeComps_nLL) + sum(r$Srv_caal_nLL),
               dat$ss3$likelihoods["Age_comp", "values"] * (1 + length(dat$obs_ages) * c_age), tolerance = 5e-4)
  n_idx <- sum(dat$UseSrvIdx)
  expect_equal(sum(r$SrvIdx_nLL), dat$ss3$likelihoods["Survey", "values"] + n_idx * 0.5 * log(2 * pi), tolerance = 5e-3)

  # The catchability prior is the assessment's normal on the log scale
  expect_equal(as.numeric(sum(r$srv_q_nLL)) - 0.5 * log(2 * pi * dat$q$prior_sd^2),
               dat$ss3$likelihoods["Parm_priors", "values"], tolerance = 1e-3)

  # The recruitment penalty: SS3 writes sum(dev^2) / (2 sigma^2) + sum(b) log(sigma)
  # over the early and main deviations together. SPoRC's two arrays carry the normal
  # constants, the main deviations' (1 - b/2) log(sigma) adjustment, and one log(sigma)
  # per initial age, so the crosswalk subtracts exactly those. A shared deviation has
  # to be counted once for this to close; counting it per region or per sex fails it.
  sig <- dat$rec$sigmaR; n_main <- length(main_yrs); n_early <- length(early_yrs)
  b_main <- as.vector(r$bias_ramp)[seq_len(n_main)]; b_early <- as.vector(r$init_bias_ramp)[seq_len(n_early)]
  rec_net <- sum(r$Rec_nLL) + sum(r$Init_Rec_nLL) - (n_main + n_early) * 0.5 * log(2 * pi) +
    (0.5 * sum(b_main) + sum(b_early) - n_early) * log(sig)
  expect_equal(rec_net, dat$ss3$likelihoods["Recruitment", "values"], tolerance = 1e-4)
})
