# Self-validating bridge test. The expectations are the 2024 eastern Bering Sea
# Pacific cod assessment's own reported quantities (Stock Synthesis Model 24.1),
# shipped in sgl_rg_ebs_pcod_data$ss3. Every parameter is set to the assessment's
# maximum likelihood estimate and the model is evaluated there without
# optimizing, so a failure means time-varying growth, the cohort propagation of
# size at age, the Richards curve, length-based selectivity, the length bin map
# or the selection-weighted weight at age no longer reproduces the assessment at
# a point where it is known to. Do not loosen the tolerances to make this pass.
# The assessment's report carries six significant digits, which is the floor
# under every comparison here.
#
# Three convention differences are deliberate and are compared after restating
# SPoRC's quantity in the assessment's terms, never by loosening a tolerance:
#   1. SPoRC's single-sex spawning biomass is the female share at an even sex
#      ratio; a one-sex Stock Synthesis model's spawning output counts every
#      mature fish, so the two differ by exactly two.
#   2. SPoRC's deviation penalties carry the full normal constant
#      (log(sigma) + 0.5 log(2 pi) per deviation) where Stock Synthesis carries
#      only the bias-ramp share of log(sigma). With sigmaR and the deviation
#      standard deviations fixed, as they are here, that is a constant.
#   3. SPoRC applies half the bias ramp to the log(sigmaR) term of the
#      recruitment penalty where Stock Synthesis applies the whole ramp.

library(SPoRC)
library(testthat)
data("sgl_rg_ebs_pcod_data")

test_that("EBS Pacific cod bridges to the 2024 Stock Synthesis assessment at its own estimate", {

  dat <- sgl_rg_ebs_pcod_data
  yrs <- dat$years; n_yrs <- length(yrs)
  ages <- dat$ages; n_ages <- length(ages)
  s3 <- dat$ss3

  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)
  obj <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  r <- obj$rep

  pct <- function(a, b) max(abs(100 * (a - b) / b), na.rm = TRUE)
  yr_row <- function(m, y) { rr <- as.integer(rownames(m)); m[as.character(max(rr[rr <= y])), ] }

  # Growth ------------------------------------------------------------------
  # Size at age is carried cohort by cohort from 2000, so a year late in the
  # series only reproduces if every earlier year's increment did.
  for(y in c(1977, 1999, 2000, 2001, 2012, 2023)) {
    iy <- match(y, yrs)
    expect_lt(pct(r$WAA[1, 1, iy, 1, , 1], dat$wtatage[["0"]][as.character(y), ]), 1e-2)      # start-of-year weight
    expect_lt(pct(r$WAA_fish[1, 1, iy, 1, , 1, 1], dat$wtatage[["1"]][as.character(y), ]), 1e-2) # selection-weighted
    # the survey index is in numbers, so its weight at age enters nothing; this is
    # the population weight at the survey's mid-year timing, which tests the key there
    expect_lt(pct(r$WAA_srv[1, 1, iy, 1, , 1, 1], dat$wtatage[["-1"]][as.character(y), ]), 1e-2)
  }
  # the year-by-year growth parameters, deviations applied on the bounded logit scale
  gy <- s3$growth_by_year
  for(y in c(1977, 2000, 2010, 2024)) {
    iy <- match(y, yrs); g <- gy[gy$Yr == y, ]
    expect_lt(pct(r$growth_pars_y[1, 1, iy, 1, 1], g$L1), 1e-3)  # length at the young reference age
    expect_lt(pct(r$growth_pars_y[1, 1, iy, 3, 1], g$K), 1e-3)   # growth rate
    expect_lt(pct(r$growth_pars_y[1, 1, iy, 2, 1], g$L2), 1e-3)  # asymptote, not varying
    expect_lt(pct(r$growth_pars_y[1, 1, iy, 6, 1], g$rho), 1e-3) # Richards coefficient, not varying
  }

  # Selectivity -------------------------------------------------------------
  # Length-based double normal on the population bins, anchored at the first
  # data bin, with two fishery blocks and annual deviations on the survey width
  fish_l <- srv_l <- fish_a <- srv_a <- 0
  for(y in yrs) {
    iy <- match(y, yrs)
    fish_l <- max(fish_l, abs(r$fish_sel_l[1, iy, , 1, 1] - yr_row(s3$lsel[[1]], y)))
    srv_l  <- max(srv_l,  abs(r$srv_sel_l[1, iy, , 1, 1]  - yr_row(s3$lsel[[2]], y)))
    fish_a <- max(fish_a, abs(r$fish_sel[1, 1, iy, 1, , 1, 1] - yr_row(s3$asel2[[1]], y)))
    srv_a  <- max(srv_a,  abs(r$srv_sel[1, 1, iy, 1, , 1, 1]  - yr_row(s3$asel2[[2]], y)))
  }
  expect_lt(fish_l, 1e-5); expect_lt(srv_l, 1e-5)
  # folding the length selectivity to age goes through the age-length key, so
  # these also test the key in every year
  expect_lt(fish_a, 1e-5); expect_lt(srv_a, 1e-5)

  # Population --------------------------------------------------------------
  naa <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
  expect_lt(pct(naa, s3$NAA), 1e-2)
  expect_lt(pct(r$Rec[1, 1, ], s3$Rec), 1e-2)
  # SPoRC's single-sex spawning biomass is the female share of the population
  expect_lt(pct(2 * r$SSB[1, 1, 1:n_yrs], s3$SSB), 1e-2)
  expect_lt(pct(apply(naa * r$WAA[1, 1, 1:n_yrs, 1, , 1], 1, sum), s3$Bio_all), 1e-2)
  expect_lt(pct(r$PredCatch[1, 1, , 1, 1], s3$dead_B), 1e-2)
  cp <- s3$cpue
  expect_lt(pct(r$PredSrvIdx[1, 1, match(cp$Yr, yrs), 1, 1], cp$Exp), 1e-2)

  # Expected compositions ---------------------------------------------------
  # The fishery's are formed at mid season on the season-long catch with the
  # length selectivity applied at length, the survey's at its own timing, and
  # both are mapped from the model's 121 population bins onto the 24 data bins
  LBM <- dat$LenBinMap
  exp_len <- function(v) { w <- as.vector(v %*% LBM); w / sum(w) }
  ld <- s3$lendbase; ad <- s3$agedbase
  fl <- sl <- sa <- 0
  for(y in unique(ld$Yr[ld$Fleet == 1])) {
    e <- exp_len(r$CAL[1, 1, match(y, yrs), 1, , 1, 1])
    fl <- max(fl, abs(e - ld$Exp[ld$Yr == y & ld$Fleet == 1]))
  }
  for(y in unique(ld$Yr[ld$Fleet == 2])) {
    e <- exp_len(r$SrvIAL[1, 1, match(y, yrs), 1, , 1, 1])
    sl <- max(sl, abs(e - ld$Exp[ld$Yr == y & ld$Fleet == 2]))
  }
  for(y in unique(ad$Yr[ad$Fleet == 2])) {
    iy <- match(y, yrs)
    e <- as.vector(r$SrvIAA[1, 1, iy, 1, , 1, 1] %*% dat$AgeingError[iy, , ]); e <- e / sum(e)
    sa <- max(sa, abs(e - ad$Exp[ad$Yr == y & ad$Fleet == 2]))
  }
  expect_lt(fl, 1e-4); expect_lt(sl, 1e-4); expect_lt(sa, 1e-4)

  # Likelihood components ---------------------------------------------------
  # Each is restated in the assessment's convention by removing the normal
  # constants SPoRC carries and Stock Synthesis does not.
  sigmaR <- dat$rec$sigmaR; lsr <- log(sigmaR); c2pi <- 0.5 * log(2 * pi)
  L <- s3$likelihoods; lbf <- s3$likelihoods_by_fleet
  by_fleet <- function(lab, col) lbf[[col]][lbf$Label == lab]

  est_r <- which(!is.na(obj$data$map_ln_RecDevs[1, 1, ]))
  est_i <- which(!is.na(obj$data$map_ln_InitDevs[1, 1, , 1]))
  use_i <- obj$data$init_devs_pen_use[1, 1, , 1]
  ramp_r <- sum(r$bias_ramp[est_r]); ramp_i <- sum(r$init_bias_ramp[est_i] * use_i[est_i])

  # recruitment: SPoRC's half ramp on log(sigmaR) restated as the assessment's whole ramp
  rec <- sum(r$Rec_nLL) - length(est_r) * c2pi - 0.5 * ramp_r * lsr +
    sum(r$Init_Rec_nLL) - sum(use_i[est_i]) * (lsr + c2pi) + (ramp_r + ramp_i) * lsr
  expect_lt(pct(rec, L["Recruitment", "values"]), 1e-2)

  n_idx <- sum(obj$data$UseSrvIdx)
  expect_lt(pct(sum(r$SrvIdx_nLL) - n_idx * c2pi, L["Survey", "values"]), 1e-2)
  expect_lt(pct(sum(r$FishLenComps_nLL), by_fleet("Length_like", "fishery")), 1e-2)
  expect_lt(pct(sum(r$SrvLenComps_nLL), by_fleet("Length_like", "survey")), 1e-2)
  expect_lt(pct(sum(r$SrvAgeComps_nLL), by_fleet("Age_like", "survey")), 1e-2)

  # the initial equilibrium recruitment's offset from the recruitment level
  expect_lt(pct(r$rinit_nLL - (log(obj$data$rinit_pen_sd) + c2pi), L["InitEQ_Regime", "values"]), 1e-2)

  # the two growth deviation series and the survey selectivity deviations
  sg <- c(dat$growth$dev_sd[["L1"]], dat$growth$dev_sd[["K"]]); n_gd <- length(dat$growth$dev_years$L1)
  expect_lt(pct(r$growth_tv_nLL - n_gd * sum(log(sg) + c2pi), sum(s3$parm_devs$Like_devs[1:2])), 1e-2)
  n_sd <- sum(!is.na(obj$data$map_ln_srvsel_devs))
  expect_lt(pct(r$sel_nLL - n_sd * (log(dat$mle$sel$survey$dev_sd) + c2pi), s3$parm_devs$Like_devs[3]), 1e-2)

  # catch is fit essentially exactly on both sides, so compare the kernel to zero
  sc <- exp(as.vector(input_list$par$ln_sigmaC)); obs <- as.vector(dat$ObsCatch); ok <- !is.na(obs)
  expect_lt(sum(r$Catch_nLL) - sum(log(sc[ok]) + c2pi), 1e-4)

  # the gradient at the assessment's estimate is finite everywhere; fishing
  # mortality is conditioned on the catch in the assessment rather than
  # estimated, so its deviations carry the residual that the soft catch
  # likelihood stands in for and the gradient is not expected to be zero
  expect_true(all(is.finite(obj$gr(obj$par))))
})


test_that("a survey can take its weight at age on the fish it selects", {

  # The assessment itself reads the population weight at age for its index, so
  # this is not a bridge comparison: it checks that srv_waa_selected reaches the
  # survey's weight at age the same way fish_waa_selected reaches the fishery's,
  # on a model that has everything the option needs (length-based survey
  # selectivity and weight at age derived from growth).
  dat <- sgl_rg_ebs_pcod_data
  base_input <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)

  sel_input <- base_input
  sel_input$data$srv_waa_selected <- 1

  r0 <- fit_model(base_input$data, base_input$par, base_input$map, do_optim = FALSE, silent = TRUE)$rep
  r1 <- fit_model(sel_input$data, sel_input$par, sel_input$map, do_optim = FALSE, silent = TRUE)$rep

  w_mid <- dat$wtlen[1] * dat$lens^dat$wtlen[2]
  n_yrs <- length(dat$years)

  for(y in c(1, round(n_yrs / 2), n_yrs)) {

    key <- r1$SizeAgeTrans_srv[1, 1, y, 1, , , 1, 1]
    sel <- r1$srv_sel_l[1, y, , 1, 1]
    hand <- as.vector(t(key) %*% (sel * w_mid)) / as.vector(t(key) %*% sel)

    expect_equal(r1$WAA_srv[1, 1, y, 1, , 1, 1], hand, tolerance = 1e-10)
    # and it is a different quantity from the population weight the default reads
    expect_false(isTRUE(all.equal(r1$WAA_srv[1, 1, y, 1, , 1, 1], r0$WAA_srv[1, 1, y, 1, , 1, 1])))

  } # end y loop

  # the fishery's weight is untouched by the survey's flag
  expect_equal(r1$WAA_fish, r0$WAA_fish)
})
