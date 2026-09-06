# Self-validating bridge test: the expectations are the 2025 West Coast sablefish assessment's own
# quantities in sgl_rg_wc_sablefish_data$ss3, evaluated at its estimate without optimizing.
#
# SS3's report file has six significant digits, the floor under every comparison here.

library(SPoRC)
library(testthat)
data("sgl_rg_wc_sablefish_data")

test_that("West Coast sablefish bridges to the 2025 Stock Synthesis assessment at its own estimate", {

  dat <- sgl_rg_wc_sablefish_data
  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_srv <- dat$n_srv_fleets

  input_list <- seed_wc_sablefish_mle(build_wc_sablefish_input(dat), dat)
  obj <- fit_model(input_list$data, input_list$par, input_list$map,
                   do_optim = FALSE, silent = TRUE)
  r <- obj$rep

  pct <- function(a, b) max(abs(100 * (a - b) / b))

  # Population dynamics ----
  expect_lt(pct(as.vector(r$NAA[1, 1, 1:n_yrs, 1, , ]), as.vector(dat$ss3$NAA)), 1e-2)
  expect_lt(pct(as.vector(r$SSB[1, 1, 1:n_yrs]), dat$ss3$SSB), 1e-2)
  expect_lt(pct(as.vector(r$Rec[1, 1, 1:n_yrs]), dat$ss3$Rec), 1e-2)
  totb_jan1 <- sapply(1:n_yrs, function(y) sum(r$NAA[1, 1, y, 1, , ] * dat$WAA[1, 1, y, 1, , ]))
  expect_lt(pct(totb_jan1, dat$ss3$Bio_all), 1e-2)

  # The bias ramp is built from the four breakpoints in deviation index space, and
  # a wrong one is silent in every other quantity, so it is checked outright
  expect_equal(as.vector(r$bias_ramp), dat$mle$bias_adj, tolerance = 1e-8)

  # Predicted observations ----
  i_catch <- dat$UseCatch[1, , 1, 1:6] == 1
  expect_lt(pct(as.vector(r$PredCatch[1, 1, , 1, 1:6])[i_catch], as.vector(dat$ss3$pred_catch)[i_catch]), 1e-2)
  # fleets 1-4 hold the trawl survey indices; fleet 5 is compositions only and
  # fleet 6 the recruitment index, checked separately below
  for(sf in 1:4) {
    ci <- dat$ss3$pred_idx[dat$ss3$pred_idx$Fleet == dat$srv_src[sf], ]
    expect_lt(pct(r$PredSrvIdx[1, 1, match(ci$Yr, yrs), 1, sf], ci$Exp), 1e-2)
  } # end sf loop

  # Selectivity. The assessment's parameters go in as starting values on SPoRC's own
  # double normal, so its whole surface should come back, blocks, mirrored fleets and
  # male offsets included.
  sel_fish <- expand_wc_sablefish_sel(dat$fish_sel_blocks_ss3, n_yrs, n_ages, dat$n_sexes)
  sel_srv <- expand_wc_sablefish_sel(dat$srv_sel_blocks_ss3, n_yrs, n_ages, dat$n_sexes)
  expect_lt(max(abs(r$fish_sel[1, 1, 1:n_yrs, 1, , , ] - sel_fish[1, 1, , 1, , , ])), 1e-5)
  # the recruitment index fleet observes the deviations and reads no curve
  i_sel_srv <- which(dat$srv_src != 11)
  expect_lt(max(abs(r$srv_sel[1, 1, 1:n_yrs, 1, , , i_sel_srv] - sel_srv[1, 1, , 1, , , i_sel_srv])), 1e-5)

  # Expected age compositions, formed the way the likelihood forms them and put through
  # the composition constant so they compare with the assessment's own table
  exp_comp <- function(p_f, p_m, joint) {
    if(joint) {
      p <- c(p_f, p_m) / sum(p_f + p_m)
      e <- c(as.vector(p[1:n_ages] %*% dat$AgeingError), as.vector(p[n_ages + 1:n_ages] %*% dat$AgeingError))
    } else e <- as.vector(((p_f + p_m) / sum(p_f + p_m)) %*% dat$AgeingError)
    e <- e / sum(e)
    (e + dat$addtocomp) / sum(e + dat$addtocomp)
  }
  comp_gap <- function(src, sexcode, numbers) {
    max(sapply(seq_along(src), function(f) {
      db <- dat$ss3$agedbase[dat$ss3$agedbase$Fleet == src[f] & dat$ss3$agedbase$Sexes == sexcode[f], ]
      if(nrow(db) == 0) return(0) # a fleet with no compositions, such as the recruitment index
      max(sapply(sort(unique(db$Yr)), function(y) {
        e <- exp_comp(numbers(match(y, yrs), 1, f), numbers(match(y, yrs), 2, f), sexcode[f] == 3)
        dby <- db[db$Yr == y, ]
        max(abs(e - dby$Exp[order(dby$Sex, dby$Bin)]))
      }))
    }))
  }
  expect_lt(comp_gap(dat$fish_src, dat$fish_sex, function(y, s, f) r$CAA[1, 1, y, 1, , s, f]), 1e-5)
  expect_lt(comp_gap(dat$srv_src, dat$srv_sex, function(y, s, f) r$SrvIAA[1, 1, y, 1, , s, f]), 1e-5)

  # Likelihoods. SPoRC evaluates proper densities where the assessment drops
  # normalizing constants, so each Gaussian block is compared net of a closed-form
  # offset. A change to any of those constants shows up here.
  lc <- function(sigma, n) n * (log(sigma) + 0.5 * log(2 * pi))
  n_est_dev <- sum(yrs %in% dat$yrs_rec_est)

  expect_equal(sum(r$Catch_nLL) - lc(dat$sigmaC, sum(dat$UseCatch == 1, na.rm = TRUE)),
               dat$ss3$lik$catch, tolerance = 1e-4)
  # fleets 1-4 are the trawl surveys and 6 the recruitment index; fleet 5
  # has compositions only. The assessment reports all five in one row.
  for(sf in c(1:4, 6)) {
    expect_equal(sum(r$SrvIdx_nLL[, , , sf]) - 0.5 * log(2 * pi) * sum(dat$UseSrvIdx[, , , sf]),
                 dat$ss3$lik$index[min(sf, 5)], tolerance = 1e-3)
  } # end sf loop
  # the trawl fleet's two composition data sources are one fleet in the assessment
  expect_equal(sum(r$FishAgeComps_nLL[, , , , c(1, dat$n_fish_fleets)]), dat$ss3$lik$age[1], tolerance = 1e-3)
  for(f in 2:6) expect_equal(sum(r$FishAgeComps_nLL[, , , , f]), dat$ss3$lik$age[f], tolerance = 1e-3)
  for(sf in 1:3) expect_equal(sum(r$SrvAgeComps_nLL[, , , , sf]), dat$ss3$lik$age[6 + sf], tolerance = 1e-3)
  expect_equal(sum(r$SrvAgeComps_nLL[, , , , 4:5]), dat$ss3$lik$age[10], tolerance = 1e-3)
  expect_equal(sum(r$Rec_nLL) - n_est_dev * 0.5 * log(2 * pi) - 0.5 * sum(dat$mle$bias_adj) * log(dat$sigmaR),
               dat$ss3$lik$recruitment + dat$ss3$lik$forecast_recruitment - sum(dat$mle$bias_adj) * log(dat$sigmaR),
               tolerance = 1e-3)
  expect_equal(r$M_nLL - lc(dat$M_prior$sd, 1), dat$ss3$lik$priors, tolerance = 1e-4)
  # the recruitment index reads the deviations rather than the population, so
  # its predicted values should be the assessment's own recruitment deviations
  # times its catchability
  i_ri <- match(dat$rec_idx$yr, yrs)
  expect_equal(as.vector(r$RecDev_anom[1, 1, i_ri]), dat$mle$recdev[i_ri], tolerance = 1e-8)
  expect_equal(as.vector(r$PredSrvIdx[1, 1, i_ri, 1, n_srv]),
               (dat$ss3$pred_idx$Exp[dat$ss3$pred_idx$Fleet == 11]), tolerance = 1e-5)

  # and the objective itself, net of every constant the assessment omits
  const <- lc(dat$sigmaC, sum(dat$UseCatch == 1, na.rm = TRUE)) +
    0.5 * log(2 * pi) * sum(dat$UseSrvIdx) + n_est_dev * 0.5 * log(2 * pi) -
    0.5 * sum(dat$mle$bias_adj) * log(dat$sigmaR) + lc(dat$M_prior$sd, 1)
  expect_equal(obj$fn(obj$par) - const, dat$mle$objective, tolerance = 1e-6)

  # terms the assessment does not have should not be contributing
  for(quant_name in c("Fmort_nLL", "sel_nLL", "srv_q_nLL", "fish_q_nLL", "h_nLL", "R0_nLL",
              "SR_pen_nLL", "Init_Rec_nLL", "FishIdx_nLL")) {
    expect_equal(sum(r[[quant_name]]), 0, tolerance = 1e-10)
  } # end quant_name loop
})
