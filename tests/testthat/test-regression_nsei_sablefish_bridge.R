# checks SPoRC against the 2025 NSEI sablefish TMB model at its saved estimate, without optimizing
# the fishery cpue and mark-recapture are fishery indices here (mean numbers over the year), tmb predicts
# them in august, so those two are checked against their equation and everything else against tmb
library(SPoRC)
library(testthat)
data("sgl_rg_nsei_sablefish_data")

test_that("NSEI sablefish reproduces the 2025 TMB assessment at its saved estimate", {

  tmb_data <- sgl_rg_nsei_sablefish_data$inputs
  tmb <- sgl_rg_nsei_sablefish_data$tmb_saved
  n_yrs <- tmb_data$nyr

  # setup model and fix values at the tmb estimate
  input_list <- suppressWarnings(build_nsei_sablefish_input(tmb_data = tmb_data, tmb_par = tmb$par))
  input_list <- seed_nsei_sablefish_tmb(input_list = input_list, tmb_par = tmb$par)

  # keep the estimation map, mapping everything off drops the deviation penalties
  obj <- fit_model(data = input_list$data,
                   parameters = input_list$par,
                   mapping = input_list$map,
                   random = NULL,
                   do_optim = FALSE,
                   silent = TRUE)
  sporc_rep <- obj$rep

  # Compare dynamics (tmb arrays are male first)
  expect_equal(as.vector(sporc_rep$NAA[1,1,1:n_yrs,1,,]), as.vector(tmb$N[,,2:1]), tolerance = 1e-8)
  expect_equal(as.vector(sporc_rep$Rec), tmb$pred_rec, tolerance = 1e-8)
  expect_equal(as.vector(sporc_rep$SSB), tmb$ssb, tolerance = 1e-8)
  expect_equal(sporc_rep$Fmort[1,,1,1], tmb$Fmort, tolerance = 1e-8)
  expect_equal(as.vector(sporc_rep$ZAA[1,1,,1,,]), as.vector(tmb$Z[,,2:1]), tolerance = 1e-8)
  expect_equal(as.vector(sporc_rep$fish_sel[1,1,1:n_yrs,1,,,1]), as.vector(tmb$fsh_slx[,,2:1]), tolerance = 1e-8)
  expect_equal(as.vector(sporc_rep$srv_sel[1,1,1:n_yrs,1,,,1]), as.vector(tmb$srv_slx[,,2:1]), tolerance = 1e-8)

  # the mark-recapture fleet doesn't fish and has the fishery selectivity
  expect_equal(max(sporc_rep$Fmort[1,,1,2]), 0)
  expect_equal(sporc_rep$fish_sel[,,,,,,2], sporc_rep$fish_sel[,,,,,,1])

  # Compare catch (kg in SPoRC, tonnes in tmb) and the longline survey index
  expect_equal(sporc_rep$PredCatch[1,1,,1,1] / 1e3, tmb$pred_landed, tolerance = 1e-8)
  expect_equal(sporc_rep$PredSrvIdx[1,1,tmb_data$yrs_srv_cpue + 1,1,1], tmb$pred_srv_cpue, tolerance = 1e-8)

  # Fishery cpue and mark-recapture, q times the mean numbers over the year times fishery selectivity
  # and retention (times fishery weight at age for the fishery cpue)
  fsh_cpue_yrs <- tmb_data$yrs_fsh_cpue + 1
  mr_yrs <- tmb_data$yrs_mr + 1
  pred_fsh_cpue <- rep(0, length(fsh_cpue_yrs))
  pred_mr <- rep(0, length(mr_yrs))
  for(i in 1:length(fsh_cpue_yrs)) {
    y <- fsh_cpue_yrs[i]
    mean_naa <- sporc_rep$NAA[1,1,y,1,,] * (1 - exp(-sporc_rep$ZAA[1,1,y,1,,])) / sporc_rep$ZAA[1,1,y,1,,]
    pred_fsh_cpue[i] <- sporc_rep$fish_q[1,y,1] * sum(mean_naa * sporc_rep$fish_sel[1,1,y,1,,,1] * sporc_rep$ret_sel[1,1,y,1,,,1] * input_list$data$WAA_fish[1,1,y,1,,,1])
  } # end i loop
  for(i in 1:length(mr_yrs)) {
    y <- mr_yrs[i]
    mean_naa <- sporc_rep$NAA[1,1,y,1,,] * (1 - exp(-sporc_rep$ZAA[1,1,y,1,,])) / sporc_rep$ZAA[1,1,y,1,,]
    pred_mr[i] <- sporc_rep$fish_q[1,y,2] * sum(mean_naa * sporc_rep$fish_sel[1,1,y,1,,,2] * sporc_rep$ret_sel[1,1,y,1,,,2])
  } # end i loop
  expect_equal(sporc_rep$PredFishIdx[1,1,fsh_cpue_yrs,1,1], pred_fsh_cpue, tolerance = 1e-10)
  expect_equal(sporc_rep$PredFishIdx[1,1,mr_yrs,1,2], pred_mr, tolerance = 1e-10)

  # the mean over the year is 2 - 3% above the august numbers tmb uses
  expect_true(all(pred_mr / tmb$pred_mr > 1.02 & pred_mr / tmb$pred_mr < 1.03))

  # Compare age comps
  # tmb ageing error rows don't sum to 1 and tmb doesn't renormalize after, so renormalize tmb here
  fsh_age_yrs <- tmb_data$yrs_fsh_age + 1
  srv_age_yrs <- tmb_data$yrs_srv_age + 1
  sporc_fsh_age <- array(0, dim = dim(tmb$pred_fsh_age))
  sporc_srv_age <- array(0, dim = dim(tmb$pred_srv_age))
  for(i in 1:length(fsh_age_yrs)) {
    caa <- rowSums(sporc_rep$CAA[1,1,fsh_age_yrs[i],1,,,1]) # sexes combined
    caa <- as.vector((caa / sum(caa)) %*% tmb_data$ageing_error) # ageing error
    sporc_fsh_age[i,] <- caa / sum(caa)
  } # end i loop
  for(i in 1:length(srv_age_yrs)) {
    iaa <- rowSums(sporc_rep$SrvIAA[1,1,srv_age_yrs[i],1,,,1])
    iaa <- as.vector((iaa / sum(iaa)) %*% tmb_data$ageing_error)
    sporc_srv_age[i,] <- iaa / sum(iaa)
  } # end i loop
  expect_equal(sporc_fsh_age, tmb$pred_fsh_age / rowSums(tmb$pred_fsh_age), tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(sporc_srv_age, tmb$pred_srv_age / rowSums(tmb$pred_srv_age), tolerance = 1e-8, ignore_attr = TRUE)

  # Compare length comps (by sex)
  sporc_fsh_len <- array(0, dim = dim(tmb$pred_fsh_len))
  sporc_srv_len <- array(0, dim = dim(tmb$pred_srv_len))
  for(s in 1:2) {
    for(i in 1:length(tmb_data$yrs_fsh_len)) {
      cal <- sporc_rep$CAL[1,1,tmb_data$yrs_fsh_len[i] + 1,1,,s,1]
      sporc_fsh_len[i,,3 - s] <- cal / sum(cal)
    } # end i loop
    for(i in 1:length(tmb_data$yrs_srv_len)) {
      ial <- sporc_rep$SrvIAL[1,1,tmb_data$yrs_srv_len[i] + 1,1,,s,1]
      sporc_srv_len[i,,3 - s] <- ial / sum(ial)
    } # end i loop
  } # end s loop
  expect_equal(sporc_fsh_len, tmb$pred_fsh_len, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(sporc_srv_len, tmb$pred_srv_len, tolerance = 1e-4, ignore_attr = TRUE)

  # Compare comp likelihoods, tmb multinomial (1e-5 inside each log) on the SPoRC predictions
  tmb_mltnml <- function(obs, pred, n) {
    obs <- obs / rowSums(obs)
    sum(n * rowSums((obs + 1e-5) * (log(obs + 1e-5) - log(pred + 1e-5))))
  } # end function
  expect_equal(sum(sporc_rep$FishAgeComps_nLL), tmb_mltnml(tmb_data$data_fsh_age, sporc_fsh_age, tmb_data$effn_fsh_age), tolerance = 1e-8)
  expect_equal(sum(sporc_rep$SrvAgeComps_nLL), tmb_mltnml(tmb_data$data_srv_age, sporc_srv_age, tmb_data$effn_srv_age), tolerance = 1e-8)
  fsh_len_like <- 0
  srv_len_like <- 0
  for(k in 1:2) {
    fsh_len_like <- fsh_len_like + tmb_mltnml(tmb_data$data_fsh_len[,,k], sporc_fsh_len[,,k], tmb_data$effn_fsh_len[,1,k])
    srv_len_like <- srv_len_like + tmb_mltnml(tmb_data$data_srv_len[,,k], sporc_srv_len[,,k], tmb_data$effn_srv_len[,1,k])
  } # end k loop
  expect_equal(sum(sporc_rep$FishLenComps_nLL), fsh_len_like, tolerance = 1e-8)
  expect_equal(sum(sporc_rep$SrvLenComps_nLL), srv_len_like, tolerance = 1e-8)

  # Compare other likelihoods
  # SPoRC normal densities keep log(sd) + log(2 pi) / 2, tmb sums of squares don't, so take that off
  norm_const <- function(sd) sum(log(sd) + 0.5 * log(2 * pi))
  sigmaR_sporc <- exp(tmb$par$log_sigma_r) / sqrt(2)
  n_init_devs <- length(tmb$par$log_rinit_devs)

  expect_equal(sum(sporc_rep$Catch_nLL) - norm_const(tmb_data$sigma_catch), tmb$catch_like, tolerance = 1e-6)
  expect_equal(sum(sporc_rep$SrvIdx_nLL) - norm_const(tmb_data$sigma_srv_cpue), tmb$index_like[2], tolerance = 1e-8)
  expect_equal(sum(sporc_rep$Rec_nLL) + sum(sporc_rep$Init_Rec_nLL) - norm_const(rep(sigmaR_sporc, n_yrs + n_init_devs)), tmb$rec_like, tolerance = 1e-8)
  expect_equal(sum(sporc_rep$Fmort_nLL) - norm_const(rep(sqrt(1 / (2 * tmb_data$wt_fpen)), n_yrs)), tmb$fpen, tolerance = 1e-8)
  expect_equal(sporc_rep$fish_q_nLL - norm_const(tmb_data$sigma_mr_q), tmb$priors[3], tolerance = 1e-8)

  # fishery cpue and mark-recapture, lognormal with 1e-5 inside each log
  fsh_cpue_like <- -sum(dnorm(log(tmb_data$data_fsh_cpue + 1e-5), log(pred_fsh_cpue + 1e-5), tmb_data$sigma_fsh_cpue, log = TRUE))
  mr_like <- -sum(dnorm(log(tmb_data$data_mr + 1e-5), log(pred_mr + 1e-5), tmb_data$sigma_mr, log = TRUE))
  expect_equal(sum(sporc_rep$FishIdx_nLL[,,,1]), fsh_cpue_like, tolerance = 1e-8)
  expect_equal(sum(sporc_rep$FishIdx_nLL[,,,2]), mr_like, tolerance = 1e-8)

  # Compare SPR reference points (F35 - F70 and spawning biomass per recruit)
  for(x in 1:length(tmb_data$Fxx_levels)) {
    ref_pts <- Get_Reference_Points(data = input_list$data,
                                    rep = sporc_rep,
                                    SPR_x = tmb_data$Fxx_levels[x],
                                    t_spawn = tmb_data$spawn_month,
                                    calc_rec_yrs = 1:n_yrs,
                                    type = "single_region",
                                    what = "SPR",
                                    n_avg_yrs = 1)
    expect_equal(as.vector(ref_pts$f_ref_pt), tmb$Fxx[x + 1], tolerance = 1e-8)
    expect_equal(as.vector(ref_pts$b_ref_pt / ref_pts$mean_rec), 0.5 * tmb$SBPR[x + 1], tolerance = 1e-8) # tmb sbpr is both sexes
  } # end x loop

  expect_jnLL_decomposes(obj)
})
