library(SPoRC)
library(testthat)
data("sgl_rg_ebswp_data")

test_that("Single-region EBS Pollock (nonpar selex + BTS covariance) RTMB model produces expected results", {

  ## Initialize model dimensions and data list ----
  input_list <- Setup_Mod_Dim(
    years = sgl_rg_ebswp_data$years,
    ages = sgl_rg_ebswp_data$ages,
    lens = NA,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 3,
    n_seas = sgl_rg_ebswp_data$n_seas,
    n_pop = sgl_rg_ebswp_data$n_pop,
    natal_region = sgl_rg_ebswp_data$natal_region,
    verbose = FALSE
  )

  inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

  # Recruitment ----
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(c(log(sqrt(1 / (2 * 0.1))), log(1)),
                      dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    rec_model = "bh_rec",
    steepness_h = array(inv_steepness(0.67),
                        dim = c(input_list$data$n_pop, input_list$data$n_regions)),
    h_spec = "est_shared_pop_r",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    ln_global_R0 = 10,
    t_spawn = 0.25,
    equil_init_age_strc = 2,
    use_rinit = 1,
    dont_est_recdev_last = 2,
    Use_h_prior = 1,
    h_prior = data.frame(pop = 1, region = 1, mu = 0.6, sd = 0.12)
  )

  # Fixed natural mortality ----
  fix_natmort <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                  length(input_list$data$years),
                                  length(input_list$data$ages), 1))
  fix_natmort[, , , 1, ]        <- 0.9  # age 1 M
  fix_natmort[, , , 2, ]        <- 0.45 # age 2 M
  fix_natmort[, , , -c(1, 2), ] <- 0.3  # age 3+ M

  suppressWarnings(
    input_list <- Setup_Mod_Biologicals(
      input_list = input_list,
      WAA = sgl_rg_ebswp_data$WAA,
      WAA_fish = sgl_rg_ebswp_data$WAA_fish,
      WAA_srv = sgl_rg_ebswp_data$WAA_srv,
      MatAA = sgl_rg_ebswp_data$MatAA,
      fit_lengths = 0,
      M_spec = "fix",
      Fixed_natmort = fix_natmort,
      addtocomp = 0.001,
      addtosrvidx = 0
    )
  )

  # Movement ----
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  # Catch and F ----
  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list,
      ObsCatch = sgl_rg_ebswp_data$ObsCatch,
      UseCatch = sgl_rg_ebswp_data$UseCatch,
      Use_F_pen = 1,
      sigmaC_spec = "fix",
      ln_sigmaC = array(log(0.05),
                        dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1))
    )
  )

  # Fishery index and comps ----
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sgl_rg_ebswp_data$ObsFishIdx,
    ObsFishIdx_SE = sgl_rg_ebswp_data$ObsFishIdx_SE,
    UseFishIdx = sgl_rg_ebswp_data$UseFishIdx,
    ObsFishAgeComps = sgl_rg_ebswp_data$ObsFishAgeComps,
    UseFishAgeComps = sgl_rg_ebswp_data$UseFishAgeComps,
    ISS_FishAgeComps = sgl_rg_ebswp_data$ISS_FishAgeComps,
    ObsFishLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, 1)),
    UseFishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1)),
    ISS_FishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, 1)),
    fish_idx_type = c("biom"),
    FishAgeComps_LikeType = c("Multinomial"),
    FishLenComps_LikeType = c("none"),
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  # Survey indices and comps ----
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sgl_rg_ebswp_data$ObsSrvIdx,
    ObsSrvIdx_SE = sgl_rg_ebswp_data$ObsSrvIdx_SE,
    UseSrvIdx = sgl_rg_ebswp_data$UseSrvIdx,
    ObsSrvAgeComps = sgl_rg_ebswp_data$ObsSrvAgeComps,
    ISS_SrvAgeComps = sgl_rg_ebswp_data$ISS_SrvAgeComps,
    UseSrvAgeComps = sgl_rg_ebswp_data$UseSrvAgeComps,
    ObsSrvLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, 3)),
    UseSrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, 3)),
    ISS_SrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, 3)),
    srv_idx_type = c("biom", "biom", "biom"),
    SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "Multinomial"),
    SrvLenComps_LikeType = c("none", "none", "none"),
    SrvAgeComps_Type = c(
      "agg_Year_1-terminal_Fleet_1",
      "agg_Year_1-terminal_Fleet_2",
      "none_Year_1-terminal_Fleet_3"
    ),
    SrvLenComps_Type = c(
      "none_Year_1-terminal_Fleet_1",
      "none_Year_1-terminal_Fleet_2",
      "none_Year_1-terminal_Fleet_3"
    ),
    t_srv = array(c(0.5, 0.5, 1),
                  dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_srv_fleets))
  )

  # Fishery selectivity and catchability (nonparametric, random walk) ----
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = c("rw_Fleet_1"),
    fishsel_pe_pars_spec = "fix",
    fish_sel_devs_spec = "est_all",
    fish_sel_blocks = c("none_Fleet_1"),
    fish_sel_model = c("nonpar_Fleet_1"),
    fish_q_blocks = c("none_Fleet_1"),
    fish_fixed_sel_pars_spec = c("est_all"),
    fish_q_spec = c("est_all"),
    fish_sel_nonpar_est_bins = list(list(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12:15)),
    fishsel_devs_shared_bins = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12:15)
  )

  input_list$data$map_ln_fishsel_devs[, , 13, , ] <- input_list$data$map_ln_fishsel_devs[, , 12, , ]
  input_list$data$map_ln_fishsel_devs[, , 14, , ] <- input_list$data$map_ln_fishsel_devs[, , 12, , ]
  input_list$data$map_ln_fishsel_devs[, , 15, , ] <- input_list$data$map_ln_fishsel_devs[, , 12, , ]
  input_list$map$ln_fishsel_devs <- factor(input_list$data$map_ln_fishsel_devs)

  # Survey selectivity and catchability ----
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = c("rw_Fleet_1", "rw_Fleet_2", "rw_Fleet_3"),
    srvsel_pe_pars_spec = c("fix", "fix", "fix"),
    srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2"),
    srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
    srv_sel_model = c(
      "logist1_Fleet_1",
      "nonpar_Fleet_2",
      "nonpar_Fleet_3"
    ),
    srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
    srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2"),
    srv_q_spec = c("est_all", "est_all", "est_all"),
    Use_srv_q_prior = 1,
    srv_q_prior = data.frame(
      region = 1,
      block = 1,
      fleet = 1:2,
      mu = 1,
      sd = 2
    ),
    srv_sel_nonpar_est_bins =
      list(NULL,
           list(1:2, 3, 4, 5, 6, 7, 8:15),
           list(1:2, 3, 4, 5, 6, 7, 8:15)
      )
  )

  input_list$data$map_ln_srvsel_devs[, , 10, , 2:3] <- input_list$data$map_ln_srvsel_devs[, , 10, , 2:3]
  input_list$data$map_ln_srvsel_devs[, , 11, , 2:3] <- input_list$data$map_ln_srvsel_devs[, , 10, , 2:3]
  input_list$data$map_ln_srvsel_devs[, , 12, , 2:3] <- input_list$data$map_ln_srvsel_devs[, , 10, , 2:3]
  input_list$data$map_ln_srvsel_devs[, , 13, , 2:3] <- input_list$data$map_ln_srvsel_devs[, , 10, , 2:3]
  input_list$data$map_ln_srvsel_devs[, , 14, , 2:3] <- input_list$data$map_ln_srvsel_devs[, , 10, , 2:3]
  input_list$data$map_ln_srvsel_devs[, , 15, , 2:3] <- input_list$data$map_ln_srvsel_devs[, , 10, , 2:3]
  input_list$map$ln_srvsel_devs <- factor(input_list$data$map_ln_srvsel_devs)

  # Tagging (off) ----
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  # Weighting (with selectivity smoothing penalties) ----
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_SrvAgeComps  = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_SrvLenComps  = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    fish_sel_pen_wts = list(
      smooth_bin_curve   = 10,
      smooth_yr_diff     = 3,
      smooth_mean_center = 0
    ),
    srv_sel_pen_wts = list(
      smooth_bin_curve   = 10,
      smooth_yr_diff     = 5,
      smooth_mean_center = 0
    )
  )

  # Extract lists updated by helper functions ----
  data       <- input_list$data
  parameters <- input_list$par
  mapping    <- input_list$map

  # Fixed selectivity process-error sigmas (penalized likelihood) ----
  parameters$fishsel_pe_pars[, 4, , ]    <- log(0.5)   # fishery selex variance
  parameters$srvsel_pe_pars[, 1:2, , 1]  <- log(0.5)   # survey BTS - a50 and delta variance
  parameters$srvsel_pe_pars[, 4, , 2]    <- log(0.138) # survey ATS - ato variance

  # Recruitment bias ramp switched on at fit time ----
  data$bias_year       <- rep(61, 4)
  data$do_rec_bias_ramp <- 1

  # Fit model ----
  ebswp_rtmb_model <- fit_model(
    data,
    parameters,
    mapping,
    newton_loops = 3,
    silent = TRUE,
    do_optim = TRUE
  )

  ebswp_rtmb_model$sdrep <- RTMB::sdreport(ebswp_rtmb_model)

  # ---- Comparison baseline ----
  ssb_expected_vec <- c(366.8626, 457.713, 550.0043, 735.4656, 969.1294, 1275.5534,
                        1531.884, 1776.353, 1743.4722, 1518.9423, 1163.5947, 1120.1498,
                        1073.5644, 1274.3666, 1153.179, 1096.6816, 1344.2482, 2108.8223,
                        2952.8546, 3923.7941, 3975.1037, 4334.5781, 4063.6311, 3967.1995,
                        3699.1711, 3191.4743, 2702.5179, 2177.0367, 2093.0028, 2785.6286,
                        3222.8723, 3369.467, 3226.8465, 3082.7834, 2605.3678, 2740.7255,
                        2668.5781, 2736.9726, 2572.1176, 2639.5445, 2956.2391, 2610.4811,
                        2400.8533, 1975.0247, 1518.9446, 1626.6975, 1660.7823, 1975.0049,
                        2297.6076, 2536.8032, 2467.4582, 2337.0739, 2646.5515, 3141.7698,
                        3059.7701, 3073.4181, 2473.1161, 2648.1925, 3483.6239, 3342.8881,
                        3411.0757)

  rec_expected_vec <- c(5419.7147, 20768.3937, 15161.0254, 28360.7394, 24853.5719,
                        28797.2502, 24051.6128, 14263.111, 11640.7215, 28381.4595, 20071.2335,
                        16936.36, 13823.7898, 15249.9781, 28607.8904, 68272.0631, 30755.428,
                        32631.5845, 15568.8288, 45929.6206, 12710.0416, 30335.896, 14496.5441,
                        7475.6779, 6016.1821, 12624.0896, 50221.8444, 26297.9252, 21618.8472,
                        40402.4542, 14292.7797, 10744.9609, 22490.5451, 30119.145, 15104.2283,
                        16571.7474, 25077.2388, 34729.3722, 22908.4732, 14551.9982, 6643.6772,
                        4708.9962, 12305.9918, 25202.1208, 12794.1875, 43705.8681, 20260.3778,
                        14053.801, 12287.9564, 42715.7257, 49462.3543, 22325.8386, 9373.124,
                        11084.2854, 19177.1941, 75200.1184, 17212.529, 13544.3165, 14310.2958,
                        23409.9632, 23271.4757)

  expect_equal(ebswp_rtmb_model$rep$SSB[1, 1, ], ssb_expected_vec, tolerance = 1e-5)
  expect_equal(ebswp_rtmb_model$rep$Rec[1, 1, ], rec_expected_vec, tolerance = 1e-5)
  expect_true(ebswp_rtmb_model$sdrep$pdHess)
  expect_equal(ebswp_rtmb_model$rep$jnLL, 2595.243, tolerance = 1e-5)
  expect_jnLL_decomposes(ebswp_rtmb_model)

})
