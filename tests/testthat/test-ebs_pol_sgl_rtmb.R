library(SPoRC)
library(testthat)
data("sgl_rg_ebswp_data")

test_that("Single-region EBS Pollock RTMB model produces expected results", {

  ## Initialize model dimensions and data list----
  input_list <- Setup_Mod_Dim(
    years = sgl_rg_ebswp_data$years,
    # vector of years
    ages = sgl_rg_ebswp_data$ages,
    # vector of ages
    lens = NA,
    # number of lengths
    n_regions = 1,
    # number of regions
    n_sexes = 1,
    # number of sexes
    n_fish_fleets = 1,
    # number of fishery fleets
    n_srv_fleets = 3, # number of survey fleets
    # number of seasons
    n_seas = sgl_rg_ebswp_data$n_seas,
    # Populaiton stuff
    n_pop = sgl_rg_ebswp_data$n_pop,
    natal_region = sgl_rg_ebswp_data$natal_region,
    verbose = TRUE
  )

  inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(
    input_list = input_list,

    # Model options
    do_rec_bias_ramp = 0,
    # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
    sigmaR_switch = 1,
    # when to switch from early to late sigmaR (switch in first year)
    ln_sigmaR = array(log(1), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    # Starting values for early and late sigmaR
    rec_model = "bh_rec",
    # recruitment model
    steepness_h = array(inv_steepness(0.623013), dim = c(input_list$data$n_pop, input_list$data$n_regions)),
    h_spec = "fix",
    # fixing steepness
    sigmaR_spec = "fix",
    # fix early sigmaR and late sigmaR
    init_age_strc = 1,
    ln_global_R0 = 10,
    t_spawn = 0.25,
    equil_init_age_strc = 2
    # starting value for r0
  )

  # Setup a fixed natural mortality array for use
  fix_natmort <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), 1))
  fix_natmort[,,,1,] <- 0.9 # age 1 M
  fix_natmort[,,,2,] <- 0.45 # age 2 M
  fix_natmort[,,,-c(1,2),] <- 0.3 # age 3+ M

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,

    # Data inputs
    WAA = sgl_rg_ebswp_data$WAA,
    MatAA = sgl_rg_ebswp_data$MatAA,

    # Model options
    # mean and sd for M prior
    fit_lengths = 0,
    # don't fit length compositions
    M_spec = "fix",
    # fixing natural mortality
    Fixed_natmort = fix_natmort
  )

  # Setup movement stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,

    # Data inputs
    ObsCatch = sgl_rg_ebswp_data$ObsCatch,
    UseCatch = sgl_rg_ebswp_data$UseCatch,

    # Model options
    Use_F_pen = 1,
    # whether to use f penalty, == 0 don't use, == 1 use
    sigmaC_spec = "fix",
    # fixing catch standard deviation
    ln_sigmaC = array(log(0.05), dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1))
    # starting / fixed value for catch standard deviation
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    # data inputs
    ObsFishIdx = sgl_rg_ebswp_data$ObsFishIdx,
    ObsFishIdx_SE = sgl_rg_ebswp_data$ObsFishIdx_SE,
    UseFishIdx = sgl_rg_ebswp_data$UseFishIdx,
    ObsFishAgeComps = sgl_rg_ebswp_data$ObsFishAgeComps,
    UseFishAgeComps = sgl_rg_ebswp_data$UseFishAgeComps,
    ISS_FishAgeComps = sgl_rg_ebswp_data$ISS_FishAgeComps,
    ObsFishLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, 1)),
    UseFishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1)),
    ISS_FishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, 1)),

    # Model options
    fish_idx_type = c("biom"),
    # indices for fishery
    FishAgeComps_LikeType = c("Multinomial"),
    # age comp likelihoods for fishery fleet
    FishLenComps_LikeType = c("none"),
    # length comp likelihoods for fishery
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    # age comp structure for fishery
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
    # length comp structure for fishery
  )

  # Setup survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,

    # data inputs
    ObsSrvIdx = sgl_rg_ebswp_data$ObsSrvIdx,
    ObsSrvIdx_SE = sgl_rg_ebswp_data$ObsSrvIdx_SE,
    UseSrvIdx = sgl_rg_ebswp_data$UseSrvIdx,
    ObsSrvAgeComps = sgl_rg_ebswp_data$ObsSrvAgeComps,
    ISS_SrvAgeComps = sgl_rg_ebswp_data$ISS_SrvAgeComps,
    UseSrvAgeComps = sgl_rg_ebswp_data$UseSrvAgeComps,
    ObsSrvLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, 3)),
    UseSrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, 3)),
    ISS_SrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, 3)),

    # Model options
    srv_idx_type = c("biom", "biom", "biom"),
    # abundance and biomass for survey fleet 1, 2, and 3
    SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "Multinomial"),
    # survey age composition likelihood for survey fleet 1, 2, and 3
    SrvLenComps_LikeType = c("none", "none", "none"),
    #  survey length composition likelihood for survey fleet 1, 2, and 3
    SrvAgeComps_Type = c(
      "agg_Year_1-terminal_Fleet_1",
      "agg_Year_1-terminal_Fleet_2",
      "none_Year_1-terminal_Fleet_3"
    ),
    # survey age comp type

    SrvLenComps_Type = c(
      "none_Year_1-terminal_Fleet_1",
      "none_Year_1-terminal_Fleet_2",
      "none_Year_1-terminal_Fleet_3"
    )
    # survey length comp type
  )


  # Setup fishery selectivity and catchability
  input_list <- Setup_Mod_Fishsel_and_Q(

    input_list = input_list,

    # Model options
    # fishery selectivity, whether continuous time-varying
    cont_tv_fish_sel = c("2dar1_Fleet_1"),
    fishsel_pe_pars_spec = "fix", # doing penalized likelihood for selex devs
    fish_sel_devs_spec = "est_all", # estimating all sel devs
    corr_opt_semipar = "corr_zero_y_b", # making sure 2d correaltions are 0, collapses to a simple iid case
    # fishery selectivity blocks
    fish_sel_blocks = c("none_Fleet_1"),
    # fishery selectivity form
    fish_sel_model = c("logist1_Fleet_1"),
    # fishery catchability blocks
    fish_q_blocks = c("none_Fleet_1"),
    # whether to estiamte all fixed effects for fishery selectivity
    fish_fixed_sel_pars_spec = c("est_all"),
    # whether to estiamte all fixed effects for fishery catchability
    fish_q_spec = c("est_all")
  )

  # Setup survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,

    # Model options
    # survey selectivity, whether continuous time-varying
    cont_tv_srv_sel = c("iid_Fleet_1", "2dar1_Fleet_2", "2dar1_Fleet_3"),
    srvsel_pe_pars_spec = c("fix", "fix", "fix"), # penalize survey selex devs
    srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2"), # estimating all srv selex devs
    corr_opt_semipar = c(NA, "corr_zero_y_b", "corr_zero_y_b"), # setting corelations at 0, so 2dar1 collapses to simple iid semi-parametric devs

    # survey selectivity blocks
    srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
    # survey selectivity form
    srv_sel_model = c(
      "logist1_Fleet_1",
      "logist1_Fleet_2",
      "logist1_Fleet_3"
    ),
    # survey catchability blocks
    srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
    # whether to estiamte all fixed effects for survey selectivity
    srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2"),
    # whether to estiamte all fixed effects for survey catchability
    srv_q_spec = c("est_all", "est_all", "est_all")
  )

  # Setup tagging stuff
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years),
                                       input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_srv_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years),
                                       input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_srv_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets))
  )


  # extract out lists updated with helper functions
  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # selex sigma to fix at, given penalized likelihood
  parameters$fishsel_pe_pars[,4,,] <- log(0.075) # fishery selex variance
  parameters$srvsel_pe_pars[,1:2,,1] <- log(0.075) # survey BTS - a50 and delta variance
  parameters$srvsel_pe_pars[,4,,2] <- log(0.15) # survey ATS and ato variance


  # Fit model
  ebswp_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                # random = NULL,
                                newton_loops = 3,
                                silent = TRUE
  )

  ebswp_rtmb_model$sdrep <- RTMB::sdreport(ebswp_rtmb_model)

  ssb_expected_vec <- c(604.8303, 620.3900, 650.9218, 749.7160, 864.9073, 1001.7853, 1079.4247, 1164.7102, 1120.7780, 987.6462,
                        720.1733, 667.4834, 635.9955, 747.5421, 664.7550, 658.3988, 914.7824, 1591.3212, 2307.3483, 3140.5967,
                        3256.9845, 3711.8637, 3619.9472, 3706.5877, 3581.7596, 3136.6307, 2664.5500, 2166.4447, 2150.8002, 2928.0622,
                        3416.5021, 3590.8871, 3444.9857, 3261.1167, 2758.6261, 2956.9983, 2957.1241, 3082.2228, 2882.4844, 2903.2505,
                        3196.0122, 2801.1688, 2554.6853, 2136.4215, 1699.6423, 1896.1574, 1951.9939, 2284.2466, 2612.3691, 2860.9857,
                        2775.7142, 2659.4668, 2963.2610, 3354.2383, 3063.3879, 2824.2168, 2038.8854, 2196.5457, 3155.6046, 3226.9139,
                        3437.5710)

  rec_expected_vec <- c(4336.183, 17914.462, 10405.798, 20871.329, 20385.116, 25215.438, 23145.713, 14095.170, 11009.571, 20797.617,
                        11958.594, 11396.141, 11543.669, 12878.865, 23161.020, 55626.388, 26242.962, 30956.401, 15945.294, 48113.748,
                        14181.902, 35934.771, 15256.019,  7509.865,  5909.875, 12336.732, 56068.126, 29073.899, 22745.314, 43862.523,
                        15026.768, 11208.157, 25745.114, 37052.757, 17211.965, 18061.800, 26783.645, 37348.415, 24088.436, 15601.540,
                        7200.247,  5118.634, 13881.948, 29676.666, 12726.764, 47689.853, 21942.654, 14328.258, 12741.027, 49722.905,
                        49113.069, 18013.886,  5952.349,  6010.924, 12352.300, 83426.122, 23304.418, 16428.846, 14114.836, 27634.359,
                        50430.467)

expect_equal(ebswp_rtmb_model$rep$SSB[1,1,], ssb_expected_vec, tolerance = 1e-2)
expect_equal(ebswp_rtmb_model$rep$Rec[1,1,], rec_expected_vec, tolerance = 1e-2)
expect_true(ebswp_rtmb_model$sdrep$pdHess)


})

