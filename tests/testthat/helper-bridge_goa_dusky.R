# GOA dusky rockfish bridge: the 2024 Gulf of Alaska dusky rockfish assessment
# rebuilt in SPoRC.
#
# One Setup_Mod_* call per section, in the order vignette(
# "x_goa_dusky_rockfish_case_study") walks through them, with a reason for each
# argument that follows the assessment rather than a SPoRC default.
#
# The model is one area, one sex, one season, ages 4-33 with a plus group and
# lengths 21-52 cm, years 1977-2024.
#
#   Source                     Years        Observations  Likelihood
#   Catch                      1977-2024    48            Lognormal, weighted
#   Survey biomass             1990-2023    16            Lognormal
#   Fishery age comps          2000-2022    15            Multinomial
#   Fishery length comps       1991-2023    18            Multinomial
#   Survey age comps           1990-2023    16            Multinomial
#
# dev/make_sporc_obj_figs/make_goa_dusky_bridge_figs.R and
# tests/testthat/test-regression_dusky.R both build from this file, so the case
# study figures and the pinned regression cannot drift apart: a change to the
# specification moves both or neither.

build_goa_dusky_input <- function(dat) {

  ## Model dimensions ---------------------------------------------------------
  # one area, one sex, one season, so the region, sex and season subscripts in
  # vignette("c_model_equations") all collapse to one
  input_list <- Setup_Mod_Dim(
    n_pop = dat$n_pop,
    natal_region = dat$natal_region,
    years = dat$years,
    ages = dat$mod_ages,
    lens = dat$lens,
    n_regions = dat$n_regions,
    n_sexes = dat$n_sexes,
    n_fish_fleets = dat$n_fish_fleets,
    n_srv_fleets = dat$n_srv_fleets,
    n_seas = dat$n_seas,
    verbose = FALSE,
    store_config = TRUE
  )

  ## Recruitment --------------------------------------------------------------
  # a mean with annual deviations rather than a stock recruit function, sigmaR
  # fixed at exp(-0.1068576), roughly 0.899, and spawning at the start of the
  # year. the bias ramp is switched on but every ramp year is set to the
  # terminal year, which holds the ramp at zero over the whole series and so
  # leaves recruitment uncorrected. that combination is the assessment's own
  # convention written out, not a SPoRC default
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 1,
    bias_year = rep(length(dat$years), 4),
    sigmaR_switch = 1,
    ln_sigmaR = array(-0.1068576, dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    ln_global_R0 = log(2.7),
    t_spawn = dat$spwn_month
  )

  ## Biological dynamics ------------------------------------------------------
  # natural mortality is fixed at 0.07 for every age and year, so no prior is
  # needed. length compositions are fit through a size at age transition matrix
  # supplied as data, and age compositions pass through an ageing error matrix
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = dat$waa_arr,
    MatAA = dat$mataa_arr,
    fit_lengths = 1,
    SizeAgeTrans = dat$sizeage,
    AgeingError = dat$age_error_matrix,
    M_spec = "fix",
    Fixed_natmort = dat$fix_natmort,
    addtocomp = 0.00001
  )

  ## Movement and tagging -----------------------------------------------------
  # one area, so movement is the identity and nothing is tagged. both still have
  # to be declared
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  ## Catch and fishing mortality ----------------------------------------------
  # the assessment writes its catch and F statements as weighted sums of
  # squares. a weighted sum of squares and a normal likelihood with a fixed
  # standard deviation are the same statement up to a constant, related by
  # sigma = 1 / sqrt(2 w), so both sigmas are fixed at 1 / sqrt(2) to make the
  # quadratic term an unweighted sum of squares. the year specific catch weights
  # are then applied through Wt_Catch in the weighting section
  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list,
      ObsCatch = dat$ObsCatch,
      UseCatch = dat$UseCatch,
      Use_F_pen = 1,
      sigmaC_spec = "fix",
      ln_sigmaC = array(log(sqrt(1 / 2)), dim = c(input_list$data$n_regions,
                                                  length(input_list$data$years),
                                                  input_list$data$n_seas,
                                                  input_list$data$n_fish_fleets)),
      ln_sigmaF = array(log(sqrt(1 / 2)), dim = c(input_list$data$n_regions,
                                                  input_list$data$n_seas,
                                                  input_list$data$n_fish_fleets))
    )
  )

  ## Fishery compositions -----------------------------------------------------
  # no fishery index in this assessment, only compositions, so fish_idx_type is
  # "none". both age and length compositions are aggregated over the region and
  # fit multinomially
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = dat$ObsFishIdx,
    ObsFishIdx_SE = dat$ObsFishIdx_SE,
    UseFishIdx = dat$UseFishIdx,
    ObsFishAgeComps = dat$ObsFishAgeComps,
    UseFishAgeComps = dat$UseFishAgeComps,
    ISS_FishAgeComps = dat$ISS_FishAgeComps,
    ObsFishLenComps = dat$ObsFishLenComps,
    UseFishLenComps = dat$UseFishLenComps,
    ISS_FishLenComps = dat$ISS_FishLenComps,
    fish_idx_type = c("none"),
    FishAgeComps_LikeType = c("Multinomial"),
    FishLenComps_LikeType = c("Multinomial"),
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1")
  )

  ## Survey index and compositions --------------------------------------------
  # the index is lognormal, so its standard error has to be on the log scale.
  # the assessment reports an arithmetic standard error, and the conversion is
  # the coefficient of variation, which for a lognormal is the log scale
  # standard deviation to first order. that is the only transformation applied
  # to the survey data
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = dat$ObsSrvIdx,
    ObsSrvIdx_SE = dat$ObsSrvIdx_SE / dat$ObsSrvIdx,
    UseSrvIdx = dat$UseSrvIdx,
    ObsSrvAgeComps = dat$ObsSrvAgeComps,
    ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
    UseSrvAgeComps = dat$UseSrvAgeComps,
    ObsSrvLenComps = dat$ObsSrvLenComps,
    UseSrvLenComps = dat$UseSrvLenComps,
    ISS_SrvLenComps = dat$ISS_SrvLenComps,
    srv_idx_type = c("biom"),
    SrvAgeComps_LikeType = c("Multinomial"),
    SrvLenComps_LikeType = c("Multinomial"),
    SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    SrvLenComps_Type = c("agg_Year_1-terminal_Fleet_1")
  )

  ## Fishery selectivity and catchability -------------------------------------
  # the a50 and a95 parameterisation of the logistic, which is logist2.
  # selectivity is time invariant, so there are no deviations and no process
  # error, and fishery catchability is not used because there is no fishery
  # index to scale
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = c("none_Fleet_1"),
    fish_sel_blocks = c("none_Fleet_1"),
    fish_sel_model = c("logist2_Fleet_1"),
    fish_q_blocks = c("none_Fleet_1"),
    fish_fixed_sel_pars_spec = c("est_all"),
    fish_q_spec = c("fix")
  )

  ## Survey selectivity and catchability --------------------------------------
  # the same logistic form, with catchability estimated under a lognormal prior
  # centred on 1 with a standard deviation of 1 / sqrt(5)
  srv_q_prior <- data.frame(
    region = 1,
    block = 1,
    fleet = 1,
    mu = 1,
    sd = 0.447213595
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = c("none_Fleet_1"),
    srv_sel_blocks = c("none_Fleet_1"),
    srv_sel_model = c("logist2_Fleet_1"),
    srv_q_blocks = c("none_Fleet_1"),
    srv_fixed_sel_pars_spec = c("est_all"),
    srv_q_spec = c("est_all"),
    Use_srv_q_prior = 1,
    srv_q_prior = srv_q_prior,
    t_srv = array(0, dim = c(input_list$data$n_regions,
                             input_list$data$n_seas,
                             input_list$data$n_srv_fleets))
  )

  ## Weighting ----------------------------------------------------------------
  # the early catch series was reconstructed, so the assessment down-weights it
  # relative to the observer era: 2 through 1991 and 50 after. survey length
  # compositions are present in the data object but carry a weight of zero,
  # which is how the assessment treats them
  Wt_Catch <- array(0, dim = c(dat$n_regions, length(dat$years), dat$n_seas, dat$n_fish_fleets))
  Wt_Catch[, which(dat$years %in% 1977:1991), 1, ] <- 2
  Wt_Catch[, -which(dat$years %in% 1977:1991), 1, ] <- 50

  comp_dim <- c(input_list$data$n_regions,
                length(input_list$data$years),
                input_list$data$n_seas,
                input_list$data$n_sexes,
                input_list$data$n_fish_fleets)

  srv_comp_dim <- c(input_list$data$n_regions,
                    length(input_list$data$years),
                    input_list$data$n_seas,
                    input_list$data$n_sexes,
                    input_list$data$n_srv_fleets)

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = Wt_Catch,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1.66,
    Wt_Rec = 1,
    Wt_F = 2,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = comp_dim),
    Wt_FishLenComps = array(1, dim = comp_dim),
    Wt_SrvAgeComps = array(1, dim = srv_comp_dim),
    Wt_SrvLenComps = array(0, dim = srv_comp_dim)
  )

  return(input_list)

} # end build_goa_dusky_input
