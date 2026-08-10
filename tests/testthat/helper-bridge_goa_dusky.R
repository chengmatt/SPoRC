# Shared specification for the 2024 GOA dusky rockfish bridge. The figures in
# dev/make_sporc_obj_figs/make_goa_dusky_bridge_figs.R source this file so the
# case study figures and the pinned regression test cannot drift apart: a change
# to the specification moves both or neither. tests/testthat/test-regression_dusky.R
# builds from here as well.

build_goa_dusky_input <- function(dat) {

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

  # The bias ramp is switched on but its ramp years are set to the terminal year,
  # which leaves recruitment uncorrected, matching the assessment.
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

  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

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

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = dat$ObsSrvIdx,
    ObsSrvIdx_SE = dat$ObsSrvIdx_SE / dat$ObsSrvIdx, # the assessment reports an arithmetic SE
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

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = c("none_Fleet_1"),
    fish_sel_blocks = c("none_Fleet_1"),
    fish_sel_model = c("logist2_Fleet_1"),
    fish_q_blocks = c("none_Fleet_1"),
    fish_fixed_sel_pars_spec = c("est_all"),
    fish_q_spec = c("fix")
  )

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

  # The assessment down-weights the early catch series, which was reconstructed,
  # relative to the observer era.
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
