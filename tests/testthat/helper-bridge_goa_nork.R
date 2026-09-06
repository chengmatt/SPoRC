# The 2024 GOA northern rockfish assessment (ADMB) rebuilt in SPoRC. One area, one sex, one season,
# ages 2-51 with a plus group reported over 2-45, lengths 15-45 cm, years 1961-2024.
#
#   Source                     Years        Observations  Likelihood
#   Catch                      1961-2024    64            Lognormal, weighted
#   Survey biomass             1990-2023    16            Lognormal
#   Fishery age comps          1998-2022    16            Multinomial
#   Fishery length comps       1991-2023    17            Multinomial
#   Survey age comps           1990-2023    16            Multinomial
#
# test-regression_goa_nork_bridge.R evaluates this at the ADMB estimate without optimizing and
# test-regression_goa_nork_sgl.R refits from it, so a specification change moves both or neither.

# Build the input_list for the 2024 assessment configuration.
build_goa_nork_input <- function(dat) {

  yrs <- dat$years
  n_yrs <- length(yrs)

  ## Model dimensions ---------------------------------------------------------
  # one area, one sex, one season, so the region, sex and season subscripts in
  # vignette("c_model_equations") all collapse to one
  input_list <- Setup_Mod_Dim(
    n_pop = dat$n_pop,
    years = yrs,
    ages = dat$ages,
    lens = dat$lens,
    n_regions = dat$n_regions,
    n_sexes = dat$n_sexes,
    n_fish_fleets = dat$n_fish_fleets,
    n_srv_fleets = dat$n_srv_fleets,
    n_seas = dat$n_seas,
    verbose = FALSE
  )

  ## Recruitment --------------------------------------------------------------
  # a mean with deviations in every year and no lognormal bias correction, which
  # is what the ADMB template does. do_rec_bias_ramp = 0 says so directly
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec",
    do_rec_bias_ramp = 0,
    bias_year = rep(n_yrs, 4),
    sigmaR_switch = 1,
    ln_sigmaR = array(log(dat$sigmaR), dim = c(2, dat$n_pop, dat$n_regions)),
    dont_est_recdev_last = 0,
    sigmaR_spec = "fix",
    init_age_strc = 1,
    t_spawn = 0
  )

  ## Biological dynamics ------------------------------------------------------
  # natural mortality is estimated under a lognormal prior centered on the
  # assessment's mean with its coefficient of variation. length compositions are
  # fit through a size at age transition matrix and age compositions pass
  # through an ageing error matrix, so both are declared here rather than at the
  # composition call
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = dat$WAA,
    MatAA = dat$MatAA,
    fit_lengths = 1,
    SizeAgeTrans = dat$SizeAgeTrans,
    AgeingError = dat$AgeingError,
    M_spec = "est_ln_M",
    Use_M_prior = 1,
    M_prior = data.frame(
      popblk = 1,
      regionblk = 1,
      yearblk = 1,
      ageblk = 1,
      sexblk = 1,
      mu = dat$mean_M,
      sd = dat$cv_M
    ),
    addtosrvidx = 0.00001,
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
  # sigma = 1 / sqrt(2 w), so the weights enter through ln_sigmaC rather than a
  # separate multiplier. the reconstructed early catches have a weight of 5 and
  # the observer era series 50. the F deviations take sigma = 1 / sqrt(2), with
  # their overall weight applied in the weighting section
  ln_sigmaC <- array(NA_real_, dim = c(dat$n_regions, n_yrs, dat$n_seas, dat$n_fish_fleets))
  ln_sigmaC[1, , 1, 1] <- log(sqrt(1 / (2 * dat$catch_wt)))
  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list,
      ObsCatch = dat$ObsCatch,
      UseCatch = dat$UseCatch,
      Use_F_pen = 1,
      sigmaC_spec = "fix",
      ln_sigmaC = ln_sigmaC,
      ln_sigmaF = array(log(sqrt(1 / 2)),
                        dim = c(dat$n_regions, dat$n_seas, dat$n_fish_fleets))
    )
  )

  ## Fishery compositions -----------------------------------------------------
  # no fishery index in this assessment, only compositions, so fish_idx_type is
  # "none" and the index arrays are declared empty. both age and length
  # compositions are aggregated over the region and fit multinomially
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = array(NA, dim = c(dat$n_regions, n_yrs, dat$n_seas, dat$n_fish_fleets)),
    ObsFishIdx_SE = array(NA, dim = c(dat$n_regions, n_yrs, dat$n_seas, dat$n_fish_fleets)),
    UseFishIdx = array(0, dim = c(dat$n_regions, n_yrs, dat$n_seas, dat$n_fish_fleets)),
    ObsFishAgeComps = dat$ObsFishAgeComps,
    UseFishAgeComps = dat$UseFishAgeComps,
    ISS_FishAgeComps = dat$ISS_FishAgeComps,
    ObsFishLenComps = dat$ObsFishLenComps,
    UseFishLenComps = dat$UseFishLenComps,
    ISS_FishLenComps = dat$ISS_FishLenComps,
    fish_idx_type = "none",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "Multinomial",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "agg_Year_1-terminal_Fleet_1"
  )

  ## Survey index and compositions --------------------------------------------
  # the bottom trawl survey supplies a biomass index and age compositions. the
  # index is lognormal with year specific standard errors, and the survey is fit
  # at the start of the year, which t_srv = 0 sets in the selectivity section
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = dat$ObsSrvIdx,
    ObsSrvIdx_SE = dat$ObsSrvIdx_SE,
    UseSrvIdx = dat$UseSrvIdx,
    ObsSrvAgeComps = dat$ObsSrvAgeComps,
    ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
    UseSrvAgeComps = dat$UseSrvAgeComps,
    ObsSrvLenComps = dat$ObsSrvLenComps,
    UseSrvLenComps = dat$UseSrvLenComps,
    ISS_SrvLenComps = dat$ISS_SrvLenComps,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "Multinomial",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "agg_Year_1-terminal_Fleet_1"
  )

  ## Fishery selectivity and catchability -------------------------------------
  # the a50 and a95 parameterization of the logistic, which is logist2.
  # selectivity is time invariant, so there are no deviations and no process
  # error, and fishery catchability is not used because there is no fishery
  # index to scale
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = "none_Fleet_1",
    fish_sel_blocks = "none_Fleet_1",
    fish_sel_model = "logist2_Fleet_1",
    fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix"
  )

  ## Survey selectivity and catchability --------------------------------------
  # the same logistic form, with catchability estimated under a lognormal prior
  # centered on the assessment's mean, loose enough to let the data move it.
  # t_srv = 0 puts the survey at the start of the year
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = "none_Fleet_1",
    srv_sel_blocks = "none_Fleet_1",
    srv_sel_model = "logist2_Fleet_1",
    srv_q_blocks = "none_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all",
    Use_srv_q_prior = 1,
    srv_q_prior = data.frame(
      region = 1,
      block = 1,
      fleet = 1,
      mu = dat$mean_q,
      sd = dat$cv_q
    ),
    t_srv = array(0, dim = c(dat$n_regions, dat$n_seas, dat$n_srv_fleets))
  )

  ## Weighting ----------------------------------------------------------------
  # the survey index and the F penalty hold the assessment's own fixed weights,
  # and every composition source holds the assessment's multipliers, which
  # ship in the data object
  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = dat$srv_wt,
    Wt_Rec = 1,
    Wt_F = dat$fmort_wt,
    Wt_Tagging = 0,
    Wt_FishAgeComps = dat$Wt_FishAgeComps,
    Wt_FishLenComps = dat$Wt_FishLenComps,
    Wt_SrvAgeComps = dat$Wt_SrvAgeComps,
    Wt_SrvLenComps = dat$Wt_SrvLenComps
  )
} # end build_goa_nork_input


#' Set every parameter to the assessment's maximum likelihood estimate
#'
#' Evaluating at a known point separates a specification error from an
#' optimization difference: if the population and the likelihood agree at the
#' ADMB solution, the two models are the same model. Every parameter can be
#' assigned directly, because recruitment is a mean with deviations rather than
#' a stock recruit function and nothing has to be solved by substitution.
#'
#' The initial age structure deviations were back derived from the ADMB numbers
#' at age when the data object was built, so seeding them reproduces the ADMB
#' starting conditions exactly.
seed_goa_nork_mle <- function(input_list, dat) {

  mle <- dat$mle

  ## Recruitment and the initial age structure --------------------------------
  input_list$par$ln_global_R0[] <- mle$log_mean_R
  input_list$par$ln_RecDevs[1, 1, ] <- mle$log_Rt
  input_list$par$ln_InitDevs[1, 1, , ] <- mle$init_devs

  ## Natural mortality --------------------------------------------------------
  input_list$par$ln_M[] <- log(mle$M)

  ## Fishing mortality --------------------------------------------------------
  input_list$par$ln_F_mean[] <- mle$log_mean_F
  input_list$par$ln_F_devs[1, , 1, 1] <- mle$log_Ft

  ## Selectivity and catchability ---------------------------------------------
  # SPoRC estimates the logistic parameters on the log scale
  input_list$par$fish_fixed_sel_pars[] <- log(c(mle$a50C, mle$deltaC))
  input_list$par$srv_fixed_sel_pars[] <- log(c(mle$a50S, mle$deltaS))
  input_list$par$ln_srv_q[] <- log(mle$q)

  input_list
} # end seed_goa_nork_mle
