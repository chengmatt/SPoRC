# BSAI blackspotted and rougheye rockfish bridge: the 2024 Bering Sea and
# Aleutian Islands assessment (ADMB) rebuilt in SPoRC.
#
# One Setup_Mod_* call per section, in the order vignette(
# "z_bsai_rougheye_rockfish_case_study") walks through them, with a reason for
# each argument that follows the assessment rather than a SPoRC default.
#
# The model is one area, one sex, one season, ages 3-54 with a plus group and
# ages reported over 3-45 with the last column pooling model ages 45-54, lengths
# 12-50 cm, years 1977-2024.
#
#   Source                     Years        Observations  Likelihood
#   Catch                      1977-2024    48            Lognormal, weighted
#   AI survey biomass          1991-2024    14            Lognormal
#   Fishery age comps          2004-2023    13            Multinomial
#   Fishery length comps       1979-2022    11            Multinomial
#   Survey age comps           1991-2022    13            Multinomial
#   Survey length comps        2024          1            Multinomial
#
# Both tests build from here. test-regression_rebs_bridge.R evaluates this
# configuration at the ADMB maximum likelihood estimate without optimizing, and
# test-regression_rebs_sgl.R refits from it, so a specification change cannot
# move one without the other.

# Build the input_list for the 2024 assessment configuration.
build_rebs_input <- function(dat) {

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

  ## Recruitment and the initial age structure --------------------------------
  # a mean with deviations rather than a stock recruit function, and the last
  # three years take the mean outright, which is dont_est_recdev_last = 3.
  #
  # the initial age structure is where this assessment is unusual. its
  # fyear_ac_option 3 gives the initial ages their own scalar, ln_rinit, rather
  # than reusing mean recruitment, with deviations that the ages beyond the
  # observed range share. init_age_devs_shared says which ages share which
  # deviation: ages 1-42 get their own and the remaining nine reuse the 42nd
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec",
    do_rec_bias_ramp = 1,
    bias_year = rep(n_yrs, 4),
    sigmaR_switch = 1,
    ln_sigmaR = array(log(dat$sigmaR), dim = c(2, dat$n_pop, dat$n_regions)),
    equil_init_age_strc = "stoch_shared_ages",
    init_age_devs_shared = c(1:42, rep(42, 9)),
    dont_est_recdev_last = 3,
    sigmaR_spec = "fix",
    init_age_strc = 1,
    t_spawn = (3 - 1) / 12,
    use_rinit = 1
  )

  ## Biological dynamics ------------------------------------------------------
  # natural mortality is estimated under a lognormal prior. the prior median is
  # shifted by exp(-cv^2 / 2) so that its MEAN is the assessment's 0.045, since
  # the assessment states the prior on the mean and SPoRC's is on the median.
  # length compositions are fit through a size at age transition matrix and age
  # compositions pass through an ageing error matrix
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = dat$WAA,
    MatAA = dat$MatAA,
    fit_lengths = 1,
    SizeAgeTrans = dat$SizeAgeTrans,
    AgeingError = dat$AgeingError,
    M_spec = "est_ln_M",
    Use_M_prior = 1,
    M_prior = data.frame(popblk = 1, regionblk = 1, yearblk = 1, ageblk = 1, sexblk = 1,
                         mu = dat$mean_M * exp(-dat$cv_M^2 / 2), sd = dat$cv_M),
    addtosrvidx = 1e-13,
    addtocomp = 1e-13
  )

  ## Movement and tagging -----------------------------------------------------
  # one area, so movement is the identity and nothing is tagged. both still have
  # to be declared
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  ## Catch and fishing mortality ----------------------------------------------
  # the assessment writes its catch and F penalties as sums of squares with
  # weights 50 and 0.1. a weighted sum of squares and a normal likelihood with a
  # fixed standard deviation are the same statement up to a constant, related by
  # sigma = 1 / sqrt(2 w), so the weights enter as these standard deviations
  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list,
      ObsCatch = dat$ObsCatch,
      UseCatch = dat$UseCatch,
      Use_F_pen = 1,
      sigmaC_spec = "fix",
      ln_sigmaC = array(log(sqrt(1 / (2 * 50))),
                        dim = c(dat$n_regions, n_yrs, dat$n_seas, dat$n_fish_fleets)),
      ln_sigmaF = array(log(sqrt(1 / (2 * 0.1))),
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
  # the Aleutian Islands survey supplies a biomass index with year specific
  # standard errors, age compositions, and a single year of length compositions
  # in 2024
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
  # logist1 is the a50 and slope parameterisation, which is the form this
  # assessment uses rather than the a50 and a95 form the rockfish bridges use.
  # selectivity is time invariant, and fishery catchability is not used because
  # there is no fishery index to scale
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = "none_Fleet_1",
    fish_sel_blocks = "none_Fleet_1",
    fish_sel_model = "logist1_Fleet_1",
    fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix"
  )

  ## Survey selectivity and catchability --------------------------------------
  # the same logistic form, with catchability estimated under a lognormal prior
  # whose median carries the same exp(-cv^2 / 2) shift as the M prior. the
  # survey is read at mid year
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = "none_Fleet_1",
    srv_sel_blocks = "none_Fleet_1",
    srv_sel_model = "logist1_Fleet_1",
    srv_q_blocks = "none_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all",
    Use_srv_q_prior = 1,
    srv_q_prior = data.frame(region = 1, block = 1, fleet = 1,
                             mu = dat$mean_q * exp(-dat$cv_q^2 / 2), sd = dat$cv_q),
    t_srv = array(0.5, dim = c(dat$n_regions, dat$n_seas, dat$n_srv_fleets))
  )

  ## Weighting ----------------------------------------------------------------
  # the catch and F weights already sit in their standard deviations above, so
  # only the composition weights are set here, and they are the assessment's
  # stage-2 multipliers
  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1,
    Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
    Wt_FishAgeComps = dat$Wt_FishAgeComps,
    Wt_FishLenComps = dat$Wt_FishLenComps,
    Wt_SrvAgeComps = dat$Wt_SrvAgeComps,
    Wt_SrvLenComps = dat$Wt_SrvLenComps
  )
} # end build_rebs_input


#' Set every parameter to the assessment's maximum likelihood estimate
#'
#' Evaluating at a known point separates a specification error from an
#' optimization difference. Most parameters are assigned directly; the initial
#' age structure needs a conversion, which is the block at the foot.
seed_rebs_mle <- function(input_list, dat) {

  mle <- dat$mle
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

  ## Recruitment, mortality, catchability and selectivity ---------------------
  input_list$par$ln_global_R0[] <- mle$mean_log_rec
  input_list$par$ln_rinit <- mle$log_rinit
  input_list$par$ln_M[] <- log(mle$M)
  input_list$par$ln_srv_q[] <- log(mle$q_srv)
  input_list$par$ln_F_mean[] <- mle$log_avg_fmort
  input_list$par$ln_F_devs[1, , 1, 1] <- mle$fmort_dev
  input_list$par$fish_fixed_sel_pars[] <- log(c(mle$sel_a50_fish, mle$sel_aslope_fish))
  input_list$par$srv_fixed_sel_pars[] <- log(c(mle$sel_a50_srv, mle$sel_aslope_srv))
  input_list$par$ln_RecDevs[1, 1, 1:length(mle$rec_dev)] <- mle$rec_dev

  ## Initial age structure ----------------------------------------------------
  # the assessment parameterises it as N(styr, j) = exp(log_rinit - M (j - 1) +
  # fydev_j), an absolute statement about numbers at age. SPoRC carries
  # multiplicative deviations from an equilibrium age structure, so the
  # deviations it wants are the log ratio of the two. build both and divide.
  # ages beyond the observed range reuse the last deviation, and the plus group
  # takes the geometric accumulation
  NAA_equil <- exp(mle$log_rinit) * exp(-(0:(n_ages - 1)) * mle$M)
  NAA_equil[n_ages] <- NAA_equil[n_ages - 1] * exp(-mle$M) / (1 - exp(-mle$M))
  NAA_styr <- NAA_equil
  for(j in 2:n_obs_ages) {
    NAA_styr[j] <- exp(mle$log_rinit - mle$M * (j - 1) + mle$fydev[j - 1])
  } # end j loop
  for(j in (n_obs_ages + 1):n_ages) {
    NAA_styr[j] <- exp(mle$log_rinit - mle$M * (j - 1) + mle$fydev[length(mle$fydev)])
  } # end j loop
  NAA_styr[n_ages] <- exp(mle$log_rinit - mle$M * (n_ages - 1) +
                            mle$fydev[length(mle$fydev)]) / (1 - exp(-mle$M))
  input_list$par$ln_InitDevs[1, 1, , ] <- (log(NAA_styr) - log(NAA_equil))[-1]

  input_list
} # end seed_rebs_mle
