# The 2023 BSAI northern rockfish assessment (ADMB) rebuilt in SPoRC. One area, one sex, one season,
# ages 3-45 with a plus group reported over 3-40, lengths 15-38 cm, years 1977-2023.
#
#   Source                     Years        Observations  Likelihood
#   Catch                      1977-2023    47            Lognormal, weighted
#   AI survey biomass          1991-2022    13            Lognormal
#   Fishery age comps          2000-2021    17            Multinomial
#   Fishery length comps       1977-2022    11            Multinomial
#   Survey age comps           1991-2022    13            Multinomial
#
# test-regression_bsai_nork_bridge.R evaluates this at the ADMB estimate without optimizing and
# test-regression_bsai_nork_sgl.R refits from it, so a specification change moves both or neither.

# Build the input_list for the 2023 assessment configuration.
build_bsai_nork_input <- function(dat) {

  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

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
  # a mean with deviations rather than a stock recruit function, and the last three
  # years take the mean outright.
  #
  # the initial age structure has its own scalar, ln_rinit, rather than reusing mean
  # recruitment, with the ages beyond the observed range sharing one deviation.
  #
  # the bias ramp stays flat, with the assessment's sigmaR^2 / 2 correction kept in
  # the seeds instead, which reproduces the series and its penalty at the seed point
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec",
    do_rec_bias_ramp = 1,
    bias_year = rep(n_yrs, 4),
    sigmaR_switch = 1,
    ln_sigmaR = array(log(dat$sigmaR), dim = c(2, dat$n_pop, dat$n_regions)),
    equil_init_age_strc = "stoch_shared_ages",
    init_age_devs_shared = c(1:(n_obs_ages - 1), rep(n_obs_ages - 1, n_ages - n_obs_ages)),
    dont_est_recdev_last = 3,
    sigmaR_spec = "fix",
    init_age_strc = 1,
    t_spawn = (dat$spawn_mo - 1) / 12,
    use_rinit = 1
  )

  ## Biological dynamics ------------------------------------------------------
  # natural mortality is estimated under a lognormal prior, its median shifted by
  # exp(-cv^2 / 2) so the mean is the assessment's 0.06, which is where it states the prior.
  #
  # weight at age is year varying, and the fishery has its own, which weighs the catch while
  # the population's goes to spawning biomass and the survey
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = dat$WAA,
    WAA_fish = dat$WAA_fish,
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
      mu = dat$mean_M * exp(-dat$cv_M^2 / 2),
      sd = dat$cv_M
    ),
    addtosrvidx = 1e-13,
    addtocomp = 1e-13
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
  # the assessment writes its catch and F penalties as weighted sums of squares at
  # weights 200 and 0.1, which is a normal at a fixed sigma = 1 / sqrt(2 w)
  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list,
      ObsCatch = dat$ObsCatch,
      UseCatch = dat$UseCatch,
      Use_F_pen = 1,
      sigmaC_spec = "fix",
      ln_sigmaC = array(log(sqrt(1 / (2 * dat$catch_wt))),
                        dim = c(dat$n_regions, n_yrs, dat$n_seas, dat$n_fish_fleets)),
      ln_sigmaF = array(log(sqrt(1 / (2 * dat$fmort_wt))),
                        dim = c(dat$n_regions, dat$n_seas, dat$n_fish_fleets))
    )
  )

  ## Fishery compositions -----------------------------------------------------
  # no fishery index, only age and length compositions, aggregated over the region
  # and fit multinomially
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
  # standard errors and age compositions
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
  # logist1 is the a50 and slope parameterization. selectivity is time invariant, and
  # there is no fishery index for a fishery catchability to scale
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
  # the penalty is on selectivity at age 30 itself, toward 1 at a standard deviation of
  # 0.003, not on the parameters, and the catchability prior is tight enough to hold q at 1
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = "none_Fleet_1",
    srv_sel_blocks = "none_Fleet_1",
    srv_sel_model = "logist1_Fleet_1",
    srv_q_blocks = "none_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all",
    Use_srv_q_prior = 1,
    srv_q_prior = data.frame(
      region = 1,
      block = 1,
      fleet = 1,
      mu = dat$mean_q * exp(-dat$cv_q^2 / 2),
      sd = dat$cv_q
    ),
    t_srv = array(0.5, dim = c(dat$n_regions, dat$n_seas, dat$n_srv_fleets)),
    Use_srv_selex_prior = 1,
    srv_selex_prior = data.frame(
      region = 1,
      fleet = 1,
      block = 1,
      sex = 1,
      par = which(dat$ages == 30),
      mu = 1.0,
      sd = 0.003,
      type = "value"
    )
  )

  ## Weighting ----------------------------------------------------------------
  # the catch and F weights already sit in their standard deviations above, so only the
  # composition weights are set here, the assessment's McAllister Ianelli multipliers
  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = dat$Wt_FishAgeComps,
    Wt_FishLenComps = dat$Wt_FishLenComps,
    Wt_SrvAgeComps = dat$Wt_SrvAgeComps,
    Wt_SrvLenComps = dat$Wt_SrvLenComps
  )
} # end build_bsai_nork_input


#' Set every parameter to the assessment's maximum likelihood estimate
#'
#' Two blocks here need converting rather than assigning: the recruitment bias correction
#' and the initial age structure.
seed_bsai_nork_mle <- function(input_list, dat) {

  mle <- dat$mle
  s2 <- dat$sigmaR^2 / 2
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

  ## Recruitment --------------------------------------------------------------
  # the assessment's three terminal recruits take exp(mean_log_rec + sigmaR^2 / 2) while the
  # estimated years are raw, and its deviations sum to zero.
  #
  # putting the bias correction into ln_global_R0 and shifting every seeded deviation down by
  # the same amount reproduces the series and the penalty
  input_list$par$ln_global_R0[] <- mle$mean_log_rec + s2
  input_list$par$ln_RecDevs[1, 1, ] <- mle$rec_dev - s2
  input_list$par$ln_rinit <- mle$log_rinit

  ## Mortality, catchability and selectivity ----------------------------------
  input_list$par$ln_M[] <- log(mle$M)
  input_list$par$ln_srv_q[] <- log(mle$q_srv)
  input_list$par$ln_F_mean[] <- mle$log_avg_fmort
  input_list$par$ln_F_devs[1, , 1, 1] <- mle$fmort_dev
  input_list$par$fish_fixed_sel_pars[] <- log(c(mle$sel_a50_fish, mle$sel_aslope_fish))
  input_list$par$srv_fixed_sel_pars[] <- log(c(mle$sel_a50_srv, mle$sel_aslope_srv))

  ## Initial age structure ----------------------------------------------------
  # the assessment writes it as N(styr, j) = exp(log_rinit - M (j - 1) + fydev_j), with ages
  # past the observed range sharing the last deviation and the plus group a geometric series.
  #
  # SPoRC wants multiplicative deviations from an equilibrium age structure, so build both
  # and take the log ratio
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
} # end seed_bsai_nork_mle


#' Hold the survey selectivity flat past age 30 for the seed evaluation
#'
#' The assessment evaluates its logistic over ages 3 to 30 only and holds both curves flat
#' beyond that, which `logist1` cannot express.
#'
#' The seed evaluation supplies the survey curve as a fixed input with that hold applied. A
#' refit estimates the uncapped logistic instead, a small and documented difference.
cap_bsai_nork_srv_sel <- function(data, dat) {
  n_yrs <- length(dat$years)
  srv_sel <- 1 / (1 + exp(-dat$mle$sel_aslope_srv * (dat$ages - dat$mle$sel_a50_srv)))
  srv_sel[dat$ages > 30] <- srv_sel[dat$ages == 30]
  srv_sel_input <- array(0, dim = c(dat$n_pop, dat$n_regions, n_yrs, dat$n_seas,
                                    length(dat$ages), dat$n_sexes, dat$n_srv_fleets))
  for(y in seq_len(n_yrs)) srv_sel_input[1, 1, y, 1, , 1, 1] <- srv_sel
  data$use_fixed_srv_sel <- 1
  data$srv_sel_input <- srv_sel_input
  data
} # end cap_bsai_nork_srv_sel
