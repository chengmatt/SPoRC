# Shared setup for the BSAI northern rockfish bridge tests. The bridge test
# evaluates this configuration at the ADMB maximum likelihood estimate without
# optimizing; the pinned test refits from it. Keeping the two on one builder
# means a specification change cannot move one without the other.

# Build the input_list for the 2023 assessment configuration.
build_bsai_nork_input <- function(dat) {

  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

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

  # Recruitment is a mean with deviations rather than a stock recruit function,
  # and the last three years take the mean outright. The assessment's
  # fyear_ac_option 3 gives the initial age structure its own scalar, ln_rinit,
  # with deviations that ages beyond the observed range share. The bias ramp
  # stays off: the assessment's sigmaR^2 / 2 correction is carried in the seeds
  # instead (see seed_bsai_nork_mle), which reproduces both the recruitment and
  # its penalty at the seed point.
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

  # Natural mortality is estimated under a lognormal prior. The prior median is
  # shifted by exp(-cv^2 / 2) so its mean is the assessment's 0.06. Weight at
  # age is year varying, and the fishery carries its own block, which prices
  # the catch while the population block prices spawning biomass and the
  # survey.
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
    M_prior = data.frame(popblk = 1, regionblk = 1, yearblk = 1, ageblk = 1, sexblk = 1,
                         mu = dat$mean_M * exp(-dat$cv_M^2 / 2), sd = dat$cv_M),
    addtosrvidx = 1e-13,
    addtocomp = 1e-13
  )

  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  # The assessment writes its catch and F penalties as sums of squares with
  # weights 200 and 0.1, which are the same statements as normal likelihoods
  # with these standard deviations.
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

  # There is no fishery index in this assessment, only compositions.
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

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = "none_Fleet_1",
    fish_sel_blocks = "none_Fleet_1",
    fish_sel_model = "logist1_Fleet_1",
    fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix"
  )

  # The assessment's survey selectivity constraint is not a prior on the
  # parameters: it penalises the realized selectivity at age 30 towards 1 with
  # a standard deviation of 0.003. The type value row makes that exact
  # statement, and it is load bearing: without it the survey age compositions
  # do not identify the selectivity asymptote. The catchability prior is tight
  # enough to pin q at 1.
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
    t_srv = array(0.5, dim = c(dat$n_regions, dat$n_seas, dat$n_srv_fleets)),
    Use_srv_selex_prior = 1,
    srv_selex_prior = data.frame(region = 1, fleet = 1, block = 1, sex = 1,
                                 par = which(dat$ages == 30),
                                 mu = 1.0, sd = 0.003, type = "value")
  )

  # The composition weights are the assessment's McAllister Ianelli multipliers.
  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1,
    Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
    Wt_FishAgeComps = dat$Wt_FishAgeComps,
    Wt_FishLenComps = dat$Wt_FishLenComps,
    Wt_SrvAgeComps = dat$Wt_SrvAgeComps,
    Wt_SrvLenComps = dat$Wt_SrvLenComps
  )
}

# Set every parameter to the assessment's maximum likelihood estimate.
seed_bsai_nork_mle <- function(input_list, dat) {

  mle <- dat$mle
  s2 <- dat$sigmaR^2 / 2
  n_ages <- length(dat$ages)
  n_obs_ages <- length(dat$obs_ages)

  # The assessment builds its three dev-free terminal recruits as
  # exp(mean_log_rec + sigmaR^2 / 2) while the estimated years are raw, and its
  # rec_dev is a dev_vector summing to zero. Carrying the bias correction in
  # ln_global_R0 and shifting every seeded deviation down by the same amount
  # reproduces both the recruitment series and the penalty value exactly.
  input_list$par$ln_global_R0[] <- mle$mean_log_rec + s2
  input_list$par$ln_RecDevs[1, 1, ] <- mle$rec_dev - s2
  input_list$par$ln_rinit <- mle$log_rinit
  input_list$par$ln_M[] <- log(mle$M)
  input_list$par$ln_srv_q[] <- log(mle$q_srv)
  input_list$par$ln_F_mean[] <- mle$log_avg_fmort
  input_list$par$ln_F_devs[1, , 1, 1] <- mle$fmort_dev
  input_list$par$fish_fixed_sel_pars[] <- log(c(mle$sel_a50_fish, mle$sel_aslope_fish))
  input_list$par$srv_fixed_sel_pars[] <- log(c(mle$sel_a50_srv, mle$sel_aslope_srv))

  # The assessment parameterises the initial age structure as
  # N(styr, j) = exp(log_rinit - M (j - 1) + fydev_j), with ages beyond the
  # observed range sharing the last deviation and the plus group solved as a
  # geometric series. SPoRC carries multiplicative deviations from an
  # equilibrium age structure, so the deviations are the log ratio of the two.
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
}

# The assessment evaluates its logistic selectivity over ages 3 to 30 only
# (nselages) and holds both curves at the age 30 value beyond that, which
# SPoRC's logist1 form cannot express. For the seed evaluation the survey
# curve is supplied as a fixed input with that edge hold applied, so the
# comparison isolates the likelihoods from the selectivity form. A refit
# estimates the uncapped logistic instead, which is a documented and
# negligible difference.
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
}
