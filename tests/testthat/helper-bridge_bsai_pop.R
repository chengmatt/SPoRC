# The 2024 BSAI Pacific ocean perch assessment (ADMB) rebuilt in SPoRC. One area, one sex, one season,
# ages 3-46 with a plus group reported over 3-40, lengths 15-39 cm, years 1960-2024.
#
# The most structurally involved of the five rockfish bridges: two survey fleets, a first year
# equilibrium under fixed historical F, and a bicubic spline fishery selectivity over year and age.
#
#   Source                        Years        Observations  Likelihood
#   Catch                         1960-2024    65            Lognormal, weighted
#   AI survey biomass             1991-2024    14            Lognormal
#   EBS slope survey biomass      2002-2016     6            Lognormal
#   Fishery age comps             1981-2023    22            Multinomial
#   Fishery length comps          1964-2022    32            Multinomial
#   AI survey age comps           1991-2022    13            Multinomial
#   EBS survey age comps          2002-2016     6            Multinomial
#   AI survey length comps        2024          1            Multinomial
#
# test-regression_bsai_pop_bridge.R evaluates this at the ADMB estimate without optimizing and
# test-regression_bsai_pop_sgl.R refits from it, so a specification change moves both or neither.

# Build the input_list for the 2024 assessment configuration.
build_bsai_pop_input <- function(dat) {

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
  # few years take the mean outright.
  #
  # the first year sits in equilibrium under a FIXED historical F of 0.01, which
  # the assessment holds independent of the estimated mean F. init_F_form "abs"
  # is what keeps the two separate, so the initial age structure cannot be
  # depleted by raising the mean F.
  #
  # bias_year is indexed in deviation space rather than calendar years, so this
  # range puts the whole series in the fully bias corrected limb, which centers
  # the penalty on -sigmaR^2 / 2 to match the shifted deviations the seeds have
  init_F_par <- array(log(1e-10), dim = c(dat$n_regions, dat$n_seas, dat$n_fish_fleets))
  init_F_par[1, 1, 1] <- log(dat$mle$historic_F)

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec",
    do_rec_bias_ramp = 1,
    bias_year = c(1, 1, n_yrs + 1, n_yrs + 1),
    sigmaR_switch = 1,
    ln_sigmaR = array(log(c(dat$sigmaR, dat$sigmaR)), dim = c(2, dat$n_pop, dat$n_regions)),
    dont_est_recdev_last = dat$fixedrec,
    equil_init_age_strc = "equil",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    init_F_form = "abs",
    init_F_spec = "fix",
    init_F_par = init_F_par,
    t_spawn = (dat$spawn_mo - 1) / 12,
    use_rinit = 1
  )

  ## Biological dynamics ------------------------------------------------------
  # natural mortality is estimated under a lognormal prior. the prior median is
  # shifted by exp(-cv^2 / 2) so that its MEAN is the assessment's 0.05, since
  # the assessment states the prior on the mean and SPoRC's is on the median.
  #
  # maturity is estimated inside the assessment template and absent from its
  # report object, so the fitted logistic is rebuilt when the data object is
  # made and kept fixed here
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
  # the assessment writes its catch and F penalties as sums of squares with
  # weights 500 and 0.1. a weighted sum of squares and a normal likelihood with
  # a fixed standard deviation are the same statement up to a constant, related
  # by sigma = 1 / sqrt(2 w), so the weights enter as these standard deviations
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
  # no fishery index in this assessment, only compositions, so fish_idx_type is
  # "none" and the index arrays are declared empty
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
  # two survey fleets, the Aleutian Islands bottom trawl survey and the eastern
  # Bering Sea slope survey, both fit to biomass. every argument is now a vector
  # over fleets
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
    srv_idx_type = rep("biom", dat$n_srv_fleets),
    SrvAgeComps_LikeType = rep("Multinomial", dat$n_srv_fleets),
    SrvLenComps_LikeType = rep("Multinomial", dat$n_srv_fleets),
    SrvAgeComps_Type = paste0("agg_Year_1-terminal_Fleet_", 1:dat$n_srv_fleets),
    SrvLenComps_Type = paste0("agg_Year_1-terminal_Fleet_", 1:dat$n_srv_fleets)
  )

  ## Fishery selectivity and catchability -------------------------------------
  # a bicubic spline over a 5 year by 5 age node grid, exponentiated with no
  # normalization. SelStyr and NSelBins are not cosmetic: they set the year
  # range and bin count the smoothness penalties normalize by, and they impose
  # the assessment's edge holds, flat before 1964 and flat beyond age 40
  fish_sel_spec <- paste0("bicubic_Bin_", dat$fsh_age_nodes,
                          "_Yr_", dat$fsh_yr_nodes,
                          "_Fleet_1",
                          "_SelStyr_", dat$fsh_sel_styr,
                          "_NSelBins_", dat$nselages)

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = "none_Fleet_1",
    fish_sel_blocks = "none_Fleet_1",
    fish_sel_model = fish_sel_spec,
    fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix",
    use_fixed_fish_sel = 0
  )

  ## Survey selectivity and catchability --------------------------------------
  # logistic for both fleets. only the Aleutian Islands survey has a
  # catchability prior, and supplying a single row for fleet 1 is what restricts
  # it to that fleet
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = paste0("none_Fleet_", 1:dat$n_srv_fleets),
    srv_sel_blocks = paste0("none_Fleet_", 1:dat$n_srv_fleets),
    srv_sel_model = paste0("logist1_Fleet_", 1:dat$n_srv_fleets),
    srv_q_blocks = paste0("none_Fleet_", 1:dat$n_srv_fleets),
    srv_fixed_sel_pars_spec = rep("est_all", dat$n_srv_fleets),
    srv_q_spec = rep("est_all", dat$n_srv_fleets),
    Use_srv_q_prior = 1,
    srv_q_prior = data.frame(
      region = 1,
      block = 1,
      fleet = 1,
      mu = dat$mean_q[1] * exp(-dat$cv_q[1]^2 / 2),
      sd = dat$cv_q[1]
    ),
    t_srv = array(0.5, dim = c(dat$n_regions, dat$n_seas, dat$n_srv_fleets))
  )

  ## Weighting ----------------------------------------------------------------
  # the composition weights are the assessment's McAllister Ianelli multipliers,
  # and the selectivity penalty weights are its lambdas 3 to 6 plus the mean
  # centering term the template hardcodes for bicubic fishery selectivity.
  # survey selectivity is logistic with no deviations, so it has no
  # smoothness penalty
  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = dat$lam_rec,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = dat$Wt_FishAgeComps,
    Wt_FishLenComps = dat$Wt_FishLenComps,
    Wt_SrvAgeComps = dat$Wt_SrvAgeComps,
    Wt_SrvLenComps = dat$Wt_SrvLenComps,
    fish_sel_pen_wts = dat$sel_pen_wts,
    srv_sel_pen_wts = NULL
  )
} # end build_bsai_pop_input


#' Set every parameter to the assessment's maximum likelihood estimate
#'
#' Evaluating at a known point separates a specification error from an
#' optimization difference. Two blocks here need a conversion rather than a
#' direct assignment: the recruitment bias corrections and the bicubic node
#' grid's storage order.
seed_bsai_pop_mle <- function(input_list, dat) {

  mle <- dat$mle
  s2 <- dat$sigmaR^2 / 2

  ## Recruitment and the initial equilibrium ----------------------------------
  # the assessment builds its deviation-free terminal recruits as
  # exp(mean_log_rec + sigmaR^2 / 2) while the estimated years are raw, and it
  # starts the first year equilibrium from the MEAN of the lognormal rather than
  # exp(log_rinit). with both bias corrections in ln_global_R0 and ln_rinit
  # and shifting every seeded deviation down by the same amount reproduces the
  # recruitment series, the initial age structure and the penalty value
  input_list$par$ln_global_R0[] <- mle$mean_log_rec + s2
  input_list$par$ln_rinit <- mle$log_rinit + s2
  input_list$par$ln_RecDevs[1, 1, ] <- mle$rec_dev - s2

  ## Mortality and catchability -----------------------------------------------
  input_list$par$ln_M[] <- log(mle$M)
  input_list$par$ln_F_mean[] <- mle$log_avg_fmort
  input_list$par$ln_F_devs[1, , 1, 1] <- mle$fmort_dev
  input_list$par$ln_srv_q[1, 1, 1] <- log(mle$q_srv[1])
  input_list$par$ln_srv_q[1, 1, 2] <- log(mle$q_srv[2])

  ## Survey selectivity -------------------------------------------------------
  for(sf in 1:dat$n_srv_fleets) {
    input_list$par$srv_fixed_sel_pars[1, , 1, 1, sf] <-
      c(log(mle$sel_a50_srv[sf]), log(mle$sel_aslope_srv[sf]))
  } # end sf loop

  ## Fishery selectivity, the bicubic node grid -------------------------------
  # the assessment declares its node grid with year nodes as the outer dimension
  # and writes the parameter file row major, so the file's five rows are AGE
  # nodes. SPoRC builds the surface as exp(Wyr %*% nodes %*% t(Wbin)) and wants
  # [year node x age node] flattened column major, hence the transpose
  node_admb <- matrix(
    mle$fsh_sel_par,
    nrow = dat$fsh_age_nodes,
    ncol = dat$fsh_yr_nodes,
    byrow = TRUE
  )
  input_list$par$fish_fixed_sel_pars[1, , 1, 1, 1] <- as.vector(t(node_admb))

  input_list
} # end seed_bsai_pop_mle


#' Normalize a selectivity surface by its within-year maximum
#'
#' The assessment rescales for REPORTING only: selectivity is divided by its
#' within-year maximum and F multiplied by the same factor, leaving their product
#' invariant. Its internal selectivity is the raw exponentiated bicubic surface,
#' so a like-for-like comparison normalizes SPoRC's surface the same way and
#' divides the reported F by the same factor.
norm_bsai_pop_sel_by_year <- function(m) m / apply(m, 1, max)
