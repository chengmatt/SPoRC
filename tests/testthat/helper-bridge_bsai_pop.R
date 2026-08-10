# Shared setup for the BSAI Pacific ocean perch bridge tests. The bridge test
# evaluates this configuration at the ADMB maximum likelihood estimate without
# optimizing; the pinned test refits from it. Keeping the two on one builder
# means a specification change cannot move one without the other.

# Build the input_list for the 2024 assessment configuration.
build_bsai_pop_input <- function(dat) {

  yrs <- dat$years
  n_yrs <- length(yrs)

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
  # and the last three years take the mean outright. The first year sits in
  # equilibrium under a fixed historical F of 0.01, which the assessment holds
  # independent of the estimated mean F; init_F_form "abs" is what keeps the
  # two separate, so the initial age structure cannot be depleted by raising
  # the mean F. bias_year is indexed in deviation space rather than calendar
  # years, so this range puts the whole series in the fully bias corrected
  # limb, which centres the penalty on -sigmaR^2 / 2 to match the shifted
  # deviations the seeds carry.
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

  # Natural mortality is estimated under a lognormal prior. The prior median is
  # shifted by exp(-cv^2 / 2) so its mean is the assessment's 0.05. Maturity is
  # estimated inside the assessment template and absent from its report object,
  # so the fitted logistic is rebuilt and held fixed here.
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

  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  # The assessment writes its catch and F penalties as sums of squares with
  # weights 500 and 0.1, which are the same statements as normal likelihoods
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

  # Two survey fleets, the Aleutian Islands bottom trawl survey and the eastern
  # Bering Sea slope survey, both fit to biomass.
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

  # Fishery selectivity is a bicubic spline over a 5 year by 5 age node grid,
  # exponentiated with no normalisation. SelStyr and NSelBins are not cosmetic:
  # they set the year range and bin count the smoothness penalties normalise by,
  # and they impose the assessment's edge holds, flat before 1964 and flat
  # beyond age 40.
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

  # Survey selectivity is logistic for both fleets. Only the Aleutian Islands
  # survey carries a catchability prior; supplying a single row for fleet 1 is
  # what restricts it to that fleet.
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = paste0("none_Fleet_", 1:dat$n_srv_fleets),
    srv_sel_blocks = paste0("none_Fleet_", 1:dat$n_srv_fleets),
    srv_sel_model = paste0("logist1_Fleet_", 1:dat$n_srv_fleets),
    srv_q_blocks = paste0("none_Fleet_", 1:dat$n_srv_fleets),
    srv_fixed_sel_pars_spec = rep("est_all", dat$n_srv_fleets),
    srv_q_spec = rep("est_all", dat$n_srv_fleets),
    Use_srv_q_prior = 1,
    srv_q_prior = data.frame(region = 1, block = 1, fleet = 1,
                             mu = dat$mean_q[1] * exp(-dat$cv_q[1]^2 / 2), sd = dat$cv_q[1]),
    t_srv = array(0.5, dim = c(dat$n_regions, dat$n_seas, dat$n_srv_fleets))
  )

  # The composition weights are the assessment's McAllister Ianelli
  # multipliers, and the selectivity penalty weights are its lambdas 3 to 6
  # plus the mean centering term the template hardcodes for bicubic fishery
  # selectivity. Survey selectivity is logistic with no deviations, so it
  # carries no smoothness penalty.
  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1,
    Wt_Rec = dat$lam_rec, Wt_F = 1, Wt_Tagging = 0,
    Wt_FishAgeComps = dat$Wt_FishAgeComps,
    Wt_FishLenComps = dat$Wt_FishLenComps,
    Wt_SrvAgeComps = dat$Wt_SrvAgeComps,
    Wt_SrvLenComps = dat$Wt_SrvLenComps,
    fish_sel_pen_wts = dat$sel_pen_wts,
    srv_sel_pen_wts = NULL
  )
}

# Set every parameter to the assessment's maximum likelihood estimate.
seed_bsai_pop_mle <- function(input_list, dat) {

  mle <- dat$mle
  s2 <- dat$sigmaR^2 / 2

  # The assessment builds its three dev-free terminal recruits as
  # exp(mean_log_rec + sigmaR^2 / 2) while the estimated years are raw, and it
  # starts the first year equilibrium from the mean of the lognormal rather
  # than exp(log_rinit). Carrying both bias corrections in ln_global_R0 and
  # ln_rinit and shifting every seeded deviation down by the same amount
  # reproduces the recruitment series, the initial age structure, and the
  # penalty value at the seed point.
  input_list$par$ln_global_R0[] <- mle$mean_log_rec + s2
  input_list$par$ln_rinit <- mle$log_rinit + s2
  input_list$par$ln_RecDevs[1, 1, ] <- mle$rec_dev - s2
  input_list$par$ln_M[] <- log(mle$M)
  input_list$par$ln_F_mean[] <- mle$log_avg_fmort
  input_list$par$ln_F_devs[1, , 1, 1] <- mle$fmort_dev
  input_list$par$ln_srv_q[1, 1, 1] <- log(mle$q_srv[1])
  input_list$par$ln_srv_q[1, 1, 2] <- log(mle$q_srv[2])

  for(sf in 1:dat$n_srv_fleets) {
    input_list$par$srv_fixed_sel_pars[1, , 1, 1, sf] <-
      c(log(mle$sel_a50_srv[sf]), log(mle$sel_aslope_srv[sf]))
  } # end sf loop

  # The assessment declares its bicubic node grid with year nodes as the outer
  # dimension and writes the parameter file row major, so the file's five rows
  # are age nodes. SPoRC builds the surface as exp(Wyr %*% nodes %*% t(Wbin))
  # and wants [year node x age node] flattened column major, hence the
  # transpose.
  node_admb <- matrix(mle$fsh_sel_par, nrow = dat$fsh_age_nodes,
                      ncol = dat$fsh_yr_nodes, byrow = TRUE)
  input_list$par$fish_fixed_sel_pars[1, , 1, 1, 1] <- as.vector(t(node_admb))

  input_list
}

# The assessment rescales for reporting only: selectivity is divided by its
# within-year maximum and F multiplied by the same factor, leaving their
# product invariant. Its internal selectivity is the raw exponentiated bicubic
# surface, so a like-for-like comparison normalises SPoRC's surface the same
# way and divides the reported F by the same factor.
norm_bsai_pop_sel_by_year <- function(m) m / apply(m, 1, max)
