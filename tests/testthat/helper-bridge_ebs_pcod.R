# EBS Pacific cod bridge: the 2024 eastern Bering Sea Pacific cod assessment
# (Stock Synthesis Model 24.1, Barbeaux et al. 2024) rebuilt in SPoRC.
#
# The two functions here are the whole bridge. build_ebs_pcod_input() specifies
# the model, one Setup_Mod_* call per section in the order the case study
# vignette walks through them, and seed_ebs_pcod_mle() sets every parameter to
# the assessment's own maximum likelihood estimate so the model can be checked
# before it is ever optimized. The data list comes from
# dev/pcod_bridge/R/build_pcod_data.R and is packaged as sgl_rg_ebs_pcod_data.
#
# The model is one area, one sex, one season, ages 0-20, lengths on 121 one
# centimeter population bins reported to 24 five centimeter data bins, and years
# 1977-2024.
#
#   Source                     Years        Observations  Likelihood
#   Catch                      1977-2024    48            Lognormal
#   Survey abundance (numbers) 1982-2024    42            Lognormal, extra SD
#   Fishery length comps       1977-2024    48            Multinomial
#   Survey length comps        1982-2024    19            Multinomial
#   Survey age comps           2000-2023    23            Multinomial
#
# What the assessment does that a SPoRC default would not, all of it set below:
# Richards growth carried cohort by cohort with annual deviations on the length
# at age 1.5 and on K from 2000; length based selectivity applied at length
# rather than folded to age; a fishing fleet's compositions read at mid season
# whatever month the data carry; the catch in biomass on the selection weighted
# weight at age. See dev/pcod_bridge/README.md for the ones that were measured
# rather than assumed, and vignette("ae_ebs_pacific_cod_case_study") for the
# walkthrough.

#' Build the SPoRC input list for the EBS Pacific cod bridge
#'
#' @param dat the list from build_pcod_data() with the SS3 report attached
#' @keywords internal
build_ebs_pcod_input <- function(dat) {

  yrs <- dat$years; n_yrs <- length(yrs)
  ages <- dat$ages; n_ages <- length(ages)
  n_reg <- 1; n_sex <- 1; n_fish <- 1; n_srv <- 1
  n_lens <- length(dat$lens)
  sigmaR <- dat$rec$sigmaR

  ## Model dimensions ---------------------------------------------------------
  # one area, one sex, one season, so the region, sex and season subscripts all
  # collapse to one. lengths are the population bins the assessment grows fish
  # on, not the coarser bins the compositions are reported on
  input_list <- Setup_Mod_Dim(
    years = yrs, ages = ages, lens = dat$lens,
    n_regions = n_reg, n_sexes = n_sex, n_fish_fleets = n_fish, n_srv_fleets = n_srv,
    n_seas = 1, n_pop = 1, natal_region = 1, verbose = FALSE
  )

  ## Recruitment --------------------------------------------------------------
  # recruitment is a mean with annual deviations rather than a stock recruit
  # function, age 0 entering at the start of the year. the start year's ages
  # 1-20 are set by twenty early deviations, and the plus group's deviation is
  # estimated and penalized along with them, which is equil_init_age_strc =
  # "stoch_all". a separate initial equilibrium recruitment is penalized towards
  # R0 with the assessment's standard deviation, sigmaR over the average age,
  # and initial fishing mortality is held at the estimate because SPoRC has no
  # equilibrium catch to fit it to.
  # bias ramp years are given in deviation-index space, 1 being the first model year
  ave_age <- 1 / dat$growth$M - 0.5
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec", rec_dd = "global", rec_lag = 0, SR_ref_yr = 1, t_spawn = 0,
    sigmaR_spec = "fix", ln_sigmaR = array(log(sigmaR), dim = c(2, 1, n_reg)), sigmaR_switch = 1,
    do_rec_bias_ramp = 1,
    bias_year = dat$rec$bias_years - yrs[1] + 1,
    max_bias_ramp_fct = dat$rec$max_bias_adj,
    RecDevs_spec = "est_shared_pop_r", RecDevs_pen_center = "fixed", dont_est_recdev_last = 0,
    init_age_strc = 2, equil_init_age_strc = 2,
    InitDevs_spec = "est_shared_pop_r", InitDevs_pen_center = "fixed",
    rec_region_prop_spec = NULL,
    use_rinit = 1, Use_rinit_pen = 1, rinit_pen_sd = sigmaR / ave_age,
    init_F_form = "abs", init_F_spec = "fix", init_F_par = array(log(dat$mle$init_F), dim = c(n_reg, 1, n_fish)),
    ln_global_R0 = dat$mle$ln_R0, ln_rinit = dat$mle$ln_R0 + dat$mle$regime
  )

  ## Biological dynamics ------------------------------------------------------
  # natural mortality is fixed. growth is the Richards curve, a sixth parameter
  # raising the lengths to a power, with weight at age derived from it through
  # the weight-length relationship rather than supplied as data.
  #
  # two things here are specific to this assessment. growth_tv_type = "cohort"
  # carries size at age forward cohort by cohort instead of reading each year's
  # curve, so a cohort grows by the increment the current year's parameters
  # imply from the size it has already reached, and the plus group's size is the
  # numbers weighted blend of the cohort entering it with the fish already
  # there. growth_tv_link = "logit" makes each deviation an offset inside the
  # parameter's bounds rather than a multiplier, which is how Stock Synthesis
  # writes deviations on a bounded parameter. only L1 and K vary, so the sigma
  # vector carries their standard deviations and a placeholder for the four
  # parameters that do not.
  #
  # maturity at age is taken from the assessment rather than modeled, since
  # maturity is length based and fixed there. LenBinMap maps the 121 population
  # bins onto the 24 bins the compositions are recorded on.
  g <- dat$growth
  # the deviations' standard deviations sit in the first stream of the shared
  # process error array, one slot per growth parameter; the second stream is
  # unread with no semi-parametric surface
  pe_start <- array(log(0.05), dim = c(1, n_reg, max(4, n_ages, 6), n_sex, 2))
  pe_start[1, 1, 1:6, 1, 1] <- log(c(g$dev_sd[["L1"]], 1, g$dev_sd[["K"]], 1, 1, 1))
  MatAA <- array(0, dim = c(1, n_reg, n_yrs, 1, n_ages, n_sex))
  MatAA[1, 1, , 1, , 1] <- dat$MatAA
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = NULL, MatAA = MatAA,
    fit_lengths = 1, SizeAgeTrans = NA,
    AgeingError = dat$AgeingError,
    M_spec = "fix", Fixed_natmort = array(g$M, dim = c(1, n_reg, n_yrs, n_ages, n_sex)),
    addtocomp = dat$comp$addtocomp_len, comp_const_obs = 1, addtosrvidx = 0, addtofishidx = 0,
    growth_model = "richards",
    growth_spec = "est_all", growth_fix = !g$est[c("L1", "L2", "K", "CV1", "CV2", "rho")],
    ln_growth_pars = array(log(dat$mle$growth[c("L1", "L2", "K", "CV1", "CV2", "rho")]), dim = c(1, n_reg, n_sex, 6)),
    growth_A1 = g$A1,
    # the control file records a second reference age of 999 to mean that L2 is
    # the asymptote itself
    growth_A2 = if(g$A2 == 999) "Linf" else g$A2,
    growth_len_lower = dat$lens_lower, growth_L0 = dat$lens_lower[1],
    growth_plus_group = "mixture",
    growth_tv_model = c(L1 = "iid", K = "iid"), growth_tv_years = g$dev_years,
    growth_tv_link = "logit", growth_par_bounds = g$bounds[c("L1", "L2", "K", "CV1", "CV2", "rho"), ],
    growth_pe_pars = pe_start, growth_tv_sigma_spec = "fix",
    growth_tv_type = "cohort",
    waa_model = "wt_len", wt_len_pars = dat$wtlen,
    LenBinMap = dat$LenBinMap
  )

  ## Movement and tagging -----------------------------------------------------
  # one area, so there is nothing to move and nothing tagged
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  ## Catch and fishing mortality ----------------------------------------------
  # the assessment solves fishing mortality from the catch with its hybrid
  # method, so it spends no parameters on F. SPoRC estimates a deviation per
  # year against the assessment's catch error instead, with no penalty on the
  # deviations, which fits the catch essentially exactly
  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = dat$ObsCatch, UseCatch = dat$UseCatch,
    Use_F_pen = 0, ln_F_mean_spec = "fix",
    sigmaC_spec = "fix", ln_sigmaC = array(log(dat$catch_se), dim = c(n_reg, n_yrs, 1, n_fish))
  )

  joint <- function(n) paste0("spltRjntS_Year_1-terminal_Fleet_", seq_len(n))
  none <- function(n) paste0("none_Year_1-terminal_Fleet_", seq_len(n))

  ## Fishery compositions -----------------------------------------------------
  # no fishery index and no fishery ages, lengths only. two settings matter
  # here. FishLenComps_sel = "length" applies selectivity at length and sums
  # over ages, rather than folding selectivity to age first, which keeps the
  # covariance of length and selection within an age. fish_waa_selected = 1
  # puts the catch in biomass on the selection weighted weight at age.
  #
  # t_fish is 0.5 because Stock Synthesis assigns a fishing fleet's
  # compositions to mid season whatever month the data file records, so the
  # key, the selected weight and the compositions all sit there
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    t_fish = array(0.5, dim = c(n_reg, 1, n_fish)),
    FishLenComps_sel = "length", fish_waa_selected = 1,
    ObsFishIdx = array(NA_real_, dim = c(n_reg, n_yrs, 1, n_fish)),
    ObsFishIdx_SE = array(NA_real_, dim = c(n_reg, n_yrs, 1, n_fish)),
    UseFishIdx = array(0, dim = c(n_reg, n_yrs, 1, n_fish)),
    fish_idx_type = rep("none", n_fish), FishIdx_LikeType = rep("lognormal", n_fish),
    ObsFishAgeComps = array(NA_real_, dim = c(n_reg, n_yrs, 1, length(dat$obs_ages), n_sex, n_fish)),
    UseFishAgeComps = array(0, dim = c(n_reg, n_yrs, 1, n_fish)),
    ObsFishLenComps = dat$ObsFishLenComps, UseFishLenComps = dat$UseFishLenComps, ISS_FishLenComps = dat$ISS_FishLenComps,
    FishAgeComps_LikeType = rep("none", n_fish), FishLenComps_LikeType = rep("Multinomial", n_fish),
    FishAgeComps_Type = none(n_fish), FishLenComps_Type = joint(n_fish)
  )

  ## Survey index and compositions --------------------------------------------
  # the index is in numbers, so weight at age never enters it. the extra
  # standard deviation the assessment estimates is added to the observed
  # standard errors here rather than estimated, which is the one parameter of
  # the assessment's that SPoRC folds into data. ages come through the ageing
  # error definitions, lengths on the data bins, all read at mid year
  t_srv <- array(dat$t_srv, dim = c(n_reg, 1, n_srv))
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = dat$ObsSrvIdx, ObsSrvIdx_SE = dat$ObsSrvIdx_SE, UseSrvIdx = dat$UseSrvIdx,
    # The assessment's extra survey standard deviation, carried as a parameter
    # started at its own estimate. Seeded evaluation is identical to adding it to
    # the standard errors by hand; the free fit estimates it instead.
    sigmaSrvIdx_spec = "est_additive", ln_sigmaSrvIdx = log(dat$mle$extra_sd),
    srv_idx_type = rep("abd", n_srv), SrvIdx_LikeType = rep("lognormal", n_srv),
    ObsSrvAgeComps = dat$ObsSrvAgeComps, UseSrvAgeComps = dat$UseSrvAgeComps, ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
    ObsSrvLenComps = dat$ObsSrvLenComps, UseSrvLenComps = dat$UseSrvLenComps, ISS_SrvLenComps = dat$ISS_SrvLenComps,
    SrvAgeComps_LikeType = rep("Multinomial", n_srv), SrvLenComps_LikeType = rep("Multinomial", n_srv),
    SrvAgeComps_Type = joint(n_srv), SrvLenComps_Type = joint(n_srv),
    SrvLenComps_sel = "length"
  )

  ## Fishery selectivity ------------------------------------------------------
  # a double normal at length on the population bins, in two blocks. the
  # 1977-1989 block replaces the peak and the ascending width and leaves the
  # rest of the base block's parameters, so only two of the six are estimated
  # per block. dbnrml_raw leaves the ascending limb as a raw Gaussian rather
  # than rescaling it, and dbnrml_startbin anchors the limb at the first DATA
  # bin, bins below taking (b / b_start)^2 times the selectivity there
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_selex_type = "length",
    fish_sel_model = c("dbnrml_Fleet_1_Block_1", "dbnrml_Fleet_1_Block_2"),
    cont_tv_fish_sel = paste0("none_Fleet_", seq_len(n_fish)),
    # block years are given in year-index space, 1 being the first model year
    fish_sel_blocks = c(paste0("Block_1_Year_1-", dat$sel$blocks[[4]][2] - yrs[1] + 1, "_Fleet_1"),
                        paste0("Block_2_Year_", dat$sel$blocks[[4]][2] - yrs[1] + 2, "-terminal_Fleet_1")),
    fish_q_blocks = paste0("none_Fleet_", seq_len(n_fish)),
    fish_fixed_sel_pars_spec = rep("est_all", n_fish),
    fish_sel_dbnrml_raw = matrix(c(1, 0), n_fish, 2, byrow = TRUE),
    fish_sel_dbnrml_startbin = rep(dat$startbin, n_fish),
    fish_q_spec = rep("fix", n_fish)
  )

  ## Survey selectivity -------------------------------------------------------
  # the same double normal with independent annual deviations on the ascending
  # width from 1982, which is cont_tv_srv_sel = "iid" with the process error
  # standard deviation held at the assessment's value. catchability is estimated
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_selex_type = "length",
    srv_sel_model = paste0("dbnrml_Fleet_", seq_len(n_srv)),
    cont_tv_srv_sel = paste0("iid_Fleet_", seq_len(n_srv)),
    srv_sel_blocks = paste0("none_Fleet_", seq_len(n_srv)),
    srv_q_blocks = paste0("none_Fleet_", seq_len(n_srv)),
    srv_fixed_sel_pars_spec = rep("est_all", n_srv),
    srv_sel_devs_spec = rep("est_all", n_srv),
    srvsel_pe_pars_spec = rep("fix", n_srv),
    srv_sel_dbnrml_startbin = rep(dat$startbin, n_srv),
    srv_q_spec = rep("est_all", n_srv),
    t_srv = t_srv
  )

  ## Weighting ----------------------------------------------------------------
  # one Francis weight per composition type. Stock Synthesis adds a constant to
  # every observed and expected proportion and then renormalizes, so each
  # likelihood is the SPoRC value scaled by 1 + n_bins * constant; dividing the
  # weight by that factor absorbs it
  va <- dat$var_adj
  w_len_fish <- va$value[va$factor == 4 & va$fleet == dat$fish_fleets[1]]
  w_len_srv <- va$value[va$factor == 4 & va$fleet == dat$srv_fleets[1]]
  w_age_srv <- va$value[va$factor == 5 & va$fleet == dat$srv_fleets[1]]
  n_obs_lens <- ncol(dat$LenBinMap); n_obs_ages <- length(dat$obs_ages)
  per_fleet <- function(w, n_fl) { arr <- array(1, dim = c(n_reg, n_yrs, 1, n_sex, n_fl)); for(f in seq_len(n_fl)) arr[, , , , f] <- w; arr }
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
    Wt_FishAgeComps = per_fleet(0, n_fish),
    Wt_FishLenComps = per_fleet(w_len_fish / (1 + n_obs_lens * dat$comp$addtocomp_len), n_fish),
    Wt_SrvAgeComps = per_fleet(w_age_srv / (1 + n_obs_ages * dat$comp$addtocomp_age), n_srv),
    Wt_SrvLenComps = per_fleet(w_len_srv / (1 + n_obs_lens * dat$comp$addtocomp_len), n_srv)
  )

  input_list
} # end build_ebs_pcod_input


#' Set every parameter to the assessment's maximum likelihood estimate
#'
#' Seeding the model at the assessment's own estimate is what makes the bridge
#' checkable: every reported quantity and every likelihood component can be
#' compared before the optimizer is allowed to move anything. Each block below
#' also fixes the map, since which parameters are estimated is part of the
#' specification.
#'
#' @keywords internal
seed_ebs_pcod_mle <- function(input_list, dat) {

  yrs <- dat$years; n_yrs <- length(yrs)
  ages <- dat$ages; n_ages <- length(ages)
  sigmaR <- dat$rec$sigmaR

  ## Recruitment level --------------------------------------------------------
  input_list$par$ln_global_R0[] <- dat$mle$ln_R0
  input_list$par$ln_rinit[] <- dat$mle$ln_R0 + dat$mle$regime

  ## Recruitment deviations ---------------------------------------------------
  # SPoRC's deviation is the log of recruitment over the level, which carries
  # the bias correction inside it, so the assessment's realized ratio can be
  # used directly. the main deviations 1977-2022 are estimated; the two later
  # years carry the ramp's adjustment only and are held
  rr <- dat$mle$rec_ratio[as.character(yrs)]
  input_list$par$ln_RecDevs[] <- log(rr)
  main_yrs <- as.integer(names(dat$mle$main_recdev))
  map_rec <- array(NA_real_, dim = dim(input_list$par$ln_RecDevs))
  map_rec[1, 1, match(main_yrs, yrs)] <- seq_along(main_yrs)
  input_list$map$ln_RecDevs <- factor(map_rec)

  ## Initial age structure ----------------------------------------------------
  # the deviation for year styr - a lands on age a, less its own bias
  # correction, since the ramp is defined on calendar years. all twenty are
  # estimated and penalized, which is what the assessment does and what
  # equil_init_age_strc = "stoch_all" gives. under "stoch_plus_grp" SPoRC holds
  # the accumulator age's deviation at zero and leaves it out of the initial age
  # penalty, and overriding the map alone would estimate it without penalizing
  # it, which is a parameter with neither a penalty nor data behind it
  early_yrs <- as.integer(names(dat$mle$early_recdev))
  early_adj <- dat$mle$early_recdev - 0.5 * dat$mle$biasadj[as.character(early_yrs)] * sigmaR^2
  input_list$par$ln_InitDevs[] <- 0
  map_init <- array(NA_real_, dim = dim(input_list$par$ln_InitDevs))
  for(i in seq_along(early_yrs)) {
    a <- yrs[1] - early_yrs[i]
    input_list$par$ln_InitDevs[1, 1, a, 1] <- early_adj[i]
    map_init[1, 1, a, 1] <- i
  } # end i loop
  input_list$map$ln_InitDevs <- factor(map_init)

  ## Fishing mortality --------------------------------------------------------
  # a fixed mean plus a free deviation per year, which is how SPoRC carries an F
  # series the assessment solved rather than estimated
  lf <- log(dat$mle$Fmort[as.character(yrs)])
  input_list$par$ln_F_mean[] <- mean(lf)
  input_list$par$ln_F_devs[] <- 0
  input_list$par$ln_F_devs[1, , 1, 1] <- lf - mean(lf)
  map_F <- array(NA_real_, dim = dim(input_list$par$ln_F_devs))
  map_F[1, , 1, 1] <- seq_len(n_yrs)
  input_list$map$ln_F_devs <- factor(map_F)

  ## Growth -------------------------------------------------------------------
  # the assessment reports unit normal deviations, so each is multiplied by its
  # own standard deviation to give the offset on the logit scale SPoRC expects.
  # slot 1 is L1 and slot 3 is K
  input_list$par$ln_growth_pars[1, 1, 1, ] <- log(dat$mle$growth[c("L1", "L2", "K", "CV1", "CV2", "rho")])
  dev_yrs <- match(dat$growth$dev_years$L1, yrs)
  input_list$par$ln_growth_devs[1, 1, dev_yrs, 1, 1] <- dat$mle$growth_devs$L1 * dat$growth$dev_sd[["L1"]]
  input_list$par$ln_growth_devs[1, 1, match(dat$growth$dev_years$K, yrs), 3, 1] <- dat$mle$growth_devs$K * dat$growth$dev_sd[["K"]]

  ## Fishery selectivity ------------------------------------------------------
  # the 1977-1989 block takes the base block's six parameters with the peak and
  # the ascending width replaced, and only those two are estimated in each block
  fs <- dat$mle$sel$fishery
  block1 <- fs$base; block1[c(1, 3)] <- fs$block1
  input_list$par$fish_fixed_sel_pars[1, , 1, 1, 1] <- block1
  input_list$par$fish_fixed_sel_pars[1, , 2, 1, 1] <- fs$base
  map_fish <- array(NA_real_, dim = dim(input_list$par$fish_fixed_sel_pars))
  map_fish[1, 1, 1, 1, 1] <- 1; map_fish[1, 3, 1, 1, 1] <- 2; map_fish[1, 1, 2, 1, 1] <- 3; map_fish[1, 3, 2, 1, 1] <- 4
  input_list$map$fish_fixed_sel_pars <- factor(map_fish)

  ## Survey selectivity -------------------------------------------------------
  # the same two of six estimated, plus deviations on the ascending width from
  # 1982, again the assessment's unit normal deviations times their standard
  # deviation. the data copy of the map is set too, since the process error
  # penalty reads it to know which deviations it scores
  sv <- dat$mle$sel$survey
  input_list$par$srv_fixed_sel_pars[1, , 1, 1, 1] <- sv$base
  map_srv <- array(NA_real_, dim = dim(input_list$par$srv_fixed_sel_pars))
  map_srv[1, 1, 1, 1, 1] <- 1; map_srv[1, 3, 1, 1, 1] <- 2
  input_list$map$srv_fixed_sel_pars <- factor(map_srv)

  srv_p3 <- grep("SizeSel_P_3_survey", rownames(dat$sel$size))
  dev_start <- as.integer(dat$sel$size[srv_p3, "dev_minyr"])
  dev_end <- as.integer(dat$sel$size[srv_p3, "dev_maxyr"])
  sel_dev_yrs <- match(dev_start:dev_end, yrs)
  input_list$par$ln_srvsel_devs[] <- 0
  input_list$par$ln_srvsel_devs[1, sel_dev_yrs, 3, 1, 1] <- dat$mle$srv_sel_devs * sv$dev_sd
  map_sdev <- array(NA_real_, dim = dim(input_list$par$ln_srvsel_devs))
  map_sdev[1, sel_dev_yrs, 3, 1, 1] <- seq_along(sel_dev_yrs)
  input_list$map$ln_srvsel_devs <- factor(map_sdev)
  input_list$data$map_ln_srvsel_devs <- map_sdev
  input_list$par$srvsel_pe_pars[] <- 0
  input_list$par$srvsel_pe_pars[1, 3, 1, 1] <- log(sv$dev_sd)

  ## Catchability -------------------------------------------------------------
  input_list$par$ln_srv_q[] <- dat$mle$ln_q

  input_list
} # end seed_ebs_pcod_mle
