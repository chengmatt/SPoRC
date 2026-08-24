# West Coast sablefish bridge: the 2025 West Coast sablefish assessment (Stock
# Synthesis 3) rebuilt in SPoRC.
#
# One Setup_Mod_* call per section, in the order vignette(
# "ac_wc_sablefish_case_study") walks through them, with a reason for each
# argument that follows the assessment rather than a SPoRC default.
#
# The model is one area, two sexes, one season, ages 0-70 with compositions
# binned to 0-50, years 1890-2024. Recruitment is age 0 off the same year's
# spawning biomass on a Beverton-Holt curve with steepness fixed at 0.75 and
# sigmaR 1.4. Six catch fleets, four trawl surveys, and an index built on the
# recruitment deviations themselves.
#
#   Source                                Years        Observations  Likelihood
#   Catch, three gears and their discards 1890-2024    620           Lognormal, CV 0.01
#   Four trawl survey indices             1980-2024     35           Lognormal
#   Recruitment index                     2020-2024      5           Normal on the deviation
#   Fishery age comps, six fleets         1983-2024    198           Multinomial
#   Survey age comps, four surveys        1983-2024     52           Multinomial
#
# test-regression_wc_sablefish_bridge.R evaluates this configuration at the
# Stock Synthesis maximum likelihood estimate without optimizing, and the case
# study figures refit from it, so a specification change cannot move one without
# the other.

#' Expand the assessment's block selectivity into SPoRC's fixed input array
#'
#' The assessment stores selectivity as its distinct block rows plus the year
#' each block starts. SPoRC takes a year by age by sex by fleet array, so each
#' year is looked up against the block it falls in.
expand_wc_sablefish_sel <- function(blocks, n_yrs, n_ages, n_sexes) {
  n_fleets <- length(blocks)
  out <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, n_fleets))
  for(f in seq_len(n_fleets)) {
    blk <- findInterval(seq_len(n_yrs), blocks[[f]]$blk_yr)
    out[1,1,,1,,,f] <- blocks[[f]]$sel[blk, , , drop = FALSE]
  } # end f loop
  out
} # end expand_wc_sablefish_sel

# Build the input_list for the 2025 assessment configuration.
build_wc_sablefish_input <- function(dat) {

  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_fish <- dat$n_fish_fleets
  n_srv <- dat$n_srv_fleets

  ## Model dimensions ---------------------------------------------------------
  # one area and one season, so those subscripts collapse, but the sex subscript
  # is real: growth, selectivity and spawning all differ by sex
  input_list <- Setup_Mod_Dim(
    years = yrs,
    ages = dat$ages,
    lens = dat$lens,
    n_regions = dat$n_regions,
    n_sexes = dat$n_sexes,
    n_fish_fleets = n_fish,
    n_srv_fleets = n_srv,
    n_seas = dat$n_seas,
    n_pop = dat$n_pop,
    natal_region = dat$natal_region,
    verbose = FALSE
  )

  ## Recruitment --------------------------------------------------------------
  # age 0 off the year's OWN spawning biomass, which is rec_lag = 0, on a
  # Beverton-Holt curve at fixed steepness, started from the unfished
  # equilibrium. the bias ramp is indexed in deviation space rather than in
  # calendar years. InitDevs are fixed because the assessment starts in 1890 at
  # equilibrium and estimates no initial age deviations
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "bh_rec",
    rec_lag = 0,
    SR_ref_yr = 1,
    t_spawn = dat$t_spawn,
    h_spec = "fix",
    steepness_h = array(stats::qlogis((dat$steepness - 0.2) / 0.8), dim = c(1, 1)),
    sigmaR_spec = "fix",
    ln_sigmaR = array(log(dat$sigmaR), dim = c(2, 1, 1)),
    sigmaR_switch = 1,
    do_rec_bias_ramp = 1,
    bias_year = dat$bias_year,
    max_bias_ramp_fct = dat$max_bias_adj,
    RecDevs_pen_center = "fixed",
    dont_est_recdev_last = 0,
    init_age_strc = 2,
    equil_init_age_strc = 0,
    InitDevs_spec = "fix",
    ln_global_R0 = dat$mle$ln_R0
  )

  ## Biological dynamics ------------------------------------------------------
  # one natural mortality for both sexes under the assessment's lognormal prior,
  # empirical weight and fecundity at age rather than a growth curve, and the
  # composition constant added on both sides of the multinomial. fit_lengths = 0
  # because this assessment fits ages only
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = dat$WAA,
    WAA_fish = array(dat$WAA, dim = c(dim(dat$WAA), n_fish)),
    WAA_srv = array(dat$WAA, dim = c(dim(dat$WAA), n_srv)),
    MatAA = dat$MatAA,
    fit_lengths = 0,
    AgeingError = dat$AgeingError,
    M_spec = "est_ln_M",
    Use_M_prior = 1,
    M_prior = dat$M_prior,
    addtocomp = dat$addtocomp,
    comp_const_obs = 1,
    addtosrvidx = 0,
    addtofishidx = 0
  )

  ## Movement and tagging -----------------------------------------------------
  # one area, so movement is the identity and nothing is tagged
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  ## Catch and fishing mortality ----------------------------------------------
  # the assessment solves fishing mortality from the catch rather than
  # estimating it, so there is no penalty on it and no mean to estimate. SPoRC
  # fits it against the assessment's own CV of 0.01
  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = dat$ObsCatch,
    UseCatch = dat$UseCatch,
    Use_F_pen = 0,
    ln_F_mean_spec = "fix",
    sigmaC_spec = "fix",
    ln_sigmaC = array(log(dat$sigmaC), dim = c(1, n_yrs, 1, n_fish))
  )

  ## Fishery compositions -----------------------------------------------------
  # ages only, no index and no lengths. a fleet whose data are sexed takes the
  # split-region joint-sex form and an unsexed fleet takes the aggregated one,
  # which is what dat$fish_sex records
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = array(NA_real_, dim = c(1, n_yrs, 1, n_fish)),
    ObsFishIdx_SE = array(NA_real_, dim = c(1, n_yrs, 1, n_fish)),
    UseFishIdx = array(0, dim = c(1, n_yrs, 1, n_fish)),
    fish_idx_type = rep("none", n_fish),
    FishIdx_LikeType = rep("lognormal", n_fish),
    ObsFishAgeComps = dat$ObsFishAgeComps,
    UseFishAgeComps = dat$UseFishAgeComps,
    ISS_FishAgeComps = dat$ISS_FishAgeComps,
    ObsFishLenComps = array(NA_real_, dim = c(1, n_yrs, 1, 1, dat$n_sexes, n_fish)),
    UseFishLenComps = array(0, dim = c(1, n_yrs, 1, n_fish)),
    ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, dat$n_sexes, n_fish)),
    FishAgeComps_LikeType = rep("Multinomial", n_fish),
    FishLenComps_LikeType = rep("none", n_fish),
    FishAgeComps_Type = paste0(ifelse(dat$fish_sex == 3, "spltRjntS", "agg"),
                               "_Year_1-terminal_Fleet_", seq_len(n_fish)),
    FishLenComps_Type = paste0("none_Year_1-terminal_Fleet_", seq_len(n_fish))
  )

  ## Survey index and compositions --------------------------------------------
  # six survey fleets doing three different jobs
  t_srv <- array(dat$t_srv, dim = c(1, 1, n_srv))

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = dat$ObsSrvIdx,
    ObsSrvIdx_SE = dat$ObsSrvIdx_SE,
    UseSrvIdx = dat$UseSrvIdx,
    # fleets 1-4 are the trawl surveys, 5 carries the unsexed compositions of
    # the last of them and no index, and 6 is the recruitment index, which
    # observes the deviations themselves under a normal likelihood
    srv_idx_type = c(rep("biom", 4), "none", "recdev"),
    SrvIdx_LikeType = c(rep("lognormal", n_srv - 1), "normal"),
    ObsSrvAgeComps = dat$ObsSrvAgeComps,
    UseSrvAgeComps = dat$UseSrvAgeComps,
    ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
    ObsSrvLenComps = array(NA_real_, dim = c(1, n_yrs, 1, 1, dat$n_sexes, n_srv)),
    UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
    ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, dat$n_sexes, n_srv)),
    SrvAgeComps_LikeType = rep("Multinomial", n_srv),
    SrvLenComps_LikeType = rep("none", n_srv),
    SrvAgeComps_Type = paste0(ifelse(dat$srv_sex == 3, "spltRjntS", "agg"),
                              "_Year_1-terminal_Fleet_", seq_len(n_srv)),
    SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", seq_len(n_srv)),
    t_srv = t_srv
  )

  ## Fishery selectivity ------------------------------------------------------
  # every fleet is on the age based double normal. the trawl fleet and the two
  # hook and line fleets carry time blocks, the pot fleet mirrors hook and line
  # including its male offsets, and each composition twin mirrors its parent.
  # blocks are named by the years the assessment's own surfaces change, which is
  # what its block patterns come to
  blk_string <- function(blocks, fleet) {
    st <- blocks[[fleet]]$blk_yr
    en <- c(st[-1] - 1, n_yrs)
    paste0("Block_", seq_along(st), "_Year_", st, "-", ifelse(en == n_yrs, "terminal", en), "_Fleet_", fleet)
  }

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = paste0("dbnrml_Fleet_", seq_len(n_fish)),
    cont_tv_fish_sel = paste0("none_Fleet_", seq_len(n_fish)),
    fish_sel_blocks = unlist(lapply(seq_len(n_fish), function(f) blk_string(dat$fish_sel_blocks_ss3, f))),
    fish_q_blocks = paste0("none_Fleet_", seq_len(n_fish)),
    # the pot fleet mirrors hook and line, and the unsexed trawl twin mirrors trawl
    fish_fixed_sel_pars_spec = replace(rep("est_all", n_fish), c(3, n_fish), c("est_shared_f_2", "est_shared_f_1")),
    # male parameters are offsets on the female's; the hook and line and pot
    # fleets additionally carry their own apical selectivity, which the male
    # limbs are built up to
    fish_sel_sex_offset = replace(rep("par", n_fish), c(2, 3), "par_apical"),
    fish_q_spec = rep("fix", n_fish)
  )

  ## Survey selectivity and catchability --------------------------------------
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = paste0("dbnrml_Fleet_", seq_len(n_srv)),
    cont_tv_srv_sel = paste0("none_Fleet_", seq_len(n_srv)),
    srv_sel_blocks = unlist(lapply(seq_len(n_srv), function(sf) blk_string(dat$srv_sel_blocks_ss3, sf))),
    srv_q_blocks = paste0("none_Fleet_", seq_len(n_srv)),
    # the unsexed twin shares the curve of the survey it belongs to, and the
    # recruitment index reads no curve at all
    srv_fixed_sel_pars_spec = c(rep("est_all", 4), "est_shared_f_4", "fix"),
    srv_sel_sex_offset = rep("par", n_srv),
    # catchability floats in the assessment, which is the same optimum as
    # estimating it; the unsexed twin carries no index, and the recruitment
    # index has a catchability of its own
    srv_q_spec = c(rep("est_all", 4), "fix", "est_all"),
    t_srv = t_srv
  )

  ## Weighting ----------------------------------------------------------------
  # the assessment applies no variance adjustments here, so every source carries
  # a weight of one
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, dat$n_sexes, n_fish)),
    Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, dat$n_sexes, n_fish)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, dat$n_sexes, n_srv)),
    Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, dat$n_sexes, n_srv))
  )

  ## Age 0 spawners in the equilibrium year -----------------------------------
  # spawning biomass is formed before the year's recruits settle, so the age 0
  # cell is empty in every year the model runs. the unfished equilibrium is
  # different: it already holds R0 / 2 at age 0, and that is what the first
  # year's spawning biomass and S0 are built from. setup refuses a non-zero
  # maturity at the recruit age under rec_lag = 0, so this cell is set on the
  # data list afterwards
  input_list$data$MatAA[1,1,1,1,1,1] <- dat$mat_age0_yr1

  input_list
} # end build_wc_sablefish_input


#' Set every parameter to the assessment's maximum likelihood estimate
#'
#' Also maps off the deviations the assessment does not estimate, since which
#' parameters are free is part of the specification.
seed_wc_sablefish_mle <- function(input_list, dat) {

  yrs <- dat$years
  n_yrs <- length(yrs)
  n_fish <- dat$n_fish_fleets
  n_srv <- dat$n_srv_fleets

  ## Recruitment level, mortality and catchability -----------------------------
  input_list$par$ln_global_R0[] <- dat$mle$ln_R0
  input_list$par$ln_M[] <- log(dat$mle$M)
  input_list$par$ln_srv_q[1, 1, seq_along(dat$mle$ln_srv_q)] <- dat$mle$ln_srv_q
  # the recruitment index's catchability is on the natural scale in the
  # assessment and on the log scale here
  input_list$par$ln_srv_q[1, 1, dat$n_srv_fleets] <- log(dat$mle$q_rec_idx)

  ## Recruitment deviations ---------------------------------------------------
  # SPoRC's deviation is the assessment's less the bias correction. the early
  # years sit at a phase the assessment never estimates, so they stay at zero
  input_list$par$ln_RecDevs[1, 1, ] <- dat$mle$recdev - 0.5 * dat$mle$bias_adj * dat$sigmaR^2
  map_rec <- rep(NA_real_, n_yrs)
  est <- yrs %in% dat$yrs_rec_est
  map_rec[est] <- seq_len(sum(est))
  input_list$map$ln_RecDevs <- factor(map_rec)

  ## Fishing mortality --------------------------------------------------------
  # a fixed mean plus deviations per fleet. closures carry no deviation, and the
  # trawl fleet's composition twin sits at a rate low enough to leave the numbers
  # at age untouched
  for(f in seq_len(ncol(dat$mle$Fmort))) {
    has <- dat$UseCatch[1, , 1, f] == 1
    lf <- log(dat$mle$Fmort[has, f])
    input_list$par$ln_F_mean[1, 1, f] <- mean(lf)
    input_list$par$ln_F_devs[1, has, 1, f] <- lf - mean(lf)
  } # end f loop
  input_list$par$ln_F_mean[1, 1, n_fish] <- log(1e-12)
  map_F <- array(NA_real_, dim = dim(input_list$par$ln_F_devs))
  map_F[dat$UseCatch == 1] <- seq_len(sum(dat$UseCatch == 1))
  input_list$map$ln_F_devs <- factor(map_F)

  ## Selectivity, the female curves -------------------------------------------
  # every one of the six double normal parameters is on the same scale in SPoRC
  # as in the assessment, so its values go straight in. a parameter is estimated
  # where the assessment estimates it, and blocks drawing on the same underlying
  # parameter share one level, which is what src_id records
  lev <- 0
  map_fish <- array(NA_real_, dim = dim(input_list$par$fish_fixed_sel_pars))
  map_srv <- array(NA_real_, dim = dim(input_list$par$srv_fixed_sel_pars))
  seen <- character(0)
  put <- function(map, tab, f) {
    for(b in seq_len(ncol(tab$pars))) {
      for(k in 1:6) {
        if(!tab$est[k, b]) next
        id <- tab$src_id[k, b]
        if(!id %in% seen) { seen[[length(seen) + 1]] <<- id; lev <<- lev + 1 }
        map[1, k, b, 1, f] <- match(id, seen) + 0
      } # end k loop
    } # end b loop
    map
  }
  for(f in c(1, 2, 4, 5, 6)) { # the pot fleet mirrors hook and line and is copied below
    tab <- dat$sel_fish[[f]]
    input_list$par$fish_fixed_sel_pars[1, , seq_len(ncol(tab$pars)), 1, f] <- tab$pars
    map_fish <- put(map_fish, tab, f)
  } # end f loop
  for(sf in seq_len(4)) { # 5 mirrors 4 and 6 reads no curve
    tab <- dat$sel_srv[[sf]]
    input_list$par$srv_fixed_sel_pars[1, , seq_len(ncol(tab$pars)), 1, sf] <- tab$pars
    map_srv <- put(map_srv, tab, sf)
  } # end sf loop

  ## Selectivity, the hook and line male curve --------------------------------
  # written exactly as the assessment writes it: an offset on the peak in bins,
  # an offset on the selectivity at the last bin, and its own apical selectivity
  # that the two limbs are built up to. the assessment gives the male ascending
  # and descending widths no offset, and reads the female's parameter for the
  # selectivity at the first bin, so those three offsets stay at zero
  A <- dat$sel_male$value[["scale"]]
  tab2 <- dat$sel_fish[[2]]
  n_blk2 <- ncol(tab2$pars)
  for(b in seq_len(n_blk2)) {
    input_list$par$fish_fixed_sel_pars[1, 1, b, 2, 2] <- dat$sel_male$value[["peak"]]
    input_list$par$fish_fixed_sel_pars[1, 6, b, 2, 2] <- dat$sel_male$value[["final"]]
  } # end b loop
  if(dat$sel_male$est[["peak"]]) { lev <- lev + 1; map_fish[1, 1, seq_len(n_blk2), 2, 2] <- lev }
  if(dat$sel_male$est[["final"]]) { lev <- lev + 1; map_fish[1, 6, seq_len(n_blk2), 2, 2] <- lev }

  input_list$par$ln_fishsel_sex_scale[] <- 0
  input_list$par$ln_fishsel_sex_scale[1, seq_len(n_blk2), 2, 2:3] <- log(A)
  map_scale <- array(NA_real_, dim = dim(input_list$par$ln_fishsel_sex_scale))
  if(dat$sel_male$est[["scale"]]) map_scale[1, seq_len(n_blk2), 2, 2:3] <- 1

  ## Mirrored fleets ----------------------------------------------------------
  # a mirrored fleet carries the same values AND the same map levels, so it is
  # the same parameter rather than a copy that could drift
  input_list$par$fish_fixed_sel_pars[1, , , , 3] <- input_list$par$fish_fixed_sel_pars[1, , , , 2]
  map_fish[1, , , , 3] <- map_fish[1, , , , 2]
  input_list$par$fish_fixed_sel_pars[1, , , , n_fish] <- input_list$par$fish_fixed_sel_pars[1, , , , 1]
  map_fish[1, , , , n_fish] <- map_fish[1, , , , 1]
  input_list$par$srv_fixed_sel_pars[1, , , , 5] <- input_list$par$srv_fixed_sel_pars[1, , , , 4]
  map_srv[1, , , , 5] <- map_srv[1, , , , 4]

  input_list$map$fish_fixed_sel_pars <- factor(map_fish)
  input_list$map$srv_fixed_sel_pars <- factor(map_srv)
  input_list$map$ln_fishsel_sex_scale <- factor(map_scale)

  input_list
} # end seed_wc_sablefish_mle
