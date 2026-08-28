# GOA rex sole bridge: Model 25.1 of the 2025 Gulf of Alaska rex sole assessment
# (Stock Synthesis 3) rebuilt in SPoRC and evaluated at the assessment's own
# maximum likelihood estimate.
#
# One Setup_Mod_* call per section, in the order vignette(
# "ad_goa_rex_sole_case_study") walks through them, with a reason for each
# argument that follows the assessment rather than a SPoRC default. The data
# list comes from dev/rex_bridge/R/build_rex_data.R.
#
# This is the bridge that motivated SPoRC's parametric growth module and its
# conditional age-at-length likelihood: growth is estimated separately in each
# area and sex from the survey's conditional age-at-length data.
#
# The model is two areas (Western-Central = 1, Eastern = 2), two sexes, one
# season, ages 0-20 with ages observed from 1, lengths 29 two centimeter bins
# from 10 to 66 cm, years 1982-2024, sigmaR 0.6.
#
#   Source                          Years        Observations  Likelihood
#   Catch, Western-Central          1982-2024    43            Lognormal, CV 0.01
#   Survey biomass, both areas      1990-2023    31            Lognormal
#   Fishery length comps            1982-2024    21            Multinomial
#   Fishery age comps               1992-2022    18            Multinomial
#   Survey length comps, both       1990-2023    31            Multinomial
#   Survey conditional age-at-length 1993-2019   482 rows      Multinomial
#
# The Eastern area carries no catch, so its population is shaped entirely by the
# survey.

#' Build the SPoRC input list for the rex sole bridge
#'
#' @param dat the list from build_rex_data() with the SS3 report attached
#' @keywords internal
build_goa_rex_input <- function(dat) {

  yrs <- dat$years; n_yrs <- length(yrs)
  ages <- dat$ages; n_ages <- length(ages)
  n_reg <- dat$n_regions; n_sex <- dat$n_sexes
  n_fish <- dat$n_fish_fleets; n_srv <- dat$n_srv_fleets
  n_lens <- length(dat$lens)

  ## Model dimensions ---------------------------------------------------------
  # the assessment's areas are separate growth patterns sharing one recruitment
  # series with nothing moving between them, so SPoRC carries them as two
  # regions of one population with the identity movement matrix. the regional
  # subscripts are real here and the movement ones are not
  input_list <- Setup_Mod_Dim(
    years = yrs, ages = ages, lens = dat$lens,
    n_regions = n_reg, n_sexes = n_sex, n_fish_fleets = n_fish, n_srv_fleets = n_srv,
    n_seas = 1, n_pop = 1, natal_region = 1, verbose = FALSE
  )

  ## Recruitment --------------------------------------------------------------
  # age 0 at the start of the year, mean recruitment apportioned across the areas
  # by an estimated logit, the bias ramp on the deviations, and start-year ages
  # 1-17 set by early deviations shared across areas and sexes.
  # bias ramp years are given in deviation-index space, 1 being the first model
  # year, and each early deviation reads the ramp at its own birth year
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec", rec_dd = "global", rec_lag = 0, SR_ref_yr = 1, t_spawn = 0,
    sigmaR_spec = "fix", ln_sigmaR = array(log(dat$rec$sigmaR), dim = c(2, 1, n_reg)), sigmaR_switch = 1,
    do_rec_bias_ramp = 1,
    bias_year = dat$rec$bias_years - yrs[1] + 1,
    max_bias_ramp_fct = dat$rec$max_bias_adj,
    RecDevs_spec = "est_shared_pop_r", RecDevs_pen_center = "fixed", dont_est_recdev_last = 0,
    init_age_strc = 2, equil_init_age_strc = 1,
    InitDevs_spec = "est_shared_pop_r", InitDevs_sex_spec = "est_shared_s", InitDevs_pen_center = "fixed",
    rec_region_prop_spec = NULL,
    ln_global_R0 = dat$mle$ln_R0
  )

  ## Biological dynamics ------------------------------------------------------
  # natural mortality is fixed, and maturity is logistic on age for females with
  # none below the first mature age.
  #
  # growth is the Schnute-Francis form with a linear start: mean length at age is
  # L1 at A1 and L2 at A2, von Bertalanffy between them, and linear from the
  # first length bin's lower edge at age 0 up to L1. the coefficient of variation
  # is linear in mean length between the two reference ages and flat outside
  # them. the plus group is not the curve at age 20 but an exponentially weighted
  # mixture of the ages it holds, the rule the assessment inherits from SS3.24,
  # which growth_plus_group = "mixture" applies. weight at age is the age-length
  # key times weight at the bin midpoints, which is waa_model = "wt_len"
  MatAA <- array(0, dim = c(1, n_reg, n_yrs, 1, n_ages, n_sex))
  mat_f <- 1 / (1 + exp(dat$mat$slope * (ages - dat$mat$a50)))
  mat_f[ages < dat$mat$first_mature_age] <- 0
  for(r in 1:n_reg) for(y in 1:n_yrs) MatAA[1, r, y, 1, , 1] <- mat_f
  growth_start <- dat$mle$growth # [1, area, sex, 5]
  wl <- array(NA_real_, dim = c(1, n_reg, n_sex, 2))
  for(r in 1:n_reg) { wl[1, r, 1, ] <- dat$wtlen$fem; if(n_sex > 1) wl[1, r, 2, ] <- dat$wtlen$mal }

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = NULL, MatAA = MatAA,
    fit_lengths = 1, SizeAgeTrans = NA,
    AgeingError = dat$AgeingError,
    M_spec = "fix", Fixed_natmort = array(dat$growth[[1]]$M, dim = c(1, n_reg, n_yrs, n_ages, n_sex)),
    addtocomp = dat$comp$addtocomp_age, comp_const_obs = 1, addtosrvidx = 0, addtofishidx = 0,
    growth_model = "vb_schnute", growth_spec = "est_all",
    ln_growth_pars = log(growth_start),
    growth_A1 = dat$growth_A1, growth_A2 = dat$growth_A2,
    growth_len_lower = dat$lens_lower, growth_L0 = dat$lens_lower[1],
    growth_plus_group = "mixture",
    waa_model = "wt_len", wt_len_pars = wl
  )

  ## Movement and tagging -----------------------------------------------------
  # nothing moves between the two areas, so movement is the identity, and
  # nothing is tagged
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  ## Catch and fishing mortality ----------------------------------------------
  # the assessment solves fishing mortality from the catch with its hybrid
  # method, so it spends no parameters on F. SPoRC estimates a deviation per year
  # against a tight catch error with no penalty on the deviations
  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = dat$ObsCatch, UseCatch = dat$UseCatch,
    Use_F_pen = 0, ln_F_mean_spec = "fix",
    sigmaC_spec = "fix", ln_sigmaC = array(log(dat$catch_se_value), dim = c(n_reg, n_yrs, 1, n_fish))
  )

  joint <- function(n) paste0("spltRjntS_Year_1-terminal_Fleet_", seq_len(n))
  none <- function(n) paste0("none_Year_1-terminal_Fleet_", seq_len(n))

  ## Fishery compositions -----------------------------------------------------
  # no fishery index; joint-sex marginal ages and lengths, and no conditional
  # age-at-length. t_fish = 0.5 reads the fishery's key and weight at mid season,
  # where the assessment reads them
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    t_fish = array(0.5, dim = c(n_reg, 1, n_fish)),
    ObsFishIdx = array(NA_real_, dim = c(n_reg, n_yrs, 1, n_fish)),
    ObsFishIdx_SE = array(NA_real_, dim = c(n_reg, n_yrs, 1, n_fish)),
    UseFishIdx = array(0, dim = c(n_reg, n_yrs, 1, n_fish)),
    fish_idx_type = rep("none", n_fish), FishIdx_LikeType = rep("lognormal", n_fish),
    ObsFishAgeComps = dat$ObsFishAgeComps, UseFishAgeComps = dat$UseFishAgeComps, ISS_FishAgeComps = dat$ISS_FishAgeComps,
    ObsFishLenComps = dat$ObsFishLenComps, UseFishLenComps = dat$UseFishLenComps, ISS_FishLenComps = dat$ISS_FishLenComps,
    FishAgeComps_LikeType = rep("Multinomial", n_fish), FishLenComps_LikeType = rep("Multinomial", n_fish),
    FishAgeComps_Type = joint(n_fish), FishLenComps_Type = joint(n_fish)
  )

  ## Survey indices and compositions ------------------------------------------
  # biomass indices and joint-sex lengths, plus the conditional age-at-length
  # that carries the age information. a CAAL row is an age composition WITHIN a
  # length bin: each length bin of a survey year holds the ages of the otoliths
  # read from that bin, fit as its own multinomial with the number aged as its
  # sample size, one sex per row. the assessment carries the marginal survey ages
  # as ghosts, so they are read in but not fit
  t_srv <- array(rep(dat$t_srv, each = n_reg), dim = c(n_reg, 1, n_srv))
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = dat$ObsSrvIdx, ObsSrvIdx_SE = dat$ObsSrvIdx_SE, UseSrvIdx = dat$UseSrvIdx,
    srv_idx_type = rep("biom", n_srv), SrvIdx_LikeType = rep("lognormal", n_srv),
    ObsSrvAgeComps = dat$ObsSrvAgeComps, UseSrvAgeComps = array(0, dim = dim(dat$UseSrvAgeComps)), ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
    ObsSrvLenComps = dat$ObsSrvLenComps, UseSrvLenComps = dat$UseSrvLenComps, ISS_SrvLenComps = dat$ISS_SrvLenComps,
    SrvAgeComps_LikeType = rep("none", n_srv), SrvLenComps_LikeType = rep("Multinomial", n_srv),
    SrvAgeComps_Type = none(n_srv), SrvLenComps_Type = joint(n_srv),
    ObsSrv_caal = dat$ObsSrv_caal, UseSrv_caal = dat$UseSrv_caal, ISS_Srv_caal = dat$ISS_Srv_caal,
    Srv_caal_LikeType = rep("Multinomial", n_srv),
    Srv_caal_Type = paste0("spltRspltS_Year_1-terminal_Fleet_", seq_len(n_srv)),
    t_srv = t_srv
  )

  ## Fishery selectivity ------------------------------------------------------
  # an age based double normal with male offsets on the parameters, which is
  # fish_sel_sex_offset = "par"
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = paste0("dbnrml_Fleet_", seq_len(n_fish)),
    cont_tv_fish_sel = paste0("none_Fleet_", seq_len(n_fish)),
    fish_sel_blocks = paste0("none_Fleet_", seq_len(n_fish)),
    fish_q_blocks = paste0("none_Fleet_", seq_len(n_fish)),
    fish_fixed_sel_pars_spec = rep("est_all", n_fish),
    fish_sel_sex_offset = rep("par", n_fish),
    # the assessment leaves the ascending limb unanchored (p5 = -999) and takes the
    # descending limb to one (p6 = 999)
    fish_sel_dbnrml_raw = matrix(c(1, 0), n_fish, 2, byrow = TRUE),
    fish_q_spec = rep("fix", n_fish)
  )

  ## Survey selectivity and catchability --------------------------------------
  # the same double normal with male offsets. the Western-Central survey
  # estimates its catchability under the assessment's normal prior on the log
  # scale, and the Eastern survey MIRRORS it, which the seeding expresses by
  # giving both cells of ln_srv_q one map level
  q_prior <- data.frame(region = 1, fleet = 1, block = 1, mu = exp(dat$q$prior_mean), sd = dat$q$prior_sd)
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = paste0("dbnrml_Fleet_", seq_len(n_srv)),
    cont_tv_srv_sel = paste0("none_Fleet_", seq_len(n_srv)),
    srv_sel_blocks = paste0("none_Fleet_", seq_len(n_srv)),
    srv_q_blocks = paste0("none_Fleet_", seq_len(n_srv)),
    srv_fixed_sel_pars_spec = rep("est_all", n_srv),
    srv_sel_sex_offset = rep("par", n_srv),
    srv_sel_dbnrml_raw = matrix(c(1, 0), n_srv, 2, byrow = TRUE),
    srv_q_spec = rep("est_all", n_srv),
    Use_srv_q_prior = 1, srv_q_prior = q_prior,
    t_srv = t_srv
  )

  ## Weighting ----------------------------------------------------------------
  # Francis weights, one per fleet, for lengths and ages. the conditional
  # age-at-length takes the age weight and needs the length dimension as an extra
  # axis, which is what the `extra` argument builds
  wl_f <- dat$var_adj_len[dat$fish_fleets]; wa_f <- dat$var_adj_age[dat$fish_fleets]
  wl_s <- dat$var_adj_len[dat$srv_fleets]; wa_s <- dat$var_adj_age[dat$srv_fleets]
  per_fleet <- function(w, n_fl, extra = NULL) {
    d <- c(n_reg, n_yrs, 1, if(!is.null(extra)) extra, n_sex, n_fl)
    arr <- array(1, dim = d)
    for(f in seq_len(n_fl)) if(is.null(extra)) arr[, , , , f] <- w[f] else arr[, , , , , f] <- w[f]
    arr
  }
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
    Wt_FishAgeComps = per_fleet(wa_f, n_fish), Wt_FishLenComps = per_fleet(wl_f, n_fish),
    Wt_SrvAgeComps = per_fleet(rep(1, n_srv), n_srv), Wt_SrvLenComps = per_fleet(wl_s, n_srv),
    Wt_Srv_caal = per_fleet(wa_s, n_srv, extra = n_lens)
  )

  input_list
} # end build_goa_rex_input


#' Set every parameter to the assessment's maximum likelihood estimate
#'
#' Seeding at the assessment's own estimate is what makes the bridge checkable:
#' every reported quantity and every likelihood component can be compared before
#' the optimizer is allowed to move anything. Each block also fixes the map,
#' since which parameters are estimated is part of the specification.
#'
#' @keywords internal
seed_goa_rex_mle <- function(input_list, dat) {

  yrs <- dat$years; n_yrs <- length(yrs)
  ages <- dat$ages; n_ages <- length(ages)
  n_reg <- dat$n_regions; n_sex <- dat$n_sexes
  n_fish <- dat$n_fish_fleets; n_srv <- dat$n_srv_fleets
  sigmaR <- dat$rec$sigmaR

  ## Recruitment level and apportionment --------------------------------------
  input_list$par$ln_global_R0[] <- dat$mle$ln_R0
  # the first area is the reference of the multinomial logit
  input_list$par$rec_region_prop_pars[1, ] <- dat$mle$rec_dist_area2

  ## Recruitment deviations ---------------------------------------------------
  # SPoRC's deviation is the assessment's less its bias correction. main
  # deviations run through 2022; the two later years carry none and are mapped
  # off. the series is shared across areas, so every region cell takes the same
  # map level and is penalized once
  main_yrs <- as.integer(names(dat$mle$main_recdev))
  dev_adj <- dat$mle$main_recdev - 0.5 * dat$mle$biasadj[as.character(main_yrs)] * sigmaR^2
  input_list$par$ln_RecDevs[] <- 0
  for(r in 1:n_reg) input_list$par$ln_RecDevs[1, r, match(main_yrs, yrs)] <- dev_adj
  map_rec <- array(NA_real_, dim = dim(input_list$par$ln_RecDevs))
  for(r in 1:n_reg) map_rec[1, r, match(main_yrs, yrs)] <- seq_along(main_yrs)
  input_list$map$ln_RecDevs <- factor(map_rec)

  ## Initial age structure ----------------------------------------------------
  # the deviation for year styr - a lands on age a, again less its own bias
  # correction, since the ramp is defined on calendar years. older ages carry
  # none
  early_yrs <- as.integer(names(dat$mle$early_recdev))
  early_adj <- dat$mle$early_recdev - 0.5 * dat$mle$biasadj[as.character(early_yrs)] * sigmaR^2
  input_list$par$ln_InitDevs[] <- 0
  map_init <- array(NA_real_, dim = dim(input_list$par$ln_InitDevs))
  for(i in seq_along(early_yrs)) {
    a <- yrs[1] - early_yrs[i] # the age this deviation sets in the start year
    input_list$par$ln_InitDevs[1, , a, ] <- early_adj[i]
    map_init[1, , a, ] <- i
  } # end i loop
  input_list$map$ln_InitDevs <- factor(map_init)

  ## Fishing mortality --------------------------------------------------------
  # a fixed mean plus deviations in the area with the fishery. the Eastern area
  # carries no fishery, so its mean sits at a rate that removes nothing
  lf <- log(dat$mle$Fmort[as.character(yrs)])
  input_list$par$ln_F_mean[] <- log(1e-12)
  input_list$par$ln_F_mean[1, 1, 1] <- mean(lf)
  input_list$par$ln_F_devs[] <- 0
  input_list$par$ln_F_devs[1, , 1, 1] <- lf - mean(lf)
  map_F <- array(NA_real_, dim = dim(input_list$par$ln_F_devs))
  map_F[1, , 1, 1] <- seq_len(n_yrs)
  input_list$map$ln_F_devs <- factor(map_F)

  ## Growth -------------------------------------------------------------------
  # every growth parameter is estimated, per area and sex
  input_list$par$ln_growth_pars[] <- log(dat$mle$growth)

  ## Selectivity --------------------------------------------------------------
  # female parameters straight in, male offsets in the second sex's slots: peak
  # in bins, ascending width on the log scale. the assessment gives no offset to
  # the plateau, the selectivity at the first bin, or the last bin here, so those
  # three slots stay at zero and unmapped
  put_sel <- function(par, map, tab, f, lev) {
    for(r in 1:n_reg) {
      par[r, , 1, 1, f] <- tab$female
      if(n_sex > 1) par[r, , 1, 2, f] <- c(tab$male[["Peak"]], 0, tab$male[["Ascend"]], tab$male[["Descend"]], 0, tab$male[["Final"]])
    }
    for(k in 1:6) if(tab$female_est[k]) { lev <- lev + 1; map[, k, 1, 1, f] <- lev }
    if(n_sex > 1) for(k in c(1, 3, 4, 6)) {
      nm <- c("Peak", NA, "Ascend", "Descend", NA, "Final")[k]
      if(isTRUE(tab$male_est[[nm]])) { lev <- lev + 1; map[, k, 1, 2, f] <- lev }
    }
    list(par = par, map = map, lev = lev)
  }
  map_fish <- array(NA_real_, dim = dim(input_list$par$fish_fixed_sel_pars)); lev <- 0
  for(f in seq_len(n_fish)) {
    out <- put_sel(input_list$par$fish_fixed_sel_pars, map_fish, dat$mle$sel[[dat$fish_fleets[f]]], f, lev)
    input_list$par$fish_fixed_sel_pars <- out$par; map_fish <- out$map; lev <- out$lev
  } # end f loop
  map_srv <- array(NA_real_, dim = dim(input_list$par$srv_fixed_sel_pars)); lev <- 0
  for(sf in seq_len(n_srv)) {
    out <- put_sel(input_list$par$srv_fixed_sel_pars, map_srv, dat$mle$sel[[dat$srv_fleets[sf]]], sf, lev)
    input_list$par$srv_fixed_sel_pars <- out$par; map_srv <- out$map; lev <- out$lev
  } # end sf loop
  input_list$map$fish_fixed_sel_pars <- factor(map_fish)
  input_list$map$srv_fixed_sel_pars <- factor(map_srv)

  ## Catchability -------------------------------------------------------------
  # the assessment mirrors the Eastern survey's q on the Western-Central
  # survey's, which the map expresses by giving the two cells one level
  input_list$par$ln_srv_q[] <- dat$mle$ln_q
  map_q <- array(NA_real_, dim = dim(input_list$par$ln_srv_q))
  map_q[1, 1, 1] <- 1 # Western-Central survey in area 1
  map_q[2, 1, 2] <- 1 # Eastern survey in area 2, the same parameter
  input_list$map$ln_srv_q <- factor(map_q)

  input_list
} # end seed_goa_rex_mle
