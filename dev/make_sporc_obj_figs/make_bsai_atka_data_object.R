# Purpose: Build sgl_rg_bsai_atka_data, the container for the BSAI Atka mackerel
#          case study. Everything the case study and its regression tests need
#          lives in the object: model inputs, the AMAK maximum likelihood
#          estimate used as a seed, and the AMAK output the bridge is compared
#          against.
# Creator: Matthew LH. Cheng
# Date Created: 8/18/26
#
# Sources, all under dev/dev_data/atka_m16.0b_2024:
#   am2024.dat   2024 assessment input file
#   amak.dat     control file, carries the selectivity weights and the phases
#   amak.par     AMAK parameter file, read at full precision
#   For_R.rep    AMAK report object, used to VALIDATE the reconstruction rather
#                than as the comparison target; it prints six significant digits
#   amak.tpl     the template, for reference
#
# Atka differs from the rockfish case studies in four ways the object has to
# carry: fishery selectivity is a non-parametric log scale surface with a
# separate coefficient vector in all but one year, the assessment fits a
# Beverton-Holt curve as a penalty on free recruitment rather than generating
# recruitment from it, the survey standardizes selectivity over the ages
# catchability is defined on rather than over all ages, and both selectivity
# shape weights arrive in the control file as standard deviations and are
# converted to precisions before the likelihood sees them.

library(here)
source(here("dev", "make_sporc_obj_figs", "helper-amak.R"))

src <- function(f) here("dev", "dev_data", "atka_m16.0b_2024", f)
dat <- read_amak_dat(src("am2024.dat"))
ctl <- read_amak_ctl(src("amak.dat"), dat)
apar <- read_amak_par(src("amak.par"))
arep <- read_amak_R_rep(src("For_R.rep"))

# The comparison target, rebuilt from amak.par at full precision -------------
dyn <- amak_dynamics(dat, ctl, apar)
aobj <- amak_objective(dat, ctl, apar, dyn)

# Guard: the reconstruction has to reproduce every reported component before
# anything downstream is entitled to trust it.
stopifnot(max(abs(as.numeric(aobj$comps) - arep$Like_Comp) /
                pmax(abs(arep$Like_Comp), 1e-8)) < 1e-4)
stopifnot(abs(aobj$comps[["total"]] - apar$objective) < 1e-6)

yrs <- dat$years
n_yrs <- dat$nyrs
n_ages <- dat$nages

i_fsh_age <- match(dat$yrs_fsh_age[1, ], yrs)
i_srv <- match(dat$yrs_ind[1, ], yrs)
i_srv_age <- match(dat$yrs_ind_age[1, ], yrs)

# Model inputs ---------------------------------------------------------------
ObsCatch <- array(dat$catch_bio[1, ], dim = c(1, n_yrs, 1, 1))
UseCatch <- array(1, dim = c(1, n_yrs, 1, 1))

ObsFishAgeComps <- array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1))
ObsFishAgeComps[1, i_fsh_age, 1, , 1, 1] <- dat$oac_fsh
UseFishAgeComps <- array(0, dim = c(1, n_yrs, 1, 1))
UseFishAgeComps[1, i_fsh_age, 1, 1] <- 1
ISS_FishAgeComps <- array(0, dim = c(1, n_yrs, 1, 1, 1))
ISS_FishAgeComps[1, i_fsh_age, 1, 1, 1] <- dat$n_sample_fsh_age[1, ]

# AMAK converts the arithmetic scale survey standard error to a log scale one
# before the likelihood ever sees it, so the converted value is what belongs in
# the object.
ObsSrvIdx <- array(NA_real_, dim = c(1, n_yrs, 1, 1))
ObsSrvIdx_SE <- array(NA_real_, dim = c(1, n_yrs, 1, 1))
UseSrvIdx <- array(0, dim = c(1, n_yrs, 1, 1))
ObsSrvIdx[1, i_srv, 1, 1] <- dat$obs_ind[1, ]
ObsSrvIdx_SE[1, i_srv, 1, 1] <- dat$obs_lse_ind[1, ]
UseSrvIdx[1, i_srv, 1, 1] <- 1

ObsSrvAgeComps <- array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1))
ObsSrvAgeComps[1, i_srv_age, 1, , 1, 1] <- dat$oac_ind
UseSrvAgeComps <- array(0, dim = c(1, n_yrs, 1, 1))
UseSrvAgeComps[1, i_srv_age, 1, 1] <- 1
ISS_SrvAgeComps <- array(0, dim = c(1, n_yrs, 1, 1, 1))
ISS_SrvAgeComps[1, i_srv_age, 1, 1, 1] <- dat$n_sample_ind_age[1, ]

# Population weight at age is a single vector in AMAK; the fishery and the
# survey each carry an annual matrix. The maturity ogive goes in RAW: amak.tpl
# halves an ogive that reaches one, and SPoRC halves spawning biomass itself in
# a single sex model, so passing the halved vector would halve it twice.
WAA <- array(rep(dat$wt_pop, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
WAA_fish <- array(dat$wt_fsh, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
WAA_srv <- array(dat$wt_ind, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
MatAA <- array(rep(dat$maturity_in, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
Fixed_natmort <- array(apar$Mest, dim = c(1, 1, n_yrs, n_ages, 1))

# AMAK forms expected compositions as `age_err %*% p`, which puts the observed
# age axis on the ROWS. SPoRC forms them as `p %*% AgeingError`, so the matrix
# transposes on the way in.
AgeingError <- t(dat$age_err)

# Selectivity configuration --------------------------------------------------
# n_sel_ch is incremented and styr prepended on read, so "46 changes" is 47
# blocks, and the terminal year shares the last one. nselages edge holds the
# oldest estimated coefficient to the plus group BEFORE the standardization,
# which is a bin grouping rather than a separate parameter.
n_blk_fsh <- ctl$n_sel_ch_fsh
blk_yr_fsh <- match(ctl$yrs_sel_ch_fsh, yrs)
blk_yr_srv <- match(ctl$yrs_sel_ch_ind, yrs)

fish_sel_blocks <- c(
  paste0("Block_", seq_len(n_blk_fsh - 1), "_Year_",
         blk_yr_fsh[-n_blk_fsh], "-", blk_yr_fsh[-1] - 1, "_Fleet_1"),
  paste0("Block_", n_blk_fsh, "_Year_", blk_yr_fsh[n_blk_fsh], "-terminal_Fleet_1")
)

# Every shape weight is a precision by the time the likelihood sees it. The
# curvature sigma is converted at amak.tpl:948, the dome sigma on read at
# amak.tpl:610, and the random walk sigma inside Sel_Like. Only the change years
# carry a weight; the years in between inherit their block's coefficients and
# are not penalized again.
zero_but <- function(idx, value) { w <- rep(0, n_yrs); w[idx] <- value; w }
fish_sel_pen_wts <- list(list(
  smooth_bin_curve = zero_but(blk_yr_fsh, ctl$curv_pen_fsh),
  smooth_bin_diff = 0,
  smooth_yr_diff = zero_but(blk_yr_fsh[-1], 0.5 / ctl$sel_sigma_fsh[-1]^2),
  smooth_yr_curve = 0,
  smooth_dome = zero_but(blk_yr_fsh, 0.5 / ctl$seldec_pen_fsh),
  smooth_mean_center = 0,
  normalize = FALSE,
  bin_range = list(smooth_dome = c(ctl$seldecage - 1, ctl$nselages_fsh)),
  yr_diff_ref = NULL
))
srv_sel_pen_wts <- list(list(
  smooth_bin_curve = zero_but(blk_yr_srv, ctl$curv_pen_ind),
  smooth_bin_diff = 0, smooth_yr_diff = 0, smooth_yr_curve = 0,
  smooth_dome = zero_but(blk_yr_srv, 0.5 / ctl$seldec_pen_ind),
  smooth_mean_center = 0,
  normalize = FALSE,
  bin_range = list(smooth_dome = c(ctl$seldecage - 1, ctl$nselages_ind)),
  yr_diff_ref = NULL
))

# obj_fun += 20 * log(mean(exp(coefficients)))^2 for every block, which AMAK
# folds straight into the objective rather than into its reported selectivity
# component.
fish_selex_penalty <- data.frame(region = 1, fleet = 1, block = seq_len(n_blk_fsh),
                                 sex = 1, wt = 20)
fish_selex_penalty$par <- rep(list(seq_len(ctl$nselages_fsh)), n_blk_fsh)
srv_selex_penalty <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 20)
srv_selex_penalty$par <- list(seq_len(ctl$nselages_ind))

# The survey standardizes over the ages catchability is defined on while SPoRC
# standardizes over all of them. That is a constant rescaling of selectivity and
# the reciprocal rescaling of catchability, so the PRIOR MEAN has to move with
# it or the same prior statement lands on a different number.
p_ind <- c(apar$log_selcoffs_ind_1, apar$log_selcoffs_ind_1[ctl$nselages_ind])
q_sel_rescale <- mean(exp(p_ind[ctl$q_age_min_idx:ctl$q_age_max_idx])) / mean(exp(p_ind))

sgl_rg_bsai_atka_data <- list(
  # Dimensions
  years = yrs, ages = dat$ages, lens = NA,
  n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
  n_seas = 1, n_pop = 1, natal_region = 1,

  # Observation years, for indexing figures and fits
  yrs_fsh_age = dat$yrs_fsh_age[1, ], yrs_srv = dat$yrs_ind[1, ],
  yrs_srv_age = dat$yrs_ind_age[1, ],

  # Data
  ObsCatch = ObsCatch, UseCatch = UseCatch,
  ObsFishAgeComps = ObsFishAgeComps, UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,
  ObsSrvIdx = ObsSrvIdx, ObsSrvIdx_SE = ObsSrvIdx_SE, UseSrvIdx = UseSrvIdx,
  ObsSrvAgeComps = ObsSrvAgeComps, UseSrvAgeComps = UseSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,

  # Biologicals
  WAA = WAA, WAA_fish = WAA_fish, WAA_srv = WAA_srv, MatAA = MatAA,
  AgeingError = AgeingError, Fixed_natmort = Fixed_natmort,
  t_spawn = dat$spmo_frac, t_srv = dat$ind_month_frac[1],
  sigmaC = dat$catch_bio_lsd[1, 1],
  offset_fsh = dat$offset_fsh, offset_ind = dat$offset_ind,

  # Selectivity and catchability configuration
  n_blk_fsh = n_blk_fsh, blk_yr_fsh = blk_yr_fsh,
  fish_sel_blocks = fish_sel_blocks,
  nselages_fsh = ctl$nselages_fsh, nselages_srv = ctl$nselages_ind,
  fish_sel_pen_wts = fish_sel_pen_wts, srv_sel_pen_wts = srv_sel_pen_wts,
  fish_selex_penalty = fish_selex_penalty, srv_selex_penalty = srv_selex_penalty,
  q_sel_rescale = q_sel_rescale,
  # The ages the index selectivity standardization averages over (amak.tpl:1947)
  q_age_min = ctl$q_age_min_idx, q_age_max = ctl$q_age_max_idx,
  srv_q_prior = data.frame(region = 1, fleet = 1, block = 1,
                           mu = ctl$qprior, sd = ctl$cvqprior),

  # Recruitment configuration
  steepness = apar$steepness, sigmaR = exp(apar$log_sigmar),
  styr_rec_est = ctl$styr_rec_est, endyr_rec_est = ctl$endyr_rec_est,
  nrecs_est = ctl$nrecs_est, styr_rec = dat$styr_rec,

  # AMAK maximum likelihood estimate, the seed
  mle = list(
    mean_log_rec = apar$mean_log_rec, log_Rzero = apar$log_Rzero,
    log_sigmar = apar$log_sigmar, Mest = apar$Mest,
    rec_dev = apar$rec_dev, fmort = apar$fmort,
    log_selcoffs_fsh = apar$log_selcoffs_fsh_1,
    log_selcoffs_ind = apar$log_selcoffs_ind_1,
    log_q_ind = apar$log_q_ind_1,
    objective = apar$objective, max_grad = apar$max_grad
  ),

  # AMAK output at full precision, the comparison target
  amak = list(
    NAA = dyn$natage[1:n_yrs, ], SSB = as.numeric(dyn$Sp_Biom[as.character(yrs)]),
    Rec = dyn$natage[1:n_yrs, 1], FAA = dyn$F, ZAA = dyn$Z, CAA = dyn$catage,
    sel_fsh = dyn$sel_fsh, sel_ind = dyn$sel_ind,
    pred_catch = dyn$pred_catch, pred_ind = dyn$pred_ind,
    eac_fsh = dyn$eac_fsh, eac_ind = dyn$eac_ind,
    rec_dev = dyn$rec_dev, chi = aobj$chi, Bzero = dyn$Bzero, phizero = dyn$phizero,
    Like_Comp = aobj$comps, sel_fsh_parts = aobj$sel_fsh_parts,
    sel_ind_parts = aobj$sel_ind_parts, rec_parts = aobj$rec_parts,
    SSQRec = aobj$SSQRec, avgsel_fsh = aobj$avgsel_fsh, avgsel_ind = aobj$avgsel_ind,
    Like_Comp_reported = arep$Like_Comp
  )
)

usethis::use_data(sgl_rg_bsai_atka_data, internal = FALSE, overwrite = TRUE)
