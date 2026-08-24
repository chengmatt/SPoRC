# Purpose: Bridge the 2025 GOA rex sole assessment (Model 25.1) to SPoRC from the
#          packaged data object, and render the case study figures. The model is
#          specified the way the assessment specifies it: one population
#          recruiting into two areas by an estimated apportionment with no
#          movement between them, mean recruitment with the bias ramp and 17
#          early deviations setting the start year's ages 1-17, growth estimated
#          per area and sex in Schnute's form with weight at age derived from it,
#          age based double normal selectivity with male offsets, survey
#          catchability mirrored between the two surveys, and the survey's
#          conditional age-at-length carrying the age information.
#
#          Three stages: set every parameter to the assessment's maximum
#          likelihood estimate and check the population and the growth there,
#          optimize, then compare.
#
#          The specification below is written out in full rather than sourced,
#          so the script reads on its own. The regression test builds the same
#          model from tests/testthat/helper-bridge_goa_rex.R, and the check at
#          the foot of this script compares the two so they cannot drift apart
#          silently.
# Creator: Matthew LH. Cheng
# Date Created: 8/22/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

label <- "2025 GOA Rex Sole Assessment"
area_names <- c("Western-Central GOA", "Eastern GOA")
sex_names <- c("female", "male")
fleet_labels <- c("Fishery", "Western-Central survey", "Eastern survey")

dat <- SPoRC::mlt_rg_goa_rex_data
yrs <- dat$years
n_yrs <- length(yrs)
ages <- dat$ages
n_ages <- length(ages)
n_reg <- dat$n_regions
n_sex <- dat$n_sexes
n_fish <- dat$n_fish_fleets
n_srv <- dat$n_srv_fleets
n_lens <- length(dat$lens)
mle <- dat$mle

# composition type strings, one per fleet
joint_type <- function(n) paste0("spltRjntS_Year_1-terminal_Fleet_", seq_len(n))
split_type <- function(n) paste0("spltRspltS_Year_1-terminal_Fleet_", seq_len(n))
no_type <- function(n) paste0("none_Year_1-terminal_Fleet_", seq_len(n))

# Setup ------------------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  years = yrs, ages = ages, lens = dat$lens,
  n_regions = n_reg, n_sexes = n_sex,
  n_fish_fleets = n_fish, n_srv_fleets = n_srv,
  n_seas = 1, n_pop = 1, natal_region = 1, verbose = FALSE
)

# Recruitment: age 0 at the start of the year, mean recruitment apportioned across
# the two areas by an estimated logit, the bias ramp on the deviations, and the
# start year's ages 1-17 set by early deviations shared across areas and sexes.
# The ramp's breakpoints are given in deviation-index space, where 1 is the first
# model year, so the two that precede the model start come out negative.
input_list <- Setup_Mod_Rec(
  input_list = input_list,
  rec_model = "mean_rec", rec_dd = "global", rec_lag = 0, SR_ref_yr = 1, t_spawn = 0,
  sigmaR_spec = "fix", sigmaR_switch = 1,
  ln_sigmaR = array(log(dat$rec$sigmaR), dim = c(2, 1, n_reg)),
  do_rec_bias_ramp = 1,
  bias_year = dat$rec$bias_years - yrs[1] + 1,
  max_bias_ramp_fct = dat$rec$max_bias_adj,
  RecDevs_spec = "est_shared_pop_r", RecDevs_pen_center = "fixed", dont_est_recdev_last = 0,
  init_age_strc = 2, equil_init_age_strc = 1,
  InitDevs_spec = "est_shared_pop_r", InitDevs_sex_spec = "est_shared_s",
  InitDevs_pen_center = "fixed",
  rec_region_prop_spec = NULL,
  ln_global_R0 = mle$ln_R0
)

# Maturity: logistic on age for females, with none below the first mature age.
# Males carry none, since spawning biomass is female only.
MatAA <- array(0, dim = c(1, n_reg, n_yrs, 1, n_ages, n_sex))
mat_f <- 1 / (1 + exp(dat$mat$slope * (ages - dat$mat$a50)))
mat_f[ages < dat$mat$first_mature_age] <- 0
for(r in 1:n_reg) for(y in 1:n_yrs) MatAA[1, r, y, 1, , 1] <- mat_f

# Weight-length parameters, one pair per area and sex
wt_len <- array(NA_real_, dim = c(1, n_reg, n_sex, 2))
for(r in 1:n_reg) {
  wt_len[1, r, 1, ] <- dat$wtlen$fem
  if(n_sex > 1) wt_len[1, r, 2, ] <- dat$wtlen$mal
} # end r loop

# Biology: natural mortality fixed, growth estimated per area and sex with weight
# at age derived from the key rather than read as data, the assessment's ageing
# error matrix, and its composition constant on both sides of the multinomial.
input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA = NULL, MatAA = MatAA,
  fit_lengths = 1, SizeAgeTrans = NA,
  AgeingError = dat$AgeingError,
  M_spec = "fix",
  Fixed_natmort = array(dat$growth[[1]]$M, dim = c(1, n_reg, n_yrs, n_ages, n_sex)),
  addtocomp = dat$comp$addtocomp_age, comp_const_obs = 1,
  addtosrvidx = 0, addtofishidx = 0,
  growth_model = "vb_schnute", growth_pars = mle$growth, growth_spec = "est_all",
  growth_A1 = dat$growth_A1, growth_A2 = dat$growth_A2,
  growth_len_lower = dat$lens_lower, growth_L0 = dat$lens_lower[1],
  growth_plus_group = "mixture",
  waa_model = "wt_len", wt_len_pars = wt_len
)

# Nothing moves between the areas, and there are no tags
input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                 Fixed_Movement = NA, do_recruits_move = 0)

input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

# The assessment solves fishing mortality from the catch under its hybrid method,
# so here it is a free parameter per year with no penalty, held to the catch by a
# tight lognormal error
input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,
  ObsCatch = dat$ObsCatch, UseCatch = dat$UseCatch,
  Use_F_pen = 0, ln_F_mean_spec = "fix",
  sigmaC_spec = "fix",
  ln_sigmaC = array(log(dat$catch_se_value), dim = c(n_reg, n_yrs, 1, n_fish))
)

# Fishery: no index, joint-sex marginal ages and lengths, no conditional
# age-at-length. Its key and weight at age are read at the mid season, which is
# where the assessment reads the catch's.
input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  t_fish = array(0.5, dim = c(n_reg, 1, n_fish)),
  ObsFishIdx = array(NA_real_, dim = c(n_reg, n_yrs, 1, n_fish)),
  ObsFishIdx_SE = array(NA_real_, dim = c(n_reg, n_yrs, 1, n_fish)),
  UseFishIdx = array(0, dim = c(n_reg, n_yrs, 1, n_fish)),
  fish_idx_type = rep("none", n_fish),
  FishIdx_LikeType = rep("lognormal", n_fish),
  ObsFishAgeComps = dat$ObsFishAgeComps, UseFishAgeComps = dat$UseFishAgeComps,
  ISS_FishAgeComps = dat$ISS_FishAgeComps,
  ObsFishLenComps = dat$ObsFishLenComps, UseFishLenComps = dat$UseFishLenComps,
  ISS_FishLenComps = dat$ISS_FishLenComps,
  FishAgeComps_LikeType = rep("Multinomial", n_fish),
  FishLenComps_LikeType = rep("Multinomial", n_fish),
  FishAgeComps_Type = joint_type(n_fish),
  FishLenComps_Type = joint_type(n_fish)
)

# Surveys: a biomass index in each area, joint-sex lengths, and conditional
# age-at-length split by sex, since each row holds one sex's otoliths. The
# assessment carries the surveys' marginal ages as ghost fleets, so they go in
# but are not fit.
t_srv <- array(rep(dat$t_srv, each = n_reg), dim = c(n_reg, 1, n_srv))
input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,
  t_srv = t_srv,
  ObsSrvIdx = dat$ObsSrvIdx, ObsSrvIdx_SE = dat$ObsSrvIdx_SE, UseSrvIdx = dat$UseSrvIdx,
  srv_idx_type = rep("biom", n_srv),
  SrvIdx_LikeType = rep("lognormal", n_srv),
  ObsSrvAgeComps = dat$ObsSrvAgeComps,
  UseSrvAgeComps = array(0, dim = dim(dat$UseSrvAgeComps)),
  ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
  ObsSrvLenComps = dat$ObsSrvLenComps, UseSrvLenComps = dat$UseSrvLenComps,
  ISS_SrvLenComps = dat$ISS_SrvLenComps,
  SrvAgeComps_LikeType = rep("none", n_srv),
  SrvLenComps_LikeType = rep("Multinomial", n_srv),
  SrvAgeComps_Type = no_type(n_srv),
  SrvLenComps_Type = joint_type(n_srv),
  ObsSrv_caal = dat$ObsSrv_caal, UseSrv_caal = dat$UseSrv_caal,
  ISS_Srv_caal = dat$ISS_Srv_caal,
  Srv_caal_LikeType = rep("Multinomial", n_srv),
  Srv_caal_Type = split_type(n_srv)
)

# Selectivity is the age based double normal everywhere, time invariant, with the
# male curve a parameter offset from the female one. The assessment leaves the
# ascending limb's starting height unanchored and takes the descending limb to
# one, which the raw specification carries.
input_list <- Setup_Mod_Fishsel_and_Q(
  input_list = input_list,
  fish_sel_model = paste0("dbnrml_Fleet_", seq_len(n_fish)),
  cont_tv_fish_sel = paste0("none_Fleet_", seq_len(n_fish)),
  fish_sel_blocks = paste0("none_Fleet_", seq_len(n_fish)),
  fish_q_blocks = paste0("none_Fleet_", seq_len(n_fish)),
  fish_fixed_sel_pars_spec = rep("est_all", n_fish),
  fish_sel_sex_offset = rep("par", n_fish),
  fish_sel_dbnrml_raw = matrix(c(1, 0), n_fish, 2, byrow = TRUE),
  fish_q_spec = rep("fix", n_fish)
)

# Catchability: the Western-Central survey estimates its own under the
# assessment's normal prior on the log scale, and the Eastern survey mirrors it,
# which the seeding below expresses through the map.
q_prior <- data.frame(region = 1, fleet = 1, block = 1,
                      mu = exp(dat$q$prior_mean), sd = dat$q$prior_sd)

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  t_srv = t_srv,
  srv_sel_model = paste0("dbnrml_Fleet_", seq_len(n_srv)),
  cont_tv_srv_sel = paste0("none_Fleet_", seq_len(n_srv)),
  srv_sel_blocks = paste0("none_Fleet_", seq_len(n_srv)),
  srv_q_blocks = paste0("none_Fleet_", seq_len(n_srv)),
  srv_fixed_sel_pars_spec = rep("est_all", n_srv),
  srv_sel_sex_offset = rep("par", n_srv),
  srv_sel_dbnrml_raw = matrix(c(1, 0), n_srv, 2, byrow = TRUE),
  srv_q_spec = rep("est_all", n_srv),
  Use_srv_q_prior = 1, srv_q_prior = q_prior
)

# Francis weights, one per fleet on the lengths and one on the ages, with the
# conditional age-at-length taking each survey's age weight
per_fleet_wt <- function(w, n_fleets, extra = NULL) {
  d <- c(n_reg, n_yrs, 1, if(!is.null(extra)) extra, n_sex, n_fleets)
  arr <- array(1, dim = d)
  for(f in seq_len(n_fleets)) {
    if(is.null(extra)) arr[, , , , f] <- w[f] else arr[, , , , , f] <- w[f]
  } # end f loop
  arr
}
wt_len_fish <- dat$var_adj_len[dat$fish_fleets]
wt_age_fish <- dat$var_adj_age[dat$fish_fleets]
wt_len_srv <- dat$var_adj_len[dat$srv_fleets]
wt_age_srv <- dat$var_adj_age[dat$srv_fleets]

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
  Wt_FishAgeComps = per_fleet_wt(wt_age_fish, n_fish),
  Wt_FishLenComps = per_fleet_wt(wt_len_fish, n_fish),
  Wt_SrvAgeComps = per_fleet_wt(rep(1, n_srv), n_srv),
  Wt_SrvLenComps = per_fleet_wt(wt_len_srv, n_srv),
  Wt_Srv_caal = per_fleet_wt(wt_age_srv, n_srv, extra = n_lens)
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Stage 1: start at the assessment's estimate ----------------------------------
parameters$ln_global_R0[] <- mle$ln_R0
# regional apportionment, with the first area the reference of the logit
parameters$rec_region_prop_pars[1, ] <- mle$rec_dist_area2

# SPoRC's deviation is the assessment's less its bias correction. The main
# deviations run through 2022; the two later years carry none and are mapped off.
main_yrs <- as.integer(names(mle$main_recdev))
main_dev <- mle$main_recdev - 0.5 * mle$biasadj[as.character(main_yrs)] * dat$rec$sigmaR^2
parameters$ln_RecDevs[] <- 0
map_rec <- array(NA_real_, dim = dim(parameters$ln_RecDevs))
for(r in 1:n_reg) {
  parameters$ln_RecDevs[1, r, match(main_yrs, yrs)] <- main_dev
  map_rec[1, r, match(main_yrs, yrs)] <- seq_along(main_yrs)
} # end r loop
mapping$ln_RecDevs <- factor(map_rec)

# The early deviations set the start year's ages 1-17: the deviation belonging to
# year styr - a lands on age a, again less its bias correction. Older ages carry
# none, and the deviation is shared across areas and sexes.
early_yrs <- as.integer(names(mle$early_recdev))
early_dev <- mle$early_recdev - 0.5 * mle$biasadj[as.character(early_yrs)] * dat$rec$sigmaR^2
parameters$ln_InitDevs[] <- 0
map_init <- array(NA_real_, dim = dim(parameters$ln_InitDevs))
for(i in seq_along(early_yrs)) {
  a <- yrs[1] - early_yrs[i] # the age this deviation sets in the start year
  parameters$ln_InitDevs[1, , a, ] <- early_dev[i]
  map_init[1, , a, ] <- i
} # end i loop
mapping$ln_InitDevs <- factor(map_init)

# Fishing mortality as a fixed mean plus deviations in the area that has the
# fishery. The other area carries none, so its mean sits at a rate that removes
# nothing and its deviations stay mapped off.
log_f <- log(mle$Fmort[as.character(yrs)])
parameters$ln_F_mean[] <- log(1e-12)
parameters$ln_F_mean[1, 1, 1] <- mean(log_f)
parameters$ln_F_devs[] <- 0
parameters$ln_F_devs[1, , 1, 1] <- log_f - mean(log_f)
map_F <- array(NA_real_, dim = dim(parameters$ln_F_devs))
map_F[1, , 1, 1] <- seq_len(n_yrs)
mapping$ln_F_devs <- factor(map_F)

# Growth at the assessment's values, every parameter estimated
parameters$ln_growth_pars[] <- log(mle$growth)

# Selectivity: the female parameters go straight into the first sex's slots and
# the male offsets into the second's. The assessment offsets only the peak, the
# two widths and the selectivity at the last bin, so the plateau and the
# first-bin parameter carry no male offset. Which parameters it estimates is
# carried in the same table, and the map gives each estimated one its own level.
sel_par_names <- c("Peak", "Plateau", "Ascend", "Descend", "First", "Final")
male_offset_pars <- c(1, 3, 4, 6)

seed_selex <- function(par, map, tab, f, level) {
  male <- rep(0, length(sel_par_names))
  if(n_sex > 1) male[male_offset_pars] <- unlist(tab$male[sel_par_names[male_offset_pars]])
  for(r in 1:n_reg) {
    par[r, , 1, 1, f] <- tab$female
    if(n_sex > 1) par[r, , 1, 2, f] <- male
  } # end r loop
  for(k in seq_along(sel_par_names)) {
    if(isTRUE(tab$female_est[k])) { level <- level + 1; map[, k, 1, 1, f] <- level }
  } # end k loop
  if(n_sex > 1) for(k in male_offset_pars) {
    if(isTRUE(tab$male_est[[sel_par_names[k]]])) { level <- level + 1; map[, k, 1, 2, f] <- level }
  } # end k loop
  list(par = par, map = map, level = level)
}

map_fish_sel <- array(NA_real_, dim = dim(parameters$fish_fixed_sel_pars))
level <- 0
for(f in seq_len(n_fish)) {
  out <- seed_selex(parameters$fish_fixed_sel_pars, map_fish_sel, mle$sel[[dat$fish_fleets[f]]], f, level)
  parameters$fish_fixed_sel_pars <- out$par; map_fish_sel <- out$map; level <- out$level
} # end f loop
mapping$fish_fixed_sel_pars <- factor(map_fish_sel)

map_srv_sel <- array(NA_real_, dim = dim(parameters$srv_fixed_sel_pars))
level <- 0
for(sf in seq_len(n_srv)) {
  out <- seed_selex(parameters$srv_fixed_sel_pars, map_srv_sel, mle$sel[[dat$srv_fleets[sf]]], sf, level)
  parameters$srv_fixed_sel_pars <- out$par; map_srv_sel <- out$map; level <- out$level
} # end sf loop
mapping$srv_fixed_sel_pars <- factor(map_srv_sel)

# Catchability: the assessment mirrors the Eastern survey's q on the
# Western-Central survey's, which the map expresses by giving the two cells one
# level so they hold one parameter.
parameters$ln_srv_q[] <- mle$ln_q
map_q <- array(NA_real_, dim = dim(parameters$ln_srv_q))
map_q[1, 1, 1] <- 1 # Western-Central survey, area 1
map_q[2, 1, 2] <- 1 # Eastern survey, area 2, the same parameter
mapping$ln_srv_q <- factor(map_q)

seed <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
r <- seed$rep

# Stage 1 comparison. The assessment's report carries six significant digits, so
# agreement is judged against that rather than machine precision.
pd_max <- function(a, b) max(abs(100 * (a - b) / b), na.rm = TRUE)

naa_pct <- ssb_pct <- rec_pct <- totb_pct <- numeric(n_reg)
for(a in 1:n_reg) {
  naa_pct[a] <- pd_max(r$NAA[1, a, 1:n_yrs, 1, , ], dat$ss3$NAA[a, , , ])
  ssb_pct[a] <- pd_max(r$SSB[1, a, 1:n_yrs], dat$ss3$SSB[, a])
  rec_pct[a] <- pd_max(r$Rec[1, a, 1:n_yrs], dat$ss3$Rec[, a])
  totb_jan1 <- sapply(1:n_yrs, function(y) sum(r$NAA[1, a, y, 1, , ] * r$WAA[1, a, y, 1, , ]))
  totb_pct[a] <- pd_max(totb_jan1, dat$ss3$Bio_all[, a])
} # end a loop

catch_pct <- pd_max(r$PredCatch[1, 1, , 1, 1], dat$ss3$dead_B)

srv_idx_pct <- numeric(n_srv)
for(sf in 1:n_srv) {
  cpue <- dat$ss3$cpue[dat$ss3$cpue$Fleet == sf + 1, ]
  srv_idx_pct[sf] <- pd_max(r$PredSrvIdx[1, sf, match(cpue$Yr, yrs), 1, sf], cpue$Exp)
} # end sf loop

# Growth, read at the survey's timing, which is the mid season here
laa_pct <- sd_pct <- waa_pct <- array(NA_real_, dim = c(n_reg, n_sex))
for(a in 1:n_reg) for(s in 1:n_sex) {
  growth_ss3 <- dat$ss3$growth[[a]][[s]]
  laa_pct[a, s] <- pd_max(r$mean_LAA_srv[1, a, 1, 1, , s, 1], growth_ss3$Len_Mid)
  sd_pct[a, s] <- pd_max(r$sd_LAA_srv[1, a, 1, 1, , s, 1], growth_ss3$SD_Mid)
  waa_pct[a, s] <- pd_max(r$WAA_fish[1, a, 1, 1, , s, 1], growth_ss3$Wt_Mid)
} # end a and s loops

cat("=== Stage 1: at the assessment's estimate (its report carries six significant digits) ===\n")
stage1 <- data.frame(
  quantity = c("numbers at age", "spawning biomass", "recruitment", "total biomass (Jan 1)",
               "predicted catch", "survey indices", "mean length at age",
               "SD of length at age", "weight at age"),
  max_pct = c(max(naa_pct), max(ssb_pct), max(rec_pct), max(totb_pct),
              catch_pct, max(srv_idx_pct), max(laa_pct), max(sd_pct), max(waa_pct))
)
print(stage1, row.names = FALSE, digits = 3)

# Selectivity against the assessment's own curves, every fleet and sex
sel_gap <- function(rep) {
  gaps <- numeric(0)
  for(s in 1:n_sex) gaps <- c(gaps, max(abs(rep$fish_sel[1, 1, 1, 1, , s, 1] - dat$ss3$sel[[1]][[s]])))
  for(sf in 1:n_srv) for(s in 1:n_sex) {
    gaps <- c(gaps, max(abs(rep$srv_sel[1, sf, 1, 1, , s, sf] - dat$ss3$sel[[sf + 1]][[s]])))
  } # end sf and s loops
  max(gaps)
}
cat("selectivity against the assessment's own curves, max absolute difference:",
    signif(sel_gap(r), 3), "\n")

# The age-length key the survey compositions read, against the assessment's own
alk_ss3 <- dat$ss3$ALK[rev(seq_len(dim(dat$ss3$ALK)[1])), , "Seas: 1 Sub_Seas: 2 Morph: 1"]
alk_sporc <- r$SizeAgeTrans_srv[1, 1, 1, 1, , , 1, 1]
cat("mid-season age-length key, max absolute difference:",
    signif(max(abs(alk_sporc - alk_ss3)), 3), "\n")

# Stage 2: optimize ------------------------------------------------------------
est <- fit_model(data, parameters, mapping, do_optim = TRUE, newton_loops = 3, silent = TRUE)
est$sdrep <- tryCatch(RTMB::sdreport(est, hessian.fixed = est$he(est$optim$par)),
                      error = function(e) NULL)
pl <- est$env$parList(est$optim$par)

cat(sprintf("\n=== Stage 2: %d parameters, objective %.6f, fallen %.4f from the assessment's estimate, max |gradient| %.2e ===\n",
            length(est$par), est$optim$objective,
            seed$fn(seed$par) - est$optim$objective, max(abs(est$gr(est$optim$par)))))
cat("selectivity after fitting, max absolute difference from the assessment's:",
    signif(sel_gap(est$rep), 3), "\n")

# Stage 3: compare -------------------------------------------------------------
ssb_area <- est$rep$SSB[1, , 1:n_yrs] # [area, year]
rec_area <- est$rep$Rec[1, , 1:n_yrs]
ssb_tot <- colSums(ssb_area)
rec_tot <- colSums(rec_area)
ss3_ssb_tot <- rowSums(dat$ss3$SSB)
ss3_rec_tot <- rowSums(dat$ss3$Rec)

# Standard errors on the natural scale. The log quantities are reported in array
# order (population, area, year), and the leading population dimension is one, so
# the entries lay out area fastest and then year. Total spawning biomass has its
# own reported entry; total recruitment is formed from the covariance of the two
# areas' log recruitments, since the areas are correlated.
se_from <- function(nm, n_lead) {
  s <- est$sdrep$sd[names(est$sdrep$value) == nm]
  matrix(s, nrow = n_lead)[, 1:n_yrs, drop = FALSE]
}
ssb_se_area <- rec_se_area <- array(NA_real_, dim = c(n_reg, n_yrs))
ssb_se_tot <- rec_se_tot <- rep(NA_real_, n_yrs)

if(!is.null(est$sdrep)) {
  ssb_se_area <- se_from("log_SSB", n_reg) * ssb_area
  rec_se_area <- se_from("log_Rec", n_reg) * rec_area
  ssb_se_tot <- as.vector(se_from("log_Aggregated_SSB", 1)) * ssb_tot
  if(!is.null(est$sdrep$cov)) {
    i_rec <- which(names(est$sdrep$value) == "log_Rec")
    V <- est$sdrep$cov[i_rec, i_rec]
    entry <- function(a, y) 1 + (a - 1) + n_reg * (y - 1) # position of (1, a, y) in log_Rec
    for(y in 1:n_yrs) {
      ii <- sapply(1:n_reg, entry, y = y)
      R <- rec_area[, y]
      rec_se_tot[y] <- sqrt(as.numeric(t(R) %*% V[ii, ii] %*% R))
    } # end y loop
  }
}

cat("\n=== Stage 3: optimized SPoRC against the assessment ===\n")
cmp_rows <- list(
  bridge_cmp("Spawning biomass, total", ssb_tot, ss3_ssb_tot, signed = TRUE),
  bridge_cmp("Recruitment, total", rec_tot, ss3_rec_tot, signed = TRUE)
)
for(a in 1:n_reg) {
  cmp_rows <- c(cmp_rows, list(
    bridge_cmp(paste0("Spawning biomass, ", area_names[a]), ssb_area[a, ], dat$ss3$SSB[, a], signed = TRUE),
    bridge_cmp(paste0("Recruitment, ", area_names[a]), rec_area[a, ], dat$ss3$Rec[, a], signed = TRUE)
  ))
} # end a loop
print(do.call(rbind, cmp_rows), row.names = FALSE, digits = 4)

# Estimated parameters against the assessment's, growth one row per area and sex
growth_est <- exp(pl$ln_growth_pars)
par_names <- c("ln R0", "catchability", "area 2 apportionment (logit)")
par_sporc <- c(pl$ln_global_R0[1], exp(pl$ln_srv_q[1, 1, 1]), pl$rec_region_prop_pars[1, 1])
par_ss3 <- c(mle$ln_R0, exp(mle$ln_q), mle$rec_dist_area2)
for(a in 1:n_reg) for(s in 1:n_sex) {
  par_names <- c(par_names, paste0(c("L1", "L2", "K", "CV1", "CV2"), ", ", area_names[a], " ", sex_names[s]))
  par_sporc <- c(par_sporc, growth_est[1, a, s, ])
  par_ss3 <- c(par_ss3, mle$growth[1, a, s, ])
} # end a and s loops
print(data.frame(quantity = par_names, SPoRC = par_sporc, assessment = par_ss3),
      row.names = FALSE, digits = 6)

# Figures ----------------------------------------------------------------------
# Total stock, the standard bridge figure
p_ts <- bridge_ts_figure(yrs = yrs, ssb = ssb_tot, rec = rec_tot,
                         admb_ssb = ss3_ssb_tot, admb_rec = ss3_rec_tot,
                         label = label, ssb_se = ssb_se_tot, rec_se = rec_se_tot)
ggplot2::ggsave(here("vignettes", "figures", "ad_goa_rex_ts_comparison.png"), p_ts,
                width = 17, height = 9, dpi = 200)

# By area: the same layout with the two areas side by side
ts_rows <- list()
pd_rows <- list()
for(a in 1:n_reg) {
  ts_rows <- c(ts_rows, list(
    data.frame(Year = yrs, value = ssb_area[a, ], se = ssb_se_area[a, ],
               Par = "Spawning Biomass", Area = area_names[a], type = "SPoRC"),
    data.frame(Year = yrs, value = rec_area[a, ], se = rec_se_area[a, ],
               Par = "Recruitment", Area = area_names[a], type = "SPoRC"),
    data.frame(Year = yrs, value = dat$ss3$SSB[, a], se = NA,
               Par = "Spawning Biomass", Area = area_names[a], type = label),
    data.frame(Year = yrs, value = dat$ss3$Rec[, a], se = NA,
               Par = "Recruitment", Area = area_names[a], type = label)
  ))
  pd_rows <- c(pd_rows, list(
    data.frame(Year = yrs, value = 100 * (ssb_area[a, ] - dat$ss3$SSB[, a]) / dat$ss3$SSB[, a],
               Par = "Spawning Biomass", Area = area_names[a]),
    data.frame(Year = yrs, value = 100 * (rec_area[a, ] - dat$ss3$Rec[, a]) / dat$ss3$Rec[, a],
               Par = "Recruitment", Area = area_names[a])
  ))
} # end a loop

ts_df <- dplyr::bind_rows(ts_rows)
ts_df$Area <- factor(ts_df$Area, levels = area_names)
pd_df <- dplyr::bind_rows(pd_rows)
pd_df$Area <- factor(pd_df$Area, levels = area_names)

p_area <- ggplot2::ggplot(ts_df, ggplot2::aes(x = Year, y = value,
                                              ymin = value - 1.96 * se, ymax = value + 1.96 * se,
                                              color = type, fill = type)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(alpha = 0.3, color = NA) +
  ggplot2::facet_grid(Par ~ Area, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  ggplot2::theme_bw(base_size = 20) +
  ggplot2::coord_cartesian(ylim = c(0, NA)) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2), fill = ggplot2::guide_legend(nrow = 2)) +
  ggplot2::labs(x = "Year", y = "Value", color = "Type", fill = "Type")

p_pd_area <- ggplot2::ggplot(pd_df, ggplot2::aes(x = Year, y = value)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
  ggplot2::geom_line(color = "#E69F00") +
  ggplot2::geom_point(size = 2, color = "#E69F00") +
  ggplot2::facet_grid(Par ~ Area, scales = "free_y") +
  ggplot2::theme_bw(base_size = 20) +
  ggplot2::labs(x = "Year", y = "SPoRC vs SS3 (%)")

ggplot2::ggsave(here("vignettes", "figures", "ad_goa_rex_ts_comparison_area.png"),
                patchwork::wrap_plots(p_area, p_pd_area, ncol = 2, widths = c(1.25, 1)),
                width = 20, height = 9, dpi = 200)

# Selectivity in the terminal year, one panel per fleet and sex
sel_rows <- list()
for(s in 1:n_sex) {
  sel_rows <- c(sel_rows, list(bridge_sel_rows(
    ages, est$rep$fish_sel[1, 1, n_yrs, 1, , s, 1], dat$ss3$sel[[1]][[s]],
    paste0(fleet_labels[1], ", ", sex_names[s]), label)))
} # end s loop
for(sf in 1:n_srv) for(s in 1:n_sex) {
  sel_rows <- c(sel_rows, list(bridge_sel_rows(
    ages, est$rep$srv_sel[1, sf, n_yrs, 1, , s, sf], dat$ss3$sel[[sf + 1]][[s]],
    paste0(fleet_labels[sf + 1], ", ", sex_names[s]), label)))
} # end sf and s loops

sel_df <- dplyr::bind_rows(sel_rows)
sel_df$Gear <- factor(sel_df$Gear, levels = unlist(lapply(fleet_labels, function(f) paste0(f, ", ", sex_names))))
ggplot2::ggsave(here("vignettes", "figures", "ad_goa_rex_sel_comparison.png"),
                bridge_sel_figure(sel_df, nrow = 3), width = 12, height = 12, dpi = 170)

# Growth: mean length at age at the survey's timing, per area and sex, estimated
# here against the assessment's own curves
growth_rows <- list()
for(a in 1:n_reg) for(s in 1:n_sex) {
  growth_rows <- c(growth_rows, list(bridge_sel_rows(
    ages, est$rep$mean_LAA_srv[1, a, 1, 1, , s, 1], dat$ss3$growth[[a]][[s]]$Len_Mid,
    paste0(area_names[a], ", ", sex_names[s]), label)))
} # end a and s loops

growth_df <- dplyr::bind_rows(growth_rows)
growth_df$Gear <- factor(growth_df$Gear, levels = unlist(lapply(area_names, function(a) paste0(a, ", ", sex_names))))
ggplot2::ggsave(here("vignettes", "figures", "ad_goa_rex_growth_comparison.png"),
                bridge_sel_figure(growth_df, ylab = "Mean length at age (cm)", nrow = 2),
                width = 12, height = 9, dpi = 170)

saveRDS(list(est_rep = est$rep, optim = est$optim, parList = pl,
             sdrep_summary = if(is.null(est$sdrep)) NULL else summary(est$sdrep)),
        here("dev", "rex_bridge", "output", "rex_sporc_refit.rds"))

# The regression test builds this same model from the test helper. The two
# specifications are written out separately, so compare them rather than trusting
# that they still agree.
helper_path <- here("tests", "testthat", "helper-bridge_goa_rex.R")
if(file.exists(helper_path)) {
  source(helper_path, local = TRUE)
  ref <- seed_goa_rex_mle(suppressWarnings(suppressMessages(build_goa_rex_input(dat))), dat)
  drifted <- c(data = !isTRUE(all.equal(ref$data, data)),
               par = !isTRUE(all.equal(ref$par, parameters)),
               map = !isTRUE(all.equal(ref$map, mapping)))
  if(any(drifted)) {
    warning("this script and helper-bridge_goa_rex.R have drifted apart in: ",
            paste(names(drifted)[drifted], collapse = ", "), call. = FALSE)
  } else {
    cat("\nspecification matches helper-bridge_goa_rex.R\n")
  }
}

cat("\nfigures written to vignettes/figures\n")
