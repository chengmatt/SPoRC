# Purpose: Build sgl_rg_wc_sablefish_data, the inputs, SS3 seed, and SS3 output for the West Coast sablefish case study
# Creator: Matthew LH. Cheng
# Date Created: 8/21/26
#
# Sources, the 2025 assessment's Stock Synthesis 3.30.23 run:
#   2025_sablefish_dat.ss   data file
#   2025_sablefish_ctl.ss   control file, holds the blocks, the variance
#                           adjustments and the phases
#   wtatage.ss              empirical weight and fecundity at age
#   ss3.par                 parameter file, read at twelve significant digits;
#                           Report.sso and r4ss's parameter table have six
#   Report.sso              read through r4ss, the comparison target
#
# West Coast sablefish differs from the Alaska case studies: age-0 spawners in year one, an additive
# extra survey sd, a recruitment-deviation index, duplicate fleets, and SS3 double normal selectivity

library(here)
library(dplyr)
library(r4ss)

ss3_dir <- "/Users/matthewcheng/Downloads/2025_sablefish_model_files"

replist <- SS_output(dir = ss3_dir, verbose = FALSE, printstats = FALSE, covar = FALSE)
inputs <- SS_read(dir = ss3_dir, verbose = FALSE)
pars_ss3 <- SS_readpar_3.30(
  parfile = file.path(ss3_dir, "ss3.par"),
  datsource = file.path(ss3_dir, "2025_sablefish_dat.ss"),
  ctlsource = file.path(ss3_dir, "2025_sablefish_ctl.ss"),
  verbose = FALSE
)
obj_ss3 <- as.numeric(sub(".*Objective function value = ([-0-9.eE+]+).*", "\\1",
                          readLines(file.path(ss3_dir, "ss3.par"), n = 1)))

yrs <- inputs$dat$styr:inputs$dat$endyr
n_yrs <- length(yrs)
ages <- 0:inputs$dat$Nages
n_ages <- length(ages)
age_cols <- as.character(ages)
obs_ages <- inputs$dat$agebin_vector
n_obs_ages <- length(obs_ages)
n_sexes <- 2
fleet_names <- replist$FleetNames

# Fleets ----------------------------------------------------------------------
# SS3 fleets 1-6 are catch (three gears plus discards), 7-10 trawl surveys, 11 the recruitment index.
# trawl fishery and bottom trawl survey each carry two comp sources, so each gets a duplicate fleet
n_fish <- 7
n_srv <- 6
fish_src <- c(1:6, 1) # the Stock Synthesis fleet each SPoRC fishery fleet draws from
fish_sex <- c(3, 3, 3, 0, 0, 0, 0) # 3 = sexed compositions, 0 = unsexed
# surveys 1-4 are the trawl surveys, 5 the unsexed composition twin of the last
# of them, and 6 the recruitment index, which observes the deviations directly
srv_src <- c(7:10, 10, 11)
srv_sex <- c(3, 3, 3, 3, 0, 0)

# Biologicals -----------------------------------------------------------------
# wtatage.ss fleet -1 is the population weight at age and fleet -2 the female
# fecundity at age; every fleet specific weight equals the population's
wt <- replist$wtatage
get_wt <- function(fleet, sex) {
  df <- wt[wt$fleet == fleet & wt$sex == sex, ]
  unname(as.matrix(df[match(yrs, df$year), age_cols]))
}
stopifnot(all(sapply(0:11, function(f) isTRUE(all.equal(get_wt(f, 1), get_wt(-1, 1))))))

WAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes))
for(s in 1:n_sexes) WAA[1,1,,1,,s] <- get_wt(-1, s)

# maturity is fecundity over weight, so spawning biomass is numbers times fecundity. age 0 fish are
# spawners in the equilibrium year only, since later years form SSB before recruits settle
MatAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes))
MatAA[1,1,,1,,1] <- get_wt(-2, 1) / get_wt(-1, 1)
mat_age0_yr1 <- MatAA[1,1,1,1,1,1]
MatAA[1,1,,1,1,1] <- 0

# the ageing error key is printed with the observed bins descending, one block per sex; transposed
# it maps true ages onto observed bins, the last of which accumulates
AgeingError <- t(replist$AAK[1, n_obs_ages:1, ])
stopifnot(max(abs(rowSums(AgeingError) - 1)) < 1e-5)
AgeingError <- AgeingError / rowSums(AgeingError)

# Catch -----------------------------------------------------------------------
# years without a catch record are closures. the duplicate trawl fleet has no catch observation
# at all, which leaves its fishing mortality estimated rather than forced to zero
ObsCatch <- UseCatch <- array(0, dim = c(1, n_yrs, 1, n_fish))
for(f in 1:6) {
  cf <- inputs$dat$catch %>% dplyr::filter(fleet == f, year %in% yrs)
  ObsCatch[1, match(cf$year, yrs), 1, f] <- cf$catch
  UseCatch[1, match(cf$year, yrs), 1, f] <- 1
} # end f loop
ObsCatch[1,,1,n_fish] <- NA

# Survey indices --------------------------------------------------------------
cpue <- replist$cpue
extra_sd <- pars_ss3$Q_parms[grep("^Q_extraSD", rownames(pars_ss3$Q_parms)), "ESTIM"]
stopifnot(max(abs(cpue$SE[cpue$Fleet %in% 7:10] - cpue$SE_input[cpue$Fleet %in% 7:10] -
                    extra_sd[match(cpue$Fleet[cpue$Fleet %in% 7:10], 7:10)])) < 1e-6)

ObsSrvIdx <- ObsSrvIdx_SE <- array(NA_real_, dim = c(1, n_yrs, 1, n_srv))
UseSrvIdx <- array(0, dim = c(1, n_yrs, 1, n_srv))
for(sf in c(1:4, 6)) {
  ci <- cpue %>% dplyr::filter(Fleet == srv_src[sf])
  ObsSrvIdx[1, match(ci$Yr, yrs), 1, sf] <- ci$Obs
  ObsSrvIdx_SE[1, match(ci$Yr, yrs), 1, sf] <- ci$SE
  UseSrvIdx[1, match(ci$Yr, yrs), 1, sf] <- 1
} # end sf loop

# Age compositions ------------------------------------------------------------
# rows with a negative fleet are excluded by the assessment. sample size is the input size times
# the control file's variance adjustment, over 1 + nbins * comp constant: SS3's multinomial scale
addtocomp <- 1e-3
ac <- inputs$dat$agecomp %>% dplyr::filter(fleet > 0)
va <- replist$Age_Comp_Fit_Summary
var_adj <- stats::setNames(va$Curr_Var_Adj, va$Fleet)

build_comps <- function(src, sexcode, n_fleets) {
  Obs <- array(0, dim = c(1, n_yrs, 1, n_obs_ages, n_sexes, n_fleets))
  Use <- array(0, dim = c(1, n_yrs, 1, n_fleets))
  ISS <- array(0, dim = c(1, n_yrs, 1, n_sexes, n_fleets))
  for(f in seq_along(src)) {
    rows <- ac %>% dplyr::filter(fleet == src[f], sex == sexcode[f])
    n_bins <- if(sexcode[f] == 3) 2 * n_obs_ages else n_obs_ages
    for(i in seq_len(nrow(rows))) {
      y <- match(rows$year[i], yrs)
      Obs[1,y,1,,1,f] <- as.numeric(rows[i, paste0("f", obs_ages)])
      if(sexcode[f] == 3) Obs[1,y,1,,2,f] <- as.numeric(rows[i, paste0("m", obs_ages)])
      Use[1,y,1,f] <- 1
      ISS[1,y,1,,f] <- rows$Nsamp[i] * var_adj[as.character(src[f])] / (1 + n_bins * addtocomp)
    } # end i loop
  } # end f loop
  list(Obs = Obs, Use = Use, ISS = ISS)
}
fish_comps <- build_comps(fish_src, fish_sex, n_fish)
srv_comps <- build_comps(srv_src, srv_sex, n_srv)

# Selectivity -----------------------------------------------------------------
# every fleet is age based double normal, three with time blocks, hook and line with male offsets,
# pot mirroring it. the realized surfaces go in fixed, stored as distinct block rows plus start year
asel <- replist$ageselex %>% dplyr::filter(Factor == "Asel2", Yr %in% yrs)
get_sel <- function(fleet) {
  out <- array(0, dim = c(n_yrs, n_ages, n_sexes))
  for(s in 1:n_sexes) {
    df <- asel %>% dplyr::filter(Fleet == fleet, Sex == s)
    out[,,s] <- as.matrix(df[match(yrs, df$Yr), age_cols])
  } # end s loop
  out
}
# distinct rows, indexed by the year they start
compress_sel <- function(src) {
  out <- lapply(src, function(fleet) {
    s <- get_sel(fleet)
    key <- apply(matrix(s, nrow = n_yrs), 1, paste, collapse = "|")
    starts <- which(c(TRUE, key[-1] != key[-n_yrs]))
    list(blk_yr = starts, sel = s[starts, , , drop = FALSE])
  })
  out
}
fish_sel_blocks_ss3 <- compress_sel(fish_src)
srv_sel_blocks_ss3 <- compress_sel(srv_src[1:5])
# the recruitment index observes the deviations, so no curve of its own is ever
# read; a flat one keeps the array dimensions right
srv_sel_blocks_ss3[[6]] <- list(blk_yr = 1, sel = array(1, dim = c(1, n_ages, n_sexes)))

# the selectivity parameters behind those surfaces, so the case study can estimate the curves.
# recording which row supplies each block lets the estimation share one parameter across blocks
sel_par_table <- function(src, blk_yr) {
  base_rows <- inputs$ctl$age_selex_parms
  tv_rows <- inputs$ctl$age_selex_parms_tv
  nm_f <- gsub(" +$", "", fleet_names[src])
  n_blk <- length(blk_yr)
  pars <- matrix(NA_real_, 6, n_blk)
  est <- matrix(FALSE, 6, n_blk)
  src_id <- matrix(NA_character_, 6, n_blk)
  for(k in 1:6) {
    lab <- sprintf("AgeSel_P_%d_%s(%d)", k, nm_f, src)
    base_val <- pars_ss3$S_parms[lab, "ESTIM"]
    base_ph <- base_rows[lab, "PHASE"]
    blkpat <- base_rows[lab, "Block"]
    # the years each sub-block of this parameter's pattern covers, and the row
    # supplying them
    covered <- rep(NA_character_, n_yrs)
    if(!is.na(blkpat) && blkpat > 0) {
      design <- matrix(inputs$ctl$Block_Design[[blkpat]], nrow = 2)
      for(j in seq_len(ncol(design))) {
        rep_lab <- sprintf("%s_BLK%drepl_%d", lab, blkpat, design[1, j])
        if(rep_lab %in% rownames(pars_ss3$S_parms))
          covered[which(yrs >= design[1, j] & yrs <= design[2, j])] <- rep_lab
      } # end j loop
    }
    for(b in seq_len(n_blk)) {
      lab_b <- covered[blk_yr[b]]
      if(is.na(lab_b)) {
        pars[k, b] <- base_val; est[k, b] <- base_ph > 0; src_id[k, b] <- lab
      } else {
        pars[k, b] <- pars_ss3$S_parms[lab_b, "ESTIM"]
        est[k, b] <- tv_rows[lab_b, "PHASE"] > 0
        src_id[k, b] <- lab_b
      }
    } # end b loop
  } # end k loop
  list(pars = pars, est = est, src_id = src_id)
}
sel_fish <- lapply(seq_along(fish_src), function(f) {
  # the pot fleet mirrors hook and line and has no parameter rows of its own
  if(fish_src[f] == 3) return(NULL)
  sel_par_table(fish_src[f], fish_sel_blocks_ss3[[f]]$blk_yr)
})
sel_fish[[3]] <- sel_fish[[2]]
sel_srv <- lapply(seq_len(5), function(sf) sel_par_table(srv_src[sf], srv_sel_blocks_ss3[[sf]]$blk_yr))
# nothing to estimate for the recruitment index's curve, which is never read
sel_srv[[6]] <- list(
  pars = matrix(0, 6, 1),
  est = matrix(FALSE, 6, 1),
  src_id = matrix("recruitment index, unused", 6, 1)
)

# The hook and line fleet's male parameters are offsets on the female's, the
# last of them a scale on the whole curve. The pot fleet mirrors all of it.
male_lab <- sprintf("AgeSel_PMalOff_%d_Hook_and_Line(2)", 1:5)
sel_male <- list(value = stats::setNames(pars_ss3$S_parms[male_lab, "ESTIM"],
                                         c("peak", "ascend", "descend", "final", "scale")),
                 est = stats::setNames(inputs$ctl$age_selex_parms[male_lab, "PHASE"] > 0,
                                       c("peak", "ascend", "descend", "final", "scale")))

# Fishing mortality, recruitment and the bias ramp -----------------------------
ts <- replist$timeseries %>% dplyr::filter(Yr %in% yrs) %>% dplyr::arrange(Yr)
F_ss3 <- sapply(1:6, function(f) ts[[paste0("F:_", f)]])
rec_tab <- replist$recruit %>% dplyr::filter(Yr %in% yrs) %>% dplyr::arrange(Yr)
bias_adj <- rec_tab$biasadjuster
recdev <- numeric(n_yrs)
recdev[match(pars_ss3$recdev2[, "year"], yrs)] <- pars_ss3$recdev2[, "recdev"]
recdev[match(max(yrs), yrs)] <- pars_ss3$recdev_forecast[pars_ss3$recdev_forecast[, "year"] == max(yrs), "recdev"]

# The recruitment index is a likelihood on the deviations themselves
ri <- cpue %>% dplyr::filter(Fleet == 11)

# Numbers at age, expected compositions and the likelihood table ----------------
nat <- replist$natage %>% dplyr::filter(`Beg/Mid` == "B", Yr %in% yrs)
NAA_ss3 <- array(0, dim = c(n_yrs, n_ages, n_sexes))
for(s in 1:n_sexes) {
  df <- nat %>% dplyr::filter(Sex == s)
  NAA_ss3[,,s] <- as.matrix(df[match(yrs, df$Yr), age_cols])
} # end s loop

agedb <- replist$agedbase %>% dplyr::select(Fleet, Sexes, Yr, Sex, Bin, Obs, Exp, Nsamp_in, Nsamp_adj)
lbf <- replist$likelihoods_by_fleet
lik_ss3 <- function(label, fleet) as.numeric(lbf[lbf$Label == label, fleet_names[fleet]])
dq <- replist$derived_quants

sgl_rg_wc_sablefish_data <- list(

  # Dimensions
  years = yrs, ages = ages, obs_ages = obs_ages, lens = NA,
  n_regions = 1, n_sexes = n_sexes, n_fish_fleets = n_fish, n_srv_fleets = n_srv,
  n_seas = 1, n_pop = 1, natal_region = 1,
  fleet_names = fleet_names, fish_src = fish_src, fish_sex = fish_sex,
  srv_src = srv_src, srv_sex = srv_sex,

  # Data
  ObsCatch = ObsCatch, UseCatch = UseCatch,
  ObsSrvIdx = ObsSrvIdx, ObsSrvIdx_SE = ObsSrvIdx_SE, UseSrvIdx = UseSrvIdx,
  ObsFishAgeComps = fish_comps$Obs, UseFishAgeComps = fish_comps$Use, ISS_FishAgeComps = fish_comps$ISS,
  ObsSrvAgeComps = srv_comps$Obs, UseSrvAgeComps = srv_comps$Use, ISS_SrvAgeComps = srv_comps$ISS,
  WAA = WAA, MatAA = MatAA, mat_age0_yr1 = mat_age0_yr1, AgeingError = AgeingError,

  # Settings the assessment fixes
  addtocomp = addtocomp, sigmaC = 0.01, t_srv = 0.5, t_spawn = 0,
  sigmaR = pars_ss3$SR_parms["SR_sigmaR", "ESTIM"],
  steepness = pars_ss3$SR_parms["SR_BH_steep", "ESTIM"],
  # the four breakpoints of the assessment's bias ramp, in deviation index
  # space rather than calendar years, which is how SPoRC reads them
  bias_year = {
    brk <- c(inputs$ctl$last_early_yr_nobias_adj, inputs$ctl$first_yr_fullbias_adj,
             inputs$ctl$last_yr_fullbias_adj, inputs$ctl$first_recent_yr_nobias_adj)
    stopifnot(length(brk) == 4, !any(is.na(brk)))
    brk - yrs[1] + 1
  },
  max_bias_adj = inputs$ctl$max_bias_adj,
  yrs_rec_est = min(pars_ss3$recdev2[, "year"]):max(yrs),
  M_prior = data.frame(
    popblk = 1,
    regionblk = 1,
    yearblk = 1,
    ageblk = 1,
    sexblk = 1,
    mu = exp(-2.631),
    sd = 0.31
  ),
  fish_sel_blocks_ss3 = fish_sel_blocks_ss3, srv_sel_blocks_ss3 = srv_sel_blocks_ss3,
  sel_fish = sel_fish, sel_srv = sel_srv, sel_male = sel_male,
  rec_idx = data.frame(yr = ri$Yr, obs = ri$Obs, se = ri$SE),

  # The assessment's maximum likelihood estimate, at twelve significant digits
  mle = list(
    ln_R0 = pars_ss3$SR_parms["SR_LN(R0)", "ESTIM"],
    M = pars_ss3$MG_parms["NatM_p_1_Fem_GP_1", "ESTIM"],
    ln_srv_q = unname(pars_ss3$Q_parms[grep("^LnQ_base", rownames(pars_ss3$Q_parms)), "ESTIM"]),
    extra_sd = unname(extra_sd),
    q_rec_idx = pars_ss3$Q_parms["Q_base_Recruitment_Index(11)", "ESTIM"],
    recdev = recdev,
    bias_adj = bias_adj,
    Fmort = F_ss3,
    sel_pars = pars_ss3$S_parms,
    objective = obj_ss3
  ),

  # The assessment's output, the comparison target
  ss3 = list(
    SSB = ts$SpawnBio,
    Rec = ts$Recruit_0,
    Bio_all = ts$Bio_all,
    NAA = NAA_ss3,
    SSB_sd = dq$StdDev[match(paste0("SSB_", yrs), dq$Label)],
    Rec_sd = dq$StdDev[match(paste0("Recr_", yrs), dq$Label)],
    pred_catch = sapply(1:6, function(f) {
      cc <- replist$catch %>% dplyr::filter(Fleet == f, Yr %in% yrs) %>% dplyr::arrange(Yr)
      cc$Exp[match(yrs, cc$Yr)] }),
    pred_idx = cpue %>% dplyr::filter(Fleet %in% 7:11) %>% dplyr::select(Fleet, Yr, Obs, Exp, SE),
    agedbase = agedb,
    lik = list(catch = sum(sapply(1:6, function(f) lik_ss3("Catch_like", f))),
               index = sapply(7:11, function(f) lik_ss3("Surv_like", f)),
               age = sapply(1:10, function(f) lik_ss3("Age_like", f)),
               recruitment = replist$likelihoods_used["Recruitment", "values"],
               forecast_recruitment = replist$likelihoods_used["Forecast_Recruitment", "values"],
               priors = replist$likelihoods_used["Parm_priors", "values"],
               total = replist$likelihoods_used["TOTAL", "values"])
  )
)

usethis::use_data(sgl_rg_wc_sablefish_data, internal = FALSE, overwrite = TRUE, compress = "xz")
