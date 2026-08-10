# Purpose: Build sgl_rg_bsai_nork_data, the container for the BSAI northern
#          rockfish case study. Everything the case study and its regression
#          tests need lives in the object: model inputs, the ADMB maximum
#          likelihood estimate used as a seed, and the ADMB output the bridge
#          is compared against.
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26
#
# Sources, all under dev/dev_data/bsai_nork23:
#   nork23.dat         2023 assessment input file
#   northern.ctl       control file, carries the stage 1 haul counts
#   compweights.ctl    applied McAllister Ianelli composition weights
#   nork23.par         ADMB parameter file, read at full precision
#   nork23.rdat        ADMB report object, the comparison target

library(here)

src <- function(f) here("dev", "dev_data", "bsai_nork23", f)
lines <- readLines(src("nork23.dat"))
ctl <- readLines(src("northern.ctl"))
par_lines <- readLines(src("nork23.par"))
mod <- dget(src("nork23.rdat"))

read_nums <- function(lns) {
  txt <- gsub("#.*", "", lns)
  vals <- as.numeric(unlist(strsplit(trimws(paste(txt, collapse = " ")), "[[:space:]]+")))
  vals[!is.na(vals)]
}
read_matrix <- function(start, n, ncol) {
  m <- matrix(NA_real_, nrow = n, ncol = ncol)
  for(r in seq_len(n)) m[r, ] <- read_nums(lines[start + r - 1])
  m
}
# Comp rows in the .dat are interleaved with single "#year" comment lines, so
# the data sit on every other line.
read_comp_stack <- function(start, n, ncol) {
  m <- matrix(NA_real_, nrow = n, ncol = ncol)
  dl <- start
  for(r in seq_len(n)) {
    m[r, ] <- read_nums(lines[dl])
    dl <- dl + 2L
  } # end r loop
  m
}
get_par <- function(name) {
  i <- grep(paste0("^# ", name, ":"), par_lines)
  as.numeric(unlist(strsplit(trimws(par_lines[i + 1]), "[[:space:]]+")))
}

# Dimensions -------------------------------------------------------------------
# The model carries 43 age classes (3 to 45) while the compositions are
# reported over 38 (3 to 40), so the ageing error matrix maps the model ages
# onto the shorter observed range.
n_pop <- 1
n_regions <- 1
n_seas <- 1
n_sexes <- 1
n_fish_fleets <- 1
n_srv_fleets <- 1

years <- 1977:2023
ages <- 3:45
obs_ages <- 3:40
lens <- 15:38

n_yrs <- length(years)
n_ages <- length(ages)
n_obs_ages <- length(obs_ages)
n_lens <- length(lens)

# Priors -----------------------------------------------------------------------
mean_q <- read_nums(lines[338])
cv_q <- read_nums(lines[341])
mean_M <- read_nums(lines[345])
cv_M <- read_nums(lines[347])
sigmaR <- read_nums(lines[334])
spawn_mo <- read_nums(lines[315])

# Biologicals ------------------------------------------------------------------
# Weight at age is genuinely year varying here, and the .dat carries separate
# population and fishery blocks that differ by up to 50 percent at young ages.
# The population block prices spawning biomass and the survey; the fishery
# block prices the catch.
pop_waa <- read_matrix(219, n_yrs, n_ages)
WAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
WAA[1, 1, , 1, , 1] <- pop_waa

fish_waa <- read_matrix(267, n_yrs, n_ages)
WAA_fish <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
WAA_fish[1, 1, , 1, , 1, 1] <- fish_waa

# Maturity is estimated inside the ADMB template from the TenBrink and Shaw
# binomial data rather than read from the .dat, so the fitted logistic is
# rebuilt from mat_beta1 and mat_beta2 in the parameter file and held fixed in
# SPoRC.
mat_beta1 <- get_par("mat_beta1")
mat_beta2 <- get_par("mat_beta2")
maa_vec <- 1 / (1 + exp(-(mat_beta1 + mat_beta2 * ages)))
MatAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
MatAA[1, 1, , 1, , 1] <- matrix(rep(maa_vec, each = n_yrs), nrow = n_yrs)

# The .dat carries the size at age matrix as [length bin x age] with each age
# column a distribution over length bins, so the columns are normalised.
size_age_raw <- read_matrix(90, n_lens, n_ages)
col_sums <- colSums(size_age_raw)
col_sums[col_sums == 0] <- 1
size_age_mat <- sweep(size_age_raw, 2, col_sums, "/")
SizeAgeTrans <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes))
for(y in seq_len(n_yrs)) SizeAgeTrans[1, 1, y, 1, , , 1] <- size_age_mat

# Ageing error is stored as [observed age x true age]. The ADMB template
# column-normalises it to P(obs | true); SPoRC multiplies predicted numbers at
# age on the right, so it needs [true x obs], and the observed plus group
# collapses rows 38 to 43 onto row 38.
age_error_raw <- read_matrix(46, n_ages, n_ages)
col_sums <- colSums(age_error_raw)
col_sums[col_sums == 0] <- 1
age_error_raw <- sweep(age_error_raw, 2, col_sums, "/")
age_error_collapsed <- rbind(
  age_error_raw[1:(n_obs_ages - 1), ],
  colSums(age_error_raw[n_obs_ages:n_ages, ])
)
AgeingError_mat <- t(age_error_collapsed)
AgeingError_mat <- AgeingError_mat / rowSums(AgeingError_mat)
AgeingError <- array(NA_real_, dim = c(n_yrs, n_ages, n_obs_ages))
for(y in seq_len(n_yrs)) AgeingError[y, , ] <- AgeingError_mat

# Catch ------------------------------------------------------------------------
catch_vec <- read_nums(lines[32])
stopifnot(length(catch_vec) == n_yrs)
ObsCatch <- array(catch_vec, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
UseCatch <- array(1L, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))

# Survey index -----------------------------------------------------------------
# The Aleutian Islands bottom trawl survey is the only index. The .dat reports
# a standard deviation on the arithmetic scale, which is carried to SPoRC as a
# lognormal coefficient of variation.
srv_yrs <- read_nums(lines[36])
srv_obs <- read_nums(lines[38])
srv_sd <- read_nums(lines[40])
srv_ind <- as.integer(years %in% srv_yrs)

ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
UseSrvIdx <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx[1, srv_ind == 1, 1, 1] <- srv_obs
ObsSrvIdx_SE[1, srv_ind == 1, 1, 1] <- srv_sd / srv_obs
UseSrvIdx[1, srv_ind == 1, 1, 1] <- 1L

# Composition helpers ----------------------------------------------------------
make_obs_comp <- function(obs_mat, ind_vec, n_bins, n_fl) {
  arr <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_bins, n_sexes, n_fl))
  iy <- 1L
  for(t in seq_len(n_yrs)) {
    if(ind_vec[t] == 1L) {
      arr[1, t, 1, , 1, 1] <- obs_mat[, iy]
      iy <- iy + 1L
    }
  } # end t loop
  arr
}
make_use_comp <- function(ind_vec, n_fl) {
  arr <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_fl))
  arr[1, , 1, 1] <- as.integer(ind_vec)
  arr
}
make_iss_comp <- function(iss_vec, ind_vec, n_fl) {
  arr <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fl))
  iy <- 1L
  for(t in seq_len(n_yrs)) {
    if(ind_vec[t] == 1L) {
      arr[1, t, 1, 1, 1] <- iss_vec[iy]
      iy <- iy + 1L
    }
  } # end t loop
  arr
}
make_wt_comp <- function(ind_vec, wt_val, n_fl) {
  arr <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fl))
  arr[1, ind_vec == 1, 1, 1, 1] <- wt_val
  arr
}

# Compositions -----------------------------------------------------------------
# The stage 1 input sample sizes are the haul counts from the control file, and
# the weights are the applied McAllister Ianelli multipliers from
# compweights.ctl, which the template folds into the sample sizes.
fish_age_stg1 <- read_nums(ctl[27])
fish_len_stg1 <- read_nums(ctl[31])
srv_age_stg1 <- read_nums(ctl[35])
comp_wts <- read_nums(readLines(src("compweights.ctl")))
wt_fish_age <- comp_wts[1]
wt_fish_len <- comp_wts[2]
wt_srv_age <- comp_wts[3]

fish_age_yrs <- read_nums(lines[117])
fish_age_ind <- as.integer(years %in% fish_age_yrs)
fish_age_mat <- read_comp_stack(120, length(fish_age_yrs), n_obs_ages)
fish_age_mat <- t(fish_age_mat / rowSums(fish_age_mat))
ObsFishAgeComps <- make_obs_comp(fish_age_mat, fish_age_ind, n_obs_ages, n_fish_fleets)
UseFishAgeComps <- make_use_comp(fish_age_ind, n_fish_fleets)
ISS_FishAgeComps <- make_iss_comp(fish_age_stg1, fish_age_ind, n_fish_fleets)
Wt_FishAgeComps <- make_wt_comp(fish_age_ind, wt_fish_age, n_fish_fleets)

fish_len_yrs <- read_nums(lines[156])
fish_len_ind <- as.integer(years %in% fish_len_yrs)
fish_len_mat <- read_comp_stack(159, length(fish_len_yrs), n_lens)
fish_len_mat <- t(fish_len_mat / rowSums(fish_len_mat))
ObsFishLenComps <- make_obs_comp(fish_len_mat, fish_len_ind, n_lens, n_fish_fleets)
UseFishLenComps <- make_use_comp(fish_len_ind, n_fish_fleets)
ISS_FishLenComps <- make_iss_comp(fish_len_stg1, fish_len_ind, n_fish_fleets)
Wt_FishLenComps <- make_wt_comp(fish_len_ind, wt_fish_len, n_fish_fleets)

srv_age_yrs <- read_nums(lines[183])
srv_age_ind <- as.integer(years %in% srv_age_yrs)
srv_age_mat <- read_comp_stack(186, length(srv_age_yrs), n_obs_ages)
srv_age_mat <- t(srv_age_mat / rowSums(srv_age_mat))
ObsSrvAgeComps <- make_obs_comp(srv_age_mat, srv_age_ind, n_obs_ages, n_srv_fleets)
UseSrvAgeComps <- make_use_comp(srv_age_ind, n_srv_fleets)
ISS_SrvAgeComps <- make_iss_comp(srv_age_stg1, srv_age_ind, n_srv_fleets)
Wt_SrvAgeComps <- make_wt_comp(srv_age_ind, wt_srv_age, n_srv_fleets)

# The survey contributes no length compositions.
ObsSrvLenComps <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets))
UseSrvLenComps <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))
Wt_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))

stopifnot(
  sum(srv_ind) == length(srv_age_yrs),
  all(abs(colSums(fish_age_mat) - 1) < 1e-6),
  all(abs(colSums(fish_len_mat) - 1) < 1e-6),
  all(abs(colSums(srv_age_mat) - 1) < 1e-6),
  all(abs(rowSums(AgeingError_mat) - 1) < 1e-6)
)

# ADMB maximum likelihood estimate ---------------------------------------------
# Read at full precision from the parameter file. rec_dev covers 1977 to 2020,
# the last three years taking mean recruitment; fydev carries the 37 initial
# age deviations, whose last value is shared by ages beyond the observed range.
# The domestic fishery selectivity parameters in the file are inactive: the
# report's fishery selectivity equals the foreign fishery logistic in every
# year.
mle <- list(
  M = get_par("M"),
  mean_log_rec = get_par("mean_log_rec"),
  log_rinit = get_par("log_rinit"),
  log_avg_fmort = get_par("log_avg_fmort"),
  q_srv = get_par("q_srv3"),
  sel_a50_fish = get_par("sel_a50_forfish"),
  sel_aslope_fish = get_par("sel_aslope_forfish"),
  sel_a50_srv = get_par("sel_a50_srv3"),
  sel_aslope_srv = get_par("sel_aslope_srv3"),
  fmort_dev = get_par("fmort_dev"),
  rec_dev = get_par("rec_dev"),
  fydev = get_par("fydev"),
  mat_beta1 = mat_beta1,
  mat_beta2 = mat_beta2
)
stopifnot(length(mle$fmort_dev) == n_yrs,
          length(mle$rec_dev) == n_yrs - 3,
          length(mle$fydev) == n_obs_ages - 1)

# ADMB output ------------------------------------------------------------------
# The comparison target. natage is reported over the 38 observed ages with the
# model's ages 41 to 45 pooled into the final column.
admb <- list(
  NAA = mod$natage,
  SSB = mod$t.series$spbiom,
  Rec = mod$t.series$a3recs,
  Fmort = mod$t.series$fmort,
  TotBiom = mod$t.series$totbiom,
  pred_srv = mod$t.series$predsrv,
  sel_fsh = if(is.null(dim(mod$selfish))) mod$selfish else mod$selfish[nrow(mod$selfish), ],
  sel_srv = mod$srvsel$srvsel,
  datalikecomp = mod$datalikecomp,
  pen_likecomp = mod$pen_likecomp
)

sgl_rg_bsai_nork_data <- list(
  years = years,
  ages = ages,
  obs_ages = obs_ages,
  lens = lens,
  n_regions = n_regions,
  n_sexes = n_sexes,
  n_fish_fleets = n_fish_fleets,
  n_srv_fleets = n_srv_fleets,
  n_seas = n_seas,
  n_pop = n_pop,
  natal_region = 1,

  WAA = WAA,
  WAA_fish = WAA_fish,
  MatAA = MatAA,
  SizeAgeTrans = SizeAgeTrans,
  AgeingError = AgeingError,

  ObsCatch = ObsCatch,
  UseCatch = UseCatch,

  ObsSrvIdx = ObsSrvIdx,
  ObsSrvIdx_SE = ObsSrvIdx_SE,
  UseSrvIdx = UseSrvIdx,

  ObsFishAgeComps = ObsFishAgeComps,
  UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,
  Wt_FishAgeComps = Wt_FishAgeComps,
  ObsFishLenComps = ObsFishLenComps,
  UseFishLenComps = UseFishLenComps,
  ISS_FishLenComps = ISS_FishLenComps,
  Wt_FishLenComps = Wt_FishLenComps,

  ObsSrvAgeComps = ObsSrvAgeComps,
  UseSrvAgeComps = UseSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,
  Wt_SrvAgeComps = Wt_SrvAgeComps,
  ObsSrvLenComps = ObsSrvLenComps,
  UseSrvLenComps = UseSrvLenComps,
  ISS_SrvLenComps = ISS_SrvLenComps,
  Wt_SrvLenComps = Wt_SrvLenComps,

  yrs_srv = srv_yrs,
  yrs_fish_age = fish_age_yrs,
  yrs_fish_len = fish_len_yrs,
  yrs_srv_age = srv_age_yrs,

  sigmaR = sigmaR,
  spawn_mo = spawn_mo,
  mean_M = mean_M,
  cv_M = cv_M,
  mean_q = mean_q,
  cv_q = cv_q,
  catch_wt = 200,
  fmort_wt = 0.1,

  mle = mle,
  admb = admb
)

usethis::use_data(sgl_rg_bsai_nork_data, internal = FALSE, overwrite = TRUE)
