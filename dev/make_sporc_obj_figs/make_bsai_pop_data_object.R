# Purpose: Build sgl_rg_bsai_pop_data, the container for the BSAI Pacific ocean
#          perch case study. Everything the case study and its regression tests
#          need lives in the object: model inputs, the ADMB maximum likelihood
#          estimate used as a seed, and the ADMB output the bridge is compared
#          against.
# Creator: Matthew LH. Cheng
# Date Created: 8/8/26
#
# Sources, all under dev/dev_data/bsai_pop24:
#   pop24.dat          2024 assessment input file
#   pop24.ctl          control file, carries the stage 1 sample sizes and lambdas
#   compweights.ctl    applied McAllister Ianelli composition weights, read by
#                      the template separately from pop24.ctl
#   pop24.par          ADMB parameter file, read at full precision
#   pop24.rdat         ADMB report object, the comparison target
#
# POP differs from the BSAI northern rockfish case study in four ways that the
# object has to carry: fishery selectivity is a bicubic spline over a 5 by 5
# node grid rather than a logistic, there are two survey fleets rather than one,
# the first year sits in equilibrium under a fixed historical F, and weight at
# age is a single vector broadcast across years.

library(here)

src <- function(f) here("dev", "dev_data", "bsai_pop24", f)
lines <- readLines(src("pop24.dat"))
ctl <- readLines(src("pop24.ctl"))
par_lines <- readLines(src("pop24.par"))
mod <- dget(src("pop24.rdat"))

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
# Parameter vectors in the .par can run over several lines, so read until the
# next comment header.
get_par <- function(name) {
  i <- grep(paste0("^# ", name, ":"), par_lines)
  if(length(i) != 1) stop("parameter '", name, "' not found uniquely in pop24.par")
  j <- i + 1L
  vals <- c()
  while(j <= length(par_lines) && !grepl("^#", par_lines[j])) {
    vals <- c(vals, scan(text = par_lines[j], quiet = TRUE))
    j <- j + 1L
  } # end while
  vals
}

# Dimensions -------------------------------------------------------------------
# The model carries 44 age classes (3 to 46) while the compositions are reported
# over 38 (3 to 40), so the ageing error matrix maps the model ages onto the
# shorter observed range. The bicubic surface is fit over the first 38 ages.
n_pop <- 1
n_regions <- 1
n_seas <- 1
n_sexes <- 1
n_fish_fleets <- 1
n_srv_fleets <- 2

styr <- read_nums(lines[4])
endyr <- read_nums(lines[8])
n_ages <- read_nums(lines[12])
n_obs_ages <- read_nums(lines[14])
nselages <- read_nums(lines[20])
recage <- read_nums(lines[22])
n_lens <- read_nums(lines[24])

n_yrs <- endyr - styr + 1
years <- styr:endyr
ages <- recage:(recage + n_ages - 1)
obs_ages <- recage:(recage + n_obs_ages - 1)
lens <- read_nums(lines[26])

stopifnot(n_yrs == 65, n_ages == 44, n_obs_ages == 38, n_lens == 25, nselages == 38)

# Priors and switches ----------------------------------------------------------
spawn_mo <- read_nums(lines[322])
sigmaR <- read_nums(lines[342])
fixedrec <- read_nums(lines[340])
mean_q <- c(read_nums(lines[425]), read_nums(lines[426]))
cv_q <- c(read_nums(lines[428]), read_nums(lines[429]))
mean_M <- read_nums(lines[440])
cv_M <- read_nums(lines[442])

# The fishery selectivity block. A bicubic spline over a 5 year by 5 age node
# grid, fit from 1964 over the first 38 ages; years before 1964 and ages beyond
# 38 are edge held by the template.
fsh_sel_styr <- read_nums(lines[357])
fsh_yr_nodes <- read_nums(lines[375])
fsh_age_nodes <- read_nums(lines[377])
stopifnot(read_nums(lines[355]) == 3, all(read_nums(lines[383]) == 1))

# Selectivity penalty lambdas from the control file. These map one to one onto
# SPoRC's smoothness penalty terms; the mean centering weight is hardcoded in
# the template rather than read from the control file.
lam_rec <- read_nums(ctl[91])
lam_catch <- read_nums(ctl[93])
lam_dome <- read_nums(ctl[95])
lam_bin_curve <- read_nums(ctl[97])
lam_yr_diff <- read_nums(ctl[99])
lam_yr_curve <- read_nums(ctl[101])
lam_mean_ctr <- 10000

# Biologicals ------------------------------------------------------------------
# Weight at age is a single vector broadcast across years. The "(ages 3 - 25)"
# comment on these lines in the .dat is stale, both rows carry all 44 ages, and
# the population and fishery vectors are identical in this file.
pop_waa <- read_nums(lines[318])
fish_waa <- read_nums(lines[320])
stopifnot(length(pop_waa) == n_ages, length(fish_waa) == n_ages)

WAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
WAA[1, 1, , 1, , 1] <- matrix(rep(pop_waa, each = n_yrs), nrow = n_yrs)

WAA_fish <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
WAA_fish[1, 1, , 1, , 1] <- matrix(rep(fish_waa, each = n_yrs), nrow = n_yrs)

# Maturity is estimated inside the ADMB template rather than read from the .dat,
# so it is absent from the report object. The fitted logistic is rebuilt from
# mat_beta1 and mat_beta2 in the parameter file and held fixed in SPoRC.
mat_beta1 <- get_par("mat_beta1")
mat_beta2 <- get_par("mat_beta2")
maa_vec <- 1 / (1 + exp(-(mat_beta1 + mat_beta2 * ages)))
MatAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
MatAA[1, 1, , 1, , 1] <- matrix(rep(maa_vec, each = n_yrs), nrow = n_yrs)

# The .dat carries the size at age matrix as [length bin x age] with each age
# column a distribution over length bins, so the columns are normalised.
size_age_raw <- read_matrix(121, n_lens, n_ages)
col_sums <- colSums(size_age_raw)
col_sums[col_sums == 0] <- 1
size_age_mat <- sweep(size_age_raw, 2, col_sums, "/")
SizeAgeTrans <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes))
for(y in seq_len(n_yrs)) SizeAgeTrans[1, 1, y, 1, , , 1] <- size_age_mat

# Ageing error is stored as [observed age x true age]. The ADMB template
# column-normalises it to P(obs | true); SPoRC multiplies predicted numbers at
# age on the right, so it needs [true x obs], and the observed plus group
# collapses rows 38 to 44 onto row 38.
age_error_raw <- read_matrix(76, n_ages, n_ages)
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

# Survey indices ---------------------------------------------------------------
# The .dat carries two survey blocks, biomass in tonnes and abundance in
# millions. The control file fits biomass and does not fit abundance, so only
# the biomass block is bridged. Fleet 1 is the Aleutian Islands bottom trawl
# survey, fleet 2 the eastern Bering Sea slope survey. The .dat reports a
# standard deviation on the arithmetic scale, carried to SPoRC as a lognormal
# coefficient of variation.
srv_nyrs <- read_nums(lines[39])
ai_srv_yrs <- read_nums(lines[42])
ebs_srv_yrs <- read_nums(lines[43])
ai_srv_obs <- read_nums(lines[45])
ebs_srv_obs <- read_nums(lines[46])
ai_srv_sd <- read_nums(lines[48])
ebs_srv_sd <- read_nums(lines[49])
stopifnot(length(ai_srv_yrs) == srv_nyrs[1], length(ebs_srv_yrs) == srv_nyrs[2])

ai_ind <- as.integer(years %in% ai_srv_yrs)
ebs_ind <- as.integer(years %in% ebs_srv_yrs)

ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
UseSrvIdx <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx[1, ai_ind == 1, 1, 1] <- ai_srv_obs
ObsSrvIdx_SE[1, ai_ind == 1, 1, 1] <- ai_srv_sd / ai_srv_obs
UseSrvIdx[1, ai_ind == 1, 1, 1] <- 1L
ObsSrvIdx[1, ebs_ind == 1, 1, 2] <- ebs_srv_obs
ObsSrvIdx_SE[1, ebs_ind == 1, 1, 2] <- ebs_srv_sd / ebs_srv_obs
UseSrvIdx[1, ebs_ind == 1, 1, 2] <- 1L

# Composition helpers ----------------------------------------------------------
# obs_mat is [n_bins x n_obs_years]; every helper fills a single fleet.
fill_obs <- function(arr, obs_mat, ind_vec, fleet) {
  iy <- 1L
  for(t in seq_len(n_yrs)) {
    if(ind_vec[t] == 1L) {
      arr[1, t, 1, , 1, fleet] <- obs_mat[, iy]
      iy <- iy + 1L
    }
  } # end t loop
  arr
}
fill_iss <- function(arr, iss_vec, ind_vec, fleet) {
  iy <- 1L
  for(t in seq_len(n_yrs)) {
    if(ind_vec[t] == 1L) {
      arr[1, t, 1, 1, fleet] <- iss_vec[iy]
      iy <- iy + 1L
    }
  } # end t loop
  arr
}

# Composition sample sizes and weights -----------------------------------------
# The ADMB objective uses applied_N = compweight * stage1_N, where the stage 1
# sample sizes are the square root sample size rows in pop24.ctl and the
# compweights live in a separate file the template opens itself. The multinomial
# nLL is linear in sample size, so ISS = stage 1 with Wt = compweight is
# equivalent to ISS = applied with Wt = 1; the former keeps the two sources
# visible and separable.
fish_age_stg1 <- read_nums(ctl[64])
fish_len_stg1 <- read_nums(ctl[70])
ai_age_stg1 <- read_nums(ctl[76])
ebs_age_stg1 <- read_nums(ctl[77])
ai_len_stg1 <- read_nums(ctl[84])

cw <- read_nums(readLines(src("compweights.ctl")))
stopifnot(length(cw) == 6)
wt_fish_age <- cw[1]
wt_fish_len <- cw[2]
wt_ai_age <- cw[3]
wt_ebs_age <- cw[4]
wt_ai_len <- cw[5]

# Compositions -----------------------------------------------------------------
fish_age_yrs <- read_nums(lines[149])
n_fish_age <- read_nums(lines[147])
fish_age_ind <- as.integer(years %in% fish_age_yrs)
stopifnot(length(fish_age_yrs) == n_fish_age, n_fish_age == length(fish_age_stg1))
fish_age_mat <- read_comp_stack(152, n_fish_age, n_obs_ages)
fish_age_mat <- t(fish_age_mat / rowSums(fish_age_mat))

# The 1964 to 1972 fishery length rows are raw counts and the 1977 onwards rows
# are proportions, so row normalising handles both.
fish_len_yrs <- read_nums(lines[198])
n_fish_len <- read_nums(lines[196])
fish_len_ind <- as.integer(years %in% fish_len_yrs)
stopifnot(length(fish_len_yrs) == n_fish_len, n_fish_len == length(fish_len_stg1))
fish_len_mat <- read_comp_stack(201, n_fish_len, n_lens)
fish_len_mat <- t(fish_len_mat / rowSums(fish_len_mat))

srv_age_nyrs <- read_nums(lines[265])
ai_age_yrs <- read_nums(lines[267])
ebs_age_yrs <- read_nums(lines[268])
n_ai_age <- srv_age_nyrs[1]
n_ebs_age <- srv_age_nyrs[2]
ai_age_ind <- as.integer(years %in% ai_age_yrs)
ebs_age_ind <- as.integer(years %in% ebs_age_yrs)
stopifnot(length(ai_age_yrs) == n_ai_age, n_ai_age == length(ai_age_stg1),
          length(ebs_age_yrs) == n_ebs_age, n_ebs_age == length(ebs_age_stg1))
ai_age_mat <- read_comp_stack(271, n_ai_age, n_obs_ages)
ai_age_mat <- t(ai_age_mat / rowSums(ai_age_mat))
ebs_age_mat <- read_comp_stack(298, n_ebs_age, n_obs_ages)
ebs_age_mat <- t(ebs_age_mat / rowSums(ebs_age_mat))

# The Aleutian Islands survey contributes a single year of length compositions
# and the eastern Bering Sea slope survey none.
srv_len_nyrs <- read_nums(lines[310])
ai_len_yrs <- read_nums(lines[312])
ai_len_ind <- as.integer(years %in% ai_len_yrs)
stopifnot(srv_len_nyrs[1] == 1, srv_len_nyrs[2] == 0)
ai_len_vec <- read_nums(lines[314])
ai_len_mat <- matrix(ai_len_vec / sum(ai_len_vec), ncol = 1)

ObsFishAgeComps <- fill_obs(array(NA_real_, c(n_regions, n_yrs, n_seas, n_obs_ages, n_sexes, n_fish_fleets)),
                            fish_age_mat, fish_age_ind, 1)
UseFishAgeComps <- array(0L, c(n_regions, n_yrs, n_seas, n_fish_fleets))
UseFishAgeComps[1, , 1, 1] <- fish_age_ind
ISS_FishAgeComps <- fill_iss(array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)),
                             fish_age_stg1, fish_age_ind, 1)
Wt_FishAgeComps <- array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
Wt_FishAgeComps[1, fish_age_ind == 1, 1, 1, 1] <- wt_fish_age

ObsFishLenComps <- fill_obs(array(NA_real_, c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets)),
                            fish_len_mat, fish_len_ind, 1)
UseFishLenComps <- array(0L, c(n_regions, n_yrs, n_seas, n_fish_fleets))
UseFishLenComps[1, , 1, 1] <- fish_len_ind
ISS_FishLenComps <- fill_iss(array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)),
                             fish_len_stg1, fish_len_ind, 1)
Wt_FishLenComps <- array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
Wt_FishLenComps[1, fish_len_ind == 1, 1, 1, 1] <- wt_fish_len

ObsSrvAgeComps <- array(NA_real_, c(n_regions, n_yrs, n_seas, n_obs_ages, n_sexes, n_srv_fleets))
ObsSrvAgeComps <- fill_obs(ObsSrvAgeComps, ai_age_mat, ai_age_ind, 1)
ObsSrvAgeComps <- fill_obs(ObsSrvAgeComps, ebs_age_mat, ebs_age_ind, 2)
UseSrvAgeComps <- array(0L, c(n_regions, n_yrs, n_seas, n_srv_fleets))
UseSrvAgeComps[1, ai_age_ind == 1, 1, 1] <- 1L
UseSrvAgeComps[1, ebs_age_ind == 1, 1, 2] <- 1L
ISS_SrvAgeComps <- array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))
ISS_SrvAgeComps <- fill_iss(ISS_SrvAgeComps, ai_age_stg1, ai_age_ind, 1)
ISS_SrvAgeComps <- fill_iss(ISS_SrvAgeComps, ebs_age_stg1, ebs_age_ind, 2)
Wt_SrvAgeComps <- array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))
Wt_SrvAgeComps[1, ai_age_ind == 1, 1, 1, 1] <- wt_ai_age
Wt_SrvAgeComps[1, ebs_age_ind == 1, 1, 1, 2] <- wt_ebs_age

ObsSrvLenComps <- array(NA_real_, c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets))
ObsSrvLenComps <- fill_obs(ObsSrvLenComps, ai_len_mat, ai_len_ind, 1)
UseSrvLenComps <- array(0L, c(n_regions, n_yrs, n_seas, n_srv_fleets))
UseSrvLenComps[1, ai_len_ind == 1, 1, 1] <- 1L
ISS_SrvLenComps <- array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))
ISS_SrvLenComps <- fill_iss(ISS_SrvLenComps, ai_len_stg1, ai_len_ind, 1)
Wt_SrvLenComps <- array(0, c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))
Wt_SrvLenComps[1, ai_len_ind == 1, 1, 1, 1] <- wt_ai_len

stopifnot(
  all(abs(colSums(fish_age_mat) - 1) < 1e-6),
  all(abs(colSums(fish_len_mat) - 1) < 1e-6),
  all(abs(colSums(ai_age_mat) - 1) < 1e-6),
  all(abs(colSums(ebs_age_mat) - 1) < 1e-6),
  all(abs(colSums(ai_len_mat) - 1) < 1e-6),
  all(abs(rowSums(AgeingError_mat) - 1) < 1e-6)
)

# ADMB maximum likelihood estimate ---------------------------------------------
# Read at full precision from the parameter file. rec_dev covers 62 years, the
# last three taking mean recruitment. historic_F is fixed at phase -1 and sets
# the first year equilibrium; log_avg_fmort is estimated and sets the F series,
# so the two are independent and SPoRC carries historic_F in its own ln_init_F
# slot. fsh_sel_par is the 5 by 5 bicubic node grid, written row major with the
# rows as age nodes.
mle <- list(
  M = exp(get_par("log_avg_M")),
  mean_log_rec = get_par("mean_log_rec"),
  log_rinit = get_par("log_rinit"),
  log_avg_fmort = get_par("log_avg_fmort"),
  historic_F = get_par("historic_F"),
  q_srv = exp(c(get_par("log_q_srv\\[1\\]"), get_par("log_q_srv\\[2\\]"))),
  sel_a50_srv = c(get_par("sel_a50_srv\\[1\\]"), get_par("sel_a50_srv\\[2\\]")),
  sel_aslope_srv = c(get_par("sel_aslope_srv\\[1\\]"), get_par("sel_aslope_srv\\[2\\]")),
  fsh_sel_par = get_par("fsh_sel_par"),
  fmort_dev = get_par("fmort_dev"),
  rec_dev = get_par("rec_dev"),
  mat_beta1 = mat_beta1,
  mat_beta2 = mat_beta2
)
stopifnot(length(mle$fmort_dev) == n_yrs,
          length(mle$rec_dev) == n_yrs - fixedrec,
          length(mle$fsh_sel_par) == fsh_yr_nodes * fsh_age_nodes)

# ADMB output ------------------------------------------------------------------
# The comparison target. natage is reported over the 38 observed ages with the
# model's ages 41 to 46 pooled into the final column. selfish and fmort are both
# rescaled for reporting: the template divides selectivity by its within-year
# maximum and multiplies F by the same factor, leaving F times selectivity
# invariant, so the raw quantities have to be recovered before comparing.
admb <- list(
  NAA = as.matrix(mod$natage),
  SSB = mod$t.series$spbiom,
  Rec = mod$t.series$a3recs,
  Fmort = mod$t.series$fmort,
  TotBiom = mod$t.series$totbiom,
  yrs = mod$t.series$year,
  sel_fsh = as.matrix(mod$selfish),
  sel_srv = mod$srvsel,
  datalikecomp = mod$datalikecomp,
  pen_likecomp = mod$pen_likecomp
)

sgl_rg_bsai_pop_data <- list(
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
  nselages = nselages,

  WAA = WAA,
  WAA_fish = WAA_fish,
  MatAA = MatAA,
  SizeAgeTrans = SizeAgeTrans,
  AgeingError = AgeingError,
  pop_waa = pop_waa,

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

  yrs_srv_ai = ai_srv_yrs,
  yrs_srv_ebs = ebs_srv_yrs,
  yrs_fish_age = fish_age_yrs,
  yrs_fish_len = fish_len_yrs,
  yrs_srv_age_ai = ai_age_yrs,
  yrs_srv_age_ebs = ebs_age_yrs,
  yrs_srv_len_ai = ai_len_yrs,

  sigmaR = sigmaR,
  fixedrec = fixedrec,
  spawn_mo = spawn_mo,
  mean_M = mean_M,
  cv_M = cv_M,
  mean_q = mean_q,
  cv_q = cv_q,

  fsh_sel_styr = fsh_sel_styr,
  fsh_yr_nodes = fsh_yr_nodes,
  fsh_age_nodes = fsh_age_nodes,

  lam_rec = lam_rec,
  catch_wt = lam_catch,
  fmort_wt = 0.1,
  sel_pen_wts = c(
    smooth_dome = lam_dome,
    smooth_bin_curve = lam_bin_curve,
    smooth_yr_diff = lam_yr_diff,
    smooth_yr_curve = lam_yr_curve,
    smooth_mean_center = lam_mean_ctr
  ),

  mle = mle,
  admb = admb
)

usethis::use_data(sgl_rg_bsai_pop_data, internal = FALSE, overwrite = TRUE)
