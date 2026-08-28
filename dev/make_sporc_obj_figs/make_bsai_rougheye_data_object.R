# Purpose: Build sgl_rg_rebs_data, the container for the BSAI blackspotted and
#          rougheye rockfish case study. Everything the case study and its
#          regression tests need lives in the object: model inputs, the ADMB
#          maximum likelihood estimate used as a seed, and the ADMB output the
#          bridge is compared against.
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26
#
# Sources, all under dev/dev_data:
#   rougheye24.dat            2024 assessment input file
#   m_24_1_reweighted.rdat    ADMB report object, the comparison target
#
# The maximum likelihood estimates below are transcribed from the
# m_24_1_reweighted ADMB parameter file, which is not redistributed here. They
# are the seed for the bridge test, which evaluates SPoRC at the assessment's
# own optimum without optimizing.

library(here)

dat_path <- here("dev", "dev_data", "rougheye24.dat")
lines <- readLines(dat_path)
mod <- dget(here("dev", "dev_data", "m_24_1_reweighted.rdat"))

read_nums <- function(lns) {
  txt <- gsub("#.*", "", lns)
  vals <- as.numeric(unlist(strsplit(trimws(paste(txt, collapse = " ")), "\\s+")))
  vals[!is.na(vals)]
}

# Dimensions -----------------------------------------------------------------
# The model carries 52 age classes (3 to 54) but the compositions are reported
# over 43 (3 to 45), so the ageing error matrix maps the model ages onto the
# shorter observed range.
n_pop <- 1
n_regions <- 1
n_seas <- 1
n_sexes <- 1
n_fish_fleets <- 1
n_srv_fleets <- 1

years <- 1977:2024
ages <- 3:54
obs_ages <- 3:45
lens <- 12:50

n_yrs <- length(years)
n_ages <- length(ages)
n_obs_ages <- length(obs_ages)
n_lens <- length(lens)

# Priors ---------------------------------------------------------------------
mean_q <- 1.0
cv_q <- 0.05
mean_M <- 0.045
cv_M <- 0.05
sigmaR <- 0.75

# Biologicals ----------------------------------------------------------------
# Weight at age is time invariant and is carried in the .dat in grams.
waa_vec <- read_nums(lines[277])
stopifnot(length(waa_vec) == n_ages)
WAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
WAA[1, 1, , 1, , 1] <- matrix(rep(waa_vec, each = n_yrs), nrow = n_yrs)

# Maturity is estimated inside the ADMB template from the binomial maturity
# data rather than read from the .dat, so the fitted ogive is transcribed from
# the assessment report and held fixed in SPoRC.
maa_vec <- c(
  0.00340888, 0.0044329,  0.00576276, 0.00748858, 0.00972618, 0.0126239,
  0.0163706,  0.0212055,  0.0274284,  0.0354115,  0.0456091,  0.0585651,
  0.0749125,  0.0953608,  0.120663,   0.151553,   0.188655,   0.232353,
  0.282646,   0.339018,   0.400358,   0.464989,   0.530817,   0.59559,
  0.657196,   0.713924,   0.764628,   0.808752,   0.846267,   0.877538,
  0.903175,   0.923911,   0.940499,   0.953652,   0.964008,   0.972118,
  0.978442,   0.983356,   0.987164,   0.99011,    0.992385,   0.99414,
  0.995492,   0.996533,   0.997335,   0.997951,   0.998425,   0.99879,
  0.99907,    0.999286,   0.999286,   0.999789
)
stopifnot(length(maa_vec) == n_ages)
MatAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
MatAA[1, 1, , 1, , 1] <- matrix(rep(maa_vec, each = n_yrs), nrow = n_yrs)

# The .dat carries the size at age matrix as [length bin x age] with each age
# column a distribution over length bins, so the columns are normalized.
sizeage_raw <- matrix(read_nums(lines[122:(122 + n_lens - 1)]),
                      nrow = n_lens, ncol = n_ages, byrow = TRUE)
col_sums <- colSums(sizeage_raw)
col_sums[col_sums == 0] <- 1
size_age_mat <- sweep(sizeage_raw, 2, col_sums, "/")
SizeAgeTrans <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes))
for(y in seq_len(n_yrs)) SizeAgeTrans[1, 1, y, 1, , , 1] <- size_age_mat

# Ageing error is stored as [observed age x true age]. The ADMB template
# column-normalizes it to P(obs | true); SPoRC multiplies predicted numbers at
# age on the right, so it needs [true x obs], and the observed plus group
# collapses rows 43 to 52 onto row 43.
age_error_raw <- matrix(read_nums(lines[69:(69 + n_ages - 1)]),
                        nrow = n_ages, ncol = n_ages, byrow = TRUE)
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

# Catch ----------------------------------------------------------------------
catch_vec <- c(
  155, 2423, 3077, 660, 595, 189, 58, 35, 10, 21, 79, 75, 381, 1619, 137,
  1181, 924, 749, 395, 816, 954, 526, 385, 280, 550, 273, 174, 185, 78, 197,
  157, 171, 184, 201, 131, 182, 303, 179, 159, 121, 191, 232, 345, 477, 412,
  341, 469, 453
)
stopifnot(length(catch_vec) == n_yrs)
ObsCatch <- array(catch_vec, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
UseCatch <- array(1L, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))

# Survey index ---------------------------------------------------------------
# The Aleutian Islands bottom trawl survey is the only index. The .dat reports
# a standard deviation on the arithmetic scale, which is carried to SPoRC as a
# lognormal coefficient of variation.
srv_yrs <- c(1991, 1994, 1997, 2000, 2002, 2004, 2006, 2010, 2012, 2014,
             2016, 2018, 2022, 2024)
srv_obs <- c(10637.8, 13414.8, 10905.1, 14217.7, 8422.9, 14385.3, 8281.9,
             8541.3, 12401.2, 4425.2, 9468.7, 9514.6, 15682.5, 24087.1)
srv_sd <- c(4957.313, 3742.261, 2431.213, 3241.826, 1793.523, 3781.592,
            2101.659, 2255.739, 4772.594, 848.230, 2471.500, 4498.506,
            7785.548, 10060.479)
srv_ind <- as.integer(years %in% srv_yrs)

ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
UseSrvIdx <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx[1, srv_ind == 1, 1, 1] <- srv_obs
ObsSrvIdx_SE[1, srv_ind == 1, 1, 1] <- srv_sd / srv_obs
UseSrvIdx[1, srv_ind == 1, 1, 1] <- 1L

# Composition helpers --------------------------------------------------------
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
# Comp rows in the .dat are interleaved with single "#year" comment lines, so
# the data sit on every other line.
read_comp_stack <- function(start, n) {
  out <- vector("list", n)
  dl <- start
  for(i in seq_len(n)) {
    out[[i]] <- read_nums(lines[dl])
    dl <- dl + 2L
  } # end i loop
  m <- do.call(rbind, out)
  t(m / rowSums(m))
}

# Compositions ---------------------------------------------------------------
# The composition weights are the stage-2 multipliers from compweights1.ctl of
# the reweighted model, applied in SPoRC as likelihood weights.
fish_age_yrs <- c(2004, 2005, 2007, 2008, 2009, 2011, 2013, 2015,
                  2017, 2019, 2020, 2021, 2023)
fish_age_ind <- as.integer(years %in% fish_age_yrs)
fish_age_iss <- c(90, 65, 83, 74, 90, 85, 160, 126, 156, 201, 224, 397, 422)
fish_age_mat <- read_comp_stack(166L, length(fish_age_yrs))
ObsFishAgeComps <- make_obs_comp(fish_age_mat, fish_age_ind, n_obs_ages, n_fish_fleets)
UseFishAgeComps <- make_use_comp(fish_age_ind, n_fish_fleets)
ISS_FishAgeComps <- make_iss_comp(fish_age_iss, fish_age_ind, n_fish_fleets)
Wt_FishAgeComps <- make_wt_comp(fish_age_ind, 0.0423816, n_fish_fleets)

fish_len_yrs <- c(1979, 1990, 1992, 1993, 2003, 2010, 2012, 2014, 2016, 2018, 2022)
fish_len_ind <- as.integer(years %in% fish_len_yrs)
fish_len_iss <- c(93, 20, 67, 39, 100, 375, 164, 213, 130, 331, 621)
fish_len_mat <- read_comp_stack(197L, length(fish_len_yrs))
ObsFishLenComps <- make_obs_comp(fish_len_mat, fish_len_ind, n_lens, n_fish_fleets)
UseFishLenComps <- make_use_comp(fish_len_ind, n_fish_fleets)
ISS_FishLenComps <- make_iss_comp(fish_len_iss, fish_len_ind, n_fish_fleets)
Wt_FishLenComps <- make_wt_comp(fish_len_ind, 0.480182, n_fish_fleets)

srv_age_yrs <- c(1991, 1994, 1997, 2000, 2002, 2004, 2006, 2010, 2012,
                 2014, 2016, 2018, 2022)
srv_age_ind <- as.integer(years %in% srv_age_yrs)
srv_age_iss <- c(23, 55, 83, 71, 66, 83, 76, 68, 84, 68, 87, 89, 120)
srv_age_mat <- read_comp_stack(228L, length(srv_age_yrs))
ObsSrvAgeComps <- make_obs_comp(srv_age_mat, srv_age_ind, n_obs_ages, n_srv_fleets)
UseSrvAgeComps <- make_use_comp(srv_age_ind, n_srv_fleets)
ISS_SrvAgeComps <- make_iss_comp(srv_age_iss, srv_age_ind, n_srv_fleets)
Wt_SrvAgeComps <- make_wt_comp(srv_age_ind, 0.251453, n_srv_fleets)

# The survey contributes length compositions in 2024 only.
srv_len_yrs <- 2024
srv_len_ind <- as.integer(years %in% srv_len_yrs)
srv_len_vec <- read_nums(lines[272])
srv_len_mat <- array(srv_len_vec / sum(srv_len_vec), dim = c(n_lens, 1))
ObsSrvLenComps <- make_obs_comp(srv_len_mat, srv_len_ind, n_lens, n_srv_fleets)
UseSrvLenComps <- make_use_comp(srv_len_ind, n_srv_fleets)
ISS_SrvLenComps <- make_iss_comp(102, srv_len_ind, n_srv_fleets)
Wt_SrvLenComps <- make_wt_comp(srv_len_ind, 0.251453, n_srv_fleets)

# ADMB maximum likelihood estimate -------------------------------------------
# Transcribed from the m_24_1_reweighted parameter file. rec_dev covers
# 1977 to 2021, the last three years taking mean recruitment; fydev carries the
# 42 initial age deviations, whose last value is shared by ages beyond the
# observed range.
mle <- list(
  M = 0.0500461712767,
  mean_log_rec = 0.322819748196,
  log_rinit = 0.279306493366,
  sel_aslope_fish = 0.631073905734,
  sel_a50_fish = 13.7067447267,
  sel_aslope_srv = 0.363751750329,
  sel_a50_srv = 15.6702662974,
  q_srv = 1.04415928660,
  log_avg_fmort = -3.71801297027,
  fmort_dev = c(
    -0.925828865129, 1.86550835452, 2.24719001278, 0.800368702796,
    0.705437966135, -0.447630916116, -1.65355359851, -2.18798744764,
    -3.46757887788, -2.75410472118, -1.45397498933, -1.52338497107,
    0.0929633073861, 1.58546767953, -0.833042086020, 1.34822351422,
    1.16934855826, 1.01283587231, 0.406109107579, 1.16995949988,
    1.39722149375, 0.861989982042, 0.585537142040, 0.294732063690,
    1.01635358090, 0.352512765671, -0.0787199309256, -0.000329546720982,
    -0.857964248745, 0.0797262076344, -0.140531397793, -0.0540104778028,
    0.0137727301986, 0.0885524539813, -0.371336372293, -0.0828582241821,
    0.388368692135, -0.196254465443, -0.379051756169, -0.717841567445,
    -0.323326511628, -0.183932519751, 0.163996835406, 0.443891746827,
    0.246108769304, -0.0122580394508, 0.220645559001, 0.0886789332426
  ),
  rec_dev = c(
    -0.0904859091640, -0.138741063984, -0.178203395920, -0.211684067338,
    -0.227507488057, -0.216711065071, -0.195681305383, -0.194097183965,
    -0.255914291937, -0.383890892956, -0.539620819849, -0.685452317016,
    -0.813249853871, -0.924926810566, -1.00993523778, -1.05373008752,
    -1.04638544493, -0.996624281592, -0.906310273195, -0.791927380800,
    -0.650717323419, -0.483845006620, -0.304766482936, -0.0245869917680,
    0.455175175123, 0.563053725847, 0.495867856346, 0.588223025050,
    0.862829994440, 0.464652124822, 0.286710065419, 0.343967919990,
    0.508969225128, 0.566912130379, 0.626396028829, 0.491285564135,
    0.331309163550, 3.00673727971, 0.225423910406, 0.415681848590,
    0.653754717105, 0.769266197507, 0.494967364147, 0.162617104448,
    0.0111945546718
  ),
  fydev = c(
    -0.0291191871166, 0.0395015505620, 0.114145911257, 0.186162806402,
    0.252753608351, 0.297765494303, 0.351725661071, 0.414548313026,
    0.485634452595, 0.512336774426, 0.535063540015, 0.505448760443,
    0.434934028997, 0.342834893935, 0.263480131673, 0.183933528324,
    0.112993162060, 0.0554465267841, 0.000579894203366, -0.0450019153374,
    -0.0829461029226, -0.115318904457, -0.140173613641, -0.158239270435,
    -0.176236957235, -0.188088349514, -0.195888516170, -0.200341057984,
    -0.202626691571, -0.203773742269, -0.203890490733, -0.202423936324,
    -0.199938964188, -0.197655501912, -0.196437138530, -0.196392885858,
    -0.191097422325, -0.187211293753, -0.181490934661, -0.175849181328,
    -0.170551980271, -1.24859499989
  )
)

# ADMB output ----------------------------------------------------------------
# The comparison target. natage is reported over the 43 observed ages with the
# model's ages 45 to 54 pooled into the final column.
admb <- list(
  NAA = mod$natage,
  SSB = mod$t.series$spbiom,
  Rec = mod$t.series$a3recs,
  Fmort = mod$t.series$fmort,
  TotBiom = mod$t.series$totbiom,
  pred_srv = mod$t.series$AI_survey_pred_bio,
  sel_fsh = mod$selfish,
  sel_srv = mod$srv_sel$AI_survey_sel,
  datalikecomp = mod$datalikecomp,
  pen_likecomp = mod$pen_likecomp,
  controlrule = mod$controlrule
)

sgl_rg_rebs_data <- list(
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
  yrs_srv_len = srv_len_yrs,

  sigmaR = sigmaR,
  mean_M = mean_M,
  cv_M = cv_M,
  mean_q = mean_q,
  cv_q = cv_q,

  mle = mle,
  admb = admb
)

usethis::use_data(sgl_rg_rebs_data, internal = FALSE, overwrite = TRUE)
