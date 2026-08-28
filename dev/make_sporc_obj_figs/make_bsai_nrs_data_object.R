# Purpose: Build sgl_rg_bsai_nrs_data, the container for the BSAI northern rock
#          sole case study. Everything the case study and its regression test
#          need lives in the object: model inputs, the flatfish model's maximum
#          likelihood estimate used as a seed, and the assessment output the
#          bridge is compared against.
# Creator: Matthew LH. Cheng
# Date Created: 8/21/26
#
# Sources, all under dev/dev_data/nrs_m24.2_2024:
#   c3.dat          2024 assessment input file for Model 24.2
#   mod.ctl         control file, carries the weights, sigmas and phases
#   fm.par          parameter file, read at full precision
#   fm.std          reported time series with standard errors
#   fm.rep          nLogPosterior, the likelihood components at the estimate
#   fm_legacy.rep   predicted catch
#   fm.tpl          the template, for reference
#
# Northern rock sole differs from the other single region case studies in four
# ways the object has to carry: the model is two sex with sex specific natural
# mortality, the initial numbers at age are free and estimated separately for
# each sex about one shared level, fishery selectivity is logistic in both
# parameters of both sexes with annual deviations and the male curve scaled by a
# constant, and the survey is read twice, once as a July biomass index and once
# as January 1 compositions.

library(here)
source(here("dev", "make_sporc_obj_figs", "helper-fm.R"))

src <- function(f) here("dev", "dev_data", "nrs_m24.2_2024", f)
dat <- read_fm_dat(src("c3.dat"))
ctl <- read_fm_ctl(src("mod.ctl"))
fpar <- read_fm_par(src("fm.par"))
fstd <- read_fm_std(src("fm.std"))
frep <- read_fm_rep(src("fm.rep"))
pred_catch <- read_fm_legacy_catch(src("fm_legacy.rep"))

yrs <- dat$years
n_yrs <- length(yrs)
n_ages <- dat$n_ages
n_sexes <- 2
n_srv <- 2                                   # index fleet and composition fleet

# Guard: the reported series have to line up with the model years, and the
# reported components have to rebuild the reported objective, before anything
# downstream is entitled to trust the parse. fm.tpl scales catch by lambda(3)
# inside the objective but reports it unscaled.
stopifnot(length(fstd$SSB$value) == n_yrs, length(fstd$pred_rec$value) == n_yrs,
          length(fpar$rec_dev) == n_yrs, length(fpar$fmort_dev) == n_yrs)
comp_total <- sum(unlist(frep)) + frep$catch * (ctl$lambda[3] - 1)
# nLogPosterior prints six significant digits, so the largest component alone
# carries five thousandths; the check is at the precision the file was written
stopifnot(abs(comp_total - fpar$objective) < 0.01)

i_fsh_age <- match(dat$yrs_fsh_age, yrs)
i_srv <- match(dat$yrs_srv, yrs)
i_srv_age <- match(dat$yrs_srv_age, yrs)

# Model inputs ---------------------------------------------------------------
ObsCatch <- array(dat$obs_catch, dim = c(1, n_yrs, 1, 1))
UseCatch <- array(1, dim = c(1, n_yrs, 1, 1))

# Compositions are joint across sexes, females then males in the file, and are
# normalized the way fm.tpl normalizes them before it builds its offset.
oac_fsh <- dat$oac_fsh / rowSums(dat$oac_fsh)
oac_srv <- dat$oac_srv / rowSums(dat$oac_srv)

ObsFishAgeComps <- array(NA_real_, dim = c(1, n_yrs, 1, n_ages, n_sexes, 1))
ObsFishAgeComps[1, i_fsh_age, 1, , 1, 1] <- oac_fsh[, 1:n_ages]
ObsFishAgeComps[1, i_fsh_age, 1, , 2, 1] <- oac_fsh[, n_ages + 1:n_ages]
UseFishAgeComps <- array(0, dim = c(1, n_yrs, 1, 1))
UseFishAgeComps[1, i_fsh_age, 1, 1] <- 1
ISS_FishAgeComps <- array(0, dim = c(1, n_yrs, 1, n_sexes, 1))
ISS_FishAgeComps[1, i_fsh_age, 1, , 1] <- dat$iss_fsh_age   # a joint composition reads the first sex's slot

# The survey index sits on fleet 1 and the compositions on fleet 2, so the two
# can be read at their own times off one shared selectivity curve.
ObsSrvIdx <- ObsSrvIdx_SE <- array(NA_real_, dim = c(1, n_yrs, 1, n_srv))
UseSrvIdx <- array(0, dim = c(1, n_yrs, 1, n_srv))
ObsSrvIdx[1, i_srv, 1, 1] <- dat$obs_srv
ObsSrvIdx_SE[1, i_srv, 1, 1] <- sqrt(log((dat$obs_se_srv / dat$obs_srv)^2 + 1))
UseSrvIdx[1, i_srv, 1, 1] <- 1

ObsSrvAgeComps <- array(NA_real_, dim = c(1, n_yrs, 1, n_ages, n_sexes, n_srv))
ObsSrvAgeComps[1, i_srv_age, 1, , 1, 2] <- oac_srv[, 1:n_ages]
ObsSrvAgeComps[1, i_srv_age, 1, , 2, 2] <- oac_srv[, n_ages + 1:n_ages]
UseSrvAgeComps <- array(0, dim = c(1, n_yrs, 1, n_srv))
UseSrvAgeComps[1, i_srv_age, 1, 2] <- 1
ISS_SrvAgeComps <- array(0, dim = c(1, n_yrs, 1, n_sexes, n_srv))
ISS_SrvAgeComps[1, i_srv_age, 1, , 2] <- dat$iss_srv_age

# Weights and maturity are year specific empirical inputs. Spawning biomass is
# female only, so the male slice of maturity is never read.
WAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes))
WAA[1, 1, , 1, , 1] <- dat$wt_pop_f
WAA[1, 1, , 1, , 2] <- dat$wt_pop_m
WAA_fish <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 1))
WAA_fish[1, 1, , 1, , 1, 1] <- dat$wt_fsh[, 1:n_ages]
WAA_fish[1, 1, , 1, , 2, 1] <- dat$wt_fsh[, n_ages + 1:n_ages]
WAA_srv <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, n_srv))
for(sf in 1:n_srv) {
  WAA_srv[1, 1, , 1, , 1, sf] <- dat$wt_srv_f
  WAA_srv[1, 1, , 1, , 2, sf] <- dat$wt_srv_m
} # end sf loop
MatAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes))
MatAA[1, 1, , 1, , 1] <- dat$maturity

# The multinomial offsets, recomputed the way fm.tpl builds them, so the case
# study can show its composition likelihoods on the assessment's own scale.
offset_fsh <- -sum(dat$iss_fsh_age * rowSums((oac_fsh + 1e-3) * log(oac_fsh + 1e-3)))
offset_srv <- -sum(dat$iss_srv_age * rowSums((oac_srv + 1e-3) * log(oac_srv + 1e-3)))

sgl_rg_bsai_nrs_data <- list(

  # Dimensions
  years = yrs, ages = 1:n_ages, lens = NA,
  n_regions = 1, n_sexes = n_sexes, n_fish_fleets = 1, n_srv_fleets = n_srv,
  n_seas = 1, n_pop = 1, natal_region = 1,
  yrs_fsh_age = dat$yrs_fsh_age, yrs_srv = dat$yrs_srv, yrs_srv_age = dat$yrs_srv_age,

  # Observations
  ObsCatch = ObsCatch, UseCatch = UseCatch,
  ObsFishAgeComps = ObsFishAgeComps, UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,
  ObsSrvIdx = ObsSrvIdx, ObsSrvIdx_SE = ObsSrvIdx_SE, UseSrvIdx = UseSrvIdx,
  ObsSrvAgeComps = ObsSrvAgeComps, UseSrvAgeComps = UseSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,

  # Biologicals
  WAA = WAA, WAA_fish = WAA_fish, WAA_srv = WAA_srv, MatAA = MatAA,
  AgeingError = diag(n_ages),
  t_spawn = (dat$spawnmo - 1) / 12,
  t_srv = c((dat$srv_mo - 1) / 12, 0),         # July index, January 1 compositions

  # Specification carried by the control file
  nselages = ctl$nselages,
  sigmaC = 1 / sqrt(2 * ctl$lambda[3]),        # lambda(3) as a lognormal sigma
  sigmaR = 1 / sqrt(2 * ctl$lambda[1]),
  a50_sigma = ctl$a50_sigma, slp_sigma = ctl$slp_sigma,
  sr_pen_yrs = ctl$styr_sr:ctl$endyr_sr, sr_pen_sigma = ctl$sigmaR_exp,
  q_prior = list(mu = ctl$q_exp, sd = ctl$q_sigma / sqrt(dat$n_srv_yrs)),  # a vector prior over the survey years
  m_prior = list(mu = ctl$m_exp, sd = ctl$m_sigma),
  offset_fsh = offset_fsh, offset_srv = offset_srv,

  # The assessment's maximum likelihood estimate, the case study's seed
  mle = fpar,

  # The assessment's own output, the comparison target
  fm = list(
    SSB = fstd$SSB$value, SSB_sd = fstd$SSB$sd,
    Rec = fstd$pred_rec$value, Rec_sd = fstd$pred_rec$sd,
    TotBiom = fstd$TotBiom$value, TotBiom_sd = fstd$TotBiom$sd,
    pred_catch = pred_catch, Like_Comp = frep,
    objective = fpar$objective, max_grad = fpar$max_grad, n_par = fpar$n_par
  )
)

usethis::use_data(sgl_rg_bsai_nrs_data, internal = FALSE, overwrite = TRUE)
