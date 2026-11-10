# Purpose: Build sgl_rg_goa_pollock_data, the inputs and Rceattle estimates of GOA pollock models 26a and 26b for the case study
# Creator: Matthew LH. Cheng
# Date Created: 11/9/26
#
# Sources, all under dev/goa_pollock_rceattle_bridge/m26a:
#   rceattle_goa_pollock_26a_2026_fit.rds   Rceattle fit of model 26a, from 01_fit_rceattle_26.R
#   rceattle_goa_pollock_26b_2026_fit.rds   Rceattle fit of model 26b, the same script with MODEL_NAME=26b

library(here)

m26_dir <- here("dev", "goa_pollock_rceattle_bridge", "m26a")
rce_fits <- list(
  m26a = readRDS(file.path(m26_dir, "rceattle_goa_pollock_26a_2026_fit.rds")),
  m26b = readRDS(file.path(m26_dir, "rceattle_goa_pollock_26b_2026_fit.rds"))
)

# the observations and biology are the same in both models, so they are read from 26a
rce_data <- rce_fits$m26a$tmb_data # what the rceattle template was handed
rce_rep <- rce_fits$m26a$quantities # report at the joint mode

# Dimensions -----------------------------------------------------------------
years <- rce_data$styr:rce_data$endyr # 1970-2026
n_yrs <- length(years)
ages <- 1:10 # age 10 is the plus group
n_ages <- length(ages)
rce_srv_flt <- c(1, 2, 3, 6) # shelikof, nmfs bottom trawl, adfg bottom trawl, summer acoustic; 4, 5 and 7 are off
n_srv <- length(rce_srv_flt)
fish_flt <- 8 # rceattle fleet code of the fishery
comp_offset <- rce_data$comp_offset # 1e-5, added to observed and expected comps
srv_names <- c("Shelikof acoustic", "NMFS bottom trawl", "ADF&G bottom trawl", "Summer acoustic")

# Biology --------------------------------------------------------------------
# weight_hat rows: 1 population, 2 spawning, 2 + fleet code for each fleet
natmort <- array(0, dim = c(1, 1, n_yrs, n_ages, 1))
WAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, 1))
MatAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, 1))
WAA_fish <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
WAA_srv <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, 1, n_srv))

for(y in 1:n_yrs) {

  natmort[1,1,y,,1] <- rce_rep$M_at_age[1,1,,y] # fixed at age, the same every year
  WAA[1,1,y,1,,1] <- rce_rep$weight_hat[2,1,,y] # spawning weight
  MatAA[1,1,y,1,,1] <- rce_data$maturity[1,1:n_ages] # no sex ratio, SPoRC halves a one sex SSB itself
  WAA_fish[1,1,y,1,,1,1] <- rce_rep$weight_hat[2 + fish_flt,1,,y] # fishery weight

  for(sf in 1:n_srv) WAA_srv[1,1,y,1,,1,sf] <- rce_rep$weight_hat[2 + rce_srv_flt[sf],1,,y] # each survey's weight

} # end y loop

# Ageing error with the young ages pooled --------------------------------------
# rceattle folds the bins below Comp_accum_young into that bin after ageing error
age_error <- rce_data$age_error[1,,]
pool_fish <- diag(n_ages)
pool_fish[1,] <- 0
pool_fish[1,2] <- 1 # fishery: observed age 1 goes into age 2
pool_shel <- diag(n_ages)
pool_shel[1:2,] <- 0
pool_shel[1:2,3] <- 1 # shelikof: observed ages 1 and 2 go into age 3

AgeingError_fish <- array(0, dim = c(n_ages, n_ages, 1))
AgeingError_fish[,,1] <- age_error %*% pool_fish
AgeingError_srv <- array(0, dim = c(n_ages, n_ages, n_srv))
AgeingError_srv[,,1] <- age_error %*% pool_shel
AgeingError_srv[,,2] <- age_error
AgeingError_srv[,,3] <- age_error
AgeingError_srv[,,4] <- age_error

n_bins_fish <- 9 # bins left after pooling
n_bins_srv <- c(8, 10, 10, 10)

# Catch ----------------------------------------------------------------------
ObsCatch <- array(NA, dim = c(1, n_yrs, 1, 1))
UseCatch <- array(0, dim = c(1, n_yrs, 1, 1))

for(i in 1:nrow(rce_data$catch_ctl)) {

  y <- rce_data$catch_ctl[i,3] - rce_data$styr + 1
  if(y > n_yrs) next # projection rows
  ObsCatch[1,y,1,1] <- rce_data$catch_obs[i,1]
  UseCatch[1,y,1,1] <- 1

} # end i loop

# Survey indices -------------------------------------------------------------
ObsSrvIdx <- array(NA, dim = c(1, n_yrs, 1, n_srv))
ObsSrvIdx_SE <- array(NA, dim = c(1, n_yrs, 1, n_srv))
UseSrvIdx <- array(0, dim = c(1, n_yrs, 1, n_srv))
t_srv <- array(0, dim = c(1, 1, n_srv))

for(i in 1:nrow(rce_data$index_ctl)) {

  sf <- match(rce_data$index_ctl[i,1], rce_srv_flt)
  y <- rce_data$index_ctl[i,3] - rce_data$styr + 1
  if(is.na(sf)) next # a fleet that is off

  ObsSrvIdx[1,y,1,sf] <- rce_data$index_obs[i,1]
  ObsSrvIdx_SE[1,y,1,sf] <- rce_data$index_obs[i,2] # log scale sd
  UseSrvIdx[1,y,1,sf] <- 1
  t_srv[1,1,sf] <- rce_data$index_n[i,1] / 12 # survey month as a fraction of the year

} # end i loop

# Age compositions -----------------------------------------------------------
# observed comps go in already pooled, and the input sample size is scaled by the sum of
# obs + offset over the fitted bins, which is how rceattle counts them
ObsFishAgeComps <- array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1))
UseFishAgeComps <- array(0, dim = c(1, n_yrs, 1, 1))
ISS_FishAgeComps <- array(0, dim = c(1, n_yrs, 1, 1, 1))
ObsSrvAgeComps <- array(0, dim = c(1, n_yrs, 1, n_ages, 1, n_srv))
UseSrvAgeComps <- array(0, dim = c(1, n_yrs, 1, n_srv))
ISS_SrvAgeComps <- array(0, dim = c(1, n_yrs, 1, 1, n_srv))

for(i in 1:nrow(rce_data$comp_ctl)) {

  flt <- rce_data$comp_ctl[i,1]
  y <- rce_data$comp_ctl[i,5] - rce_data$styr + 1
  n_samp <- rce_data$comp_n[i,2]
  obs <- rce_data$comp_obs[i,1:n_ages]

  if(flt == fish_flt) {

    ObsFishAgeComps[1,y,1,,1,1] <- as.vector(obs %*% pool_fish)
    UseFishAgeComps[1,y,1,1] <- 1
    ISS_FishAgeComps[1,y,1,1,1] <- n_samp * (1 + n_bins_fish * comp_offset)

  } # end if fishery

  sf <- match(flt, rce_srv_flt)

  if(!is.na(sf)) {

    if(sf == 1) obs <- as.vector(obs %*% pool_shel)
    ObsSrvAgeComps[1,y,1,,1,sf] <- obs
    UseSrvAgeComps[1,y,1,sf] <- 1
    ISS_SrvAgeComps[1,y,1,1,sf] <- n_samp * (1 + n_bins_srv[sf] * comp_offset)

  } # end if survey

} # end i loop

# Covariates -----------------------------------------------------------------
# Qcov is the series shelikof catchability is measured through, NA in the years without a value
shel_slot <- rce_data$linkage_re_sigma == 0 # shelikof's slots in the linkage vectors
Qcov <- rce_data$linkage_re_obs_value[shel_slot]
Qcov[rce_data$linkage_re_obs_mask[shel_slot] == 0] <- NA

# the ten series of 26b's dsem as rceattle was handed them, three logged and seven scaled
dsem_vars <- rce_fits$m26b$dsem$variables
cov_names <- dsem_vars[dsem_vars != "recdevs1"]
covariates <- data.frame(year = years)
for(cov_name in cov_names) covariates[[cov_name]] <- rce_fits$m26b$dsem$y_tj[,which(dsem_vars == cov_name)]

inputs <- list(
  years = years,
  ages = ages,
  n_yrs = n_yrs,
  n_ages = n_ages,
  n_srv_fleets = n_srv,
  srv_names = srv_names,
  WAA = WAA,
  WAA_fish = WAA_fish,
  WAA_srv = WAA_srv,
  MatAA = MatAA,
  natmort = natmort,
  AgeingError = age_error,
  AgeingError_fish = AgeingError_fish,
  AgeingError_srv = AgeingError_srv,
  t_spawn = rce_data$spawn_month / 12,
  ObsCatch = ObsCatch,
  UseCatch = UseCatch,
  catch_sd = rce_data$catch_obs[1,2], # 0.05 on the log scale, every year
  ObsSrvIdx = ObsSrvIdx,
  ObsSrvIdx_SE = ObsSrvIdx_SE,
  UseSrvIdx = UseSrvIdx,
  t_srv = t_srv,
  ObsFishAgeComps = ObsFishAgeComps,
  UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,
  ObsSrvAgeComps = ObsSrvAgeComps,
  UseSrvAgeComps = UseSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,
  comp_offset = comp_offset,
  Qcov = Qcov,
  covariates = covariates,
  srv_q_prior_mean = exp(rce_data$index_log_q_prior[2]), # nmfs bottom trawl, on q at the reference age
  srv_q_prior_sd = exp(rce_fits$m26a$estimated_params$index_q_log_sd[2]),
  fish_sel_ref_age = 7, # rceattle divides each curve by its value at this age
  srv_sel_ref_age = c(3, 10, 10, 1)
)

# Rceattle's estimates ---------------------------------------------------------
# kept in rceattle's own terms, the case study converts them where the two models are written differently
rce_est <- list()

for(model_name in names(rce_fits)) {

  rce <- rce_fits[[model_name]]
  fit_data <- rce$tmb_data
  fit_par <- rce$estimated_params # parameters at the joint mode
  fit_rep <- rce$quantities

  log_slp <- fit_par$log_sel_slp[,,1] # log slope, [limb, fleet], limb 1 ascending and 2 descending
  inf_age <- fit_par$sel_inf[,,1] # age at 50%, natural scale
  re_all <- fit_rep$beta_linkage_re_all # every linkage deviation, in groups

  # fishery curve, and the annual offsets of its ascending limb
  fish_sel_pars <- c(inf_age[1,fish_flt], log_slp[1,fish_flt], inf_age[2,fish_flt], log_slp[2,fish_flt])
  names(fish_sel_pars) <- c("asc_age50", "asc_log_slope", "desc_age50", "desc_log_slope")
  fish_sel_walk <- cbind(age50 = re_all[fit_data$linkage_re_sigma == 3], # added to the age at 50%
                         log_slope = re_all[fit_data$linkage_re_sigma == 2]) # added to the log slope

  # survey curves, the descending limb for the two acoustic surveys and the ascending one for the trawl surveys
  srv_sel_pars <- rbind(c(inf_age[2,1], log_slp[2,1]),
                        c(inf_age[1,2], log_slp[1,2]),
                        c(inf_age[1,3], log_slp[1,3]),
                        c(inf_age[2,6], log_slp[2,6]))
  dimnames(srv_sel_pars) <- list(srv_names, c("age50", "log_slope"))

  # predicted indices in SPoRC's layout
  index_hat <- array(NA, dim = c(n_yrs, n_srv))
  for(i in 1:nrow(fit_data$index_ctl)) {

    sf <- match(fit_data$index_ctl[i,1], rce_srv_flt)
    y <- fit_data$index_ctl[i,3] - fit_data$styr + 1
    if(!is.na(sf)) index_hat[y,sf] <- fit_rep$index_hat[i]

  } # end i loop

  est <- list(
    ssb = fit_rep$ssb[1,1:n_yrs],
    rec = fit_rep$R[1,1:n_yrs],
    ln_R0 = fit_par$rec_pars[1,1],
    sigmaR = fit_rep$R_sd,
    rec_dev = log(fit_rep$R[1,1:n_yrs] / fit_rep$R_init), # log recruitment over mean recruitment
    log_F = fit_par$log_F[fish_flt,1:n_yrs], # at the fishery's reference age
    fish_sel_pars = fish_sel_pars,
    fish_sel_walk = fish_sel_walk,
    srv_sel_pars = srv_sel_pars,
    log_q = fit_par$index_log_q[rce_srv_flt], # at each survey's reference age
    q_dev = t(fit_rep$q_linkage_offset[rce_srv_flt,1:n_yrs]), # annual log scale offsets, [year, survey]
    q_beta = fit_par$beta_linkage_obs[1], # effect of the latent Qcov series on shelikof log q
    q_sigma = exp(fit_par$log_sigma_linkage[1]), # marginal sd of that latent series
    q_rho = 2 / (1 + exp(-2 * fit_par$trans_rho_linkage[1])) - 1, # and its correlation
    qcov_latent = re_all[shel_slot],
    adfg_q_sd = exp(fit_par$log_sigma_linkage[2]), # 0.05, fixed
    log_theta_fish = fit_par$comp_weights[fish_flt] + log(1 + n_bins_fish * comp_offset), # shifted as the sample size is
    log_theta_srv = fit_par$comp_weights[rce_srv_flt] + log(1 + n_bins_srv * comp_offset),
    index_hat = index_hat,
    nll = rce$opt$objective, # marginal negative log likelihood
    n_fixed = length(rce$opt$par)
  )

  # 26b: the dsem's coefficients by arrow name, the covariate means, and the covariate cells
  if(model_name == "m26b") {

    sem_full <- rce$dsem$sem_full
    dsem_beta <- fit_par$dsem_beta_z
    names(dsem_beta) <- sem_full$name[match(seq_along(dsem_beta), as.integer(sem_full$parameter))]
    dsem_mu <- fit_par$dsem_mu_j
    names(dsem_mu) <- dsem_vars
    dsem_x <- fit_par$dsem_x_tj
    colnames(dsem_x) <- dsem_vars

    est$dsem_beta <- dsem_beta
    est$dsem_mu <- dsem_mu[cov_names]
    est$dsem_x <- dsem_x[,cov_names]

  } # end if 26b

  rce_est[[model_name]] <- est

} # end model_name loop

sgl_rg_goa_pollock_data <- list(inputs = inputs, m26a = rce_est$m26a, m26b = rce_est$m26b)

usethis::use_data(sgl_rg_goa_pollock_data, internal = FALSE, overwrite = TRUE)
