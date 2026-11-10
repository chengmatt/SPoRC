# Purpose: Fit Rceattle's GOA pollock models 26a and 26b in SPoRC and render the case study figures
# Creator: Matthew LH. Cheng
# Date Created: 11/9/26

library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

fits <- list() # what the figures and the case study's numbers read, by model

# Build and fit, the case study's code for each model ---------------------------

for(model_name in c("26a", "26b")) {

  d <- sgl_rg_goa_pollock_data$inputs
  rce <- sgl_rg_goa_pollock_data[[paste0("m", model_name)]] # Rceattle's estimates for that model
  years <- d$years
  ages <- d$ages
  n_yrs <- d$n_yrs
  n_ages <- d$n_ages
  n_srv <- d$n_srv_fleets

  input_list <- Setup_Mod_Dim(
    n_pop = 1,
    years = years,
    ages = ages,
    lens = NA,
    n_regions = 1,
    n_sexes = 1,
    n_seas = 1,
    n_fish_fleets = 1,
    n_srv_fleets = n_srv,
    bias_correct_pe = "none", 
    bias_correct_oe = 1,
    verbose = FALSE
  )

  rec_devs_model <- "iid"
  sigmaR_spec <- "est_shared_all"
  if(model_name == "26b") {
    rec_devs_model <- "dsem" # the arrows give the deviations their distribution
    sigmaR_spec <- "fix" # the standard deviation of recruitment is a line of the arrows
  }

  input_list <- Setup_Mod_Rec(
    input_list,
    rec_model = "mean_rec",
    RecDevs_model = rec_devs_model,
    sigmaR_spec = sigmaR_spec,
    sigmaR_switch = 1,
    do_rec_bias_ramp = 0,
    init_age_strc = "scalar_no_move",
    equil_init_age_strc = "stoch_shared_ages",
    init_age_devs_shared = c(1, 1, 1, 1, 1, 1, 1, 1, 1), # all init ages share 1 dev
    t_spawn = d$t_spawn
  )

  input_list <- Setup_Mod_Biologicals(
    input_list,
    WAA = d$WAA,
    WAA_fish = d$WAA_fish,
    WAA_srv = d$WAA_srv,
    MatAA = d$MatAA,
    AgeingError = d$AgeingError,
    AgeingError_fish = d$AgeingError_fish,
    AgeingError_srv = d$AgeingError_srv,
    fit_lengths = 0,
    M_spec = "fix",
    Fixed_natmort = d$natmort
  )

  input_list <- Setup_Mod_Catch_and_F(
    input_list,
    ObsCatch = d$ObsCatch,
    UseCatch = d$UseCatch,
    Use_F_pen = 0,
    ln_F_mean_spec = "fix",
    sigmaC_spec = "fix",
    ln_sigmaC = array(log(d$catch_sd), dim = c(1, n_yrs, 1, 1))
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list,
    ObsFishIdx = array(NA, dim = c(1, n_yrs, 1, 1)),
    ObsFishIdx_SE = array(NA, dim = c(1, n_yrs, 1, 1)),
    UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
    ObsFishAgeComps = d$ObsFishAgeComps,
    UseFishAgeComps = d$UseFishAgeComps,
    ISS_FishAgeComps = d$ISS_FishAgeComps,
    ObsFishLenComps = array(NA, dim = c(1, n_yrs, 1, 1, 1, 1)),
    UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    fish_idx_type = "none",
    FishAgeComps_LikeType = "Dirichlet-Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1",
    FishAgeComps_bins = list(2:10)
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list,
    ObsSrvIdx = d$ObsSrvIdx,
    ObsSrvIdx_SE = d$ObsSrvIdx_SE,
    UseSrvIdx = d$UseSrvIdx,
    ObsSrvAgeComps = d$ObsSrvAgeComps,
    UseSrvAgeComps = d$UseSrvAgeComps,
    ISS_SrvAgeComps = d$ISS_SrvAgeComps,
    ObsSrvLenComps = array(NA, dim = c(1, n_yrs, 1, 1, 1, n_srv)),
    UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
    ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
    srv_idx_type = c("biom", "biom", "biom", "biom"),
    SrvAgeComps_LikeType = c("Dirichlet-Multinomial", "Dirichlet-Multinomial", "Dirichlet-Multinomial", "Dirichlet-Multinomial"),
    SrvLenComps_LikeType = c("none", "none", "none", "none"),
    SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2",
                         "agg_Year_1-terminal_Fleet_3", "agg_Year_1-terminal_Fleet_4"),
    SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2",
                         "none_Year_1-terminal_Fleet_3", "none_Year_1-terminal_Fleet_4"),
    SrvAgeComps_bins = list(3:10, NULL, NULL, NULL)
  )

  fish_selex_prior <- data.frame(
    region = 1,
    fleet = 1,
    block = 1,
    sex = 1,
    par = c(1, 2, 3, 4), # ascending age at 50%, ascending slope, descending age at 50%, descending slope
    mu = c(0, exp(-1), 10, exp(-1)),
    sd = c(3, 1.5, 3, 1.5),
    type = c("par_normal", "par", "par_normal", "par")
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list,
    fish_sel_model = "dbllogist1_Fleet_1",
    fish_sel_blocks = "none_Fleet_1",
    fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix",
    cont_tv_fish_sel = "rw_Fleet_1",
    fish_sel_devs_pars = list(1:2),
    fish_sel_devs_form = list(c("additive", "multiplicative", "multiplicative", "multiplicative")),
    fish_sel_devs_spec = "est_all",
    fishsel_pe_pars_spec = "fix",
    fishsel_dont_est_dev_first = 1,
    Use_fish_selex_prior = 1,
    fish_selex_prior = fish_selex_prior
  )

  # figure out survey's curve at Rceattle's reference age
  srv_limb <- c(1, -1, -1, 1) # 1 for a descending curve, -1 for an ascending one
  srv_ref_value <- 1 / (1 + exp(srv_limb * exp(rce$srv_sel_pars[, "log_slope"]) * (d$srv_sel_ref_age - rce$srv_sel_pars[, "age50"])))

  srv_q_prior <- data.frame(
    region = 1,
    block = 1,
    fleet = 2,
    mu = d$srv_q_prior_mean / srv_ref_value[2],
    sd = d$srv_q_prior_sd
  )

  srv_selex_prior <- data.frame(
    region = 1,
    block = 1,
    sex = 1,
    fleet = c(1, 1, 2, 2, 3, 3, 4, 4),
    par = c(1, 2, 1, 2, 1, 2, 1, 2), # age at 50% then slope
    mu = c(10, exp(-1), 0, exp(-1), 0, exp(-1), 10, exp(-1)),
    sd = c(3, 1.5, 3, 1.5, 3, 1.5, 3, 1.5),
    type = c("par_normal", "par", "par_normal", "par", "par_normal", "par", "par_normal", "par")
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list,
    srv_sel_model = c("desclogist1_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3", "desclogist1_Fleet_4"),
    srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4"),
    srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4"),
    srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all", "est_all"),
    srv_q_spec = c("est_all", "est_all", "est_all", "est_all"),
    cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4"),
    srv_sel_zero_bins = list(1:2, NULL, NULL, NULL),
    srv_q_model = c("dsem", "none", "dsem", "none"),
    sigma_srv_q_spec = "fix",
    srv_q_rho_spec = "fix",
    srv_q_re_years = list(NULL, NULL, 2:n_yrs, NULL),
    Use_srv_q_prior = 1,
    srv_q_prior = srv_q_prior,
    Use_srv_selex_prior = 1,
    srv_selex_prior = srv_selex_prior,
    t_srv = d$t_srv
  )

  dsem_arrows <- c(
    "Qcov -> Qcov, 1, rho_q", # the series Qcov measures, correlated from year to year
    "Qcov <-> Qcov, 0, sd_q", # and its innovation sd
    "Qcov -> srv_q_Region_1_Fleet_1, 0, beta_q", # Shelikof log q deviation is beta_q times that series
    "srv_q_Region_1_Fleet_1 <-> srv_q_Region_1_Fleet_1, 0, NA, 0", # exactly, with no noise of its own
    "srv_q_Region_1_Fleet_3 -> srv_q_Region_1_Fleet_3, 1, NA, 1", # ADF&G log q deviation, a random walk
    "srv_q_Region_1_Fleet_3 <-> srv_q_Region_1_Fleet_3, 0, NA, 0.05" # with a step sd fixed at 0.05
  )
  dsem_data <- data.frame(year = years, Qcov = d$Qcov)
  dsem_family <- c(Qcov = "gaussian_fixed_sd") # Qcov measures its series with a known error
  dsem_fixed_sd <- data.frame(year = years, Qcov = 0.02) # the sd of that error
  dsem_mu_spec <- "zero" # mean stays at zero

  if(model_name == "26b") {

    dsem_arrows <- c(
      dsem_arrows,
      "Spr_SST -> Spr_SST, 1, AR_SST",
      "Wind_NS -> Wind_NS, 1, AR_Wind",
      "Sum_LCopepod_Shelikof -> Sum_LCopepod_Shelikof, 1, AR_Cope",
      "Sum_Euph_Kodiak -> Sum_Euph_Kodiak, 1, AR_Euph",
      "Spr_Larvae_Shelikof -> Spr_Larvae_Shelikof, 1, AR_larvae",
      "Sum_OffYOY_Cond_Shelikof -> Sum_OffYOY_Cond_Shelikof, 1, AR_CondOffYOY",
      "Sum_OffYOY_Shelikof -> Sum_OffYOY_Shelikof, 1, AR_OffYOY",
      "Sum_NearYOY_Kodiak -> Sum_NearYOY_Kodiak, 1, AR_NearYOY",
      "Sum_Juv_EuphDiet -> Sum_Juv_EuphDiet, 1, AR_JuvEuphDiet",
      "Fal_Adult_Cond_Fishery -> Fal_Adult_Cond_Fishery, 1, AR_CondAd",
      "rec <-> rec, 0, sigmaR, 1",
      "Sum_Euph_Kodiak -> Fal_Adult_Cond_Fishery, 0, Euph_to_CondAd",
      "Fal_Adult_Cond_Fishery -> Spr_Larvae_Shelikof, 1, CondAd_to_Larvae",
      "Wind_NS -> Spr_Larvae_Shelikof, 0, Wind_to_Larvae",
      "Spr_Larvae_Shelikof -> Sum_OffYOY_Shelikof, 0, Larvae_to_OffYOY",
      "Sum_OffYOY_Shelikof <-> Sum_NearYOY_Kodiak, 0, OffYOY_to_NearYOY",
      "Sum_OffYOY_Shelikof -> rec, 1, offYOY_to_R",
      "Spr_SST -> rec, 1, SST_to_R",
      "Spr_SST -> Sum_OffYOY_Cond_Shelikof, 0, SST_to_CondYOY",
      "Sum_OffYOY_Cond_Shelikof -> rec, 1, CondYOY_to_R"
    )
    dsem_data <- cbind(dsem_data, d$covariates[, -1])
    dsem_family[names(d$covariates)[-1]] <- "fixed" # the ten series are known in the years they were observed
    dsem_mu_spec <- names(d$covariates)[-1] # the ten series whose means are estimated

  } # end if 26b

  input_list <- Setup_Mod_DSEM(
    input_list,
    dsem_arrows = dsem_arrows,
    dsem_data = dsem_data,
    dsem_family = dsem_family,
    dsem_fixed_sd = dsem_fixed_sd,
    dsem_mu_spec = dsem_mu_spec
  )

  input_list <- Setup_Mod_Weighting(
    input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_Init_Rec = 0,
    Wt_F = 1,
    Wt_Tagging = 0,
    addtocomp = d$comp_offset,
    comp_const_obs = 1,
    addtosrvidx = 0,
    addtofishidx = 0
  )

  input_list$par$ln_global_R0 <- rce$ln_R0
  input_list$par$ln_sigmaR[] <- log(rce$sigmaR)
  input_list$par$ln_RecDevs[1, 1, ] <- rce$rec_dev
  input_list$par$ln_InitDevs[1, 1, , 1] <- rce$rec_dev[1]

  # fishery curve, and the annual deviations of its ascending limb
  fish_sel <- rce$fish_sel_pars
  fish_walk <- rce$fish_sel_walk
  input_list$par$fish_fixed_sel_pars[1, , 1, 1, 1] <- c(log(fish_sel["asc_age50"]), fish_sel["asc_log_slope"],
                                                        log(fish_sel["desc_age50"]), fish_sel["desc_log_slope"])
  input_list$par$ln_fishsel_devs[1, , 1, 1, 1] <- fish_walk[, "age50"] # added to the age at 50%
  input_list$par$ln_fishsel_devs[1, , 2, 1, 1] <- fish_walk[, "log_slope"] # added to the log slope
  input_list$par$fishsel_pe_pars[1, 1, 1, 1] <- log(0.2) # step sd of the age at 50%, in years of age
  input_list$par$fishsel_pe_pars[1, 2, 1, 1] <- log(0.05) # step sd of the log slope

  ref_age <- d$fish_sel_ref_age
  asc_ref <- 1 / (1 + exp(-exp(fish_sel["asc_log_slope"] + fish_walk[, "log_slope"]) * (ref_age - fish_sel["asc_age50"] - fish_walk[, "age50"])))
  desc_ref <- 1 / (1 + exp(exp(fish_sel["desc_log_slope"]) * (ref_age - fish_sel["desc_age50"])))
  input_list$par$ln_F_mean[] <- 0
  input_list$par$ln_F_devs[1, , 1, 1] <- rce$log_F - log(asc_ref * desc_ref)

  input_list$par$srv_fixed_sel_pars[1, 1, 1, 1, ] <- log(rce$srv_sel_pars[, "age50"])
  input_list$par$srv_fixed_sel_pars[1, 2, 1, 1, ] <- rce$srv_sel_pars[, "log_slope"]
  input_list$par$ln_srv_q[1, 1, ] <- rce$log_q - log(srv_ref_value)
  input_list$par$ln_srv_q_devs[1, , ] <- rce$q_dev

  input_list$par$ln_FishAge_theta_agg[] <- rce$log_theta_fish
  input_list$par$ln_SrvAge_theta_agg[] <- rce$log_theta_srv

  dsem_start <- c(rho_q = rce$q_rho,
                  sd_q = rce$q_sigma * sqrt(1 - rce$q_rho^2),
                  beta_q = rce$q_beta)
  if(model_name == "26b") dsem_start <- c(dsem_start, rce$dsem_beta)

  dsem_model <- input_list$data$dsem_model
  input_list$par$dsem_beta[] <- dsem_start[dsem_model$beta_names]
  input_list$par$ln_dsem_sd[] <- log(abs(dsem_start[dsem_model$ln_sd_names]))
  input_list$par$dsem_x[, match("Qcov", dsem_model$variables)] <- rce$qcov_latent # the series Qcov measures

  if(model_name == "26b") {
    cov_names <- names(d$covariates)[-1]
    input_list$par$dsem_mu[match(cov_names, dsem_model$variables)] <- rce$dsem_mu[cov_names]
  }

  re_names <- c("ln_RecDevs", "ln_srv_q_devs", "dsem_x")
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = re_names,
                   do_optim = TRUE, newton_loops = 3, silent = TRUE)
  sd_rep <- RTMB::sdreport(fit)

  # reported series with their standard errors, and the estimates the case study quotes
  par_list <- fit$env$parList()
  dsem_est <- c(as.vector(par_list$dsem_beta), exp(as.vector(par_list$ln_dsem_sd)))
  names(dsem_est) <- c(dsem_model$beta_names, dsem_model$ln_sd_names)
  dsem_mu_est <- as.vector(par_list$dsem_mu)
  names(dsem_mu_est) <- dsem_model$variables
  fixed_se <- summary(sd_rep, "fixed")[, "Std. Error"]
  dsem_se <- fixed_se[names(fixed_se) == "dsem_beta"] # standard error of each path, in the order of the arrows
  names(dsem_se) <- dsem_model$beta_names

  fits[[model_name]] <- list(
    input_list = input_list, # data, starting values and map, for anything that refits the model
    ssb = fit$rep$SSB[1, 1, ],
    rec = fit$rep$Rec[1, 1, ],
    ssb_se = bridge_se(sd_rep, "log_SSB", fit$rep$SSB[1, 1, ]),
    rec_se = bridge_se(sd_rep, "log_Rec", fit$rep$Rec[1, 1, ]),
    ln_srv_q = par_list$ln_srv_q[1, 1, ], # log q before the annual deviations
    srv_q = fit$rep$srv_q[1, , ], # catchability by year and survey
    srv_sel = fit$rep$srv_sel[1, 1, , 1, , 1, ], # selectivity by year, age and survey
    dsem_est = dsem_est,
    dsem_se = dsem_se,
    dsem_mu_est = dsem_mu_est,
    sigmaR = exp(par_list$ln_sigmaR[1]),
    init_dev = par_list$ln_InitDevs[1, 1, 1, 1], # the deviation ages 2 to 10 share in 1970
    rec_dev_first = par_list$ln_RecDevs[1, 1, 1],
    nll = fit$optim$objective,
    max_grad = max(abs(fit$gr(fit$optim$par))),
    n_fixed = length(fit$optim$par),
    n_random = length(fit$env$random),
    pd_hess = sd_rep$pdHess
  )

} # end model_name loop

saveRDS(fits, here("dev", "dev_output", "goa_pollock_bridge.rds"))

# Spawning biomass and recruitment against Rceattle ------------------------------

for(model_name in c("26a", "26b")) {

  rce <- sgl_rg_goa_pollock_data[[paste0("m", model_name)]]
  ggplot2::ggsave(
    here("vignettes", "figures", paste0("al_goa_pollock_ts_", model_name, ".png")),
    bridge_ts_figure(
      yrs = years,
      ssb = fits[[model_name]]$ssb,
      rec = fits[[model_name]]$rec,
      admb_ssb = rce$ssb,
      admb_rec = rce$rec,
      label = paste0("Rceattle (Model ", model_name, ")"),
      ssb_se = fits[[model_name]]$ssb_se,
      rec_se = fits[[model_name]]$rec_se,
      ref_name = "Rceattle"
    ),
    width = 17,
    height = 9,
    dpi = 150
  )

} # end model_name loop

# Catchability of the two surveys whose q changes by year, model 26a --------------

# both models at the survey's reference age, which is where Rceattle's q applies
rce <- sgl_rg_goa_pollock_data$m26a
fit_a <- fits[["26a"]]
q_surveys <- d$srv_names[c(1, 3)]
q_df <- NULL

for(sf in c(1, 3)) {

  ref_sel <- fit_a$srv_sel[, d$srv_sel_ref_age[sf], sf] # SPoRC's curve at the reference age, by year
  q_df <- rbind(
    q_df,
    data.frame(Year = years, value = fit_a$srv_q[, sf] * ref_sel, Survey = d$srv_names[sf], type = "SPoRC"),
    data.frame(Year = years, value = exp(rce$log_q[sf] + rce$q_dev[, sf]), Survey = d$srv_names[sf], type = "Rceattle (Model 26a)")
  )

} # end sf loop
q_df$Survey <- factor(q_df$Survey, levels = q_surveys)

# shelikof q the covariate alone implies, the base q times exp(beta_q times Qcov)
qcov_df <- data.frame(
  Year = years,
  value = exp(fit_a$ln_srv_q[1] + fit_a$dsem_est["beta_q"] * d$Qcov) * fit_a$srv_sel[, d$srv_sel_ref_age[1], 1],
  Survey = factor(q_surveys[1], levels = q_surveys)
)
qcov_df <- qcov_df[!is.na(qcov_df$value), ]

p_q <- ggplot(q_df, aes(Year, value)) +
  geom_line(aes(color = type, linetype = type), linewidth = 1) +
  geom_point(data = qcov_df, size = 2.5, shape = 1) +
  facet_wrap(~Survey, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Catchability", color = "Type", linetype = "Type")
ggplot2::ggsave(here("vignettes", "figures", "al_goa_pollock_q.png"), p_q, width = 16, height = 7, dpi = 150)

# The dsem of model 26b: its paths as a graph, and each estimate against Rceattle's ---------

fit_b <- fits[["26b"]]
rce <- sgl_rg_goa_pollock_data$m26b

# the series of the recruitment dsem, where each sits in the graph, and the name its arrows use for its ar1 path
dag_nodes <- data.frame(
  series = c("Sum_Euph_Kodiak", "Fal_Adult_Cond_Fishery", "Wind_NS", "Spr_Larvae_Shelikof", "Sum_OffYOY_Shelikof",
             "Sum_NearYOY_Kodiak", "Spr_SST", "Sum_OffYOY_Cond_Shelikof", "rec", "Sum_LCopepod_Shelikof", "Sum_Juv_EuphDiet"),
  label = c("Summer euphausiids", "Fall adult condition", "Spring wind", "Spring larvae", "Offshore YOY",
            "Nearshore YOY", "Spring SST", "YOY condition", "Recruitment",
            "Summer copepods", "Juvenile euphausiid diet"),
  x = c(0, 2, 2, 4, 6, 6, 4, 6, 8.5, 0, 0),
  y = c(3, 3, 4.2, 3, 3, 4.2, 0.7, 1.9, 2.2, 1.7, 0.7),
  ar_name = c("AR_Euph", "AR_CondAd", "AR_Wind", "AR_larvae", "AR_OffYOY", "AR_NearYOY", "AR_SST", "AR_CondOffYOY", NA,
              "AR_Cope", "AR_JuvEuphDiet")
)

# each node reads its series, then its own year to year path, or for recruitment the sd the paths leave
dag_nodes$text <- paste0(dag_nodes$label, "\nAR ", formatC(fit_b$dsem_est[dag_nodes$ar_name], format = "f", digits = 2))
dag_nodes$text[dag_nodes$series == "rec"] <- paste0("Recruitment\nsd ", formatC(fit_b$dsem_est["sigmaR"], format = "f", digits = 2))

# the paths between series, by arrow name. lag 1 is an effect a year later, the last is a covariance
dag_edges <- data.frame(
  name = c("Euph_to_CondAd", "CondAd_to_Larvae", "Wind_to_Larvae", "Larvae_to_OffYOY", "offYOY_to_R", "SST_to_R",
           "SST_to_CondYOY", "CondYOY_to_R", "OffYOY_to_NearYOY"),
  from = c("Sum_Euph_Kodiak", "Fal_Adult_Cond_Fishery", "Wind_NS", "Spr_Larvae_Shelikof", "Sum_OffYOY_Shelikof", "Spr_SST",
           "Spr_SST", "Sum_OffYOY_Cond_Shelikof", "Sum_OffYOY_Shelikof"),
  to = c("Fal_Adult_Cond_Fishery", "Spr_Larvae_Shelikof", "Spr_Larvae_Shelikof", "Sum_OffYOY_Shelikof", "rec", "rec",
         "Sum_OffYOY_Cond_Shelikof", "rec", "Sum_NearYOY_Kodiak"),
  lag = c(0, 1, 0, 0, 1, 1, 0, 1, 0),
  covariance = c(0, 0, 0, 0, 0, 0, 0, 0, 1),
  text_at = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.5, 0.5, 0.5) # how far along its line each estimate is written
)
dag_edges$est <- fit_b$dsem_est[dag_edges$name]
dag_edges$se <- fit_b$dsem_se[dag_edges$name]
dag_edges$rceattle <- rce$dsem_beta[dag_edges$name]

# end points, pulled back to the edge of each node's box so the arrow heads show
from_row <- match(dag_edges$from, dag_nodes$series)
to_row <- match(dag_edges$to, dag_nodes$series)
box_half_width <- 0.031 * nchar(dag_nodes$label) + 0.2 # in x units, wider for a longer label
box_half_height <- 0.3
dx <- dag_nodes$x[to_row] - dag_nodes$x[from_row]
dy <- dag_nodes$y[to_row] - dag_nodes$y[from_row]
shrink_from <- pmin(box_half_width[from_row] / abs(dx), box_half_height / abs(dy)) # fraction of the line inside the first box
shrink_to <- pmin(box_half_width[to_row] / abs(dx), box_half_height / abs(dy)) # and inside the second
dag_edges$x <- dag_nodes$x[from_row] + shrink_from * dx
dag_edges$y <- dag_nodes$y[from_row] + shrink_from * dy
dag_edges$xend <- dag_nodes$x[to_row] - shrink_to * dx
dag_edges$yend <- dag_nodes$y[to_row] - shrink_to * dy

# each estimate sits beside the middle of its line, a little way off it
edge_length <- sqrt(dx^2 + dy^2)
dag_edges$x_text <- dag_edges$x + dag_edges$text_at * (dag_edges$xend - dag_edges$x) - 0.2 * dy / edge_length
dag_edges$y_text <- dag_edges$y + dag_edges$text_at * (dag_edges$yend - dag_edges$y) + 0.2 * abs(dx) / edge_length
dag_edges$Effect <- ifelse(dag_edges$covariance == 1, "Covariance", ifelse(dag_edges$est > 0, "Positive", "Negative"))
dag_edges$Timing <- ifelse(dag_edges$lag == 1, "A year later", "Same year")
dag_edges$text <- formatC(dag_edges$est, format = "f", digits = 2)

path_rows <- dag_edges$covariance == 0
p_dag <- ggplot() +
  geom_segment(data = dag_edges[path_rows,], aes(x = x, y = y, xend = xend, yend = yend, color = Effect, linetype = Timing),
               linewidth = 1, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), key_glyph = "path") +
  geom_segment(data = dag_edges[!path_rows,], aes(x = x, y = y, xend = xend, yend = yend, color = Effect),
               linewidth = 1, arrow = arrow(length = unit(0.3, "cm"), ends = "both"), key_glyph = "path") +
  geom_text(data = dag_edges, aes(x = x_text, y = y_text, label = text), size = 5.5) +
  geom_label(data = dag_nodes, aes(x = x, y = y, label = text), size = 5, lineheight = 0.9, label.padding = unit(0.5, "lines")) +
  scale_color_manual(values = c(Positive = "#2166AC", Negative = "#B2182B", Covariance = "grey45")) +
  scale_linetype_manual(values = c("Same year" = "solid", "A year later" = "dashed")) +
  coord_cartesian(xlim = c(-1, 9.4), ylim = c(0.3, 4.6)) +
  theme_void(base_size = 20) +
  theme(legend.position = "top", legend.key.width = unit(1.6, "cm"), plot.margin = margin(15, 15, 15, 15)) +
  labs(color = "Path", linetype = "Timing")
ggplot2::ggsave(here("vignettes", "figures", "al_goa_pollock_dag_26b.png"), p_dag, width = 16, height = 9, dpi = 150, bg = "white")

# every path with its 95% interval, the paths between series and each series' own year to year path
ar_rows <- !is.na(dag_nodes$ar_name)
effect_df <- rbind(
  data.frame(panel = "Paths between series",
             path = paste0(dag_nodes$label[from_row], ifelse(dag_edges$covariance == 1, " with ", " to "), dag_nodes$label[to_row]),
             est = dag_edges$est, se = dag_edges$se, rceattle = dag_edges$rceattle),
  data.frame(panel = "Year to year path of each series (AR1)",
             path = dag_nodes$label[ar_rows],
             est = fit_b$dsem_est[dag_nodes$ar_name[ar_rows]], se = fit_b$dsem_se[dag_nodes$ar_name[ar_rows]],
             rceattle = rce$dsem_beta[dag_nodes$ar_name[ar_rows]])
)
effect_df$path <- factor(effect_df$path, levels = rev(unique(effect_df$path)))

p_effect <- ggplot(effect_df, aes(y = path)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_linerange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), color = "#E69F00", linewidth = 1) +
  geom_point(aes(x = est, shape = "SPoRC"), color = "#E69F00", size = 4) +
  geom_point(aes(x = rceattle, shape = "Rceattle (Model 26b)"), color = "black", size = 4, stroke = 1.2) +
  scale_shape_manual(values = c("SPoRC" = 16, "Rceattle (Model 26b)" = 1)) +
  facet_wrap(~panel, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "Estimate", y = NULL, shape = "Type")
ggplot2::ggsave(here("vignettes", "figures", "al_goa_pollock_effects_26b.png"), p_effect, width = 14, height = 11, dpi = 150)
