# Purpose: To demonstrate a full assessment workflow using Dusky Rockfish as an example
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 8/18/25


# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(PBSmodelling)
devtools::load_all(here("R"))

# load in data
data("sgl_rg_dusky_data")

# Setup Model -------------------------------------------------------------

input_list <- Setup_Mod_Dim(
  years = sgl_rg_dusky_data$years,
  # vector of years
  ages = sgl_rg_dusky_data$mod_ages,
  # vector of ages
  lens = sgl_rg_dusky_data$lens,
  # number of lengths
  n_regions = sgl_rg_dusky_data$n_regions,
  # number of regions
  n_sexes = sgl_rg_dusky_data$n_sexes,
  # number of sexes
  n_fish_fleets = sgl_rg_dusky_data$n_fish_fleets,
  # number of fishery fleets
  n_srv_fleets = sgl_rg_dusky_data$n_srv_fleets, # number of survey fleets
  n_seas = sgl_rg_dusky_data$n_seas, # number of seasons
  # Population stuff
  n_pop = sgl_rg_dusky_data$n_pop,
  natal_region = sgl_rg_dusky_data$natal_region,
  verbose = TRUE # whether to output messages
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(
  input_list = input_list,

  # Model options
  # Doing bias ramp, but basically setting it so that no lognormal bias correction happens (as in the dusky model)
  do_rec_bias_ramp = 1,
  bias_year = rep(length(sgl_rg_dusky_data$years), 4),
  # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
  sigmaR_switch = 1,
  # when to switch from early to late sigmaR (switch in first year)
  ln_sigmaR = array(-0.1068576, dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
  # Starting values for early and late sigmaR
  rec_model = "mean_rec",
  sigmaR_spec = "fix",
  # fix early sigmaR and late sigmaR
  # recruitment sex ratio
  init_age_strc = 1, # geometric series to derive initial age structure
  ln_global_R0 = log(2.7), # starting value for mean_rec
  t_spawn = sgl_rg_dusky_data$spwn_month
)

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,

  # Data inputs
  WAA = sgl_rg_dusky_data$waa_arr,
  MatAA = sgl_rg_dusky_data$mataa_arr,

  # Model options
  # fit lengths
  fit_lengths = 1,
  SizeAgeTrans = sgl_rg_dusky_data$sizeage,
  AgeingError = sgl_rg_dusky_data$age_error_matrix,
  M_spec = "fix",
  # fixing natural mortality
  Fixed_natmort = sgl_rg_dusky_data$fix_natmort,
  addtocomp = 0.00001
)

# Setup movement stuff (using defaults for other stuff)
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = NA,
  do_recruits_move = 0
)

input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)


input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,

  # Data inputs
  ObsCatch = sgl_rg_dusky_data$ObsCatch,
  UseCatch = sgl_rg_dusky_data$UseCatch,

  # Model options
  Use_F_pen = 1,
  # whether to use f penalty, == 0 don't use, == 1 use
  sigmaC_spec = "fix",
  Catch_Constant = 0.00001,

  # Fixing sigma C and F
  ln_sigmaC = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
  ln_sigmaF = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,

  # data inputs
  ObsFishIdx = sgl_rg_dusky_data$ObsFishIdx, # fishery index
  ObsFishIdx_SE = sgl_rg_dusky_data$ObsFishIdx_SE, # standard errors
  UseFishIdx = sgl_rg_dusky_data$UseFishIdx, # whether fishery indices are used
  ObsFishAgeComps = sgl_rg_dusky_data$ObsFishAgeComps, # observed fishery ages
  UseFishAgeComps = sgl_rg_dusky_data$UseFishAgeComps, # whether fishery ages are used
  ISS_FishAgeComps = sgl_rg_dusky_data$ISS_FishAgeComps, # input sample size for fishery ages
  ObsFishLenComps = sgl_rg_dusky_data$ObsFishLenComps, # observed fishery lengths
  UseFishLenComps = sgl_rg_dusky_data$UseFishLenComps, # whether fishery lengths are used
  ISS_FishLenComps = sgl_rg_dusky_data$ISS_FishLenComps, # input sample size for fishery lengths

  # Model options
  fish_idx_type = c("none"),
  # indices for fishery
  FishAgeComps_LikeType = c("Multinomial"),
  # age comp likelihoods for fishery fleet
  FishLenComps_LikeType = c("Multinomial"),
  # length comp likelihoods for fishery
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
  # age comp structure for fishery
  FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1")
  # length comp structure for fishery
)

# Setup survey indices and compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,

  # data inputs
  ObsSrvIdx = sgl_rg_dusky_data$ObsSrvIdx, # observed survey index
  ObsSrvIdx_SE = sgl_rg_dusky_data$ObsSrvIdx_SE / sgl_rg_dusky_data$ObsSrvIdx, # lognormal SD
  UseSrvIdx = sgl_rg_dusky_data$UseSrvIdx, # whether survey indices are used
  ObsSrvAgeComps = sgl_rg_dusky_data$ObsSrvAgeComps, # observed survey ages
  ISS_SrvAgeComps = sgl_rg_dusky_data$ISS_SrvAgeComps, # input sample size for survey ages
  UseSrvAgeComps = sgl_rg_dusky_data$UseSrvAgeComps, # whether survey ages are used
  ObsSrvLenComps = sgl_rg_dusky_data$ObsSrvLenComps, # observed survey lengths
  UseSrvLenComps = sgl_rg_dusky_data$UseSrvLenComps, # whether survey lengths are used
  ISS_SrvLenComps = sgl_rg_dusky_data$ISS_SrvLenComps, # input sample size for survey lengths

  # Model options
  srv_idx_type = c("biom"),
  # abundance and biomass for survey fleet 1
  SrvAgeComps_LikeType = c("Multinomial"),
  # survey age composition likelihood for survey fleet 1
  SrvLenComps_LikeType = c("Multinomial"),
  #  survey length composition likelihood for survey fleet 1
  SrvAgeComps_Type = c(
    "agg_Year_1-terminal_Fleet_1"
  ),
  # survey age comp type

  SrvLenComps_Type = c(
    "agg_Year_1-terminal_Fleet_1"
  )
  # survey length comp type
)


input_list <- Setup_Mod_Fishsel_and_Q(

  input_list = input_list,

  # Model options
  # fishery selectivity, whether continuous time-varying
  cont_tv_fish_sel = c("none_Fleet_1"),
  # fishery selectivity blocks
  fish_sel_blocks = c("none_Fleet_1"),
  # fishery selectivity form
  fish_sel_model = c("logist2_Fleet_1"),
  # fishery catchability blocks
  fish_q_blocks = c("none_Fleet_1"),
  # whether to estiamte all fixed effects for fishery selectivity
  fish_fixed_sel_pars_spec = c("est_all"),
  # whether to estiamte all fixed effects for fishery catchability
  fish_q_spec = c("fix")
)

# Setup survey selectivity and catchability
# Set up prior for survey catchability
srv_q_prior <- data.frame(
  region = 1,
  block = 1,
  fleet = 1,
  mu = 1,
  sd = 0.447213595
)

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,

  # Model options
  # survey selectivity, whether continuous time-varying
  cont_tv_srv_sel = c("none_Fleet_1"),
  # survey selectivity blocks
  srv_sel_blocks = c("none_Fleet_1"),
  # survey selectivity form
  srv_sel_model = c("logist2_Fleet_1"),
  # survey catchability blocks
  srv_q_blocks = c("none_Fleet_1"),
  # whether to estiamte all fixed effects for survey selectivity
  srv_fixed_sel_pars_spec = c("est_all"),
  # whether to estiamte all fixed effects for survey catchability
  srv_q_spec = c("est_all"),
  Use_srv_q_prior = 1,
  # Use catchability prior
  srv_q_prior = srv_q_prior,
  # survey timing
  t_srv = array(0, dim = c(input_list$data$n_regions,
                           input_list$data$n_seas,
                           input_list$data$n_srv_fleets))
)

# catch weigthing for duskies
Wt_Catch <- array(0, dim = c(sgl_rg_dusky_data$n_regions, length(sgl_rg_dusky_data$years), input_list$data$n_seas, sgl_rg_dusky_data$n_fish_fleets))
Wt_Catch[,which(sgl_rg_dusky_data$years %in% 1977:1991),,] <- 2
Wt_Catch[,-which(sgl_rg_dusky_data$years %in% 1977:1991),,] <- 50

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = Wt_Catch,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1.66,
  Wt_Rec = 1,
  Wt_F = 2,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_seas,
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_seas,
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_seas,
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets)),
  Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_seas,
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets))
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Optimize Values (without Francis) ---------------------------------------------------------

# Fit model
nofrancis_model <- fit_model(data,
                             parameters,
                             mapping,
                             random = NULL,
                             newton_loops = 3,
                             silent = FALSE
)

nofrancis_model$sdrep <- RTMB::sdreport(nofrancis_model) # get standard error report


# Optimize Values (with Francis) ------------------------------------------
francis_data <- data # redefine data list for francis (replacing data weights)

# run francis reweighting
francis_runs <- run_francis(francis_data,
                             parameters,
                             mapping,
                             n_francis_iter = 10,
                             newton_loops =  3)

# get obj
francis_model <- francis_runs$obj
francis_data <- francis_runs$obj$data
francis_model$sdrep <- RTMB::sdreport(francis_model) # get standard error report

# Check Model -------------------------------------------------------------

# check no francis model
post_optim_sanity_checks(nofrancis_model$sdrep,
                         nofrancis_model$rep,
                         gradient_tol = 1e-4,
                         se_tol = 5,
                         corr_tol = 0.95)

# check francis model
post_optim_sanity_checks(francis_model$sdrep,
                         francis_model$rep,
                         gradient_tol = 1e-4,
                         se_tol = 5,
                         corr_tol = 0.95)


# Compare Time Series -----------------------------------------------------

ts_plot <- get_ts_plot(rep = list(francis_model$rep, nofrancis_model$rep),
                       sd_rep = list(francis_model$sdrep, nofrancis_model$sdrep),
                       model_names = c("With Francis", "No Francis"))

png(here("vignettes", "figures", "o_ts_comb.png"), width = 1000, height = 800)
ts_plot[[1]]
dev.off()

# Catch Fits --------------------------------------------------------------

png(here("vignettes", "figures", "o_catch_fits.png"), width = 1000, height = 800)
get_catch_fits_plot(list(data, francis_data),
                    list(nofrancis_model$rep, francis_model$rep),
                    c("No Francis", "With Francis"))
dev.off()

# Index Fits --------------------------------------------------------------

png(here("vignettes", "figures", "o_idx_fits.png"), width = 1000, height = 800)
get_idx_fits_plot(list(data, francis_data),
                  list(nofrancis_model$rep, francis_model$rep),
                  c("No Francis", "With Francis"))
dev.off()

# Composition Fits --------------------------------------------------------

# Extract observed and expected compositions
comp_prop <- get_comp_prop(francis_data,
                           francis_model$rep,
                           age_labels = 4:30,
                           len_labels = 21:52,
                           year_labels = 1977:2024)


### Fishery Ages ------------------------------------------------------------

# get one step ahead fishery ages
fishages <- get_osa(obs_mat = comp_prop$Obs_FishAge_mat, # observed fishery age compositions
                    exp_mat = comp_prop$Pred_FishAge_mat, # predicted fishery age compositions
                    N = francis_data$ISS_FishAgeComps[1,which(francis_data$UseFishAgeComps[,,1,1] == 1),1,1,1] * # input sample size
                      unique(francis_data$Wt_FishAgeComps[1,which(francis_data$UseFishAgeComps[,,1,1] == 1),1,,1]), # francis weight
                    years = which(francis_data$UseFishAgeComps[,,1,1] == 1), # years with fishery ages
                    fleet = 1, # fleet
                    bins = 4:30, # age bins
                    seas = 1, # seasons
                    comp_type = 0, # composition type (age-specific)
                    bin_label = "Ages" # bin labels
                    )

# plot OSA residuals
resid_plot <- SPoRC::plot_resids(osa_results = fishages)

# Get Aggregated Plot
fishage_agg <- comp_prop$Fishery_Ages %>%
  group_by(Age, Fleet, Seas) %>%
  summarize(obs = mean(obs), pred = mean(pred)) %>%
  filter(Fleet == 1) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
  geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
  theme_bw(base_size = 15) +
  labs(x = 'Age', y = 'Proportions')

png(here("vignettes", "figures", "o_fishage_fits.png"), width = 1000, height = 800)
cowplot::plot_grid(resid_plot[[2]], # Bubble
                   cowplot::plot_grid(resid_plot[[1]], fishage_agg, nrow = 1), # QQ and Aggregated
                   ncol = 1, rel_heights = c(0.6, 0.4))
dev.off()

### Fishery Lengths ---------------------------------------------------------

# get one step ahead fishery lengths
fishlens <- get_osa(obs_mat = comp_prop$Obs_FishLen_mat, # observed fishery length compositions
                    exp_mat = comp_prop$Pred_FishLen_mat, # predicted fishery length compositions
                    N = francis_data$ISS_FishLenComps[1,which(francis_data$UseFishLenComps[,,1,1] == 1),1,1,1] * # input sample size
                      unique(francis_data$Wt_FishLenComps[1,which(francis_data$UseFishLenComps[,,1,1] == 1),1,,1]), # francis weight
                    years = which(francis_data$UseFishLenComps[,,1,1] == 1), # years with fishery ages
                    fleet = 1, # fleet
                    bins = 21:52, # age bins
                    seas = 1, # seasons
                    comp_type = 0, # composition type (age-specific)
                    bin_label = "Lengths" # bin labels
)

# plot OSA residuals
resid_plot <- SPoRC::plot_resids(osa_results = fishlens)

# Get Aggregated Plot
fishlen_agg <- comp_prop$Fishery_Lens %>%
  group_by(Len, Fleet, Seas) %>%
  summarize(obs = mean(obs), pred = mean(pred)) %>%
  filter(Fleet == 1) %>%
  ggplot() +
  geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
  geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
  theme_bw(base_size = 15) +
  labs(x = 'Lengths', y = 'Proportions')

png(here("vignettes", "figures", "o_fishlen_fits.png"), width = 1000, height = 800)
cowplot::plot_grid(resid_plot[[2]], # Bubble
                   cowplot::plot_grid(resid_plot[[1]], fishlen_agg, nrow = 1), # QQ and Aggregated
                   ncol = 1, rel_heights = c(0.6, 0.4))
dev.off()

### Survey Ages -------------------------------------------------------------

# get one step ahead survey ages
srvages <- get_osa(obs_mat = comp_prop$Obs_SrvAge_mat, # observed survey age compositions
                   exp_mat = comp_prop$Pred_SrvAge_mat, # predicted survey age compositions
                   N = francis_data$ISS_SrvAgeComps[1,which(francis_data$UseSrvAgeComps[,,1,1] == 1),1,1,1] * # input sample size
                     unique(francis_data$Wt_SrvAgeComps[1,which(francis_data$UseSrvAgeComps[,,1,1] == 1),1,,1]), # francis weight
                   years = which(francis_data$UseSrvAgeComps[,,1,1] == 1), # years with Srvery ages
                   fleet = 1, # fleet
                   bins = 4:30, # age bins
                   seas = 1, # seasons
                   comp_type = 0, # composition type (age-specific)
                   bin_label = "Ages" # bin labels
)

# plot OSA residuals
resid_plot <- SPoRC::plot_resids(osa_results = srvages)

# Get Aggregated Plot
srvage_agg <- comp_prop$Survey_Ages %>%
  group_by(Age, Fleet, Seas) %>%
  summarize(obs = mean(obs), pred = mean(pred)) %>%
  filter(Fleet == 1) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
  geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
  theme_bw(base_size = 15) +
  labs(x = 'Age', y = 'Proportions')

png(here("vignettes", "figures", "o_srvage_fits.png"), width = 1000, height = 800)
cowplot::plot_grid(resid_plot[[2]], # Bubble
                   cowplot::plot_grid(resid_plot[[1]], srvage_agg, nrow = 1), # QQ and Aggregated
                   ncol = 1, rel_heights = c(0.6, 0.4))
dev.off()

# Retrospectives ----------------------------------------------------------

# do retrospective w/ francis
francis_retro <- do_retrospective(
  n_retro = 10, # number of peels
  data = francis_data, # data list (francis data)
  parameters = parameters, # parameters list
  mapping = mapping, # mapping list
  random = NULL, # random effects
  do_par = TRUE, # whether to parallellize
  n_cores = 8,  # number of cores to use
  do_francis = TRUE,  # whether to do francis within a given retrospective peel
  n_francis_iter = 10 # number of francis iterations to run within a given retrospective peel
)

# default plots
francis_retro_plot <- get_retrospective_plot(francis_retro, 4)

png(here("vignettes", "figures", "o_francis_rel_retro.png"), width = 1000, height = 800)
francis_retro_plot[[1]]
dev.off()

png(here("vignettes", "figures", "o_francis_abs_retro.png"), width = 1000, height = 800)
francis_retro_plot[[2]]
dev.off()

png(here("vignettes", "figures", "o_francis_squid_retro.png"), width = 1000, height = 800)
francis_retro_plot[[3]]
dev.off()

# do retrospective without francis
nofrancis_retro <- do_retrospective(
  n_retro = 10, # number of peels
  data = data, # data list (not francis data)
  parameters = parameters, # parameters list
  mapping = mapping, # mapping list
  random = NULL, # random effects
  do_par = TRUE, # whether to parallellize
  n_cores = 8,  # number of cores to use
  do_francis = FALSE,  # whether to do francis within a given retrospective peel
)

# default plots
nofrancis_retro_plot <- get_retrospective_plot(nofrancis_retro, Rec_Age = 4)

png(here("vignettes", "figures", "o_nofrancis_rel_retro.png"), width = 1000, height = 800)
nofrancis_retro_plot[[1]]
dev.off()

png(here("vignettes", "figures", "o_nofrancis_abs_retro.png"), width = 1000, height = 800)
nofrancis_retro_plot[[2]]
dev.off()

png(here("vignettes", "figures", "o_nofrancis_squid_retro.png"), width = 1000, height = 800)
nofrancis_retro_plot[[3]]
dev.off()


# Likelihood Profiles -----------------------------------------------------

# do likelihood profile with francis weights
francis_meanrec_prof <- do_likelihood_profile(
  francis_data, # francis data list
  parameters, # parameter list
  mapping, # mapping list
  random = NULL, # random effects
  what = 'ln_global_R0', # parameter to profile
  min_val = log(francis_model$rep$R0) * 0.1,  # min values to profile across
  max_val = log(francis_model$rep$R0) * 2,  # max values to profile across
  inc = 0.1, # increment for min and max values to profile across
  do_par =  TRUE, # whether to parrallelize
  n_cores = 8 # number of cores
)

# summarize profile
francis_mean_rec_profile <- francis_meanrec_prof$agg_nLL %>%
  mutate(Summarized_Type = case_when(
    str_detect(type, "Pen|Prior") ~ "Other",
    str_detect(type, "Len") ~ "Length Comps",
    str_detect(type, "Age") ~ "Age Comps",
    str_detect(type, "Idx") ~ "Indices",
    str_detect(type, "Catch") ~ "Catch",
    str_detect(type, "jnLL") ~ "jnLL",
  )) %>%
  filter(value != 0) %>%
  group_by(Summarized_Type, prof_val) %>%
  summarize(value = sum(value), .groups = "drop") %>%
  group_by(Summarized_Type) %>%
  mutate(value = value - min(value))

# likelihood profile on R0
png(here("vignettes", "figures", "o_meanrec_profile.png"), width = 1000, height = 800)
ggplot(francis_mean_rec_profile, aes(x = prof_val, y = value, color = Summarized_Type)) +
  geom_line(lwd = 1.3) +
  geom_vline(xintercept = log(francis_model$rep$R0), lty = 2) +
  labs(x = 'Log Mean Recruitment', y = "Scaled nLL", color = "Type") +
  theme_bw(base_size = 15)
dev.off()

# Jitter ------------------------------------------------------------------

# get jitter results
jitter_res <- SPoRC::do_jitter(data = francis_data, # francis data list
                               parameters = parameters, # parameter list
                               mapping = mapping, # mapping list
                               random = NULL, # random effects
                               sd = 0.5, # standard deviation for jitter
                               n_jitter = 50, # number of jitters
                               n_newton_loops = 3, # newton loops to od
                               do_par = TRUE, # whether to parrallelize
                               n_cores = 8 # number of cores to use
)

# get proportion converged
prop_converged <- jitter_res %>%
  filter(Year == 1, Type == 'Recruitment') %>%
  summarize(prop_conv = sum(Hessian) / length(Hessian))

# get final model results
final_mod <- reshape2::melt(francis_model$rep$SSB) %>% rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  mutate(Type = 'SSB') %>%
  bind_rows(reshape2::melt(francis_model$rep$Rec) %>%
              rename(Pop = Var1, Region = Var2, Year = Var3) %>% mutate(Type = 'Recruitment'))

# comparison of SSB and recruitment
png(here("vignettes", "figures", "o_jiter_ts.png"), width = 1000, height = 800)
ggplot() +
  geom_line(jitter_res, mapping = aes(x = Year + 1976, y = value, group = jitter, color = Hessian), lwd = 1) +
  geom_line(final_mod, mapping = aes(x = Year + 1976, y = value), color = "black", lwd = 1.3 , lty = 2) +
  facet_grid(Type~Region, scales = 'free',
             labeller = labeller(Region = function(x) paste0("Region ", x),
                                 Type = c("Recruitment" = "Age 2 Recruitment (millions)", "SSB" = 'SSB (kt)'))) +
  labs(x = "Year", y = "Value") +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("red", 'grey')) +
  geom_text(data = jitter_res %>% filter(Type == 'SSB', Year == 1, jitter == 1),
            aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
            hjust = 1.1, vjust = 1.9, size = 6, color = "black")
dev.off()

# compare jitter of max gradient and hessian PD
png(here("vignettes", "figures", "o_jitter_res.png"), width = 1000, height = 800)
ggplot(jitter_res, aes(x = jitter, y = jnLL, color = Max_Gradient, shape = Hessian)) +
  geom_point(size = 5, alpha = 0.3) +
  geom_hline(yintercept = min(francis_model$rep$jnLL), lty = 2, size = 2, color = "blue") +
  facet_wrap(~Hessian, labeller = labeller(
    Hessian = c("FALSE" = "non-PD Hessian", "TRUE" = 'PD Hessian')
  )) +
  scale_color_viridis_c() +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barwidth = 15, barheight = 0.5)) +
  labs(x = 'Jitter') +
  geom_text(data = jitter_res %>% filter(Hessian == TRUE, Year == 1, jitter == 1),
            aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
            hjust = 1.1, vjust = 1.9, size = 6, color = "black")
dev.off()

# MCMC --------------------------------------------------------------------

# load in adnuts for MCMC
# remotes::install_github("noaa-afsc/SparseNUTS")
library(SparseNUTS)

# run MCMC
mcmc <- sample_snuts(francis_model, num_samples = 1e3, control = list(adapt_delta = 0.99))

# Check MCMC summary
diag_df <- mcmc$monitor

png(here("vignettes", "figures", "o_mcmc_pairs_trace.png"), width = 1000, height = 800)
pairs(mcmc)
dev.off()

png(here("vignettes", "figures", "o_mcmc_dens.png"), width = 1000, height = 800)
plot_marginals(mcmc, pars = 'ln_srv_q')
dev.off()

png(here("vignettes", "figures", "o_mcmc_mle_comp.png"), width = 1000, height = 800)
plot_uncertainties(mcmc)
dev.off()

# get mcmc time series plots
mcmc_ts_plot <- get_model_rep_from_mcmc(rtmb_obj = francis_model, adnuts_obj = mcmc, what = c("SSB", "Rec"), n_cores = 4)
saveRDS(mcmc_ts_plot, here("dev", "dev_output", "1_Region_Model_Dusky_Rockfish", "mcmc_results.RDS"))

png(here("vignettes", "figures", "o_mcmc_ssb_plot.png"), width = 500, height = 500)
# summarize results
ssb_summry <- mcmc_ts_plot$SSB %>%
  group_by(Var2, Var3) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975))
ggplot() +
  geom_line(ssb_summry, mapping = aes(x = Var3, y = median)) +
  geom_ribbon(ssb_summry, mapping = aes(x = Var3, y = median, ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line(reshape2::melt(francis_model$rep$SSB), mapping = aes(x = Var3, y = value), col = 'red', lwd = 1.3, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Year', y = 'Spawning Stock Biomass') +
  theme_bw(base_size = 15)
dev.off()

png(here("vignettes", "figures", "o_mcmc_rec_plot.png"), width = 500, height = 500)
# summarize results
rec_summry <- mcmc_ts_plot$Rec %>%
  group_by(Var2, Var3) %>%
  summarize(median = median(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975))
ggplot() +
  geom_line(rec_summry, mapping = aes(x = Var3, y = median)) +
  geom_ribbon(rec_summry, mapping = aes(x = Var3, y = median, ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line(reshape2::melt(francis_model$rep$Rec), mapping = aes(x = Var3, y = value), col = 'red', lwd = 1.3, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Year', y = 'Recruitment') +
  theme_bw(base_size = 15)
dev.off()

# Reference Points and Projections -------------------------------------------------------------
# Notes: To match assessment reference points exactly, need to make a few changes.
# t_spawn = 2/12
# calc_rec_st_yr = 5
### Reference Points --------------------------------------------------------

# get reference points
spr_35 <- Get_Reference_Points(data = francis_data,
                               rep = francis_model$rep,
                               SPR_x = 0.35,
                               t_spawn = 0,
                               sex_ratio_f = 0.5,
                               type = "single_region",
                               what = 'SPR',
                               calc_rec_st_yr = 3,
                               rec_age = 4)

spr_40 <- Get_Reference_Points(data = francis_data,
                               rep = francis_model$rep,
                               type = "single_region",
                               what = 'SPR',
                               SPR_x = 0.4,
                               t_spawn = 0,
                               sex_ratio_f = 0.5,
                               calc_rec_st_yr = 3, rec_age = 4)

spr_60 <- Get_Reference_Points(data = francis_data,
                               rep = francis_model$rep,
                               type = "single_region",
                               what = 'SPR',
                               SPR_x = 0.6,
                               t_spawn = 0,
                               sex_ratio_f = 0.5,
                               calc_rec_st_yr = 3, rec_age = 4)

# Extract reference points
b40 <- spr_40$b_ref_pt
b60 <- spr_60$b_ref_pt
b35 <- spr_35$b_ref_pt
f40 <- spr_40$f_ref_pt
f35 <- spr_35$f_ref_pt
f60 <- spr_60$f_ref_pt


### Projections -------------------------------------------------------------
# Define HCR to use for projections
HCR_function <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp # define stock status
  # If stock status is > 1
  if(stock_status >= 1) f <- frp
  # If stock status is between brp and alpha
  if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
  # If stock status is less than alpha
  if(stock_status < alpha) f <- 0
  return(f)
}

# define projection parameters
n_sims <- 1e3
t_spawn <- 0
n_proj_yrs <- 25
n_regions <- 1
n_ages <- length(francis_data$ages)
n_fish_fleets <- 1
n_sexes <- 1
n_seas <- 1
n_pop <- 1
do_recruits_move <- 0
terminal_NAA <- array(francis_model$rep$NAA[,,length(francis_data$years),,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
terminal_NAA0 <- array(francis_model$rep$NAA0[,,length(francis_data$years),,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
WAA <- array(rep(francis_data$WAA[,,length(francis_data$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(francis_data$WAA[,,length(francis_data$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # weight at age for fishery
MatAA <- array(rep(francis_data$MatAA[,,length(francis_data$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(francis_model$rep$fish_sel[,,length(francis_data$years),,,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- array(rep(francis_model$rep$Movement[,,,length(francis_data$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # movement - not used
terminal_F <- array(francis_model$rep$Fmort[,length(francis_data$years),,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F
natmort <- array(rep(francis_model$rep$natmort[,,length(francis_data$years),,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_ages, n_sexes)) # natural mortaility
recruitment <- array(francis_model$rep$Rec[,,3:(length(francis_data$years) - 4)], dim = c(n_pop, n_regions, length(3:(length(francis_data$years) - 4)))) # recruitment from years 3 - terminal (corresponds to 1979)
sexratio <- array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)) # recruitment sex ratio
rec_seas_prop <- array(1, dim = c(1,1)) # recruitment seasonal apportionment

# Define the F used for each scenario (Based on BSAI Intro Report - Alaska Scenarios)
proj_inputs <- list(
  # Scenario 1 - Using HCR to adjust maxFABC
  list(f_ref_pt = array(f40, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b40, dim = c(n_pop, n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 2 - Using HCR to adjust maxFABC based on last year's value (constant fraction - author specified F)
  list(f_ref_pt = array(f40 * (f40 / 0.091), dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b40, dim = c(n_pop, n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 3 - Using an F input of last 5 years average F, and
  list(f_ref_pt = array(mean(francis_model$rep$Fmort[1, 43:47,,]), dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = NULL,
       fmort_opt = 'Input'
  ),
  # Scenario 4 - Using F60
  list(f_ref_pt = array(f60, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = NULL,
       fmort_opt = 'Input'
  ),
  # Scenario 5 - F is set at 0
  list(f_ref_pt = array(0, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = NULL,
       fmort_opt = 'Input'
  ),
  # Scenario 6 - Using HCR to adjust FOFL
  list(f_ref_pt = array(f35, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b35, dim = c(n_pop, n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 7 - Using HCR to adjust FABC in first 2 projection years, and then later years are adjusting FOFL
  list(f_ref_pt = array(c(rep(f40, 2), rep(f35, n_proj_yrs - 2)), dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(c(rep(b40, 2), rep(b35, n_proj_yrs - 2)), dim = c(n_pop, n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  )
)

# store outputs
all_scenarios_f <- array(0, dim = c(n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
all_scenarios_ssb <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
all_scenarios_catch <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets, n_sims, length(proj_inputs)))

set.seed(123)
for (i in seq_along(proj_inputs)) {
  for (sim in 1:n_sims) {

    # do population projection
    out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                    n_regions = n_regions,
                                    n_ages = n_ages,
                                    n_sexes = n_sexes,
                                    n_pop = n_pop,
                                    sexratio = sexratio,
                                    n_fish_fleets = n_fish_fleets,
                                    do_recruits_move = do_recruits_move,
                                    recruitment = recruitment,
                                    terminal_NAA = terminal_NAA,
                                    terminal_NAA0 = terminal_NAA0,
                                    terminal_F = terminal_F,
                                    natmort = natmort,
                                    rec_seas_prop = rec_seas_prop,
                                    WAA = WAA,
                                    WAA_fish = WAA_fish,
                                    MatAA = MatAA,
                                    fish_sel = fish_sel,
                                    Movement = Movement,
                                    f_ref_pt = proj_inputs[[i]]$f_ref_pt,
                                    b_ref_pt = proj_inputs[[i]]$b_ref_pt,
                                    HCR_function = HCR_function,
                                    recruitment_opt = "inv_gauss",
                                    fmort_opt = proj_inputs[[i]]$fmort_opt,
                                    t_spawn = t_spawn
    )

    all_scenarios_ssb[,,,sim,i] <- out$proj_SSB
    all_scenarios_catch[,,,,,sim,i] <- out$proj_Catch
    all_scenarios_f[,,sim,i] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used
  } # end sim loop
  print(i)
} # end i loop

#### SSB Projections ---------------------------------------------------------

# Get historical SSB
historical <- reshape2::melt(array(rep(francis_model$rep$SSB, n_sims),
                                   dim = c(n_pop, n_regions, length(francis_data$years), n_sims))) %>%
  mutate(Year = Var3 + 1976,
         Scenario = "FABC (F40)",  # or change to match the scenarios you're plotting
         Type = "Historical") %>%
  rename(Pop = Var1, Region = Var2, Simulation = Var4, SSB = value)

# Get all scenario projections
scenarios <- reshape2::melt(all_scenarios_ssb) %>%
  mutate(Year = Var3 + 2023,
         Scenario = case_when(
           Var5 == 1 ~ "S1: FABC (F40)",
           Var5 == 2 ~ "S2: FABC Ratio",
           Var5 == 3 ~ "S3: F Last 5 Years",
           Var5 == 4 ~ "S4: F60 SPR",
           Var5 == 5 ~ "S5: No Fishing",
           Var5 == 6 ~ "S6: FOFL",
           Var5 == 7 ~ "S7: FABC -> FOFL",
           TRUE ~ paste("Scenario", Var4)
         ),
         Type = "Projection") %>%
  rename(Pop = Var1, Region = Var2, Simulation = Var4, SSB = value)

# expand historical SSB for plotting
scenarios_unique <- unique(scenarios$Scenario)
historical_expanded <- historical[rep(1:nrow(historical), times = length(scenarios_unique)), ]
historical_expanded$Scenario <- rep(scenarios_unique, each = nrow(historical))

# combine
combined_ssb <- bind_rows(historical_expanded, scenarios)

# Plot
png(here("vignettes", "figures", "o_projall_ssb.png"), width = 1000, height = 800)
combined_ssb %>%
  ggplot(aes(x = Year, y = SSB, group = interaction(Scenario, Simulation), color = Type)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~Scenario, scales = 'free') +
  geom_hline(yintercept = b40, lty = 2) +
  scale_color_manual(values = c("Historical" = "black", "Projection" = "blue")) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')
dev.off()

png(here("vignettes", "figures", "o_projzoom_ssb.png"), width = 1000, height = 800)
combined_ssb %>%
  filter(Year > 2024) %>%
  group_by(Year, Scenario, Type) %>%
  summarize(lwr = quantile(SSB, 0.025),
            upr = quantile(SSB, 0.975),
            SSB = mean(SSB)) %>%
  ggplot(aes(x = Year, y = SSB, ymin = lwr, ymax = upr)) +
  geom_line(alpha = 1, lwd = 1.3) +
  geom_ribbon(color = NA, alpha = 0.3) +
  facet_wrap(~Scenario) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_hline(yintercept = b40, lty = 2) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')
dev.off()

combined_ssb %>%
  filter(Year > 2024) %>%
  group_by(Year, Scenario, Type) %>%
  summarize(lwr = quantile(SSB, 0.025),
            upr = quantile(SSB, 0.975),
            SSB = mean(SSB)) %>%
  select(-c(Type, lwr, upr)) %>%
  pivot_wider(names_from = Scenario, values_from = SSB)

#### Catch Projections -------------------------------------------------------

# Get historical catch
historical <- reshape2::melt(array(rep(francis_data$ObsCatch, n_sims),
                                   dim = c(n_pop, n_regions, length(francis_data$years), francis_data$n_fish_fleets, n_sims))) %>%
  mutate(Year = Var3 + 1976,
         Scenario = "FABC (F40)",  # or change to match the scenarios you're plotting
         Type = "Historical") %>%
  rename(Pop = Var1, Region = Var2, Simulation = Var5, Fleet = Var4, Catch = value) %>%
  select(-Var3)

historical$Catch[is.na(historical$Catch)] <- 0

# Get all scenario projections
scenarios <- reshape2::melt(all_scenarios_catch) %>%
  mutate(Year = Var3 + 2023,
         Scenario = case_when(
           Var7 == 1 ~ "S1: FABC (F40)",
           Var7 == 2 ~ "S2: FABC Ratio",
           Var7 == 3 ~ "S3: F Last 5 Years",
           Var7 == 4 ~ "S4: F60 SPR",
           Var7 == 5 ~ "S5: No Fishing",
           Var7 == 6 ~ "S6: FOFL",
           Var7 == 7 ~ "S7: FABC -> FOFL",
           TRUE ~ paste("Scenario", Var6)
         ),
         Type = "Projection") %>%
  rename(Pop = Var1, Region = Var2, Simulation = Var6, Catch = value, Fleet = Var5, Seas = Var4) %>%
  select(-c(Var3, Var7))

# expand historical SSB for plotting
scenarios_unique <- unique(scenarios$Scenario)
historical_expanded <- historical[rep(1:nrow(historical), times = length(scenarios_unique)), ]
historical_expanded$Scenario <- rep(scenarios_unique, each = nrow(historical))

# combine
combined_cat <- bind_rows(historical_expanded, scenarios)

# Plot
png(here("vignettes", "figures", "o_projall_catch.png"), width = 1000, height = 800)
combined_cat %>%
  group_by(Year, Scenario, Simulation, Type, Region) %>%
  summarize(Catch = sum(Catch)) %>%
  ggplot(aes(x = Year, y = Catch, group = interaction(Scenario, Simulation), color = Type)) +
  geom_line(alpha = 1) +
  facet_wrap(~Scenario) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = c("Historical" = "black", "Projection" = "blue")) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')
dev.off()

png(here("vignettes", "figures", "o_projzoom_catch.png"), width = 1000, height = 800)
combined_cat %>%
  filter(Year > 2024) %>%
  group_by(Year, Scenario, Simulation, Type, Region) %>%
  summarize(Catch = sum(Catch)) %>%
  group_by(Year, Scenario, Type) %>%
  summarize(lwr = quantile(Catch, 0.025),
            upr = quantile(Catch, 0.975),
            Catch = mean(Catch)) %>%
  ggplot(aes(x = Year, y = Catch, ymin = lwr, ymax = upr)) +
  geom_line(alpha = 1, lwd = 1.3) +
  geom_ribbon(color = NA, alpha = 0.3) +
  facet_wrap(~Scenario) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')
dev.off()

# Catch advice under F40
png(here("vignettes", "figures", "o_projzoom_catch_f40.png"), width = 1000, height = 800)
combined_cat %>%
  filter(Year > 2024, Scenario == "S1: FABC (F40)") %>%
  group_by(Year, Scenario, Simulation, Type, Region) %>%
  summarize(Catch = sum(Catch)) %>%
  group_by(Year, Scenario, Type) %>%
  summarize(lwr = quantile(Catch, 0.025),
            upr = quantile(Catch, 0.975),
            Catch = mean(Catch)) %>%
  ggplot(aes(x = Year, y = Catch, ymin = lwr, ymax = upr)) +
  geom_line(alpha = 1, lwd = 1.3) +
  geom_ribbon(color = NA, alpha = 0.3) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')
dev.off()

combined_cat %>%
  filter(Year > 2024) %>%
  group_by(Year, Scenario, Simulation, Type, Region) %>%
  summarize(Catch = sum(Catch)) %>%
  group_by(Year, Scenario, Type) %>%
  summarize(lwr = quantile(Catch, 0.025),
            upr = quantile(Catch, 0.975),
            Catch = mean(Catch)) %>%
  select(-c(Type, lwr, upr)) %>%
  pivot_wider(names_from = Scenario, values_from = Catch)

#### F Projections -----------------------------------------------------------

# Get historical catch
historical <- reshape2::melt(array(rep(as.vector(apply(francis_model$rep$Fmort, c(1,2), sum)), n_sims),
                                   dim = c(n_regions, length(francis_data$years), n_sims))) %>%
  mutate(Year = Var2 + 1976,
         Scenario = "FABC (F40)",  # or change to match the scenarios you're plotting
         Type = "Historical") %>%
  rename(Region = Var1, Simulation = Var3, Fmort = value)

# Get all scenario projections
scenarios <- reshape2::melt(all_scenarios_f) %>%
  mutate(Year = Var2 + 2023,
         Scenario = case_when(
           Var4 == 1 ~ "S1: FABC (F40)",
           Var4 == 2 ~ "S2: FABC Ratio",
           Var4 == 3 ~ "S3: F Last 5 Years",
           Var4 == 4 ~ "S4: F60 SPR",
           Var4 == 5 ~ "S5: No Fishing",
           Var4 == 6 ~ "S6: FOFL",
           Var4 == 7 ~ "S7: FABC -> FOFL",
           TRUE ~ paste("Scenario", Var4)
         ),
         Type = "Projection") %>%
  rename(Region = Var1, Simulation = Var3, Fmort = value)

# expand historical SSB for plotting
scenarios_unique <- unique(scenarios$Scenario)
historical_expanded <- historical[rep(1:nrow(historical), times = length(scenarios_unique)), ]
historical_expanded$Scenario <- rep(scenarios_unique, each = nrow(historical))

# combine
combined_fmort <- bind_rows(historical_expanded, scenarios)

# Plot
png(here("vignettes", "figures", "o_projall_f.png"), width = 1000, height = 800)
combined_fmort %>%
  ggplot(aes(x = Year, y = Fmort, group = interaction(Scenario, Simulation), color = Type)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~Scenario, scales = 'free') +
  scale_color_manual(values = c("Historical" = "black", "Projection" = "blue")) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none') +
  labs(y = 'Fully Selected Fishing Mortality')
dev.off()

png(here("vignettes", "figures", "o_projzoom_f.png"), width = 1000, height = 800)
combined_fmort %>%
  filter(Year > 2024) %>%
  group_by(Year, Scenario, Type) %>%
  summarize(lwr = quantile(Fmort, 0.025),
            upr = quantile(Fmort, 0.975),
            Fmort = mean(Fmort)) %>%
  ggplot(aes(x = Year, y = Fmort, ymin = lwr, ymax = upr)) +
  geom_line(alpha = 1, lwd = 1.3) +
  geom_ribbon(color = NA, alpha = 0.3) +
  facet_wrap(~Scenario) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none') +
  labs(y = 'Fully Selected Fishing Mortality')
dev.off()

combined_fmort %>%
  filter(Year > 2024) %>%
  group_by(Year, Scenario, Simulation, Type, Region) %>%
  group_by(Year, Scenario, Type) %>%
  summarize(lwr = quantile(Fmort, 0.025),
            upr = quantile(Fmort, 0.975),
            Fmort = mean(Fmort)) %>%
  select(-c(Type, lwr, upr)) %>%
  pivot_wider(names_from = Scenario, values_from = Fmort)

