# Purpose: To demonstrate the use of SPoRC in a single region context for BSAI pollock
# Creator: Matthew LH. Cheng
# Date Created: 5/1/25

# Set up ------------------------------------------------------------------
# unloadNamespace("SPoRC")
library(here)
library(tidyverse)
library(ebswp)
devtools::load_all(here("R"))

# Prepare Data and Inputs -------------------------------------------------

## Initialize model dimensions and data list----
input_list <- Setup_Mod_Dim(
  years = sgl_rg_ebswp_data$years,
  # vector of years
  ages = sgl_rg_ebswp_data$ages,
  # vector of ages
  lens = NA,
  # number of lengths
  n_regions = 1,
  # number of regions
  n_sexes = 1,
  # number of sexes
  n_fish_fleets = 1,
  # number of fishery fleets
  n_srv_fleets = 3, # number of survey fleets
  # number of seasons
  n_seas = sgl_rg_ebswp_data$n_seas,
  verbose = TRUE
)

inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(
  input_list = input_list,

  # Model options
  do_rec_bias_ramp = 0,
  # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
  sigmaR_switch = 1,
  # when to switch from early to late sigmaR (switch in first year)
  ln_sigmaR = log(c(1, 1)),
  # Starting values for early and late sigmaR
  rec_model = "bh_rec",
  # recruitment model
  steepness_h = inv_steepness(0.623013),
  h_spec = "fix",
  # fixing steepness
  sigmaR_spec = "fix",
  # fix early sigmaR and late sigmaR
  init_age_strc = 1,
  ln_global_R0 = 10,
  t_spawn = 0.25,
  equil_init_age_strc = 2
  # starting value for r0
)

# Setup a fixed natural mortality array for use
fix_natmort <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), 1))
fix_natmort[,,1,] <- 0.9 # age 1 M
fix_natmort[,,2,] <- 0.45 # age 2 M
fix_natmort[,,-c(1,2),] <- 0.3 # age 3+ M

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,

  # Data inputs
  WAA = sgl_rg_ebswp_data$WAA,
  MatAA = sgl_rg_ebswp_data$MatAA,

  # Model options
  # mean and sd for M prior
  fit_lengths = 0,
  # don't fit length compositions
  M_spec = "fix",
  # fixing natural mortality
  Fixed_natmort = fix_natmort
)

# Setup movement stuff (using defaults for other stuff)
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = NA,
  do_recruits_move = 0
)

input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,

  # Data inputs
  ObsCatch = sgl_rg_ebswp_data$ObsCatch,
  UseCatch = sgl_rg_ebswp_data$UseCatch,

  # Model options
  Use_F_pen = 1,
  # whether to use f penalty, == 0 don't use, == 1 use
  sigmaC_spec = "fix",
  # fixing catch standard deviation
  ln_sigmaC = array(log(0.05), dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1))
  # starting / fixed value for catch standard deviation
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  # data inputs
  ObsFishIdx = sgl_rg_ebswp_data$ObsFishIdx,
  ObsFishIdx_SE = sgl_rg_ebswp_data$ObsFishIdx_SE,
  UseFishIdx = sgl_rg_ebswp_data$UseFishIdx,
  ObsFishAgeComps = sgl_rg_ebswp_data$ObsFishAgeComps,
  UseFishAgeComps = sgl_rg_ebswp_data$UseFishAgeComps,
  ISS_FishAgeComps = sgl_rg_ebswp_data$ISS_FishAgeComps,
  ObsFishLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, 1)),
  UseFishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1)),
  ISS_FishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, 1)),

  # Model options
  fish_idx_type = c("biom"),
  # indices for fishery
  FishAgeComps_LikeType = c("Multinomial"),
  # age comp likelihoods for fishery fleet
  FishLenComps_LikeType = c("none"),
  # length comp likelihoods for fishery
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
  # age comp structure for fishery
  FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  # length comp structure for fishery
)

# Setup survey indices and compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,

  # data inputs
  ObsSrvIdx = sgl_rg_ebswp_data$ObsSrvIdx,
  ObsSrvIdx_SE = sgl_rg_ebswp_data$ObsSrvIdx_SE,
  UseSrvIdx = sgl_rg_ebswp_data$UseSrvIdx,
  ObsSrvAgeComps = sgl_rg_ebswp_data$ObsSrvAgeComps,
  ISS_SrvAgeComps = sgl_rg_ebswp_data$ISS_SrvAgeComps,
  UseSrvAgeComps = sgl_rg_ebswp_data$UseSrvAgeComps,
  ObsSrvLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, 3)),
  UseSrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, 3)),
  ISS_SrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, 3)),

  # Model options
  srv_idx_type = c("biom", "biom", "biom"),
  # abundance and biomass for survey fleet 1, 2, and 3
  SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "Multinomial"),
  # survey age composition likelihood for survey fleet 1, 2, and 3
  SrvLenComps_LikeType = c("none", "none", "none"),
  #  survey length composition likelihood for survey fleet 1, 2, and 3
  SrvAgeComps_Type = c(
    "agg_Year_1-terminal_Fleet_1",
    "agg_Year_1-terminal_Fleet_2",
    "none_Year_1-terminal_Fleet_3"
  ),
  # survey age comp type

  SrvLenComps_Type = c(
    "none_Year_1-terminal_Fleet_1",
    "none_Year_1-terminal_Fleet_2",
    "none_Year_1-terminal_Fleet_3"
  )
  # survey length comp type
)


# Setup fishery selectivity and catchability
input_list <- Setup_Mod_Fishsel_and_Q(

  input_list = input_list,

  # Model options
  # fishery selectivity, whether continuous time-varying
  cont_tv_fish_sel = c("2dar1_Fleet_1"),
  fishsel_pe_pars_spec = "fix", # doing penalized likelihood for selex devs
  fish_sel_devs_spec = "est_all", # estimating all sel devs
  corr_opt_semipar = "corr_zero_y_b", # making sure 2d correaltions are 0, collapses to a simple iid case
  # fishery selectivity blocks
  fish_sel_blocks = c("none_Fleet_1"),
  # fishery selectivity form
  fish_sel_model = c("logist1_Fleet_1"),
  # fishery catchability blocks
  fish_q_blocks = c("none_Fleet_1"),
  # whether to estiamte all fixed effects for fishery selectivity
  fish_fixed_sel_pars_spec = c("est_all"),
  # whether to estiamte all fixed effects for fishery catchability
  fish_q_spec = c("est_all")
)

# Setup survey selectivity and catchability
input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,

  # Model options
  # survey selectivity, whether continuous time-varying
  cont_tv_srv_sel = c("iid_Fleet_1", "2dar1_Fleet_2", "2dar1_Fleet_3"),
  srvsel_pe_pars_spec = c("fix", "fix", "fix"), # penalize survey selex devs
  srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2"), # estimating all srv selex devs
  corr_opt_semipar = c(NA, "corr_zero_y_b", "corr_zero_y_b"), # setting corelations at 0, so 2dar1 collapses to simple iid semi-parametric devs

  # survey selectivity blocks
  srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  # survey selectivity form
  srv_sel_model = c(
    "logist1_Fleet_1",
    "logist1_Fleet_2",
    "logist1_Fleet_3"
  ),
  # survey catchability blocks
  srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  # whether to estiamte all fixed effects for survey selectivity
  srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2"),
  # whether to estiamte all fixed effects for survey catchability
  srv_q_spec = c("est_all", "est_all", "est_all")
)

# Setup tagging stuff
input_list <- Setup_Mod_Tagging(input_list = input_list, UseTagging = 0)

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_seas,
                                     input_list$data$n_sexes,
                                     input_list$data$n_srv_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_seas,
                                     input_list$data$n_sexes,
                                     input_list$data$n_srv_fleets)),
  Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_seas,
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets)),
  Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_seas,
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets))
)


# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# selex sigma to fix at, given penalized likelihood
parameters$fishsel_pe_pars[,4,,] <- log(0.075) # fishery selex variance
parameters$srvsel_pe_pars[,1:2,,1] <- log(0.075) # survey BTS - a50 and delta variance
parameters$srvsel_pe_pars[,4,,2] <- log(0.15) # survey ATS and ato variance

# Fit model
ebswp_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              newton_loops = 3,
                              silent = FALSE, do_optim = T
                              )

ebswp_rtmb_model$sdrep <- RTMB::sdreport(ebswp_rtmb_model)

sdrep <- ebswp_rtmb_model$sdrep
rep <- ebswp_rtmb_model$rep

# Plot! -------------------------------------------------------------------
# Get recruitment
rec_series <- reshape2::melt((ebswp_rtmb_model$rep$Rec)) %>%
  mutate(se = ebswp_rtmb_model$sdrep$sd[names(ebswp_rtmb_model$sdrep$value) == 'log(Rec)'] * t(ebswp_rtmb_model$rep$Rec))
rec_series$Par <- "Recruitment"

# Get SSB time-series
ssb_series <- reshape2::melt((ebswp_rtmb_model$rep$SSB)) %>%
  mutate(se = ebswp_rtmb_model$sdrep$sd[names(ebswp_rtmb_model$sdrep$value) == 'log(SSB)'] * t(ebswp_rtmb_model$rep$SSB))
ssb_series$Par <- "Spawning Biomass"

# bind together
ts_df <- rbind(ssb_series,rec_series) %>%
  dplyr::rename(Region = Var1, Year = Var2) %>%
  dplyr::mutate(Year = Year + 1963, type = 'SPoRC')

# Get actual assessment results
ssb_ass <- data.frame(
  Region = 1,
  Year = 1964:2024,
  value = sgl_rg_ebswp_data$SSB[,2],
  se = sgl_rg_ebswp_data$SSB[,3],
  Par = 'Spawning Biomass',
  type = '2024 Pollock Assessment'
)

# recruitment
rec_ass <- data.frame(
  Region = 1,
  Year = 1964:2024,
  value = sgl_rg_ebswp_data$R[,2],
  se = sgl_rg_ebswp_data$R[,3],
  Par = 'Recruitment',
  type = '2024 Pollock Assessment'
)

# bind
ts_df <- ts_df %>% bind_rows(ssb_ass, rec_ass)

# plot!
png(here("vignettes", "figures", "f_ebs_pol_ts_comparison.png"), width = 1000, height = 500)
ggplot(ts_df, aes(x = Year, y = value, ymin = value - (1.96 * se),
                  ymax = value + (1.96 * se), color = type, fill = type)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~Par, scales = 'free') +
  geom_ribbon(alpha = 0.3, color = NA) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(y = "Value")  +
  theme_bw(base_size = 13) +
  ylim(0, NA) +
  labs(x = 'Year', y = 'Value', color = 'Type', fill = 'Type')
dev.off()

png(here("vignettes", "figures", "f_ebs_pol_fishsel_comparison.png"), width = 1000, height = 500)
reshape2::melt(ebswp_rtmb_model$rep$fish_sel) %>%
  mutate(value = value/max(value)) %>%
  rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Fleet = Var5) %>%
  filter(Age %in% 3:11) %>%
  ggplot(aes(x = Year + 1963, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Age) +
  theme_bw(base_size = 15) +
  labs(x = 'Year', y = 'Relative Selectivity')
dev.off()

