# Purpose: To demonstrate the use of SPoRC in a single region context for BSAI pollock
# Creator: Matthew LH. Cheng
# Date Created: 5/1/25

# Set up ------------------------------------------------------------------
# unloadNamespace("SPoRC")
library(here)
library(tidyverse)
library(RTMB)
# library(SPoRC)
library(ebswp)
devtools::load_all(here("R"))

# Read in data
in_dat <- ebswp::read_dat(here::here("dev", "dev_data", "pm_24_db.dat"))
rep <- ebswp::read_rep(here::here("dev", "dev_data", "pm.rep"))

# Reformat data list to fit SPoRC
ebswp_SPoRC_data <- list()
# Dimensions
ebswp_SPoRC_data$years <- seq(in_dat$styr, in_dat$endyr)
ebswp_SPoRC_data$ages <- 1:in_dat$nages

# Biologicals
ebswp_SPoRC_data$WAA <- array(0, dim = c(1, length(ebswp_SPoRC_data$years), length(ebswp_SPoRC_data$ages), 1))
ebswp_SPoRC_data$MatAA <- array(0, dim = c(1, length(ebswp_SPoRC_data$years), length(ebswp_SPoRC_data$ages), 1))
for(a in ebswp_SPoRC_data$ages) {
  for(y in 1:length(ebswp_SPoRC_data$years )) {
    ebswp_SPoRC_data$WAA[,y,a,] <- in_dat$wt_ssb[y,a]
    ebswp_SPoRC_data$MatAA[,y,a,] <- in_dat$p_mature[a]
  } # end y loop
} # end a loop

ebswp_SPoRC_data$AgeingError <- as.matrix(in_dat$age_err)

# Catch stuff
ebswp_SPoRC_data$ObsCatch <- array(in_dat$obs_catch, dim = c(1, length(in_dat$obs_catch), 1))
ebswp_SPoRC_data$UseCatch <- array(1, dim = c(1, length(in_dat$obs_catch), 1))

# Fishery Indices and Compositions
# set up observed age compositions for fishery
yrs <- in_dat$styr:in_dat$endyr
idx <- c(1:61)[seq(1964, 2024, 1) %in% in_dat$yrs_fsh_data]
oac_tmp <- (array(NA_real_, dim = c(1, nrow(in_dat$wt_fsh), ncol(in_dat$wt_fsh), 1, 1)))
for (i in 1:length(idx)) oac_tmp[1, idx[i], 1:15, 1, 1] <- in_dat$oac_fsh_data[i, 1:15]

# set up fishery age comp indicators
use_tmp <- (array(0, dim = c(1, nrow(in_dat$wt_fsh), 1)))
for (i in 1:length(idx))  use_tmp[1, idx[i], 1] <- 1

# set up fishery ISS
iss_tmp <- (array(0, dim = c(1, nrow(in_dat$wt_fsh), 1, 1)))
for (i in 1:length(idx)) iss_tmp[1, idx[i], 1, 1] <- in_dat$sam_fsh[i]

ebswp_SPoRC_data$ObsFishAgeComps <- oac_tmp
ebswp_SPoRC_data$UseFishAgeComps <- use_tmp
ebswp_SPoRC_data$ISS_FishAgeComps <- iss_tmp

# Survey Indicies and Compositions
# Define array dimensions to set up data inputs for survey indices and age compositions
srv_tmp <- (array(NA_real_, dim = c(1, nrow(in_dat$wt_fsh), 3)))
srvse_tmp <- (array(NA_real_, dim = c(1, nrow(in_dat$wt_fsh), 3)))
oac_tmp <- (array(NA_real_, dim = c(1, nrow(in_dat$wt_fsh), ncol(in_dat$wt_fsh), 1, 3)))
use_tmp <- (array(0, dim = c(1, nrow(in_dat$wt_fsh), 3)))
iss_tmp <- (array(0, dim = c(1, nrow(in_dat$wt_fsh), 1, 3)))
use_oac <- (array(0, dim = c(1, nrow(in_dat$wt_fsh), 3)))

# Set up bottom trawl survey data
idx <- c(1:61)[seq(1964, 2024, 1) %in% in_dat$yrs_bts]
for (i in 1:length(idx)) {
  srv_tmp[1, idx[i], 1] <- in_dat$ob_bts[i]
  srvse_tmp[1, idx[i], 1] <- in_dat$ob_bts_std[i]
  use_tmp[1, idx[i], 1] <- 1
  use_oac[1, idx[i], 1] <- 1
  oac_tmp[1, idx[i], 1:15, 1, 1] <- in_dat$oac_bts[i, 1:15]
  iss_tmp[1, idx[i], 1, 1] <- in_dat$sam_bts[i]
}

# Set up ATS survey data
idx <- c(1:61)[seq(1964, 2024, 1) %in% in_dat$yrs_ats]
for (i in 1:length(idx)) {
  srv_tmp[1, idx[i], 2] <- in_dat$ob_ats[i]
  srvse_tmp[1, idx[i], 2] <- in_dat$ob_ats_std[i]
  use_tmp[1, idx[i], 2] <- 1
  use_oac[1, idx[i], 2] <- 1
  oac_tmp[1, idx[i], 1:15, 1, 2] <- in_dat$oac_ats[i, 1:15]
  iss_tmp[1, idx[i], 1, 2] <- in_dat$sam_ats[i]
}

# Set up AVO survey data
idx <- c(1:61)[seq(1964, 2024, 1) %in% in_dat$yrs_avo]
for (i in 1:length(idx)) {
  srv_tmp[1, idx[i], 3] <- in_dat$ob_avo[i]
  srvse_tmp[1, idx[i], 3] <- in_dat$ob_avo_std[i]
  use_tmp[1, idx[i], 3] <- 1
  use_oac[1, idx[i], 3] <- 0
}

ebswp_SPoRC_data$ObsSrvIdx <- srv_tmp
ebswp_SPoRC_data$ObsSrvIdx_SE <- srvse_tmp / srv_tmp
ebswp_SPoRC_data$UseSrvIdx <- use_tmp
ebswp_SPoRC_data$ObsSrvAgeComps <- oac_tmp
ebswp_SPoRC_data$ISS_SrvAgeComps <- iss_tmp
ebswp_SPoRC_data$UseSrvAgeComps <- use_oac

# Setup fishery index
ebswp_SPoRC_data$ObsFishIdx <- array(NA, dim = c(1, length(ebswp_SPoRC_data$years), 1))
ebswp_SPoRC_data$ObsFishIdx_SE <- array(NA, dim = c(1, length(ebswp_SPoRC_data$years), 1))
ebswp_SPoRC_data$UseFishIdx <- array(0, dim = c(1, length(ebswp_SPoRC_data$years), 1))
ebswp_SPoRC_data$ObsFishIdx[,2:13,1] <- in_dat$obs_cpue
ebswp_SPoRC_data$ObsFishIdx_SE[,2:13,1] <- in_dat$obs_cpue_std / in_dat$obs_cpue
ebswp_SPoRC_data$UseFishIdx[,2:13,1] <- 1

# Save report values from actual assessment
ebswp_SPoRC_data$SSB <- rep$SSB
ebswp_SPoRC_data$R <- rep$R

# write out data
sgl_rg_ebswp_data <- ebswp_SPoRC_data
usethis::use_data(sgl_rg_ebswp_data, internal = FALSE, overwrite = TRUE)

# Prepare Data and Inputs -------------------------------------------------
{
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
    ln_sigmaC = array(log(0.05), dim = c(1, length(input_list$data$years), 1))
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
    ObsFishLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), length(input_list$data$lens), 1, 1)),
    UseFishLenComps = array(0, dim = c(1, length(input_list$data$years), 1)),
    ISS_FishLenComps = NULL,

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
    ObsSrvLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), length(input_list$data$lens), 1, 3)),
    UseSrvLenComps = array(0, dim = c(1, length(input_list$data$years), 3)),
    ISS_SrvLenComps = NULL,

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
                                       input_list$data$n_sexes,
                                       input_list$data$n_srv_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years),
                                       input_list$data$n_sexes,
                                       input_list$data$n_srv_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets))
  )
}

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
                              # random = NULL,
                              newton_loops = 3,
                              silent = FALSE
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

