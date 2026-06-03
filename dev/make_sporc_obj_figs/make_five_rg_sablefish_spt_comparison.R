# Purpose: To bridge to the spatial assessment in Cheng and Marsh et al. 2025  using SPoRC
# Creator: Matthew LH. Cheng
# Date Created: 2/5/25

# Set up ------------------------------------------------------------------
unloadNamespace("SPoRC")
library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)
devtools::load_all(here("R"))

data("mlt_rg_sable_data")

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(mlt_rg_sable_data$years),
                            # vector of years (1 - 62)
                            ages = 1:length(mlt_rg_sable_data$ages),
                            # vector of ages (1 - 30)
                            lens = mlt_rg_sable_data$lens,
                            # number of lengths (41 - 99)
                            n_regions = mlt_rg_sable_data$n_regions,
                            # number of regions (5)
                            n_sexes = mlt_rg_sable_data$n_sexes,
                            # number of sexes (2)
                            n_fish_fleets = mlt_rg_sable_data$n_fish_fleets,
                            # number of fishery fleet (2)
                            n_srv_fleets = mlt_rg_sable_data$n_srv_fleets,
                            # number of survey fleets (2)
                            n_seas = mlt_rg_sable_data$n_seas,
                            # number of seasons
                            n_pop = mlt_rg_sable_data$n_pop,
                            natal_region = mlt_rg_sable_data$natal_region,
                            # population stuff
                            verbose = TRUE
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            do_rec_bias_ramp = 0, # not using bias ramp
                            sigmaR_switch = 16, # switch to using late sigma in year 16
                            dont_est_recdev_last = 1, # don't estimate last rec dev
                            # Model options
                            rec_model = "mean_rec", # recruitment model
                            sigmaR_spec = "fix", # fixing
                            InitDevs_spec = "est_shared_r",
                            # initial deviations are shared across regions,
                            # but recruitment deviations are region specific
                            ln_sigmaR = array(log(c(0.4, 1.2)), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
                            # values to fix sigmaR at, or starting values
                            ln_global_R0 = log(20),
                            use_rec_region_prop_prior = 1,
                            rec_region_prop_prior = data.frame(pop = 1, alpha = I(list(rep(3, input_list$data$n_regions)))),
                            # starting value for global R0
                            rec_region_prop_pars = array(c(0.2, 0.2, 0.2, 0.2),
                                                         dim = c(input_list$data$n_pop, input_list$data$n_regions - 1))
                            # starting value for R0 proportions in multinomial logit space
)

# Setup biological stuff (using defaults for other stuff)
input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = mlt_rg_sable_data$WAA, # weight at age
                                    MatAA = mlt_rg_sable_data$MatAA, # maturity at age
                                    AgeingError = mlt_rg_sable_data$AgeingError,
                                    # ageing error matrix
                                    fit_lengths = 1, # fitting lengths
                                    SizeAgeTrans = mlt_rg_sable_data$SizeAgeTrans,
                                    # size age transition matrix
                                    M_spec = "fix", # fix natural mortality
                                    Fixed_natmort = array(0.104884, dim = c(mlt_rg_sable_data$n_pop,
                                                                            mlt_rg_sable_data$n_regions,
                                                                            length(mlt_rg_sable_data$years),
                                                                            length(mlt_rg_sable_data$ages),
                                                                            mlt_rg_sable_data$n_sexes))
                                    # values to fix natural mortality at
)

# setting up movement parameterization
Movement_prior <- expand.grid(
  pop = 1, # populations
  region_from = 1:5, # regions
  year = 1, # penalize first year since no blocks
  seas = 1,
  age = c(6,7,16), # age blocks
  sex = 1, # sex
  alpha = I(list(rep(3, 5))) # prior alpha to each row
)

input_list <- Setup_Mod_Movement(input_list = input_list,
                                 # Model options
                                 Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)),
                                 # estimating movement in 3 age blocks
                                 # (ages 1-6, ages 7-15, ages 16-30)
                                 Movement_popblk_spec = 'constant', # population-invariant movement
                                 Movement_yearblk_spec = "constant", # time-invariant movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 Movement_seasblk_spec = 'constant', # seasonal blocks
                                 do_recruits_move = 0, # recruits do not move
                                 use_fixed_movement = 0, # estimating movement
                                 Use_Movement_Prior = 1, # priors used for movement
                                 Movement_prior = Movement_prior
                                 # vague prior to penalize movement away from the extremes
)

# setting up tagging parameterization

# setup tagging priors
tag_prior <- data.frame(
  region = 1,
  block = c(1,2),
  fleet = 1,
  mu = NA, # no mean, since symmetric beta
  sd = 5, # sd = 5
  type = 0 # symmetric beta
)

input_list <- Setup_Mod_Tagging(input_list = input_list,
                                use_conv_fish_tagging = c(1,0), # using tagging data for fixed gear
                                conv_tag_max_liberty = 15, # maximum number of years to track a cohort

                                # Data Inputs
                                conv_tag_release_indicator = mlt_rg_sable_data$conv_tag_release_indicator,
                                # tag release indicator (first col = tag region,
                                # second col = tag year),
                                # total number of rows = number of tagged cohorts
                                conv_tagged_fish = mlt_rg_sable_data$conv_tagged_fish, # Released fish
                                # dimensioned by total number of tagged cohorts, (implicitly
                                # tracks the release year and region), pop, age, and sex
                                obs_conv_tag_fish_recap = mlt_rg_sable_data$obs_conv_tag_fish_recap,
                                # dimensioned by max tag liberty, tagged cohorts, pop, regions,
                                # ages, and sexes

                                # Model options
                                conv_fish_tag_like = "Multinomial_Release", # Negative Binomial
                                conv_tag_mixing_period = 2, # Don't fit tagging until release year + 1
                                conv_tag_t_tagging = 0.5, # tagging happens midway through the year,
                                # movement does not occur within that year
                                use_conv_tag_fishrep_prior = 1, # tag reporting rate priors are used
                                conv_tag_fishrep_prior = tag_prior,
                                conv_tag_age_pool = as.list(1:30), # whether or
                                # not to pool tagging data when fitting (for computational cost)
                                conv_tag_sex_pool = list(c(1:2)), # whether or not to pool
                                # sex-specific data when fitting
                                init_conv_tag_mort_spec = "fix", # fixing initial tag mortality
                                conv_tag_shed_spec = "fix", # fixing chronic shedding
                                conv_tagrep_spec = "est_shared_r_f", # tag reporting rates are
                                # not region specific
                                # Time blocks for tag reporting rates
                                conv_tag_fish_reporting_blocks = c(
                                  apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
                                    paste0("Block_1_Year_1-35_Region_", x[1], "_Fleet_", x[2])),
                                  apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
                                    paste0("Block_2_Year_36-terminal_Region_", x[1], "_Fleet_", x[2]))
                                ),
                                conv_fish_tag_attr = 'p_a_s',
                                # Specify starting values or fixing values
                                ln_init_conv_tag_mort = log(0.1), # fixing initial tag mortality
                                ln_conv_tag_shed = log(0.02),  # fixing tag shedding
                                ln_conv_fish_tag_theta = log(0.5),
                                # starting value for tagging overdispersion
                                conv_tag_fish_reporting_pars = array(log(0.2 / (1-0.2)), dim = c(input_list$data$n_regions, 2, input_list$data$n_fish_fleets))
                                # starting values for tag reporting pars

)

# setting up catch data
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = mlt_rg_sable_data$ObsCatch,
                                    UseCatch = mlt_rg_sable_data$UseCatch,
                                    # Model options
                                    Use_F_pen = 1,
                                    # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix',
                                    ln_sigmaC =
                                      array(log(0.05), dim = c(input_list$data$n_regions,
                                                               length(input_list$data$years),
                                                               input_list$data$n_seas,
                                                               input_list$data$n_fish_fleets)),
                                    # fixing catch sd at small value
                                    ln_F_mean = array(-2, dim = c(input_list$data$n_regions,
                                                                  input_list$data$n_seas,
                                                                  input_list$data$n_fish_fleets))
                                    # some starting values for fishing mortality
)

# Fishery Indices and Compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = mlt_rg_sable_data$ObsFishIdx,
                                          ObsFishIdx_SE = mlt_rg_sable_data$ObsFishIdx_SE,
                                          UseFishIdx =  mlt_rg_sable_data$UseFishIdx,
                                          ObsFishAgeComps = mlt_rg_sable_data$ObsFishAgeComps,
                                          UseFishAgeComps = mlt_rg_sable_data$UseFishAgeComps,
                                          ISS_FishAgeComps = mlt_rg_sable_data$ISS_FishAgeComps,
                                          ObsFishLenComps = mlt_rg_sable_data$ObsFishLenComps,
                                          UseFishLenComps = mlt_rg_sable_data$UseFishLenComps,
                                          ISS_FishLenComps = mlt_rg_sable_data$ISS_FishLenComps,

                                          # Model options
                                          fish_idx_type = c("none", "none"),
                                          # fishery indices not used
                                          FishAgeComps_LikeType =
                                            c("Multinomial", "none"),
                                          # age comp likelihoods for fishery fleet 1 and 2
                                          FishLenComps_LikeType =
                                            c("Multinomial", "Multinomial"),
                                          # length comp likelihoods for fishery fleet 1 and 2
                                          FishAgeComps_Type =
                                            c("spltRjntS_Year_1-terminal_Fleet_1",
                                              "none_Year_1-terminal_Fleet_2"),
                                          # age comp structure for fishery fleet 1 and 2
                                          FishLenComps_Type =
                                            c("spltRjntS_Year_1-terminal_Fleet_1",
                                              "spltRjntS_Year_1-terminal_Fleet_2")
                                          # length comp structure for fishery fleet 1 and 2
)

# Survey Indices and Compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         # data inputs
                                         ObsSrvIdx = mlt_rg_sable_data$ObsSrvIdx,
                                         ObsSrvIdx_SE = mlt_rg_sable_data$ObsSrvIdx_SE,
                                         UseSrvIdx =  mlt_rg_sable_data$UseSrvIdx,
                                         ObsSrvAgeComps = mlt_rg_sable_data$ObsSrvAgeComps,
                                         ISS_SrvAgeComps = mlt_rg_sable_data$ISS_SrvAgeComps,
                                         UseSrvAgeComps = mlt_rg_sable_data$UseSrvAgeComps,
                                         ObsSrvLenComps = mlt_rg_sable_data$ObsSrvLenComps,
                                         UseSrvLenComps = mlt_rg_sable_data$UseSrvLenComps,
                                         ISS_SrvLenComps = mlt_rg_sable_data$ISS_SrvLenComps,

                                         # Model options
                                         srv_idx_type = c("abd", "abd"),
                                         # abundance and biomass for survey fleet 1 and 2
                                         SrvAgeComps_LikeType =
                                           c("Multinomial", "Multinomial"),
                                         # survey age composition likelihood for survey fleet
                                         # 1, and 2
                                         SrvLenComps_LikeType =
                                           c("none", "none"),
                                         #  no length compositions used for survey
                                         SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                                                              "spltRjntS_Year_1-terminal_Fleet_2"),
                                         # survey age comp type
                                         SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                                                              "none_Year_1-terminal_Fleet_2")
                                         # survey length comp type
)

# Fishery Selectivity and Catchability

# defining priors
sex_par <- expand.grid(sex = 1:2, par = 1:2)
fleet_blocks <- data.frame(
  fleet = c(1, 2),
  block = 1
)

# merge together (note that unlike the operational assessment, selectivity
# blocks are reduced from 3 to 2)
fish_selex_structure <- merge(fleet_blocks, sex_par)

# Merge to get all valid combinations
fish_selex_structure <- merge(fleet_blocks, sex_par) %>%
  dplyr::filter(!(fleet == 1 & block == 1 & sex == 2 & par == 2)) %>%              # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING
  dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 1))                  # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
fish_selex_prior <- cbind(
  region = 1,
  fish_selex_structure,
  mu = 2,                                                                      # All selex means = 1 (means should be defined in normal space)
  sd = 3                                                                       # All selex sd = 5
)

fish_selex_prior_tf <- fish_selex_prior %>%                                    # set tighter selex prior for TF
  dplyr::filter((fleet == 2 & par == 1)) %>%
  dplyr::mutate(mu = 2, sd = 1) %>%
  dplyr::full_join(fish_selex_prior %>%  dplyr::filter(!(fleet == 2 & par == 1 )))

fish_selex_prior_tf <- fish_selex_prior_tf %>%                                    # set tighter selex prior for TF
  dplyr::filter((fleet == 2 & par == 2)) %>%
  dplyr::mutate(mu = 5, sd = 2) %>%
  dplyr::full_join(fish_selex_prior_tf %>%  dplyr::filter(!(fleet == 2 & par == 2)))

input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      # fishery selectivity, whether continuous time-varying

                                      # fishery selectivity blocks
                                      fish_sel_blocks =
                                        c("none_Fleet_1",
                                          "none_Fleet_2"),
                                      # no blocks for trawl fishery

                                      # fishery selectivity form
                                      fish_sel_model =
                                        c("logist1_Fleet_1", "gamma_Fleet_2"),

                                      # fishery catchability blocks
                                      fish_q_blocks =
                                        c("none_Fleet_1", "none_Fleet_2"),
                                      # no blocks since q is not estimated

                                      # sharing fishery selex parameters
                                      fish_fixed_sel_pars =
                                        c("est_shared_r", "est_shared_r"),

                                      # whether to estimate all fixed effects
                                      # for fishery catchability
                                      fish_q_spec =
                                        c("fix", "fix"),
                                      Use_fish_selex_prior = 1,
                                      fish_selex_prior = fish_selex_prior
)


# setup survey selectivity
# Define sex and parameter combinations
sex_par <- expand.grid(sex = 1:2, par = 1:2)

# Define valid fleet-block combinations (only estimating domestic and jp LLS)
fleet_blocks <- data.frame(
  fleet = c(1, 2),
  block = c(1, 1)
)

# Merge to get all valid combinations
srv_selex_structure <- merge(fleet_blocks, sex_par)

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
srv_selex_prior <- cbind(
  region = 1,
  srv_selex_structure,
  mu = 1,
  sd = 5
) %>%
  filter(!(fleet == 2 & par == 2 & sex == 2)) %>%
  mutate(mu = ifelse(fleet == 2, 2, mu),
         sd = ifelse(fleet == 2, 3, sd))

input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"
                                       ),

                                     # survey selectivity blocks
                                     srv_sel_blocks =                          # survey selectivity time blocks if not TV specified above for a given fleet
                                       c("none_Fleet_1",
                                         "none_Fleet_2"                        # No blocks for JPN LLS
                                       ),

                                     # survey selectivity form
                                     srv_sel_model =
                                       c("logist1_Fleet_1",
                                         "logist1_Fleet_2"
                                       ),

                                     # survey catchability blocks
                                     srv_q_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"
                                       ),

                                     # whether to estiamte all fixed effects
                                     # for survey selectivity and later
                                     # modify to fix/share parameters
                                     srv_fixed_sel_pars_spec =
                                       c("est_shared_r",
                                         "est_shared_r"
                                       ),

                                     # whether to estiamte all
                                     # fixed effects for survey catchability
                                     # spatially-invariant q
                                     srv_q_spec =
                                       c("est_shared_r",
                                         "est_shared_r"
                                       ),
                                     Use_srv_selex_prior = 1,
                                     srv_selex_prior = srv_selex_prior
)

# set up model weighting stuff
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  Wt_Catch = 1,
                                  Wt_FishIdx = 1,
                                  Wt_SrvIdx = 1,
                                  Wt_Rec = 1,
                                  Wt_F = 1,
                                  Wt_Tagging = 0.5,
                                  # Composition model weighting
                                  Wt_FishAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_FishLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_SrvAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets)),
                                  Wt_SrvLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets))
)

# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map


# Additional Model Specifications -----------------------------------------

# Survey Ages (~100 total across all regions)
data$ISS_SrvAgeComps[] <- 20

# Fishery Ages
data$ISS_FishAgeComps[1,,,,] <- 25  # BS
data$ISS_FishAgeComps[2,,,,] <- 20  # AI
data$ISS_FishAgeComps[3,,,,] <- 14  # WGOA
data$ISS_FishAgeComps[4,,,,] <- 18  # CGOA
data$ISS_FishAgeComps[5,,,,] <- 18  # EGOA

# Fishery Lengths - Fixed Gear
data$ISS_FishLenComps[1,,,,] <- 12  # BS
data$ISS_FishLenComps[2,,,,] <- 12  # AI
data$ISS_FishLenComps[3,,,,] <-  7  # WGOA
data$ISS_FishLenComps[4,,,,] <-  7  # CGOA
data$ISS_FishLenComps[5,,,,] <-  7  # EGOA

# Fishery Lengths - Trawl Gear
data$ISS_FishLenComps[1,,,,] <- 24  # BS
data$ISS_FishLenComps[2,,,,] <- 8  # AI
data$ISS_FishLenComps[3,,,,] <-  4  # WGOA
data$ISS_FishLenComps[4,,,,] <-  4  # CGOA
data$ISS_FishLenComps[5,,,,] <-  4  # EGOA

# Map off early delta for fishery
map_fish_fixed <- array(mapping$fish_fixed_sel_pars, dim = dim(parameters$fish_fixed_sel_pars))
map_fish_fixed[,2,1,2,1]  <- map_fish_fixed[,2,1,1,1] # share deltas

# Map off bmax for trawl females
map_fish_fixed[,1,1,2,2]  <- map_fish_fixed[,1,1,1,2] # share deltas
mapping$fish_fixed_sel_pars <- factor(map_fish_fixed)

# Map off delta for JP LLS
map_srv_fixed <- array(mapping$srv_fixed_sel_pars, dim = dim(parameters$srv_fixed_sel_pars))
map_srv_fixed[,2,1,2,2]  <- map_srv_fixed[,2,1,1,2] # share deltas
mapping$srv_fixed_sel_pars <- factor(map_srv_fixed)

# Some starting values to help out the model
parameters$srv_fixed_sel_pars[] <- log(5)
parameters$fish_fixed_sel_pars[,,,,1] <- log(5) # fixed gear
parameters$fish_fixed_sel_pars[,,,,2] <- log(8) # trawl gear

# make AD model function
# Fit model
st <- Sys.time()
sabie_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = FALSE
)
en <- Sys.time()
print(en - st)

sd_rep <- sdreport(sabie_rtmb_model)
rep <- sabie_rtmb_model$rep

saveRDS(data, here("dev", "dev_output", "5_Region_Model_Sablefish", "data.RDS"))
saveRDS(sd_rep, here("dev", "dev_output", "5_Region_Model_Sablefish", "sd_rep.RDS"))
saveRDS(rep, here("dev", "dev_output", "5_Region_Model_Sablefish", "rep.RDS"))

# write out report file
mlt_rg_sable_rep <- rep
usethis::use_data(mlt_rg_sable_rep, internal = FALSE, overwrite = TRUE, compress = 'xz')

# Vignette Plotshow
rec_series <- reshape2::melt((sabie_rtmb_model$rep$Rec))
rec_series$Par <- "Recruitment"

# Get SSB time-series
ssb_series <- reshape2::melt((sabie_rtmb_model$rep$SSB))
ssb_series$Par <- "Spawning Biomass"

ts_df <- rbind(ssb_series,rec_series) # bind together

# Do some data munging here
ts_df <- ts_df %>% dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS',
    Region == 2 ~ 'AI',
    Region == 3 ~ 'WGOA',
    Region == 4 ~ 'CGOA',
    Region == 5 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
  Year = Year + 1959)

png(here("vignettes", "figures", "g_ts_comparison.png"), width = 1000, height = 800)
ggplot(ts_df, aes(x = Year, y = value, color = Region)) +
  geom_line(size = 1.3) +
  facet_grid(Par~Region, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  labs(y = "Value")  +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none')
dev.off()

get_idx_fits_plot(list(data), list(rep), 1)
get_selex_plot(list(rep), 1)

