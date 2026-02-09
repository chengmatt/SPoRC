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

data("three_rg_sable_data")

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(three_rg_sable_data$years), # vector of years (1 - 62)
                            ages = 1:30, # vector of ages (1 - 30)
                            lens = three_rg_sable_data$lens, # number of lengths (41 - 99)
                            n_regions = three_rg_sable_data$n_regions, # number of regions (5)
                            n_sexes = three_rg_sable_data$n_sexes, # number of sexes (2)
                            n_fish_fleets = three_rg_sable_data$n_fish_fleets, # number of fishery fleet (2)
                            n_srv_fleets = three_rg_sable_data$n_srv_fleets, # number of survey fleets (2)
                            verbose = TRUE
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            do_rec_bias_ramp = 0, # not using bias ramp
                            sigmaR_switch = 16, # switch to using late sigma in year 16
                            dont_est_recdev_last = 1, # don't estimate last rec dev

                            # Model options
                            rec_model = "mean_rec", # recruitment model
                            h_spec = 'fix',
                            sigmaR_spec = "fix", # fixing
                            InitDevs_spec = "est_shared_r", # initial deviations are shared across regions,
                            # but recruitment deviations are region specific
                            ln_sigmaR = log(c(0.4, 1.2)), # values to fix sigmaR at, or starting values
                            ln_global_R0 = log(20),
                            # starting value for global R0
                            R0_prop = array(c(0.2, 0.2),
                                            dim = c(input_list$data$n_regions - 1))
                            # starting value for R0 proportions in multinomial logit space
)

# Setup biological stuff (using defaults for other stuff)
input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = three_rg_sable_data$WAA, # weight at age
                                    MatAA = three_rg_sable_data$MatAA, # maturity at age
                                    AgeingError = three_rg_sable_data$AgeingError, # ageing error matrix
                                    fit_lengths = 1, # fitting lengths
                                    SizeAgeTrans = three_rg_sable_data$SizeAgeTrans, # size age transition matrix
                                    M_spec = "fix", # fix natural mortality
                                    Fixed_natmort = array(0.104884, dim = c(three_rg_sable_data$n_regions,
                                                                            length(three_rg_sable_data$years),
                                                                            length(three_rg_sable_data$ages),
                                                                            three_rg_sable_data$n_sexes))
                                    # value to fix natural mortality at
)

# setting up movement parameterization
Movement_prior <- expand.grid(
  region_from = 1:3, # regions
  year = 1, # penalize first year since no blocks
  age = c(6,7,16), # age blocks
  sex = 1, # sex
  alpha = I(list(rep(2.5, 3))) # prior alpha to each row
)

input_list <- Setup_Mod_Movement(input_list = input_list,
                                 # Model options
                                 Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)), # estimating movement in 3 age blocks
                                 # (ages 1-6, ages 7-15, ages 16-30)
                                 Movement_yearblk_spec = "constant", # time-invariant for movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 do_recruits_move = 0, # recruits do not move
                                 use_fixed_movement = 0, # estimating movement
                                 Use_Movement_Prior = 1, # priors used for movement
                                 Movement_prior = Movement_prior, # vague prior to penalize movement away from the extremes
                                 cont_vary_movement = 'none'
)

# setting up tagging parameterization
# setup tagging priors
tag_prior <- data.frame(
  region = 1,
  block = c(1,2),
  mu = NA, # no mean, since symmetric beta
  sd = 5, # sd = 5
  type = 0 # symmetric beta
)

input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 1, # using tagging data
                                max_tag_liberty = 15, # maximum number of years to track a cohort

                                # Data Inputs
                                tag_release_indicator = three_rg_sable_data$tag_release_indicator,
                                # tag release indicator (first col = tag region, second col = tag year),
                                # total number of rows = number of tagged cohorts
                                Tagged_Fish = three_rg_sable_data$Tagged_Fish, # Released fish
                                # dimensioned by total number of tagged cohorts, (implicitly
                                # tracks the release year and region), age, and sex
                                Obs_Tag_Recap = three_rg_sable_data$Obs_Tag_Recap,
                                # dimensioned by max tag liberty, tagged cohorts, regions,
                                # ages, and sexes

                                # Model options
                                Tag_LikeType = "NegBin", # Negative Binomial
                                mixing_period = 2, # Don't fit tagging until release year + 1
                                t_tagging = 0.5, # tagging happens midway through the year,
                                # movement does not occur within that year
                                tag_selex = "SexSp_AllFleet", # tagging recapture selectivity
                                # is a weighted average of fishery selectivity of two fleets
                                tag_natmort = "AgeSp_SexSp", # tagging natural mortality is
                                # age and sex-specific
                                Use_TagRep_Prior = 1, # tag reporting rate priors are used
                                TagRep_Prior = tag_prior,
                                move_age_tag_pool = list(c(1:6), c(7:15), c(16:30)), # whether or
                                # not to pool tagging data when fitting (for computational cost)
                                move_sex_tag_pool = list(c(1:2)), # whether or not to pool
                                # sex-specific data when fitting
                                Init_Tag_Mort_spec = "fix", # fixing initial tag mortality
                                Tag_Shed_spec = "fix", # fixing chronic shedding
                                TagRep_spec = "est_shared_r", # tag reporting rates are not region specific
                                # Time blocks for tag reporting rates
                                Tag_Reporting_blocks = c(
                                  paste("Block_1_Year_1-35_Region_", c(1:input_list$data$n_regions), sep = ''),
                                  paste("Block_2_Year_36-terminal_Region_", c(1:input_list$data$n_regions), sep = '')
                                ),

                                # Specify starting values or fixing values
                                ln_Init_Tag_Mort = log(0.1), # fixing initial tag mortality
                                ln_Tag_Shed = log(0.02),  # fixing tag shedding
                                ln_tag_theta = log(0.5), # starting value for tagging overdispersion
                                Tag_Reporting_Pars = array(log(0.2 / (1-0.2)), # starting values for tag reporting pars
                                                           dim = c(input_list$data$n_regions, 2))
)


# setting up catch data
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = three_rg_sable_data$ObsCatch,
                                    UseCatch = three_rg_sable_data$UseCatch,
                                    # Model options
                                    Use_F_pen = 1,
                                    # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix',
                                    ln_sigmaC =
                                      array(log(0.05), dim = c(input_list$data$n_regions,
                                                               length(input_list$data$years),
                                                               input_list$data$n_fish_fleets)),
                                    # fixing catch sd at small value
                                    ln_F_mean = array(-2, dim = c(input_list$data$n_regions,
                                                                  input_list$data$n_fish_fleets))
                                    # some starting values for fishing mortality
)

# Fishery Indices and Compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = three_rg_sable_data$ObsFishIdx,
                                          ObsFishIdx_SE = three_rg_sable_data$ObsFishIdx_SE,
                                          UseFishIdx =  three_rg_sable_data$UseFishIdx,
                                          ObsFishAgeComps = three_rg_sable_data$ObsFishAgeComps,
                                          UseFishAgeComps = three_rg_sable_data$UseFishAgeComps,
                                          ISS_FishAgeComps = three_rg_sable_data$ISS_FishAgeComps,
                                          ObsFishLenComps = three_rg_sable_data$ObsFishLenComps,
                                          UseFishLenComps = three_rg_sable_data$UseFishLenComps,
                                          ISS_FishLenComps = three_rg_sable_data$ISS_FishLenComps,

                                          # Model options
                                          fish_idx_type = c("none", "none"),
                                          # fishery indices not used
                                          FishAgeComps_LikeType = c("Multinomial", "none"),
                                          # age comp likelihoods for fishery fleet 1 and 2
                                          FishLenComps_LikeType = c("Multinomial", "Multinomial"),
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
                                         ObsSrvIdx = three_rg_sable_data$ObsSrvIdx,
                                         ObsSrvIdx_SE = three_rg_sable_data$ObsSrvIdx_SE,
                                         UseSrvIdx =  three_rg_sable_data$UseSrvIdx,
                                         ObsSrvAgeComps = three_rg_sable_data$ObsSrvAgeComps,
                                         ISS_SrvAgeComps = three_rg_sable_data$ISS_SrvAgeComps,
                                         UseSrvAgeComps = three_rg_sable_data$UseSrvAgeComps,
                                         ObsSrvLenComps = three_rg_sable_data$ObsSrvLenComps,
                                         UseSrvLenComps = three_rg_sable_data$UseSrvLenComps,
                                         ISS_SrvLenComps = three_rg_sable_data$ISS_SrvLenComps,

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
)

# Fishery Selectivity and Catchability
input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      # fishery selectivity, whether continuous time-varying

                                      # fishery selectivity blocks
                                      fish_sel_blocks =
                                        c("Block_1_Year_1-56_Fleet_1",
                                          # block 1, fishery ll selex
                                          "Block_2_Year_57-terminal_Fleet_1",
                                          # block 3 fishery ll selex
                                          "none_Fleet_2"),
                                      # no blocks for trawl fishery

                                      # fishery selectivity form
                                      fish_sel_model =
                                        c("logist1_Fleet_1",
                                          "gamma_Fleet_2"),

                                      # fishery catchability blocks
                                      fish_q_blocks =
                                        c("none_Fleet_1",
                                          "none_Fleet_2"),
                                      # no blocks since q is not estimated

                                      # whether to estimate all fixed effects
                                      # for fishery selectivity and later modify
                                      # to fix and share parameters
                                      fish_fixed_sel_pars_spec =
                                        c("est_all", "est_all"),

                                      # whether to estimate all fixed effects
                                      # for fishery catchability
                                      fish_q_spec =
                                        c("fix", "fix")
                                      # fix fishery q since not used
)

# Custom parameter sharing for fishery selectivity
map_ln_fish_fixed_sel_pars <- input_list$par$ln_fish_fixed_sel_pars # mapping fishery selectivity

# Fixed gear fleet, unique parameters for each sex (time block 1)
map_ln_fish_fixed_sel_pars[,1,1,1,1] <- 1 # a50, female, time block 1, fixed gear
map_ln_fish_fixed_sel_pars[,2,1,1,1] <- 2 # delta, female, time block 1, fixed gear (shared with time block 2 and sex)
map_ln_fish_fixed_sel_pars[,1,1,2,1] <- 3 # a50, male, time block 1, fixed gear
map_ln_fish_fixed_sel_pars[,2,1,2,1] <- 4 # delta, male, time block 1, fixed gear (shared with time block 2 and sex)

# time block 2, fixed gear fishery
map_ln_fish_fixed_sel_pars[,1,2,1,1] <- 5 # a50, female, time block 2, fixed gear
map_ln_fish_fixed_sel_pars[,2,2,1,1] <- 2 # delta, female, time block 2, fixed gear (shared with time block 1 and sex)
map_ln_fish_fixed_sel_pars[,1,2,2,1] <- 6 # a50, male, time block 2, fixed gear
map_ln_fish_fixed_sel_pars[,2,2,2,1] <- 4 # delta, male, time block 2, fixed gear (shared with time block 1 and sex)

# time block 1 and 2, trawl gear fishery
map_ln_fish_fixed_sel_pars[,1,1,1,2] <- 7 # amax, female, time block 1, trawl gear
map_ln_fish_fixed_sel_pars[,2,1,1,2] <- 8 # delta, female, time block 1, trawl gear (shared by sex)
map_ln_fish_fixed_sel_pars[,1,1,2,2] <- 9 # amax, male, time block 1, trawl gear
map_ln_fish_fixed_sel_pars[,2,1,2,2] <- 8 # delta, male, time block 1, trawl gear (shared by sex)
map_ln_fish_fixed_sel_pars[,,2,,2] <- NA # no parameters estimated for time block 2 trawl gear

input_list$map$ln_fish_fixed_sel_pars <- factor(map_ln_fish_fixed_sel_pars) # input into map list
input_list$par$ln_fish_fixed_sel_pars[] <- log(0.1) # some more inforamtive starting values

# Survey Selectivity and Catchability
input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"),

                                     # survey selectivity blocks
                                     srv_sel_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"
                                       ), # no blocks for jp and domestic survey

                                     # survey selectivity form
                                     srv_sel_model =
                                       c("logist1_Fleet_1",
                                         "logist1_Fleet_2"),

                                     # survey catchability blocks
                                     srv_q_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"),

                                     # whether to estiamte all fixed effects
                                     # for survey selectivity and later
                                     # modify to fix/share parameters
                                     srv_fixed_sel_pars_spec =
                                       c("est_all",
                                         "est_all"),

                                     # whether to estiamte all
                                     # fixed effects for survey catchability
                                     # spatially-invariant q
                                     srv_q_spec =
                                       c("est_shared_r", "est_shared_r"),

                                     # Starting values for survey catchability
                                     ln_srv_q = array(8.75,
                                                      dim = c(input_list$data$n_regions, 1,
                                                              input_list$data$n_srv_fleets))
)

# Custom mapping survey selectivity stuff
map_ln_srv_fixed_sel_pars <- input_list$par$ln_srv_fixed_sel_pars # set up mapping factor stuff

# Coop survey (japanese)
map_ln_srv_fixed_sel_pars[,1,1,1,1] <- 1 # a50, coop survey, time block 1, female
map_ln_srv_fixed_sel_pars[,2,1,1,1] <- 2 # delta, coop survey, time block 1, female (sharing with domestic survey)
map_ln_srv_fixed_sel_pars[,1,1,2,1] <- 3 # a50, coop survey, time block 1, male
map_ln_srv_fixed_sel_pars[,2,1,2,1] <- 2 # delta, coop survey, time block 1, male (sharing with domestic survey)

# domestic survey
map_ln_srv_fixed_sel_pars[,1,1,1,2] <- 5 # a50, domestic survey, time block 1, female
map_ln_srv_fixed_sel_pars[,2,1,1,2] <- 2 # delta, domestic survey, time block 1, female (sharing with coop survey)
map_ln_srv_fixed_sel_pars[,1,1,2,2] <- 6 # a50, domestic survey, time block 1, male
map_ln_srv_fixed_sel_pars[,2,1,2,2] <- 2 # delta, domestic survey, time block 1, male (sharing with coop survey)

input_list$map$ln_srv_fixed_sel_pars <- factor(map_ln_srv_fixed_sel_pars)  # input into map list
input_list$par$ln_srv_fixed_sel_pars[] <- log(0.1) # some more informative starting values


# set up model weighting stuff
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  Wt_Catch = 1,
                                  Wt_FishIdx = 1,
                                  Wt_SrvIdx = 1,
                                  Wt_Rec = 1,
                                  Wt_F = 1,
                                  Wt_Tagging = 1,
                                  # Composition model weighting
                                  Wt_FishAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_FishLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_SrvAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets)),
                                  Wt_SrvLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets))
)

# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Pre-specify ISS
data$ISS_FishAgeComps[] <- 50
data$ISS_FishLenComps[] <- 40
data$ISS_SrvAgeComps[] <- 60

# make AD model function
obj <- RTMB::MakeADFun(SPoRC:::cmb(SPoRC:::SPoRC_rtmb, data), parameters = parameters,
                       map = mapping, silent = FALSE)

# Now, optimize the function
optim <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

# newton steps
try_improve <- tryCatch(expr =
                          for(i in 1:5) {
                            g = as.numeric(obj$gr(optim$par))
                            h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                            optim$par = optim$par - solve(h,g)
                            optim$objective = obj$fn(optim$par)
                          }
                        , error = function(e){e}, warning = function(w){w})

obj$sd_rep <- sdreport(obj)
obj$rep <- obj$report(obj$env$last.par.best)

sd_rep <- obj$sd_rep
rep <- obj$rep

saveRDS(data, here("dev", "dev_output", "3_Region_Model_Sablefish", "data.RDS"))
saveRDS(input_list, here("dev", "dev_output", "3_Region_Model_Sablefish", "input_list.RDS"))
saveRDS(sd_rep, here("dev", "dev_output", "3_Region_Model_Sablefish", "sd_rep.RDS"))
saveRDS(rep, here("dev", "dev_output", "3_Region_Model_Sablefish", "rep.RDS"))

# Quick plots for inspection
reshape2::melt(rep$Rec) %>%
  ggplot(aes(x = Var2 + 1959, y = value)) +
  geom_line() +
  facet_wrap(~Var1) +
  labs(x = 'Year', y = 'Recruitment')

reshape2::melt(rep$SSB) %>%
  ggplot(aes(x = Var2, y = value)) +
  geom_line() +
  ylim(0, NA) +
  facet_wrap(~Var1) +
  labs(x = 'Year', y = 'SSB')

get_idx_fits(data, rep, 1960:2021) %>%
  filter(obs != 0) %>%
  ggplot() +
  geom_pointrange(aes(x = Year, y = obs, ymin = lci, ymax = uci)) +
  geom_line(aes(x = Year, y = value)) +
  facet_grid(Region~Fleet, scales = 'free_x')

# Fishing Mortality
f_plot_df <- reshape2::melt(rep$Fmort) %>%
  dplyr::rename(Region = Var1, Year = Var2, Type = Var3) %>%
  dplyr::mutate(Region = paste("Region", Region),
                Type = paste("Fleet", Type, "F"),
                lwr = NA,
                upr = NA,
                se = NA)

ggplot2::ggplot(f_plot_df,
                ggplot2::aes(x = Year, y = value)) +
  ggplot2::geom_line(lwd = 0.9) +
  ggplot2::facet_grid(Type~Region, scales = 'free') +
  ggplot2::labs(x = 'Year', y = 'Value') +
  ggplot2::coord_cartesian(ylim = c(0,NA)) +
  theme_sablefish()

