# Purpose: To bridge to the spatial assessment in Cheng and Marsh et al. 2025  using SPoRC
# Creator: Matthew LH. Cheng
# Date Created: 2/5/25

# install.packages("devtools")
# # install.packages("TMB")
# install.packages("RTMB")
# TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
# remotes::install_github("fishfollower/compResidual/compResidual")
# devtools::install_github("chengmatt/SPoRC", dependencies = c("Depends", "Imports"))

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)
library(pbapply)
library(doParallel)


# load in libraries
library(igraph)
library(splines2)
library(mgcv)

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
                            n_seas = three_rg_sable_data$n_seas,
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
                                Tag_LikeType = "Multinomial_Release", # Negative Binomial
                                mixing_period = 2, # Don't fit tagging until release year + 1
                                t_tagging = 0.5, # tagging happens midway through the year,
                                # movement does not occur within that year
                                tag_selex = "SexSp_AllFleet", # tagging recapture selectivity
                                # is a weighted average of fishery selectivity of two fleets
                                tag_natmort = "AgeSp_SexSp", # tagging natural mortality is
                                # age and sex-specific
                                Use_TagRep_Prior = 1, # tag reporting rate priors are used
                                TagRep_Prior = tag_prior,
                                move_age_tag_pool = as.list(1:30), # whether or
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
                                                               three_rg_sable_data$n_seas,
                                                               input_list$data$n_fish_fleets))
                                    # fixing catch sd at small value
)

# ln_F_mean cannot be passed through ... because R partially matches the name to
# the ln_F_mean_spec formal, so the starting value is assigned post-hoc.
input_list$par$ln_F_mean[] <- -2

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

                                      # whether to estimate all fixed effects for fishery
                                      # selectivity, modified below to fix or share
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
input_list$par$ln_fish_fixed_sel_pars[] <- log(6) # some more inforamtive starting values

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

                                     # whether to estimate all fixed effects for survey
                                     # selectivity, modified below to fix or share
                                     srv_fixed_sel_pars_spec =
                                       c("est_all",
                                         "est_all"),

                                     # whether to estimate all fixed effects for survey
                                     # catchability, spatially-invariant q
                                     srv_q_spec =
                                       c("est_shared_r", "est_shared_r"),

                                     # Starting values for survey catchability
                                     ln_srv_q = array(8.75,
                                                      dim = c(input_list$data$n_regions, 1,
                                                              input_list$data$n_srv_fleets)),
                                     t_srv = array(0.5, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_srv_fleets))
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
input_list$par$ln_srv_fixed_sel_pars[] <- log(3) # some more informative starting values


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
                                                     length(input_list$data$years), input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_FishLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),  input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_SrvAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),  input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets)),
                                  Wt_SrvLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),  input_list$data$n_seas,
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets))
)


# Setup Movement EMs ------------------------------------------------------


### CTMC Stuff --------------------------------------------------------------
# CTMC stuff
# make adjacency matrix
adjacency <- as_adjacency_matrix(
  make_graph(
    ~ 1 - 2,
    2 - 3,
    1 - 3
  )
)

# make ctmc data
ctmc_data <- expand.grid(
  regions = 1:input_list$data$n_regions,
  years = input_list$data$years,
  seas = 1:input_list$data$n_seas,
  ages = input_list$data$ages,
  sexes = 1:input_list$data$n_sexes
)

diffusion_formula = ~1
preference_formula = ~ 0 + factor(regions):splines2::bSpline(ages, df = 4, intercept = TRUE)


### Model Scenarios ---------------------------------------------------------
model_df <- data.frame(
  move_type = c(0, 0, 1),
  ctmc_diffusion_bounds = c(NA, NA, 0),
  stringsAsFactors = FALSE
)

# Add in movement age block specification
model_df$age_blk_spec <- list(
  "constant",
  list(c(1:6), c(7:15), c(16:30)),
  "constant"
)

# Setup Parrallel ---------------------------------------------------------

# Initialize parallel cluster
cl <- makeCluster(nrow(model_df))
registerDoParallel(cl)

# Load packages and export objects to cluster
clusterEvalQ(cl, {
  library(SPoRC)
  library(RTMB)
  library(here)
})

clusterExport(cl,
              c("input_list", "ctmc_data", "adjacency", "diffusion_formula",
                "preference_formula", "model_df"),
              envir = environment())

clusterEvalQ(cl, {
  devtools::load_all(here("R"))
})


# Run Parallel -----------------------------------------------------------

# Run parallel loop across movement model scenarios
model <- pblapply(1:nrow(model_df), function(i) {

  # Setup movement model
  input_list <- SPoRC:::Setup_Mod_Movement(
    input_list = input_list,
    move_type = model_df$move_type[i],
    do_recruits_move = 0,
    ctmc_move_dat = ctmc_data,
    Movement_ageblk_spec = model_df$age_blk_spec[[i]],
    adjacency_mat = adjacency,
    area_r = c(0.745, 0.145, 0.11),
    diffusion_formula = diffusion_formula,
    preference_formula = preference_formula,
    ctmc_diffusion_bounds = model_df$ctmc_diffusion_bounds[i],
    Movement_prior = 0,
  )

  # Extract updated lists
  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map
  input_list$data$Use_Movement_Prior <- 0
  parameters$log_move_diffusion_pars <- log(0.15) # mess a bit with starting values
  parameters$ln_fish_fixed_sel_pars[,,,,2] <- log(10) # some more inforamtive starting values
  parameters$ln_srv_fixed_sel_pars[,,,,] <- log(0.5) # some more inforamtive starting values

  # fit model
  obj <- fit_model(
    data = data,
    parameters = parameters,
    mapping = mapping,
    random = NULL,
    silent = F,
    newton_loops = 3
  )

  obj$gr()

  # Generate reports
  obj$sd_rep <- sdreport(obj)

  plot(obj$rep$srv_sel[1,1,,1,1])

  return(obj)
}, cl = cl)

# Stop cluster
stopCluster(cl)

saveRDS(model, here("dev", "dev_output", "CTMC", "movement_model_param.RDS"))

# Results ---------------------------------------------------------------
model <- readRDS((here("dev", "dev_output", "CTMC", "movement_model_param.RDS")))
model_names <- c('unstrc_mark_const', 'unstrc_mark_age', 'ctmc_age', 'ctmc_age_bnd')
model_colors <- c("#C73E3A", "#2A9D8F", "#E9B44C", "#9B59B6")

### SSB ---------------------------------------------------------------------
ts_res <- bind_rows(
  lapply(seq_along(model), function(i) {
    bind_rows(
      # SSB data
      reshape2::melt(model[[i]]$rep$SSB) %>%
        mutate(model = model_names[i], type = 'Spawning Stock Biomass (kt)') %>%
        bind_cols(se = model[[i]]$sd_rep$sd[names(model[[i]]$sd_rep$value) == 'log(SSB)']),
      # Recruitment data
      reshape2::melt(model[[i]]$rep$Rec) %>%
        mutate(model = model_names[i], type = 'Age-2 Recruitment (millions)') %>%
        bind_cols(se = model[[i]]$sd_rep$sd[names(model[[i]]$sd_rep$value) == 'log(Rec)'])
    )
  })
) %>%
  rename(Region = Var1, Year = Var2) %>%
  mutate(
    Region = case_when(
      Region == 1 ~ "BS+AI+WGOA",
      Region == 2 ~ "CGOA",
      Region == 3 ~ "EGOA"
    )
  )

# remove slipstream diffussion from plot
ts_plot <- ggplot(ts_res %>% filter(model != 'ctmc_age_bnd') %>%
                    mutate(type = factor(type, levels = c("Spawning Stock Biomass (kt)", "Age-2 Recruitment (millions)"))),
                   aes(x = Year + 1959, y = value, color = model, ymin = value * exp(-1.96 * se),
                                ymax = value * exp(1.96 * se), lty = model, fill = model)) +
  geom_ribbon(alpha = 0.25, color = NA, lty = NA) +
  geom_line(lwd = 1) +
  facet_grid(type~Region) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 20) +
  theme(legend.position = 'top') +
  labs(x = 'Year', y = 'Value', color = 'Model', lty = 'Model', fill = 'Model')

ggsave(
  here("dev", "spatial_movement_ctmc", "ts_sab.png"),
  ts_plot,
  height = 10, width = 15
)


### Movement ----------------------------------------------------------------
move_res <- reshape2::melt(model[[1]]$rep$Movement) %>% mutate(model = model_names[1]) %>%
  bind_rows(
    reshape2::melt(model[[2]]$rep$Movement) %>% mutate(model = model_names[2]),
    reshape2::melt(model[[3]]$rep$Movement) %>% mutate(model = model_names[3]),
    reshape2::melt(model[[4]]$rep$Movement) %>% mutate(model = model_names[4])
  ) %>%
  mutate(
    from = case_when(
      from == 1 ~ "from BS+AI+WGOA",
      from == 2 ~ "from CGOA",
      from == 3 ~ "from EGOA",
    ),
    to = case_when(
      to == 1 ~ "to BS+AI+WGOA",
      to == 2 ~ "to CGOA",
      to == 3 ~ "to EGOA",
    ))

movement_plot <- ggplot(move_res %>% filter(ages != 1, model != 'ctmc_age_bnd'), aes(x = ages, y = value, color = model, linetype = model)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = model_colors) +
  ggh4x::facet_grid2(to~from, scales = 'free_y', independent = 'y') +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 20) +
  theme(legend.position = 'top') +
  labs(x = 'Ages', y = 'Fraction of Individuals Moving', color = 'Model', lty = 'Model',)

ggsave(
  here("dev", "spatial_movement_ctmc", "move_sab.png"),
  movement_plot,
  height = 10, width = 13
)

### AIC ---------------------------------------------------------------------

marg_AIC(model[[1]]$optim)
marg_AIC(model[[2]]$optim)
marg_AIC(model[[3]]$optim)
# marg_AIC(model[[4]]$optim)


# SPoRC::get_nLL_plot(
#   list(data, data, data, data),
#   rep = list(model[[1]]$rep, model[[2]]$rep,
#              model[[3]]$rep, model[[4]]$rep),
#   model_names
# )

