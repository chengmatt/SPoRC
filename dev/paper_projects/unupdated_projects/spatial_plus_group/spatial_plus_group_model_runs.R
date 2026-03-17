# Purpose: To demonstrate implications of plus group initialization for spatial assessments
# Creator: Matthew LH. Cheng
# Date Created: 2/5/25

# Set up ------------------------------------------------------------------
unloadNamespace("SPoRC")
library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)
library(bench)
library(foreach)
library(doParallel)
library(patchwork)
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

# setting up movement parameterization
# vague prior to penalize movement away from the extremes
Movement_prior <- expand.grid(
  region_from = 1:3, # regions
  year = 1, # years
  seas = 1, # seasons
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
                                                               input_list$data$n_seas,
                                                               input_list$data$n_fish_fleets)),
                                    # fixing catch sd at small value
                                    ln_F_mean = array(-2, dim = c(input_list$data$n_regions, input_list$data$n_seas,
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


# Fit Model ---------------------------------------------------------------

# Fit model
init_age_strc <- c(0,1,2,3) # ext_ages, scalar w/ no movement, matrix, scalar with movement

# Initialize parallel cluster
cl <- makeCluster(length(init_age_strc))
registerDoParallel(cl)

# load in sporc and other objects to cluster
clusterEvalQ(cl,
 {library(SPoRC)
  library(here)})

clusterExport(cl, c("init_age_strc", "data", "parameters", "mapping"), envir = environment())

# run parrallel loop
model <- parLapply(cl, 1:length(init_age_strc), function(i) {
  devtools::load_all(here("R"))

  # iterate age structure init options
  data$init_age_strc <- init_age_strc[i]

  obj <- SPoRC::fit_model(
    data,
    parameters,
    mapping,
    random = NULL,
    newton_loops = 3
  )

  # get sd report
  obj$sd_rep <- sdreport(obj)

  obj
})


# save model results
saveRDS(model, here("dev", "paper_projects", "spatial_plus_group", "model_results.RDS"))

# Plots -------------------------------------------------------------------

model <- readRDS(here("dev", "paper_projects", "spatial_plus_group", "model_results.RDS"))
# loop to store results
rep_list <- list()
sd_list <- list()

for(i in 1:length(init_age_strc)) {
  rep_list[[i]] <- model[[i]]$rep
  sd_list[[i]] <- model[[i]]$sd_rep
}

plots <- get_ts_plot(
  rep_list,
  sd_list,
  model_names = c('ext_ages', 'no_move_all', 'matrix', "no_move_plus"),
  do_ci = F
)

# get initial numbers
derive_sbl_init_naa <- function(init_age_strc, rep, sd) {

  # dims
  n_regions <- 3
  n_ages <- 30
  n_sexes <- 2
  n_fish_fleets <- 2
  n_seas <- 1

  out <- Get_Init_NAA(
    init_age_strc = init_age_strc, # initial age structure
    init_iter = n_ages * 5, # if init_age_strc == 0, number of iterations to run
    n_regions = n_regions, # regions
    n_sexes = n_sexes, # sexes
    n_ages = n_ages, # ages
    natmort = array(rep$natmort[,1,,], dim = c(n_regions, n_ages, n_sexes)), # natural mortality in first year
    init_F = 0, # initial F applied (0 for unfished)
    fish_sel = array(rep$fish_sel[,1,,,], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)), # fishery selectivity in first year
    R0_r = rep$R0 * rep$Rec_trans_prop, # regional mean or virgin recruitment
    sexratio = array(rep$sexratio[,1,], dim = c(n_regions, n_sexes)), # sex ratio in first year
    Movement = array(rep$Movement[,,1,,,], dim = c(n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
    do_recruits_move = 0, # whether recruits move
    # seasonal stuff
    n_seas = 1,
    seasdur = 1,
    ln_InitDevs = array(sd$par.fixed[names(sd$par.fixed) == 'ln_InitDevs'], dim = c(n_regions,n_ages-1)) # initial deviations
  )
  return(out)
}


init_naa_sbl <- reshape2::melt(
  derive_sbl_init_naa(init_age_strc = 0, rep = rep_list[[1]], sd = sd_list[[1]])
) %>%
  mutate(Model = 'ext_ages') %>%
  bind_rows(

    reshape2::melt(
      derive_sbl_init_naa(init_age_strc = 1, rep = rep_list[[2]], sd = sd_list[[2]])
    ) %>% mutate(Model = 'no_move_all'),

    reshape2::melt(
      derive_sbl_init_naa(init_age_strc = 2, rep = rep_list[[3]], sd = sd_list[[3]])
    ) %>% mutate(Model = 'matrix'),

    reshape2::melt(
      derive_sbl_init_naa(init_age_strc = 3, rep = rep_list[[4]], sd = sd_list[[4]])
    ) %>% mutate(Model = 'no_move_plus')
  ) %>%
  rename(Region = Var1, Age = Var2)

init_naa_sbl <- init_naa_sbl %>%
  mutate(Region = case_when(
    Region == 1 ~ "BS+AI+WGOA",
    Region == 2 ~ "CGOA",
    Region == 3 ~ "EGOA"
  ),
  Region = factor(Region, levels = c("BS+AI+WGOA", "CGOA", "EGOA")),
  Model = factor(Model, levels = c("no_move_all", "no_move_plus", "ext_ages", "matrix")))

## Bar Plot ----------------------------------------------------------------
# Inset theme
inset_theme <- theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

main <- init_naa_sbl %>%
  ggplot(aes(x = Age, y = value, lty = Model, color = Model)) +
  geom_line(lwd = 1.3) +
  facet_wrap(~Region) +
  labs(fill = 'Model', y = 'Equilibrium Abundance') +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none')

# Function to create inset plot for a specific region
create_inset <- function(region_val) {
  init_naa_sbl %>%
    filter(Age == 30, Region == region_val) %>%
    ggplot(aes(x = Age, y = value, fill = Model)) +
    geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
    ggthemes::scale_fill_colorblind() +
    labs(x = NULL, y = NULL) +
    inset_theme
}

inset_region1 <- create_inset("BS+AI+WGOA")
inset_region2 <- create_inset("CGOA")
inset_region3 <- create_inset("EGOA")

# Combine with insets
abd_plot <- main +
  inset_element(inset_region1, left = 0.15, bottom = 0.70, right = 0.3, top = 0.95) +
  inset_element(inset_region2, left = 0.48, bottom = 0.70, right = 0.63, top = 0.95) +
  inset_element(inset_region3, left = 0.83, bottom = 0.70, right = 0.98, top = 0.95)


## Time Series -------------------------------------------------------------

# get ssb estimates
ssb_data <- plots[[4]]$data

# make relative to matrix and do some munging
ssb_data <- ssb_data %>%
  select(Region, Year, value, Model) %>%
  group_by(Region, Year) %>%
  mutate(
    matrix_value = value[Model == "matrix"],
    rel_diff = (value - matrix_value) / matrix_value,
  ) %>%
  ungroup() %>%
  mutate(Region = case_when(
    Region == "Region 1" ~ "BS+AI+WGOA",
    Region == "Region 2" ~ "CGOA",
    Region == "Region 3" ~ "EGOA"
  ),
  Region = factor(Region, levels = c("BS+AI+WGOA", "CGOA", "EGOA")),
  Model = factor(Model, levels = c("no_move_all", "no_move_plus", "ext_ages", "matrix"))) %>%
  group_by(Region, Model) %>%
  mutate(Depletion = value / value[which.min(Year)]) %>%
  group_by(Region, Year) %>%
  mutate(
    dep_matrix_value = Depletion[Model == "matrix"],
    dep_rel_diff = (Depletion - dep_matrix_value) / dep_matrix_value,
  ) %>%
  ungroup()

# spawning stock biomass
ssb_plot <- ggplot(ssb_data, aes(x = Year + 1959, y = value, color = Model, lty = Model)) +
  geom_line(lwd = 1.3) +
  facet_wrap(~Region) +
  coord_cartesian(ylim = c(0,NA)) +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top') +
  labs(x = 'Year', y = 'Spawning Stock Biomass (kt)')

# depletion
dep_plot <- ggplot(ssb_data, aes(x = Year + 1959, y = Depletion, color = Model, lty = Model)) +
  geom_line(lwd = 1.3) +
  facet_wrap(~Region) +
  coord_cartesian(ylim = c(0,NA)) +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none') +
  labs(x = 'Year', y = 'Depletion')

# combined plot
comb <- cowplot::plot_grid(ssb_plot, dep_plot, abd_plot, ncol = 1, align = 'v', axis = 'l', labels = c('A', 'B', 'C'), label_size = 30)

ggsave(
  here("dev", "paper_projects", "spatial_plus_group", "sablefish_plus_group_plot.png"),
  comb,
  width = 13, height = 13
)

# Condensed plot
ggsave(
  here("dev", "paper_projects", "spatial_plus_group", "sablefish_plus_group_condensed_plot.png"),
  cowplot::plot_grid(ssb_plot +
                       theme_bw(base_size = 25) +
                       theme(legend.position = 'top'),
                     # Combine with insets
                      main + theme_bw(base_size = 25) + theme(legend.position = 'none') +
                       inset_element(inset_region1, left = 0.15, bottom = 0.70, right = 0.3, top = 0.95) +
                       inset_element(inset_region2, left = 0.48, bottom = 0.70, right = 0.63, top = 0.95) +
                       inset_element(inset_region3, left = 0.83, bottom = 0.70, right = 0.98, top = 0.95),
                     ncol = 1, align = 'v', axis = 'l'), width = 12, height = 12
)


## Recruitment and Movement -------------------------------------------------------------

# get rec estimates
rec_data <- plots[[3]]$data

# make relative to matrix and do some munging
rec_data <- rec_data %>%
  select(Region, Year, value, Model) %>%
  group_by(Region, Year) %>%
  mutate(
    matrix_value = value[Model == "matrix"],
    rel_diff = (value - matrix_value) / matrix_value,
  ) %>%
  ungroup() %>%
  mutate(Region = case_when(
    Region == "Region 1" ~ "BS+AI+WGOA",
    Region == "Region 2" ~ "CGOA",
    Region == "Region 3" ~ "EGOA"
  ),
  Region = factor(Region, levels = c("BS+AI+WGOA", "CGOA", "EGOA")),
  Model = factor(Model, levels = c("no_move_all", "no_move_plus", "ext_ages", "matrix"))) %>%
  group_by(Region, Model) %>%
  mutate(Depletion = value / value[which.min(Year)]) %>%
  group_by(Region, Year) %>%
  mutate(
    dep_matrix_value = Depletion[Model == "matrix"],
    dep_rel_diff = (Depletion - dep_matrix_value) / dep_matrix_value,
  ) %>%
  ungroup()

# recruitment
rec_plot <- ggplot(rec_data %>% filter(Model == 'matrix'),
                   aes(x = Year + 1959, y = value)) +
  geom_line(lwd = 1.3) +
  facet_wrap(~Region) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 18) +
  labs(x = 'Year', y = 'Age-2 Recruitment (millions)')

# movement
move_plot <- reshape2::melt(model[[4]]$rep$Movement) %>%
  rename(From = from, To = to, Year = years, Seas = seas,
         Age = ages, Sex = sexes) %>%
  filter(Age %in% c(3, 10, 30), Year == 1) %>%
  mutate(
    Age = case_when(Age == 3 ~ "Age Block 1 (Young)",
                    Age == 10 ~ "Age Block 2 (Intermediate)",
                    Age == 30 ~ "Age Block 3 (Old)"),
    From = case_when(
      From == 1 ~ "BS + AI + WGOA",
      From == 2 ~ "CGOA",
      From == 3 ~ "EGOA"
    ),
    To = case_when(
      To == 1 ~ "BS + AI + WGOA",
      To == 2 ~ "CGOA",
      To == 3 ~ "EGOA"
    )
  ) %>%
  ggplot(aes(x = From, y = To, fill = value, label = round(value, 2))) +
  geom_tile(alpha = 0.5) +
  geom_text(size = 7) +
  scale_fill_viridis_c() +
  theme_bw(base_size = 18) +
  facet_wrap(~Age) +
  theme(legend.position = 'top',
        legend.key.width = unit(1, 'cm')) +
  labs(x = "From", y = 'To', fill = 'Movement Probability')

rec_move_plot <- cowplot::plot_grid(
  rec_plot, move_plot, ncol = 1, align = 'v', axis = 'l',
  labels = c("A", 'B'), label_size = 30, hjust = -2
)

ggsave(
  here("dev", "paper_projects", "spatial_plus_group", "sablefish_rec_move_plot.png"),
  rec_move_plot,
  width = 15, height = 15
)

