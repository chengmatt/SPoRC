# Purpose: Fit the five-region spatial Alaska sablefish model and render the spatial case study figures
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26
#
# no operational counterpart, so nothing here is a bridge. slow: most of the runtime is MakeADFun
# tape construction, which scales with tag cohorts times maximum tag liberty. cohorts are release EVENTS

library(here)
library(dplyr)
library(ggplot2)
devtools::load_all(here())

dat <- SPoRC::mlt_rg_sable_data
yrs <- dat$years
n_yrs <- length(dat$years)
n_ages <- length(dat$ages)
n_regions <- dat$n_regions

region_lab <- c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')
region_lvl <- c("BS", "AI", "WGOA", "CGOA", "EGOA")

# Setup ----------------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  years = 1:length(dat$years),
  ages = 1:length(dat$ages),
  lens = dat$lens,
  n_regions = dat$n_regions,
  n_sexes = dat$n_sexes,
  n_fish_fleets = dat$n_fish_fleets,
  n_srv_fleets = dat$n_srv_fleets,
  n_pop = dat$n_pop,
  verbose = FALSE
)

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  rec_model = "mean_rec",
  do_rec_bias_ramp = 0,
  sigmaR_switch = 16,
  dont_est_recdev_last = 1,
  sigmaR_spec = "fix",
  # initial deviations shared across regions: estimating both these and regional recruitment
  # deviations leaves the initial condition underdetermined, since nothing before year 1 separates them
  InitDevs_spec = "est_shared_r",
  ln_sigmaR = array(log(c(0.4, 1.2)), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
  ln_global_R0 = log(20),
  # alpha = 3 everywhere puts the Dirichlet mode at EQUAL recruitment across
  # regions with a concentration of sum(alpha) - n_r = 10 pseudo-observations.
  # This is informative, not vague; alpha = 1 would be uniform over the simplex.
  use_rec_region_prop_prior = 1,
  rec_region_prop_prior = data.frame(pop = 1, alpha = I(list(rep(3, input_list$data$n_regions)))),
  rec_region_prop_pars = array(c(0.2, 0.2, 0.2, 0.2),
                               dim = c(input_list$data$n_pop, input_list$data$n_regions - 1))
)

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA = dat$WAA,
  MatAA = dat$MatAA,
  AgeingError = dat$AgeingError,
  fit_lengths = 1,
  SizeAgeTrans = dat$SizeAgeTrans,
  # M and movement both remove fish from a region, so they are confounded once
  # movement is estimated. Fixing M removes a confounding the data cannot solve.
  M_spec = "fix",
  Fixed_natmort = array(0.104884, dim = c(dat$n_pop, dat$n_regions, n_yrs, n_ages, dat$n_sexes))
)

# Movement -------------------------------------------------------------------
# one prior row per origin region, age block (one age standing for each block) and sex.
# I() keeps expand.grid from splitting the alpha list across rows
Movement_prior <- expand.grid(
  pop = 1,
  region_from = 1:5,
  year = 1,
  seas = 1,
  age = c(6, 7, 16),
  sex = 1,
  alpha = I(list(rep(3, 5)))
)

input_list <- Setup_Mod_Movement(
  input_list = input_list,
  # juveniles, maturing fish, adults
  Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)),
  Movement_yearblk_spec = "constant",
  Movement_sexblk_spec = "constant",
  # recruits moving is confounded with recruitment apportionment: both determine
  # where age-1 fish are found
  do_recruits_move = 0,
  use_fixed_movement = 0,
  Use_Movement_Prior = 1,
  Movement_prior = Movement_prior
)

# Tagging --------------------------------------------------------------------
tag_prior <- data.frame(
  region = 1,
  block = c(1, 2),
  fleet = 1,
  mu = NA,  # symmetric beta has no mean
  sd = 5,
  type = 0
)

input_list <- Setup_Mod_Tagging(
  input_list = input_list,
  use_conv_fish_tagging = c(1, 0),
  # dominant driver of tape build cost
  conv_tag_max_liberty = 15,
  conv_tag_release_indicator = dat$conv_tag_release_indicator,
  conv_tagged_fish = dat$conv_tagged_fish,
  obs_conv_tag_fish_recap = dat$obs_conv_tag_fish_recap,
  conv_fish_tag_like = "Multinomial_Release",
  # newly released fish have not mixed and would read as extreme residency
  conv_tag_mixing_period = 2,
  conv_tag_t_tagging = 0.5,
  use_conv_tag_fishrep_prior = 1,
  conv_tag_fishrep_prior = tag_prior,
  conv_tag_age_pool = as.list(1:30),
  conv_tag_sex_pool = list(c(1:2)),
  # both confounded with natural mortality, which is itself fixed
  init_conv_tag_mort_spec = "fix",
  conv_tag_shed_spec = "fix",
  conv_tagrep_spec = "est_shared_r_f",
  conv_tag_fish_reporting_blocks = c(
    apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
      paste0("Block_1_Year_1-35_Region_", x[1], "_Fleet_", x[2])),
    apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
      paste0("Block_2_Year_36-terminal_Region_", x[1], "_Fleet_", x[2]))
  ),
  conv_fish_tag_attr = 'p_a_s',
  ln_init_conv_tag_mort = log(0.1),
  ln_conv_tag_shed = log(0.02),
  ln_conv_fish_tag_theta = log(0.5),
  conv_tag_fish_reporting_pars = array(log(0.2 / (1 - 0.2)),
                                       dim = c(input_list$data$n_regions, 2,
                                               input_list$data$n_fish_fleets))
)

input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,
  ObsCatch = dat$ObsCatch,
  UseCatch = dat$UseCatch,
  Use_F_pen = 1,
  sigmaC_spec = 'fix',
  ln_sigmaC = array(log(0.05), dim = c(input_list$data$n_regions, n_yrs,
                                       input_list$data$n_seas,
                                       input_list$data$n_fish_fleets))
)

# ln_F_mean cannot be passed through ... because R partially matches the name to
# the ln_F_mean_spec formal, so the starting value is assigned post-hoc.
input_list$par$ln_F_mean[] <- -2

# spltRjntS: each region's compositions sum to one across both sexes, so they hold sex-ratio
# information within a region but say nothing about relative abundance between regions
input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  ObsFishIdx = dat$ObsFishIdx,
  ObsFishIdx_SE = dat$ObsFishIdx_SE,
  UseFishIdx = dat$UseFishIdx,
  ObsFishAgeComps = dat$ObsFishAgeComps,
  UseFishAgeComps = dat$UseFishAgeComps,
  ISS_FishAgeComps = dat$ISS_FishAgeComps,
  ObsFishLenComps = dat$ObsFishLenComps,
  UseFishLenComps = dat$UseFishLenComps,
  ISS_FishLenComps = dat$ISS_FishLenComps,
  fish_idx_type = c("none", "none"),
  FishAgeComps_LikeType = c("Multinomial", "none"),
  FishLenComps_LikeType = c("Multinomial", "Multinomial"),
  FishAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                        "none_Year_1-terminal_Fleet_2"),
  FishLenComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                        "spltRjntS_Year_1-terminal_Fleet_2")
)

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,
  ObsSrvIdx = dat$ObsSrvIdx,
  ObsSrvIdx_SE = dat$ObsSrvIdx_SE,
  UseSrvIdx = dat$UseSrvIdx,
  ObsSrvAgeComps = dat$ObsSrvAgeComps,
  ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
  UseSrvAgeComps = dat$UseSrvAgeComps,
  ObsSrvLenComps = dat$ObsSrvLenComps,
  UseSrvLenComps = dat$UseSrvLenComps,
  ISS_SrvLenComps = dat$ISS_SrvLenComps,
  srv_idx_type = c("abd", "abd"),
  SrvAgeComps_LikeType = c("Multinomial", "Multinomial"),
  SrvLenComps_LikeType = c("none", "none"),
  SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                       "spltRjntS_Year_1-terminal_Fleet_2"),
  SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                       "none_Year_1-terminal_Fleet_2"),
  t_srv = array(0.5, dim = c(input_list$data$n_regions,
                             input_list$data$n_seas,
                             input_list$data$n_srv_fleets))
)

# Selectivity ----------------------------------------------------------------
# Prior rows MUST match the mapping: a prior on a mapped-off parameter adds a
# constant to the objective and a gradient to nothing.
sex_par <- expand.grid(sex = 1:2, par = 1:2)
fleet_blocks <- data.frame(fleet = c(1, 2), block = 1)

fish_selex_structure <- merge(fleet_blocks, sex_par) %>%
  dplyr::filter(!(fleet == 1 & block == 1 & sex == 2 & par == 2)) %>%
  dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 1))

fish_selex_prior <- cbind(region = 1, fish_selex_structure, mu = 2, sd = 3)

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list = input_list,
  cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
  fish_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),
  fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  # spatially-invariant selectivity: regional differences in the compositions are
  # attributed to age structure, which movement and apportionment produce
  fish_fixed_sel_pars = c("est_shared_r", "est_shared_r"),
  fish_q_spec = c("fix", "fix"),
  Use_fish_selex_prior = 1,
  fish_selex_prior = fish_selex_prior
)

map_fish_fixed <- array(input_list$map$fish_fixed_sel_pars, dim = dim(input_list$par$fish_fixed_sel_pars))
map_fish_fixed[, 2, 1, 2, 1] <- map_fish_fixed[, 2, 1, 1, 1] # fixed-gear delta across sexes
map_fish_fixed[, 1, 1, 2, 2] <- map_fish_fixed[, 1, 1, 1, 2] # trawl bmax across sexes
input_list$map$fish_fixed_sel_pars <- factor(map_fish_fixed)

srv_selex_prior <- cbind(
  region = 1,
  merge(data.frame(fleet = c(1, 2), block = c(1, 1)), sex_par),
  mu = 1,
  sd = 5
) %>%
  dplyr::filter(!(fleet == 2 & par == 2 & sex == 2)) %>%
  dplyr::mutate(mu = ifelse(fleet == 2, 2, mu),
                sd = ifelse(fleet == 2, 3, sd))

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2"),
  srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  srv_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2"),
  srv_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  srv_fixed_sel_pars_spec = c("est_shared_r", "est_shared_r"),
  srv_q_spec = c("est_shared_r", "est_shared_r"),
  Use_srv_selex_prior = 1,
  srv_selex_prior = srv_selex_prior
)

map_srv_fixed <- array(input_list$map$srv_fixed_sel_pars, dim = dim(input_list$par$srv_fixed_sel_pars))
map_srv_fixed[, 2, 1, 2, 2] <- map_srv_fixed[, 2, 1, 1, 2] # JP longline delta across sexes
input_list$map$srv_fixed_sel_pars <- factor(map_srv_fixed)

# Weighting ------------------------------------------------------------------
comp_wt <- function(n_fleets) {
  array(1, dim = c(input_list$data$n_regions, n_yrs, input_list$data$n_seas,
                   input_list$data$n_sexes, n_fleets))
} # end comp_wt

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  # tag data comprise a very large number of cells, so their nominal likelihood
  # contribution outruns their information content
  Wt_Tagging = 0.5,
  Wt_FishAgeComps = comp_wt(input_list$data$n_fish_fleets),
  Wt_FishLenComps = comp_wt(input_list$data$n_fish_fleets),
  Wt_SrvAgeComps = comp_wt(input_list$data$n_srv_fleets),
  Wt_SrvLenComps = comp_wt(input_list$data$n_srv_fleets)
)

# Sample sizes total ~100 ACROSS regions, not per region. Splitting comps by
# region multiplies cells by n_r without adding any sampling.
input_list$data$ISS_SrvAgeComps[] <- 20
input_list$data$ISS_FishAgeComps[1, , , , ] <- 25  # BS
input_list$data$ISS_FishAgeComps[2, , , , ] <- 20  # AI
input_list$data$ISS_FishAgeComps[3, , , , ] <- 14  # WGOA
input_list$data$ISS_FishAgeComps[4, , , , ] <- 18  # CGOA
input_list$data$ISS_FishAgeComps[5, , , , ] <- 18  # EGOA
input_list$data$ISS_FishLenComps[1, , , , ] <- 12  # BS
input_list$data$ISS_FishLenComps[2, , , , ] <- 12  # AI
input_list$data$ISS_FishLenComps[3, , , , ] <-  7  # WGOA
input_list$data$ISS_FishLenComps[4, , , , ] <-  7  # CGOA
input_list$data$ISS_FishLenComps[5, , , , ] <-  7  # EGOA

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Fit ------------------------------------------------------------------------
est <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 5, silent = FALSE)
est$sd_rep <- RTMB::sdreport(est)

cat("\n=== Fit ===\n")
cat("free parameters:", length(est$par), "\n")
cat("final jnLL:", est$optim$objective, "  max |gradient|:", max(abs(est$gr(est$optim$par))), "\n")

rep <- est$rep

# Regional trajectories ------------------------------------------------------
ts_df <- bind_rows(
  reshape2::melt(rep$SSB) %>% dplyr::mutate(Par = "Spawning Biomass"),
  reshape2::melt(rep$Rec) %>% dplyr::mutate(Par = "Recruitment")
) %>%
  dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  dplyr::mutate(Region = factor(region_lab[Region], levels = region_lvl),
                Year = Year + 1959)

p_ts <- ggplot2::ggplot(ts_df, ggplot2::aes(x = Year, y = value, color = Region)) +
  ggplot2::geom_line(linewidth = 1.3) +
  ggplot2::facet_grid(Par ~ Region, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  ggplot2::labs(x = "Year", y = "Value") +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(legend.position = 'none')

ggplot2::ggsave(
  here("vignettes", "figures", "g_ts_comparison.png"),
  p_ts,
  width = 12,
  height = 5.5,
  dpi = 150
)

# Movement -------------------------------------------------------------------
# Movement is [pop, from, to, year, season, age, sex]; year and sex are constant
# blocks here, so the first index of each stands for all of them.
move_df <- reshape2::melt(rep$Movement[1, , , 1, 1, , 1]) %>%
  dplyr::rename(From = Var1, To = Var2, Age = Var3) %>%
  dplyr::filter(Age %in% c(6, 7, 16)) %>%
  dplyr::mutate(dplyr::across(c(From, To), ~ factor(region_lab[.x], levels = region_lvl)),
                Age = factor(Age, levels = c(6, 7, 16),
                             labels = c("Ages 1-6", "Ages 7-15", "Ages 16-30")))

p_move <- ggplot2::ggplot(move_df, ggplot2::aes(x = To, y = From, fill = value)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 3.2) +
  ggplot2::facet_wrap(~Age) +
  ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
  ggplot2::labs(x = "To region", y = "From region", fill = "Movement\nprobability") +
  ggplot2::theme_bw(base_size = 13)

ggplot2::ggsave(
  here("vignettes", "figures", "g_movement.png"),
  p_move,
  width = 12,
  height = 4.5,
  dpi = 150
)

cat("\n=== Regional recruitment apportionment ===\n")
print(data.frame(Region = region_lvl,
                 proportion = as.vector(rep$Rec_prop[1, ])),
      row.names = FALSE, digits = 3)

cat("\n=== Residency by age block (diagonal of the movement matrix) ===\n")
print(move_df %>%
        dplyr::filter(From == To) %>%
        dplyr::select(Age, Region = From, residency = value) %>%
        tidyr::pivot_wider(names_from = Age, values_from = residency),
      row.names = FALSE, digits = 3)

saveRDS(list(
  rep = rep,
  sdrep = est$sd_rep,
  opt = est$optim,
  data = data,
  parameters = parameters,
  mapping = mapping
),
        here("dev", "dev_output", "mlt_rg_sablefish_fit.rds"))
