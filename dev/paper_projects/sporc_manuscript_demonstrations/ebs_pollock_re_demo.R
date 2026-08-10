# Purpose: To demonstrate the use of SPoRC in a single region context for BSAI pollock
# Creator: Matthew LH. Cheng
# Date Created: 5/1/25

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)
library(ebswp)
devtools::load_all(here("R"))
data("sgl_rg_ebswp_data")

# Prepare Data and Inputs -------------------------------------------------
#' Setup Single-Region Population Model for EBS Pollock
#'
#' Constructs a single-region population model input list, tailored to the
#' \code{sgl_rg_ebswp_data} dataset. This function initializes dimensions,
#' recruitment, natural mortality, biologicals, movement, tagging, catch,
#' fishery indices and compositions, survey indices and compositions,
#' selectivity and catchability, and component weighting.
#'
#' @param cont_tv_fish_sel Character vector. Whether to estimate continuous
#' time-varying fishery selectivit (see \code{\link{Setup_Mod_FishIdx_and_Comps}}).
#' @param fishsel_pe_pars_spec Character vector. Specification for penalized
#' likelihood parameters for fishery selectivity deviations (see \code{\link{Setup_Mod_FishIdx_and_Comps}}).
#' @param corr_opt_semipar Character vector. Correlation options for
#' semi-parametric selectivity (see \code{\link{Setup_Mod_FishIdx_and_Comps}}).
#' @param fish_sel_devs_spec Character vector. Specification of fishery
#' selectivity deviations to be estimated (see \code{\link{Setup_Mod_FishIdx_and_Comps}}).
#'
#' @details
#' The function relies on the global dataset \code{sgl_rg_ebswp_data} to provide
#' years, ages, weight-at-age, maturity-at-age, observed catches, fishery and
#' survey indices, age and length compositions, and other required inputs.
#'
#' The model is configured for:
#' - one region
#' - one sex
#' - one fishery fleet
#' - four survey fleets
#'
#' Survey fleet 1 is the bottom trawl survey, 2 the acoustic trawl survey, 3 the
#' acoustic vessel of opportunity index, and 4 the acoustic survey's age 1
#' abundance, which the assessment fits as an index of its own.
#'
#' Natural mortality (\eqn{M}) is fixed at:
#' - 0.9 for age-1,
#' - 0.45 for age-2,
#' - 0.3 for age-3+.
#'
#' Recruitment is modeled using a Beverton–Holt stock–recruitment function
#' with fixed steepness.
#'
#' Selectivity and catchability are parameterized separately for fishery and
#' survey fleets, with user control over random effects and correlation
#' structure for fishery selectivity deviations.
#'
#' @return A fully specified model input list to pass onto subsequent model fitting functions.
#'
#' @seealso
#' \code{\link{Setup_Mod_Dim}}, \code{\link{Setup_Mod_Rec}},
#' \code{\link{Setup_Mod_Biologicals}}, \code{\link{Setup_Mod_Movement}},
#' \code{\link{Setup_Mod_Tagging}}, \code{\link{Setup_Mod_Catch_and_F}},
#' \code{\link{Setup_Mod_FishIdx_and_Comps}},
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}},
#' \code{\link{Setup_Mod_Fishsel_and_Q}}, \code{\link{Setup_Mod_Srvsel_and_Q}},
#' \code{\link{Setup_Mod_Weighting}}
#'
pol_model <- function(cont_tv_fish_sel,
                      fishsel_pe_pars_spec,
                      corr_opt_semipar,
                      fish_sel_devs_spec
) {

  n_srv <- sgl_rg_ebswp_data$n_srv_fleets

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
    n_srv_fleets = n_srv, # number of survey fleets
    n_pop = 1, # number of populations
    verbose = FALSE
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
    ln_sigmaR = array(log(c(5, 1)), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    # Starting values for early and late sigmaR
    rec_model = "bh_rec",
    # recruitment model
    steepness_h = array(inv_steepness(0.623013), dim = c(input_list$data$n_pop, input_list$data$n_regions)),
    h_spec = "fix",
    # fixing steepness
    sigmaR_spec = "fix",
    # fix early sigmaR and late sigmaR
    init_age_strc = 1,
    ln_global_R0 = 10,
    t_spawn = 0.25,
    equil_init_age_strc = 2
  )

  # Setup a fixed natural mortality array for use
  fix_natmort <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), 1))
  fix_natmort[,,,1,] <- 0.9 # age 1 M
  fix_natmort[,,,2,] <- 0.45 # age 2 M
  fix_natmort[,,,-c(1,2),] <- 0.3 # age 3+ M

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,

    # Data inputs
    WAA = sgl_rg_ebswp_data$WAA,
    # the assessment carries a separate weight at age matrix for the fishery
    # and for each survey index
    WAA_fish = sgl_rg_ebswp_data$WAA_fish,
    WAA_srv = sgl_rg_ebswp_data$WAA_srv,
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

  # Setup tagging stuff
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

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
    ObsSrvLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, n_srv)),
    UseSrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, n_srv)),
    ISS_SrvLenComps = NULL,

    # Model options
    srv_idx_type = c("biom", "biom", "biom", "abd"),
    # biomass for survey fleets 1 to 3, abundance for survey fleet 4
    srv_idx_ages = list(NULL, NULL, NULL, 1),
    # fleet 4 is the acoustic survey's age 1 abundance, so it sees age 1 only
    SrvIdx_LikeType = c("mvn", "lognormal", "normal", "lognormal"),
    # index likelihood for survey fleet 1, 2, 3, and 4
    SrvIdx_Cov = list(sgl_rg_ebswp_data$SrvIdx_Cov, NULL, NULL, NULL),
    # the bottom trawl index is fit with a full covariance matrix
    SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "none", "none"),
    # survey age composition likelihood for survey fleet 1, 2, 3, and 4
    SrvAgeComps_bins = list(NULL, 2:15, NULL, NULL),
    # the acoustic compositions are normalised over ages 2-15 only
    SrvLenComps_LikeType = rep("none", n_srv),
    #  survey length composition likelihood for survey fleet 1, 2, 3, and 4
    SrvAgeComps_Type = c(
      "agg_Year_1-terminal_Fleet_1",
      "agg_Year_1-terminal_Fleet_2",
      "none_Year_1-terminal_Fleet_3",
      "none_Year_1-terminal_Fleet_4"
    ),
    # survey age comp type

    SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:n_srv),
    # survey length comp type
    t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv))
    # fraction of the year elapsed when each survey occurs
  )


  # Setup fishery selectivity and catchability
  input_list <- Setup_Mod_Fishsel_and_Q(

    input_list = input_list,

    # Model options (NOTE: Iterating Different Fishery Selectivity Random Effects Here!)
    cont_tv_fish_sel = cont_tv_fish_sel,  # fishery selectivity, whether continuous time-varying
    fishsel_pe_pars_spec = fishsel_pe_pars_spec, # doing penalized likelihood for selex devs
    fish_sel_devs_spec = fish_sel_devs_spec, # estimating all sel devs
    corr_opt_semipar = corr_opt_semipar, # correlation options

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
    # survey selectivity blocks
    srv_sel_blocks = paste0("none_Fleet_", 1:n_srv),
    # survey selectivity form
    srv_sel_model = paste0("logist1_Fleet_", 1:n_srv),
    # survey catchability blocks
    srv_q_blocks = paste0("none_Fleet_", 1:n_srv),
    # whether to estiamte all fixed effects for survey selectivity. The vessel
    # of opportunity index shares the acoustic survey's, and fleet 4 sees age 1
    # only, where selectivity is absorbed into catchability and so is fixed.
    srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2", "fix"),
    # whether to estiamte all fixed effects for survey catchability
    srv_q_spec = rep("est_all", n_srv),
    t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv))
  )

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets))
  )

  return(input_list)
}


# model storage
models <- list()

# models to iterate through
pol_model_var <- data.frame(
  cont_tv_fish_sel = c("none_Fleet_1", "iid_Fleet_1", "rw_Fleet_1", "2dar1_Fleet_1", "3dcond_Fleet_1", "3dcond_Fleet_1"),
  fishsel_pe_pars_spec = c("none", rep("est_all", 5)),
  fish_sel_devs_spec = c("none", rep("est_all", 5)),
  corr_opt_semipar = c(rep(NA, 5), "corr_zero_y_b_c")
  )

# loop through models
for(i in 1:nrow(pol_model_var)) {

  # set up random effects
  if(str_detect(pol_model_var$cont_tv_fish_sel[i], "none")) random <- NULL
  else random <- "ln_fishsel_devs"

  # get input list
  input_list <- pol_model(cont_tv_fish_sel = pol_model_var$cont_tv_fish_sel[i],
                          fishsel_pe_pars_spec = pol_model_var$fishsel_pe_pars_spec[i],
                          fish_sel_devs_spec = pol_model_var$fish_sel_devs_spec[i],
                          corr_opt_semipar = if(is.na(pol_model_var$corr_opt_semipar[i])) NULL else pol_model_var$corr_opt_semipar[i]
                          )

  # extract out lists updated with helper functions
  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # Fit model
  ebswp_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = random,
                                newton_loops = 3,
                                silent = FALSE
  )

  ebswp_rtmb_model$sdrep <- RTMB::sdreport(ebswp_rtmb_model)
  sdrep <- ebswp_rtmb_model$sdrep
  rep <- ebswp_rtmb_model$rep
  models[[i]] <- ebswp_rtmb_model

}

# save models
saveRDS(models, here("dev", "paper_projects",  "sporc_manuscript_demonstrations", "pollock_re_comparison", "pol_models_re.RDS"))

# Plot! -------------------------------------------------------------------
models <- readRDS(here("dev", "paper_projects",  "sporc_manuscript_demonstrations", "pollock_re_comparison", "pol_models_re.RDS"))
model_names <- c("constant", "iid_p", "rw_p", "2dar1_sp", "3dgmrf_sp", "iid_sp")
fishsel_all_df <- data.frame() # empty dataframe to bind to
ts_all_df <- data.frame() # empty dataframe to bind to
for(i in 1:length(models)) {

  # Get recruitment time-series
  rec_series <- reshape2::melt((models[[i]]$rep$Rec)) %>%
    mutate(se = models[[i]]$sdrep$sd[names(models[[i]]$sdrep$value) == 'log_Rec'] * as.vector(models[[i]]$rep$Rec))
  rec_series$Par <- "Recruitment"
  rec_series$Model <- model_names[i]

  # Get SSB time-series
  ssb_series <- reshape2::melt((models[[i]]$rep$SSB)) %>%
    mutate(se = models[[i]]$sdrep$sd[names(models[[i]]$sdrep$value) == 'log_SSB'] * as.vector(models[[i]]$rep$SSB))
  ssb_series$Par <- "Spawning Stock Biomass"
  ssb_series$Model <- model_names[i]

  # Get fishery selectivity estimates
  fishsel_df <- reshape2::melt(models[[i]]$rep$fish_sel) %>%
    rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Age = Var5, Sex = Var6, Fleet = Var7) %>%
    group_by(Year, Region, Sex, Fleet) %>%
    mutate(value = value/max(value),
           Year = Year + 1963)
  fishsel_df$Model <- model_names[i]

  # bind together
  ts_df <- rbind(ssb_series,rec_series) %>%
    dplyr::rename(Region = Var2, Year = Var3) %>%
    dplyr::mutate(Year = Year + 1963)

  ts_all_df <- rbind(ts_all_df, ts_df)
  fishsel_all_df <- rbind(fishsel_df, fishsel_all_df)

} # end i loop

# Refactor to order models
fishsel_all_df <- fishsel_all_df %>% mutate(Model = factor(Model, levels = c("constant", "iid_p", "rw_p", "iid_sp", "2dar1_sp", "3dgmrf_sp")))
ts_all_df <- ts_all_df %>% mutate(Model = factor(Model, levels = c("constant", "iid_p", "rw_p", "iid_sp", "2dar1_sp", "3dgmrf_sp")))

cols <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7") # colors

# time series plot
ts_plot <- ggplot(ts_all_df, aes(x = Year, y = value,
                                 ymin = value - (1.96 * se), ymax = value + (1.96 * se),
                                 color = Model, fill = Model, lty = Model)) +
  geom_line(lwd = 1.3) +
  facet_wrap(~Par, scales = 'free', ncol = 1) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(y = "Value")  +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top',
        legend.background = element_blank()) +
  ylim(0, NA) +
  labs(x = 'Year', y = 'Value', color = 'Model', fill = 'Model')

# selectivity plot
selex_plot <- ggplot(fishsel_all_df, aes(x = Year, y = Age, fill = value)) +
  geom_tile() +
  scale_fill_continuous(palette = "magma") +
  facet_wrap(~Model, scales = 'free') +
  theme_bw(base_size = 15) +
  labs(x = 'Year', y = 'Age', fill = 'Relative Selectivity') +
  theme(legend.position = "top",
        legend.key.size = unit(1, "cm"),
        legend.key.height = unit(0.2, "cm"))

png(here("vignettes", "figures", "n_ebs_pol_re.png"), width = 800, height = 1e3)
ts_plot
dev.off()

png(here("vignettes", "figures", "n_fishsel_re.png"), width = 1000, height = 500)
selex_plot
dev.off()

ggsave(
  here("dev", "paper_projects",  "sporc_manuscript_demonstrations", "figs", "ebs_pol_randef.png"),
  cowplot::plot_grid(
    ts_plot,
    selex_plot,
    labels = c('A', 'B'),
    label_size = 30,
    rel_widths = c(0.4, 0.6), align = 'h', vjust = 4
  ),
  height = 8, width = 15
)

# Presentation plots
# time series plot
ts_plot_pres <- ggplot(ts_all_df %>% filter(Model %in% c('constant', 'rw_p', '2dar1_sp', '3dgmrf_sp')),
                       aes(x = Year, y = value, ymin = value - (1.96 * se), ymax = value + (1.96 * se),
                                 color = Model, fill = Model, lty = Model)) +
  geom_line(lwd = 1.3) +
  facet_wrap(~Par, scales = 'free', ncol = 1) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(y = "Value")  +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top',
        legend.background = element_blank()) +
  ylim(0, NA) +
  labs(x = 'Year', y = 'Value', color = 'Model', fill = 'Model')

# selectivity plot
selex_plot_pres <- ggplot(fishsel_all_df %>% filter(Model %in% c('constant', 'rw_p', '2dar1_sp', '3dgmrf_sp')),
                     aes(x = Year, y = Age, fill = value)) +
  geom_tile() +
  scale_fill_continuous(palette = "magma") +
  facet_wrap(~Model, scales = 'free') +
  theme_bw(base_size = 15) +
  labs(x = 'Year', y = 'Age', fill = 'Relative Selectivity') +
  theme(legend.position = "top",
        legend.key.size = unit(1, "cm"),
        legend.key.height = unit(0.2, "cm"))

ggsave(
  here("dev", "paper_projects",  "sporc_manuscript_demonstrations", "figs", "ebs_pol_ts.png"),
  ts_plot_pres,
  height = 9, width = 8
)

ggsave(
  here("dev", "paper_projects",  "sporc_manuscript_demonstrations", "figs", "ebs_pol_selex.png"),
  selex_plot_pres,
  height = 6, width = 7
)

ggsave(
  here("dev", "paper_projects",  "sporc_manuscript_demonstrations", "figs", "ebs_pol_selex_zoom.png"),
  ggplot(fishsel_all_df %>% filter(Model %in% c('constant', 'rw_p', '2dar1_sp', '3dgmrf_sp'), Year == 2005) %>%
           group_by(Model) %>% mutate(value = value / max(value)),
         aes(x = Age, y = value, color = Model)) +
    geom_line(lwd = 1) +
    scale_color_manual(values = cols) +

    theme_bw(base_size = 14) +
    labs(x = 'Age', y = 'Relative Selectivity') +
    theme(legend.position = "top",
          legend.key.size = unit(1, "cm"),
          legend.key.height = unit(0.2, "cm")),
  height = 6, width = 6
)
