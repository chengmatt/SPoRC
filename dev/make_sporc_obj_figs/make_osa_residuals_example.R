# Purpose: Demonstrate external (post-hoc) and internal (model-based) OSA
#          residual diagnostics via get_osa()/plot_resids()
# Creator: Matthew LH. Cheng (UAF-CFOS)

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(cowplot)
devtools::load_all(here("R"))

# External vs. internal OSA on the same real dataset (single-region GOA Dusky Rockfish)
data("sgl_rg_dusky_data")

# Builds the dusky input_list; `do_internal_comp_osa`
setup_dusky <- function(do_internal_comp_osa = FALSE) {

  input_list <- Setup_Mod_Dim(
    years = sgl_rg_dusky_data$years,
    ages = sgl_rg_dusky_data$mod_ages,
    lens = sgl_rg_dusky_data$lens,
    n_regions = sgl_rg_dusky_data$n_regions,
    n_sexes = sgl_rg_dusky_data$n_sexes,
    n_fish_fleets = sgl_rg_dusky_data$n_fish_fleets,
    n_srv_fleets = sgl_rg_dusky_data$n_srv_fleets,
    n_seas = sgl_rg_dusky_data$n_seas,
    n_pop = sgl_rg_dusky_data$n_pop,
    natal_region = sgl_rg_dusky_data$natal_region,
    verbose = FALSE,
    do_internal_comp_osa = do_internal_comp_osa
  )

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 1,
    bias_year = rep(length(sgl_rg_dusky_data$years), 4),
    sigmaR_switch = 1,
    ln_sigmaR = array(-0.1068576, dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    ln_global_R0 = log(2.7),
    t_spawn = sgl_rg_dusky_data$spwn_month
  )

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = sgl_rg_dusky_data$waa_arr,
    MatAA = sgl_rg_dusky_data$mataa_arr,
    fit_lengths = 1,
    SizeAgeTrans = sgl_rg_dusky_data$sizeage,
    AgeingError = sgl_rg_dusky_data$age_error_matrix,
    M_spec = "fix",
    Fixed_natmort = sgl_rg_dusky_data$fix_natmort,
    addtocomp = 0.00001
  )

  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sgl_rg_dusky_data$ObsCatch,
    UseCatch = sgl_rg_dusky_data$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    Catch_Constant = 0.00001,
    ln_sigmaC = array(log(sqrt(1 / (2 * c(rep(2, 15), rep(50, 33))) )), dim = c(input_list$data$n_regions, length(input_list$data$years),
                                                                                input_list$data$n_seas, input_list$data$n_fish_fleets)),
    ln_sigmaF = array(log(sqrt(1 / 4)), dim = c(input_list$data$n_regions, input_list$data$n_seas,
                                                input_list$data$n_fish_fleets))
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sgl_rg_dusky_data$ObsFishIdx,
    ObsFishIdx_SE = sgl_rg_dusky_data$ObsFishIdx_SE,
    UseFishIdx = sgl_rg_dusky_data$UseFishIdx,
    ObsFishAgeComps = sgl_rg_dusky_data$ObsFishAgeComps,
    UseFishAgeComps = sgl_rg_dusky_data$UseFishAgeComps,
    ISS_FishAgeComps = sgl_rg_dusky_data$ISS_FishAgeComps,
    ObsFishLenComps = sgl_rg_dusky_data$ObsFishLenComps,
    UseFishLenComps = sgl_rg_dusky_data$UseFishLenComps,
    ISS_FishLenComps = sgl_rg_dusky_data$ISS_FishLenComps,
    fish_idx_type = c("none"),
    FishAgeComps_LikeType = c("Multinomial"),
    FishLenComps_LikeType = c("Multinomial"),
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sgl_rg_dusky_data$ObsSrvIdx,
    ObsSrvIdx_SE = (sgl_rg_dusky_data$ObsSrvIdx_SE / sgl_rg_dusky_data$ObsSrvIdx) / sqrt(1.66),
    UseSrvIdx = sgl_rg_dusky_data$UseSrvIdx,
    ObsSrvAgeComps = sgl_rg_dusky_data$ObsSrvAgeComps,
    ISS_SrvAgeComps = sgl_rg_dusky_data$ISS_SrvAgeComps,
    UseSrvAgeComps = sgl_rg_dusky_data$UseSrvAgeComps,
    ObsSrvLenComps = sgl_rg_dusky_data$ObsSrvLenComps,
    UseSrvLenComps = sgl_rg_dusky_data$UseSrvLenComps,
    ISS_SrvLenComps = sgl_rg_dusky_data$ISS_SrvLenComps,
    srv_idx_type = c("biom"),
    SrvAgeComps_LikeType = c("Multinomial"),
    SrvLenComps_LikeType = c("Multinomial"),
    SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    SrvLenComps_Type = c("agg_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    cont_tv_fish_sel = c("none_Fleet_1"),
    fish_sel_blocks = c("none_Fleet_1"),
    fish_sel_model = c("logist2_Fleet_1"),
    fish_q_blocks = c("none_Fleet_1"),
    fish_fixed_sel_pars_spec = c("est_all"),
    fish_q_spec = c("fix")
  )

  srv_q_prior <- data.frame(region = 1, block = 1, fleet = 1, mu = 1, sd = 0.447213595)

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    cont_tv_srv_sel = c("none_Fleet_1"),
    srv_sel_blocks = c("none_Fleet_1"),
    srv_sel_model = c("logist2_Fleet_1"),
    srv_q_blocks = c("none_Fleet_1"),
    srv_fixed_sel_pars_spec = c("est_all"),
    srv_q_spec = c("est_all"),
    Use_srv_q_prior = 1,
    srv_q_prior = srv_q_prior,
    t_srv = array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_srv_fleets))
  )

  input_list <- Setup_Mod_Weighting(
    input_list
  )


  input_list
}

# External (post-hoc) OSA ----------------------------------------------
input_list_ext <- setup_dusky(do_internal_comp_osa = FALSE)
model_ext <- fit_model(input_list_ext$data, input_list_ext$par, input_list_ext$map,
                       random = NULL, newton_loops = 3, silent = TRUE)
obs_age_bins <- 4:30
comp_prop <- get_comp_prop(input_list_ext$data, model_ext$rep,
                           age_labels = obs_age_bins,
                           len_labels = sgl_rg_dusky_data$lens,
                           year_labels = sgl_rg_dusky_data$years)

fishages_ext <- get_osa(
  obs_mat = comp_prop$Obs_FishAge_mat,
  exp_mat = comp_prop$Pred_FishAge_mat,
  N = input_list_ext$data$ISS_FishAgeComps[1, which(input_list_ext$data$UseFishAgeComps[,,1,1] == 1), 1, 1, 1] *
    unique(input_list_ext$data$Wt_FishAgeComps[1, which(input_list_ext$data$UseFishAgeComps[,,1,1] == 1), 1, , 1]),
  years = which(input_list_ext$data$UseFishAgeComps[,,1,1] == 1),
  fleet = 1,
  bins = obs_age_bins,
  seas = 1,
  comp_type = 0,
  bin_label = "Ages"
)
resid_ext <- plot_resids(fishages_ext)

png(here("vignettes", "figures", "u_external_comp.png"), width = 1000, height = 500)
cowplot::plot_grid(resid_ext[[1]], resid_ext[[2]], ncol = 2)
dev.off()

# Internal (model-based) OSA --------------------------------------------
input_list_int <- setup_dusky(do_internal_comp_osa = TRUE)
model_int <- fit_model(input_list_int$data, input_list_int$par, input_list_int$map,
                       random = NULL, newton_loops = 3, silent = TRUE)

fishages_int <- get_osa(model = model_int, data = input_list_int$data,
                        comp_source = "FishAge", family = "discrete",
                        bins = obs_age_bins, bin_label = "Ages")
resid_int <- plot_resids(fishages_int)

png(here("vignettes", "figures", "u_internal_comp.png"), width = 1000, height = 500)
cowplot::plot_grid(resid_int[[1]], resid_int[[2]], ncol = 2)
dev.off()

resid_catch  <- plot_resids(get_osa(model = model_int, data = input_list_int$data, index_source = "Catch"))
resid_srvidx <- plot_resids(get_osa(model = model_int, data = input_list_int$data, index_source = "SrvIdx"))

png(here("vignettes", "figures", "u_internal_index.png"), width = 1000, height = 800)
cowplot::plot_grid(
  resid_catch[[1]], resid_catch[[2]],
  resid_srvidx[[1]], resid_srvidx[[2]],
  ncol = 2, nrow = 2
)
dev.off()

# Population-specific composition, index, and tagging OSA
set.seed(123)
sim_list <- Setup_Sim_Dim(
  n_sims = 1, n_yrs = 10, n_regions = 2, n_ages = 8, n_lens = NULL,
  n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 2, n_pop = 3,
  natal_region = c(1, 1, 2) # pops 1 & 2 natal to region 1; pop 3 natal to region 2
)

sim_list <- Setup_Sim_Containers(sim_list)

sim_list <- Setup_Sim_Fishing(
  sim_list = sim_list,
  fish_sel_input = replicate(
    n = sim_list$n_sims,
    {
      arr <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                               sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))
      for (r in 1:sim_list$n_regions) for (y in 1:sim_list$n_yrs) for (s in 1:sim_list$n_sexes)
        for (p in 1:sim_list$n_pop) for (seas in 1:sim_list$n_seas)
          arr[p, r, y, seas, , s, 1] <- 1 / (1 + exp(-1.5 * (1:sim_list$n_ages - 3)))
      arr
    }
  ),
  Fmort_input = {
    n <- sim_list$n_yrs * sim_list$n_seas * sim_list$n_sims * sim_list$n_fish_fleets
    t <- seq(0, 2 * pi, length.out = n)
    arr <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                             sim_list$n_fish_fleets, sim_list$n_sims))
    arr[1, , , , ] <- 0.15 * exp(sin(t) + rnorm(n, 0, 0.1))
    arr[2, , , , ] <- 0.05 * exp(-sin(t) + rnorm(n, 0, 0.1))
    arr
  },
  ISS_FishAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                        sim_list$n_fish_fleets, sim_list$n_sims)),
  ISS_FishAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                    sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                    sim_list$n_fish_fleets, sim_list$n_sims)),
  ln_sigmaC = array(log(0.01), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
  ln_sigmaC_pop = array(log(0.01), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
)

sim_list <- Setup_Sim_Survey(
  sim_list = sim_list,
  srv_sel_input = replicate(
    n = sim_list$n_sims,
    {
      arr <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                               sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))
      for (r in 1:sim_list$n_regions) for (y in 1:sim_list$n_yrs) for (s in 1:sim_list$n_sexes)
        for (p in 1:sim_list$n_pop) for (seas in 1:sim_list$n_seas)
          arr[p, r, y, seas, , s, 1] <- 1 / (1 + exp(-1 * (1:sim_list$n_ages - 2.5)))
      arr
    }
  ),
  ISS_SrvAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                       sim_list$n_srv_fleets, sim_list$n_sims)),
  ISS_SrvAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                   sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                   sim_list$n_srv_fleets, sim_list$n_sims)),
  ObsSrvIdx_SE = array(0.15, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets)),
  ObsSrvIdx_pop_SE = array(0.15, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))
)

sim_list <- Setup_Sim_Biologicals(
  sim_list = sim_list,
  natmort_input = array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages,
                                     sim_list$n_sexes, sim_list$n_sims)),
  WAA_input = replicate(
    n = sim_list$n_sims,
    array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas, times = sim_list$n_sexes),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))
  ),
  WAA_fish_input = replicate(
    n = sim_list$n_sims,
    array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas, times = sim_list$n_sexes * sim_list$n_fish_fleets),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))
  ),
  WAA_srv_input = replicate(
    n = sim_list$n_sims,
    array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas, times = sim_list$n_sexes * sim_list$n_srv_fleets),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))
  ),
  MatAA_input = replicate(
    n = sim_list$n_sims,
    array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas, times = sim_list$n_sexes),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))
  )
)

sim_list <- Setup_Sim_Tagging(
  sim_list = sim_list, use_conv_fish_tagging = 1,
  n_tags = 1e3, conv_tag_max_liberty = 8, conv_fish_tag_like = "Poisson"
)

sim_list$Movement <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
                                      sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                                      sim_list$n_sexes, sim_list$n_sims))
stay_prob <- c(0.7, 0.3, 0.7)
disperse_prob <- (1 - stay_prob) / (sim_list$n_regions - 1)
non_natal_rate <- 0.15
for (p in seq_len(sim_list$n_pop)) {
  nr <- sim_list$natal_region[p]
  for (r_from in seq_len(sim_list$n_regions)) {
    for (r_to in seq_len(sim_list$n_regions)) { # Season 1: dispersal
      prob <- if (r_to == r_from) stay_prob[p] else disperse_prob[p]
      sim_list$Movement[p, r_from, r_to, , 1, , , ] <- prob
    }
    for (r_to in seq_len(sim_list$n_regions)) { # Season 2: natal return + straying
      prob <- if (r_to == nr) 1 - non_natal_rate * (sim_list$n_regions - 1) else non_natal_rate
      sim_list$Movement[p, r_from, r_to, , 2, , , ] <- prob
    }
  }
}

sim_list <- Setup_Sim_Rec(
  sim_list = sim_list,
  R0_input = {
    arr <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
    arr[1, 1, , ] <- 7; arr[2, 1, , ] <- 7; arr[3, 2, , ] <- 7
    arr
  },
  ln_sigmaR = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
  init_age_strc = "matrix", recruitment_opt = "bh_rec", rec_dd = "local",
  spawn_seas = 2, t_spawn = 0.5,
  h_input = {
    arr <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
    arr[1, 1, , ] <- 0.75; arr[2, 1, , ] <- 0.75; arr[3, 2, , ] <- 0.75
    arr
  }
)

sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_yrs, sim = 1)

input_pop <- Setup_Mod_Dim(
  years = 1:sim_obj$n_years, ages = 1:sim_obj$n_ages, lens = sim_obj$n_lens,
  n_regions = sim_obj$n_regions, n_sexes = sim_obj$n_sexes,
  n_fish_fleets = sim_obj$n_fish_fleets, n_srv_fleets = sim_obj$n_srv_fleets,
  n_seas = sim_obj$n_seas, n_pop = sim_obj$n_pop, seasdur = sim_obj$seasdur,
  natal_region = c(1, 1, 2), verbose = FALSE,
  do_internal_comp_osa = TRUE, do_internal_conv_tag_osa = TRUE
)

input_pop <- Setup_Mod_Rec(
  input_list = input_pop, do_rec_bias_ramp = 0, sigmaR_switch = 1,
  init_age_strc = "matrix", equil_init_age_strc = "stoch_all",
  spawn_seas = sim_obj$spawn_seas, t_spawn = sim_obj$t_spawn,
  rec_model = "bh_rec", sigmaR_spec = "fix", rec_dd = "local",
  InitDevs_spec = "est_shared_r", RecDevs_spec = "est_shared_r",
  sexratio_spec = "fix", rec_region_prop_spec = "no_dispersal", h_spec = "fix",
  steepness_h = {
    arr <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions))
    arr[1, 1] <- qlogis((0.75 - 0.2) / 0.8); arr[2, 1] <- qlogis((0.75 - 0.2) / 0.8); arr[3, 2] <- qlogis((0.75 - 0.2) / 0.8)
    arr
  },
  ln_sigmaR = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
  ln_global_R0 = array(c(log(7), log(5), log(10)), dim = sim_list$n_pop)
)

input_pop <- Setup_Mod_Biologicals(
  input_list = input_pop, WAA = sim_data$WAA, MatAA = sim_data$MatAA,
  WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv, fit_lengths = 0,
  AgeingError = sim_data$AgeingError, M_spec = "est_ln_M"
)

input_pop <- Setup_Mod_Tagging(
  input_list = input_pop, use_conv_fish_tagging = 1,
  conv_tagged_fish = sim_data$conv_tagged_fish_attr,
  conv_tag_max_liberty = dim(sim_data$obs_conv_tag_fish_recap)[1],
  obs_conv_tag_fish_recap = sim_data$obs_conv_tag_fish_recap,
  conv_fish_tag_like = "Poisson", init_conv_tag_mort_spec = "fix",
  conv_tag_shed_spec = "fix", conv_tagrep_spec = "est_shared_r",
  conv_fish_tag_attr = "p_a_s", conv_tag_release_indicator = sim_data$conv_tag_release_indicator
)

input_pop <- Setup_Mod_Movement(
  input_list = input_pop, do_recruits_move = 0, use_fixed_movement = 0,
  Movement_popblk_spec = list(1, 2, 3), Movement_seasblk_spec = list(1, 2)
)

input_pop <- Setup_Mod_Catch_and_F(
  input_list = input_pop, ObsCatch = sim_data$ObsCatch,
  UseCatch = array(0, dim = dim(sim_data$UseCatch)),
  ObsCatch_pop = sim_data$ObsCatch_pop, UseCatch_pop = sim_data$UseCatch_pop,
  Use_F_pen = 1, sigmaC_spec = "fix", sigmaC_pop_spec = "fix",
  ln_sigmaC = sim_data$ln_sigmaC, ln_sigmaC_pop = sim_data$ln_sigmaC_pop,
  ln_sigmaF = array(log(1), dim = c(input_pop$data$n_regions, input_pop$data$n_seas, input_pop$data$n_fish_fleets))
)

input_pop <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_pop, ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
  UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
  ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = sim_data$ObsFishLenComps,
  UseFishAgeComps = array(0, dim = dim(sim_data$UseFishAgeComps)), UseFishLenComps = sim_data$UseFishLenComps,
  ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = sim_data$ISS_FishLenComps,
  ObsFishAgeComps_pop = sim_data$ObsFishAgeComps_pop, UseFishAgeComps_pop = sim_data$UseFishAgeComps_pop,
  ISS_FishAgeComps_pop = sim_data$ISS_FishAgeComps_pop,
  FishAgeComps_pop_LikeType = c("Multinomial"), FishAgeComps_pop_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
  fish_idx_type = "none", FishAgeComps_LikeType = c("Multinomial"), FishLenComps_LikeType = c("none"),
  FishAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"), FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
)

input_pop <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_pop, ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
  UseSrvIdx = array(0, dim = dim(sim_data$UseSrvIdx)),
  ObsSrvIdx_pop = sim_data$ObsSrvIdx_pop, ObsSrvIdx_pop_SE = sim_data$ObsSrvIdx_pop_SE,
  UseSrvIdx_pop = array(1, dim = dim(sim_data$UseSrvIdx_pop)),
  ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = sim_data$ObsSrvLenComps,
  UseSrvAgeComps = array(0, dim = dim(sim_data$UseSrvAgeComps)), UseSrvLenComps = sim_data$UseSrvLenComps,
  ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
  ObsSrvAgeComps_pop = sim_data$ObsSrvAgeComps_pop, UseSrvAgeComps_pop = sim_data$UseSrvAgeComps_pop,
  ISS_SrvAgeComps_pop = sim_data$ISS_SrvAgeComps_pop,
  SrvAgeComps_pop_LikeType = c("Multinomial"), SrvAgeComps_pop_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
  srv_idx_type = c("biom"), SrvAgeComps_LikeType = c("Multinomial"), SrvLenComps_LikeType = c("none"),
  SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"), SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
)

input_pop <- Setup_Mod_Fishsel_and_Q(
  input_list = input_pop, fish_sel_model = c("logist1_Fleet_1"),
  fish_fixed_sel_pars_spec = c("est_shared_r"), fish_q_spec = c("fix")
)

input_pop <- Setup_Mod_Srvsel_and_Q(
  input_list = input_pop, srv_sel_model = c("logist1_Fleet_1"),
  srv_fixed_sel_pars_spec = c("est_shared_r"), srv_q_spec = c("est_shared_r")
)

input_pop <- Setup_Mod_Weighting(
  input_list = input_pop, Wt_Catch = 1, Wt_Catch_pop = 1, Wt_SrvIdx_pop = 1,
  Wt_FishIdx = 1, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 1,
  Wt_FishAgeComps = array(0, dim = c(input_pop$data$n_regions, length(input_pop$data$years),
                                     input_pop$data$n_seas, input_pop$data$n_sexes, input_pop$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(0, dim = c(input_pop$data$n_regions, length(input_pop$data$years),
                                    input_pop$data$n_seas, input_pop$data$n_sexes, input_pop$data$n_srv_fleets)),
  Wt_FishAgeComps_pop = array(1, dim = c(input_pop$data$n_pop, input_pop$data$n_regions, length(input_pop$data$years),
                                         input_pop$data$n_seas, input_pop$data$n_sexes, input_pop$data$n_fish_fleets)),
  Wt_SrvAgeComps_pop = array(1, dim = c(input_pop$data$n_pop, input_pop$data$n_regions, length(input_pop$data$years),
                                        input_pop$data$n_seas, input_pop$data$n_sexes, input_pop$data$n_srv_fleets))
)

model_pop <- fit_model(input_pop$data, input_pop$par, input_pop$map,
                       random = NULL, newton_loops = 3, silent = TRUE)

# Population-specific composition OSA (pop = TRUE, joint-sex) ----------

comp_pop <- get_osa(model = model_pop, data = input_pop$data, comp_source = "FishAge",
                    pop = TRUE, family = "discrete", bins = input_pop$data$ages, bin_label = "Ages")
resid_comp_pop <- plot_resids(comp_pop)

png(here("vignettes", "figures", "u_internal_comp_pop.png"), width = 1600, height = 600)
cowplot::plot_grid(resid_comp_pop[[1]], resid_comp_pop[[2]], ncol = 2)
dev.off()

# Conventional tagging OSA (tag = TRUE) ---------------------------------

tag_osa <- get_osa(model = model_pop, data = input_pop$data, tag = TRUE)
resid_tag <- plot_resids(tag_osa)

png(here("vignettes", "figures", "u_internal_tag.png"), width = 1600, height = 600)
cowplot::plot_grid(resid_tag[[1]], resid_tag[[2]], ncol = 2)
dev.off()

# Population-specific index OSA (pop = TRUE) ----------------------------

resid_catch_pop  <- plot_resids(get_osa(model = model_pop, data = input_pop$data, index_source = "Catch", pop = TRUE))
resid_srvidx_pop <- plot_resids(get_osa(model = model_pop, data = input_pop$data, index_source = "SrvIdx", pop = TRUE))

png(here("vignettes", "figures", "u_internal_index_pop.png"), width = 1600, height = 900)
cowplot::plot_grid(
  resid_catch_pop[[1]], resid_catch_pop[[2]],
  resid_srvidx_pop[[1]], resid_srvidx_pop[[2]],
  ncol = 2, nrow = 2
)
dev.off()
