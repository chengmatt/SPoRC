# Purpose: Set up the 2024 Alaska sablefish assessment in SPoRC and render the
#          single region case study figures. Everything is specified the way the
#          sablefish assessment specifies it: mean recruitment with a Methot and
#          Taylor bias ramp, early/late sigmaR, a historical F that is a fixed
#          proportion of the fixed-gear mean F, blocked logistic and gamma
#          selectivity, and the assessment's own Francis weights.
#
#          This is a CLOSE bridge, not an exact one. tem.tpl carries several
#          conventions SPoRC deliberately does not reproduce (sum-to-zero
#          deviation vectors, a recruitment penalty centred at zero, an
#          uncorrected terminal recruitment, and a survey index reweighted by an
#          observed sex ratio). The measured gap is reported at the bottom and is
#          what the vignette's difference table quotes.
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26

library(here)
library(dplyr)
library(ggplot2)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

dat <- SPoRC::sgl_rg_sable_data
yrs <- 1960:2024
n_yrs <- length(dat$years)
n_ages <- length(dat$ages)

# Setup ----------------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  years = 1:length(dat$years),
  ages = 1:length(dat$ages),
  lens = seq(41, 99, 2),
  n_regions = 1,
  n_sexes = dat$n_sexes,
  n_fish_fleets = dat$n_fish_fleets,
  n_srv_fleets = dat$n_srv_fleets,
  n_pop = dat$n_pop,
  verbose = FALSE
)

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  rec_model = "mean_rec",
  do_rec_bias_ramp = 1,
  # bias_year is in DEV-INDEX space, not calendar years. Passing calendar years
  # leaves every range empty and silently gives bias_ramp = 0 throughout.
  bias_year = c(length(1960:1979),
                length(1960:1989),
                (length(1960:2023) - 5),
                length(1960:2024) - 2) + 1,
  sigmaR_switch = as.integer(length(1960:1975)),
  ln_sigmaR = array(log(c(0.4, 1.2)), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
  sigmaR_spec = "fix_early_est_late",
  dont_est_recdev_last = 1,
  init_age_strc = 1,
  # tem.tpl builds hist_hal_F as hist_hal_prop * exp(log_avg_F_fish1), a fixed
  # proportion of the fixed-gear mean F, which is exactly init_F_form = "prop".
  # Assessments carrying an INDEPENDENT historical F need "abs" instead.
  init_F_form = "prop",
  init_F_spec = "fix",
  init_F_par = array(stats::qlogis(0.1), dim = c(input_list$data$n_regions,
                                                 input_list$data$n_seas,
                                                 input_list$data$n_fish_fleets))
)

fixed_natmort <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                  n_yrs, n_ages, input_list$data$n_sexes))
fixed_natmort[, , , , 1] <- 0.1134156 # female M
fixed_natmort[, , , , 2] <- 0.1052175 # male M

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA = dat$WAA,
  MatAA = dat$MatAA,
  AgeingError = as.matrix(dat$age_error),
  SizeAgeTrans = dat$SizeAgeTrans,
  Use_M_prior = 0,
  fit_lengths = 1,
  # tem.tpl estimates a single logm with a male offset mdelta. SPoRC has no
  # offset parameterisation, so both are fixed at the 2024 estimates rather than
  # left as an uncontrolled difference.
  M_spec = "fix",
  Fixed_natmort = fixed_natmort
)

input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = NA,
  do_recruits_move = 0
)

input_list <- Setup_Mod_Tagging(
  input_list = input_list,
  use_conv_fish_tagging = rep(0, input_list$data$n_fish_fleets)
)

# tem.tpl writes catch and the F penalty as unweighted sums of squares
# (wt_ssqcatch * norm2(...) and wt_fmort_reg * norm2(...)). sigma = 1/sqrt(2)
# makes 1/(2 sigma^2) = 1, so SPoRC's Gaussian collapses to the same sum of
# squares and the assessment's lambdas transfer over unchanged.
input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,
  ObsCatch = dat$ObsCatch,
  UseCatch = dat$UseCatch,
  Use_F_pen = 1,
  sigmaC_spec = 'fix',
  sigmaF_spec = "fix",
  ln_sigmaC = array(log(sqrt(1 / 2)), dim = c(input_list$data$n_regions, n_yrs,
                                              input_list$data$n_seas,
                                              input_list$data$n_fish_fleets)),
  ln_sigmaF = array(log(sqrt(1 / 2)), dim = c(input_list$data$n_regions,
                                              input_list$data$n_seas,
                                              input_list$data$n_fish_fleets))
)

# tem.tpl divides the index residual by the CV (obs_se/obs_biom), not the SE.
# Passing raw SEs weights the index by a factor of obs too many.
input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  ObsFishIdx = dat$ObsFishIdx,
  ObsFishIdx_SE = dat$ObsFishIdx_SE / dat$ObsFishIdx,
  UseFishIdx = dat$UseFishIdx,
  ObsFishAgeComps = dat$ObsFishAgeComps,
  UseFishAgeComps = dat$UseFishAgeComps,
  ISS_FishAgeComps = dat$ISS_FishAgeComps,
  ObsFishLenComps = dat$ObsFishLenComps,
  UseFishLenComps = dat$UseFishLenComps,
  ISS_FishLenComps = dat$ISS_FishLenComps,
  fish_idx_type = c("biom", "none"),
  FishAgeComps_LikeType = c("Multinomial", "none"),
  FishLenComps_LikeType = c("Multinomial", "Multinomial"),
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1",
                        "none_Year_1-terminal_Fleet_2"),
  FishLenComps_Type = c("spltRspltS_Year_1-terminal_Fleet_1",
                        "spltRspltS_Year_1-terminal_Fleet_2")
)

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,
  ObsSrvIdx = dat$ObsSrvIdx,
  ObsSrvIdx_SE = dat$ObsSrvIdx_SE / dat$ObsSrvIdx,
  UseSrvIdx = dat$UseSrvIdx,
  ObsSrvAgeComps = dat$ObsSrvAgeComps,
  ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
  UseSrvAgeComps = dat$UseSrvAgeComps,
  ObsSrvLenComps = dat$ObsSrvLenComps,
  UseSrvLenComps = dat$UseSrvLenComps,
  ISS_SrvLenComps = dat$ISS_SrvLenComps,
  srv_idx_type = c("abd", "biom", "abd"),
  SrvAgeComps_LikeType = c("Multinomial", "none", "Multinomial"),
  SrvLenComps_LikeType = c("Multinomial", "Multinomial", "Multinomial"),
  SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1",
                       "none_Year_1-terminal_Fleet_2",
                       "agg_Year_1-terminal_Fleet_3"),
  SrvLenComps_Type = c("spltRspltS_Year_1-terminal_Fleet_1",
                       "spltRspltS_Year_1-terminal_Fleet_2",
                       "spltRspltS_Year_1-terminal_Fleet_3")
)

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list = input_list,
  cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
  # pre-IFQ, IFQ, and recent blocks for the fixed-gear fleet
  fish_sel_blocks = c("Block_1_Year_1-35_Fleet_1",
                      "Block_2_Year_36-56_Fleet_1",
                      "Block_3_Year_57-terminal_Fleet_1",
                      "none_Fleet_2"),
  fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),
  fish_q_blocks = c("Block_1_Year_1-35_Fleet_1",
                    "Block_2_Year_36-56_Fleet_1",
                    "Block_3_Year_57-terminal_Fleet_1",
                    "none_Fleet_2"),
  fish_fixed_sel_pars_spec = c("est_all", "est_all"),
  fish_q_spec = c("est_all", "fix")
)

# share the slope across sexes in the early fixed-gear block
input_list$map$fish_fixed_sel_pars <- factor(c(1:7, 2, 8:11, rep(12:13, 3), rep(c(14, 13), 3)))

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_sel_blocks = c("Block_1_Year_1-56_Fleet_1",
                     "Block_2_Year_57-terminal_Fleet_1",
                     "none_Fleet_2",
                     "none_Fleet_3"),
  srv_sel_model = c("logist1_Fleet_1", "exponential_Fleet_2", "logist1_Fleet_3"),
  srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all"),
  srv_q_spec = c("est_all", "est_all", "est_all"),
  # tem.tpl predicts survey indices from mid-year survival (S_f_mid / S_m_mid)
  t_srv = array(0.5, dim = c(input_list$data$n_regions,
                             input_list$data$n_seas,
                             input_list$data$n_srv_fleets))
)

# Longline survey slopes (indices 2 and 5) are shared across time blocks and with
# the cooperative Japanese survey, which estimates nothing of its own. The trawl
# survey's power function carries a single parameter per sex (indices 7, 8).
input_list$map$srv_fixed_sel_pars <-
  factor(c(1:3, 2, 4:6, 5, rep(7, 4),
           rep(8, 4), rep(c(NA, 2), 2), rep(c(NA, 5), 2)))

# Cooperative Japanese survey ended in 1994 and its logistic is held at the 2024
# values rather than estimated from the remaining data.
input_list$par$srv_fixed_sel_pars[1, , , 1, 3] <- c(0.980660760456, 0.9295241)
input_list$par$srv_fixed_sel_pars[1, , , 2, 3] <- c(1.22224502478, 0.8821623)

# Weighting ------------------------------------------------------------------
# Converged Francis weights from the 2024 assessment.
dim_comps <- function(n_fleets) {
  array(NA, dim = c(input_list$data$n_regions, n_yrs, input_list$data$n_seas,
                    input_list$data$n_sexes, n_fleets))
}

Wt_FishAgeComps <- dim_comps(input_list$data$n_fish_fleets)
Wt_FishAgeComps[1, , , 1, 1] <- 0.826107286513784

Wt_FishLenComps <- dim_comps(input_list$data$n_fish_fleets)
Wt_FishLenComps[1, , , 1, 1] <- 4.1837057381917
Wt_FishLenComps[1, , , 2, 1] <- 4.26969350917589
Wt_FishLenComps[1, , , 1, 2] <- 0.316485920691651
Wt_FishLenComps[1, , , 2, 2] <- 0.229396580680981

Wt_SrvAgeComps <- dim_comps(input_list$data$n_srv_fleets)
Wt_SrvAgeComps[1, , , 1, 1] <- 3.79224544725927
Wt_SrvAgeComps[1, , , 1, 3] <- 1.31681114024037

Wt_SrvLenComps <- dim_comps(input_list$data$n_srv_fleets)
Wt_SrvLenComps[1, , , 1, 1] <- 1.43792019016567
Wt_SrvLenComps[1, , , 2, 1] <- 1.07053763450712
Wt_SrvLenComps[1, , , 1, 2] <- 0.670883273592302
Wt_SrvLenComps[1, , , 2, 2] <- 0.465207132450763
Wt_SrvLenComps[1, , , 1, 3] <- 1.27772810174693
Wt_SrvLenComps[1, , , 2, 3] <- 0.857519546948587

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 50,
  Wt_FishIdx = 0.448,
  Wt_SrvIdx = 0.448,
  Wt_Rec = 1.5,
  Wt_F = 0.1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = Wt_FishAgeComps,
  Wt_FishLenComps = Wt_FishLenComps,
  Wt_SrvAgeComps = Wt_SrvAgeComps,
  Wt_SrvLenComps = Wt_SrvLenComps
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Fit ------------------------------------------------------------------------
est <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 3, silent = TRUE)
est$sd_rep <- RTMB::sdreport(est)

cat("=== Fit ===\n")
cat("free parameters:", length(est$par), "\n")
cat("final jnLL:", est$optim$objective, "  max |gradient|:", max(abs(est$gr(est$optim$par))), "\n")

# Compare --------------------------------------------------------------------
sdr <- est$sd_rep
rep <- est$rep

ssb <- as.vector(rep$SSB[1, 1, 1:n_yrs])
recr <- as.vector(rep$Rec[1, 1, 1:n_yrs])

ggplot2::ggsave(here("vignettes", "figures", "e_ts_comparison.png"),
                bridge_ts_figure(yrs = yrs, ssb = ssb, rec = recr,
                                 admb_ssb = dat$admb_spbiom, admb_rec = dat$admb_recr,
                                 label = "2024 Sablefish Assessment",
                                 ssb_se = bridge_se(sdr, "log_SSB", ssb, exact = TRUE),
                                 rec_se = bridge_se(sdr, "log_Rec", recr, exact = TRUE),
                                 legend_nrow = NULL),
                width = 17, height = 9, dpi = 150)

cat("\n=== SPoRC against the 2024 sablefish assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, dat$admb_spbiom, signed = TRUE),
            bridge_cmp("Recruitment", recr, dat$admb_recr, signed = TRUE)),
      row.names = FALSE, digits = 4)

# The terminal year is the largest single discrepancy by construction: tem.tpl
# builds terminal recruitment as exp(log_mean_rec) with NO bias correction, while
# SPoRC applies the ramp uniformly. Reported separately so it is not read as a
# population dynamics error.
cat("\nterminal-year recruitment: SPoRC", recr[n_yrs],
    " ADMB", dat$admb_recr[n_yrs],
    " ratio", recr[n_yrs] / dat$admb_recr[n_yrs], "\n")
cat("expected ratio from the missing bias correction, exp(sigmaR^2/2):",
    exp(1.2^2 / 2), "\n")

cat("\nexcluding the terminal year:\n")
print(rbind(bridge_cmp("SSB", ssb[-n_yrs], dat$admb_spbiom[-n_yrs], signed = TRUE),
            bridge_cmp("Recruitment", recr[-n_yrs], dat$admb_recr[-n_yrs], signed = TRUE)),
      row.names = FALSE, digits = 4)

saveRDS(list(rep = rep, sdrep = est$sd_rep, opt = est$optim,
             data = data, parameters = parameters, mapping = mapping),
        here("dev", "dev_output", "sgl_rg_sablefish_bridge.rds"))
