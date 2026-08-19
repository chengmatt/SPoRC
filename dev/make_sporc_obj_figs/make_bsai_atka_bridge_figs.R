# Purpose: Bridge the 2024 BSAI Atka mackerel assessment (AMAK Model 16.0b) to
#          SPoRC and render the case study figures. The model is specified the
#          way AMAK specifies it, including AMAK's own 1967 start year, so the
#          period AMAK calls pre-model is carried as model years under a true
#          closure rather than as an initial age structure with deviations.
#          Recruitment is a mean with deviations and the Beverton-Holt curve is
#          fitted as a penalty on the residual, taking its scale from the
#          initial equilibrium so one parameter plays AMAK's log_Rzero in both
#          places it appears.
#
#          Three stages: set every parameter to the AMAK maximum likelihood
#          estimate and check the objective there, optimize, then compare.
# Creator: Matthew LH. Cheng
# Date Created: 8/18/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

dat <- SPoRC::sgl_rg_bsai_atka_data
n_pre <- 10                                        # 1967 to 1976
yrs <- (min(dat$years) - n_pre):max(dat$years)     # AMAK's own model years
n_yrs <- length(yrs)
n_ages <- length(dat$ages)
n_obs <- length(dat$years)
i_amak <- which(yrs %in% dat$years)                # the years AMAK reports on

i_srv <- which(yrs %in% dat$yrs_srv)
inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

# The assessment's arrays run 1977 to 2024. Put ten years in front of whichever
# margin holds the years, so they run 1967 to 2024.
pad <- function(x, fill = NULL) {
  m <- which(dim(x) == n_obs)
  pre <- if(is.null(fill)) abind::asub(x, rep(1, n_pre), m, drop = FALSE)
         else array(fill, dim = replace(dim(x), m, n_pre))
  unname(abind::abind(pre, x, along = m))
}
# Setup ----------------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  years = yrs,
  ages = dat$ages,
  lens = NA,
  n_regions = dat$n_regions,
  n_sexes = dat$n_sexes,
  n_fish_fleets = dat$n_fish_fleets,
  n_srv_fleets = dat$n_srv_fleets,
  n_seas = dat$n_seas,
  n_pop = dat$n_pop,
  natal_region = dat$natal_region,
  verbose = FALSE
)

# AMAK's precision on each recruitment deviation is 1 everywhere from
# rec_like(2), plus 0.5/sigmaR^2 on the tails from rec_like(4). SPoRC charges
# Wt_Rec * eps^2 / (2 * sigma^2), so with sigma at AMAK's sigmaR the weights are
# 2*sigmaR^2 and 2*sigmaR^2 + 1. Both tails sit inside one deviation vector.
tail_yrs <- c(dat$styr_rec:dat$styr_rec_est, dat$endyr_rec_est:max(yrs))

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  # Recruitment is a mean with deviations; the curve is scored, not used.
  rec_model = "mean_rec",
  rec_lag = 1,
  SR_ref_yr = 1,
  sr_penalty = "bh",
  # rec_like(1), with AMAK's m_sigmar^2/2 inflation folded into the sigma
  sr_pen_sigma = dat$sigmaR / sqrt(1 + 1 / (2 * dat$nrecs_est)),
  # AMAK's window in full. The pre-model years supply 1977's lagged spawning
  # biomass, which is what a 1977 start could not do.
  sr_pen_yrs = dat$styr_rec_est:dat$endyr_rec_est,
  # One parameter plays AMAK's log_Rzero in both roles, the 1967 equilibrium
  # and the curve's scale.
  sr_R0_spec = "rinit",
  # Steepness is fixed at 0.8 and carried on a logit bounded to (0.2, 1), so the
  # fixed value goes in transformed. Passing 0.8 directly is swallowed by the
  # starting value mechanism and leaves the default h = 0.6.
  steepness_h = array(inv_steepness(dat$steepness), dim = c(1, 1)),
  h_spec = "fix",
  # do_rec_bias_ramp = 0 does NOT centre the penalty on zero: it sets the ramp
  # to one throughout, which centres on -sigmaR^2/2. Breaks past the last year
  # are what leave the offset out.
  do_rec_bias_ramp = 1,
  bias_year = rep(n_yrs + 1, 4),
  sigmaR_switch = 1,
  sigmaR_spec = "fix",
  ln_sigmaR = array(rep(log(dat$sigmaR), 2), dim = c(2, 1, 1)),
  # AMAK's 1967 age structure is the Rzero equilibrium with nothing on it.
  InitDevs_spec = "fix",
  init_age_strc = 1,
  equil_init_age_strc = 2,
  ln_global_R0 = dat$mle$mean_log_rec,   # under mean_rec this IS mean recruitment
  ln_rinit = dat$mle$log_Rzero,
  use_rinit = 1,
  t_spawn = dat$t_spawn,
  # AMAK estimates every deviation from 1967 to 2024 and restricts only its
  # stock recruit LIKELIHOOD; dont_est_recdev_last would delete them instead.
  dont_est_recdev_last = 0,
  # rec_like(2) is carried by Rec_nLL directly.
  Use_rec_level_pen = 0
)

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA = pad(dat$WAA),
  WAA_fish = pad(dat$WAA_fish),
  WAA_srv = pad(dat$WAA_srv),
  MatAA = pad(dat$MatAA),
  AgeingError = dat$AgeingError,
  fit_lengths = 0,
  M_spec = "fix",
  Fixed_natmort = pad(dat$Fixed_natmort),
  addtocomp = 1e-3,
  comp_const_obs = 1,
  addtosrvidx = 0,
  addtofishidx = 0
)

input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                 Fixed_Movement = NA, do_recruits_move = 0)

# ObsCatch = 0 with UseCatch = 0 is a true closure: F is forced to zero and the
# deviation leaves the parameter vector. NA would be a missing observation
# instead, which still estimates an F.
input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,
  ObsCatch = pad(dat$ObsCatch, 0),
  UseCatch = pad(dat$UseCatch, 0),
  # AMAK has no log_avg_F: fmort is the annual rate itself, so fishing mortality
  # is free annual log F. ln_F_mean_spec = "fix" is that parameterization, and it
  # pairs with Use_F_pen = 0 here because AMAK's only F statement penalizes the
  # realized F at age surface toward 0.2 rather than the deviations, which no
  # SPoRC term expresses.
  ln_F_mean_spec = "fix",
  Use_F_pen = 0,
  sigmaC_spec = "fix",
  ln_sigmaC = array(log(dat$sigmaC), dim = c(1, n_yrs, 1, 1))
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  ObsFishIdx = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
  ObsFishIdx_SE = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
  UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
  ObsFishAgeComps = pad(dat$ObsFishAgeComps, 0),
  UseFishAgeComps = pad(dat$UseFishAgeComps, 0),
  ISS_FishAgeComps = pad(dat$ISS_FishAgeComps, 0),
  ObsFishLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), 1, 1)),
  UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
  ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
  fish_idx_type = c("biom"),
  FishIdx_LikeType = c("lognormal"),
  FishAgeComps_LikeType = c("Multinomial"),
  FishLenComps_LikeType = c("none"),
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
  FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
)

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,
  ObsSrvIdx = pad(dat$ObsSrvIdx, NA_real_),
  ObsSrvIdx_SE = pad(dat$ObsSrvIdx_SE, NA_real_),
  UseSrvIdx = pad(dat$UseSrvIdx, 0),
  ObsSrvAgeComps = pad(dat$ObsSrvAgeComps, 0),
  UseSrvAgeComps = pad(dat$UseSrvAgeComps, 0),
  ISS_SrvAgeComps = pad(dat$ISS_SrvAgeComps, 0),
  ObsSrvLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), 1, 1)),
  UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
  ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
  srv_idx_type = c("biom"),
  SrvIdx_LikeType = c("lognormal"),
  SrvAgeComps_LikeType = c("Multinomial"),
  SrvLenComps_LikeType = c("none"),
  SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
  SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1"),
  t_srv = array(dat$t_srv, dim = c(1, 1, 1))
)

# Selectivity ----------------------------------------------------------------
# The oldest estimated coefficient is held flat to the plus group before the
# standardization, which is a bin grouping rather than a separate parameter.
bin_groups <- function(nsel) c(as.list(seq_len(nsel - 1)), list(nsel:n_ages))

# The closed years ride on AMAK's first selectivity block, which never touches
# the dynamics because F is zero there.
blk <- c(paste0("Block_1_Year_1-", n_pre + 1, "_Fleet_1"),
         paste0("Block_", 2:(dat$n_blk_fsh - 1), "_Year_", (n_pre + 2):(n_yrs - 2), "-",
                (n_pre + 2):(n_yrs - 2), "_Fleet_1"),
         paste0("Block_", dat$n_blk_fsh, "_Year_", n_yrs - 1, "-terminal_Fleet_1"))
stopifnot(length(blk) == dat$n_blk_fsh)

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list = input_list,
  cont_tv_fish_sel = c("none_Fleet_1"),
  fish_sel_blocks = blk,
  fish_sel_model = c("nonparlog_Fleet_1"),
  fish_sel_nonpar_est_bins = list(rep(list(bin_groups(dat$nselages_fsh)), dat$n_blk_fsh)),
  fish_fixed_sel_pars_spec = c("est_all"),
  fish_q_blocks = c("none_Fleet_1"),
  fish_q_spec = c("fix"),
  use_fixed_fish_sel = 0,
  Use_fish_selex_penalty = 1,
  fish_selex_penalty = dat$fish_selex_penalty
)

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  cont_tv_srv_sel = c("none_Fleet_1"),
  srv_sel_blocks = c("none_Fleet_1"),
  srv_sel_model = c("nonparlog_Fleet_1"),
  srv_sel_nonpar_est_bins = list(list(bin_groups(dat$nselages_srv))),
  # AMAK standardizes the index selectivity over the ages catchability is
  # defined against rather than over every age (amak.tpl:1947). The difference
  # is a constant that catchability would absorb were the q prior uninformative;
  # with a prior it is not free, so the window has to match.
  srv_sel_norm_bins = list(dat$q_age_min:dat$q_age_max),
  srv_fixed_sel_pars_spec = c("est_all"),
  srv_q_blocks = c("none_Fleet_1"),
  srv_q_spec = c("est_all"),
  srv_q_type = c("est"),
  use_fixed_srv_sel = 0,
  Use_srv_selex_penalty = 1,
  srv_selex_penalty = dat$srv_selex_penalty,
  Use_srv_q_prior = 1,
  srv_q_prior = dat$srv_q_prior,
  t_srv = array(dat$t_srv, dim = c(1, 1, 1))
)

input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

# The smoothness weights are per year, so the closed years carry none: block 1
# would otherwise be charged its curvature and dome penalty eleven times over.
fw <- dat$fish_sel_pen_wts
for(nm in c("smooth_bin_curve", "smooth_yr_diff", "smooth_dome")) {
  fw[[1]][[nm]] <- c(rep(0, n_pre), dat$fish_sel_pen_wts[[1]][[nm]])
} # end nm loop
fw[[1]]$smooth_yr_diff[n_pre + 1] <- 0

sw <- dat$srv_sel_pen_wts
for(nm in c("smooth_bin_curve", "smooth_yr_diff", "smooth_dome")) {
  if(length(sw[[1]][[nm]]) == n_obs) sw[[1]][[nm]] <- c(rep(0, n_pre), dat$srv_sel_pen_wts[[1]][[nm]])
} # end nm loop

Wt_Rec <- array(2 * dat$sigmaR^2, dim = c(1, 1, n_yrs))
Wt_Rec[1, 1, which(yrs %in% tail_yrs)] <- 2 * dat$sigmaR^2 + 1

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1,
  Wt_Rec = Wt_Rec, Wt_Init_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
  Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  fish_sel_pen_wts = fw,
  srv_sel_pen_wts = sw
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Seed every parameter at the AMAK maximum likelihood estimate ----------------
parameters$ln_F_devs[1, n_pre + seq_len(n_obs), 1, 1] <- log(dat$mle$fmort)
for(b in seq_len(dat$n_blk_fsh)) {
  parameters$fish_fixed_sel_pars[1, , b, 1, 1] <-
    c(dat$mle$log_selcoffs_fsh[b, ], dat$mle$log_selcoffs_fsh[b, dat$nselages_fsh])
} # end b loop
parameters$srv_fixed_sel_pars[1, , 1, 1, 1] <-
  c(dat$mle$log_selcoffs_ind, dat$mle$log_selcoffs_ind[dat$nselages_srv])
# The deviations are AMAK's rec_dev outright, over AMAK's own year range.
parameters$ln_RecDevs[1, 1, ] <- dat$amak$rec_dev[as.character(yrs)]

# Catchability is whatever puts the predicted index on AMAK's, so it is added to
# the CURRENT value rather than set to AMAK's.
pass1 <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
parameters$ln_srv_q[1, 1, 1] <- parameters$ln_srv_q[1, 1, 1] +
  log(mean(dat$amak$pred_ind / pass1$rep$PredSrvIdx[1, 1, i_srv, 1, 1]))

seed <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
r0 <- seed$rep

pd <- function(a, b) 100 * max(abs(as.vector(a) - as.vector(b)) / abs(as.vector(b)))
cat("\n=== Stage 1: at the AMAK MLE ===\n")
cat("numbers at age            max pct diff:", pd(r0$NAA[1, 1, i_amak, 1, , 1], dat$amak$NAA), "\n")
cat("fishing mortality at age  max pct diff:", pd(r0$tot_FAA[1, 1, i_amak, 1, , 1, 1], dat$amak$FAA), "\n")
cat("catch at age              max pct diff:", pd(r0$CAA[1, 1, i_amak, 1, , 1, 1], dat$amak$CAA), "\n")
cat("spawning biomass          max pct diff:", pd(r0$SSB[1, 1, i_amak], dat$amak$SSB), "\n")
cat("recruitment               max pct diff:", pd(r0$Rec[1, 1, i_amak], dat$amak$Rec), "\n")
cat("predicted catch           max pct diff:", pd(r0$PredCatch[1, 1, i_amak, 1, 1], dat$amak$pred_catch), "\n")
cat("predicted survey index    max pct diff:", pd(r0$PredSrvIdx[1, 1, i_srv, 1, 1], dat$amak$pred_ind), "\n")
cat("fishery selectivity       max pct diff:", pd(r0$fish_sel[1, 1, i_amak, 1, , 1, 1], dat$amak$sel_fsh), "\n")
cat("survey selectivity        max pct diff:",
    pd(r0$srv_sel[1, 1, 1, 1, , 1, 1] / r0$srv_sel[1, 1, 1, 1, 1, 1, 1],
       dat$amak$sel_ind[1, ] / dat$amak$sel_ind[1, 1]), "\n")
# The 1967 equilibrium is AMAK's Get_Bzero, so its unfished spawning biomass is
# comparable now that the curve's scale is the same parameter.
cat("unfished spawning biomass:", format(r0$SSB[1, 1, 1], digits = 10),
    "vs AMAK", format(dat$amak$Bzero, digits = 10), "\n")
# Under mean_rec the deviations are AMAK's rec_dev, not its stock recruit
# residual, so that is what to check.
cat("recruitment deviations against AMAK's rec_dev, max abs:",
    max(abs(r0$ln_RecDevs[1, 1, ] - dat$amak$rec_dev[as.character(yrs)])), "\n")
sig_sr <- dat$sigmaR / sqrt(1 + 1 / (2 * dat$nrecs_est))
n_pen <- sum(data$sr_pen_yrs)
cat("stock recruit penalty, quadratic part:",
    format(sum(r0$SR_pen_nLL) - n_pen * (0.5 * log(2 * pi) + log(sig_sr)), digits = 10),
    "vs AMAK", format(dat$amak$rec_parts[2] - dat$nrecs_est * log(dat$sigmaR), digits = 10), "\n")

# Likelihood crosswalk. SPoRC writes each component as a proper density while
# AMAK drops normalising constants, so a like for like comparison subtracts
# exactly the constants AMAK omits.
lg2pi <- log(2 * pi)
const <- c(catch = n_obs * (0.5 * lg2pi + log(dat$sigmaC)),
           fish_age = 0,
           srv_idx = sum(0.5 * lg2pi + log(dat$ObsSrvIdx_SE[1, which(dat$years %in% dat$yrs_srv), 1, 1])),
           srv_age = 0,
           sel = 0,
           srv_q = 0.5 * lg2pi + log(dat$srv_q_prior$sd))
sel_amak <- sum(dat$amak$sel_fsh_parts) + sum(dat$amak$sel_ind_parts) +
  20 * sum(dat$amak$avgsel_fsh^2) + 20 * dat$amak$avgsel_ind^2

lik <- data.frame(
  component = c("Catch", "Fishery age compositions", "Survey index",
                "Survey age compositions", "Selectivity penalties",
                "Survey catchability prior"),
  sporc = c(sum(r0$Catch_nLL), sum(r0$FishAgeComps_nLL), sum(r0$SrvIdx_nLL),
            sum(r0$SrvAgeComps_nLL), r0$sel_nLL, r0$srv_q_nLL),
  amak = c(dat$amak$Like_Comp[["catch_like"]], dat$amak$Like_Comp[["age_like_fsh"]],
           dat$amak$Like_Comp[["ind_like"]], dat$amak$Like_Comp[["age_like_ind"]],
           sel_amak, dat$amak$Like_Comp[["post_priors_indq"]])
)
lik$sporc_less_constant <- lik$sporc - as.numeric(const)
lik$difference <- lik$sporc_less_constant - lik$amak
cat("\n=== Likelihood, SPoRC less its density constants against AMAK ===\n")
print(lik[, c("component", "sporc_less_constant", "amak", "difference")],
      row.names = FALSE, digits = 9)

# Stage 2: optimize -----------------------------------------------------------
est <- fit_model(data, parameters, mapping, do_optim = TRUE, newton_loops = 3, silent = TRUE)
cat("\n=== Stage 2: optimized ===\n")
cat("free parameters:", length(est$par), "\n")
cat("final jnLL:", est$optim$objective, "  max |gradient|:",
    max(abs(est$gr(est$optim$par))), "\n")
cat("jnLL at the AMAK MLE:", seed$rep$jnLL,
    "  improvement available:", seed$rep$jnLL - est$rep$jnLL, "\n")
g <- as.vector(seed$gr(seed$par)); blk <- names(seed$par)
gd <- data.frame(block = unique(blk),
                 max_abs_grad = sapply(unique(blk), function(k) max(abs(g[blk == k]))))
cat("\ngradient at the AMAK MLE, by parameter block:\n")
print(gd[order(-gd$max_abs_grad), ], row.names = FALSE, digits = 4)
est$sdrep <- RTMB::sdreport(est, hessian.fixed = est$he())

# Stage 3: compare ------------------------------------------------------------
sdr <- est$sdrep
rep <- est$rep
ssb <- as.vector(rep$SSB[1, 1, i_amak])
rec <- as.vector(rep$Rec[1, 1, i_amak])

ggplot2::ggsave(here("vignettes", "figures", "ab_bsai_atka_ts_comparison.png"),
                bridge_ts_figure(yrs = dat$years, ssb = ssb, rec = rec,
                                 admb_ssb = dat$amak$SSB, admb_rec = dat$amak$Rec,
                                 label = "2024 Atka Mackerel Assessment",
                                 ssb_se = bridge_se(sdr, "log_SSB", rep$SSB[1, 1, ], exact = TRUE)[i_amak],
                                 rec_se = bridge_se(sdr, "log_Rec", rep$Rec[1, 1, ], exact = TRUE)[i_amak],
                                 legend_nrow = NULL),
                width = 17, height = 9, dpi = 150)

# Fishery selectivity changes in all but one year, so a handful of years carries
# the shape change; the full surface is what the bridge test checks.
sel_yrs <- c(1977, 1990, 2005, 2024)
sel_df <- bind_rows(lapply(sel_yrs, function(y) {
  bridge_sel_rows(ages = dat$ages,
                  sporc = rep$fish_sel[1, 1, which(yrs == y), 1, , 1, 1],
                  admb = dat$amak$sel_fsh[which(dat$years == y), ],
                  gear = paste("Fishery", y),
                  label = "2024 Atka Mackerel Assessment", facet_by = "Gear")
})) %>%
  bind_rows(bridge_sel_rows(
    ages = dat$ages,
    # Both are standardized over AMAK's catchability ages now, so this only
    # puts AMAK's reported curve on the same footing.
    sporc = rep$srv_sel[1, 1, 1, 1, , 1, 1],
    admb = dat$amak$sel_ind[1, ] / mean(dat$amak$sel_ind[1, dat$q_age_min:dat$q_age_max]),
    gear = "Survey", label = "2024 Atka Mackerel Assessment", facet_by = "Gear"))

ggplot2::ggsave(here("vignettes", "figures", "ab_bsai_atka_sel_comparison.png"),
                bridge_sel_figure(sel_df, facet_by = "Gear", nrow = 1),
                width = 17, height = 6, dpi = 150)

cat("\n=== Stage 3: optimized SPoRC against the AMAK assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, dat$amak$SSB),
            bridge_cmp("Recruitment", rec, dat$amak$Rec)),
      row.names = FALSE, digits = 4)

saveRDS(list(rep = rep, sdrep = est$sdrep, opt = est$optim, seed_rep = r0, lik = lik,
             data = data, parameters = parameters, mapping = mapping, yrs = yrs, i_amak = i_amak),
        here("dev", "dev_output", "bsai_atka_bridge.rds"))
