# Purpose: Bridge the 2024 BSAI blackspotted and rougheye rockfish assessment to SPoRC and render its figures
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("tests", "testthat", "helper-bridge_rebs.R"))
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

label <- "2024 BSAI Rougheye Assessment"

dat <- SPoRC::sgl_rg_rebs_data
yrs <- dat$years
n_yrs <- length(yrs)
n_ages <- length(dat$ages)
n_obs_ages <- length(dat$obs_ages)

input_list <- seed_rebs_mle(build_rebs_input(dat), dat)

# Stage 1: evaluate at the assessment's MLE ------------------------------------
obj <- fit_model(input_list$data, input_list$par, input_list$map,
                 do_optim = FALSE, silent = TRUE)
r <- obj$rep

cat("=== Stage 1: at the assessment's MLE, no optimization ===\n")
cat("fishery selectivity max pct diff:",
    100 * max(abs(as.vector(r$fish_sel[1, 1, 1, 1, 1:n_obs_ages, 1, 1]) /
                    as.vector(dat$admb$sel_fsh[1, ]) - 1)), "\n")
cat("survey selectivity  max pct diff:",
    100 * max(abs(as.vector(r$srv_sel[1, 1, 1, 1, 1:n_obs_ages, 1, 1]) /
                    as.vector(dat$admb$sel_srv) - 1)), "\n")

# the assessment reports numbers at age over 43 observed ages with ages 45 to 54 pooled, so SPoRC's
# 52 ages are pooled the same way. the comparison follows the bridge test's windows
n_est <- length(dat$mle$rec_dev)
naa <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
naa_pooled <- cbind(naa[, 1:(n_obs_ages - 1)], rowSums(naa[, n_obs_ages:n_ages]))
cat("N (ages 4+)  max pct diff:",
    100 * max(abs(naa_pooled[, -(1:3)] / dat$admb$NAA[, -(1:3)] - 1)), "\n")
cat("SSB          max pct diff:",
    100 * max(abs(as.vector(r$SSB)[1:n_yrs] / dat$admb$SSB - 1)), "\n")
cat("F            max pct diff:",
    100 * max(abs(as.vector(r$Fmort)[1:n_yrs] / dat$admb$Fmort - 1)), "\n")
cat("terminal recruitment ratio (should be exp(-sigmaR^2/2) =",
    exp(-dat$sigmaR^2 / 2), "):",
    unique(round(as.vector(r$Rec)[(n_est + 1):n_yrs] /
                   dat$admb$Rec[(n_est + 1):n_yrs], 6)), "\n")
cat("jnLL at the assessment MLE:", obj$fn(obj$par), "\n")
cat("max |gradient| there      :", max(abs(obj$gr(obj$par))), "\n")

# Stage 2: optimize -----------------------------------------------------------
est <- fit_model(
  input_list$data,
  input_list$par,
  input_list$map,
  random = NULL,
  newton_loops = 3,
  silent = TRUE
)
est$sdrep <- RTMB::sdreport(est)
cat("\n=== Stage 2: optimized ===\n")
cat("free parameters:", length(est$optim$par), "\n")
cat("final jnLL:", est$optim$objective, "  max |gradient|:",
    max(abs(est$gr(est$optim$par))), "\n")
cat("pdHess:", est$sdrep$pdHess, "\n")

# Stage 3: compare ------------------------------------------------------------
sdr <- est$sdrep
rep <- est$rep

ssb <- as.vector(rep$SSB)[1:n_yrs]
rec <- as.vector(rep$Rec)[1:n_yrs]

# the two series overplot, so the difference gets its own panel. the three terminal recruitments
# sit exactly exp(-sigmaR^2/2) low by the documented convention difference
ggplot2::ggsave(
  here("vignettes", "figures", "z_rebs_ts_comparison.png"),
  bridge_ts_figure(
    yrs,
    ssb,
    rec,
    dat$admb$SSB,
    dat$admb$Rec,
    label,
    ssb_se = bridge_se(sdr, "log_SSB", ssb),
    rec_se = bridge_se(sdr, "log_Rec", rec),
    mark_year = yrs[n_est] + 0.5
  ),
  width = 17,
  height = 9,
  dpi = 150
)

# Selectivity is time invariant, so a single curve per fleet suffices.
sel_df <- bind_rows(
  bridge_sel_rows(dat$obs_ages, rep$fish_sel[1, 1, 1, 1, 1:n_obs_ages, 1, 1],
                  dat$admb$sel_fsh[1, ], "Fishery", label),
  bridge_sel_rows(dat$obs_ages, rep$srv_sel[1, 1, 1, 1, 1:n_obs_ages, 1, 1],
                  dat$admb$sel_srv, "AI Survey", label)
)

ggplot2::ggsave(
  here("vignettes", "figures", "z_rebs_sel_comparison.png"),
  bridge_sel_figure(sel_df, base_size = 15),
  width = 12,
  height = 5,
  dpi = 150
)

cat("\n=== Stage 3: optimized SPoRC against the assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, dat$admb$SSB),
            bridge_cmp("Recruitment (estimated years)", rec[1:n_est], dat$admb$Rec[1:n_est]),
            bridge_cmp("Recruitment (terminal 3)", rec[(n_est + 1):n_yrs],
                       dat$admb$Rec[(n_est + 1):n_yrs])),
      row.names = FALSE, digits = 4)
