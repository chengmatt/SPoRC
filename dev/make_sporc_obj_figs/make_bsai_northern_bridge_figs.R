# Purpose: Bridge the 2023 BSAI northern rockfish assessment to SPoRC and render
#          the case study figures. The configuration is not restated here: it is
#          sourced from the bridge test helper, so the figures, the bridge test,
#          and the pinned regression test all run the same specification and a
#          change to one cannot silently leave the others behind.
#
#          Three stages: set every parameter to the assessment's maximum
#          likelihood estimate and check the population and likelihood there,
#          optimize, then compare. The bridge stage supplies the assessment's
#          age 30 edge hold on survey selectivity as a fixed input, which
#          SPoRC's logist1 form cannot express; the refit estimates the uncapped
#          logistic instead. That is the only specification difference.
# Creator: Matthew LH. Cheng
# Date Created: 8/8/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("tests", "testthat", "helper-bridge_bsai_nork.R"))
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

label <- "2023 BSAI Northern Rockfish Assessment"

dat <- SPoRC::sgl_rg_bsai_nork_data
yrs <- dat$years
n_yrs <- length(yrs)
n_ages <- length(dat$ages)
n_obs_ages <- length(dat$obs_ages)

input_list <- seed_bsai_nork_mle(build_bsai_nork_input(dat), dat)
bridge_list <- input_list
bridge_list$data <- cap_bsai_nork_srv_sel(bridge_list$data, dat)

# Stage 1: evaluate at the assessment's MLE ------------------------------------
obj <- fit_model(bridge_list$data, bridge_list$par, bridge_list$map,
                 do_optim = FALSE, silent = TRUE)
r <- obj$rep

naa <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
naa_obs <- cbind(naa[, 1:(n_obs_ages - 1)], rowSums(naa[, n_obs_ages:n_ages]))

cat("=== Stage 1: population at the assessment MLE ===\n")
cat("fishery selectivity max pct diff:",
    100 * max(abs(as.vector(r$fish_sel[1, 1, 1, 1, seq_len(n_obs_ages), 1, 1]) - dat$admb$sel_fsh) / dat$admb$sel_fsh), "\n")
cat("survey selectivity  max pct diff:",
    100 * max(abs(as.vector(r$srv_sel[1, 1, 1, 1, seq_len(n_obs_ages), 1, 1]) - dat$admb$sel_srv) / dat$admb$sel_srv), "\n")
cat("N     max pct diff:", 100 * max(abs(naa_obs - dat$admb$NAA) / dat$admb$NAA), "\n")
cat("SSB   max pct diff:", 100 * max(abs(as.vector(r$SSB)[1:n_yrs] - dat$admb$SSB) / dat$admb$SSB), "\n")
cat("Rec   max pct diff:", 100 * max(abs(as.vector(r$Rec)[1:n_yrs] - dat$admb$Rec) / dat$admb$Rec), "\n")
cat("Fmort max pct diff:", 100 * max(abs(as.vector(r$Fmort) - dat$admb$Fmort) / dat$admb$Fmort), "\n")

# The likelihood crosswalk. SPoRC writes each component as a proper density
# while the assessment drops normalising constants, so each comparison subtracts
# exactly the constants the assessment omits.
c2pi <- 0.5 * log(2 * pi)
n_catch_obs <- sum(dat$UseCatch)
n_recdev <- length(dat$mle$rec_dev)
n_fydev <- length(dat$mle$fydev)

srv_like <- sum(r$SrvIdx_nLL) - sum(c2pi + log(dat$ObsSrvIdx_SE[dat$UseSrvIdx == 1]))
f_pen <- as.vector(sum(r$Fmort_nLL) - n_catch_obs * (c2pi + as.vector(bridge_list$par$ln_sigmaF)))
rec_like <- sum(r$Rec_nLL) + sum(r$Init_Rec_nLL) - (n_recdev + n_fydev) * c2pi
sel_pri <- sum(r$sel_nLL) - log(sqrt(2 * pi) * 0.003)
catch_ssq <- sum(as.vector(r$Catch_nLL) - (c2pi + as.vector(bridge_list$par$ln_sigmaC)))

crosswalk <- data.frame(
  component = c("catch ssq", "survey index", "fishery age comps", "fishery length comps",
                "survey age comps", "recruitment", "F regularity", "survey selectivity prior",
                "M prior", "q prior"),
  SPoRC = c(catch_ssq, srv_like, sum(r$FishAgeComps_nLL), sum(r$FishLenComps_nLL),
            sum(r$SrvAgeComps_nLL), rec_like, f_pen, sel_pri,
            r$M_nLL - log(sqrt(2 * pi) * dat$cv_M),
            r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q)),
  assessment = c(dat$admb$datalikecomp[["catch.like"]],
                 dat$admb$datalikecomp[["aisurvlike"]],
                 dat$admb$datalikecomp[["fish.unbiased.ac"]],
                 dat$admb$datalikecomp[["fish.lc"]],
                 dat$admb$datalikecomp[["aisrv.ac"]],
                 dat$admb$pen_likecomp[["reclike"]],
                 dat$admb$pen_likecomp[["Fmortpen"]],
                 dat$admb$pen_likecomp[["prior_sel"]],
                 dat$admb$pen_likecomp[["prior_m"]],
                 dat$admb$pen_likecomp[["prior_q"]])
)

cat("\n=== Stage 1: likelihood components at the assessment MLE ===\n")
print(crosswalk, row.names = FALSE, digits = 6)
cat("like for like total:", sum(crosswalk$SPoRC), " assessment:",
    dat$admb$datalikecomp[["obj_fun"]] - dat$admb$datalikecomp[["mat_like"]], "\n")

# Stage 2: optimize -------------------------------------------------------------
est <- fit_model(input_list$data, input_list$par, input_list$map,
                 random = NULL, newton_loops = 3, silent = TRUE)
est$sdrep <- RTMB::sdreport(est)

cat("\n=== Stage 2: optimized ===\n")
cat("free parameters:", length(est$optim$par), "\n")
cat("final jnLL:", est$optim$objective, "  max |gradient|:", max(abs(est$gr(est$optim$par))), "\n")
cat("pdHess:", est$sdrep$pdHess, "\n")

# Stage 3: compare --------------------------------------------------------------
sdr <- est$sdrep
rep <- est$rep

ssb <- as.vector(rep$SSB)[1:n_yrs]
rec <- as.vector(rep$Rec)[1:n_yrs]

# The three terminal recruits are deviation free and the two models build them on
# different conventions, so the boundary is marked rather than left to read as
# drift.
ggplot2::ggsave(here("vignettes", "figures", "y_bsai_nork_ts_comparison.png"),
                bridge_ts_figure(yrs, ssb, rec, dat$admb$SSB, dat$admb$Rec, label,
                                 ssb_se = bridge_se(sdr, "log_SSB", ssb),
                                 rec_se = bridge_se(sdr, "log_Rec", rec),
                                 mark_year = yrs[n_yrs - 3] + 0.5),
                width = 17, height = 9, dpi = 150)

# Selectivity is time invariant, so a single curve per gear carries everything.
# The survey panel is where the refit and the assessment part company, because
# the assessment holds its curve flat past age 30 and the logistic does not.
sel_df <- bind_rows(
  bridge_sel_rows(dat$obs_ages, rep$fish_sel[1, 1, 1, 1, seq_len(n_obs_ages), 1, 1],
                  dat$admb$sel_fsh, "Fishery", label),
  bridge_sel_rows(dat$obs_ages, rep$srv_sel[1, 1, 1, 1, seq_len(n_obs_ages), 1, 1],
                  dat$admb$sel_srv, "Survey", label)
)

ggplot2::ggsave(here("vignettes", "figures", "y_bsai_nork_sel_comparison.png"),
                bridge_sel_figure(sel_df), width = 12, height = 7, dpi = 150)

est_yrs <- seq_len(n_yrs - 3)
cat("\n=== Stage 3: optimized SPoRC against the assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, dat$admb$SSB),
            bridge_cmp("Recruitment (estimated years)", rec[est_yrs], dat$admb$Rec[est_yrs])),
      row.names = FALSE, digits = 4)

# The terminal recruits differ by a fixed factor because the assessment builds
# them as the mean of the lognormal and SPoRC as the median.
term_yrs <- seq(n_yrs - 2, n_yrs)
cat("terminal recruit ratio (ADMB / SPoRC):",
    unique(round(dat$admb$Rec[term_yrs] / rec[term_yrs], 6)),
    "  exp(sigmaR^2 / 2):", round(exp(dat$sigmaR^2 / 2), 6), "\n")
