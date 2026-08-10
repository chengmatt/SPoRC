# Purpose: Bridge the 2024 GOA northern rockfish assessment to SPoRC and render
#          the case study figures. The configuration is not restated here: it is
#          sourced from the bridge test helper, so the figures, the bridge test,
#          and the pinned regression test all run the same specification and a
#          change to one cannot silently leave the others behind.
#
#          Three stages: set every parameter to the assessment's maximum
#          likelihood estimate and check the population and likelihood there,
#          optimize, then compare.
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("tests", "testthat", "helper-bridge_goa_nork.R"))
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

label <- "2024 Northern Rockfish Assessment"

dat <- SPoRC::sgl_rg_goa_nork_data
yrs <- dat$years
n_yrs <- length(yrs)
n_ages <- length(dat$ages)

input_list <- seed_goa_nork_mle(build_goa_nork_input(dat), dat)

# Stage 1: evaluate at the assessment's MLE ------------------------------------
obj <- fit_model(input_list$data, input_list$par, input_list$map,
                 do_optim = FALSE, silent = TRUE)
r <- obj$rep

cat("=== Stage 1: population at the assessment MLE ===\n")
cat("fishery selectivity max pct diff:",
    100 * max(abs(as.vector(r$fish_sel[1, 1, 1, 1, , 1, 1]) - dat$admb$sel[, 1]) / dat$admb$sel[, 1]), "\n")
cat("survey selectivity  max pct diff:",
    100 * max(abs(as.vector(r$srv_sel[1, 1, 1, 1, , 1, 1]) - dat$admb$sel[, 2]) / dat$admb$sel[, 2]), "\n")
cat("N   max pct diff:",
    100 * max(abs(r$NAA[1, 1, 1:n_yrs, 1, , 1] - t(dat$admb$NAA)) / t(dat$admb$NAA)), "\n")
cat("SSB max pct diff:",
    100 * max(abs(as.vector(r$SSB)[1:n_yrs] - dat$admb$SSB) / dat$admb$SSB), "\n")
cat("Rec max pct diff:",
    100 * max(abs(as.vector(r$Rec)[1:n_yrs] - dat$admb$Rec) / dat$admb$Rec), "\n")

# The likelihood crosswalk. SPoRC writes each component as a proper density
# while the assessment drops normalising constants, so each comparison subtracts
# exactly the constants the assessment omits.
ln_sigmaC <- input_list$par$ln_sigmaC
n_catch_obs <- sum(dat$UseCatch)
n_srv_obs <- sum(dat$UseSrvIdx)
n_recdev <- length(dat$mle$log_Rt)
n_initdev_pen <- n_ages - 2
s <- dat$sigmaR

crosswalk <- data.frame(
  component = c("catch ssq", "survey index", "fishery age comps", "fishery length comps",
                "survey age comps", "recruitment", "F regularity", "M prior", "q prior"),
  SPoRC = c(
    sum(as.vector(r$Catch_nLL) - (0.5 * log(2 * pi) + as.vector(ln_sigmaC))),
    dat$srv_wt * (sum(r$SrvIdx_nLL) - n_srv_obs * 0.5 * log(2 * pi)),
    sum(r$FishAgeComps_nLL),
    sum(r$FishLenComps_nLL),
    sum(r$SrvAgeComps_nLL),
    sum(r$Rec_nLL) - n_recdev * (0.5 * log(2 * pi) + log(s)) +
      sum(r$Init_Rec_nLL) - n_initdev_pen * (0.5 * log(2 * pi) + log(s)),
    dat$fmort_wt * (sum(r$Fmort_nLL) - n_catch_obs * (0.5 * log(2 * pi) + as.vector(input_list$par$ln_sigmaF))),
    r$M_nLL - log(sqrt(2 * pi) * dat$cv_M),
    r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q)
  ),
  assessment = c(dat$admb$ssqcatch, dat$admb$like_srv, dat$admb$like_fish_age,
                 dat$admb$like_fish_size, dat$admb$like_srv_age, dat$admb$like_rec,
                 dat$admb$f_regularity, dat$admb$nll_M, dat$admb$nll_q)
)

cat("\n=== Stage 1: likelihood components at the assessment MLE ===\n")
print(crosswalk, row.names = FALSE, digits = 6)

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

ggplot2::ggsave(here("vignettes", "figures", "w_goa_nork_ts_comparison.png"),
                bridge_ts_figure(yrs, ssb, rec, dat$admb$SSB, dat$admb$Rec, label,
                                 ssb_se = bridge_se(sdr, "log_SSB", ssb),
                                 rec_se = bridge_se(sdr, "log_Rec", rec)),
                width = 17, height = 9, dpi = 150)

# Selectivity is time invariant, so a single curve per gear carries everything.
sel_df <- bind_rows(
  bridge_sel_rows(dat$ages, rep$fish_sel[1, 1, 1, 1, , 1, 1], dat$admb$sel[, 1], "Fishery", label),
  bridge_sel_rows(dat$ages, rep$srv_sel[1, 1, 1, 1, , 1, 1], dat$admb$sel[, 2], "Survey", label)
)

ggplot2::ggsave(here("vignettes", "figures", "w_goa_nork_sel_comparison.png"),
                bridge_sel_figure(sel_df), width = 12, height = 6, dpi = 150)

cat("\n=== Stage 3: optimized SPoRC against the assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, dat$admb$SSB),
            bridge_cmp("Recruitment", rec, dat$admb$Rec)),
      row.names = FALSE, digits = 4)
