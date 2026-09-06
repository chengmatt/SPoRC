# Purpose: Bridge the 2024 BSAI Pacific ocean perch assessment to SPoRC and render its figures
# Creator: Matthew LH. Cheng
# Date Created: 8/8/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("tests", "testthat", "helper-bridge_bsai_pop.R"))
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

label <- "2024 BSAI Pacific Ocean Perch Assessment"

dat <- SPoRC::sgl_rg_bsai_pop_data
yrs <- dat$years
n_yrs <- length(yrs)
n_obs_ages <- length(dat$obs_ages)
nsel <- dat$nselages
yr_ind <- match(yrs, dat$admb$yrs)

input_list <- seed_bsai_pop_mle(build_bsai_pop_input(dat), dat)

# Stage 1: evaluate at the assessment's MLE ------------------------------------
obj <- fit_model(input_list$data, input_list$par, input_list$map,
                 do_optim = FALSE, silent = TRUE)
r <- obj$rep

# the assessment normalizes selectivity to a maximum of one within each year and multiplies F by
# the same factor, so both sides are put on the reported convention before comparing
sel_bridge <- r$fish_sel[1, 1, 1:n_yrs, 1, 1:nsel, 1, 1]
sel_bridge_norm <- norm_bsai_pop_sel_by_year(sel_bridge)
naa_all <- r$NAA[1, 1, 1:n_yrs, 1, , 1]

cat("=== Stage 1: population at the assessment MLE ===\n")
cat("fishery selectivity surface max pct diff:",
    100 * max(abs(sel_bridge_norm - dat$admb$sel_fsh) / dat$admb$sel_fsh), "\n")
cat("N     max pct diff:",
    100 * max(abs(as.vector(r$NAA[1, 1, 1, 1, 1:(n_obs_ages - 1), 1]) - dat$admb$NAA[1, 1:(n_obs_ages - 1)]) /
                dat$admb$NAA[1, 1:(n_obs_ages - 1)]), "\n")
cat("SSB   max pct diff:", 100 * max(abs(as.vector(r$SSB)[1:n_yrs] - dat$admb$SSB[yr_ind]) / dat$admb$SSB[yr_ind]), "\n")
cat("Rec   max pct diff:", 100 * max(abs(as.vector(r$Rec)[1:n_yrs] - dat$admb$Rec[yr_ind]) / dat$admb$Rec[yr_ind]), "\n")
cat("Total biomass max pct diff:",
    100 * max(abs(as.vector(naa_all %*% dat$pop_waa) - dat$admb$TotBiom[yr_ind]) / dat$admb$TotBiom[yr_ind]), "\n")

# the likelihood crosswalk. SPoRC writes each component as a proper density while the
# assessment drops normalizing constants, so each comparison subtracts those constants
c2pi <- 0.5 * log(2 * pi)
dlc <- unlist(dat$admb$datalikecomp)
plc <- unlist(dat$admb$pen_likecomp)

srv_like <- function(fl) {
  keep <- dat$UseSrvIdx[, , , fl] == 1
  sum(r$SrvIdx_nLL[, , , fl][keep]) - sum(c2pi + log(dat$ObsSrvIdx_SE[, , , fl][keep]))
}

n_recdev <- dim(input_list$par$ln_RecDevs)[3]
f_pen <- as.vector(sum(r$Fmort_nLL) - n_yrs * (c2pi + as.vector(input_list$par$ln_sigmaF)))
rec_like <- sum(r$Rec_nLL) + sum(r$Init_Rec_nLL) - n_recdev * c2pi + n_recdev * 0.5 * log(dat$sigmaR)
catch_ssq <- sum(as.vector(r$Catch_nLL) - (c2pi + as.vector(input_list$par$ln_sigmaC)))

crosswalk <- data.frame(
  component = c("catch ssq", "AI survey index", "EBS survey index", "fishery age comps",
                "fishery length comps", "AI survey age comps", "EBS survey age comps",
                "AI survey length comps", "recruitment", "F regularity",
                "selectivity penalty", "M prior", "AI q prior"),
  SPoRC = c(catch_ssq, srv_like(1), srv_like(2), sum(r$FishAgeComps_nLL),
            sum(r$FishLenComps_nLL), sum(r$SrvAgeComps_nLL[, , , , 1]),
            sum(r$SrvAgeComps_nLL[, , , , 2]), sum(r$SrvLenComps_nLL[, , , , 1]),
            rec_like, f_pen, sum(r$sel_nLL),
            r$M_nLL - log(sqrt(2 * pi) * dat$cv_M),
            r$srv_q_nLL - log(sqrt(2 * pi) * dat$cv_q[1])),
  assessment = c(dlc[["catch.like"]], dlc[["AI_survey.biom"]], dlc[["EBS_survey.biom"]],
                 dlc[["fish.ac"]], dlc[["fish.lc"]], dlc[["AI_survey.sac"]],
                 dlc[["EBS_survey.sac"]], dlc[["AI_survey.slc"]], plc[["reclike"]],
                 plc[["Fmortpen"]], dlc[["selpen"]], plc[["prior_M_bin_1"]],
                 plc[["AI_survey_prior_q_bin_1"]])
)

cat("\n=== Stage 1: likelihood components at the assessment MLE ===\n")
print(crosswalk, row.names = FALSE, digits = 6)
cat("like for like total:", sum(crosswalk$SPoRC), " assessment:",
    dlc[["obj_fun"]] - dlc[["mat_like"]], "\n")

# Stage 2: optimize -------------------------------------------------------------
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
cat("final jnLL:", est$optim$objective, "  max |gradient|:", max(abs(est$gr(est$optim$par))), "\n")
cat("pdHess:", est$sdrep$pdHess, "\n")

# Stage 3: compare --------------------------------------------------------------
sdr <- est$sdrep
rep <- est$rep

ssb <- as.vector(rep$SSB)[1:n_yrs]
rec <- as.vector(rep$Rec)[1:n_yrs]

# The terminal recruits are deviation free and hold the recruitment level shift,
# so the boundary is marked rather than left to read as a real change.
ggplot2::ggsave(
  here("vignettes", "figures", "aa_bsai_pop_ts_comparison.png"),
  bridge_ts_figure(
                  yrs,
                  ssb,
                  rec,
                  dat$admb$SSB[yr_ind],
                  dat$admb$Rec[yr_ind],
                  label,
                  ssb_se = bridge_se(sdr, "log_SSB", ssb),
                  rec_se = bridge_se(sdr, "log_Rec", rec),
                  mark_year = yrs[n_yrs - dat$fixedrec] + 0.5
                ),
  width = 17,
  height = 9,
  dpi = 150
)

# fishery selectivity is a bicubic surface over year and age, flat before 1964 and past the last
# node. survey selectivity is logistic and time invariant in each of the two surveys
sel_fit <- norm_bsai_pop_sel_by_year(rep$fish_sel[1, 1, 1:n_yrs, 1, 1:nsel, 1, 1])
show_yrs <- c(1960, 1980, 2000, yrs[n_yrs])

fsh_df <- bind_rows(
  lapply(show_yrs, function(y) {
    i <- match(y, yrs)
    bridge_sel_rows(dat$obs_ages, sel_fit[i, ], dat$admb$sel_fsh[i, ], y, label, facet_by = "Year")
  })
)

p_fsh <- bridge_sel_figure(
  fsh_df,
  facet_by = "Year",
  ylab = "Fishery selectivity",
  base_size = 18,
  nrow = 1
)

srv_logistic <- function(fl) {
  1 / (1 + exp(-dat$mle$sel_aslope_srv[fl] * (dat$obs_ages - dat$mle$sel_a50_srv[fl])))
}
srv_names <- c("Aleutian Islands survey", "Eastern Bering Sea survey")

srv_df <- bind_rows(
  lapply(seq_len(dat$n_srv_fleets), function(fl) {
    bridge_sel_rows(dat$obs_ages, rep$srv_sel[1, 1, 1, 1, seq_len(n_obs_ages), 1, fl],
                    srv_logistic(fl), srv_names[fl], label)
  })
)

p_srv <- bridge_sel_figure(
  srv_df,
  ylab = "Survey selectivity",
  base_size = 18,
  legend = "none",
  nrow = 1
)

ggplot2::ggsave(
  here("vignettes", "figures", "aa_bsai_pop_sel_comparison.png"),
  patchwork::wrap_plots(p_fsh, p_srv, nrow = 2, heights = c(1, 0.9)),
  width = 15,
  height = 10,
  dpi = 150
)

est_yrs <- seq_len(n_yrs - dat$fixedrec)
cat("\n=== Stage 3: optimized SPoRC against the assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, dat$admb$SSB[yr_ind]),
            bridge_cmp("Recruitment (estimated years)", rec[est_yrs], dat$admb$Rec[yr_ind][est_yrs])),
      row.names = FALSE, digits = 4)

# the terminal recruits are deviation free and differ by the shift in recruitment level, which the
# assessment's sum-to-zero dev_vector cannot make and SPoRC's free deviations can
fit_par <- est$env$parList(est$env$last.par.best)
term_yrs <- seq(n_yrs - dat$fixedrec + 1, n_yrs)
r0_shift <- (dat$mle$mean_log_rec + dat$sigmaR^2 / 2) - as.vector(fit_par$ln_global_R0)[1]
cat("terminal recruit ratio (ADMB / SPoRC):",
    unique(round(dat$admb$Rec[yr_ind][term_yrs] / rec[term_yrs], 6)),
    "  exp(R0 shift):", round(exp(r0_shift), 6), "\n")
