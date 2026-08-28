# Purpose: Bridge the 2025 West Coast sablefish assessment to SPoRC and render
#          the case study figures. The configuration is not restated here: it is
#          sourced from the bridge test helper, so the figures, the bridge test
#          and the case study all run the same specification and a change to one
#          cannot silently leave the others behind.
#
#          Three stages: set every parameter to the assessment's maximum
#          likelihood estimate and check the population and the likelihoods
#          there, optimize, then compare. Selectivity is estimated, on SPoRC's
#          own double normal with the assessment's time blocks, mirrored fleets
#          and male offsets, so the assessment's selectivity parameters go in
#          as starting values rather than its curves going in as data.
# Creator: Matthew LH. Cheng
# Date Created: 8/21/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("tests", "testthat", "helper-bridge_wc_sablefish.R"))
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

label <- "2025 West Coast Sablefish Assessment"

dat <- SPoRC::sgl_rg_wc_sablefish_data
yrs <- dat$years
n_yrs <- length(yrs)
n_ages <- length(dat$ages)
n_obs_ages <- length(dat$obs_ages)
n_fish <- dat$n_fish_fleets
n_srv <- dat$n_srv_fleets

input_list <- seed_wc_sablefish_mle(build_wc_sablefish_input(dat), dat)

# the assessment's own selectivity surfaces, for comparison only
sel_fish_ss3 <- expand_wc_sablefish_sel(dat$fish_sel_blocks_ss3, n_yrs, n_ages, dat$n_sexes)
sel_srv_ss3 <- expand_wc_sablefish_sel(dat$srv_sel_blocks_ss3, n_yrs, n_ages, dat$n_sexes)
# the recruitment index fleet reads no curve, so it is left out of this
i_sel_srv <- which(dat$srv_src != 11)
sel_gap <- function(rep) c(fishery = max(abs(rep$fish_sel[1, 1, 1:n_yrs, 1, , , ] - sel_fish_ss3[1, 1, , 1, , , ])),
                           survey = max(abs(rep$srv_sel[1, 1, 1:n_yrs, 1, , , i_sel_srv] - sel_srv_ss3[1, 1, , 1, , , i_sel_srv])))

# Stage 1: evaluate at the assessment's estimate -------------------------------
obj <- fit_model(input_list$data, input_list$par, input_list$map,
                 do_optim = FALSE, silent = TRUE)
r <- obj$rep

pd <- function(a, b) 100 * (a - b) / b
pd_max <- function(a, b) max(abs(pd(a, b)))

naa_sp <- r$NAA[1, 1, 1:n_yrs, 1, , ]
totb_jan1 <- sapply(1:n_yrs, function(y) sum(naa_sp[y, , ] * dat$WAA[1, 1, y, 1, , ]))
i_catch <- dat$UseCatch[1, , 1, 1:6] == 1

cat("=== Stage 1: at the assessment's estimate (its report carries six significant digits) ===\n")
stage1 <- data.frame(
  quantity = c("numbers at age", "spawning biomass", "recruitment", "total biomass (Jan 1)",
               "predicted catch", "survey indices"),
  max_pct = c(pd_max(as.vector(naa_sp), as.vector(dat$ss3$NAA)),
              pd_max(as.vector(r$SSB[1, 1, 1:n_yrs]), dat$ss3$SSB),
              pd_max(as.vector(r$Rec[1, 1, 1:n_yrs]), dat$ss3$Rec),
              pd_max(totb_jan1, dat$ss3$Bio_all),
              pd_max(as.vector(r$PredCatch[1, 1, , 1, 1:6])[i_catch], as.vector(dat$ss3$pred_catch)[i_catch]),
              pd_max(unlist(lapply(1:4, function(sf) {
                       ci <- dat$ss3$pred_idx %>% dplyr::filter(Fleet == dat$srv_src[sf])
                       r$PredSrvIdx[1, 1, match(ci$Yr, yrs), 1, sf] })),
                     unlist(lapply(1:4, function(sf) {
                       (dat$ss3$pred_idx %>% dplyr::filter(Fleet == dat$srv_src[sf]))$Exp })))))
print(stage1, row.names = FALSE, digits = 3)
cat("bias ramp against the assessment's own, max absolute difference:",
    signif(max(abs(r$bias_ramp - dat$mle$bias_adj)), 3), "\n")
cat("selectivity against the assessment's own surfaces, max absolute difference: fishery",
    signif(sel_gap(r)[["fishery"]], 3), " survey", signif(sel_gap(r)[["survey"]], 3), "\n")

# Expected compositions, built the way the likelihood forms them and put through
# the composition constant so they compare with the assessment's own table
exp_comp <- function(p_f, p_m, joint) {
  if(joint) {
    p <- c(p_f, p_m) / sum(p_f + p_m)
    e <- c(as.vector(p[1:n_ages] %*% dat$AgeingError), as.vector(p[n_ages + 1:n_ages] %*% dat$AgeingError))
  } else e <- as.vector(((p_f + p_m) / sum(p_f + p_m)) %*% dat$AgeingError)
  e <- e / sum(e)
  (e + dat$addtocomp) / sum(e + dat$addtocomp)
}
comp_gap <- function(src, sexcode, numbers) {
  sapply(seq_along(src), function(f) {
    db <- dat$ss3$agedbase %>% dplyr::filter(Fleet == src[f], Sexes == sexcode[f])
    yy <- sort(unique(db$Yr))
    if(length(yy) == 0) return(0) # a fleet with no compositions, such as the recruitment index
    max(sapply(yy, function(y) {
      e <- exp_comp(numbers(match(y, yrs), 1, f), numbers(match(y, yrs), 2, f), sexcode[f] == 3)
      max(abs(e - (db %>% dplyr::filter(Yr == y) %>% dplyr::arrange(Sex, Bin))$Exp))
    }))
  })
}
cat("expected age compositions against the assessment's, max absolute difference:",
    signif(max(c(comp_gap(dat$fish_src, dat$fish_sex, function(y, s, f) r$CAA[1, 1, y, 1, , s, f]),
                 comp_gap(dat$srv_src, dat$srv_sex, function(y, s, f) r$SrvIAA[1, 1, y, 1, , s, f]))), 3), "\n")

# Likelihoods. SPoRC writes each component as a proper density and the
# assessment drops normalizing constants, so the comparison subtracts exactly
# the constants it omits.
lc <- function(sigma, n) n * (log(sigma) + 0.5 * log(2 * pi))
n_est_dev <- sum(yrs %in% dat$yrs_rec_est)
lik <- data.frame(
  component = c("catch", paste0("index: ", dat$fleet_names[c(7:10, 11)]),
                paste0("ages: ", dat$fleet_names[1:6]), paste0("ages: ", dat$fleet_names[7:10]),
                "recruitment deviations", "natural mortality prior"),
  SPoRC = c(sum(r$Catch_nLL) - lc(dat$sigmaC, sum(dat$UseCatch == 1, na.rm = TRUE)),
            sapply(c(1:4, n_srv), function(sf) sum(r$SrvIdx_nLL[, , , sf]) - 0.5 * log(2 * pi) * sum(dat$UseSrvIdx[, , , sf])),
            sum(r$FishAgeComps_nLL[, , , , c(1, 7)]), sapply(2:6, function(f) sum(r$FishAgeComps_nLL[, , , , f])),
            sapply(1:3, function(sf) sum(r$SrvAgeComps_nLL[, , , , sf])), sum(r$SrvAgeComps_nLL[, , , , 4:5]),
            sum(r$Rec_nLL) - n_est_dev * 0.5 * log(2 * pi) - 0.5 * sum(dat$mle$bias_adj) * log(dat$sigmaR),
            r$M_nLL - lc(dat$M_prior$sd, 1)),
  assessment = c(dat$ss3$lik$catch, dat$ss3$lik$index[1:5], dat$ss3$lik$age,
                 dat$ss3$lik$recruitment + dat$ss3$lik$forecast_recruitment - sum(dat$mle$bias_adj) * log(dat$sigmaR),
                 dat$ss3$lik$priors))
lik$difference <- lik$SPoRC - lik$assessment
cat("\n=== Stage 1: likelihoods, SPoRC net of the constants the assessment omits ===\n")
print(lik, row.names = FALSE, digits = 8)

const <- lc(dat$sigmaC, sum(dat$UseCatch == 1, na.rm = TRUE)) + 0.5 * log(2 * pi) * sum(dat$UseSrvIdx) +
  n_est_dev * 0.5 * log(2 * pi) - 0.5 * sum(dat$mle$bias_adj) * log(dat$sigmaR) + lc(dat$M_prior$sd, 1)
cat(sprintf("\ntotal objective: SPoRC %.8f net of constants against the assessment's %.8f (difference %.2e)\n",
            obj$fn(obj$par) - const, dat$mle$objective, obj$fn(obj$par) - const - dat$mle$objective))

# The assessment solves fishing mortality from the catch where SPoRC estimates
# it, so SPoRC's partials at its estimate still carry the dependence of F on
# every other parameter. The comparable quantity profiles them out.
H <- as.matrix(obj$he(obj$par)); g <- as.vector(obj$gr(obj$par)); nms <- names(obj$par)
iF <- which(nms == "ln_F_devs"); iO <- which(nms != "ln_F_devs")
g_prof <- as.vector(g[iO] - H[iO, iF] %*% solve(H[iF, iF], g[iF]))
cat("\nmax |gradient| by parameter block, profiled over the fishing mortality deviations:\n")
print(round(tapply(abs(g_prof), nms[iO], max), 5))

# Stage 2: optimize ------------------------------------------------------------
est <- fit_model(input_list$data, input_list$par, input_list$map,
                 do_optim = TRUE, newton_loops = 3, silent = TRUE)
est$sdrep <- tryCatch(RTMB::sdreport(est, hessian.fixed = est$he(est$optim$par)), error = function(e) NULL)
pl <- est$env$parList(est$optim$par)
cat(sprintf("\n=== Stage 2: %d parameters (%d of them selectivity), objective %.6f, fallen %.4e from the assessment's estimate, max |gradient| %.2e ===\n",
            length(est$par), sum(names(est$par) %in% c("fish_fixed_sel_pars", "srv_fixed_sel_pars", "ln_fishsel_sex_scale")),
            est$optim$objective, obj$fn(obj$par) - est$optim$objective, max(abs(est$gr(est$optim$par)))))
cat("selectivity after fitting, max absolute difference from the assessment's: fishery",
    signif(sel_gap(est$rep)[["fishery"]], 3), " survey", signif(sel_gap(est$rep)[["survey"]], 3), "\n")

# Stage 3: compare -------------------------------------------------------------
ssb <- as.vector(est$rep$SSB[1, 1, 1:n_yrs])
rec <- as.vector(est$rep$Rec[1, 1, 1:n_yrs])
ssb_se <- if(is.null(est$sdrep)) NA else bridge_se(est$sdrep, "log_SSB", ssb, exact = TRUE)
rec_se <- if(is.null(est$sdrep)) NA else bridge_se(est$sdrep, "log_Rec", rec, exact = TRUE)

cat("\n=== Stage 3: optimized SPoRC against the assessment ===\n")
print(rbind(bridge_cmp("Spawning biomass", ssb, dat$ss3$SSB, signed = TRUE),
            bridge_cmp("Recruitment", rec, dat$ss3$Rec, signed = TRUE)), row.names = FALSE, digits = 4)
print(data.frame(
  quantity = c("ln R0", "natural mortality", paste0("catchability: ", dat$fleet_names[7:10]),
               "catchability: recruitment index"),
  SPoRC = c(pl$ln_global_R0[1], exp(pl$ln_M[1]), exp(pl$ln_srv_q[1, 1, c(1:4, n_srv)])),
  assessment = c(dat$mle$ln_R0, dat$mle$M, exp(dat$mle$ln_srv_q), dat$mle$q_rec_idx)),
  row.names = FALSE, digits = 7)

ggplot2::ggsave(here("vignettes", "figures", "ac_wc_sablefish_ts_comparison.png"),
                bridge_ts_figure(yrs, ssb, rec, dat$ss3$SSB, dat$ss3$Rec, label,
                                 ssb_se = ssb_se, rec_se = rec_se),
                width = 17, height = 9, dpi = 200)

# Selectivity, one panel per gear in the terminal year. Each gear's two sexes
# are drawn together, since the hook and line and pot fleets carry male offsets.
sel_df <- dplyr::bind_rows(lapply(1:6, function(f) {
  dplyr::bind_rows(
    bridge_sel_rows(dat$ages, est$rep$fish_sel[1, 1, n_yrs, 1, , 1, f], sel_fish_ss3[1, 1, n_yrs, 1, , 1, f],
                    paste0(gsub("_", " ", dat$fleet_names[f]), ", female"), label),
    bridge_sel_rows(dat$ages, est$rep$fish_sel[1, 1, n_yrs, 1, , 2, f], sel_fish_ss3[1, 1, n_yrs, 1, , 2, f],
                    paste0(gsub("_", " ", dat$fleet_names[f]), ", male"), label))
}))
sel_df <- dplyr::bind_rows(sel_df, dplyr::bind_rows(lapply(1:4, function(sf) {
  bridge_sel_rows(dat$ages, est$rep$srv_sel[1, 1, n_yrs, 1, , 1, sf], sel_srv_ss3[1, 1, n_yrs, 1, , 1, sf],
                  paste0(gsub("_", " ", dat$fleet_names[dat$srv_src[sf]]), ", female"), label)
})))
sel_df <- sel_df %>% dplyr::filter(Age <= 30)

ggplot2::ggsave(here("vignettes", "figures", "ac_wc_sablefish_sel_comparison.png"),
                bridge_sel_figure(sel_df, nrow = 4), width = 15, height = 12, dpi = 170)

cat("\nfigures written to vignettes/figures\n")
