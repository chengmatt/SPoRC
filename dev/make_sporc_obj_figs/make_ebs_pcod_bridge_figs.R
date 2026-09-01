# Purpose: Bridge the 2024 EBS Pacific cod assessment (Stock Synthesis Model
#          24.1) to SPoRC from the packaged data object, and render the case
#          study figures. The model is specified the way the assessment
#          specifies it: one area and one sex, ages 0-20 with age 0 recruiting
#          at the start of the year, mean recruitment with the bias ramp and
#          twenty early deviations setting the start year's ages, a separate
#          initial equilibrium recruitment penalized toward the recruitment
#          level, Richards growth with annual deviations on the length at age
#          1.5 and on K from 2000 carried cohort by cohort, length based double
#          normal selectivity with two fishery blocks and annual deviations on
#          the survey's ascending width, and compositions recorded on 24 five
#          centimeter bins mapped from the model's 121 population bins.
#
#          Three stages: set every parameter to the assessment's maximum
#          likelihood estimate and check the population, the growth and the
#          likelihood components there, optimize, then compare.
#
#          The model is built through the same helper the regression test uses,
#          tests/testthat/helper-bridge_ebs_pcod.R, so the two cannot drift
#          apart. That helper is the readable specification; see also
#          vignette("ae_ebs_pacific_cod_case_study").
# Creator: Matthew LH. Cheng
# Date Created: 8/23/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))
source(here("tests", "testthat", "helper-bridge_ebs_pcod.R"))

label <- "2024 EBS Pacific Cod Assessment"

dat <- SPoRC::sgl_rg_ebs_pcod_data
yrs <- dat$years; n_yrs <- length(yrs)
ages <- dat$ages; n_ages <- length(ages)
s3 <- dat$ss3
LBM <- dat$LenBinMap
dat_lens <- dat$dat_lens_lower + 2.5 # midpoints of the composition bins

## Stage 1: evaluate at the assessment's estimate -----------------------------
input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)
at_mle <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
r <- at_mle$rep

library(Matrix)
image(at_mle$env$spHess(random = T))

pct <- function(a, b) 100 * max(abs(a / b - 1), na.rm = TRUE)
yr_row <- function(m, y) { rr <- as.integer(rownames(m)); m[as.character(max(rr[rr <= y])), ] }

cat("== At the assessment's estimate ==\n")
cat("numbers at age      ", sprintf("%.4g%%", pct(r$NAA[1,1,1:n_yrs,1,,1], s3$NAA)), "\n")
cat("spawning biomass    ", sprintf("%.4g%%", pct(2 * r$SSB[1,1,1:n_yrs], s3$SSB)), "\n")
cat("recruitment         ", sprintf("%.4g%%", pct(r$Rec[1,1,], s3$Rec)), "\n")
cat("catch               ", sprintf("%.4g%%", pct(r$PredCatch[1,1,,1,1], s3$dead_B)), "\n")
cat("survey index        ", sprintf("%.4g%%", pct(r$PredSrvIdx[1,1,match(s3$cpue$Yr, yrs),1,1], s3$cpue$Exp)), "\n")

## Stage 2: refit -------------------------------------------------------------
refit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = TRUE,
                   silent = F, newton_loops = 3)
# fit_model returns the RTMB object without an sdreport, so it is taken here.
# the analytical Hessian is passed rather than letting sdreport call optimHess,
# which reports a non positive definite Hessian on models where it is fine
sd_rep <- RTMB::sdreport(refit, hessian.fixed = refit$he(refit$env$last.par.best))
cat("pdHess", sd_rep$pdHess, "\n")
cat("\n== Refit ==\n")
cat("objective drop", at_mle$rep$jnLL - refit$rep$jnLL,
    " max |gradient|", signif(max(abs(refit$gr(refit$env$last.par.best))), 3), "\n")

## Stage 3: figures -----------------------------------------------------------
if(!dir.exists(here("vignettes", "figures"))) dir.create(here("vignettes", "figures"), recursive = TRUE)

## Spawning biomass and recruitment -------------------------------------------
# At the estimate and refitted. SPoRC's single
# sex spawning biomass is the female share of the population, so it is doubled
# to sit on the assessment's own definition.
ssb_mle <- 2 * as.vector(r$SSB[1,1,1:n_yrs]); rec_mle <- as.vector(r$Rec[1,1,])
ssb_ref <- 2 * as.vector(refit$rep$SSB[1,1,1:n_yrs]); rec_ref <- as.vector(refit$rep$Rec[1,1,])

# the first panel carries no interval: it is the model evaluated at the
# assessment's estimate rather than fitted, so the curvature there is not an
# uncertainty statement about anything SPoRC estimated
p_ts <- bridge_ts_figure(yrs = yrs, ssb = ssb_mle, rec = rec_mle,
                         admb_ssb = s3$SSB, admb_rec = s3$Rec, label = label,
                         ref_name = "Stock Synthesis")
ggplot2::ggsave(here("vignettes", "figures", "ae_ebs_pcod_ts_comparison.png"), p_ts,
                width = 15, height = 10, dpi = 170)

# the refit does, from the sdreport. bridge_se reads the log scale standard
# deviation and returns it on the natural scale, so passing the doubled spawning
# biomass gives the doubled interval with it
p_ts_ref <- bridge_ts_figure(yrs = yrs, ssb = ssb_ref, rec = rec_ref,
                             admb_ssb = s3$SSB, admb_rec = s3$Rec, label = label,
                             ref_name = "Stock Synthesis",
                             ssb_se = bridge_se(sd_rep, "log_SSB", ssb_ref),
                             rec_se = bridge_se(sd_rep, "log_Rec", rec_ref))
ggplot2::ggsave(here("vignettes", "figures", "ae_ebs_pcod_ts_refit.png"), p_ts_ref,
                width = 15, height = 10, dpi = 170)

## Selectivity at length ------------------------------------------------------
# The fishery's two blocks and the survey in years its ascending width moved.
# Both are on the population's length bins, not the coarser data bins.
sel_rows <- function(y, fleet, gear) {
  sporc <- if(fleet == 1) r$fish_sel_l[1, match(y, yrs), , 1, 1] else r$srv_sel_l[1, match(y, yrs), , 1, 1]
  ss3 <- yr_row(s3$lsel[[fleet]], y)
  dplyr::bind_rows(
    data.frame(Age = dat$lens, value = as.vector(sporc), Gear = gear, type = "SPoRC"),
    data.frame(Age = dat$lens, value = as.vector(ss3), Gear = gear, type = label)
  )
}
sel_df <- dplyr::bind_rows(
  sel_rows(1977, 1, "Fishery, 1977-1989"),
  sel_rows(1990, 1, "Fishery, 1990-2024"),
  sel_rows(1982, 2, "Survey, 1982"),
  sel_rows(2000, 2, "Survey, 2000"),
  sel_rows(2012, 2, "Survey, 2012"),
  sel_rows(2024, 2, "Survey, 2024")
)
ggplot2::ggsave(here("vignettes", "figures", "ae_ebs_pcod_sel_comparison.png"),
                bridge_sel_figure(sel_df, nrow = 2) + ggplot2::labs(x = "Length (cm)"),
                width = 14, height = 9, dpi = 170)

## Mean length at age ---------------------------------------------------------
# Years spanning the time-varying period. Under cohort growth a late year only
# reproduces if every earlier year's increment did, so 2023 is the strong test.
gs <- s3$growthseries
growth_df <- dplyr::bind_rows(lapply(c(1977, 2000, 2005, 2010, 2016, 2023), function(y) {
  iy <- match(y, yrs)
  dplyr::bind_rows(
    data.frame(Age = ages, value = r$mean_LAA_spawn[1,1,iy,1,,1], Gear = paste("Year", y), type = "SPoRC"),
    data.frame(Age = ages, value = as.numeric(gs[gs$Yr == y, as.character(ages)]), Gear = paste("Year", y), type = label)
  )
}))
ggplot2::ggsave(here("vignettes", "figures", "ae_ebs_pcod_growth_comparison.png"),
                bridge_sel_figure(growth_df, ylab = "Mean length at age (cm)", nrow = 2),
                width = 14, height = 9, dpi = 170)

## The two varying growth parameters ------------------------------------------
# Through the series, against the assessment's own realized values
gy <- s3$growth_by_year
tv_df <- dplyr::bind_rows(
  data.frame(Year = yrs, value = r$growth_pars_y[1,1,,1,1], Par = "Length at age 1.5 (cm)", type = "SPoRC"),
  data.frame(Year = yrs, value = gy$L1[match(yrs, gy$Yr)], Par = "Length at age 1.5 (cm)", type = label),
  data.frame(Year = yrs, value = r$growth_pars_y[1,1,,3,1], Par = "Richards K", type = "SPoRC"),
  data.frame(Year = yrs, value = gy$K[match(yrs, gy$Yr)], Par = "Richards K", type = label)
)
p_tv <- ggplot2::ggplot(tv_df, ggplot2::aes(x = Year, y = value, color = type, linetype = type)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::facet_wrap(~Par, scales = "free_y", ncol = 1) +
  ggthemes::scale_color_colorblind() +
  ggplot2::theme_bw(base_size = 20) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::labs(x = "Year", y = "Value", color = "Type", linetype = "Type")
ggplot2::ggsave(here("vignettes", "figures", "ae_ebs_pcod_growth_tv.png"), p_tv,
                width = 12, height = 9, dpi = 170)

## Expected length compositions -----------------------------------------------
# On the bins the data are recorded on. The fishery's are formed at mid season
# on the season-long catch with the length selectivity applied at length.
exp_len <- function(v) { w <- as.vector(v %*% LBM); w / sum(w) }
ld <- s3$lendbase
comp_rows <- function(y, fleet, gear) {
  if(!any(ld$Yr == y & ld$Fleet == fleet)) return(NULL)
  sporc <- if(fleet == 1) exp_len(r$CAL[1,1,match(y, yrs),1,,1,1]) else exp_len(r$SrvIAL[1,1,match(y, yrs),1,,1,1])
  dplyr::bind_rows(
    data.frame(Age = dat_lens, value = sporc, Gear = gear, type = "SPoRC"),
    data.frame(Age = dat_lens, value = ld$Exp[ld$Yr == y & ld$Fleet == fleet], Gear = gear, type = label)
  )
}
# the survey's length compositions are fit in 1982-1999 and 2024 only
comp_df <- dplyr::bind_rows(
  comp_rows(1990, 1, "Fishery, 1990"), comp_rows(2005, 1, "Fishery, 2005"), comp_rows(2024, 1, "Fishery, 2024"),
  comp_rows(1985, 2, "Survey, 1985"), comp_rows(1995, 2, "Survey, 1995"), comp_rows(2024, 2, "Survey, 2024")
)
ggplot2::ggsave(here("vignettes", "figures", "ae_ebs_pcod_comp_comparison.png"),
                bridge_sel_figure(comp_df, ylab = "Proportion", nrow = 2) + ggplot2::labs(x = "Length (cm)"),
                width = 15, height = 9, dpi = 170)

cat("\nfigures written to vignettes/figures/ae_ebs_pcod_*.png\n")
