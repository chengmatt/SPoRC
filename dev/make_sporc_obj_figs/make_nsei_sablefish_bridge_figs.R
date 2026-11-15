# Purpose: To fit the NSEI sablefish TMB bridge in SPoRC and make the case study figures
# Creator: Matthew LH. Cheng
# Date: 11/13/26

# Setup -------------------------------------------------------------------

library(here)
library(ggplot2)
library(patchwork)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))
source(here("tests", "testthat", "helper-bridge_nsei_sablefish.R"))

# read in tmb data, parameters and reported values
tmb_data <- sgl_rg_nsei_sablefish_data$inputs
tmb_saved <- sgl_rg_nsei_sablefish_data$tmb_saved
tmb_converged <- sgl_rg_nsei_sablefish_data$tmb_converged
years <- 1975:2025
ages <- 2:31
n_yrs <- length(years)
label <- "TMB (converged)"

# Setup Model -------------------------------------------------------------

input_list <- build_nsei_sablefish_input(tmb_data = tmb_data, tmb_par = tmb_saved$par)
input_list <- seed_nsei_sablefish_tmb(input_list = input_list, tmb_par = tmb_saved$par)

# fit model, no newton steps since the 2017 - 2025 survey selectivity is knife edge
nsei_model <- fit_model(input_list$data,
                        input_list$par,
                        input_list$map,
                        random = NULL,
                        newton_loops = 0,
                        silent = TRUE
)

# get standard error report at the optimum, with the exact hessian
nsei_model$sdrep <- RTMB::sdreport(nsei_model,
                                   par.fixed = nsei_model$optim$par,
                                   hessian.fixed = nsei_model$he(nsei_model$optim$par))
sporc_sd <- summary(nsei_model$sdrep, "report")
log_ssb_sd <- sporc_sd[rownames(sporc_sd) == "log_SSB", 2]
log_rec_sd <- sporc_sd[rownames(sporc_sd) == "log_Rec", 2]

ssb <- nsei_model$rep$SSB[1,1,] / 1e3 # tonnes
rec <- nsei_model$rep$Rec[1,1,] / 1e6 # millions

# SSB and Recruitment -----------------------------------------------------

ssb_recruitment_plot <- bridge_ts_figure(yrs = years,
                                         ssb = ssb,
                                         rec = rec,
                                         admb_ssb = tmb_converged$ssb / 1e3,
                                         admb_rec = tmb_converged$pred_rec / 1e6,
                                         label = label,
                                         ssb_se = ssb * log_ssb_sd, # delta method, natural scale
                                         rec_se = rec * log_rec_sd,
                                         ref_name = "TMB"
)
ggsave(here("vignettes", "figures", "am_nsei_sablefish_ts.png"), ssb_recruitment_plot, width = 17, height = 9, dpi = 150)

# Index Fits --------------------------------------------------------------

srv_cpue_yrs <- tmb_data$yrs_srv_cpue + 1 # tmb years start at 0
fsh_cpue_yrs <- tmb_data$yrs_fsh_cpue + 1
mr_yrs <- tmb_data$yrs_mr + 1
index_names <- c("Longline survey CPUE", "Fishery CPUE", "Mark-recapture (millions)")

# observed indices
idx_obs <- rbind(
  data.frame(Year = years[srv_cpue_yrs], Index = index_names[1], Obs = tmb_data$data_srv_cpue, SE = tmb_data$sigma_srv_cpue),
  data.frame(Year = years[fsh_cpue_yrs], Index = index_names[2], Obs = tmb_data$data_fsh_cpue, SE = tmb_data$sigma_fsh_cpue),
  data.frame(Year = years[mr_yrs], Index = index_names[3], Obs = tmb_data$data_mr, SE = tmb_data$sigma_mr)
)
idx_obs$lwr <- idx_obs$Obs * exp(-1.96 * idx_obs$SE) # 95% lognormal interval
idx_obs$upr <- idx_obs$Obs * exp(1.96 * idx_obs$SE)

# predicted indices (fishery cpue and mark-recapture are fishery indices on fleets 1 and 2)
idx_pred <- rbind(
  data.frame(Year = years[srv_cpue_yrs], Index = index_names[1],
             SPoRC = nsei_model$rep$PredSrvIdx[1,1,srv_cpue_yrs,1,1], TMB = tmb_converged$pred_srv_cpue),
  data.frame(Year = years[fsh_cpue_yrs], Index = index_names[2],
             SPoRC = nsei_model$rep$PredFishIdx[1,1,fsh_cpue_yrs,1,1], TMB = tmb_converged$pred_fsh_cpue),
  data.frame(Year = years[mr_yrs], Index = index_names[3],
             SPoRC = nsei_model$rep$PredFishIdx[1,1,mr_yrs,1,2], TMB = tmb_converged$pred_mr)
)
idx_pred$pct_diff <- (idx_pred$SPoRC - idx_pred$TMB) / idx_pred$TMB * 100

idx_long <- reshape2::melt(idx_pred[, c("Year", "Index", "SPoRC", "TMB")], id.vars = c("Year", "Index"), variable.name = "Type")
levels(idx_long$Type) <- c("SPoRC", label)
idx_obs$Index <- factor(idx_obs$Index, levels = index_names)
idx_long$Index <- factor(idx_long$Index, levels = index_names)
idx_pred$Index <- factor(idx_pred$Index, levels = index_names)

idx_plot <- ggplot() +
  geom_pointrange(data = idx_obs, aes(x = Year, y = Obs, ymin = lwr, ymax = upr), color = "grey55", size = 0.4) +
  geom_line(data = idx_long, aes(x = Year, y = value, color = Type, linetype = Type), linewidth = 1) +
  facet_wrap(~Index, scales = "free_y", ncol = 1) +
  ggthemes::scale_color_colorblind() +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Index")

# percent difference
idx_pd_plot <- ggplot(idx_pred, aes(x = Year, y = pct_diff)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(color = "#E69F00") +
  geom_point(size = 2, color = "#E69F00") +
  facet_wrap(~Index, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 20) +
  labs(x = "Year", y = "SPoRC vs TMB (%)")

ggsave(here("vignettes", "figures", "am_nsei_sablefish_indices.png"),
       idx_plot + idx_pd_plot + plot_layout(widths = c(1.25, 1)),
       width = 17, height = 12, dpi = 150)

# Survey Selectivity ------------------------------------------------------

# longline survey selectivity by block and sex
blk_labels <- c("1975-1999", "2000-2016", "2017-2025")
blk_first_yr <- c(1, 26, 43)
sex_labels <- c("female", "male")
sel_df <- NULL
for(b in 1:3) {
  for(s in 1:2) {
    a50 <- exp(tmb_converged$par$log_srv_slx_pars[b,1,3 - s]) # tmb arrays are male first
    slope <- exp(tmb_converged$par$log_srv_slx_pars[b,2,3 - s])
    tmb_sel <- 1 / (1 + exp(-slope * (ages - 2 - a50))) # tmb logistic is in age - 2
    sel_df <- rbind(sel_df, bridge_sel_rows(ages = ages,
                                            sporc = nsei_model$rep$srv_sel[1,1,blk_first_yr[b],1,,s,1],
                                            admb = tmb_sel,
                                            gear = paste(blk_labels[b], sex_labels[s], sep = ", "),
                                            label = label))
  } # end s loop
} # end b loop
sel_df <- sel_df[sel_df$Age <= 15, ]
sel_df$Gear <- factor(sel_df$Gear, levels = unique(sel_df$Gear))

ggsave(here("vignettes", "figures", "am_nsei_sablefish_srv_sel.png"), bridge_sel_figure(sel_df, nrow = 3), width = 14, height = 12, dpi = 150)

# Uncertainty -------------------------------------------------------------

# CV of SSB and recruitment (sd of log in SPoRC, sd / estimate in tmb)
cv_df <- rbind(data.frame(Year = years, value = log_ssb_sd, Par = "Spawning Biomass", Type = "SPoRC"),
               data.frame(Year = years, value = tmb_converged$ssb_sd / tmb_converged$ssb, Par = "Spawning Biomass", Type = label),
               data.frame(Year = years, value = log_rec_sd, Par = "Recruitment", Type = "SPoRC"),
               data.frame(Year = years, value = tmb_converged$rec_sd / tmb_converged$pred_rec, Par = "Recruitment", Type = label)
)

cv_plot <- ggplot(cv_df, aes(x = Year, y = value, color = Type, linetype = Type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Par, scales = "free_y", ncol = 1) +
  ggthemes::scale_color_colorblind() +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Coefficient of variation")
ggsave(here("vignettes", "figures", "am_nsei_sablefish_cv.png"), cv_plot, width = 12, height = 9, dpi = 150)
