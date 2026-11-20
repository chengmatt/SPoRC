# Purpose: To fit the EBS pollock Model 26.0 (Rceattle) bridge in SPoRC and make the case study figures
# Creator: Matthew LH. Cheng
# Date: 11/19/26

# Setup -------------------------------------------------------------------

library(here)
library(ggplot2)
library(patchwork)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))
source(here("tests", "testthat", "helper-bridge_ebs_pollock_m26.R"))

inputs <- sgl_rg_ebswp_m26_data$inputs
rce <- sgl_rg_ebswp_m26_data$rceattle
years <- inputs$years
n_yrs <- length(years)
label <- "Rceattle (Model 26.0)"

# Setup Model -------------------------------------------------------------

input_list <- build_ebs_pollock_m26_input(inputs = inputs, rce = rce)
input_list <- seed_ebs_pollock_m26_rceattle(input_list = input_list, rce = rce)

ebs_pollock_model <- fit_model(input_list$data,
                               input_list$par,
                               input_list$map,
                               newton_loops = 2,
                               silent = TRUE
)
ebs_pollock_model$sdrep <- RTMB::sdreport(ebs_pollock_model)
sporc_sd <- summary(ebs_pollock_model$sdrep, "report")
log_ssb_sd <- sporc_sd[rownames(sporc_sd) == "log_SSB", 2]
log_rec_sd <- sporc_sd[rownames(sporc_sd) == "log_Rec", 2]

rep <- ebs_pollock_model$rep
ssb <- as.vector(rep$SSB[1,1,]) # thousand tons
rec <- as.vector(rep$Rec[1,1,]) # millions

# numbers the case study quotes
ssb_pct <- 100 * (ssb / rce$ssb - 1)
rec_pct <- 100 * (rec / rce$rec - 1)
length(ebs_pollock_model$par) # number of parameters
ebs_pollock_model$optim$objective # objective
max(abs(ebs_pollock_model$gr(ebs_pollock_model$optim$par))) # gradient
ebs_pollock_model$sdrep$pdHess # positive definite hessian
years[which.max(abs(ssb_pct))]; max(abs(ssb_pct)) # SSB, all years
max(abs(ssb_pct[years >= 1980 & years < 2000])) # 1980 - 1999
max(abs(ssb_pct[years >= 2000])) # from 2000
median(ssb_pct) # median
years[which.max(abs(rec_pct))]; max(abs(rec_pct)) # recruitment, all years
max(abs(rec_pct[years >= 1980 & years < 2000])) # 1980 - 1999
max(abs(rec_pct[years >= 2000])) # from 2000
median(rec_pct) # median

# SSB and Recruitment -----------------------------------------------------

ssb_recruitment_plot <- bridge_ts_figure(yrs = years,
                                         ssb = ssb,
                                         rec = rec,
                                         admb_ssb = rce$ssb,
                                         admb_rec = rce$rec,
                                         label = label,
                                         ssb_se = ssb * log_ssb_sd, # delta method, natural scale
                                         rec_se = rec * log_rec_sd,
                                         legend_nrow = NULL,
                                         ref_name = "Rceattle"
)
ggsave(here("vignettes", "figures", "f_ebs_pol_m26_ts.png"), ssb_recruitment_plot, width = 17, height = 9, dpi = 150)

# Fishery CPUE ------------------------------------------------------------

cpue_yrs <- which(inputs$UseFishIdx[1,,1,1] == 1)
cpue_df <- rbind(
  data.frame(Year = years[cpue_yrs], value = rep$PredFishIdx[1,1,cpue_yrs,1,1], type = "SPoRC"),
  data.frame(Year = years[cpue_yrs], value = rce$cpue_hat[cpue_yrs], type = label)
)
cpue_obs <- data.frame(Year = years[cpue_yrs],
                       value = inputs$ObsFishIdx[1,cpue_yrs,1,1],
                       se = inputs$ObsFishIdx_SE[1,cpue_yrs,1,1])

cpue_plot <- ggplot() +
  geom_pointrange(data = cpue_obs, aes(x = Year, y = value, ymin = value - 1.96 * se, ymax = value + 1.96 * se), color = "grey40") +
  geom_line(data = cpue_df, aes(x = Year, y = value, color = type, linetype = type), linewidth = 1.2) +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 20) +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Year", y = "Fishery CPUE", color = "Type", linetype = "Type")
ggsave(here("vignettes", "figures", "f_ebs_pol_m26_cpue.png"), cpue_plot, width = 12, height = 7, dpi = 150)

# Selectivity -------------------------------------------------------------

sel_ages <- c(1, 2, 3, 4, 6, 8)
fleet_names <- c("Fishery", "Bottom trawl", "Acoustic trawl")
sporc_sel <- list(rep$fish_sel[1,1,1:n_yrs,1,,1,1], rep$srv_sel[1,1,1:n_yrs,1,,1,1], rep$srv_sel[1,1,1:n_yrs,1,,1,2])
rce_sel <- list(rce$fish_sel, rce$bts_sel, rce$ats_sel)
first_year <- c(1964, 1982, 1994) # first year each fleet has data

sel_df <- NULL
for(f in 1:3) {
  keep <- years >= first_year[f]
  for(a in sel_ages) {
    if(f == 3 && a == 1) next # acoustic age 1 is not used here, and rceattle sets it to zero
    sel_df <- rbind(sel_df,
                    data.frame(Year = years[keep], Age = paste("Age", a), Fleet = fleet_names[f], value = sporc_sel[[f]][keep,a], type = "SPoRC"),
                    data.frame(Year = years[keep], Age = paste("Age", a), Fleet = fleet_names[f], value = rce_sel[[f]][keep,a], type = label))
  } # end a loop
} # end f loop
sel_df$Age <- factor(sel_df$Age, levels = paste("Age", sel_ages))
sel_df$Fleet <- factor(sel_df$Fleet, levels = fleet_names)

sel_plot <- ggplot(sel_df, aes(x = Year, y = value, color = type, linetype = type)) +
  geom_line(linewidth = 0.9) +
  facet_grid(Fleet ~ Age, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 16) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Selectivity", color = "Type", linetype = "Type")
ggsave(here("vignettes", "figures", "f_ebs_pol_m26_sel.png"), sel_plot, width = 17, height = 9, dpi = 150)
