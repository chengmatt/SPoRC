# Purpose: To plot simulation results comparing CTMC and unstructured markov
# Creator: Matthew LH. Cheng
# 10/31/25

# setup -------------------------------------------------------------------
library(here)
library(tidyverse)

# load in simulation data
const <- readRDS(here("dev", "dev_output", "CTMC", "const.RDS"))
age_move <- readRDS(here("dev", "dev_output", "CTMC", "age_move.RDS"))
time_move <- readRDS(here("dev", "dev_output", "CTMC", "time_move.RDS"))
timeage_move <- readRDS(here("dev", "dev_output", "CTMC", "timeage_move.RDS"))
age_env_move <- readRDS(here("dev", "dev_output", "CTMC", "age_env_move.RDS"))
sim_data <- list(const, age_move, time_move, timeage_move, age_env_move)

# load in simulation results
const_res <- readRDS(here("dev", "dev_output", "CTMC", "const_sim_results.RDS"))
age_move_res <- readRDS(here("dev", "dev_output", "CTMC", "age_move_sim_results.RDS"))
time_move_res <- readRDS(here("dev", "dev_output", "CTMC", "time_move_sim_results.RDS"))
timeage_move_res <- readRDS(here("dev", "dev_output", "CTMC", "timeage_move_sim_results.RDS"))
age_env_move_rec <- readRDS(here("dev", "dev_output", "CTMC", "age_env_move_sim_results.RDS"))
sim_model_res <- list(const_res, age_move_res, time_move_res, timeage_move_res, age_env_move_rec)

# scenario names
scenario <- c("const", "age", "time", "timeage", "age_env")

# extract results ---------------------------------------------------------

# empty containers
ssb_df <- data.frame()
move_df <- data.frame()

for(i in 1:length(sim_model_res)) {

  tmp_res <- sim_model_res[[i]] # get simulation results
  tmp_sim_dat <- sim_data[[i]] # get simulation data

  for(j in 1:length(tmp_res)) {

    # extract out report
    tmp_rep <- tmp_res[[j]]$model$rep

    # get ssb results
    tmp_ssb <- reshape2::melt(tmp_rep$SSB) %>%
      mutate(sim_id = tmp_res[[j]]$sim_id,
             model_id = tmp_res[[j]]$model_id,
             pd = tmp_res[[j]]$pdHess,
             grad = max(abs(tmp_res[[j]]$model$sd_rep$gradient.fixed)),
             scenario = scenario[i]) %>%
      left_join(
        reshape2::melt(tmp_sim_dat$SSB[,,tmp_res[[j]]$sim_id]) %>%
          rename(true = value),
        by = c("Var1", "Var2")
      )

    # get movement results
    tmp_move <- reshape2::melt(tmp_rep$Movement) %>%
      mutate(sim_id = tmp_res[[j]]$sim_id,
             model_id = tmp_res[[j]]$model_id,
             pd = tmp_res[[j]]$pdHess,
             grad = max(abs(tmp_res[[j]]$model$sd_rep$gradient.fixed)),
             scenario = scenario[i]) %>%
      left_join(
        reshape2::melt(tmp_sim_dat$Movement[,,,,1,tmp_res[[j]]$sim_id]) %>%
          rename(from = Var1, to = Var2, years = Var3, ages = Var4, true = value),
        by = c("from", "to", "years", "ages")
      )

    ssb_df <- rbind(ssb_df, tmp_ssb)
    move_df <- rbind(tmp_move, move_df)
    print(j)

  } # end j loop
} # end i loop

# Convergence -------------------------------------------------------------

# Figure out convergence
conv_rates <- ssb_df %>%
  filter(Var1 == 1, Var2 == 1, (pd == T & grad < 0.01)) %>%
  select(model_id, sim_id, pd, scenario) %>%
  group_by(model_id, scenario) %>%
  summarize(num_pd = sum(pd),
            prop_pd = num_pd / 100) %>%
  mutate(scenario = factor(scenario, levels = c("const", "age", "time", "timeage", "age_env")))

# get convergence rates for plotting
conv_rates_plot <- conv_rates %>%
  ungroup() %>%
  mutate(model_id = ifelse(model_id == 1, "unstrc_mark", "ctmc")) %>%
  mutate(region = 'Region 1') # add for plotting

# Data Munging ------------------------------------------------------------

# Munging ssb dataframe
ssb_sum_df <- ssb_df %>%
  filter((pd == T)) %>%
  rename(region = Var1, years = Var2) %>%
  mutate(re = (value - true) / true,
         region = paste("Region", region),
         model_id = ifelse(model_id == 1, "unstrc_mark", "ctmc")
         ) %>%
  group_by(model_id, region, years, scenario) %>%
  summarize(med = median(re),
            lwr_95 = quantile(re, 0.025),
            upr_95 = quantile(re, 0.975),
            lwr_75 = quantile(re, 0.125),
            upr_75 = quantile(re, 0.875)) %>%
  mutate(scenario = factor(scenario, levels = c("const", "age", "time", "timeage", "age_env")))

# Munging movement dataframe
move_summ_short_df <- move_df %>%
  filter(ages != 1, (pd == T)) %>%
  # Filter by scenario-specific conditions
  filter(
    (scenario == 'const' & ages == 2 & years == 1) | # filter to only one set of movement pars
      (scenario == 'age' & years == 1) | # filter to only ages
      (scenario == 'time' & ages == 2) | # filter to only years
      (scenario %in% c('timeage', "age_env") ) # filter to both ages and years
  ) %>%
  mutate( se = (value - true)^2 ) %>%
  mutate(
    from = paste("From region", from),
    to = paste("To region", to)
  ) %>%
  group_by(scenario, model_id) %>%
  summarize(
    rmse = sqrt(mean(se)),
    lwr_95_rmse = quantile(sqrt(se), 0.025),
    upr_95_rmse = quantile(sqrt(se), 0.975),
    lwr_75_rmse = quantile(sqrt(se), 0.125),
    upr_75_rmse = quantile(sqrt(se), 0.875),
    .groups = 'drop'
  ) %>%
  mutate(model_id = ifelse(model_id == 1, "unstrc_mark", "ctmc"),
         scenario = factor(scenario, levels = c("const", "age", "time", "timeage", "age_env"))) %>%
  left_join(conv_rates_plot, by = c('model_id', "scenario"))

# movement filtered to converged
move_conv_df <- move_df %>%
  filter(ages != 1, (pd == T)) %>%
  # Filter by scenario-specific conditions
  filter(
    (scenario == 'const' & ages == 2 & years == 1) | # filter to only one set of movement pars
      (scenario == 'age' & years == 1) | # filter to only ages
      (scenario == 'time' & ages == 2) | # filter to only years
      (scenario %in% c('timeage', "age_env") ) # filter to both ages and years
  )

# Get summary in addition to individual runs
move_conv_idv_df <- move_conv_df %>%
  filter(from == 3, to == 3) %>%
  group_by(scenario, ages, years, model_id) %>%
  mutate(median = median(value),
         lwr_95 = quantile((value), 0.025),
         upr_95 = quantile((value), 0.975),
         lwr_75 = quantile((value), 0.125),
         upr_75 = quantile((value), 0.875)) %>%
  mutate(model_id = ifelse(model_id == 1, "unstrc_mark", "ctmc"))


# SSB Plots -------------------------------------------------------------------
model_colors <- c("#C73E3A", "#2A9D8F")

ssb_re_plot <- ssb_sum_df %>%
  ggplot(aes(x = years, y = med, fill = factor(model_id), color = factor(model_id), lty = factor(model_id))) +
  geom_line() +
  geom_text(conv_rates_plot %>% filter(model_id == 'unstrc_mark') %>% mutate(scenario = factor(scenario, levels = c("const", "age", "time", "timeage", "age_env"))),
            mapping = aes(x = 3, y = -0.9, label = paste(num_pd, "%", sep ='')),
            show.legend = FALSE) +
  geom_text(conv_rates_plot %>% filter(model_id == 'ctmc') %>% mutate(scenario = factor(scenario, levels = c("const", "age", "time", "timeage", "age_env"))),
            mapping = aes(x = 3, y = -0.7, label = paste(num_pd, "%", sep ='')),
            show.legend = FALSE) +
  geom_ribbon(aes(x = years, y = med, ymin = lwr_75, ymax = upr_75, fill = factor(model_id)), color = NA, alpha = 0.25) +
  geom_ribbon(aes(x = years, y = med, ymin = lwr_95, ymax = upr_95, fill = factor(model_id)), color = NA, alpha = 0.2) +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, col = 'black', lty = 2, lwd = 1.3) +
  facet_grid(scenario~region) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.085, 0.95),
        legend.background = element_blank()) +
  labs(x = 'Year', y = 'Relative Error in Spawning Stock Biomass', fill = 'Model', color = 'Model', lty = 'Model')

ggsave(
  here("dev", "spatial_movement_ctmc", "ssb_re.png"),
  ssb_re_plot,
  width = 14, height = 15
)


# Movement Plots ----------------------------------------------------------
# RMSE for Movement
move_rmse_plot <- ggplot() +
  geom_pointrange(move_summ_short_df, mapping = aes(x = scenario, y = rmse, ymin = lwr_95_rmse, ymax = upr_95_rmse, color = model_id, alpha = num_pd),
                  position = position_dodge(width = 0.45), linewidth = 1) +
  geom_pointrange(move_summ_short_df, mapping = aes(x = scenario, y = rmse, ymin = lwr_75_rmse, ymax = upr_75_rmse, color = model_id, alpha = num_pd),
                  position = position_dodge(width = 0.45), linewidth = 3, size = 2) +
  geom_errorbar(move_summ_short_df, mapping = aes(x = scenario, ymin = lwr_95_rmse, ymax = upr_95_rmse, color = model_id, alpha = num_pd),
                position = position_dodge(width = 0.45), width = 0.25, linewidth = 1) +
  scale_color_manual(values = model_colors) +
  scale_alpha_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  labs(x = 'Scenario', y = 'Root Mean Squared Error in Movement Estimates', color = 'Model', alpha = "Convergence (%)") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.11, 0.725),
        legend.background = element_blank(),
        legend.key.height = unit(1, "cm"))

# Constant movement
const_est_plot <- ggplot() +
  geom_jitter(
    move_conv_idv_df %>% filter(scenario == 'const', years == 1, ages == 2),
    mapping = aes(x = model_id, y = value, color = model_id), position = position_jitter(width = 0.05, height = 0.05),
    alpha = 0.075
  ) +
  geom_point(
    move_conv_idv_df %>% filter(sim_id == 1, scenario == 'const', years == 1, ages == 2),
    mapping = aes(x = model_id, y = median),
    shape = 18, color = 'black', size = 5, stroke = 1
  ) +
  scale_color_manual(values = model_colors) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = 'Model', y = 'Movement', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none',
        legend.background = element_blank())


# Age varying movement
age_est_plot <- ggplot() +
  # Median line
  geom_line(move_conv_idv_df %>% filter(sim_id == 1, scenario == 'age', years == 1),
            mapping = aes(x = ages, y = median, color = factor(model_id), group = interaction(factor(model_id), sim_id)), lty = 2, lwd = 1.1) +
  # Individual simulations line
  geom_line(move_conv_idv_df %>% filter(sim_id %in% seq(1,100,5), scenario == 'age', years == 1), lwd = 0.5,
            mapping = aes(x = ages, y = value, color = factor(model_id), group = interaction(factor(model_id), sim_id)), alpha = 0.15) +
  # Shading 95%
  geom_ribbon(move_conv_idv_df %>% filter(sim_id == 1, scenario == 'age', years == 1),
              mapping = aes(x = ages, y = value, ymin = lwr_95, ymax = upr_95, fill = factor(model_id)), color = NA, alpha = 0.15) +
  # True
  geom_line(move_conv_df %>% filter(scenario == 'age', from == 3, to == 3, years == 1) ,
            mapping = aes(x = ages, y = true), lty = 2, color = 'black', lwd = 1.1) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(x = 'Ages', y = 'Movement', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.16, 0.15),
        legend.background = element_blank()
  )


# Time varying movement
time_est_plot <- ggplot() +
  # Median line
  geom_line(move_conv_idv_df %>% filter(sim_id == 1, scenario == 'time', ages == 2),
            mapping = aes(x = years, y = median, color = factor(model_id), group = interaction(factor(model_id), sim_id)), lty = 2, lwd = 1.1) +
  # Individual simulations line
  geom_line(move_conv_idv_df %>% filter(sim_id %in% seq(1,100,5), scenario == 'time', ages == 2), lwd = 0.5,
            mapping = aes(x = years, y = value, color = factor(model_id), group = interaction(factor(model_id), sim_id)), alpha = 0.15) +
  # Shading 95%
  geom_ribbon(move_conv_idv_df %>% filter(sim_id == 1, scenario == 'time', ages == 2),
              mapping = aes(x = years, y = value, ymin = lwr_95, ymax = upr_95, fill = factor(model_id)), color = NA, alpha = 0.15) +
  # True
  geom_line(move_conv_df %>% filter(scenario == 'time', from == 3, to == 3, ages == 2) ,
            mapping = aes(x = years, y = true), lty = 2, color = 'black', lwd = 1.1) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(x = 'Years', y = 'Movement', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.16, 0.15),
        legend.background = element_blank()
        )

timeage_est_plot <- ggplot() +
  geom_tile(move_conv_idv_df %>% filter(sim_id == 97, scenario == 'timeage', ages != 1) %>% mutate(model_id = 'true'),
            mapping = aes(x = years, y = ages, fill = true), alpha = 0.75) +
  geom_tile(move_conv_idv_df %>% filter(sim_id == 97, scenario == 'timeage', ages != 1),
            mapping = aes(x = years, y = ages, fill = median), alpha = 0.75) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0.5, limits = c(0,1)) +
  labs(x = 'Years', y = 'Ages', fill = 'Movement') +
  facet_wrap(~factor(model_id, levels = c("true", 'ctmc', 'unstrc_mark')), ncol = 1) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.2, 0.125),
        legend.background = element_blank()
        )

age_env_est_plot <- ggplot() +
  geom_tile(move_conv_idv_df %>% filter(sim_id == 10, scenario == 'age_env', ages != 1) %>% mutate(model_id = 'true'),
            mapping = aes(x = years, y = ages, fill = true), alpha = 0.75) +
  geom_tile(move_conv_idv_df %>% filter(sim_id == 10, scenario == 'age_env', ages != 1),
            mapping = aes(x = years, y = ages, fill = median), alpha = 0.75) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0.5, limits = c(0,1)) +
  labs(x = 'Years', y = 'Ages', fill = 'Movement') +
  facet_wrap(~factor(model_id, levels = c("true", 'ctmc', 'unstrc_mark')), ncol = 1) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.2, 0.125),
        legend.background = element_blank()
  )

# combine plots
one <- cowplot::plot_grid(
  const_est_plot, age_est_plot, time_est_plot, nrow = 1,
  rel_widths = c(0.4, 0.6, 0.6),
  labels = c("B", 'C', 'D'),
  label_size = 30, hjust = -0.3
)

two <- cowplot::plot_grid(
  timeage_est_plot, age_env_est_plot,
  labels = c('E', 'F'),
  label_size = 30
)

three <- cowplot::plot_grid(
  move_rmse_plot, two,
  rel_widths = c(0.5,0.5),
  labels = c('A', ''),
  label_size = 30
)

four <- cowplot::plot_grid(
  three, one, ncol = 1,
  rel_heights = c(0.7, 0.3)
)

ggsave(
  here("dev", "spatial_movement_ctmc", "move_sim_error.png"),
  four,
  width = 25, height = 16
)

