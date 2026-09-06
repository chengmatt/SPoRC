# Purpose: To demonstrate simulation movement scenarios and trajectories
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 10/30/25


# Setup -------------------------------------------------------------------

library(igraph)
library(SPoRC)
library(here)
library(tidyverse)

# read in simulations
const <- readRDS(here("dev", "dev_output", "CTMC", "const.RDS"))
age_move <- readRDS(here("dev", "dev_output", "CTMC", "age_move.RDS"))
time_move <- readRDS(here("dev", "dev_output", "CTMC", "time_move.RDS"))
timeage_move <- readRDS(here("dev", "dev_output", "CTMC", "timeage_move.RDS"))
age_env_move <- readRDS(here("dev", "dev_output", "CTMC", "age_env_move.RDS"))


# SSB plot ----------------------------------------------------------------
model_colors <- c("#C73E3A", "#2A9D8F", "#E9B44C", "#9B59B6")

ssb_sim <- reshape2::melt(const$SSB) %>% mutate(om = "const") %>%
  bind_rows(
    reshape2::melt(age_move$SSB) %>% mutate(om = "age"),
    reshape2::melt(time_move$SSB) %>% mutate(om = "time"),
    reshape2::melt(timeage_move$SSB) %>% mutate(om = "timeage"),
    reshape2::melt(age_env_move$SSB) %>% mutate(om = "age_env")
  ) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3) %>%
  mutate(Region = paste("Region", Region),
         om = factor(om, levels = c('const', 'age', 'time', 'timeage', "age_env")))

# get SSB summary
ssb_sim_sum <- ssb_sim %>%
  group_by(Region, Year, om) %>%
  summarize(
    median = median(value),
    lwr_95 = quantile(value, 0.025),
    upr_95 = quantile(value, 0.975),
    lwr_75 = quantile(value, 0.125),
    upr_75 = quantile(value, 0.875)
  )

ssb_plot <- ssb_sim_sum %>%
  ggplot(aes(x = Year, y = median, ymin = lwr_95, ymax = upr_95)) +
  geom_line(aes(x = Year, y = median), lwd = 0.9 ) +
  geom_ribbon(aes(x = Year, y = median, ymin = lwr_95, ymax = upr_95), color = NA, alpha = 0.25) +
  geom_ribbon(aes(x = Year, y = median, ymin = lwr_75, ymax = upr_75), color = NA, alpha = 0.3) +
  facet_grid(om~Region) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none') +
  labs(x = 'Year', y = 'Spawning Stock Biomass')

ggsave(
  here("dev", "spatial_movement_ctmc", "ssb_plot_sim.png"),
  ssb_plot,
  height = 10, width = 7.5
)


# Movement plot ---------------------------------------------------------------
colors <- c("#C73E3A", "#2A9D8F", "#E9B44C", "#9B59B6")

# constant movement
constmove_df <- reshape2::melt(const$Movement) %>%
  rename(from = Var1, to = Var2, year = Var3, age = Var4, sex = Var5, sim = Var6) %>%
  mutate(from = paste("Region", from),
         to = paste("Region", to))

const_move_plot <- constmove_df %>%
  filter(sex == 1, year == 1, sim == 1) %>%
  ggplot(aes(x = age, y = value, color = factor(to))) +
  geom_line(lwd = 1.3) +
  facet_wrap(~from) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 20) +
  labs(x = 'Age', y = 'Movement', color = '', lty = '') +
  theme(legend.position = c(0.08, 0.5), legend.background = element_blank())

# Time plots
timemove_df <- reshape2::melt(time_move$Movement) %>%
  rename(from = Var1, to = Var2, year = Var3, age = Var4, sex = Var5, sim = Var6) %>%
  mutate(from = paste("From Region", from), to = paste("To Region", to))

timemove_plot <- timemove_df %>%
  filter(sex == 1, age == 15, sim == 1) %>%
  ggplot(aes(x = year, y = value, color = factor(to))) +
  geom_line(lwd = 1.3) +
  facet_wrap(~from) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 20) +
  labs(x = 'Year', y = 'Movement', color = '', lty = '') +
  theme(legend.position = c(0.1, 0.875), legend.background = element_blank())

# Age plots
agemove_df <- reshape2::melt(age_move$Movement) %>%
  rename(from = Var1, to = Var2, year = Var3, age = Var4, sex = Var5, sim = Var6) %>%
  mutate(from = paste("From Region", from), to = paste("To Region", to))

agemove_plot <- agemove_df %>%
  filter(sex == 1, age != 1, sim == 1, year == 30) %>%
  ggplot(aes(x = age, y = value, color = factor(to))) +
  geom_line(lwd = 1.3) +
  facet_wrap(~from) +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw(base_size = 20) +
  labs(x = 'Age', y = 'Movement', color = '', lty = '') +
  theme(legend.position = c(0.1, 0.875), legend.background = element_blank())

# time age plots
timeage_move_df <- reshape2::melt(timeage_move$Movement) %>%
  rename(from = Var1, to = Var2, year = Var3, age = Var4, sex = Var5, sim = Var6) %>%
  mutate(from = paste("From Region", from), to = paste("To Region", to))

timeage_plot <- timeage_move_df %>%
  filter(sim == 1) %>%
  ggplot(aes(x = year, y = age, z = value)) +
  geom_raster(aes(fill = value), interpolate = F, alpha = 0.45) +
  geom_contour(color = "black", alpha = 0.5, bins = 8, linewidth = 1) +
  facet_grid(to~from) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0.5) +
  theme_bw(base_size = 20) +
  labs(x = 'Year', y = 'Age', fill = 'Movement') +
  theme(legend.position = c(0.1, 0.875), legend.background = element_blank())

gaussian_dome <- function(temp, optimal_temp, width) {
  exp(-((temp - optimal_temp)^2) / (2 * width^2))
}

# Environmental plots
age_envmove_df <- reshape2::melt(age_env_move$Movement) %>%
  rename(from = Var1, to = Var2, year = Var3, age = Var4, sex = Var5, sim = Var6) %>%
  mutate(from = paste("From Region", from), to = paste("To Region", to))

age_env_plot <- age_envmove_df %>%
  filter(sex == 1, age == 15, sim == 1) %>%
  ggplot(aes(x = year, y = value, color = factor(to))) +
  geom_line(lwd = 1.3) +
  facet_wrap(~from) +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 20) +
  labs(x = 'Year', y = 'Movement', color = '', lty = '') +
  theme(legend.position = c(0.1, 0.875), legend.background = element_blank())

temp_pref <- data.frame(
  temp = seq(-3, 3, 0.01), pref = gaussian_dome(seq(-3, 3, 0.01), 0, 1)
)

pref_plot <- ggplot(temp_pref, aes(x = temp, y = pref * 0.85)) +
  geom_line(lwd = 1.3) +
  theme_bw(base_size = 20) +
  labs(x = 'Covariate Value', y = 'Preference', color = '', lty = '') +
  theme(legend.position = c(0.1, 0.875), legend.background = element_blank())

temp_effect <- as.matrix(read.csv(here("dev", "dev_output", "CTMC", "age_env_temp.csv")))
temp_df <- reshape2::melt(temp_effect) %>%
  rename(regions = Var1, years = Var2, temp = value) %>%
  filter(years != "X") %>%
  mutate(years = as.numeric(str_remove(years, 'V')),
         taxis = gaussian_dome(temp, 0, 1))

cov_plot <- ggplot(temp_df, aes(x = years, y = taxis, color = paste("Region", factor(regions)))) +
  geom_line(lwd = 1.3) +
  scale_color_manual(values = colors)  +
  theme_bw(base_size = 20) +
  labs(x = 'Year', y = 'Environmental Taxis (Preference)', color = '', lty = '') +
  theme(legend.position = c(0.15, 0.89), legend.background = element_blank())

# Time age varying movement plot
timeage_move_plots <- cowplot::plot_grid(
  const_move_plot,
  agemove_plot,
  timemove_plot,
  age_env_plot,
  ncol = 1,
  labels = c("A", "B", "C", 'D'),
  label_size = 30,
  hjust = -1.5
)


# combine with constant movement plot
move_plot <- cowplot::plot_grid(

  cov_plot,
  labels = c('', 'E'),
  label_size = 30,
  rel_widths = c(0.55, 0.45)
)

ggsave(
  here("dev", "spatial_movement_ctmc", "move_plot_sim.png"),
  move_plot,
  height = 13, width = 20
)


# Demographics ------------------------------------------------------------
colors <- c("#C73E3A", "#2A9D8F", "#E9B44C", "#9B59B6")

# natural mortality
natmort_plot <- reshape2::melt(
  const$natmort
) %>%
  ggplot(aes(x = Var3, y = value)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Age', y = 'Instantaneous Natural Mortality Rate') +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')

# weight at age
waa_plot <- reshape2::melt(
  const$WAA
) %>%
  ggplot(aes(x = Var3, y = value)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Age', y = 'Weight') +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')

# fishery selex
fishsel_plot <- reshape2::melt(
  const$fish_sel
) %>%
  ggplot(aes(x = Var3, y = value)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Age', y = 'Fishery Selectivity') +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')

# survey selex
srvsel_plot <- reshape2::melt(
  const$srv_sel
) %>%
  ggplot(aes(x = Var3, y = value)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Age', y = 'Survey Selectivity') +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')


ggsave(
  here("dev", "spatial_movement_ctmc", "demographic_sim.png"),
  cowplot::plot_grid(
    natmort_plot, waa_plot, fishsel_plot, srvsel_plot,
    labels = c("A", 'B', 'C', 'D'), label_size = 20
  ),
  height = 13, width = 15,
)

# Recruitment -------------------------------------------------------------

rec_plot <- reshape2::melt(
  const$Rec
) %>% mutate(type = 'const') %>%
  bind_rows(
    reshape2::melt(
      age_move$Rec
    ) %>% mutate(type = 'age'),
    reshape2::melt(
      time_move$Rec
    ) %>% mutate(type = 'time'),
    reshape2::melt(
      timeage_move$Rec
    ) %>% mutate(type = 'timeage'),
    reshape2::melt(
      age_env_move$Rec
    ) %>% mutate(type = 'age_env')
  ) %>%
  mutate(type = factor(type, levels = c('const', 'age', 'time', 'timeage', 'age_env'))) %>%
  ggplot(aes(x = Var2, y = value, group = Var3)) +
  geom_line(alpha = 0.25) +
  facet_grid(type~paste('Region', Var1)) +
  labs(x = 'Year', y = 'Recruitment') +
  theme_bw(base_size = 15)

ggsave(
  here("dev", "spatial_movement_ctmc", "rec_sim.png"),
  rec_plot,
  height = 13, width = 15,
)

# Fishing Mortality -------------------------------------------------------

f_plot <- reshape2::melt(
  const$Fmort
) %>% mutate(type = 'const') %>%
  bind_rows(
    reshape2::melt(
      age_move$Fmort
    ) %>% mutate(type = 'age'),
    reshape2::melt(
      time_move$Fmort
    ) %>% mutate(type = 'time'),
    reshape2::melt(
      timeage_move$Fmort
    ) %>% mutate(type = 'timeage'),
    reshape2::melt(
      age_env_move$Fmort
    ) %>% mutate(type = 'age_env')
  ) %>%
  mutate(type = factor(type, levels = c('const', 'age', 'time', 'timeage', 'age_env'))) %>%
  ggplot(aes(x = Var2, y = value, group = Var4)) +
  geom_line(alpha = 0.25) +
  facet_grid(type~paste('Region', Var1)) +
  labs(x = 'Year', y = 'Instantaneous Fishing Mortality Rate') +
  theme_bw(base_size = 15)

ggsave(
  here("dev", "spatial_movement_ctmc", "f_sim.png"),
  f_plot,
  height = 13, width = 15,
)

