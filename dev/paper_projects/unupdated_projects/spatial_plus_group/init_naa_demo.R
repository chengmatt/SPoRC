# Purpose: To set up different demographic scenarios and evaluate implications
# of alterantive methods to initialize population
# Creator: Matthew LH. Cheng
# Date Created: 10/24/25


# setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)
library(cowplot)
library(patchwork)
devtools::load_all(here("R"))


# setup demographic scenarios ---------------------------------------------

# dimensions
n_regions = 2
n_ages = 15
n_sexes = 1
n_fish_fleets = 1
n_seas = 1

# variables held constant
global_R0 = 1
init_F = 0
fish_sel = array(1, dim = c(n_regions, n_ages, n_sexes, n_fish_fleets))
sex_ratio = array(1, dim = c(n_regions, n_sexes))
do_recruits_move = 0
ln_InitDevs = array(0, dim = c(n_regions, n_ages - 1))

# initial age structure scenarios
init_age_strc <- c(0,1,2,3) # ext_ages, scalar w/ no movement, matrix, scalar with movement

# Movement scenarios
# No movement
# Random movement
# One way movement
# ontogenetic movement
move_list = list(

  # no movement
  no = array(diag(1, n_regions), dim = c(n_regions, n_regions, n_seas,  n_ages, n_sexes)),

  # random movement
  rand = array(1/n_regions, dim = c(n_regions, n_regions, n_seas, n_ages, n_sexes)),

  # one way movement
  oneway = {
    rate = 0.3
    arr = array(0, dim = c(n_regions, n_regions, n_seas, n_ages, n_sexes))
    arr[1,1,,,] = rate
    arr[1,2,,,] = 1 - rate
    arr[2,2,,,] = 1 - rate
    arr[2,1,,,] = rate
    arr
  },

  # ontogenetic movement
  onto = {

    # logistic parameters
    midpoint = n_ages / 2
    steepness = 0.1
    max_leave_rate = 0.8

    # get array
    arr = array(0, dim = c(n_regions, n_regions, n_seas, n_ages, n_sexes))

    # get age-specific leaving rates using logistic
    ages = 1:n_ages
    leave_rates = max_leave_rate / (1 + exp(-steepness * (ages - midpoint)))

    # get age-specific movement from region 1 to region 2
    for(age in 1:n_ages) {
      arr[1, 1, , age, ] = 1 - leave_rates[age]  # stay in region 1
      arr[1, 2, , age, ] = leave_rates[age]      # move to region 2
      # fish in region 2 stay put
      arr[2, 1, , age, ] = 1 - max_leave_rate
      arr[2, 2, , age, ] = max_leave_rate
    }
    arr
  }

)

# Natural Mortality scenarios
# high, medium and low
natmort_list = list(
  eq = array(0.25, dim = c(n_regions, n_ages, n_sexes)),
  uneq = array(c(0.2, 0.25), dim = c(n_regions, n_ages, n_sexes))
)

# regional R0
# equal, source-area
r0r_list = list(
  eq = global_R0 * c(0.5, 0.5),
  src = global_R0 * c(0.7, 0.3)
)

# factorial design combining all scenarios
factorial_scenarios = expand.grid(
  movement = names(move_list),
  natmort = names(natmort_list),
  r0r = names(r0r_list),
  stringsAsFactors = FALSE
)

# loop through to populate scenarios
scenario_list = vector("list", nrow(factorial_scenarios))
for(i in 1:nrow(factorial_scenarios)) {
  scenario_list[[i]] = list(
    name = paste(factorial_scenarios$movement[i],
                 factorial_scenarios$natmort[i],
                 factorial_scenarios$r0r[i],
                 sep = "_"),
    movement = move_list[[factorial_scenarios$movement[i]]],
    natmort = natmort_list[[factorial_scenarios$natmort[i]]],
    r0r = r0r_list[[factorial_scenarios$r0r[i]]]
  )
}
names(scenario_list) = sapply(scenario_list, function(x) x$name) # name



# Run Scenarios -----------------------------------------------------------
init_naa_results = data.frame() # storage

for(i in 1:length(scenario_list)) {

  # get scenario
  tmp_scenario = scenario_list[[i]]

  for(j in 1:length(init_age_strc)) {
    tmp_naa = Get_Init_NAA(
      init_age_strc = init_age_strc[j], # initial age structure
      init_iter = n_ages * 5,
      n_regions = n_regions,
      n_sexes = n_sexes,
      n_ages = n_ages,
      natmort = tmp_scenario$natmort,
      init_F = 0,
      fish_sel = fish_sel,
      R0_r = tmp_scenario$r0r,
      sexratio = sex_ratio,
      Movement = tmp_scenario$movement,
      do_recruits_move = do_recruits_move,
      ln_InitDevs = ln_InitDevs,
      n_seas = 1, seasdur = 1
    )

    # bind results
    tmp_naa_df = reshape2::melt(tmp_naa)
    tmp_naa_df$scenario = tmp_scenario$name
    tmp_naa_df$init_age_strc = init_age_strc[j]
    init_naa_results = rbind(tmp_naa_df, init_naa_results)
  } # end j loop

} # end i loop

# rename
init_naa_results = init_naa_results %>% rename(Region = Var1, Age = Var2, Sex = Var3)

# parse out scenario
init_naa_results <- init_naa_results %>%
  separate(scenario, into = c("movement", "natmort", "r0r"),  sep = "_", remove = FALSE)

# rename init_age_strc
init_naa_results = init_naa_results %>%
  mutate(init_age_strc =
           case_when(
             init_age_strc == 0 ~ "ext_ages",
             init_age_strc == 1 ~ "no_move_all",
             init_age_strc == 2 ~ "matrix",
             init_age_strc == 3 ~ "no_move_plus"
           ),
         init_age_strc = factor(init_age_strc, levels = c("no_move_all", "no_move_plus",
                                                          "ext_ages", "matrix")),
         movement = factor(movement, levels = c("no", "rand", "oneway", "onto"),
                           labels = paste("move:", c("no", "rand", "oneway", "onto"))),
         natmort = factor(natmort, levels = c("eq", 'uneq', "higheq", 'highuneq',
                                              "eqage", 'uneqage', "higheqage", 'highuneqage'),
                          labels = paste("natmort:", c("eq", 'uneq', "higheq", 'highuneq',
                                                       "eqage", 'uneqage', "higheqage", 'highuneqage'))))

# compute relative difference compared to matrix method
init_naa_results <- init_naa_results %>%
  group_by(scenario, Region, Age, Sex) %>%
  mutate(
    matrix_value = value[init_age_strc == "matrix"],
    rel_diff = (value - matrix_value) / matrix_value,
  ) %>%
  ungroup()

# Plot --------------------------------------------------------------------

### Region 1 ----------------------------------------------------------------

# Main plot
main <- init_naa_results %>%
  filter(Region == 1, r0r == 'src') %>%
  ggplot(aes(x = Age, y = value, color = init_age_strc, lty = init_age_strc)) +
  geom_line(lwd = 1.3) +
  ggh4x::facet_grid2(movement ~ natmort, scales = 'free_y', independent = 'y') +
  labs(lty = 'Method', color = 'Method', lty = "Method", y = 'Equilibrium Abundance (Region 1)') +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none')

# First panel (move: no, natmort: eq) - top left
first <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: no', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Second panel (move: no, natmort: uneq) - top right
second <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: no', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Third panel (move: rand, natmort: eq)
third <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: rand', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Fourth panel (move: rand, natmort: uneq)
fourth <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: rand', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Fifth panel (move: oneway, natmort: eq)
fifth <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: oneway', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Sixth panel (move: oneway, natmort: uneq)
sixth <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: oneway', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Seventh panel (move: onto, natmort: eq)
seventh <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: onto', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Eighth panel (move: onto, natmort: uneq)
eighth <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: onto', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Combine all insets
region_1 <- main +
  inset_element(first, left = 0.28, bottom = 0.85, right = 0.43, top = 0.99) +    # First panel
  inset_element(second, left = 0.81, bottom = 0.85, right = 0.96, top = 0.99) +   # Top-right panel
  inset_element(third, left = 0.28, bottom = 0.59, right = 0.43, top = 0.73) +    # Second row left
  inset_element(fourth, left = 0.81, bottom = 0.59, right = 0.96, top = 0.73) +   # Second row right
  inset_element(fifth, left = 0.28, bottom = 0.33, right = 0.43, top = 0.47) +    # Third row left
  inset_element(sixth, left = 0.81, bottom = 0.33, right = 0.96, top = 0.47) +    # Third row right
  inset_element(seventh, left = 0.28, bottom = 0.08, right = 0.43, top = 0.22) +  # Bottom row left
  inset_element(eighth, left = 0.81, bottom = 0.08, right = 0.96, top = 0.22)     # Bottom row right



### Region 2 ----------------------------------------------------------------

# Main plot
main <- init_naa_results %>%
  filter(Region == 2, r0r == 'src') %>%
  ggplot(aes(x = Age, y = value, color = init_age_strc, lty = init_age_strc)) +
  geom_line(lwd = 1.3) +
  ggh4x::facet_grid2(movement ~ natmort, scales = 'free_y', independent = 'y') +
  labs(lty = 'Method', color = 'Method', lty = "Method", y = 'Equilibrium Abundance (Region 2)') +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none')

# First panel (move: no, natmort: eq) - top left
first <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: no', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Second panel (move: no, natmort: uneq) - top right
second <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: no', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Third panel (move: rand, natmort: eq)
third <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: rand', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Fourth panel (move: rand, natmort: uneq)
fourth <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: rand', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Fifth panel (move: oneway, natmort: eq)
fifth <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: oneway', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Sixth panel (move: oneway, natmort: uneq)
sixth <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: oneway', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Seventh panel (move: onto, natmort: eq)
seventh <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: onto', natmort == 'natmort: eq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Eighth panel (move: onto, natmort: uneq)
eighth <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: onto', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Combine all insets
region_2 <- main +
  inset_element(first, left = 0.28, bottom = 0.85, right = 0.43, top = 0.99) +    # First panel
  inset_element(second, left = 0.81, bottom = 0.85, right = 0.96, top = 0.99) +   # Top-right panel
  inset_element(third, left = 0.28, bottom = 0.59, right = 0.43, top = 0.73) +    # Second row left
  inset_element(fourth, left = 0.81, bottom = 0.59, right = 0.96, top = 0.73) +   # Second row right
  inset_element(fifth, left = 0.28, bottom = 0.33, right = 0.43, top = 0.47) +    # Third row left
  inset_element(sixth, left = 0.81, bottom = 0.33, right = 0.96, top = 0.47) +    # Third row right
  inset_element(seventh, left = 0.28, bottom = 0.08, right = 0.43, top = 0.22) +  # Bottom row left
  inset_element(eighth, left = 0.81, bottom = 0.08, right = 0.96, top = 0.22)     # Bottom row right

# get legend
legend <- cowplot::get_legend(main + theme(legend.position = 'top') + labs(fill = 'Method'))
comb <- cowplot::plot_grid(region_1, region_2, labels = c("A", "B"), label_size = 30)
comb2 <- cowplot::plot_grid(legend, comb, ncol = 1, rel_heights = c(0.05, 0.95))

ggsave(
  here("dev", "paper_projects", "spatial_plus_group", "plus_group_demo_plot.png"),
  comb2,
  width = 13, height = 10
)

# Condensed Plot ----------------------------------------------------------
main_condensed <- init_naa_results %>%
  filter(r0r == 'src', movement %in% c("move: oneway", "move: onto"), natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, color = init_age_strc, lty = init_age_strc)) +
  geom_line(lwd = 1.3) +
  ggh4x::facet_grid2(paste("Region", Region) ~ movement, scales = 'free_y', independent = 'y') +
  labs(lty = 'Method', color = 'Method', lty = "Method", y = 'Equilibrium Abundance') +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top')

# First panel (move: oneway, natmort: uneq)
one <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: oneway', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Second panel (move: onto, natmort: uneq)
two <- init_naa_results %>%
  filter(Region == 1, r0r == 'src', Age == 15,
         movement == 'move: onto', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Third panel (move: oneway, natmort: uneq)
three <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: oneway', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Fourth panel (move: onto, natmort: uneq)
four <- init_naa_results %>%
  filter(Region == 2, r0r == 'src', Age == 15,
         movement == 'move: onto', natmort == 'natmort: uneq') %>%
  ggplot(aes(x = Age, y = value, fill = init_age_strc)) +
  geom_col(position = position_dodge(), color = 'black', alpha = 0.85, lwd = 0.1) +
  ggthemes::scale_fill_colorblind() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())


# Combine all insets
condensed <- main_condensed + theme_bw(base_size = 20) + theme(legend.position = 'top') +
  inset_element(one, left = 0.28, bottom = 0.85, right = 0.43, top = 0.99) +
  inset_element(two, left = 0.81, bottom = 0.85, right = 0.96, top = 0.99) +
  inset_element(three, left = 0.28, bottom = 0.33, right = 0.43, top = 0.47) +
  inset_element(four, left = 0.81, bottom = 0.33, right = 0.96, top = 0.47)

ggsave(
  here("dev", "paper_projects", "spatial_plus_group", "plus_group_condensed_plot.png"),
  condensed,
  width = 10, height = 7
)


# Movement Scenarios ------------------------------------------------------
movement_scenarios <- reshape2::melt(move_list) %>%
  rename(from = Var1, to = Var2, seas = Var3, age = Var4, sex = Var5, scenario = L1) %>%
  mutate(scenario = paste("move:", scenario))

# No Movement
nomove <- ggplot(movement_scenarios %>%
                   filter(age == 1, sex == 1, scenario == 'move: no'),
                 aes(x = paste("Region", factor(from)), y = paste("Region", factor(to)), fill = value, label = value)) +
  geom_tile(alpha = 0.5) +
  geom_text(size = 10) +
  scale_fill_viridis_c() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top',
        legend.key.width = unit(1, 'cm')) +
  labs(x = "From", y = 'To', fill = 'Movement Probability')

# Random Movement
randmove <- ggplot(movement_scenarios %>%
                     filter(age == 1, sex == 1, scenario == 'move: rand'),
                   aes(x = paste("Region", factor(from)), y = paste("Region", factor(to)), fill = value, label = value)) +
  geom_tile(alpha = 0.5) +
  geom_text(size = 10) +
  scale_fill_viridis_c() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top',
        legend.key.width = unit(1, 'cm')) +
  labs(x = "From", y = 'To', fill = 'Movement Probability')

# One way Movement
onewaymove <- ggplot(movement_scenarios %>%
                   filter(age == 1, sex == 1, scenario == 'move: oneway'),
                 aes(x = paste("Region", factor(from)), y = paste("Region", factor(to)), fill = value, label = value)) +
  geom_tile(alpha = 0.5) +
  geom_text(size = 10) +
  scale_fill_viridis_c() +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top',
        legend.key.width = unit(1, 'cm')) +
  labs(x = "From", y = 'To', fill = 'Movement Probability')

# Onotgenetic Movement
ontomove <- ggplot(movement_scenarios %>%
                     filter(sex == 1, scenario == 'move: onto'),
                   aes(x = age, y = value,
                       color = interaction(paste("From R", factor(from)), paste("To R", factor(to)), sep = ' '))) +
  geom_line(lwd = 1.3) +
  scale_fill_viridis_c() +
  ylim(0,1) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top') +
  labs(x = "Age", y = 'Movement Probability', fill = 'Movement Probability', color = '')

ggsave(
  here("dev", "paper_projects", "spatial_plus_group", "movement_scenarios.png"),
  cowplot::plot_grid(
    nomove, randmove, onewaymove, ontomove, align = 'hv',
    labels = c('A', 'B', 'C', 'D'), label_size = 30, hjust = -3, vjust = 3
  ),
  width = 18.5, height = 13
)

