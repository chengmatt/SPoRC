# Purpose: To demonstrate the CTMC model conceptually as a 1d chain (only adjacent regions are connected)
# Creator: Matthew LH. Cheng
# Date 10/30/25


# Setup -------------------------------------------------------------------

library(igraph)
library(SPoRC)
library(tidyverse)
library(here)


# CTMC Model --------------------------------------------------------------
create_1d_chain <- function(
  n_cells = 25,
  n_ages = 5,
  diffusion = 3,
  age_pref_coeff,
  t_interval = c(1, 3, 5),
  grid_start = c(1, 3, 25)
) {

  # create containers
  diffusion_mat <- adjacency <- array(0, c(n_cells, n_cells))
  movement <- m <- taxis <- array(0, dim = c(n_cells, n_cells, n_ages))
  age_pref <- pref <- array(0, dim = c(n_cells, n_ages))
  movement_t <- array(0, dim = c(n_cells, n_cells, n_ages, length(t_interval))) # movement by time interval
  distribution <- array(0, dim = c(n_cells, n_ages, length(t_interval) + 1, length(grid_start)))

  # Create spatial cells
  regions <- 1:n_cells

  # create 1d adjacency matrix w/ reflective boundary
  for(i in 1:n_cells) {
    if(i > 1) adjacency[i-1, i] <- 1
    if(i < n_cells) adjacency[i+1, i] <- 1
  }

  # Get diffusion (constant by age)
  diffusion_mat[adjacency == 1] <- diffusion
  diag(diffusion_mat) <- -colSums(diffusion_mat) # make cols sum to 0 to conserve abundance

  # Get depth preference
  depth <- 1 + 5 * sin(pi * (1:n_cells) / (n_cells + 1))

  for(a in 1:n_ages) {

    # get age preference
    age_pref[,a] <- age_pref_coeff[a] * 1:n_cells

    # get total preference
    pref[,a] <- depth + age_pref[,a]

    # get taxis
    for(i in 1:n_cells) {
      for(j in 1:n_cells) {
        if(adjacency[i,j] == 1) taxis[i,j,a] <- pref[i,a] - pref[j,a] # get linear preference
      } # end i
    } # end j
    diag(taxis[,,a]) <- -colSums(taxis[,,a])

    # get instantaneous movement rates
    m[,,a] <- diffusion_mat + taxis[,,a]

    # get movement fractions
    movement[,,a] <- as.matrix(Matrix::expm(m[,,a]))

    # get movement across time
    for(t in 1:length(t_interval)) {
      movement_t[,,a,t] <- as.matrix(Matrix::expm(t_interval[t] * m[,,a])) # movement by time

      # get distribution after a given time interval
      for(z in 1:length(grid_start)) {
        tmp_dist <- rep(0, n_cells)
        tmp_dist[grid_start[z]] <- 1 # where the cohort starts
        distribution[,a,t,z] <- movement_t[,,a,t] %*% tmp_dist
      } # end z loop
    } # end t loop

    # get stationary movement
    eigen_vals <- eigen(movement[,,a])
    for(z in 1:length(grid_start)) {
      distribution[,a, length(t_interval) + 1, z] <- Re(eigen_vals$vectors[,1]) / sum(Re(eigen_vals$vectors[,1]))
    } # end z loop

  } # end age loop

  # munge to dataframes
  pref_df <- reshape2::melt(pref) %>% rename(cell = Var1, age = Var2)
  age_pref_df <- reshape2::melt(age_pref) %>% rename(cell = Var1, age = Var2)
  depth_df <- data.frame(cell = 1:n_cells, value = depth)
  diffusion_df <- reshape2::melt(diffusion_mat) %>% rename(to = Var1, from = Var2)
  taxis_df <- reshape2::melt(taxis) %>% rename(to = Var1, from = Var2, age = Var3)
  m_df <- reshape2::melt(m) %>% rename(to = Var1, from = Var2, age = Var3)
  movement_df <- reshape2::melt(movement) %>% rename(to = Var1, from = Var2, age = Var3)
  distribution_df <- reshape2::melt(distribution) %>% rename(cell = Var1, age = Var2, time = Var3, grid = Var4)

  return(list(
    total_pref = pref_df,
    age_pref = age_pref_df,
    depth_pref = depth_df,
    diffusion = diffusion_df,
    taxis = taxis_df,
    m_rate = m_df,
    distribution = distribution_df,
    movement = movement_df
  ))

} # end function

create_1d_chain <- function(n_cells = 25,
                            n_ages = 5,
                            diffusion = 3,
                            age_pref_coeff,
                            t_interval = c(1, 3, 5),
                            grid_start = c(1, 3, 25)) {
  # create containers
  diffusion_mat <- adjacency <- array(0, c(n_cells, n_cells))
  movement <- m <- taxis <- array(0, dim = c(n_cells, n_cells, n_ages))
  age_pref <- pref <- array(0, dim = c(n_cells, n_ages))
  movement_t <- array(0, dim = c(n_cells, n_cells, n_ages, length(t_interval)))
  distribution <- array(0, dim = c(n_cells, n_ages, length(t_interval) + 1, length(grid_start)))

  regions <- 1:n_cells

  # create 1d adjacency matrix w/ reflective boundary
  # Convention: adjacency[i,j] = 1 means you can move FROM i TO j
  for(i in 1:n_cells) {
    if(i > 1) adjacency[i, i-1] <- 1      # from i to i-1
    if(i < n_cells) adjacency[i, i+1] <- 1 # from i to i+1
  }

  # Get diffusion (constant by age)
  diffusion_mat[adjacency == 1] <- diffusion
  diag(diffusion_mat) <- -rowSums(diffusion_mat) # rows sum to 0

  # Get depth preference
  depth <- 1 + 5 * sin(pi * (1:n_cells) / (n_cells + 1))

  for(a in 1:n_ages) {
    age_pref[,a] <- age_pref_coeff[a] * 1:n_cells
    pref[,a] <- depth + age_pref[,a]

    # get taxis
    # taxis[i,j] = preference-driven rate from i to j
    for(i in 1:n_cells) {
      for(j in 1:n_cells) {
        if(adjacency[i,j] == 1) {
          # Move from i to j based on preference difference
          taxis[i,j,a] <- pref[j,a] - pref[i,a]
        }
      }
    }
    diag(taxis[,,a]) <- -rowSums(taxis[,,a]) # rows sum to 0

    m[,,a] <- diffusion_mat + taxis[,,a]
    movement[,,a] <- as.matrix(Matrix::expm(m[,,a]))

    for(t in 1:length(t_interval)) {
      movement_t[,,a,t] <- as.matrix(Matrix::expm(t_interval[t] * m[,,a]))
      for(z in 1:length(grid_start)) {
        tmp_dist <- rep(0, n_cells)
        tmp_dist[grid_start[z]] <- 1
        distribution[,a,t,z] <- t(tmp_dist) %*% movement_t[,,a,t]
      }
    }

    eigen_vals <- eigen(t(movement[,,a]))
    for(z in 1:length(grid_start)) {
      distribution[,a, length(t_interval) + 1, z] <- Re(eigen_vals$vectors[,1]) / sum(Re(eigen_vals$vectors[,1]))
    }
  }

  pref_df <- reshape2::melt(pref) %>% rename(cell = Var1, age = Var2)
  age_pref_df <- reshape2::melt(age_pref) %>% rename(cell = Var1, age = Var2)
  depth_df <- data.frame(cell = 1:n_cells, value = depth)
  diffusion_df <- reshape2::melt(diffusion_mat) %>% rename(from = Var1, to = Var2)
  taxis_df <- reshape2::melt(taxis) %>% rename(from = Var1, to = Var2, age = Var3)
  m_df <- reshape2::melt(m) %>% rename(from = Var1, to = Var2, age = Var3)
  movement_df <- reshape2::melt(movement) %>% rename(from = Var1, to = Var2, age = Var3)
  distribution_df <- reshape2::melt(distribution) %>% rename(cell = Var1, age = Var2, time = Var3, grid = Var4)

  return(list(
    total_pref = pref_df,
    age_pref = age_pref_df,
    depth_pref = depth_df,
    diffusion = diffusion_df,
    taxis = taxis_df,
    m_rate = m_df,
    distribution = distribution_df,
    movement = movement_df
  ))
}
# Visualize CTMC ----------------------------------------------------------
ctmc_list <- create_1d_chain(
  n_cells = 25,
  n_ages = 3,
  diffusion = 3,
  age_pref_coeff = seq(0, 0.5, length.out = 3),
  t_interval = c(1, 5, 10),
  grid_start = c(1, 13, 25)
)

# data munging
# rename distribution
ctmc_list$distribution <- ctmc_list$distribution %>%
  mutate(time = case_when(
    time == 1 ~ "t = 1",
    time == 2 ~ "t = 3",
    time == 3 ~ "t = 5",
    time == 4 ~ "t = Inf"
  ),
  grid = case_when(
    grid == 1 ~ "g = 1",
    grid == 2 ~ "g = 13",
    grid == 3 ~ "g = 25"
  ),
  age = paste("age =", age))

# rename ages
ctmc_list$taxis <- ctmc_list$taxis %>%
  mutate(age = paste("age =", age))


# Plots -------------------------------------------------------------------
colors <- c("#C73E3A", "#2A9D8F", "#E9B44C", "black")

# Cohort distribution
dist_plot <- ggplot(ctmc_list$distribution %>% filter(grid != 'g = 13'),
                    aes(x = cell, y = value, color = factor(time))) +
  geom_step(lwd = 1) +
  scale_color_manual(values = c("blue3", "steelblue", "#fc9272", "#de2d26")) +
  facet_grid(grid~age) +
  labs(x = 'Grid cell (g)', y = "Distribution", color = "Time")+
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.05, 0.35),
        legend.background = element_blank(),
        legend.key.size = unit(3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

# Preference functions
totalpref_plot <- ggplot(ctmc_list$total_pref, aes(x = cell, y = value, color = factor(age))) +
  geom_step(lwd = 1) +
  scale_color_manual(values = colors) +
  labs(x = 'Grid cell (g)', y = "Total Preference (Depth + Ages)", color = "Ages") +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.13, 0.81),
        legend.background = element_blank(),
        legend.key.size = unit(3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

# Depth Preference
depthpref_plot <- ggplot(ctmc_list$depth_pref, aes(x = cell, y = value)) +
  geom_step() +
  labs(x = 'Grid cell (g)', y = "Depth Preference") +
  theme_bw(base_size = 15)

# Age Preference
agepref_plot <- ggplot(ctmc_list$age_pref, aes(x = cell, y = value, color = factor(age))) +
  geom_step(lwd = 1) +
  scale_color_manual(values = colors) +
  labs(x = 'Grid cell (g)', y = "Age Preference", color = 'Ages') +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.13, 0.81),
        legend.background = element_blank(),
        legend.key.size = unit(3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

# Diffusion plot
diffusion_plot <- ggplot(ctmc_list$diffusion, aes(x = from, y = to, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red",
                       midpoint = 0) +
  annotate("text", x = 6.5, y = 23.5, label = 'Diffusion Matrix', size = 5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.9, 0.3),
        legend.background = element_blank(),
        legend.key.size = unit(3, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)) +
  labs(fill = '', x = 'From grid cell (g)', y = "To grid cell (g)")

# Taxis plot
taxis_plot <- ggplot(ctmc_list$taxis, aes(x = from, y = to, fill = value)) +
  geom_tile() +
  geom_text(
    data = ctmc_list$taxis %>% filter(age == unique(age)[1]) %>% slice(1),
    aes(x = 5.5, y = 23.5, label = "Taxis Matrix"),
    size = 5,
    inherit.aes = FALSE
  )  +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red",
                       midpoint = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~age, ncol = 1) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.9, 0.085),
        legend.background = element_blank(),
        legend.key.size = unit(3, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)) +
  labs(fill = '', x = 'From grid cell (g)', y = "To grid cell (g)")


# Combine Plots -----------------------------------------------------------
# Combine preference plots
pref_plots <- cowplot::plot_grid(
  totalpref_plot,
  agepref_plot,
  depthpref_plot,
  nrow = 1,
  labels = c("C", "D", "E"),
  label_size = 20
)

# matrix plots
rate_plots <- cowplot::plot_grid(
  diffusion_plot,
  taxis_plot,
  ncol = 1,
  rel_heights = c(0.25,0.75),
  labels = c('A', "B"),
  label_size = 20
)

# combine dist and preference plots
left_plots <- cowplot::plot_grid(
  pref_plots,
  dist_plot,
  ncol = 1,
  rel_heights = c(0.35, 0.65),
  labels = c("", "F"),
  label_size = 20
)

ggsave(
  here("dev", "paper_projects", "spatial_movement_ctmc", "ctmc_conceptual.png"),
  cowplot::plot_grid(
    rate_plots, left_plots, rel_widths = c(0.25, 0.75)
  ),
  width = 17, height = 13
)

