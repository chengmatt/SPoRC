# purpose: to create a 3-box spatial model with different movement parameterizations

library(here)
library(SPoRC)
library(here)
devtools::load_all(here("R"))

# Create movement matrix via CTMC with logistic preference for region 3
create_ctmc_model <- function(
  n_cells = 3,
  n_ages = 15,
  n_time = 30,
  n_sexes = 1,
  diffusion = 3,
  age_pref_coeff,
  use_age_preference = TRUE,
  use_time_preference = FALSE,
  use_env_preference = FALSE,
  env_effects = NULL,
  logistic_params = list(L = 1.5, k = 0.3, x0 = 3),
  sine_params = list(amplitude = 1, frequency = 1, phase = 0),
  from_cell_weights = NULL
) {
  # create containers
  diffusion_mat <- adjacency <- array(0, c(n_cells, n_cells))
  movement <- array(0, dim = c(n_cells, n_cells, n_time, 1, n_ages, n_sexes))
  m <- taxis <- array(0, dim = c(n_cells, n_cells, n_ages, n_time))
  age_pref <- array(0, dim = c(n_cells, n_ages))
  time_pref <- array(0, dim = c(n_cells, n_time))
  pref <- array(0, dim = c(n_cells, n_ages, n_time))

  # Create spatial cells
  regions <- 1:n_cells

  # create adjacency matrix with all regions connected
  adjacency[] <- 1
  diag(adjacency) <- 0

  # Get diffusion (constant by age and time)
  diffusion_mat[adjacency == 1] <- diffusion
  diag(diffusion_mat) <- -colSums(diffusion_mat)

  # Logistic function for region 3 preference (age-based)
  logistic_preference <- function(age, L, k, x0) {
    L / (1 + exp(-k * (age - x0)))
  }

  # Sine wave function for region 3 preference (time-based)
  sine_preference <- function(time, amplitude, frequency, phase) {
    # Scaled to range [0, amplitude]
    amplitude * (sin(2 * pi * frequency * time / n_time + phase) + 1) / 2
  }

  # Temperature effect
  gaussian_dome <- function(temp, optimal_temp, width) {
    exp(-((temp - optimal_temp)^2) / (2 * width^2))
  }

  # Calculate age preference if enabled
  if(use_age_preference) {
    for(a in 1:n_ages) {
      region3_strength <- logistic_preference(a,
                                              L = logistic_params$L,
                                              k = logistic_params$k,
                                              x0 = logistic_params$x0)

      for(i in 1:n_cells) {
        if(i == 3) {
          age_pref[i, a] <- region3_strength * age_pref_coeff[a]
        } else {
          age_pref[i, a] <- (1 - region3_strength) * age_pref_coeff[a] * (4 - i) / 2
        }
      }
    }
  }

  # Calculate time preference if enabled
  if(use_time_preference) {
    for(t in 1:n_time) {
      region3_strength <- sine_preference(t,
                                          amplitude = sine_params$amplitude,
                                          frequency = sine_params$frequency,
                                          phase = sine_params$phase)

      for(i in 1:n_cells) {
        if(i == 3) {
          time_pref[i, t] <- region3_strength
        } else {
          time_pref[i, t] <- (1 - region3_strength) * (4 - i) / 2
        }
      }
    }
  }

  if(use_env_preference) {
    env_pref <- array(0, dim = c(n_cells, n_ages, n_time))
    for(i in 1:n_cells) {
      for(t in 1:n_time) {
        temp_response <- gaussian_dome(env_effects$temp_array[i, t],
                                       optimal_temp = env_effects$optimal_temp,
                                       width = env_effects$temp_width) * 0.85
        env_pref[i, , t] <- temp_response  # Same across ages, varies by cell/time
      }
    }
  }

  # Combine preferences based on what's enabled
  for(a in 1:n_ages) {
    for(t in 1:n_time) {
      if(use_age_preference && use_time_preference && use_env_preference) {
        # All three enabled
        pref[, a, t] <- age_pref[, a] + time_pref[, t] + env_pref[, a, t]
      } else if(use_age_preference && use_time_preference) {
        # Age and time only
        pref[, a, t] <- age_pref[, a] + time_pref[, t]
      } else if(use_age_preference && use_env_preference) {
        # Age and env only
        pref[, a, t] <- age_pref[, a] + env_pref[, a, t]
      } else if(use_time_preference && use_env_preference) {
        # Time and env only
        pref[, a, t] <- time_pref[, t] + env_pref[, a, t]
      } else if(use_age_preference) {
        # Age only
        pref[, a, t] <- age_pref[, a]
      } else if(use_time_preference) {
        # Time only
        pref[, a, t] <- time_pref[, t]
      } else if(use_env_preference) {
        # Env only
        pref[, a, t] <- env_pref[, a, t]
      } else {
        # None enabled: no preference (uniform)
        pref[, a, t] <- 0
      }
    }
  }

  # Calculate movement matrices for each age and time combination
  for(a in 1:n_ages) {
    for(t in 1:n_time) {
      # get taxis for this age-time combination
      taxis_temp <- array(0, dim = c(n_cells, n_cells))

      for(i in 1:n_cells) {
        for(j in 1:n_cells) {
          if(adjacency[i,j] == 1) {
            taxis_temp[i, j] <- from_cell_weights[j] * (pref[i, a, t] - pref[j, a, t])
          }
        }
      }
      diag(taxis_temp) <- -colSums(taxis_temp)

      # get instantaneous movement rates
      m[,,a,t] <- diffusion_mat + taxis_temp

      # get movement fractions - replicate across sex dimension
      movement_age_time <- t(as.matrix(Matrix::expm(m[,,a,t]))) # from, to

      for(s in 1:n_sexes) {
        movement[,,t,1,a,s] <- movement_age_time
      }
    }
  }

  return(list(
    movement_mat = movement,
    age_pref = age_pref,
    time_pref = time_pref,
    combined_pref = pref,
    temp = env_effects$temp_array
  ))

}

# Get movement matrices here

# Age preference only
model_age <- create_ctmc_model(
  n_cells = 3,
  n_ages = 15,
  n_time = 30,
  n_sexes = 2,
  diffusion = 1,
  age_pref_coeff = rep(1, 15),
  use_age_preference = TRUE,
  use_time_preference = FALSE,
  logistic_params = list(L = 0.9, k = 0.5, x0 = 2),
  from_cell_weights = c(0.05, 0.9, 1)
)

plot(model_age$movement_mat[3,3,1,,1], ylim = c(0,1))
lines(model_age$movement_mat[3,2,1,,1])
lines(model_age$movement_mat[3,1,1,,1])

sum(model_age$movement_mat > 1)
sum(model_age$movement_mat < 0)

# Time preference only
model_time <- create_ctmc_model(
  n_cells = 3,
  n_ages = 15,
  n_time = 30,
  n_sexes = 2,
  diffusion = 1,
  age_pref_coeff = rep(1, 15),
  use_age_preference = FALSE,
  use_time_preference = TRUE,
  sine_params = list(amplitude = 1.05, frequency = 1.5, phase = 0),
  from_cell_weights = c(0.1, 0.35, 0.55)
)

plot(model_time$movement_mat[2,2,,2,1], ylim = c(0,1))
lines(model_time$movement_mat[3,3,,2,1])
lines(model_time$movement_mat[2,1,,2,1])

sum(model_time$movement_mat > 1)
sum(model_time$movement_mat < 0)

# both age and time preference
model_both <- create_ctmc_model(
  n_cells = 3,
  n_ages = 15,
  n_time = 30,
  n_sexes = 2,
  diffusion = 1,
  age_pref_coeff = rep(1, 15),
  use_age_preference = TRUE,
  use_time_preference = TRUE,
  logistic_params = list(L = 1, k = 2, x0 = 0.5),
  sine_params = list(amplitude = 0.5, frequency = 1.5, phase = 0),
  from_cell_weights = c(0.1, 0.1, 1)
)

plot(model_both$movement_mat[2,2,,2,1], ylim = c(0,1))
lines(model_both$movement_mat[2,1,,2,1])
lines(model_both$movement_mat[2,3,,2,1])

plot(model_both$movement_mat[2,2,1,,1], ylim = c(0,1))
lines(model_both$movement_mat[2,3,1,,1])
lines(model_both$movement_mat[2,1,1,,1])

sum(model_both$movement_mat > 1)
sum(model_both$movement_mat < 0)

# age and environmental effect
env_effects <- list()
env_effects$temp_array <- array(0, dim = c(n_cells, n_time))
env_effects$temp_array[1,] <- seq(-3, 2, length.out = n_time)
env_effects$temp_array[2,] <- seq(-2, 3, length.out = n_time)
env_effects$temp_array[3,] <- seq(-1, 1, length.out = n_time)
env_effects$optimal_temp <- 0
env_effects$temp_width <- 1

model_env <- create_ctmc_model(
  n_cells = 3,
  n_ages = 15,
  n_time = 30,
  n_sexes = 2,
  diffusion = 1,
  age_pref_coeff = rep(1, 15),
  use_age_preference = TRUE,
  use_time_preference = FALSE,
  use_env_preference = TRUE,
  env_effects = env_effects,
  logistic_params = list(L = 0.6, k = 0.5, x0 = 2),
  from_cell_weights = c(0.05, 0.9, 1)
)

plot(model_env$movement_mat[1,1,,2,1], ylim = c(0,1))
lines(model_env$movement_mat[1,2,,2,1])
lines(model_env$movement_mat[2,3,,2,1])

plot(model_env$movement_mat[2,2,1,,1], ylim = c(0,1))
lines(model_env$movement_mat[2,3,1,,1])
lines(model_env$movement_mat[2,1,1,,1])

sum(model_env$movement_mat > 1)
sum(model_env$movement_mat < 0)


# Setup Operating Model ---------------------------------------------------
### Setup Model Dimensions --------------------------------------------------
sim_list <- Setup_Sim_Dim(n_sims = 100, # number of simulations
                          n_yrs = 30, # number of years
                          n_regions = 3,  # number of regions
                          n_ages = 15, # number of ages
                          n_lens = NULL, # number of lengths
                          n_sexes = 1, # number of sexes
                          n_fish_fleets = 1, # number of fishery fleets
                          n_srv_fleets = 1, # number of survey fleets
                          n_seas = 1 # number of seasons
                          )

### Setup Simulation Containers ---------------------------------------------
sim_list <- Setup_Sim_Containers(sim_list)

### Setup Fishing Processes -------------------------------------------------
sim_list <- Setup_Sim_Fishing(sim_list = sim_list, # update simulate list
                              # Logistic selectivity
                              fish_sel_input = replicate(
                                n = sim_list$n_sims,
                                array(rep(1 / (1 + exp(-1 * ((1:sim_list$n_ages) - 5))), each = sim_list$n_yrs * sim_list$n_regions),
                                      dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                                              sim_list$n_sexes, sim_list$n_fish_fleets))
                              )
)

# Two way trip fishing mortality pattern
sim_list$Fmort[1,,,1,] <- replicate(sim_list$n_sims, c(seq(0.01, 0.15, length.out = sim_list$n_yrs/2), seq(0.15, 0.05, length.out = sim_list$n_yrs/2)))
sim_list$Fmort[2,,,1,] <- replicate(sim_list$n_sims, c(seq(0.01, 0.1, length.out = sim_list$n_yrs/2), seq(0.1, 0.05, length.out = sim_list$n_yrs/2)))
sim_list$Fmort[3,,,1,] <- replicate(sim_list$n_sims, c(seq(0.01, 0.05, length.out = sim_list$n_yrs/2), seq(0.05, 0.05, length.out = sim_list$n_yrs/2)))
sim_list$Fmort <- sim_list$Fmort * exp(rnorm(length(sim_list$Fmort), 0, 0.1)) # add some process error

### Setup Survey Process ----------------------------------------------------
sim_list <- Setup_Sim_Survey(
  sim_list = sim_list,
  # Logistic selectivity
  srv_sel_input = replicate(
    n = sim_list$n_sims,
    array(rep(1 / (1 + exp(-0.85 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs * sim_list$n_regions),
          dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages,
                  sim_list$n_sexes, sim_list$n_srv))
  ),
  ObsSrvIdx_SE = array(0.15, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))
)

### Setup Biological Dynamics -----------------------------------------------
sim_list <- Setup_Sim_Biologicals(
  sim_list = sim_list, # simualtion list
  natmort_input = replicate(n = sim_list$n_sims, array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs,
                                                                    sim_list$n_ages, sim_list$n_sexes))), # natural mortality
  WAA_input = replicate(n = sim_list$n_sims, array(rep(3 / (1 + exp(-0.85 * ((1:sim_list$n_ages) - 4))), each = sim_list$n_yrs * sim_list$n_regions),
                                                   dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))), # weight at age
  WAA_fish_input = replicate(n = sim_list$n_sims, array(rep(3 / (1 + exp(-0.85 * ((1:sim_list$n_ages) - 4))), each = sim_list$n_yrs * sim_list$n_regions),
                                                        dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))), # fishery weight at age
  WAA_srv_input = replicate(n = sim_list$n_sims, array(rep(3 / (1 + exp(-0.85 * ((1:sim_list$n_ages) - 4))), each = sim_list$n_yrs * sim_list$n_regions),
                                                       dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))), # survey weight at age
  MatAA_input = replicate(n = sim_list$n_sims, array(rep(1 / (1 + exp(-0.85 * ((1:sim_list$n_ages) - 7))), each = sim_list$n_yrs * sim_list$n_regions),
                                                     dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))) # maturity at age
)

### Setup Recruitment Processes ---------------------------------------------
sim_list <- Setup_Sim_Rec(
  sim_list = sim_list,
  R0_input = replicate(n = sim_list$n_sims, expr = array(c(10, 8, 5), dim = c(sim_list$n_regions, sim_list$n_yrs))), # R0
  ln_sigmaR = log(c(1, 1)),
  recruitment_opt = 'mean_rec',
  init_age_strc = 2,
  do_recruits_move = 0,
  rec_dd = "local",
  init_dd = 'local'
)


### Setup Tagging and Movement -----------------------------------------------------------
sim_list <- Setup_Sim_Tagging(
  sim_list = sim_list, # simulation list
  UseTagging = 1,
  n_tags = 200,
  max_liberty = sim_list$n_ages,
  tag_like = 2,
  tag_release_indicator = expand.grid(regions = 1:sim_list$n_regions, tag_years = seq(1, sim_list$n_yrs, 1)), # tags every 3 years
  t_tagging = 0.5,
  Tag_Reporting_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
  ln_Init_Tag_Mort = log(1e-100),
  ln_Tag_Shed = log(1e-100)
)


# Movement Scenarios ------------------------------------------------------
### Constant ----------------------------------------------------------------
const_move <- array(
  rep(c(
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8
  ),
  length.out = prod(c(sim_list$n_regions, sim_list$n_regions,
                      sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))),
  dim = c(sim_list$n_regions, sim_list$n_regions,
          sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes)
)

# input movement in
sim_list$Movement <- replicate(sim_list$n_sims, const_move)
const_move <- Simulate_Pop_Static(sim_list = sim_list, output_path = here("dev", "dev_output", "CTMC", "const.RDS")) # get simulated datasets

### Age Varying Movement ----------------------------------------------------
sim_list$Movement <- replicate(sim_list$n_sims, model_age$movement_mat)

age_move <- Simulate_Pop_Static(
  sim_list = sim_list,
  output_path = here("dev", "dev_output", "CTMC", "age_move.RDS")
)

### Time Varying Movement ---------------------------------------------------
# Input into simulation
sim_list$Movement <- replicate(sim_list$n_sims, model_time$movement_mat)
Simulate_Pop_Static(
  sim_list = sim_list,
  output_path = here("dev", "dev_output", "CTMC", "time_move.RDS")
)

### Age and Time Varying Movement ---------------------------------------------------
# Input into simulation
sim_list$Movement <- replicate(sim_list$n_sims, model_both$movement_mat)

# Input to simulation
Simulate_Pop_Static(
  sim_list = sim_list,
  output_path = here("dev", "dev_output", "CTMC", "timeage_move.RDS")
)


### Age and Env Effect ------------------------------------------------------
# Input into simulation
sim_list$Movement <- replicate(sim_list$n_sims, model_env$movement_mat)

# Input to simulation
Simulate_Pop_Static(
  sim_list = sim_list,
  output_path = here("dev", "dev_output", "CTMC", "age_env_move.RDS")
)

# write out temp effect
write.csv(model_env$temp, here("dev", "dev_output", "CTMC", "age_env_temp.csv"))

