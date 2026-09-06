# Purpose: To compare between unstructured markov and ctmc via simulation
# Creator: Matthew LH. Cheng
# UAF CFOS 10/28/25


# setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)
library(doParallel)
library(doSNOW)
library(foreach)
library(igraph)

# source basic EM setup
source(here("dev", "spatial_movement_ctmc", "setup_basic_spatial_em.R"))

# Run Simulation for Constant Movement ----------------------------------------------------------
n_cores <- parallel::detectCores() - 3
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# load in sim data
const <- readRDS(here("dev", "dev_output", "CTMC", "const.RDS"))

### Setup CTMC Components ---------------------------------------------------

# Setup adjacency matrix
adjacency <- as_adjacency_matrix(
  make_graph(
    ~ 1 - 2,
    2 - 3,
    1 - 3
  )
)

# make ctmc data
ctmc_data <- expand.grid(
  regions = 1:const$n_regions,
  years = 1:const$n_years,
  seas = 1:const$n_seas,
  ages = 1:const$n_ages,
  sexes = 1:const$n_sexes
)

# setup formulas for CTMC
diffusion_formula = ~0 + factor(regions) # constant diffusion
preference_formula = ~0

# setup model grid
model_df <- data.frame(
  move_type = c(0, 1),
  ctmc_diffusion_bounds = c(NA, 0),
  stringsAsFactors = FALSE
)

# movement age block
model_df$ageblock <- list(
  "constant",
  "constant"
)

# movement year blocks
model_df$yearblock <- list(
  "constant",
  "constant"
)


# parameter grid to loop through
param_grid <- expand.grid(
  sim_id = 1:const$n_sims,
  model_id = 1:nrow(model_df)
)


### Run Model ---------------------------------------------------------------

# setup progress par
n_tasks <- nrow(param_grid)
pb <- txtProgressBar(max = n_tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results_list <- foreach(
  row_idx = 1:n_tasks,
  .options.snow = opts,
  .packages = c("SPoRC", "tidyverse", "here")
) %dopar% {
                          devtools::load_all(here("R")) # load SPoRC in
                          fit_single_model(
                            row_idx = row_idx,
                            sim_obj = const,
                            model_df = model_df,
                            ctmc_data = ctmc_data,
                            adjacency = adjacency,
                            diffusion_formula = diffusion_formula,
                            preference_formula = preference_formula
                          )
                        }

close(pb)
stopCluster(cl)

saveRDS(results_list, here("dev", "dev_output", "CTMC", "const_sim_results.RDS")) # save results
beepr::beep(4)

# Run Simulation for Age Varying Movement ----------------------------------------------------------
n_cores <- parallel::detectCores() - 3
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# load in sim data
age_move <- readRDS(here("dev", "dev_output", "CTMC", "age_move.RDS"))

### Setup CTMC Components ---------------------------------------------------

# Setup adjacency matrix
adjacency <- as_adjacency_matrix(
  make_graph(
    ~ 1 - 2,
    2 - 3,
    1 - 3
  )
)
# make ctmc data
ctmc_data <- expand.grid(
  regions = 1:age_move$n_regions,
  years = 1:age_move$n_years,
  seas = 1:age_move$n_seas,
  ages = 1:age_move$n_ages,
  sexes = 1:age_move$n_sexes
)

# setup formulas for CTMC
diffusion_formula = ~1
preference_formula = ~ 0 + factor(regions):splines2::bSpline(ages, df = 4, intercept = TRUE)

# setup model grid
model_df <- data.frame(
  move_type = c(0, 1),
  ctmc_diffusion_bounds = c(NA, 0),
  stringsAsFactors = FALSE
)

# movement age block
model_df$ageblock <- list(
  list(1:5, 6:10, 11:15),
  "constant"
)

# movement year blocks
model_df$yearblock <- list(
  "constant",
  "constant"
)

# parameter grid to loop through
param_grid <- expand.grid(
  sim_id = 1:age_move$n_sims,
  model_id = 1:nrow(model_df)
)


### Run Model ---------------------------------------------------------------

# setup progress par
n_tasks <- nrow(param_grid)
pb <- txtProgressBar(max = n_tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results_list <- foreach(
  row_idx = 1:n_tasks,
  .options.snow = opts,
  .packages = c("SPoRC", "tidyverse", "here")
) %dopar% {
                          devtools::load_all(here("R")) # load SPoRC in
                          fit_single_model(
                            row_idx = row_idx,
                            sim_obj = age_move,
                            model_df = model_df,
                            ctmc_data = ctmc_data,
                            adjacency = adjacency,
                            diffusion_formula = diffusion_formula,
                            preference_formula = preference_formula
                          )
                        }

close(pb)
stopCluster(cl)

saveRDS(results_list, here("dev", "dev_output", "CTMC", "age_move_sim_results.RDS")) # save results
beepr::beep(4)

# Run Simulation for Time Varying Movement ----------------------------------------------------------
n_cores <- parallel::detectCores() - 3
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# load in sim data
time_move <- readRDS(here("dev", "dev_output", "CTMC", "time_move.RDS"))

### Setup CTMC Components ---------------------------------------------------

# Setup adjacency matrix
adjacency <- as_adjacency_matrix(
  make_graph(
    ~ 1 - 2,
    2 - 3,
    1 - 3
  )
)

# make ctmc data
ctmc_data <- expand.grid(
  regions = 1:time_move$n_regions,
  years = 1:time_move$n_years,
  seas = 1:time_move$n_seas,
  ages = 1:time_move$n_ages,
  sexes = 1:time_move$n_sexes
)

# setup formulas for CTMC
diffusion_formula = ~1
preference_formula = ~ 0 + factor(regions):splines2::bSpline(years, df = 5, intercept = TRUE) # increase df to allow more wiggliness

# setup model grid
model_df <- data.frame(
  move_type = c(0, 1),
  ctmc_diffusion_bounds = c(NA, 0),
  stringsAsFactors = FALSE
)

# movement age block
model_df$ageblock <- list(
  "constant",
  "constant"
)

# movement year blocks
model_df$yearblock <- list(
  list(1:6, 7:11, 12:17, 18:23, 24:30),
  "constant"
)

# parameter grid to loop through
param_grid <- expand.grid(
  sim_id = 1:time_move$n_sims,
  model_id = 1:nrow(model_df)
)


### Run Model ---------------------------------------------------------------

# setup progress par
n_tasks <- nrow(param_grid)
pb <- txtProgressBar(max = n_tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results_list <- foreach(
  row_idx = 1:n_tasks,
  .options.snow = opts,
  .packages = c("SPoRC", "tidyverse", "here")
) %dopar% {
                          devtools::load_all(here("R")) # load SPoRC in
                          fit_single_model(
                            row_idx = row_idx,
                            sim_obj = time_move,
                            model_df = model_df,
                            ctmc_data = ctmc_data,
                            adjacency = adjacency,
                            diffusion_formula = diffusion_formula,
                            preference_formula = preference_formula
                          )
                        }

close(pb)
stopCluster(cl)

saveRDS(results_list, here("dev", "dev_output", "CTMC", "time_move_sim_results.RDS")) # save results
beepr::beep(4)

# Run Simulation for Time and Age Varying Movement ----------------------------------------------------------
n_cores <- parallel::detectCores() - 3
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# load in sim data
timeage_move <- readRDS(here("dev", "dev_output", "CTMC", "timeage_move.RDS"))

### Setup CTMC Components ---------------------------------------------------

# Setup adjacency matrix
adjacency <- as_adjacency_matrix(
  make_graph(
    ~ 1 - 2,
    2 - 3,
    1 - 3
  )
)

# make ctmc data
ctmc_data <- expand.grid(
  regions = 1:timeage_move$n_regions,
  years = 1:timeage_move$n_years,
  seas = 1:timeage_move$n_seas,
  ages = 1:timeage_move$n_ages,
  sexes = 1:timeage_move$n_sexes
)

# setup formulas for CTMC
diffusion_formula = ~1
preference_formula = ~ 0 + factor(regions):splines2::bSpline(years, df = 5, intercept = TRUE):splines2::bSpline(ages, df = 4, intercept = TRUE)

# setup model grid
model_df <- data.frame(
  move_type = c(0, 1),
  ctmc_diffusion_bounds = c(NA, 0),
  stringsAsFactors = FALSE
)

# movement age block
model_df$ageblock <- list(
  list(1:5, 6:10, 11:15),
  "constant"
)

# movement year blocks
model_df$yearblock <- list(
  list(1:10, 11:20, 21:30),
  "constant"
)

# parameter grid to loop through
param_grid <- expand.grid(
  sim_id = 1:timeage_move$n_sims,
  model_id = 1:nrow(model_df)
)

### Run Model ---------------------------------------------------------------

# setup progress par
n_tasks <- nrow(param_grid)
pb <- txtProgressBar(max = n_tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results_list <- foreach(
  row_idx = 1:n_tasks,
  .options.snow = opts,
  .packages = c("SPoRC", "tidyverse", "here")
) %dopar% {
                          devtools::load_all(here("R")) # load SPoRC in
                          fit_single_model(
                            row_idx = row_idx,
                            sim_obj = timeage_move,
                            model_df = model_df,
                            ctmc_data = ctmc_data,
                            adjacency = adjacency,
                            diffusion_formula = diffusion_formula,
                            preference_formula = preference_formula
                          )
                        }

close(pb)
stopCluster(cl)

saveRDS(results_list, here("dev", "dev_output", "CTMC", "timeage_move_sim_results.RDS")) # save results
beepr::beep(4)

# Run Simulation for Age Varying Movement and Environmental Effect ----------------------------------------------------------
n_cores <- parallel::detectCores() - 3
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# load in sim data
age_env_move <- readRDS(here("dev", "dev_output", "CTMC", "age_env_move.RDS"))
temp_effect <- as.matrix(read.csv(here("dev", "dev_output", "CTMC", "age_env_temp.csv")))
temp_df <- reshape2::melt(temp_effect) %>%
  rename(regions = Var1, years = Var2, temp = value) %>%
  filter(years != "X") %>%
  mutate(years = as.numeric(str_remove(years, 'V')))

### Setup CTMC Components ---------------------------------------------------

# Setup adjacency matrix
adjacency <- as_adjacency_matrix(
  make_graph(
    ~ 1 - 2,
    2 - 3,
    1 - 3
  )
)

# make ctmc data
ctmc_data <- expand.grid(
  regions = 1:age_env_move$n_regions,
  years = 1:age_env_move$n_years,
  seas = 1:age_env_move$n_seas,
  ages = 1:age_env_move$n_ages,
  sexes = 1:age_env_move$n_sexes
)

ctmc_data <- ctmc_data %>%
  left_join(temp_df, by = c("regions", 'years'))

# setup formulas for CTMC
diffusion_formula = ~1
preference_formula = ~ 0 + factor(regions):splines2::bSpline(temp, df = 5, intercept = TRUE):splines2::bSpline(ages, df = 4, intercept = TRUE)

# setup model grid
model_df <- data.frame(
  move_type = c(0, 1),
  ctmc_diffusion_bounds = c(NA, 0),
  stringsAsFactors = FALSE
)

# movement age block
model_df$ageblock <- list(
  list(1:5, 6:10, 11:15),
  "constant"
)

# movement year blocks
model_df$yearblock <- list(
  list(1:10, 11:20, 21:30),
  "constant"
)

# parameter grid to loop through
param_grid <- expand.grid(
  sim_id = 1:timeage_move$n_sims,
  model_id = 1:nrow(model_df)
)



### Run Model ---------------------------------------------------------------

# setup progress par
n_tasks <- nrow(param_grid)
pb <- txtProgressBar(max = n_tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results_list <- foreach(
  row_idx = 1:n_tasks,
  .options.snow = opts,
  .packages = c("SPoRC", "tidyverse", "here")
) %dopar% {
                          devtools::load_all(here("R")) # load SPoRC in
                          fit_single_model(
                            row_idx = row_idx,
                            sim_obj = age_env_move,
                            model_df = model_df,
                            ctmc_data = ctmc_data,
                            adjacency = adjacency,
                            diffusion_formula = diffusion_formula,
                            preference_formula = preference_formula
                          )
                        }

close(pb)
stopCluster(cl)

saveRDS(results_list, here("dev", "dev_output", "CTMC", "age_env_move_sim_results.RDS")) # save results
beepr::beep(4)







