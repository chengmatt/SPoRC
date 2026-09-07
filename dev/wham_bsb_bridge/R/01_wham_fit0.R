# Purpose: rebuild the multi-WHAM black sea bass fit_0 report at its published MLE
# Date Created: 9/7/26

library(wham)
library(here)

bsb_dir <- "/Users/matthewcheng/Downloads/multi_wham_bsb-main"
out_dir <- here("dev", "wham_bsb_bridge", "output")

# wham 2.1.0.9009/9010 leaves selpars_ini at the asap row count when n_selblocks is
# overridden, so the bounds comparison in set_selectivity goes non-conformable
sel_src <- paste(readLines("/Users/matthewcheng/Downloads/wham-master/R/set_selectivity.R", warn = FALSE), collapse = "\n")

anchor <- "  if(!is.null(selectivity$n_selblocks)){ #override asap structure\n    data$n_selblocks <- selectivity$n_selblocks\n  }"

trim <- paste0(anchor, "\n",
               "  if(!is.null(selpars_ini) && NROW(selpars_ini) != data$n_selblocks) {\n",
               "    selpars_ini <- selpars_ini[seq_len(data$n_selblocks), , drop = FALSE]\n",
               "    if(!is.null(estimate_selpars)) estimate_selpars <- estimate_selpars[seq_len(data$n_selblocks), , drop = FALSE]\n",
               "  }")

patch_env <- new.env(parent = asNamespace("wham"))
eval(parse(text = sub(anchor, trim, sel_src, fixed = TRUE)), envir = patch_env)
patched_set_selectivity <- get("set_selectivity", envir = patch_env)
environment(patched_set_selectivity) <- asNamespace("wham")
utils::assignInNamespace("set_selectivity", patched_set_selectivity, ns = "wham")

asap <- read_asap3_dat(file.path(bsb_dir, "data", c("north.dat", "south.dat")))
temp <- prepare_wham_input(asap)

# 11 seasons, one month each apart from a two month june-july step
seasons <- c(rep(1, 5), 2, rep(1, 5)) / 12

basic_info <- list(region_names = c("North", "South"),
                   stock_names = paste0("BSB_", c("North", "South")))

basic_info$fracyr_seasons <- seasons
basic_info$NAA_where <- array(1, dim = c(2, 2, 8))
basic_info$NAA_where[1, 2, 1] <- 0 # north age 1 never starts a year in the south
basic_info$NAA_where[2, 1, ] <- 0  # south stock never occupies the north
basic_info$XSPR_R_avg_yrs <- which(temp$years > 1999)
basic_info$XSPR_R_opt <- 2

# bottom temperature, centered, ar1 process, no recruitment link in fit_0
north_bt <- read.csv(file.path(bsb_dir, "data", "bsb_bt_temp_nmab_1959-2022.csv"))
south_bt <- read.csv(file.path(bsb_dir, "data", "bsb_bt_temp_smab_1959-2022.csv"))

ecov <- list(label = c("North_BT", "South_BT"))
ecov$mean <- cbind(north_bt[, "mean"], south_bt[, "mean"])
ecov$mean <- t(t(ecov$mean) - apply(ecov$mean, 2, mean))
ecov$logsigma <- log(cbind(north_bt[, "se"], south_bt[, "se"]))
ecov$year <- north_bt[, "year"]
ecov$use_obs <- matrix(1, NROW(ecov$mean), NCOL(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- apply(ecov$mean, 2, mean)
ecov$recruitment_how <- matrix("none", 2, 2)

# state-space rec+1 with 2dar1; north fish sitting in the south are near scaa
NAA_re <- list(sigma = list("rec+1", "rec+1"),
               cor = list("2dar1", "2dar1"),
               N1_model = rep("equilibrium", 2))

NAA_re$sigma_vals <- array(1, dim = c(2, 2, temp$data$n_ages))
NAA_re$sigma_vals[1, 2, 2:temp$data$n_ages] <- 0.05
NAA_re$sigma_map <- array(NA, dim = c(2, 2, temp$data$n_ages))
NAA_re$sigma_map[1, 1:2, 1] <- 1
NAA_re$sigma_map[2, 2, 1] <- 3
NAA_re$sigma_map[1, 1, 2:temp$data$n_ages] <- 2
NAA_re$sigma_map[2, 2, 2:temp$data$n_ages] <- 4

cor_map <- array(NA, dim = dim(temp$par$trans_NAA_rho))
cor_map[1, 1, 1:3] <- 1:3
cor_map[2, 2, 1:3] <- 4:6
NAA_re$cor_map <- cor_map

# only the north stock moves, and it must return north before spawning
move <- list(stock_move = c(TRUE, FALSE), separable = TRUE)
move$must_move <- array(0, dim = c(2, length(seasons), 2))
move$must_move[1, 5, 2] <- 1
move$can_move <- array(0, dim = c(2, length(seasons), 2, 2))
move$can_move[1, 1:4, 2, 1] <- 1
move$can_move[1, 7:11, 1, 2] <- 1
move$can_move[1, 5, 2, ] <- 1

mus <- array(0, dim = c(2, length(seasons), 2, 1))
mus[1, 1:11, 1, 1] <- 0.02214863
mus[1, 1:11, 2, 1] <- 0.3130358
move$mean_vals <- mus
move$mean_model <- matrix("stock_constant", 2, 1)
move$use_prior <- array(0, dim = c(2, length(seasons), 2, 1))
move$use_prior[1, 1, 1, 1] <- 1
move$use_prior[1, 1, 2, 1] <- 1
move$prior_sigma <- array(0, dim = c(2, length(seasons), 2, 1))
move$prior_sigma[1, 1, 1, 1] <- 0.2
move$prior_sigma[1, 1, 2, 1] <- 0.2

# 8 blocks, four fleets then four indices
sel <- list(n_selblocks = 8,
            model = rep(c("age-specific", "logistic", "age-specific"), c(2, 2, 4)))

sel$initial_pars <- list(
  rep(c(0.5, 1), c(3, 5)),    # north commercial
  rep(c(0.5, 1), c(6, 2)),    # north recreational
  c(5, 1),                    # south commercial
  c(5, 1),                    # south recreational
  rep(c(0.5, 1, 1), c(1, 1, 6)), # north rec cpa
  rep(c(0.5, 1), c(4, 4)),    # north vast
  rep(c(0.5, 1, 1), c(2, 4, 2)), # south rec cpa
  rep(c(0.5, 1), c(1, 7))     # south vast
)

sel$fix_pars <- list(4:8, 7:8, NULL, NULL, 2:8, 5:8, 3:8, 2:8)
sel$re <- rep(c("2dar1", "2dar1", "none", "ar1_y", "2dar1", "none"), c(1, 1, 2, 1, 1, 2))

# comps get a flat input sample size because the fits use d-m and logistic-normal
catch_Neff <- temp$data$catch_Neff
catch_Neff[] <- 1000
catch_info <- list(catch_Neff = catch_Neff)
fleet_blocks <- temp$data$selblock_pointer_fleets
fleet_blocks[] <- rep(1:4, each = NROW(fleet_blocks))
catch_info$selblock_pointer_fleets <- fleet_blocks

index_Neff <- temp$data$index_Neff
index_Neff[] <- 1000
index_info <- list(index_Neff = index_Neff)
index_blocks <- temp$data$selblock_pointer_indices
index_blocks[] <- rep(5:8, each = NROW(index_blocks))
index_info$selblock_pointer_indices <- index_blocks
index_info$initial_index_sd_scale <- c(5, 1, 5, 1)
index_info$map_index_sd_scale <- c(1, NA, 2, NA)

age_comp <- list(
  fleets = c("dir-mult", "logistic-normal-miss0", "logistic-normal-ar1-miss0", "logistic-normal-ar1-miss0"),
  indices = c("logistic-normal-miss0", "dir-mult", "logistic-normal-ar1-miss0", "logistic-normal-ar1-miss0")
)

input_0 <- prepare_wham_input(asap,
                              selectivity = sel,
                              NAA_re = NAA_re,
                              basic_info = basic_info,
                              move = move,
                              ecov = ecov,
                              catch_info = catch_info,
                              index_info = index_info,
                              age_comp = age_comp)

input_0$fleet_names <- paste0(rep(c("North_", "South_"), each = 2), temp$fleet_names)
input_0$index_names <- paste0(rep(c("North_", "South_"), c(2, 2)), temp$index_names)

# the paper ships parameter lists but no fitted objects, so paste the mle back in
parLists <- readRDS(file.path(bsb_dir, "results", "parLists_no_M_re.RDS"))
par0 <- parLists[[1]]

# Ecov_obs_sigma_par is the one shape mismatch and is inert, obs sigma comes from logsigma
same_shape <- vapply(names(par0),
                     function(par_name) identical(dim(input_0$par[[par_name]]), dim(par0[[par_name]])),
                     logical(1))

input_0$par[names(par0)[same_shape]] <- par0[names(par0)[same_shape]]

fit0 <- fit_wham(input_0,
                 do.fit = FALSE,
                 do.brps = FALSE,
                 do.sdrep = FALSE,
                 do.osa = FALSE,
                 do.retro = FALSE,
                 MakeADFun.silent = TRUE)

saveRDS(list(input = input_0,
             rep = fit0$report(),
             parList = par0,
             nll = fit0$fn(),
             max_grad = max(abs(fit0$gr()))),
        file.path(out_dir, "01_wham_fit0.rds"))
