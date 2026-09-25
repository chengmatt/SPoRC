# Checks the general setup: series names from any deviation array, cells read back exactly, the same
# answer as the recruitment-only path (1e-8), and the map mirror switch.

data("sgl_rg_dusky_data")

test_that("series names and cells come out of the parameter array", {

  input_list <- sweep_input() # three regions, and movement, so both processes enumerate
  n_yrs <- length(input_list$data$years)
  link <- get_dsem_link(input_list, c("rec", "move"), c("x -> rec_Pop_1_Region_2, 0, b", "move_Pop_1_From_3_To_2_Seas_1_Age_4_Sex_2 -> x, 1, p"), n_yrs)

  expect_equal(link$offered[1:3], paste0("rec_Pop_1_Region_", 1:3))
  expect_equal(link$offered[4], "move_Pop_1_From_1_To_1_Seas_1_Age_1_Sex_1")
  expect_equal(length(link$offered), 3 + prod(dim(input_list$par$move_devs)[-4])) # every dim but the year

  expect_equal(link$name, c("rec_Pop_1_Region_2", "move_Pop_1_From_3_To_2_Seas_1_Age_4_Sex_2"))
  expect_equal(link$par, c("ln_RecDevs", "move_devs"))
  expect_equal(link$penalty_reads_map, c(TRUE, TRUE)) # both penalties skip a blanked cell

  # the cells are where those years sit in the flattened array
  input_list$par$ln_RecDevs[1,2,] <- 100 + seq_len(n_yrs)
  input_list$par$move_devs[1,3,2,,1,4,2] <- 200 + seq_len(n_yrs)
  expect_equal(as.numeric(input_list$par$ln_RecDevs[link$cell[[1]]]), 100 + seq_len(n_yrs))
  expect_equal(as.numeric(input_list$par$move_devs[link$cell[[2]]]), 200 + seq_len(n_yrs))

})

test_that("every process in the table reads its mirror, so every one can be linked", {

  # the switch blanks map_<parameter>, so a process is linkable only if its penalty honours that.
  # sort(unique(map_sel_devs)) drops NA for growth; NAA and movement each skip a blanked cell
  reads <- vapply(dsem_process_table(), function(e) e$penalty_reads_map, logical(1))
  expect_true(all(reads))
  expect_equal(vapply(dsem_process_table(), function(e) e$label, ""), c("rec", "growth", "growth_semipar", "NAA", "move", "fish_q", "srv_q"))

})

test_that("an empty array is skipped and growth does not answer for growth_semipar", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)

  # dusky is one region, so move_devs has no destination and the process offers nothing
  expect_equal(dim(input_list$par$move_devs)[3], 0L)
  expect_equal(get_dsem_link(input_list, c("rec", "move"), "env -> rec, 0, b", 48)$offered, "rec")

  # the two growth arrays share a prefix, so the name match has to stop at the underscore
  input_list$par$ln_growth_devs <- array(0, dim = c(1, 1, 48, 1, 1))
  input_list$par$ln_growth_semipar_devs <- array(0, dim = c(1, 1, 48, 1, 1))
  expect_equal(get_dsem_link(input_list, c("growth", "growth_semipar"), "env -> growth, 0, b", 48)$name, "growth")
  expect_equal(get_dsem_link(input_list, c("growth", "growth_semipar"), "env -> growth_semipar, 0, b", 48)$name, "growth_semipar")

})

test_that("an unknown process, a clashing covariate name and a missing mirror are refused, and movement links", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  env <- data.frame(year = input_list$data$years, env = rnorm(length(input_list$data$years)))

  expect_error(get_dsem_link(input_list, "biomass", "env -> rec, 0, b", 48), "cannot be linked")
  expect_error(Setup_Mod_DSEM(input_list, c("env -> env, 1, r", "env <-> env, 0, s"), env), "No deviation series")
  expect_error(Setup_Mod_DSEM(input_list, c("rec -> rec, 1, r", "rec <-> rec, 0, s"),
                                    data.frame(year = input_list$data$years, rec = 1)), "name of a deviation series")

  # growth deviations exist but this model keeps no mirror for them, so the penalty cannot come off
  input_list$par$ln_growth_devs <- array(0, dim = c(1, 1, 48, 1, 1))
  expect_error(Setup_Mod_DSEM(input_list, c("env -> growth, 0, b", "env <-> env, 0, s", "growth <-> growth, 0, sg"), env,
                                    dsem_processes = "growth"), "no map mirror")

  # a series whose deviations are all fixed by the map has nothing to describe
  fixed <- input_list; fixed$map$ln_RecDevs <- factor(rep(NA, length(input_list$par$ln_RecDevs)))
  expect_error(Setup_Mod_DSEM(fixed, c("env -> rec, 0, b", "env <-> env, 0, s", "rec <-> rec, 0, sr"), env), "nothing for the dsem to describe")

  # a movement series links now that Get_move_PE_loglik skips a blanked cell. the sweep model maps
  # every movement deviation off, so estimate one series first, otherwise there is nothing to describe
  sweep <- sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA)) # estimated, since fixed movement never reads its deviations
  n_yrs <- length(sweep$data$years)
  ms <- "move_Pop_1_From_1_To_1_Seas_1_Age_1_Sex_1"
  expect_error(Setup_Mod_DSEM(sweep, c(paste0("x -> ", ms, ", 0, b"), "x <-> x, 0, s", paste0(ms, " <-> ", ms, ", 0, sm")),
                              data.frame(year = sweep$data$years, x = rnorm(n_yrs)), dsem_processes = "move"), "nothing for the dsem to describe")
  map_move <- array(as.integer(sweep$map$move_devs), dim = dim(sweep$par$move_devs))
  map_move[1,1,1,,1,1,1] <- seq_len(n_yrs)
  sweep$map$move_devs <- factor(map_move)
  sweep$data$map_move_devs <- array(as.numeric(sweep$map$move_devs), dim = dim(sweep$par$move_devs))
  linked <- Setup_Mod_DSEM(sweep, c(paste0("x -> ", ms, ", 0, b"), "x <-> x, 0, s", paste0(ms, " <-> ", ms, ", 0, sm")),
                           data.frame(year = sweep$data$years, x = rnorm(n_yrs)), dsem_processes = "move")
  expect_equal(linked$data$dsem_link_par, "move_devs")
  expect_equal(length(linked$data$dsem_link_row[[1]]), n_yrs)

})

test_that("the general setup gives the recruitment-only setup's answer on dusky", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  set.seed(20)
  input_list$par$ln_RecDevs[] <- rnorm(length(input_list$par$ln_RecDevs), 0, 0.6)
  env <- data.frame(year = input_list$data$years, env = rnorm(length(input_list$data$years)))
  sigmaR <- exp(input_list$par$ln_sigmaR[2,1,1])
  arrows <- c(sprintf("rec <-> rec, 0, NA, %.17g", sigmaR), "env <-> env, 0, NA, 1")

  value_of <- function(model, il) {
    data <- apply_dsem_link_switch(sync_dev_map_data(il$data, il$map))
    obj <- RTMB::MakeADFun(cmb(model, data), il$par, map = il$map, silent = TRUE)
    obj$fn(obj$par)
  }

  old <- Setup_Mod_DSEM(input_list, arrows, env, dsem_mu_spec = "fix")
  new <- Setup_Mod_DSEM(input_list, arrows, env, dsem_mu_spec = "fix")
  env_nLL <- -sum(dnorm(env$env, mean(env$env), 1, log = TRUE)) # the unlinked covariate's own density

  expect_equal(value_of(SPoRC_rtmb, new), value_of(SPoRC_rtmb, old), tolerance = 1e-8)
  expect_equal(value_of(SPoRC_rtmb, new) - env_nLL,
               value_of(SPoRC_rtmb, input_list), tolerance = 1e-8) # reproduces SPoRC's own penalty

})

test_that("two linked series share one grid and the gather differentiates", {

  input_list <- sweep_input()
  set.seed(21)
  input_list$par$ln_RecDevs[] <- rnorm(length(input_list$par$ln_RecDevs), 0, 0.4)
  n_yrs <- length(input_list$data$years)
  cov_data <- data.frame(year = input_list$data$years, sst = rnorm(n_yrs))
  arrows <- c("sst -> rec_Pop_1_Region_2, 0, b_2", "sst -> rec_Pop_1_Region_3, 1, b_3",
              "sst -> sst, 1, rho", "sst <-> sst, 0, s_sst", "rec_Pop_1_Region_2 <-> rec_Pop_1_Region_2, 0, s_2",
              "rec_Pop_1_Region_3 <-> rec_Pop_1_Region_3, 0, s_3")

  both <- Setup_Mod_DSEM(input_list, arrows, cov_data, dsem_mu_spec = "fix")
  expect_equal(both$data$dsem_var_names, c("sst", "rec_Pop_1_Region_2", "rec_Pop_1_Region_3"))
  expect_equal(both$data$dsem_link_par, c("ln_RecDevs", "ln_RecDevs"))

  # a linked series sits under its process' own level parameter, so its dsem mean stays at zero
  expect_equal(both$par$dsem_mu[2:3], c(0, 0))
  expect_true(all(is.na(as.integer(both$map$dsem_mu)[2:3])))

  # and its cells are not estimated twice, once here and once in the parameter array
  map_x <- matrix(as.integer(both$map$dsem_x), nrow = both$data$dsem_n_grid_yrs)
  expect_true(all(is.na(map_x[, 2:3])))

  data <- apply_dsem_link_switch(sync_dev_map_data(both$data, both$map))
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), both$par, map = both$map, silent = TRUE)
  g <- as.vector(obj$gr(obj$par))
  fd <- numeric(length(obj$par))
  for(i in seq_along(obj$par)) {
    h <- 1e-4 * max(1, abs(obj$par[i])); p1 <- p2 <- obj$par
    p1[i] <- p1[i] + h; p2[i] <- p2[i] - h
    fd[i] <- (obj$fn(p1) - obj$fn(p2)) / (2 * h)
  }
  # the floor here is the objective's own size, about 6e6, not the gradient: central differences
  # bottom out at 1.6e-6 at this step and get worse either side of it
  expect_lt(max(abs(g - fd) / (1 + abs(g))), 1e-5)

})

test_that("the switch blanks the linked cells and leaves the rest penalized", {

  input_list <- sweep_input() # three regions, so some recruitment series stay unlinked
  n_yrs <- length(input_list$data$years)
  cov_data <- data.frame(year = input_list$data$years, sst = rnorm(n_yrs))
  arrows <- c("sst -> rec_Pop_1_Region_2, 0, b", "sst -> sst, 1, rho", "sst <-> sst, 0, s_sst", "rec_Pop_1_Region_2 <-> rec_Pop_1_Region_2, 0, s_rec")

  linked <- Setup_Mod_DSEM(input_list, arrows, cov_data, dsem_mu_spec = "fix")
  expect_equal(linked$data$dsem_link_par, "ln_RecDevs")

  mirror <- apply_dsem_link_switch(sync_dev_map_data(linked$data, linked$map))$map_ln_RecDevs
  expect_true(all(is.na(mirror[1,2,]))) # the linked region loses its penalty
  expect_false(any(is.na(mirror[1,1,]))) # the others keep it
  expect_false(any(is.na(mirror[1,3,])))

})
