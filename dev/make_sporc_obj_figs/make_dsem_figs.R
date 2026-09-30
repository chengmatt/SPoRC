# Purpose: Render the figures for vignette("ag_dsem"), the dsem fits, self tests, cross test and non-stationary fits
# Creator: Matthew LH. Cheng
# Date: 9/30/26

library(here)
library(RTMB)
devtools::load_all(here())

fig_dir <- here("vignettes", "figures")
if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
n_cores <- min(10, future::availableCores())

# the vignette's chunk sizes in inches, so the panels look the same as they did when knitted
fig_open <- function(name, width, height) png(file.path(fig_dir, name), width = width, height = height, units = "in", res = 150)

# Dusky base configuration --------------------------------------------------
# the penalized fit reproduces the packaged one, which has no bias correction. the random
# effect fits take the full lognormal correction, so they read a base that asks for it
data("dusky_rtmb_model")
base <- list(data = dusky_rtmb_model$data, par = dusky_rtmb_model$parameters, map = dusky_rtmb_model$mapping, verbose = FALSE, store_config = FALSE)
re_base <- base
re_base$data$bias_correct_pe <- 1

# setup recruitment for dsem, penalized likelihood with sigmaR fixed
pen_base <- Setup_Mod_Rec(
  input_list = base,
  sigmaR_switch = 1,
  ln_sigmaR = array(-0.1068576, dim = c(2, base$data$n_pop, base$data$n_regions)),
  rec_model = "mean_rec",
  init_age_strc = 1,
  ln_global_R0 = log(2.7),
  t_spawn = base$data$t_spawn,
  RecDevs_model = "dsem" # set recdev model as dsem
)

pen_recdev <- Setup_Mod_DSEM(
  pen_base,
  dsem_arrows = "
  rec <-> rec, 0, NA, 0.8986536
  ",
  dsem_data = NULL
)

pen_recdev_mod <- fit_model(
  pen_recdev$data,
  pen_recdev$par,
  pen_recdev$map,
  NULL, silent = TRUE
)
pen_recdev_mod$sd_rep <- sdreport(pen_recdev_mod)

cat("== Penalized dsem against the packaged fit ==\n")
cat("jnLL dsem", pen_recdev_mod$rep$jnLL, " packaged", dusky_rtmb_model$rep$jnLL, "\n")

fig_open("ag_dsem_pen_recdevs.png", 7, 5)
plot(pen_recdev_mod$rep$dsem_x_grid, ylab = 'RecDev', xlab = 'Year')
lines(as.vector(dusky_rtmb_model$rep$ln_RecDevs))
dev.off()

# Random effect recruitment, dsem against iid -------------------------------
# random effect deviations take the full lognormal correction, so no ramp
dsem_base <- Setup_Mod_Rec(
  input_list = re_base,
  do_rec_bias_ramp = 0,
  sigmaR_switch = 1,
  ln_sigmaR = array(-0.1068576, dim = c(2, re_base$data$n_pop, re_base$data$n_regions)),
  rec_model = "mean_rec",
  init_age_strc = 1,
  ln_global_R0 = log(2.7),
  t_spawn = re_base$data$t_spawn,
  RecDevs_model = "dsem"
)

# iid state space with sigmaR estimated through the arrows
re_dsem_recdev <- Setup_Mod_DSEM(
  dsem_base,
  dsem_arrows = "
  rec <-> rec, 0, sd_rec, 0.8986536
  ",
  dsem_data = NULL
)

re_dsem_recdev_mod <- fit_model(
  re_dsem_recdev$data,
  re_dsem_recdev$par,
  re_dsem_recdev$map,
  "ln_RecDevs", silent = TRUE
)
re_dsem_recdev_mod$sd_rep <- sdreport(re_dsem_recdev_mod)

# the same model without the dsem, one sigmaR across the early and late period
nodsem_base <- Setup_Mod_Rec(
  input_list = re_base,
  do_rec_bias_ramp = 0,
  sigmaR_switch = 1,
  ln_sigmaR = array(-0.1068576, dim = c(2, re_base$data$n_pop, re_base$data$n_regions)),
  rec_model = "mean_rec",
  init_age_strc = 1,
  ln_global_R0 = log(2.7),
  t_spawn = re_base$data$t_spawn,
  RecDevs_model = "iid",
  sigmaR_spec = "est_shared_all"
)

re_nodsem_recdev_mod <- fit_model(
  nodsem_base$data,
  nodsem_base$par,
  nodsem_base$map,
  "ln_RecDevs", silent = TRUE
)
re_nodsem_recdev_mod$sd_rep <- sdreport(re_nodsem_recdev_mod)

cat("\n== Random effect iid recruitment, dsem against no dsem ==\n")
cat("jnLL dsem", re_dsem_recdev_mod$rep$jnLL, " no dsem", re_nodsem_recdev_mod$rep$jnLL, "\n")

fig_open("ag_dsem_re_recdevs.png", 7, 5)
plot(re_dsem_recdev_mod$rep$dsem_x_grid, ylab = 'RecDev', xlab = 'Year')
lines(as.vector(re_nodsem_recdev_mod$rep$ln_RecDevs))
dev.off()

# AR(1) recruitment through the arrows --------------------------------------
re_ar1_dsem_recdev <- Setup_Mod_DSEM(
  dsem_base,
  dsem_arrows = "
  rec -> rec, 1, ar1_rho, 0.5
  rec <-> rec, 0, sd_rec, 0.8986536
  ",
  dsem_data = NULL
)

re_ar1_dsem_recdev_mod <- fit_model(
  re_ar1_dsem_recdev$data,
  re_ar1_dsem_recdev$par,
  re_ar1_dsem_recdev$map,
  "ln_RecDevs", silent = TRUE
)
re_ar1_dsem_recdev_mod$sd_rep <- sdreport(re_ar1_dsem_recdev_mod)

cat("\n== AR(1) recruitment ==\n")
ar1_tab <- summary(re_ar1_dsem_recdev_mod$sd_rep, "fixed")
print(ar1_tab[rownames(ar1_tab) %in% c("dsem_beta", "ln_dsem_sd"), , drop = FALSE])

fig_open("ag_dsem_ar1_recdevs.png", 7, 5)
plot(re_ar1_dsem_recdev_mod$rep$dsem_x_grid, ylab = 'RecDev', xlab = 'Year', col = 'blue', ylim = c(-2, 3), type = 'l')
lines(re_dsem_recdev_mod$rep$dsem_x_grid, col = 'red')
dev.off()

# Mediation: two contrived covariates ---------------------------------------
# y is the packaged deviations plus noise, x is y plus noise, so x reaches recruitment only through y
set.seed(777)
rec_dev <- as.vector(dusky_rtmb_model$rep$ln_RecDevs[1,1,])
y <- as.vector(scale(rec_dev + rnorm(length(rec_dev), 0, 0.9)))
x <- as.vector(scale(y + rnorm(length(y), 0, 0.5)))
cov_df <- data.frame(year = base$data$years, x = x, y = y)

cat("\n== Mediation covariates ==\n")
print(c(cor_xy = cor(x, y),
        cor_x_rec = cor(x, rec_dev),
        cor_y_rec = cor(y, rec_dev),
        partial_x = coef(lm(rec_dev ~ x + y))["x"]))

# both covariates are standardized, so their means are fixed at zero
fit_dsem_cov <- function(arrows) {
  il <- Setup_Mod_DSEM(
    dsem_base,
    dsem_arrows = arrows,
    dsem_data = cov_df,
    dsem_mu_spec = c(x = 0, y = 0)
  )
  mod <- fit_model(il$data, il$par, il$map, "ln_RecDevs", silent = TRUE)
  mod$sd_rep <- sdreport(mod)
  mod$arrow_names <- c(il$data$dsem_model$beta_names, il$data$dsem_model$ln_sd_names)
  mod
} # end fit_dsem_cov

# x regressed on its own, with the mediator left out of recruitment
x_only <- fit_dsem_cov("
  x <-> x, 0, sd_x, 1
  x -> y, 0, b_xy, 0.9
  y <-> y, 0, sd_y, 0.5
  x -> rec, 0, b_x_rec, 0.3
  rec <-> rec, 0, sd_rec, 0.8986536
")

# both paths in, so the fit can identify a direct effect of x if one is there
both <- fit_dsem_cov("
  x <-> x, 0, sd_x, 1
  x -> y, 0, b_xy, 0.9
  y <-> y, 0, sd_y, 0.5
  x -> rec, 0, b_x_rec, 0
  y -> rec, 0, b_y_rec, 0.3
  rec <-> rec, 0, sd_rec, 0.8986536
")

cat("\n== Mediation arrow estimates ==\n")
med_fits <- list(x_only = x_only, both = both)
for(nm in names(med_fits)) {
  m <- med_fits[[nm]]
  tab <- summary(m$sd_rep, "fixed")
  tab <- tab[rownames(tab) %in% c("dsem_beta", "ln_dsem_sd"), , drop = FALSE]
  rownames(tab) <- m$arrow_names
  cat(nm, "\n")
  print(round(tab, 4))
} # end nm loop

# Self tests ----------------------------------------------------------------
# the three panels each self test draws: median relative error in SSB and recruitment by year,
# and the relative error of every arrow parameter across replicates
plot_self_test <- function(st, name) {

  fig_open(name, 11, 4.5)
  par(mfrow = c(1,3), mar = c(5, 5.5, 1, 1), mgp = c(3.2, 0.9, 0),
      cex.lab = 1.5, cex.axis = 1.25, las = 1)

  plot(apply((st$SSB - st$truth$SSB) / st$truth$SSB, 3, median),
       ylim = c(-0.5, 0.5), xlab = 'Year', ylab = 'Relative Error in SSB',
       type = 'l', lwd = 2)
  abline(h = 0, lty = 1)

  plot(apply((st$Rec - st$truth$Rec) / st$truth$Rec, 3, median),
       ylim = c(-0.5, 0.5), xlab = 'Year', ylab = 'Relative Error in Recruitment',
       type = 'l', lwd = 2)
  abline(h = 0, lty = 1)

  # the sd rows come back on the log scale, so they are compared on the natural one
  beta_mat <- (st$dsem_beta - st$truth$dsem_beta) / st$truth$dsem_beta
  sd_mat <- (exp(st$ln_dsem_sd) - exp(st$truth$ln_dsem_sd)) / exp(st$truth$ln_dsem_sd)
  dsem_par_mat <- rbind(beta_mat, sd_mat)
  rownames(dsem_par_mat) <- c("b_xy", "b_x", "b_y", "sd_x", "sd_y", "sd_rec")
  stripchart(as.data.frame(t(dsem_par_mat)),
             vertical = TRUE, method = "jitter", pch = 16, cex = 1.1, xaxt = "n", xlab = "",
             ylab = "Relative error", ylim = c(-3, 3))
  axis(1, at = 1:nrow(dsem_par_mat), labels = rownames(dsem_par_mat), las = 2)
  abline(h = 0, lty = 2)

  dev.off()

} # end plot_self_test

# conditional: the states are kept at their estimates and only the observations are redrawn
self_test_cond <- simulation_self_test(
  both$data,
  both$parameters,
  both$mapping,
  random = 'ln_RecDevs',
  rep = both$rep,
  sd_rep = both$sd_rep,
  n_sims = 10,
  newton_loops = 2,
  n_cores = n_cores,
  what = c("SSB", "Rec", "dsem_x_grid"),
  what_par = c("dsem_beta", "ln_dsem_sd"),
  do_par = TRUE,
  sim_type = 'conditional'
)
plot_self_test(self_test_cond, "ag_dsem_self_test_cond.png")

# conditional with near-zero observation error, so what fails to return is not sampling noise
self_test_cond_perfect <- simulation_self_test(
  both$data,
  both$parameters,
  both$mapping,
  random = 'ln_RecDevs',
  rep = both$rep,
  sd_rep = both$sd_rep,
  n_sims = 10,
  newton_loops = 2,
  n_cores = n_cores,
  what = c("SSB", "Rec", "dsem_x_grid"),
  what_par = c("dsem_beta", "ln_dsem_sd"),
  do_par = TRUE,
  perfect_data = TRUE,
  sim_type = 'conditional'
)
plot_self_test(self_test_cond_perfect, "ag_dsem_self_test_cond_perfect.png")

# joint: the parameters are drawn from the joint precision as well
both$sd_rep <- sdreport(both, getJointPrecision = TRUE)
self_test_joint <- simulation_self_test(
  both$data,
  both$parameters,
  both$mapping,
  random = 'ln_RecDevs',
  obj = both,
  rep = both$rep,
  sd_rep = both$sd_rep,
  n_sims = 10,
  newton_loops = 2,
  n_cores = n_cores,
  what = c("SSB", "Rec", "dsem_x_grid"),
  what_par = c("dsem_beta", "ln_dsem_sd"),
  do_par = TRUE,
  sim_type = 'joint'
)
plot_self_test(self_test_joint, "ag_dsem_self_test_joint.png")

# Cross test: operating model -----------------------------------------------
# the recruitment and initial age deviations are centered on minus half their variance here,
# and the estimating models below say the same thing
sim_list <- Setup_Sim_Dim(
  n_sims        = 30,  # number of simulations
  n_yrs         = 30,  # number of years
  n_regions     = 1,   # single region
  n_ages        = 10,  # number of ages
  n_lens        = NULL,# no length structure
  n_sexes       = 1,   # single sex
  n_fish_fleets = 1,   # one fishery fleet
  n_srv_fleets  = 1,   # one survey fleet
  n_pop         = 1,   # number of pops
  bias_correct_pe = "rec"
)

sim_list <- Setup_Sim_Containers(sim_list)

# thin age data, so the year class strengths are not already fixed by the compositions
sim_list <- Setup_Sim_Fishing(sim_list = sim_list,
                              fish_sel_input = replicate(
                                n = sim_list$n_sims,
                                array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 5))), each = sim_list$n_yrs),
                                      dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                                              sim_list$n_sexes, sim_list$n_fish_fleets))
                              ),
                              ISS_FishAgeComps = array(25, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                                                   sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
)

sim_list <- Setup_Sim_Biologicals(
  sim_list = sim_list,
  natmort_input = replicate(n = sim_list$n_sims, array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                                                                    sim_list$n_ages, sim_list$n_sexes))), # natural mortality
  WAA_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                   dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))), # weight at age
  WAA_fish_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))), # fishery weight at age
  WAA_srv_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                       dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))), # survey weight at age
  MatAA_input = replicate(n = sim_list$n_sims, array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                     dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))) # maturity at age
)

# the population starts in equilibrium, so the initial age deviations are zero rather than drawn at sigmaR
sim_list <- Setup_Sim_Rec(
  sim_list = sim_list,
  R0_input = replicate(n = sim_list$n_sims, expr = array(5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs))), # R0
  ln_sigmaR = array(log(1), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
  recruitment_opt = 'mean_rec',
  ln_InitDevs_input = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_ages - 1,
                                       sim_list$n_sexes, sim_list$n_sims))
)

# no tagging and no movement
sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
sim_list$Movement <- array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

# a believable survey, so a change in catchability is worth mistaking for abundance
sim_list <- Setup_Sim_Survey(
  sim_list = sim_list,
  srv_sel_input = replicate(
    n = sim_list$n_sims,
    array(rep(1 / (1 + exp(-1 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                  sim_list$n_sexes, sim_list$n_srv))
  ),
  ObsSrvIdx_SE = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets)),
  ISS_SrvAgeComps = array(25, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                      sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
)

# catchability deviations are given to the arrows instead of a process error of their own
sim_list <- Setup_Sim_q_devs(sim_list, srv_q_model = "dsem")

# a contrived climate linked map: a basin index warms the shelf and moves the prey field, those drive
# recruitment, temperature also shifts availability, and availability and recruitment drive catchability
arrows <- c("pdo -> pdo, 1, rho_pdo",
            "pdo <-> pdo, 0, sd_pdo",
            "pdo -> temp, 0, b_pdo_temp",
            "temp -> temp, 1, rho_temp",
            "temp <-> temp, 0, sd_temp",
            "pdo -> prey, 0, b_pdo_prey",
            "prey <-> prey, 0, sd_prey",
            "temp -> rec, 0, b_temp_rec",
            "prey -> rec, 0, b_prey_rec",
            "rec <-> rec, 0, sd_rec",
            "temp -> cond, 0, b_temp_cond",
            "cond <-> cond, 0, sd_cond",
            "cond -> srv_q, 0, b_cond_q",
            "rec -> srv_q, 0, b_rec_q",
            "srv_q <-> srv_q, 0, sd_q")

cov_names <- c("pdo", "temp", "prey", "cond") # the four covariates
env_obs_sd <- 0.2 # each covariate is observed with error

# true values for the causal map
om_values <- c(rho_pdo = 0.9, sd_pdo = 0.45, b_pdo_temp = 0.7, rho_temp = 0.4, sd_temp = 0.4,
               b_pdo_prey = -0.5, sd_prey = 0.5, b_temp_rec = 0.45, b_prey_rec = 0.35, sd_rec = 0.35,
               b_temp_cond = 0.6, sd_cond = 0.4, b_cond_q = 0.5, b_rec_q = 0.3, sd_q = 0.2)

# dsem time-series, on two processes at once
sim_list <- Setup_Sim_DSEM(sim_list,
                           dsem_arrows = arrows,
                           dsem_values = om_values,
                           dsem_processes = c("rec", "srv_q"),
                           dsem_cov_obs_sd = setNames(rep(env_obs_sd, length(cov_names)), cov_names))

set.seed(123)
sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL) # get simulated datasets

# the drawn grid, six columns in the order the arrows introduce them
fig_open("ag_dsem_sim_series.png", 11, 7)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 1.5, 1), cex.lab = 1.2, cex.axis = 1.1)
panels <- list(c("PDO", "1"), c("Temperature", "2"), c("Prey", "3"), c("Availability", "5"),
               c("Rec dev", "4"), c("Q dev", "6"))
for(pan in panels) {
  for(i in 1:sim_obj$n_sims) {
    y <- sim_obj$dsem_x_sim[, as.integer(pan[2]), i]
    if(i == 1) plot(y, type = 'l', xlab = 'Year', ylab = pan[1]) else lines(y)
  } # end i loop
} # end pan loop
dev.off()

# Cross test: estimation model ----------------------------------------------
# settings match the operating model; spec says which arrows the EM believes and which processes it declares
setup_em <- function(sim_obj, sim, spec) {

  # extract simulation data for the terminal year and this replicate
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = sim)

  input_list <- Setup_Mod_Dim(
    years = 1:sim_obj$n_years,
    ages = 1:sim_obj$n_ages,
    lens = sim_obj$n_lens,
    n_regions = sim_obj$n_regions,
    n_sexes = sim_obj$n_sexes,
    n_fish_fleets = sim_obj$n_fish_fleets,
    n_srv_fleets = sim_obj$n_srv_fleets,
    n_pop = sim_obj$n_pop,
    bias_correct_pe = "rec", # the same switch the operating model was given
    verbose = FALSE
  )

  # a declared recruitment leaves sigmaR unread; the arrows take the deviations, or they stay iid
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(1), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    rec_model = "mean_rec",
    sigmaR_spec = if(spec$rec) "fix" else "fix_early_est_late",
    init_age_strc = "matrix",
    equil_init_age_strc = "stoch_all",
    ln_global_R0 = log(5),
    RecDevs_model = if(spec$rec) "dsem" else "iid"
  )

  # the operating model started in equilibrium, so the initial age deviations are known
  input_list$par$ln_InitDevs[] <- 0
  input_list$map$ln_InitDevs <- factor(rep(NA, length(input_list$par$ln_InitDevs)))

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = sim_data$WAA,
    MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0, # not fitting lengths
    AgeingError = sim_data$AgeingError,
    M_spec = "fix", # fixing natural mortality
    Fixed_natmort = array(0.3, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),  length(input_list$data$ages), input_list$data$n_sexes))
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sim_data$ObsCatch,
    UseCatch = sim_data$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = sim_data$UseFishIdx,
    ObsFishAgeComps = sim_data$ObsFishAgeComps,
    ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ISS_FishLenComps = sim_data$ISS_FishLenComps,
    fish_idx_type = c("biom"),
    FishAgeComps_LikeType = c("Multinomial"),
    FishLenComps_LikeType = c("none"),
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
    ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    srv_idx_type = c("biom"),
    SrvAgeComps_LikeType = c("Multinomial"),
    SrvLenComps_LikeType = c("none"),
    SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = c("logist2_Fleet_1"),
    fish_fixed_sel_pars_spec = c("est_all"),
    fish_q_spec = "est_all"
  )

  # "dsem" gives the catchability deviations to the arrows, "none" keeps catchability constant
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = c("logist2_Fleet_1"),
    srv_fixed_sel_pars_spec = c("est_all"),
    srv_q_spec = c("est_all"),
    srv_q_model = spec$q,
    sigma_srv_q_spec = if(spec$q == "dsem") "fix" else "est_all"
  )

  if(!is.null(spec$arrows)) {

    # the observed covariates, with their measurement error, are all the estimating model sees
    cov_df <- data.frame(year = 1:sim_obj$n_years)
    for(cov in spec$covs) cov_df[[cov]] <- sim_obj$dsem_cov_obs_sim[, match(cov, cov_names), sim]

    input_list <- Setup_Mod_DSEM(
      input_list,
      dsem_data = cov_df,
      dsem_arrows = spec$arrows,
      dsem_processes = c(if(spec$rec) "rec", if(spec$q == "dsem") "srv_q"),
      dsem_family = setNames(rep("normal", length(spec$covs)), spec$covs) # estimate obs error too
    )

  } # end if arrows

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets))
  )

  return(input_list)

} # end setup_em

# the arrows both estimating models share
shared <- c("pdo -> pdo, 1, rho_pdo",
            "pdo <-> pdo, 0, sd_pdo",
            "temp -> temp, 1, rho_temp",
            "pdo -> temp, 0, b_pdo_temp",
            "temp <-> temp, 0, sd_temp",
            "pdo -> prey, 0, b_pdo_prey",
            "prey <-> prey, 0, sd_prey",
            "temp -> rec, 0, b_temp_rec",
            "prey -> rec, 0, b_prey_rec",
            "rec <-> rec, 0, sd_rec",
            "temp -> cond, 0, b_temp_cond",
            "cond <-> cond, 0, sd_cond")

# the map the data were drawn under, and the same map with catchability constant
specs <- list(
  dsem_q = list(covs = cov_names, rec = TRUE, q = "dsem",
                arrows = c(shared,
                           "cond -> srv_q, 0, b_cond_q",
                           "rec -> srv_q, 0, b_rec_q",
                           "srv_q <-> srv_q, 0, sd_q")),
  const_q = list(covs = cov_names, rec = TRUE, q = "none", arrows = shared)
)

n_fit <- 1
n_yrs_em <- sim_obj$n_years

res <- NULL
ssb_ts <- array(NA_real_, dim = c(n_yrs_em, length(specs), n_fit), dimnames = list(NULL, names(specs), NULL))
q_ts <- ssb_ts

for(sim in 1:n_fit) {
  for(name in names(specs)) {

    spec <- specs[[name]]
    input_list <- setup_em(sim_obj, sim, spec)

    # what is random depends on whether catchability deviations exist
    has_q_dev <- !is.null(input_list$map$ln_srv_q_devs) && any(!is.na(input_list$map$ln_srv_q_devs))
    rand <- c("ln_RecDevs", if(has_q_dev) "ln_srv_q_devs", "dsem_x")
    model <- fit_model(input_list$data, input_list$par, input_list$map, random = rand, newton_loops = 1, silent = TRUE)

    par_list <- model$env$parList()
    ssb_ts[, name, sim] <- model$rep$SSB[1,1,]
    q_ts[, name, sim] <- model$rep$srv_q[1,,1]

    res <- rbind(res, data.frame(
      sim = sim,
      spec = name,
      max_grad = max(abs(model$gr(model$env$last.par.best[-model$env$random]))),
      betas = paste(round(par_list$dsem_beta, 2), collapse = " "),
      sd_q = if(spec$q == "dsem") exp(par_list$ln_dsem_sd[match("sd_q", input_list$data$dsem_model$ln_sd_names)]) else NA, # truth 0.2
      ssb_mean_abs_err = mean(abs(model$rep$SSB[1,1,] / sim_obj$SSB[1,1,,sim] - 1)),
      q_mean_abs_err = mean(abs(model$rep$srv_q[1,,1] / sim_obj$srv_q[1,,1,sim] - 1)),
      rec_dev_rmse = sqrt(mean((par_list$ln_RecDevs[1,1,] - sim_obj$ln_RecDevs[1,1,,sim])^2))
    ))

  } # end name loop
} # end sim loop

cat("\n== Cross test ==\n")
print(res)

# the error of every replicate against its own operating model values, then the median over them
median_err <- function(ts, truth) {
  err <- sweep(ts, c(1, 3), truth, "/") - 1 # [year, spec, replicate]
  apply(err, c(1, 2), median, na.rm = TRUE)
} # end median_err

ssb_err <- median_err(ssb_ts, sim_obj$SSB[1,1,,1:n_fit])
q_err <- median_err(q_ts, sim_obj$srv_q[1,,1,1:n_fit])

fig_open("ag_dsem_cross_test.png", 11, 9)
par(mfrow = c(2, 1), mar = c(4.5, 4.8, 1.5, 1), cex.lab = 1.3, cex.axis = 1.2)
cols <- c(dsem_q = "#000000", const_q = "#CC79A7")
yrs <- 1:n_yrs_em

plot(range(yrs), range(ssb_err, na.rm = TRUE), type = 'n', xlab = 'Year', ylab = 'Median relative error in SSB')
for(name in names(specs)) lines(yrs, ssb_err[, name], col = cols[name], lwd = 2.5)
abline(h = 0, lty = 2)
legend('topright', names(specs), col = cols, lwd = 2.5, bty = 'n', cex = 1.1)

plot(range(yrs), range(q_err, na.rm = TRUE), type = 'n', xlab = 'Year', ylab = 'Median relative error in survey q')
for(name in names(specs)) lines(yrs, q_err[, name], col = cols[name], lwd = 2.5)
abline(h = 0, lty = 2)
dev.off()

# Non-stationary: an effect that weakens over time --------------------------
# a contrived covariate that tracks recruitment before 2000 and is noise after, observed from 1990
years <- base$data$years
n_yrs <- length(years)
eps_std <- as.vector(scale(dusky_rtmb_model$rep$ln_RecDevs[1,1,])) # the packaged deviations, standardized

set.seed(17)
env_fade <- 0.7 * eps_std * (years < 2000) + rnorm(n_yrs, 0, 0.7)
env_fade[years < 1990] <- NA # the index begins in 1990, so the earlier years are states
env_fade <- as.vector(scale(env_fade))
fade_df <- data.frame(year = years, env = env_fade)

cat("\n== Fading covariate ==\n")
print(c(observed_years = sum(!is.na(env_fade)), cor_before_2000 = cor(env_fade[years %in% 1990:1999], eps_std[years %in% 1990:1999])))

# the path coefficient estimates and standard errors, named by arrow
dsem_pars <- function(fit, input_list) {
  fixed_names <- names(fit$sd_rep$par.fixed)
  est <- fit$sd_rep$par.fixed[fixed_names == "dsem_beta"]
  se <- sqrt(diag(fit$sd_rep$cov.fixed))[fixed_names == "dsem_beta"]
  names(est) <- names(se) <- input_list$data$dsem_model$beta_names
  round(rbind(est = est, se = se), 3)
} # end dsem_pars

# one coefficient for the whole series
flat_arrows <- c("env -> env, 1, rho_env",
                 "env <-> env, 0, sd_env",
                 "env -> rec, 0, b_env",
                 "rec <-> rec, 0, sd_rec")

flat_il <- Setup_Mod_DSEM(dsem_base, dsem_arrows = flat_arrows, dsem_data = fade_df)
flat_fit <- fit_model(flat_il$data, flat_il$par, flat_il$map, random = c("ln_RecDevs", "dsem_x"),
                      newton_loops = 2, silent = TRUE)
flat_fit$sd_rep <- sdreport(flat_fit)
print(dsem_pars(flat_fit, flat_il))

# the same effect as a series of its own, with its autocorrelation and innovation sd kept fixed
slope_arrows <- c("env -> env, 1, rho_env",
                  "env <-> env, 0, sd_env",
                  "env -> rec, 0, b_env",         # the effect of the covariate is the b_env series
                  "b_env -> b_env, 1, NA, 0.9",   # how fast it may move, kept fixed
                  "b_env <-> b_env, 0, NA, 0.15", # how far it may move each year, kept fixed
                  "rec <-> rec, 0, sd_rec")

slope_il <- Setup_Mod_DSEM(dsem_base,
                           dsem_arrows = slope_arrows,
                           dsem_data = cbind(fade_df, b_env = NA_real_), # unmeasured, so every year is a state
                           dsem_mu_spec = "b_env")

# start the effect flat at the single coefficient, as any state space fit is started
b_flat <- flat_fit$env$parList()$dsem_beta[match("b_env", flat_il$data$dsem_model$beta_names)]
slope_il$par$dsem_x[, match("b_env", slope_il$data$dsem_model$variables)] <- b_flat
slope_il$par$dsem_mu[match("b_env", slope_il$data$dsem_model$variables)] <- b_flat

slope_fit <- fit_model(slope_il$data, slope_il$par, slope_il$map, random = c("ln_RecDevs", "dsem_x"),
                       newton_loops = 2, silent = TRUE)
slope_fit$sd_rep <- sdreport(slope_fit)
print(dsem_pars(slope_fit, slope_il))

# the effect by year, and its standard error from the random effect part of the report
b_yr <- slope_fit$rep$dsem_x_grid[, match("b_env", slope_il$data$dsem_model$variables)]
in_b_col <- rep(slope_il$data$dsem_model$variables == "b_env", each = n_yrs) # the grid runs year within series
estimated <- !is.na(as.integer(slope_il$map$dsem_x))
is_x <- names(slope_fit$sd_rep$par.random) == "dsem_x"
b_se <- sqrt(slope_fit$sd_rep$diag.cov.random[is_x])[in_b_col[estimated]]

print(c(one_coefficient = flat_fit$optim$objective, effect_by_year = slope_fit$optim$objective))
print(c(mean_before_2000 = mean(b_yr[years < 2000]), mean_after_2000 = mean(b_yr[years >= 2000])))
cat("mean se on the effect", mean(b_se), "\n")

fig_open("ag_dsem_effect_by_year.png", 11, 4.5)
par(mfrow = c(1, 2), mar = c(5, 5.5, 1, 1), mgp = c(3.2, 0.9, 0), cex.lab = 1.4, cex.axis = 1.2, las = 1)

plot(years, env_fade, type = 'l', xlab = 'Year', ylab = 'Environmental covariate')
points(years, env_fade, pch = 16, cex = 0.8)

plot(years, b_yr, type = 'n', ylim = range(c(b_yr - 1.96 * b_se, b_yr + 1.96 * b_se)),
     xlab = 'Year', ylab = 'Effect on log recruitment')
polygon(c(years, rev(years)), c(b_yr - 1.96 * b_se, rev(b_yr + 1.96 * b_se)),
        col = '#00000022', border = NA)
lines(years, b_yr, lwd = 2.5)
abline(h = b_flat, lty = 2, lwd = 2)
legend('topright', c('Effect by year', 'One coefficient'), lwd = c(2.5, 2), lty = c(1, 2), bty = 'n', cex = 1.1)
dev.off()

# Non-stationary: a dome shaped recruitment response ------------------------
# the square of the covariate is a derived series, and its mean is subtracted through a constant series
set.seed(31)
dome_raw <- sqrt(pmax(0, (stats::quantile(eps_std, 0.95) - eps_std) / 2)) * sample(c(-1, 1), n_yrs, TRUE)
env_dome <- as.vector(scale(dome_raw + rnorm(n_yrs, 0, 0.25))) # recruitment is highest at intermediate values
sq_mean <- mean(env_dome^2)

dome_df <- data.frame(year = years,
                      one = 1,             # a constant series, so an arrow out of it is an intercept
                      env = env_dome,
                      env_sq = NA_real_,   # the square, which the arrows work out
                      env_ctr = NA_real_)  # and the square less its mean

dome_arrows <- c(
  # 1. a series that is 1 in every year. every series needs an sd line, and since this one is
  # data in every year the line only adds a constant, so it is kept fixed rather than estimated
  "one <-> one, 0, NA, 1",

  # 2. the covariate itself, an ar1 with its own innovation sd, both estimated
  "env -> env, 1, rho_env",
  "env <-> env, 0, sd_env",

  # 3. env_sq = env * env. the coefficient field names env, so in year t the coefficient is env_t
  # and env_sq_t = env_t * env_t. an sd of 0 says env_sq has no innovation, so it is only the square
  "env -> env_sq, 0, env",
  "env_sq <-> env_sq, 0, NA, 0",

  # 4. env_ctr = env_sq - mean(env_sq). two fixed arrows into env_ctr: the square with a coefficient
  # of 1, and the constant series with a coefficient of minus the mean. an sd of 0 again, so it is only that sum
  "env_sq -> env_ctr, 0, NA, 1",
  paste0("one -> env_ctr, 0, NA, ", -round(sq_mean, 4)),
  "env_ctr <-> env_ctr, 0, NA, 0",

  # 5. rec = b_lin * env + b_quad * env_ctr + innovation. the two coefficients and the sd are estimated
  "env -> rec, 0, b_lin",
  "env_ctr -> rec, 0, b_quad",
  "rec <-> rec, 0, sd_rec"
)

dome_il <- Setup_Mod_DSEM(dsem_base, dsem_arrows = dome_arrows, dsem_data = dome_df)
dome_fit <- fit_model(dome_il$data, dome_il$par, dome_il$map, random = c("ln_RecDevs", "dsem_x"),
                      newton_loops = 2, silent = TRUE)
dome_fit$sd_rep <- sdreport(dome_fit)

cat("\n== Dome shaped response ==\n")
print(dsem_pars(dome_fit, dome_il))

# the derived series read back off the grid, and the covariate value where recruitment is highest
dome_grid <- dome_fit$rep$dsem_x_grid
colnames(dome_grid) <- dome_il$data$dsem_model$variables
cat("max centered square error", max(abs(dome_grid[,'env_ctr'] - (dome_grid[,'env']^2 - sq_mean))), "\n")
dome_beta <- dsem_pars(dome_fit, dome_il)["est",]
cat("covariate at peak recruitment", -dome_beta[["b_lin"]] / (2 * dome_beta[["b_quad"]]), "\n")
cat("deviation in 2000", dome_grid[years == 2000, 'rec'], " packaged", dusky_rtmb_model$rep$ln_RecDevs[1,1,years == 2000], "\n")
cat("deviation range", range(dome_grid[,'rec']), " packaged", range(dusky_rtmb_model$rep$ln_RecDevs), "\n")

fig_open("ag_dsem_dome.png", 7, 5)
par(mar = c(5, 5.5, 1, 1), mgp = c(3.2, 0.9, 0), cex.lab = 1.4, cex.axis = 1.2, las = 1)
env_seq <- seq(min(env_dome), max(env_dome), length.out = 100)
plot(env_dome, dome_grid[,'rec'], pch = 16, xlab = 'Environmental covariate', ylab = 'Recruitment deviation')
lines(env_seq, dome_beta[["b_lin"]] * env_seq + dome_beta[["b_quad"]] * (env_seq^2 - sq_mean), lwd = 2.5)
dev.off()
