# Purpose: qq panels of the dsem covariate OSA residuals under a correct model, one per family
# Creator: Matthew LH. Cheng (UAF-CFOS)

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(cowplot)
devtools::load_all(here("R"))
source(here("tests", "testthat", "helper-sweep_setup.R"))

n_reps <- 30
rho <- 0.5 # the covariate's autoregression
state_sd <- 0.5 # and its innovation sd, both kept fixed
known_sd <- 0.25 # the measurement sd gaussian_fixed_sd is told
random_pars <- c("ln_RecDevs", "dsem_x")

# a small model, since the residual code is the same whatever the assessment around it
input_list <- sweep_input(dims = list(n_yrs = 15, n_ages = 5, n_regions = 1, n_sexes = 1,
                                      n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 1, n_pop = 1))
yrs <- input_list$data$years
n_yrs <- length(yrs)

# rec has to appear but nothing joins it to env, so env alone is a draw from the model's own prior
arrows <- c(paste0("env -> env, 1, NA, ", rho), paste0("env <-> env, 0, NA, ", state_sd),
            "rec <-> rec, 0, NA, 0.9")

# the spread each family reads, and where its state sits on the link scale
obs_spread <- c(normal = 0.3, gaussian_fixed_sd = 0, lognormal = 0.3, gamma = 0.4,
                bernoulli = 0, poisson = 0, tweedie = 1.2)
offset <- c(normal = 0, gaussian_fixed_sd = 0, lognormal = 0, gamma = log(5),
            bernoulli = 0, poisson = log(8), tweedie = log(3))

# Replicates --------------------------------------------------------------

res_all <- data.frame()
t_start <- Sys.time()

for(family in names(obs_spread)) {

  fam_code <- unname(dsem_family_codes()[family])
  link_code <- dsem_default_link(fam_code)
  osa_fam <- if(family %in% c("bernoulli", "poisson", "tweedie")) family else "continuous"

  for(rep_i in seq_len(n_reps)) {

    # the dsem's own prior: year one is the innovation alone, then the ar1 runs forward
    set.seed(1000 * match(family, names(obs_spread)) + rep_i)
    x <- numeric(n_yrs)
    x[1] <- rnorm(1, 0, state_sd)
    for(t in 2:n_yrs) x[t] <- rho * x[t - 1] + rnorm(1, 0, state_sd)
    x <- x + offset[[family]]

    # the operating model's own draw, so the residual is checked against the density it was fit with
    obs <- as.numeric(draw_dsem_cov_obs(array(x, dim = c(n_yrs, 1, 1)), fam_code, link_code, obs_spread[[family]], 1.5,
                                        if(family == "gaussian_fixed_sd") rep(known_sd, n_yrs) else NULL))

    # the mean pinned where the state was drawn, the spread at the value it was drawn with, no refit
    d <- suppressMessages(Setup_Mod_DSEM(input_list, arrows, data.frame(year = yrs, env = obs),
                                         dsem_mu_spec = c(env = offset[[family]]), dsem_family = c(env = family),
                                         dsem_fixed_sd = if(family == "gaussian_fixed_sd") data.frame(year = yrs, env = rep(known_sd, n_yrs)) else NULL))
    if(obs_spread[[family]] > 0) d$par$ln_dsem_obs_sd[] <- log(obs_spread[[family]])
    d$map$ln_dsem_obs_sd <- factor(rep(NA, length(d$par$ln_dsem_obs_sd)))
    fit <- suppressWarnings(fit_model(d$data, d$par, d$map, random = random_pars, do_optim = FALSE, silent = TRUE))

    # a randomisation seed per replicate, or every discrete residual takes the same uniform draw
    osa_rep <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE, family = osa_fam, seed = 50000 + rep_i))
    res_all <- rbind(res_all, data.frame(family = family, rep = rep_i, osa_rep$res))

  } # end rep_i loop

  cat(format(Sys.time(), "%H:%M:%S"), family, "done, minutes",
      round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "\n")

} # end family loop

write.csv(res_all, here("dev", "make_sporc_obj_figs", "osa_dsem_validation_residuals.csv"), row.names = FALSE)

# Figure ------------------------------------------------------------------

# plot_resids draws the qq panel and annotates it with the sd, so each family's pooled residuals go
# through it as one data source, with the family on the x axis in place of a title
qq_family <- function(which_family) {
  pooled <- res_all[res_all$family == which_family & is.finite(res_all$resid), ]
  plot_resids(list(res = pooled))[[1]] + ggplot2::labs(x = which_family)
}

png(here("vignettes", "figures", "u_internal_dsem_validation.png"), width = 1400, height = 800)
cowplot::plot_grid(plotlist = lapply(names(obs_spread), qq_family), ncol = 4)
dev.off()
