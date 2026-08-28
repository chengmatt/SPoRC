# Purpose: Render the figures for vignette("af_growth_options"), the guide to
#          choosing the growth options in Setup_Mod_Biologicals(). Each figure
#          isolates one decision and shows what the switch does, built from the
#          same functions the objective function calls, so nothing here is a
#          sketch of the model's behavior.
#
#          Four figures:
#            1. the Schnute curve dissected, and what growth_model changes
#            2. the spread at age and the age-length key it produces
#            3. time variation, growth_tv_type curve against cohort, on the
#               EBS Pacific cod model
#            4. a semi-parametric surface recovered from simulated data, against
#               the parametric curve fit to the same data
#
#          Figures 3 and 4 reuse the bridge and self-test helpers, so the
#          vignette cannot drift from the regression tests. The reference ages
#          and the accumulator age in figure 1 are chosen so the switch being
#          demonstrated is visible: a species that is already at its asymptote
#          by the second reference age would show nothing.
# Creator: Matthew LH. Cheng
# Date Created: 8/24/26

library(here)
library(dplyr)
library(ggplot2)
devtools::load_all(here())

fig_dir <- here("vignettes", "figures")
thm <- ggplot2::theme_bw(base_size = 17) +
  ggplot2::theme(legend.position = "top",
                 legend.title = ggplot2::element_text(size = 14),
                 legend.text = ggplot2::element_text(size = 13),
                 strip.background = ggplot2::element_rect(fill = "gray92"))
lvl <- function(x, ...) factor(x, levels = c(...))

## 1. The curve --------------------------------------------------------------
# One set of parameters throughout, in the Schnute form the module reads: the
# length at a young reference age, the length at an old one, and the rate
# joining them, with the linear phase below A1 anchored at L0.

L0 <- 5; L1 <- 20; L2 <- 70; K <- 0.22; CV1 <- 0.14; CV2 <- 0.06
A1 <- 2; A2 <- 20
ages <- 0:25
lower <- seq(2, 92, by = 2)

crv <- get_laa_curve(ages, L0, L1, L2, K, CV1, CV2, A1, A2)

# (a) the pieces of the curve, with the interval one standard deviation either
# side of the mean, which is the spread the age-length key is built from
anat <- data.frame(Age = ages, L = crv$L, lo = crv$L - crv$sd, hi = crv$L + crv$sd)
p_anat <- ggplot(anat, aes(Age, L)) +
  annotate("rect", xmin = 0, xmax = A1, ymin = -Inf, ymax = Inf, fill = "gray90") +
  annotate("text", x = A1 / 2, y = 88, label = "linear\nphase", size = 4.2, lineheight = 0.9) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray60", alpha = 0.35) +
  geom_hline(yintercept = crv$Linf, linetype = "dashed", color = "gray30") +
  annotate("text", x = 21.5, y = crv$Linf + 4, label = "L[infinity]", parse = TRUE, size = 5) +
  geom_line(linewidth = 1.1) +
  geom_point(data = data.frame(Age = c(0, A1, A2), L = c(L0, L1, L2)), size = 3.2) +
  annotate("text", x = c(0, A1, A2) + 1.5, y = c(L0, L1, L2) - 5,
           label = c("L[0]", "L[1]~at~A[1]", "L[2]~at~A[2]"), parse = TRUE, size = 5) +
  ylim(0, 95) + labs(x = "Age", y = "Mean length") + thm

# (b) growth_model: the Richards coefficient bends the approach to the
# asymptote, and rho of one is exactly the von Bertalanffy curve
rho_df <- bind_rows(lapply(c(0.4, 1, 2.5), function(rr) {
  cc <- get_laa_curve(ages, L0, L1, L2, K, CV1, CV2, A1, A2, rho = rr)
  data.frame(Age = ages, L = cc$L,
             Option = if(rr == 1) "vb_schnute (rho = 1)" else paste0("richards, rho = ", rr))
}))
rho_df$Option <- lvl(rho_df$Option, "vb_schnute (rho = 1)", "richards, rho = 0.4", "richards, rho = 2.5")
p_rho <- ggplot(rho_df, aes(Age, L, color = Option, linetype = Option)) +
  geom_line(linewidth = 1.1) + ggthemes::scale_color_colorblind() +
  labs(x = "Age", y = "Mean length", color = "growth_model", linetype = "growth_model") +
  thm + theme(legend.direction = "vertical")

# (c) growth_A2: the same L2 read as the length at a reference age, from which
# the asymptote is solved, or read as the asymptote itself. The gap is the whole
# point of the switch, so the reference age here is young enough for the curve
# to still be climbing at it
A2_y <- 8
crv_y <- get_laa_curve(ages, L0, L1, L2, K, CV1, CV2, A1, A2_y)
a2_df <- bind_rows(
  data.frame(Age = ages, L = crv_y$L, Option = "growth_A2 = 8"),
  data.frame(Age = ages, L = get_laa_curve(ages, L0, L1, L2, K, CV1, CV2, A1, max(ages),
                                           L2_asymptote = 1)$L,
             Option = 'growth_A2 = "Linf"')
)
p_a2 <- ggplot(a2_df, aes(Age, L, color = Option, linetype = Option)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = L2, color = "gray40", linetype = "dotted") +
  geom_point(data = data.frame(Age = A2_y, L = L2), inherit.aes = FALSE,
             aes(Age, L), color = "gray20", size = 3) +
  annotate("text", x = 20, y = L2 - 4.5, label = "L[2]==70", parse = TRUE, size = 5) +
  ggthemes::scale_color_colorblind() +
  labs(x = "Age", y = "Mean length", color = "L2 is read as",
       linetype = "L2 is read as") +
  thm + theme(legend.direction = "vertical")

# (d) growth_plus_group: the accumulator age holds every older fish, so its mean
# size is a mixture of them rather than the curve's value there. An accumulator
# age well short of the asymptote is where the two answers separate
acc <- 8
crv_acc <- get_laa_curve(0:acc, L0, L1, L2, K, CV1, CV2, A1, A2)
pg <- crv_acc$L; pg[acc + 1] <- plus_group_size(crv_acc$L[acc + 1], crv_acc$Linf, acc)
pg_df <- bind_rows(
  data.frame(Age = 0:acc, L = crv_acc$L, Option = 'growth_plus_group = "curve"'),
  data.frame(Age = 0:acc, L = pg, Option = 'growth_plus_group = "mixture"')
)
p_pg <- ggplot(pg_df, aes(Age, L, color = Option, linetype = Option)) +
  geom_hline(yintercept = crv_acc$Linf, linetype = "dashed", color = "gray30") +
  annotate("text", x = 1.4, y = crv_acc$Linf - 3, label = "L[infinity]", parse = TRUE, size = 5) +
  geom_line(linewidth = 1.1) +
  geom_point(data = subset(pg_df, Age == acc), size = 3.4) +
  annotate("segment", x = acc - 0.9, xend = acc - 0.9, y = crv_acc$L[acc + 1], yend = pg[acc + 1],
           arrow = grid::arrow(ends = "both", length = grid::unit(0.14, "cm")), color = "gray20") +
  annotate("text", x = acc - 1.2, y = mean(c(crv_acc$L[acc + 1], pg[acc + 1])),
           label = sprintf("%.1f", pg[acc + 1] - crv_acc$L[acc + 1]), hjust = 1, size = 4.6) +
  ggthemes::scale_color_colorblind() +
  labs(x = "Age (8 accumulates)", y = "Mean length", color = "Plus group", linetype = "Plus group") +
  thm + theme(legend.direction = "vertical")

ggsave(file.path(fig_dir, "af_growth_curve.png"),
       patchwork::wrap_plots(p_anat, p_rho, p_a2, p_pg, ncol = 2),
       width = 14, height = 11, dpi = 150)

## 2. The spread and the key --------------------------------------------------
# The curve gives a mean length at age; the age-length key needs a spread around
# it, and these three switches are what set it.

# (a) growth_cv_type: the CV runs from CV1 to CV2 between the reference ages,
# interpolated either on mean length or on age
cv_df <- bind_rows(
  data.frame(Age = ages, cv = crv$cv, Option = 'growth_cv_type = "len"'),
  data.frame(Age = ages, cv = get_laa_curve(ages, L0, L1, L2, K, CV1, CV2, A1, A2, cv_type = 1)$cv,
             Option = 'growth_cv_type = "age"')
)
p_cv <- ggplot(cv_df, aes(Age, cv, color = Option, linetype = Option)) +
  geom_line(linewidth = 1.1) + ggthemes::scale_color_colorblind() +
  labs(x = "Age", y = "CV of length at age", color = NULL, linetype = NULL) +
  thm + theme(legend.direction = "vertical")

# (b) growth_sd_type: the same two parameters read as CVs, which scale the mean,
# or as standard deviations, which do not
sd_df <- bind_rows(
  data.frame(Age = ages, sd = crv$sd, Option = 'growth_sd_type = "cv"'),
  data.frame(Age = ages, sd = get_laa_curve(ages, L0, L1, L2, K, CV1 = 3, CV2 = 5, A1, A2, sd_type = 1)$sd,
             Option = 'growth_sd_type = "sd"')
)
p_sd <- ggplot(sd_df, aes(Age, sd, color = Option, linetype = Option)) +
  geom_line(linewidth = 1.1) + ggthemes::scale_color_colorblind() +
  labs(x = "Age", y = "SD of length at age", color = NULL, linetype = NULL) +
  thm + theme(legend.direction = "vertical")

# (c) growth_dist: the distribution the key bins. The two forms are close at a
# modest CV, so this panel is drawn at a wide one, which is where the lognormal's
# skew and its floor at zero actually change the key.
#
# The module always hands get_alk() the spread that growth_sd_type produced,
# sd = cv * L under "cv" and the parameter itself under "sd". A log-scale
# transform needs a log-scale spread, so the pairings compared here are the two
# that are dimensionally coherent: normal with "cv", lognormal with "sd". A log
# SD of 0.35 is very nearly a CV of 0.35, so the two are on equal footing.
cv_wide <- 0.35
crv_w <- get_laa_curve(ages, L0, L1, L2, K, CV1 = cv_wide, CV2 = cv_wide, A1, A2)
crv_ws <- get_laa_curve(ages, L0, L1, L2, K, CV1 = cv_wide, CV2 = cv_wide, A1, A2, sd_type = 1)
mid <- lower + 1
key_wn <- get_alk(lower, crv_w$L, crv_w$sd, dist = 0)   # sd_type = "cv"
key_wl <- get_alk(lower, crv_ws$L, crv_ws$sd, dist = 1) # sd_type = "sd", log scale
# ages young enough that the whole distribution sits inside the bins, so the
# comparison is not confused by the tail accumulating in the last bin
dist_df <- bind_rows(lapply(c(1, 4), function(a) {
  ia <- match(a, ages)
  bind_rows(data.frame(Length = mid, p = key_wn[, ia], Option = 'normal, sd_type = "cv"'),
            data.frame(Length = mid, p = key_wl[, ia], Option = 'lognormal, sd_type = "sd"')) %>%
    mutate(Age = paste("Age", a))
}))
dist_df$Age <- lvl(dist_df$Age, "Age 1", "Age 4")
dist_df$Option <- lvl(dist_df$Option, 'normal, sd_type = "cv"', 'lognormal, sd_type = "sd"')
p_dist <- ggplot(dist_df, aes(Length, p, color = Option, linetype = Option)) +
  geom_line(linewidth = 1.1) + facet_wrap(~Age, scales = "free", ncol = 1) +
  ggthemes::scale_color_colorblind() +
  labs(x = "Length", y = "Proportion at length", color = "Spread of 0.35",
       linetype = "Spread of 0.35") +
  thm + theme(legend.direction = "vertical")

# (d) the key itself, which is what the length compositions and the conditional
# age-at-length rows are read through. On a square root scale, since the age
# zero column is nearly a spike and would otherwise flatten everything else
key_n <- get_alk(lower, crv$L, crv$sd, dist = 0)
key_df <- expand.grid(Length = mid, Age = ages) %>% mutate(p = as.vector(key_n))
p_key <- ggplot(key_df, aes(Age, Length, fill = p)) +
  geom_raster() +
  scale_fill_viridis_c(option = "mako", direction = -1, trans = "sqrt",
                       breaks = c(0.01, 0.1, 0.3, 0.6)) +
  labs(x = "Age", y = "Length", fill = "P(length | age)") +
  thm + theme(legend.key.width = unit(1.5, "cm"))

ggsave(file.path(fig_dir, "af_growth_spread.png"),
       patchwork::wrap_plots(p_cv, p_sd, p_dist, p_key, ncol = 2),
       width = 14, height = 11, dpi = 150)

## 3. Time variation ----------------------------------------------------------
# The EBS Pacific cod model, which puts independent annual deviations on the
# length at age 1.5 and on K from 2000 and carries size at age forward cohort by
# cohort.

source(here("tests", "testthat", "helper-bridge_ebs_pcod.R"))
pdat <- SPoRC::sgl_rg_ebs_pcod_data
pyrs <- pdat$years; pages <- pdat$ages

p_inp <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(pdat))), pdat)
coh <- fit_model(p_inp$data, p_inp$par, p_inp$map, do_optim = FALSE, silent = TRUE)$rep

# (a) the two parameters that vary, as the assessment's own estimate realizes
# them year by year
tv_df <- bind_rows(
  data.frame(Year = pyrs, value = coh$growth_pars_y[1, 1, , 1, 1], Par = "Length at age 1.5 (cm)"),
  data.frame(Year = pyrs, value = coh$growth_pars_y[1, 1, , 3, 1], Par = "Richards K")
)
p_pars <- ggplot(tv_df, aes(Year, value)) +
  geom_line(linewidth = 1) + facet_wrap(~Par, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "Value in effect") + thm

# (b) what curve and cohort do differently. The assessment's own deviations are
# small enough that the two readings almost coincide, so this panel exaggerates
# them into a single sustained drop in K from 2000 to show the mechanism: the
# same deviations, read two ways. Not the assessment's estimate.
step_inp <- p_inp
kdev <- array(0, dim = dim(step_inp$par$ln_growth_devs))
kdev[1, 1, match(2000:2024, pyrs), 3, 1] <- -1.2 # logit link; K falls from 0.115 to 0.037
step_inp$par$ln_growth_devs <- kdev

step_coh <- fit_model(step_inp$data, step_inp$par, step_inp$map, do_optim = FALSE, silent = TRUE)$rep
step_inp$data$growth_tv_type <- 0 # read every year off that year's own curve
step_cur <- fit_model(step_inp$data, step_inp$par, step_inp$map, do_optim = FALSE, silent = TRUE)$rep

laa_df <- bind_rows(lapply(c(1999, 2003, 2012), function(y) {
  iy <- match(y, pyrs)
  bind_rows(
    data.frame(Age = pages, L = step_cur$mean_LAA_spawn[1, 1, iy, 1, , 1], Option = 'growth_tv_type = "curve"'),
    data.frame(Age = pages, L = step_coh$mean_LAA_spawn[1, 1, iy, 1, , 1], Option = 'growth_tv_type = "cohort"')
  ) %>% mutate(Year = paste("Year", y))
}))
laa_df$Year <- lvl(laa_df$Year, "Year 1999", "Year 2003", "Year 2012")
p_laa <- ggplot(laa_df, aes(Age, L, color = Option, linetype = Option)) +
  geom_line(linewidth = 1.1) + facet_wrap(~Year, nrow = 1) +
  ggthemes::scale_color_colorblind() +
  labs(x = "Age", y = "Mean length at age (cm)", color = NULL, linetype = NULL) + thm

ggsave(file.path(fig_dir, "af_growth_tv.png"),
       patchwork::wrap_plots(p_pars, p_laa, ncol = 1, heights = c(0.75, 1)),
       width = 13, height = 10, dpi = 150)

## 4. The semi-parametric surface ---------------------------------------------
# Simulated data whose mean length at age moves by a known year-by-age surface.
# One model estimates the parametric curve alone, the other adds a 2D AR(1)
# surface on top of it, and both see the same length compositions and
# conditional age-at-length.

source(here("tests", "testthat", "helper-selftest_growth_semipar.R"))
sp_sim <- semipar_simulate(seed = 11)
sp_in <- suppressWarnings(semipar_input("2dar1", obs = sp_sim$obs))
sp_in0 <- suppressWarnings(semipar_input("none", obs = sp_sim$obs))
sp_fit <- fit_model(sp_in$data, sp_in$par, sp_in$map, random = NULL, silent = TRUE,
                    do_optim = TRUE, newton_loops = 2)
sp_flat <- fit_model(sp_in0$data, sp_in0$par, sp_in0$map, random = NULL, silent = TRUE,
                     do_optim = TRUE, newton_loops = 2)

sp_ages <- 1:spcfg$n_ages; sp_yrs <- 1:spcfg$n_yrs
sp_est <- sp_fit$rep$mean_LAA_srv[1, 1, , 1, , 1, 1]
sp_par <- sp_flat$rep$mean_LAA_srv[1, 1, , 1, , 1, 1]

# (a) the surface that was simulated, as a proportional departure from the curve
surf_df <- expand.grid(Year = sp_yrs, Age = sp_ages) %>% mutate(dev = as.vector(sp_sim$devs))
p_surf <- ggplot(surf_df, aes(Year, Age, fill = dev)) +
  geom_raster() +
  scale_fill_gradient2(low = "#3B4CC0", mid = "gray95", high = "#B40426") +
  scale_y_continuous(breaks = seq(2, 12, 2)) +
  labs(x = "Year", y = "Age", fill = "Departure from the curve,\nlog mean length at age") +
  thm + theme(legend.key.width = unit(1.5, "cm"))

# (b) mean length at age recovered, against the truth, in three years
sp_laa <- bind_rows(lapply(c(5, 15, 27), function(y) {
  bind_rows(
    data.frame(Age = sp_ages, L = sp_sim$mean_LAA[y, ], Option = "Simulated truth"),
    data.frame(Age = sp_ages, L = sp_par[y, ], Option = 'growth_semipar = "none"'),
    data.frame(Age = sp_ages, L = sp_est[y, ], Option = 'growth_semipar = "2dar1"')
  ) %>% mutate(Year = paste("Year", y))
}))
sp_laa$Year <- lvl(sp_laa$Year, "Year 5", "Year 15", "Year 27")
sp_laa$Option <- lvl(sp_laa$Option, "Simulated truth", 'growth_semipar = "none"', 'growth_semipar = "2dar1"')
p_sp_laa <- ggplot(sp_laa, aes(Age, L, color = Option, linetype = Option)) +
  geom_line(linewidth = 1.1) + facet_wrap(~Year, nrow = 1) +
  ggthemes::scale_color_colorblind() +
  labs(x = "Age", y = "Mean length at age", color = NULL, linetype = NULL) + thm

ggsave(file.path(fig_dir, "af_growth_semipar.png"),
       patchwork::wrap_plots(p_surf, p_sp_laa, ncol = 1, heights = c(1, 0.85)),
       width = 13, height = 10, dpi = 150)

cat("semi-parametric max relative error in mean length at age:",
    sprintf("%.2f%%", 100 * max(abs(sp_est / sp_sim$mean_LAA - 1))), "\n")
cat("parametric only:",
    sprintf("%.2f%%", 100 * max(abs(sp_par / sp_sim$mean_LAA - 1))), "\n")
cat("K in 2012 under the illustrative step:",
    sprintf("%.3f (base %.3f)", step_coh$growth_pars_y[1, 1, match(2012, pyrs), 3, 1],
            coh$growth_pars_y[1, 1, 1, 3, 1]), "\n")
