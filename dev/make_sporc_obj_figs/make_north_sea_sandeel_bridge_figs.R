# Purpose: Bridge the smsR assessment for North Sea sandeel in area 1r to SPoRC
#          and render the case study figures. The model is specified the way
#          smsR specifies it: one region and one sex, ages 0-4, two half-year
#          seasons with age 0 recruiting into the second, mean recruitment with
#          free initial numbers at age, and fishing mortality driven by the
#          observed effort series rather than estimated year by year.
# Creator: Matthew LH. Cheng
# Date Created: 8/26/26

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

bundle <- readRDS(here("dev", "dev_data", "north_sea_sandeel_1r.rds"))
dat <- bundle$data
ref <- bundle$reference
lhs <- dat$lhs                               # weights, maturity and natural mortality

# Dimensions -----------------------------------------------------------------
yrs <- 1983:2021
ages <- 0:4
n_yrs <- length(yrs)
n_ages <- length(ages)
n_seas <- 2
n_regions <- 1
n_pop <- 1
n_sexes <- 1
n_fish <- 2                                  # one fishery fleet per season
n_srv <- 2                                   # the Dredge and the RTM

# The Dredge runs in season 2 over ages 0-1, the RTM in season 1 over ages 1-3,
# each entering part way through its season.
srv_ages <- list(0:1, 1:3)
srv_seas <- c(2, 1)
srv_t <- c(0.75, 0)

# Biologicals ----------------------------------------------------------------
# smsR carries its biologicals as [age, year, season]. SPoRC wants
# [pop, region, year, season, age, sex], the same numbers with three dimensions
# of length one, and one more for fleet on the weights a fleet sees.
WAA      <- array(0, dim = c(1, 1, n_yrs, n_seas, n_ages, 1))
MatAA    <- array(0, dim = c(1, 1, n_yrs, n_seas, n_ages, 1))
WAA_fish <- array(0, dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 2))
WAA_srv  <- array(0, dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 2))

for(y in 1:n_yrs) {
  for(seas in 1:n_seas) {
    WAA[1, 1, y, seas, , 1]   = lhs$west[, y, seas]   # stock weight at age
    MatAA[1, 1, y, seas, , 1] = lhs$mat[, y, seas]    # proportion mature at age
    for(f in 1:2) WAA_fish[1, 1, y, seas, , 1, f] = lhs$weca[, y, seas]  # catch weight
    for(k in 1:2) WAA_srv[1, 1, y, seas, , 1, k]  = lhs$west[, y, seas]
  } # end seas loop
} # end y loop

# smsR carries natural mortality by season; SPoRC scales an annual rate by season
# duration, so the annual sum preserves the annual trajectory exactly.
M_sms <- lhs$M[, 1:n_yrs, ]
natmort <- array(0, dim = c(1, 1, n_yrs, n_ages, 1))
for(y in 1:n_yrs) {
  annual_M = M_sms[, y, 1] + M_sms[, y, 2]
  annual_M[1] = M_sms[1, y, 2] / 0.5   # age 0 is only present in season 2
  natmort[1, 1, y, , 1] = annual_M
} # end y loop

# Catch at age ---------------------------------------------------------------
# Ages 1-4 only, with nocatch naming the year and season cells that were fished.
# Fleet f fishes season f alone, so a season needing its own observation error
# gets it through the fleet rather than through an extra key dimension.
catch_obs <- dat$Catch
nocatch <- as.matrix(dat$nocatch)
ObsCatchAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))
UseCatchAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))
for(y in 1:n_yrs) {
  for(s in 1:n_seas) {
    for(a in 2:n_ages) {
      obs <- catch_obs[a, y, s]
      if(!is.na(obs) && obs > 0 && nocatch[y, s] == 1) {
        ObsCatchAA[1, y, s, a, 1, s] <- obs
        UseCatchAA[1, y, s, a, 1, s] <- 1
      } # end fished cell
    } # end a loop
  } # end s loop
} # end y loop

# Survey index at age --------------------------------------------------------
# The aggregated survey stream stays empty: each survey is fit age by age.
survey_obs <- dat$survey
ObsSrvIdx <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_srv))
UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
ObsSrvIdxAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))
UseSrvIdxAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))
for(k in 1:n_srv) {
  for(a in srv_ages[[k]]) {
    age_row <- which(ages == a)
    for(y in 1:n_yrs) {
      obs <- survey_obs[age_row, y, k]
      if(!is.na(obs) && obs > 0) {
        ObsSrvIdxAA[1, y, srv_seas[k], age_row, 1, k] <- obs
        UseSrvIdxAA[1, y, srv_seas[k], age_row, 1, k] <- 1
      } # end observed cell
    } # end y loop
  } # end a loop
} # end k loop

# Observation error keyed by age, as smsR groups it. Catchability at age is not
# a separate parameter: it is the survey selectivity, through the "nonparfree"
# form, whose values carry the height of the curve as well as its shape.
srv_sd_key <- array(NA_integer_, dim = c(n_ages, n_sexes, n_srv))
# The Dredge's age-0 standard deviation is held rather than estimated. Age 0
# appears in one place only, that index, and the recruitment deviations are free,
# so they can fit it exactly: the residual goes to zero, log sigma runs to
# negative infinity, and the likelihood is unbounded. The optimizer reports
# convergence from inside that hole, at an objective well below the real one.
# Held at smsR's own value, every starting value from -5 to -15 reaches the same
# optimum; estimated, only a narrow band of starting values does.
srv_sd_key[1, 1, 1] <- NA                    # Dredge age 0, held (see above)
srv_sd_key[2, 1, 1] <- 1                     # Dredge age 1
srv_sd_key[2, 1, 2] <- 2                     # RTM age 1
srv_sd_key[3:4, 1, 2] <- 3                   # RTM ages 2 and 3 share

srv_sd_start <- array(log(0.5), dim = c(n_ages, n_sexes, n_srv))
srv_sd_start[1, 1, 1] <- log(0.4052)         # the held value, smsR's own

# smsR groups its catch standard deviations as ages 1-2 and ages 3 and above,
# separately by season. With a fleet per season that is the plain [age, fleet] key.
sdc_key <- array(c(NA, 1L, 1L, 2L, 2L,
                   NA, 3L, 3L, 4L, 4L), dim = c(n_ages, n_sexes, n_fish))

# Setup ----------------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  n_pop = n_pop, years = yrs, ages = ages, lens = NA,
  n_regions = n_regions, n_sexes = n_sexes, n_seas = n_seas, seasdur = c(0.5, 0.5),
  n_fish_fleets = n_fish, n_srv_fleets = n_srv, verbose = FALSE)

# Recruitment enters in season 2, which is what fixed_rec_seas_prop states. smsR
# estimates the first year's numbers at age as free parameters rather than as
# departures from an equilibrium, which is init_age_strc "free": ln_InitDevs are
# then the initial log numbers themselves. It applies no penalty to them, and
# because those deviations are log-numbers rather than log-ratios a penalty here
# would act as a prior on initial abundance, so equil_init_age_strc is "equil".
input_list <- Setup_Mod_Rec(
  input_list, rec_model = "mean_rec", rec_lag = 1,
  spawn_seas = 1, t_spawn = 0, use_fixed_rec_seas_prop = 1,
  fixed_rec_seas_prop = matrix(c(0, 1), nrow = n_pop, ncol = n_seas),
  sigmaR_spec = "fix", do_rec_bias_ramp = 0,
  init_age_strc = "free", equil_init_age_strc = "equil",
  ln_global_R0 = log(2e8))
input_list$par$ln_sigmaR[] <- log(4.0128)    # nllfactor 0.05 at smsR's sigmaR 0.897

input_list <- Setup_Mod_Biologicals(
  input_list, WAA = WAA, WAA_fish = WAA_fish, WAA_srv = WAA_srv,
  MatAA = MatAA * 2,                           # cancels SPoRC's single-sex SSB halving
  fit_lengths = 0, M_spec = "fix", Fixed_natmort = natmort)

input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                 Fixed_Movement = NA, do_recruits_move = 0)

# The aggregated catch stream is switched off entirely: this fleet fits catch at
# age, and supplying both is an error rather than a warning.
suppressWarnings(
  input_list <- Setup_Mod_Catch_and_F(
    input_list,
    ObsCatch = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    UseCatch = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    ObsCatchAA = ObsCatchAA, UseCatchAA = UseCatchAA,
    sigmaCAA_key = sdc_key, sigmaCAA_spec = "est",
    ln_sigmaCAA = array(log(0.5), dim = c(n_ages, n_sexes, n_fish)),
    catch_units = array("abd", dim = c(n_fish)),
    Use_F_pen = 0, sigmaC_spec = "fix", sigmaF_spec = "fix"))

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list,
  ObsFishIdx = array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  UseFishIdx = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ObsFishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish)),
  UseFishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ISS_FishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  ObsFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_fish)),
  UseFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ISS_FishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  fish_idx_type = c("none", "none"),
  FishAgeComps_LikeType = c("none", "none"),
  FishLenComps_LikeType = c("none", "none"),
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2"),
  FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2"))

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list,
  ObsSrvIdx = ObsSrvIdx, ObsSrvIdx_SE = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
  UseSrvIdx = UseSrvIdx,
  ObsSrvIdxAA = ObsSrvIdxAA, UseSrvIdxAA = UseSrvIdxAA,
  sigmaSrvIdxAA_key = srv_sd_key, sigmaSrvIdxAA_spec = "est",
  ln_sigmaSrvIdxAA = srv_sd_start,
  ObsSrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv)),
  UseSrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
  ISS_SrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
  ObsSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_srv)),
  UseSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
  ISS_SrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
  srv_idx_type = c("abd", "abd"),
  SrvAgeComps_LikeType = c("none", "none"),
  SrvLenComps_LikeType = c("none", "none"),
  SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2"),
  SrvLenComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2"))

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list, cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
  fish_sel_blocks = c("Block_1_Year_1-16_Fleet_1", "Block_2_Year_17-terminal_Fleet_1",
                      "Block_1_Year_1-16_Fleet_2", "Block_2_Year_17-terminal_Fleet_2"),
  fish_sel_model = c("nonparfree_Fleet_1", "nonparfree_Fleet_2"),
  fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  fish_sel_nonpar_est_bins = rep(list(list(list(1, 2, 3, c(4, 5)), list(1, 2, 3, c(4, 5)))), 2),
  fish_fixed_sel_pars_spec = c("est_all", "est_all"),
  fish_q_spec = c("fix", "fix"))

t_srv <- array(0, dim = c(n_regions, n_seas, n_srv))
for(k in 1:n_srv) t_srv[1, srv_seas[k], k] <- srv_t[k]
input_list <- Setup_Mod_Srvsel_and_Q(
  input_list, cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2"),
  srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  srv_sel_model = c("nonparfree_Fleet_1", "nonparfree_Fleet_2"),
  srv_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  # one value per age each survey observes: the Dredge ages 0-1, the RTM ages
  # 1-3. Ages left out get no parameter, which is what an unobserved age wants.
  srv_sel_nonpar_est_bins = list(list(list(1, 2)), list(list(2, 3, 4))),
  srv_fixed_sel_pars_spec = c("est_all", "est_all"),
  srv_q_spec = c("est_all", "est_all"), t_srv = t_srv)
input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

input_list <- Setup_Mod_Weighting(
  input_list, Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  Wt_FishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  Wt_SrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
  Wt_SrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)))

# Fishery selectivity is one value per age, per block, shared by the two season
# fleets because smsR's age pattern does not vary by season. Three things are not
# estimated: age 0, which is never fished, and block one's top group, which is the
# reference the rest of the curve is measured against. Holding that reference at
# one is what keeps the level of fishing in ln_F_mean rather than splitting it
# between the two.
#
#             age 0    age 1    age 2    age 3    age 4
#   block 1     .        1        2      <-- reference, held at 1 -->
#   block 2     .        3        4        5        5
#
sel_par <- rbind(block_1 = c(NA, 1L, 2L, NA, NA),
                 block_2 = c(NA, 3L, 4L, 5L, 5L))

# fish_fixed_sel_pars is [region, age, block, sex, fleet]. Under "nonparfree" a
# value of 0 is a selectivity of one, and -10 is a selectivity of zero.
sel_start <- input_list$par$fish_fixed_sel_pars
sel_map <- array(NA_integer_, dim = dim(sel_start))
sel_start[1, , , 1, ] <- 0
sel_start[1, 1, , 1, ] <- -10                # age 0, never fished
for(b in 1:2) {
  for(f in 1:n_fish) sel_map[1, , b, 1, f] <- sel_par[b, ]
} # end b loop

input_list$par$fish_fixed_sel_pars <- sel_start
input_list$map$fish_fixed_sel_pars <- factor(sel_map)

# Effort-driven fishing mortality, and the free initial numbers at age. An at-age
# catch likelihood with free annual fishing mortality is close to saturated, so
# the effort series is what closes that ridge.
effort <- dat$effort
dev_dim <- dim(input_list$par$ln_F_devs)
ln_F_devs <- array(0, dim = dev_dim)
for(y in 1:dev_dim[2]) {
  for(s in 1:dev_dim[3]) ln_F_devs[1, y, s, ] <- log(effort[y, s])
} # end y loop
input_list$par$ln_F_devs <- ln_F_devs
input_list$map$ln_F_devs <- factor(array(NA_integer_, dim = dev_dim))
input_list$par$ln_F_mean[] <- -5
# One initial numbers-at-age parameter per age, all estimated, and this has to be
# said by hand. equil_init_age_strc "equil" neither penalizes nor estimates the
# deviations, and "stoch_all" does both; under init_age_strc "free" what is
# wanted is the combination neither offers, estimated and unpenalized, because
# the deviations are the numbers themselves rather than departures from an
# equilibrium. Setting the map here is one way; init_devs_pen_use[] <- 0 with
# "stoch_all" is the other, and gives the same fit. Either has to come after the
# setup calls, which would otherwise overwrite it.
input_list$map$ln_InitDevs <- factor(seq_along(input_list$par$ln_InitDevs))

# Reference series -----------------------------------------------------------
# smsR carries fishing mortality and numbers by season, so they are summed or
# read at season one to line up with what SPoRC reports.
ref_F <- t(apply(ref$F0[, 1:n_yrs, ], c(1, 2), sum))
ref_N <- t(ref$Nsave[, 1:n_yrs, 1])
qs <- c("Recruitment", "Spawning Biomass", "Fbar (ages 1-2)", "Numbers at Age 4")

# The two series overplot at bridge accuracy, so the percent difference gets its
# own column and its own axis.
sandeel_fig <- function(rep_obj, file) {
  sp_N <- rep_obj$NAA[1,1,1:n_yrs,1,,1]
  sp_F <- apply(rep_obj$tot_FAA[1,1,1:n_yrs,,,1,], c(1,3), sum)
  d <- bind_rows(
    data.frame(Year = yrs, smsR = ref$Rsave[1:n_yrs],    SPoRC = as.numeric(rep_obj$Rec)[1:n_yrs], q = qs[1]),
    data.frame(Year = yrs, smsR = ref$SSB[1:n_yrs],      SPoRC = as.numeric(rep_obj$SSB)[1:n_yrs], q = qs[2]),
    data.frame(Year = yrs, smsR = rowMeans(ref_F[,2:3]), SPoRC = rowMeans(sp_F[,2:3]),             q = qs[3]),
    data.frame(Year = yrs, smsR = ref_N[,5],             SPoRC = sp_N[,5],                         q = qs[4])
  ) %>% mutate(q = factor(q, levels = qs))

  lvl <- d %>% pivot_longer(c(smsR, SPoRC), names_to = "Type", values_to = "Value") %>%
    mutate(Type = factor(ifelse(Type == "smsR", "smsR (sandeel area 1r)", "SPoRC"),
                         levels = c("smsR (sandeel area 1r)", "SPoRC")))
  cols <- c("smsR (sandeel area 1r)" = "black", "SPoRC" = "#E69F00")
  base <- theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  p1 <- ggplot(lvl, aes(Year, Value, color = Type)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.1) +
    facet_wrap(~q, scales = "free_y", ncol = 1) +
    scale_color_manual(values = cols) + scale_y_continuous(labels = scales::comma) +
    labs(x = "Year", y = "Value", color = "Type") + base
  p2 <- ggplot(d, aes(Year, 100 * (SPoRC - smsR) / smsR)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_line(color = "#E69F00", linewidth = 0.6) + geom_point(color = "#E69F00", size = 1.1) +
    facet_wrap(~q, scales = "free_y", ncol = 1) +
    labs(x = "Year", y = "SPoRC vs smsR (%)") + base

  # file = NULL reports the comparison without writing a figure, which is what
  # the seeded stage wants: it is a check, not something the vignette shows.
  if(!is.null(file)) ggplot2::ggsave(here("vignettes", "figures", file), p1 | p2,
                                     width = 13, height = 10, dpi = 150)
  bind_rows(lapply(qs, function(k) {
    rows <- d[d$q == k, ]
    bridge_cmp(k, rows$SPoRC, rows$smsR)
  }))
} # end sandeel_fig

# Stage one: seeded at smsR's maximum likelihood estimate ---------------------
F0 <- ref$F0[, 1:n_yrs, ]
blk <- c(rep(1L, 16), rep(2L, 23))           # matches fish_sel_blocks
seeded <- input_list

pat <- sapply(1:2, function(b) {
  y <- which(blk == b & F0[3, , 1] > 1e-10)[1]
  F0[, y, 1] / max(F0[, y, 1])
})
fsp_seed <- seeded$par$fish_fixed_sel_pars
for(b in 1:2) for(f in 1:n_fish) {
  fsp_seed[1, 1, b, 1, f] <- log(1e-12)     # age 0 unfished, to machine zero
  fsp_seed[1, 2:5, b, 1, f] <- log(pmax(pat[2:5, b], 1e-12) / pat[4, 1])
} # end b loop
seeded$par$fish_fixed_sel_pars <- fsp_seed

seeded$par$ln_F_mean[] <- 0
seeded$map$ln_F_mean <- factor(array(NA_integer_, dim = dim(seeded$par$ln_F_mean)))
seeded$par$ln_F_devs[] <- log(1e-12)
seeded$map$ln_F_devs <- factor(array(NA_integer_, dim = dim(seeded$par$ln_F_devs)))
probe <- fit_model(seeded$data, seeded$par, seeded$map, do_optim = FALSE, silent = TRUE)
selmax <- apply(probe$rep$fish_sel[1, 1, , 1, , 1, 1], 1, max)

lnFd <- array(log(1e-12), dim = dim(seeded$par$ln_F_devs))
for(y in 1:n_yrs) for(s in 1:n_seas) {
  mx <- max(F0[, y, s])
  lnFd[1, y, s, s] <- if(mx > 1e-10) log(mx / selmax[y]) else log(1e-12)
} # end y loop
seeded$par$ln_F_devs <- lnFd
seeded$map$ln_F_devs <- factor(array(NA_integer_, dim = dim(lnFd)))

lnR0 <- mean(log(ref$Rsave[1:n_yrs]))
seeded$par$ln_global_R0 <- lnR0
seeded$map$ln_global_R0 <- factor(NA)
seeded$par$ln_RecDevs[] <- log(ref$Rsave[1:n_yrs]) - lnR0
seeded$map$ln_RecDevs <- factor(rep(NA, length(seeded$par$ln_RecDevs)))

# solve the initial deviations so year one matches smsR exactly
seeded$par$ln_InitDevs[] <- 0
s0 <- fit_model(seeded$data, seeded$par, seeded$map, do_optim = FALSE, silent = TRUE)
base_N <- s0$rep$NAA[1, 1, 1, 1, , 1]
seeded$par$ln_InitDevs[] <- log(ref$Nsave[, 1, 1][2:5] / base_N[2:5])

seed <- fit_model(seeded$data, seeded$par, seeded$map, do_optim = FALSE, silent = TRUE)

sp_F <- apply(seed$rep$tot_FAA[1,1,1:n_yrs,,,1,], c(1,2,3), sum)
sms_F <- aperm(F0, c(2,3,1))
rel <- abs(sp_F - sms_F) / pmax(sms_F, 1e-12)
rel[sms_F < 1e-10] <- 0
message("seeded at smsR's MLE, not optimized")
message("  fishing mortality at age, every cell: ", signif(100 * max(rel), 3), " %")
print(sandeel_fig(seed$rep, NULL))

# Stage two: fit freely -------------------------------------------------------
est <- fit_model(input_list$data, input_list$par, input_list$map,
                 random = NULL, newton_loops = 3, silent = TRUE)
for(i in 1:5) est$optim <- stats::nlminb(est$optim$par, est$fn, est$gr,
                    control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
message("\nfree fit: ", length(est$optim$par), " parameters, jnLL ",
        round(est$optim$objective, 3), ", max |gradient| ",
        signif(max(abs(est$gr(est$optim$par))), 2))
r <- est$report(est$env$last.par.best)
print(sandeel_fig(r, "ah_sandeel_ts_refit.png"))

sp_N <- r$NAA[1, 1, 1:n_yrs, 1, , 1]
message("\nnumbers at age, free fit")
for(a in 2:5) message("  age ", a - 1,
                      "  max ", round(max(abs(100*(sp_N[,a]-ref_N[,a])/ref_N[,a])), 2),
                      " %  cor ", round(stats::cor(log(ref_N[,a]), log(sp_N[,a])), 4))
