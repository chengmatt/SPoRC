# Purpose: To demonstrate the collapsibility of SPoRC using Alaska sablefish
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 6/11/25


# Setup -------------------------------------------------------------------

library(here)
library(ggplot2)
library(tidyverse)
library(SPoRC)
library(geomtextpath)
devtools::load_all(here("R"))

# Read in model outputs
five_dat <- readRDS(here("dev", "dev_output", "5_Region_Model_Sablefish", "data.RDS"))
three_dat <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "data.RDS"))
sgl_dat <- readRDS(here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "data.RDS"))

five_rep <- readRDS(here("dev", "dev_output", "5_Region_Model_Sablefish", "rep.RDS"))
three_rep <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "rep.RDS"))
sgl_rep <- readRDS(here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "rep.RDS"))

five_sdrep <- readRDS(here("dev", "dev_output", "5_Region_Model_Sablefish", "sd_rep.RDS"))
three_sdrep <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "sd_rep.RDS"))
sgl_sdrep <- readRDS(here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "sd_rep.RDS"))

# Projections -------------------------------------------------------------

# Define HCR to use
HCR_function <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp # define stock status
  # If stock status is > 1
  if(stock_status >= 1) f <- frp
  # If stock status is between brp and alpha
  if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
  # If stock status is less than alpha
  if(stock_status < alpha) f <- 0
  return(f)
}

# single area model
sgl_ref_pt <- Get_Reference_Points(data = sgl_dat,
                                   rep = sgl_rep,
                                   SPR_x = 0.4,
                                   type = 'single_region',
                                   what = 'SPR',
                                   calc_rec_st_yr = 20,
                                   rec_age = 2)

# set up quantities to use in projection function
n_sims <- 1
t_spawn <- 0 # spawn timing
n_proj_yrs <- 30 # number of projection years
n_regions <- 1 # number of regions
n_ages <- length(sgl_dat$ages) # number of ages
n_sexes <- sgl_dat$n_sexes # number of sexes
n_fish_fleets <- sgl_dat$n_fish_fleets # number of fishery fleets
n_seas <- 1
n_pop <- 1
do_recruits_move <- 0 # recruits don't move
terminal_NAA <- array(sgl_rep$NAA[,,length(sgl_dat$years),,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
terminal_NAA0 <- array(sgl_rep$NAA0[,,length(sgl_dat$years),,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
WAA <- array(rep(sgl_dat$WAA[,,length(sgl_dat$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(sgl_dat$WAA[,,length(sgl_dat$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # weight at age fishery
MatAA <- array(rep(sgl_dat$MatAA[,,length(sgl_dat$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(sgl_rep$fish_sel[,length(sgl_dat$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- array(rep(sgl_rep$Movement[,,,length(sgl_dat$years),,,], each = n_proj_yrs), dim = c(n_pop, n_regions, n_regions, n_seas, n_proj_yrs, n_ages, n_sexes)) # movement
terminal_F <- array(sgl_rep$Fmort[,length(sgl_dat$years),,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F
natmort <- array(sgl_rep$natmort[,,length(sgl_dat$years),,], dim = c(n_pop, n_regions, n_proj_yrs, n_ages, n_sexes)) # natural mortality
recruitment <- array(sgl_rep$Rec[,,20:(length(sgl_dat$years) - 2)], dim = c(n_pop, n_regions, length(20:(length(sgl_dat$years) - 2)))) # recruitment values to use for mean recruitment calculations or inverse gaussian parameterization
sexratio <- array(0.5, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)) # recruitment sex ratio

# storage
sgl_f_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
sgl_ssb_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
sgl_catch_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims))

# do population projection
for(sim in 1:n_sims) {

  # do projection
  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  sexratio = sexratio,
                                  n_fish_fleets = n_fish_fleets,
                                  do_recruits_move = do_recruits_move,
                                  recruitment = recruitment,
                                  terminal_NAA = terminal_NAA,
                                  terminal_NAA0 = terminal_NAA0,
                                  terminal_F = terminal_F,
                                  natmort = natmort,
                                  WAA = WAA,
                                  n_pop = n_pop,
                                  WAA_fish = WAA_fish,
                                  MatAA = MatAA,
                                  fish_sel = fish_sel,
                                  Movement = Movement,
                                  f_ref_pt = array(sgl_ref_pt$f_ref_pt, dim = c(sgl_dat$n_regions, n_proj_yrs)),
                                  b_ref_pt = array(sgl_ref_pt$b_ref_pt, dim = c(sgl_dat$n_pop, sgl_dat$n_regions, n_proj_yrs)),
                                  HCR_function = HCR_function,
                                  recruitment_opt = "mean_rec",
                                  fmort_opt = "Input",
                                  t_spawn = t_spawn
  )

  sgl_ssb_proj[,,sim] <- out$proj_SSB
  sgl_catch_proj[,,,sim] <- out$proj_Catch
  sgl_f_proj[,,sim] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used

}

# Three area model
three_ref_pt <- Get_Reference_Points(data = three_dat,
                                     rep = three_rep,
                                     SPR_x = 0.4,
                                     type = 'multi_region',
                                     what = 'global_SPR',
                                     calc_rec_st_yr = 20,
                                     rec_age = 2)

# quantities to use in projection
n_sims <- 1
t_spawn <- 0
sexratio <- 0.5
n_proj_yrs <- 30
n_regions <- 3
n_ages <- length(three_dat$ages)
n_sexes <- three_dat$n_sexes
n_fish_fleets <- 2
do_recruits_move <- 0
n_seas <- 1
terminal_NAA <- array(three_rep$NAA[,,length(three_dat$years),,,], dim = c(1, n_regions, n_seas, n_ages, n_sexes))
terminal_NAA0 <- array(three_rep$NAA0[,,length(three_dat$years),,,], dim = c(1, n_regions, n_seas, n_ages, n_sexes))
WAA <- array(rep(three_dat$WAA[,,length(three_dat$years),,,], each = n_proj_yrs), dim = c(1, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(three_dat$WAA[,,length(sgl_dat$years),,,], each = n_proj_yrs), dim = c(1, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # weight at age
MatAA <- array(rep(three_dat$MatAA[,,length(three_dat$years),,,], each = n_proj_yrs), dim = c(1, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(three_rep$fish_sel[,length(three_dat$years),,,], each = n_proj_yrs), dim = c( n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- abind::abind(replicate(n_proj_yrs,  three_rep$Movement[,,,length(three_dat$years),,,,drop = FALSE],  simplify = FALSE), along = 4) # movement
terminal_F <- array(three_rep$Fmort[,length(three_dat$years),,], dim = c(n_regions, n_seas, n_fish_fleets))
natmort <- array(three_rep$natmort[,,length(three_dat$years),,], dim = c(1,n_regions, n_proj_yrs, n_ages, n_sexes))
recruitment <- array(three_rep$Rec[,,20:(length(three_dat$years)  - 2)], dim = c(1,n_regions, length(20:(length(three_dat$years) - 2))))
sexratio <- array(0.5, dim = c(1, n_regions, n_proj_yrs, n_sexes))

# storage
three_f_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
three_ssb_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
three_catch_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims))

bh_rec_opt <- list(
  rec_dd = 0,
  rec_lag = three_dat$rec_lag,
  R0 = three_rep$R0,
  h = array(three_rep$h_trans, dim = c(1, 3)),
  rec_region_prop = three_rep$rec_region_prop,
  WAA = array(three_dat$WAA[,,1,,,1], dim = c(1,3,1,30)),
  MatAA = array(three_dat$MatAA[,,1,,,1], dim = c(1,3,1,30)),
  SSB = three_rep$SSB,
  Movement = array(Movement[,,,1,,,1], dim = c(1, 3, 3, 1, 30)),
  do_recruits_move = three_dat$do_recruits_move,
  # t_spawn = three_dat$t_spawn,
  sex_ratio_f = array(0.5, dim = c(1, n_regions)),
  sgl_seas_spawning_movement = three_dat$sgl_seas_spawning_movement[,,,1,,],
  stray_rate = array(1, dim = c(1)),
  # rec_seas_prop = rec_seas_prop,
  natmort = array(natmort[,,1,,1], dim = c(1, 3, 30))

)

# do population projection
for(sim in 1:n_sims) {

  # do projection
  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                         n_regions = n_regions,
                                         n_ages = n_ages,
                                         n_sexes = n_sexes,
                                         sexratio = sexratio,
                                         n_fish_fleets = n_fish_fleets,
                                         do_recruits_move = do_recruits_move,
                                         recruitment = recruitment,
                                         terminal_NAA = terminal_NAA,
                                         terminal_NAA0 = terminal_NAA0,
                                         terminal_F = terminal_F,
                                         natmort = natmort,
                                         WAA = WAA,
                                         n_pop = 1,
                                         WAA_fish = WAA_fish,
                                         MatAA = MatAA,
                                         fish_sel = fish_sel,
                                         Movement = Movement,
                                         f_ref_pt = array(three_ref_pt$f_ref_pt, dim = c(three_dat$n_regions, n_proj_yrs)),
                                         b_ref_pt = array(three_ref_pt$b_ref_pt, dim = c(1, three_dat$n_regions, n_proj_yrs)),
                                         HCR_function = HCR_function,
                                         recruitment_opt = "mean_rec",
                                         fmort_opt = "HCR",
                                         t_spawn = t_spawn,
                                         bh_rec_opt = bh_rec_opt,
                                         rec_seas_prop = array(1, dim = c(1,1))
  )

  three_ssb_proj[,,sim] <- out$proj_SSB
  three_catch_proj[,,,sim] <- out$proj_Catch
  three_f_proj[,,sim] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used

}


# five area model
five_ref_pt <- Get_Reference_Points(data = five_dat,
                                    rep = five_rep,
                                    SPR_x = 0.4,
                                    type = 'multi_region',
                                    what = 'global_SPR',
                                    calc_rec_st_yr = 20,
                                    rec_age = 2)

# quantities to use in projection
n_sims <- 1
t_spawn <- 0
sexratio <- 0.5
n_proj_yrs <- 30
n_regions <- 5
n_ages <- length(five_dat$ages)
n_sexes <- five_dat$n_sexes
n_fish_fleets <- 2
do_recruits_move <- 0
n_seas <- 1
terminal_NAA <- array(five_rep$NAA[,,length(five_dat$years),,,], dim = c(1, n_regions, n_seas, n_ages, n_sexes))
terminal_NAA0 <- array(five_rep$NAA0[,,length(five_dat$years),,,], dim = c(1, n_regions, n_seas, n_ages, n_sexes))
WAA <- array(rep(five_dat$WAA[,,length(five_dat$years),,,], each = n_proj_yrs), dim = c(1, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(five_dat$WAA[,,length(sgl_dat$years),,,], each = n_proj_yrs), dim = c(1, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # weight at age
MatAA <- array(rep(five_dat$MatAA[,,length(five_dat$years),,,], each = n_proj_yrs), dim = c(1, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(five_rep$fish_sel[,length(five_dat$years),,,], each = n_proj_yrs), dim = c( n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- abind::abind(replicate(n_proj_yrs,  five_rep$Movement[,,,length(five_dat$years),,,,drop = FALSE],  simplify = FALSE), along = 4) # movement
terminal_F <- array(five_rep$Fmort[,length(five_dat$years),,], dim = c(n_regions, n_seas, n_fish_fleets))
natmort <- array(five_rep$natmort[,,length(five_dat$years),,], dim = c(1,n_regions, n_proj_yrs, n_ages, n_sexes))
recruitment <- array(five_rep$Rec[,,20:(length(five_dat$years)  - 2 )], dim = c(1,n_regions, length(20:(length(five_dat$years) - 2))))
sexratio <- array(0.5, dim = c(1, n_regions, n_proj_yrs, n_sexes))

# storage
five_f_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
five_ssb_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
five_catch_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims))

# do population projection
for(sim in 1:n_sims) {

  # do projection
  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                         n_regions = n_regions,
                                         n_ages = n_ages,
                                         n_sexes = n_sexes,
                                         sexratio = sexratio,
                                         n_fish_fleets = n_fish_fleets,
                                         do_recruits_move = do_recruits_move,
                                         recruitment = recruitment,
                                         terminal_NAA = terminal_NAA,
                                         terminal_NAA0 = terminal_NAA0,
                                         terminal_F = terminal_F,
                                         natmort = natmort,
                                         WAA = WAA,
                                         WAA_fish = WAA_fish,
                                         MatAA = MatAA,
                                         n_pop = n_pop,
                                         fish_sel = fish_sel,
                                         Movement = Movement,
                                         f_ref_pt = array(five_ref_pt$f_ref_pt, dim = c(five_dat$n_regions, n_proj_yrs)),
                                         b_ref_pt = array(five_ref_pt$b_ref_pt, dim = c(1, five_dat$n_regions, n_proj_yrs)),
                                         HCR_function = HCR_function,
                                         recruitment_opt = "mean_rec",
                                         fmort_opt = "HCR",
                                         t_spawn = t_spawn
  )

  five_ssb_proj[,,sim] <- out$proj_SSB
  five_catch_proj[,,,sim] <- out$proj_Catch
  five_f_proj[,,sim] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used

}


# Combine Projections -----------------------------------------------------

# bind reports
all_rg_ssb <- reshape2::melt(sgl_rep$Aggregated_SSB) %>%
  mutate(Type = "Single-Region",
         log_se = sgl_sdrep$sd[names(sgl_sdrep$value) == 'log_Aggregated_SSB'],
         Year = 1960:2021) %>%
  bind_rows(
    reshape2::melt(three_rep$Aggregated_SSB) %>%
      mutate(Type = 'Three-Region',
             log_se = three_sdrep$sd[names(three_sdrep$value) == 'log_Aggregated_SSB'],
             Year = 1960:2021),
    reshape2::melt(five_rep$Aggregated_SSB) %>%
      mutate(Type = 'Five-Region',
             log_se = five_sdrep$sd[names(five_sdrep$value) == 'log_Aggregated_SSB'],
             Year = 1960:2021)
  ) %>%
  mutate(Type = factor(Type, levels = c("Single-Region", "Three-Region", "Five-Region")))

# Five region SSB data
five_rg_ssb <- reshape2::melt(five_rep$SSB) %>%
  rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  left_join(data.frame(Region = 1:5, b40 = as.vector(five_ref_pt$b_ref_pt)), by = 'Region') %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS',
    Region == 2 ~ 'AI',
    Region == 3 ~ 'WGOA',
    Region == 4 ~ 'CGOA',
    Region == 5 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = five_sdrep$sd[names(five_sdrep$value) == 'log_SSB'],
         Year = Year + 1959)

# Three region SSB data
three_rg_ssb <- reshape2::melt(three_rep$SSB) %>%
  rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  left_join(data.frame(Region = 1:3, b40 = as.vector(three_ref_pt$b_ref_pt)), by = 'Region') %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS + AI + WGOA',
    Region == 2 ~ 'CGOA',
    Region == 3 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS + AI + WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = three_sdrep$sd[names(three_sdrep$value) == 'log_SSB'],
         Year = Year + 1959)

# bind five region historical estimates
five_rg_ssb <- five_rg_ssb %>% mutate(Sim = 1, Type = 'Historical') %>%
  bind_rows(reshape2::melt(five_ssb_proj) %>%
              rename(Region = Var1, Year = Var2, Sim = Var3) %>%
              dplyr::mutate(Region = dplyr::case_when(
                Region == 1 ~ 'BS',
                Region == 2 ~ 'AI',
                Region == 3 ~ 'WGOA',
                Region == 4 ~ 'CGOA',
                Region == 5 ~ 'EGOA'
              ),
              Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
              log_se = NA,
              Year = Year + 2020,
              Type = 'Projection'))

# bind three region historical estimates
three_rg_ssb <- three_rg_ssb %>% mutate(Sim = 1, Type = 'Historical') %>%
  bind_rows(reshape2::melt(three_ssb_proj) %>%
              rename(Region = Var1, Year = Var2, Sim = Var3) %>%
              dplyr::mutate(Region = dplyr::case_when(
                Region == 1 ~ 'BS + AI + WGOA',
                Region == 2 ~ 'CGOA',
                Region == 3 ~ 'EGOA'
              ),
              Region = factor(Region, levels = c("BS + AI + WGOA", "CGOA", "EGOA")),
              log_se = NA,
              Year = Year + 2020,
              Type = 'Projection'))

# bind aggregated historical estimates
all_rg_ssb <- all_rg_ssb %>% mutate(Sim = 1, Type_Period = 'Historical') %>%
  bind_rows(

    reshape2::melt(sgl_ssb_proj) %>%
      rename(Region = Var1, Year = Var2, Sim = Var3) %>%
      dplyr::mutate(,
                    log_se = NA,
                    Year = Year + 2020,
                    Type_Period = 'Projection',
                    Type = 'Single-Region') %>%
      group_by(Year, Sim, Type_Period, Type) %>%
      summarize(value = sum(value)),

    reshape2::melt(three_ssb_proj) %>%
      rename(Region = Var1, Year = Var2, Sim = Var3) %>%
      dplyr::mutate(,
                    log_se = NA,
                    Year = Year + 2020,
                    Type_Period = 'Projection',
                    Type = 'Three-Region') %>%
      group_by(Year, Sim, Type_Period, Type) %>%
      summarize(value = sum(value)),

    reshape2::melt(five_ssb_proj) %>%
      rename(Region = Var1, Year = Var2, Sim = Var3) %>%
      dplyr::mutate(,
                    log_se = NA,
                    Year = Year + 2020,
                    Type_Period = 'Projection',
                    Type = 'Five-Region') %>%
      group_by(Year, Sim, Type_Period, Type) %>%
      summarize(value = sum(value))

  ) %>%
  left_join(data.frame(Type = c("Single-Region", "Three-Region", "Five-Region"),
                       b40 = c(sgl_ref_pt$b_ref_pt, sum(three_ref_pt$b_ref_pt), sum(five_ref_pt$b_ref_pt)))) %>%
  mutate(Type = factor(Type, levels = c("Single-Region", "Three-Region", "Five-Region")))

# Three region Rec data
three_rg_rec <- reshape2::melt(three_rep$Rec) %>%
  rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS + AI + WGOA',
    Region == 2 ~ 'CGOA',
    Region == 3 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS + AI + WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = three_sdrep$sd[names(three_sdrep$value) == 'log_Rec'],
         Year = Year + 1959)

# Five region rec data
five_rg_rec <- reshape2::melt(five_rep$Rec) %>%
  rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS',
    Region == 2 ~ 'AI',
    Region == 3 ~ 'WGOA',
    Region == 4 ~ 'CGOA',
    Region == 5 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = five_sdrep$sd[names(five_sdrep$value) == 'log_Rec'],
         Year = Year + 1959)

# Plots -------------------------------------------------------------------

## Regional Biomass --------------------------------------------------------
# five region plot
five_ssb_plot <- ggplot() +
  geom_line(five_rg_ssb %>% filter(Type == 'Historical'),
            mapping = aes(x = Year, y = value), color = '#00a473', lty = 1, lwd = 1.3) +
  geom_ribbon(five_rg_ssb %>% filter(Type == 'Historical'),
              mapping = aes(x = Year, y = value, ymin = exp(log(value) - 1.96 * log_se),
                            ymax = exp(log(value) + 1.96 * log_se)), alpha = 0.25, color = NA, fill = '#00a473') +
  geom_line(five_rg_ssb %>% filter(Type == 'Projection'),
            mapping = aes(x = Year, y = value, group = Sim), color = '#00a473', lty = 1, lwd = 1.3, alpha = 1) +
  geom_hline(five_rg_ssb, mapping = aes(yintercept = b40), lty = 2, lwd = 1, color = '#00a473') +
  geom_vline(xintercept = 2021, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~Region, nrow = 5) +
  labs(x = 'Year', y = '') +
  theme_bw(base_size = 25)

# three region plot
three_ssb_plot <- ggplot() +
  geom_line(three_rg_ssb %>% filter(Type == 'Historical'),
            mapping = aes(x = Year, y = value), color = '#ef5f10', lty = 1, lwd = 1.3) +
  geom_ribbon(three_rg_ssb %>% filter(Type == 'Historical'),
              mapping = aes(x = Year, y = value, ymin = exp(log(value) - 1.96 * log_se),
                            ymax = exp(log(value) + 1.96 * log_se)), alpha = 0.25, color = NA, fill = '#ef5f10') +
  geom_line(three_rg_ssb %>% filter(Type == 'Projection'),
            mapping = aes(x = Year, y = value, group = Sim), color = '#ef5f10', lty = 1, lwd = 1.3, alpha = 1) +
  geom_hline(three_rg_ssb, mapping = aes(yintercept = b40), lty = 2, lwd = 1, color = '#ef5f10') +
  geom_vline(xintercept = 2021, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~Region, nrow = 3) +
  labs(x = 'Year', y = 'Spawning Stock Biomass ') +
  theme_bw(base_size = 25)

# aggregated plot
all_ssb_plot <- ggplot() +
  geom_line(all_rg_ssb %>% filter(Type_Period == 'Historical'),
            mapping = aes(x = Year, y = value, color = Type), lty = 1, lwd = 1.3) +
  geom_ribbon(all_rg_ssb %>% filter(Type_Period == 'Historical'),
              mapping = aes(x = Year, y = value, ymin = exp(log(value) - 1.96 * log_se),
                            ymax = exp(log(value) + 1.96 * log_se), fill = Type), alpha = 0.25, color = NA) +
  geom_line(all_rg_ssb %>% filter(Type_Period == 'Projection'),
            mapping = aes(x = Year, y = value, group = interaction(Sim, Type), color = Type), lty = 1, lwd = 1.3, alpha = 1) +
  geom_hline(all_rg_ssb, mapping = aes(yintercept = b40, color = Type), lty = 2, lwd = 1) +
  geom_vline(xintercept = 2021, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = c("#7e79b0", "#ef5f10", "#00a473")) +
  scale_fill_manual(values = c("#7e79b0", "#ef5f10", "#00a473")) +
  labs(x = 'Year', y = 'Spawning Stock Biomass ', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 25) +
  theme(legend.position = c(0.895, 0.125),
        legend.background = element_blank())


# combine plots
multi_ssb_plots <- cowplot::plot_grid(three_ssb_plot, five_ssb_plot, ncol = 2,
                                      labels = c("B", "C"), label_size = 30, label_x = 0.03)

ggsave(
  here("dev", "paper_projects", "sporc_manuscript_demonstrations", "figs", "spatial_proj.png"),
  print(
    cowplot::plot_grid(all_ssb_plot, multi_ssb_plots, rel_heights = c(0.6, 1), nrow = 2,
                       labels = c("A", ""), label_size = 30)
  ),
  width = 14, height = 20
)

ggsave(
  here("dev", "paper_projects", "sporc_manuscript_demonstrations", "figs", "spatial_proj_presentation.png"),
  all_ssb_plot + theme(legend.position = c(0.15, 0.125)),
  width = 10, height = 8
)

ggsave(
  here("dev", "paper_projects", "sporc_manuscript_demonstrations", "figs", "three_area_spatial_proj_presentation.png"),
  three_ssb_plot,
  width = 5, height = 8
)

## Movement ----------------------------------------------------------------
n_t <- 60 # 60 forward projections
NAA_three <- array(0, dim = c(three_dat$n_regions, n_t)) # numbers at age container
NAA_three[,1] <- three_rep$rec_region_prop # input recruitment proportions for first time step
NAA_five <- array(0, dim = c(five_dat$n_regions, n_t)) # numbers at age container
NAA_five[,1] <- five_rep$rec_region_prop # input recruitment proportions for first time step

for(t in 2:n_t) {

  # three region stationary movement
  if(t <= dim(three_rep$Movement)[6]) NAA_three[,t] <- NAA_three[,t-1] %*% three_rep$Movement[1,,,1,1,t,1]
  else NAA_three[,t] <- NAA_three[,t-1] %*% three_rep$Movement[1,,,1,1,30,1]

  # five region stationary movement
  if(t <= dim(five_rep$Movement)[6]) NAA_five[,t] <- NAA_five[,t-1] %*% five_rep$Movement[1,,,1,1,t,1]
  else NAA_five[,t] <- NAA_five[,t-1] %*% five_rep$Movement[1,,,1,1,30,1]

}

# format arrays
NAA_three_df <- reshape2::melt(NAA_three) %>%
  rename(Time = Var2, Region = Var1) %>%
  mutate(
    Region = case_when(
      Region == 1 ~ "BS + AI + WGOA",
      Region == 2 ~ "CGOA",
      Region == 3 ~ "EGOA"
    )
  )


# format arrays
NAA_five_df <- reshape2::melt(NAA_five) %>%
  rename(Time = Var2, Region = Var1) %>%
  mutate(
    Region = case_when(
      Region == 1 ~ "BS",
      Region == 2 ~ "AI",
      Region == 3 ~ "WGOA",
      Region == 4 ~ "CGOA",
      Region == 5 ~ "EGOA"
    ),
    Region = factor(Region,
                    levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))
  )

# plot!
cohort_three_area_plot <- ggplot(NAA_three_df, aes(x = Time, y = value)) +
  geom_line(color = "#ef5f10", lwd = 1.1) +
  facet_wrap( ~ Region, ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = 'Time', y = "Proportion of Cohort") +
  theme_bw(base_size = 18)

cohort_five_area_plot <- ggplot(NAA_five_df, aes(x = Time, y = value)) +
  geom_line(color = "#00a473", lwd = 1.1) +
  facet_wrap( ~ Region, ncol = 1) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = 'Time', y = "Proportion of Cohort") +
  theme_bw(base_size = 18)

rec_three_area_plot <- ggplot(three_rg_rec, aes(x = Year, y = value,
                                                ymin = exp(log(value) - 1.96 * log_se),
                                                ymax = exp(log(value) + 1.96 * log_se) )) +
  geom_line(color = "#ef5f10") +
  geom_ribbon(alpha = 0.25, fill = "#ef5f10") +
  facet_wrap( ~ Region, ncol = 1) +
  labs(x = 'Year', y = 'Age-2 Recruitment') +
  theme_bw(base_size = 18)

rec_five_area_plot <- ggplot(five_rg_rec, aes(x = Year, y = value,
                                                ymin = exp(log(value) - 1.96 * log_se),
                                                ymax = exp(log(value) + 1.96 * log_se) )) +
  geom_line(color = "#00a473") +
  coord_cartesian(ylim = c(0,100))+
  geom_ribbon(alpha = 0.25, fill = "#00a473") +
  facet_wrap( ~ Region, ncol = 1) +
  labs(x = 'Year', y = 'Age-2 Recruitment') +
  theme_bw(base_size = 18)

ggsave(
  here("dev", "paper_projects", "sporc_manuscript_demonstrations", "figs", "move_proj.png"),
  print(
    cowplot::plot_grid(rec_three_area_plot, cohort_three_area_plot, rec_five_area_plot,
                       cohort_five_area_plot,
                       ncol = 2, labels = c("A", "B", "C", "D"), label_size = 25, label_x = 0.)
  ),
  width = 10, height = 16
)

ggsave(
  here("dev", "paper_projects", "sporc_manuscript_demonstrations", "figs", "rec_presentation.png"),
  rec_three_area_plot + theme(legend.position = c(0.15, 0.125)),
  width = 5, height = 5
)


## HCR Comparison ----------------------------------------------------------

# One Region HCR
one_hcr_df <- expand.grid(j = 1, i = seq(0, 300, 0.01)) %>%
  mutate(
    frp = mapply(function(j) sgl_ref_pt$f_ref_pt[j], j),
    brp = mapply(function(j) sgl_ref_pt$b_ref_pt[j], j),
    F = mapply(function(i, j) {
      HCR_function(x = i, frp = sgl_ref_pt$f_ref_pt[j], brp = sgl_ref_pt$b_ref_pt[j])
    }, i, j),
    F_F40 = F / frp,
    SSB_B40 = i / brp
  ) %>%
  mutate(j = 'BS + AI + WGOA + CGOA + EGOA')

# Three Region HCR
three_hcr_df <- expand.grid(j = 1:3, i = seq(0, 100, 0.01)) %>%
  mutate(
    frp = mapply(function(j) three_ref_pt$f_ref_pt[j], j),
    brp = mapply(function(j) three_ref_pt$b_ref_pt[j], j),
    F = mapply(function(i, j) {
      HCR_function(x = i, frp = three_ref_pt$f_ref_pt[j], brp = three_ref_pt$b_ref_pt[j])
    }, i, j),
    F_F40 = F / frp,
    SSB_B40 = i / brp
  ) %>%
  dplyr::mutate(j = dplyr::case_when(
    j == 1 ~ 'BS + AI + WGOA',
    j == 2 ~ 'CGOA',
    j == 3 ~ 'EGOA'
  ),
  j = factor(j, levels = c("BS + AI + WGOA", "CGOA", "EGOA")))

# Five Region HCR
five_hcr_df <- expand.grid(j = 1:5, i = seq(0, 50, 0.01)) %>%
  mutate(
    frp = mapply(function(j) five_ref_pt$f_ref_pt[j], j),
    brp = mapply(function(j) five_ref_pt$b_ref_pt[j], j),
    F = mapply(function(i, j) {
      HCR_function(x = i, frp = five_ref_pt$f_ref_pt[j], brp = five_ref_pt$b_ref_pt[j])
    }, i, j),
    F_F40 = F / frp,
    SSB_B40 = i / brp
  ) %>%
  dplyr::mutate(j = dplyr::case_when(
    j == 1 ~ 'BS',
    j == 2 ~ 'AI',
    j == 3 ~ 'WGOA',
    j == 4 ~ 'CGOA',
    j == 5 ~ 'EGOA'
  ),
  j = factor(j, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")))

one_hcr_plot <- ggplot(one_hcr_df, aes(x = i, y = F, color = factor(j))) +
  geom_textpath(aes(label = j), color = '#7e79b0', lwd = 1.3, size = 3.5, hjust = 0.3, vjust = 0.5, textcolour = 'black') +
  facet_wrap(~"Single-Region", nrow = 1) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none') +
  labs(x = '', y = "")

three_hcr_plot <- ggplot(three_hcr_df, aes(x = i, y = F, color = factor(j), lty = factor(j))) +
  geom_textpath(aes(label = j), color = '#ef5f10', lwd = 1.3, size = 3.5, hjust = 0.25, vjust = 0.5, textcolour = 'black') +
  facet_wrap(~"Three-Region", nrow = 1) +
  labs(x = "SSB", y = "") +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none') +
  labs(x = '', y = bquote(Instantaneous~Fishing~Mortality~(y^-1)))

five_hcr_plot <- ggplot(five_hcr_df, aes(x = i, y = F, group = factor(j), lty = factor(j)) ) +
  geom_textpath(
    aes(label = j), color = '#00a473', lwd = 1.3, size = 3.5, hjust = 0.4, vjust = 0.6, textcolour = 'black') +
  facet_wrap(~"Five-Region", nrow = 1) +
  labs( x = "Spawning Stock Biomass", y = '') +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

# combine plots
multi_hcr_plots <- cowplot::plot_grid(one_hcr_plot, three_hcr_plot, five_hcr_plot, ncol = 1, align = 'hv')

ggsave(
  here("dev", "paper_projects", "sporc_manuscript_demonstrations", "figs", "spatial_hcr.png"),
  multi_hcr_plots,
  height = 13, width = 8
)

ggsave(
  here("dev", "paper_projects", "sporc_manuscript_demonstrations", "figs", "spatial_hcr_presentation.png"),
  three_hcr_plot +facet_null(),
  height = 5, width = 5
)

