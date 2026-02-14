# Purpose: Example script for deriving reference points and projections
# Creator: Matthew LH. Cheng
# Date 6/29/25


# Setup -------------------------------------------------------------------

library(SPoRC)
library(here)
library(tidyverse)

# load in datasets
data("sgl_rg_sable_rep") # read in single region report
data("sgl_rg_sable_data") # read in single region data
data("mlt_rg_sable_rep") # read in multi region report
data("mlt_rg_sable_data") # read in multi region data

# Single Region SPR -------------------------------------------------------

# single area model
sgl_ref_pt <- SPoRC::Get_Reference_Points(data = sgl_rg_sable_data, # data file
                                          rep = sgl_rg_sable_rep, # report file
                                          SPR_x = 0.4, # spr rate
                                          type = 'single_region', # single region reference point
                                          what = 'SPR', # SPR reference point
                                          calc_rec_st_yr = 20, # first year to calculate mean recruitment
                                          rec_age = 2,  # exclues the last rec_age years when computing mean recruitment
                                          t_spawn = 0,  # spawn timing
                                          sex_ratio_f = 0.5 # recruitment sex-ratio for females
)
sgl_ref_pt$f_ref_pt # F40
sgl_ref_pt$b_ref_pt # B40


# Multi Region SPR --------------------------------------------------------
## Independent -------------------------------------------------------------
# multi region model with independent SPR
mlt_ref_pt_indp <- Get_Reference_Points(data = mlt_rg_sable_data, # data file
                                               rep = mlt_rg_sable_rep, # report file
                                               SPR_x = 0.4, # spr rate
                                               type = 'multi_region', # multi region reference point
                                               what = 'independent_SPR', # SPR reference point
                                               calc_rec_st_yr = 20, # first year to calculate mean recruitment
                                               rec_age = 2,  # exclues the last rec_age years when computing mean recruitment
                                               t_spawn = 0,  # spawn timing
                                               sex_ratio_f = rep(0.5, mlt_rg_sable_data$n_regions) # recruitment sex-ratio for females
)
mlt_ref_pt_indp$f_ref_pt # F40
mlt_ref_pt_indp$b_ref_pt # B40

## Global ------------------------------------------------------------------
# multi region model with global SPR
mlt_ref_pt_global <- Get_Reference_Points(data = mlt_rg_sable_data, # data file
                                          rep = mlt_rg_sable_rep, # report file
                                          SPR_x = 0.4, # spr rate
                                          type = 'multi_region', # multi region reference point
                                          what = 'global_SPR', # SPR reference point
                                          calc_rec_st_yr = 20, # first year to calculate mean recruitment
                                          rec_age = 2,  # exclues the last rec_age years when computing mean recruitment
                                          t_spawn = 0,  # spawn timing
                                          sex_ratio_f = rep(0.5, mlt_rg_sable_data$n_regions) # recruitment sex-ratio for females
)
mlt_ref_pt_global$f_ref_pt # F40
mlt_ref_pt_global$b_ref_pt # B40


# Catch Projections -------------------------------------------------------
## Single Region -----------------------------------------------------------

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

# Create a tibble for plotting
hcr_df <- tibble(
  i = 1:200,
  SSB_B40 = i / sgl_ref_pt$b_ref_pt,
  F = sapply(i, function(x) {
    HCR_function(x = x, frp = sgl_ref_pt$f_ref_pt, brp = sgl_ref_pt$b_ref_pt)
  })
)

# Plot
png(here("vignettes", "figures", "i_threshold_hcr.png"), width = 1000, height = 800)
ggplot(hcr_df, aes(x = SSB_B40, y = F)) +
  geom_line(color = "steelblue", size = 1) +
  labs(x = "SSB / B40", y = "F") +
  theme_bw(base_size = 13)
dev.off()

# Setup necessary inputs
t_spawn <- 0 # spawn timing
n_proj_yrs <- 15 # number of projection years
n_regions <- 1 # number of regions
n_ages <- length(sgl_rg_sable_data$ages) # number of ages
n_sexes <- sgl_rg_sable_data$n_sexes # number of sexes
n_fish_fleets <- sgl_rg_sable_data$n_fish_fleets # number of fishery fleets
n_seas <- 1
do_recruits_move <- 0 # recruits don't move
terminal_NAA <- array(sgl_rg_sable_rep$NAA[,length(sgl_rg_sable_data$years),,,], dim = c(n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
terminal_NAA0 <- array(sgl_rg_sable_rep$NAA0[,length(sgl_rg_sable_data$years),,,], dim = c(n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
WAA <- array(rep(sgl_rg_sable_data$WAA[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(sgl_rg_sable_data$WAA[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # weight at age fishery
MatAA <- array(rep(sgl_rg_sable_data$MatAA[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(sgl_rg_sable_rep$fish_sel[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- array(rep(sgl_rg_sable_rep$Movement[,,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_regions, n_seas, n_proj_yrs, n_ages, n_sexes)) # movement
terminal_F <- array(sgl_rg_sable_rep$Fmort[,length(sgl_rg_sable_data$years),,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F
natmort <- array(sgl_rg_sable_rep$natmort[,length(sgl_rg_sable_data$years),,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # natural mortality
recruitment <- array(sgl_rg_sable_rep$Rec[,20:(length(sgl_rg_sable_data$years) - 2)], dim = c(n_regions, length(20:length(sgl_rg_sable_data$years) - 2))) # recruitment values to use for mean recruitment calculations or inverse gaussian parameterization
sexratio <- array(0.5, dim = c(n_regions, n_proj_yrs, n_sexes)) # recruitment sex ratio

# Define reference points to use in HCR
f_ref_pt = array(sgl_ref_pt$f_ref_pt, dim = c(n_regions, n_proj_yrs))
b_ref_pt = array(sgl_ref_pt$b_ref_pt, dim = c(n_regions, n_proj_yrs))

# do population projection
out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs, # Number of projection years
                                n_regions = n_regions, # number of regions
                                n_ages = n_ages, # number of ages
                                n_sexes = n_sexes, # number of sexes
                                sexratio = sexratio, # sex ratio for recruitment
                                n_fish_fleets = n_fish_fleets, # number of fishery fleets
                                do_recruits_move = do_recruits_move, # whether recruits move (not used since single area)
                                recruitment = recruitment, # recruitment values to use for mean recruitment
                                terminal_NAA = terminal_NAA, # terminal numbers at age
                                terminal_NAA0 = terminal_NAA0,
                                terminal_F = terminal_F, # terminal F
                                natmort = natmort, # natural mortality values to use in projection
                                WAA = WAA, # weight at age values to use in projection
                                WAA_fish = WAA_fish, # weight at age values to use for fishery catch
                                MatAA = MatAA, # maturity at age values to use in projection
                                fish_sel = fish_sel, # fishery selectivity values to use in projection
                                Movement = Movement, # movement values (not used since single area)
                                f_ref_pt = f_ref_pt, # fishery reference points (f40)
                                b_ref_pt = b_ref_pt, # biological reference points (b40)
                                HCR_function = HCR_function, # threshold control rule defined above
                                recruitment_opt = "mean_rec", # recruitment assumption utilizes the mean recruits for the supplied recruitment estimates
                                fmort_opt = "HCR", # Fishing mortality in projection years are determined using a HCR
                                t_spawn = t_spawn # Spawn timing
)

# ssb
combined_ssb <- c(sgl_rg_sable_rep$SSB[1, -65], out$proj_SSB[1,]) # removing terminal year becauase repeated in projection calculations
years <- 1960:(2023 + n_proj_yrs)

ssb_df <- tibble(
  Year = years,
  SSB = combined_ssb
)

# Plot
png(here("vignettes", "figures", "i_sgl_proj_ssb.png"), width = 1000, height = 800)
ggplot(ssb_df, aes(x = Year, y = SSB)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed") + # projection start
  scale_y_continuous(limits = c(0, 300)) +
  labs(x = "Year", y = "SSB (kt)") +
  theme_bw(base_size = 13)
dev.off()

# catch
combined_catch <- c(
  rowSums(sgl_rg_sable_rep$PredCatch[1, -65, 1, ]), # removing terminal year becauase repeated in projection calculations
  rowSums(out$proj_Catch[1, , 1, ])
)

years <- 1960:(2023 + n_proj_yrs)

catch_df <- tibble(
  Year = years,
  Catch = combined_catch
)

# Plot
png(here("vignettes", "figures", "i_sgl_proj_catch.png"), width = 1000, height = 800)
ggplot(catch_df, aes(x = Year, y = Catch)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed") +  # projection start
  labs(x = "Year", y = "Catch (kt)") +
  theme_bw(base_size = 13)
dev.off()
sum(out$proj_Catch[1,2,1,]) # Catch advice in terminal year + 1


## Multi Region ------------------------------------------------------------

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

# Build a dataframe by looping over j and i scalars
hcr_df <- expand.grid(
  j = 1:5,
  i = 1:50
) %>%
  mutate(
    frp = mapply(function(j) mlt_ref_pt_indp$f_ref_pt[j], j),
    brp = mapply(function(j) mlt_ref_pt_indp$b_ref_pt[j], j),
    F = mapply(function(i, j) {
      HCR_function(x = i, frp = mlt_ref_pt_indp$f_ref_pt[j], brp = mlt_ref_pt_indp$b_ref_pt[j])
    }, i, j),
    SSB_B40 = i / brp
  )

png(here("vignettes", "figures", "i_mlt_threshold_hcr.png"), width = 1000, height = 800)
ggplot(hcr_df, aes(x = SSB_B40, y = F, color = factor(j))) +
  geom_line(lwd = 1.3) +
  facet_wrap(~j, scales = 'free') +
  labs(x = "SSB / B40", y = "F",  color = 'Region') +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none')
dev.off()

# Setup necessary inputs
t_spawn <- 0 # spawn timing
n_proj_yrs <- 15 # number of projection years
n_regions <- 5 # number of regions
n_ages <- length(mlt_rg_sable_data$ages) # number of ages
n_sexes <- mlt_rg_sable_data$n_sexes # number of sexes
n_fish_fleets <- mlt_rg_sable_data$n_fish_fleets # number of fishery fleets
do_recruits_move <- 0 # recruits don't move
n_seas <- 1
terminal_NAA <- array(mlt_rg_sable_rep$NAA[,length(mlt_rg_sable_data$years),,,], dim = c(n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
terminal_NAA0 <- array(mlt_rg_sable_rep$NAA0[,length(mlt_rg_sable_data$years),,,], dim = c(n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
WAA <- array(rep(mlt_rg_sable_data$WAA[,length(mlt_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(mlt_rg_sable_data$WAA[,length(mlt_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # weight at age for fishery
MatAA <- array(rep(mlt_rg_sable_data$MatAA[,length(mlt_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(mlt_rg_sable_rep$fish_sel[,length(mlt_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- abind::abind(replicate(n_proj_yrs,  mlt_rg_sable_rep$Movement[,,length(mlt_rg_sable_data$years),,,,drop = FALSE],  simplify = FALSE),along = 3) # movement
terminal_F <- array(mlt_rg_sable_rep$Fmort[,length(mlt_rg_sable_data$years),,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F
natmort <- array(mlt_rg_sable_rep$natmort[,length(mlt_rg_sable_data$years),,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # natural mortality
recruitment <- array(mlt_rg_sable_rep$Rec[,20:(length(mlt_rg_sable_data$years) - 2)], dim = c(n_regions, length(20:length(mlt_rg_sable_data$years) - 2))) # recruitment values to use for mean recruitment calculations or inverse gaussian parameterization
sexratio <- array(0.5, dim = c(n_regions, n_proj_yrs, n_sexes)) # recruitment sex ratio


# Define independent SPR reference points to use in HCR
f_ref_pt_indp = array(mlt_ref_pt_indp$f_ref_pt, dim = c(n_regions, n_proj_yrs))
b_ref_pt_indp = array(mlt_ref_pt_indp$b_ref_pt, dim = c(n_regions, n_proj_yrs))

# do population projection
out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs, # Number of projection years
                                n_regions = n_regions, # number of regions
                                n_ages = n_ages, # number of ages
                                n_sexes = n_sexes, # number of sexes
                                sexratio = sexratio, # sex ratio for recruitment
                                n_fish_fleets = n_fish_fleets, # number of fishery fleets
                                do_recruits_move = do_recruits_move, # whether recruits move (not used since single area)
                                recruitment = recruitment, # recruitment values to use for mean recruitment
                                terminal_NAA = terminal_NAA, # terminal numbers at age
                                terminal_NAA0 = terminal_NAA0,
                                terminal_F = terminal_F, # terminal F
                                natmort = natmort, # natural mortality values to use in projection
                                WAA = WAA, # weight at age values to use in projection
                                WAA_fish = WAA_fish, # weight at age values for fishery catch
                                MatAA = MatAA, # maturity at age values to use in projection
                                fish_sel = fish_sel, # fishery selectivity values to use in projection
                                Movement = Movement, # movement values (not used since single area)
                                f_ref_pt = f_ref_pt_indp, # fishery reference points (f40)
                                b_ref_pt = b_ref_pt_indp, # biological reference points (b40)
                                HCR_function = HCR_function, # threshold control rule defined above
                                recruitment_opt = "mean_rec", # recruitment assumption utilizes the mean recruits for the supplied recruitment estimates
                                fmort_opt = "HCR", # Fishing mortality in projection years are determined using a HCR
                                t_spawn = t_spawn # Spawn timing
)

combined_ssb <- cbind(mlt_rg_sable_rep$SSB[,-62], out$proj_SSB[,]) # removing terminal year becauase repeated in projection calculations
combined_ssb_df <- reshape2::melt(combined_ssb) %>%
  rename(Region = Var1, Year = Var2, SSB = value)

# Plot
png(here("vignettes", "figures", "i_mlt_proj_ssb.png"), width = 1000, height = 800)
ggplot(combined_ssb_df, aes(x = Year + 1959, y = SSB, color = factor(Region))) +
  geom_line(size = 1) +
  geom_vline(xintercept = 2021, linetype = "dashed") + # projection start
  facet_wrap(~Region) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "SSB (kt)") +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none')
dev.off()

combined_catch <- cbind(apply(mlt_rg_sable_rep$PredCatch, c(1,2), sum), apply(out$proj_Catch, c(1,2), sum))
combined_catch_df <- reshape2::melt(combined_catch) %>%
  rename(Region = Var1, Year = Var2, Catch = value)

# Plot
png(here("vignettes", "figures", "i_mlt_proj_catch.png"), width = 1000, height = 800)
ggplot(combined_catch_df, aes(x = Year + 1959, y = Catch, color = factor(Region))) +
  geom_line(size = 1) +
  geom_vline(xintercept = 2021, linetype = "dashed") + # projection start
  facet_wrap(~Region) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "Catch") +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none')
dev.off()

rowSums(out$proj_Catch[,2,1,]) # Catch advice by region in terminal year + 1


# Stochastic Projections --------------------------------------------------

# Setup necessary inputs
n_sims <- 1e3 # number of simulations to conduct
t_spawn <- 0 # spawn timing
n_proj_yrs <- 15 # number of projection years
n_regions <- 1 # number of regions
n_ages <- length(sgl_rg_sable_data$ages) # number of ages
n_sexes <- sgl_rg_sable_data$n_sexes # number of sexes
n_fish_fleets <- sgl_rg_sable_data$n_fish_fleets # number of fishery fleets
do_recruits_move <- 0 # recruits don't move
n_seas <- 1
terminal_NAA <- array(sgl_rg_sable_rep$NAA[,length(sgl_rg_sable_data$years),,,], dim = c(n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
terminal_NAA0 <- array(sgl_rg_sable_rep$NAA0[,length(sgl_rg_sable_data$years),,,], dim = c(n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
WAA <- array(rep(sgl_rg_sable_data$WAA[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(sgl_rg_sable_data$WAA[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # weight at age fpr fishery catch
MatAA <- array(rep(sgl_rg_sable_data$MatAA[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(sgl_rg_sable_rep$fish_sel[,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- array(rep(sgl_rg_sable_rep$Movement[,,length(sgl_rg_sable_data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # movement
terminal_F <- array(sgl_rg_sable_rep$Fmort[,length(sgl_rg_sable_data$years),,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F
natmort <- array(sgl_rg_sable_rep$natmort[,length(sgl_rg_sable_data$years),,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # natural mortality
recruitment <- array(sgl_rg_sable_rep$Rec[,20:(length(sgl_rg_sable_data$years) - 2)], dim = c(n_regions, length(20:length(sgl_rg_sable_data$years) - 2))) # recruitment values to use for mean recruitment calculations or inverse gaussian parameterization
sexratio <- array(0.5, dim = c(n_regions, n_proj_yrs, n_sexes)) # recruitment sex ratio

# Define the F used for each scenario
proj_inputs <- list(
  # Scenario 1 - Using HCR to adjust f40
  list(f_ref_pt = array(sgl_ref_pt$f_ref_pt, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(sgl_ref_pt$b_ref_pt, dim = c(n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 2 - F is set at 0
  list(f_ref_pt = array(0, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = NULL,
       fmort_opt = 'Input'
  )
)

# store outputs
all_scenarios_f <- array(0, dim = c(n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
all_scenarios_ssb <- array(0, dim = c(n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
all_scenarios_catch <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims, length(proj_inputs)))

set.seed(123)
for (i in seq_along(proj_inputs)) {
  for (sim in 1:n_sims) {

    # do population projection
    out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs, # number of projection years
                                    n_regions = n_regions, # number of regions
                                    n_ages = n_ages, # number of ages
                                    n_sexes = n_sexes, # number of sexes
                                    sexratio = sexratio, # sex ratio
                                    n_fish_fleets = n_fish_fleets, # number of fleets
                                    do_recruits_move = do_recruits_move, # whether recruits move
                                    recruitment = recruitment, # recruitment values to use to parameterize inverse gaussian
                                    terminal_NAA = terminal_NAA, # terminal numbers at age
                                    terminal_NAA0 = terminal_NAA, # unfished terminal numbers at age
                                    terminal_F = terminal_F, # terminal fishing mortality at age
                                    natmort = natmort, # natural mortality
                                    WAA = WAA, # weight at age
                                    WAA_fish = WAA_fish, # weight at age for fishery catch
                                    MatAA = MatAA, # maturity at age
                                    fish_sel = fish_sel, # fishery selectivity
                                    Movement = Movement, # movement
                                    f_ref_pt = proj_inputs[[i]]$f_ref_pt, # fishing mortality to use for projection
                                    b_ref_pt = proj_inputs[[i]]$b_ref_pt, # biological reference point to use in HCR
                                    HCR_function = HCR_function, # harvest control rule function
                                    recruitment_opt = "inv_gauss", # stochastic simulation for recruitment
                                    fmort_opt = proj_inputs[[i]]$fmort_opt, # fishing mortality option (either HCR or a user input)
                                    t_spawn = t_spawn # spawn timing
    )

    # store results
    all_scenarios_ssb[,,sim,i] <- out$proj_SSB
    all_scenarios_catch[,,,sim,i] <- out$proj_Catch
    all_scenarios_f[,,sim,i] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used

  } # end sim loop
  print(i)
} # end i loop

# Get historical SSB
historical <- reshape2::melt(array(rep(sgl_rg_sable_rep$SSB, n_sims),
                                   dim = c(n_regions, length(sgl_rg_sable_data$years), n_sims))) %>%
  mutate(Year = Var2 + 1959,
         Scenario = "FABC (F40)",  # or change to match the scenarios you're plotting
         Type = "Historical") %>%
  rename(Region = Var1, Simulation = Var3, SSB = value)

# Get all scenario projections
scenarios <- reshape2::melt(all_scenarios_ssb) %>%
  mutate(Year = Var2 + 2023,
         Scenario = case_when(
           Var4 == 1 ~ "S1: FABC (F40)",
           Var4 == 2 ~ "S2: F = 0"
         ),
         Type = "Projection") %>%
  rename(Region = Var1, Simulation = Var3, SSB = value)

# expand historical SSB for plotting
scenarios_unique <- unique(scenarios$Scenario)
historical_expanded <- historical[rep(1:nrow(historical), times = length(scenarios_unique)), ]
historical_expanded$Scenario <- rep(scenarios_unique, each = nrow(historical))

# combine
combined_ssb <- bind_rows(historical_expanded, scenarios)

# Plot
png(here("vignettes", "figures", "i_sgl_stoch_proj_ssb.png"), width = 1000, height = 800)
combined_ssb %>%
  ggplot(aes(x = Year, y = SSB, group = interaction(Scenario, Simulation), color = Type)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~Scenario, scales = 'free') +
  geom_hline(yintercept = sgl_ref_pt$b_ref_pt, lty = 2) + # b40
  geom_vline(xintercept = 2024, lty = 2) + # projection start
  scale_color_manual(values = c("Historical" = "black", "Projection" = "blue")) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none')
dev.off()

# Get historical catch
historical <- reshape2::melt(array(rep(sgl_rg_sable_data$ObsCatch, n_sims),
                                   dim = c(n_regions, length(sgl_rg_sable_data$years), sgl_rg_sable_data$n_fish_fleets, n_sims))) %>%
  mutate(Year = Var2 + 1959,
         Scenario = "FABC (F40)",  # or change to match the scenarios you're plotting
         Type = "Historical") %>%
  rename(Region = Var1, Simulation = Var4, Fleet = Var3, Catch = value) %>%
  select(-Var2)

historical$Catch[is.na(historical$Catch)] <- 0

# Get all scenario projections
scenarios <- reshape2::melt(all_scenarios_catch) %>%
  mutate(Year = Var2 + 2023,
         Scenario = case_when(
           Var5 == 1 ~ "S1: FABC (F40)",
           Var5 == 2 ~ "S2: F = 0"
         ),
         Type = "Projection") %>%
  rename(Region = Var1, Simulation = Var4, Catch = value, Fleet = Var3) %>%
  select(-c(Var2, Var5))

# expand historical SSB for plotting
scenarios_unique <- unique(scenarios$Scenario)
historical_expanded <- historical[rep(1:nrow(historical), times = length(scenarios_unique)), ]
historical_expanded$Scenario <- rep(scenarios_unique, each = nrow(historical))

# combine
combined_cat <- bind_rows(historical_expanded, scenarios)

# Plot
png(here("vignettes", "figures", "i_sgl_stoch_proj_catch.png"), width = 1000, height = 800)
combined_cat %>%
  group_by(Year, Scenario, Simulation, Type, Region) %>%
  summarize(Catch = sum(Catch)) %>%
  ggplot(aes(x = Year, y = Catch, group = interaction(Scenario, Simulation), color = Type)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~Scenario) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = c("Historical" = "black", "Projection" = "blue")) +
  theme_bw() +
  theme(legend.position = 'none')
dev.off()
