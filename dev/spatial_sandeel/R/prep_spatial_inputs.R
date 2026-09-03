# Purpose: clean the 3 region North Sea sandeel catch and survey at age into SPoRC arrays
# Date Created: 9/3/26

# read stuff in
library(here)
raw <- utils::read.csv(here("dev", "spatial_sandeel", "inputs", "spatial_catches.csv"), stringsAsFactors = FALSE)
srv_raw <- utils::read.csv(here("dev", "spatial_sandeel", "inputs", "spatial_survey.csv"), stringsAsFactors = FALSE)
eff_raw <- utils::read.csv(here("dev", "spatial_sandeel", "inputs", "Standardized_effort_days_2025.csv"), stringsAsFactors = FALSE)

# fxn to help read in the .in files
read_in <- function(file, year_first) {
  x <- readLines(here("dev", "spatial_sandeel", "inputs", file), warn = FALSE)
  x <- x[!grepl("^\\s*#", x) & nzchar(trimws(x))]
  lab <- strsplit(trimws(sub(".*#\\s*", "", x)), "\\s+")
  yr <- as.integer(sapply(lab, `[`, if(year_first) 1 else 2))
  sn <- as.integer(sapply(lab, `[`, if(year_first) 2 else 1))
  val <- t(sapply(x, function(z) scan(text = sub("#.*", "", z), quiet = TRUE)))
  list(year = yr, seas = sn, val = unname(val))
} # end read_in

# get biologicals
bio <- list(M = read_in("natmor.in", FALSE),
            mat = read_in("propmat.in", TRUE),
            weca = read_in("weca.in", TRUE),
            west = read_in("west.in", TRUE))

# setup --------------------------------------------------------------------
end_year <- 2025 # end year
apply_nocatch <- TRUE   # drop the year and season cells the zero catch file names
fit_age0_catch <- FALSE # whether to map off selex at age-0 for catch
rec_seas <- 2  # recruitment season
effort_model <- "Mixed Annual Weekly"
effort_col <- "newE" # E times the area scale

# dimensions ------------------------------------------------------------------
regions <- c("EU1", "EU2", "UK")
ages <- 0:4
n_ages <- length(ages)
n_seas <- 2
seasdur <- c(0.5, 0.5)
n_pop <- 1
n_sexes <- 1
n_regions <- length(regions)
first_year <- 1983
yrs <- first_year:end_year
n_yrs <- length(yrs)

# one fishery fleet per region and season (3 regions x 2 season = 6 fleets)
fish_region <- rep(1:n_regions, each = n_seas)
fish_seas <- rep(1:n_seas, times = n_regions)
fish_names <- paste0(regions[fish_region], "_s", fish_seas)
n_fish <- length(fish_names)
fish_fleet <- matrix(1:n_fish, nrow = n_regions, ncol = n_seas, byrow = TRUE) # indexing for fleet used later

# survey setup stuff
srv_names <- c("EU1", "EU2", "UK")
srv_seas <- c(2, 2, 2)
srv_t <- c(0.75, 0.75, 0.75)
srv_type <- c("spltRspltS", "spltRspltS", "spltRspltS")
n_srv <- length(srv_names)

# catch at age ----------------------------------------------------------------
cols <- paste0("N", ages)
raw$r <- match(raw$area, regions)
raw$yi <- match(raw$year, yrs)

ObsCatchAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))
UseCatchAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))

# 0 = no catch that year and season, area wide
zc <- readLines(here("dev", "spatial_sandeel", "inputs", "zero_catch_year_season.in"), warn = FALSE)
zc <- zc[!grepl("^\\s*#", zc) & nzchar(trimws(zc))]
zc_val <- unname(t(sapply(zc, function(z) scan(text = sub("#.*", "", z), quiet = TRUE))))

# read by position, row 2 is mislabelled 1974 where 1984 belongs
nocatch <- zc_val[seq_len(n_yrs), , drop = FALSE]

# age 0 only exists from the recruitment season
fit_age <- function(a, s) {
  if(a >= 2) return(TRUE)
  fit_age0_catch && s >= rec_seas
} # end fit_age

# a fleet fishes one region and one season, so its other cells stay empty
for(i in seq_len(nrow(raw))) {

  y <- raw$yi[i]
  if(is.na(y)) next # outside model window

  r <- raw$r[i]
  s <- raw$season[i]
  f <- fish_fleet[r, s]

  for(a in seq_len(n_ages)) {

    obs <- raw[i, cols[a]]
    ObsCatchAA[r, y, s, a, 1, f] <- obs

    # fit if positive, age is available, and not masked out
    if(obs > 0 && fit_age(a, s) && (!apply_nocatch || nocatch[y, s] == 1)) {
      UseCatchAA[r, y, s, a, 1, f] <- 1
    }

  } # end a loop
} # end i loop

# aggregated catch off, fleets fit the at age form
ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))
UseCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))

# survey index at age ---------------------------------------------------------
ObsSrvIdxAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))
UseSrvIdxAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))

srv_raw$r <- match(srv_raw$area, regions)
srv_raw$yi <- match(srv_raw$year, yrs)
srv_raw$ai <- match(srv_raw$age, ages)

# fleet k is area k and sees region k only
for(i in seq_len(nrow(srv_raw))) {
  y <- srv_raw$yi[i]
  if(is.na(y)) next  # outside the model window
  if(srv_raw$index[i] <= 0) next
  r <- srv_raw$r[i]
  ObsSrvIdxAA[r, y, srv_seas[r], srv_raw$ai[i], 1, r] <- srv_raw$index[i]
  UseSrvIdxAA[r, y, srv_seas[r], srv_raw$ai[i], 1, r] <- 1

} # end i loop

ObsSrvIdx <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_srv))
UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))

# biologicals -----------------------------------------------------------------
WAA      <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
MatAA    <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
WAA_fish <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))
WAA_srv  <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))

# func to help get a row of a .in file for one year and season
bio_row <- function(b, year, seas) b$val[which(b$year == year & b$seas == seas)[1], ]

for(y in 1:n_yrs) {
  for(seas in 1:n_seas) {

    west_y <- bio_row(bio$west, yrs[y], seas)
    weca_y <- bio_row(bio$weca, yrs[y], seas)
    mat_y <- bio_row(bio$mat, yrs[y], seas)

    for(r in 1:n_regions) {
      WAA[1, r, y, seas, , 1] <- west_y  # stock weight at age
      MatAA[1, r, y, seas, , 1] <- mat_y # proportion mature at age

      for(f in 1:n_fish) WAA_fish[1, r, y, seas, , 1, f] <- weca_y  # catch weight
      for(k in 1:n_srv) WAA_srv[1, r, y, seas, , 1, k] <- west_y

    } # end r loop
  } # end seas loop
} # end y loop

# SPoRC scales an annual rate by seasdur, so use the season sum to approximate for now (will add feature later ... )
natmort <- array(0, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes))

for(y in 1:n_yrs) {
  M_s1 <- bio_row(bio$M, yrs[y], 1)
  M_s2 <- bio_row(bio$M, yrs[y], 2)
  annual_M <- M_s1 + M_s2
  annual_M[1] <- M_s2[1] / seasdur[2]   # age 0 is only present in season 2
  for(r in 1:n_regions) natmort[1, r, y, , 1] <- annual_M

} # end y loop

# effort ----------------------------------------------------------------------
effort <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas))
effort_status <- array(NA, dim = c(n_regions, n_yrs, n_seas))

eff_raw$r <- match(eff_raw$area, regions)
eff_raw$yi <- match(eff_raw$year, 1983:2025)

for(i in seq_len(nrow(eff_raw))) {

  y <- eff_raw$yi[i]
  if(is.na(y)) next                    # outside the model window

  r <- eff_raw$r[i]
  s <- eff_raw$season[i]

  if(eff_raw$model[i] == effort_model) {
    effort[r, y, s] <- eff_raw[i, effort_col]
    effort_status[r, y, s] <- 1
  }

  if(eff_raw$model[i] == "Landings, no effort") effort_status[r, y, s] <- 2
  if(eff_raw$model[i] == "No landings, no effort") effort_status[r, y, s] <- 3

} # end i loop

# observation error keys ------------------------------------------------------

# setup obs error keys fishery - NA excludes, same numbers = share
sigmaCAA_key <- array(NA, dim = c(n_ages, n_sexes, n_fish))
counter <- 0 # setup counter

for(f in 1:n_fish) {
  # age 0 only in the fleets that fish it
  if(fit_age0_catch && fish_seas[f] >= rec_seas) {
    counter <- counter + 1
    sigmaCAA_key[1, 1, f] <- counter
  }
  counter <- counter + 1
  sigmaCAA_key[2:3, 1, f] <- counter # share age 1 and 2
  counter <- counter + 1
  sigmaCAA_key[4:5, 1, f] <- counter # share age 3 and above

} # end f loop

# setup obs error key survey
sigmaSrvIdxAA_key <- array(NA, dim = c(n_ages, n_sexes, n_srv))
for(k in 1:n_srv) {
  sigmaSrvIdxAA_key[1:2, 1, k] <- (2 * (k - 1)) + 1   # ages 0 and 1 share
  sigmaSrvIdxAA_key[3, 1, k] <- (2 * (k - 1)) + 2     # age 2 on its own
} # end k loop

ln_sigmaSrvIdxAA <- array(log(0.5), dim = c(n_ages, n_sexes, n_srv))

# save shit ------------------------------------------------------------------------
dims <- list(
  n_pop = n_pop,
  n_regions = n_regions,
  regions = regions,
  years = yrs,
  ages = ages,
  n_ages = n_ages,
  n_yrs = n_yrs,
  n_seas = n_seas,
  seasdur = seasdur,
  n_sexes = n_sexes,
  n_fish_fleets = n_fish,
  n_srv_fleets = n_srv,
  fish_names = fish_names,
  fish_region = fish_region,
  fish_seas = fish_seas,
  fish_fleet = fish_fleet,
  srv_names = srv_names,
  srv_seas = srv_seas,
  srv_t = srv_t,
  srv_type = srv_type
)

spatial <- list(
  dims = dims,
  ObsCatchAA = ObsCatchAA,
  UseCatchAA = UseCatchAA,
  ObsCatch = ObsCatch,
  UseCatch = UseCatch,
  ObsSrvIdxAA = ObsSrvIdxAA,
  UseSrvIdxAA = UseSrvIdxAA,
  ObsSrvIdx = ObsSrvIdx,
  UseSrvIdx = UseSrvIdx,
  WAA = WAA,
  MatAA = MatAA,
  WAA_fish = WAA_fish,
  WAA_srv = WAA_srv,
  natmort = natmort,
  effort = effort,
  effort_status = effort_status,
  nocatch = nocatch,
  sigmaCAA_key = sigmaCAA_key,
  sigmaSrvIdxAA_key = sigmaSrvIdxAA_key,
  ln_sigmaSrvIdxAA = ln_sigmaSrvIdxAA
)

saveRDS(spatial, here("dev", "spatial_sandeel", "inputs", "north_sea_sandeel_3r.rds"))
