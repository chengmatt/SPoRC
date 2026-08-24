# Build the SPoRC data list for the EBS Pacific cod bridge.
#
# Source: Stock Synthesis Model 24.1 of the 2024 eastern Bering Sea Pacific cod
# assessment (Barbeaux et al. 2024), NOVEMBER_MODELS/APPENDIX_2.3_2024_MODELS/
# Model_24.1 in the AFSC EBS_PCOD repository.
#
# Two functions, run in order:
#
#   build_pcod_data(run_dir)          parses the SS3 input files with r4ss into
#                                     SPoRC arrays and returns the dat list
#   add_pcod_ss3_report(dat, report)  attaches the run's reported quantities as
#                                     dat$ss3, which is what the regression test
#                                     compares against
#
# Nothing is typed in by hand; every number comes out of the SS3 files, so
# rerunning against a new model version picks the changes up. The packaged
# result is data/sgl_rg_ebs_pcod_data.rda.
#
# The assessment repository ships inputs only. Report.sso and wtatage.ss_new
# have to be produced by running SS3 in a copy of the Model 24.1 folder with
# init_values_src set to 1 and -maxfn 0 -nohess, which evaluates at ss.par
# without estimating and reproduces the assessment's own reported objective.
# dev/pcod_bridge/README.md has the build steps.
#
# What the SS3 files say, carried into the list as it stands:
#   ages 0-20 in the population, observed age bins 0-12 with 12 accumulating
#   two ageing error definitions, a biased reader through 2007, unbiased after
#   121 population length bins, compositions on 24 bins of 5 cm from 4.5
#   one sex, one area, spawning and settlement at the start of the year
#   the survey at month 7, the fishery length compositions at month 1
#
# The last of those is a trap. SS3 reads a fishing fleet's compositions at mid
# season whatever month the data file records, so month 1 is not the timing the
# model uses, and the bridge helper sets t_fish = 0.5 instead. See
# vignette("ae_ebs_pacific_cod_case_study") for what the numbers here become.

library(r4ss)

pcod_paths <- function(run_dir) {
  list(dat = file.path(run_dir, "BSPcod24_OCT_5cm_NB.dat"), ctl = file.path(run_dir, "Model_24.1.ctl"),
       starter = file.path(run_dir, "starter.ss"), forecast = file.path(run_dir, "forecast.ss"),
       par = file.path(run_dir, "ss.par"), dir = run_dir)
}

#' Stock Synthesis's population-to-data length bin map
#'
#' Each population bin goes to the data bin whose range holds its lower edge;
#' bins below the first data bin go into the first, bins at or above the last
#' data bin into the last.
#'
#' @param pop_lower lower edges of the population bins
#' @param dat_lower lower edges of the data bins
#' @return matrix [n_pop_bins x n_dat_bins] of zeros and ones
ss3_len_bin_map <- function(pop_lower, dat_lower) {
  map <- matrix(0, length(pop_lower), length(dat_lower))
  for(z in seq_along(pop_lower)) {
    k <- findInterval(pop_lower[z], dat_lower) # 0 below the first data bin
    map[z, max(1, k)] <- 1
  }
  map
}

#' Parse the SS3 inputs into SPoRC arrays
#'
#' @param run_dir directory holding the Model 24.1 files
#' @return a named list, the `dat` object the bridge helper builds the model from
build_pcod_data <- function(run_dir) {

  p <- pcod_paths(run_dir)
  d <- SS_readdat(p$dat, verbose = FALSE)
  c <- SS_readctl(p$ctl, datlist = d, verbose = FALSE)

  years <- d$styr:d$endyr
  n_yrs <- length(years)
  ages <- 0:d$Nages # population ages, 0 through the accumulator
  n_ages <- length(ages)
  obs_ages <- d$agebin_vector # 0..12
  n_obs_ages <- length(obs_ages)
  pop_lower <- d$lbin_vector_pop # 121 lower edges
  n_lens <- length(pop_lower)
  widths <- diff(pop_lower)
  pop_mid <- pop_lower + c(widths, widths[n_lens - 1]) / 2 # SS3's len_bins_m
  dat_lower <- d$lbin_vector # 24 lower edges
  n_obs_lens <- length(dat_lower)
  LenBinMap <- ss3_len_bin_map(pop_lower, dat_lower)
  startbin <- which(pop_lower >= dat_lower[1])[1] # first population bin in the data range
  n_sexes <- d$Nsexes
  n_regions <- d$N_areas
  fleetinfo <- d$fleetinfo
  fish_fleets <- which(fleetinfo$type == 1)
  srv_fleets <- which(fleetinfo$type == 3)
  n_fish <- length(fish_fleets); n_srv <- length(srv_fleets)
  fleet_area <- fleetinfo$area

  yidx <- function(y) match(y, years)

  # Catch (biomass) and the equilibrium catch -----------------------------------
  ObsCatch <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_fish))
  UseCatch <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))
  catch_se <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_fish))
  for(i in seq_len(nrow(d$catch))) {
    ci <- d$catch[i, ]
    if(ci$year < d$styr) next
    f <- match(ci$fleet, fish_fleets); r <- fleet_area[ci$fleet]
    ObsCatch[r, yidx(ci$year), 1, f] <- ci$catch
    UseCatch[r, yidx(ci$year), 1, f] <- 1
    catch_se[r, yidx(ci$year), 1, f] <- ci$catch_se
  }
  init_equil_catch <- d$catch$catch[d$catch$year < d$styr]
  init_equil_catch_se <- d$catch$catch_se[d$catch$year < d$styr]

  # Survey index (numbers, lognormal); the -2 fleet code in 2020 is a skipped row ----
  ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_srv))
  ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_srv))
  UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))
  cpue <- d$CPUE
  idx_fleet_col <- grep("^index$|^fleet$|FltSvy", names(cpue), value = TRUE)[1]
  for(i in seq_len(nrow(cpue))) {
    fcode <- cpue[[idx_fleet_col]][i]; sf <- match(abs(fcode), srv_fleets); r <- fleet_area[abs(fcode)]
    ObsSrvIdx[r, yidx(cpue$year[i]), 1, sf] <- cpue$obs[i]
    ObsSrvIdx_SE[r, yidx(cpue$year[i]), 1, sf] <- cpue$se_log[i]
    UseSrvIdx[r, yidx(cpue$year[i]), 1, sf] <- as.numeric(fcode > 0 & cpue$obs[i] > 0)
  }
  srv_month <- unique(cpue$month[cpue[[idx_fleet_col]] > 0])
  t_srv <- (srv_month - 1) / 12 # fraction of the year at the survey

  # Length compositions on the 24 data bins; negative fleet codes are ghosts ------
  lc <- d$lencomp
  lc_fleet_col <- grep("^fleet$|FltSvy", names(lc), value = TRUE)[1]
  lc_n_col <- grep("^Nsamp$", names(lc), value = TRUE)[1]
  lc_vals <- function(row) as.numeric(row[(which(names(lc) == lc_n_col) + 1):ncol(lc)])
  ObsFishLenComps <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_obs_lens, n_sexes, n_fish))
  ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fish))
  UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))
  ObsSrvLenComps <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_obs_lens, n_sexes, n_srv))
  ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_srv))
  UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))
  fish_len_month <- NA; srv_len_month <- NA
  for(i in seq_len(nrow(lc))) {
    fcode <- lc[[lc_fleet_col]][i]; fl <- abs(fcode); y <- yidx(abs(lc$year[i])); r <- fleet_area[fl]
    use_row <- as.numeric(fcode > 0 & lc$year[i] > 0)
    if(is.na(y)) next
    v <- lc_vals(lc[i, ])[1:n_obs_lens]
    if(fl %in% fish_fleets) {
      f <- match(fl, fish_fleets); fish_len_month <- lc$month[i]
      # a ghost row never overwrites a fitted one
      if(use_row == 1 || UseFishLenComps[r, y, 1, f] == 0) {
        ObsFishLenComps[r, y, 1, , 1, f] <- v
        ISS_FishLenComps[r, y, 1, , f] <- lc[[lc_n_col]][i]
        UseFishLenComps[r, y, 1, f] <- max(UseFishLenComps[r, y, 1, f], use_row)
      }
    } else {
      sf <- match(fl, srv_fleets); if(use_row == 1) srv_len_month <- lc$month[i]
      if(use_row == 1 || UseSrvLenComps[r, y, 1, sf] == 0) {
        ObsSrvLenComps[r, y, 1, , 1, sf] <- v
        ISS_SrvLenComps[r, y, 1, , sf] <- lc[[lc_n_col]][i]
        UseSrvLenComps[r, y, 1, sf] <- max(UseSrvLenComps[r, y, 1, sf], use_row)
      }
    }
  }
  t_fish_comp <- (fish_len_month - 1) / 12

  # Age compositions: the survey's marginal ages (full length range) with their
  # ageing error definition by year; the conditional age-at-length rows are ghosts
  ac <- d$agecomp
  ac_fleet_col <- grep("^fleet$|FltSvy", names(ac), value = TRUE)[1]
  ac_n_col <- grep("^Nsamp$", names(ac), value = TRUE)[1]
  ac_vals <- function(row) as.numeric(row[(which(names(ac) == ac_n_col) + 1):ncol(ac)])
  ObsSrvAgeComps <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_obs_ages, n_sexes, n_srv))
  ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_srv))
  UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))
  ageerr_by_year <- rep(NA_integer_, n_yrs)
  for(i in seq_len(nrow(ac))) {
    fcode <- ac[[ac_fleet_col]][i]; fl <- abs(fcode); y <- yidx(abs(ac$year[i])); r <- fleet_area[fl]
    use <- as.numeric(fcode > 0 & ac$year[i] > 0)
    if(!(fl %in% srv_fleets)) next
    if(ac$Lbin_lo[i] != min(ac$Lbin_lo) || ac$Lbin_hi[i] != max(ac$Lbin_hi)) next # a conditional row
    sf <- match(fl, srv_fleets)
    if(use == 1 || UseSrvAgeComps[r, y, 1, sf] == 0) {
      ObsSrvAgeComps[r, y, 1, , 1, sf] <- ac_vals(ac[i, ])[1:n_obs_ages]
      ISS_SrvAgeComps[r, y, 1, , sf] <- ac[[ac_n_col]][i]
      UseSrvAgeComps[r, y, 1, sf] <- max(UseSrvAgeComps[r, y, 1, sf], use)
      if(use == 1) ageerr_by_year[y] <- ac$ageerr[i]
    }
  }

  # Ageing error: the two definitions as SS3 forms them (get_age_age), the
  # normal CDF of the read age at the integer edges of observed bins 2..n with
  # the given mean (a biased reader's expected age, or true age + 0.5) and SD,
  # tails into the end bins
  ae <- d$ageerror
  n_defs <- nrow(ae) / 2
  AgeingError_defs <- vector("list", n_defs)
  for(k in seq_len(n_defs)) {
    ae_mean <- as.numeric(ae[2 * k - 1, ]); ae_sd <- as.numeric(ae[2 * k, ])
    ae_mean[ae_mean < 0] <- ages[ae_mean < 0] + 0.5
    M <- matrix(0, n_ages, n_obs_ages)
    edges <- c(-Inf, obs_ages[-1], Inf)
    for(a in seq_len(n_ages)) M[a, ] <- diff(pnorm(edges, ae_mean[a], ae_sd[a]))
    AgeingError_defs[[k]] <- M
  }
  # years without an age composition take the last definition in use
  last_def <- n_defs
  AgeingError <- array(0, dim = c(n_yrs, n_ages, n_obs_ages))
  for(y in seq_len(n_yrs)) {
    if(!is.na(ageerr_by_year[y])) last_def <- ageerr_by_year[y]
    AgeingError[y, , ] <- AgeingError_defs[[last_def]]
  }

  # Biology and growth from the control file ------------------------------------
  mg <- c$MG_parms
  grow <- list(
    M = mg[grep("^NatM", rownames(mg)), "INIT"],
    L1 = mg[grep("L_at_Amin", rownames(mg)), "INIT"], L2 = mg[grep("L_at_Amax", rownames(mg)), "INIT"],
    K = mg[grep("VonBert_K", rownames(mg)), "INIT"], rho = mg[grep("Richards", rownames(mg)), "INIT"],
    CV1 = mg[grep("CV_young", rownames(mg)), "INIT"], CV2 = mg[grep("CV_old", rownames(mg)), "INIT"],
    bounds = rbind(L1 = as.numeric(mg[grep("L_at_Amin", rownames(mg)), c("LO", "HI")]),
                   L2 = as.numeric(mg[grep("L_at_Amax", rownames(mg)), c("LO", "HI")]),
                   K = as.numeric(mg[grep("VonBert_K", rownames(mg)), c("LO", "HI")]),
                   CV1 = as.numeric(mg[grep("CV_young", rownames(mg)), c("LO", "HI")]),
                   CV2 = as.numeric(mg[grep("CV_old", rownames(mg)), c("LO", "HI")]),
                   rho = as.numeric(mg[grep("Richards", rownames(mg)), c("LO", "HI")])),
    est = c(L1 = TRUE, L2 = TRUE, K = TRUE, CV1 = FALSE, CV2 = FALSE, rho = TRUE),
    dev_se = c(L1 = mg[grep("L_at_Amin", rownames(mg)), "dev_link"] > 0, K = mg[grep("VonBert_K", rownames(mg)), "dev_link"] > 0),
    dev_years = list(L1 = mg[grep("L_at_Amin", rownames(mg)), "dev_minyr"]:mg[grep("L_at_Amin", rownames(mg)), "dev_maxyr"],
                     K = mg[grep("VonBert_K", rownames(mg)), "dev_minyr"]:mg[grep("VonBert_K", rownames(mg)), "dev_maxyr"]),
    A1 = c$Growth_Age_for_L1, A2 = c$Growth_Age_for_L2
  )
  grow$est["rho"] <- mg[grep("Richards", rownames(mg)), "PHASE"] > 0
  tv <- c$MG_parms_tv
  grow$dev_sd <- c(L1 = tv[grep("L_at_Amin.*dev_se", rownames(tv)), "INIT"],
                   K = tv[grep("(VonBert_K|Richards).*dev_se", rownames(tv))[1], "INIT"])
  wtlen <- mg[grep("^Wtlen_[12]", rownames(mg)), "INIT"][1:2]
  mat <- list(L50 = mg[grep("^Mat50%", rownames(mg)), "INIT"], slope = mg[grep("^Mat_slope", rownames(mg)), "INIT"])

  sr <- c$SR_parms
  rec <- list(ln_R0 = sr["SR_LN(R0)", "INIT"], h = sr["SR_BH_steep", "INIT"], sigmaR = sr["SR_sigmaR", "INIT"],
              regime = c$SR_parms_tv[1, "INIT"],
              main_first = c$MainRdevYrFirst, main_last = c$MainRdevYrLast,
              early_start = c$recdev_early_start,
              bias_years = c(c$last_early_yr_nobias_adj, c$first_yr_fullbias_adj, c$last_yr_fullbias_adj, c$first_recent_yr_nobias_adj),
              max_bias_adj = c$max_bias_adj)

  sel <- list(size = c$size_selex_parms, size_tv = c$size_selex_parms_tv, size_types = c$size_selex_types,
              blocks = c$Block_Design, blocks_per_pattern = c$blocks_per_pattern)
  q <- list(options = c$Q_options, parms = c$Q_parms)
  var_adj <- c$Variance_adjustment_list
  init_F <- c$init_F

  list(
    source = run_dir, years = years, ages = ages, obs_ages = obs_ages,
    lens_lower = pop_lower, lens = pop_mid, dat_lens_lower = dat_lower, LenBinMap = LenBinMap, startbin = startbin,
    n_regions = n_regions, n_sexes = n_sexes, n_fish_fleets = n_fish, n_srv_fleets = n_srv,
    fleetnames = fleetinfo$fleetname, fish_fleets = fish_fleets, srv_fleets = srv_fleets, fleet_area = fleet_area,
    spawn_month = d$spawn_month, n_subseas = d$Nsubseasons, t_srv = t_srv, t_fish_comp = t_fish_comp,
    ObsCatch = ObsCatch, UseCatch = UseCatch, catch_se = catch_se,
    init_equil_catch = init_equil_catch, init_equil_catch_se = init_equil_catch_se, init_F = init_F,
    ObsSrvIdx = ObsSrvIdx, ObsSrvIdx_SE = ObsSrvIdx_SE, UseSrvIdx = UseSrvIdx,
    ObsFishLenComps = ObsFishLenComps, ISS_FishLenComps = ISS_FishLenComps, UseFishLenComps = UseFishLenComps,
    ObsSrvLenComps = ObsSrvLenComps, ISS_SrvLenComps = ISS_SrvLenComps, UseSrvLenComps = UseSrvLenComps,
    ObsSrvAgeComps = ObsSrvAgeComps, ISS_SrvAgeComps = ISS_SrvAgeComps, UseSrvAgeComps = UseSrvAgeComps,
    AgeingError = AgeingError, AgeingError_defs = AgeingError_defs, ageerr_by_year = ageerr_by_year, ageerror_raw = ae,
    growth = grow, wtlen = wtlen, mat = mat, rec = rec,
    sel = sel, q = q, var_adj = var_adj,
    comp = list(addtocomp_len = d$len_info$addtocomp[1], addtocomp_age = d$age_info$addtocomp[1]),
    ctl = c, dat_raw = d
  )
}

if(sys.nframe() == 0) {
  run_dir <- "/Users/matthewcheng/Downloads/EBS_PCOD-main/2024_ASSESSMENT/NOVEMBER_MODELS/APPENDIX_2.3_2024_MODELS/Model_24.1"
  dat <- build_pcod_data(run_dir)
  out <- file.path("dev/pcod_bridge/output/pcod_bridge_data.rds")
  saveRDS(dat, out)
  cat("wrote", out, "\n")
  cat("length comp years fish/srv:", sum(dat$UseFishLenComps), sum(dat$UseSrvLenComps), " survey age years:", sum(dat$UseSrvAgeComps), " index years:", sum(dat$UseSrvIdx), "\n")
  cat("startbin", dat$startbin, " t_fish_comp", dat$t_fish_comp, " t_srv", dat$t_srv, "\n")
  str(dat$growth[c("L1", "L2", "K", "rho", "CV1", "CV2", "dev_sd", "A1", "A2")])
}


#' Attach the SS3 report quantities the bridge test compares against
#'
#' Read SS3's year-by-year weight and fecundity at age
#'
#' `wtatage.ss_new` carries one row per year and "fleet", where fleet 0 is the
#' population weight at the start of the season, -1 the mid-season weight, -2
#' the fecundity (maturity times weight at the spawning time), and a positive
#' fleet its own selection-weighted body weight.
#'
#' @param path wtatage.ss_new written by the run
#' @param years model years
#' @param ages model ages
#' @return list of matrices [year x age], one per fleet code present
read_ss3_wtatage <- function(path, years, ages) {
  ln <- readLines(path)
  ln <- ln[grepl("^[0-9]", ln)]
  m <- do.call(rbind, lapply(strsplit(trimws(sub("#.*$", "", ln)), " +"), as.numeric))
  m <- m[m[, 1] %in% years, , drop = FALSE]
  out <- list()
  for(fl in unique(m[, 6])) {
    sub <- m[m[, 6] == fl, , drop = FALSE]
    tab <- matrix(NA_real_, length(years), length(ages), dimnames = list(years, ages))
    tab[match(sub[, 1], years), ] <- sub[, 7:(6 + length(ages))]
    out[[as.character(fl)]] <- tab
  }
  out
}

#' Attach the SS3 report quantities the bridge test compares against
#'
#' @param dat the list from build_pcod_data()
#' @param report the r4ss::SS_output() list for the same run
#' @param par_file path to the ss.par the report was evaluated at (12 digits)
#' @param wtatage_file path to wtatage.ss_new, read for maturity at age
#' @return dat with a `$ss3` element and an `$mle` element of parameter values
add_pcod_ss3_report <- function(dat, report, par_file = NULL, wtatage_file = NULL) {

  r <- report
  years <- dat$years; n_yrs <- length(years); n_ages <- length(dat$ages)
  p <- r$parameters
  pv <- function(lab) { v <- p$Value[p$Label == lab]; if(length(v) != 1) stop("parameter not found: ", lab); v }

  # the par file carries twelve digits against the report's six
  pf <- if(!is.null(par_file)) readLines(par_file) else NULL
  par_block <- function(name) {
    if(is.null(pf)) return(NULL)
    i <- which(pf == paste0("# ", name, ":"))
    if(length(i) != 1) return(NULL)
    as.numeric(strsplit(trimws(pf[i + 1]), " +")[[1]])
  }
  mgp <- sapply(1:19, function(k) par_block(paste0("MGparm[", k, "]")))
  srp <- sapply(1:6, function(k) par_block(paste0("SR_parm[", k, "]")))
  selp <- sapply(1:16, function(k) par_block(paste0("selparm[", k, "]")))

  growth_mle <- c(L1 = mgp[2], L2 = mgp[3], K = mgp[4], CV1 = mgp[6], CV2 = mgp[7], rho = mgp[5])
  sel_mle <- list(
    fishery = list(base = selp[1:6], block1 = selp[13:14]), # peak and ascending width replaced in 1977-1989
    survey = list(base = selp[7:12], dev_sd = selp[15])
  )
  rec <- r$recruit
  main <- p[grep("^Main_RecrDev_", p$Label), ]
  early <- p[grep("^Early_InitAge_", p$Label), ]
  early_yr <- years[1] - as.integer(sub("Early_InitAge_", "", early$Label))
  mle <- list(
    ln_R0 = srp[1], sigmaR = srp[3], regime = srp[6],
    ln_q = par_block("Q_parm[1]"), extra_sd = par_block("Q_parm[2]"),
    init_F = par_block("init_F[1]"),
    growth = growth_mle,
    growth_devs = list(L1 = par_block("parm_dev[1]"), K = par_block("parm_dev[2]")),
    srv_sel_devs = par_block("parm_dev[3]"),
    sel = sel_mle,
    main_recdev = setNames(par_block("recdev1"), seq(dat$rec$main_first, dat$rec$main_last)),
    early_recdev = setNames(par_block("recdev_early"), seq(years[1] + dat$rec$early_start, years[1] - 1)),
    # recruitment as realized, so the seeds need no convention about where the bias term sits
    rec_ratio = setNames(rec$pred_recr / rec$exp_recr, rec$Yr),
    biasadj = setNames(rec$biasadjuster, rec$Yr),
    # the fleet's own fishing mortality (the rate multiplying selectivity), not
    # the report's annual_F, which is an exploitation-rate summary
    Fmort = setNames(r$exploitation[[dat$fleetnames[dat$fish_fleets[1]]]], r$exploitation$Yr)
  )

  # Derived quantities ----------------------------------------------------------
  na <- r$natage
  NAA <- NAA_mid <- array(NA_real_, dim = c(n_yrs, n_ages))
  rows <- na[na$Area == 1 & na$Sex == 1 & na[["Beg/Mid"]] == "B" & na$Era == "TIME", ]
  NAA[match(rows$Yr, years), ] <- as.matrix(rows[, as.character(dat$ages)])
  rows <- na[na$Area == 1 & na$Sex == 1 & na[["Beg/Mid"]] == "M" & na$Era == "TIME", ]
  NAA_mid[match(rows$Yr, years), ] <- as.matrix(rows[, as.character(dat$ages)])
  ts <- r$timeseries
  t <- ts[ts$Area == 1 & ts$Era == "TIME", ]
  SSB <- t$SpawnBio[match(years, t$Yr)]; Rec <- t$Recruit_0[match(years, t$Yr)]; Bio_all <- t$Bio_all[match(years, t$Yr)]
  dead_B <- t[["dead(B):_1"]][match(years, t$Yr)]
  init <- ts[ts$Era == "INIT", ]

  mg <- r$MGparmAdj
  growth_by_year <- data.frame(Yr = mg$Yr, L1 = mg$L_at_Amin_Fem_GP_1, L2 = mg$L_at_Amax_Fem_GP_1, K = mg$VonBert_K_Fem_GP_1, rho = mg$Richards_Fem_GP_1)
  gs <- r$growthseries
  as <- r$ageselex
  asel2 <- lapply(1:2, function(f) { e <- as[as$Fleet == f & as$Factor == "Asel2", ]; m <- as.matrix(e[, as.character(dat$ages)]); rownames(m) <- e$Yr; m })
  bodywt <- { e <- as[as$Fleet == 1 & as$Factor == "bodywt", ]; m <- as.matrix(e[, as.character(dat$ages)]); rownames(m) <- e$Yr; m }
  ss <- r$sizeselex
  lsel <- lapply(1:2, function(f) { e <- ss[ss$Fleet == f & ss$Factor == "Lsel", ]; m <- as.matrix(e[, as.character(dat$lens)]); rownames(m) <- e$Yr; m })

  # Maturity at age, as the share of the weight at age that spawns, taken
  # straight from SS3's own fecundity and start-of-season weight. Maturity is
  # length based and fixed in the assessment, so this is data to SPoRC while
  # weight at age stays derived from the estimated growth.
  MatAA <- fecundity <- popwt_beg <- NULL
  if(!is.null(wtatage_file) && file.exists(wtatage_file)) {
    wa <- read_ss3_wtatage(wtatage_file, years, dat$ages)
    fecundity <- wa[["-2"]]; popwt_beg <- wa[["0"]]
    MatAA <- fecundity / popwt_beg
    MatAA[!is.finite(MatAA)] <- 0
    dat$wtatage <- wa
    dat$MatAA <- MatAA
  }

  dat$ss3 <- list(
    likelihoods = r$likelihoods_used, likelihoods_by_fleet = r$likelihoods_by_fleet, parm_devs = r$Parm_devs_detail,
    NAA = NAA, NAA_mid = NAA_mid, SSB = SSB, Rec = Rec, Bio_all = Bio_all, dead_B = dead_B,
    SSB_virgin = r$derived_quants$Value[r$derived_quants$Label == "SSB_Virgin"],
    SSB_initial = r$derived_quants$Value[r$derived_quants$Label == "SSB_Initial"],
    init = init[, c("Yr", "Bio_all", "SpawnBio", "Recruit_0")],
    cpue = r$cpue[, intersect(c("Fleet", "Yr", "Obs", "Exp", "Calc_Q", "SE", "SE_input", "Like"), names(r$cpue))],
    growth_by_year = growth_by_year, growthseries = gs, endgrowth = r$endgrowth,
    asel2 = asel2, bodywt = bodywt, lsel = lsel, AAK = r$AAK,
    agedbase = r$agedbase[, intersect(c("Yr", "Fleet", "Bin", "Obs", "Exp", "Nsamp_adj", "Nsamp_in", "Like"), names(r$agedbase))],
    lendbase = r$lendbase[, intersect(c("Yr", "Fleet", "Bin", "Obs", "Exp", "Nsamp_adj", "Nsamp_in", "Like"), names(r$lendbase))],
    recruit = rec[, intersect(c("Yr", "era", "exp_recr", "bias_adjusted", "pred_recr", "SpawnBio", "dev", "biasadjuster"), names(rec))],
    version = r$SS_version
  )
  dat$mle <- mle
  dat
}

if(sys.nframe() == 0) {
  rep_file <- "dev/pcod_bridge/output/pcod_ss3_report.rds"
  if(file.exists(rep_file)) {
    wt_file <- "dev/pcod_bridge/output/wtatage.ss_new"
    dat <- add_pcod_ss3_report(dat, readRDS(rep_file), par_file = file.path(run_dir, "ss.par"), wtatage_file = wt_file)
    saveRDS(dat, out)
    cat("attached SS3 report; ln_R0", dat$mle$ln_R0, " regime", dat$mle$regime, " ln_q", dat$mle$ln_q, " init_F", dat$mle$init_F, "\n")
    cat("growth MLE:", dat$mle$growth, "\n")
    cat("fishery sel base:", dat$mle$sel$fishery$base, " block:", dat$mle$sel$fishery$block1, "\n")
    cat("survey sel:", dat$mle$sel$survey$base, "\n")
    if(!is.null(dat$MatAA)) cat("maturity at age (1977, ages 0-6):", signif(dat$MatAA[1, 1:7], 4), "\n")
  }
}
