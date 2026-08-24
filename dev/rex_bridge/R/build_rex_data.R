# Build the SPoRC-ready data list for the GOA rex sole bridge (SS3 Model 25.1,
# runs/2025_models/run_14_asfor7_ageing_error in Carey McGilliard's goa_rex repo).
#
# Everything is parsed with r4ss from the SS3 input files, so nothing here is typed
# in by hand. Report.sso, when available, is read by add_rex_ss3_report() and
# attached as dat$ss3 with the quantities the regression test compares against,
# following the West Coast sablefish bridge (tests/testthat/helper-bridge_wc_sablefish.R).
#
# SS3 conventions carried into the list:
#   ages 0-20 in the population, observed age bins 1-20 (age 0 maps into bin 1
#   through the ageing error matrix), length bins 9-65 by 2 with population bins
#   equal to data bins, two areas with no movement and one growth pattern per area,
#   all compositions at month 7 (mid year), spawning at month 1.

library(r4ss)

rex_paths <- function(run_dir) {
  list(dat = file.path(run_dir, "goa_rex.dat"), ctl = file.path(run_dir, "goa_rex.ctl"),
       starter = file.path(run_dir, "starter.ss"), forecast = file.path(run_dir, "forecast.ss"),
       dir = run_dir)
}

#' Parse the SS3 inputs into SPoRC arrays
#'
#' @param run_dir directory holding goa_rex.dat and goa_rex.ctl
#' @return a named list, the `dat` object the bridge helper builds the model from
build_rex_data <- function(run_dir) {

  p <- rex_paths(run_dir)
  d <- SS_readdat(p$dat, verbose = FALSE)
  c <- SS_readctl(p$ctl, datlist = d, verbose = FALSE)

  years <- d$styr:d$endyr
  n_yrs <- length(years)
  ages <- 0:d$Nages                     # population ages, 0 through the accumulator
  n_ages <- length(ages)
  obs_ages <- d$agebin_vector            # 1..20
  n_obs_ages <- length(obs_ages)
  lens_lower <- d$lbin_vector            # lower edges, 9..65 by 2
  n_lens <- length(lens_lower)
  lens_mid <- lens_lower + d$binwidth / 2
  n_sexes <- d$Nsexes
  n_regions <- d$N_areas
  fleetinfo <- d$fleetinfo
  fish_fleets <- which(fleetinfo$type == 1)
  srv_fleets <- which(fleetinfo$type == 3)
  n_fish <- length(fish_fleets); n_srv <- length(srv_fleets)
  fleet_area <- fleetinfo$area

  yidx <- function(y) match(y, years)

  # Catch (biomass, one fishery in area 1) --------------------------------------
  ObsCatch <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_fish))
  UseCatch <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))
  catch_se <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_fish))
  for(i in seq_len(nrow(d$catch))) {
    ci <- d$catch[i, ]
    if(ci$year < d$styr) next # the -999 equilibrium catch row
    f <- match(ci$fleet, fish_fleets); r <- fleet_area[ci$fleet]
    ObsCatch[r, yidx(ci$year), 1, f] <- ci$catch
    UseCatch[r, yidx(ci$year), 1, f] <- 1
    catch_se[r, yidx(ci$year), 1, f] <- ci$catch_se
  }
  init_equil_catch <- d$catch$catch[d$catch$year < d$styr]

  # Survey indices (biomass, lognormal); negative fleet codes are excluded years ----
  ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_srv))
  ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_srv))
  UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))
  cpue <- d$CPUE
  idx_fleet_col <- grep("^index$|^fleet$|FltSvy", names(cpue), value = TRUE)[1]
  for(i in seq_len(nrow(cpue))) {
    fcode <- cpue[[idx_fleet_col]][i]; sf <- match(abs(fcode), srv_fleets); r <- fleet_area[abs(fcode)]
    ObsSrvIdx[r, yidx(cpue$year[i]), 1, sf] <- cpue$obs[i]
    ObsSrvIdx_SE[r, yidx(cpue$year[i]), 1, sf] <- cpue$se_log[i]
    UseSrvIdx[r, yidx(cpue$year[i]), 1, sf] <- as.numeric(fcode > 0)
  }
  srv_month <- fleetinfo$surveytiming[srv_fleets]
  t_srv <- (srv_month - 1) / 12 # fraction of the year at the survey

  # Length compositions: joint sex rows (female then male), proportions ----------
  lc <- d$lencomp
  lc_fleet_col <- grep("^fleet$|FltSvy", names(lc), value = TRUE)[1]
  lc_sex_col <- grep("^sex$|Gender", names(lc), value = TRUE)[1]
  lc_n_col <- grep("^Nsamp$", names(lc), value = TRUE)[1]
  lc_vals <- function(row) as.numeric(row[(which(names(lc) == lc_n_col) + 1):ncol(lc)])
  ObsFishLenComps <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_fish))
  ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fish))
  UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))
  ObsSrvLenComps <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_srv))
  ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_srv))
  UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))
  len_sex <- list(fish = rep(NA, n_fish), srv = rep(NA, n_srv))
  for(i in seq_len(nrow(lc))) {
    # a negative fleet or a negative year is a row SS3 carries but does not fit
    fcode <- lc[[lc_fleet_col]][i]; fl <- abs(fcode); y <- yidx(abs(lc$year[i])); r <- fleet_area[fl]
    use_row <- as.numeric(fcode > 0 & lc$year[i] > 0)
    if(is.na(y)) next # a year outside the model (the fishery carries three 1977 rows)
    v <- lc_vals(lc[i, ]); sx <- lc[[lc_sex_col]][i]
    vals <- if(sx == 3) matrix(v, nrow = n_lens, ncol = n_sexes) else matrix(v[1:n_lens], nrow = n_lens, ncol = 1)
    if(fl %in% fish_fleets) {
      f <- match(fl, fish_fleets); len_sex$fish[f] <- sx
      ObsFishLenComps[r, y, 1, , , f] <- if(sx == 3) vals else cbind(vals, 0)
      ISS_FishLenComps[r, y, 1, , f] <- lc[[lc_n_col]][i]
      UseFishLenComps[r, y, 1, f] <- use_row
    } else {
      sf <- match(fl, srv_fleets); len_sex$srv[sf] <- sx
      ObsSrvLenComps[r, y, 1, , , sf] <- if(sx == 3) vals else cbind(vals, 0)
      ISS_SrvLenComps[r, y, 1, , sf] <- lc[[lc_n_col]][i]
      UseSrvLenComps[r, y, 1, sf] <- use_row
    }
  }

  # Age compositions: fishery marginal (joint sex, Lbin -1) and survey CAAL -----
  # (sex-specific rows, one length bin each, counts of otoliths). Ghost fleets
  # (negative) are stored but never used.
  ac <- d$agecomp
  ac_fleet_col <- grep("^fleet$|FltSvy", names(ac), value = TRUE)[1]
  ac_sex_col <- grep("^sex$|Gender", names(ac), value = TRUE)[1]
  ac_n_col <- grep("^Nsamp$", names(ac), value = TRUE)[1]
  ac_vals <- function(row) as.numeric(row[(which(names(ac) == ac_n_col) + 1):ncol(ac)])
  ObsFishAgeComps <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_obs_ages, n_sexes, n_fish))
  ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fish))
  UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))
  ObsSrvAgeComps <- array(NA_real_, dim = c(n_regions, n_yrs, 1, n_obs_ages, n_sexes, n_srv))
  ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_srv))
  UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))
  ObsSrv_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_obs_ages, n_sexes, n_srv))
  ISS_Srv_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_srv))
  UseSrv_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_srv))
  ObsFish_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_obs_ages, n_sexes, n_fish))
  ISS_Fish_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_fish))
  UseFish_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_fish))
  age_sex <- list(fish = rep(NA, n_fish), srv = rep(NA, n_srv))
  for(i in seq_len(nrow(ac))) {
    fcode <- ac[[ac_fleet_col]][i]; fl <- abs(fcode); y <- yidx(abs(ac$year[i])); r <- fleet_area[fl]
    use <- as.numeric(fcode > 0 & ac$year[i] > 0)
    v <- ac_vals(ac[i, ]); sx <- ac[[ac_sex_col]][i]; nn <- ac[[ac_n_col]][i]
    lo <- ac$Lbin_lo[i]
    if(lo < 0) {
      # marginal age composition across all lengths
      vals <- if(sx == 3) matrix(v, nrow = n_obs_ages, ncol = n_sexes) else matrix(v[1:n_obs_ages], nrow = n_obs_ages, ncol = 1)
      if(fl %in% fish_fleets) {
        f <- match(fl, fish_fleets); age_sex$fish[f] <- sx
        ObsFishAgeComps[r, y, 1, , , f] <- if(sx == 3) vals else cbind(vals, 0)
        ISS_FishAgeComps[r, y, 1, , f] <- nn; UseFishAgeComps[r, y, 1, f] <- use
      } else {
        sf <- match(fl, srv_fleets); age_sex$srv[sf] <- sx
        ObsSrvAgeComps[r, y, 1, , , sf] <- if(sx == 3) vals else cbind(vals, 0)
        ISS_SrvAgeComps[r, y, 1, , sf] <- nn; UseSrvAgeComps[r, y, 1, sf] <- use
      }
    } else {
      # conditional age at length: one length bin (Lbin_lo == Lbin_hi, lengths), one sex
      if(ac$Lbin_hi[i] != lo) stop("CAAL row spans more than one length bin; not handled")
      l <- match(lo, lens_lower)
      if(is.na(l)) stop("CAAL Lbin_lo not a length bin edge: ", lo)
      # sex codes 1 (female) / 2 (male): counts sit in that sex's block
      counts <- if(sx == 1) v[1:n_obs_ages] else if(sx == 2) v[n_obs_ages + 1:n_obs_ages] else stop("joint-sex CAAL row not handled")
      if(fl %in% fish_fleets) {
        f <- match(fl, fish_fleets)
        ObsFish_caal[r, y, 1, l, , sx, f] <- counts
        ISS_Fish_caal[r, y, 1, l, sx, f] <- nn
        UseFish_caal[r, y, 1, l, f] <- max(UseFish_caal[r, y, 1, l, f], use)
      } else {
        sf <- match(fl, srv_fleets)
        ObsSrv_caal[r, y, 1, l, , sx, sf] <- counts
        ISS_Srv_caal[r, y, 1, l, sx, sf] <- nn
        UseSrv_caal[r, y, 1, l, sf] <- max(UseSrv_caal[r, y, 1, l, sf], use)
      }
    }
  }

  # Ageing error: SS3 gives mean (-1 = unbiased, read as true age + 0.5) and SD of the
  # observed age for each true age 0..20. SS3's get_age_age() takes the normal CDF at
  # the lower edges of observed bins 2..n (the integer ages), so bin b holds
  # [age_b, age_b+1), with the lower tail in the first bin and the upper in the last.
  ae_sd <- as.numeric(d$ageerror[2, ])
  ae_mean <- as.numeric(d$ageerror[1, ]); ae_mean[ae_mean < 0] <- ages[ae_mean < 0] + 0.5
  AgeingError <- matrix(0, n_ages, n_obs_ages)
  edges <- c(-Inf, obs_ages[-1], Inf)
  for(a in seq_len(n_ages)) {
    cdf <- pnorm(edges, ae_mean[a], ae_sd[a])
    AgeingError[a, ] <- diff(cdf)
  }

  # Biology and growth from the control file ------------------------------------
  mg <- c$MG_parms
  getp <- function(pattern) mg[grep(pattern, rownames(mg)), ]
  growth <- list()
  for(s in 1:n_sexes) for(gp in 1:n_regions) {
    sx <- c("Fem", "Mal")[s]
    rows <- mg[grep(paste0("_", sx, "_GP_", gp, "$"), rownames(mg)), ]
    growth[[paste0("s", s, "_gp", gp)]] <- list(
      M = rows["NatM_p_1_" %+% sx %+% "_GP_" %+% gp, "INIT"],
      L1 = rows[grep("L_at_Amin", rownames(rows)), "INIT"], L2 = rows[grep("L_at_Amax", rownames(rows)), "INIT"],
      K = rows[grep("VonBert_K", rownames(rows)), "INIT"], CV1 = rows[grep("CV_young", rownames(rows)), "INIT"],
      CV2 = rows[grep("CV_old", rownames(rows)), "INIT"],
      est = rows[grep("L_at_Amin|L_at_Amax|VonBert_K|CV_young|CV_old", rownames(rows)), "PHASE"] > 0)
  }
  wtlen <- list(
    fem = mg[grep("^Wtlen_[12]_Fem", rownames(mg)), "INIT"][1:2],
    mal = mg[grep("^Wtlen_[12]_Mal", rownames(mg)), "INIT"][1:2]
  )
  # maturity is read from the first growth pattern; the assessment gives every pattern the same curve
  mat <- list(a50 = mg[grep("^Mat50%_Fem", rownames(mg))[1], "INIT"], slope = mg[grep("^Mat_slope_Fem", rownames(mg))[1], "INIT"],
              first_mature_age = c$First_Mature_Age)
  rec_dist <- mg[grep("^RecrDist", rownames(mg)), c("INIT", "PHASE")]

  sr <- c$SR_parms
  rec <- list(ln_R0 = sr["SR_LN(R0)", "INIT"], h = sr["SR_BH_steep", "INIT"], sigmaR = sr["SR_sigmaR", "INIT"],
              main_first = c$MainRdevYrFirst, main_last = c$MainRdevYrLast,
              early_start = c$recdev_early_start, early_phase = c$recdev_early_phase,
              bias_years = c(c$last_early_yr_nobias_adj, c$first_yr_fullbias_adj, c$last_yr_fullbias_adj, c$first_recent_yr_nobias_adj),
              max_bias_adj = c$max_bias_adj)

  sel <- list(age = c$age_selex_parms, age_types = c$age_selex_types, size_types = c$size_selex_types)
  q <- list(options = c$Q_options, parms = c$Q_parms)
  var_adj <- c$Variance_adjustment_list

  list(
    source = run_dir, years = years, ages = ages, obs_ages = obs_ages, lens_lower = lens_lower, lens = lens_mid,
    n_regions = n_regions, n_sexes = n_sexes, n_fish_fleets = n_fish, n_srv_fleets = n_srv,
    fleetnames = fleetinfo$fleetname, fish_fleets = fish_fleets, srv_fleets = srv_fleets, fleet_area = fleet_area,
    spawn_month = d$spawn_month, n_subseas = d$Nsubseasons, t_srv = t_srv,
    ObsCatch = ObsCatch, UseCatch = UseCatch, catch_se = catch_se, init_equil_catch = init_equil_catch,
    ObsSrvIdx = ObsSrvIdx, ObsSrvIdx_SE = ObsSrvIdx_SE, UseSrvIdx = UseSrvIdx,
    ObsFishLenComps = ObsFishLenComps, ISS_FishLenComps = ISS_FishLenComps, UseFishLenComps = UseFishLenComps,
    ObsSrvLenComps = ObsSrvLenComps, ISS_SrvLenComps = ISS_SrvLenComps, UseSrvLenComps = UseSrvLenComps,
    ObsFishAgeComps = ObsFishAgeComps, ISS_FishAgeComps = ISS_FishAgeComps, UseFishAgeComps = UseFishAgeComps,
    ObsSrvAgeComps = ObsSrvAgeComps, ISS_SrvAgeComps = ISS_SrvAgeComps, UseSrvAgeComps = UseSrvAgeComps,
    ObsFish_caal = ObsFish_caal, ISS_Fish_caal = ISS_Fish_caal, UseFish_caal = UseFish_caal,
    ObsSrv_caal = ObsSrv_caal, ISS_Srv_caal = ISS_Srv_caal, UseSrv_caal = UseSrv_caal,
    len_sex = len_sex, age_sex = age_sex,
    AgeingError = AgeingError, ageerror_raw = d$ageerror,
    growth = growth, wtlen = wtlen, mat = mat, rec_dist = rec_dist, rec = rec,
    sel = sel, q = q, var_adj = var_adj,
    comp = list(addtocomp_len = d$len_info$addtocomp[1], addtocomp_age = d$age_info$addtocomp[1],
                minsamplesize = d$len_info$minsamplesize[1], combine_M_F_age = d$age_info$combine_M_F[1]),
    ctl = c, dat_raw = d
  )
}

`%+%` <- function(a, b) paste0(a, b)

if(sys.nframe() == 0) {
  run_dir <- "/Users/matthewcheng/Downloads/goa_rex-main/runs/2025_models/run_14_asfor7_ageing_error"
  dat <- build_rex_data(run_dir)
  out <- file.path("dev/rex_bridge/output/rex_bridge_data.rds")
  saveRDS(dat, out)
  cat("wrote", out, "\n")
  cat("CAAL survey rows used:", sum(dat$UseSrv_caal), " fishery marginal age years:", sum(dat$UseFishAgeComps), "\n")
  cat("length comp years fish/srv:", sum(dat$UseFishLenComps), sum(dat$UseSrvLenComps), "\n")
  str(dat$growth[[1]])
}


#' Attach the SS3 report quantities the bridge test compares against
#'
#' @param dat the list from build_rex_data()
#' @param report the r4ss::SS_output() list for the same run
#' @return dat with a `$ss3` element and an `$mle` element of parameter values
add_rex_ss3_report <- function(dat, report) {

  r <- report
  years <- dat$years; n_yrs <- length(years); n_ages <- length(dat$ages); n_sexes <- dat$n_sexes; n_regions <- dat$n_regions
  p <- r$parameters

  # Parameter values at the MLE -------------------------------------------------
  pv <- function(lab) { v <- p$Value[p$Label == lab]; if(length(v) != 1) stop("parameter not found: ", lab); v }
  growth_mle <- array(NA_real_, dim = c(1, n_regions, n_sexes, 5), dimnames = list(NULL, NULL, c("Fem", "Mal")[1:n_sexes], c("L1", "L2", "K", "CV1", "CV2")))
  for(s in 1:n_sexes) for(gp in 1:n_regions) {
    sx <- c("Fem", "Mal")[s]
    growth_mle[1, gp, s, ] <- c(pv(paste0("L_at_Amin_", sx, "_GP_", gp)), pv(paste0("L_at_Amax_", sx, "_GP_", gp)),
                                pv(paste0("VonBert_K_", sx, "_GP_", gp)), pv(paste0("CV_young_", sx, "_GP_", gp)), pv(paste0("CV_old_", sx, "_GP_", gp)))
  }
  sel_mle <- list()
  for(f in seq_along(dat$fleetnames)) {
    nm <- dat$fleetnames[f]
    rows <- p[grep(paste0("Age_DblN_.*_", nm, "\\(", f, "\\)$"), p$Label), ]
    male <- p[grep(paste0("AgeSel_", f, "Male_.*_", nm, "$"), p$Label), ]
    fnames <- sub(paste0("_", nm, ".*"), "", sub("Age_DblN_", "", rows$Label))
    mnames <- sub(paste0("_", nm, "$"), "", sub(paste0("AgeSel_", f, "Male_"), "", male$Label))
    sel_mle[[nm]] <- list(female = setNames(rows$Value, fnames), female_est = setNames(!is.na(rows$Active_Cnt), fnames),
                          male = setNames(male$Value, mnames), male_est = setNames(!is.na(male$Active_Cnt), mnames))
  }
  rec <- r$recruit
  main <- p[grep("^Main_RecrDev_", p$Label), ]
  # early deviations are labelled by the age they set in the start year: Early_InitAge_n is the year styr - n
  early <- p[grep("^Early_InitAge_|^Early_RecrDev_", p$Label), ]
  early_yr <- ifelse(grepl("InitAge", early$Label), years[1] - as.integer(sub("Early_InitAge_", "", early$Label)), as.integer(sub("Early_RecrDev_", "", early$Label)))
  mle <- list(
    ln_R0 = pv("SR_LN(R0)"),
    rec_dist_area2 = pv("RecrDist_GP_2_area_2_month_1"),
    ln_q = pv(paste0("LnQ_base_", dat$fleetnames[2], "(2)")),
    growth = growth_mle,
    sel = sel_mle,
    main_recdev = setNames(main$Value, as.integer(sub("Main_RecrDev_", "", main$Label))),
    early_recdev = setNames(early$Value, early_yr)[order(early_yr)],
    biasadj = setNames(rec$biasadjuster, rec$Yr),
    Fmort = setNames(r$exploitation$Fishery, r$exploitation$Yr),
    sigmaR = dat$rec$sigmaR
  )

  # Derived quantities ----------------------------------------------------------
  na <- r$natage
  NAA <- array(NA_real_, dim = c(n_regions, n_yrs, n_ages, n_sexes))
  for(a in 1:n_regions) for(s in 1:n_sexes) {
    rows <- na[na$Area == a & na$Bio_Pattern == a & na$Sex == s & na[["Beg/Mid"]] == "B" & na$Era == "TIME", ]
    NAA[a, match(rows$Yr, years), , s] <- as.matrix(rows[, as.character(dat$ages)])
  }
  ts <- r$timeseries
  SSB <- sapply(1:n_regions, function(a) { t <- ts[ts$Area == a & ts$Era == "TIME", ]; t$SpawnBio[match(years, t$Yr)] })
  Rec <- sapply(1:n_regions, function(a) { t <- ts[ts$Area == a & ts$Era == "TIME", ]; t$Recruit_0[match(years, t$Yr)] })
  Bio_all <- sapply(1:n_regions, function(a) { t <- ts[ts$Area == a & ts$Era == "TIME", ]; t$Bio_all[match(years, t$Yr)] })
  dead_B <- { t <- ts[ts$Area == 1 & ts$Era == "TIME", ]; t[["dead(B):_1"]][match(years, t$Yr)] }

  eg <- r$endgrowth
  growth_tab <- lapply(1:n_regions, function(gp) lapply(1:n_sexes, function(s) {
    e <- eg[eg$Bio_Pattern == gp & eg$Sex == s & eg$Platoon == 1, ]
    e[order(e$Age_Beg), c("Age_Beg", "Len_Beg", "SD_Beg", "Len_Mid", "SD_Mid", "Wt_Beg", "Wt_Mid", "Age_Mat")]
  }))

  as <- r$ageselex
  sel_tab <- lapply(seq_along(dat$fleetnames), function(f) lapply(1:n_sexes, function(s) {
    e <- as[as$Fleet == f & as$Factor == "Asel" & as$Sex == s, ]
    e <- e[e$Yr == max(e$Yr[e$Yr <= years[1]]), ]
    as.numeric(e[1, as.character(dat$ages)])
  }))

  list_ss3 <- list(
    likelihoods = r$likelihoods_used, NAA = NAA, SSB = SSB, Rec = Rec, Bio_all = Bio_all, dead_B = dead_B,
    SSB_virgin = r$derived_quants$Value[r$derived_quants$Label == "SSB_Virgin"],
    cpue = r$cpue[, intersect(c("Fleet", "Yr", "Obs", "Exp", "Calc_Q", "SE", "Like"), names(r$cpue))],
    growth = growth_tab, ALK = r$ALK, sel = sel_tab,
    condbase = r$condbase[, intersect(c("Yr", "Fleet", "Sex", "Lbin_lo", "Bin", "Obs", "Exp", "Nsamp_adj", "Nsamp_in", "Like"), names(r$condbase))],
    agedbase = r$agedbase[, intersect(c("Yr", "Fleet", "Sexes", "Sex", "Bin", "Obs", "Exp", "Nsamp_adj", "Nsamp_in", "Like"), names(r$agedbase))],
    lendbase = r$lendbase[, intersect(c("Yr", "Fleet", "Sexes", "Sex", "Bin", "Obs", "Exp", "Nsamp_adj", "Nsamp_in", "Like"), names(r$lendbase))],
    recruit = rec[, intersect(c("Yr", "era", "exp_recr", "bias_adjusted", "pred_recr", "SpawnBio", "dev", "biasadjuster"), names(rec))],
    version = r$SS_version
  )
  dat$ss3 <- list_ss3
  dat$mle <- mle
  dat
}

if(sys.nframe() == 0) {
  rep_file <- "dev/rex_bridge/output/rex_ss3_report.rds"
  if(file.exists(rep_file)) {
    dat <- add_rex_ss3_report(dat, readRDS(rep_file))
    saveRDS(dat, out)
    cat("attached SS3 report; ln_R0", dat$mle$ln_R0, " ln_q", dat$mle$ln_q, " early devs", length(dat$mle$early_recdev), " main devs", length(dat$mle$main_recdev), "\n")
    cat("growth MLE (female, area 1):", dat$mle$growth[1, 1, 1, ], "\n")
    cat("fishery female sel:", dat$mle$sel[[1]]$female, "\n"); cat("fishery male offsets:", dat$mle$sel[[1]]$male, "\n")
  }
}
