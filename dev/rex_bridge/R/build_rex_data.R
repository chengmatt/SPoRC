# Purpose: Build SS3 inputs for GOA rex sole into SPoRC arrays
# Date Created: 8/24/26

rex_paths <- function(run_dir) {
  list(dat = file.path(run_dir, "goa_rex.dat"), ctl = file.path(run_dir, "goa_rex.ctl"),
       starter = file.path(run_dir, "starter.ss"), forecast = file.path(run_dir, "forecast.ss"),
       dir = run_dir)
}

build_rex_data <- function(run_dir) {

  # read the data and control files
  p <- rex_paths(run_dir)
  d <- r4ss::SS_readdat(p$dat, verbose = FALSE)
  ctl <- r4ss::SS_readctl(p$ctl, datlist = d, verbose = FALSE)

  # model years
  years <- d$styr:d$endyr
  n_yrs <- length(years)

  # model ages and the ages compositions are reported at
  ages <- 0:d$Nages
  n_ages <- length(ages)
  obs_ages <- d$agebin_vector
  n_obs_ages <- length(obs_ages)

  # length bins, lower edge and midpoint in cm
  lens_lower <- d$lbin_vector
  n_lens <- length(lens_lower)
  lens_mid <- lens_lower + d$binwidth / 2

  # sexes and areas
  n_sexes <- d$Nsexes
  n_regions <- d$N_areas

  # fleet type 1 is a fishery and type 3 is a survey
  fleetinfo <- d$fleetinfo
  fish_fleets <- which(fleetinfo$type == 1)
  srv_fleets <- which(fleetinfo$type == 3)
  n_fish <- length(fish_fleets)
  n_srv <- length(srv_fleets)
  fleet_area <- fleetinfo$area

  # position of a calendar year in the model years
  yidx <- function(y) match(y, years)

  # Catch ----------------------------------------------------------------------

  # catch is in metric tons, one fishery in area 1
  ObsCatch <- array(NA, dim = c(n_regions, n_yrs, 1, n_fish))
  UseCatch <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))
  catch_se <- array(NA, dim = c(n_regions, n_yrs, 1, n_fish))

  for(i in seq_len(nrow(d$catch))) {

    ci <- d$catch[i, ]

    # the -999 row is equilibrium catch, not a model year
    if(ci$year < d$styr) next

    f <- match(ci$fleet, fish_fleets)
    r <- fleet_area[ci$fleet]

    ObsCatch[r, yidx(ci$year), 1, f] <- ci$catch
    UseCatch[r, yidx(ci$year), 1, f] <- 1
    catch_se[r, yidx(ci$year), 1, f] <- ci$catch_se

  } # end i loop

  init_equil_catch <- d$catch$catch[d$catch$year < d$styr]

  # Survey indices -------------------------------------------------------------

  # indices are survey biomass fit lognormally
  ObsSrvIdx <- array(NA, dim = c(n_regions, n_yrs, 1, n_srv))
  ObsSrvIdx_SE <- array(NA, dim = c(n_regions, n_yrs, 1, n_srv))
  UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))

  cpue <- d$CPUE
  idx_fleet_col <- grep("^index$|^fleet$|FltSvy", names(cpue), value = TRUE)[1]

  for(i in seq_len(nrow(cpue))) {

    # a negative fleet code marks a year SS3 holds out of the fit
    fcode <- cpue[[idx_fleet_col]][i]
    sf <- match(abs(fcode), srv_fleets)
    r <- fleet_area[abs(fcode)]

    ObsSrvIdx[r, yidx(cpue$year[i]), 1, sf] <- cpue$obs[i]
    ObsSrvIdx_SE[r, yidx(cpue$year[i]), 1, sf] <- cpue$se_log[i]
    UseSrvIdx[r, yidx(cpue$year[i]), 1, sf] <- as.numeric(fcode > 0)

  } # end i loop

  # fraction of the year elapsed at the survey
  srv_month <- fleetinfo$surveytiming[srv_fleets]
  t_srv <- (srv_month - 1) / 12

  # Length compositions --------------------------------------------------------

  # rows are proportions at length, a joint sex row runs females then males
  lc <- d$lencomp
  lc_fleet_col <- grep("^fleet$|FltSvy", names(lc), value = TRUE)[1]
  lc_sex_col <- grep("^sex$|Gender", names(lc), value = TRUE)[1]
  lc_n_col <- grep("^Nsamp$", names(lc), value = TRUE)[1]
  lc_vals <- function(row) as.numeric(row[(which(names(lc) == lc_n_col) + 1):ncol(lc)])

  ObsFishLenComps <- array(NA, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_fish))
  ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fish))
  UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))

  ObsSrvLenComps <- array(NA, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_srv))
  ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_srv))
  UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))

  len_sex <- list(fish = rep(NA, n_fish), srv = rep(NA, n_srv))

  for(i in seq_len(nrow(lc))) {

    # a negative fleet code or a negative year is a row SS3 stores but does not fit
    fcode <- lc[[lc_fleet_col]][i]
    fl <- abs(fcode)
    y <- yidx(abs(lc$year[i]))
    r <- fleet_area[fl]
    use_row <- as.numeric(fcode > 0 & lc$year[i] > 0)

    # year outside the model, the fishery carries three 1977 rows
    if(is.na(y)) next

    # sex code 3 is a joint sex row, otherwise only the first block of bins is filled
    v <- lc_vals(lc[i, ])
    sx <- lc[[lc_sex_col]][i]

    if(sx == 3) vals <- matrix(v, nrow = n_lens, ncol = n_sexes)
    else vals <- matrix(v[1:n_lens], nrow = n_lens, ncol = 1)

    if(fl %in% fish_fleets) {

      f <- match(fl, fish_fleets)
      len_sex$fish[f] <- sx

      if(sx == 3) ObsFishLenComps[r, y, 1, , , f] <- vals
      else ObsFishLenComps[r, y, 1, , , f] <- cbind(vals, 0)

      ISS_FishLenComps[r, y, 1, , f] <- lc[[lc_n_col]][i]
      UseFishLenComps[r, y, 1, f] <- use_row

    } else {

      sf <- match(fl, srv_fleets)
      len_sex$srv[sf] <- sx

      if(sx == 3) ObsSrvLenComps[r, y, 1, , , sf] <- vals
      else ObsSrvLenComps[r, y, 1, , , sf] <- cbind(vals, 0)

      ISS_SrvLenComps[r, y, 1, , sf] <- lc[[lc_n_col]][i]
      UseSrvLenComps[r, y, 1, sf] <- use_row

    } # end if fishery or survey fleet

  } # end i loop

  # Age compositions -----------------------------------------------------------

  # the fishery gives marginal ages (joint sex, Lbin -1) and the surveys give conditional age at
  # length (one sex and one length bin per row, otolith counts). ghost fleets are stored, never fit
  ac <- d$agecomp
  ac_fleet_col <- grep("^fleet$|FltSvy", names(ac), value = TRUE)[1]
  ac_sex_col <- grep("^sex$|Gender", names(ac), value = TRUE)[1]
  ac_n_col <- grep("^Nsamp$", names(ac), value = TRUE)[1]
  ac_vals <- function(row) as.numeric(row[(which(names(ac) == ac_n_col) + 1):ncol(ac)])

  ObsFishAgeComps <- array(NA, dim = c(n_regions, n_yrs, 1, n_obs_ages, n_sexes, n_fish))
  ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fish))
  UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_fish))

  ObsSrvAgeComps <- array(NA, dim = c(n_regions, n_yrs, 1, n_obs_ages, n_sexes, n_srv))
  ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_srv))
  UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, 1, n_srv))

  ObsFish_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_obs_ages, n_sexes, n_fish))
  ISS_Fish_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_fish))
  UseFish_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_fish))

  ObsSrv_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_obs_ages, n_sexes, n_srv))
  ISS_Srv_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_sexes, n_srv))
  UseSrv_caal <- array(0, dim = c(n_regions, n_yrs, 1, n_lens, n_srv))

  age_sex <- list(fish = rep(NA, n_fish), srv = rep(NA, n_srv))

  for(i in seq_len(nrow(ac))) {

    # a negative fleet code or a negative year is a row SS3 stores but does not fit
    fcode <- ac[[ac_fleet_col]][i]
    fl <- abs(fcode)
    y <- yidx(abs(ac$year[i]))
    r <- fleet_area[fl]
    use <- as.numeric(fcode > 0 & ac$year[i] > 0)

    v <- ac_vals(ac[i, ])
    sx <- ac[[ac_sex_col]][i]
    nn <- ac[[ac_n_col]][i]
    lo <- ac$Lbin_lo[i]

    # marginal ages, pooled across all length bins
    if(lo < 0) {

      if(sx == 3) vals <- matrix(v, nrow = n_obs_ages, ncol = n_sexes)
      else vals <- matrix(v[1:n_obs_ages], nrow = n_obs_ages, ncol = 1)

      if(fl %in% fish_fleets) {

        f <- match(fl, fish_fleets)
        age_sex$fish[f] <- sx

        if(sx == 3) ObsFishAgeComps[r, y, 1, , , f] <- vals
        else ObsFishAgeComps[r, y, 1, , , f] <- cbind(vals, 0)

        ISS_FishAgeComps[r, y, 1, , f] <- nn
        UseFishAgeComps[r, y, 1, f] <- use

      } else {

        sf <- match(fl, srv_fleets)
        age_sex$srv[sf] <- sx

        if(sx == 3) ObsSrvAgeComps[r, y, 1, , , sf] <- vals
        else ObsSrvAgeComps[r, y, 1, , , sf] <- cbind(vals, 0)

        ISS_SrvAgeComps[r, y, 1, , sf] <- nn
        UseSrvAgeComps[r, y, 1, sf] <- use

      } # end if fishery or survey fleet

    # conditional age at length, one length bin (Lbin_lo equals Lbin_hi, in cm) and one sex
    } else {

      if(ac$Lbin_hi[i] != lo) stop("CAAL row spans more than one length bin; not handled")

      l <- match(lo, lens_lower)
      if(is.na(l)) stop("CAAL Lbin_lo not a length bin edge: ", lo)

      # sex code 1 is female and 2 is male, otolith counts sit in that sex's block of the row
      if(sx == 1) counts <- v[1:n_obs_ages]
      else if(sx == 2) counts <- v[n_obs_ages + 1:n_obs_ages]
      else stop("joint-sex CAAL row not handled")

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

      } # end if fishery or survey fleet

    } # end if marginal ages or conditional age at length

  } # end i loop

  # Ageing error ---------------------------------------------------------------

  # SS3 stores the mean and standard deviation of the observed age at each true age 0 to 20, and a
  # mean of -1 means unbiased, read as true age plus 0.5
  ae_sd <- as.numeric(d$ageerror[2, ])
  ae_mean <- as.numeric(d$ageerror[1, ])
  ae_mean[ae_mean < 0] <- ages[ae_mean < 0] + 0.5

  # observed bin b spans [age_b, age_b + 1), with the lower tail in the first bin and the upper in the last
  AgeingError <- matrix(0, n_ages, n_obs_ages)
  edges <- c(-Inf, obs_ages[-1], Inf)

  for(a in seq_len(n_ages)) {

    cdf <- stats::pnorm(edges, ae_mean[a], ae_sd[a])
    AgeingError[a, ] <- diff(cdf)

  } # end a loop

  # Biology and growth ---------------------------------------------------------

  # growth and natural mortality are read per sex and growth pattern, lengths in cm and M per year
  mg <- ctl$MG_parms
  growth <- list()

  for(s in 1:n_sexes) {
    for(gp in 1:n_regions) {

      sx <- c("Fem", "Mal")[s]
      rows <- mg[grep(paste0("_", sx, "_GP_", gp, "$"), rownames(mg)), ]

      growth[[paste0("s", s, "_gp", gp)]] <- list(
        M = rows[paste0("NatM_p_1_", sx, "_GP_", gp), "INIT"],
        L1 = rows[grep("L_at_Amin", rownames(rows)), "INIT"],
        L2 = rows[grep("L_at_Amax", rownames(rows)), "INIT"],
        K = rows[grep("VonBert_K", rownames(rows)), "INIT"],
        CV1 = rows[grep("CV_young", rownames(rows)), "INIT"],
        CV2 = rows[grep("CV_old", rownames(rows)), "INIT"],
        est = rows[grep("L_at_Amin|L_at_Amax|VonBert_K|CV_young|CV_old", rownames(rows)), "PHASE"] > 0)

    } # end gp loop
  } # end s loop

  # weight length coefficients, weight in kg and length in cm
  wtlen <- list(
    fem = mg[grep("^Wtlen_[12]_Fem", rownames(mg)), "INIT"][1:2],
    mal = mg[grep("^Wtlen_[12]_Mal", rownames(mg)), "INIT"][1:2]
  )

  # maturity is read from the first growth pattern, the assessment gives every pattern the same curve
  mat <- list(a50 = mg[grep("^Mat50%_Fem", rownames(mg))[1], "INIT"],
              slope = mg[grep("^Mat_slope_Fem", rownames(mg))[1], "INIT"],
              first_mature_age = ctl$First_Mature_Age)

  rec_dist <- mg[grep("^RecrDist", rownames(mg)), c("INIT", "PHASE")]

  # stock recruit relationship and the recruitment deviation year blocks
  sr <- ctl$SR_parms
  rec <- list(ln_R0 = sr["SR_LN(R0)", "INIT"], h = sr["SR_BH_steep", "INIT"], sigmaR = sr["SR_sigmaR", "INIT"],
              main_first = ctl$MainRdevYrFirst, main_last = ctl$MainRdevYrLast,
              early_start = ctl$recdev_early_start, early_phase = ctl$recdev_early_phase,
              bias_years = c(ctl$last_early_yr_nobias_adj, ctl$first_yr_fullbias_adj, ctl$last_yr_fullbias_adj, ctl$first_recent_yr_nobias_adj),
              max_bias_adj = ctl$max_bias_adj)

  # selectivity, catchability, and the francis variance adjustments
  sel <- list(age = ctl$age_selex_parms, age_types = ctl$age_selex_types, size_types = ctl$size_selex_types)
  q <- list(options = ctl$Q_options, parms = ctl$Q_parms)
  var_adj <- ctl$Variance_adjustment_list

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
    ctl = ctl, dat_raw = d
  )

} # end function

# attach the SS3 report quantities the bridge compares against
add_rex_ss3_report <- function(dat, report) {

  r <- report
  p <- r$parameters

  years <- dat$years
  n_yrs <- length(years)
  n_ages <- length(dat$ages)
  n_sexes <- dat$n_sexes
  n_regions <- dat$n_regions

  # Parameter values at the maximum likelihood estimate -------------------------

  pv <- function(lab) {
    v <- p$Value[p$Label == lab]
    if(length(v) != 1) stop("parameter not found: ", lab)
    v
  }

  growth_mle <- array(NA, dim = c(1, n_regions, n_sexes, 5),
                      dimnames = list(NULL, NULL, c("Fem", "Mal")[1:n_sexes], c("L1", "L2", "K", "CV1", "CV2")))

  for(s in 1:n_sexes) {
    for(gp in 1:n_regions) {

      sx <- c("Fem", "Mal")[s]

      growth_mle[1, gp, s, ] <- c(pv(paste0("L_at_Amin_", sx, "_GP_", gp)),
                                  pv(paste0("L_at_Amax_", sx, "_GP_", gp)),
                                  pv(paste0("VonBert_K_", sx, "_GP_", gp)),
                                  pv(paste0("CV_young_", sx, "_GP_", gp)),
                                  pv(paste0("CV_old_", sx, "_GP_", gp)))

    } # end gp loop
  } # end s loop

  # female selectivity parameters and the male offsets from them
  sel_mle <- list()

  for(f in seq_along(dat$fleetnames)) {

    nm <- dat$fleetnames[f]
    rows <- p[grep(paste0("Age_DblN_.*_", nm, "\\(", f, "\\)$"), p$Label), ]
    male <- p[grep(paste0("AgeSel_", f, "Male_.*_", nm, "$"), p$Label), ]

    fnames <- sub(paste0("_", nm, ".*"), "", sub("Age_DblN_", "", rows$Label))
    mnames <- sub(paste0("_", nm, "$"), "", sub(paste0("AgeSel_", f, "Male_"), "", male$Label))

    sel_mle[[nm]] <- list(female = stats::setNames(rows$Value, fnames),
                          female_est = stats::setNames(!is.na(rows$Active_Cnt), fnames),
                          male = stats::setNames(male$Value, mnames),
                          male_est = stats::setNames(!is.na(male$Active_Cnt), mnames))

  } # end f loop

  # early deviations are labeled by the age they set in the start year, so Early_InitAge_n is year styr minus n
  rec <- r$recruit
  main <- p[grep("^Main_RecrDev_", p$Label), ]
  early <- p[grep("^Early_InitAge_|^Early_RecrDev_", p$Label), ]
  early_yr <- ifelse(grepl("InitAge", early$Label), years[1] - as.integer(sub("Early_InitAge_", "", early$Label)),
                     as.integer(sub("Early_RecrDev_", "", early$Label)))

  mle <- list(
    ln_R0 = pv("SR_LN(R0)"),
    rec_dist_area2 = pv("RecrDist_GP_2_area_2_month_1"),
    ln_q = pv(paste0("LnQ_base_", dat$fleetnames[2], "(2)")),
    growth = growth_mle,
    sel = sel_mle,
    main_recdev = stats::setNames(main$Value, as.integer(sub("Main_RecrDev_", "", main$Label))),
    early_recdev = stats::setNames(early$Value, early_yr)[order(early_yr)],
    biasadj = stats::setNames(rec$biasadjuster, rec$Yr),
    Fmort = stats::setNames(r$exploitation$Fishery, r$exploitation$Yr),
    sigmaR = dat$rec$sigmaR
  )

  # Derived quantities ----------------------------------------------------------

  # numbers at age at the start of the year
  na <- r$natage
  NAA <- array(NA, dim = c(n_regions, n_yrs, n_ages, n_sexes))

  for(a in 1:n_regions) {
    for(s in 1:n_sexes) {

      rows <- na[na$Area == a & na$Bio_Pattern == a & na$Sex == s & na[["Beg/Mid"]] == "B" & na$Era == "TIME", ]
      NAA[a, match(rows$Yr, years), , s] <- as.matrix(rows[, as.character(dat$ages)])

    } # end s loop
  } # end a loop

  # spawning biomass, age 0 recruits, and total biomass by area
  ts <- r$timeseries
  SSB <- Rec <- Bio_all <- matrix(NA, n_yrs, n_regions)

  for(a in 1:n_regions) {

    tsa <- ts[ts$Area == a & ts$Era == "TIME", ]
    yr <- match(years, tsa$Yr)

    SSB[, a] <- tsa$SpawnBio[yr]
    Rec[, a] <- tsa$Recruit_0[yr]
    Bio_all[, a] <- tsa$Bio_all[yr]

  } # end a loop

  # dead catch in biomass, the fishery sits in area 1
  ts1 <- ts[ts$Area == 1 & ts$Era == "TIME", ]
  dead_B <- ts1[["dead(B):_1"]][match(years, ts1$Yr)]

  # length, weight, and maturity at age by growth pattern and sex
  eg <- r$endgrowth
  growth_tab <- vector("list", n_regions)

  for(gp in 1:n_regions) {

    growth_tab[[gp]] <- vector("list", n_sexes)

    for(s in 1:n_sexes) {

      e <- eg[eg$Bio_Pattern == gp & eg$Sex == s & eg$Platoon == 1, ]
      growth_tab[[gp]][[s]] <- e[order(e$Age_Beg), c("Age_Beg", "Len_Beg", "SD_Beg", "Len_Mid", "SD_Mid", "Wt_Beg", "Wt_Mid", "Age_Mat")]

    } # end s loop
  } # end gp loop

  # selectivity at age in the first model year, selectivity is time invariant here
  agesel <- r$ageselex
  sel_tab <- vector("list", length(dat$fleetnames))

  for(f in seq_along(dat$fleetnames)) {

    sel_tab[[f]] <- vector("list", n_sexes)

    for(s in 1:n_sexes) {

      e <- agesel[agesel$Fleet == f & agesel$Factor == "Asel" & agesel$Sex == s, ]
      e <- e[e$Yr == max(e$Yr[e$Yr <= years[1]]), ]
      sel_tab[[f]][[s]] <- as.numeric(e[1, as.character(dat$ages)])

    } # end s loop
  } # end f loop

  dat$ss3 <- list(
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

  dat$mle <- mle
  dat

} # end function
