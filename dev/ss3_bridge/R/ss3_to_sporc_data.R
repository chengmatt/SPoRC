# Purpose: Parse generally ... any SS3 run directory into SPoRC-shaped arrays
# Date Created: 9/2/26
#
#

# Files ------------------------------------------------------------------------

#' Resolve the SS3 files in a run directory
#'
#' starter.ss names the data and control files, so the run directory is the only
#' argument a caller supplies. wtatage.ss is optional and only read when the
#' control file asks for it.
ss3_paths <- function(run_dir) {

  starter_path <- file.path(run_dir, "starter.ss")
  if(!file.exists(starter_path)) stop("no starter.ss in ", run_dir)
  starter <- r4ss::SS_readstarter(starter_path, verbose = FALSE)

  # SS3 accepts wtatage.ss or wtatage.ss_new, preferring whichever is present
  wtatage <- file.path(run_dir, c("wtatage.ss", "wtatage.ss_new"))
  wtatage <- wtatage[file.exists(wtatage)][1]

  list(
    dir = run_dir,
    starter = starter_path,
    dat = file.path(run_dir, starter$datfile),
    ctl = file.path(run_dir, starter$ctlfile),
    forecast = file.path(run_dir, "forecast.ss"),
    wtatage = wtatage,
    report = file.path(run_dir, "Report.sso")
  )

} # end function

#' Bin midpoints from a vector of lower edges
#'
#' The last bin has no upper edge in the file, so it reuses the width of the one
#' before it.
bin_midpoints <- function(lower) {

  if(length(lower) == 0) return(numeric(0))
  if(length(lower) == 1) return(lower)

  w <- diff(lower)
  lower + c(w, w[length(w)]) / 2

} # end function

# Time blocks ------------------------------------------------------------------

#' Expand one SS3 block pattern into model year indices
#'
#' Block_Design holds begin and end years in pairs. Returns one integer vector of
#' year positions per block, clipped to the model years, so a terminal block that
#' runs past endyr for the forecast contributes only its modeled part.
ss3_block_years <- function(ctl, pattern, years) {

  design <- ctl$Block_Design[[pattern]]
  begins <- design[seq(1, length(design), by = 2)]
  ends <- design[seq(2, length(design), by = 2)]

  lapply(seq_along(begins), function(b) which(years >= begins[b] & years <= ends[b]))

} # end function

#' Apply a parameter's time blocks to give one value per model year
#'
#' base is the parameter's INIT, tv is the vector of block values in file order,
#' and fxn is the SS3 block function code. Years outside every block keep base.
ss3_apply_blocks <- function(base, tv, fxn, blk_years, n_yrs) {

  out <- rep(base, n_yrs)

  for(b in seq_along(blk_years)) {

    yi <- blk_years[[b]]
    if(length(yi) == 0) next

    # 0 multiplies the base, 1 adds to it, 2 replaces it, 3 accumulates from the
    # value the previous block left behind
    if(fxn == 0) out[yi] <- base * exp(tv[b])
    if(fxn == 1) out[yi] <- base + tv[b]
    if(fxn == 2) out[yi] <- tv[b]
    if(fxn == 3) out[yi] <- out[min(yi) - 1] + tv[b]

  } # end b loop

  out

} # end function

#' Pull the block values SS3 stored for one base parameter
#'
#' Time-varying rows are named <base>_BLK<pattern><fxn tag>_<block start year> and
#' are written in block order, so file order is the block order.
ss3_tv_rows <- function(tv_table, base_name, pattern) {

  if(is.null(tv_table) || nrow(tv_table) == 0) return(numeric(0))
  hits <- grep(paste0("^", base_name, "_BLK", pattern), rownames(tv_table))
  tv_table[hits, "INIT"]

} # end function

# Natural mortality ------------------------------------------------------------

#' Build fixed natural mortality at [pop, region, year, age, sex]
#'
#' Covers natM_type 0 (one rate), 1 (breakpoints interpolated over age) and 3
#' (one rate per age), each optionally with time blocks. Growth patterns map
#' onto regions, which holds whenever an SS3 area owns one pattern.
ss3_natmort <- function(ctl, ages, years, n_pop, n_regions, n_sexes) {

  if(!ctl$natM_type %in% c(0, 1, 3)) {
    stop("natM_type ", ctl$natM_type, " is not handled; only 0, 1 and 3 are")
  }

  n_ages <- length(ages)
  n_yrs <- length(years)
  mg <- ctl$MG_parms
  M <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes))

  # number of M parameters per sex and growth pattern
  n_M <- switch(as.character(ctl$natM_type), "0" = 1, "1" = ctl$N_natMparms, "3" = n_ages)

  for(s in seq_len(n_sexes)) {
    for(gp in seq_len(n_regions)) {

      sx <- c("Fem", "Mal")[s]

      # one column of per-year values per M parameter, blocks applied first
      by_parm <- matrix(NA, n_yrs, n_M)

      for(k in seq_len(n_M)) {

        nm <- paste0("NatM_p_", k, "_", sx, "_GP_", gp)
        if(!nm %in% rownames(mg)) stop("control file has no parameter ", nm)

        base <- mg[nm, "INIT"]
        pattern <- mg[nm, "Block"]

        if(pattern > 0) {
          tv <- ss3_tv_rows(ctl$MG_parms_tv, nm, pattern)
          by_parm[, k] <- ss3_apply_blocks(base, tv, mg[nm, "Block_Fxn"],
                                           ss3_block_years(ctl, pattern, years), n_yrs)
        } else {
          by_parm[, k] <- base
        } # end if else

      } # end k loop

      # spread the parameters over ages: type 3 is already one per age, type 1
      # interpolates between the breakpoint ages and holds the end values flat
      for(y in seq_len(n_yrs)) {

        if(ctl$natM_type == 0) M_age <- rep(by_parm[y, 1], n_ages)
        if(ctl$natM_type == 3) M_age <- by_parm[y, ]
        if(ctl$natM_type == 1) {
          M_age <- stats::approx(ctl$M_ageBreakPoints, by_parm[y, ], xout = ages,
                                 rule = 2)$y
        } # end if

        M[, gp, y, , s] <- rep(M_age, each = n_pop)

      } # end y loop

    } # end gp loop
  } # end s loop

  M

} # end function

# Ageing error -----------------------------------------------------------------

#' Convert SS3 ageing error definitions into [true age, observed age] matrices
#'
#' SS3 stores the mean and standard deviation of observed age at each true age; a
#' mean of -1 means unbiased and is read as true age plus a half year. Observed
#' bin b spans [age_b, age_b + 1), with the tails accumulated into the end bins.
ss3_ageing_error <- function(dat, ages, obs_ages) {

  n_defs <- dat$N_ageerror_definitions
  out <- array(0, dim = c(n_defs, length(ages), length(obs_ages)))
  edges <- c(-Inf, obs_ages[-1], Inf)

  for(k in seq_len(n_defs)) {

    rows <- dat$ageerror[(2 * k - 1):(2 * k), ]
    ae_mean <- as.numeric(rows[1, ])
    ae_sd <- as.numeric(rows[2, ])
    ae_mean[ae_mean < 0] <- ages[ae_mean < 0] + 0.5

    for(a in seq_along(ages)) out[k, a, ] <- diff(stats::pnorm(edges, ae_mean[a], ae_sd[a]))

  } # end k loop

  out

} # end function

# Growth ------------------------------------------------------------------------

#' Length at age from the SS3 von Bertalanffy parameters
#'
#' L1 is the length at age A1 and L2 the length at age A2. When A2 is the 999
#' the L2 flag value is already Linf; otherwise Linf is solved from the two anchors.
#' Growth is linear below A1, which SS3 holds at the L1 value.
ss3_laa <- function(ages, L1, L2, K, A1, A2) {

  Linf <- if(A2 == 999) L2 else L1 + (L2 - L1) / (1 - exp(-K * (A2 - A1)))
  LAA <- Linf + (L1 - Linf) * exp(-K * (ages - A1))
  LAA[ages < A1] <- L1
  LAA

} # end function

#' Standard deviation of length at age
#'
#' CV_Growth_Pattern selects whether the two parameters are coefficients of
#' variation or standard deviations, and whether they are interpolated over
#' length at age or over age. Interpolation is linear and kept flat outside the
#' A1 to A2 anchors.
ss3_sd_laa <- function(ages, LAA, CV1, CV2, A1, A2, pattern) {

  if(!pattern %in% 0:3) stop("CV_Growth_Pattern ", pattern, " is not handled")

  # anchors are lengths for the even patterns and ages for the odd ones
  A2_use <- if(A2 == 999) max(ages) else A2
  x <- if(pattern %in% c(0, 2)) LAA else ages
  x1 <- if(pattern %in% c(0, 2)) LAA[which.min(abs(ages - A1))] else A1
  x2 <- if(pattern %in% c(0, 2)) LAA[which.min(abs(ages - A2_use))] else A2_use

  frac <- if(x2 == x1) rep(0, length(x)) else pmin(1, pmax(0, (x - x1) / (x2 - x1)))
  v <- CV1 + (CV2 - CV1) * frac

  # patterns 0 and 1 give a coefficient of variation, 2 and 3 give the deviation
  if(pattern %in% c(0, 1)) v * LAA else v

} # end function

#' Maturity at age from the SS3 maturity options
#'
#' Option 2 is logistic in age and reads straight off. Option 1 is logistic in
#' length, so it is integrated over the normal distribution of length at age.
#' Options 3 and 4 read a matrix this parser does not have.
ss3_maturity <- function(ctl, mg, sx, ages, LAA, sd_LAA, lens_lower) {

  if(ctl$maturity_option == 2) {

    a50 <- mg[paste0("Mat50%_", sx, "_GP_1"), "INIT"]
    slope <- mg[paste0("Mat_slope_", sx, "_GP_1"), "INIT"]
    mat <- 1 / (1 + exp(slope * (ages - a50)))

  } else if(ctl$maturity_option == 1) {

    l50 <- mg[paste0("Mat50%_", sx, "_GP_1"), "INIT"]
    slope <- mg[paste0("Mat_slope_", sx, "_GP_1"), "INIT"]
    if(length(lens_lower) == 0) stop("maturity_option 1 needs population length bins")

    # weight maturity at length by the share of each age sitting in that bin
    mids <- bin_midpoints(lens_lower)
    edges <- c(-Inf, lens_lower[-1], Inf)
    mat_len <- 1 / (1 + exp(slope * (mids - l50)))

    mat <- vapply(seq_along(ages), function(a) {
      pl <- diff(stats::pnorm(edges, LAA[a], sd_LAA[a]))
      sum(pl * mat_len)
    }, numeric(1))

  } else {

    stop("maturity_option ", ctl$maturity_option, " is not handled; only 1 and 2 are")

  } # end if else

  mat[ages < ctl$First_Mature_Age] <- 0
  mat

} # end function

# Weight and maturity at age ---------------------------------------------------

#' Weight at age and maturity at age on the SPoRC axes
#'
#' When the control file sets EmpiricalWAA, SS3 takes both from wtatage.ss:
#' fleet 0 is the population weight at the start of the season, fleet -1 the
#' same at mid season, fleet -2 maturity times fecundity, and fleets 1..N their
#' own weights. Otherwise both are rebuilt from the growth, weight-length and
#' maturity parameters in the control file.
#'
#' Fleet -2 is already a weight, so it is NOT SPoRC's MatAA. SPoRC forms
#' spawning biomass as \code{NAA * WAA * MatAA}, which would multiply weight in
#' twice. Dividing fecundity by the population weight recovers the maturity
#' ogive that makes the product right.
ss3_waa <- function(dat, ctl, wtatage_path, ages, years, n_pop, n_regions, n_seas,
                    n_sexes, n_fish, n_srv, fish_fleets, srv_fleets, lens_lower,
                    waa_fallback = FALSE) {

  n_ages <- length(ages)
  n_yrs <- length(years)
  dims <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)

  WAA <- array(NA, dim = dims)
  MatAA <- array(NA, dim = dims)
  WAA_fish <- array(NA, dim = c(dims, n_fish))
  WAA_srv <- array(NA, dim = c(dims, n_srv))

  # a missing wtatage.ss is a different model, not a missing detail, so it stops
  # unless the caller has said in so many words to rebuild weights from growth
  empirical <- isTRUE(ctl$EmpiricalWAA == 1) && !is.na(wtatage_path)

  if(isTRUE(ctl$EmpiricalWAA == 1) && is.na(wtatage_path)) {
    if(!waa_fallback) {
      stop("the control file sets EmpiricalWAA but no wtatage.ss is in the run directory; ",
           "pass waa_fallback = TRUE to rebuild weight and maturity at age from the ",
           "growth parameters instead, which will not reproduce the assessment")
    }
    warning("no wtatage.ss found; weight and maturity at age rebuilt from growth parameters")
  } # end if

  if(empirical) {

    w <- r4ss::SS_readwtatage(wtatage_path, verbose = FALSE)
    age_cols <- as.character(ages)

    # a row's year is negative when SS3 should hold it constant from that year on
    fill <- function(fleet_code, sex) {

      sub <- w[w$fleet == fleet_code & w$sex == sex, , drop = FALSE]
      if(nrow(sub) == 0) return(NULL)

      vals <- matrix(NA, n_yrs, n_ages)

      for(i in seq_len(nrow(sub))) {

        yr <- abs(sub$year[i])
        yi <- if(sub$year[i] < 0) which(years >= yr) else match(yr, years)
        yi <- yi[!is.na(yi)]
        if(length(yi) == 0) next
        vals[yi, ] <- rep(as.numeric(sub[i, age_cols]), each = length(yi))

      } # end i loop

      vals

    } # end function

    for(s in seq_len(n_sexes)) {

      pop <- fill(0, s)
      fec <- fill(-2, s)
      if(is.null(pop) || is.null(fec)) stop("wtatage.ss is missing fleet 0 or -2 for sex ", s)

      # the ogive that returns fecundity when SPoRC multiplies it back by weight;
      # an age with no weight has no spawning output either
      mat <- ifelse(pop > 0, fec / pop, 0)

      for(y in seq_len(n_yrs)) {
        WAA[, , y, , , s] <- pop[y, ]
        MatAA[, , y, , , s] <- mat[y, ]
      } # end y loop

      # a fleet with no row of its own inherits the population weights
      for(f in seq_along(fish_fleets)) {
        v <- fill(fish_fleets[f], s)
        if(is.null(v)) v <- pop
        for(y in seq_len(n_yrs)) WAA_fish[, , y, , , s, f] <- v[y, ]
      } # end f loop

      for(sf in seq_along(srv_fleets)) {
        v <- fill(srv_fleets[sf], s)
        if(is.null(v)) v <- pop
        for(y in seq_len(n_yrs)) WAA_srv[, , y, , , s, sf] <- v[y, ]
      } # end sf loop

    } # end s loop

    # maturity was recovered from fecundity rather than read as an ogive
    mat_is_fecundity <- TRUE

  } else {

    mg <- ctl$MG_parms

    for(s in seq_len(n_sexes)) {

      sx <- c("Fem", "Mal")[s]
      A1 <- ctl$Growth_Age_for_L1
      A2 <- ctl$Growth_Age_for_L2
      L1 <- mg[paste0("L_at_Amin_", sx, "_GP_1"), "INIT"]
      L2 <- mg[paste0("L_at_Amax_", sx, "_GP_1"), "INIT"]
      K <- mg[paste0("VonBert_K_", sx, "_GP_1"), "INIT"]
      CV1 <- mg[paste0("CV_young_", sx, "_GP_1"), "INIT"]
      CV2 <- mg[paste0("CV_old_", sx, "_GP_1"), "INIT"]
      wa <- mg[paste0("Wtlen_1_", sx, "_GP_1"), "INIT"]
      wb <- mg[paste0("Wtlen_2_", sx, "_GP_1"), "INIT"]

      LAA <- ss3_laa(ages, L1, L2, K, A1, A2)
      sd_LAA <- ss3_sd_laa(ages, LAA, CV1, CV2, A1, A2, ctl$CV_Growth_Pattern)
      W <- wa * LAA^wb
      mat <- ss3_maturity(ctl, mg, sx, ages, LAA, sd_LAA, lens_lower)

      WAA[, , , , , s] <- rep(W, each = n_pop * n_regions * n_yrs * n_seas)
      MatAA[, , , , , s] <- rep(mat, each = n_pop * n_regions * n_yrs * n_seas)
      for(f in seq_len(n_fish)) WAA_fish[, , , , , s, f] <- WAA[, , , , , s]
      for(sf in seq_len(n_srv)) WAA_srv[, , , , , s, sf] <- WAA[, , , , , s]

    } # end s loop

    # here maturity is a proportion, so spawning biomass is maturity times weight
    mat_is_fecundity <- FALSE

  } # end if else

  list(
    WAA = WAA,
    MatAA = MatAA,
    WAA_fish = WAA_fish,
    WAA_srv = WAA_srv,
    mat_is_fecundity = mat_is_fecundity
  )

} # end function

# Main entry point -------------------------------------------------------------

#' Parse an SS3 run directory into SPoRC arrays
#'
#' @param run_dir Directory holding starter.ss and the files it names.
#' @param waa_fallback Rebuild weight and maturity at age from the growth
#'   parameters when the control file asks for a wtatage.ss that is not present.
#'   Off by default because the result is a different model.
#' @return A list of SPoRC-shaped observation arrays and the parsed SS3 biology,
#'   selectivity, catchability and recruitment tables a bridge needs.
ss3_to_sporc_data <- function(run_dir, waa_fallback = FALSE) {

  p <- ss3_paths(run_dir)
  d <- r4ss::SS_readdat(p$dat, verbose = FALSE)
  ctl <- r4ss::SS_readctl(p$ctl, datlist = d, verbose = FALSE)

  # Dimensions -----------------------------------------------------------------

  years <- d$styr:d$endyr
  n_yrs <- length(years)

  # SS3 counts ages from 0 and treats Nages as the accumulator age
  ages <- 0:d$Nages
  n_ages <- length(ages)
  obs_ages <- d$agebin_vector
  n_obs_ages <- length(obs_ages)

  n_seas <- d$nseas
  n_sexes <- d$Nsexes
  n_regions <- d$N_areas
  n_pop <- 1

  # SS3 keeps two length axes: lbin_vector is what compositions are reported on and what SPoRC
  # models, lbin_vector_pop the finer internal grid. bins are given as lower edges
  lens_lower <- if(length(d$lbin_vector) > 0) as.numeric(d$lbin_vector) else numeric(0)
  n_lens <- length(lens_lower)

  # midpoints come from the bin edges rather than binwidth, which a model reading its population
  # bins as a vector does not have; a NULL binwidth would collapse the length axis silently
  lens_mid <- bin_midpoints(lens_lower)
  pop_lens_lower <- if(length(d$lbin_vector_pop) > 0) d$lbin_vector_pop else lens_lower

  # type 1 is a fishery, 2 is bycatch only, 3 is a survey and 4 is ignored
  fleetinfo <- d$fleetinfo
  fish_fleets <- which(fleetinfo$type %in% c(1, 2))
  srv_fleets <- which(fleetinfo$type == 3)
  n_fish <- length(fish_fleets)
  n_srv <- length(srv_fleets)
  fleet_area <- fleetinfo$area

  # positions along the model axes; SS3 addresses these by value, SPoRC by index
  yidx <- function(y) match(y, years)
  sidx <- function(s) pmax(1, pmin(n_seas, s))

  # in a two-sex model every composition row has one block of bins per sex: code 3 fills both,
  # 1 and 2 fill their own, and 0 is a combined row in the first block
  comp_blocks <- function(v, n_bins, code) {

    out <- matrix(0, n_bins, n_sexes)
    if(n_sexes == 1) { out[, 1] <- v[seq_len(n_bins)]; return(out) }

    if(code == 3) {
      out[, 1] <- v[seq_len(n_bins)]
      out[, 2] <- v[n_bins + seq_len(n_bins)]
    } else {
      slot <- if(code == 2) 2 else 1
      out[, slot] <- v[(slot - 1) * n_bins + seq_len(n_bins)]
    } # end if else

    out

  } # end function

  # which sex blocks a row actually has data for
  sex_slots <- function(code) if(code == 3) seq_len(n_sexes) else if(code == 2) 2 else 1

  # SS3 holds a row out of the likelihood by negating either its year or its fleet
  is_fit <- function(year, fleet) year > 0 & fleet > 0

  # Catch ----------------------------------------------------------------------

  ObsCatch <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish))
  UseCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))
  catch_se <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish))

  for(i in seq_len(nrow(d$catch))) {

    ci <- d$catch[i, ]

    # rows before styr hold the equilibrium catch, not a model year
    if(ci$year < d$styr) next

    f <- match(ci$fleet, fish_fleets)
    if(is.na(f)) next
    r <- fleet_area[ci$fleet]

    ObsCatch[r, yidx(ci$year), sidx(ci$seas), f] <- ci$catch
    UseCatch[r, yidx(ci$year), sidx(ci$seas), f] <- as.numeric(ci$catch > 0)
    catch_se[r, yidx(ci$year), sidx(ci$seas), f] <- ci$catch_se

  } # end i loop

  eq <- d$catch[d$catch$year < d$styr, , drop = FALSE]
  init_equil_catch <- stats::setNames(eq$catch, eq$fleet)

  # units of catch are 1 for biomass and 2 for numbers
  catch_units <- fleetinfo$units[fish_fleets]

  # Indices --------------------------------------------------------------------

  ObsSrvIdx <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_srv))
  ObsSrvIdx_SE <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_srv))
  UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))

  ObsFishIdx <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish))
  ObsFishIdx_SE <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish))
  UseFishIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))

  cpue <- d$CPUE
  idx_col <- grep("^index$|^fleet$|FltSvy", names(cpue), value = TRUE)[1]

  for(i in seq_len(nrow(cpue))) {

    fl <- abs(cpue[[idx_col]][i])
    yr <- abs(cpue$year[i])
    yi <- yidx(yr)
    if(is.na(yi)) next

    r <- fleet_area[fl]
    se <- sidx(match(cpue$month[i], sort(unique(cpue$month))))
    use <- as.numeric(is_fit(cpue$year[i], cpue[[idx_col]][i]))
    if(n_seas == 1) se <- 1

    if(fl %in% srv_fleets) {
      sf <- match(fl, srv_fleets)
      ObsSrvIdx[r, yi, se, sf] <- cpue$obs[i]
      ObsSrvIdx_SE[r, yi, se, sf] <- cpue$se_log[i]
      UseSrvIdx[r, yi, se, sf] <- use
    } else {
      f <- match(fl, fish_fleets)
      if(is.na(f)) next
      ObsFishIdx[r, yi, se, f] <- cpue$obs[i]
      ObsFishIdx_SE[r, yi, se, f] <- cpue$se_log[i]
      UseFishIdx[r, yi, se, f] <- use
    } # end if else

  } # end i loop

  # CPUEinfo gives units 0 for numbers and 1 for biomass, which is what SPoRC's
  # idx_type switch wants, and an error type of 0 for lognormal
  cpue_info <- d$CPUEinfo
  ci_col <- function(nm) cpue_info[[grep(paste0("^", nm, "$"), names(cpue_info), ignore.case = TRUE)[1]]]
  ci_fleet <- ci_col("fleet")
  idx_units <- ci_col("units")[match(seq_len(nrow(fleetinfo)), ci_fleet)]
  idx_errtype <- ci_col("errtype")[match(seq_len(nrow(fleetinfo)), ci_fleet)]
  srv_idx_type <- idx_units[srv_fleets]
  fish_idx_type <- idx_units[fish_fleets]

  # fraction of the year elapsed when each survey is taken, read off each observation row.
  # fleetinfo$surveytiming is only a fallback for a fleet with no rows and is routinely stale
  srv_month <- rep(NA_real_, n_srv)

  for(k in seq_along(srv_fleets)) {

    # a model can have no age data at all, so the age table may be absent
    ac_month <- if(NROW(d$agecomp) > 0) d$agecomp$month[abs(d$agecomp$fleet) == srv_fleets[k]] else numeric(0)
    m <- c(cpue$month[abs(cpue[[idx_col]]) == srv_fleets[k]], ac_month)
    srv_month[k] <- if(length(m) > 0) stats::median(m) else max(fleetinfo$surveytiming[srv_fleets[k]], 1)

  } # end k loop

  t_srv <- (srv_month - 1) / 12

  # Compositions ---------------------------------------------------------------

  # a comp row's numbers start immediately after its sample size column
  comp_vals <- function(tab, row, n_col) as.numeric(tab[row, (which(names(tab) == n_col) + 1):ncol(tab)])

  ObsFishAgeComps <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_obs_ages, n_sexes, n_fish))
  ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish))
  UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))

  ObsSrvAgeComps <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_obs_ages, n_sexes, n_srv))
  ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv))
  UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))

  ObsFish_caal <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_obs_ages, n_sexes, n_fish))
  ISS_Fish_caal <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish))
  UseFish_caal <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_fish))

  ObsSrv_caal <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_obs_ages, n_sexes, n_srv))
  ISS_Srv_caal <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv))
  UseSrv_caal <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_srv))

  # which ageing error definition each fleet's age data uses, by year: SS3 reads
  # it per row, and a fleet can switch definitions partway through a series
  ageerr_fish <- matrix(NA_integer_, n_yrs, n_fish)
  ageerr_srv <- matrix(NA_integer_, n_yrs, n_srv)
  ageerr_fish_caal <- matrix(NA_integer_, n_yrs, n_fish)
  ageerr_srv_caal <- matrix(NA_integer_, n_yrs, n_srv)

  ac <- d$agecomp

  for(i in seq_len(NROW(ac))) {

    fl <- abs(ac$fleet[i])
    yr <- abs(ac$year[i])
    yi <- yidx(yr)
    if(is.na(yi)) next

    r <- fleet_area[fl]
    se <- if(n_seas == 1) 1 else sidx(ceiling(ac$month[i] / (12 / n_seas)))
    use <- as.numeric(is_fit(ac$year[i], ac$fleet[i]))
    slots <- sex_slots(ac$sex[i])
    vals <- comp_blocks(comp_vals(ac, i, "Nsamp"), n_obs_ages, ac$sex[i])

    # SS3 writes a marginal age composition with either a negative Lbin_lo or a range spanning
    # every bin. anything else is conditional age at length, one length bin per row
    spans_all <- ac$Lbin_lo[i] <= min(lens_lower) && ac$Lbin_hi[i] >= max(lens_lower)
    marginal <- ac$Lbin_lo[i] < 0 || spans_all

    if(!marginal) {
      if(ac$Lbin_hi[i] != ac$Lbin_lo[i]) {
        stop("age comp row ", i, " spans length bins ", ac$Lbin_lo[i], " to ", ac$Lbin_hi[i],
             ", which is neither one bin nor the whole range")
      }
      li <- match(ac$Lbin_lo[i], lens_lower)
      if(is.na(li)) stop("CAAL Lbin_lo ", ac$Lbin_lo[i], " is not a composition length bin edge")
    } # end if

    if(fl %in% srv_fleets) {

      sf <- match(fl, srv_fleets)
      if(use == 1) { if(marginal) ageerr_srv[yi, sf] <- abs(ac$ageerr[i]) else ageerr_srv_caal[yi, sf] <- abs(ac$ageerr[i]) }

      if(marginal) {
        ObsSrvAgeComps[r, yi, se, , , sf] <- vals
        ISS_SrvAgeComps[r, yi, se, slots, sf] <- ac$Nsamp[i]
        UseSrvAgeComps[r, yi, se, sf] <- use
      } else {
        # each sex gets its own row for a length bin, so write only this row's
        # block or the second row would zero out the first
        ObsSrv_caal[r, yi, se, li, , slots, sf] <- vals[, slots]
        ISS_Srv_caal[r, yi, se, li, slots, sf] <- ac$Nsamp[i]
        UseSrv_caal[r, yi, se, li, sf] <- max(UseSrv_caal[r, yi, se, li, sf], use)
      } # end if else

    } else {

      f <- match(fl, fish_fleets)
      if(is.na(f)) next
      if(use == 1) { if(marginal) ageerr_fish[yi, f] <- abs(ac$ageerr[i]) else ageerr_fish_caal[yi, f] <- abs(ac$ageerr[i]) }

      if(marginal) {
        ObsFishAgeComps[r, yi, se, , , f] <- vals
        ISS_FishAgeComps[r, yi, se, slots, f] <- ac$Nsamp[i]
        UseFishAgeComps[r, yi, se, f] <- use
      } else {
        ObsFish_caal[r, yi, se, li, , slots, f] <- vals[, slots]
        ISS_Fish_caal[r, yi, se, li, slots, f] <- ac$Nsamp[i]
        UseFish_caal[r, yi, se, li, f] <- max(UseFish_caal[r, yi, se, li, f], use)
      } # end if else

    } # end if else

  } # end i loop

  ObsFishLenComps <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish))
  ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish))
  UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))

  ObsSrvLenComps <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv))
  ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv))
  UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))

  lc <- d$lencomp

  for(i in seq_len(NROW(lc))) {

    fl <- abs(lc$fleet[i])
    yr <- abs(lc$year[i])
    yi <- yidx(yr)
    if(is.na(yi)) next

    r <- fleet_area[fl]
    se <- if(n_seas == 1) 1 else sidx(ceiling(lc$month[i] / (12 / n_seas)))
    use <- as.numeric(is_fit(lc$year[i], lc$fleet[i]))
    slots <- sex_slots(lc$sex[i])

    vals <- comp_blocks(comp_vals(lc, i, "Nsamp"), n_lens, lc$sex[i])

    if(fl %in% srv_fleets) {
      sf <- match(fl, srv_fleets)
      ObsSrvLenComps[r, yi, se, , , sf] <- vals
      ISS_SrvLenComps[r, yi, se, slots, sf] <- lc$Nsamp[i]
      UseSrvLenComps[r, yi, se, sf] <- use
    } else {
      f <- match(fl, fish_fleets)
      if(is.na(f)) next
      ObsFishLenComps[r, yi, se, , , f] <- vals
      ISS_FishLenComps[r, yi, se, slots, f] <- lc$Nsamp[i]
      UseFishLenComps[r, yi, se, f] <- use
    } # end if else

  } # end i loop

  # Biology --------------------------------------------------------------------

  ae_defs <- ss3_ageing_error(d, ages, obs_ages)
  natmort <- ss3_natmort(ctl, ages, years, n_pop, n_regions, n_sexes)
  waa <- ss3_waa(d, ctl, p$wtatage, ages, years, n_pop, n_regions, n_seas, n_sexes,
                 n_fish, n_srv, fish_fleets, srv_fleets, pop_lens_lower, waa_fallback)

  # Setup_Mod_Biologicals takes [n_ages, n_obs_ages, n_fleets] or the time-varying form with a
  # leading year dim. emit the time-varying one; a year with no age data uses the first definition
  ae_by_fleet <- function(which_def) {

    n_fl <- ncol(which_def)
    out <- array(0, dim = c(n_yrs, length(ages), n_obs_ages, n_fl))

    for(k in seq_len(n_fl)) {

      # a year with no age data holds the most recent definition forward, so a
      # switch partway through a series holds for every later year
      last_def <- 1

      for(y in seq_len(n_yrs)) {
        if(!is.na(which_def[y, k])) last_def <- which_def[y, k]
        out[y, , , k] <- ae_defs[last_def, , ]
      } # end y loop

    } # end k loop

    out

  } # end function

  mg <- ctl$MG_parms

  growth <- lapply(seq_len(n_sexes), function(s) {
    sx <- c("Fem", "Mal")[s]
    rows <- mg[grep(paste0("_", sx, "_GP_1$"), rownames(mg)), , drop = FALSE]
    list(L1 = rows[grep("L_at_Amin", rownames(rows)), "INIT"],
         L2 = rows[grep("L_at_Amax", rownames(rows)), "INIT"],
         K = rows[grep("VonBert_K", rownames(rows)), "INIT"],
         CV1 = rows[grep("CV_young", rownames(rows)), "INIT"],
         CV2 = rows[grep("CV_old", rownames(rows)), "INIT"],
         wtlen = rows[grep("^Wtlen_[12]_", rownames(rows)), "INIT"],
         est = rows[grep("L_at_Amin|L_at_Amax|VonBert_K|CV_young|CV_old", rownames(rows)), "PHASE"] > 0)
  })

  # Parameters -----------------------------------------------------------------

  sr <- ctl$SR_parms
  sr_val <- function(nm) if(nm %in% rownames(sr)) sr[nm, "INIT"] else NA

  rec <- list(
    SR_function = ctl$SR_function,
    ln_R0 = sr_val("SR_LN(R0)"),
    h = sr_val("SR_BH_steep"),
    sigmaR = sr_val("SR_sigmaR"),
    autocorr = sr_val("SR_autocorr"),
    main_first = ctl$MainRdevYrFirst,
    main_last = ctl$MainRdevYrLast,
    early_start = ctl$recdev_early_start,
    early_phase = ctl$recdev_early_phase,
    bias_years = c(ctl$last_early_yr_nobias_adj, ctl$first_yr_fullbias_adj,
                             ctl$last_yr_fullbias_adj, ctl$first_recent_yr_nobias_adj),
    max_bias_adj = ctl$max_bias_adj
  )

  fmort <- list(method = ctl$F_Method, maxF = ctl$maxF, init_F = ctl$init_F)

  sel <- list(
    age_types = ctl$age_selex_types,
    age = ctl$age_selex_parms,
    size_types = ctl$size_selex_types,
    size = ctl$size_selex_parms,
    age_tv = ctl$age_selex_parms_tv,
    size_tv = ctl$size_selex_parms_tv
  )

  q <- list(options = ctl$Q_options, parms = ctl$Q_parms)

  # comp error is 0 for multinomial and 1 for Dirichlet-multinomial, with
  # ParmSelect naming the row of dirichlet_parms a fleet uses
  comp <- list(
    age_info = d$age_info,
    len_info = d$len_info,
    dirichlet = ctl$dirichlet_parms,
    var_adj = ctl$Variance_adjustment_list
  )

  list(
    source = run_dir,
    years = years, ages = ages, obs_ages = obs_ages,
    # lens is the composition grid; pop_lens is SS3's finer internal grid, which
    # is the model's length axis in assessments where the two differ
    lens_lower = lens_lower, lens = lens_mid,
    pop_lens_lower = pop_lens_lower, pop_lens = bin_midpoints(pop_lens_lower),
    n_pop = n_pop, n_regions = n_regions, n_seas = n_seas, n_sexes = n_sexes,
    n_fish_fleets = n_fish, n_srv_fleets = n_srv,
    fleetnames = fleetinfo$fleetname, fish_fleets = fish_fleets, srv_fleets = srv_fleets,
    fleet_area = fleet_area, spawn_month = d$spawn_month, t_srv = t_srv,
    srv_idx_type = srv_idx_type, fish_idx_type = fish_idx_type, idx_errtype = idx_errtype,
    catch_units = catch_units,

    ObsCatch = ObsCatch, UseCatch = UseCatch, catch_se = catch_se,
    init_equil_catch = init_equil_catch,
    ObsSrvIdx = ObsSrvIdx, ObsSrvIdx_SE = ObsSrvIdx_SE, UseSrvIdx = UseSrvIdx,
    ObsFishIdx = ObsFishIdx, ObsFishIdx_SE = ObsFishIdx_SE, UseFishIdx = UseFishIdx,
    ObsFishAgeComps = ObsFishAgeComps, ISS_FishAgeComps = ISS_FishAgeComps, UseFishAgeComps = UseFishAgeComps,
    ObsSrvAgeComps = ObsSrvAgeComps, ISS_SrvAgeComps = ISS_SrvAgeComps, UseSrvAgeComps = UseSrvAgeComps,
    ObsFishLenComps = ObsFishLenComps, ISS_FishLenComps = ISS_FishLenComps, UseFishLenComps = UseFishLenComps,
    ObsSrvLenComps = ObsSrvLenComps, ISS_SrvLenComps = ISS_SrvLenComps, UseSrvLenComps = UseSrvLenComps,
    ObsFish_caal = ObsFish_caal, ISS_Fish_caal = ISS_Fish_caal, UseFish_caal = UseFish_caal,
    ObsSrv_caal = ObsSrv_caal, ISS_Srv_caal = ISS_Srv_caal, UseSrv_caal = UseSrv_caal,

    WAA = waa$WAA, MatAA = waa$MatAA, WAA_fish = waa$WAA_fish, WAA_srv = waa$WAA_srv,
    mat_is_fecundity = waa$mat_is_fecundity,
    Fixed_natmort = natmort,
    AgeingError = ae_defs[1, , ], AgeingError_defs = ae_defs,
    AgeingError_fish = ae_by_fleet(ageerr_fish), AgeingError_srv = ae_by_fleet(ageerr_srv),
    AgeingError_fish_caal = ae_by_fleet(ageerr_fish_caal),
    AgeingError_srv_caal = ae_by_fleet(ageerr_srv_caal),
    ageerr_fish = ageerr_fish, ageerr_srv = ageerr_srv,
    ageerr_fish_caal = ageerr_fish_caal, ageerr_srv_caal = ageerr_srv_caal,

    growth = growth, rec = rec, fmort = fmort, sel = sel, q = q, comp = comp,
    ctl = ctl, dat_raw = d
  )

} # end function
