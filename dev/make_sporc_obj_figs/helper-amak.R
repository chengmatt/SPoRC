# Shared readers and reference reconstruction for the AMAK assessments bridged
# in this directory. Two jobs live here.
#
# The first is parsing. AMAK's .dat files carry comments that are decorative
# rather than structural, so the only reliable contract is the read order in
# amak.tpl's DATA_SECTION, and the readers below are positional for that reason.
#
# The second is the comparison target itself. For_R.rep prints six significant
# digits, which puts a 1e-4 percent floor under everything and hides whether a
# residual is structural or rounding. amak_dynamics and amak_objective rebuild
# the population and every objective component from amak.par at full precision
# instead, and reproduce For_R.rep's $Like_Comp to its printed precision, which
# is what makes them safe to compare a bridge against.
#
# Sourced by make_bsai_atka_data_object.R.
# Creator: Matthew LH. Cheng
# Date Created: 8/17/26

# Strip comments and blank lines, returning one numeric stream ---------------
amak_tokens <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("#.*$", "", lines)
  toks <- unlist(strsplit(trimws(lines), "[ \t]+"))
  toks[nzchar(toks)]
}

# A cursor over that stream. AMAK reads strictly in order, so a stateful
# reader mirrors the TPL line for line and makes a misread obvious.
amak_reader <- function(path, numeric_only = TRUE) {
  toks <- amak_tokens(path)
  i <- 0
  list(
    num = function(n = 1) {
      out <- as.numeric(toks[(i + 1):(i + n)])
      i <<- i + n
      out
    },
    int = function(n = 1) {
      out <- as.integer(toks[(i + 1):(i + n)])
      i <<- i + n
      out
    },
    chr = function(n = 1) {
      out <- toks[(i + 1):(i + n)]
      i <<- i + n
      out
    },
    mat = function(nrow, ncol) {
      out <- matrix(as.numeric(toks[(i + 1):(i + nrow * ncol)]), nrow = nrow, ncol = ncol, byrow = TRUE)
      i <<- i + nrow * ncol
      out
    },
    left = function() length(toks) - i,
    pos = function() i
  )
}

#' Read an AMAK data file (am2024.dat)
read_amak_dat <- function(path) {
  r <- amak_reader(path)
  d <- list()

  d$styr <- r$int()
  d$endyr <- r$int()
  d$rec_age <- r$int()
  d$oldest_age <- r$int()
  d$nages <- d$oldest_age - d$rec_age + 1
  d$nyrs <- d$endyr - d$styr + 1
  d$years <- d$styr:d$endyr
  d$ages <- d$rec_age:d$oldest_age

  d$nlength <- r$int()
  d$len_bins <- r$num(d$nlength)

  d$nfsh <- r$int()
  d$fshname <- r$chr(1)

  d$catch_bio <- r$mat(d$nfsh, d$nyrs)
  d$catch_bio_sd <- r$mat(d$nfsh, d$nyrs)

  d$nyrs_fsh_age <- r$int(d$nfsh)
  d$nyrs_fsh_length <- r$int(d$nfsh)
  d$yrs_fsh_age <- r$mat(d$nfsh, d$nyrs_fsh_age)
  d$yrs_fsh_length <- if(sum(d$nyrs_fsh_length) > 0) r$mat(d$nfsh, d$nyrs_fsh_length) else matrix(numeric(0), d$nfsh, 0)
  d$n_sample_fsh_age <- r$mat(d$nfsh, d$nyrs_fsh_age)
  d$n_sample_fsh_length <- if(sum(d$nyrs_fsh_length) > 0) r$mat(d$nfsh, d$nyrs_fsh_length) else matrix(numeric(0), d$nfsh, 0)
  d$oac_fsh <- r$mat(d$nyrs_fsh_age, d$nages) # nfsh == 1 throughout this model
  d$olc_fsh <- if(sum(d$nyrs_fsh_length) > 0) r$mat(d$nyrs_fsh_length, d$nlength) else matrix(numeric(0), 0, d$nlength)
  d$wt_fsh <- r$mat(d$nyrs, d$nages)

  d$nind <- r$int()
  d$indname <- r$chr(1)
  d$nyrs_ind <- r$int(d$nind)
  d$yrs_ind <- r$mat(d$nind, d$nyrs_ind)
  d$mo_ind <- r$num(d$nind)
  d$obs_ind <- r$mat(d$nind, d$nyrs_ind)
  d$obs_se_ind <- r$mat(d$nind, d$nyrs_ind)

  d$nyrs_ind_age <- r$int(d$nind)
  d$nyrs_ind_length <- r$int(d$nind)
  d$yrs_ind_age <- r$mat(d$nind, d$nyrs_ind_age)
  d$yrs_ind_length <- if(sum(d$nyrs_ind_length) > 0) r$mat(d$nind, d$nyrs_ind_length) else matrix(numeric(0), d$nind, 0)
  d$n_sample_ind_age <- r$mat(d$nind, d$nyrs_ind_age)
  d$n_sample_ind_length <- if(sum(d$nyrs_ind_length) > 0) r$mat(d$nind, d$nyrs_ind_length) else matrix(numeric(0), d$nind, 0)
  d$oac_ind <- r$mat(d$nyrs_ind_age, d$nages)
  d$olc_ind <- if(sum(d$nyrs_ind_length) > 0) r$mat(d$nyrs_ind_length, d$nlength) else matrix(numeric(0), 0, d$nlength)
  d$wt_ind <- r$mat(d$nyrs, d$nages)

  d$wt_pop <- r$num(d$nages)
  d$maturity_in <- r$num(d$nages)
  # amak.tpl:290 halves an ogive that reaches one, which is how a female
  # fraction is folded into a single-sex model.
  d$maturity <- if(max(d$maturity_in) > 0.9) d$maturity_in / 2 else d$maturity_in
  d$wt_mature <- d$wt_pop * d$maturity

  d$spawnmo <- r$num()
  d$spmo_frac <- (d$spawnmo - 1) / 12
  d$ind_month_frac <- (d$mo_ind - 1) / 12

  d$age_err <- r$mat(d$nages, d$nages)

  # Derived quantities the TPL builds in its LOCAL_CALCS
  d$styr_rec <- d$styr - d$nages + 1
  d$styr_sp <- d$styr_rec - d$rec_age - 1
  d$catch_bio_lva <- log(d$catch_bio_sd^2 + 1)
  d$catch_bio_lsd <- sqrt(d$catch_bio_lva)
  d$obs_lse_ind <- sqrt(log((d$obs_se_ind / d$obs_ind)^2 + 1))

  # Composition data are normalized to sum to one, and the multinomial offset
  # is built from the normalized observations.
  d$oac_fsh <- d$oac_fsh / rowSums(d$oac_fsh)
  d$oac_ind <- d$oac_ind / rowSums(d$oac_ind)
  d$offset_fsh <- -sum(d$n_sample_fsh_age[1, ] * rowSums((d$oac_fsh + 0.001) * log(d$oac_fsh + 0.001)))
  d$offset_ind <- -sum(d$n_sample_ind_age[1, ] * rowSums((d$oac_ind + 0.001) * log(d$oac_ind + 0.001)))

  if(r$left() != 0) warning("read_amak_dat: ", r$left(), " tokens left unread in ", basename(path))
  d
}

#' Read an AMAK control file (amak.dat)
#'
#' Only the branches Model 16.0b actually takes are implemented: selectivity
#' option 1 (age coefficients) for both the fishery and the index. Any other
#' option stops rather than silently misreading the stream that follows.
read_amak_ctl <- function(path, dat) {
  toks <- amak_tokens(path)
  # The first two tokens are the data file name and the model name.
  i <- 2
  nxt <- function(n = 1) { out <- as.numeric(toks[(i + 1):(i + n)]); i <<- i + n; out }

  ct <- list()
  ct$model_name <- toks[2]
  nfsh_and_ind <- dat$nfsh + dat$nind
  ct$sel_map <- matrix(nxt(2 * nfsh_and_ind), nrow = 2, byrow = TRUE)

  ct$SrType <- nxt()
  ct$use_age_err <- nxt()
  ct$retro <- nxt()
  ct$steepnessprior <- nxt(); ct$cvsteepnessprior <- nxt(); ct$phase_srec <- nxt()
  ct$sigmarprior <- nxt(); ct$cvsigmarprior <- nxt(); ct$phase_sigmar <- nxt()
  ct$styr_rec_est <- nxt(); ct$endyr_rec_est <- nxt()
  ct$Linfprior <- nxt(); ct$cvLinfprior <- nxt(); ct$phase_Linf <- nxt()
  ct$kprior <- nxt(); ct$cvkprior <- nxt(); ct$phase_k <- nxt()
  ct$Loprior <- nxt(); ct$cvLoprior <- nxt(); ct$phase_Lo <- nxt()
  ct$sdageprior <- nxt(); ct$cvsdageprior <- nxt(); ct$phase_sdage <- nxt()
  ct$natmortprior <- nxt(); ct$cvnatmortprior <- nxt(); ct$phase_M <- nxt()

  ct$npars_Mage <- nxt()
  if(ct$npars_Mage > 0) { ct$ages_M_changes <- nxt(ct$npars_Mage); ct$Mage_in <- nxt(ct$npars_Mage) }
  ct$phase_Mage <- nxt()
  ct$phase_rw_M <- nxt()
  ct$npars_rw_M <- nxt()
  if(ct$npars_rw_M > 0) { ct$yrs_rw_M <- nxt(ct$npars_rw_M); ct$sigma_rw_M <- nxt(ct$npars_rw_M) }

  ct$qprior <- nxt(dat$nind); ct$cvqprior <- nxt(dat$nind); ct$phase_q <- nxt(dat$nind)
  ct$q_power_prior <- nxt(dat$nind); ct$cvq_power_prior <- nxt(dat$nind); ct$phase_q_power <- nxt(dat$nind)
  ct$phase_rw_q <- nxt(dat$nind)
  ct$npars_rw_q <- nxt(dat$nind)
  if(sum(ct$npars_rw_q) > 0) { ct$yrs_rw_q <- nxt(sum(ct$npars_rw_q)); ct$sigma_rw_q <- nxt(sum(ct$npars_rw_q)) }
  ct$q_age_min <- nxt(dat$nind); ct$q_age_max <- nxt(dat$nind)
  # Mapped from age to age index, as at amak.tpl:449
  ct$q_age_min_idx <- ct$q_age_min - dat$rec_age + 1
  ct$q_age_max_idx <- ct$q_age_max - dat$rec_age + 1

  ct$cv_catchbiomass <- nxt()
  ct$catchbiomass_pen <- 1 / (2 * ct$cv_catchbiomass^2)
  ct$nproj_yrs <- nxt()

  ct$fsh_sel_opt <- nxt()
  if(ct$fsh_sel_opt != 1) stop("read_amak_ctl only implements fishery selectivity option 1 (coefficients).")
  ct$nselages_fsh <- nxt()
  ct$phase_sel_fsh <- nxt()
  # Both shape weights are entered as standard deviations and converted to
  # precisions, but at different points in the file: the curvature sigma at
  # amak.tpl:948 (only when the coefficients are estimated) and the dome sigma
  # on read at amak.tpl:610. The control file comments call them "sigma", the
  # likelihood uses them as weights.
  ct$curv_sigma_fsh <- nxt()
  ct$curv_pen_fsh <- if(ct$phase_sel_fsh > 0) 1 / (2 * ct$curv_sigma_fsh^2) else ct$curv_sigma_fsh
  ct$seldec_pen_fsh <- nxt()^2
  ct$n_sel_ch_fsh <- nxt() + 1 # styr is always the first block
  ct$yrs_sel_ch_fsh <- c(dat$styr, nxt(ct$n_sel_ch_fsh - 1))
  ct$sel_sigma_fsh <- c(NA, nxt(ct$n_sel_ch_fsh - 1))
  ct$sel_fsh_init <- nxt(dat$nages)

  ct$ind_sel_opt <- nxt()
  if(ct$ind_sel_opt != 1) stop("read_amak_ctl only implements index selectivity option 1 (coefficients).")
  ct$nselages_ind <- nxt()
  ct$phase_sel_ind <- nxt()
  ct$curv_sigma_ind <- nxt()
  ct$curv_pen_ind <- if(ct$phase_sel_ind > 0) 1 / (2 * ct$curv_sigma_ind^2) else ct$curv_sigma_ind
  ct$seldec_pen_ind <- nxt()^2
  ct$n_sel_ch_ind <- nxt() + 1
  ct$yrs_sel_ch_ind <- c(dat$styr, if(ct$n_sel_ch_ind > 1) nxt(ct$n_sel_ch_ind - 1) else numeric(0))
  ct$sel_sigma_ind <- c(NA, if(ct$n_sel_ch_ind > 1) nxt(ct$n_sel_ch_ind - 1) else numeric(0))
  ct$sel_ind_init <- nxt(dat$nages)

  ct$test <- nxt()
  if(ct$test != 123456789) warning("read_amak_ctl: trailing test value is ", ct$test, ", not 123456789. The control file may have been misread.")

  ct$seldecage <- floor(dat$nages / 2) # int(nages/2), amak.tpl:492
  ct$nrecs_est <- ct$endyr_rec_est - ct$styr_rec_est + 1
  ct
}

#' Read an AMAK parameter file (amak.par)
read_amak_par <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr <- lines[1]
  out <- list()
  out$objective <- as.numeric(sub(".*Objective function value *= *([-0-9.eE+]+).*", "\\1", hdr))
  out$n_pars <- as.numeric(sub(".*Number of parameters *= *([0-9]+).*", "\\1", hdr))
  out$max_grad <- as.numeric(sub(".*Maximum gradient component *= *([-0-9.eE+]+).*", "\\1", hdr))

  nm_idx <- grep("^#", lines)
  nm_idx <- nm_idx[nm_idx > 1]
  for(k in seq_along(nm_idx)) {
    start <- nm_idx[k] + 1
    end <- if(k < length(nm_idx)) nm_idx[k + 1] - 1 else length(lines)
    nm <- gsub("^# *|:.*$", "", lines[nm_idx[k]])
    nm <- gsub("\\[|\\]", "_", nm)
    nm <- sub("_$", "", nm)
    body <- lines[start:end]
    vals <- lapply(body, function(x) as.numeric(unlist(strsplit(trimws(x), "[ \t]+"))))
    vals <- vals[lengths(vals) > 0]
    out[[nm]] <- if(length(vals) == 1) vals[[1]] else do.call(rbind, vals)
  }
  out
}

#' Read For_R.rep, which AMAK writes as `$name` followed by its values
read_amak_R_rep <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr <- grep("^\\s*\\$", lines)
  out <- list()
  for(k in seq_along(hdr)) {
    nm <- sub("^\\s*\\$", "", trimws(lines[hdr[k]]))
    start <- hdr[k] + 1
    end <- if(k < length(hdr)) hdr[k + 1] - 1 else length(lines)
    if(end < start) { out[[nm]] <- numeric(0); next }
    body <- lines[start:end]
    vals <- lapply(body, function(x) suppressWarnings(as.numeric(unlist(strsplit(trimws(x), "[ \t]+")))))
    vals <- vals[lengths(vals) > 0]
    if(length(vals) == 0) { out[[nm]] <- numeric(0); next }
    # A block whose rows are ragged, or whose entries are not numeric, is a
    # name table rather than a matrix; keep those as raw text.
    if(any(sapply(vals, function(v) any(is.na(v))))) {
      out[[nm]] <- trimws(body[nzchar(trimws(body))])
    } else if(length(vals) == 1) {
      out[[nm]] <- vals[[1]]
    } else if(length(unique(lengths(vals))) == 1) {
      out[[nm]] <- do.call(rbind, vals)
    } else {
      out[[nm]] <- vals
    }
  }
  out
}


# Selectivity ---------------------------------------------------------------
# Both streams use option 1: free coefficients over nselages bins, the oldest
# bin held flat to the plus group, then standardized. The fishery standardizes
# over every age; the index standardizes over the ages q is defined on, which
# is a pure rescaling that catchability absorbs.
amak_selex <- function(coffs, nselages, nages, norm_bins = NULL) {
  ls <- c(coffs[1:nselages], rep(coffs[nselages], nages - nselages))
  nb <- if(is.null(norm_bins)) seq_len(nages) else norm_bins
  ls - log(mean(exp(ls[nb])))
}

# Map each model year onto its selectivity block. AMAK advances the block
# pointer on the change year and holds it thereafter, so years past the last
# change share the last block.
amak_block_of_year <- function(years, yrs_ch) findInterval(years, yrs_ch)

#' Rebuild AMAK's population dynamics from its parameter vector
amak_dynamics <- function(dat, ctl, par) {

  nages <- dat$nages; nyrs <- dat$nyrs; years <- dat$years
  styr <- dat$styr; endyr <- dat$endyr; styr_rec <- dat$styr_rec

  M <- matrix(par$Mest, nyrs, nages)
  survtmp <- exp(-M[1, ])

  # Selectivity surfaces
  blk_fsh <- amak_block_of_year(years, ctl$yrs_sel_ch_fsh)
  log_sel_fsh <- t(sapply(blk_fsh, function(b) amak_selex(par$log_selcoffs_fsh_1[b, ], ctl$nselages_fsh, nages)))
  sel_fsh <- exp(log_sel_fsh)

  log_sel_ind_1 <- amak_selex(par$log_selcoffs_ind_1, ctl$nselages_ind, nages,
                              norm_bins = ctl$q_age_min_idx:ctl$q_age_max_idx)
  log_sel_ind <- matrix(log_sel_ind_1, nyrs, nages, byrow = TRUE)
  sel_ind <- exp(log_sel_ind)

  # Mortality
  Fmat <- par$fmort * sel_fsh
  Z <- M + Fmat
  S <- exp(-Z)

  # Equilibrium and pre-model recruitment (Get_Bzero)
  Rzero <- exp(par$log_Rzero)
  n_pre <- styr - styr_rec + 1 # 1967..1977 inclusive
  natagetmp <- matrix(0, n_pre, nages)
  natagetmp[1, 1] <- Rzero
  for(j in 2:nages) natagetmp[1, j] <- natagetmp[1, j - 1] * survtmp[j - 1]
  natagetmp[1, nages] <- natagetmp[1, nages] / (1 - survtmp[nages])

  wt_mature <- dat$wt_mature
  Bzero <- sum(wt_mature * survtmp^dat$spmo_frac * natagetmp[1, ])
  phizero <- Bzero / Rzero

  h <- par$steepness
  if(ctl$SrType != 2) stop("amak_dynamics only implements SrType 2 (Beverton-Holt).")
  alpha <- Bzero * (1 - (h - 0.2) / (0.8 * h)) / Rzero
  beta <- (5 * h - 1) / (4 * h * Rzero)

  # Sp_Biom is indexed from styr_sp; hold it in a named vector on calendar years
  sp_years <- dat$styr_sp:(endyr + 1)
  Sp_Biom <- setNames(rep(0, length(sp_years)), sp_years)
  Sp_Biom[as.character(dat$styr_sp:(styr_rec - 1))] <- Bzero

  rec_dev <- setNames(par$rec_dev, styr_rec:endyr)
  for(i in 1:(n_pre - 1)) {
    yr <- styr_rec + i - 1
    Sp_Biom[as.character(yr)] <- sum(natagetmp[i, ] * survtmp^dat$spmo_frac * wt_mature)
    natagetmp[i, 1] <- exp(rec_dev[as.character(yr)] + par$mean_log_rec)
    natagetmp[i + 1, 2:nages] <- natagetmp[i, 1:(nages - 1)] * survtmp[1:(nages - 1)]
    natagetmp[i + 1, nages] <- natagetmp[i + 1, nages] + natagetmp[i, nages] * survtmp[nages]
  }
  natagetmp[n_pre, 1] <- exp(rec_dev[as.character(styr)] + par$mean_log_rec)
  mod_rec_pre <- setNames(natagetmp[, 1], styr_rec:styr)

  # Forward dynamics
  natage <- matrix(0, nyrs + 1, nages)
  natage[1, ] <- natagetmp[n_pre, ]
  Sp_Biom[as.character(styr)] <- sum(natage[1, ] * survtmp^dat$spmo_frac * wt_mature)
  for(y in 2:nyrs) natage[y, 1] <- exp(par$mean_log_rec + rec_dev[as.character(years[y])])

  catage <- matrix(0, nyrs, nages)
  pred_catch <- rep(0, nyrs)
  for(y in 1:nyrs) {
    natage[y + 1, 2:nages] <- natage[y, 1:(nages - 1)] * S[y, 1:(nages - 1)]
    natage[y + 1, nages] <- natage[y + 1, nages] + natage[y, nages] * S[y, nages]
    catage[y, ] <- Fmat[y, ] / Z[y, ] * (1 - S[y, ]) * natage[y, ]
    pred_catch[y] <- sum(catage[y, ] * dat$wt_fsh[y, ])
    Sp_Biom[as.character(years[y])] <- sum(natage[y, ] * S[y, ]^dat$spmo_frac * wt_mature)
  }

  mod_rec <- setNames(c(mod_rec_pre[-length(mod_rec_pre)], natage[1:nyrs, 1]), styr_rec:endyr)

  # Predicted index and compositions
  i_ind <- match(dat$yrs_ind[1, ], years)
  frac <- dat$ind_month_frac[1]
  q <- exp(par$log_q_ind_1)
  q_power <- exp(par$log_q_power_ind_1)
  pred_ind <- q * sapply(i_ind, function(y) sum(natage[y, ] * S[y, ]^frac * sel_ind[y, ] * dat$wt_ind[y, ]))^q_power

  # AMAK's `age_err * p` is a matrix-vector product over the SECOND index, so
  # the observed-age axis is the ROW axis. Anything downstream that wants
  # `p %*% A` needs the transpose.
  ae <- dat$age_err
  i_fsh_age <- match(dat$yrs_fsh_age[1, ], years)
  eac_fsh <- t(sapply(i_fsh_age, function(y) {
    p <- catage[y, ] / sum(catage[y, ])
    e <- as.vector(ae %*% p)
    e / sum(e)
  }))

  i_ind_age <- match(dat$yrs_ind_age[1, ], years)
  eac_ind <- t(sapply(i_ind_age, function(y) {
    tn <- S[y, ]^frac * sel_ind[y, ] * natage[y, ]
    p <- tn / sum(tn)
    as.vector(ae %*% p)
  }))

  # Stock-recruit prediction: pred_rec(i) = SRecruit(Sp_Biom(i - rec_age))
  rec_years <- styr_rec:endyr
  pred_rec <- setNames(sapply(rec_years, function(i) {
    s <- Sp_Biom[as.character(i - dat$rec_age)]
    s / (alpha + beta * s)
  }), rec_years)

  list(M = M, sel_fsh = sel_fsh, sel_ind = sel_ind, log_sel_fsh = log_sel_fsh,
       log_sel_ind = log_sel_ind, blk_fsh = blk_fsh, F = Fmat, Z = Z, S = S,
       natage = natage, catage = catage, pred_catch = pred_catch,
       Sp_Biom = Sp_Biom, Bzero = Bzero, phizero = phizero, Rzero = Rzero,
       alpha = alpha, beta = beta, natagetmp = natagetmp,
       mod_rec = mod_rec, pred_rec = pred_rec, rec_dev = rec_dev,
       pred_ind = pred_ind, eac_fsh = eac_fsh, eac_ind = eac_ind,
       i_ind = i_ind, i_fsh_age = i_fsh_age, i_ind_age = i_ind_age,
       totbiom = as.vector(natage[1:nyrs, ] %*% dat$wt_pop))
}

#' Rebuild every AMAK objective function component
amak_objective <- function(dat, ctl, par, dyn) {

  nages <- dat$nages; nyrs <- dat$nyrs; years <- dat$years

  # 1 catch: lognormal on a log-scale variance derived from the CV
  catch_like <- sum(0.5 * (log(dat$catch_bio[1, ] + 1e-4) - log(dyn$pred_catch + 1e-4))^2 / dat$catch_bio_lva[1, ])

  # 2 fishery age comps: ADMB multinomial with the 0.001 constant on both sides
  age_like_fsh <- -sum(dat$n_sample_fsh_age[1, ] * rowSums((dat$oac_fsh + 0.001) * log(dyn$eac_fsh + 0.001))) - dat$offset_fsh

  # 4 fishery selectivity penalties, evaluated only in the block change years
  d1 <- function(x) diff(x)
  iyr_fsh <- match(ctl$yrs_sel_ch_fsh, years)
  sel_curve_fsh <- sum(sapply(iyr_fsh, function(y) ctl$curv_pen_fsh * sum(d1(d1(dyn$log_sel_fsh[y, ]))^2)))
  sel_rw_fsh <- sum(sapply(iyr_fsh[-1], function(y) 0.5 * sum((dyn$log_sel_fsh[y - 1, ] - dyn$log_sel_fsh[y, ])^2) / ctl$sel_sigma_fsh[2]^2))
  sel_dome_fsh <- sum(sapply(iyr_fsh, function(y) {
    dd <- dyn$log_sel_fsh[y, ctl$seldecage:ctl$nselages_fsh - 1] - dyn$log_sel_fsh[y, ctl$seldecage:ctl$nselages_fsh]
    0.5 * sum(pmax(dd, 0)^2) / ctl$seldec_pen_fsh
  }))
  sel_like_fsh <- sel_curve_fsh + sel_rw_fsh + sel_dome_fsh

  # 5 index: lognormal about the log-scale SE derived from the arithmetic SE
  ind_like <- sum((log(dat$obs_ind[1, ]) - log(dyn$pred_ind))^2 / (2 * dat$obs_lse_ind[1, ]^2))

  # 6 index age comps
  age_like_ind <- -sum(dat$n_sample_ind_age[1, ] * rowSums((dat$oac_ind + 0.001) * log(dyn$eac_ind + 0.001))) - dat$offset_ind

  # 8 index selectivity penalties, evaluated only in its single change year
  iyr_ind <- match(ctl$yrs_sel_ch_ind, years)
  sel_curve_ind <- sum(sapply(iyr_ind, function(y) ctl$curv_pen_ind * sum(d1(d1(dyn$log_sel_ind[y, ]))^2)))
  sel_dome_ind <- sum(sapply(iyr_ind, function(y) {
    dd <- dyn$log_sel_ind[y, ctl$seldecage:ctl$nselages_ind - 1] - dyn$log_sel_ind[y, ctl$seldecage:ctl$nselages_ind]
    0.5 * sum(pmax(dd, 0)^2) / ctl$seldec_pen_ind
  }))
  sel_like_ind <- sel_curve_ind + sel_dome_ind

  # 9 recruitment
  sigmar <- exp(par$log_sigmar)
  sigmarsq <- sigmar^2
  est_yrs <- as.character(ctl$styr_rec_est:ctl$endyr_rec_est)
  chi <- log(dyn$mod_rec[est_yrs]) - log(dyn$pred_rec[est_yrs])
  SSQRec <- sum(chi^2)
  m_sigmarsq <- SSQRec / ctl$nrecs_est
  rec_like_1 <- (SSQRec + m_sigmarsq / 2) / (2 * sigmarsq) + ctl$nrecs_est * par$log_sigmar
  rec_like_2 <- sum(par$rec_dev^2)
  n_fut <- length(par$rec_dev_future)
  rec_like_3 <- sum(par$rec_dev_future^2) / (2 * sigmarsq) + n_fut * log(sigmar)
  # ADMB's inclusive index ranges make the sum-of-squares carry one more
  # deviation than the log-sigma count on each tail.
  early <- as.character(dat$styr_rec:ctl$styr_rec_est)
  late <- as.character(ctl$endyr_rec_est:dat$endyr)
  rec_like_4 <- 0.5 * sum(dyn$rec_dev[early]^2) / sigmarsq + (ctl$styr_rec_est - dat$styr_rec) * log(sigmar) +
    0.5 * sum(dyn$rec_dev[late]^2) / sigmarsq + (dat$endyr - ctl$endyr_rec_est) * log(sigmar)
  rec_like <- rec_like_1 + rec_like_2 + rec_like_3 + rec_like_4

  # 10 F penalty
  fpen <- 1e-4 * sum((dyn$F - 0.2)^2)

  # 11 catchability prior
  post_priors_indq <- (log(exp(par$log_q_ind_1) / ctl$qprior))^2 / (2 * ctl$cvqprior^2)

  # 12 other priors: only sigmaR is active at this configuration
  post_priors <- (log(sigmar / ctl$sigmarprior))^2 / (2 * ctl$cvsigmarprior^2)

  # 13 terms AMAK adds straight onto obj_fun rather than into obj_comps
  avgsel_fsh <- apply(par$log_selcoffs_fsh_1[, 1:ctl$nselages_fsh, drop = FALSE], 1, function(p) log(mean(exp(p))))
  avgsel_ind <- log(mean(exp(par$log_selcoffs_ind_1[1:ctl$nselages_ind])))
  residual <- 0.5 * (par$log_Rzero - par$mean_log_rec)^2 + 20 * sum(avgsel_fsh^2) + 20 * sum(avgsel_ind^2)

  comps <- c(catch_like = catch_like, age_like_fsh = age_like_fsh, length_like_fsh = 0,
             sel_like_fsh = sel_like_fsh, ind_like = ind_like, age_like_ind = age_like_ind,
             length_like_ind = 0, sel_like_ind = sel_like_ind, rec_like = rec_like,
             fpen = fpen, post_priors_indq = post_priors_indq, post_priors = post_priors,
             residual = residual)

  list(comps = c(comps, total = sum(comps)),
       sel_fsh_parts = c(curve = sel_curve_fsh, rw = sel_rw_fsh, dome = sel_dome_fsh),
       sel_ind_parts = c(curve = sel_curve_ind, rw = 0, dome = sel_dome_ind),
       rec_parts = c(sigmar = sigmar, rec_like_1, rec_like_2, rec_like_3, rec_like_4),
       chi = chi, SSQRec = SSQRec, m_sigmar = sqrt(m_sigmarsq),
       avgsel_fsh = avgsel_fsh, avgsel_ind = avgsel_ind)
}
