# Readers for the flatfish model (fm.tpl) files behind the BSAI northern rock sole case study, kept
# apart from the case study script so parsing and specification stay separable (cf. helper-amak.R).
# Creator: Matthew LH. Cheng

# fm.par and writeinput.log both write a "# name" line followed by values, so
# one reader covers both. The par file's names end in a colon, the log's do not.
fm_grab <- function(txt, nm, sep = ":") {
  i <- grep(paste0("^# ", nm, sep, "?$"), txt)
  stopifnot(length(i) == 1)
  j <- i + 1
  v <- numeric(0)
  while(j <= length(txt) && !grepl("^#", txt[j])) {
    v <- c(v, scan(text = txt[j], quiet = TRUE))
    j <- j + 1
  } # end while
  v
} # end fm_grab

# fm.tpl reads its data file as one token data source in a fixed order, so the file is parsed the
# same way here. the 123456 check at the end confirms the data source stayed aligned
read_fm_dat <- function(path) {

  tok <- scan(text = gsub("#.*$", "", readLines(path)), quiet = TRUE)
  ptr <- 0
  take <- function(n) { v <- tok[(ptr + 1):(ptr + n)]; ptr <<- ptr + n; v }

  out <- list()
  out$styr <- take(1); out$endyr <- take(1); out$n_ages <- take(1)
  out$lw_pars <- take(4)
  out$n_fsh <- take(1)
  stopifnot(out$n_fsh == 1)
  out$years <- out$styr:out$endyr
  n_yrs <- length(out$years); n_ages <- out$n_ages

  out$obs_catch <- take(n_yrs)
  n_fsh_age_c <- take(1); out$n_fsh_age <- take(1)
  stopifnot(n_fsh_age_c == 0)                       # 24.2 fits split sex only
  out$yrs_fsh_age <- take(out$n_fsh_age)
  out$oac_fsh <- matrix(take(out$n_fsh_age * 2 * n_ages), nrow = out$n_fsh_age, byrow = TRUE)
  out$wt_fsh <- matrix(take(n_yrs * 2 * n_ages), nrow = n_yrs, byrow = TRUE)
  out$iss_fsh_age <- take(out$n_fsh_age)

  n_srv <- take(1)
  stopifnot(n_srv == 1)
  out$n_srv_yrs <- take(1)
  out$yrs_srv <- take(out$n_srv_yrs)
  out$mo_srv <- take(1)
  out$obs_srv <- take(out$n_srv_yrs)
  out$obs_se_srv <- take(out$n_srv_yrs)

  n_srv_age_c <- take(1); out$n_srv_age <- take(1)
  stopifnot(n_srv_age_c == 0)
  out$yrs_srv_age <- take(out$n_srv_age)
  out$iss_srv_age <- take(out$n_srv_age)
  out$oac_srv <- matrix(take(out$n_srv_age * 2 * n_ages), nrow = out$n_srv_age, byrow = TRUE)

  out$wt_srv_f <- matrix(take(n_yrs * n_ages), nrow = n_yrs, byrow = TRUE)
  out$wt_srv_m <- matrix(take(n_yrs * n_ages), nrow = n_yrs, byrow = TRUE)
  out$wt_pop_f <- matrix(take(n_yrs * n_ages), nrow = n_yrs, byrow = TRUE)
  out$wt_pop_m <- matrix(take(n_yrs * n_ages), nrow = n_yrs, byrow = TRUE)
  out$maturity <- matrix(take(n_yrs * n_ages), nrow = n_yrs, byrow = TRUE)

  out$init_age_comp <- take(1)                      # 1: initial ages free of recruitment
  n_env_cov <- take(1)
  stopifnot(n_env_cov == 0)                         # no catchability covariates in 24.2
  out$spawnmo <- take(1)
  out$srv_mo <- take(1)
  out$n_wts <- take((2012 - 1982 + 1) * n_ages)     # weight sample sizes, unused under empirical weights
  out$growth_cov <- take(n_yrs)
  stopifnot(take(1) == 123456)

  out
} # end read_fm_dat

# The control file is a bare sequence of values with names in comments, so it is
# read positionally. Only the entries the case study uses are named.
read_fm_ctl <- function(path) {
  v <- scan(text = gsub("#.*$", "", readLines(path)), quiet = TRUE)
  list(
    growth_option = v[1],
    nselages = v[34],
    lambda = v[35:44], # lambda(4) is a phase, not a weight
    styr_sr = v[45],
    endyr_sr = v[46],
    a50_sigma = v[26],
    slp_sigma = v[27],
    q_exp = v[28],
    q_sigma = v[29],
    m_exp = v[30],
    m_sigma = v[31],
    sigmaR_exp = v[32],
    sigmaR_sigma = v[33]
  )
} # end read_fm_ctl

# Every estimated parameter, read at full precision.
read_fm_par <- function(path) {
  txt <- readLines(path)
  g <- function(nm) fm_grab(txt, nm)
  # the header line names its values in words, so each is taken by its label
  after <- function(lab) as.numeric(sub(paste0("^.*", lab, " *= *([0-9.eE+-]+).*$"), "\\1", txt[1]))
  list(
    ln_q = g("ln_q_srv"),
    M_f = g("natmort_f"),
    M_m = g("natmort_m"),
    mean_log_rec = g("mean_log_rec"),
    rec_dev = g("rec_dev"),
    mean_log_init = g("mean_log_init"),
    init_dev_f = g("init_dev_f"),
    init_dev_m = g("init_dev_m"),
    log_avg_fmort = g("log_avg_fmort"),
    fmort_dev = g("fmort_dev"),
    sel_slope_fsh_f = g("sel_slope_fsh_f"),
    sel50_fsh_f = g("sel50_fsh_f"),
    slope_devs_f = g("sel_slope_fsh_devs_f"),
    sel50_devs_f = g("sel50_fsh_devs_f"),
    sel_slope_fsh_m = g("sel_slope_fsh_m"),
    sel50_fsh_m = g("sel50_fsh_m"),
    slope_devs_m = g("sel_slope_fsh_devs_m"),
    sel50_devs_m = g("sel50_fsh_devs_m"),
    male_sel_offset = g("male_sel_offset"),
    sel_slope_srv = g("sel_slope_srv"),
    sel50_srv = g("sel50_srv"),
    sel_slope_srv_m = g("sel_slope_srv_m"),
    sel50_srv_m = g("sel50_srv_m"),
    R_logalpha = g("R_logalpha"),
    R_logbeta = g("R_logbeta"),
    n_par = after("Number of parameters"),
    objective = after("Objective function value"),
    max_grad = after("Maximum gradient component")
  )
} # end read_fm_par

# fm.std holds the reported time series with their standard errors.
read_fm_std <- function(path) {
  x <- utils::read.table(text = readLines(path)[-1], col.names = c("index", "name", "value", "sd"))
  pick <- function(nm) list(value = x$value[x$name == nm], sd = x$sd[x$name == nm])
  list(SSB = pick("SSB"), TotBiom = pick("TotBiom"), pred_rec = pick("pred_rec"))
} # end read_fm_std

# fm.rep holds nLogPosterior, the likelihood components at the estimate. The
# report is written once per phase, so the last block is the converged one.
read_fm_rep <- function(path) {
  txt <- readLines(path)
  v <- scan(text = txt[tail(grep("^nLogPosterior", txt), 1) + 1], quiet = TRUE)
  list(
    wt_like = v[1:3],
    wt_fut_like = v[4],
    wt_msy_like = v[5],
    init_like = v[6],
    srv = v[7],
    catch = v[8],
    fsh_age = v[9],
    srv_age = v[10],
    rec = v[11],
    init = v[12],
    sr = v[13],
    rec_fut = v[14],
    sel_slope = v[15],
    sel_a50 = v[16],
    q_prior = v[17],
    sigmaR_prior = v[18],
    m_prior = v[19],
    fpen = v[20]
  )
} # end read_fm_rep

# The predicted catch series from the legacy report, likewise the last block.
read_fm_legacy_catch <- function(path) {
  txt <- readLines(path)
  scan(text = txt[tail(grep("predicted_catch_bomass", txt), 1) + 1], quiet = TRUE)
} # end read_fm_legacy_catch
