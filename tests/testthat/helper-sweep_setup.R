# A model built for sweeping the option surface rather than for testing one feature.
#
# Every dimension has a distinct extent, so an index walking the wrong dim lands somewhere visible, and every
# data source is on, so an option is never dead for want of data. Both would otherwise read as "no change".
#
# The stages take override lists rather than dots, so a sweep can reach one argument of one stage.

# Distinct extents for every dimension that has more than one cell. Sharing an
# extent between two dimensions is what lets a transposed or mis-strided index
# return the right shape, so these stay pairwise distinct even when a smaller
# model would run faster.
sweep_dims <- list(
  n_yrs = 13,
  n_ages = 7,
  n_regions = 3,
  n_sexes = 2,
  n_fish_fleets = 5,
  n_srv_fleets = 1,
  n_seas = 1,
  n_pop = 1,
  natal_region = NA
)

#' Build the sweep input list
#'
#' @param dims Named list overriding \code{sweep_dims}.
#' @param rec,biol,move,tag,catch,fishidx,srvidx,fishsel,srvsel,wt Named lists of
#'   arguments overriding the defaults passed to the corresponding
#'   \code{Setup_Mod_*} stage.
#' @param stop_after Name of a stage, or \code{NULL}. When given, the build
#'   returns as soon as that stage completes, which lets a sweep exercise a stage
#'   whose arguments make a later stage unbuildable.
#'
#' @keywords internal
sweep_input <- function(
  dims = list(),
  rec = list(),
  biol = list(),
  move = list(),
  tag = list(),
  catch = list(),
  fishidx = list(),
  srvidx = list(),
  fishsel = list(),
  srvsel = list(),
  wt = list(),
  stop_after = NULL
) {

  d <- utils::modifyList(sweep_dims, dims)
  n_yrs <- d$n_yrs; n_ages <- d$n_ages; n_regions <- d$n_regions
  n_sexes <- d$n_sexes; n_seas <- d$n_seas; n_pop <- d$n_pop
  n_fish <- d$n_fish_fleets; n_srv <- d$n_srv_fleets
  ff <- seq_len(n_fish); sf <- seq_len(n_srv)

  # region by year by season by age by sex, with a leading population dim, is
  # the biological array layout every stage reads
  biol_d <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)
  logistic <- function(k, x0, scale = 1) scale / (1 + exp(-k * (seq_len(n_ages) - x0)))

  dim_args <- list(
    n_pop = n_pop,
    years = seq_len(n_yrs),
    ages = seq_len(n_ages),
    lens = NA,
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_seas = n_seas,
    n_fish_fleets = n_fish,
    n_srv_fleets = n_srv,
    verbose = FALSE
  )
  # passing natal_region = NA is not the same as leaving it out: the NA reaches the
  # data list and the objective indexes on it, so it is only supplied when set
  if(!all(is.na(d$natal_region))) dim_args$natal_region <- d$natal_region
  il <- do.call(Setup_Mod_Dim, dim_args)
  if(identical(stop_after, "dim")) return(il)

  il <- do.call(Setup_Mod_Rec, utils::modifyList(list(
    input_list = il, rec_model = "mean_rec", sigmaR_spec = "fix", do_rec_bias_ramp = 0,
    init_age_strc = 1,
    # one unfished recruitment per population: a single value is read position by
    # position further in and indexes past its end once there is more than one
    ln_global_R0 = rep(log(1e6), n_pop)
  ), rec))
  if(identical(stop_after, "rec")) return(il)

  il <- suppressWarnings(do.call(Setup_Mod_Biologicals, utils::modifyList(list(
    input_list = il,
    WAA = array(rep(logistic(1, 3, 5), each = prod(biol_d[1:4])), dim = biol_d),
    WAA_fish = array(rep(logistic(1, 3, 5), each = prod(biol_d[1:4])), dim = c(biol_d, n_fish)),
    WAA_srv = array(rep(logistic(1, 3, 5), each = prod(biol_d[1:4])), dim = c(biol_d, n_srv)),
    MatAA = array(rep(logistic(1, 3), each = prod(biol_d[1:4])), dim = biol_d),
    fit_lengths = 0,
    M_spec = "fix",
    Fixed_natmort = array(0.2, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes))
  ), biol)))
  if(identical(stop_after, "biol")) return(il)

  il <- do.call(Setup_Mod_Movement, utils::modifyList(list(
    input_list = il,
    use_fixed_movement = 1,
    do_recruits_move = 0,
    Fixed_Movement = if(n_regions == 1) NA else
      array(1 / n_regions, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes))
  ), move))
  if(identical(stop_after, "move")) return(il)

  il <- do.call(Setup_Mod_Tagging, utils::modifyList(list(
    input_list = il, use_conv_fish_tagging = 0
  ), tag))
  if(identical(stop_after, "tag")) return(il)

  il <- suppressWarnings(do.call(Setup_Mod_Catch_and_F, utils::modifyList(list(
    input_list = il,
    ObsCatch = array(1e4, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    UseCatch = array(1, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    sigmaC_spec = "fix",
    sigmaF_spec = "fix"
  ), catch)))
  if(identical(stop_after, "catch")) return(il)

  il <- suppressWarnings(do.call(Setup_Mod_FishIdx_and_Comps, utils::modifyList(list(
    input_list = il,
    ObsFishIdx = array(1e5, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    ObsFishIdx_SE = array(0.2, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    UseFishIdx = array(1, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    ObsFishAgeComps = array(1 / n_ages, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish)),
    UseFishAgeComps = array(1, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    ISS_FishAgeComps = array(100, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
    ObsFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_fish)),
    UseFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
    ISS_FishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
    fish_idx_type = rep("biom", n_fish),
    FishAgeComps_LikeType = rep("Multinomial", n_fish),
    FishLenComps_LikeType = rep("none", n_fish),
    FishAgeComps_Type = paste0("agg_Year_1-terminal_Fleet_", ff),
    FishLenComps_Type = paste0("none_Year_1-terminal_Fleet_", ff)
  ), fishidx)))
  if(identical(stop_after, "fishidx")) return(il)

  il <- suppressWarnings(do.call(Setup_Mod_SrvIdx_and_Comps, utils::modifyList(list(
    input_list = il,
    ObsSrvIdx = array(1e5, dim = c(n_regions, n_yrs, n_seas, n_srv)),
    ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_yrs, n_seas, n_srv)),
    UseSrvIdx = array(1, dim = c(n_regions, n_yrs, n_seas, n_srv)),
    ObsSrvAgeComps = array(1 / n_ages, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv)),
    UseSrvAgeComps = array(1, dim = c(n_regions, n_yrs, n_seas, n_srv)),
    ISS_SrvAgeComps = array(100, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
    ObsSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_srv)),
    UseSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
    ISS_SrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
    srv_idx_type = rep("abd", n_srv),
    SrvAgeComps_LikeType = rep("Multinomial", n_srv),
    SrvLenComps_LikeType = rep("none", n_srv),
    SrvAgeComps_Type = paste0("agg_Year_1-terminal_Fleet_", sf),
    SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", sf)
  ), srvidx)))
  if(identical(stop_after, "srvidx")) return(il)

  il <- do.call(Setup_Mod_Fishsel_and_Q, utils::modifyList(list(
    input_list = il,
    cont_tv_fish_sel = paste0("none_Fleet_", ff),
    fish_sel_blocks = paste0("none_Fleet_", ff),
    fish_sel_model = paste0("logist1_Fleet_", ff),
    fish_q_blocks = paste0("none_Fleet_", ff),
    fish_fixed_sel_pars_spec = rep("est_all", n_fish),
    fish_q_spec = rep("est_all", n_fish)
  ), fishsel))
  if(identical(stop_after, "fishsel")) return(il)

  il <- do.call(Setup_Mod_Srvsel_and_Q, utils::modifyList(list(
    input_list = il,
    cont_tv_srv_sel = paste0("none_Fleet_", sf),
    srv_sel_blocks = paste0("none_Fleet_", sf),
    srv_sel_model = paste0("logist1_Fleet_", sf),
    srv_q_blocks = paste0("none_Fleet_", sf),
    srv_fixed_sel_pars_spec = rep("est_all", n_srv),
    srv_q_spec = rep("est_all", n_srv)
  ), srvsel))
  if(identical(stop_after, "srvsel")) return(il)

  do.call(Setup_Mod_Weighting, utils::modifyList(list(
    input_list = il,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
    Wt_FishLenComps = array(1, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
    Wt_SrvAgeComps = array(1, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
    Wt_SrvLenComps = array(1, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv))
  ), wt))
}


# Which override slot of sweep_input feeds each Setup_Mod_* stage, so a sweep can
# go from an argument name to the call that receives it.
sweep_stage_slot <- c(
  Setup_Mod_Rec              = "rec",
  Setup_Mod_Biologicals      = "biol",
  Setup_Mod_Movement         = "move",
  Setup_Mod_Tagging          = "tag",
  Setup_Mod_Catch_and_F      = "catch",
  Setup_Mod_FishIdx_and_Comps = "fishidx",
  Setup_Mod_SrvIdx_and_Comps = "srvidx",
  Setup_Mod_Fishsel_and_Q    = "fishsel",
  Setup_Mod_Srvsel_and_Q     = "srvsel"
)


#' Build with one argument overridden
#'
#' @param stage Name of the \code{Setup_Mod_*} function receiving the argument.
#' @param arg Argument name.
#' @param value Value to pass.
#' @param dims Dimension overrides.
#' @param extra Further overrides for the same stage, as a named list.
#'
#' @return The input list, or a condition object if the build failed.
#'
#' @keywords internal
sweep_build_with <- function(
  stage,
  arg,
  value,
  dims = list(),
  extra = list(),
  other = list(),
  full = FALSE
) {
  slot <- sweep_stage_slot[[stage]]
  d <- utils::modifyList(sweep_dims, dims)
  # a spec with no live configuration to apply arrives as NULL
  if(is.null(extra)) extra <- list()
  if(is.null(other)) other <- list()

  attempt <- function(v) {
    # a sweep that only reads the map stops at the stage under test; one that has
    # to evaluate the objective needs the remaining stages too
    args <- c(list(dims = dims, stop_after = if(full) NULL else slot), other)
    args[[slot]] <- utils::modifyList(stats::setNames(list(v), arg), extra)
    # warnings are suppressed rather than caught: several specs warn about weakly
    # informed parameters and still build, and a caught warning would abort the
    # build and read as a rejected spec
    tryCatch(suppressWarnings(do.call(sweep_input, args)), error = function(e) e)
  }

  res <- attempt(value)

  # Some specs take one value per fleet. Only a few say so; the rest reject the
  # scalar by complaining about its value, or fall over on an NA comparison, so
  # the retry cannot be driven off the message and is simply always tried.
  if(inherits(res, "condition") && length(value) == 1) {
    n <- if(grepl("srv", stage, ignore.case = TRUE)) d$n_srv_fleets else d$n_fish_fleets
    if(n > 1) {
      retry <- attempt(rep(value, n))
      # a length complaint says nothing about the value, so the fleet-length
      # attempt is the informative one whether it succeeded or failed for some
      # further reason
      if(!inherits(retry, "condition") ||
         grepl("not length|length of n_|needs to have a length|entr(y|ies) for [0-9]+ fleet",
               conditionMessage(res))) res <- retry
    }
  }

  res
}


#' Whether an error is the model declining a spec its configuration cannot have
#'
#' Some specs are legal only alongside another option: sharing over selectivity
#' bins needs a bin-indexed deviation form, for instance. At any one test setup
#' configuration those specs are legitimately unbuildable, which is the check
#' working rather than a broken spec.
#'
#' @keywords internal
sweep_is_config_refusal <- function(e) {
  if(!inherits(e, "condition")) return(FALSE)
  msg <- conditionMessage(e)
  grepl(paste0("only supported when|only available when|requires|not compatible|cannot be used|",
               # the family that names another argument the caller must set too
               "supply |Pair it with|it needs |so it needs |but no starting value|",
               "please provide|must also"),
        msg)
}


#' Put a bare option value into the form its argument requires
#'
#' Several options are given as \code{"<value>_Year_<a>-<b>_Fleet_<f>"} or
#' \code{"<value>_Fleet_<f>"}, but name only the \code{<value>} half in the error
#' listing their legal settings. A sweep reading that listing therefore holds a
#' value the argument will reject on its shape rather than its content. The
#' required shape is recovered the same way the values were: by asking, and
#' reading the answer.
#'
#' @param stage,arg Stage function and argument.
#' @param value Bare value from \code{sweep_legal_specs}.
#' @param dims Dimension overrides.
#'
#' @return A character vector the argument will accept, one entry per fleet where
#'   the argument is per-fleet.
#'
#' @keywords internal
sweep_format_value <- function(stage, arg, value, dims = list()) {
  d <- utils::modifyList(sweep_dims, dims)
  n <- if(grepl("Srv|srv", arg)) d$n_srv_fleets else d$n_fish_fleets

  res <- sweep_build_with(stage, arg, value, dims = dims)
  if(!inherits(res, "condition")) return(value)
  msg <- conditionMessage(res)

  if(grepl("Year_.*Fleet_|Fleet_x", msg) && grepl("Year", msg))
    return(paste0(value, "_Year_1-terminal_Fleet_", seq_len(n)))
  if(grepl("Fleet_x|_Fleet_", msg))
    return(paste0(value, "_Fleet_", seq_len(n)))
  value
}


#' Whether an error is the identifiability guard declining a spec
#'
#' \code{check_spec_map_identifiable} refuses a spec that leaves an observation
#' error parameter with too few observations to estimate. Which specs it refuses
#' depends on the model's dimensions, so at any one test setup size some legal specs
#' are legitimately unbuildable. That is the guard working, not a broken spec, and
#' the sweeps treat it as a skip.
#'
#' @keywords internal
sweep_is_identifiability_refusal <- function(e) {
  inherits(e, "condition") &&
    grepl("informed by fewer than|cannot be estimated|which creates", conditionMessage(e))
}


#' Legal values of a spec argument, read off the package's own error message
#'
#' Every spec argument validates its input against a list it names in the error
#' it raises, so the surface is discovered by asking rather than by keeping a
#' second copy of it here that would fall out of step.
#'
#' @param stage,arg Stage function name and argument name.
#' @param dims Dimension overrides.
#'
#' @return Character vector of legal values, or \code{character(0)} when the
#'   argument does not validate itself in a form this can read.
#'
#' @keywords internal
sweep_legal_specs <- function(stage, arg, dims = list(), extra = list()) {
  d <- utils::modifyList(sweep_dims, dims)
  n <- if(grepl("Srv|srv", arg)) d$n_srv_fleets else d$n_fish_fleets

  # Arguments differ in whether they check a value's length or its content first.
  # One that checks length first answers a scalar probe with a length complaint
  # and never reaches the listing, so both shapes are probed and whichever names
  # values is the one read.
  probes <- list("__not_a_spec__")
  if(n > 1) probes[[2]] <- rep("__not_a_spec__", n)

  msg <- ""
  for(probe in probes) {
    res <- sweep_build_with(stage, arg, probe, dims = dims, extra = extra)
    if(!inherits(res, "condition")) next
    cand <- gsub("\n", " ", conditionMessage(res))
    if(grepl("one of|Valid options|Should be|must be", cand)) { msg <- cand; break }
    if(!nzchar(msg)) msg <- cand
  }
  if(!nzchar(msg)) return(character(0))

  # The message forms the package uses to name what it will accept. They differ
  # by author rather than by meaning, so all of them are read here instead of
  # keeping a second copy of the surface that would fall out of step with the first.
  patterns <- c("(?<=Should be one of these: ).*", "(?<=Valid options: ).*",
                "(?<=Must be one of: ).*", "(?<=Must be one of these: ).*",
                "(?<=Should be one of: ).*", "(?<=one of these: ).*",
                "(?<=Should be either ).*", "(?<=Should be ).*", "(?<=must be ).*")
  m <- character(0)
  for(pat in patterns) {
    m <- regmatches(msg, regexpr(pat, msg, perl = TRUE))
    if(length(m) > 0) break
  }
  if(length(m) == 0) return(character(0))

  # the list ends at the first sentence break; prose after it describes the
  # arguments rather than naming more values
  txt <- sub("[.!] .*$", "", m[1])
  txt <- sub("[.!]$", "", txt)
  parts <- trimws(unlist(strsplit(txt, ",| or |/")))

  # a trailing family name ("est_shared_f_# (where # is fleet number)") is a
  # pattern rather than a value, and quoted or bracketed prose is not a value
  parts <- gsub("^['\"`]|['\"`]$", "", parts)
  parts <- parts[!grepl("#|\\(|\\)|[[:space:]]", parts)]
  parts <- parts[nzchar(parts) & parts != "NULL" & parts != "the"]
  unique(parts)
}


#' A comparable digest of what a build produced
#'
#' Records the estimation structure (how many free parameters each block has and
#' how they are grouped) and the model configuration separately, so a sweep can
#' tell an option that changes what is estimated from one that changes what is
#' computed, and both from one that changes nothing.
#'
#' @param il An input list from \code{sweep_input}.
#'
#' @keywords internal
sweep_signature <- function(il) {
  map_sig <- lapply(il$map, function(m) {
    f <- as.integer(as.factor(as.character(m)))
    list(n_free = length(unique(stats::na.omit(f))), grouping = paste(f, collapse = ","))
  })
  par_sig <- vapply(il$par, function(p) paste0(paste(dim(p), collapse = "x"), ":", length(p)),
                    character(1))
  # the data list holds the switches the objective reads, so a configuration
  # change with no parameter change still shows up here
  data_sig <- vapply(il$data, function(x)
    paste(utils::capture.output(str(x, give.attr = FALSE, vec.len = 10)), collapse = "|"),
    character(1))

  list(map = map_sig, par = par_sig, data = data_sig)
}


#' Which blocks differ between two signatures
#'
#' @keywords internal
sweep_diff <- function(a, b) {
  cmp <- function(x, y) {
    keys <- union(names(x), names(y))
    keys[!vapply(keys, function(k) identical(x[[k]], y[[k]]), logical(1))]
  }
  list(map = cmp(a$map, b$map), par = cmp(a$par, b$par), data = cmp(a$data, b$data))
}

#' TRUE when two signatures are identical in every respect
#'
#' @keywords internal
sweep_identical <- function(a, b) {
  d <- sweep_diff(a, b)
  length(d$map) == 0 && length(d$par) == 0 && length(d$data) == 0
}


# ---------------------------------------------------------------------------
# Which configuration makes a spec live
#
# Several specs govern a parameter that only exists once some other option is
# switched on: the age-correlation specs need an at-age data source with a
# non-iid correlation, the selectivity process-error specs need a time-varying
# selectivity form. Outside that configuration the spec has nothing to map, and
# every one of its values builds the same model.
#
# Recording the configuration turns "this option looks dead" into a question with
# an answer. A spec listed here must be live in the configuration named; a spec
# still inert there is wired to nothing.
# ---------------------------------------------------------------------------

#' Extra arguments that put a spec's parameter into the model
#'
#' @param arg Spec argument name.
#' @param dims Dimension overrides, so the arrays match the test setup being swept.
#'
#' @return Named list of extra arguments for the same stage, or \code{NULL} when
#'   the spec needs no special configuration.
#'
#' @keywords internal
sweep_live_config <- function(arg, dims = list()) {
  d <- utils::modifyList(sweep_dims, dims)
  d <- utils::modifyList(d, sweep_live_dims(arg))
  aa_dim <- function(n_fleets) c(d$n_regions, d$n_yrs, d$n_seas, d$n_ages, d$n_sexes, n_fleets)
  agg_dim <- function(n_fleets) c(d$n_regions, d$n_yrs, d$n_seas, n_fleets)

  # An at-age data source replaces its aggregated counterpart rather than joining it,
  # so the aggregate is switched off wherever the at-age form is switched on. The
  # type is set sex-split because the default sums over sexes, which a model
  # with two sexes of observations cannot do.
  aa_stream <- function(prefix, n_fleets, corr_name, agg_use, se = FALSE, extra = list()) {
    out <- list()
    out[[paste0("Obs", prefix)]] <- array(100, dim = aa_dim(n_fleets))
    out[[paste0("Use", prefix)]] <- array(1, dim = aa_dim(n_fleets))
    if(se) out[[paste0("Obs", prefix, "_SE")]] <- array(0.2, dim = aa_dim(n_fleets))
    out[[paste0(prefix, "_Type")]] <- "spltRspltS"
    if(!is.null(agg_use)) out[[agg_use]] <- array(0, dim = agg_dim(n_fleets))
    # an iid correlation has no parameter to map, so the data source is switched on
    # with the simplest form that does have one
    if(!is.null(corr_name)) out[[corr_name]] <- "1dar1"
    utils::modifyList(out, extra)
  }

  n_f <- d$n_fish_fleets; n_s <- d$n_srv_fleets

  switch(arg,
    rho_catch_spec      = aa_stream("CatchAA", n_f, "AgeObsCorr_catch", "UseCatch"),
    rho_discard_spec    = aa_stream("DiscardAA", n_f, "AgeObsCorr_discard", "UseDiscard",
                                    extra = list(discard_units = rep("abd", n_f))),
    sigmaDAA_spec       = aa_stream("DiscardAA", n_f, NULL, "UseDiscard",
                                    extra = list(discard_units = rep("abd", n_f))),
    sigmaCAA_spec       = aa_stream("CatchAA", n_f, NULL, "UseCatch"),
    rho_srv_idx_spec    = aa_stream("SrvIdxAA", n_s, "AgeObsCorr_srv_idx", "UseSrvIdx", se = TRUE),
    sigmaSrvIdxAA_spec  = aa_stream("SrvIdxAA", n_s, NULL, "UseSrvIdx", se = TRUE),
    # Time-varying selectivity needs BOTH the process-error spec and the deviation
    # spec; supplying one alone passes the NULL guard and fails downstream.
    srvsel_pe_pars_spec  = list(cont_tv_srv_sel = paste0("iid_Fleet_", seq_len(n_s)),
                                srv_sel_devs_spec = rep("est_all", n_s)),
    srv_sel_devs_spec    = list(cont_tv_srv_sel = paste0("iid_Fleet_", seq_len(n_s)),
                                srvsel_pe_pars_spec = rep("est_all", n_s)),
    fishsel_pe_pars_spec = list(cont_tv_fish_sel = paste0("iid_Fleet_", seq_len(n_f)),
                                fish_sel_devs_spec = rep("est_all", n_f)),
    fish_sel_devs_spec   = list(cont_tv_fish_sel = paste0("iid_Fleet_", seq_len(n_f)),
                                fishsel_pe_pars_spec = rep("est_all", n_f)),
    # retention selectivity is fixed by default, so an est variant is refused
    # until the fixed form is switched off
    ret_fixed_sel_pars_spec = list(use_fixed_ret_sel = rep(0, n_f),
                                   ret_sel_model = paste0("logist1_Fleet_", seq_len(n_f))),
    # switching the fixed retention curve off makes the retained parameters
    # something to estimate, which its own spec then has to agree with
    ret_sel_devs_spec    = list(use_fixed_ret_sel = rep(0, n_f),
                                ret_sel_model = paste0("logist1_Fleet_", seq_len(n_f)),
                                ret_fixed_sel_pars_spec = rep("est_all", n_f),
                                cont_tv_ret_sel = paste0("iid_Fleet_", seq_len(n_f)),
                                retsel_pe_pars_spec = rep("est_all", n_f)),
    retsel_pe_pars_spec  = list(use_fixed_ret_sel = rep(0, n_f),
                                ret_sel_model = paste0("logist1_Fleet_", seq_len(n_f)),
                                ret_fixed_sel_pars_spec = rep("est_all", n_f),
                                cont_tv_ret_sel = paste0("iid_Fleet_", seq_len(n_f)),
                                ret_sel_devs_spec = rep("est_all", n_f)),
    NULL)
}


#' Dimension overrides a spec needs to have anything to map
#'
#' The population-specific data sources exist only in a model with more than one
#' population, so their specs are swept at two populations rather than at the
#' test setup's default of one.
#'
#' @keywords internal
sweep_live_dims <- function(arg) {
  if(grepl("_pop_spec$", arg)) list(n_pop = 2, n_regions = 2, natal_region = c(1, 2)) else list()
}


#' Settings on stages earlier than the spec's own that its configuration needs
#'
#' A population-specific data source needs a model that has populations, and a
#' model with more than one population only recruits under local density
#' dependence. That setting belongs to the recruitment stage, several stages
#' before the spec being swept.
#'
#' @keywords internal
sweep_live_other <- function(arg) {
  if(grepl("_pop_spec$", arg)) list(rec = list(rec_dd = "local")) else list()
}


# Extent of the dimension each abbreviation names, so a spec that shares only over
# dimensions this model has one of can be recognized as legitimately equal to
# est_all rather than reported as dead wiring.
sweep_abbrev_extent <- function(stage, dims = list()) {
  d <- utils::modifyList(sweep_dims, dims)
  n_f <- if(grepl("Srv", stage)) d$n_srv_fleets else d$n_fish_fleets
  c(
    p = d$n_pop,
    pop = d$n_pop,
    r = d$n_regions,
    y = d$n_yrs,
    seas = d$n_seas,
    s = d$n_sexes,
    x = d$n_sexes,
    f = n_f
  )
}

#' Whether a sharing spec collapses only dimensions this model has one of
#'
#' @keywords internal
sweep_spec_is_degenerate <- function(stage, value, dims = list()) {
  if(!grepl("^est_shared_", value)) return(FALSE)
  parts <- strsplit(sub("^est_shared_", "", value), "_")[[1]]
  ext <- sweep_abbrev_extent(stage, dims)
  # an abbreviation this table does not have (selectivity bins, say) cannot be
  # ruled degenerate, so the spec is treated as meaningful
  if(!all(parts %in% names(ext))) return(FALSE)
  all(ext[parts] < 2)
}
