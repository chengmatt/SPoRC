# Stage 1 of 3: model setup
#
#' Set Up A Projection Inside The Estimation Model
#'
#' Projects the population forward inside the model, so that \code{sdreport} gives standard
#' errors for projected spawning biomass, catch, fishing mortality, numbers at age and stock
#' status, with the uncertainty of the fit and of future recruitment in them. The projection
#' is the same one \code{\link{Do_Population_Projection}} runs after a fit, but fed the
#' model's own estimates rather than copies of them: the terminal year's numbers at age and
#' F, the demographics averaged over the terminal years, and mean recruitment. Projection
#' year 1 is the terminal assessment year run again, so year 2 is the first year after the
#' assessment and the year catch advice is read from.
#'
#' \code{\link{fit_model}} fits the model with the projection switched off, since the
#' projection adds nothing to the likelihood, then rebuilds it at the estimates with the
#' projection on.
#'
#' @param input_list List with \code{$data}, \code{$par} and \code{$map} from the earlier
#'   setup calls.
#' @param n_proj_yrs Integer. Number of years to project, counting the terminal assessment
#'   year, which is run again as projection year 1. 0 (default) leaves the projection off
#'   and the model as it was.
#' @param f_opt How fishing mortality is set in the projection years. \code{"terminal"}
#'   (default) fishes every year at the terminal year's F. \code{"input"} fishes at the F
#'   given in \code{f_input}. \code{"hcr"} runs a harvest control rule on reference points
#'   given in \code{f_ref_pt} and \code{b_ref_pt}, treated as known numbers. \code{"refpt"}
#'   runs the rule on reference points the model solves for itself from its own estimates,
#'   so the uncertainty in F40, B40 or MSY is part of the projection's standard errors.
#' @param f_input Array \code{[n_regions, n_proj_yrs]} of the F to fish at in each region
#'   and projection year, read under \code{f_opt = "input"}. Year 1 replays the terminal
#'   year, so column 1 is not read.
#' @param catch_input Retained catch to take in each region and projection year, an array
#'   \code{[n_regions, n_proj_yrs]}, or \code{[n_regions, n_proj_yrs, n_seas]} for catch
#'   by season, in the units of \code{proj_Catch}. The model finds the F that takes it, so
#'   the catch has no uncertainty and the uncertainty moves into that F and the projected
#'   biomass. A year left \code{NA} in every region is fished under \code{f_opt}, so known
#'   catch in the first years can be followed by the rule; a year is given in every region
#'   or in none. Year 1 replays the terminal year, so column 1 is not read.
#'   \code{F_constraint} applies only in the years left to \code{f_opt}. Default
#'   \code{NULL}, no catch given.
#' @param f_ref_pt The F reference point the rule reads under \code{f_opt = "hcr"}: the
#'   target F, one value per region, or \code{[n_regions, n_proj_yrs]} to vary it by year.
#'   A fixed number here, so its own uncertainty is not in the standard errors;
#'   \code{f_opt = "refpt"} is the way to include it.
#' @param b_ref_pt The biomass reference point the rule reads under \code{f_opt = "hcr"}:
#'   one value per population and region, or \code{[n_pop, n_regions, n_proj_yrs]} to vary
#'   it by year. A fixed number here, as \code{f_ref_pt} is.
#' @param hcr_function The harvest control rule: a function of \code{x}, the spawning biomass
#'   the rule reads, \code{frp}, the F reference point, and \code{brp}, the biomass reference
#'   point, returning the F to apply, as \code{\link{Do_Population_Projection}} takes it.
#'   Default \code{\link{hcr_threshold}}: the target F above \code{brp}, no fishing below
#'   0.05 times \code{brp}, and a straight line between. For another cutoff, pass
#'   \code{function(x, frp, brp) hcr_threshold(x, frp, brp, alpha = 0.1)}. A rule of your
#'   own runs inside the model, where the spawning biomass it
#'   is given is a quantity the model is differentiating, and such a quantity cannot be
#'   compared with \code{>} or \code{<}. Compare \code{ad_value(x)}, the plain number, to
#'   decide which part of the rule applies, and use \code{x} itself in the arithmetic that
#'   follows, as \code{hcr_threshold} does.
#' @param hcr_scope Which spawning biomass a region's rule reads. \code{"region"} (default),
#'   its own, against its own biomass reference point. \code{"global"}, the total over
#'   regions against the reference points summed over regions, so every region shares one
#'   stock status.
#' @param F_constraint A limit on the catch the rule's or input F would take.
#'   \code{\link{catch_cap}} lowers F wherever the catch would exceed a cap;
#'   \code{\link{catch_stability}} keeps the catch within a band of the previous year's,
#'   lowering or raising F; a list of them is applied in the order given, so put the cap
#'   last when nothing may exceed it. \code{NULL} (default) applies none.
#' @param refpt_type Under \code{f_opt = "refpt"}, \code{"single_region"} (default) or
#'   \code{"multi_region"}, named as in \code{\link{Get_Reference_Points}}.
#' @param refpt_what Under \code{f_opt = "refpt"}, which reference point the model solves,
#'   named as in \code{\link{Get_Reference_Points}}: \code{"SPR"} (default) or \code{"MSY"}
#'   for one region, and \code{"independent_SPR"}, \code{"independent_MSY"},
#'   \code{"global_SPR"}, \code{"global_MSY"} or \code{"local_MSY"} for several. The MSY
#'   options need a model fitted with a stock-recruit curve.
#' @param SPR_x The target spawning potential ratio for the SPR reference points, for
#'   example 0.4 for F40 and B40.
#' @param sex_ratio_f Array \code{[n_pop, n_regions]}. Female share of recruits the
#'   per-recruit calculation behind the reference points starts from. Default \code{NULL}
#'   reads it from the model: in a two-sex model the recruitment sex ratio averaged over the
#'   last \code{n_avg_yrs} years, so an estimated sex ratio and its uncertainty reach the
#'   reference points; in a one-sex model 0.5, since the model halves its spawning biomass.
#'   A value given here is used instead, for instance to match a post hoc
#'   \code{\link{Get_Reference_Points}} run.
#' @param is_discard_fleet Integer vector \code{[n_fish_fleets]}, 1 for a fleet whose catch
#'   is discarded rather than landed, so the MSY options do not count it as yield. Default
#'   all zero.
#' @param local_bh_msy_newton_steps Integer. Steps the local MSY solve takes for the
#'   equilibrium recruitment at each trial F. Default 6, as in
#'   \code{\link{Get_Reference_Points}}.
#' @param refpt_newton_steps Integer. The reference points are solved numerically first,
#'   then refined with this many Newton steps that the model differentiates through, which
#'   is what gives a reference point its standard error. Default 2, which is enough. 0
#'   leaves the reference points as fixed numbers with no standard error.
#' @param refpt_log_f_start Numeric. Where the numerical solve for a reference point starts,
#'   on the log scale. Default \code{log(0.1)}.
#' @param n_avg_yrs Integer. How many terminal years weight at age, maturity, natural
#'   mortality, sex ratio and straying are averaged over, for the projection and for the
#'   reference points alike. Default 1, the terminal year alone.
#' @param sel_opt Selectivity in the projection years. \code{"average"} (default) repeats
#'   its average over the last \code{n_avg_yrs}. \code{"propagate"} uses the years the
#'   model projected itself, which needs \code{n_proj_yrs_devs} in
#'   \code{\link{Setup_Mod_Dim}} and a time-varying selectivity, so its projected
#'   deviations and their uncertainty come along.
#' @param move_opt Movement in the projection years, with the same two choices as
#'   \code{sel_opt}: \code{"average"} (default) over the last \code{n_avg_yrs}, or
#'   \code{"propagate"} from the years the model projected itself.
#' @param growth_opt Weight at age in the projection years, when a growth model derives it
#'   (\code{waa_model = "wt_len"} in \code{\link{Setup_Mod_Biologicals}}). \code{"average"}
#'   (default) repeats its average over the last \code{n_avg_yrs}. \code{"propagate"} uses
#'   the years the model projected growth itself, which needs \code{n_proj_yrs_devs} in
#'   \code{\link{Setup_Mod_Dim}}, time-varying growth, and \code{growth_tv_type = "curve"},
#'   so the projected growth deviations and their uncertainty come along. Weight at age
#'   read from data has nothing to propagate and is always averaged.
#' @param rec_opt Recruitment in the projection years, always mean recruitment over the
#'   years \code{calc_rec_yrs} picks, times a deviation. \code{"propagate"} (default) uses
#'   the deviations the model estimated for the projection years, which needs
#'   \code{n_proj_yrs_devs} in \code{\link{Setup_Mod_Dim}}. With no data behind them they
#'   sit at zero, so the projection runs at the mean, but each has the recruitment process
#'   error as its standard error, so the projected standard errors widen as a forecast's
#'   should. \code{"mean"} applies no deviation: recruitment is the mean, and its
#'   uncertainty is only that of the mean, which is small. \code{"terminal"} applies the
#'   terminal year's deviation in every projection year, one value shared by all of them,
#'   so the projected cohorts rise or fall together rather than averaging out.
#' @param naa_opt The numbers at age state of a state-space model in the projection years.
#'   \code{"propagate"} (default) applies the state innovations the model estimated for the
#'   projection years, which needs \code{n_proj_yrs_devs} in \code{\link{Setup_Mod_Dim}};
#'   they sit at their forecast, zero under \code{"iid"}, and bring the state's process
#'   error into the standard errors. \code{"terminal"} advances the terminal year's state
#'   with no further innovations. Ignored by a model without a state.
#' @param calc_rec_st_yr Integer. First model year of the window mean recruitment is
#'   averaged over, as in \code{\link{Get_Reference_Points}}; the window is
#'   \code{calc_rec_st_yr:(n_years - rec_age)}. Default 1, the first model year.
#' @param rec_age Integer. Years left off the end of that window, the most recent years
#'   whose recruitment is still poorly known. Default 1.
#' @param calc_rec_yrs Integer vector of model years to average recruitment over instead
#'   of the window, for a mean that runs to the terminal year or skips years. Model year
#'   indices, 1 for the first year, not calendar years. Default \code{NULL}, the window.
#'
#' @return \code{input_list} with the settings stored in \code{$data} as \code{proj_*} fields.
#' @export Setup_Mod_Projection
#' @family Reference Points and Projections
#' @examples
#' \dontrun{
#' # six years under the threshold rule on post hoc F40 and B40, capped at 25 kt
#' input_list <- Setup_Mod_Projection(input_list, n_proj_yrs = 6, f_opt = "hcr",
#'                                    f_ref_pt = refpts$f_ref_pt, b_ref_pt = refpts$b_ref_pt,
#'                                    F_constraint = catch_cap(25))
#' # 20 kt taken in each of the next two years, then the same rule
#' catch <- array(NA, dim = c(1, 6))
#' catch[1, 2:3] <- 20
#' input_list <- Setup_Mod_Projection(input_list, n_proj_yrs = 6, f_opt = "hcr",
#'                                    f_ref_pt = refpts$f_ref_pt, b_ref_pt = refpts$b_ref_pt,
#'                                    catch_input = catch)
#' fit <- fit_model(input_list$data, input_list$par, input_list$map, do_sdrep = TRUE)
#' fit$sd_rep # log_proj_SSB, log_proj_Catch, log_proj_F, log_proj_Total_Biom
#' }
Setup_Mod_Projection <- function(input_list,
                                 n_proj_yrs = 0,
                                 f_opt = "terminal",
                                 f_input = NULL,
                                 catch_input = NULL,
                                 f_ref_pt = NULL,
                                 b_ref_pt = NULL,
                                 hcr_function = hcr_threshold,
                                 hcr_scope = "region",
                                 F_constraint = NULL,
                                 refpt_type = "single_region",
                                 refpt_what = "SPR",
                                 SPR_x = NULL,
                                 sex_ratio_f = NULL,
                                 is_discard_fleet = NULL,
                                 local_bh_msy_newton_steps = 6,
                                 refpt_newton_steps = 2,
                                 refpt_log_f_start = log(0.1),
                                 n_avg_yrs = 1,
                                 sel_opt = "average",
                                 move_opt = "average",
                                 growth_opt = "average",
                                 rec_opt = "propagate",
                                 naa_opt = "propagate",
                                 calc_rec_st_yr = 1,
                                 rec_age = 1,
                                 calc_rec_yrs = NULL) {

  messages_list <<- character(0)
  data <- input_list$data
  n_regions <- data$n_regions
  n_pop <- data$n_pop
  n_yrs <- length(data$years)

  # Options -----------------------------------------------------------------

  f_code <- convert_to_numeric(f_opt, list(terminal = 1, input = 2, hcr = 3, refpt = 4))
  scope_code <- convert_to_numeric(hcr_scope, list(region = 1, global = 2))
  refpt_type_code <- convert_to_numeric(refpt_type, list(single_region = 1, multi_region = 2))
  refpt_what_code <- convert_to_numeric(refpt_what, list(SPR = 1, MSY = 2, independent_SPR = 3, independent_MSY = 4, global_SPR = 5, global_MSY = 6, local_MSY = 7))
  sel_code <- convert_to_numeric(sel_opt, list(average = 1, propagate = 2))
  move_code <- convert_to_numeric(move_opt, list(average = 1, propagate = 2))
  growth_code <- convert_to_numeric(growth_opt, list(average = 1, propagate = 2))
  rec_code <- convert_to_numeric(rec_opt, list(propagate = 1, mean = 2, terminal = 3))
  naa_code <- convert_to_numeric(naa_opt, list(propagate = 1, terminal = 2))

  # Input Validation --------------------------------------------------------

  if(!is.numeric(n_proj_yrs) || length(n_proj_yrs) != 1 || n_proj_yrs < 0 || n_proj_yrs != round(n_proj_yrs)) stop("n_proj_yrs must be a single non-negative integer")

  if(n_proj_yrs > 0) {

    # growth can only be propagated where a growth model builds it over the projection years
    if(growth_code == 2) {
      if(is.null(data$derive_waa) || data$derive_waa != 1) stop("growth_opt = 'propagate' needs weight at age derived by a growth model (waa_model = 'wt_len' in Setup_Mod_Biologicals); weight at age read from data has nothing to propagate.")
      if(!is.null(data$growth_tv_type) && data$growth_tv_type == 1) stop("growth_opt = 'propagate' needs growth_tv_type = 'curve' in Setup_Mod_Biologicals, since cohort growth is not built over the projection years.")
    }

    # a process propagated from the fit needs its projected deviations in the fit
    n_proj_devs <- if(is.null(data$n_proj_yrs_devs)) 0 else data$n_proj_yrs_devs
    propagated <- c(if(rec_code == 1) "rec_opt", if(sel_code == 2) "sel_opt", if(move_code == 2) "move_opt", if(growth_code == 2) "growth_opt")
    if(length(propagated) > 0 && n_proj_devs < n_proj_yrs - 1) stop("A projection of ", n_proj_yrs, " years propagates ", paste(propagated, collapse = ", "), " from the fit, which needs deviations for ", n_proj_yrs - 1, " years past the fit, but Setup_Mod_Dim was given n_proj_yrs_devs = ", n_proj_devs, ". Set n_proj_yrs_devs to at least ", n_proj_yrs - 1, ", or project those at the mean, the terminal year or an average instead.")
    has_state <- !is.null(input_list$par$ln_NAA) && !is.null(data$n_est_naa_re) && data$n_est_naa_re > 0
    if(naa_code == 1 && has_state && dim(input_list$par$ln_NAA)[3] < n_yrs + n_proj_yrs - 1) stop("naa_opt = 'propagate' needs the numbers at age state over ", n_proj_yrs - 1, " years past the fit, but it has ", dim(input_list$par$ln_NAA)[3] - n_yrs, ". Set n_proj_yrs_devs in Setup_Mod_Dim to at least ", n_proj_yrs - 1, ", or naa_opt = 'terminal'.")
    if(n_avg_yrs < 1 || n_avg_yrs > n_yrs) stop("n_avg_yrs must be between 1 and the number of model years")

    # F by year
    if(f_code == 2) {
      if(is.null(f_input)) stop("f_opt = 'input' needs f_input, an array [n_regions, n_proj_yrs]")
      if(!identical(as.integer(dim(f_input)), as.integer(c(n_regions, n_proj_yrs)))) stop("f_input should be dimensioned [", n_regions, ", ", n_proj_yrs, "], but is [", paste(dim(f_input), collapse = ", "), "]")
      if(any(!is.finite(f_input)) || any(f_input < 0)) stop("f_input holds negative or non-finite values")
    }

    # catch by year, solved for the F that takes it; a year left NA in every region is fished under f_opt
    if(!is.null(catch_input)) {
      if(is.null(dim(catch_input)) || !length(dim(catch_input)) %in% c(2, 3)) stop("catch_input should be an array [n_regions, n_proj_yrs], or [n_regions, n_proj_yrs, n_seas] for catch by season")
      want <- if(length(dim(catch_input)) == 3) c(n_regions, n_proj_yrs, data$n_seas) else c(n_regions, n_proj_yrs)
      if(!identical(as.integer(dim(catch_input)), as.integer(want))) stop("catch_input should be dimensioned [", paste(want, collapse = ", "), "], but is [", paste(dim(catch_input), collapse = ", "), "]")
      if(any(is.nan(catch_input)) || any(!is.na(catch_input) & (catch_input < 0 | is.infinite(catch_input)))) stop("catch_input holds negative or non-finite catch. Use NA to leave a year to f_opt, and 0 for no fishing.")
      n_set <- apply(!is.na(catch_input), 2, sum)[-1] # year 1 replays the terminal year and is not read
      n_cell <- length(catch_input) / n_proj_yrs
      part <- which(n_set > 0 & n_set < n_cell) + 1
      if(length(part) > 0) stop("catch_input is only partly given in projection year(s) ", paste(part, collapse = ", "), ". Give a year's catch in every region", if(length(want) == 3) " and season" else "", ", or leave the whole year NA to fish it under f_opt.")
      if(all(n_set == 0)) stop("catch_input gives no catch from projection year 2 on, the first year after the assessment.")
    }

    # the rule and what it reads
    if(f_code == 3) {
      if(is.null(f_ref_pt) || is.null(b_ref_pt)) stop("f_opt = 'hcr' needs f_ref_pt [n_regions] and b_ref_pt [n_pop, n_regions], or either by projection year as well")
      if(!length(f_ref_pt) %in% c(n_regions, n_regions * n_proj_yrs)) stop("f_ref_pt should have ", n_regions, " values, one per region, or ", n_regions * n_proj_yrs, " by region and projection year, not ", length(f_ref_pt))
      if(!length(b_ref_pt) %in% c(n_pop * n_regions, n_pop * n_regions * n_proj_yrs)) stop("b_ref_pt should have ", n_pop * n_regions, " values, by population and region, or ", n_pop * n_regions * n_proj_yrs, " by projection year as well, not ", length(b_ref_pt))
      if(any(!is.finite(f_ref_pt)) || any(f_ref_pt < 0) || any(!is.finite(b_ref_pt)) || any(b_ref_pt <= 0)) stop("f_ref_pt should be finite and non-negative, b_ref_pt finite and positive")
      if(!is.null(hcr_function)) {
        if(!is.function(hcr_function)) stop("hcr_function should be a function of x, frp and brp")
        if(!all(c("x", "frp", "brp") %in% names(formals(hcr_function)))) stop("hcr_function must take arguments x, frp and brp")
      }
    }

    # the rule on reference points solved inside the model
    if(f_code == 4) {
      if(!is.null(f_ref_pt) || !is.null(b_ref_pt)) stop("f_opt = 'refpt' solves the reference points inside the model, so f_ref_pt and b_ref_pt are not read. Leave them NULL, or give them under f_opt = 'hcr'.")
      if(refpt_type_code == 1 && n_regions > 1) stop("refpt_type = 'single_region' with ", n_regions, " regions. Use 'multi_region' with one of independent_SPR, independent_MSY, global_SPR, global_MSY or local_MSY.")
      if(refpt_type_code == 1 && !refpt_what_code %in% c(1, 2)) stop("refpt_type = 'single_region' takes refpt_what 'SPR' or 'MSY'")
      if(refpt_type_code == 2 && !refpt_what_code %in% 3:7) stop("refpt_type = 'multi_region' takes refpt_what 'independent_SPR', 'independent_MSY', 'global_SPR', 'global_MSY' or 'local_MSY'")
      if(refpt_what_code %in% c(1, 3, 5) && (is.null(SPR_x) || SPR_x <= 0 || SPR_x >= 1)) stop("The SPR reference points need SPR_x in (0, 1)")
      check_msy_rec_model(refpt_what, if(is.null(data$rec_model)) 1 else data$rec_model, if(is.null(data$sr_penalty)) 0 else data$sr_penalty)
      if(refpt_what_code == 6 && n_pop > 1) stop("global_MSY needs one population. Use independent_SPR, independent_MSY, global_SPR or local_MSY.")
      if(!is.null(sex_ratio_f) && length(sex_ratio_f) != n_pop * n_regions) stop("sex_ratio_f should be [n_pop, n_regions]")
      if(!is.null(is_discard_fleet) && length(is_discard_fleet) != data$n_fish_fleets) stop("is_discard_fleet should have one value per fishery fleet")
      if(refpt_newton_steps < 0 || local_bh_msy_newton_steps < 1) stop("refpt_newton_steps must be 0 or more and local_bh_msy_newton_steps 1 or more")
      if(!is.null(hcr_function)) {
        if(!is.function(hcr_function)) stop("hcr_function should be a function of x, frp and brp")
        if(!all(c("x", "frp", "brp") %in% names(formals(hcr_function)))) stop("hcr_function must take arguments x, frp and brp")
      }
    }

    # constraints on F
    if(!is.null(F_constraint)) {
      if(is.function(F_constraint)) F_constraint <- list(F_constraint)
      if(!is.list(F_constraint) || length(F_constraint) == 0) stop("F_constraint should be a function (F, state), or a list of them applied in order. See catch_cap() and catch_stability().")
      for(i in seq_along(F_constraint)) if(!is.function(F_constraint[[i]])) stop("F_constraint[[", i, "]] is not a function. Each constraint is a function (F, state); see catch_cap() and catch_stability().")
    }

    # mean recruitment years, as Get_Reference_Points reads them
    if(is.null(calc_rec_yrs)) calc_rec_yrs <- calc_rec_st_yr:(n_yrs - rec_age)
    if(length(calc_rec_yrs) == 0 || any(!calc_rec_yrs %in% 1:n_yrs)) stop("Recruitment years (calc_rec_yrs, or calc_rec_st_yr:(n_years - rec_age)) must be model years between 1 and ", n_yrs, ", not calendar years.")

  } # end checks that apply when the projection is on

  # Defaults for Unsupplied Inputs ------------------------------------------

  if(is.null(f_input)) f_input <- array(0, dim = c(n_regions, max(n_proj_yrs, 1)))
  if(is.null(f_ref_pt)) f_ref_pt <- array(0, dim = c(n_regions, max(n_proj_yrs, 1)))
  if(is.null(b_ref_pt)) b_ref_pt <- array(1, dim = c(n_pop, n_regions, max(n_proj_yrs, 1)))
  if(is.null(calc_rec_yrs)) calc_rec_yrs <- 1:n_yrs
  if(is.null(is_discard_fleet)) is_discard_fleet <- array(0, dim = data$n_fish_fleets)
  if(is.null(SPR_x)) SPR_x <- 0

  # year 1 replays the terminal year, so its catch is never read
  if(n_proj_yrs == 0) catch_input <- NULL
  if(!is.null(catch_input)) {
    catch_input <- array(as.numeric(catch_input), dim = dim(catch_input))
    if(length(dim(catch_input)) == 3) catch_input[,1,] <- NA else catch_input[,1] <- NA
  }

  # one value per region is repeated over the projection years
  f_ref_pt <- array(f_ref_pt, dim = c(n_regions, max(n_proj_yrs, 1)))
  b_ref_pt <- array(b_ref_pt, dim = c(n_pop, n_regions, max(n_proj_yrs, 1)))

  # Populate Data List ------------------------------------------------------

  input_list$data$proj_n_yrs <- n_proj_yrs
  input_list$data$proj_f_opt <- f_code
  input_list$data$proj_f_input <- f_input
  input_list$data["proj_catch_input"] <- list(catch_input) # NULL has to be stored, not dropped
  input_list$data$proj_f_ref_pt <- f_ref_pt
  input_list$data$proj_b_ref_pt <- b_ref_pt
  input_list$data["proj_hcr_function"] <- list(hcr_function) # NULL has to be stored, not dropped
  input_list$data$proj_hcr_scope <- scope_code
  input_list$data["proj_F_constraint"] <- list(F_constraint) # NULL has to be stored, not dropped
  input_list$data$proj_refpt_type <- refpt_type_code
  input_list$data$proj_refpt_what <- refpt_what_code
  input_list$data$proj_spr_x <- SPR_x
  input_list$data["proj_sex_ratio_f"] <- list(if(is.null(sex_ratio_f)) NULL else array(sex_ratio_f, dim = c(n_pop, n_regions))) # NULL reads the model's own
  input_list$data$proj_is_discard_fleet <- is_discard_fleet
  input_list$data$proj_local_bh_msy_newton_steps <- local_bh_msy_newton_steps
  input_list$data$proj_refpt_newton_steps <- refpt_newton_steps
  input_list$data$proj_refpt_log_f_start <- refpt_log_f_start
  input_list$data$proj_n_avg_yrs <- n_avg_yrs
  input_list$data$proj_sel_opt <- sel_code
  input_list$data$proj_move_opt <- move_code
  input_list$data$proj_growth_opt <- growth_code
  input_list$data$proj_rec_opt <- rec_code
  input_list$data$proj_naa_opt <- naa_code
  input_list$data$proj_calc_rec_yrs <- calc_rec_yrs

  # Print Messages ----------------------------------------------------------

  if(n_proj_yrs == 0) collect_message("Projection inside the model is off")
  else {
    collect_message("Projection inside the model: ", n_proj_yrs, " years, fishing mortality from ", f_opt)
    if(f_code %in% c(3, 4)) collect_message("Rule: ", if(is.null(hcr_function) || identical(hcr_function, hcr_threshold)) "the threshold rule, hcr_threshold" else "given", ", reading ", hcr_scope, " spawning biomass")
    if(f_code == 4) collect_message("Reference point ", refpt_what, " (", refpt_type, ") solved inside the model, refined with ", refpt_newton_steps, " Newton steps; female share of recruits ",
                                    if(is.null(sex_ratio_f)) "from the model" else "given")
    if(!is.null(catch_input)) {
      catch_yrs <- which(apply(!is.na(catch_input), 2, any))
      collect_message("Catch given in projection year(s) ", paste(catch_yrs, collapse = ", "), ", with the F that takes it solved inside the model",
                      if(length(catch_yrs) < n_proj_yrs - 1) paste0("; the other years fished under ", f_opt) else "")
    }
    if(!is.null(F_constraint)) collect_message("Constraints on F: ", length(F_constraint), " applied in order")
    collect_message("Demographics averaged over the last ", n_avg_yrs, " year(s); selectivity ", sel_opt, ", movement ", move_opt, ", growth ", growth_opt)
    collect_message("Projected recruitment: mean over model years ", min(calc_rec_yrs), " to ", max(calc_rec_yrs),
                    c(" times the fit's projected deviations", " alone", " times the terminal year's deviation")[rec_code])
    if(has_state) collect_message("Numbers at age state: ", c("the fit's projected innovations", "the terminal state advanced with none")[naa_code])
  }

  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

} # end Setup_Mod_Projection
