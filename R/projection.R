# Stage 3 of 3: post fit
#
# Forward projection off a fitted model under a specified catch or fishing mortality, with optional
# stochastic recruitment. Short horizon advice, not the closed loop simulation in sim_closed_loop.R.

#' Do Population Projections
#'
#' Projects the population forward under a recruitment and a fishing mortality
#' scenario, starting from the terminal assessment year and advancing numbers at
#' age over \code{[population x region x year x season x age x sex]} through
#' recruitment, movement, mortality, ageing and a harvest control rule. Recruitment
#' is generated annually and spread over seasons by \code{rec_seas_prop}.
#'
#' @param n_proj_yrs Integer. Number of projection years.
#' @param n_pop Integer. Number of populations, which may exceed the number of
#'   regions under natal homing.
#' @param n_regions Integer. Number of spatial regions.
#' @param n_ages Integer. Number of age classes including the plus group.
#' @param n_sexes Integer. Number of sexes.
#' @param sexratio Array `[n_pop, n_regions, n_proj_yrs, n_sexes]` allocating
#'   projected recruits by sex.
#' @param n_fish_fleets Integer. Number of fishing fleets.
#' @param do_recruits_move Integer (0 or 1). Whether age-1 recruits move. Default 0.
#' @param rec_seas_prop Array `[n_pop, n_seas]` of the share of annual recruitment
#'   entering in each season, summing to 1 for each population.
#' @param recruitment Array `[n_pop, n_regions, n_yrs]` of historical recruitment,
#'   used to condition the stochastic options.
#' @param terminal_NAA Array `[n_pop, n_regions, n_seas, n_ages, n_sexes]` of fished
#'   numbers at age in the terminal assessment year.
#' @param terminal_NAA0 As \code{terminal_NAA}, unfished.
#' @param terminal_F Array `[n_regions, n_seas, n_fish_fleets]`. Sets F in
#'   projection year 1 and the seasonal F ratios used in later years.
#' @param natmort Natural mortality, a rate per year in each season, either
#'   `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]` or the same without
#'   the season dim, which is expanded across seasons. Scaled internally by season
#'   duration.
#' @param WAA Array `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]` of
#'   weight-at-age, used for spawning biomass.
#' @param WAA_fish As \code{WAA} with a trailing `n_fish_fleets` dim, used for catch
#'   biomass.
#' @param MatAA Array dimensioned like \code{WAA}, maturity-at-age.
#' @param fish_sel Array `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes,
#'   n_fish_fleets]` of fishery selectivity-at-age.
#' @param Movement Array `[n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages,
#'   n_sexes]` of seasonal movement transition matrices.
#' @param Mrate Array dimensioned like \code{Movement} holding the generator rather
#'   than the realized fractions. Only read when `move_timing` is 1 or 2, so `NULL`
#'   (default) is valid under `move_timing = 0`.
#' @param expm_nsub Integer controlling how the matrix exponential is evaluated
#'   under `move_timing = 2`: `0` uses `Matrix::expm`, `n >= 1` uses `n` implicit
#'   backward Euler substeps. See [mat_exp()].
#' @param move_timing When movement happens relative to mortality within a season.
#'   `0` (default) moves then kills, `1` kills then moves, and `2` runs the two
#'   continuously, which also switches catch at age to the spatial Baranov form on
#'   season-integrated abundance. Must match the timing the reference points were
#'   derived under.
#' @param sgl_seas_spawning_movement Array `[n_pop, n_regions, n_regions,
#'   n_proj_yrs, n_ages, n_sexes]` redistributing fish to natal grounds before SSB.
#'   Only read when `n_seas = 1` and `n_pop > 1`, so `NULL` (default) is otherwise
#'   valid.
#' @param stray_rate Array `[n_pop, n_proj_yrs]` used when accumulating effective
#'   SSB across populations. Only read when `n_pop > 1`, so `NULL` (default) is
#'   otherwise valid.
#' @param f_ref_pt Array `[n_regions, n_proj_yrs]` of the fishing mortality
#'   reference point or fixed input F, depending on `fmort_opt`.
#' @param b_ref_pt Array `[n_pop, n_regions, n_proj_yrs]` of the biomass reference
#'   point used in the control rule.
#' @param HCR_function Harvest control rule taking `x` (SSB), `frp` and `brp`. A
#'   rule that also declares a `state` argument (or `...`) is handed this year's
#'   population as a named list, holding `y`, `r`, `NAA`, `SSB`, `Total_Biom` and
#'   `Catch`, so it can be written on more than spawning biomass. Rules without it
#'   are called as before and no state is assembled.
#' @param recruitment_opt Recruitment scenario: `"inv_gauss"`, `"mean_rec"`,
#'   `"zero"`, or `"bh_rec"`.
#' @param fmort_opt Fishing mortality scenario: `"HCR"`, `"HCR_global"`, `"Input"`,
#'   or `"Catch"`, which solves each year's F so realized catch matches
#'   `catch_input` and leaves every other quantity untouched.
#' @param catch_input Catch targets in biomass, read under `fmort_opt = "Catch"`.
#'   Either `[n_regions, n_proj_yrs]` of annual targets or `[n_regions,
#'   n_proj_yrs, n_seas]` of seasonal ones, and the shape decides what is solved.
#'   `catch_input[r, y]` is the catch removed in projection year `y`, indexed as
#'   `proj_Catch` is rather than with the one year lag `f_ref_pt` uses. Targets are
#'   totals over populations and fleets, and over seasons in the annual case; a
#'   target of 0 sets `F = 0` without a solve.
#'
#'   Annual targets solve one F per region, split over seasons at the terminal year
#'   shares. Seasonal targets solve an F per region and season, with the fleet split
#'   within a season still at terminal year ratios, so fleet-specific targets are
#'   not supported either way. A season the terminal year did not fish has no fleet
#'   split to inherit and can take no catch, which is an error.
#'
#'   A year set to `NA` falls back to `catch_fallback_opt`, which is the usual shape
#'   of catch advice. `NA` and `0` are different things. A year must be all target
#'   or all `NA` across regions and seasons; a partly specified year is an error.
#'   Column 1 is read only when `catch_terminal_yr = TRUE`.
#' @param catch_fallback_opt Which rule sets F in the years `catch_input` leaves
#'   `NA`: `"HCR"`, `"HCR_global"` or `"Input"`. Defaults to `"HCR"` under
#'   `fmort_opt = "Catch"` and to `fmort_opt` itself otherwise, where it is unused.
#'   Mind the indexing when mixing the two: `catch_input[r, y]` is the catch taken
#'   in year `y`, while `f_ref_pt[r, y]` sets F in year `y + 1`.
#' @param catch_terminal_yr Logical. Whether projection year 1, which replays the
#'   terminal assessment year, is solved against its catch target rather than fished
#'   at `terminal_F`. Default `FALSE`. `TRUE` suits the common case where the
#'   terminal year is not yet complete, but it overrides the F the assessment
#'   estimated and so changes the numbers entering year 2. With `n_seas > 1` the
#'   terminal year takes all its seasons from `terminal_NAA`, so only the last
#'   season's F feeds year 2.
#' @param catch_f_max Upper bound on the F searched under `fmort_opt = "Catch"`.
#'   Default 5. An unreachable target caps F here, undershoots, and warns with the
#'   regions named.
#' @param catch_tol Relative catch tolerance for the F solver. Default 1e-6.
#' @param catch_max_iter Maximum solver iterations per projection year. Default 100.
#' @param t_spawn Fraction of the spawning season elapsed before spawning.
#' @param srr_opt Named list of inputs for deterministic stock-recruit recruitment
#'   under `recruitment_opt = "bh_rec"` or `"ricker_rec"`, passed straight to
#'   \code{\link{Get_Det_Recruitment}} and holding every argument that function
#'   needs. Formerly `bh_rec_opt`. The arrays are \code{R0} \code{[n_pop]}, \code{h}
#'   and \code{rec_region_prop} \code{[n_pop, n_regions]}, \code{rec_seas_prop}
#'   \code{[n_pop, n_seas]}, \code{SSB} \code{[n_pop, n_regions, n_yrs]},
#'   \code{WAA}, \code{MatAA} and \code{natmort} \code{[n_pop, n_regions, n_seas,
#'   n_ages]} (natmort also accepted without the season dim), \code{Movement}
#'   \code{[n_pop, n_regions, n_regions, n_seas, n_ages]},
#'   \code{sgl_seas_spawning_movement} \code{[n_pop, n_regions, n_regions,
#'   n_ages]}, \code{stray_rate} \code{[n_pop]}, \code{init_F} \code{[n_regions,
#'   n_seas, n_fish_fleets]}, \code{fish_sel} and \code{ret_sel} \code{[n_pop,
#'   n_regions, n_seas, n_ages, n_fish_fleets]}, \code{dmr} \code{[n_regions,
#'   n_seas, n_fish_fleets]} and \code{sex_ratio_f} \code{[n_pop, n_regions]}. The
#'   scalars are \code{rec_dd}, \code{rec_lag}, \code{n_pop}, \code{n_regions},
#'   \code{n_ages}, \code{n_seas}, \code{spawn_seas}, \code{seasdur},
#'   \code{t_spawn} and \code{do_recruits_move}. Spawning biomass is built
#'   internally by appending projected SSB to \code{srr_opt$SSB}.
#'
#'   \code{srr_opt$rec_lag = 1} computes each year's recruitment up front from the
#'   prior year's SSB, as \code{"inv_gauss"} and \code{"mean_rec"} do.
#'   \code{srr_opt$rec_lag = 0} computes it from the year's own SSB once
#'   \code{spawn_seas} is reached, and inserts recruits no earlier than that season,
#'   so \code{rec_seas_prop} must be zero before it. Reference points and the
#'   seasonal SBPR are unaffected either way.
#' @param bh_rec_opt Deprecated former name of \code{srr_opt}; supplying it warns
#'   and forwards, and supplying both is an error.
#' @param n_seas Integer. Number of seasons. Default 1.
#' @param seasdur Numeric vector `[n_seas]` of season durations as fractions of a
#'   year.
#' @param spawn_seas Integer spawning season index.
#' @param natal_region Integer vector `[n_pop]` of each population's natal region.
#'   Only read when `n_pop > 1`, so `NULL` (default) is otherwise valid.
#' @param dmr Array \code{[n_regions, n_seas, n_fish_fleets]} of discard mortality
#'   rate. Default \code{0}, which with \code{ret_sel = 1} means a fleet discards
#'   nothing.
#' @param ret_sel Array \code{[n_pop, n_regions, n_proj_yrs, n_seas, n_ages,
#'   n_sexes, n_fish_fleets]} of retention selectivity-at-age. Default \code{1},
#'   full retention.
#' @param rec_devs Optional array \code{[n_pop, n_regions, n_proj_yrs]} of
#'   multiplicative deviations applied to whatever recruitment
#'   \code{recruitment_opt} produces, so a deterministic option becomes stochastic
#'   under deviations the caller draws. \code{NULL} (default) leaves recruitment as
#'   the option gives it. Year 1 generates no recruitment, so its slice is never
#'   read. Drawing outside is what lets replicates share recruitment across
#'   management procedures, and what lets the projection be differentiated with
#'   respect to the rule with the deviations kept fixed.
#'
#' @return A named list of projected quantities. Year index 1 is the terminal
#'   assessment year replayed, so year 2 is the first projected year and catch
#'   advice for terminal year + 1 is read from index 2. \code{proj_NAA},
#'   \code{proj_NAA0}, \code{proj_F} and \code{proj_F_seas} fill their trailing
#'   \code{n_proj_yrs + 1} year slot; \code{proj_ZAA}, \code{proj_ret_FAA} and
#'   \code{proj_disc_FAA} leave it at 0.
#'
#' \describe{
#'   \item{\code{proj_F}}{`[n_regions, n_proj_yrs + 1]`. Annual F by region, summed
#'     over seasons and fleets. The trailing column holds the F the rule or input
#'     would apply in the year after the projection, and stays 0 under
#'     `fmort_opt = "Catch"`.}
#'   \item{\code{proj_F_seas}}{`[n_regions, n_proj_yrs + 1, n_seas]`. The same F by
#'     season, so `rowSums(proj_F_seas[, y, ])` recovers `proj_F[, y]`. The only
#'     place the answer lives under seasonal catch targets.}
#'   \item{\code{proj_ret_FAA}}{`[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages,
#'     n_sexes, n_fish_fleets]`. Retained fishing mortality-at-age, the component
#'     that lands catch.}
#'   \item{\code{proj_disc_FAA}}{Dimensioned as `proj_ret_FAA`. Discard fishing
#'     mortality-at-age, set by `ret_sel` and `dmr`. Total F at age is the sum of
#'     the two.}
#'   \item{\code{proj_Catch}}{`[n_pop, n_regions, n_proj_yrs, n_seas,
#'     n_fish_fleets]`. Retained catch in biomass, built from `proj_ret_FAA`, and
#'     the quantity `catch_input` is matched against.}
#'   \item{\code{proj_SSB}}{`[n_pop, n_regions, n_proj_yrs]`. Female spawning
#'     biomass, accumulated in `spawn_seas` with the `t_spawn` correction. Halved
#'     when `n_sexes = 1`.}
#'   \item{\code{proj_eff_SSB}}{`[n_pop, n_proj_yrs]`. Effective spawning biomass at
#'     each population's natal region, with cross-population terms scaled by
#'     `stray_rate`. Equal to SSB summed across regions when `n_pop = 1`.}
#'   \item{\code{proj_Total_Biom}}{`[n_pop, n_regions, n_proj_yrs]`. Total biomass
#'     over all ages and sexes, at the same point in the season as `proj_SSB` and
#'     on the estimation model's definition, so the series continues without a
#'     discontinuity at the terminal year.}
#'   \item{\code{proj_Dynamic_SSB0}}{`[n_pop, n_regions, n_proj_yrs]`. Spawning
#'     biomass under the same realized recruitment but no fishing, for dynamic
#'     depletion.}
#'   \item{\code{proj_NAA}}{`[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages,
#'     n_sexes]`. Fished numbers at age at the start of each season, before that
#'     season's mortality and ageing. Movement has already been applied under
#'     `move_timing = 0` and not under 1 or 2. The trailing slot holds the numbers
#'     passed into the year after the projection.}
#'   \item{\code{proj_NAA0}}{Dimensioned as `proj_NAA`, decremented by natural
#'     mortality alone.}
#'   \item{\code{proj_ZAA}}{Dimensioned as `proj_NAA`. Total mortality-at-age for
#'     the season: natural mortality scaled by season duration plus retained and
#'     discard F summed over fleets.}
#'   \item{\code{proj_catch_resid}}{Shaped like `catch_input`. Relative miss on each
#'     target, `(realized - target) / target`, and `NA` for years with no target.
#'     Should be at or below `catch_tol` wherever the solve converged, and is worth
#'     checking directly rather than relying on warnings.}
#' }
#'
#' @details
#' A projection year generates and allocates recruitment, builds F at age from the
#' annual F, the terminal year's seasonal ratios and selectivity, moves fish each
#' season, applies within-season mortality and ages the survivors at the end of the
#' final season, computes spawning biomass in \code{spawn_seas} with a mid-season
#' correction (after spawning movement under natal homing with one season), takes
#' catch by the Baranov equation, and sets next year's F from the control rule or
#' input.
#'
#' Under \code{srr_opt$rec_lag == 0} the recruitment and spawning steps are
#' reordered within \code{spawn_seas}: movement runs first, spawning biomass is
#' computed from the survivors alone, that SSB generates this year's recruitment,
#' and the recruits are inserted immediately before mortality and ageing. Year 1
#' holds the terminal state forward with no recruitment event.
#'
#' Under \code{fmort_opt = "Catch"} the F step moves to the front of the following
#' year, since the F that lands a target depends on that year's own numbers at age
#' rather than the previous year's spawning biomass. The year is run at trial F
#' values until realized catch matches the target, then run once more at the
#' accepted F and committed; no demographic input is modified. Regions are solved
#' jointly, because movement (and, under \code{move_timing = 2}, the
#' season-integrated abundance) makes each region's catch depend on every other
#' region's F. Seasonal targets are swept forward one season at a time, which is
#' exact because a season's catch depends only on the F in that season and earlier
#' ones. The catch solved against is retained catch, so a fleet that discards
#' exerts more total F than the target implies.
#'
#' Spawning biomass is multiplied by 0.5 when \code{n_sexes = 1}, and movement is
#' skipped when \code{n_regions = 1}.
#'
#' @section Differentiating through the projection:
#'
#' The projection uses RTMB's replacement operators, so it can be taped with
#' \code{\link[RTMB]{MakeTape}} or \code{\link[RTMB]{MakeADFun}} and optimized with
#' an exact gradient, giving an F schedule solved against an objective rather than
#' scanned over a grid.
#'
#' Two options are refused on AD types, since neither has a derivative and both
#' would otherwise return a wrong gradient rather than an error:
#' \code{recruitment_opt = "inv_gauss"} draws at random, and
#' \code{fmort_opt = "Catch"} inverts the target numerically. Tape under
#' \code{"mean_rec"}, \code{"bh_rec"} or \code{"ricker_rec"} with
#' \code{fmort_opt = "Input"}. A control rule that branches on stock status stops
#' on its own, since comparing an AD spawning biomass raises an error inside RTMB,
#' so optimizing through a rule means writing a smooth one.
#'
#' @export Do_Population_Projection
#' @family Reference Points and Projections
#' @import abind
Do_Population_Projection <- function(
  n_proj_yrs = 2,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  sexratio,
  n_fish_fleets,
  do_recruits_move = 0,
  recruitment,
  terminal_NAA,
  terminal_NAA0,
  terminal_F,
  dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
  natmort,
  natal_region = NULL,
  WAA,
  WAA_fish,
  MatAA,
  fish_sel,
  ret_sel = array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
  Movement,
  sgl_seas_spawning_movement = NULL,
  stray_rate = NULL,
  f_ref_pt = NULL,
  b_ref_pt = NULL,
  HCR_function = NULL,
  recruitment_opt = "inv_gauss",
  fmort_opt = 'HCR',
  catch_input = NULL,
  catch_fallback_opt = if(fmort_opt == "Catch") "HCR" else fmort_opt,
  catch_terminal_yr = FALSE,
  catch_f_max = 5,
  catch_tol = 1e-6,
  catch_max_iter = 100,
  t_spawn,
  srr_opt = NULL,
  bh_rec_opt = NULL,
  n_seas = 1,
  seasdur = rep(1 / n_seas, n_seas),
  spawn_seas = 1,
  rec_seas_prop = {
    rec_seas_prop = array(0, dim = c(n_pop, n_seas))
    rec_seas_prop[] <- 1 / n_seas
    rec_seas_prop
  }, Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0,
  rec_devs = NULL
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # srr_opt was bh_rec_opt when Beverton-Holt was the only stock-recruit curve.
  # It now has either curve, so the bh_ prefix is wrong rather than redundant.
  if(!is.null(bh_rec_opt)) {
    if(!is.null(srr_opt)) stop("Supply either srr_opt or the deprecated bh_rec_opt, not both.")
    warning("'bh_rec_opt' is deprecated and will be removed; use 'srr_opt'. It now holds the Ricker as well, so the bh_ prefix no longer describes it.", call. = FALSE)
    srr_opt <- bh_rec_opt
  }

  # Error Checking ----------------------------------------------------------

  if(!recruitment_opt %in% c("inv_gauss", "mean_rec", "zero", "bh_rec", "ricker_rec")) stop("Recruitment options are not specified correctly! Should be inv_gauss, mean_rec, zero, bh_rec, or ricker_rec")
  if(!fmort_opt %in% c("HCR", "Input", "HCR_global", "Catch")) stop("Fishing Mortality options are not specified correctly! Should be HCR, Input, HCR_global, or Catch")
  if(!catch_fallback_opt %in% c("HCR", "Input", "HCR_global")) stop("Catch fallback options are not specified correctly! Should be HCR, Input, or HCR_global")

  if(!is.null(rec_devs)) {
    want <- c(n_pop, n_regions, n_proj_yrs)
    if(!identical(as.integer(dim(rec_devs)), as.integer(want))) stop(paste0("rec_devs should be dimensioned [", paste(want, collapse = ", "), "], but is [", paste(dim(rec_devs), collapse = ", "), "]."))
    if(any(!is.finite(rec_devs)) || any(rec_devs < 0)) stop("rec_devs holds negative or non-finite values. They multiply recruitment, so they should be positive.")
  }

  # Taping the projection turns its inputs into AD types. Two options cannot be
  # differentiated through, and both would give a silently wrong gradient rather
  # than an error, so they are refused here instead.
  if(any(vapply(list(f_ref_pt, catch_input, terminal_NAA, terminal_F, fish_sel, natmort, WAA),
                inherits, logical(1), "advector"))) {
    if(recruitment_opt == "inv_gauss") stop("recruitment_opt = 'inv_gauss' draws recruitment at random, which has no derivative. Tape the projection under 'mean_rec', 'bh_rec' or 'ricker_rec'.")
    if(fmort_opt == "Catch") stop("fmort_opt = 'Catch' inverts the catch target with a numerical solve, and the F it returns carries no derivative. Differentiating through it needs implicit differentiation; tape the projection under fmort_opt = 'Input' instead.")
  }

  # Backwards compatibility for seasonal natural mortality ...
  natmort <- expand_natmort_seasons(natmort, n_seas)
  if(!is.null(srr_opt) && !is.null(srr_opt$natmort))
    srr_opt$natmort <- expand_natmort_seasons(srr_opt$natmort, n_seas, seas_dim = 3, n_dim = 4)


  # Set up for catch stuff
  # which years in the projection are driven by a catch target
  catch_yr_targeted <- rep(FALSE, n_proj_yrs)
  # whether targets have a season dimension, set from catch_input's shape below
  catch_seasonal <- FALSE
  # Setup fmort rule
  fmort_rule <- if(fmort_opt == "Catch") catch_fallback_opt else fmort_opt

  if(fmort_opt == "Catch") {

    # check dimensions and input
    if(is.null(catch_input)) stop("fmort_opt = 'Catch' requires catch_input, an array [n_regions, n_proj_yrs] or [n_regions, n_proj_yrs, n_seas] of catch targets.")
    if(is.null(dim(catch_input)) || length(dim(catch_input)) == 1) { # accept a vector and reshape it
      if(length(catch_input) != n_regions * n_proj_yrs) stop(paste0("catch_input has ", length(catch_input), " values but n_regions * n_proj_yrs = ", n_regions * n_proj_yrs, "."))
      catch_input <- array(catch_input, dim = c(n_regions, n_proj_yrs))
    }

    # error checking for seasonal catches
    catch_seasonal <- length(dim(catch_input)) == 3
    want <- if(catch_seasonal) c(n_regions, n_proj_yrs, n_seas) else c(n_regions, n_proj_yrs)
    if(!identical(as.integer(dim(catch_input)), as.integer(want))) stop(paste0("catch_input should be dimensioned [", paste(want, collapse = ", "), "], but is [", paste(dim(catch_input), collapse = ", "), "]."))
    # NA means no target that year, which is different from a target of 0 (no fishing). Anything else has to be a usable catch.
    if(any(catch_input < 0 | is.nan(catch_input) | is.infinite(catch_input), na.rm = TRUE)) stop("catch_input holds negative or non-finite catch targets. Use NA to leave a year to the fallback rule, and 0 to ask for no fishing.")
    if(catch_f_max <= 0) stop("catch_f_max should be a positive upper bound on the F searched.")

    # Check to see if any missing values mid-season
    n_set <- apply(!is.na(catch_input), 2, sum)
    n_cell <- length(catch_input) / n_proj_yrs
    part <- which(n_set > 0 & n_set < n_cell)
    if(length(part) > 0) stop(paste0("catch_input is only partly specified in projection year(s) ", paste(part, collapse = ", "), ". Give every region", if(catch_seasonal) " and season" else "", " in a year a target, or set them all NA to leave that year to catch_fallback_opt."))
    catch_yr_targeted <- n_set == n_cell
    if(!catch_terminal_yr) catch_yr_targeted[1] <- FALSE # year 1 replays the terminal assessment year
    if(!any(catch_yr_targeted)) stop("fmort_opt = 'Catch' but no projection year has a catch target. Note that year 1 is only solved when catch_terminal_yr = TRUE.")

    # Year 1 always falls back to terminal_F rather than to catch_fallback_opt,
    # so only years 2 onward can call on the fallback rule. Check its inputs are present in the fxn
    if(n_proj_yrs > 1 && any(!catch_yr_targeted[2:n_proj_yrs])) {
      if(is.null(f_ref_pt)) stop(paste0("catch_input leaves projection year(s) ", paste(which(!catch_yr_targeted[2:n_proj_yrs]) + 1, collapse = ", "), " to catch_fallback_opt = '", catch_fallback_opt, "', which needs f_ref_pt."))
      if(catch_fallback_opt %in% c("HCR", "HCR_global") && (is.null(HCR_function) || is.null(b_ref_pt))) stop(paste0("catch_fallback_opt = '", catch_fallback_opt, "' needs HCR_function and b_ref_pt for the projection years catch_input leaves NA."))
    }

    # A season the terminal year did not fish has no fleet selectivity split to use, so no F can be apportioned into it and no catch can be taken there.
    if(catch_seasonal) {
      asked <- apply(array(catch_input[,catch_yr_targeted,, drop = FALSE], dim = c(n_regions, sum(catch_yr_targeted), n_seas)), c(1,3), max)
      dead <- which(apply(terminal_F, c(1,2), sum) == 0 & asked > 0, arr.ind = TRUE)
      if(nrow(dead) > 0) stop(paste0("catch_input asks for catch in region ", dead[1,1], ", season ", dead[1,2], ", but terminal_F is 0 there, so there is no fleet split to apportion F with."))
    }
  }

  # error checking for bh_opt
  if(recruitment_opt %in% c("bh_rec", "ricker_rec")) {
    required_fields <- c("rec_dd", "rec_lag", "R0", "h", "rec_region_prop",
                         "WAA", "MatAA", "natmort", "SSB", "Movement",
                         "sex_ratio_f", "stray_rate", "fish_sel", "ret_sel", "dmr", "init_F")
    diff <- setdiff(required_fields, names(srr_opt)) # find difference
    if(length(diff) > 0) stop(paste("srr_opt is missing the following required fields:", paste(diff)))
  }

  # Define Containers -------------------------------------------------------

  # Get splits by season and region and split up terminal F
  seas_share <- array(0, dim = c(n_regions, n_seas))
  fratio_fleet <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  for(r in 1:n_regions) {
    for(seas in 1:n_seas) {
      seas_tot <- sum(terminal_F[r,seas,])
      seas_share[r,seas] <- seas_tot / sum(terminal_F[r,,])
      # A season the terminal year did not fish has no fleet split to inherit,
      # and gets a zero share anyway, so leave the split at zero rather than 0/0.
      if(seas_tot > 0) for(f in 1:n_fish_fleets) fratio_fleet[r,seas,f] <- terminal_F[r,seas,f] / seas_tot
    } # end seas loop
  } # end r loop

  proj_NAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_NAA0 <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_ZAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_tot_FAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_ret_FAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_disc_FAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_CAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_DAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_Catch <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets))
  proj_SSB <- array(0, dim = c(n_pop, n_regions, n_proj_yrs))
  proj_eff_SSB <- array(0, dim = c(n_pop, n_proj_yrs))
  proj_Total_Biom <- array(0, dim = c(n_pop, n_regions, n_proj_yrs))
  proj_Dynamic_SSB0 <- array(0, dim = c(n_pop, n_regions, n_proj_yrs))
  proj_F <- array(0, dim = c(n_regions, n_proj_yrs + 1))
  proj_F_seas <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_seas))
  proj_catch_resid <- array(NA_real_, dim = if(catch_seasonal) c(n_regions, n_proj_yrs, n_seas) else c(n_regions, n_proj_yrs)) # relative catch miss, only filled under fmort_opt = 'Catch'
  tmp_rec <- NULL # this year's recruitment, generated below and handed to the season loop

  # Start Projection --------------------------------------------------------

  # Input terminal year assessment at age
  proj_NAA[,,1,,,] <- terminal_NAA
  proj_NAA0[,,1,,,] <- terminal_NAA0

  # the two stock-recruit options share the per-recruit calculation, the lag and the apportionment,
  # and differ only in the curve Get_Det_Recruitment evaluates. absent srr_opt means Beverton-Holt
  if(!is.null(srr_opt)) srr_opt$rec_model <- if(recruitment_opt == "ricker_rec") 2 else 1

  # Flag for age-0 stock-recruit recruitment
  age0_rec <- recruitment_opt %in% c("bh_rec", "ricker_rec") && !is.null(srr_opt) && srr_opt$rec_lag == 0

  # Arguments for run_proj_yr used below
  proj_args <- list(
    n_pop = n_pop,
    n_regions = n_regions,
    n_ages = n_ages,
    n_sexes = n_sexes,
    n_seas = n_seas,
    n_fish_fleets = n_fish_fleets,
    fratio_fleet = fratio_fleet,
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    dmr = dmr,
    natmort = natmort,
    seasdur = seasdur,
    Movement = Movement,
    Mrate = Mrate,
    move_timing = move_timing,
    expm_nsub = expm_nsub,
    do_recruits_move = do_recruits_move,
    WAA = WAA,
    MatAA = MatAA,
    WAA_fish = WAA_fish,
    t_spawn = t_spawn,
    spawn_seas = spawn_seas,
    sgl_seas_spawning_movement = sgl_seas_spawning_movement,
    natal_region = natal_region,
    stray_rate = stray_rate,
    sexratio = sexratio,
    rec_seas_prop = rec_seas_prop,
    age0_rec = age0_rec,
    srr_opt = srr_opt,
    rec_devs = rec_devs
  )

  for(y in 1:n_proj_yrs) {

    # use terminal F in the first year (subsequent years use F derived from reference points and HCR)
    if(y == 1) proj_F[,y] <- rowSums(terminal_F)

    # Recruitment Processes (rec_lag != 0, or non-BH recruitment) -------------

    # For age0_rec, recruitment for the year is instead generated inline once
    # spawn_seas is reached within the season loop below.
    if(y > 1 && !age0_rec) {

      # Get annual recruitment
      tmp_rec <- switch(recruitment_opt,

                        "inv_gauss" = { # if inverse gaussian
                          sapply(1:n_regions, function(r)
                            sapply(1:n_pop, function(p)
                              rinvgauss_rec(1, recruitment[p, r, ])
                            )
                          )
                        },

                        "mean_rec" = { # if mean recruitment
                          sapply(1:n_regions, function(r)
                            sapply(1:n_pop, function(p)
                              mean(recruitment[p, r, ])
                            )
                          )
                        },

                        "zero" = { # if zero recruitment
                          array(0, dim = c(n_pop, n_regions))
                        },

                        "bh_rec" = , # both stock-recruit options land here
                        "ricker_rec" = { # Beverton-Holt or Ricker, per srr_opt$rec_model
                          Get_Det_Recruitment(recruitment_model = srr_opt$rec_model,
                                              rec_dd = srr_opt$rec_dd,
                                              n_pop = n_pop,
                                              sgl_seas_spawning_movement = srr_opt$sgl_seas_spawning_movement,
                                              natal_region = natal_region,
                                              y = y + dim(srr_opt$SSB)[3],
                                              rec_lag = srr_opt$rec_lag,
                                              R0 = srr_opt$R0,
                                              rec_region_prop = srr_opt$rec_region_prop,
                                              rec_seas_prop = rec_seas_prop,
                                              h = srr_opt$h,
                                              n_regions = n_regions,
                                              n_ages = n_ages,
                                              WAA = srr_opt$WAA,
                                              MatAA = srr_opt$MatAA,
                                              n_seas = n_seas,
                                              seasdur = seasdur,
                                              spawn_seas = spawn_seas,
                                              natmort = srr_opt$natmort,
                                              SSB_vals = bind_proj_SSB(srr_opt$SSB, proj_SSB),
                                              Movement = srr_opt$Movement,
                                              # SSB0 behind the stock recruit curve has to use the same movement
                                              # sequencing as the projection itself, so forward both of these.
                                              Mrate = srr_opt$Mrate,
                                              stray_rate = srr_opt$stray_rate,
                                              do_recruits_move = do_recruits_move,
                                              t_spawn = t_spawn,
                                              sexratio_f = srr_opt$sex_ratio_f,
                                              init_F = srr_opt$init_F,
                                              n_fish_fleets = n_fish_fleets,
                                              fish_sel = srr_opt$fish_sel,
                                              ret_sel = srr_opt$ret_sel,
                                              dmr = srr_opt$dmr,
                                              move_timing = move_timing,
                                              expm_nsub = expm_nsub)
                        }
      )

      # coerce into array
      tmp_rec <- array(tmp_rec, dim = c(n_pop, n_regions))
      if(!is.null(rec_devs)) tmp_rec <- tmp_rec * array(rec_devs[,,y], dim = c(n_pop, n_regions))

      # Apply recruitment to projected proj_NAA
      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          tmp <- tmp_rec[p,r] * sexratio[p,r,y,] * rec_seas_prop[p,1]
          proj_NAA[p,r,y,1,1,] <- proj_NAA0[p,r,y,1,1,]  <- tmp
        } # end r loop
      } # end p loop

    } # if y > 1

    # Grab state arguments
    state <- list(proj_NAA = proj_NAA,
                 proj_NAA0 = proj_NAA0,
                 proj_ZAA = proj_ZAA,
                 proj_ret_FAA = proj_ret_FAA,
                 proj_disc_FAA = proj_disc_FAA,
                 proj_tot_FAA = proj_tot_FAA,
                 proj_CAA = proj_CAA,
                 proj_DAA = proj_DAA,
                 proj_Catch = proj_Catch,
                 proj_SSB = proj_SSB,
                 proj_Dynamic_SSB0 = proj_Dynamic_SSB0,
                 proj_eff_SSB = proj_eff_SSB,
                 proj_Total_Biom = proj_Total_Biom)

    # Solve This Year's F By Region And Season -------------------------------
    # Solve for catch here since need full year abundance to know what catch is
    if(fmort_opt == 'Catch' && catch_yr_targeted[y]) { # note that thius overwrites the terminal F value provided if catch_terminal_yr is set TRUE

      # Solve for catch to F
      catch_solve <- solve_proj_year_F(y = y,
                                       target = if(catch_seasonal) array(catch_input[,y,], dim = c(n_regions, n_seas)) else catch_input[,y],
                                       seasonal = catch_seasonal,
                                       seas_share = seas_share,
                                       f_start = if(y > 1) array(proj_F_seas[,y - 1,], dim = c(n_regions, n_seas)) else apply(terminal_F, c(1,2), sum),
                                       state = state,
                                       tmp_rec = tmp_rec,
                                       proj_args = proj_args,
                                       catch_f_max = catch_f_max,
                                       catch_tol = catch_tol,
                                       catch_max_iter = catch_max_iter)

      F_y <- catch_solve$F_y
      proj_F[,y] <- rowSums(F_y) # annual total
      if(catch_seasonal) proj_catch_resid[,y,] <- catch_solve$resid else proj_catch_resid[,y] <- catch_solve$resid

    } else {

      # Distribute annual F to seasonal rates using seasonal splits determined above
      F_y <- array(proj_F[,y] * seas_share, dim = c(n_regions, n_seas))

    } # end if catch target

    proj_F_seas[,y,] <- F_y

    # Run the projection with the new F determined
    state <- do.call(run_proj_year, c(list(y = y, F_y = F_y, tmp_rec = tmp_rec), state, proj_args))

    proj_NAA <- state$proj_NAA
    proj_NAA0 <- state$proj_NAA0
    proj_ZAA <- state$proj_ZAA
    proj_ret_FAA <- state$proj_ret_FAA
    proj_disc_FAA <- state$proj_disc_FAA
    proj_tot_FAA <- state$proj_tot_FAA
    proj_CAA <- state$proj_CAA
    proj_Catch <- state$proj_Catch
    proj_SSB <- state$proj_SSB
    proj_Dynamic_SSB0 <- state$proj_Dynamic_SSB0
    proj_eff_SSB <- state$proj_eff_SSB
    proj_Total_Biom <- state$proj_Total_Biom


    # compute F for next year. fmort_rule is fmort_opt itself except under Catch, where it is
    # catch_fallback_opt and only runs when next year needs it
    if(fmort_opt != 'Catch' || (y + 1 <= n_proj_yrs && !catch_yr_targeted[y + 1])) {

      # A rule that declares a state argument is handed this year's numbers and
      # biomass, so a policy can read more than spawning biomass alone: mean
      # weight, age structure, last year's catch. Rules without one are called
      # exactly as before, so the state is only assembled when it is wanted.
      hcr_state <- if(fmort_rule %in% c("HCR", "HCR_global") &&
                      any(c("state", "...") %in% names(formals(HCR_function)))) {
        list(y = y,
             NAA = array(proj_NAA[,,y,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
             SSB = array(proj_SSB[,,y], dim = c(n_pop, n_regions)),
             Total_Biom = array(proj_Total_Biom[,,y], dim = c(n_pop, n_regions)),
             Catch = array(proj_Catch[,,y,,], dim = c(n_pop, n_regions, n_seas, n_fish_fleets)))
      } else NULL

      for(r in 1:n_regions) {

        # Project F using HCR and reference points -----------------------------------------------------
        if(fmort_rule == 'HCR') {
          if(is.null(hcr_state))
            proj_F[r,y + 1] <- HCR_function(x = sum(proj_SSB[,r,y]),
                                          frp = f_ref_pt[r,y],
                                          brp = sum(b_ref_pt[,r,y]))
          else
            proj_F[r,y + 1] <- HCR_function(x = sum(proj_SSB[,r,y]),
                                          frp = f_ref_pt[r,y],
                                          brp = sum(b_ref_pt[,r,y]),
                                          state = c(hcr_state, list(r = r)))
        }

        if(fmort_rule == 'HCR_global') {
          if(is.null(hcr_state))
            proj_F[r,y + 1] <- HCR_function(x = sum(proj_SSB[,,y]),
                                          frp = f_ref_pt[r,y],
                                          brp = sum(b_ref_pt[,,y]))
          else
            proj_F[r,y + 1] <- HCR_function(x = sum(proj_SSB[,,y]),
                                          frp = f_ref_pt[r,y],
                                          brp = sum(b_ref_pt[,,y]),
                                          state = c(hcr_state, list(r = r)))
        }

        # Project F using User Inputs ---------------------------------------------
        if(fmort_rule == 'Input') proj_F[r,y + 1] <- f_ref_pt[r,y]

      } # end r loop
    } # end if the year needs an F rule

  } # end y loop

  # The year loop never reaches n_proj_yrs + 1, but the HCR and Input rules leave
  # an F there, so give it the same seasonal split as the rest of proj_F_seas.
  proj_F_seas[,n_proj_yrs + 1,] <- array(proj_F[,n_proj_yrs + 1] * seas_share, dim = c(n_regions, n_seas))

  return(list(proj_F = proj_F,
              proj_ret_FAA = proj_ret_FAA,
              proj_disc_FAA = proj_disc_FAA,
              proj_Catch = proj_Catch,
              proj_CAA = proj_CAA,
              proj_SSB = proj_SSB,
              proj_eff_SSB = proj_eff_SSB,
              proj_Total_Biom = proj_Total_Biom,
              proj_Dynamic_SSB0 = proj_Dynamic_SSB0,
              proj_NAA = proj_NAA,
              proj_NAA0 = proj_NAA0,
              proj_ZAA = proj_ZAA,
              proj_F_seas = proj_F_seas,
              proj_catch_resid = proj_catch_resid)
  )

} # end function




# Catch Targeted Projection Helpers -------------------------------------------
#
# fmort_opt = "Catch" inverts a catch target back to fishing mortality, which has no closed form.
# run_proj_year() replays the year at trial F values, as a pure function so no trial leaks out.
#
# Call order, outermost first, once per projection year:
#
#   Do_Population_Projection()
#    +- solve_proj_year_F()        splits the year into blocks to solve
#        +- solve_proj_F_catch()   solves ONE block: one F per region
#            +- proj_catch_at_F()      what catch does a trial F give?
#            |   +- build_proj_F()         assembles the trial F matrix
#            |   +- run_proj_year()        replays the season loop
#            +- proj_target_catch()    reduces that catch to what the target is on
#            +- proj_log_catch_resid() the residual nleqslv is handed
#
# Every trial F matrix is built the same way, in build_proj_F():
#
#   F_y[r, seas] = F_base[r, seas] + F_reg[r] * seas_profile[r, seas]
#
#   F_reg        the one number per region being solved for
#   seas_profile how that number is split across seasons
#   F_base       F already settled and kept fixed
#
# Annual targets are one block for the whole year: seas_profile is the terminal year's seasonal
# shares, and the target is read against catch summed over seasons. Seasonal targets are one block
# per season, swept forward, which is exact because a season's catch depends only on the F in that
# season and earlier ones. Within a block, regions solve together, since movement makes each
# region's catch depend on every other region's F: one free region bisects, several use nleqslv.

#' Join Historical And Projected Spawning Biomass
#'
#' The stock-recruit curve reads spawning biomass from one array spanning assessment
#' and projection years. \code{abind()} drops the AD class, which would take
#' recruitment off the tape with no error raised, so the two are copied into a
#' container here instead.
#'
#' @param hist Array \code{[n_pop, n_regions, n_hist_yrs]} of assessment spawning biomass.
#' @param proj Array \code{[n_pop, n_regions, n_proj_yrs]} of projected spawning biomass.
#' @return Array \code{[n_pop, n_regions, n_hist_yrs + n_proj_yrs]}.
#' @keywords internal
#' @noRd
bind_proj_SSB <- function(hist, proj) {

  "[<-" <- RTMB::ADoverload("[<-")

  n_hist <- dim(hist)[3]
  n_proj <- dim(proj)[3]
  out <- array(0, dim = c(dim(hist)[1], dim(hist)[2], n_hist + n_proj))
  out[,,1:n_hist] <- hist
  out[,,(n_hist + 1):(n_hist + n_proj)] <- proj

  return(out)
}


#' Run One Projection Year At A Given Fishing Mortality
#'
#' Advances the projection through every season of year \code{y} at the fishing
#' mortality supplied, returning the updated state. Split out of
#' \code{\link{Do_Population_Projection}} so that a catch target can be solved
#' for by replaying the year, which requires the year to be reproducible from
#' its arguments alone.
#'
#' @param y Integer. Projection year to run.
#' @param F_y Numeric matrix \code{[n_regions, n_seas]}. Total fishing mortality
#'   by region and season, before the fleet split in \code{fratio_fleet}.
#' @param tmp_rec Numeric array \code{[n_pop, n_regions]} or \code{NULL}. This
#'   year's recruitment when it is already known. Ignored (and regenerated
#'   internally) when \code{age0_rec} is \code{TRUE}.
#' @param fratio_fleet Array \code{[n_regions, n_seas, n_fish_fleets]}. Fleet
#'   split of F within a season, summing to 1 across fleets, or all 0 for a
#'   season the terminal year did not fish.
#' @param age0_rec Logical. Whether recruitment is age-0 Beverton-Holt, in which
#'   case it is generated inside the season loop from this year's own SSB.
#' @param proj_NAA,proj_NAA0,proj_ZAA,proj_ret_FAA,proj_disc_FAA,proj_tot_FAA,proj_CAA,proj_DAA,proj_Catch,proj_SSB,proj_Dynamic_SSB0,proj_eff_SSB
#'   The projection arrays, as built in \code{Do_Population_Projection}.
#' @param n_pop,n_regions,n_ages,n_sexes,n_seas,n_fish_fleets,fish_sel,ret_sel,dmr,natmort,seasdur,Movement,Mrate,move_timing,expm_nsub,do_recruits_move,WAA,MatAA,WAA_fish,t_spawn,spawn_seas,sgl_seas_spawning_movement,natal_region,stray_rate,sexratio,rec_seas_prop,srr_opt,rec_devs
#'   Static projection inputs, documented in \code{\link{Do_Population_Projection}}.
#'
#' @return A named list holding the same eleven arrays, advanced through year
#'   \code{y}.
#'
#' @keywords internal
#' @noRd
run_proj_year <- function(y,
                          F_y,
                          tmp_rec,
                          proj_NAA, proj_NAA0, proj_ZAA,
                          proj_ret_FAA, proj_disc_FAA, proj_tot_FAA,
                          proj_CAA, proj_DAA, proj_Catch,
                          proj_SSB, proj_Dynamic_SSB0, proj_eff_SSB, proj_Total_Biom,
                          n_pop, n_regions, n_ages, n_sexes, n_seas, n_fish_fleets,
                          fratio_fleet, fish_sel, ret_sel, dmr, natmort, seasdur,
                          Movement, Mrate, move_timing, do_recruits_move,
                          WAA, MatAA, WAA_fish, t_spawn, spawn_seas,
                          sgl_seas_spawning_movement, natal_region, stray_rate,
                          sexratio, rec_seas_prop, age0_rec, srr_opt,
                          expm_nsub = 0, rec_devs = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  for(seas in 1:n_seas) {

    # insert seasonal recruits already known from earlier this year.
    # under age0_rec spawn_seas generates and inserts its own share below
    if(y > 1 && (if(age0_rec) seas > spawn_seas else seas > 1)) {
      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            proj_NAA[p,r,y,seas,1,s]  = proj_NAA[p,r,y,seas,1,s]  + tmp_rec[p,r] * rec_seas_prop[p,seas] * sexratio[p,r,y,s]
            proj_NAA0[p,r,y,seas,1,s] = proj_NAA0[p,r,y,seas,1,s] + tmp_rec[p,r] * rec_seas_prop[p,seas] * sexratio[p,r,y,s]
          } # end s loop
        } # end r loop
      } # end p loop
    } # end if

    # Construct Mortality Processes -------------------------------------------
    for(r in 1:n_regions) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          for(f in 1:n_fish_fleets) {
            # get fishing mortality at age
            for(p in 1:n_pop) {
              proj_ret_FAA[p,r,y,seas,a,s,f] <- F_y[r,seas] * fratio_fleet[r,seas,f] * fish_sel[p,r,y,seas,a,s,f] * ret_sel[p,r,y,seas,a,s,f] # retained F
              proj_disc_FAA[p,r,y,seas,a,s,f] <- F_y[r,seas] * fratio_fleet[r,seas,f] * fish_sel[p,r,y,seas,a,s,f] * (1 - ret_sel[p,r,y,seas,a,s,f]) * dmr[r,seas,f] # discarded F
              proj_tot_FAA[p,r,y,seas,a,s,f] <- proj_ret_FAA[p,r,y,seas,a,s,f] + proj_disc_FAA[p,r,y,seas,a,s,f] # total F
            } # end p loop
          } # end f loop

          # Get Total Mortality at Age
          for(p in 1:n_pop) {
            proj_ZAA[p,r,y,seas,a,s] <- (natmort[p,r,y,seas,a,s] * seasdur[seas]) + sum(proj_tot_FAA[p,r,y,seas,a,s,])
          }

        } # end s loop
      } # end a loop
    }

    # Movement Processes ------------------------------------------------------
    # Only apply movement if more than 1 region, or if y > 1 (because terminal proj_NAA already has movement applied).
    # Under move_timing 1 and 2 movement is deferred to the mortality/ageing step below.
    if(n_regions > 1 && y > 1 && move_timing == 0) {
      for(p in 1:n_pop) {
        # Recruits don't move
        if(do_recruits_move == 0) {
          # Apply movement after ageing processes - start movement at age 2
          for(a in 2:n_ages) for(s in 1:n_sexes) proj_NAA[p,,y,seas,a,s] = t(proj_NAA[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # fished
          for(a in 2:n_ages) for(s in 1:n_sexes) proj_NAA0[p,,y,seas,a,s] = t(proj_NAA0[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # unfished
        } # end if recruits don't move
        # Recruits move here
        if(do_recruits_move == 1) {
          for(a in 1:n_ages) for(s in 1:n_sexes) proj_NAA[p,,y,seas,a,s] = t(proj_NAA[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # fished
          for(a in 1:n_ages) for(s in 1:n_sexes) proj_NAA0[p,,y,seas,a,s] = t(proj_NAA0[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # unfished
        }
      } # end p loop
    } # only compute if spatial

    # Derive Biomass + Recruitment (age0_rec only) ------------------------------
    # SSB is fully determined by the survivors here, so generate this year's recruitment from
    # it and insert the spawn_seas share before mortality and ageing run below
    if(age0_rec && seas == spawn_seas) {

      biom <- derive_proj_biom(y, seas, proj_NAA, proj_NAA0, WAA, MatAA, proj_ZAA, natmort, t_spawn, seasdur,
                              n_seas, n_pop, n_regions, n_ages, n_sexes,
                              sgl_seas_spawning_movement, natal_region, stray_rate,
                              Movement, Mrate, move_timing, do_recruits_move, expm_nsub = expm_nsub)
      proj_SSB[,, y] <- biom$SSB_y
      proj_Dynamic_SSB0[,,y] <- biom$Dynamic_SSB0_y
      proj_eff_SSB[,y] <- biom$eff_SSB_y
      proj_Total_Biom[,,y] <- biom$Total_Biom_y

      if(y > 1) {

        tmp_rec <- Get_Det_Recruitment(recruitment_model = srr_opt$rec_model,
                                       rec_dd = srr_opt$rec_dd,
                                       n_pop = n_pop,
                                       sgl_seas_spawning_movement = srr_opt$sgl_seas_spawning_movement,
                                       natal_region = natal_region,
                                       y = y + dim(srr_opt$SSB)[3],
                                       rec_lag = srr_opt$rec_lag,
                                       R0 = srr_opt$R0,
                                       rec_region_prop = srr_opt$rec_region_prop,
                                       rec_seas_prop = rec_seas_prop,
                                       h = srr_opt$h,
                                       n_regions = n_regions,
                                       n_ages = n_ages,
                                       WAA = srr_opt$WAA,
                                       MatAA = srr_opt$MatAA,
                                       n_seas = n_seas,
                                       seasdur = seasdur,
                                       spawn_seas = spawn_seas,
                                       natmort = srr_opt$natmort,
                                       SSB_vals = bind_proj_SSB(srr_opt$SSB, proj_SSB),
                                       Movement = srr_opt$Movement,
                                       # SSB0 behind the stock recruit curve has to use the same movement
                                       # sequencing as the projection itself, so forward both of these.
                                       Mrate = srr_opt$Mrate,
                                       stray_rate = srr_opt$stray_rate,
                                       do_recruits_move = do_recruits_move,
                                       t_spawn = t_spawn,
                                       sexratio_f = srr_opt$sex_ratio_f,
                                       init_F = srr_opt$init_F,
                                       n_fish_fleets = n_fish_fleets,
                                       fish_sel = srr_opt$fish_sel,
                                       ret_sel = srr_opt$ret_sel,
                                       dmr = srr_opt$dmr,
                                       move_timing = move_timing,
                                       expm_nsub = expm_nsub)

        tmp_rec <- array(tmp_rec, dim = c(n_pop, n_regions))
        if(!is.null(rec_devs)) tmp_rec <- tmp_rec * array(rec_devs[,,y], dim = c(n_pop, n_regions))

        for(p in 1:n_pop) {
          for(r in 1:n_regions) {
            proj_NAA[p,r,y,spawn_seas,1,]  <- proj_NAA[p,r,y,spawn_seas,1,]  + tmp_rec[p,r] * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,]
            proj_NAA0[p,r,y,spawn_seas,1,] <- proj_NAA0[p,r,y,spawn_seas,1,] + tmp_rec[p,r] * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,]
          } # end r loop
        } # end p loop

        # recruits just inserted missed this season's movement step, which had to run before
        # SSB was knowable. catch age 1 up when recruits are supposed to move from birth

        # Only needed under move_timing == 0; under timings 1 and 2 these recruits are
        # picked up by the end-of-season transition below.
        if(do_recruits_move == 1 && n_regions > 1 && move_timing == 0) {
          for(p in 1:n_pop) {
            for(s in 1:n_sexes) proj_NAA[p,,y,seas,1,s] = t(proj_NAA[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
            for(s in 1:n_sexes) proj_NAA0[p,,y,seas,1,s] = t(proj_NAA0[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
          } # end p loop
        }

      } # end if y > 1

    } # end if age0_rec && seas == spawn_seas

    # Movement (timing 1 and 2), Mortality and Ageing --------------------------
    # Post-season state at every age, before the ageing shift. Under move_timing == 0
    # movement was applied above so this reduces to the original elementwise survival.
    if(move_timing == 0 || n_regions == 1) {
      pstep_NAA <- array(proj_NAA[,,y,seas,1:n_ages,] * exp(-proj_ZAA[,,y,seas,1:n_ages,]),
                         dim = c(n_pop, n_regions, n_ages, n_sexes))
      pstep_NAA0 <- array(proj_NAA0[,,y,seas,1:n_ages,] * exp(-natmort[,,y,seas,1:n_ages,] * seasdur[seas]),
                          dim = c(n_pop, n_regions, n_ages, n_sexes))
    } else {

      pstep_NAA <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
      pstep_NAA0 <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))

      # Advance fish throughout the season
      for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
        moves <- (do_recruits_move == 1 || a > 1)
        Mv <- if(moves) Movement[p,,,y,seas,a,s] else diag(n_regions)
        Qv <- if(moves) Mrate[p,,,y,seas,a,s] else matrix(0, n_regions, n_regions)
        pstep_NAA[p,,a,s] <- advance_seas(proj_NAA[p,,y,seas,a,s], Mv, proj_ZAA[p,,y,seas,a,s],
                                          Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
        pstep_NAA0[p,,a,s] <- advance_seas(proj_NAA0[p,,y,seas,a,s], Mv, natmort[p,,y,seas,a,s] * seasdur[seas],
                                           Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
      }
    }

    # Input fish into seasonal containers / fish at the end of the season / year
    if(seas < n_seas && y > 1) { # within season mortality
      proj_NAA[,,y,seas + 1,1:n_ages,] = pstep_NAA
      proj_NAA0[,,y,seas + 1,1:n_ages,] = pstep_NAA0
    } else { # age advancement
      # age advancement and enter into first season of next year
      proj_NAA[,,y + 1,1,2:n_ages,] = pstep_NAA[,,1:(n_ages - 1),] # Exponential mortality for individuals not in plus group
      proj_NAA[,,y + 1,1,n_ages,] = proj_NAA[,,y + 1,1,n_ages,] + pstep_NAA[,,n_ages,] # Acuumulate plus group
      proj_NAA0[,,y + 1,1,2:n_ages,] = pstep_NAA0[,,1:(n_ages - 1),] # Exponential mortality for individuals not in plus group
      proj_NAA0[,,y + 1,1,n_ages,] = proj_NAA0[,,y + 1,1,n_ages,] + pstep_NAA0[,,n_ages,] # Acuumulate plus group
    }

    # Derive Biomass (age0_rec: already computed above, before mortality/ageing),
    if(seas == spawn_seas && !age0_rec) {
      biom <- derive_proj_biom(y, seas, proj_NAA, proj_NAA0, WAA, MatAA, proj_ZAA, natmort, t_spawn, seasdur,
                              n_seas, n_pop, n_regions, n_ages, n_sexes,
                              sgl_seas_spawning_movement, natal_region, stray_rate,
                              Movement, Mrate, move_timing, do_recruits_move, expm_nsub = expm_nsub)
      proj_SSB[,, y] <- biom$SSB_y
      proj_Dynamic_SSB0[,,y] <- biom$Dynamic_SSB0_y
      proj_eff_SSB[,y] <- biom$eff_SSB_y
      proj_Total_Biom[,,y] <- biom$Total_Biom_y
    } # calculate biomass


    # Season-integrated abundance for the spatial Baranov under continuous movement.
    # Computed once per season across all regions, since the integral couples them.
    if(move_timing == 2) {
      proj_NAA_int <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
      for(p in 1:n_pop) {
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            proj_NAA_int[p,,a,s] <- integrate_seas_abundance(proj_NAA[p,,y,seas,a,s], proj_ZAA[p,,y,seas,a,s],
                                                            Mrate[p,,,y,seas,a,s], seasdur[seas], expm_nsub = expm_nsub)
          } # end s loop
        } # end a loop
      } # end p loop
    }

    # Derive Catches ----------------------------------------------------------
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(f in 1:n_fish_fleets) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              if(move_timing == 2) {
                # Spatial Baranov: fish redistribute among regions while dying, so catch
                # uses the season-integrated abundance rather than N (1 - exp(-Z)) / Z
                proj_CAA[p,r,y,seas,a,s,f] <- proj_ret_FAA[p,r,y,seas,a,s,f] * proj_NAA_int[p,r,a,s]
                proj_DAA[p,r,y,seas,a,s,f] <- proj_disc_FAA[p,r,y,seas,a,s,f] * proj_NAA_int[p,r,a,s]
              } else {
                # Get catch and discards at age with Baranov's
                proj_CAA[p,r,y,seas,a,s,f] <- (proj_ret_FAA[p,r,y,seas,a,s,f] / proj_ZAA[p,r,y,seas,a,s]) *
                  proj_NAA[p,r,y,seas,a,s] * (1 - exp(-proj_ZAA[p,r,y,seas,a,s]))
                proj_DAA[p,r,y,seas,a,s,f] <- (proj_disc_FAA[p,r,y,seas,a,s,f] / proj_ZAA[p,r,y,seas,a,s]) *
                  proj_NAA[p,r,y,seas,a,s] * (1 - exp(-proj_ZAA[p,r,y,seas,a,s]))
              }
            } # end s loop
          } # end a loop

          # Get total catch
          proj_Catch[p,r,y,seas,f] <- sum(proj_CAA[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f])

        } # end f loop
      } # end r loop
    } # end p loop

  } # end seas loop

  return(list(proj_NAA = proj_NAA,
              proj_NAA0 = proj_NAA0,
              proj_ZAA = proj_ZAA,
              proj_ret_FAA = proj_ret_FAA,
              proj_disc_FAA = proj_disc_FAA,
              proj_tot_FAA = proj_tot_FAA,
              proj_CAA = proj_CAA,
              proj_DAA = proj_DAA,
              proj_Catch = proj_Catch,
              proj_SSB = proj_SSB,
              proj_Dynamic_SSB0 = proj_Dynamic_SSB0,
              proj_eff_SSB = proj_eff_SSB,
              proj_Total_Biom = proj_Total_Biom))

} # end run_proj_year


#' Assemble A Trial F Matrix From The One F Per Region Being Solved For
#'
#' Builds \code{F_base + F_reg * seas_profile}. See the section header above for
#' what the three terms are and how the annual and seasonal cases fill them in.
#'
#' @param F_reg Numeric vector \code{[n_regions]}. The F being solved for.
#' @param F_base,seas_profile Numeric matrices \code{[n_regions, n_seas]}.
#' @return Numeric matrix \code{[n_regions, n_seas]} of total F, for \code{run_proj_year}.
#' @keywords internal
#' @noRd
build_proj_F <- function(F_reg, F_base, seas_profile) {
  return(F_base + as.vector(F_reg) * seas_profile) # note that F_reg goes down columns, so F_reg[r] scales row r
}


#' Regional And Seasonal Catch Produced By A Trial F
#'
#' @param F_y Numeric matrix \code{[n_regions, n_seas]}. Trial fishing mortality.
#' @param y Integer. Projection year.
#' @param state Named list of the mutable projection arrays.
#' @param tmp_rec This year's recruitment, passed through to \code{run_proj_year}.
#' @param proj_args Named list of the static \code{run_proj_year} arguments.
#' @return Numeric matrix \code{[n_regions, n_seas]} of catch in biomass, summed
#'   over populations and fleets.
#' @keywords internal
#' @noRd
proj_catch_at_F <- function(F_y, y, state, tmp_rec, proj_args) {

  yr <- do.call(run_proj_year, c(list(y = y, F_y = F_y, tmp_rec = tmp_rec), state, proj_args))

  catch_mat <- array(0, dim = c(proj_args$n_regions, proj_args$n_seas))
  for(r in 1:proj_args$n_regions) {
    for(seas in 1:proj_args$n_seas) catch_mat[r,seas] <- sum(yr$proj_Catch[,r,y,seas,])
  } # end r loop

  return(catch_mat)
}


#' Reduce A Catch Matrix To The Quantity A Target Is Set On
#'
#' @param catch_mat Numeric matrix \code{[n_regions, n_seas]}.
#' @param target_seas Integer season the target applies to, or \code{NULL} for
#'   an annual target, which sums across seasons.
#' @return Numeric vector \code{[n_regions]}.
#' @keywords internal
#' @noRd
proj_target_catch <- function(catch_mat, target_seas) {
  if(is.null(target_seas)) return(rowSums(catch_mat))
  return(catch_mat[, target_seas])
}


#' Log Scale Catch Residual For The Joint Regional Solve
#'
#' Residuals and unknowns both sit on the log scale: F stays positive with no
#' constraints to enforce, and a residual in log catch is a relative catch error,
#' which is the tolerance the caller specifies.
#'
#' @param theta Numeric vector. log F for the free regions.
#' @param F_reg_fixed Numeric vector \code{[n_regions]}. F for the regions not
#'   being solved (zero targets, or regions already capped at the F bound).
#' @param free Integer vector. Indices of the regions being solved.
#' @param target Numeric vector \code{[n_regions]}. Catch targets.
#' @param seas_profile,F_base Passed to \code{build_proj_F}.
#' @param target_seas Passed to \code{proj_target_catch}.
#' @param y,state,tmp_rec,proj_args Passed to \code{proj_catch_at_F}.
#' @param catch_f_max Numeric. Upper bound on F.
#' @return Numeric vector, one residual per free region.
#' @keywords internal
#' @noRd
proj_log_catch_resid <- function(theta, F_reg_fixed, free, target, seas_profile, F_base,
                                 target_seas, y, state, tmp_rec, proj_args, catch_f_max) {

  F_reg <- F_reg_fixed
  F_reg[free] <- pmin(exp(theta), catch_f_max)
  catch_mat <- proj_catch_at_F(build_proj_F(F_reg, F_base, seas_profile), y, state, tmp_rec, proj_args)
  realized_catch <- proj_target_catch(catch_mat, target_seas)[free]

  return(log(pmax(realized_catch, 1e-12)) - log(target[free]))
}


#' Solve One Block Of Fishing Mortalities Against A Catch Target
#'
#' Finds the one F per region that makes realized catch match \code{target}.
#' Catch rises monotonically with each region's F, so a block with one free region
#' bisects, which needs no start value and lets the bracket double as a
#' feasibility check. Regions in a block are coupled by movement, so a block with
#' several free regions solves jointly instead.
#'
#' @param y Integer. Projection year.
#' @param target Numeric vector \code{[n_regions]}. Catch targets; 0 means no
#'   fishing rather than something to solve.
#' @param seas_profile,F_base Numeric matrices \code{[n_regions, n_seas]} placing
#'   \code{F_reg} into the year, see \code{build_proj_F}.
#' @param target_seas Integer or \code{NULL}, see \code{proj_target_catch}.
#' @param state,tmp_rec,proj_args Projection state and inputs.
#' @param f_start Numeric vector \code{[n_regions]}. Starting values for the joint
#'   solve, normally the previous year's F.
#' @param catch_f_max,catch_tol,catch_max_iter Solver settings, documented in
#'   \code{\link{Do_Population_Projection}}.
#' @param label Character. Names what failed in warning messages.
#'
#' @return Named list with \code{F_reg}, the solved F per region, and
#'   \code{resid}, the relative miss on each target.
#' @keywords internal
#' @noRd
solve_proj_F_catch <- function(y, target, seas_profile, F_base, target_seas,
                               state, tmp_rec, proj_args, f_start,
                               catch_f_max, catch_tol, catch_max_iter, label) {

  n_regions <- proj_args$n_regions
  F_reg <- rep(0, n_regions)
  capped <- rep(FALSE, n_regions)
  free <- which(target > 0) # a zero target is F = 0, not something to solve for

  if(length(free) > 0) {

    # Bounding the catch
    F_reg_cap <- F_reg
    F_reg_cap[free] <- catch_f_max
    cap_catch <- proj_target_catch(proj_catch_at_F(build_proj_F(F_reg_cap, F_base, seas_profile),
                                                    y, state, tmp_rec, proj_args), target_seas)
    infeas <- free[cap_catch[free] < target[free]]

    if(length(infeas) > 0) {
      warning(paste0("Catch target for ", label, " is not reachable in region(s) ",
                     paste(infeas, collapse = ", "), " at the F bound catch_f_max = ",
                     catch_f_max, ". F is capped there and the target is undershot."))
      F_reg[infeas] <- catch_f_max
      capped[infeas] <- TRUE
      free <- setdiff(free, infeas) # anything left solves against the capped regions
    }
  }

  # One region only: bisect the bracket already shown to contain the root
  if(length(free) == 1) {
    lb <- 0
    ub <- catch_f_max
    for(i in seq_len(catch_max_iter)) {
      F_reg[free] <- (lb + ub) / 2
      realized_catch <- proj_target_catch(proj_catch_at_F(build_proj_F(F_reg, F_base, seas_profile),
                                                y, state, tmp_rec, proj_args), target_seas)[free]
      if(abs(realized_catch - target[free]) <= catch_tol * target[free]) break
      if(realized_catch < target[free]) lb <- F_reg[free] else ub <- F_reg[free]
    } # end i loop
  }

  # Several regions: joint solve on the log F scale
  if(length(free) > 1) {

    F_reg_start <- F_reg
    F_reg_start[free] <- pmin(pmax(f_start[free], 1e-4), catch_f_max)

    # get starting point
    c0 <- proj_target_catch(proj_catch_at_F(build_proj_F(F_reg_start, F_base, seas_profile),
                                             y, state, tmp_rec, proj_args), target_seas)
    scaling <- ifelse(c0[free] > 0, target[free] / c0[free], 1)
    F_reg_start[free] <- pmin(pmax(F_reg_start[free] * scaling, 1e-8), catch_f_max)

    # solve for F
    solve_out <- nleqslv::nleqslv(
      log(F_reg_start[free]),
      proj_log_catch_resid, # function to be optimized across (computes the projection cycle)
      F_reg_fixed = F_reg,
      free = free,
      target = target,
      seas_profile = seas_profile,
      F_base = F_base,
      target_seas = target_seas,
      y = y,
      state = state,
      tmp_rec = tmp_rec,
      proj_args = proj_args,
      catch_f_max = catch_f_max,
      control = list(ftol = catch_tol, xtol = 1e-10,  maxit = catch_max_iter)
    )
    F_reg[free] <- pmin(exp(solve_out$x), catch_f_max)
  }

  # Relative miss on the F actually being returned, handed back to the caller
  realized_catch <- proj_target_catch(proj_catch_at_F(build_proj_F(F_reg, F_base, seas_profile),
                                            y, state, tmp_rec, proj_args), target_seas)
  resid <- rep(0, n_regions)
  pos <- target > 0
  resid[pos] <- (realized_catch[pos] - target[pos]) / target[pos]

  missed <- which(pos & !capped & abs(resid) > max(catch_tol, 1e-4))
  if(length(missed) > 0) {
    warning(paste0("Catch target for ", label, " did not converge in region(s) ",
                   paste(missed, collapse = ", "), ". Largest relative catch error is ",
                   signif(max(abs(resid[missed])), 3), "."))
  }

  return(list(F_reg = F_reg, resid = resid))

} # end solve_proj_F_catch


#' Solve A Projection Year's Fishing Mortality Against Its Catch Target
#'
#' Breaks the year into blocks and hands each to \code{solve_proj_F_catch}: one
#' block for an annual target, one per season for seasonal ones, swept forward
#' with each solved F kept in \code{F_base}. See the section header above for why
#' the sweep is exact.
#'
#' @param y Integer. Projection year.
#' @param target Numeric \code{[n_regions]} for annual targets, \code{[n_regions,
#'   n_seas]} for seasonal ones.
#' @param seasonal Logical. Whether targets are seasonal.
#' @param seas_share Numeric matrix \code{[n_regions, n_seas]}. Terminal year
#'   seasonal shares of annual F, the profile for annual targets.
#' @param f_start Numeric matrix \code{[n_regions, n_seas]}. Previous year's F,
#'   used to start the joint solves.
#' @param state,tmp_rec,proj_args Projection state and inputs.
#' @param catch_f_max,catch_tol,catch_max_iter Solver settings.
#'
#' @return Named list with \code{F_y} \code{[n_regions, n_seas]} and
#'   \code{resid}, shaped like \code{target}.
#' @keywords internal
#' @noRd
solve_proj_year_F <- function(y, target, seasonal, seas_share, f_start,
                              state, tmp_rec, proj_args,
                              catch_f_max, catch_tol, catch_max_iter) {

  n_regions <- proj_args$n_regions
  n_seas <- proj_args$n_seas
  F_base <- array(0, dim = c(n_regions, n_seas))

  # Annual targets: one annual F per region, split at the terminal year splits at the end
  if(!seasonal) {
    sol <- solve_proj_F_catch(
      y = y,
      target = target,
      seas_profile = seas_share,
      F_base = F_base,
      target_seas = NULL,
      state = state,
      tmp_rec = tmp_rec,
      proj_args = proj_args,
      f_start = rowSums(f_start),
      catch_f_max = catch_f_max,
      catch_tol = catch_tol,
      catch_max_iter = catch_max_iter,
      label = paste0("projection year ", y)
    )
    return(list(F_y = build_proj_F(sol$F_reg, F_base, seas_share), resid = sol$resid))
  }

  # Seasonal targets
  resid <- array(0, dim = c(n_regions, n_seas))
  for(seas in 1:n_seas) {
    seas_profile <- array(0, dim = c(n_regions, n_seas))
    seas_profile[,seas] <- 1 # this season's value is just this season's F
    sol <- solve_proj_F_catch(
      y = y,
      target = target[,seas],
      seas_profile = seas_profile,
      F_base = F_base,
      target_seas = seas,
      state = state,
      tmp_rec = tmp_rec,
      proj_args = proj_args,
      f_start = f_start[,seas],
      catch_f_max = catch_f_max,
      catch_tol = catch_tol,
      catch_max_iter = catch_max_iter,
      label = paste0("projection year ", y, ", season ", seas)
    )
    F_base[,seas] <- sol$F_reg # settled, and passed into the next season's solve
    resid[,seas] <- sol$resid
  } # end seas loop

  return(list(F_y = F_base, resid = resid))

} # end solve_proj_year_F
