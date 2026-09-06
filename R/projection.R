# Stage 3 of 3: post fit
#
# Forward projection off a fitted model under a specified catch or fishing mortality, with optional
# stochastic recruitment. Short horizon advice, not the closed loop simulation in sim_closed_loop.R.

#' Do Population Projections
#'
#' Projects population dynamics forward in time under alternative recruitment
#' and fishing mortality scenarios. The model initializes from terminal
#' assessment quantities and advances numbers-at-age through recruitment,
#' seasonal movement, mortality, ageing, and harvest control rules across
#' multiple seasons and years.
#'
#' Population dynamics are tracked over
#' \code{[population x region x year x season x age x sex]}. Recruitment is
#' generated annually and then distributed across seasons using
#' \code{rec_seas_prop}, allowing intra-annual timing of recruitment within
#' the first age class.
#'
#' @param n_proj_yrs Integer. Number of projection years.
#' @param n_pop Integer. Number of populations (may exceed regions when
#'   natal homing is modeled).
#' @param n_regions Integer. Number of spatial regions.
#' @param n_ages Integer. Number of age classes including the plus group.
#' @param n_sexes Integer. Number of sexes.
#' @param sexratio Array `[n_pop, n_regions, n_proj_yrs, n_sexes]`.
#'   Recruitment sex ratio used to allocate projected recruits by sex.
#' @param n_fish_fleets Integer. Number of fishing fleets.
#' @param do_recruits_move Integer (0 or 1). Whether age-1 recruits are
#'   subject to movement. Default = 0.
#' @param rec_seas_prop Array `[n_pop, n_seas]`. Proportion of annual
#'   recruitment entering in each season. Must sum to 1 across seasons for
#'   each population.
#' @param recruitment Array `[n_pop, n_regions, n_yrs]`. Historical
#'   recruitment used to condition stochastic projection options.
#' @param terminal_NAA Array `[n_pop, n_regions, n_seas, n_ages, n_sexes]`.
#'   Fished numbers-at-age in the terminal assessment year.
#' @param terminal_NAA0 Array `[n_pop, n_regions, n_seas, n_ages, n_sexes]`.
#'   Unfished numbers-at-age in the terminal assessment year.
#' @param terminal_F Array `[n_regions, n_seas, n_fish_fleets]`. Terminal
#'   fishing mortality; sets F in projection year 1 and defines the seasonal
#'   F ratios applied in subsequent years.
#' @param natmort Array of natural mortality, a rate per year in each season.
#'   Either `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]` or the same
#'   without the season dim, which is expanded across seasons.
#'   Annual natural mortality-at-age, scaled internally by season duration.
#' @param WAA Array `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]`.
#'   Weight-at-age used in spawning biomass calculations.
#' @param WAA_fish Array
#'   `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]`.
#'   Fishery weight-at-age used in catch biomass calculations.
#' @param MatAA Array `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]`.
#'   Maturity-at-age.
#' @param fish_sel Array `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]`.
#'   Fishery selectivity-at-age.
#' @param Movement Array
#'   `[n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]`.
#'   Seasonal movement transition matrices.
#' @param Mrate Array dimensioned like `Movement`, holding the instantaneous
#'   movement rates (the generator) rather than the realized transition
#'   fractions. Only read when `move_timing` is 1 or 2. `NULL` (default) is
#'   valid for `move_timing = 0`, where movement is applied as a transition
#'   matrix and no generator is needed.
#' @param expm_nsub Integer controlling how the matrix exponential is evaluated under
#'   `move_timing = 2`: `0` uses `Matrix::expm`, `n >= 1` uses `n` implicit backward
#'   Euler substeps. See [mat_exp()].
#' @param move_timing Integer. When movement happens relative to mortality
#'   within a season. `0` (default) applies movement first and mortality
#'   afterwards; `1` applies mortality first and movement afterwards; `2`
#'   runs the two continuously and simultaneously, which also switches
#'   catch-at-age to the spatial Baranov form built on season-integrated
#'   abundance. Must match the timing used to derive the reference points the
#'   projection is run against.
#' @param sgl_seas_spawning_movement Array
#'   `[n_pop, n_regions, n_regions, n_proj_yrs, n_ages, n_sexes]`.
#'   Spawning movement matrix applied when `n_seas = 1` and `n_pop > 1`
#'   to redistribute fish to natal grounds prior to SSB calculation. Only read in
#'   that case, so `NULL` (the default) is valid otherwise.
#' @param stray_rate Array `[n_pop, n_proj_yrs]`. Per-population stray rate
#'   used when accumulating effective SSB contributions across populations. Only
#'   read when `n_pop > 1`, so `NULL` (the default) is valid otherwise.
#' @param f_ref_pt Array `[n_regions, n_proj_yrs]`. Fishing mortality
#'   reference point (e.g., F_MSY) or fixed input F, depending on
#'   `fmort_opt`.
#' @param b_ref_pt Array `[n_pop, n_regions, n_proj_yrs]`. Biomass reference
#'   point used in harvest control rules.
#' @param HCR_function Function. Harvest control rule with arguments `x`
#'   (SSB), `frp` (F reference point), and `brp` (B reference point). A rule that
#'   also declares a `state` argument (or `...`) is handed this year's population
#'   as a named list, so a rule can be written on more than spawning biomass:
#'   `y` (projection year), `r` (region the F is being set for), `NAA`, `SSB`,
#'   `Total_Biom` and `Catch`. Mean weight, age structure and last year's catch
#'   all follow from those. Rules that do not declare it are called exactly as
#'   before and the state is not assembled.
#' @param recruitment_opt Character. Recruitment scenario:
#'   `"inv_gauss"`, `"mean_rec"`, `"zero"`, or `"bh_rec"`.
#' @param fmort_opt Character. Fishing mortality scenario:
#'   `"HCR"`, `"HCR_global"`, `"Input"`, or `"Catch"`. `"Catch"` solves each
#'   projection year's fishing mortality so that realized catch matches
#'   `catch_input`, leaving every other model quantity untouched.
#' @param catch_input Catch targets in biomass, used when `fmort_opt = "Catch"`
#'   and ignored otherwise. Either an array `[n_regions, n_proj_yrs]` of annual
#'   targets or `[n_regions, n_proj_yrs, n_seas]` of seasonal ones; which shape
#'   is supplied decides what gets solved. `catch_input[r, y]` is the catch
#'   removed from region `r` during projection year `y`, indexed the same way
#'   `proj_Catch` is in the returned list rather than with the one year lag
#'   `f_ref_pt` does. Targets are totals over populations and fleets, and over
#'   seasons too in the annual case. A target of 0 sets `F = 0` there without a
#'   solve.
#'
#'   With annual targets one annual F per region is solved for and split over
#'   seasons at the terminal year seasonal shares, exactly as the other
#'   `fmort_opt` settings do. With seasonal targets that constraint is released
#'   and a separate F is solved per region and season; the split across fleets
#'   within a season still stays at terminal year ratios, so fleet specific
#'   targets are not supported either way. A season the terminal year did not
#'   fish has no fleet split to inherit and so can take no catch, which is an
#'   error rather than a silent zero.
#'
#'   Not every projection year has to have a target. Set a year to `NA` and it
#'   falls back to `catch_fallback_opt` instead, which is the usual shape of catch
#'   advice: a year or two of agreed catch followed by the harvest control rule.
#'   `NA` (no target, use the fallback) and `0` (a target of no fishing) are
#'   different things. A year has to be all target or all `NA` across regions and
#'   seasons, since splitting one annual fallback F across only some seasons has
#'   no defensible reading; a partly specified year is an error.
#'
#'   Column 1 is only used when `catch_terminal_yr = TRUE`; see that argument.
#' @param catch_fallback_opt Character. Which rule sets F in the projection years
#'   `catch_input` leaves `NA`: `"HCR"`, `"HCR_global"`, or `"Input"`. Defaults to
#'   `"HCR"` under `fmort_opt = "Catch"` and to `fmort_opt` itself otherwise,
#'   where it is unused.
#'   Ignored unless `fmort_opt = "Catch"`, and its usual inputs (`f_ref_pt`,
#'   `b_ref_pt`, `HCR_function`) are only needed if some year actually falls back.
#'   Mind the indexing difference when mixing the two: `catch_input[r, y]` is the
#'   catch taken in year `y`, but `f_ref_pt[r, y]` sets F in year `y + 1`, which
#'   is the lag the HCR and Input options have always kept.
#' @param catch_terminal_yr Logical. Whether projection year 1, which replays the
#'   terminal assessment year, is also solved against its catch target rather
#'   than fished at `terminal_F`. Default `FALSE`. Set it `TRUE` for the common
#'   assessment case where the terminal year's catch is itself a projection
#'   because the year is not yet complete. Note that this overrides the F the
#'   assessment estimated for that year, and so changes the numbers-at-age
#'   entering year 2. Note also that with `n_seas > 1` the terminal year takes
#'   all its seasons from `terminal_NAA` rather than propagating them, so only
#'   the last season's F feeds year 2; the earlier seasons still take their
#'   catch, but do not otherwise propagate.
#' @param catch_f_max Numeric. Upper bound on the F searched when
#'   `fmort_opt = "Catch"`. Default 5. A target that cannot be taken even at this
#'   F is unreachable, in which case F is capped here, the target is undershot,
#'   and a warning names the regions involved.
#' @param catch_tol Numeric. Relative catch tolerance for the F solver.
#'   Default 1e-6.
#' @param catch_max_iter Integer. Maximum solver iterations per projection year.
#'   Default 100.
#' @param t_spawn Numeric scalar. Fraction of the spawning season elapsed
#'   before spawning; used for mid-season SSB calculations.
#' @param srr_opt Named list of inputs for deterministic stock-recruit
#'   recruitment when `recruitment_opt` is `"bh_rec"` or `"ricker_rec"`.
#'   The curve itself is taken from `recruitment_opt`, so the same list serves
#'   both. Formerly `bh_rec_opt`.
#' @param bh_rec_opt Deprecated. Former name of `srr_opt`; supplying it warns
#'   and forwards. Supplying both is an error. This list is passed
#'   directly to \code{\link{Get_Det_Recruitment}} and must contain all
#'   required arguments for that function.
#'
#'   Required elements and their expected dimensions include:
#'   \describe{
#'     \item{\code{R0}}{Numeric vector \code{[n_pop]}. Unfished recruitment.}
#'     \item{\code{h}}{Numeric array \code{[n_pop, n_regions]}. Steepness.}
#'     \item{\code{rec_region_prop}}{Numeric array
#'       \code{[n_pop, n_regions]}. Recruitment allocation across regions
#'       (sums to 1 across regions).}
#'     \item{\code{rec_seas_prop}}{Numeric array
#'       \code{[n_pop, n_seas]}. Seasonal recruitment proportions
#'       (sums to 1 across seasons).}
#'     \item{\code{SSB}}{Numeric array
#'       \code{[n_pop, n_regions, n_yrs]}. Historical spawning biomass,
#'       to which projected SSB is appended internally.}
#'     \item{\code{WAA}}{Array
#'       \code{[n_pop, n_regions, n_seas, n_ages]}. Weight-at-age.}
#'     \item{\code{MatAA}}{Array
#'       \code{[n_pop, n_regions, n_seas, n_ages]}. Maturity-at-age.}
#'     \item{\code{natmort}}{Array
#'       \code{[n_pop, n_regions, n_seas, n_ages]}. Natural mortality, a rate per
#'       year in each season. Also accepted without the season dim.}
#'     \item{\code{Movement}}{Array
#'       \code{[n_pop, n_regions, n_regions, n_seas, n_ages]}. Movement
#'       transition matrices.}
#'     \item{\code{sgl_seas_spawning_movement}}{Array
#'       \code{[n_pop, n_regions, n_regions, n_ages]}. Spawning movement
#'       (single-season case).}
#'     \item{\code{stray_rate}}{Numeric vector \code{[n_pop]}. Straying rates.}
#'     \item{\code{init_F}}{Array
#'       \code{[n_regions, n_seas, n_fish_fleets]}. Initial fishing mortality.}
#'     \item{\code{fish_sel}}{Array
#'       \code{[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]}. Total selectivity.}
#'     \item{\code{ret_sel}}{Array
#'       \code{[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]}. Retention selectivity.}
#'     \item{\code{dmr}}{Array
#'       \code{[n_regions, n_seas, n_fish_fleets]}. Discard mortality rates.}
#'     \item{\code{sex_ratio_f}}{Numeric array
#'       \code{[n_pop, n_regions]}. Female recruitment proportion.}
#'   }
#'
#'   Additional scalar inputs include \code{rec_dd}, \code{rec_lag},
#'   \code{n_pop}, \code{n_regions}, \code{n_ages}, \code{n_seas},
#'   \code{spawn_seas}, \code{seasdur}, \code{t_spawn}, and
#'   \code{do_recruits_move}.
#'
#'   Spawning biomass used in recruitment is constructed internally by
#'   combining \code{srr_opt$SSB} with projected SSB values during the
#'   simulation.
#'
#'   \code{srr_opt$rec_lag = 1} is the classic lagged case: each
#'   projection year's recruitment is computed up front from the prior
#'   year's SSB, exactly as \code{recruitment_opt = "inv_gauss"}/
#'   \code{"mean_rec"} are. \code{srr_opt$rec_lag = 0} is age-0
#'   recruitment: recruitment for year \code{y} is computed from year
#'   \code{y}'s own SSB once \code{spawn_seas} is reached within that year's
#'   season loop, and is inserted no earlier than \code{spawn_seas}
#'   (\code{rec_seas_prop} must be zero for every season before
#'   \code{spawn_seas} in that case). Reference points and the seasonal SBPR
#'   calculation used to get \code{srr_opt$WAA}/\code{MatAA}/etc. are
#'   unaffected by this choice, \code{rec_lag} only changes which year's
#'   SSB feeds the Beverton-Holt curve, not the per-recruit math itself.
#'
#' @param n_seas Integer. Number of seasons. Default = 1.
#' @param seasdur Numeric vector `[n_seas]`. Duration of each season as a
#'   fraction of the year.
#' @param spawn_seas Integer. Spawning season index.
#' @param natal_region Integer vector `[n_pop]`. Natal region for each
#'   population. Only read when `n_pop > 1`, so `NULL` (the default) is valid for
#'   single-population models.
#' @param dmr Array \code{[n_regions, n_seas, n_fish_fleets]}. Discard mortality rate.
#'   Default behavior is no discard mortality (\code{dmr = 0}). When combined with
#'   \code{ret_sel = 1}, this implies no discarding within a given fleet (all catch is retained).
#' @param ret_sel Array \code{[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]}. Retention
#'   selectivity-at-age. Default behavior corresponds to full retention (\code{ret_sel = 1}),
#'   meaning all captured fish are retained unless otherwise specified.
#'
#'
#' @return A named list of projected quantities. Year index 1 is the terminal
#'   assessment year replayed, so year 2 is the first projected year
#'   and catch advice for terminal year + 1 is read from index 2.
#'
#'   Several arrays have a trailing `n_proj_yrs + 1` year slot, which is used
#'   inconsistently and is noted per element below. In short: `proj_NAA`,
#'   `proj_NAA0`, `proj_F` and `proj_F_seas` fill it, and `proj_ZAA`,
#'   `proj_ret_FAA` and `proj_disc_FAA` leave it at 0.
#'
#' \describe{
#'   \item{\code{proj_F}}{Array `[n_regions, n_proj_yrs + 1]`. Annual fishing
#'     mortality by region, summed over seasons and fleets. The trailing column
#'     holds the F the harvest control rule or input would apply in the year
#'     after the projection; it stays 0 under `fmort_opt = "Catch"`, where there
#'     is no further year to solve a target for.}
#'   \item{\code{proj_F_seas}}{Array `[n_regions, n_proj_yrs + 1, n_seas]`. The
#'     same fishing mortality broken out by season, so
#'     `rowSums(proj_F_seas[, y, ])` recovers `proj_F[, y]` for every `y`,
#'     including the trailing column. This is the only
#'     place the answer lives when seasonal catch targets are used, since an
#'     annual total cannot represent a seasonal solve.}
#'   \item{\code{proj_ret_FAA}}{Array
#'     `[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes, n_fish_fleets]`.
#'     Retained fishing mortality-at-age, i.e. the component that generates
#'     landed catch. Only years `1:n_proj_yrs` are filled; the trailing year slot
#'     stays 0.}
#'   \item{\code{proj_disc_FAA}}{Array dimensioned as `proj_ret_FAA`, and filled
#'     over the same years. Discard fishing mortality-at-age, i.e. the component
#'     killed but not landed, set by `ret_sel` and `dmr`. Total fishing
#'     mortality-at-age is the sum of the two.}
#'   \item{\code{proj_Catch}}{Array
#'     `[n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets]`. Retained catch in
#'     biomass. Built from `proj_ret_FAA`, so discard mortality acts on the
#'     population but is not counted here, and this is the quantity
#'     `catch_input` is matched against.}
#'   \item{\code{proj_SSB}}{Array `[n_pop, n_regions, n_proj_yrs]`. Female
#'     spawning biomass, accumulated in `spawn_seas` with the `t_spawn` mortality
#'     correction. Halved when `n_sexes = 1`.}
#'   \item{\code{proj_eff_SSB}}{Array `[n_pop, n_proj_yrs]`. Effective spawning
#'     biomass at each population's natal region, aggregating contributions from
#'     every population with cross-population terms scaled by `stray_rate`. Equal
#'     to spawning biomass summed across regions when `n_pop = 1`.}
#'   \item{\code{proj_Total_Biom}}{Array `[n_pop, n_regions, n_proj_yrs]`. Total
#'     biomass over all ages and both sexes, accumulated at the same point in the
#'     season as `proj_SSB` and using the same definition the estimation model
#'     uses for `Total_Biom`, so the projected series continues the estimated one
#'     without a discontinuity at the terminal year.}
#'   \item{\code{proj_Dynamic_SSB0}}{Array `[n_pop, n_regions, n_proj_yrs]`.
#'     Spawning biomass the population would have kept under the same realized
#'     recruitment but no fishing, for dynamic depletion.}
#'   \item{\code{proj_NAA}}{Array
#'     `[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes]`. Fished
#'     numbers-at-age kept at the start of each season, before that season's
#'     mortality and ageing. Under the default `move_timing = 0` movement has
#'     already been applied at this point; under `move_timing` 1 and 2 it has
#'     not, since movement is deferred into the mortality step. The trailing year
#'     slot is filled, and holds the numbers passed into the year after the
#'     projection ends.}
#'   \item{\code{proj_NAA0}}{Array dimensioned as `proj_NAA`. The unfished
#'     counterpart, decremented by natural mortality alone.}
#'   \item{\code{proj_ZAA}}{Array
#'     `[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes]`. Total
#'     mortality-at-age for the season: natural mortality scaled by season
#'     duration, plus retained and discard fishing mortality summed over fleets.
#'     Only years `1:n_proj_yrs` are filled; the trailing year slot stays 0.}
#'   \item{\code{proj_catch_resid}}{Array shaped like `catch_input`:
#'     `[n_regions, n_proj_yrs]` for annual targets, `[n_regions, n_proj_yrs, n_seas]`
#'     for seasonal ones. Relative miss on each catch target,
#'     `(realized - target) / target`, and `NA` for years with no target
#'     (including every year when `fmort_opt != "Catch"`). Should be at or below
#'     `catch_tol` wherever the solve converged, and is worth checking directly
#'     rather than relying on warnings alone.}
#' }
#'
#' @details
#' Each projection year proceeds as follows when
#' \code{recruitment_opt != "bh_rec"} or \code{srr_opt$rec_lag != 0}
#' (the classic case):
#' \enumerate{
#'   \item Annual recruitment is generated and allocated across regions and
#'   sexes. Seasonal recruitment is then distributed within the first age
#'   class using \code{rec_seas_prop}, with additional recruits entering in
#'   seasons \code{seas > 1}.
#'   \item Fishing mortality-at-age is constructed from annual F, seasonal
#'   F ratios derived from the terminal year, and selectivity.
#'   \item Movement is applied at each seasonal step via transition matrices.
#'   Age-1 movement is optional via \code{do_recruits_move}.
#'   \item Within-season mortality is applied using exponential decay. At the
#'   end of the final season, individuals age forward and the plus group
#'   accumulates survivors.
#'   \item Spawning biomass is computed in \code{spawn_seas} using a
#'   mid-season mortality correction. For natal homing models with a single
#'   season, spawning movement is applied prior to SSB calculation.
#'   \item Catch is calculated using the Baranov equation and aggregated to
#'   biomass using fishery-specific weights.
#'   \item Fishing mortality for the next year is updated via the specified
#'   harvest control rule or fixed input.
#' }
#'
#' When \code{srr_opt$rec_lag == 0} (age-0 recruitment), steps 1 and 5
#' above are reordered within \code{spawn_seas}: movement is applied first,
#' spawning biomass is computed from the survivor population alone (no new
#' recruits exist yet), that SSB is used to generate this year's
#' recruitment, and only then are the recruits inserted (no earlier than
#' \code{spawn_seas}) - immediately before mortality/ageing runs for that
#' season, so the new cohort is advanced exactly like any other
#' seasonal recruit pulse. Years \code{y > 1} generate recruitment this way;
#' year 1 holds the supplied terminal assessment state forward with no new
#' recruitment event, matching the classic case.
#'
#' Under \code{fmort_opt = "Catch"} step 7 moves to the front of the following
#' year instead: the F that lands a catch target depends on that year's own
#' numbers-at-age, not on the previous year's spawning biomass, so it cannot be
#' set at the end of the previous year the way the HCR and Input options are.
#' The year is run repeatedly at trial F values until realized catch matches the
#' target, then run once more at the accepted F and committed. Steps 1 to 6 are
#' otherwise unchanged, and no demographic input is modified: only F moves.
#' Regions are solved jointly rather than one at a time, because between-season
#' movement (and, under \code{move_timing = 2}, the season-integrated abundance)
#' makes each region's catch depend on the F set in every other region. Seasonal
#' targets are instead swept forward one season at a time, which is exact because
#' a season's catch depends only on the F in that season and earlier ones.
#'
#' Note that the catch solved against is retained catch, matching
#' \code{proj_Catch}. Discard mortality still acts on the population through
#' \code{dmr} and \code{ret_sel}, so a fleet that discards will exert more total
#' F than the target alone implies.
#'
#' Effective spawning biomass at each population's natal region aggregates
#' contributions from all populations, with cross-population contributions
#' scaled by \code{stray_rate} and normalized by the number of populations
#' in each natal region.
#'
#' When \code{n_sexes = 1}, spawning biomass is multiplied by 0.5. When
#' \code{n_regions = 1}, movement is skipped.
#'
#' @section Differentiating through the projection:
#'
#' The projection uses RTMB's replacement operators, so it can be taped with
#' \code{\link[RTMB]{MakeTape}} or \code{\link[RTMB]{MakeADFun}} and handed to an
#' optimizer with an exact gradient. That gives an F schedule solved against an
#' objective rather than scanned over a grid, and the delta method on projected
#' quantities from a fitted model's parameters.
#'
#' Two options are refused when the inputs are AD types, because neither has a
#' derivative and both would otherwise return a gradient that is wrong rather
#' than an error. \code{recruitment_opt = "inv_gauss"} draws recruitment at
#' random, and \code{fmort_opt = "Catch"} inverts the catch target with a
#' numerical solve, so the F it hands back has no derivative. Tape under
#' \code{"mean_rec"}, \code{"bh_rec"} or \code{"ricker_rec"} with
#' \code{fmort_opt = "Input"}.
#'
#' A control rule stops on its own rather than being refused here. The usual
#' threshold rule branches on stock status, and a comparison against an AD
#' spawning biomass raises an error inside RTMB. Optimizing through a control
#' rule means writing a smooth one.
#'
#' @param rec_devs Optional array \code{[n_pop, n_regions, n_proj_yrs]} of
#'   multiplicative deviations applied to whatever recruitment
#'   \code{recruitment_opt} produces, so a deterministic option becomes a
#'   stochastic one under deviations the caller draws and owns. Defaults to
#'   \code{NULL}, which leaves recruitment as the option alone gives it.
#'   Projection year 1 replays the terminal assessment year and generates no
#'   recruitment, so its slice is never read. Drawing the deviations outside is
#'   what makes a set of replicates share the same recruitment across management
#'   procedures, and what lets the projection be differentiated with respect to
#'   the rule while the deviations are kept fixed.
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
                                       f_start = if(y > 1) array(proj_F_seas[,y-1,], dim = c(n_regions, n_seas)) else apply(terminal_F, c(1,2), sum),
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
    if(fmort_opt != 'Catch' || (y + 1 <= n_proj_yrs && !catch_yr_targeted[y+1])) {

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
          proj_F[r,y+1] <- HCR_function(x = sum(proj_SSB[,r,y]),
                                        frp = f_ref_pt[r,y],
                                        brp = sum(b_ref_pt[,r,y]))
        else
          proj_F[r,y+1] <- HCR_function(x = sum(proj_SSB[,r,y]),
                                        frp = f_ref_pt[r,y],
                                        brp = sum(b_ref_pt[,r,y]),
                                        state = c(hcr_state, list(r = r)))
      }

      if(fmort_rule == 'HCR_global') {
        if(is.null(hcr_state))
          proj_F[r,y+1] <- HCR_function(x = sum(proj_SSB[,,y]),
                                        frp = f_ref_pt[r,y],
                                        brp = sum(b_ref_pt[,,y]))
        else
          proj_F[r,y+1] <- HCR_function(x = sum(proj_SSB[,,y]),
                                        frp = f_ref_pt[r,y],
                                        brp = sum(b_ref_pt[,,y]),
                                        state = c(hcr_state, list(r = r)))
      }

      # Project F using User Inputs ---------------------------------------------
      if(fmort_rule == 'Input') proj_F[r,y+1] <- f_ref_pt[r,y]

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
          proj_NAA[,,y,seas+1,1:n_ages,] = pstep_NAA
          proj_NAA0[,,y,seas+1,1:n_ages,] = pstep_NAA0
        } else { # age advancement
          # age advancement and enter into first season of next year
          proj_NAA[,,y+1,1,2:n_ages,] = pstep_NAA[,,1:(n_ages-1),] # Exponential mortality for individuals not in plus group
          proj_NAA[,,y+1,1,n_ages,] = proj_NAA[,,y+1,1,n_ages,] + pstep_NAA[,,n_ages,] # Acuumulate plus group
          proj_NAA0[,,y+1,1,2:n_ages,] = pstep_NAA0[,,1:(n_ages-1),] # Exponential mortality for individuals not in plus group
          proj_NAA0[,,y+1,1,n_ages,] = proj_NAA0[,,y+1,1,n_ages,] + pstep_NAA0[,,n_ages,] # Acuumulate plus group
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
