# Projection inputs for the reference point invariants, written once so a test can say which fishing
# mortality it wants and nothing else.
#
# Everything is kept at its terminal-year value, which is the constant-biology, constant-F setting the
# equilibrium reference points are defined in and what makes them comparable to the projection.

#' Project the packaged sablefish stock at a constant fishing mortality
#'
#' @param f Fishing mortality to apply in every projection year.
#' @param n_proj_yrs Years to project. Long enough to reach equilibrium.
#' @param rec_window How many years back from the terminal year mean recruitment
#'   is taken over. Clipped to the series length.
#' @param rec_yrs Explicit years to average recruitment over, overriding
#'   \code{rec_window}. Comparing against a reference point means matching the
#'   years that reference point was scaled by, which is
#'   \code{calc_rec_st_yr:(n_yrs - rec_age)}; a different window is a different
#'   mean recruitment and so a different biomass scale.
#' @param data,rep The model to project. Defaults to the packaged sablefish
#'   assessment; any fitted model can be passed, which is what lets the
#'   cross-check run at configurations the packaged one does not cover.
#' @param recruitment_opt \code{"mean_rec"} or \code{"bh_rec"}. Yield has an
#'   interior maximum in F only under a stock-recruit curve; with mean
#'   recruitment it rises without bound, so MSY needs \code{"bh_rec"}.
#' @param bh_rec_opt The deprecated spelling of \code{srr_opt}, kept only so
#'   the deprecation shim itself can be tested. Leave it \code{NULL}.
#' @param srr_opt Beverton-Holt settings, required when
#'   \code{recruitment_opt = "bh_rec"}.
#'
#' @return The list from \code{Do_Population_Projection}.
#'
#' @keywords internal
project_at_F <- function(
  f,
  n_proj_yrs = 300,
  recruitment_opt = "mean_rec",
  srr_opt = NULL,
  data = NULL,
  rep = NULL,
  rec_window = 45L,
  rec_yrs = NULL,
  bh_rec_opt = NULL,
  rec_devs = NULL,
  fmort_opt = "Input",
  HCR_function = function(x, frp, brp, alpha = 0.05) frp
) {

  d <- if(is.null(data)) sgl_rg_sable_data else data
  rp <- if(is.null(rep)) sgl_rg_sable_rep else rep
  n_yrs <- length(d$years)
  n_regions <- d$n_regions; n_ages <- length(d$ages); n_sexes <- d$n_sexes
  n_fish_fleets <- d$n_fish_fleets; n_seas <- d$n_seas; n_pop <- d$n_pop

  biol_d <- c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)

  # A slice taken with drop = TRUE collapses whichever dims happen to be
  # length one, and which those are depends on the model. Keeping every dim
  # and flattening makes the copy positional, so the same code holds a one-season
  # single-sex model and a seasonal sexed one.
  hold_last_year <- function(arr, dims) {
    idx <- rep(list(bquote()), length(dim(arr)))
    idx[[3]] <- n_yrs
    slice <- as.vector(do.call(`[`, c(list(arr), idx, list(drop = FALSE))))
    out <- array(0, dim = dims)
    n_per_year <- prod(dims) / n_proj_yrs
    # a source without a fleet dim (weight at age, which the packaged data
    # has once rather than per fleet) is repeated across the fleets
    if(length(slice) && n_per_year %% length(slice) == 0)
      slice <- rep(slice, n_per_year / length(slice))
    if(length(slice) != n_per_year)
      stop("terminal-year slice has ", length(slice), " cells where the projection wants ",
           n_per_year, "; the array layouts have diverged")
    for(y in seq_len(n_proj_yrs)) {
      i <- rep(list(bquote()), length(dims)); i[[3]] <- y
      # pinned to base: while a tape is being built RTMB's replacement operator is
      # in scope, and it does not take the empty-symbol form do.call needs here
      out <- do.call(base::`[<-`, c(list(out), i, list(slice)))
    }
    out
  }

  WAA <- hold_last_year(d$WAA, biol_d)
  MatAA <- hold_last_year(d$MatAA, biol_d)
  # weight at age in the fishery falls back to the population's where a model
  # does not have a separate one
  WAA_fish <- hold_last_year(if(is.null(dim(d$WAA_fish))) d$WAA else d$WAA_fish,
                             c(biol_d, n_fish_fleets))
  fish_sel <- hold_last_year(rp$fish_sel, c(biol_d, n_fish_fleets))
  ret_sel <- hold_last_year(rp$ret_sel, c(biol_d, n_fish_fleets))

  # M has the same layout as WAA, seasons included
  natmort <- hold_last_year(rp$natmort, biol_d)

  Do_Population_Projection(
    n_proj_yrs = n_proj_yrs, n_regions = n_regions, n_ages = n_ages,
    n_sexes = n_sexes, n_pop = n_pop, n_fish_fleets = n_fish_fleets,
    # the seasonal structure has to be handed over explicitly: left out, the
    # projection assumes a single season and rejects the arrays built for more
    n_seas = n_seas, seasdur = d$seasdur, spawn_seas = d$spawn_seas,
    natal_region = d$natal_region,
    sexratio = array(0.5, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)),
    do_recruits_move = 0,
    recruitment = local({
      # the years mean recruitment is taken over. Written as a window back from
      # the terminal year rather than a fixed start, so a short series does not
      # produce a decreasing sequence and silently index backwards
      yrs <- if(!is.null(rec_yrs)) rec_yrs else
        seq.int(max(1L, n_yrs - rec_window), max(1L, n_yrs - 2L))
      array(rp$Rec[, , yrs], dim = c(n_pop, n_regions, length(yrs)))
    }),
    terminal_NAA = array(rp$NAA[, , n_yrs, , , ],
                         dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
    terminal_NAA0 = array(rp$NAA0[, , n_yrs, , , ],
                          dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
    terminal_F = array(rp$Fmort[, n_yrs, , ], dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(rp$dmr[, n_yrs, , ], dim = c(n_regions, n_seas, n_fish_fleets)),
    natmort = natmort, WAA = WAA, WAA_fish = WAA_fish, MatAA = MatAA,
    fish_sel = fish_sel, ret_sel = ret_sel,
    Movement = array(1, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)),
    # under fmort_opt = "Input" the projection fishes at f_ref_pt directly, so
    # this is how a constant rate is applied
    f_ref_pt = array(f, dim = c(n_regions, n_proj_yrs)),
    # a rule reads this, so it has to be a real reference point under fmort_opt
    # "HCR" and is left at zero for the constant-F path the reference point
    # tests use
    b_ref_pt = array(if(fmort_opt == "Input") 0 else sum(rp$SSB[,,n_yrs]) * 0.4 / n_pop,
                     dim = c(n_pop, n_regions, n_proj_yrs)),
    HCR_function = HCR_function,
    recruitment_opt = recruitment_opt, fmort_opt = fmort_opt,
    rec_devs = rec_devs,
    t_spawn = 0, srr_opt = srr_opt, bh_rec_opt = bh_rec_opt)
}

#' Project the packaged sablefish stock under a harvest control rule
#'
#' Wraps \code{project_at_F} in \code{fmort_opt = "HCR"} so a test can hand over
#' a rule and nothing else.
#'
#' @param rule The control rule, taking \code{x}, \code{frp}, \code{brp} and
#'   optionally \code{state}.
#' @param f The F reference point the rule is given.
#' @param n_proj_yrs Years to project.
#' @keywords internal
project_with_HCR <- function(rule, f = 0.06, n_proj_yrs = 8) {
  project_at_F(f = f, n_proj_yrs = n_proj_yrs,
               fmort_opt = "HCR", HCR_function = rule)
}

# Projected quantities are laid out with year on the third dim. The final
# projection year is read one year short of the end: the F rules set the next
# year's rate at the end of the current one, so the last slot is never fished.
proj_year_total <- function(arr, offset = 1) {
  y <- dim(arr)[3] - offset
  idx <- rep(list(bquote()), length(dim(arr)))
  idx[[3]] <- y
  sum(do.call(`[`, c(list(arr), idx, list(drop = FALSE))))
}

#' Equilibrium catch reached by a constant-F projection
#'
#' @keywords internal
equilibrium_catch <- function(out) proj_year_total(out$proj_Catch)

#' Equilibrium spawning biomass reached by a constant-F projection
#'
#' @keywords internal
equilibrium_ssb <- function(out) proj_year_total(out$proj_SSB)
