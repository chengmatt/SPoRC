# Stage 1 of 3: model setup
#
# Recruitment inputs: stock recruit model, steepness, recruitment and initial age deviations, and the
# spatial, seasonal and sex apportionment of recruits. Setup_Mod_Rec estimates, Setup_Sim_Rec simulates.

#' Set up recruitment dynamics for the operating model simulation
#'
#' Sets the stock-recruit form and density dependence, \eqn{R_0}, steepness, sex
#' ratio, the recruitment and initial age deviations, seasonal allocation, spawn
#' timing and the equilibrium initialization. Call after \code{\link{Setup_Sim_Dim}}.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}}.
#' @param recruitment_opt Recruitment model, default \code{"bh_rec"}.
#'   \code{0}/\code{"mean_rec"} has no stock-recruit relationship,
#'   \code{1}/\code{"bh_rec"} is Beverton-Holt and requires \code{rec_dd = "local"}
#'   when \code{n_pop > 1}, and \code{999}/\code{"resample_from_input"} resamples
#'   years from \code{Rec_input}, using historical years as they are and sampling
#'   projection years with replacement so spatial covariance within a year is kept.
#' @param rec_dd Density dependence for the stock-recruit relationship, default
#'   \code{"global"}. \code{0}/\code{"local"} gives each region its own \eqn{R_0} and
#'   steepness, required when \code{n_pop > 1} under \code{"bh_rec"};
#'   \code{1}/\code{"global"} pools across regions.
#' @param init_dd Density dependence for equilibrium initialization, same options as
#'   \code{rec_dd}. Default \code{"global"}.
#' @param R0_input Unfished equilibrium recruitment array \code{[n_pop x n_regions x
#'   n_yrs x n_sims]}. Default \code{15}.
#' @param h_input Steepness array \code{[n_pop x n_regions x n_yrs x n_sims]}, values
#'   in \eqn{(0.2, 1)}. Default \code{0.8}.
#' @param sexratio_input Proportion of recruits per sex, array \code{[n_pop x
#'   n_regions x n_yrs x n_sexes x n_sims]}. Default \code{1} for one sex, \code{0.5}
#'   each for two.
#' @param RecDevs_model Process error the deviations are drawn under. \code{"iid"}
#'   (default) is independent, \code{"rw"} a random walk, \code{"ar1"} reverts toward
#'   zero at rate \code{RecDevs_rho}. Year one is drawn at \code{ln_sigmaR} under the
#'   first two and from \code{ln_sigmaR / sqrt(1 - RecDevs_rho^2)} under \code{"ar1"}.
#'   Matches \code{\link{Setup_Mod_Rec}}.
#' @param RecDevs_rho Matrix \code{[n_pop x n_regions]} of AR1 correlations in
#'   \eqn{(-1, 1)}. Only read under \code{RecDevs_model = "ar1"}. Default zero.
#' @param ln_sigmaR Log-scale sd of the recruitment deviations, array \code{[2 x n_pop
#'   x n_regions]}, index 1 for \code{ln_InitDevs} and 2 for \code{ln_RecDevs}.
#'   Default \code{log(1)}.
#' @param rec_seas_prop_input Seasonal allocation of annual recruitment, array
#'   \code{[n_pop x n_seas x n_sims]} summing to 1 across seasons. Default all in
#'   season 1. Must be zero before \code{spawn_seas} when \code{rec_lag = 0} and
#'   \code{spawn_seas > 1}.
#' @param spawn_seas Integer index of the spawning season. Default \code{1}.
#' @param use_rinit Integer (0/1). Whether \code{rinit_input} initializes the
#'   population separately from \code{R0_input}. Under \code{0} (default)
#'   \code{rinit_input} is ignored.
#' @param rinit_input Equilibrium recruitment used for initialization under
#'   \code{use_rinit = 1}, array \code{[n_pop x n_regions x n_sims]}. Default \code{15}.
#' @param t_spawn Spawn timing as a fraction of \code{spawn_seas} elapsed before
#'   spawning. \code{0} (default) spawns before any mortality, \code{1} after all of it.
#' @param rec_lag Integer seasons between spawning and recruitment. \code{1} (default)
#'   uses SSB from that many seasons prior, in any season. \code{0} is age-0
#'   recruitment on the same year's SSB, so recruits may only enter in
#'   \code{spawn_seas} or later and \code{rec_seas_prop_input} must be zero before it.
#' @param init_age_strc Equilibrium initialization method, default \code{2}.
#'   \code{0}/\code{"iterative"} iterates forward to approximate equilibrium,
#'   \code{1}/\code{"scalar_no_move"} is a scalar geometric series without movement,
#'   \code{2}/\code{"matrix"} is the matrix series with movement, and
#'   \code{3}/\code{"scalar_plus_only"} moves only the plus group.
#'   \code{4}/\code{"free"} projects no equilibrium: \code{ln_InitDevs} are the initial
#'   log numbers-at-age for ages 2 and above, apportioned by sex ratio, so the initial
#'   condition ignores \code{init_F_par} and \code{ln_rinit} and any penalty becomes a
#'   prior on initial abundance.
#' @param do_recruits_move Integer flag. \code{0} (default) starts movement at age 2,
#'   \code{1} moves recruits from age 1.
#' @param stray_rate_input Natal-homing stray rate array \code{[n_pop x n_yrs x
#'   n_sims]}, the proportion straying from the natal region at spawning. Default
#'   \code{0}.
#' @param Rec_input External recruitment array \code{[n_pop x n_regions x n_yrs x
#'   n_sims]}, required under \code{recruitment_opt = "resample_from_input"}.
#'   Projection years beyond its length are resampled from historical years with
#'   replacement. Default \code{NULL}.
#' @param InitDevs_sex_spec How initial age deviations are drawn across sexes when
#'   \code{ln_InitDevs_input} is not supplied. \code{"est_shared_s"} (default) draws
#'   one curve per population or region for every sex, \code{"est_all"} draws each sex
#'   its own. Names match \code{\link{Setup_Mod_Rec}}.
#' @param ln_InitDevs_input Optional log-scale initial age deviations, either
#'   \code{[n_pop x n_regions x (n_ages - 1) x n_sims]} for one shared curve or
#'   \code{[n_pop x n_regions x (n_ages - 1) x n_sexes x n_sims]} for one per sex. The
#'   \code{n_ages - 1} dim excludes the reference age. \code{NULL} (default) draws one
#'   shared curve per population and region; pass zeros to start in equilibrium.
#' @param rec_bias_correct Integer. \code{1} (default) draws the recruitment and
#'   initial age deviations as mean-one lognormal multipliers centered at
#'   \eqn{-\sigma^2/2}, matching an estimation model with the bias ramp on; \code{0}
#'   centers them at zero. A linked cell follows the same switch under the arrows.
#' @param SR_ref_yr Integer year index supplying the biological inputs to unfished
#'   spawning biomass per recruit, and so the curve's scale. Matches the estimation
#'   model's \code{SR_ref_yr}. Default \code{1}.
#'
#' @return \code{sim_list} with \code{$recruitment_opt}, \code{$rec_dd},
#'   \code{$init_dd}, \code{$R0}, \code{$h}, \code{$sexratio}, \code{$ln_sigmaR},
#'   \code{$rec_seas_prop}, \code{$spawn_seas}, \code{$t_spawn}, \code{$rec_lag},
#'   \code{$init_age_strc}, \code{$do_recruits_move}, \code{$move_age},
#'   \code{$stray_rate}, and optionally \code{$Rec_input} and
#'   \code{$ln_InitDevs_input}. Character codes are converted to integers before
#'   storage.
#'
#' @export Setup_Sim_Rec
#' @family Simulation Setup
Setup_Sim_Rec <- function(
  sim_list,
  do_recruits_move = 0,
  sexratio_input = array(if(sim_list$n_sexes == 1) 1 else 0.5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_sims)),
  R0_input = array(15, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
  rinit_input = array(15, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims)),
  h_input = array(0.8, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
  stray_rate_input = array(0, dim = c(sim_list$n_pop, sim_list$n_yrs, sim_list$n_sims)),
  ln_sigmaR = array(log(1), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
  rec_seas_prop_input = {
    rec_seas_prop = array(0, dim = c(sim_list$n_pop, sim_list$n_seas, sim_list$n_sims))
    rec_seas_prop[, 1, ] <- 1
    rec_seas_prop
  },
  recruitment_opt = 'bh_rec',
  rec_dd = 'global',
  init_dd = 'global',
  use_rinit = 0,
  init_age_strc = 2,
  spawn_seas = 1,
  t_spawn = 0,
  rec_lag = 1,
  SR_ref_yr = 1,
  Rec_input = NULL,
  ln_InitDevs_input = NULL,
  InitDevs_sex_spec = "est_shared_s",
  RecDevs_model = "iid",
  RecDevs_rho = array(0, dim = c(sim_list$n_pop, sim_list$n_regions)),
  rec_bias_correct = 1
) {

  if(rec_dd == 'global' && sim_list$n_pop > 1 && recruitment_opt == 'bh_rec') stop("Invalid recruitment density-dependence option! When n_pop > 1 and recruitment_opt == 'bh_rec', rec_dd must be local (0).")
  if(rec_lag < 0) stop("rec_lag cannot be negative!")

  # Convert character inputs to numeric codes
  recruitment_opt <- convert_to_numeric(recruitment_opt, list(mean_rec = 0, bh_rec = 1, ricker_rec = 2, resample_from_input = 999))
  rec_dd <- convert_to_numeric(rec_dd, list(local = 0, global = 1))
  init_dd <- convert_to_numeric(init_dd, list(local = 0, global = 1))
  init_age_strc <- convert_to_numeric(init_age_strc, list(iterative = 0, scalar_no_move = 1, matrix = 2, scalar_plus_only = 3, free = 4))

  check_sim_dimensions(
    sexratio_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_sexes = sim_list$n_sexes,
    n_sims = sim_list$n_sims,
    n_pop = sim_list$n_pop,
    what = "sexratio_input"
  )
  check_sim_dimensions(
    R0_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_sims = sim_list$n_sims,
    n_pop = sim_list$n_pop,
    what = "R0_input"
  )
  check_sim_dimensions(
    rinit_input,
    n_regions = sim_list$n_regions,
    n_sims = sim_list$n_sims,
    n_pop = sim_list$n_pop,
    what = "rinit_input"
  )
  check_sim_dimensions(
    h_input,
    n_regions = sim_list$n_regions,
    n_years = sim_list$n_yrs,
    n_sims  = sim_list$n_sims,
    n_pop = sim_list$n_pop,
    what = "h_input"
  )
  check_sim_dimensions(
    stray_rate_input,
    n_years = sim_list$n_yrs,
    n_sims  = sim_list$n_sims,
    n_pop = sim_list$n_pop,
    what = "stray_rate_input"
  )
  check_sim_dimensions(
    rec_seas_prop_input,
    n_seas = sim_list$n_seas,
    n_sims  = sim_list$n_sims,
    n_pop = sim_list$n_pop,
    what = "rec_seas_prop_input"
  )
  if(!is.null(ln_InitDevs_input)) {
    check_sim_dimensions(
      ln_InitDevs_input,
      n_regions = sim_list$n_regions,
      n_ages = sim_list$n_ages,
      n_sexes = sim_list$n_sexes,
      n_sims = sim_list$n_sims,
      n_pop = sim_list$n_pop,
      what = "ln_InitDevs_input"
    )
    # one shared curve (no sex dimension) broadcasts across sexes
    if(length(dim(ln_InitDevs_input)) == 4) {
      tmp_init_input <- array(0, dim = c(dim(ln_InitDevs_input)[1:3], sim_list$n_sexes, sim_list$n_sims))
      for(s in 1:sim_list$n_sexes) tmp_init_input[,,,s,] <- ln_InitDevs_input
      ln_InitDevs_input <- tmp_init_input
    }
  }

  # age-0 recruits produced by this year's spawning cannot appear before spawn_seas in the same
  # year, since SSB and recruitment are not known yet at that point
  if(rec_lag == 0 && spawn_seas > 1 && any(rec_seas_prop_input[, seq_len(spawn_seas - 1), , drop = FALSE] != 0)) {
    stop("rec_lag = 0 requires rec_seas_prop_input to be zero in every season before spawn_seas (age-0 recruits can't predate the spawning event that produced them).")
  }

  # Recruitment options
  sim_list$do_recruits_move <- do_recruits_move
  if(sim_list$do_recruits_move == 0) sim_list$move_age <- 2 else sim_list$move_age <- 1 # what age to start movement of individuals

  if(recruitment_opt == 999) {
    if(is.null(Rec_input)) stop("Recruitment input is NULL, but future recruitment is specified to be resampled!")
    rec_input_yrs <- dim(Rec_input)[3] # get years from Rec_input
    tmp_Rec_input <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
    # loop through simulations to resample years
    for(i in 1:sim_list$n_sims) {
      tmp_Rec_input[,,1:rec_input_yrs,i] <- Rec_input[,,,i]
      resampled_years <- sample(1:rec_input_yrs, length(tmp_Rec_input[1,1,-c(1:rec_input_yrs),i]), TRUE)
      tmp_Rec_input[,,-c(1:rec_input_yrs),i] <- Rec_input[,,resampled_years,i]
    } # end i loop
    Rec_input <- tmp_Rec_input # overwrite
  } # resampling

  # Output recruitment stuff into environment
  sim_list$recruitment_opt <- recruitment_opt
  sim_list$rec_dd <- rec_dd
  sim_list$init_dd <- init_dd
  sim_list$h <- h_input
  sim_list$R0 <- R0_input
  sim_list$use_rinit <- use_rinit
  sim_list$rinit <- rinit_input
  sim_list$sexratio <- sexratio_input
  sim_list$rec_lag <- rec_lag
  # every input to unfished spawning biomass per recruit is taken at this year, matching
  # the estimation model's SR_ref_yr so a self test compares like with like
  if(length(SR_ref_yr) != 1 || SR_ref_yr < 1 || SR_ref_yr > sim_list$n_yrs)
    stop("SR_ref_yr must be a single year index between 1 and ", sim_list$n_yrs, ".")
  sim_list$SR_ref_yr <- as.integer(SR_ref_yr)
  sim_list$ln_sigmaR <- ln_sigmaR
  if(!RecDevs_model %in% c("iid", "rw", "ar1")) stop("RecDevs_model incorrectly specified. Must be one of 'iid', 'rw', or 'ar1'")
  sim_list$RecDevs_model <- match(RecDevs_model, c("iid", "rw", "ar1")) # 1 = iid, 2 = rw, 3 = ar1
  if(RecDevs_model == "ar1" && any(abs(RecDevs_rho) >= 1)) stop("RecDevs_rho must be inside (-1, 1). An ar1 at 1 or beyond has no stationary variance, so year one of the series has nothing to be drawn from.")
  sim_list$RecDevs_rho <- array(RecDevs_rho, dim = c(sim_list$n_pop, sim_list$n_regions))
  sim_list$t_spawn <- t_spawn
  sim_list$rec_seas_prop <- rec_seas_prop_input
  sim_list$init_age_strc <- init_age_strc
  if(!InitDevs_sex_spec %in% c("est_shared_s", "est_all")) stop("InitDevs_sex_spec must be est_shared_s or est_all")
  if(InitDevs_sex_spec == "est_all" && sim_list$n_sexes == 1) stop("InitDevs_sex_spec = 'est_all' draws a curve per sex, so it needs n_sexes > 1 (with one sex, est_shared_s is the same thing)")
  sim_list$InitDevs_sex_spec <- InitDevs_sex_spec
  sim_list$spawn_seas <- spawn_seas
  sim_list$stray_rate <- stray_rate_input
  if(!is.null(Rec_input)) sim_list$Rec_input <- Rec_input
  if(!is.null(ln_InitDevs_input)) sim_list$ln_InitDevs_input <- ln_InitDevs_input


  # recruitment bias correction stuff
  if(!rec_bias_correct %in% c(0, 1)) stop("rec_bias_correct must be 0 or 1")
  sim_list$rec_bias_correct <- rec_bias_correct
  if(is.null(ln_InitDevs_input)) message("Setup_Sim_Rec: ln_InitDevs_input is NULL, so every replicate draws its own initial age deviations from N(",
                                         if(rec_bias_correct == 1) "-sigma^2/2" else "0", ", sigma) at the early ln_sigmaR. Pass zeros for a population that starts in equilibrium, or a fit's deviations to condition on it.")

  return(sim_list)

}

#' Map recruitment variability (sigma_R) parameters
#'
#' Builds the factor map for \code{ln_sigmaR} \code{[2 x n_pop x n_regions]}, where
#' the first index is the initial deviation period and the second the annual one.
#' The recruitment penalty reads each region's own slot, and a population's
#' non-natal region slots are mapped off under \code{rec_region_prop_spec = 1} with
#' several populations. Called by \code{\link{Setup_Mod_Rec}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#' @param sigmaR_spec Estimation structure for \code{ln_sigmaR}. \code{"est_all"}
#'   estimates each period, population and region separately. \code{"est_shared_r"}
#'   estimates per period and population, shared across regions.
#'   \code{"est_shared_all"} gives one value for everything.
#'   \code{"fix_early_est_late"} holds the initial period at its starting value and
#'   estimates the annual one per population and region. \code{"fix"} holds every
#'   value.
#'
#' @return \code{input_list} with \code{$map$ln_sigmaR} set to a factor vector of
#'   length \code{prod(dim(par$ln_sigmaR))}. Active parameters take sequential
#'   integers, fixed ones are \code{NA}.
#'
#' @keywords internal
do_sigmaR_mapping <- function(input_list, sigmaR_spec) {

  # Define valid sigmaR options
  valid_options <- c("est_all", "est_shared_r", "est_shared_all", "fix_early_est_late", "fix")

  # Checking to see if valid options
  if (!sigmaR_spec %in% valid_options) stop("Invalid sigmaR_spec. Must be one of: ", paste(valid_options, collapse = ", "))

  dims <- c(period = 2, pop = input_list$data$n_pop, region = input_list$data$n_regions)

  if(sigmaR_spec == 'fix') {
    input_list$map$ln_sigmaR <- factor(rep(NA, prod(dims)))
  } else if(sigmaR_spec == 'est_shared_all') {
    input_list$map$ln_sigmaR <- factor(rep(1, prod(dims)))
  } else {
    # "est_all" and "fix_early_est_late" resolve every period, population and
    # region; "est_shared_r" shares across regions within a period and population
    share_over <- if(sigmaR_spec == 'est_shared_r') "region" else character(0)
    map_sigmaR <- build_pe_map(dims, share_over = share_over)

    # a population confined to its natal region has no deviations elsewhere, so
    # those slots have no information and stay off
    if(isTRUE(input_list$data$rec_region_prop_spec == 1) && input_list$data$n_pop > 1) {
      for(p in 1:input_list$data$n_pop) map_sigmaR[, p, -input_list$data$natal_region[p]] <- NA
    } # end if

    if(sigmaR_spec == 'fix_early_est_late') map_sigmaR[1,,] <- NA # fix early period, keep late period estimated
    input_list$map$ln_sigmaR <- factor(as.vector(map_sigmaR))
  }

  collect_message("Recruitment Variability is specified as: ", sigmaR_spec)

  return(input_list)
}

#' The recruitment bias ramp, year by year
#'
#' The lognormal bias correction a recruitment deviation's penalty is centered
#' on, \eqn{-b_t \sigma^2 / 2}, as a factor \eqn{b_t} per estimated year.
#' \code{do_rec_bias_ramp = 0} gives the full correction in every year
#' (\eqn{b_t = 1}); \code{1} ramps it up, holds it and ramps it down over the
#' four \code{bias_year} indices, scaled by \code{max_bias_ramp_fct}, so a
#' ramp whose years all sit at the last year is zero everywhere. Used by the
#' objective and by the setup checks that ask whether the correction touches a
#' given year.
#'
#' @param do_rec_bias_ramp Integer, \code{0} or \code{1}.
#' @param bias_year Integer vector of the four ramp years, as deviation indices.
#' @param n_est_rec_devs Number of estimated recruitment deviation years.
#' @param max_bias_ramp_fct Scale of the ramp at its plateau.
#'
#' @return Numeric vector of length \code{n_est_rec_devs}.
#'
#' @keywords internal
get_rec_bias_ramp <- function(do_rec_bias_ramp, bias_year, n_est_rec_devs, max_bias_ramp_fct = 1) {

  if(do_rec_bias_ramp == 0) return(rep(1, n_est_rec_devs)) # the full correction every year

  ramp_yrs <- 1:n_est_rec_devs
  bias_ramp <- rep(0, n_est_rec_devs)
  range1 <- which(ramp_yrs >= bias_year[1] & ramp_yrs < bias_year[2]) # ascending limb
  range2 <- which(ramp_yrs >= bias_year[2] & ramp_yrs < bias_year[3]) # full correction
  range3 <- which(ramp_yrs >= bias_year[3] & ramp_yrs < bias_year[4]) # descending limb
  if(length(range1) > 0) bias_ramp[range1] <- (ramp_yrs[range1] - bias_year[1]) / (bias_year[2] - bias_year[1])
  if(length(range2) > 0) bias_ramp[range2] <- 1
  if(length(range3) > 0) bias_ramp[range3] <- 1 - ((ramp_yrs[range3] - bias_year[3]) / (bias_year[4] - bias_year[3]))

  return(bias_ramp * max_bias_ramp_fct)

} # end function

#' Map initial age-structure deviation parameters
#'
#' Builds the factor map for \code{ln_InitDevs} \code{[n_pop x n_regions x
#' (n_ages - 1) x n_sexes]}, the log-scale deviations from the equilibrium initial
#' age structure. Population, region and age are resolved on a single-sex slice and
#' then expanded across sexes by \code{InitDevs_sex_spec}. Called by
#' \code{\link{Setup_Mod_Rec}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#'   Requires \code{$data$equil_init_age_strc}, \code{$data$rec_region_prop_spec},
#'   \code{$data$natal_region} and \code{$data$rec_dd}.
#' @param InitDevs_spec Sharing structure for \code{ln_InitDevs}.
#'   \code{"est_shared_pop_r"} gives one set of age deviations across every
#'   population and region, required when \code{rec_dd = "global"} and
#'   \code{n_regions > 1}. \code{"est_shared_r"} gives one set per population,
#'   shared across its regions. \code{"fix"} holds every deviation at zero.
#'   \code{NULL} (default) estimates all independently, which is not permitted when
#'   \code{rec_region_prop_spec = 1} and \code{n_pop > 1}.
#' @param rec_dd Density dependence inherited from \code{\link{Setup_Mod_Rec}}.
#'   \code{"global"} restricts \code{InitDevs_spec} to \code{"est_shared_r"} or
#'   \code{"est_shared_pop_r"} when \code{n_regions > 1}.
#' @param init_age_devs_shared Integer vector of length \code{n_ages - 1} giving the
#'   factor level of each age position; positions sharing a value share one parameter.
#'   Read under \code{equil_init_age_strc = 3}, and respected by
#'   \code{InitDevs_spec = "est_shared_r"} (per population, with a population offset)
#'   and \code{"est_shared_pop_r"} (globally, no offset). \code{c(1:42, rep(42, 9))}
#'   gives 42 free parameters for a 52-age model with 43 data ages. Default
#'   \code{NULL}.
#' @param InitDevs_sex_spec \code{"est_shared_s"} (default) maps every sex onto one
#'   age curve penalized once; \code{"est_all"} offsets the factor levels per sex so
#'   each has its own, each penalized. Also builds \code{data$init_devs_pen_use},
#'   which flags one penalized copy of every estimated parameter.
#'
#' @return \code{input_list} with \code{$map$ln_InitDevs} set to a factor vector of
#'   length \code{prod(dim(par$ln_InitDevs))}. Active parameters take sequential
#'   integers; plus-group slots under \code{equil_init_age_strc = 1} and non-natal
#'   region slots under \code{rec_region_prop_spec = 1} are \code{NA}, and their
#'   starting values are reset to \code{0}.
#'
#' @keywords internal
do_InitDevs_mapping <- function(input_list, InitDevs_spec, rec_dd, init_age_devs_shared, InitDevs_sex_spec = "est_shared_s") {

  # validate if stoch_shared_ages
  if(input_list$data$equil_init_age_strc == 3) {
    if(is.null(init_age_devs_shared)) stop("init_age_devs_shared is NULL, but equil_init_age_strc is stoch_shared_ages!")
    if(length(init_age_devs_shared) != (length(input_list$data$ages) - 1)) stop("init_age_devs_shared must have length n_ages - 1 = ", n_age_dim, " but has length ", length(init_age_devs_shared))
  }

  # sexes either share one age curve (est_shared_s) or have their own (est_all). the logic below is
  # written on a single-sex slice, so the sex dim is peeled off here and reattached at the end
  if(!InitDevs_sex_spec %in% c("est_shared_s", "est_all")) stop("InitDevs_sex_spec must be est_shared_s or est_all")
  if(InitDevs_sex_spec == "est_all" && input_list$data$n_sexes == 1) stop("InitDevs_sex_spec = 'est_all' estimates a curve per sex, so it requires n_sexes > 1 (with one sex, est_shared_s is the same thing)")
  par_full <- input_list$par$ln_InitDevs
  input_list$par$ln_InitDevs <- array(par_full[,,,1], dim = dim(par_full)[1:3])

  # code 4 estimates the same cells code 2 does and penalizes none of them, so the two share a mapping and only get_recruitment_penalty tells them apart
  est_all_ages <- input_list$data$equil_init_age_strc %in% c(2, 4)
  if(est_all_ages && isTRUE(input_list$data$use_rinit == 1)) collect_message("use_rinit = 1 with a deviation estimated on every initial age: ln_rinit and the level of ln_InitDevs are separated only by the initial deviation penalty (ln_sigmaR[1]). Where the initial structure is known to be in equilibrium, map ln_InitDevs off at zero.")
  all_ages_msg <- if(input_list$data$equil_init_age_strc == 4)
    "Initial age deviations are estimated for all ages including the plus group, and none of them are penalized."
  else "Initial age deviations are stochastic and estimated for all ages, including the plus group"

  # Initial age deviations (equilibrium)
  if(input_list$data$equil_init_age_strc == 0) {
    input_list$par$ln_InitDevs <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$ages) - 1)) # override starting values if previously specified
    input_list$map$ln_InitDevs <- factor(rep(NA, length(input_list$par$ln_InitDevs))) # set mapping
    collect_message("Initial Age Structure is specified to be in equilibrium. No initial age deviations are estimated.")
  }

  # Initial age deviations (stochastic for all ages, including plus group)
  if(!is.null(InitDevs_spec)) {

    # Validate options
    if(!is.null(rec_dd) && rec_dd == 'global' && !InitDevs_spec %in% c("est_shared_r", "est_shared_pop_r") && input_list$data$n_regions > 1) {
      stop("Please specify a valid initial age deviations option for global recruitment density dependence (should be est_shared_r or est_shared_pop_r)!")
    }

    if(!InitDevs_spec %in% c("est_shared_pop_r", "est_shared_r", "fix")) stop("Please specify a valid initial deviations option. These include: fix, est_shared_r, est_shared_pop_r. Conversely, leave at NULL to estimate all initial deviations.")
    else collect_message("Initial Deviations is stochastic and specified as: ", InitDevs_spec)

    # set up mapping for initial age deviations
    map_InitDevs <- input_list$par$ln_InitDevs

    # Fix all initial deviations
    if(InitDevs_spec == "fix") input_list$map$ln_InitDevs <- factor(rep(NA, prod(dim(map_InitDevs))))

    # Share across regions and populations
    if(InitDevs_spec == 'est_shared_pop_r') {

      # share parameters, but no stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 1) {

        # get indices
        n_ages <- dim(input_list$par$ln_InitDevs)[3] - 1
        # each populaiton and reigon has the same indices
        map_InitDevs[,,-dim(input_list$par$ln_InitDevs)[3]] <- rep(1:n_ages, each = input_list$data$n_regions * input_list$data$n_pop)
        map_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- NA  # NA for plus group
        input_list$par$ln_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- 0
        collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
      }

      # share parameters across regions, with deviations on the plus group as well
      if(est_all_ages) {

        # get indices
        n_ages_all <- dim(input_list$par$ln_InitDevs)[3]
        map_InitDevs[] <- rep(1:n_ages_all, each = input_list$data$n_regions * input_list$data$n_pop)

        collect_message(all_ages_msg)
      }

      # share parameters across regions, with stochastic deviations on user defined ages
      if(input_list$data$equil_init_age_strc == 3) {
        # get indices
        for(p in 1:input_list$data$n_pop) for(r in 1:input_list$data$n_regions) map_InitDevs[p, r, ] <- init_age_devs_shared
        collect_message("Initial age deviations are stochastic for user-defined age sharing structure.")
      }

      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
    }

    # Share across regions and estimate for each population
    if(InitDevs_spec == "est_shared_r") {

      # share parameters, but no stochastic deviations on plus group
      if(input_list$data$equil_init_age_strc == 1) {

        # get indices
        n_ages <- dim(input_list$par$ln_InitDevs)[3] - 1
        n_region <- dim(input_list$par$ln_InitDevs)[2]

        for(p in 1:input_list$data$n_pop) {
          age_indices <- (1:n_ages) + (p - 1) * n_ages # get age indices
          # each region gets the same index for a given age (repeat each index across regions)
          map_InitDevs[p,,-dim(input_list$par$ln_InitDevs)[3]] <- matrix(rep(age_indices, each = n_region), nrow = n_region)
          map_InitDevs[p,,dim(input_list$par$ln_InitDevs)[3]] <- NA  # NA for plus group
          input_list$par$ln_InitDevs[p,,dim(input_list$par$ln_InitDevs)[3]] <- 0
        } # end p loop

        collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
      }

      # share parameters across regions, with deviations on the plus group as well
      if(est_all_ages) {

        # get indices
        n_ages_all <- dim(input_list$par$ln_InitDevs)[3]
        n_region <- dim(input_list$par$ln_InitDevs)[2]

        for(p in 1:input_list$data$n_pop) {
          age_indices <- (1:n_ages_all) + (p - 1) * n_ages_all
          # each region gets the same index for a given age (repeat each index across regions)
          map_InitDevs[p,,] <- matrix(rep(age_indices, each = n_region), nrow = n_region)
        } # end p loop

        collect_message(all_ages_msg)
      }

      # Share parameters across used-defined ages
      if(input_list$data$equil_init_age_strc == 3) {
        for(p in 1:input_list$data$n_pop) {
          pop_offset <- (p - 1) * max(init_age_devs_shared, na.rm = TRUE)
          for(r in 1:input_list$data$n_regions)
            map_InitDevs[p, r, ] <- ifelse(is.na(init_age_devs_shared), NA, init_age_devs_shared + pop_offset)
        }
        collect_message("Initial age deviations are stochastic for user-defined age sharing structure.")
      }

      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
    } # end if

  } else { # If NULL, then estimating age deviations across all dimensions

    if(input_list$data$n_pop > 1 && input_list$data$rec_region_prop_spec == 1)
      stop("Can't estimate initial age deviations for all populations and regions if no recruitment dispersal is occuring within a given region! Please specify est_shared_r or est_shared_pop_r instead!")

    map_InitDevs <- input_list$par$ln_InitDevs # set up mapping for initial age deviations

    if(input_list$data$equil_init_age_strc == 1) { # estimating all deviations across all dimensions, except for plus group
      map_InitDevs[,,-dim(input_list$par$ln_InitDevs)[3]] <- seq_along(map_InitDevs[,,-dim(input_list$par$ln_InitDevs)[3]]) # don't estimate plus group
      map_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- NA # NA for plus group
      input_list$par$ln_InitDevs[,,dim(input_list$par$ln_InitDevs)[3]] <- 0 # reset plus group starting value to 0
      input_list$map$ln_InitDevs <- factor(map_InitDevs) # input into map
      collect_message("Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.")
    }

    # Plus group and estimating deviations for all dimensions
    if(est_all_ages) {
      input_list$map$ln_InitDevs <- factor(seq_along(map_InitDevs)) # input into map
      collect_message(all_ages_msg)
    }

    # User-defined age sharing, estimated independently across all pops and regions
    if(input_list$data$equil_init_age_strc == 3) {
      max_age_idx <- max(init_age_devs_shared, na.rm = TRUE)
      for(p in 1:input_list$data$n_pop) {
        for(r in 1:input_list$data$n_regions) {
          pr_offset <- ((p - 1) * input_list$data$n_regions + (r - 1)) * max_age_idx
          map_InitDevs[p, r, ] <- ifelse(is.na(init_age_devs_shared), NA, init_age_devs_shared + pr_offset)
        }
      }
      input_list$map$ln_InitDevs <- factor(map_InitDevs)
      collect_message("Initial Age Deviations estimated independently per pop/region with user-defined age sharing structure.")
    }


  }

  # When no_dispersal, non-natal regions have no recruitment so their
  # deviations are structurally zero. Fix them regardless of InitDevs_spec.
  if(input_list$data$rec_region_prop_spec == 1 && input_list$data$n_pop > 1) {

    # extract mapping
    map_tmp <- as.integer(input_list$map$ln_InitDevs)
    dim(map_tmp) <- dim(input_list$par$ln_InitDevs)

    for(p in seq_len(input_list$data$n_pop)) {
      for(r in seq_len(input_list$data$n_regions)) {
        if(r != input_list$data$natal_region[p]) {
          input_list$par$ln_InitDevs[p, r, ] <- 0  # fix starting value
          map_tmp[p, r, ] <- NA                     # turn off estimation
        }
      }
    }

    # Re-index non-NA values sequentially (1, 2, 3, ...)
    non_na <- !is.na(map_tmp)
    map_tmp[non_na] <- as.integer(factor(map_tmp[non_na]))
    input_list$map$ln_InitDevs <- factor(map_tmp)
    collect_message("No dispersal: initial age deviations for non-natal regions fixed to 0 and not estimated.")
  }

  # reattach the sex dim. est_shared_s repeats the single-sex factor levels so every sex reads one
  # parameter set, est_all offsets them. init_devs_pen_use flags one penalized copy per parameter
  map3 <- as.integer(input_list$map$ln_InitDevs)
  dim(map3) <- dim(input_list$par$ln_InitDevs)
  par3 <- input_list$par$ln_InitDevs
  n_sexes_id <- input_list$data$n_sexes
  dims4 <- c(dim(par3), n_sexes_id)

  map4 <- array(NA_real_, dim = dims4)
  par4 <- array(0, dim = dims4)
  pen_use <- array(0, dim = dims4)
  n_free <- if(all(is.na(map3))) 0 else max(map3, na.rm = TRUE)
  for(s in 1:n_sexes_id) {
    if(InitDevs_sex_spec == "est_shared_s" || s == 1) map4[,,,s] <- map3
    else map4[,,,s] <- map3 + (s - 1) * n_free
    # under sharing every sex reads sex 1's slice, so starting values mirror it too
    par4[,,,s] <- if(InitDevs_sex_spec == "est_shared_s") par3 else array(par_full[,,,s], dim = dim(par3))
    par4[,,,s][is.na(map3)] <- par3[is.na(map3)] # fixed cells hold the (reset) single-sex value
    if(InitDevs_sex_spec == "est_all" || s == 1) pen_use[,,,s][!is.na(map3)] <- 1
  } # end s loop

  input_list$par$ln_InitDevs <- par4
  input_list$map$ln_InitDevs <- factor(map4)
  input_list$data$init_devs_pen_use <- pen_use
  if(InitDevs_sex_spec == "est_all") collect_message("Initial age deviations are sex-specific (one age curve per sex).")

  return(input_list)
}

#' Map annual recruitment deviation parameters
#'
#' Builds the factor map for \code{ln_RecDevs} \code{[n_pop x n_regions x
#' n_years]}, the log-scale annual recruitment deviations. When
#' \code{rec_region_prop_spec = 1} and \code{n_pop > 1}, non-natal regions get no
#' recruitment, so their deviations are fixed to \code{NA} whatever
#' \code{RecDevs_spec} asks for and the rest are re-numbered. Called by
#' \code{\link{Setup_Mod_Rec}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#'   Requires \code{$data$rec_region_prop_spec}, \code{$data$natal_region},
#'   \code{$data$rec_dd} and \code{$data$n_pop}.
#' @param RecDevs_spec Sharing structure for \code{ln_RecDevs}.
#'   \code{"est_shared_r"} gives one deviation series per population, shared across
#'   its regions. \code{"est_shared_pop_r"} gives one series across every
#'   population and region, required when \code{rec_dd = "global"} and
#'   \code{n_regions > 1}. \code{"fix"} holds every deviation at zero. \code{NULL}
#'   estimates all independently, which is not permitted when
#'   \code{rec_region_prop_spec = 1} and \code{n_pop > 1}.
#' @param rec_dd Density dependence inherited from \code{\link{Setup_Mod_Rec}}.
#'   \code{"global"} restricts \code{RecDevs_spec} to \code{"est_shared_r"} or
#'   \code{"est_shared_pop_r"} when \code{n_regions > 1}.
#'
#' @return \code{input_list} with \code{$map$ln_RecDevs} set to a factor vector of
#'   length \code{prod(dim(par$ln_RecDevs))}. Active parameters take sequential
#'   integers; non-natal region slots and fixed deviations are \code{NA}, and their
#'   starting values are reset to \code{0}.
#'
#' @seealso \code{\link{do_InitDevs_mapping}}, which shares the same options.
#'
#' @keywords internal
do_RecDevs_mapping <- function(input_list, RecDevs_spec, rec_dd, dont_pen_recdev_first = 0) {

  map_RecDevs <- input_list$par$ln_RecDevs # set up mapping for recruitment deviations

  # Recruitment deviations
  if(!is.null(RecDevs_spec)) {

    # Validate options
    if(!is.null(rec_dd) && rec_dd == 'global' && !RecDevs_spec %in% c("est_shared_r", "est_shared_pop_r") && input_list$data$n_regions > 1) stop("Please specify a valid recruitment deviations option for global recruitment density dependence (should be est_shared_r or est_shared_pop_r)!")
    if(!RecDevs_spec %in% c("est_shared_pop_r", "est_shared_r", "fix"))  stop("Please specify a valid recruitment deviations option. These include: fix, est_shared_r, est_shared_pop_r. Conversely, leave at NULL to estimate all recruitment deviations.")

    # Share across regions and estimate by population
    if(RecDevs_spec == "est_shared_r") {

      # get indices
      n_yrs <- dim(input_list$par$ln_RecDevs)[3]
      n_region <- dim(input_list$par$ln_RecDevs)[2]

      for(p in 1:input_list$data$n_pop) {
        yr_indices <- (1:n_yrs) + (p - 1) * n_yrs # get age indices
        # each region gets the same index for a given age (repeat each index across regions)
        map_RecDevs[p,,] <- matrix(rep(yr_indices, each = n_region), nrow = n_region)
      } # end p loop

      input_list$map$ln_RecDevs <- factor(map_RecDevs)
    } # end if

    # Share across regions and populations
    if(RecDevs_spec == 'est_shared_pop_r') {
      # get indices
      n_yrs <- dim(input_list$par$ln_RecDevs)[3]
      map_RecDevs[] <- rep(1:n_yrs, each = input_list$data$n_regions * input_list$data$n_pop)
      input_list$map$ln_RecDevs <- factor(map_RecDevs)
    }

    # Fix all recruitment deviations
    if(RecDevs_spec == "fix") input_list$map$ln_RecDevs <- factor(rep(NA, prod(dim(map_RecDevs))))

    # print message
    collect_message("Recruitment Deviations is specified as: ", RecDevs_spec)

  } else { # if NULL, estimating all dimensions

    if(input_list$data$n_pop > 1 && input_list$data$rec_region_prop_spec == 1)
      stop("Can't estimate recruitment eviations for all populations and regions if no recruitment dispersal is occuring within a given region! Please specify est_shared_r or est_shared_pop_r instead!")

    input_list$map$ln_RecDevs <- factor(seq_along(map_RecDevs)) # input into mapping
    collect_message("Recruitment Deviations is estimated for all dimensions")
  }

  # When no_dispersal, non-natal regions have no recruitment so their
  # deviations are structurally zero. Fix them regardless of RecDevs_spec.
  if(input_list$data$rec_region_prop_spec == 1 && input_list$data$n_pop > 1) {

    # extract mapping
    map_tmp <- as.integer(input_list$map$ln_RecDevs)
    dim(map_tmp) <- dim(input_list$par$ln_RecDevs)

    for(p in seq_len(input_list$data$n_pop)) {
      for(r in seq_len(input_list$data$n_regions)) {
        if(r != input_list$data$natal_region[p]) {
          input_list$par$ln_RecDevs[p, r, ] <- 0  # fix starting value
          map_tmp[p, r, ] <- NA                     # turn off estimation
        }
      }
    }

    # Re-index non-NA values sequentially (1, 2, 3, ...)
    non_na <- !is.na(map_tmp)
    map_tmp[non_na] <- as.integer(factor(map_tmp[non_na]))
    input_list$map$ln_RecDevs <- factor(map_tmp)
    collect_message("No dispersal: Recruitment deviations for non-natal regions fixed to 0 and not estimated.")
  }

  # mirror the deviation map into the data list so the recruitment penalty keys on the cells that
  # are estimated. a deviation mapped off by hand after setup is neither estimated nor penalized
  input_list$data$map_ln_RecDevs <- array(as.numeric(input_list$map$ln_RecDevs),
                                          dim = dim(input_list$par$ln_RecDevs))

  # the first years can belong to the initial condition rather than to the recruitment process.
  # dropping them from the mirror alone leaves them estimated but takes their penalty away
  n_dev_yrs <- dim(input_list$par$ln_RecDevs)[3]
  if(length(dont_pen_recdev_first) != 1 || is.na(dont_pen_recdev_first) || dont_pen_recdev_first %% 1 != 0 || dont_pen_recdev_first < 0)
    stop("dont_pen_recdev_first is '", paste(dont_pen_recdev_first, collapse = ", "), "'. Give a whole number of ",
         "leading years to leave out of the recruitment penalty, or 0 to penalize every year.")

  if(dont_pen_recdev_first >= n_dev_yrs)
    stop("dont_pen_recdev_first is ", dont_pen_recdev_first, " but there are only ", n_dev_yrs,
         " years of recruitment deviations. Leaving every year out would leave the deviations with ",
         "no process error at all, so at least one year has to stay in the penalty.")

  input_list$data$dont_pen_recdev_first <- as.integer(dont_pen_recdev_first)

  if(dont_pen_recdev_first > 0) {
    input_list$data$map_ln_RecDevs[,,seq_len(dont_pen_recdev_first)] <- NA
    collect_message("Recruitment deviations for the first ", dont_pen_recdev_first,
                    " year(s) are estimated but left out of the penalty, so they belong to the initial condition")
  }
  # do the initial age deviations too, so that deviations shared across regions
  # or sexes through the map split one penalty rather than being counted per cell
  if(!is.null(input_list$map$ln_InitDevs))
    input_list$data$map_ln_InitDevs <- array(as.numeric(input_list$map$ln_InitDevs),
                                             dim = dim(input_list$par$ln_InitDevs))

  return(input_list)
}

#' Map AR1 correlation parameter for recruitment deviations
#'
#' Constructs the \code{RecDevs_rho} factor map. \code{RecDevs_rho} is only
#' read when \code{RecDevs_model = "ar1"} (see \code{\link{Setup_Mod_Rec}});
#' under any other \code{RecDevs_model} every \code{RecDevs_rho} parameter is
#' mapped to \code{NA} regardless of \code{RecDevs_rho_spec}, since the
#' recruitment penalty never reads it.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param RecDevs_rho_spec Character string controlling the sharing and
#'   estimation structure for \code{RecDevs_rho}: one of \code{"est_all"},
#'   \code{"est_shared_pop"}, \code{"est_shared_r"},
#'   \code{"est_shared_pop_r"}, or \code{"fix"}.
#'
#' @return The input \code{input_list} with \code{$map$RecDevs_rho} set to a
#'   factor vector of length \code{prod(dim(par$RecDevs_rho))}.
#'
#' @keywords internal
do_RecDevs_rho_mapping <- function(input_list, RecDevs_rho_spec) {

  dims <- c(pop = input_list$data$n_pop,
            region = input_list$data$n_regions)

  if(input_list$data$RecDevs_model != 3) { # only ar1 reads RecDevs_rho
    input_list$map$RecDevs_rho <- factor(rep(NA, prod(dims)))
  } else {
    input_list$map$RecDevs_rho <- build_shared_spec_map(
      dims = dims,
      spec = RecDevs_rho_spec,
      dim_abbrev = c(pop = "pop", r = "region")
    )
  }

  collect_message("RecDevs_rho is specified as: ", RecDevs_rho_spec)

  return(input_list)
}

#' Map Beverton-Holt steepness parameters
#'
#' Builds the factor map for \code{steepness_h} \code{[n_pop x n_regions]}. Every
#' element is \code{NA} under \code{rec_model = 0}, where steepness has no role.
#' Called by \code{\link{Setup_Mod_Rec}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#'   Requires \code{$data$rec_model}, \code{$data$n_pop}, \code{$data$n_regions} and
#'   \code{$data$rec_dd}.
#' @param h_spec Sharing structure for \code{steepness_h}.
#'   \code{"est_shared_pop_r"} gives one value across every population and region,
#'   required when \code{rec_dd = "global"} and \code{n_regions > 1}.
#'   \code{"est_shared_r"} gives one per population. \code{"fix"} holds every value
#'   at its starting value. \code{NULL} estimates by population when
#'   \code{n_pop > 1} and by region when \code{n_pop = 1}, and is not permitted
#'   under global density dependence with \code{n_regions > 1}.
#' @param rec_dd Density dependence inherited from \code{\link{Setup_Mod_Rec}}.
#'   \code{"global"} restricts \code{h_spec} to \code{"est_shared_r"},
#'   \code{"est_shared_pop_r"} or \code{"fix"}.
#'
#' @return \code{input_list} with \code{$map$steepness_h} set to a factor vector of
#'   length \code{prod(dim(par$steepness_h))}. Active parameters take sequential
#'   integers, unused ones are \code{NA}.
#'
#' @keywords internal
do_h_mapping <- function(input_list, h_spec, rec_dd) {

  # Validate h_spec given rec_dd context
  if(input_list$data$rec_model != 0 && !is.null(rec_dd) && rec_dd == "global") {
    if(is.null(h_spec)) {
      stop("When rec_dd == `global` (global density dependence), h_spec cannot be NULL. ",
           "Steepness must be shared across the global SR relationship: use 'est_shared_pop_r' or 'est_shared_r', or 'fix'.")
    }
    if(!h_spec %in% c("est_shared_pop_r", "est_shared_r", "fix")) {
      stop("When rec_dd == `global` (global density dependence), h_spec must be ",
           "'est_shared_pop_r', 'est_shared_r', or 'fix'.")
    }
  }

  # mean recruitment with no stock-recruit penalty has no curve, so steepness stays off. with a
  # penalty there is a curve and steepness has to be reachable, so fall through to h_spec
  if(input_list$data$rec_model == 0 && input_list$data$sr_penalty == 0) {
    input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))
  } else if(!is.null(h_spec)) {

    # Validate options
    if(!h_spec %in% c("est_shared_pop_r", "est_shared_r", "fix")) stop("Please specify a valid steepness option. These include: est_shared_pop_r, fix, est_shared_r. Conversely, leave at NULL to estimate all steepness values.")

    # Share across populations and regions and estimate
    if(h_spec == "est_shared_pop_r") input_list$map$steepness_h <- factor(rep(1, length(input_list$par$steepness_h)))

    # Share across regions but estimate for each population
    if(h_spec == "est_shared_r") input_list$map$steepness_h <- factor(rep(1:input_list$data$n_pop, times = input_list$data$n_regions))

    # Fix all steepness values
    if(h_spec == "fix") input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))

    collect_message("Steepness is specified as: ", h_spec) # output message

  } else {
    # if a stock-recruit curve is used and estimating all steepness parameters
    if(input_list$data$rec_model %in% c(1, 2)) {
      # estimate steepness for all populations
      if(input_list$data$n_pop > 1) input_list$map$steepness_h <- factor(rep(1:input_list$data$n_pop, times = input_list$data$n_regions)) # estimating all steepness parameters by  population
      if(input_list$data$n_pop == 1) input_list$map$steepness_h <- factor(1:input_list$data$n_regions) # estimating all steepness parameters by region
    }
    collect_message("Steepness is estimated for all relavant dimensions")
  }
  return(input_list)
}

#' Map sex ratio parameters
#'
#' Builds the factor map for \code{sexratio_pars} \code{[n_pop x n_regions x
#' n_sexratio_blocks]}, the proportion of recruits assigned to the first sex, where
#' \code{n_sexratio_blocks} is the largest number of time blocks across population
#' and region in \code{$data$sexratio_blocks}. Called by \code{\link{Setup_Mod_Rec}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#'   Requires \code{$data$n_sexes}, \code{$data$n_pop}, \code{$data$n_regions},
#'   \code{$data$sexratio_blocks} and \code{$data$rec_region_prop_spec}.
#' @param sexratio_spec Estimation structure for \code{sexratio_pars}.
#'   \code{"est_all"} gives one parameter per population, region and block, and is
#'   not permitted when \code{rec_region_prop_spec = 1} and \code{n_pop > 1}.
#'   \code{"est_shared_r"} gives one per population and block, shared across its
#'   regions. \code{"est_shared_pop_r"} gives one per block across every population
#'   and region, and requires identical block structures. \code{"fix"} holds every
#'   parameter at its starting value and is required when \code{n_sexes = 1}.
#'
#' @return \code{input_list} with \code{$map$sexratio_pars} set to a factor vector
#'   of length \code{prod(dim(par$sexratio_pars))}. Active parameters take
#'   sequential integers, fixed or invalid cells are \code{NA}.
#'
#' @keywords internal
do_sexratio_pars_mapping <- function(input_list, sexratio_spec) {

  # Initialize arrays and counters
  map_sexratio <- input_list$par$sexratio_pars
  map_sexratio[] <- NA
  sexratio_counter <- 1

  # Validate inputs here
  if(!sexratio_spec %in% c("est_all", "est_shared_pop_r", "est_shared_r", "fix")) stop("Sex Ratio Specificaiton is not correctly specified. Needs to be fix, est_all, est_shared_pop_r, or est_shared_r")
  if(input_list$data$n_sexes == 1 && sexratio_spec != 'fix') stop('Sex Ratio is being estiamted, but there is only 1 sex!')

  # Validate whether blocking structure is appropriate
  if(sexratio_spec == 'est_shared_pop_r') {
    ref_blocks <- input_list$data$sexratio_blocks[1,1,]
    for(pp in 1:input_list$data$n_pop) {
      for(rr in 1:input_list$data$n_regions) {
        if(!identical(as.vector(input_list$data$sexratio_blocks[pp,rr,]), as.vector(ref_blocks))) {
          stop("est_shared_pop_r requires consistent sex ratio block structure across all populations and regions.")
        }
      } # end rr loop
    } # end pp loop
  } # end if

  # if we want to fix
  if(sexratio_spec == 'fix') map_sexratio[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(r in 1:input_list$data$n_regions) {

      # Get number of sex ratio rate blocks
      sexratio_blocks_tmp <- unique(as.vector(input_list$data$sexratio_blocks[p,r,]))

      for(b in seq_along(sexratio_blocks_tmp)) {

        # Estimate for all regions
        if(sexratio_spec == 'est_all') {
          if(input_list$data$n_pop > 1 && input_list$data$rec_region_prop_spec == 1)
            stop("Can't estimate recruitment sex ratio for all populations and regions if no recruitment dispersal is occuring within a given region! Please specify est_shared_r or est_shared_pop_r instead!")
          map_sexratio[p,r,b] <- sexratio_counter
          sexratio_counter <- sexratio_counter + 1
        }

        # Estimate but share sex ratio across regions
        if(sexratio_spec == 'est_shared_r' && r == 1) {
          for(rr in 1:input_list$data$n_regions) {
            # only assign if this value exists for this region
            if(sexratio_blocks_tmp[b] %in% input_list$data$sexratio_blocks[p,rr,]) {
              map_sexratio[p,rr, b] <- sexratio_counter
            } # end if
          } # end rr loop
          sexratio_counter <- sexratio_counter + 1
        }

        if(sexratio_spec == 'est_shared_pop_r' && p == 1 && r == 1) {
          for(pp in 1:input_list$data$n_pop) {
            for(rr in 1:input_list$data$n_regions) {
              # only assign if this block exists for this pop/region combo
              if(sexratio_blocks_tmp[b] %in% input_list$data$sexratio_blocks[pp,rr,]) {
                map_sexratio[pp,rr,b] <- sexratio_counter
              }
            } # end rr loop
          } # end pp loop
          sexratio_counter <- sexratio_counter + 1
        }

      } # end b loop
    } # end r loop
  } # end p loop

  collect_message("Sex ratio is specified as: ", sexratio_spec)

  # input sex ratio rates into mapping list
  input_list$map$sexratio_pars <- factor(map_sexratio) # sex ratio rates

  return(input_list)
}

#' Map recruitment regional apportionment parameters
#'
#' Builds the factor map for \code{rec_region_prop_pars} \code{[n_pop x
#' (n_regions - 1)]}, the logit-scale share of recruits per region, under a softmax
#' with one reference region omitted. Both the parameter and its map are
#' \code{NULL} when \code{n_regions = 1}. Called by \code{\link{Setup_Mod_Rec}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#'   Requires \code{$data$n_pop}, \code{$data$n_regions} and
#'   \code{$data$natal_region}.
#' @param rec_region_prop_spec Dispersal structure. \code{"no_dispersal"} assigns
#'   recruits entirely to their natal region, overwriting the starting values with
#'   \code{-20} for non-natal regions and \code{+20} for the natal region when
#'   \code{natal_region > 1}, and mapping every element to \code{NA}. \code{NULL}
#'   estimates all independently. Both require \code{n_regions > 1}.
#'
#' @return \code{input_list} with \code{$map$rec_region_prop_pars} set to a factor
#'   vector of length \code{n_pop * (n_regions - 1)}, or \code{NULL} when
#'   \code{n_regions = 1}. Starting values are overwritten under
#'   \code{"no_dispersal"}.
#'
#' @keywords internal
do_rec_region_prop_mapping <- function(input_list, rec_region_prop_spec) {

  # Validate spec options
  valid_specs <- c("no_dispersal")
  if(!is.null(rec_region_prop_spec) && !rec_region_prop_spec %in% valid_specs) {
    stop("Invalid rec_region_prop_spec: '", rec_region_prop_spec, "'. Valid options are: ", paste(valid_specs, collapse = ", "), ", or NULL to estimate all.")
  }

  # with one region there is no apportionment to fix
  if(!is.null(rec_region_prop_spec) && rec_region_prop_spec == "no_dispersal" && input_list$data$n_regions == 1) stop("'no_dispersal' is only valid when n_regions > 1. With a single region there is no recruitment apportionment to fix.")

  # par is [n_pop, n_regions-1]
  if(!is.null(rec_region_prop_spec) && rec_region_prop_spec == 'no_dispersal') {
    par_mat <- matrix(-20, nrow = input_list$data$n_pop, ncol = input_list$data$n_regions - 1)
    natal_region <- input_list$data$natal_region
    for(p in seq_len(input_list$data$n_pop)) {
      if(natal_region[p] > 1) par_mat[p, natal_region[p] - 1] <- 20
    }
    input_list$par$rec_region_prop_pars <- par_mat # fix values
    input_list$map$rec_region_prop_pars <- factor(rep(NA, length(par_mat)))  # fix all
  }

  # estimate all recruitment propostions if n_regions > 1
  if(is.null(rec_region_prop_spec) && input_list$data$n_regions > 1) input_list$map$rec_region_prop_pars <- factor(seq_along(input_list$par$rec_region_prop_pars))

  # single region - not even a parameter
  if(input_list$data$n_regions == 1) {
    input_list$par$rec_region_prop_pars <- NULL
    input_list$map$rec_region_prop_pars <- NULL
  }

  return(input_list)
}

#' Map stray rate parameters
#'
#' Builds the factor map for \code{stray_rate_pars} \code{[n_pop x
#' max_stray_blocks]}, the logit-scale stray rates, where \code{max_stray_blocks} is
#' the largest number of time blocks across populations in
#' \code{$data$stray_rate_blocks}. Every parameter is fixed when \code{n_pop = 1} or
#' \code{use_fixed_stray_rate = 1}, whatever \code{stray_rate_spec} asks for. Called
#' by \code{\link{Setup_Mod_Rec}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#'   Requires \code{$data$n_pop}, \code{$data$stray_rate_blocks} and
#'   \code{$data$use_fixed_stray_rate}.
#' @param stray_rate_spec Estimation structure. \code{"fix"} holds every parameter
#'   at its starting value, \code{"est_all"} estimates per population and block, and
#'   \code{"est_shared_pop"} gives one parameter per block shared across
#'   populations, which requires identical block structures and errors otherwise.
#'
#' @return \code{input_list} with \code{$map$stray_rate_pars} set to a factor vector
#'   of length \code{prod(dim(par$stray_rate_pars))}. Active parameters take
#'   sequential integers, fixed ones are \code{NA}.
#'
#' @keywords internal
do_stray_rate_mapping <- function(input_list, stray_rate_spec) {

  map_stray      <- input_list$par$stray_rate_pars
  map_stray[]    <- NA
  stray_counter  <- 1

  valid_specs <- c("fix", "est_all", "est_shared_pop")
  if (!stray_rate_spec %in% valid_specs)
    stop("Invalid stray_rate_spec. Must be one of: ", paste(valid_specs, collapse = ", "))

  # Not applicable for single population
  if (input_list$data$n_pop == 1) {
    input_list$map$stray_rate_pars <- factor(map_stray)
    collect_message("Stray rates fixed (n_pop == 1, straying not applicable).")
    return(input_list)
  }

  # Externally fixed, pars not used by objective function
  if (input_list$data$use_fixed_stray_rate == 1) {
    input_list$map$stray_rate_pars <- factor(map_stray)
    collect_message("Stray rates are externally fixed (use_fixed_stray_rate == 1).")
    return(input_list)
  }

  if (stray_rate_spec == "fix") {
    input_list$map$stray_rate_pars <- factor(map_stray)
    collect_message("Stray rates fixed at starting values.")
    return(input_list)
  }

  # Validate block consistency for shared estimation
  if (stray_rate_spec == "est_shared_pop") {
    ref_blks <- sort(unique(as.vector(input_list$data$stray_rate_blocks[1, ])))
    for (p in seq_len(input_list$data$n_pop)) {
      if (!identical(sort(unique(as.vector(input_list$data$stray_rate_blocks[p, ]))), ref_blks))
        stop("est_shared_pop requires identical stray rate block structure across all populations.")
    }
  }

  for (p in seq_len(input_list$data$n_pop)) {

    blks <- unique(as.vector(input_list$data$stray_rate_blocks[p, ]))

    for (b in blks) {

      if (stray_rate_spec == "est_all") {
        map_stray[p, b] <- stray_counter
        stray_counter    <- stray_counter + 1
      }

      if (stray_rate_spec == "est_shared_pop" && p == 1) {
        for (pp in seq_len(input_list$data$n_pop)) map_stray[pp, b] <- stray_counter
        stray_counter <- stray_counter + 1
      }
    }
  }

  input_list$map$stray_rate_pars <- factor(map_stray)
  collect_message("Stray rates specified as: ", stray_rate_spec)
  return(input_list)
}

#' Map recruitment seasonal apportionment parameters
#'
#' Builds the factor map for \code{rec_seas_prop_pars} \code{[n_pop x (n_seas - 1)]},
#' the logit-scale share of annual recruitment per season, under a softmax with one
#' reference season omitted. Called by \code{\link{Setup_Mod_Rec}}. When
#' \code{rec_lag = 0} and \code{spawn_seas > 1} the seasons before \code{spawn_seas}
#' are fixed at zero by a restricted softmax in the model, so the trailing unused
#' columns are forced to \code{NA} whatever \code{rec_seas_prop_spec} asks for.
#'
#' @param input_list Named list with \code{$data}, \code{$par} and \code{$map}.
#'   Requires \code{$data$n_pop}, \code{$data$n_seas} and
#'   \code{$data$use_fixed_rec_seas_prop}.
#' @param rec_seas_prop_spec Seasonal apportionment structure.
#'   \code{"est_shared_pop"} estimates one set of \code{n_seas - 1} parameters for
#'   every population and is only valid when \code{n_seas > 1}. \code{"fix"} holds
#'   every parameter at its starting value. \code{NULL} estimates all
#'   \code{n_pop x (n_seas - 1)} independently. The first and last reset
#'   \code{use_fixed_rec_seas_prop} to \code{0} with a warning if it was \code{1}.
#'
#' @return \code{input_list} with \code{$map$rec_seas_prop_pars} set to a factor
#'   vector of length \code{n_pop * (n_seas - 1)}, or \code{NULL} when
#'   \code{n_seas = 1}. Both \code{$par$rec_seas_prop_pars} and its map are
#'   \code{NULL} in that case. \code{$data$use_fixed_rec_seas_prop} may change as a
#'   side effect.
#'
#' @keywords internal
do_rec_seas_prop_mapping <- function(input_list, rec_seas_prop_spec) {

  # Validate spec options
  valid_specs <- c("fix", "est_shared_pop")
  if(!is.null(rec_seas_prop_spec) && !rec_seas_prop_spec %in% valid_specs) {
    stop("Invalid rec_seas_prop_spec: '", rec_seas_prop_spec, "'. Valid options are: ", paste(valid_specs, collapse = ", "), ", or NULL to estimate all.")
  }

  if((is.null(rec_seas_prop_spec) || rec_seas_prop_spec == 'est_shared_pop') && input_list$data$use_fixed_rec_seas_prop == 1) {
    input_list$data$use_fixed_rec_seas_prop <- 0
    warning("Recruitment seasonal apportionment is specified as estimated, but use_fixed_rec_seas_prop == 1 (fixed). Changing to use_fixed_rec_seas_prop == 0.")
  }

  # estimating recruitment seasonal apporitonment is only valid for seasonal models
  if(!is.null(rec_seas_prop_spec) && rec_seas_prop_spec == 'est_shared_pop' && input_list$data$n_seas == 1)
    stop("Estimating recruitment seasonal apportionment is only applicable for seasonal models. ")

  # estimate all recruitment seasonal proportions if n_seas > 1
  if(is.null(rec_seas_prop_spec)) {
    input_list$map$rec_seas_prop_pars <- factor(seq_along(input_list$par$rec_seas_prop_pars))
  } else if(rec_seas_prop_spec == 'est_shared_pop') { # estimate recruitment seasonal proportions but share across populations
    counter <- 1
    tmp_map = input_list$par$rec_seas_prop_pars
    for(seas in 1:(input_list$data$n_seas - 1)) {
      tmp_map[,seas] <- counter
      counter + 1
    }
    input_list$map$rec_seas_prop_pars <- factor(tmp_map)
  } else if(rec_seas_prop_spec == 'fix') input_list$map$rec_seas_prop_pars <- factor(rep(NA, length(input_list$par$rec_seas_prop_pars)))

  # single seas - not even a parameter
  if(input_list$data$n_seas == 1) {
    input_list$par$rec_seas_prop_pars <- NULL
    input_list$map$rec_seas_prop_pars <- NULL
  }

  # under age-0 recruitment with spawning after season 1, seasons before spawn_seas are fixed at
  # zero, so map the structurally unused columns of rec_seas_prop_pars to NA
  if(!is.null(input_list$map$rec_seas_prop_pars) &&
     input_list$data$rec_lag == 0 && input_list$data$spawn_seas > 1) {
    n_allowed <- input_list$data$n_seas - input_list$data$spawn_seas + 1
    n_unused <- (input_list$data$n_seas - 1) - (n_allowed - 1) # trailing unused columns per population
    if(n_unused > 0) {
      tmp_map <- matrix(as.numeric(as.character(input_list$map$rec_seas_prop_pars)),
                         nrow = input_list$data$n_pop)
      tmp_map[, (n_allowed):(input_list$data$n_seas - 1)] <- NA
      input_list$map$rec_seas_prop_pars <- factor(tmp_map)
    }
  }

  return(input_list)
}

#' Set up the recruitment module and associated processes
#'
#' Sets the stock-recruit form and density dependence, steepness and its prior,
#' \eqn{\sigma_R}, annual and initial deviations, regional and seasonal
#' apportionment, spawning movement, stray rates, sex ratio, the equilibrium
#' initialization and the bias ramp. Call after \code{\link{Setup_Mod_Dim}} and
#' \code{\link{Setup_Mod_Biologicals}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}. Dimensions must already be set in \code{$data}.
#' @param rec_model Character, required. \code{"mean_rec"} (steepness fixed and not
#'   estimated), \code{"bh_rec"}, or \code{"ricker_rec"}. Steepness is not
#'   interchangeable between the last two, see \code{\link{Get_Det_Recruitment}}.
#' @param rec_dd Density dependence. \code{"local"} is one stock-recruit relationship
#'   per population, required when \code{n_pop > 1}. \code{"global"} (default) pools
#'   across regions and restricts \code{h_spec}, \code{RecDevs_spec} and
#'   \code{InitDevs_spec} to shared or fixed options when \code{n_regions > 1}.
#' @param SR_ref_yr Integer year index supplying every input to unfished spawning
#'   biomass per recruit, and so to \code{S0} and the curve's scale: weight-at-age,
#'   maturity, natural mortality, movement, stray rate, sex ratio, and what enters
#'   through \code{init_F}. \code{R0} is the exception and is always the year's own
#'   value. Default \code{1}. Ignored under \code{rec_model = "mean_rec"}.
#' @param rec_lag Integer lag in seasons between spawning biomass and recruitment.
#'   \code{1} (default) uses SSB from that many seasons prior. \code{0} is age-0
#'   recruitment on the same year's SSB, so recruits may only enter in
#'   \code{spawn_seas} or later.
#' @param sigmaR_spec Character. Estimation structure for \eqn{\sigma_R}, stored in
#'   \code{ln_sigmaR} \code{[2 x n_pop x n_regions]} with index 1 the initial period
#'   and 2 the annual period. Default \code{"est_all"}, see
#'   \code{\link{do_sigmaR_mapping}}.
#' @param sigmaR_switch Integer year index at which \eqn{\sigma_R} switches from the
#'   early to the late value. \eqn{\leq 1} (default) uses one value throughout.
#' @param sr_penalty Character. \code{"none"} (default), \code{"bh"} or
#'   \code{"ricker"}. Only valid under \code{rec_model = "mean_rec"}. Fits the curve
#'   as a likelihood on \eqn{\log R_y - \log\widehat{R}_y} without letting it
#'   generate recruitment.
#' @param sr_pen_sigma Numeric standard deviation of that residual.
#' @param sr_pen_yrs Years the stock-recruit penalty applies over, or \code{NULL}
#'   (default) for every year with a lagged spawning biomass. Naming a year without
#'   one is an error.
#' @param sr_R0_spec Character. \code{"shared"} (default) takes the curve's scale
#'   from \code{ln_global_R0}, \code{"est"} gives it its own \code{ln_sr_R0}, and
#'   \code{"rinit"} takes it from \code{ln_rinit} and requires \code{use_rinit = 1}.
#' @param Use_rec_level_pen Integer (0/1). Whether the log recruitment series itself
#'   is penalized, separately from the deviation penalty. Default \code{0}.
#' @param rec_level_pen_sigma Numeric standard deviation of that penalty. A sum of
#'   squares with weight \eqn{w} corresponds to \eqn{1/\sqrt{2w}}. Default \code{1}.
#' @param rec_level_pen_center \code{"own_mean"} (default) centers on the mean of the
#'   log recruitment series, \code{"fixed"} centers on zero.
#' @param rec_level_pen_yrs Years the penalty applies over, or \code{NULL} (default)
#'   for every year.
#' @param Use_init_sex_pen Integer (0/1). Whether each later sex's initial age
#'   deviations are tied to the first sex's by a Gaussian on their difference.
#'   Requires \code{n_sexes > 1} and \code{InitDevs_sex_spec = "est_all"}. Enters the
#'   objective unweighted. Default \code{0}.
#' @param init_sex_pen_sigma Numeric standard deviation of that tie. Default \code{1}.
#' @param RecDevs_pen_center,InitDevs_pen_center Where the recruitment and initial age
#'   deviation penalties are centered. \code{"fixed"} (default) centers on zero or the
#'   bias-corrected \eqn{-\sigma_R^2/2}; \code{"own_mean"} centers on the estimated
#'   deviations' own mean, leaving their level free to be set elsewhere. Cannot be
#'   combined with \code{do_rec_bias_ramp = 1}.
#' @param RecDevs_spec Character or \code{NULL}. Sharing structure for
#'   \code{ln_RecDevs} \code{[n_pop x n_regions x n_years]}. Default \code{NULL}
#'   (all independent), see \code{\link{do_RecDevs_mapping}}.
#' @param RecDevs_model Process error on \code{ln_RecDevs}. \code{"iid"} (default) is
#'   independent about the prior mean, \code{"rw"} centers each deviation on the
#'   previous one with a diffuse first year, \code{"ar1"} reverts toward zero at rate
#'   \code{RecDevs_rho} with a stationary first year, and \code{"dsem"} takes the
#'   density from the arrows in \code{\link{Setup_Mod_DSEM}}. Steps span estimated
#'   years, not calendar years. \code{"rw"} and \code{"ar1"} are refused alongside
#'   \code{do_rec_bias_ramp = 1} or \code{RecDevs_pen_center = "own_mean"};
#'   \code{"dsem"} reads \code{sigmaR} off the arrows, so \code{sigmaR_spec} other
#'   than \code{"fix"}, \code{dont_est_recdev_last > 0} and a nonzero ramp are refused.
#' @param RecDevs_rho_spec Sharing structure for \code{RecDevs_rho} \code{[n_pop x
#'   n_regions]}: \code{"est_all"}, \code{"est_shared_pop"}, \code{"est_shared_r"},
#'   \code{"est_shared_pop_r"} or \code{"fix"} (default). Only read under
#'   \code{RecDevs_model = "ar1"}, see \code{\link{do_RecDevs_rho_mapping}}.
#' @param RecDevs_rw_init_sigma Standard deviation given to year one of a random walk,
#'   which sets the level of the series. Default \code{5}; \code{NA} instead starts the
#'   walk at zero under its own sigma. Only read under \code{RecDevs_model = "rw"}.
#' @param ln_global_R0_spec \code{"est"} (default) or \code{"fix"}. \code{"fix"} maps
#'   \code{ln_global_R0} off at its starting value, so the deviations hold log
#'   recruitment outright. Under \code{rec_model = "mean_rec"} a random walk with
#'   \code{dont_pen_recdev_first >= 1} leaves the level unidentified and is refused; a
#'   walk with the first year still penalized is accepted with a warning.
#' @param dont_pen_recdev_first Integer. How many leading years of recruitment
#'   deviations are estimated but left out of the penalty. \code{0} (default) penalizes
#'   every year.
#' @param dont_est_recdev_last Non-negative integer. Terminal years for which
#'   recruitment deviations are not estimated. Forced to \code{0} when
#'   \code{n_proj_yrs_devs > 0}, refused under \code{RecDevs_model = "dsem"}. Default
#'   \code{0}.
#' @param init_age_strc Initialization method, default \code{2}.
#'   \code{0}/\code{"iterative"} iterates to approximate equilibrium,
#'   \code{1}/\code{"scalar_no_move"} is a scalar geometric series without movement,
#'   \code{2}/\code{"matrix"} is the matrix series with movement, and
#'   \code{3}/\code{"scalar_plus_only"} moves only the plus group; all four treat
#'   \code{ln_InitDevs} as multiplicative deviations from the equilibrium.
#'   \code{4}/\code{"free"} projects no equilibrium: ages 2 and older are
#'   \code{exp(ln_InitDevs)} apportioned by sex ratio, so the deviations are on the
#'   scale of numbers and \code{equil_init_age_strc} becomes a prior on log abundance.
#' @param equil_init_age_strc Plus-group treatment during stochastic initialization,
#'   default \code{1}. \code{0}/\code{"equil"} estimates no \code{ln_InitDevs},
#'   \code{1}/\code{"stoch_no_plus"} estimates every age but the plus group,
#'   \code{2}/\code{"stoch_all"} estimates every age, \code{3}/\code{"stoch_shared_ages"}
#'   shares ages through \code{init_age_devs_shared} (non-\code{NULL} required; the plus
#'   group is not fixed automatically), and \code{4}/\code{"stoch_all_no_pen"} estimates
#'   every age and penalizes none, for pairing with \code{init_age_strc = "free"}.
#' @param InitDevs_spec Character or \code{NULL}. Sharing structure for
#'   \code{ln_InitDevs} \code{[n_pop x n_regions x (n_ages - 1) x n_sexes]}. Default
#'   \code{NULL} (all independent), see \code{\link{do_InitDevs_mapping}}.
#' @param InitDevs_sex_spec \code{"est_shared_s"} (default) estimates one initial age
#'   deviation curve read by every sex; \code{"est_all"} gives each sex its own,
#'   pooling the level across sexes under an \code{"own_mean"}
#'   \code{InitDevs_pen_center}. Requires \code{n_sexes > 1}.
#' @param init_F_prop Numeric array \code{[n_regions x n_seas x n_fish_fleets]}.
#'   Legacy interface. A non-zero value without \code{init_F_par} is converted to
#'   \code{ln_init_F = log(init_F_prop)} with \code{init_F_form = "prop"}. Prefer
#'   \code{init_F_par}. Default zero.
#' @param init_F_form What \code{init_F_par} means. \code{"prop"} (default) gives
#'   \code{init_F = exp(ln_init_F) * exp(ln_F_mean)}, a proportion of the estimated
#'   mean F. \code{"abs"} gives \code{init_F = exp(ln_init_F)}, independent of
#'   \code{ln_F_mean}. Use \code{"abs"} when bridging an assessment with a separate
#'   historical F; under \code{"prop"} the two collapse into one parameter and catch
#'   constrains only their product.
#' @param init_F_spec \code{"fix"} (default) or \code{"est"}, whether
#'   \code{init_F_par} is estimated. Sets only the mapping, so it combines freely with
#'   \code{init_F_form}. Refused under \code{init_age_strc = "free"}. The value comes
#'   from \code{init_F_par} \code{[n_regions x n_seas x n_fish_fleets]}, passed through
#'   \code{...}, on the logit scale under \code{"prop"} and the log scale under
#'   \code{"abs"}.
#' @param rec_region_prop_spec Character or \code{NULL}. Regional recruitment dispersal
#'   structure, default \code{NULL} (all estimated freely). Stored as
#'   \code{$data$rec_region_prop_spec}, \code{0} = full dispersal, \code{1} = none. See
#'   \code{\link{do_rec_region_prop_mapping}}.
#' @param use_rec_region_prop_prior Integer (0/1). Dirichlet priors on regional
#'   recruitment proportions. Not valid when \code{n_regions = 1}. Default \code{0}.
#' @param rec_region_prop_prior Data frame of Dirichlet concentrations with columns
#'   \code{pop} and \code{alpha}, a list-column of length-\code{n_regions} vectors.
#'   Default \code{NULL}.
#' @param rec_seas_prop_spec Character or \code{NULL}. Seasonal apportionment
#'   structure, default \code{"fix"}. See \code{\link{do_rec_seas_prop_mapping}}.
#' @param use_fixed_rec_seas_prop Integer (0/1). Whether \code{fixed_rec_seas_prop} is
#'   used. Reset to \code{0} with a warning if \code{rec_seas_prop_spec} estimates.
#'   Default \code{1}.
#' @param fixed_rec_seas_prop Array \code{[n_pop x n_seas]} of fixed seasonal
#'   proportions, default all recruitment in season 1. Must be zero before
#'   \code{spawn_seas} when \code{rec_lag = 0} and \code{spawn_seas > 1}.
#' @param use_rec_seas_prop_prior Integer (0/1). Dirichlet priors on seasonal
#'   proportions. Not valid when \code{n_seas = 1}. Evaluated only over
#'   \code{spawn_seas:n_seas} when \code{rec_lag = 0} and \code{spawn_seas > 1}.
#'   Default \code{0}.
#' @param rec_seas_prop_prior Data frame of Dirichlet concentrations with columns
#'   \code{pop} and \code{alpha}. Default \code{NULL}.
#' @param h_spec Character or \code{NULL}. Sharing structure for \code{steepness_h}
#'   \code{[n_pop x n_regions]}, on a logit scale bounded to \eqn{(0.2, 1)}. Default
#'   \code{NULL}, ignored under \code{rec_model = "mean_rec"}. See
#'   \code{\link{do_h_mapping}}.
#' @param Use_h_prior Integer (0/1). Normal priors on steepness. Default \code{0}.
#' @param h_prior Data frame with columns \code{pop}, \code{region}, \code{mu} and
#'   \code{sd}. Default \code{NULL}.
#' @param spawn_seas Integer season index in which spawning occurs. Default \code{1}.
#' @param t_spawn Numeric fraction of the spawning season elapsed before spawning.
#'   \code{0} (default) spawns before any mortality, \code{1} after all of it.
#' @param sgl_seas_spawning_movement Array \code{[n_pop x n_regions x n_regions x
#'   n_years x n_ages x n_sexes]}, each \code{[p, , r, y, a, s]} slice row-stochastic.
#'   \code{NA} (default) assumes complete natal homing and builds the array internally.
#' @param use_fixed_stray_rate Integer (0/1). Whether stray rates come from
#'   \code{fixed_stray_rate} rather than being estimated. Default \code{1}.
#' @param fixed_stray_rate Array \code{[n_pop x n_years]} of stray rates in
#'   \eqn{[0, 1]}. Default \code{0}.
#' @param stray_rate_spec Estimation structure for \code{stray_rate_pars} \code{[n_pop
#'   x max_stray_blocks]} on the logit scale. \code{"fix"} (default) holds every value,
#'   \code{"est_all"} estimates per population and block, and \code{"est_shared_pop"}
#'   gives one parameter per block shared across populations, which requires identical
#'   block structures. Ignored when \code{use_fixed_stray_rate = 1} or \code{n_pop = 1}.
#' @param stray_rate_blocks Character vector of length \code{n_pop}, either
#'   \code{"none_Pop_x"} for one block or \code{"Block_k_Year_a-b_Pop_x"}, with
#'   \code{"terminal"} allowed as the end year. Default one block each. Stray rate is
#'   generally unidentifiable from fisheries data alone, so use
#'   \code{use_stray_rate_prior} whenever \code{stray_rate_spec != "fix"}.
#' @param use_stray_rate_prior Integer (0/1). Beta priors on estimated stray rates.
#'   Only relevant when \code{use_fixed_stray_rate = 0} and \code{n_pop > 1}; an error
#'   is raised alongside \code{use_fixed_stray_rate = 1}. Default \code{0}.
#' @param stray_rate_prior Data frame with columns \code{pop}, \code{block},
#'   \code{mu} in \eqn{(0,1)} and \code{sd}, one row per population and block.
#'   Default \code{NULL}.
#' @param sexratio_spec Estimation structure for \code{sexratio_pars} \code{[n_pop x
#'   n_regions x n_blocks]}. Default \code{"fix"}, which is required when
#'   \code{n_sexes = 1}. See \code{\link{do_sexratio_pars_mapping}}.
#' @param sexratio_blocks Character vector, one entry per population and region, either
#'   \code{"none_Pop_x_Region_x"} or \code{"Block_k_Year_a-b_Pop_x_Region_x"}, with
#'   \code{"terminal"} allowed as the end year. Default one block each.
#' @param do_rec_bias_ramp Integer (0/1). Whether a bias ramp is applied to
#'   \code{ln_RecDevs}. Under \code{0} every penalty is centered on the full
#'   \eqn{-\sigma_R^2/2}, under \code{1} the center follows the ramp. Default \code{0}.
#' @param bias_year Numeric calendar year at which the ramp reaches its maximum
#'   correction. Default \code{NA}.
#' @param max_bias_ramp_fct Numeric in \eqn{[0, 1]}, the maximum correction applied at
#'   \code{bias_year}. Default \code{1}.
#' @param ... Optional named starting values: \code{ln_global_R0} \code{[n_pop]},
#'   \code{ln_rinit} \code{[n_pop]}, \code{rec_region_prop_pars} \code{[n_pop x
#'   (n_regions - 1)]}, \code{rec_seas_prop_pars} \code{[n_pop x (n_seas - 1)]},
#'   \code{steepness_h} \code{[n_pop x n_regions]} on the bounded logit scale,
#'   \code{ln_InitDevs} \code{[n_pop x n_regions x (n_ages - 1) x n_sexes]} (a 3-D
#'   array is expanded across sexes), \code{ln_RecDevs} \code{[n_pop x n_regions x
#'   n_years]}, \code{ln_sigmaR} \code{[2 x n_pop x n_regions]}, \code{sexratio_pars}
#'   \code{[n_pop x n_regions x n_blocks]}.
#' @param use_rinit Integer (0/1). Whether \code{ln_rinit} initializes the population
#'   separately from \code{ln_global_R0}. Under \code{0} (default) \code{ln_rinit} is
#'   fixed and \code{ln_global_R0} does both jobs.
#' @param init_age_devs_shared Integer vector of length \code{n_ages - 1} giving the
#'   factor level of each age position for \code{ln_InitDevs}; positions sharing a
#'   value share one parameter. Read under \code{equil_init_age_strc = 3}, and also
#'   respected by \code{InitDevs_spec = "est_shared_r"} (applied per population, with a
#'   population offset) and \code{"est_shared_pop_r"} (applied globally, no offset).
#'   \code{c(1:42, rep(42, 9))} gives 42 free parameters for a 52-age model with 43
#'   data ages. Default \code{NULL}.
#' @param Use_rinit_pen Integer (0/1). Whether \eqn{\log(R_{init} / R_0)} is penalized
#'   under \code{use_rinit = 1}. An equilibrium recruitment stands for an average of
#'   several years, so \eqn{\sigma_R / (1 / M - 0.5)} is a reasonable sd. Default 0.
#' @param rinit_pen_sd Standard deviation of that penalty, log scale. Default 1.
#' @param R0_blocks Character vector of time blocks for \code{R0}, one per population,
#'   in the selectivity block vocabulary: \code{"none_Pop_<p>"} or
#'   \code{"Block_<b>_Year_<a>-<e>_Pop_<p>"} on 1-based year indices, with
#'   \code{"terminal"} allowed. Under \code{"mean_rec"} a block is a productivity
#'   regime; under a stock-recruit form it makes the curve time-varying, so \code{S0},
#'   depletion and any reference point step at the boundary. Default \code{NULL}.
#' @param R0_ref_block Integer, the block whose \code{R0} is used wherever a single
#'   value is needed: the initial age structure, the regional apportionment, the
#'   \code{R0} prior, the \code{ln_rinit} penalty, and the stock-recruit scale under
#'   \code{sr_R0_spec = "shared"}. Default 1.
#' @param use_r0_prior Integer (0/1). Lognormal prior on \code{R0}. Default 0.
#' @param r0_prior Data frame with columns \code{pop}, \code{mu} on the natural scale
#'   and \code{sd} on the log scale. Required when \code{use_r0_prior = 1}.
#'
#' @return \code{input_list} with recruitment fields set in \code{$data} and
#'   \code{$par}, and maps built in \code{$map} for \code{rec_region_prop_pars},
#'   \code{rec_seas_prop_pars}, \code{ln_sigmaR}, \code{ln_InitDevs},
#'   \code{ln_RecDevs}, \code{RecDevs_rho}, \code{steepness_h}, \code{sexratio_pars}
#'   and \code{stray_rate_pars}. Character codes for \code{init_age_strc} and
#'   \code{equil_init_age_strc} are converted to integers before storage.
#'
#' @family Model Setup
#' @export Setup_Mod_Rec
Setup_Mod_Rec <- function(input_list,
                          rec_model,
                          rec_dd = "global",
                          rec_lag = 1,
                          SR_ref_yr = 1,
                          Use_h_prior = 0,
                          h_prior = NULL,
                          rec_region_prop_spec = NULL,
                          use_rec_region_prop_prior = 0,
                          rec_region_prop_prior = NULL,
                          rec_seas_prop_spec = "fix",
                          use_rec_seas_prop_prior = 0,
                          rec_seas_prop_prior = NULL,
                          use_fixed_rec_seas_prop = 1,
                          fixed_rec_seas_prop = {
                            rec_seas_prop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_seas))
                            rec_seas_prop[, 1] <- 1
                            rec_seas_prop
                          },
                          do_rec_bias_ramp = 0,
                          bias_year = NA,
                          max_bias_ramp_fct = 1,
                          sigmaR_switch = 1,
                          dont_est_recdev_last = 0,
                          dont_pen_recdev_first = 0,
                          init_age_strc = 2,
                          equil_init_age_strc = 1,
                          init_F_prop = array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets)),
                          init_F_form = "prop",
                          init_F_spec = "fix",
                          sigmaR_spec = "est_all",
                          InitDevs_spec = NULL,
                          InitDevs_sex_spec = "est_shared_s",
                          RecDevs_spec = NULL,
                          RecDevs_model = "iid",
                          RecDevs_rho_spec = "fix",
                          RecDevs_rw_init_sigma = 5,
                          RecDevs_pen_center = "fixed",
                          Use_rec_level_pen = 0,
                          rec_level_pen_sigma = 1,
                          rec_level_pen_center = "own_mean",
                          rec_level_pen_yrs = NULL,
                          Use_init_sex_pen = 0,
                          init_sex_pen_sigma = 1,
                          sr_penalty = "none",
                          sr_pen_sigma = 1,
                          sr_pen_yrs = NULL,
                          sr_R0_spec = "shared",
                          InitDevs_pen_center = "fixed",
                          h_spec = NULL,
                          sgl_seas_spawning_movement = NA,
                          t_spawn = 0,
                          stray_rate_spec = "fix",
                          stray_rate_blocks = paste0("none_Pop_", seq_len(input_list$data$n_pop)),
                          use_fixed_stray_rate = if(stray_rate_spec != 'fix') 0 else 1,
                          fixed_stray_rate = array(0, dim = c(input_list$data$n_pop, length(input_list$data$years))),
                          use_stray_rate_prior = 0,
                          stray_rate_prior = NULL,
                          spawn_seas = 1,
                          sexratio_spec = 'fix',
                          sexratio_blocks = {
                            grid <- expand.grid(region = 1:input_list$data$n_regions, pop = 1:input_list$data$n_pop)
                            blks <- paste0("none_Pop_", grid$pop, "_Region_", grid$region)
                            blks
                          },
                          use_rinit = 0,
                          init_age_devs_shared = NULL,
                          R0_blocks = NULL,
                          R0_ref_block = 1,
                          use_r0_prior = 0,
                          r0_prior = NULL,
                          Use_rinit_pen = 0,
                          rinit_pen_sd = 1,
                          ...,
                          ln_global_R0_spec = "est"
                          ) {

  # define stuff so not overidden
  sigmaR_given <- !missing(sigmaR_spec)
  dont_est_recdev_last_given <- dont_est_recdev_last

  messages_list <<- character(0) # nolint: object_usage_linter.
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Rec <- mget(names(formals()))[-1]

  # Convert character inputs to numeric codes for init_age_strc and equil_init_age_strc
  init_age_strc <- convert_to_numeric(init_age_strc, list(iterative = 0, scalar_no_move = 1, matrix = 2, scalar_plus_only = 3, free = 4))
  equil_init_age_strc <- convert_to_numeric(equil_init_age_strc, list(equil = 0, stoch_no_plus = 1, stoch_all = 2, stoch_shared_ages = 3, stoch_all_no_pen = 4))

  # Setting up the initial fishing mortality ----------------------------------------
  if(!(length(init_F_form) == 1 && init_F_form %in% c("prop", "abs"))) stop("init_F_form must be 'prop' (a proportion of the mean F) or 'abs' (an absolute F), but was: ", paste(init_F_form, collapse = ", "))
  if(!(length(init_F_spec) == 1 && init_F_spec %in% c("fix", "est"))) stop("init_F_spec must be 'fix' or 'est', but was: ", paste(init_F_spec, collapse = ", "))

  init_F_form_num <- convert_to_numeric(init_F_form, list(prop = 0, abs = 1))
  init_F_est      <- convert_to_numeric(init_F_spec, list(fix = 0, est = 1))

  # a free initial age structure is exp(ln_InitDevs), so no equilibrium is projected
  # under an initial F and init_F_par never reaches the objective
  if(init_F_est == 1 && init_age_strc == 4)
    stop("init_F_spec = 'est' with init_age_strc = 'free' leaves init_F_par unidentified: the free ",
         "initialization takes the numbers at age 2 and older as exp(ln_InitDevs), so no ",
         "equilibrium is projected under an initial F and init_F_par never reaches the objective. ",
         "Use init_F_spec = 'fix' here, or an equilibrium init_age_strc if the initial F is meant to ",
         "set the age structure.")

  init_F_dim <- c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets)
  init_F_off <- if(init_F_form_num == 0) stats::qlogis(1e-10) else log(1e-100)

  if("init_F_par" %in% names(starting_values)) {
    input_list$par$init_F_par <- array(starting_values$init_F_par, dim = init_F_dim)
  } else if(!all(init_F_prop == 0)) {
    # For backwards compatibility, because init_F_prop was a data array of proportions of mean F.
    if(!all(dim(init_F_prop) == init_F_dim)) stop("init_F_prop must have dimensions [n_regions, n_seas, n_fish_fleets] = [", paste(init_F_dim, collapse = ", "), "]")
    if(any(init_F_prop < 0)) stop("init_F_prop cannot be negative")
    if(init_F_form_num != 0) stop("init_F_prop is a proportion of the mean F, so it requires init_F_form = 'prop'. Supply init_F_par instead to specify an absolute initialization F.")
    if(any(init_F_prop >= 1)) stop("init_F_prop must be < 1: it is now transformed with plogis, which bounds the proportion to (0, 1). Use init_F_form = 'abs' if the initialization F should exceed the mean F.")
    input_list$par$init_F_par <- array(stats::qlogis(pmax(init_F_prop, 1e-10)), dim = init_F_dim)
  } else {
    input_list$par$init_F_par <- array(init_F_off, dim = init_F_dim)
    if(init_F_est == 1) stop("init_F_spec = 'est' but no starting value was given; supply init_F_par (estimating from an effectively-zero start is degenerate)")
  }
  input_list$map$init_F_par <- factor(if(init_F_est == 1) seq_len(prod(init_F_dim)) else rep(NA_integer_, prod(init_F_dim)))

  collect_message("Initialization F is ", ifelse(init_F_est == 1, "estimated", "fixed"), " as ",
                  ifelse(init_F_form_num == 0, "a proportion of the mean F (moves with ln_F_mean)", "an absolute F (independent of ln_F_mean)"), ".")

  # Recruitment Model Type and Options --------------------------------------

  # Recruitment model
  rec_model_map <- list(mean_rec = 0, bh_rec = 1, ricker_rec = 2)
  if (!rec_model %in% names(rec_model_map)) stop("Invalid recruitment model. Use 'mean_rec', 'bh_rec', or 'ricker_rec'")
  rec_model_val <- rec_model_map[[rec_model]]
  collect_message("Recruitment is specified as: ", rec_model)

  # Recruitment density dependence
  if (!is.null(rec_dd)) {
    rec_dd_map <- list(local = 0, global = 1)
    if (!rec_dd %in% names(rec_dd_map)) stop("Invalid rec_dd. Use 'local' or 'global'")
    rec_dd_val <- rec_dd_map[[rec_dd]]
    collect_message("Recruitment Density Dependence is specified as: ", rec_dd)
  } else {
    rec_dd_val <- ifelse(rec_model == "mean_rec", 999, 1)
  }

  # Recruitment lag
  if(rec_model != "mean_rec") collect_message("Recruitment and SSB lag is specified as: ", rec_lag)
  if(rec_lag < 0) stop("rec_lag cannot be negative!")

  # Reference year for unfished spawning biomass per recruit
  if(length(SR_ref_yr) != 1 || SR_ref_yr < 1 || SR_ref_yr > length(input_list$data$years))
    stop("SR_ref_yr must be a single year index between 1 and ", length(input_list$data$years), ".")
  if(rec_model != "mean_rec")
    collect_message("Unfished spawning biomass per recruit uses biologicals from year ",
                    input_list$data$years[SR_ref_yr], ".")

  # age-0 recruits cannot appear before spawn_seas in the same year. fixed seasonal proportions
  # enforce that here; estimated ones get it from the restricted softmax instead
  if(rec_lag == 0 && spawn_seas > 1 && use_fixed_rec_seas_prop == 1 &&
     any(fixed_rec_seas_prop[, seq_len(spawn_seas - 1), drop = FALSE] != 0)) {
    stop("rec_lag = 0 requires fixed_rec_seas_prop to be zero in every season before spawn_seas (age-0 recruits can't predate the spawning event that produced them).")
  }

  # Recruitment regional proportion prior
  if(!use_rec_region_prop_prior %in% c(0,1)) stop("use_rec_region_prop_prior must be 0 or 1")
  if(use_rec_region_prop_prior == 1 && input_list$data$n_regions == 1) stop("Priors should not be applied to recruitment regional proportions when n_regions = 1.")
  if(use_rec_region_prop_prior == 1) {
    required_cols <- c("pop", "alpha")
    missing_cols <- setdiff(required_cols, names(rec_region_prop_prior))
    if (length(missing_cols) > 0) stop("rec_region_prop_prior is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  collect_message("Recruitment regional proportion priors are: ", ifelse(use_rec_region_prop_prior == 1, "Used", "Not Used"))

  # recruitment seasonal priors
  if(!use_rec_seas_prop_prior %in% c(0,1)) stop("use_rec_seas_prop_prior must be 0 or 1")
  if(use_rec_seas_prop_prior == 1 && input_list$data$n_seas == 1) stop("Priors should not be applied to recruitment seasonal proportions when n_seass = 1.")
  if(use_rec_seas_prop_prior == 1) {
    required_cols <- c("pop", "alpha")
    missing_cols <- setdiff(required_cols, names(rec_seas_prop_prior))
    if (length(missing_cols) > 0) stop("rec_seas_prop_prior is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  collect_message("Recruitment seasonal proportion priors are: ", ifelse(use_rec_seas_prop_prior == 1, "Used", "Not Used"))
  collect_message("Recruitment seasonal proportions is: ", ifelse(is.null(rec_seas_prop_spec), "estimated for all dimensions", rec_seas_prop_spec))

  # Checking that rec_dd is local when n_pop > 1
  if(input_list$data$n_pop > 1 && rec_dd != 'local') stop("When n_pop > 1, rec_dd must be local!")


  # Spawning Movement -------------------------------------------------------
  if(is.na(sum(sgl_seas_spawning_movement))) {
    natal_region <- input_list$data$natal_region
    arr <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))
    for(p in seq_len(input_list$data$n_pop)) arr[p, , natal_region[p], , , ] <- 1
    tmp_sgl_seas_spawning_movement <- arr
    if(input_list$data$n_pop > 1 && input_list$data$n_seas == 1 && rec_model_val %in% c(1, 2)) collect_message("Using 100% natal homing rate.")
  } else {
    check_data_dimensions(
      sgl_seas_spawning_movement,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_ages = length(input_list$data$ages),
      n_sexes = input_list$data$n_sexes,
      what = 'sgl_seas_spawning_movement'
    )
    tmp_sgl_seas_spawning_movement <- sgl_seas_spawning_movement
    if(input_list$data$n_pop > 1 && input_list$data$n_seas == 1 && rec_model_val %in% c(1, 2)) collect_message("Using user input natal homing rate.")
  }


  # Straying Rates ----------------------------------------------------------
  stray_rate_blocks_mat <- array(NA, dim = c(input_list$data$n_pop,
                                             length(input_list$data$years)))

  for (i in seq_along(stray_rate_blocks)) {

    # parse
    tmp     <- stray_rate_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    if (!tmp_vec[1] %in% c("none", "Block"))
      stop("stray_rate_blocks not correctly specified. ",
           "Use 'none_Pop_x' or 'Block_k_Year_a-b_Pop_x'.")

    # if none
    if (tmp_vec[1] == "none") {
      pop <- as.numeric(tmp_vec[3])
      stray_rate_blocks_mat[pop, ] <- 1
    }

    # if blocks
    if (tmp_vec[1] == "Block") {
      block_val <- as.numeric(tmp_vec[2])
      pop        <- as.numeric(tmp_vec[6])
      if (!str_detect(tmp, "terminal")) {
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        yrs        <- year_range[1]:year_range[2]
      } else {
        yrs <- as.numeric(unlist(strsplit(tmp_vec[4], "-"))[1]):length(input_list$data$years)
      }
      stray_rate_blocks_mat[pop, yrs] <- block_val
    }
  }

  for (p in seq_len(input_list$data$n_pop))
    collect_message("Stray rates for population ", p, " specified with ", length(unique(stray_rate_blocks_mat[p, ])), " block(s).")

  # validate priors
  if (!use_stray_rate_prior %in% c(0, 1)) stop("use_stray_rate_prior must be 0 or 1")
  if (use_stray_rate_prior == 1 && input_list$data$n_pop == 1)
    stop("Stray rate priors are not applicable when n_pop == 1.")
  if (use_stray_rate_prior == 1 && use_fixed_stray_rate == 1)
    stop("use_stray_rate_prior == 1 but use_fixed_stray_rate == 1 - stray_rate_pars are not estimated so a prior has no effect.")
  if (use_stray_rate_prior == 1) {
    required_cols <- c("pop", "block", "mu", "sd")
    missing_cols  <- setdiff(required_cols, names(stray_rate_prior))
    if (length(missing_cols) > 0)
      stop("stray_rate_prior is missing columns: ", paste(missing_cols, collapse = ", "))
    if (any(stray_rate_prior$mu <= 0 | stray_rate_prior$mu >= 1))
      stop("stray_rate_prior$mu must be in (0, 1).")
  }
  collect_message("Stray rate prior is: ", ifelse(use_stray_rate_prior == 1, "Used", "Not Used"))

  # Steepness Settings ------------------------------------------------------
  if (rec_model %in% c("bh_rec", "ricker_rec")) {
    if (!Use_h_prior %in% c(0, 1)) stop("Use_h_prior must be 0 or 1")
    if (Use_h_prior == 1) {
      required_cols <- c("pop", "region", "mu", "sd")
      missing_cols <- setdiff(required_cols, names(h_prior))
      if (length(missing_cols) > 0) stop("h_prior is missing columns: ", paste(missing_cols, collapse = ", "))
    }
    collect_message("Steepness priors are: ", ifelse(Use_h_prior == 1, "Used", "Not Used"))
  }


  # Sex Ratio Options ---------------------------------------------
  sexratio_blocks_mat <- array(NA, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years)))

  for(i in seq_along(sexratio_blocks)) {

    # Extract out components from list
    tmp <- sexratio_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    if(!tmp_vec[1] %in% c("none", "Block")) stop("Sex Ratio Blocks not correctly specified. This should be either none_Pop_x_Region_x or Block_x_Year_x-y_Pop_x_Region_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      pop <- as.numeric(tmp_vec[3]) # get pop index
      region <- as.numeric(tmp_vec[5]) # get region index
      sexratio_blocks_mat[pop, region,] <- 1 # input sex ratio time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      pop <- as.numeric(tmp_vec[6]) # get pop value
      region <- as.numeric(tmp_vec[8]) # get region value

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      sexratio_blocks_mat[pop,region,years] <- block_val # input sex ratio time block
    }

  } # end i loop

  for(p in 1:input_list$data$n_pop) for(r in 1:input_list$data$n_regions)
    collect_message("Sex Ratios specified with ", length(unique(sexratio_blocks_mat[p,r,])), " block for population ", p, " and region ", r)

  # Input Validation --------------------------------------------------------

  # Helper function
  check_in <- function(x, valid, name) {
    if (!x %in% valid) stop(name, " must be one of: ", paste(valid, collapse = ", "))
  }

  # Validation
  check_in(do_rec_bias_ramp, 0:1, "do_rec_bias_ramp")
  check_in(init_age_strc, 0:4, "init_age_strc")
  if(!is.numeric(sigmaR_switch)) stop("sigmaR_switch must be numeric")
  if(max_bias_ramp_fct > 1 || max_bias_ramp_fct < 0) stop("max_bias_ramp_fct must be between 0 and 1")

  # print messages
  collect_message("Recruitment Bias Ramp is: ", ifelse(do_rec_bias_ramp == 0, "Off", 'On'))
  init_age_methods <- c("Iterated", "No Movement and Scalar Geometric Series", "Movement and Matrix Geometric Series", "Movement but Scalar Geometric Series for plus group", "Free (deviations are the numbers at age, no equilibrium)")
  collect_message("Initial Age Structure is: ", init_age_methods[init_age_strc + 1])
  if(sigmaR_switch > 1) collect_message("Sigma R switches from an early period value to a late period value at year: ", sigmaR_switch)
  collect_message("Recruitment deviations for ", ifelse(dont_est_recdev_last == 0, "every year are estimated", paste("terminal year not estimated -", dont_est_recdev_last)))
  if(dont_est_recdev_last != 0 && input_list$data$n_proj_yrs_devs != 0) {
    collect_message(
      "Recruitment deviations were specified to not be estimated for the last ",
      dont_est_recdev_last,
      " years, but n_proj_yrs_devs != 0. Because projected deviations are still computed (penalized toward the mean), those `unestimated` years are stil effectively estimated. Setting dont_est_recdev_last to 0."
    )
    dont_est_recdev_last <- 0 # overwrite at 0
  }

  # Populate Data List ------------------------------------------------------

  # # input variables into data list
  collect_message("Spawning season occurs in season ", spawn_seas)
  input_list$data$spawn_seas <- spawn_seas
  input_list$data$sgl_seas_spawning_movement <- tmp_sgl_seas_spawning_movement
  input_list$data$rec_model <- rec_model_val
  input_list$data$rec_dd <- rec_dd_val
  input_list$data$rec_lag <- rec_lag
  input_list$data$SR_ref_yr <- SR_ref_yr
  input_list$data$Use_h_prior <- Use_h_prior
  input_list$data$h_prior <- h_prior
  input_list$data$do_rec_bias_ramp <- do_rec_bias_ramp
  input_list$data$bias_year <- bias_year
  input_list$data$sigmaR_switch <- sigmaR_switch
  input_list$data$init_age_strc <- init_age_strc
  if(!RecDevs_pen_center %in% c("fixed", "own_mean")) stop("RecDevs_pen_center must be fixed or own_mean")
  if(!InitDevs_pen_center %in% c("fixed", "own_mean")) stop("InitDevs_pen_center must be fixed or own_mean")
  if(RecDevs_pen_center == "own_mean" && do_rec_bias_ramp == 1) stop("RecDevs_pen_center = own_mean estimates the deviations' mean from the deviations themselves, which leaves the bias ramp's -sigma^2/2 offset meaningless. Use one or the other.")
  input_list$data$RecDevs_pen_center <- convert_to_numeric(RecDevs_pen_center, list(fixed = 0, own_mean = 1))

  # RecDevs_model checking
  if(!RecDevs_model %in% c("iid", "rw", "ar1", "dsem")) stop("RecDevs_model incorrectly specified. Must be one of 'iid', 'rw', 'ar1', or 'dsem'")
  else collect_message("RecDevs_model is specified as: ", RecDevs_model)

  # "dsem" hands the deviations' density to Setup_Mod_DSEM, and their sd with it - so uses dsem sigma and devs rather than those here
  if(RecDevs_model == "dsem") {
    if(sigmaR_given && sigmaR_spec != "fix") stop("RecDevs_model = 'dsem' reads sigmaR off the arrows' recruitment sd line, so ln_sigmaR is not read and cannot be estimated. Leave sigmaR_spec out (it is set to 'fix') or set it to 'fix'.")
    if(dont_est_recdev_last_given > 0) stop("RecDevs_model = 'dsem' describes a recruitment deviation in every year, so dont_est_recdev_last must be 0. The arrows then describe the terminal years too; use RecDevs_model = 'iid' to leave them out.")
    ramp_here <- get_rec_bias_ramp(do_rec_bias_ramp, bias_year, length(input_list$data$years), max_bias_ramp_fct)
    if(do_rec_bias_ramp == 1 && any(ramp_here != 0)) stop("RecDevs_model = 'dsem' makes the recruitment deviations random effects under the arrows, and a random effect takes the full lognormal correction or none: the bias ramp is a device for penalized deviations and has nothing to act on. Set do_rec_bias_ramp = 0 for the full correction, or move bias_year past the last year for none.")
    sigmaR_spec <- "fix"
    input_list$data$dsem_declared <- union(input_list$data$dsem_declared, "rec")
    collect_message("RecDevs_model = 'dsem': the recruitment deviations' density and sd come from Setup_Mod_DSEM. sigmaR is read off the arrows' recruitment sd line (the initial age deviations read it too) and ln_sigmaR is not read. The linked deviations take ", if(any(ramp_here != 0)) "the full lognormal correction, minus half their variance under the arrows, so R0 scales mean recruitment as before." else "no lognormal correction, as the penalty takes none here.")
  }

  if(RecDevs_model %in% c("rw", "ar1") && do_rec_bias_ramp == 1)
    stop("RecDevs_model = '", RecDevs_model, "' centers each deviation on the previous one, so the bias ramp's -sigma^2/2 offset about zero does not apply. Set do_rec_bias_ramp = 0, or use RecDevs_model = 'iid'.")

  if(RecDevs_model %in% c("rw", "ar1") && RecDevs_pen_center == "own_mean")
    stop("RecDevs_model = '", RecDevs_model, "' centers each deviation on the previous one, so there is no single mean for RecDevs_pen_center = 'own_mean' to estimate. Use RecDevs_pen_center = 'fixed'.")

  if(RecDevs_model %in% c("rw", "ar1") && sigmaR_spec == "fix")
    warning("RecDevs_model = '", RecDevs_model, "' but sigmaR_spec = 'fix'; the process error standard deviation (ln_sigmaR) driving the ", RecDevs_model, " process is not being estimated. This may be intentional (e.g. fixing sigma at a known value), but if not, consider estimating ln_sigmaR via sigmaR_spec.")

  if(RecDevs_model == "ar1" && RecDevs_rho_spec == "fix")
    warning("RecDevs_model = 'ar1' but RecDevs_rho_spec = 'fix'; the AR1 correlation parameter (RecDevs_rho) is not being estimated. This may be intentional (e.g. fixing rho at a known value), but if not, consider estimating RecDevs_rho via RecDevs_rho_spec.")

  if(!(length(RecDevs_rw_init_sigma) == 1 && (is.na(RecDevs_rw_init_sigma) || RecDevs_rw_init_sigma > 0)))
    stop("RecDevs_rw_init_sigma must be a single positive number, or NA to start the walk at zero under its own sigma")

  if(!ln_global_R0_spec %in% c("est", "fix")) stop("ln_global_R0_spec must be est or fix")
  collect_message("ln_global_R0 is specified as: ", ln_global_R0_spec)

  # under mean recruitment log R is ln_global_R0 plus a deviation. a walk penalizes only the change
  # between deviations, so it never reads their level, and dropping the first year's term takes away
  # the one thing that did. adding a constant to R0 and taking it off every deviation is then exactly
  # flat: the fit converges, the hessian is singular and every standard error comes back NA
  if(rec_model == "mean_rec" && RecDevs_model == "rw" && ln_global_R0_spec == "est" &&
     dont_pen_recdev_first >= 1)
    stop("rec_model = 'mean_rec' with RecDevs_model = 'rw' and dont_pen_recdev_first = ", dont_pen_recdev_first,
         " leaves ln_global_R0 and the recruitment deviations mutually unidentified: the walk reads only ",
         "the change between deviations and the first year no longer reads their level, so the two trade ",
         "along an exactly flat direction. It needs ln_global_R0_spec = 'fix', which is how SAM writes it, ",
         "with the deviations holding recruitment outright. Set dont_pen_recdev_first = 0 to keep the ",
         "level readable instead.")

  # with the first year still in the penalty the level is readable, but only through that one term,
  # so it is a weak statement rather than an unidentified one
  if(rec_model == "mean_rec" && RecDevs_model == "rw" && ln_global_R0_spec == "est" &&
     dont_pen_recdev_first == 0)
    warning("rec_model = 'mean_rec' with RecDevs_model = 'rw' leaves the level of log recruitment ",
            "readable only through the first year's deviation, whose standard deviation is ",
            "RecDevs_rw_init_sigma (", RecDevs_rw_init_sigma, "). ln_global_R0 is then weakly identified ",
            "and its standard error will come back near that value. Consider ln_global_R0_spec = 'fix'.")

  # dsem needs a nonzero code, or the deviation never reaches recruitment
  # the dsem sets those cells to NA in map_ln_RecDevs, so the penalty doesn't use them and the dsem supplies their density
  input_list$data$RecDevs_model <- match(if(RecDevs_model == "dsem") "iid" else RecDevs_model, c("iid", "rw", "ar1")) # 1 = iid, 2 = rw, 3 = ar1
  input_list$data$RecDevs_rw_init_sigma <- RecDevs_rw_init_sigma
  if(!Use_rec_level_pen %in% c(0,1)) stop("Use_rec_level_pen must be 0 or 1")
  if(!rec_level_pen_center %in% c("fixed", "own_mean")) stop("rec_level_pen_center must be fixed or own_mean")
  input_list$data$Use_rec_level_pen <- Use_rec_level_pen
  input_list$data$ln_sigma_rec_level <- log(rec_level_pen_sigma)
  input_list$data$rec_level_pen_center <- convert_to_numeric(rec_level_pen_center, list(fixed = 0, own_mean = 1))
  input_list$data$rec_level_pen_yrs <- if(is.null(rec_level_pen_yrs)) rep(1, length(input_list$data$years)) else as.numeric(input_list$data$years %in% rec_level_pen_yrs)
  if(Use_rec_level_pen == 1) collect_message("A recruitment level penalty is applied, centered on: ", rec_level_pen_center)
  if(!Use_init_sex_pen %in% c(0,1)) stop("Use_init_sex_pen must be 0 or 1")
  if(Use_init_sex_pen == 1 && (input_list$data$n_sexes == 1 || InitDevs_sex_spec != "est_all"))
    stop("Use_init_sex_pen ties each sex's initial age deviations to the first sex's, so it needs n_sexes > 1 and InitDevs_sex_spec = 'est_all'. Under est_shared_s the sexes already read one curve and the tie is identically zero.")
  input_list$data$Use_init_sex_pen <- Use_init_sex_pen
  input_list$data$ln_sigma_init_sex <- log(init_sex_pen_sigma)
  if(Use_init_sex_pen == 1) collect_message("Each sex's initial age deviations are tied to the first sex's at sigma ", init_sex_pen_sigma)
  if(!sr_penalty %in% c("none", "bh", "ricker")) stop("sr_penalty must be none, bh, or ricker")
  if(!sr_R0_spec %in% c("shared", "est", "rinit")) stop("sr_R0_spec must be shared, est, or rinit")
  if(sr_R0_spec == "rinit" && use_rinit != 1)
    stop("sr_R0_spec = 'rinit' takes the curve's scale from ln_rinit, so it needs use_rinit = 1. With use_rinit = 0 the initial age structure is built from ln_global_R0 and 'shared' is the same thing.")
  if(sr_penalty != "none" && input_list$data$rec_model != 0)
    stop("sr_penalty is only valid with rec_model = 'mean_rec'. Under bh_rec or ricker_rec the stock-recruit curve already generates recruitment, and penalizing the residual as well would penalize it twice.")
  input_list$data$sr_penalty <- convert_to_numeric(sr_penalty, list(none = 0, bh = 1, ricker = 2))
  input_list$data$sr_R0_spec <- convert_to_numeric(sr_R0_spec, list(shared = 0, est = 1, rinit = 2))
  input_list$data$ln_sigma_sr_pen <- log(sr_pen_sigma)
  # a penalty year needs a spawning biomass behind it. the first rec_lag years have none and fall
  # back to the fished equilibrium, so they are dropped by default and rejected if asked for
  sr_yr_ok <- seq_along(input_list$data$years) > rec_lag
  input_list$data$sr_pen_yrs <- if(is.null(sr_pen_yrs)) as.numeric(sr_yr_ok) else as.numeric(input_list$data$years %in% sr_pen_yrs)
  if(sr_penalty != "none" && any(input_list$data$sr_pen_yrs == 1 & !sr_yr_ok))
    stop("sr_pen_yrs includes the first ", rec_lag, " year(s) of the model, which have no lagged spawning biomass. The stock-recruit prediction there falls back to the fished equilibrium, so the residual would not be a stock-recruit residual. Drop those years.")
  if(sr_penalty != "none") collect_message("Recruitment is a mean with deviations; a ", sr_penalty, " curve is fitted as a penalty on the residual, with its scale ", sr_R0_spec)
  input_list$data$InitDevs_pen_center <- convert_to_numeric(InitDevs_pen_center, list(fixed = 0, own_mean = 1))
  collect_message("Recruitment deviation penalty is centered on: ", RecDevs_pen_center)
  input_list$data$init_F_prop <- init_F_prop
  input_list$data$init_F_form <- init_F_form_num
  input_list$data$t_spawn <- t_spawn
  input_list$data$equil_init_age_strc <- equil_init_age_strc
  input_list$data$max_bias_ramp_fct <- max_bias_ramp_fct
  input_list$data$use_rec_region_prop_prior <- use_rec_region_prop_prior
  input_list$data$rec_region_prop_spec <- if(is.null(rec_region_prop_spec)) 0 else 1 # 0 = Full dispersal, 1 = no dispersal
  input_list$data$rec_region_prop_prior <- rec_region_prop_prior
  input_list$data$use_fixed_stray_rate <- use_fixed_stray_rate
  input_list$data$fixed_stray_rate     <- fixed_stray_rate
  input_list$data$stray_rate_blocks    <- stray_rate_blocks_mat
  input_list$data$sexratio_blocks <- sexratio_blocks_mat
  input_list$data$use_fixed_rec_seas_prop <- use_fixed_rec_seas_prop
  input_list$data$fixed_rec_seas_prop <- fixed_rec_seas_prop
  input_list$data$use_rec_seas_prop_prior <- use_rec_seas_prop_prior
  input_list$data$rec_seas_prop_prior <- rec_seas_prop_prior
  input_list$data$use_stray_rate_prior <- use_stray_rate_prior
  input_list$data$stray_rate_prior     <- stray_rate_prior
  input_list$data$use_rinit <- use_rinit
  if(!Use_rinit_pen %in% c(0, 1)) stop("Use_rinit_pen must be 0 or 1")
  if(Use_rinit_pen == 1 && use_rinit != 1) stop("Use_rinit_pen = 1 penalizes ln_rinit against ln_global_R0, so it needs use_rinit = 1")
  if(Use_rinit_pen == 1 && (length(rinit_pen_sd) != 1 || rinit_pen_sd <= 0)) stop("rinit_pen_sd must be a single positive number")
  input_list$data$Use_rinit_pen <- Use_rinit_pen
  input_list$data$rinit_pen_sd <- rinit_pen_sd
  if(Use_rinit_pen == 1) collect_message("Initial recruitment offset from R0 is penalized with sd ", rinit_pen_sd)
  input_list$data$init_age_devs_shared <- init_age_devs_shared
  input_list$data$use_r0_prior <- use_r0_prior
  input_list$data$r0_prior     <- r0_prior

  # Populate Parameter List -------------------------------------------------

  # global R0, one per population; use_starting_value catches a wrong-length starting value here.
  # time blocks add a column, so a model with no blocks keeps ln_global_R0 a length-n_pop vector
  n_yrs_r0 <- length(input_list$data$years)
  n_pop_r0 <- input_list$data$n_pop
  if(is.null(R0_blocks)) R0_blocks <- paste0("none_Pop_", seq_len(n_pop_r0))
  R0_blocks_arr <- array(NA, dim = c(1, n_yrs_r0, n_pop_r0))
  for(str in R0_blocks) {
    v <- unlist(strsplit(str, "_"))
    if(!v[1] %in% c("none", "Block")) stop("R0_blocks must be none_Pop_p or Block_b_Year_a-e_Pop_p")
    if(v[1] == "none") R0_blocks_arr[, , as.numeric(v[3])] <- 1
    if(v[1] == "Block") {
      pp <- as.numeric(v[6])
      bv <- as.numeric(v[2])
      rng <- unlist(strsplit(v[4], "-"))
      yy <- as.numeric(rng[1]):(if(rng[2] == "terminal") n_yrs_r0 else as.numeric(rng[2]))
      R0_blocks_arr[, yy, pp] <- bv
    }
  } # end str loop
  if(any(is.na(R0_blocks_arr))) stop("R0_blocks leaves some years unassigned for at least one population")
  n_R0_blks <- max(1, max(R0_blocks_arr))
  if(R0_ref_block < 1 || R0_ref_block > n_R0_blks) stop("R0_ref_block must lie within 1..", n_R0_blks)
  input_list$data$R0_blocks <- R0_blocks_arr
  input_list$data$R0_ref_block <- as.integer(R0_ref_block)

  input_list$par$ln_global_R0 <- array(log(15), dim = c(input_list$data$n_pop, n_R0_blks))
  # a starting value written for the old shape is a plain length-n_pop vector; give every
  # block that value rather than refusing it
  if(!is.null(starting_values$ln_global_R0) && is.null(dim(starting_values$ln_global_R0)) &&
     length(starting_values$ln_global_R0) == n_pop_r0 && n_R0_blks > 1)
    starting_values$ln_global_R0 <- matrix(rep(starting_values$ln_global_R0, n_R0_blks), n_pop_r0, n_R0_blks)
  # only reshape a value whose LENGTH already matches; anything else must still reach
  # use_starting_value so its own guard names the argument and the expected size
  if(!is.null(starting_values$ln_global_R0) && n_R0_blks == 1 &&
     length(starting_values$ln_global_R0) == n_pop_r0)
    starting_values$ln_global_R0 <- array(as.vector(starting_values$ln_global_R0), dim = c(n_pop_r0, 1L))
  input_list$par$ln_global_R0 <- use_starting_value(input_list$par$ln_global_R0, starting_values, "ln_global_R0")
  # a starting value supplied as a plain vector loses the block dimension, and every
  # downstream reader indexes ln_global_R0[pop, block]
  input_list$par$ln_global_R0 <- array(as.vector(input_list$par$ln_global_R0), dim = c(n_pop_r0, n_R0_blks))

  # The stock-recruit curve's own scale, used only when sr_R0_spec = "est".
  # Mapped off otherwise so it never enters the parameter vector by accident.
  input_list$par$ln_sr_R0 <- array(input_list$par$ln_global_R0[, R0_ref_block], dim = c(input_list$data$n_pop))
  input_list$par$ln_sr_R0 <- use_starting_value(input_list$par$ln_sr_R0, starting_values, "ln_sr_R0")
  input_list$map$ln_sr_R0 <- if(input_list$data$sr_R0_spec == 1 && input_list$data$sr_penalty > 0) {
    factor(seq_along(input_list$par$ln_sr_R0))
  } else factor(rep(NA, length(input_list$par$ln_sr_R0)))

  # Global Initial R0
  input_list$par$ln_rinit <- array(log(15), dim = c(input_list$data$n_pop))
  input_list$par$ln_rinit <- use_starting_value(input_list$par$ln_rinit, starting_values, "ln_rinit")
  if (use_rinit == 0) input_list$map$ln_rinit <- factor(rep(NA, input_list$data$n_pop))

  # R0 regional proportion (not availiable when n_regions == 1; altered in do_rec_region_prop_mapping)
  input_list$par$rec_region_prop_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions - 1))
  input_list$par$rec_region_prop_pars <- use_starting_value(input_list$par$rec_region_prop_pars, starting_values, "rec_region_prop_pars")

  # R0 seasonal proportion (not availiable when n_seas == 1; altered in do_rec_seas_prop_mapping)
  input_list$par$rec_seas_prop_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_seas - 1))
  input_list$par$rec_seas_prop_pars <- use_starting_value(input_list$par$rec_seas_prop_pars, starting_values, "rec_seas_prop_pars")

  # Steepness in bounded logit space (0.2 and 1)
  input_list$par$steepness_h <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions))
  input_list$par$steepness_h <- use_starting_value(input_list$par$steepness_h, starting_values, "steepness_h")

  # Initial age deviations, one age curve per sex; a 3-D starting value array
  # without the sex dimension is broadcast across sexes
  if("ln_InitDevs" %in% names(starting_values)) {
    input_list$par$ln_InitDevs <- starting_values$ln_InitDevs
    if(length(dim(input_list$par$ln_InitDevs)) == 3) input_list$par$ln_InitDevs <- array(rep(input_list$par$ln_InitDevs, input_list$data$n_sexes), dim = c(dim(input_list$par$ln_InitDevs), input_list$data$n_sexes))
  }
  else input_list$par$ln_InitDevs <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$ages) - 1, input_list$data$n_sexes))

  # Recruitment deviations
  input_list$par$ln_RecDevs <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years) - dont_est_recdev_last + input_list$data$n_proj_yrs_devs))
  input_list$par$ln_RecDevs <- use_starting_value(input_list$par$ln_RecDevs, starting_values, "ln_RecDevs")

  # Recruitment variability
  input_list$par$ln_sigmaR <- array(log(1), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)) # (early period 1st element, late period 2nd element)
  input_list$par$ln_sigmaR <- use_starting_value(input_list$par$ln_sigmaR, starting_values, "ln_sigmaR")

  # sexratio parameters
  max_sexratio_blks <- max(apply(input_list$data$sexratio_blocks, c(1,2), FUN = function(x) length(unique(x)))) # figure out maximum number ofsex ratio blocks for each region
  input_list$par$sexratio_pars <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, max_sexratio_blks)) # specified at 0.5 in inverse logit space
  input_list$par$sexratio_pars <- use_starting_value(input_list$par$sexratio_pars, starting_values, "sexratio_pars")

  # Stray rate parameters (logit scale)
  max_stray_blks <- if (input_list$data$n_pop > 1) max(apply(stray_rate_blocks_mat, 1, function(x) length(unique(x))))
  else 1

  # ar1 correlation for recruitment deviations (only read when RecDevs_model = 'ar1')
  input_list$par$RecDevs_rho <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions))
  input_list$par$RecDevs_rho <- use_starting_value(input_list$par$RecDevs_rho, starting_values, "RecDevs_rho")

  max_stray_blks <- if (input_list$data$n_pop > 1)  max(apply(stray_rate_blocks_mat, 1, function(x) length(unique(x)))) else 1
  input_list$par$stray_rate_pars <- array(0,  dim = c(input_list$data$n_pop, max_stray_blks))
  input_list$par$stray_rate_pars <- use_starting_value(input_list$par$stray_rate_pars, starting_values, "stray_rate_pars")

  # Mapping Options -----------------------------------------------------------

  input_list <- do_rec_region_prop_mapping(input_list, rec_region_prop_spec) # Recruitment regional proportion mapping
  input_list <- do_rec_seas_prop_mapping(input_list, rec_seas_prop_spec) # Recruitment seasonal proportion mapping
  input_list <- do_sigmaR_mapping(input_list, sigmaR_spec) # sigmaR mapping
  input_list <- do_InitDevs_mapping(input_list, InitDevs_spec, rec_dd, init_age_devs_shared, InitDevs_sex_spec) # InitDevs mapping
  # the shared-subset penalty (equil_init_age_strc == 3) reads this in the model,
  # so it goes into data as well as the map
  if(!is.null(init_age_devs_shared)) input_list$data$init_age_devs_shared <- init_age_devs_shared
  input_list <- do_RecDevs_mapping(input_list, RecDevs_spec, rec_dd, dont_pen_recdev_first) # RevDevs mapping
  input_list <- do_RecDevs_rho_mapping(input_list, RecDevs_rho_spec) # recruitment deviation ar1 correlation mapping

  # unfished recruitment, estimated per population and R0 block or fixed at its starting value
  n_R0_par <- length(input_list$par$ln_global_R0)
  input_list$map$ln_global_R0 <- factor(if(ln_global_R0_spec == "est") seq_len(n_R0_par) else rep(NA_integer_, n_R0_par))
  input_list <- do_h_mapping(input_list, h_spec, rec_dd) # steepness mapping
  input_list <- do_sexratio_pars_mapping(input_list, sexratio_spec) # sex ratio parameters
  input_list <- do_stray_rate_mapping(input_list, stray_rate_spec) # stray rates

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}
