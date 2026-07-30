# Stage 3 of 3: post fit
#
# Forward projection off a fitted model under a specified catch or fishing
# mortality, with optional stochastic recruitment. This is the short horizon
# projection used for advice, not the closed loop simulation in
# sim_closed_loop.R.

#' Do Population Projections
#'
#' Projects population dynamics forward in time under alternative recruitment
#' and fishing mortality scenarios. The model initializes from terminal
#' assessment quantities and advances numbers-at-age through recruitment,
#' seasonal movement, mortality, ageing, and harvest control rules across
#' multiple seasons and years.
#'
#' Population dynamics are tracked at full resolution over
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
#' @param natmort Array `[n_pop, n_regions, n_proj_yrs, n_ages, n_sexes]`.
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
#'   to redistribute fish to natal grounds prior to SSB calculation.
#' @param stray_rate Array `[n_pop, n_proj_yrs]`. Per-population stray rate
#'   used when accumulating effective SSB contributions across populations.
#' @param f_ref_pt Array `[n_regions, n_proj_yrs]`. Fishing mortality
#'   reference point (e.g., F_MSY) or fixed input F, depending on
#'   `fmort_opt`.
#' @param b_ref_pt Array `[n_pop, n_regions, n_proj_yrs]`. Biomass reference
#'   point used in harvest control rules.
#' @param HCR_function Function. Harvest control rule with arguments `x`
#'   (SSB), `frp` (F reference point), and `brp` (B reference point).
#' @param recruitment_opt Character. Recruitment scenario:
#'   `"inv_gauss"`, `"mean_rec"`, `"zero"`, or `"bh_rec"`.
#' @param fmort_opt Character. Fishing mortality scenario:
#'   `"HCR"`, `"HCR_global"`, or `"Input"`.
#' @param t_spawn Numeric scalar. Fraction of the spawning season elapsed
#'   before spawning; used for mid-season SSB calculations.
#' @param bh_rec_opt Named list of inputs for deterministic Beverton–Holt
#'   recruitment when `recruitment_opt = "bh_rec"`. This list is passed
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
#'       \code{[n_pop, n_regions, n_ages]}. Natural mortality.}
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
#'   combining \code{bh_rec_opt$SSB} with projected SSB values during the
#'   simulation.
#'
#'   \code{bh_rec_opt$rec_lag = 1} is the classic lagged case: each
#'   projection year's recruitment is computed up front from the prior
#'   year's SSB, exactly as \code{recruitment_opt = "inv_gauss"}/
#'   \code{"mean_rec"} are. \code{bh_rec_opt$rec_lag = 0} is age-0
#'   recruitment: recruitment for year \code{y} is computed from year
#'   \code{y}'s own SSB once \code{spawn_seas} is reached within that year's
#'   season loop, and is inserted no earlier than \code{spawn_seas}
#'   (\code{rec_seas_prop} must be zero for every season before
#'   \code{spawn_seas} in that case). Reference points and the seasonal SBPR
#'   calculation used to get \code{bh_rec_opt$WAA}/\code{MatAA}/etc. are
#'   unaffected by this choice -- \code{rec_lag} only changes which year's
#'   SSB feeds the Beverton-Holt curve, not the per-recruit math itself.
#'
#' @param n_seas Integer. Number of seasons. Default = 1.
#' @param seasdur Numeric vector `[n_seas]`. Duration of each season as a
#'   fraction of the year.
#' @param spawn_seas Integer. Spawning season index.
#' @param natal_region Integer vector `[n_pop]`. Natal region for each
#'   population.
#' @param dmr Array \code{[n_regions, n_seas, n_fish_fleets]}. Discard mortality rate.
#'   Default behavior is no discard mortality (\code{dmr = 0}). When combined with
#'   \code{ret_sel = 1}, this implies no discarding within a given fleet (all catch is retained).
#' @param ret_sel Array \code{[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]}. Retention
#'   selectivity-at-age. Default behavior corresponds to full retention (\code{ret_sel = 1}),
#'   meaning all captured fish are retained unless otherwise specified.
#'
#' @return A named list containing projected fishing mortality, catch,
#'   spawning biomass, effective spawning biomass, dynamic unfished biomass,
#'   and numbers-at-age for fished and unfished states.
#'
#' @details
#' Each projection year proceeds as follows when
#' \code{recruitment_opt != "bh_rec"} or \code{bh_rec_opt$rec_lag != 0}
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
#' When \code{bh_rec_opt$rec_lag == 0} (age-0 recruitment), steps 1 and 5
#' above are reordered within \code{spawn_seas}: movement is applied first,
#' spawning biomass is computed from the survivor population alone (no new
#' recruits exist yet), that SSB is used to generate this year's
#' recruitment, and only then are the recruits inserted (no earlier than
#' \code{spawn_seas}) - immediately before mortality/ageing runs for that
#' season, so the new cohort is carried forward exactly like any other
#' seasonal recruit pulse. Years \code{y > 1} generate recruitment this way;
#' year 1 carries the supplied terminal assessment state forward with no new
#' recruitment event, matching the classic case.
#'
#' Effective spawning biomass at each population's natal region aggregates
#' contributions from all populations, with cross-population contributions
#' scaled by \code{stray_rate} and normalised by the number of populations
#' in each natal region.
#'
#' When \code{n_sexes = 1}, spawning biomass is multiplied by 0.5. When
#' \code{n_regions = 1}, movement is skipped.
#'
#' @export
#' @family Reference Points and Projections
#' @import abind abind
Do_Population_Projection <- function(n_proj_yrs = 2,
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
                                     natal_region,
                                     WAA,
                                     WAA_fish,
                                     MatAA,
                                     fish_sel,
                                     ret_sel = array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
                                     Movement,
                                     sgl_seas_spawning_movement,
                                     stray_rate,
                                     f_ref_pt = NULL,
                                     b_ref_pt = NULL,
                                     HCR_function = NULL,
                                     recruitment_opt = "inv_gauss",
                                     fmort_opt = 'HCR',
                                     t_spawn,
                                     bh_rec_opt = NULL,
                                     n_seas = 1,
                                     seasdur = rep(1 / n_seas, n_seas),
                                     spawn_seas = 1,
                                     rec_seas_prop = {
                                       rec_seas_prop = array(0, dim = c(n_pop, n_seas))
                                       rec_seas_prop[] <- 1 / n_seas
                                       rec_seas_prop
                                     },
                                     Mrate = NULL,
                                     move_timing = 0
) {


  # Error Checking ----------------------------------------------------------

  if(!recruitment_opt %in% c("inv_gauss", "mean_rec", "zero", "bh_rec")) stop("Recruitment options are not specified correctly! Should be inv_gauss, mean_rec, zero, or bh_rec")
  if(!fmort_opt %in% c("HCR", "Input", "HCR_global")) stop("Fishing Mortality options are not specified correctly! Should be HCR, Input, HCR_global")
  if(recruitment_opt == "bh_rec") {
    required_fields <- c("rec_dd", "rec_lag", "R0", "h", "rec_region_prop",
                         "WAA", "MatAA", "natmort", "SSB", "Movement",
                         "sex_ratio_f", "stray_rate", "fish_sel", "ret_sel", "dmr", "init_F")
    diff <- setdiff(required_fields, names(bh_rec_opt)) # find difference
    if(length(diff) > 0) stop(paste("bh_rec_opt is missing the following required fields:", paste(diff)))
  }

  # Define Containers -------------------------------------------------------
  fratio <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  for(r in 1:n_regions) for(seas in 1:n_seas) for(f in 1:n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
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
  proj_Dynamic_SSB0 <- array(0, dim = c(n_pop, n_regions, n_proj_yrs))
  proj_F <- array(0, dim = c(n_regions, n_proj_yrs + 1))

  # Start Projection --------------------------------------------------------
  # Input terminal year assessment at age
  proj_NAA[,,1,,,] <- terminal_NAA
  proj_NAA0[,,1,,,] <- terminal_NAA0

  # Age-0 (rec_lag = 0) BH recruitment: this year's own SSB determines this
  # year's recruitment, which isn't known until spawn_seas is reached within
  # the season loop - see the "rec_lag == 0" handling below, mirroring the
  # equivalent restructuring in SPoRC_rtmb.R and Simulate_Population.R.
  age0_bh <- recruitment_opt == "bh_rec" && !is.null(bh_rec_opt) && bh_rec_opt$rec_lag == 0

  for(y in 1:n_proj_yrs) {

    # use terminal F in the first year (subsequent years use F derived from reference points and HCR)
    if(y == 1) proj_F[,y] <- rowSums(terminal_F)

    # Recruitment Processes (rec_lag != 0, or non-BH recruitment) -------------
    # For age0_bh, recruitment for the year is instead generated inline once
    # spawn_seas is reached within the season loop below.
    if(y > 1 && !age0_bh) {

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

                        "bh_rec" = { # if beverton holt recruitment
                          Get_Det_Recruitment(recruitment_model = 1,
                                              rec_dd = bh_rec_opt$rec_dd,
                                              n_pop = n_pop,
                                              sgl_seas_spawning_movement = bh_rec_opt$sgl_seas_spawning_movement,
                                              natal_region = natal_region,
                                              y = y + dim(bh_rec_opt$SSB)[3],
                                              rec_lag = bh_rec_opt$rec_lag,
                                              R0 = bh_rec_opt$R0,
                                              rec_region_prop = bh_rec_opt$rec_region_prop,
                                              rec_seas_prop = rec_seas_prop,
                                              h = bh_rec_opt$h,
                                              n_regions = n_regions,
                                              n_ages = n_ages,
                                              WAA = bh_rec_opt$WAA,
                                              MatAA = bh_rec_opt$MatAA,
                                              n_seas = n_seas,
                                              seasdur = seasdur,
                                              spawn_seas = spawn_seas,
                                              natmort = bh_rec_opt$natmort,
                                              SSB_vals = abind::abind(bh_rec_opt$SSB, proj_SSB, along = 3),
                                              Movement = bh_rec_opt$Movement,
                                              # SSB0 behind the stock recruit curve has to use the same movement
                                              # sequencing as the projection itself, so forward both of these.
                                              Mrate = bh_rec_opt$Mrate,
                                              stray_rate = bh_rec_opt$stray_rate,
                                              do_recruits_move = do_recruits_move,
                                              t_spawn = t_spawn,
                                              sexratio_f = bh_rec_opt$sex_ratio_f,
                                              init_F = bh_rec_opt$init_F,
                                              n_fish_fleets = n_fish_fleets,
                                              fish_sel = bh_rec_opt$fish_sel,
                                              ret_sel = bh_rec_opt$ret_sel,
                                              dmr = bh_rec_opt$dmr,
                                              move_timing = move_timing
                          )
                        }
      )

      # coerce into array
      tmp_rec <- array(tmp_rec, dim = c(n_pop, n_regions))

      # Apply recruitment to projected proj_NAA
      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          tmp <- tmp_rec[p,r] * sexratio[p,r,y,] * rec_seas_prop[p,1]
          proj_NAA[p,r,y,1,1,] <- proj_NAA0[p,r,y,1,1,]  <- tmp
        } # end r loop
      } # end p loop

    } # if y > 1

    for(seas in 1:n_seas) {

      # Insert seasonal recruits already known from earlier this year:
      # - rec_lag != 0 (or non-BH recruitment): the year's recruitment is
      #   already known (computed above), so any season past the first gets
      #   its share here, as before.
      # - age0_bh (rec_lag == 0): recruitment isn't known until spawn_seas is
      #   reached (below), so only seasons strictly after spawn_seas are
      #   handled here; spawn_seas itself generates and inserts its own share.
      if(y > 1 && (if(age0_bh) seas > spawn_seas else seas > 1)) {
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
                proj_ret_FAA[p,r,y,seas,a,s,f] <- proj_F[r,y] * fratio[r,seas,f] * fish_sel[p,r,y,seas,a,s,f] * ret_sel[p,r,y,seas,a,s,f] # retained F
                proj_disc_FAA[p,r,y,seas,a,s,f] <- proj_F[r,y] * fratio[r,seas,f] * fish_sel[p,r,y,seas,a,s,f] * (1 - ret_sel[p,r,y,seas,a,s,f]) * dmr[r,seas,f] # discarded F
                proj_tot_FAA[p,r,y,seas,a,s,f] <- proj_ret_FAA[p,r,y,seas,a,s,f] + proj_disc_FAA[p,r,y,seas,a,s,f] # total F
              } # end p loop
            } # end f loop

            # Get Total Mortality at Age
            for(p in 1:n_pop) {
              proj_ZAA[p,r,y,seas,a,s] <- (natmort[p,r,y,a,s] * seasdur[seas]) + sum(proj_tot_FAA[p,r,y,seas,a,s,])
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

      # Derive Biomass + Recruitment (age0_bh only) ------------------------------
      # This year's SSB is now fully determined by the survivor population
      # (age-0 recruits have not been produced yet, and couldn't affect SSB
      # even if they had - rec_lag == 0 requires MatAA == 0 at the recruit
      # age). Compute it now, generate this year's recruitment from it, and
      # insert the spawn_seas share BEFORE mortality/ageing runs below, so the
      # new cohort is carried forward exactly like any other seasonal recruit
      # pulse.
      if(age0_bh && seas == spawn_seas) {

        biom <- derive_proj_biom(y, seas, proj_NAA, proj_NAA0, WAA, MatAA, proj_ZAA, natmort, t_spawn, seasdur,
                                n_seas, n_pop, n_regions, n_ages, n_sexes,
                                sgl_seas_spawning_movement, natal_region, stray_rate,
                                Movement, Mrate, move_timing, do_recruits_move)
        proj_SSB[,, y] <- biom$SSB_y
        proj_Dynamic_SSB0[,,y] <- biom$Dynamic_SSB0_y
        proj_eff_SSB[,y] <- biom$eff_SSB_y

        if(y > 1) {

          tmp_rec <- Get_Det_Recruitment(recruitment_model = 1,
                                         rec_dd = bh_rec_opt$rec_dd,
                                         n_pop = n_pop,
                                         sgl_seas_spawning_movement = bh_rec_opt$sgl_seas_spawning_movement,
                                         natal_region = natal_region,
                                         y = y + dim(bh_rec_opt$SSB)[3],
                                         rec_lag = bh_rec_opt$rec_lag,
                                         R0 = bh_rec_opt$R0,
                                         rec_region_prop = bh_rec_opt$rec_region_prop,
                                         rec_seas_prop = rec_seas_prop,
                                         h = bh_rec_opt$h,
                                         n_regions = n_regions,
                                         n_ages = n_ages,
                                         WAA = bh_rec_opt$WAA,
                                         MatAA = bh_rec_opt$MatAA,
                                         n_seas = n_seas,
                                         seasdur = seasdur,
                                         spawn_seas = spawn_seas,
                                         natmort = bh_rec_opt$natmort,
                                         SSB_vals = abind::abind(bh_rec_opt$SSB, proj_SSB, along = 3),
                                         Movement = bh_rec_opt$Movement,
                                         # SSB0 behind the stock recruit curve has to use the same movement
                                         # sequencing as the projection itself, so forward both of these.
                                         Mrate = bh_rec_opt$Mrate,
                                         stray_rate = bh_rec_opt$stray_rate,
                                         do_recruits_move = do_recruits_move,
                                         t_spawn = t_spawn,
                                         sexratio_f = bh_rec_opt$sex_ratio_f,
                                         init_F = bh_rec_opt$init_F,
                                         n_fish_fleets = n_fish_fleets,
                                         fish_sel = bh_rec_opt$fish_sel,
                                         ret_sel = bh_rec_opt$ret_sel,
                                         dmr = bh_rec_opt$dmr,
                                         move_timing = move_timing
          )
          tmp_rec <- array(tmp_rec, dim = c(n_pop, n_regions))

          for(p in 1:n_pop) {
            for(r in 1:n_regions) {
              proj_NAA[p,r,y,spawn_seas,1,]  <- proj_NAA[p,r,y,spawn_seas,1,]  + tmp_rec[p,r] * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,]
              proj_NAA0[p,r,y,spawn_seas,1,] <- proj_NAA0[p,r,y,spawn_seas,1,] + tmp_rec[p,r] * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,]
            } # end r loop
          } # end p loop

          # Recruits just inserted above missed this season's movement step
          # (which already ran, since it had to happen before this year's
          # SSB, and hence recruitment, was knowable). Catch age index 1 up
          # to the rest of the cohort when recruits are supposed to move from
          # birth.
          # Only needed under move_timing == 0; under timings 1 and 2 these recruits are
          # picked up by the end-of-season transition below.
          if(do_recruits_move == 1 && n_regions > 1 && move_timing == 0) {
            for(p in 1:n_pop) {
              for(s in 1:n_sexes) proj_NAA[p,,y,seas,1,s] = t(proj_NAA[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
              for(s in 1:n_sexes) proj_NAA0[p,,y,seas,1,s] = t(proj_NAA0[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
            } # end p loop
          }

        } # end if y > 1

      } # end if age0_bh && seas == spawn_seas

      # Movement (timing 1 and 2), Mortality and Ageing --------------------------
      # Post-season state at every age, before the ageing shift. Under move_timing == 0
      # movement was applied above so this reduces to the original elementwise survival.
      if(move_timing == 0 || n_regions == 1) {
        # array() guards against R dropping length-1 pop/region/sex dimensions
        pstep_NAA <- array(proj_NAA[,,y,seas,1:n_ages,] * exp(-proj_ZAA[,,y,seas,1:n_ages,]),
                           dim = c(n_pop, n_regions, n_ages, n_sexes))
        pstep_NAA0 <- array(proj_NAA0[,,y,seas,1:n_ages,] * exp(-natmort[,,y,1:n_ages,] * seasdur[seas]),
                            dim = c(n_pop, n_regions, n_ages, n_sexes))
      } else {
        pstep_NAA <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        pstep_NAA0 <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
          moves <- (do_recruits_move == 1 || a > 1)
          Mv <- if(moves) Movement[p,,,y,seas,a,s] else diag(n_regions)
          Qv <- if(moves) Mrate[p,,,y,seas,a,s] else matrix(0, n_regions, n_regions)
          pstep_NAA[p,,a,s] <- advance_seas(proj_NAA[p,,y,seas,a,s], Mv, proj_ZAA[p,,y,seas,a,s],
                                            Qv, seasdur[seas], move_timing)
          pstep_NAA0[p,,a,s] <- advance_seas(proj_NAA0[p,,y,seas,a,s], Mv, natmort[p,,y,a,s] * seasdur[seas],
                                             Qv, seasdur[seas], move_timing)
        }
      }

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

      # Derive Biomass (age0_bh: already computed above, before mortality/ageing) --
      if(seas == spawn_seas && !age0_bh) {
        biom <- derive_proj_biom(y, seas, proj_NAA, proj_NAA0, WAA, MatAA, proj_ZAA, natmort, t_spawn, seasdur,
                                n_seas, n_pop, n_regions, n_ages, n_sexes,
                                sgl_seas_spawning_movement, natal_region, stray_rate,
                                Movement, Mrate, move_timing, do_recruits_move)
        proj_SSB[,, y] <- biom$SSB_y
        proj_Dynamic_SSB0[,,y] <- biom$Dynamic_SSB0_y
        proj_eff_SSB[,y] <- biom$eff_SSB_y
      } # calculate biomass


      # Season-integrated abundance for the spatial Baranov under continuous movement.
      # Computed once per season across all regions, since the integral couples them.
      if(move_timing == 2) {
        proj_NAA_int <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        for(p in 1:n_pop) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              proj_NAA_int[p,,a,s] <- integrate_seas_abundance(proj_NAA[p,,y,seas,a,s], proj_ZAA[p,,y,seas,a,s],
                                                              Mrate[p,,,y,seas,a,s], seasdur[seas])
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
                } else {
                  # Get catch at age with Baranov's
                  proj_CAA[p,r,y,seas,a,s,f] <- (proj_ret_FAA[p,r,y,seas,a,s,f] / proj_ZAA[p,r,y,seas,a,s]) *
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

    # compute F for next year
    for(r in 1:n_regions) {

      # Project F using HCR and reference points -----------------------------------------------------
      if(fmort_opt == 'HCR') {
        proj_F[r,y+1] <- HCR_function(x = sum(proj_SSB[,r,y]),
                                      frp = f_ref_pt[r,y],
                                      brp = sum(b_ref_pt[,r,y]))
      }

      if(fmort_opt == 'HCR_global') {
        proj_F[r,y+1] <- HCR_function(x = sum(proj_SSB[,,y]),
                                      frp = f_ref_pt[r,y],
                                      brp = sum(b_ref_pt[,,y]))
      }

      # Project F using User Inputs ---------------------------------------------
      if(fmort_opt == 'Input') proj_F[r,y+1] <- f_ref_pt[r,y]

    } # end r loop

  } # end y loop

  return(list(proj_F = proj_F,
              proj_ret_FAA = proj_ret_FAA,
              proj_disc_FAA = proj_disc_FAA,
              proj_Catch = proj_Catch,
              proj_SSB = proj_SSB,
              proj_eff_SSB = proj_eff_SSB,
              proj_Dynamic_SSB0 = proj_Dynamic_SSB0,
              proj_NAA = proj_NAA,
              proj_NAA0 = proj_NAA0,
              proj_ZAA = proj_ZAA)
  )

} # end function

