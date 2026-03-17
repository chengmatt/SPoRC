#' Do Population Projections
#'
#' Projects population dynamics forward in time under alternative recruitment
#' and fishing mortality scenarios. Initializes from terminal assessment
#' quantities and advances numbers-at-age through survival, movement,
#' recruitment, ageing, and harvest control rules across seasons and years.
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
#'   recruitment entering in each season.
#' @param recruitment Array `[n_pop, n_regions, n_yrs]`. Historical
#'   recruitment used to condition stochastic projection options.
#' @param terminal_NAA Array `[n_pop, n_regions, n_seas, n_ages, n_sexes]`.
#'   Fished numbers-at-age in the terminal assessment year.
#' @param terminal_NAA0 Array `[n_pop, n_regions, n_seas, n_ages, n_sexes]`.
#'   Unfished numbers-at-age in the terminal assessment year.
#' @param terminal_F Array `[n_regions, n_seas, n_fish_fleets]`. Terminal
#'   fishing mortality; sets F in projection year 1 and defines the seasonal
#'   F ratio applied in subsequent years.
#' @param natmort Array `[n_pop, n_regions, n_proj_yrs, n_ages, n_sexes]`.
#'   Annual natural mortality-at-age, scaled internally by season duration.
#' @param WAA Array `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]`.
#'   Weight-at-age used in spawning biomass calculations.
#' @param WAA_fish Array
#'   `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]`.
#'   Fishery weight-at-age used in catch biomass calculations.
#' @param MatAA Array `[n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]`.
#'   Maturity-at-age.
#' @param fish_sel Array `[n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets]`.
#'   Fishery selectivity-at-age.
#' @param Movement Array
#'   `[n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes]`.
#'   Seasonal movement transition matrices. Entry `[p, r1, r2, y, s, a, sx]`
#'   is the probability of a fish in population `p` moving from region `r1`
#'   to region `r2`.
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
#'   point. For `"HCR_global"`, values are summed internally across
#'   populations and regions.
#' @param HCR_function Function. Harvest control rule with arguments `x`
#'   (SSB), `frp` (F reference point), and `brp` (B reference point).
#'   Required when `fmort_opt` is `"HCR"` or `"HCR_global"`.
#' @param recruitment_opt Character. Recruitment scenario. One of:
#'   \describe{
#'     \item{`"inv_gauss"`}{Stochastic draws from an inverse Gaussian
#'       distribution conditioned on historical recruitment.}
#'     \item{`"mean_rec"`}{Mean historical recruitment by population and
#'       region.}
#'     \item{`"zero"`}{No recruitment.}
#'     \item{`"bh_rec"`}{Deterministic Beverton-Holt recruitment. Requires
#'       `bh_rec_opt`.}
#'   }
#' @param fmort_opt Character. Fishing mortality scenario. One of:
#'   \describe{
#'     \item{`"HCR"`}{Region-specific harvest control rule applied to
#'       regional SSB.}
#'     \item{`"HCR_global"`}{Harvest control rule applied to global SSB
#'       (summed across all populations and regions).}
#'     \item{`"Input"`}{Fixed F supplied directly via `f_ref_pt`.}
#'   }
#' @param t_spawn Numeric scalar. Fraction of the spawning season elapsed
#'   before spawning; used in the mid-season SSB calculation.
#' @param bh_rec_opt Named list. Required when `recruitment_opt = "bh_rec"`.
#'   Must contain:
#'   \describe{
#'     \item{`rec_dd`}{Density dependence structure: 1 = global SSB,
#'       0 = local SSB.}
#'     \item{`rec_lag`}{Recruitment lag in years.}
#'     \item{`R0`}{Unfished equilibrium recruitment.}
#'     \item{`h`}{Steepness by region.}
#'     \item{`rec_region_prop`}{Proportional recruitment apportionment
#'       by region.}
#'     \item{`WAA`}{Reference weight-at-age, indexing the first year.}
#'     \item{`MatAA`}{Reference maturity-at-age, indexing the first year.}
#'     \item{`natmort`}{Reference natural mortality, indexing the first year.}
#'     \item{`SSB`}{Historical SSB array `[n_pop, n_regions, n_yrs]`.
#'       Projection SSBs are appended internally.}
#'     \item{`Movement`}{Movement array for recruitment distribution,
#'       indexing the first year.}
#'     \item{`sgl_seas_spawning_movement`}{Single-season spawning movement
#'       array, indexing the first year.}
#'     \item{`stray_rate`}{Stray rates by population, indexing the first
#'       year.}
#'   }
#' @param n_seas Integer. Number of seasons. Default = 1.
#' @param seasdur Numeric vector `[n_seas]`. Duration of each season as a
#'   fraction of the year. Default = equal fractions.
#' @param spawn_seas Integer. Spawning season index. Default = 1.
#' @param init_F Numeric vector `[n_seas]`. Initial seasonal F values used
#'   when deriving Beverton-Holt equilibrium quantities.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{`proj_F`}{Array `[n_regions, n_proj_yrs + 1]`. Annual fishing
#'       mortality. Year 1 holds terminal F; subsequent years are HCR- or
#'       input-derived.}
#'     \item{`proj_Catch`}{Array
#'       `[n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets]`. Catch
#'       biomass by population, region, year, season, and fleet.}
#'     \item{`proj_SSB`}{Array `[n_pop, n_regions, n_proj_yrs]`. Spawning
#'       stock biomass.}
#'     \item{`proj_eff_SSB`}{Array `[n_pop, n_proj_yrs]`. Effective SSB at
#'       each population's natal region, accumulating straying contributions
#'       scaled by `stray_rate`.}
#'     \item{`proj_Dynamic_SSB0`}{Array `[n_pop, n_regions, n_proj_yrs]`.
#'       Dynamic unfished spawning biomass.}
#'     \item{`proj_NAA`}{Array
#'       `[n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes]`.
#'       Fished numbers-at-age. Index 1 holds terminal values.}
#'     \item{`proj_NAA0`}{Same dimensions as `proj_NAA`. Unfished
#'       numbers-at-age.}
#'     \item{`proj_ZAA`}{Same dimensions as `proj_NAA`. Total
#'       mortality-at-age.}
#'   }
#'
#' @details
#' Each projection year proceeds through the following steps: (1) recruitment
#' is generated and allocated to populations, regions, sexes, and seasons;
#' (2) fishing mortality-at-age is assembled from the current-year F, the
#' terminal-year seasonal F ratio, and fishery selectivity; (3) numbers-at-age
#' are redistributed among regions via the movement transition matrices
#' (recruits optionally excluded); (4) within-season survival is applied via
#' exponential decay, and at the end of the final season the population ages
#' by one year with plus-group accumulation; (5) spawning biomass is
#' calculated in `spawn_seas` with a mid-season survival correction; and
#' (6) catch is derived from the Baranov equation. Fishing mortality for the
#' following year is then set by the HCR or taken directly from `f_ref_pt`.
#'
#' Single-sex models receive a 0.5 multiplier applied to SSB. When
#' `n_regions = 1`, movement steps are skipped. When `n_seas = 1` and
#' `n_pop > 1`, `sgl_seas_spawning_movement` redistributes fish to natal
#' spawning grounds before SSB is calculated.
#'
#' @export Do_Population_Projection
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
                                     natmort,
                                     natal_region,
                                     WAA,
                                     WAA_fish,
                                     MatAA,
                                     fish_sel,
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
                                     init_F = rep(0, n_seas)
                                     ) {


# Error Checking ----------------------------------------------------------

  if(!recruitment_opt %in% c("inv_gauss", "mean_rec", "zero", "bh_rec")) stop("Recruitment options are not specified correctly! Should be inv_gauss, mean_rec, zero, or bh_rec")
  if(!fmort_opt %in% c("HCR", "Input", "HCR_global")) stop("Fishing Mortality options are not specified correctly! Should be HCR, Input, HCR_global")
  if(recruitment_opt == "bh_rec") {
    required_fields <- c("rec_dd", "rec_lag", "R0", "h", "rec_region_prop",
                         "WAA", "MatAA", "natmort", "SSB", "Movement",
                         "sex_ratio_f", "stray_rate")
    diff <- setdiff(required_fields, names(bh_rec_opt)) # find difference
    if(length(diff) > 0) stop(paste("bh_rec_opt is missing the following required fields:", paste(diff)))
  }

# Define Containers -------------------------------------------------------
  fratio <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  for(r in 1:n_regions) for(seas in 1:n_seas) for(f in 1:n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
  proj_NAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_NAA0 <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_ZAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_FAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_CAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_Catch <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets))
  proj_SSB <- array(0, dim = c(n_pop, n_regions, n_proj_yrs))
  proj_eff_SSB <- array(0, dim = c(n_pop, n_proj_yrs))
  proj_Dynamic_SSB0 <- array(0, dim = c(n_pop, n_regions, n_proj_yrs))
  proj_F <- array(0, dim = c(n_regions, n_proj_yrs + 1))

# Start Projection --------------------------------------------------------
  # Input terminal year assessment at age
  proj_NAA[,,1,,,] <- terminal_NAA
  proj_NAA0[,,1,,,] <- terminal_NAA0

  for(y in 1:n_proj_yrs) {

    # use terminal F in the first year (subsequent years use F derived from reference points and HCR)
    if(y == 1) proj_F[,y] <- rowSums(terminal_F)

    # Recruitment Processes ---------------------------------------------------
    if(y > 1) {

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
                                              stray_rate = bh_rec_opt$stray_rate,
                                              do_recruits_move = do_recruits_move,
                                              t_spawn = t_spawn,
                                              sexratio_f = bh_rec_opt$sex_ratio_f,
                                              init_F = init_F,
                                              fish_sel = array(fish_sel[,1,,1,1], dim = c(n_regions, n_ages))
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

      # Insert seasonal recruits at seas > 1
      if(seas > 1 && y > 1) {
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
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          for(f in 1:n_fish_fleets) {
            # get fishing mortality at age
            proj_FAA[,y,seas,a,s,f] <- proj_F[,y] * fratio[,seas,f] * fish_sel[,y,a,s,f]
          } # end f loop

          # Get Total Mortality at Age
          for(p in 1:n_pop) {
            proj_ZAA[p,,y,seas,a,s] <- (natmort[p,,y,a,s] * seasdur[seas]) +
              apply(proj_FAA[,y,seas,a,s,,drop = FALSE], 1, sum) # M and sum F across fleets
          } # end p loop

        } # end s loop
      } # end a loop

      # Movement Processes ------------------------------------------------------
      # Only apply movement if more than 1 region, or if y > 1 (because terminal proj_NAA already has movement applied)
      if(n_regions > 1 && y > 1) {
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

      # Mortality and Ageing ----------------------------------------------------
      if(seas < n_seas && y > 1) { # within season mortality
        proj_NAA[,,y,seas+1,1:n_ages,] = proj_NAA[,,y,seas,1:n_ages,] * exp(-proj_ZAA[,,y,seas,1:n_ages,])
        proj_NAA0[,,y,seas+1,1:n_ages,] = proj_NAA0[,,y,seas,1:n_ages,] * exp(-natmort[,,y,1:n_ages,] * seasdur[seas])
      } else { # age advancement
        # age advancement and enter into first season of next year
        proj_NAA[,,y+1,1,2:n_ages,] = proj_NAA[,,y,n_seas,1:(n_ages-1),] * exp(-proj_ZAA[,,y,seas,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group
        proj_NAA[,,y+1,1,n_ages,] = proj_NAA[,,y+1,1,n_ages,] + proj_NAA[,,y,seas,n_ages,] * exp(-proj_ZAA[,,y,seas,n_ages,]) # Acuumulate plus group
        proj_NAA0[,,y+1,1,2:n_ages,] = proj_NAA0[,,y,n_seas,1:(n_ages-1),] * exp(-natmort[,,y,1:(n_ages-1),] * seasdur[n_seas]) # Exponential mortality for individuals not in plus group
        proj_NAA0[,,y+1,1,n_ages,] = proj_NAA0[,,y+1,1,n_ages,] + proj_NAA0[,,y,n_seas,n_ages,] * exp(-natmort[,,y,n_ages,] * seasdur[n_seas]) # Acuumulate plus group
      }

      # Derive Biomass ----------------------------------------------------------
      if(seas == spawn_seas) {

        # Get proj_NAA for spawning
        tmp_NAA_spawn = proj_NAA[,,y,spawn_seas,,, drop = FALSE]
        tmp_NAA0_spawn = proj_NAA0[,,y,spawn_seas,,, drop = FALSE]

        # If we we are natal homing with 1 season
        if(n_seas == 1 && n_pop > 1) {
          # Get proj_NAA during spawning
          for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
            tmp_NAA_spawn[p,,1,1,a,s] = tmp_NAA_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
            tmp_NAA0_spawn[p,,1,1,a,s] = tmp_NAA0_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
          } # end s loop
        }

        # get SSB
        proj_SSB[,, y] = apply(tmp_NAA_spawn[,, 1, 1, , 1,drop = FALSE] *
                               WAA[,, y, spawn_seas, , 1,drop = FALSE] *
                               MatAA[,, y, spawn_seas, , 1,drop = FALSE] *
                               exp(-proj_ZAA[,, y, spawn_seas, , 1,drop = FALSE] * t_spawn), c(1,2), sum)

        # Get dynamic B0
        SSB0_array = tmp_NAA0_spawn[,, 1, 1, , 1,drop = FALSE] *  WAA[,,  y, spawn_seas, , 1, drop = FALSE] * MatAA[,,y, spawn_seas, , 1, drop = FALSE]
        mort_spawn = exp(-natmort[,, y, , 1, drop = FALSE] * t_spawn * seasdur[spawn_seas])
        mort_spawn = array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
        proj_Dynamic_SSB0[,,y] = apply(SSB0_array * mort_spawn, c(1,2), sum) # Dynamic B0

        if(n_sexes == 1) { # If single sex model, multiply SSB calculations by 0.5
          proj_SSB[,,y] = proj_SSB[,,y] * 0.5
          proj_Dynamic_SSB0[,,y] = proj_Dynamic_SSB0[,,y] * 0.5
        }

        # Accumulate effective SSB at each population's natal region across all source populations to capture stray contributions
        if(n_pop > 1) {
          for(p2 in 1:n_pop) {
            for(p in 1:n_pop) {
              if(p == p2) {
                # Own population contribution - no stray scaling
                proj_eff_SSB[p2, y] = proj_eff_SSB[p2, y] + proj_SSB[p, natal_region[p2], y]
              } else {
                # Cross-population contribution scaled by stray_rate
                proj_eff_SSB[p2, y] = proj_eff_SSB[p2, y] + stray_rate[p, y] * proj_SSB[p, natal_region[p2], y]
              }
            } # end p loop
          } # end p2 loop
        } else proj_eff_SSB[1,y] = sum(proj_SSB[1,,y])

      } # calculate biomass

      # Derive Catches ----------------------------------------------------------
      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          for(f in 1:n_fish_fleets) {
            for(a in 1:n_ages) {
              for(s in 1:n_sexes) {
                # Get catch at age with Baranov's
                proj_CAA[p,r,y,seas,a,s,f] <- (proj_FAA[r,y,seas,a,s,f] / proj_ZAA[p,r,y,seas,a,s]) *
                  proj_NAA[p,r,y,seas,a,s] * (1 - exp(-proj_ZAA[p,r,y,seas,a,s]))
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
              proj_Catch = proj_Catch,
              proj_SSB = proj_SSB,
              proj_eff_SSB = proj_eff_SSB,
              proj_Dynamic_SSB0 = proj_Dynamic_SSB0,
              proj_NAA = proj_NAA,
              proj_NAA0 = proj_NAA0,
              proj_ZAA = proj_ZAA)
  )

} # end function

