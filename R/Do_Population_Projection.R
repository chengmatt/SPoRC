#' Do Population Projections
#'
#' Projects population dynamics forward in time under alternative recruitment
#' and fishing mortality scenarios. The projection initializes from terminal
#' assessment quantities and advances numbers-at-age by season and year,
#' applying fishing mortality, natural mortality, recruitment, movement,
#' ageing, and harvest control rules.
#'
#' @param n_proj_yrs Integer. Number of projection years.
#' @param n_pop Integer. Number of populations (may exceed regions when
#'   natal homing or spawning structure is modeled).
#' @param n_regions Integer. Number of spatial regions.
#' @param n_ages Integer. Number of age classes (including plus group).
#' @param n_sexes Integer. Number of sexes.
#' @param sexratio Recruitment sex ratio array used to allocate projected
#'   recruits by sex after recruitment is generated.
#'   Dimensioned [n_pop, n_regions, n_proj_yrs, n_sexes].
#' @param n_fish_fleets Integer. Number of fishing fleets.
#' @param do_recruits_move Integer (0 or 1). Whether age-1 recruits are
#'   subject to movement.
#' @param recruitment Historical recruitment matrix used to condition
#'   projections. Dimensioned [n_regions, n_yrs].
#' @param terminal_NAA Terminal fished numbers-at-age.
#'   Dimensioned [n_pop, n_regions, n_seas, n_ages, n_sexes].
#' @param terminal_NAA0 Terminal unfished numbers-at-age.
#'   Same dimensions as `terminal_NAA`.
#' @param terminal_F Terminal fishing mortality.
#'   Dimensioned [n_regions, n_seas, n_fish_fleets].
#' @param natmort Natural mortality-at-age.
#'   Dimensioned [n_pop, n_regions, n_proj_yrs, n_ages, n_sexes].
#' @param WAA Weight-at-age used in spawning biomass calculations.
#'   Dimensioned [n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes].
#' @param WAA_fish Fishery weight-at-age.
#'   Dimensioned [n_pop, n_regions, n_proj_yrs, n_seas,
#'   n_ages, n_sexes, n_fish_fleets].
#' @param MatAA Maturity-at-age.
#'   Dimensioned [n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes].
#' @param fish_sel Fishery selectivity-at-age.
#'   Dimensioned [n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets].
#' @param Movement Seasonal movement transition matrices.
#'   Dimensioned [n_pop, n_regions, n_regions,
#'   n_proj_yrs, n_seas, n_ages, n_sexes].
#' @param f_ref_pt Fishing mortality reference point or projected F input.
#'   Dimensioned [n_regions, n_proj_yrs].
#' @param b_ref_pt Biomass reference point.
#'   For `"HCR"`: dimensioned [n_pop, n_regions, n_proj_yrs].
#'   For `"HCR_global"`: region-specific values are summed internally.
#' @param HCR_function Harvest control rule function. Must accept:
#'   \describe{
#'     \item{x}{Spawning stock biomass (SSB)}
#'     \item{frp}{Fishing mortality reference point}
#'     \item{brp}{Biomass reference point}
#'   }
#' @param recruitment_opt Recruitment option. One of:
#' \describe{
#'   \item{"inv_gauss"}{Draw recruitment from inverse Gaussian
#'   distribution conditioned on supplied values.}
#'   \item{"mean_rec"}{Use mean historical recruitment by region.}
#'   \item{"zero"}{Set recruitment to zero.}
#'   \item{"bh_rec"}{Apply deterministic Beverton–Holt recruitment.}
#' }
#' @param fmort_opt Fishing mortality option. One of:
#' \describe{
#'   \item{"HCR"}{Apply region-specific harvest control rule.}
#'   \item{"HCR_global"}{Apply harvest control rule using global SSB
#'   and summed biomass reference points.}
#'   \item{"Input"}{Use values in `f_ref_pt` directly.}
#' }
#' @param t_spawn Fraction of the spawning season elapsed prior to spawning,
#'   used in SSB calculations.
#' @param bh_rec_opt List required when `recruitment_opt = "bh_rec"`.
#'   Must contain:
#' \describe{
#'   \item{rec_dd}{Density dependence indicator (1 = global, 0 = local).}
#'   \item{rec_lag}{Recruitment lag in years.}
#'   \item{R0}{Unfished recruitment.}
#'   \item{h}{Steepness (vector by region).}
#'   \item{Rec_Prop}{Recruitment apportionment by region.}
#'   \item{WAA}{Reference WAA array for recruitment calculation. (should index the first year)}
#'   \item{MatAA}{Reference maturity array. (should index the first year)}
#'   \item{natmort}{Reference natural mortality array. (should use the first year)}
#'   \item{SSB}{Historical SSB array [n_pop, n_regions, n_yrs].}
#'   \item{Movement}{Movement array for recruitment distribution. (should index the first year)}
#'   \item{do_recruits_move}{Recruit movement indicator.}
#'   \item{t_spawn}{Spawning timing fraction.}
#'   \item{sex_ratio_f}{Female recruitment proportion used internally in the
#'     Beverton–Holt recruitment calculation. (should index the first year)}
#'   \item{sgl_seas_spawning_movement}{Single-season spawning movement array,
#'   if applicable. (should index the first year)}
#'   \item{stray_rate}{Vector of stray rates by population. (should index the first year)}
#' }
#' @param n_seas Integer. Number of seasons. Default = 1.
#' @param seasdur Vector of seasonal durations (length = n_seas).
#'   Default = rep(1 / n_seas, n_seas).
#' @param spawn_seas Integer. Spawning season index. Default = 1.
#' @param init_F Vector of initial seasonal F values (length = n_seas)
#'   used when deriving Beverton–Holt recruitment.
#' @param sgl_seas_spawning_movement Single-season spawning movement matrix
#'   used when `n_seas = 1` and `n_pop > 1`. Dimensioned [n_pop, n_regions, n_regions, n_proj_yrs, n_ages, n_sexes].
#'
#' @param stray_rate Array of stray rates dimensioned by [n_pop, n_proj_yrs]
#' @returns A list containing:
#' \describe{
#'   \item{proj_F}{Annual fishing mortality
#'   [n_regions, n_proj_yrs + 1].}
#'   \item{proj_Catch}{Seasonal catch biomass
#'   [n_pop, n_regions, n_proj_yrs, n_seas, n_fish_fleets].}
#'   \item{proj_SSB}{Projected spawning stock biomass
#'   [n_pop, n_regions, n_proj_yrs].}
#'   \item{proj_SSB}{Projected effective spawning stock biomass (by population)
#'   [n_pop, n_proj_yrs].}
#'   \item{proj_Dynamic_SSB0}{Projected dynamic unfished SSB
#'   [n_pop, n_regions, n_proj_yrs].}
#'   \item{proj_NAA}{Projected fished numbers-at-age
#'   [n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes].}
#'   \item{proj_NAA0}{Projected unfished numbers-at-age
#'   (same dimensions as `proj_NAA`).}
#'   \item{proj_ZAA}{Total mortality-at-age
#'   [n_pop, n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes].}
#' }
#'
#' @details
#' Population dynamics are projected forward by season and year using
#' exponential survival, seasonal movement, recruitment, and ageing.
#'
#' Let:
#' \itemize{
#'   \item \eqn{N_{p,r,a,s,y}} = numbers-at-age
#'   \item \eqn{F_{r,a,s,y}} = fishing mortality-at-age
#'   \item \eqn{M_{p,r,a,y}} = natural mortality
#'   \item \eqn{Z = M + \sum_f F_f} = total mortality
#' }
#'
#' **Within-season survival**
#'
#' \deqn{
#' N_{p,r,a,s+1,y} =
#' N_{p,r,a,s,y} \exp(-Z_{p,r,a,s,y})
#' }
#'
#' **Age advancement (end of final season)**
#'
#' \deqn{
#' N_{p,r,a+1,1,y+1} =
#' N_{p,r,a,n_{seas},y}
#' \exp(-Z_{p,r,a,n_{seas},y})
#' }
#'
#' The plus group accumulates survivors:
#'
#' \deqn{
#' N_{A^+,y+1} =
#' N_{A-1,y}e^{-Z_{A-1,y}} +
#' N_{A^+,y}e^{-Z_{A^+,y}}
#' }
#'
#' **Spawning stock biomass**
#'
#' Spawning biomass is calculated in `spawn_seas` as:
#'
#' \deqn{
#' SSB_{p,r,y} =
#' \sum_a N_{p,r,a,spawn,y}
#' W_{p,r,a,y}
#' Mat_{p,r,a,y}
#' \exp(-Z_{p,r,a,spawn,y} t_{spawn})
#' }
#'
#' where \eqn{t_{spawn}} is the fraction of the season elapsed before spawning.
#'
#' **Catch (Baranov equation)**
#'
#' \deqn{
#' C_{p,r,a,s,f,y} =
#' \frac{F_{r,a,s,f,y}}{Z_{p,r,a,s,y}}
#' N_{p,r,a,s,y}
#' \left(1 - e^{-Z_{p,r,a,s,y}}\right)
#' }
#'
#' **Harvest control rule**
#'
#' Fishing mortality in year \eqn{y+1} is derived as:
#'
#' \deqn{
#' F_{r,y+1} =
#' HCR(SSB, F_{ref}, B_{ref})
#' }
#'
#' either regionally (`"HCR"`) or using global SSB (`"HCR_global"`).
#'
#' Recruitment is inserted at age 1 in season 1 according to
#' `recruitment_opt` and, when `"bh_rec"` is selected, follows the
#' Beverton–Holt formulation described in `Get_Det_Recruitment()`.
#'
#' @export Do_Population_Projection
#' @family Reference Points and Projections
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
                                     init_F = rep(0, n_seas)
                                     ) {


# Error Checking ----------------------------------------------------------

  if(!recruitment_opt %in% c("inv_gauss", "mean_rec", "zero", "bh_rec")) stop("Recruitment options are not specified correctly! Should be inv_gauss, mean_rec, zero, or bh_rec")
  if(!fmort_opt %in% c("HCR", "Input", "HCR_global")) stop("Fishing Mortality options are not specified correctly! Should be HCR, Input, HCR_global")
  if(recruitment_opt == "bh_rec") {
    required_fields <- c("rec_dd", "rec_lag", "R0", "h", "Rec_Prop",
                         "WAA", "MatAA", "natmort", "SSB", "Movement", "do_recruits_move", "t_spawn",
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

    for(seas in 1:n_seas) {

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
              apply(proj_FAA[,y,seas,a,s,,drop = FALSE], c(1:5), sum) # M and sum F across fleets
          } # end p loop

        } # end s loop
      } # end a loop

      # Recruitment Processes ---------------------------------------------------
      if(y > 1 && seas == 1) {

        # Get recruitment
        tmp_rec <- switch(recruitment_opt,

                          "inv_gauss" = { # if inverse gaussian
                            sapply(1:n_regions, function(r) rinvgauss_rec(1, recruitment[r,]))
                          },

                          "mean_rec" = { # if mean recruitment
                            sapply(1:n_regions, function(r) mean(recruitment[r,]))
                          },

                          "zero" = { # if zero recruitment
                            rep(0, n_regions)
                          },

                          "bh_rec" = { # if beverton holt recruitment
                            Get_Det_Recruitment(recruitment_model = 1,
                                                rec_dd = bh_rec_opt$rec_dd,
                                                n_pop = n_pop,
                                                sgl_seas_spawning_movement = bh_rec_opt$sgl_seas_spawning_movement,
                                                natal_region = natal_region,
                                                y = y + dim(bh_rec_opt$SSB)[2],
                                                rec_lag = bh_rec_opt$rec_lag,
                                                R0 = bh_rec_opt$R0,
                                                Rec_Prop = bh_rec_opt$Rec_Prop,
                                                h = bh_rec_opt$h,
                                                n_regions = n_regions,
                                                n_ages = n_ages,
                                                WAA = bh_rec_opt$WAA,
                                                MatAA = bh_rec_opt$MatAA,
                                                n_seas = n_seas,
                                                seasdur = seasdur,
                                                spawn_seas = spawn_seas,
                                                natmort = bh_rec_opt$natmort,
                                                SSB_vals = cbind(bh_rec_opt$SSB, proj_SSB),
                                                Movement = bh_rec_opt$Movement,
                                                stray_rate = bh_rec_opt$stray_rate,
                                                do_recruits_move = bh_rec_opt$do_recruits_move,
                                                t_spawn = bh_rec_opt$t_spawn,
                                                sex_ratio_f = bh_rec_opt$sex_ratio_f,
                                                init_F = init_F,
                                                fish_sel = array(fish_sel[,1,,1,1], dim = c(n_regions, n_ages))
                            )
                          }
        )

        # Apply recruitment to projected NAA
        for(p in 1:n_pop) {
          for(r in 1:n_regions) {
            if(n_sexes == 2) tmp <- tmp_rec[p,r] * sexratio[p,r,y,]
            if(n_sexes == 1) tmp <- tmp_rec[p,r]
            proj_NAA[p,r,y,1,1,] <- proj_NAA0[p,r,y,1,1,]  <- tmp
          } # end r loop
        } # end p loop

      }

      # Movement Processes ------------------------------------------------------
      # Only apply movement if more than 1 region, or if y > 1 (because terminal NAA already has movement applied)
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
      if(seas < n_seas) { # within season mortality
        proj_NAA[,,y,seas+1,1:n_ages,] = proj_NAA[,,y,seas,1:n_ages,] * exp(-proj_ZAA[,,y,seas,1:n_ages,])
        proj_NAA0[,,y,seas+1,1:n_ages,] = proj_NAA0[,,y,seas,1:n_ages,] * exp(-natmort[,,y,1:n_ages,] * seasdur[seas])
      } else { # age advancement
        # Fished
        proj_NAA[,,y+1,1,2:n_ages,] = proj_NAA[,,y,n_seas,1:(n_ages-1),] * exp(-proj_ZAA[,,y,seas,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group
        proj_NAA[,,y+1,1,n_ages,] = proj_NAA[,,y+1,1,n_ages,] + proj_NAA[,,y,seas,n_ages,] * exp(-proj_ZAA[,,y,seas,n_ages,]) # Acuumulate plus group
        # Unfished
        proj_NAA0[,,y+1,1,2:n_ages,] = proj_NAA0[,,y,n_seas,1:(n_ages-1),] * exp(-natmort[,,y,1:(n_ages-1),] * seasdur[n_seas]) # Exponential mortality for individuals not in plus group
        proj_NAA0[,,y+1,1,n_ages,] = proj_NAA0[,,y+1,1,n_ages,] + proj_NAA0[,,y,n_seas,n_ages,] * exp(-natmort[,,y,n_ages,] * seasdur[n_seas]) # Acuumulate plus group
      }

      # Derive Biomass ----------------------------------------------------------
      if(seas == spawn_seas) {

        # Get NAA for spawning
        tmp_NAA_spawn = proj_NAA[,,y,spawn_seas,,, drop = FALSE]
        tmp_NAA0_spawn = proj_NAA0[,,y,spawn_seas,,, drop = FALSE]

        # If we we are natal homing with 1 season
        if(n_seas == 1 && n_pop > 1) {
          # Get NAA during spawning
          for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
            tmp_NAA_spawn[p,,1,1,a,s] = tmp_NAA_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
            tmp_NAA0_spawn[p,,1,1,a,s] = tmp_NAA0_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
          } # end s loop
        }

        # get SSB
        proj_SSB[,, y] = apply(tmp_NAA_spawn *
                            WAA[,, y, spawn_seas, , 1,drop = FALSE] *
                            MatAA[,, y, spawn_seas, , 1,drop = FALSE] *
                            exp(-proj_ZAA[,, y, spawn_seas, , 1,drop = FALSE] * t_spawn), c(1,2), sum)

        # Get dynamic B0
        SSB0_array = tmp_NAA0_spawn *  WAA[,,  y, spawn_seas, , 1, drop = FALSE] * MatAA[,,y, spawn_seas, , 1, drop = FALSE]
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
        } else proj_eff_SSB[1, y] = sum(proj_SSB[1,,y])

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
                                      brp = b_ref_pt[,r,y])
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

