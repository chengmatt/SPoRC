#' Do Population Projections
#'
#' @param n_proj_yrs Number of projection years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param sexratio Array of recruitment sex ratio (n_regions, n_proj_yrs, n_sexes)
#' @param n_fish_fleets Number of fishery fleets
#' @param do_recruits_move Whether recruits move (0 == don't move, 1 == move)
#' @param recruitment Recruitment matrix dimensioned by n_regions, and n_yrs that we want to summarize across, or condition our projection on
#' @param terminal_NAA Terminal Numbers at Age dimensioned by n_regions, n_seas, n_ages, n_sexes
#' @param terminal_NAA0 Terminal Unfished Numbers at Age dimensioned by n_regions, n_seas, n_ages, n_sexes
#' @param terminal_F Terminal fishing mortality rate, dimensioned by n_regions, n_seas, n_fish_fleets
#' @param natmort Natural mortality, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes
#' @param WAA Weight at age, dimensioned by n_regions, n_proj_yrs, n_seas, n_ages, n_sexes
#' @param WAA_fish Weight at age for the fishery, dimensioned by n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets
#' @param MatAA Maturity at age, dimensioned by n_regions, n_proj_yrs, n_seas, n_ages, n_sexes
#' @param fish_sel Fishery selectivity, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets
#' @param Movement Movement, dimensioned by n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes
#' @param f_ref_pt Fishing mortality reference point dimensioned by n_regions and n_proj_yrs
#' @param b_ref_pt Biological reference point dimensioned by n_regions and n_proj_yrs
#' @param HCR_function Function describing a harvest control rule. The function should always have the following arguments: x, which represents SSB, frp, which takes inputs of fishery reference points, and brp, which takes inputs of biological reference points. Any additional arguments should be specified with defaults or hard coded / fixed within the function.
#' @param recruitment_opt Recruitment simulation option, where options are "inv_gauss", which simulates future recruitment based on the the recruitment values supplied using an inverse gaussian distribution, "mean_rec", which takes the mean of the recruitment values supplied for a given region, and "zero", which assumes that future recruitment does not occur
#' @param fmort_opt Fishing mortality option. Choices are:
#'
#' * **"HCR"** – Applies the user-supplied `HCR_function` using region-specific
#'   SSB, F reference point, and biomass reference point.
#'
#' * **"HCR_global"** – Applies the `HCR_function` using global SSB (summed
#'   across regions) and a global biomass reference point (sum of the
#'   region-specific biomass reference points). Each region's biomass reference
#'   point should be defined individually; the function performs the summation.
#'
#' * **"Input"** – Uses user-supplied projected fishing mortality values directly.
#' @param t_spawn Fraction time of spawning used to compute projected SSB
#' @param bh_rec_opt A list object containing the following arguments:
#' \describe{
#' \item{recruitment_dd}{A value (0 or 1) indicating global (1) or local density dependence (0). In the case of a single region model, either local or global will give the same results}
#' \item{rec_lag}{A value indicating the number of years lagged that a given year's SSB produces recruits}
#' \item{R0}{The virgin recruitment parameter}
#' \item{Rec_Prop}{Recruitment apportionment values. In a single region model, this should be set at a value of 1. Dimensioned by n_regions}
#' \item{h}{Steepness values for the stock recruitment curve. Dimensioned by n_regions}
#' \item{WAA}{A weight-at-age array dimensioned by n_regions, n_seas, n_ages, and n_sexes, where the reference year should utilize values from the first year}
#' \item{MatAA}{A maturity at age array dimensioned by n_regions, n_seas, n_ages, and n_sexes, where the reference year should utilize values from the first year}
#' \item{natmort}{A natural mortality at age array dimensioned by n_regions, n_ages, and n_sexes, where the reference year should utilize values from the first year}
#' \item{SSB}{All SSB values estimated from a given model, dimensioned by n_regions and n_yrs}
#' \item{spawn_seas}{Spawning season index}
#' }
#' @param init_F Vector of initial F values (n_seas) to apply for deriving beverton holt recruitment; default is set at 0.
#' @param n_seas Number of seasons. Default = 1
#' @param seasdur Vector of season durations (length n_seas). Default = 1 / n_seas
#' @param spawn_seas Spawning season index. Default = 1
#'
#' @returns A list containing projected F, catch, SSB (and dynamic unfished), and Numbers at Age (and dynamic unfished). (Objects are generally dimensioned in the following order: n_regions, n_yrs, n_ages, n_sexes, n_fleets)
#' @export Do_Population_Projection
#' @family Reference Points and Projections
Do_Population_Projection <- function(n_proj_yrs = 2,
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
    required_fields <- c("recruitment_dd", "rec_lag", "R0", "h", "Rec_Prop",
                         "WAA", "MatAA", "natmort", "SSB", "Movement", "do_recruits_move", "t_spawn",
                         "sex_ratio_f")
    diff <- setdiff(required_fields, names(bh_rec_opt)) # find difference
    if(length(diff) > 0) stop(paste("bh_rec_opt is missing the following required fields:", paste(diff)))
  }

# Define Containers -------------------------------------------------------
  fratio <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  for(r in 1:n_regions) for(seas in 1:n_seas) for(f in 1:n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
  proj_NAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_NAA0 <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_ZAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes))
  proj_FAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_CAA <- array(0, dim = c(n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  proj_Catch <- array(0, dim = c(n_regions, n_proj_yrs, n_seas, n_fish_fleets))
  proj_SSB <- array(0, dim = c(n_regions, n_proj_yrs))
  proj_Dynamic_SSB0 <- array(0, dim = c(n_regions, n_proj_yrs))
  proj_F <- array(0, dim = c(n_regions, n_proj_yrs + 1))

# Start Projection --------------------------------------------------------
  # Input terminal year assessment at age
  proj_NAA[,1,,,] <- terminal_NAA
  proj_NAA0[,1,,,] <- terminal_NAA0

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
          proj_ZAA[,y,seas,a,s] <- (natmort[,y,a,s] * seasdur[seas]) +
                                    apply(proj_FAA[,y,seas,a,s,,drop = FALSE], c(1:5), sum) # M and sum F across fleets

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
                                                recruitment_dd = bh_rec_opt$recruitment_dd,
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
                                                do_recruits_move = bh_rec_opt$do_recruits_move,
                                                t_spawn = bh_rec_opt$t_spawn,
                                                sex_ratio_f = bh_rec_opt$sex_ratio_f,
                                                init_F = init_F,
                                                fish_sel = array(fish_sel[,1,,1,1], dim = c(n_regions, n_ages))
                            )
                          }
        )

        # Apply recruitment to projected NAA
        for(r in 1:n_regions) {
          if(n_sexes == 2) tmp <- tmp_rec[r] * sexratio[r,y,]
          if(n_sexes == 1) tmp <- tmp_rec[r]
          proj_NAA[r,y,1,1,] <- proj_NAA0[r,y,1,1,]  <- tmp
        } # end r loop

      }

      # Movement Processes ------------------------------------------------------
      # Only apply movement if more than 1 region, or if y > 1 (because terminal NAA already has movement applied)
      if(n_regions > 1 && y > 1) {
        # Recruits don't move
        if(do_recruits_move == 0) {
          # Apply movement after ageing processes - start movement at age 2
          for(a in 2:n_ages) for(s in 1:n_sexes) proj_NAA[,y,seas,a,s] = t(proj_NAA[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # fished
          for(a in 2:n_ages) for(s in 1:n_sexes) proj_NAA0[,y,seas,a,s] = t(proj_NAA0[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # unfished
        } # end if recruits don't move
        # Recruits move here
        if(do_recruits_move == 1) {
          for(a in 1:n_ages) for(s in 1:n_sexes) proj_NAA[,y,seas,a,s] = t(proj_NAA[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # fished
          for(a in 1:n_ages) for(s in 1:n_sexes) proj_NAA0[,y,seas,a,s] = t(proj_NAA0[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # unfished
        }
      } # only compute if spatial

      # Mortality and Ageing ----------------------------------------------------
      if(seas < n_seas) { # within season mortality
        proj_NAA[,y,seas+1,1:n_ages,] = proj_NAA[,y,seas,1:n_ages,] * exp(-proj_ZAA[,y,seas,1:n_ages,])
        proj_NAA0[,y,seas+1,1:n_ages,] = proj_NAA0[,y,seas,1:n_ages,] * exp(-natmort[,y,1:n_ages,] * seasdur[seas])
      } else { # age advancement
        # Fished
        proj_NAA[,y+1,1,2:n_ages,] = proj_NAA[,y,n_seas,1:(n_ages-1),] * exp(-proj_ZAA[,y,seas,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group
        proj_NAA[,y+1,1,n_ages,] = proj_NAA[,y+1,1,n_ages,] + proj_NAA[,y,seas,n_ages,] * exp(-proj_ZAA[,y,seas,n_ages,]) # Acuumulate plus group
        # Unfished
        proj_NAA0[,y+1,1,2:n_ages,] = proj_NAA0[,y,n_seas,1:(n_ages-1),] * exp(-natmort[,y,1:(n_ages-1),] * seasdur[n_seas]) # Exponential mortality for individuals not in plus group
        proj_NAA0[,y+1,1,n_ages,] = proj_NAA0[,y+1,1,n_ages,] + proj_NAA0[,y,n_seas,n_ages,] * exp(-natmort[,y,n_ages,] * seasdur[n_seas]) # Acuumulate plus group
      }

      # Derive Biomass ----------------------------------------------------------
      if(seas == spawn_seas) {

        # Get SSB
        proj_SSB[,y] = apply(proj_NAA[,y,spawn_seas,,1,drop = FALSE] * exp(-proj_ZAA[,y,spawn_seas,,1,drop = FALSE] * t_spawn) * WAA[,y,spawn_seas,,1,drop = FALSE] * MatAA[,y,spawn_seas,,1,drop = FALSE], 1, sum) # Spawning Stock Biomass

        # Get dynamic B0
        SSB0_array <- proj_NAA0[,y,spawn_seas,,1,drop = FALSE] *  WAA[,y,spawn_seas,,1,drop = FALSE] * MatAA[,y,spawn_seas,,1,drop = FALSE]
        mort_spawn <- exp(-natmort[,y,,1,drop = FALSE] * t_spawn * seasdur[seas])
        mort_spawn <- array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
        proj_Dynamic_SSB0[,y] <- apply(SSB0_array * mort_spawn, 1, sum) # Dynamic B0

        if(n_sexes == 1) {  # If single sex model, multiply SSB calculations by 0.5
          proj_SSB[,y] = proj_SSB[,y] * 0.5
          proj_Dynamic_SSB0[,y] = proj_Dynamic_SSB0[,y] * 0.5
        }

      } # calculate biomass

      # Derive Catches ----------------------------------------------------------
      for(r in 1:n_regions) {
        for(f in 1:n_fish_fleets) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              # Get catch at age with Baranov's
              proj_CAA[r,y,seas,a,s,f] <- (proj_FAA[r,y,seas,a,s,f] / proj_ZAA[r,y,seas,a,s]) *
                                           proj_NAA[r,y,seas,a,s] * (1 - exp(-proj_ZAA[r,y,seas,a,s]))
            } # end s loop
          } # end a loop

          # Get total catch
          proj_Catch[r,y,seas,f] <- sum(proj_CAA[r,y,seas,,,f] * WAA_fish[r,y,seas,,,f])

        } # end f loop
      } # end r loop

    } # end seas loop

    # compute F for next year
    for(r in 1:n_regions) {
      # Project F using HCR and reference points -----------------------------------------------------
      if(fmort_opt == 'HCR') {
        proj_F[r,y+1] <- HCR_function(x = proj_SSB[r,y],
                                      frp = f_ref_pt[r,y],
                                      brp = b_ref_pt[r,y])
      }

      if(fmort_opt == 'HCR_global') {
        proj_F[r,y+1] <- HCR_function(x = sum(proj_SSB[,y]),
                                      frp = f_ref_pt[r,y],
                                      brp = sum(b_ref_pt[,y]))
      }

      # Project F using User Inputs ---------------------------------------------
      if(fmort_opt == 'Input') proj_F[r,y+1] <- f_ref_pt[r,y]

    } # end r loop

  } # end y loop

  return(list(proj_F = proj_F,
              proj_Catch = proj_Catch,
              proj_SSB = proj_SSB,
              proj_Dynamic_SSB0 = proj_Dynamic_SSB0,
              proj_NAA = proj_NAA,
              proj_NAA0 = proj_NAA0,
              proj_ZAA = proj_ZAA)
  )

} # end function

