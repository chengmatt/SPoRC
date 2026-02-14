#' Initialize Numbers-at-Age (NAA) for a Population Model
#'
#' This function generates initial numbers-at-age (NAA) for a structured population
#' model across regions, sexes, and age classes. It supports multiple initialization
#' methods, including iterative solution, scalar geometric series, and matrix
#' geometric series, optionally accounting for movement and fishing mortality. Initial
#' age deviations can also be applied.
#'
#' @param init_age_strc Integer specifying the initialization method for the age structure:
#'   - 0: Iterative solution to equilibrium
#'   - 1: Scalar geometric series solution w/o movement in any groups (no movement in all groups)
#'   - 2: Matrix geometric series solution (generalizes scalar solution with movement)
#'   - 3: Scalar geometric series solution w/o movement only in plus group (no movement in plus groups)
#' @param init_iter Integer; number of iterations to run when `init_age_strc = 0`.
#' @param n_regions Integer; number of spatial regions.
#' @param n_sexes Integer; number of sexes.
#' @param n_ages Integer; number of age classes.
#' @param natmort Array of natural mortality rates with dimensions `[regions, ages, sexes]`.
#' @param init_F Vector; initial fishing mortality applied by season (0 for unfished population).
#' @param fish_sel Array of fishery selectivity with dimensions `[regions, ages, sexes, fleets]`.
#' @param R0_r Numeric vector of recruitment values for each region.
#' @param sexratio Array `[regions, sexes]` giving the proportion of each sex in recruitment.
#' @param Movement Array `[regions, regions, season, ages, sexes]` defining movement probabilities.
#' @param do_recruits_move Integer; 0 = recruits do not move, 1 = recruits move according to `Movement`.
#' @param ln_InitDevs Array `[regions, ages-1]` of log-scale deviations for initial numbers-at-age.
#' @param n_seas Number of seasons
#' @param seasdur Fraction of time within each season
#'
#' @return Array of initial numbers-at-age with dimensions `[regions, ages, sexes]`.
#'
#' @keywords internal
Get_Init_NAA <- function(init_age_strc,
                         init_iter,
                         n_regions,
                         n_sexes,
                         n_ages,
                         n_seas,
                         seasdur,
                         natmort,
                         init_F,
                         fish_sel,
                         R0_r,
                         sexratio,
                         Movement,
                         do_recruits_move,
                         ln_InitDevs
                         ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # create containers
  Init_NAA = array(0, dim = c(n_regions, n_ages, n_sexes))
  Init_NAA_next_year = array(0, dim = c(n_regions, n_ages, n_sexes))
  NAA = array(0, dim = c(n_regions, n_ages, n_sexes))

  # Iterative Solution
  if(init_age_strc == 0) {

    # initialize age structure (starting point)
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        tmp_cumsum_Z = cumsum(natmort[r,1:(n_ages-1),s] + sum(init_F) * fish_sel[r,1:(n_ages-1),s,1])
        Init_NAA[r,,s] = c(R0_r[r] * sexratio[r,s], R0_r[r] * sexratio[r,s] * exp(-tmp_cumsum_Z))
      } # end s loop
    } # end r loop

    # Apply annual cycle and iterate to equilibrium
    for(i in 1:init_iter) {
      for(s in 1:n_sexes) {
        for(seas in 1:n_seas) {

          # recruitment happens in the first season
          if(seas == 1) Init_NAA[,1,s] = R0_r * sexratio[,s]

          # movement
          if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,seas,a,s] # recruits don't move
          if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,seas,a,s] # recruits move

          # Apply mortality
          for(r in 1:n_regions) {
            # mortality wtihin season
            if(seas < n_seas) {
              Init_NAA_next_year[r,1:n_ages,s] = Init_NAA[r,1:n_ages,s] *
                exp(-((natmort[r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
            } else {
              # ageing and mortality (advance ages in the next year)
              Init_NAA_next_year[r,2:n_ages,s] = Init_NAA[r,1:(n_ages-1),s] *
                exp(-((natmort[r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
              # accumulate plus group
              Init_NAA_next_year[r,n_ages,s] = (Init_NAA_next_year[r,n_ages,s]) + (Init_NAA[r,n_ages,s] *
                                                exp(-((natmort[r,n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,n_ages,s,1]))))
            }
          } # end r loop
          Init_NAA = Init_NAA_next_year # iterate to next cycle
        } # end seas loop

      } # end s loop
    } # end i loop
    # save result
    NAA[,,] = Init_NAA
  } # end if iterative solution

  # Scalar Geometric Series Solution (no movement in all ages)
  if(init_age_strc == 1) {
    # projection initial abundance forward
    for(i in 1:n_ages) {
      for(s in 1:n_sexes) {
        for(seas in 1:n_seas) {

          if(seas == 1) Init_NAA[,1,s] = R0_r * sexratio[,s] # initialize recruitment

          for(r in 1:n_regions) {
            # within season mortality
            if(seas < n_seas) {
              Init_NAA[r,1:n_ages,s] = Init_NAA[r,1:n_ages,s] *
                exp(-((natmort[r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
            } else {
              tmp_plus_befage = Init_NAA[r,n_ages,s] # save temporary plus group before ageing
              # ageing and mortality (age advancement)
              Init_NAA[r,2:n_ages,s] = Init_NAA[r,1:(n_ages-1),s] * exp(-((natmort[r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
              # accumulate plus group
              Init_NAA[r,n_ages,s] = (Init_NAA[r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[r,n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,n_ages,s,1]))))
            }

          } # end r loop
        } # end seas loop
      } # end s loop
    } # end i loop

    # Set up analytical solution for plus group
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        # Plus group - scalar geometric series (summing init_F to get annual Z)
        Z_penult = natmort[r,n_ages-1,s] + (sum(init_F) * fish_sel[r,n_ages-1,s,1])
        Z_plus = natmort[r,n_ages,s] + (sum(init_F) * fish_sel[r,n_ages,s,1])
        Init_NAA[r,n_ages,s] = Init_NAA[r,n_ages-1,s] * exp(-Z_penult) / (1 - exp(-Z_plus))
      } # end s loop
    } # end r loop

    # save result
    NAA = Init_NAA
  } # end if

  # Matrix Geometric Series Solution (genearlizes to scalar w/o movement)
  if(init_age_strc == 2) {

    # projection initial abundance forward
    for(i in 1:n_ages) {
      for(s in 1:n_sexes) {
        for(seas in 1:n_seas) {

          if(seas == 1) Init_NAA[,1,s] = R0_r * sexratio[,s] # initialize recruitment

          # movement
          if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,seas,a,s] # recruits don't move
          if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,seas,a,s] # recruits move

          for(r in 1:n_regions) {
            # within season mortality
            if(seas < n_seas) {
              Init_NAA[r,1:n_ages,s] = Init_NAA[r,1:n_ages,s] *
                exp(-((natmort[r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
            } else {
              tmp_plus_befage = Init_NAA[r,n_ages,s] # save temporary plus group before ageing
              # ageing and mortality (age advancement)
              Init_NAA[r,2:n_ages,s] = Init_NAA[r,1:(n_ages-1),s] * exp(-((natmort[r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
              # accumulate plus group
              Init_NAA[r,n_ages,s] = (Init_NAA[r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[r,n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,n_ages,s,1]))))
            }

          } # end r loop
        } # end seas loop
      } # end s loop
    } # end i loop

    # Set up analytical solution for plus group
    for(s in 1:n_sexes) {

      # build annual transition for penultimate and plus ages
      T_penult = diag(n_regions)
      T_plus = diag(n_regions)

      for(seas in 1:n_seas) {
        S_penult = diag(exp(-((natmort[,n_ages-1,s] * seasdur[seas]) + init_F[seas] * fish_sel[,n_ages-1,s,1])), n_regions)
        S_plus = diag(exp(-((natmort[,n_ages,s] * seasdur[seas]) + init_F[seas] * fish_sel[,n_ages,s,1])), n_regions)
        T_penult = T_penult %*% t(Movement[,,seas,n_ages-1,s]) %*% S_penult
        T_plus = T_plus %*% t(Movement[,,seas,n_ages,s]) %*% S_plus
      }

      source = T_penult %*% Init_NAA[,n_ages-1,s] # compute forward projection of penultimate age
      Init_NAA[,n_ages,s] = solve(diag(n_regions) - T_plus, source)

    } # end s loop

    # save result
    NAA = Init_NAA

  } # end matrix approach

  # Scalar approach for last age, but with movement in preceeding ages
  if(init_age_strc == 3) {

    # projection initial abundance forward
    for(i in 1:n_ages) {
      for(s in 1:n_sexes) {
        for(seas in 1:n_seas) {

          if(seas == 1) Init_NAA[,1,s] = R0_r * sexratio[,s] # initialize recruitment

          # movement
          if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,seas,a,s] # recruits don't move
          if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,seas,a,s] # recruits move

          for(r in 1:n_regions) {
            # within season mortality
            if(seas < n_seas) {
              Init_NAA[r,1:n_ages,s] = Init_NAA[r,1:n_ages,s] *
                exp(-((natmort[r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
            } else {
              tmp_plus_befage = Init_NAA[r,n_ages,s] # save temporary plus group before ageing
              # ageing and mortality (age advancement)
              Init_NAA[r,2:n_ages,s] = Init_NAA[r,1:(n_ages-1),s] * exp(-((natmort[r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
              # accumulate plus group
              Init_NAA[r,n_ages,s] = (Init_NAA[r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[r,n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,n_ages,s,1]))))
            }

          } # end r loop
        } # end seas loop
      } # end s loop
    } # end i loop

    # Set up analytical solution for plus group
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        # Plus group - scalar geometric series (summing init_F to get annual Z)
        Z_penult = natmort[r,n_ages-1,s] + (sum(init_F) * fish_sel[r,n_ages-1,s,1])
        Z_plus = natmort[r,n_ages,s] + (sum(init_F) * fish_sel[r,n_ages,s,1])
        Init_NAA[r,n_ages,s] = Init_NAA[r,n_ages-1,s] * exp(-Z_penult) / (1 - exp(-Z_plus))
      } # end s loop
    } # end r loop

    # save result
    NAA = Init_NAA
  }

  # Input r0 into first age
  NAA[,1,] = R0_r * as.vector(sexratio)

  # Apply initial age deviations
  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      NAA[r,2:n_ages,s] = NAA[r,2:n_ages,s] * exp(ln_InitDevs[r,])
    } # end s loop
  } # end r loop

  return(NAA)

}
