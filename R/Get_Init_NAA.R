#' Initialize Numbers-at-Age (NAA) for a Population Model
#'
#' This function generates initial numbers-at-age (NAA) for a structured population
#' model across populations, regions, sexes, and age classes. It supports multiple
#' initialization methods, including iterative solution, scalar geometric series, and
#' matrix geometric series, optionally accounting for movement and fishing mortality.
#' Initial age deviations can also be applied.
#'
#' @param init_age_strc Integer specifying the initialization method for the age structure:
#'   - 0: Iterative solution to equilibrium
#'   - 1: Scalar geometric series solution w/o movement in any groups (no movement in all groups)
#'   - 2: Matrix geometric series solution (generalizes scalar solution with movement)
#'   - 3: Scalar geometric series solution w/o movement only in plus group (no movement in plus groups)
#' @param init_iter Integer; number of iterations to run when `init_age_strc = 0`.
#' @param n_pop Integer; number of populations.
#' @param n_regions Integer; number of spatial regions.
#' @param n_sexes Integer; number of sexes.
#' @param n_ages Integer; number of age classes.
#' @param n_seas Integer; number of seasons.
#' @param seasdur Numeric vector of length `n_seas`; fraction of the year within each season.
#' @param natmort Array of natural mortality rates with dimensions `[pop, regions, ages, sexes]`.
#' @param init_F Numeric vector of length `n_seas`; initial fishing mortality applied per season
#'   (set to 0 for an unfished population).
#' @param fish_sel Array of fishery selectivity with dimensions `[regions, ages, sexes, fleets]`.
#'   Only fleet 1 is used during initialization.
#' @param R0_r Array of virgin recruitment values with dimensions `[pop, regions]`.
#' @param sexratio Array `[pop, regions, sexes]` giving the proportion of each sex in recruitment.
#' @param Movement Array `[pop, regions, regions, seas, ages, sexes]` defining movement probabilities.
#' @param do_recruits_move Integer; 0 = recruits do not move, 1 = recruits move according to `Movement`.
#' @param ln_InitDevs Array `[pop, regions, ages-1]` of log-scale deviations for initial
#'   numbers-at-age. Applied to ages 2 through `n_ages`.
#'
#' @return Array of initial numbers-at-age with dimensions `[pop, regions, ages, sexes]`.
#'
#' @details
#' Initial numbers-at-age are derived under equilibrium recruitment
#' \eqn{R_0} and constant mortality and movement.
#'
#' Let:
#' \itemize{
#'   \item \eqn{N_{p,r,a,s}} = numbers-at-age
#'   \item \eqn{M_{p,r,a}} = natural mortality
#'   \item \eqn{F_{r,a,s}} = fishing mortality
#'   \item \eqn{Z = M + F} = total mortality
#' }
#'
#' Recruitment at age 1 is:
#' \deqn{
#' N_{p,r,1,s=1} = R_{0,p,r} \times sexratio_{p,r}
#' }
#'
#' Within-season survival is:
#' \deqn{
#' N_{p,r,a,s+1} =
#' N_{p,r,a,s} \exp(-Z_{p,r,a,s})
#' }
#'
#' Age advancement at the end of the final season:
#' \deqn{
#' N_{p,r,a+1,1} =
#' N_{p,r,a,n_{seas}}
#' \exp(-Z_{p,r,a,n_{seas}})
#' }
#'
#' The plus group accumulates survivors:
#' \deqn{
#' N_{A^+} =
#' N_{A-1}e^{-Z_{A-1}} +
#' N_{A^+}e^{-Z_{A^+}}
#' }
#'
#' ### Scalar geometric solution (no movement)
#'
#' When movement is absent, equilibrium abundance follows:
#'
#' \deqn{
#' N_a = N_1 \exp\left(-\sum_{i=1}^{a-1} Z_i \right)
#' }
#'
#' The plus group has closed-form solution:
#'
#' \deqn{
#' N_{A^+} =
#' \frac{N_{A-1} e^{-Z_{A-1}}}
#' {1 - e^{-Z_{A^+}}}
#' }
#'
#' ### Matrix geometric solution (with movement)
#'
#' When seasonal movement occurs, annual survival is represented by
#' transition matrices:
#'
#' \deqn{
#' \mathbf{T}_a =
#' \prod_{s=1}^{n_{seas}}
#' \mathbf{M}_{a,s} \mathbf{S}_{a,s}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\mathbf{M}} = movement transition matrix
#'   \item \eqn{\mathbf{S}} = diagonal survival matrix
#' }
#'
#' The plus group equilibrium is obtained from:
#'
#' \deqn{
#' \mathbf{N}_{A^+} =
#' (\mathbf{I} - \mathbf{T}_{A^+})^{-1}
#' \mathbf{T}_{A-1}
#' \mathbf{N}_{A-1}
#' }
#'
#' The iterative method (`init_age_strc = 0`) numerically applies the
#' full annual cycle repeatedly until convergence.
#'
#' Log-scale initial age deviations are applied multiplicatively to
#' ages 2 through \eqn{A} after equilibrium is derived.
#'
#' @keywords internal
Get_Init_NAA <- function(init_age_strc,
                         init_iter,
                         n_regions,
                         n_pop,
                         n_sexes,
                         n_ages,
                         n_seas,
                         seasdur,
                         rec_seas_prop,
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
  Init_NAA = array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
  Init_NAA_next_year = array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
  NAA = array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))

  # Iterative Solution
  if(init_age_strc == 0) {

    # initialize age structure (starting point)
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          tmp_cumsum_Z = cumsum(natmort[p,r,1:(n_ages-1),s] + sum(init_F) * fish_sel[r,1:(n_ages-1),s,1])
          Init_NAA[p,r,,s] = c(R0_r[p,r] * sexratio[p,r,s] * rec_seas_prop[p,1], R0_r[p,r] * sexratio[p,r,s] * rec_seas_prop[p,1] * exp(-tmp_cumsum_Z))
        } # end s loop
      } # end r loop
    }

    # Apply annual cycle and iterate to equilibrium
    for(i in 1:init_iter) {
      for(p in 1:n_pop) {
        for(s in 1:n_sexes) {
          for(seas in 1:n_seas) {

            # recruitment in the first season
            if(seas == 1) Init_NAA[p,,1,s] = R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,1]
            else Init_NAA[p,,1,s] = Init_NAA[p,,1,s] + (R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,seas]) # recruitment not in the first season

            # movement
            if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits don't move
            if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits move

            # Apply mortality
            for(r in 1:n_regions) {

              # mortality wtihin season
              if(seas < n_seas) {
                Init_NAA_next_year[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
              } else {
                # ageing and mortality (advance ages in the next year)
                Init_NAA_next_year[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] *
                  exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
                # accumulate plus group
                Init_NAA_next_year[p,r,n_ages,s] = (Init_NAA_next_year[p,r,n_ages,s]) +
                  (Init_NAA[p,r,n_ages,s] * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,n_ages,s,1]))))
              } # end else
            } # end r loop

            Init_NAA = Init_NAA_next_year # iterate to next cycle

          } # end seas loop
        } # end s loop
      } # end p loop
    } # end i loop

    # save result
    NAA[] = Init_NAA

  } # end if iterative solution


  # Scalar Geometric Series Solution (no movement in all ages)
  if(init_age_strc == 1) {

    # projection initial abundance forward
    for(p in 1:n_pop) {
      for(i in 1:n_ages) {
        for(s in 1:n_sexes) {
          for(seas in 1:n_seas) {

            # recruitment in the first season
            if(seas == 1) Init_NAA[p,,1,s] = R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,1]
            else Init_NAA[p,,1,s] = Init_NAA[p,,1,s] + (R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,seas]) # recruitment not in the first season

            for(r in 1:n_regions) {
              # within season mortality
              if(seas < n_seas) {
                Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
              } else {
                tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                # ageing and mortality (age advancement)
                Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
                # accumulate plus group
                Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,n_ages,s,1]))))
              }

            } # end r loop
          } # end seas loop
        } # end s loop
      } # end i loop

      # Set up analytical solution for plus group
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          # Plus group - scalar geometric series (summing init_F to get annual Z)
          Z_penult = natmort[p,r,n_ages-1,s] + (sum(init_F) * fish_sel[r,n_ages-1,s,1])
          Z_plus = natmort[p,r,n_ages,s] + (sum(init_F) * fish_sel[r,n_ages,s,1])
          Init_NAA[p,r,n_ages,s] = Init_NAA[p,r,n_ages-1,s] * exp(-Z_penult) / (1 - exp(-Z_plus))
        } # end s loop
      } # end r loop

    } # end p loop

    # save result
    NAA = Init_NAA
  } # end if

  # Matrix Geometric Series Solution (genearlizes to scalar w/o movement)
  if(init_age_strc == 2) {

    # projection initial abundance forward
    for(p in 1:n_pop) {
      for(i in 1:n_ages) {
        for(s in 1:n_sexes) {
          for(seas in 1:n_seas) {

            # recruitment in the first season
            if(seas == 1) Init_NAA[p,,1,s] = R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,1]
            else Init_NAA[p,,1,s] = Init_NAA[p,,1,s] + (R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,seas]) # recruitment not in the first season

            # movement
            if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits don't move
            if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits move

            for(r in 1:n_regions) {
              # within season mortality
              if(seas < n_seas) {
                Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
              } else {
                tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                # ageing and mortality (age advancement)
                Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
                # accumulate plus group
                Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,n_ages,s,1]))))
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
          S_penult = diag(exp(-((natmort[p,,n_ages-1,s] * seasdur[seas]) + init_F[seas] * fish_sel[,n_ages-1,s,1])), n_regions)
          S_plus = diag(exp(-((natmort[p,,n_ages,s] * seasdur[seas]) + init_F[seas] * fish_sel[,n_ages,s,1])), n_regions)
          T_penult = T_penult %*% t(Movement[p,,,seas,n_ages-1,s]) %*% S_penult
          T_plus = T_plus %*% t(Movement[p,,,seas,n_ages,s]) %*% S_plus
        }

        source = T_penult %*% Init_NAA[p,,n_ages-1,s] # compute forward projection of penultimate age
        Init_NAA[p,,n_ages,s] = solve(diag(n_regions) - T_plus, source)

      } # end s loop

    } # end p loop

    # save result
    NAA = Init_NAA

  } # end matrix approach

  # Scalar approach for last age, but with movement in preceeding ages
  if(init_age_strc == 3) {

    # projection initial abundance forward
    for(p in 1:n_pop) {
      for(i in 1:n_ages) {
        for(s in 1:n_sexes) {
          for(seas in 1:n_seas) {

            # recruitment in the first season
            if(seas == 1) Init_NAA[p,,1,s] = R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,1]
            else Init_NAA[p,,1,s] = Init_NAA[p,,1,s] + (R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,seas]) # recruitment not in the first season

            # movement
            if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits don't move
            if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits move

            for(r in 1:n_regions) {
              # within season mortality
              if(seas < n_seas) {
                Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:n_ages,s,1])))
              } else {
                tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                # ageing and mortality (age advancement)
                Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + (init_F[seas] * fish_sel[r,1:(n_ages-1),s,1])))
                # accumulate plus group
                Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) +
                                                                                               (init_F[seas] * fish_sel[r,n_ages,s,1]))))
              }

            } # end r loop
          } # end seas loop
        } # end s loop
      } # end i loop

      # Set up analytical solution for plus group
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          # Plus group - scalar geometric series (summing init_F to get annual Z)
          Z_penult = natmort[p,r,n_ages-1,s] + (sum(init_F) * fish_sel[r,n_ages-1,s,1])
          Z_plus = natmort[p,r,n_ages,s] + (sum(init_F) * fish_sel[r,n_ages,s,1])
          Init_NAA[p,r,n_ages,s] = Init_NAA[p,r,n_ages-1,s] * exp(-Z_penult) / (1 - exp(-Z_plus))
        } # end s loop
      } # end r loop

    } # end p loop
    # save result
    NAA = Init_NAA
  }


  # Overwrite firsta ge and apply age deviations
  for(p in 1:n_pop) {
    # Overwrite first age
    for(s in 1:n_sexes) NAA[p,,1,s] <- R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,1]
    # Apply age deviations
    for(r in 1:n_regions) for(s1 in 1:n_sexes) NAA[p,r,2:n_ages,s1] <- NAA[p,r,2:n_ages,s1] * exp(ln_InitDevs[p,r,])
  } # end p loop


  return(NAA)

}
