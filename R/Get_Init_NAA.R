#' Initialize Numbers-at-Age (NAA) for a Population Model
#'
#' Computes equilibrium initial numbers-at-age (NAA) for a structured population
#' model across populations, regions, sexes, and ages. Several initialization
#' methods are supported, ranging from numerical iteration to analytical
#' geometric-series solutions that optionally incorporate seasonal movement.
#'
#' The resulting age structure represents an equilibrium population under
#' constant recruitment (\eqn{R_0}), mortality, fishing mortality, and movement.
#'
#' @param init_age_strc Integer specifying the initialization method:
#' \itemize{
#' \item \code{0} Iterative equilibrium solution
#' \item \code{1} Scalar geometric-series solution (no movement in any age)
#' \item \code{2} Matrix geometric-series solution (movement allowed)
#' \item \code{3} Hybrid solution: movement in ages < plus group, scalar solution for plus group
#' }
#'
#' @param init_iter Integer; number of annual iterations used when
#'   \code{init_age_strc = 0}.
#'
#' @param n_regions Integer; number of spatial regions.
#' @param n_pop Integer; number of populations.
#' @param n_sexes Integer; number of sexes.
#' @param n_ages Integer; number of age classes (including the plus group).
#' @param n_seas Integer; number of seasons per year.
#' @param n_fish_fleets Integer; number of fishing fleets.
#'
#' @param seasdur Numeric vector (\code{n_seas}) giving the fraction of the year
#'   represented by each season.
#'
#' @param rec_seas_prop Matrix (\code{n_pop x n_seas}) giving seasonal recruitment
#'   proportions.
#'
#' @param natmort Array (\code{n_pop x n_regions x n_ages x n_sexes}) of natural
#'   mortality rates.
#'
#' @param init_F Numeric array (\code{n_regions x n_seas x n_fish_fleets})
#'   giving fishing mortality applied in each region, season, and fleet during
#'   initialization. Set to zero for an unfished population.
#'
#' @param fish_sel Array
#'   (\code{n_pop x n_regions x n_seas x n_ages x n_sexes x n_fish_fleets})
#'   of total fishery selectivity at age.
#'
#' @param R0_r Matrix (\code{n_pop x n_regions}) giving unfished recruitment
#'   allocated to each region.
#'
#' @param sexratio Array (\code{n_pop x n_regions x n_sexes}) giving the
#'   proportion of recruits by sex.
#'
#' @param Movement Array
#'   (\code{n_pop x origin x destination x n_seas x n_ages x n_sexes})
#'   containing seasonal movement probabilities.
#'
#' @param do_recruits_move Integer indicator:
#' \itemize{
#' \item \code{0} Recruits do not move during their first year
#' \item \code{1} Recruits move according to the movement matrix
#' }
#'
#' @param ln_InitDevs Array (\code{n_pop x n_regions x (n_ages - 1)}) containing
#'   log-scale deviations applied to ages 2 through \eqn{A}.
#'
#' @param dmr Numeric array (\code{n_regions x n_seas x n_fish_fleets})
#'   giving discard mortality rate applied in each region, season, and fleet during
#'   initialization (first year). Set to zero for an unfished population.
#' @param ret_sel Array
#'   (\code{n_pop x n_regions x n_seas x n_ages x n_sexes x n_fish_fleets})
#'   of retained fishery selectivity at age.
#'
#' @details
#'
#' Initial numbers-at-age are derived assuming constant recruitment
#' (\eqn{R_0}) and constant mortality and movement.
#'
#' Let
#'
#' \itemize{
#' \item \eqn{N_{p,r,a,s}} denote numbers-at-age
#' \item \eqn{M_{p,r,a}} denote natural mortality
#' \item \eqn{F_{r,a,s}} denote total fishing mortality (retained + dead discards)
#' \item \eqn{Z = M + F} denote total mortality
#' }
#'
#' Recruitment at age 1 is
#'
#' \deqn{
#' N_{p,r,1,s} = R_{0,p,r} \times sexratio_{p,r,s}
#' }
#'
#' Within-season survival follows
#'
#' \deqn{
#' N_{p,r,a,s+1} =
#' N_{p,r,a,s}\exp(-Z_{p,r,a,s})
#' }
#'
#' Ages advance at the end of the final season of the year:
#'
#' \deqn{
#' N_{p,r,a+1,1} =
#' N_{p,r,a,n_{seas}}
#' \exp(-Z_{p,r,a,n_{seas}})
#' }
#'
#' The plus group accumulates survivors from the terminal age:
#'
#' \deqn{
#' N_{A^+} =
#' N_{A-1} e^{-Z_{A-1}} +
#' N_{A^+} e^{-Z_{A^+}}
#' }
#'
#' Fishing mortality at age is decomposed into retained and dead discard
#' components, summed across all fleets:
#'
#' \deqn{
#' F_{p,r,s,a} = \sum_{f=1}^{n_f} F^{init}_{r,s,f} \left[
#' sel_{p,r,s,a,f} \cdot ret_{p,r,s,a,f} +
#' sel_{p,r,s,a,f} \cdot (1 - ret_{p,r,s,a,f}) \cdot dmr_{r,s,f}
#' \right]
#' }
#'
#' where \eqn{sel} is total fishery selectivity, \eqn{ret} is retention
#' selectivity, and \eqn{dmr} is the discard mortality rate.
#'
#' ### Scalar geometric-series solution
#'
#' When movement is absent, equilibrium abundance follows
#'
#' \deqn{
#' N_a =
#' N_1 \exp\left(-\sum_{i=1}^{a-1} Z_i \right)
#' }
#'
#' The plus group has a closed-form solution
#'
#' \deqn{
#' N_{A^+} =
#' \frac{N_{A-1} e^{-Z_{A-1}}}
#' {1 - e^{-Z_{A^+}}}
#' }
#'
#' ### Matrix geometric-series solution
#'
#' When movement occurs, survival and movement are combined into
#' seasonal transition matrices:
#'
#' \deqn{
#' \mathbf{T}_a =
#' \prod_{s=1}^{n_{seas}}
#' \mathbf{M}_{a,s}\mathbf{S}_{a,s}
#' }
#'
#' where
#'
#' \itemize{
#' \item \eqn{\mathbf{M}} is the movement transition matrix
#' \item \eqn{\mathbf{S}} is a diagonal matrix of survival probabilities
#' }
#'
#' The plus group equilibrium satisfies
#'
#' \deqn{
#' \mathbf{N}_{A^+} =
#' (\mathbf{I} - \mathbf{T}_{A^+})^{-1}
#' \mathbf{T}_{A-1}
#' \mathbf{N}_{A-1}
#' }
#'
#' The iterative method (\code{init_age_strc = 0}) numerically applies
#' the full seasonal population dynamics repeatedly until the population
#' converges to equilibrium.
#'
#' After equilibrium is derived, log-scale initial age deviations
#' (\code{ln_InitDevs}) are applied multiplicatively to ages
#' \eqn{2,\dots,A}.
#'
#' @keywords internal
Get_Init_NAA <- function(init_age_strc,
                         init_iter,
                         n_regions,
                         n_pop,
                         n_sexes,
                         n_ages,
                         n_seas,
                         n_fish_fleets,
                         seasdur,
                         rec_seas_prop,
                         natmort,
                         init_F,
                         dmr,
                         fish_sel,
                         ret_sel,
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
          # retained F
          tmp_ret_F = rowSums(array(sweep(fish_sel[p,r,1,1:(n_ages-1),s,] * ret_sel[p,r,1,1:(n_ages-1),s,], 2, init_F[r,1,], "*"), dim = c(n_ages - 1, n_fish_fleets)))
          # discarded F
          tmp_disc_F = rowSums(array(sweep(fish_sel[p,r,1,1:(n_ages-1),s,] * (1 - ret_sel[p,r,1,1:(n_ages-1),s,]) * dmr[r,1,], 2, init_F[r,1,], "*"), dim = c(n_ages - 1, n_fish_fleets)))
          tmp_F = tmp_ret_F + tmp_disc_F # total F
          tmp_cumsum_Z = cumsum(natmort[p,r,1:(n_ages-1),s] + tmp_F)
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

              # get tmp F
              tmp_ret_F = rowSums(array(sweep(fish_sel[p,r,seas,1:n_ages,s,] * ret_sel[p,r,seas,1:n_ages,s,], 2, init_F[r,seas,], "*"), dim = c(n_ages, n_fish_fleets)))
              tmp_disc_F = rowSums(array(sweep(fish_sel[p,r,seas,1:n_ages,s,] * (1 - ret_sel[p,r,seas,1:n_ages,s,]) * dmr[r,seas,], 2, init_F[r,seas,], "*"), dim = c(n_ages, n_fish_fleets)))
              tmp_F = tmp_ret_F + tmp_disc_F

              # mortality wtihin season
              if(seas < n_seas) {
                Init_NAA_next_year[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + tmp_F ))
              } else {
                # ageing and mortality (advance ages in the next year)
                Init_NAA_next_year[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] *
                  exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)] ))
                # accumulate plus group
                Init_NAA_next_year[p,r,n_ages,s] = (Init_NAA_next_year[p,r,n_ages,s]) +
                  (Init_NAA[p,r,n_ages,s] * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) + tmp_F[n_ages] )))
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
              tmp_ret_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * ret_sel[p,r,seas,1:n_ages,s,], dim = c(n_ages, n_fish_fleets)),
                                        2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # retained F
              tmp_disc_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * (1 - ret_sel[p,r,seas,1:n_ages,s,]) * dmr[r,seas,], dim = c(n_ages, n_fish_fleets)),
                                         2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # discarded F
              tmp_F = tmp_ret_F + tmp_disc_F # total F
              # within season mortality
              if(seas < n_seas) {
                Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + tmp_F ))
              } else {
                tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                # ageing and mortality (age advancement)
                Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)] ))
                # accumulate plus group
                Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) + tmp_F[n_ages] )))
              }

            } # end r loop
          } # end seas loop
        } # end s loop
      } # end i loop

      # Set up analytical solution for plus group
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          # Plus group - scalar geometric series (summing annual F across seasons and fleets)
          F_annual_penult = sum(array(init_F[r,,] * (fish_sel[p,r,,n_ages-1,s,] * ret_sel[p,r,,n_ages-1,s,] + # retained
                                      fish_sel[p,r,,n_ages-1,s,] * (1 - ret_sel[p,r,,n_ages-1,s,]) * dmr[r,,]), # discarded
                                      dim = c(n_seas, n_fish_fleets)))
          F_annual_plus = sum(array(init_F[r,,] * (fish_sel[p,r,,n_ages,s,] * ret_sel[p,r,,n_ages,s,] + # retained
                                    fish_sel[p,r,,n_ages,s,] * (1 - ret_sel[p,r,,n_ages,s,]) * dmr[r,,]), # discarded
                                    dim = c(n_seas, n_fish_fleets)))
          Z_penult = natmort[p,r,n_ages-1,s] + F_annual_penult
          Z_plus = natmort[p,r,n_ages,s] + F_annual_plus
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
              tmp_ret_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * ret_sel[p,r,seas,1:n_ages,s,], dim = c(n_ages, n_fish_fleets)),
                                        2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # retained F
              tmp_disc_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * (1 - ret_sel[p,r,seas,1:n_ages,s,]) * dmr[r,seas,], dim = c(n_ages, n_fish_fleets)),
                                         2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # discarded F
              tmp_F = tmp_ret_F + tmp_disc_F # total F
              # within season mortality
              if(seas < n_seas) {
                Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + tmp_F ))
              } else {
                tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                # ageing and mortality (age advancement)
                Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)]))
                # accumulate plus group
                Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) + tmp_F[n_ages])))
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
          F_penult = rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages-1,s,] * ret_sel[p,,seas,n_ages-1,s,] + # retained
                                   fish_sel[p,,seas,n_ages-1,s,] * (1 - ret_sel[p,,seas,n_ages-1,s,]) * dmr[,seas,]), # discarded
                                   dim = c(n_regions, n_fish_fleets)))
          F_plus = rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages,s,] * ret_sel[p,,seas,n_ages,s,] + # retained
                                 fish_sel[p,,seas,n_ages,s,] * (1 - ret_sel[p,,seas,n_ages,s,]) * dmr[,seas,]), # discarded
                                 dim = c(n_regions, n_fish_fleets)))
          S_penult = diag(exp(-((natmort[p,,n_ages-1,s] * seasdur[seas]) + F_penult)), n_regions)
          S_plus = diag(exp(-((natmort[p,,n_ages,s] * seasdur[seas]) + F_plus)), n_regions)
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
              tmp_ret_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * ret_sel[p,r,seas,1:n_ages,s,], dim = c(n_ages, n_fish_fleets)),
                                        2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # retained F
              tmp_disc_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * (1 - ret_sel[p,r,seas,1:n_ages,s,]) * dmr[r,seas,], dim = c(n_ages, n_fish_fleets)),
                                         2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # discarded F
              tmp_F = tmp_ret_F + tmp_disc_F # total F
              # within season mortality
              if(seas < n_seas) {
                Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                  exp(-((natmort[p,r,1:n_ages,s] * seasdur[seas]) + tmp_F ))
              } else {
                tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                # ageing and mortality (age advancement)
                Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)]))
                # accumulate plus group
                Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,n_ages,s] * seasdur[seas]) + tmp_F[n_ages])))
              }

            } # end r loop
          } # end seas loop
        } # end s loop
      } # end i loop

      # Set up analytical solution for plus group
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          # Plus group - scalar geometric series (summing annual F across seasons and fleets)
          F_annual_penult = sum(array(init_F[r,,] * (fish_sel[p,r,,n_ages-1,s,] * ret_sel[p,r,,n_ages-1,s,] + # retained
                                                       fish_sel[p,r,,n_ages-1,s,] * (1 - ret_sel[p,r,,n_ages-1,s,]) * dmr[r,,]), # discarded
                                      dim = c(n_seas, n_fish_fleets)))
          F_annual_plus = sum(array(init_F[r,,] * (fish_sel[p,r,,n_ages,s,] * ret_sel[p,r,,n_ages,s,] + # retained
                                                     fish_sel[p,r,,n_ages,s,] * (1 - ret_sel[p,r,,n_ages,s,]) * dmr[r,,]), # discarded
                                    dim = c(n_seas, n_fish_fleets)))
          Z_penult = natmort[p,r,n_ages-1,s] + F_annual_penult
          Z_plus = natmort[p,r,n_ages,s] + F_annual_plus
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
