# Stage 2 of 3: objective function
#
# Initial numbers at age: the equilibrium age structure the model starts from, and how the
# plus group is closed when there is movement.

#' Initial Numbers-at-Age (NAA)
#'
#' Numbers at age in the first model year, by population, region, age and sex, at the
#' equilibrium implied by constant recruitment, mortality and movement.
#' \code{init_age_strc} chooses how that equilibrium is solved, and \code{ln_InitDevs}
#' then deviates each age away from it.
#'
#' @param init_age_strc Integer, how the initial age structure is solved: \code{0} iterates
#'   the annual cycle \code{init_iter} times, \code{1} scalar geometric series with no
#'   movement at any age, \code{2} matrix geometric series with movement at every age,
#'   \code{3} movement below the plus group and a scalar series for the plus group,
#'   \code{4} no equilibrium at all, ages 2 and older are \code{exp(ln_InitDevs)}
#'   apportioned by sex ratio.
#' @param init_iter Integer, annual cycles run when \code{init_age_strc = 0}.
#' @param n_regions,n_pop,n_sexes,n_ages,n_seas,n_fish_fleets Integer dimensions.
#'   \code{n_ages} includes the plus group.
#' @param seasdur Numeric vector (\code{n_seas}) of each season's fraction of a year.
#' @param rec_seas_prop Matrix (\code{n_pop x n_seas}) of the share of annual recruitment
#'   entering in each season.
#' @param natmort Array (\code{n_pop x n_regions x n_seas x n_ages x n_sexes}) of natural
#'   mortality, a rate per year applied within each season.
#' @param natmort_annual Array (\code{n_pop x n_regions x n_ages x n_sexes}) of the annual
#'   total, the duration weighted sum over seasons, read by the steps that advance a whole
#'   year at once. Defaults to that sum.
#' @param init_F Numeric array (\code{n_regions x n_seas x n_fish_fleets}) of fully selected
#'   fishing mortality during initialization. Zero for an unfished population.
#' @param dmr Numeric array (\code{n_regions x n_seas x n_fish_fleets}) of the discard
#'   mortality rate during initialization.
#' @param fish_sel,ret_sel Arrays
#'   (\code{n_pop x n_regions x n_seas x n_ages x n_sexes x n_fish_fleets}) of total fishery
#'   selectivity at age and of the proportion of those fish retained.
#' @param R0_r Matrix (\code{n_pop x n_regions}) of unfished recruitment allocated to each region.
#' @param sexratio Array (\code{n_pop x n_regions x n_sexes}) of the proportion of recruits by sex.
#' @param Movement Array
#'   (\code{n_pop x n_regions x n_regions x n_seas x n_ages x n_sexes}) of seasonal movement
#'   probabilities, where \code{Movement[p,r,r2,,,]} is the fraction of the fish in region
#'   \code{r} that move to region \code{r2}.
#' @param do_recruits_move Integer, \code{0} recruits stay in their region for their first
#'   year, \code{1} recruits move with every other age.
#' @param ln_InitDevs Array (\code{n_pop x n_regions x (n_ages - 1) x n_sexes}) of log scale
#'   deviations for ages 2 and older. A 3-D array without the sex dimension is expanded
#'   across sexes as one shared curve. Under \code{init_age_strc = 4} these are the numbers
#'   themselves rather than multipliers on an equilibrium.
#' @param Mrate Array dimensioned like \code{Movement} of instantaneous movement rates.
#'   Required when \code{move_timing = 2}, ignored otherwise.
#' @param move_timing Integer ordering of movement and mortality within a season: \code{0}
#'   movement then mortality, \code{1} mortality then movement, \code{2} both at once. See
#'   \code{\link{build_seas_operator}}.
#' @param expm_nsub Integer, how the matrix exponential is taken under \code{move_timing = 2}:
#'   \code{0} uses \code{Matrix::expm}, \eqn{n \ge 1} the implicit backward Euler scheme. See
#'   \code{\link{mat_exp}}.
#'
#' @return Array (\code{n_pop x n_regions x n_ages x n_sexes}) of initial numbers at age.
#'
#' @keywords internal
Get_Init_NAA <- function(
  init_age_strc,
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
  natmort_annual = collapse_natmort_annual(natmort, seasdur, seas_dim = 3),
  init_F,
  dmr,
  fish_sel,
  ret_sel,
  R0_r,
  sexratio,
  Movement,
  do_recruits_move,
  ln_InitDevs,
  Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0
) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # A 3-D deviation array (the layout before the sex dimension existed, still
  # used by the simulation, which draws one shared curve) broadcasts across sexes
  if(length(dim(ln_InitDevs)) == 3) ln_InitDevs = array(rep(ln_InitDevs, n_sexes), dim = c(dim(ln_InitDevs), n_sexes))

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
          tmp_ret_F = rowSums(sweep(
            array(fish_sel[p,r,1,1:(n_ages-1),s,, drop=FALSE] * ret_sel[p,r,1,1:(n_ages-1),s,, drop=FALSE],
                  dim = c(n_ages - 1, n_fish_fleets)),
            2, as.vector(init_F[r,1,]), "*"
          ))
          # discarded F
          tmp_disc_F = rowSums(sweep(
            array(fish_sel[p,r,1,1:(n_ages-1),s,, drop=FALSE] *
                    (1 - ret_sel[p,r,1,1:(n_ages-1),s,, drop=FALSE]) *
                    dmr[r,1,],
                  dim = c(n_ages - 1, n_fish_fleets)),
            2, as.vector(init_F[r,1,]), "*"
          ))
          tmp_F = tmp_ret_F + tmp_disc_F # total F
          tmp_cumsum_Z = cumsum(natmort_annual[p,r,1:(n_ages-1),s] + tmp_F)
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
            # movement (applied here only under move_timing == 0; timings 1 and 2 fold it
            # into the seasonal transition operator below)
            if(move_timing == 0) {
              if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits don't move
              if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits move
            }
            # Apply mortality
            Z_ra = matrix(0, n_regions, n_ages) # total mortality by region and age, for the operator branch
            for(r in 1:n_regions) {
              # get tmp F
              tmp_ret_F = rowSums(sweep(
                array(fish_sel[p,r,seas,1:n_ages,s,, drop=FALSE] * ret_sel[p,r,seas,1:n_ages,s,, drop=FALSE],
                      dim = c(n_ages, n_fish_fleets)),
                2, as.vector(init_F[r,seas,]), "*"
              ))
              tmp_disc_F = rowSums(sweep(
                array(fish_sel[p,r,seas,1:n_ages,s,, drop=FALSE] * (1 - ret_sel[p,r,seas,1:n_ages,s,, drop=FALSE]) * dmr[r,seas,],
                      dim = c(n_ages, n_fish_fleets)),
                2, as.vector(init_F[r,seas,]), "*"
              ))
              tmp_F = tmp_ret_F + tmp_disc_F
              Z_ra[r,] = (natmort[p,r,seas,1:n_ages,s] * seasdur[seas]) + tmp_F
              # mortality wtihin season
              if(move_timing == 0) {
                if(seas < n_seas) {
                  Init_NAA_next_year[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                    exp(-((natmort[p,r,seas,1:n_ages,s] * seasdur[seas]) + tmp_F ))
                } else {
                  # ageing and mortality (advance ages in the next year)
                  Init_NAA_next_year[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] *
                    exp(-((natmort[p,r,seas,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)] ))
                  # accumulate plus group
                  Init_NAA_next_year[p,r,n_ages,s] = (Init_NAA_next_year[p,r,n_ages,s]) +
                    (Init_NAA[p,r,n_ages,s] * exp(-((natmort[p,r,seas,n_ages,s] * seasdur[seas]) + tmp_F[n_ages] )))
                } # end else
              } # end if move_timing == 0
            } # end r loop

            # Movement and mortality together for timings 1 and 2
            if(move_timing != 0) {
              step_ra = matrix(0, n_regions, n_ages)
              for(a in 1:n_ages) {
                moves = (do_recruits_move == 1 || a > 1)
                Mv = if(moves) Movement[p,,,seas,a,s] else diag(n_regions)
                Qv = if(moves) Mrate[p,,,seas,a,s] else matrix(0, n_regions, n_regions)
                step_ra[,a] = advance_seas(Init_NAA[p,,a,s], Mv, Z_ra[,a], Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
              } # end a loop
              if(seas < n_seas) {
                Init_NAA_next_year[p,,1:n_ages,s] = step_ra
              } else {
                Init_NAA_next_year[p,,2:n_ages,s] = step_ra[,1:(n_ages-1)]
                Init_NAA_next_year[p,,n_ages,s] = Init_NAA_next_year[p,,n_ages,s] + step_ra[,n_ages]
              } # end else
            } # end if move_timing != 0

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
                  exp(-((natmort[p,r,seas,1:n_ages,s] * seasdur[seas]) + tmp_F ))
              } else {
                tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                # ageing and mortality (age advancement)
                Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,seas,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)] ))
                # accumulate plus group
                Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,seas,n_ages,s] * seasdur[seas]) + tmp_F[n_ages] )))
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
          Z_penult = natmort_annual[p,r,n_ages-1,s] + F_annual_penult
          Z_plus = natmort_annual[p,r,n_ages,s] + F_annual_plus
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
            # movement (applied here only under move_timing == 0; timings 1 and 2 fold it
            # into the seasonal transition operator below)
            if(move_timing == 0) {
              if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits don't move
              if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits move
            }
            Z_ra = matrix(0, n_regions, n_ages) # total mortality by region and age, for the operator branch
            for(r in 1:n_regions) {
              tmp_ret_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * ret_sel[p,r,seas,1:n_ages,s,], dim = c(n_ages, n_fish_fleets)),
                                        2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # retained F
              tmp_disc_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * (1 - ret_sel[p,r,seas,1:n_ages,s,]) * dmr[r,seas,], dim = c(n_ages, n_fish_fleets)),
                                         2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # discarded F
              tmp_F = tmp_ret_F + tmp_disc_F # total F
              Z_ra[r,] = (natmort[p,r,seas,1:n_ages,s] * seasdur[seas]) + tmp_F
              # within season mortality
              if(move_timing == 0) {
                if(seas < n_seas) {
                  Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                    exp(-((natmort[p,r,seas,1:n_ages,s] * seasdur[seas]) + tmp_F ))
                } else {
                  tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                  # ageing and mortality (age advancement)
                  Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,seas,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)]))
                  # accumulate plus group
                  Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,seas,n_ages,s] * seasdur[seas]) + tmp_F[n_ages])))
                }
              } # end if move_timing == 0
            } # end r loop

            # Movement and mortality together for timings 1 and 2
            if(move_timing != 0) {
              step_ra = matrix(0, n_regions, n_ages)
              for(a in 1:n_ages) {
                moves = (do_recruits_move == 1 || a > 1)
                Mv = if(moves) Movement[p,,,seas,a,s] else diag(n_regions)
                Qv = if(moves) Mrate[p,,,seas,a,s] else matrix(0, n_regions, n_regions)
                step_ra[,a] = advance_seas(Init_NAA[p,,a,s], Mv, Z_ra[,a], Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
              } # end a loop
              if(seas < n_seas) {
                Init_NAA[p,,1:n_ages,s] = step_ra
              } else {
                Init_NAA[p,,2:n_ages,s] = step_ra[,1:(n_ages-1)] # ageing
                Init_NAA[p,,n_ages,s] = Init_NAA[p,,n_ages,s] + step_ra[,n_ages] # accumulate plus group
              }
            } # end if move_timing != 0
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
          # Column-convention seasonal operator (build_seas_operator returns row convention),
          # left-composed so that season 1 is applied first.
          Z_penult = (natmort[p,,seas,n_ages-1,s] * seasdur[seas]) + F_penult
          Z_plus_s = (natmort[p,,seas,n_ages,s] * seasdur[seas]) + F_plus
          Qp = if(is.null(Mrate)) NULL else Mrate[p,,,seas,n_ages-1,s]
          Ql = if(is.null(Mrate)) NULL else Mrate[p,,,seas,n_ages,s]
          T_penult = t(build_seas_operator(Movement[p,,,seas,n_ages-1,s], Z_penult, Qp, seasdur[seas], move_timing, expm_nsub = expm_nsub)) %*% T_penult
          T_plus = t(build_seas_operator(Movement[p,,,seas,n_ages,s], Z_plus_s, Ql, seasdur[seas], move_timing, expm_nsub = expm_nsub)) %*% T_plus
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
            # movement (applied here only under move_timing == 0; timings 1 and 2 fold it
            # into the seasonal transition operator below)
            if(move_timing == 0) {
              if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits don't move
              if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[p,,a,s] = t(Init_NAA[p,,a,s]) %*% Movement[p,,,seas,a,s] # recruits move
            }
            Z_ra = matrix(0, n_regions, n_ages) # total mortality by region and age, for the operator branch
            for(r in 1:n_regions) {
              tmp_ret_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * ret_sel[p,r,seas,1:n_ages,s,], dim = c(n_ages, n_fish_fleets)),
                                        2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # retained F
              tmp_disc_F = rowSums(sweep(array(fish_sel[p,r,seas,1:n_ages,s,] * (1 - ret_sel[p,r,seas,1:n_ages,s,]) * dmr[r,seas,], dim = c(n_ages, n_fish_fleets)),
                                         2, array(init_F[r,seas,], dim = n_fish_fleets), "*")) # discarded F
              tmp_F = tmp_ret_F + tmp_disc_F # total F
              Z_ra[r,] = (natmort[p,r,seas,1:n_ages,s] * seasdur[seas]) + tmp_F
              # within season mortality
              if(move_timing == 0) {
                if(seas < n_seas) {
                  Init_NAA[p,r,1:n_ages,s] = Init_NAA[p,r,1:n_ages,s] *
                    exp(-((natmort[p,r,seas,1:n_ages,s] * seasdur[seas]) + tmp_F ))
                } else {
                  tmp_plus_befage = Init_NAA[p,r,n_ages,s] # save temporary plus group before ageing
                  # ageing and mortality (age advancement)
                  Init_NAA[p,r,2:n_ages,s] = Init_NAA[p,r,1:(n_ages-1),s] * exp(-((natmort[p,r,seas,1:(n_ages-1),s] * seasdur[seas]) + tmp_F[1:(n_ages-1)]))
                  # accumulate plus group
                  Init_NAA[p,r,n_ages,s] = (Init_NAA[p,r,n_ages,s]) + (tmp_plus_befage * exp(-((natmort[p,r,seas,n_ages,s] * seasdur[seas]) + tmp_F[n_ages])))
                }
              } # end if move_timing == 0
            } # end r loop

            # Movement and mortality together for timings 1 and 2
            if(move_timing != 0) {
              step_ra = matrix(0, n_regions, n_ages)
              for(a in 1:n_ages) {
                moves = (do_recruits_move == 1 || a > 1)
                Mv = if(moves) Movement[p,,,seas,a,s] else diag(n_regions)
                Qv = if(moves) Mrate[p,,,seas,a,s] else matrix(0, n_regions, n_regions)
                step_ra[,a] = advance_seas(Init_NAA[p,,a,s], Mv, Z_ra[,a], Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
              } # end a loop
              if(seas < n_seas) {
                Init_NAA[p,,1:n_ages,s] = step_ra
              } else {
                Init_NAA[p,,2:n_ages,s] = step_ra[,1:(n_ages-1)] # ageing
                Init_NAA[p,,n_ages,s] = Init_NAA[p,,n_ages,s] + step_ra[,n_ages] # accumulate plus group
              }
            } # end if move_timing != 0
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
          Z_penult = natmort_annual[p,r,n_ages-1,s] + F_annual_penult
          Z_plus = natmort_annual[p,r,n_ages,s] + F_annual_plus
          Init_NAA[p,r,n_ages,s] = Init_NAA[p,r,n_ages-1,s] * exp(-Z_penult) / (1 - exp(-Z_plus))
        } # end s loop
      } # end r loop
    } # end p loop
    # save result
    NAA = Init_NAA
  }

  # free initial numbers at age: nothing is projected, so the deviations are the numbers rather
  # than multipliers on an equilibrium. seeding with the sex ratio reuses the shared step below
  if(init_age_strc == 4) {
    for(p in 1:n_pop) for(r in 1:n_regions) for(s in 1:n_sexes) NAA[p,r,2:n_ages,s] = sexratio[p,r,s]
  }

  # Overwrite first age and apply age deviations
  for(p in 1:n_pop) {
    # Overwrite first age
    for(s in 1:n_sexes) NAA[p,,1,s] <- R0_r[p,] * sexratio[p,,s] * rec_seas_prop[p,1]
    # Apply age deviations; sexes mapped to one shared parameter have identical values here
    for(r in 1:n_regions) for(s1 in 1:n_sexes) NAA[p,r,2:n_ages,s1] <- NAA[p,r,2:n_ages,s1] * exp(ln_InitDevs[p,r,,s1])
  } # end p loop

  return(NAA)
}
