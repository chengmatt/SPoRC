library(SPoRC)
library(testthat)

# Equilibrium yield from the multi-region Beverton-Holt MSY reference points, checked against
# the equilibrium catch of a long projection held at the same F.
#
# The other refpts test files check spawning biomass, which never touches CAA. Yield does, and
# the catch side has two failure modes this file guards:
#
#   Seasons. The penultimate age and the plus group take their catch in a loop over seasons.
#   Abundance has to be carried forward through that loop, or every season after the first is
#   charged against fish that should already be dead. This is invisible at n_seas == 1, which
#   is what the rest of the suite runs at, and grows with the number of seasons.
#
#   Continuous movement. Under move_timing = 2 fish redistribute among regions while they are
#   being caught, so the region-local Baranov form is invalid and the season-integrated
#   (spatial Baranov) abundance is required instead.
#
# Everything is synthetic and deterministic: no operating model, no fitting. The reference
# points and Do_Population_Projection are independent implementations of the same dynamics, so
# agreement to machine precision is a real check rather than a restatement of one of them.
#
# Region-specific M is what makes the timings distinguishable. The three move_timings agree
# exactly when Z is constant across regions, so a symmetric setup would pass at every timing
# regardless of whether continuous movement were implemented at all.

n_regions <- 2; n_ages <- 12; n_pop <- 1; n_sexes <- 1; n_fish_fleets <- 1
n_proj_yrs <- 150

# CTMC generator in row convention (Q[from, to]), rows summing to zero
Qgen <- matrix(c(-0.5, 0.5, 0.9, -0.9), nrow = n_regions, byrow = TRUE)
frac_from_Q <- function(Q, dur) t(as.matrix(Matrix::expm(t(Q) * dur)))

waa <- 5 / (1 + exp(-3 * ((1:n_ages) - 3)))
mat <- 1 / (1 + exp(-3 * ((1:n_ages) - 3)))
sel <- 1 / (1 + exp(-1.5 * ((1:n_ages) - 4)))
M_r <- c(0.15, 0.45)
rec_region_prop <- c(0.6, 0.4)

# Solve Fmsy, then project at it and return refpts yield alongside realised equilibrium catch
msy_vs_projection <- function(n_seas, move_timing) {

  seasdur <- rep(1 / n_seas, n_seas)
  spawn_seas <- 1
  # F_fract_flt is a share of the annual F across seasons and fleets, so it must sum to 1
  # the way Get_Reference_Points builds it from the terminal year's Fmort
  seas_w <- seq_len(n_seas) / sum(seq_len(n_seas))

  Mov <- array(0, dim = c(n_regions, n_regions, n_seas, n_ages))
  Mrt <- array(0, dim = c(n_regions, n_regions, n_seas, n_ages))
  for(seas in 1:n_seas) for(a in 1:n_ages) {
    Mov[,,seas,a] <- frac_from_Q(Qgen, seasdur[seas])
    Mrt[,,seas,a] <- Qgen
  }

  data_list <- list(
    t_spawn = 0, n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
    n_pop = n_pop, n_ages = n_ages, n_regions = n_regions,
    F_fract_flt = array(rep(seas_w, each = n_regions), dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr      = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    fish_sel = array(rep(sel, each = n_pop * n_regions * n_seas),
                     dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    ret_sel  = array(1, dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    natmort  = array(rep(M_r, times = n_ages), dim = c(n_regions, n_ages)),
    WAA      = array(rep(waa, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    MatAA    = array(rep(mat, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    Movement = Mov,
    Mrate    = if(move_timing == 2) Mrt else array(0, dim = dim(Mrt)),
    move_timing = move_timing,
    is_discard_fleet = array(0, dim = n_fish_fleets),
    do_recruits_move = 1,
    rec_region_prop = rec_region_prop,
    sex_ratio_f = 1,  # matches sexratio = 1 below, so per-recruit yield is directly comparable
    rec_seas_prop = rep(1 / n_seas, n_seas),
    h = 0.7, R0 = 20
  )

  fit <- SPoRC:::optim_ref_pts(SPoRC:::global_BH_Fmsy, data_list, list(log_Fmsy = log(0.1)))
  Fmsy <- fit$rep$Fmsy; Req <- fit$rep$Req

  # Lay a season/age slice out over projection years, which sit in dimension 3
  lay <- function(v, dims) aperm(array(rep(as.vector(v), n_proj_yrs), dim = c(dims, n_proj_yrs)),
                                 c(1, 2, length(dims) + 1, 3:length(dims)))
  # Movement arrays put years in dimension 4 instead, after the two region dimensions
  lay_move <- function(v) aperm(array(rep(as.vector(v), n_proj_yrs),
                                      dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes, n_proj_yrs)),
                                c(1, 2, 3, 7, 4, 5, 6))

  Movp <- array(0, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)); Mrtp <- Movp
  for(seas in 1:n_seas) for(a in 1:n_ages) {
    Movp[1,,,seas,a,1] <- Mov[,,seas,a]
    Mrtp[1,,,seas,a,1] <- Mrt[,,seas,a]
  }

  natmort <- aperm(array(rep(rep(M_r, times = n_ages), n_pop * n_sexes * n_proj_yrs),
                         dim = c(n_regions, n_ages, n_pop, n_sexes, n_proj_yrs)), c(3, 1, 5, 2, 4))

  # Any starting state converges; a roughly declining age structure just gets there sooner
  init <- array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  for(a in 1:n_ages) init[1,,,a,1] <- outer(Req * rec_region_prop * exp(-0.3 * (a - 1)), rep(1, n_seas))

  out <- Do_Population_Projection(
    n_proj_yrs = n_proj_yrs, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages,
    n_sexes = n_sexes, sexratio = array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)),
    n_fish_fleets = n_fish_fleets, do_recruits_move = 1,
    recruitment = array(rep(Req * rec_region_prop, each = n_pop), dim = c(n_pop, n_regions, 5)),
    terminal_NAA = init, terminal_NAA0 = init,
    terminal_F = array(rep(seas_w, each = n_regions), dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    natmort = natmort,
    WAA      = lay(array(rep(waa, each = n_pop*n_regions*n_seas),
                         dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
                   c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
    WAA_fish = lay(array(rep(waa, each = n_pop*n_regions*n_seas),
                         dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
                   c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
    MatAA    = lay(array(rep(mat, each = n_pop*n_regions*n_seas),
                         dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
                   c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
    fish_sel = lay(array(rep(sel, each = n_pop*n_regions*n_seas),
                         dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
                   c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
    ret_sel  = array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
    Movement = lay_move(Movp),
    Mrate    = if(move_timing == 2) lay_move(Mrtp) else NULL,
    move_timing = move_timing,
    f_ref_pt = array(Fmsy, dim = c(n_regions, n_proj_yrs)),
    b_ref_pt = array(1, dim = c(n_pop, n_regions, n_proj_yrs)),
    HCR_function = function(x, frp, brp, alpha = 0.05) frp,  # hold F at Fmsy
    recruitment_opt = "mean_rec", fmort_opt = "HCR", t_spawn = 0,
    n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
    rec_seas_prop = array(1 / n_seas, dim = c(n_pop, n_seas))
  )

  list(Fmsy = Fmsy,
       refpts_yield = as.numeric(fit$rep$Yield_r),
       proj_catch = as.numeric(apply(out$proj_Catch[,, n_proj_yrs,,, drop = FALSE], 2, sum)),
       ssb_final = as.numeric(out$proj_SSB[,, n_proj_yrs]),
       ssb_penult = as.numeric(out$proj_SSB[,, n_proj_yrs - 1]))
}

test_that("multi-region BH MSY yield matches the equilibrium catch of a projection held at Fmsy", {
  for(n_seas in c(1, 2, 4)) {
    for(move_timing in c(0, 2)) {
      res <- msy_vs_projection(n_seas, move_timing)
      lbl <- paste("n_seas =", n_seas, ", move_timing =", move_timing)

      # the projection has to have settled, or the comparison is against a transient
      expect_equal(res$ssb_penult, res$ssb_final, tolerance = 1e-10, info = lbl)

      expect_equal(res$proj_catch, res$refpts_yield, tolerance = 1e-8, info = lbl)
    }
  }
})

test_that("continuous movement changes the MSY reference points it is supposed to change", {
  # Guards against the check above passing vacuously: if the timings produced identical
  # answers, every assertion would hold whether or not continuous movement did anything.
  discrete   <- msy_vs_projection(2, 0)
  continuous <- msy_vs_projection(2, 2)

  expect_false(isTRUE(all.equal(discrete$Fmsy, continuous$Fmsy, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(discrete$refpts_yield, continuous$refpts_yield, tolerance = 1e-6)))
})
