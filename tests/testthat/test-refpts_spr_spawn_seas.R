# Spawning-per-recruit when spawning happens in season 1 of a multi-season year.
#
# The seasonal reference point tests all set spawn_seas = 2. That leaves
# spawn_seas = 1 with n_seas > 1 uncovered, and it is the configuration the
# North Sea sandeel case study uses (vignettes/ah_north_sea_sandeel_case_study.Rmd:
# n_pop 1, n_regions 1, n_sexes 1, n_seas 2, spawn_seas 1).
#
# The gap matters because of how the solver is arranged. Ages up to n_ages - 2
# capture spawning biomass inside the season loop, while the penultimate age and
# the plus group are handled in separate blocks that advance a start-of-year
# abundance forward to the spawning season under `if(spawn_seas > 1)`. When
# spawn_seas is 1 that advancement is skipped, so the two halves of the
# calculation meet in a way no test exercised.
#
# These check the solver against a cohort tracked forward season by season, a
# different arrangement of the same quantity: one continuous walk rather than a
# loop plus two closed-form tails. Two conventions are shared rather than
# rederived, because they are choices the model makes rather than consequences:
# a season's F is applied whole rather than scaled by season duration, and the
# plus group accumulates on annual rates.

spr_data <- function(spawn_seas, rec_seas_prop, n_seas = 2L, seasdur = NULL,
                     t_spawn = 0, n_ages = 5L, M = 0.8, sel = NULL, dmr = 0) {
  if(is.null(seasdur)) seasdur <- rep(1 / n_seas, n_seas)
  if(is.null(sel)) sel <- c(0, 0.4, 1, 1, 0.7)[seq_len(n_ages)]
  n_flt <- 1L
  list(
    n_pop = 1L, n_ages = n_ages, n_seas = n_seas, n_fish_fleets = n_flt,
    seasdur = seasdur, spawn_seas = spawn_seas, t_spawn = t_spawn,
    F_fract_flt = array(1, dim = c(n_seas, n_flt)),
    fish_sel = array(rep(sel, each = n_seas), dim = c(1, n_seas, n_ages, n_flt)),
    ret_sel = array(1, dim = c(1, n_seas, n_ages, n_flt)),
    dmr = array(dmr, dim = c(n_seas, n_flt)),
    natmort = array(M, dim = c(1, n_ages)),
    WAA = array(rep(seq(0.005, 0.025, length.out = n_ages), each = n_seas),
                dim = c(1, n_seas, n_ages)),
    MatAA = array(rep(c(0, 0.5, 1, 1, 1)[seq_len(n_ages)], each = n_seas),
                  dim = c(1, n_seas, n_ages)),
    sex_ratio_f = 1, rec_seas_prop = array(rec_seas_prop, dim = c(1, n_seas)),
    stray_rate = 0, natal_region = 1L, n_pop_in_region = 1, SPR_x = 0.4)
}

# Spawning biomass per recruit, walked forward one season at a time.
spr_oracle <- function(d, f) {
  M <- d$natmort[1, ]; sel <- d$fish_sel[1, 1, , 1]
  ss <- d$spawn_seas; ts <- d$t_spawn; sd_ <- d$seasdur
  n_ages <- d$n_ages; n_seas <- d$n_seas

  N <- numeric(n_ages); SB <- numeric(n_ages)
  start_of_year <- numeric(n_ages)
  N[1] <- d$rec_seas_prop[1, 1] * d$sex_ratio_f[1]
  start_of_year[1] <- N[1]

  for(a in 1:(n_ages - 1)) {
    for(s in 1:n_seas) {
      f_seas <- f * sel[a] * (1 + (1 - 1) * d$dmr[s, 1])
      z_seas <- f_seas + M[a] * sd_[s]
      if(s > 1 && a == 1) N[a] <- N[a] + d$rec_seas_prop[1, s] * d$sex_ratio_f[1]
      if(s == ss) SB[a] <- N[a] * d$WAA[1, ss, a] * d$MatAA[1, ss, a] * exp(-ts * z_seas)
      if(s < n_seas) {
        N[a] <- N[a] * exp(-z_seas)
      } else {
        N[a + 1] <- N[a] * exp(-z_seas)
        start_of_year[a + 1] <- N[a + 1]
      }
    } # end s loop
  } # end a loop

  # the plus group accumulates from the penultimate age at its start-of-year
  # abundance, on annual rates, then is carried to the spawning season
  z_penult <- M[n_ages - 1] + n_seas * f * sel[n_ages - 1]
  z_plus <- M[n_ages] + n_seas * f * sel[n_ages]
  n_plus <- start_of_year[n_ages - 1] * exp(-z_penult) / (1 - exp(-z_plus))
  if(ss > 1)
    for(s in 1:(ss - 1))
      n_plus <- n_plus * exp(-(M[n_ages] * sd_[s] + f * sel[n_ages]))
  SB[n_ages] <- n_plus * d$WAA[1, ss, n_ages] * d$MatAA[1, ss, n_ages] *
    exp(-ts * (M[n_ages] * sd_[ss] + f * sel[n_ages]))

  sum(SB)
}

spr_solver <- function(d, f) {
  obj <- RTMB::MakeADFun(function(p) single_region_SPR(p, d),
                         list(log_F_x = log(max(f, 1e-12))), silent = TRUE)
  rep <- obj$report(obj$par)
  if(f < 1e-9) rep$SB0 else rep$SB
}


test_that("spawning in season 1 of a multi-season year gives the per-recruit biomass of a cohort walked forward", {
  # The uncovered case, and the one the sandeel case study runs.
  for(rec in list(c(0, 1), c(1, 0), c(0.3, 0.7))) {
    d <- spr_data(spawn_seas = 1L, rec_seas_prop = rec)
    for(f in c(0, 0.1, 0.4)) {
      expect_equal(spr_solver(d, f), spr_oracle(d, f), tolerance = 1e-10,
                   label = sprintf("SBPR at spawn_seas 1, rec (%s), F = %g",
                                   paste(rec, collapse = ", "), f))
    } # end f loop
  } # end rec loop
})


test_that("the agreement holds for every season spawning can fall in", {
  # spawn_seas = 2 was the only seasonal case with coverage. Running all of them
  # is what distinguishes a solver that is right from one whose two halves happen
  # to cancel at a particular spawning season.
  for(n_seas in c(1L, 2L, 4L)) {
    for(ss in seq_len(n_seas)) {
      rec <- rep(0, n_seas); rec[n_seas] <- 1
      d <- spr_data(spawn_seas = ss, rec_seas_prop = rec, n_seas = n_seas)
      for(f in c(0, 0.25)) {
        expect_equal(spr_solver(d, f), spr_oracle(d, f), tolerance = 1e-10,
                     label = sprintf("SBPR at n_seas = %d, spawn_seas = %d, F = %g",
                                     n_seas, ss, f))
      } # end f loop
    } # end ss loop
  } # end n_seas loop
})


test_that("unequal season durations and spawning inside the season agree", {
  # t_spawn moves spawning to part way through its season, and the penultimate
  # age and plus group apply that correction through a different expression than
  # the season loop does. Uneven durations stop a bug that scales by 1 / n_seas
  # from passing.
  for(ss in 1:2) {
    for(ts in c(0, 0.5, 1)) {
      d <- spr_data(spawn_seas = ss, rec_seas_prop = c(0, 1),
                    seasdur = c(0.3, 0.7), t_spawn = ts)
      for(f in c(0, 0.3)) {
        expect_equal(spr_solver(d, f), spr_oracle(d, f), tolerance = 1e-10,
                     label = sprintf("SBPR at spawn_seas = %d, t_spawn = %g, F = %g",
                                     ss, ts, f))
      } # end f loop
    } # end ts loop
  } # end ss loop
})


test_that("spawning season changes the answer", {
  # Without this the agreement above could be read as both calculations being
  # insensitive to spawn_seas. Spawning earlier in the year means less mortality
  # has accrued, so per-recruit biomass is strictly higher.
  d1 <- spr_data(spawn_seas = 1L, rec_seas_prop = c(1, 0))
  d2 <- spr_data(spawn_seas = 2L, rec_seas_prop = c(1, 0))

  expect_gt(spr_solver(d1, 0), spr_solver(d2, 0))
  expect_false(isTRUE(all.equal(spr_solver(d1, 0.3), spr_solver(d2, 0.3))))
})


test_that("the plus group carries a real share of per-recruit biomass", {
  # The one place the two arrangements differed was the plus group, and it only
  # showed up because the plus group is a large share of the total here. A
  # fixture where it were negligible would pass while testing nothing.
  d <- spr_data(spawn_seas = 1L, rec_seas_prop = c(0, 1))
  obj <- RTMB::MakeADFun(function(p) single_region_SPR(p, d),
                         list(log_F_x = log(1e-12)), silent = TRUE)
  sb_age <- obj$report(obj$par)$SB_age[1, 1, ]

  expect_gt(sb_age[length(sb_age)] / sum(sb_age), 0.05)
})
