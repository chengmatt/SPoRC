# Tagged fish through two regions and two years, against a hand calculation written from the equations (1e-8),
# under each movement timing, with and without the state-space numbers at age factor. Mortality, movement and
# reporting differ by region and age, and recruits do not move.

library(SPoRC)
library(testthat)

# release_conv_tag_attr is called directly here, so messages_list is created the way the setup functions would
assign("messages_list", character(0), envir = .GlobalEnv)

# Setup ----------------------------------------------------------------------

tag_hand <- list(
  n_regions = 2,
  n_yrs = 2,
  n_ages = 3, # age 3 is the plus group
  released = c(400, 300, 200), # tags released in region 1 at ages 1 to 3
  init_tag_mort = 0.1, # tag induced mortality at release
  shed = 0.05, # tag shedding per year
  Fmort = c(0.3, 0.1), # fishing mortality by region
  reporting = c(0.6, 0.4), # reporting rate by region
  natmort = rbind(c(0.40, 0.25, 0.15), # region 1, ages 1 to 3
                  c(0.35, 0.20, 0.10)), # region 2
  # movement fractions [from, to], one matrix per age
  move = list(matrix(c(0.9, 0.1, 0.3, 0.7), nrow = 2, byrow = TRUE),
              matrix(c(0.6, 0.4, 0.2, 0.8), nrow = 2, byrow = TRUE),
              matrix(c(0.5, 0.5, 0.45, 0.55), nrow = 2, byrow = TRUE)),
  # movement rates [from, to] for continuous movement, rows summing to zero
  rate = list(matrix(c(-0.2, 0.2, 0.5, -0.5), nrow = 2, byrow = TRUE),
              matrix(c(-0.7, 0.7, 0.3, -0.3), nrow = 2, byrow = TRUE),
              matrix(c(-0.4, 0.4, 0.9, -0.9), nrow = 2, byrow = TRUE)),
  # state-space factor on numbers at age at the start of year 2 [region, age]
  naa_scalar = rbind(c(1.3, 0.8, 1.1),
                     c(0.7, 1.2, 0.9))
)

# what get_tagging_observation_model takes, for one cohort released in region 1 at the start of year 1
tag_hand_input <- function(move_timing, use_naa_scalar, released = tag_hand$released) {

  n_regions <- tag_hand$n_regions
  n_yrs <- tag_hand$n_yrs
  n_ages <- tag_hand$n_ages
  max_liberty <- 2

  Fmort <- array(0, dim = c(n_regions, n_yrs, 1, 1))
  natmort <- array(0, dim = c(1, n_regions, n_yrs, 1, n_ages, 1))
  Movement <- Mrate <- array(0, dim = c(1, n_regions, n_regions, n_yrs, 1, n_ages, 1))
  NAA_scalar <- array(1, dim = c(1, n_regions, n_yrs, 1, n_ages, 1))
  for(r in 1:n_regions) {
    Fmort[r,,1,1] <- tag_hand$Fmort[r]
    for(a in 1:n_ages) {
      natmort[1,r,,1,a,1] <- tag_hand$natmort[r,a]
      NAA_scalar[1,r,2,1,a,1] <- tag_hand$naa_scalar[r,a]
    } # end a loop
  } # end r loop
  for(a in 1:n_ages) {
    for(y in 1:n_yrs) {
      Movement[1,,,y,1,a,1] <- tag_hand$move[[a]]
      Mrate[1,,,y,1,a,1] <- tag_hand$rate[[a]]
    } # end y loop
  } # end a loop

  NAA_scalar_used <- NULL
  if(use_naa_scalar) NAA_scalar_used <- NAA_scalar

  list(
    n_fish_fleets = 1,
    n_regions = n_regions,
    n_conv_tag_cohorts = 1,
    n_yrs = n_yrs,
    n_seas = 1,
    n_pop = 1,
    n_ages = n_ages,
    n_sexes = 1,
    conv_tag_fish_reporting_blocks = array(1L, dim = c(n_regions, n_yrs, 1)),
    conv_tag_fish_reporting_pars = array(qlogis(tag_hand$reporting), dim = c(n_regions, 1, 1)),
    conv_tag_fish_reporting = array(0, dim = c(n_regions, n_yrs, 1)),
    conv_tag_release_indicator = matrix(c(1, 1, 1), nrow = 1), # region 1, year 1, season 1
    conv_tag_max_liberty = max_liberty,
    use_conv_fish_tagging = 1,
    Fmort = Fmort,
    fish_sel = array(1, dim = c(1, n_regions, n_yrs, 1, n_ages, 1, 1)), # fully selected
    ret_sel = array(1, dim = c(1, n_regions, n_yrs, 1, n_ages, 1, 1)), # fully retained
    dmr = array(1, dim = c(n_regions, n_yrs, 1, 1)),
    natmort = natmort,
    seasdur = 1,
    ln_conv_tag_shed = log(tag_hand$shed),
    conv_tag_t_tagging = 1, # released at the start of the season
    conv_tagged_fish = array(released, dim = c(1, 1, n_ages, 1)),
    conv_fish_tag_attr = "p_a_s", # releases known by age, so no apportionment
    conv_tag_release_platform = matrix(0, nrow = 1, ncol = 1),
    srv_sel = array(0, dim = c(1, n_regions, n_yrs, 1, n_ages, 1, 1)),
    NAA_bef = array(0, dim = c(1, n_regions, n_yrs, 1, n_ages, 1)),
    ln_init_conv_tag_mort = log(tag_hand$init_tag_mort),
    do_recruits_move = 0,
    Movement = Movement,
    conv_tag_fish_avail = array(0, dim = c(max_liberty + 1, 1, 1, 1, n_regions, n_ages, 1)),
    pred_conv_tag_fish_recap = array(0, dim = c(max_liberty, 1, 1, 1, n_regions, n_ages, 1, 1)),
    Mrate = Mrate,
    move_timing = move_timing,
    expm_nsub = 0,
    NAA_scalar = NAA_scalar_used
  )
}

# exp(A) and the integral of exp(A t) over one year, by their power series (A is 2 by 2 and small)
tag_hand_expm <- function(A) {
  term <- diag(nrow(A)) # A^k / k!
  expA <- integral <- matrix(0, nrow(A), ncol(A))
  for(k in 0:60) {
    expA <- expA + term
    integral <- integral + term / (k + 1) # the integral of t^k / k! over 0 to 1
    term <- term %*% A / (k + 1)
  } # end k loop
  list(expA = expA, integral = integral)
}

# tags available [years at liberty, region, age] and predicted recaptures, one age and region at a time
tag_hand_calc <- function(move_timing, use_naa_scalar, released = tag_hand$released) {

  n_regions <- tag_hand$n_regions
  n_ages <- tag_hand$n_ages
  avail <- array(0, dim = c(3, n_regions, n_ages))
  recap <- array(0, dim = c(2, n_regions, n_ages))
  avail[1,1,] <- released * exp(-tag_hand$init_tag_mort) # released in region 1, less tag induced mortality

  for(ry in 1:2) {

    year_end <- matrix(0, n_regions, n_ages) # tags alive at the end of the year, before ageing

    for(a in 1:n_ages) {

      Z <- tag_hand$natmort[,a] + tag_hand$Fmort + tag_hand$shed # total mortality by region, with shedding
      tags <- avail[ry,,a]
      moves <- a > 1 # recruits stay where they are

      # movement, then mortality and recaptures
      if(move_timing == 0) {
        if(moves) tags <- as.vector(tags %*% tag_hand$move[[a]])
        avail[ry,,a] <- tags # the model stores the tags after they move
        recap[ry,,a] <- tag_hand$reporting * tag_hand$Fmort / Z * tags * (1 - exp(-Z))
        year_end[,a] <- tags * exp(-Z)
      }

      # mortality and recaptures, then movement
      if(move_timing == 1) {
        recap[ry,,a] <- tag_hand$reporting * tag_hand$Fmort / Z * tags * (1 - exp(-Z))
        year_end[,a] <- tags * exp(-Z)
        if(moves) year_end[,a] <- as.vector(year_end[,a] %*% tag_hand$move[[a]])
      }

      # movement and mortality together, recaptures from the tags integrated over the year
      if(move_timing == 2) {
        A <- -diag(Z)
        if(moves) A <- tag_hand$rate[[a]] - diag(Z)
        series <- tag_hand_expm(A)
        recap[ry,,a] <- tag_hand$reporting * tag_hand$Fmort * as.vector(tags %*% series$integral)
        year_end[,a] <- as.vector(tags %*% series$expA)
      }

    } # end a loop

    # ageing, with the plus group taking the oldest true age and its own survivors
    avail[ry + 1,,2:n_ages] <- year_end[,1:(n_ages - 1)]
    avail[ry + 1,,n_ages] <- avail[ry + 1,,n_ages] + year_end[,n_ages]

    # tagged fish take the same state-space factor as the population at the start of year 2
    if(use_naa_scalar && ry == 1) avail[2,,2:n_ages] <- avail[2,,2:n_ages] * tag_hand$naa_scalar[,2:n_ages]

  } # end ry loop

  list(avail = avail, recap = recap)
}

# Tests ----------------------------------------------------------------------

test_that("tags available and predicted recaptures match the hand calculation under every movement timing", {

  for(move_timing in 0:2) {
    for(use_naa_scalar in c(FALSE, TRUE)) {

      model <- do.call(SPoRC:::get_tagging_observation_model, tag_hand_input(move_timing, use_naa_scalar))
      hand <- tag_hand_calc(move_timing, use_naa_scalar)
      label <- sprintf("move_timing %d, state-space factor %s", move_timing, use_naa_scalar)

      expect_equal(model$conv_tag_fish_avail[,1,1,1,,,1], hand$avail, tolerance = 1e-8, label = paste("tags available,", label))
      expect_equal(model$pred_conv_tag_fish_recap[,1,1,1,,,1,1], hand$recap, tolerance = 1e-8, label = paste("recaptures,", label))

    } # end use_naa_scalar loop
  } # end move_timing loop

})

test_that("tags released at age 1 stay in the release region for that year when recruits do not move", {

  # only age 1 released, in region 1, so a tag seen in region 2 during year 1 moved as a recruit
  for(move_timing in 0:2) {

    model <- do.call(SPoRC:::get_tagging_observation_model, tag_hand_input(move_timing, FALSE, released = c(400, 0, 0)))
    label <- sprintf("move_timing %d", move_timing)

    # recaptures in region 2 in the release year
    expect_equal(model$pred_conv_tag_fish_recap[1,1,1,1,2,1,1,1], 0, label = paste("region 2 recaptures at age 1,", label))
    expect_gt(model$pred_conv_tag_fish_recap[1,1,1,1,1,1,1,1], 10)

    # under movement before mortality the start of year 2 is stored after that year's movement, so it is not checked
    if(move_timing > 0) expect_equal(model$conv_tag_fish_avail[2,1,1,1,2,2,1], 0, label = paste("region 2 tags at the start of year 2,", label))

  } # end move_timing loop

})

test_that("the hand calculation separates the three movement timings and the state-space factor", {

  # if movement or the factor did nothing in this setup, the first test could not tell the timings apart
  timing_0 <- tag_hand_calc(0, FALSE)
  timing_1 <- tag_hand_calc(1, FALSE)
  timing_2 <- tag_hand_calc(2, FALSE)
  with_scalar <- tag_hand_calc(1, TRUE)

  expect_gt(max(abs(timing_0$recap - timing_1$recap)), 1)
  expect_gt(max(abs(timing_1$recap - timing_2$recap)), 1)
  expect_gt(max(abs(with_scalar$avail - timing_1$avail)), 1)

})

# Two Fleets -----------------------------------------------------------------

# one region, two ages, one year at liberty and two fleets that differ in every rate
tag_two_fleet <- list(
  released = c(400, 300), # tags released at ages 1 and 2
  natmort = c(0.30, 0.20),
  shed = 0.05,
  init_tag_mort = 0.1,
  Fmort = c(0.25, 0.10), # by fleet
  reporting = c(0.6, 0.3), # by fleet
  dmr = c(0.4, 0.8), # share of each fleet's discards that die
  fish_sel = cbind(c(0.5, 1.0), c(1.0, 0.6)), # [age, fleet]
  ret_sel = cbind(c(0.7, 1.0), c(0.4, 0.9)) # [age, fleet], the rest discarded
)

tag_two_fleet_input <- function(use_conv_fish_tagging) {

  fish_sel <- ret_sel <- array(0, dim = c(1, 1, 1, 1, 2, 1, 2))
  for(a in 1:2) {
    for(f in 1:2) {
      fish_sel[1,1,1,1,a,1,f] <- tag_two_fleet$fish_sel[a,f]
      ret_sel[1,1,1,1,a,1,f] <- tag_two_fleet$ret_sel[a,f]
    } # end f loop
  } # end a loop

  list(
    n_fish_fleets = 2,
    n_regions = 1,
    n_conv_tag_cohorts = 1,
    n_yrs = 1,
    n_seas = 1,
    n_pop = 1,
    n_ages = 2,
    n_sexes = 1,
    conv_tag_fish_reporting_blocks = array(1L, dim = c(1, 1, 2)),
    conv_tag_fish_reporting_pars = array(qlogis(tag_two_fleet$reporting), dim = c(1, 1, 2)),
    conv_tag_fish_reporting = array(0, dim = c(1, 1, 2)),
    conv_tag_release_indicator = matrix(c(1, 1, 1), nrow = 1), # region 1, year 1, season 1
    conv_tag_max_liberty = 1,
    use_conv_fish_tagging = use_conv_fish_tagging,
    Fmort = array(tag_two_fleet$Fmort, dim = c(1, 1, 1, 2)),
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    dmr = array(tag_two_fleet$dmr, dim = c(1, 1, 1, 2)),
    natmort = array(tag_two_fleet$natmort, dim = c(1, 1, 1, 1, 2, 1)),
    seasdur = 1,
    ln_conv_tag_shed = log(tag_two_fleet$shed),
    conv_tag_t_tagging = 1, # released at the start of the season
    conv_tagged_fish = array(tag_two_fleet$released, dim = c(1, 1, 2, 1)),
    conv_fish_tag_attr = "p_a_s",
    conv_tag_release_platform = matrix(0, nrow = 1, ncol = 1),
    srv_sel = array(0, dim = c(1, 1, 1, 1, 2, 1, 1)),
    NAA_bef = array(0, dim = c(1, 1, 1, 1, 2, 1)),
    ln_init_conv_tag_mort = log(tag_two_fleet$init_tag_mort),
    do_recruits_move = 0,
    Movement = array(1, dim = c(1, 1, 1, 1, 1, 2, 1)),
    conv_tag_fish_avail = array(0, dim = c(2, 1, 1, 1, 1, 2, 1)),
    pred_conv_tag_fish_recap = array(0, dim = c(1, 1, 1, 1, 1, 2, 1, 2))
  )
}

test_that("each fleet's recaptures use its own reporting, selectivity and retention, and its dead discards add to tag mortality", {

  for(use_conv_fish_tagging in list(c(1, 1), c(1, 0))) {

    model <- do.call(SPoRC:::get_tagging_observation_model, tag_two_fleet_input(use_conv_fish_tagging))
    label <- sprintf("fleets used for tags: %s", paste(use_conv_fish_tagging, collapse = ","))

    for(a in 1:2) {

      # total mortality on tagged fish: natural, shedding, and the fleets that recover tags
      Z <- tag_two_fleet$natmort[a] + tag_two_fleet$shed
      ret_F <- rep(0, 2)
      for(f in 1:2) {
        if(use_conv_fish_tagging[f] == 1) {
          ret_F[f] <- tag_two_fleet$Fmort[f] * tag_two_fleet$fish_sel[a,f] * tag_two_fleet$ret_sel[a,f]
          dead_discards <- tag_two_fleet$Fmort[f] * tag_two_fleet$fish_sel[a,f] * (1 - tag_two_fleet$ret_sel[a,f]) * tag_two_fleet$dmr[f]
          Z <- Z + ret_F[f] + dead_discards
        }
      } # end f loop

      tags <- tag_two_fleet$released[a] * exp(-tag_two_fleet$init_tag_mort)
      for(f in 1:2) {
        recaptures <- tag_two_fleet$reporting[f] * ret_F[f] / Z * tags * (1 - exp(-Z)) # zero for a fleet that recovers none
        expect_equal(model$pred_conv_tag_fish_recap[1,1,1,1,1,a,1,f], recaptures, tolerance = 1e-8,
                     label = sprintf("recaptures at age %d by fleet %d, %s", a, f, label))
      } # end f loop

      # tags left at the start of the next year, one age older, both ages in the plus group
      if(a == 1) survivors <- 0
      survivors <- survivors + tags * exp(-Z)

    } # end a loop

    expect_equal(model$conv_tag_fish_avail[2,1,1,1,1,2,1], survivors, tolerance = 1e-8, label = paste("tags surviving the year,", label))

  } # end use_conv_fish_tagging loop

})
