# Spawning biomass at spawning time with three populations that home and stray, against a hand calculation from the
# documented equations (1e-10), and the numbers at the spawning point under each movement timing. Two regions,
# two sexes, and inputs that differ by population, region, age and sex.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

# three populations, two of them born in region 1, so a region's strays are shared between two populations
spawn_hand <- list(
  n_pop = 3,
  n_regions = 2,
  n_ages = 3,
  n_sexes = 2,
  natal_region = c(1, 1, 2),
  stray_rate = c(0.10, 0.25, 0.40), # share of each population's spawning biomass that counts elsewhere
  t_spawn = 0.3, # spawning 30 pct of the way into the year
  seasdur = 1
)

# an array [pop, region, 1, 1, age, sex] with a different value in every cell
spawn_hand_array <- function(base, by_pop, by_region, by_age, by_sex) {
  out <- array(0, dim = c(spawn_hand$n_pop, spawn_hand$n_regions, 1, 1, spawn_hand$n_ages, spawn_hand$n_sexes))
  for(p in 1:spawn_hand$n_pop) {
    for(r in 1:spawn_hand$n_regions) {
      for(a in 1:spawn_hand$n_ages) {
        for(s in 1:spawn_hand$n_sexes) {
          out[p,r,1,1,a,s] <- base * by_pop[p] * by_region[r] * by_age[a] * by_sex[s]
        } # end s loop
      } # end a loop
    } # end r loop
  } # end p loop
  out
}

spawn_hand$NAA <- spawn_hand_array(1000, c(1, 0.6, 1.5), c(1, 0.4), c(1, 0.7, 0.9), c(1, 1.1))
spawn_hand$NAA0 <- spawn_hand_array(1600, c(1, 0.7, 1.4), c(1, 0.5), c(1, 0.8, 1.3), c(1, 1.05))
spawn_hand$WAA <- spawn_hand_array(1, c(1, 1.1, 0.9), c(1, 1.2), c(0.6, 1.4, 2.3), c(1, 0.85))
spawn_hand$MatAA <- spawn_hand_array(1, c(1, 0.9, 1), c(1, 0.95), c(0.1, 0.6, 1), c(1, 1))
spawn_hand$natmort <- spawn_hand_array(0.2, c(1, 1.2, 0.8), c(1, 1.3), c(1.5, 1, 0.7), c(1, 1.1)) # per year
spawn_hand$ZAA <- spawn_hand$natmort + spawn_hand_array(0.15, c(1, 1, 1), c(1, 0.5), c(0.3, 1, 1), c(1, 0.9)) # with fishing

# natal homing [pop, from, to, 1, age, sex]: a share of the fish outside the natal region return to it to spawn
spawn_hand$spawn_move <- array(0, dim = c(spawn_hand$n_pop, 2, 2, 1, spawn_hand$n_ages, spawn_hand$n_sexes))
for(p in 1:spawn_hand$n_pop) {
  for(a in 1:spawn_hand$n_ages) {
    for(s in 1:spawn_hand$n_sexes) {
      home <- spawn_hand$natal_region[p]
      away <- 3 - home
      returning <- 0.2 * a + 0.05 * s # older fish return more
      spawn_hand$spawn_move[p,home,home,1,a,s] <- 1
      spawn_hand$spawn_move[p,away,home,1,a,s] <- returning
      spawn_hand$spawn_move[p,away,away,1,a,s] <- 1 - returning
    } # end s loop
  } # end a loop
} # end p loop

spawn_hand_model <- function(n_sexes = spawn_hand$n_sexes) {
  sexes <- 1:n_sexes
  SPoRC:::biom_at_spawn(
    NAA_s = spawn_hand$NAA[,,,,,sexes, drop = FALSE],
    NAA0_s = spawn_hand$NAA0[,,,,,sexes, drop = FALSE],
    WAA_s = spawn_hand$WAA[,,,,,sexes, drop = FALSE],
    MatAA_s = spawn_hand$MatAA[,,,,,sexes, drop = FALSE],
    ZAA_s = spawn_hand$ZAA[,,,,,sexes, drop = FALSE],
    natmort_s = spawn_hand$natmort[,,,,,sexes, drop = FALSE],
    spawn_move_s = spawn_hand$spawn_move[,,,,,sexes, drop = FALSE],
    stray_rate_y = spawn_hand$stray_rate,
    t_spawn = spawn_hand$t_spawn,
    seasdur_seas = spawn_hand$seasdur,
    n_seas = 1,
    n_pop = spawn_hand$n_pop,
    n_regions = spawn_hand$n_regions,
    n_ages = spawn_hand$n_ages,
    n_sexes = n_sexes,
    natal_region = spawn_hand$natal_region,
    move_timing = 0
  )
}

# the same quantities one population, region, age and sex at a time
spawn_hand_calc <- function(n_sexes = spawn_hand$n_sexes) {

  n_pop <- spawn_hand$n_pop
  n_regions <- spawn_hand$n_regions
  n_ages <- spawn_hand$n_ages
  female_share <- 1
  if(n_sexes == 1) female_share <- 0.5 # one sex holds both, so half are female

  SSB <- Dynamic_SSB0 <- Total_Biom <- matrix(0, n_pop, n_regions)

  for(p in 1:n_pop) {
    for(a in 1:n_ages) {
      for(s in 1:n_sexes) {
        for(r_to in 1:n_regions) {

          # fish in each region after natal homing
          homed <- homed0 <- 0
          for(r_from in 1:n_regions) {
            homed <- homed + spawn_hand$NAA[p,r_from,1,1,a,s] * spawn_hand$spawn_move[p,r_from,r_to,1,a,s]
            homed0 <- homed0 + spawn_hand$NAA0[p,r_from,1,1,a,s] * spawn_hand$spawn_move[p,r_from,r_to,1,a,s]
          } # end r_from loop

          # mortality up to spawning, where the fish now are
          at_spawning <- homed * exp(-spawn_hand$t_spawn * spawn_hand$ZAA[p,r_to,1,1,a,s])
          at_spawning0 <- homed0 * exp(-spawn_hand$t_spawn * spawn_hand$natmort[p,r_to,1,1,a,s] * spawn_hand$seasdur)

          Total_Biom[p,r_to] <- Total_Biom[p,r_to] + at_spawning * spawn_hand$WAA[p,r_to,1,1,a,s]

          # spawning biomass is females, the first sex
          if(s == 1) {
            SSB[p,r_to] <- SSB[p,r_to] + female_share * at_spawning * spawn_hand$WAA[p,r_to,1,1,a,s] * spawn_hand$MatAA[p,r_to,1,1,a,s]
            Dynamic_SSB0[p,r_to] <- Dynamic_SSB0[p,r_to] + female_share * at_spawning0 * spawn_hand$WAA[p,r_to,1,1,a,s] * spawn_hand$MatAA[p,r_to,1,1,a,s]
          }

        } # end r_to loop
      } # end s loop
    } # end a loop
  } # end p loop

  # effective spawning biomass: a population's own in its natal region, plus the strays of the others there,
  # shared among the populations born in that region
  eff_SSB <- rep(0, n_pop)
  for(p in 1:n_pop) {
    home <- spawn_hand$natal_region[p]
    n_born_here <- sum(spawn_hand$natal_region == home)
    eff_SSB[p] <- SSB[p,home]
    for(p_other in 1:n_pop) {
      if(p_other != p) eff_SSB[p] <- eff_SSB[p] + spawn_hand$stray_rate[p_other] / n_born_here * SSB[p_other,home]
    } # end p_other loop
  } # end p loop

  list(SSB = SSB, Dynamic_SSB0 = Dynamic_SSB0, Total_Biom = Total_Biom, eff_SSB = eff_SSB)
}

# Tests ----------------------------------------------------------------------

test_that("spawning biomass, dynamic B0 and total biomass match the hand calculation after natal homing", {

  for(n_sexes in 1:2) {

    model <- spawn_hand_model(n_sexes)
    hand <- spawn_hand_calc(n_sexes)
    label <- sprintf("%d sex", n_sexes)

    expect_equal(as.numeric(model$SSB_y), as.numeric(hand$SSB), tolerance = 1e-10, label = paste("spawning biomass,", label))
    expect_equal(as.numeric(model$Dynamic_SSB0_y), as.numeric(hand$Dynamic_SSB0), tolerance = 1e-10, label = paste("dynamic B0,", label))
    expect_equal(as.numeric(model$Total_Biom_y), as.numeric(hand$Total_Biom), tolerance = 1e-10, label = paste("total biomass,", label))

  } # end n_sexes loop

})

test_that("effective spawning biomass adds the other populations' strays, shared among those born in the region", {

  model <- spawn_hand_model()
  hand <- spawn_hand_calc()

  expect_equal(as.numeric(model$eff_SSB_y), hand$eff_SSB, tolerance = 1e-10)

  # strays add to every population here, so a model that left them out would differ
  for(p in 1:spawn_hand$n_pop) expect_gt(hand$eff_SSB[p], 1.01 * hand$SSB[p,spawn_hand$natal_region[p]])

})

test_that("numbers at the spawning point follow each movement timing", {

  numbers <- c(500, 300) # start of season, by region
  Z <- c(0.45, 0.30) # total mortality over the season
  t_spawn <- 0.3
  move <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE) # movement fractions [from, to]
  rate <- matrix(c(-0.6, 0.6, 0.25, -0.25), nrow = 2, byrow = TRUE) # movement rates [from, to]

  # movement first: fish move, then die at the destination up to spawning
  moved <- c(numbers[1] * move[1,1] + numbers[2] * move[2,1],
             numbers[1] * move[1,2] + numbers[2] * move[2,2])
  expect_equal(SPoRC:::spawn_state(N = numbers, Move = move, Z = Z, Q = rate, dur = 1, t_spawn = t_spawn, move_timing = 0),
               moved * exp(-t_spawn * Z), tolerance = 1e-10)

  # movement at the end of the season: fish die where they started up to spawning
  expect_equal(SPoRC:::spawn_state(N = numbers, Move = move, Z = Z, Q = rate, dur = 1, t_spawn = t_spawn, move_timing = 1),
               numbers * exp(-t_spawn * Z), tolerance = 1e-10)

  # movement and mortality together for t_spawn of the season, by the power series of the matrix exponential
  generator <- (rate - diag(Z)) * t_spawn
  term <- diag(2)
  transition <- matrix(0, 2, 2)
  for(k in 0:40) {
    transition <- transition + term
    term <- term %*% generator / (k + 1)
  } # end k loop
  expect_equal(SPoRC:::spawn_state(N = numbers, Move = move, Z = Z, Q = rate, dur = 1, t_spawn = t_spawn, move_timing = 2),
               as.vector(numbers %*% transition), tolerance = 1e-8)

})

test_that("natal homing at the spawning point comes before the mortality under movement first, and after it otherwise", {

  numbers <- c(500, 300) # start of season, by region
  Z <- c(0.45, 0.30) # total mortality over the season
  t_spawn <- 0.3
  move <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE) # movement fractions [from, to]
  rate <- matrix(c(-0.6, 0.6, 0.25, -0.25), nrow = 2, byrow = TRUE) # movement rates [from, to]
  homing <- matrix(c(1, 0, 0.6, 0.4), nrow = 2, byrow = TRUE) # 60 pct of the fish in region 2 return to region 1 to spawn

  # movement first: fish move, home, then die where they spawn
  moved <- c(numbers[1] * move[1,1] + numbers[2] * move[2,1],
             numbers[1] * move[1,2] + numbers[2] * move[2,2])
  homed <- c(moved[1] * homing[1,1] + moved[2] * homing[2,1],
             moved[1] * homing[1,2] + moved[2] * homing[2,2])
  expect_equal(SPoRC:::spawn_state(N = numbers, Move = move, Z = Z, Q = rate, dur = 1, t_spawn = t_spawn, move_timing = 0, spawn_move = homing),
               homed * exp(-t_spawn * Z), tolerance = 1e-10)

  # movement at the end of the season: fish die where they started, then home
  survivors <- numbers * exp(-t_spawn * Z)
  homed <- c(survivors[1] * homing[1,1] + survivors[2] * homing[2,1],
             survivors[1] * homing[1,2] + survivors[2] * homing[2,2])
  expect_equal(SPoRC:::spawn_state(N = numbers, Move = move, Z = Z, Q = rate, dur = 1, t_spawn = t_spawn, move_timing = 1, spawn_move = homing),
               homed, tolerance = 1e-10)

  # continuous movement: the numbers at the spawning point, then home
  at_spawning <- SPoRC:::spawn_state(N = numbers, Move = move, Z = Z, Q = rate, dur = 1, t_spawn = t_spawn, move_timing = 2)
  homed <- c(at_spawning[1] * homing[1,1] + at_spawning[2] * homing[2,1],
             at_spawning[1] * homing[1,2] + at_spawning[2] * homing[2,2])
  expect_equal(SPoRC:::spawn_state(N = numbers, Move = move, Z = Z, Q = rate, dur = 1, t_spawn = t_spawn, move_timing = 2, spawn_move = homing),
               homed, tolerance = 1e-10)

})

test_that("recruits that do not move spawn where they are under continuous movement", {

  n_regions <- 2
  n_ages <- 3
  t_spawn <- 0.3
  numbers <- rbind(c(500, 300, 200), c(250, 400, 100)) # start of season [region, age]
  Z <- rbind(c(0.50, 0.35, 0.25), c(0.40, 0.30, 0.20)) # total mortality over the season [region, age]
  natmort <- rbind(c(0.30, 0.20, 0.15), c(0.25, 0.18, 0.12)) # per year
  WAA <- c(0.6, 1.4, 2.3) # weight at age (kg)
  MatAA <- c(0.4, 0.8, 1) # recruits are partly mature here, so where they spawn enters the answer
  rate <- matrix(c(-0.6, 0.6, 0.25, -0.25), nrow = 2, byrow = TRUE) # movement rates [from, to], every age

  cells <- c(1, n_regions, 1, 1, n_ages, 1) # [pop, region, 1, 1, age, sex]
  Mrate <- Movement <- array(0, dim = c(1, n_regions, n_regions, 1, 1, n_ages, 1))
  WAA_s <- MatAA_s <- array(0, dim = cells)
  for(a in 1:n_ages) {
    Mrate[1,,,1,1,a,1] <- rate
    Movement[1,,,1,1,a,1] <- diag(n_regions) # not read under continuous movement
    WAA_s[1,,1,1,a,1] <- WAA[a]
    MatAA_s[1,,1,1,a,1] <- MatAA[a]
  } # end a loop

  spawning <- function(do_recruits_move) {
    SPoRC:::biom_at_spawn(
      NAA_s = array(numbers, dim = cells),
      NAA0_s = array(numbers, dim = cells),
      WAA_s = WAA_s,
      MatAA_s = MatAA_s,
      ZAA_s = array(Z, dim = cells),
      natmort_s = array(natmort, dim = cells),
      spawn_move_s = NULL,
      stray_rate_y = NULL,
      t_spawn = t_spawn,
      seasdur_seas = 1,
      n_seas = 1,
      n_pop = 1,
      n_regions = n_regions,
      n_ages = n_ages,
      n_sexes = 1,
      natal_region = 1,
      Movement_s = Movement,
      Mrate_s = Mrate,
      move_timing = 2,
      do_recruits_move = do_recruits_move
    )
  }

  # numbers at the spawning point by region for one age: movement and mortality together for t_spawn of the season,
  # by the power series of the matrix exponential, or mortality alone when the age does not move
  at_spawning <- function(a, mortality, moves) {
    if(!moves) return(numbers[,a] * exp(-t_spawn * mortality[,a]))
    generator <- (rate - diag(mortality[,a])) * t_spawn
    term <- diag(n_regions)
    transition <- matrix(0, n_regions, n_regions)
    for(k in 0:40) {
      transition <- transition + term
      term <- term %*% generator / (k + 1)
    } # end k loop
    as.vector(numbers[,a] %*% transition)
  }

  for(do_recruits_move in 0:1) {

    SSB <- Dynamic_SSB0 <- rep(0, n_regions)
    for(a in 1:n_ages) {
      moves <- do_recruits_move == 1 || a > 1
      SSB <- SSB + 0.5 * at_spawning(a, Z, moves) * WAA[a] * MatAA[a] # one sex, so half are female
      Dynamic_SSB0 <- Dynamic_SSB0 + 0.5 * at_spawning(a, natmort, moves) * WAA[a] * MatAA[a]
    } # end a loop

    model <- spawning(do_recruits_move)
    label <- sprintf("do_recruits_move %d", do_recruits_move)
    expect_equal(as.numeric(model$SSB_y), SSB, tolerance = 1e-8, label = paste("spawning biomass,", label))
    expect_equal(as.numeric(model$Dynamic_SSB0_y), Dynamic_SSB0, tolerance = 1e-8, label = paste("dynamic B0,", label))

  } # end do_recruits_move loop

  # whether recruits move changes spawning biomass by region here
  expect_gt(max(abs(as.numeric(spawning(0)$SSB_y) / as.numeric(spawning(1)$SSB_y) - 1)), 0.005)

})
