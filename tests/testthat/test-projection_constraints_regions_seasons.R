# catch_cap and catch_stability in a projection with two regions and two seasons, against what each is defined to do
# (1e-4, the solver's tolerance). Half of region 1 moves to region 2 each season, so lowering one region's F raises
# the other's catch.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

constraint_hand <- list(
  n_proj_yrs = 5,
  n_regions = 2,
  n_seas = 2,
  n_ages = 4,
  terminal_F = array(c(0.20, 0.10, 0.12, 0.18), dim = c(2, 2, 1)), # year 1 F [region, season, fleet]
  f_input = matrix(c(0.6, 0.5, 0.6, 0.5, 0.7, # f_input[r, y] is region r's annual F in year y + 1
                     0.3, 0.3, 0.35, 0.3, 0.3), nrow = 2, byrow = TRUE)
)

constraint_hand_run <- function(F_constraint = NULL) {

  n_proj_yrs <- constraint_hand$n_proj_yrs
  n_regions <- constraint_hand$n_regions
  n_seas <- constraint_hand$n_seas
  n_ages <- constraint_hand$n_ages
  biol_dim <- c(1, n_regions, n_proj_yrs, n_seas, n_ages, 1)

  natmort <- WAA <- MatAA <- array(0, dim = biol_dim)
  WAA_fish <- fish_sel <- array(0, dim = c(biol_dim, 1))
  for(a in 1:n_ages) {
    natmort[,,,,a,] <- c(0.3, 0.25, 0.2, 0.15)[a] # per year
    WAA[,,,,a,] <- c(0.5, 1.1, 1.8, 2.6)[a] # weight at age (kg)
    MatAA[,,,,a,] <- c(0, 0.3, 0.8, 1)[a]
    WAA_fish[,,,,a,,] <- c(0.5, 1.1, 1.8, 2.6)[a]
    fish_sel[,,,,a,,] <- c(0.3, 0.7, 1, 1)[a]
  } # end a loop

  # movement fractions [from, to]: half of region 1 moves to region 2 each season, a tenth comes back
  Movement <- array(0, dim = c(1, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, 1))
  for(y in 1:n_proj_yrs) for(seas in 1:n_seas) for(a in 1:n_ages) Movement[1,,,y,seas,a,1] <- matrix(c(0.5, 0.5, 0.1, 0.9), nrow = 2, byrow = TRUE)

  # numbers at the start of each season of year 1 [pop, region, season, age, sex]
  terminal_NAA <- array(0, dim = c(1, n_regions, n_seas, n_ages, 1))
  terminal_NAA[1,1,1,,1] <- c(1000, 800, 600, 700)
  terminal_NAA[1,2,1,,1] <- c(400, 300, 250, 300)
  terminal_NAA[1,1,2,,1] <- c(700, 560, 430, 520)
  terminal_NAA[1,2,2,,1] <- c(600, 480, 380, 430)

  Do_Population_Projection(
    n_proj_yrs = n_proj_yrs,
    n_pop = 1,
    n_regions = n_regions,
    n_ages = n_ages,
    n_sexes = 1,
    n_fish_fleets = 1,
    sexratio = array(1, dim = c(1, n_regions, n_proj_yrs, 1)),
    do_recruits_move = 1,
    recruitment = array(c(900, 1100, 400, 500), dim = c(1, n_regions, 2)),
    terminal_NAA = terminal_NAA,
    terminal_NAA0 = terminal_NAA,
    terminal_F = constraint_hand$terminal_F,
    natmort = natmort,
    WAA = WAA,
    WAA_fish = WAA_fish,
    MatAA = MatAA,
    fish_sel = fish_sel,
    Movement = Movement,
    f_ref_pt = constraint_hand$f_input,
    b_ref_pt = array(1, dim = c(1, n_regions, n_proj_yrs)),
    recruitment_opt = "mean_rec",
    fmort_opt = "Input",
    t_spawn = 0,
    n_seas = n_seas,
    spawn_seas = 1,
    F_constraint = F_constraint
  )
}

# retained catch by region and year, summed over seasons
constraint_hand_catch <- function(out) {
  catch <- matrix(0, constraint_hand$n_regions, constraint_hand$n_proj_yrs)
  for(r in 1:constraint_hand$n_regions) for(y in 1:constraint_hand$n_proj_yrs) catch[r,y] <- sum(out$proj_Catch[,r,y,,])
  catch
}

# the projection with no constraint, which every test compares against
constraint_hand_free <- constraint_hand_run()
constraint_hand_free_catch <- constraint_hand_catch(constraint_hand_free)

# Tests ----------------------------------------------------------------------

test_that("a cap per region is met in both regions when lowering one raises the other", {

  # region 1 is cut to half its catch in year 2, and region 2's cap sits just above what the input F takes there
  cap <- matrix(1e9, constraint_hand$n_regions, constraint_hand$n_proj_yrs)
  cap[1,2] <- 0.5 * constraint_hand_free_catch[1,2]
  cap[2,2] <- 1.003 * constraint_hand_free_catch[2,2]

  out <- constraint_hand_run(F_constraint = catch_cap(cap))
  catch <- constraint_hand_catch(out)

  expect_equal(catch[1,2], cap[1,2], tolerance = 1e-4)
  expect_equal(catch[2,2], cap[2,2], tolerance = 1e-4)

  # region 2 was under its cap at the input F, so its F only moved because region 1's fish arrived
  expect_lt(out$proj_F[2,2], 0.999 * out$proj_F_rule[2,2])
  expect_equal(out$proj_F_rule[,2], constraint_hand$f_input[,1], tolerance = 1e-12)

  # year 1 is not constrained, and the later years' caps are far above the catch
  expect_equal(catch[,1], constraint_hand_free_catch[,1], tolerance = 1e-12)
  expect_equal(out$proj_F[,3:5], out$proj_F_rule[,3:5], tolerance = 1e-12)

})

test_that("a constrained F keeps the terminal year's split across seasons", {

  cap <- matrix(1e9, constraint_hand$n_regions, constraint_hand$n_proj_yrs)
  cap[1,2] <- 0.5 * constraint_hand_free_catch[1,2]
  out <- constraint_hand_run(F_constraint = catch_cap(cap))

  # region 1 is solved to its cap while region 2, never over its own, keeps the input F in both seasons
  expect_equal(constraint_hand_catch(out)[1,2], cap[1,2], tolerance = 1e-4)
  expect_equal(out$proj_F[2,2], constraint_hand$f_input[2,1], tolerance = 1e-12)

  for(r in 1:constraint_hand$n_regions) {
    terminal_share <- constraint_hand$terminal_F[r,,1] / sum(constraint_hand$terminal_F[r,,1])
    for(y in 2:constraint_hand$n_proj_yrs) {
      expect_equal(out$proj_F_seas[r,y,] / out$proj_F[r,y], terminal_share, tolerance = 1e-8,
                   label = sprintf("seasonal share of F, region %d year %d", r, y))
    } # end y loop
  } # end r loop

})

test_that("a global cap is met in total and shared among regions as the input F would have", {

  cap <- 0.6 * sum(constraint_hand_free_catch[,2])
  out <- constraint_hand_run(F_constraint = catch_cap(c(1e9, cap, 1e9, 1e9, 1e9), scope = "global"))
  catch <- constraint_hand_catch(out)

  expect_equal(sum(catch[,2]), cap, tolerance = 1e-4)
  expect_equal(catch[,2] / sum(catch[,2]), constraint_hand_free_catch[,2] / sum(constraint_hand_free_catch[,2]), tolerance = 1e-4)

})

test_that("catch stability keeps each region, or the total, within the band of last year's catch", {

  down <- 0.10
  up <- 0.05

  # the input F alone leaves the band in every year, so each year is moved to an edge
  free_change <- constraint_hand_free_catch[,2:5] / constraint_hand_free_catch[,1:4]
  expect_true(any(free_change < 1 - down) && any(free_change > 1 + up))

  # by region
  by_region <- constraint_hand_catch(constraint_hand_run(F_constraint = catch_stability(down = down, up = up)))
  change <- by_region[,2:5] / by_region[,1:4]
  expect_true(all(change > (1 - down) * (1 - 1e-4)))
  expect_true(all(change < (1 + up) * (1 + 1e-4)))
  expect_equal(change[1,1], 1 - down, tolerance = 1e-4) # region 1 would have fallen 16 pct in year 2
  expect_equal(change[2,1], 1 + up, tolerance = 1e-4) # region 2 would have doubled

  # on the total
  global <- constraint_hand_catch(constraint_hand_run(F_constraint = catch_stability(down = down, up = up, scope = "global")))
  total_change <- colSums(global)[2:5] / colSums(global)[1:4]
  expect_true(all(total_change > (1 - down) * (1 - 1e-4)))
  expect_true(all(total_change < (1 + up) * (1 + 1e-4)))

})

test_that("a constraint is given last year's catch and spawning biomass", {

  seen <- list()
  record_state <- function(F, state) {
    seen[[length(seen) + 1]] <<- list(y = state$y, SSB = state$SSB, last_catch = state$last_catch, F_rule = state$F_rule)
    F # left as the input F
  }
  out <- constraint_hand_run(F_constraint = record_state)
  catch <- constraint_hand_catch(out)

  expect_equal(length(seen), constraint_hand$n_proj_yrs - 1) # years 2 to 5
  for(i in seq_along(seen)) {
    y <- seen[[i]]$y
    expect_equal(as.numeric(seen[[i]]$SSB), as.numeric(out$proj_SSB[1,,y - 1]), tolerance = 1e-10)
    expect_equal(as.numeric(seen[[i]]$last_catch), catch[,y - 1], tolerance = 1e-10)
    expect_equal(as.numeric(seen[[i]]$F_rule), constraint_hand$f_input[,y - 1], tolerance = 1e-12)
  } # end i loop

  # spawning biomass differs between years here, so last year's cannot pass for this year's
  expect_gt(abs(out$proj_SSB[1,1,2] / out$proj_SSB[1,1,1] - 1), 0.05)

})

test_that("catch stability checks the regions again when holding one moves another out of its band", {

  # how many times the constraint solved for F in each year
  solves <- list()
  count_solves <- function(constraint) {
    function(F, state) {
      n_solves <- 0
      solve_F <- state$F_at_catch
      state$F_at_catch <- function(target, F_keep) {
        n_solves <<- n_solves + 1
        solve_F(target, F_keep)
      }
      F_out <- constraint(F, state)
      solves[[length(solves) + 1]] <<- c(y = state$y, n_solves = n_solves)
      F_out
    }
  }

  down <- 0
  up <- 0.5
  out <- suppressWarnings(constraint_hand_run(F_constraint = count_solves(catch_stability(down = down, up = up))))
  catch <- constraint_hand_catch(out)

  # in year 4 only region 1 starts outside its band, and raising its F takes region 2 below last year's catch
  expect_equal(solves[[3]][["y"]], 4)
  expect_equal(solves[[3]][["n_solves"]], 2)
  change <- catch[,4] / catch[,3]
  expect_true(all(change > (1 - down) * (1 - 1e-4)))
  expect_true(all(change < (1 + up) * (1 + 1e-4)))
  expect_equal(change[2], 1 - down, tolerance = 1e-4) # region 2 ends on its floor

})
