# Constraints on the fishing mortality a control rule or input gives the projection: a catch cap and a
# stability band, applied through F_constraint. Realized catch is checked at 1e-4 relative, the solver's
# own tolerance, and year 1 is never constrained since it replays the terminal year.

library(SPoRC)
library(testthat)

rule <- function(x, frp, brp) frp * min(1, x / brp)

# regional annual retained catch by projection year, summed over populations, seasons and fleets
annual_catch <- function(out) {
  n_regions <- dim(out$proj_Catch)[2]
  n_proj_yrs <- dim(out$proj_Catch)[3]
  m <- array(0, dim = c(n_regions, n_proj_yrs))
  for(r in 1:n_regions) for(y in 1:n_proj_yrs) m[r,y] <- sum(out$proj_Catch[,r,y,,])
  m
}

# a two region, one season, one fleet, one sex stock with mixing, so a joint solve is live
two_region_inputs <- function(n_proj_yrs = 5, n_ages = 10, f = 0.1) {
  n_pop <- 1
  n_regions <- 2
  n_seas <- 1
  n_sexes <- 1
  n_fish_fleets <- 1
  ages <- 1:n_ages
  biol_d <- c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)

  natmort <- array(0.15, dim = biol_d)
  WAA <- array(0, dim = biol_d)
  MatAA <- array(0, dim = biol_d)
  fish_sel <- array(0, dim = c(biol_d, n_fish_fleets))
  for(r in 1:n_regions) for(y in 1:n_proj_yrs) {
    WAA[1,r,y,1,,1] <- (1 - exp(-0.2 * ages))^3 * 5
    MatAA[1,r,y,1,,1] <- 1 / (1 + exp(-0.8 * (ages - 4)))
    fish_sel[1,r,y,1,,1,1] <- 1 / (1 + exp(-1.1 * (ages - 3 - 0.3 * r)))
  }
  WAA_fish <- array(WAA, dim = c(biol_d, n_fish_fleets))

  # a fifth of each region moves to the other every year
  Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  for(y in 1:n_proj_yrs) for(a in 1:n_ages) Movement[1,,,y,1,a,1] <- matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)
  sgl_move <- array(0, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_ages, n_sexes))
  for(y in 1:n_proj_yrs) for(a in 1:n_ages) sgl_move[1,,,y,a,1] <- diag(n_regions)

  terminal_NAA <- array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  for(r in 1:n_regions) terminal_NAA[1,r,1,,1] <- 1e6 * exp(-0.3 * (ages - 1)) * (1 + 0.3 * r)

  list(
    n_proj_yrs = n_proj_yrs, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages,
    n_sexes = n_sexes, n_fish_fleets = n_fish_fleets, n_seas = n_seas,
    seasdur = 1, spawn_seas = 1, t_spawn = 0, natal_region = 1, do_recruits_move = 0,
    stray_rate = array(0, dim = c(n_pop, n_proj_yrs)),
    sexratio = array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)),
    recruitment = array(3e5, dim = c(n_pop, n_regions, 10)),
    terminal_NAA = terminal_NAA, terminal_NAA0 = terminal_NAA * 1.5,
    terminal_F = array(f, dim = c(n_regions, n_seas, n_fish_fleets)),
    natmort = natmort, WAA = WAA, WAA_fish = WAA_fish, MatAA = MatAA, fish_sel = fish_sel,
    Movement = Movement, sgl_seas_spawning_movement = sgl_move,
    recruitment_opt = "mean_rec",
    # a reference biomass far below the stock, so the rule returns frp in every year
    f_ref_pt = array(f, dim = c(n_regions, n_proj_yrs)),
    b_ref_pt = array(1, dim = c(n_pop, n_regions, n_proj_yrs)),
    HCR_function = rule
  )
}

test_that("a cap that never binds changes nothing, and proj_F_rule equals proj_F", {

  base <- project_with_HCR(rule, n_proj_yrs = 8)
  capped <- project_with_HCR(rule, n_proj_yrs = 8, F_constraint = catch_cap(1e6))

  expect_equal(capped$proj_F, base$proj_F)
  expect_equal(capped$proj_Catch, base$proj_Catch)
  expect_equal(capped$proj_SSB, base$proj_SSB)
  expect_equal(capped$proj_F_rule, capped$proj_F)
  expect_equal(base$proj_F_rule, base$proj_F)
})

test_that("a binding cap holds catch at the cap and lowers F below the rule's, under a rule or an input F", {

  cap <- 30
  base <- project_with_HCR(rule, n_proj_yrs = 8)
  capped <- project_with_HCR(rule, n_proj_yrs = 8, F_constraint = catch_cap(cap))
  base_catch <- annual_catch(base)
  catch <- annual_catch(capped)

  # year 1 replays the terminal year, and the rule's catch in year 2 is well over the cap
  expect_equal(catch[,1], base_catch[,1])
  expect_gt(base_catch[1,2], cap)
  expect_true(all(catch[,2:8] <= cap * (1 + 1e-4)))

  # year 2 starts from the same numbers as the unconstrained run, so the rule's F is the same and the cap binds exactly
  expect_equal(capped$proj_F_rule[1,2], base$proj_F[1,2])
  expect_equal(catch[1,2], cap, tolerance = 1e-4)
  expect_lt(capped$proj_F[1,2], capped$proj_F_rule[1,2])

  # in every later year F is either lowered onto the cap or left at the rule's
  for(y in 3:8) {
    if(capped$proj_F[1,y] < capped$proj_F_rule[1,y]) expect_equal(catch[1,y], cap, tolerance = 1e-4)
    else expect_equal(capped$proj_F[1,y], capped$proj_F_rule[1,y])
  }

  # the same cap on a fixed input F
  input <- project_at_F(f = 0.06, n_proj_yrs = 6, F_constraint = catch_cap(cap))
  expect_true(all(annual_catch(input)[,2:6] <= cap * (1 + 1e-4)))
  expect_equal(annual_catch(input)[1,2], cap, tolerance = 1e-4)
})

test_that("a stability band keeps each year's catch within the band of last year's", {

  down <- 0.1
  up <- 0.1
  out <- project_with_HCR(rule, n_proj_yrs = 8, F_constraint = catch_stability(down = down, up = up))
  catch <- annual_catch(out)
  change <- catch[1,2:8] / catch[1,1:7] - 1

  expect_true(all(change >= -down - 1e-4 & change <= up + 1e-4))

  # the rule's catch rises by half in year 2, so the band binds there at the top edge
  expect_equal(catch[1,2], catch[1,1] * (1 + up), tolerance = 1e-4)
  expect_lt(out$proj_F[1,2], out$proj_F_rule[1,2])

  # a year where F was moved sits at an edge of the band; one where it was not is at the rule's F
  for(y in 2:8) {
    if(out$proj_F[1,y] != out$proj_F_rule[1,y]) expect_true(abs(change[y - 1] - up) < 1e-4 || abs(change[y - 1] + down) < 1e-4)
  }
})

test_that("constraints apply in the order given, so the last one is the one respected", {

  cap <- 20
  floor_only <- catch_stability(down = 0, up = Inf) # catch may not fall below last year's

  cap_last <- project_with_HCR(rule, n_proj_yrs = 4, F_constraint = list(floor_only, catch_cap(cap)))
  floor_last <- project_with_HCR(rule, n_proj_yrs = 4, F_constraint = list(catch_cap(cap), floor_only))

  # the cap wins when it comes last; the floor lifts the capped catch back to last year's when it comes last
  expect_equal(annual_catch(cap_last)[1,2], cap, tolerance = 1e-4)
  expect_equal(annual_catch(floor_last)[1,2], annual_catch(floor_last)[1,1], tolerance = 1e-4)
  expect_gt(annual_catch(floor_last)[1,2], cap)
})

test_that("the constructors invert a known catch, skip what cannot bind, and refuse bad settings", {

  # a stand-in for the projection's state: catch rises linearly in F, so a target is easy to invert
  state <- list(y = 2, F_rule = c(0.1, 0.2), last_catch = c(0, 50),
                catch_at_F = function(F) 100 * F,
                F_at_catch = function(target, F) ifelse(is.na(target), F, target / 100))
  F <- c(0.1, 0.2)

  expect_equal(catch_cap(100)(F, state), F) # catches of 10 and 20, both under
  expect_equal(catch_cap(15)(F, state), c(0.1, 0.15)) # region 2 over, lowered to 15
  expect_equal(catch_cap(24, scope = "global")(F, state), c(0.08, 0.16)) # total 30 to 24, shares kept
  expect_equal(catch_cap(c(100, 15))(F, state), c(0.1, 0.15)) # one cap per region
  cap_by_year <- matrix(c(100, 100, 5, 100), nrow = 2) # year 2: region 1 capped at 5
  expect_equal(catch_cap(cap_by_year)(F, state), c(0.05, 0.2))

  # region 1 took nothing last year so its band is open; region 2 is below 80 percent of 50 and is raised to 40
  expect_equal(catch_stability(down = 0.2, up = 0.2)(F, state), c(0.1, 0.4))
  expect_equal(catch_stability(down = 1, up = Inf)(F, state), F) # a band that cannot bind
  no_catch <- state
  no_catch$last_catch <- c(0, 0)
  expect_equal(catch_stability(scope = "global")(F, no_catch), F)

  expect_error(catch_cap(-1), "non-negative")
  expect_error(catch_cap(1, scope = "fleet"), "scope")
  expect_error(catch_cap(matrix(1, 2, 2), scope = "global"), "not a matrix")
  expect_error(catch_cap(c(1, 2, 3))(F, state), "one value")
  expect_error(catch_stability(down = 1.5), "between 0 and 1")
  expect_error(catch_stability(up = -0.1), "non-negative")
})

test_that("F_constraint is validated by the projection", {

  expect_error(project_with_HCR(rule, n_proj_yrs = 4, F_constraint = "cap"), "function")
  expect_error(project_with_HCR(rule, n_proj_yrs = 4, F_constraint = list(catch_cap(30), 2)), "not a function")
  expect_error(project_with_HCR(rule, n_proj_yrs = 4, F_constraint = function(F, state) -F), "non-negative")
  expect_error(project_with_HCR(rule, n_proj_yrs = 4, F_constraint = function(F, state) c(F, F)), "length")
})

test_that("ad_value reads values off a tape and leaves plain objects alone", {

  plain <- list(a = matrix(1:4, 2), b = NULL, c = list(d = 2))
  expect_identical(ad_value(plain), plain)
  expect_identical(ad_value(3), 3)

  seen <- new.env()
  RTMB::MakeTape(function(x) {
    "[<-" <- RTMB::ADoverload("[<-")
    m <- matrix(0, 2, 2)
    m[1,2] <- x
    v <- ad_value(list(m = m, keep = NULL, n = 5))
    seen$class <- class(v$m)
    seen$dim <- dim(v$m)
    seen$value <- v$m[1,2]
    seen$names <- names(v)
    seen$vec_names <- names(ad_value(x * c(a = 1, b = 2)))
    seen$dimnames <- dimnames(ad_value(array(x * 1:4, dim = c(2, 2), dimnames = list(c("r1", "r2"), c("c1", "c2")))))
    sum(m)
  }, 1.5)
  expect_false("advector" %in% seen$class)
  expect_equal(seen$dim, c(2, 2))
  expect_equal(seen$value, 1.5)
  expect_equal(seen$names, c("m", "keep", "n"))
  expect_equal(seen$vec_names, c("a", "b")) # names and dimnames come through with the numbers
  expect_equal(seen$dimnames, list(c("r1", "r2"), c("c1", "c2")))
})

# the projection taped over log F and a log scale on the terminal numbers, so the cap's F solve
# depends on the parameters through this year's numbers as well as through the rate
sable_constrained <- function(p, F_constraint, n_proj_yrs = 5, what = "proj_SSB") {
  rp <- sgl_rg_sable_rep
  rp$NAA <- rp$NAA * exp(p[2])
  out <- project_at_F(f = exp(p[1]), n_proj_yrs = n_proj_yrs, rep = rp, F_constraint = F_constraint)
  sum(out[[what]])
}

test_that("a binding cap tapes, and its gradient matches finite differences", {

  p0 <- c(log(0.06), 0)
  for(what in c("proj_SSB", "proj_Catch")) {
    f <- function(p) sable_constrained(p, catch_cap(30), what = what)
    tape <- RTMB::MakeTape(f, p0)
    expect_equal(tape(p0), f(p0), tolerance = 1e-10)
    expect_equal(as.vector(tape$jacobian(p0)), numDeriv::grad(f, p0), tolerance = 1e-5)
  }

  # the cap binds in every projected year, so the catch is fixed and nothing depends on the rule's F;
  # the derivative through the solved F shows up against the terminal numbers instead
  f_ssb <- function(p) sable_constrained(p, catch_cap(30), what = "proj_SSB")
  grad_ssb <- as.vector(RTMB::MakeTape(f_ssb, p0)$jacobian(p0))
  expect_lt(abs(grad_ssb[1]), 1e-8)
  expect_gt(abs(grad_ssb[2]), 1)

  # a cap that never binds leaves the gradient of the unconstrained projection
  f_loose <- function(p) sable_constrained(p, catch_cap(1e6))
  f_none <- function(p) sable_constrained(p, NULL)
  expect_equal(as.vector(RTMB::MakeTape(f_loose, p0)$jacobian(p0)),
               as.vector(RTMB::MakeTape(f_none, p0)$jacobian(p0)), tolerance = 1e-10)
})

test_that("a stability band tapes, and its gradient matches finite differences", {

  p0 <- c(log(0.06), 0)
  f <- function(p) sable_constrained(p, catch_stability(down = 0.1, up = 0.1))
  tape <- RTMB::MakeTape(f, p0)
  expect_equal(tape(p0), f(p0), tolerance = 1e-10)
  expect_equal(as.vector(tape$jacobian(p0)), numDeriv::grad(f, p0), tolerance = 1e-5)
  # the band holds every year, so the path depends on the terminal numbers and not on F
  expect_lt(abs(tape$jacobian(p0)[1,1]), 1e-6)
  expect_gt(abs(tape$jacobian(p0)[1,2]), 1e-3)
})

test_that("a global cap on two regions tapes with the shares on the tape", {

  inp <- two_region_inputs(n_proj_yrs = 4)
  base <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "HCR_global")))
  cap <- 0.8 * sum(annual_catch(base)[,2])

  # on the tape the rule reads its status as a value, the documented way to branch
  inp$HCR_function <- function(x, frp, brp) frp * min(1, ad_value(x) / brp)
  f <- function(ln_F) {
    inp$f_ref_pt <- array(exp(ln_F), dim = c(2, 4))
    out <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "HCR_global", F_constraint = catch_cap(cap, scope = "global"), catch_tol = 1e-10)))
    sum(out$proj_SSB)
  }
  ln_F0 <- log(c(0.1, 0.12))
  tape <- RTMB::MakeTape(f, ln_F0)
  expect_equal(tape(ln_F0), f(ln_F0), tolerance = 1e-10)
  expect_equal(as.vector(tape$jacobian(ln_F0)), numDeriv::grad(f, ln_F0), tolerance = 1e-5)
})

test_that("a global cap under a global rule holds the total and keeps the rule's regional shares", {

  inp <- two_region_inputs(n_proj_yrs = 5)
  base <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "HCR_global")))
  base_catch <- annual_catch(base)
  cap <- 0.8 * sum(base_catch[,2])

  capped <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "HCR_global", F_constraint = catch_cap(cap, scope = "global"))))
  catch <- annual_catch(capped)

  expect_equal(sum(catch[,2]), cap, tolerance = 1e-4)
  expect_equal(catch[,2] / sum(catch[,2]), base_catch[,2] / sum(base_catch[,2]), tolerance = 1e-4)
  expect_true(all(colSums(catch[,2:5]) <= cap * (1 + 1e-4)))
  expect_true(all(capped$proj_F[,2] < capped$proj_F_rule[,2]))

  # a regional cap under the same rule caps each region on its own, solved jointly through the mixing
  reg_cap <- 0.9 * base_catch[,2]
  reg <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "HCR_global", F_constraint = catch_cap(reg_cap))))
  expect_equal(annual_catch(reg)[,2], reg_cap, tolerance = 1e-4)

  # capping one region alone leaves the other at the rule's F, with its catch moved by the mixing only
  one <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "HCR_global", F_constraint = catch_cap(c(reg_cap[1], 1e9)))))
  expect_equal(annual_catch(one)[1,2], reg_cap[1], tolerance = 1e-4)
  expect_equal(one$proj_F[2,2], one$proj_F_rule[2,2])
})
