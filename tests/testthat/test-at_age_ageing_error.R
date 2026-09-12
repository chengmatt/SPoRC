# Ageing error on the at-age data sources: each fleet's prediction read through its own matrix onto the
# observed ages, the arrays sized by those ages, and the operating model drawing the same way.

library(SPoRC)
library(testthat)

# prediction arrays for a direct call: 1 population, 2 regions, 3 years, 1 season, 4 ages, 1 sex, 2 fleets
ae_test_arrays <- function(seed = 1) {
  set.seed(seed)
  d <- c(1, 2, 3, 1, 4, 1, 2)
  list(
    CAA = array(stats::runif(prod(d), 10, 100), dim = d),
    DAA = array(stats::runif(prod(d), 1, 10), dim = d),
    SrvIAA = array(stats::runif(prod(d), 100, 1000), dim = d),
    WAA_fish = array(stats::runif(prod(d), 0.5, 2), dim = d),
    dmr = array(0.5, dim = c(2, 3, 1, 2)),
    catch_units = c(0, 1),   # fleet 1 in numbers, fleet 2 in weight
    discard_units = c(0, 1)
  )
}

# a share p of each age read one year young and p one year old
misread_matrix <- function(n_ages, p = 0.1) {
  m <- diag(1 - 2 * p, n_ages)
  for(a in seq_len(n_ages - 1)) {
    m[a, a + 1] <- p
    m[a + 1, a] <- p
  }
  m[1, 1] <- 1 - p              # the youngest can only be read older
  m[n_ages, n_ages] <- 1 - p    # and the plus group only younger
  m
}

# one matrix per fleet, repeated over years
ae_by_fleet <- function(mats, n_yrs = 3) {
  out <- array(0, dim = c(n_yrs, nrow(mats[[1]]), ncol(mats[[1]]), length(mats)))
  for(f in seq_along(mats)) for(y in seq_len(n_yrs)) out[y, , , f] <- mats[[f]]
  out
}

# the model's quantity at each model age for one region, year and fleet, read onto the observed ages
expected_at_age <- function(arrays, source, ae, r, y, f) {
  units <- if(source == "catch") arrays$catch_units[f] else arrays$discard_units[f]
  wt <- if(source != "srv_index" && units == 1) arrays$WAA_fish[1, r, y, 1, , 1, f] else 1
  x <- switch(source,
              catch = arrays$CAA[1, r, y, 1, , 1, f] * wt,
              discard = arrays$DAA[1, r, y, 1, , 1, f] / arrays$dmr[r, y, 1, f] * wt,
              srv_index = arrays$SrvIAA[1, r, y, 1, , 1, f])
  as.vector(x %*% ae[y, , , f])
}

# direct call of the at-age likelihood with every cell fit and the two regions split
at_age_call <- function(arrays, source, ae, n_obs_ages, pop = FALSE) {
  use_dim <- c(if(pop) 1, 2, 3, 1, n_obs_ages, 1, 2)
  use <- array(1, dim = use_dim)
  rho_dim <- c(if(pop) 1, 2, 1, 2)
  SPoRC:::get_at_age_source_nLL(
    obs_t = SPoRC:::prep_at_age_obs(array(50, dim = use_dim), use, c(0, 0)),
    use = use,
    ln_sigma = array(log(0.2), dim = c(if(pop) 1, n_obs_ages, 1, 2)),
    source = source,
    pop = pop,
    arrays = arrays,
    trans_rho = array(0, dim = rho_dim),
    trans_rho_year = array(0, dim = rho_dim),
    ageing_error = ae
  )
}

# catch or index at age at a plausible level, every cell fit
aa_obs <- function(n_yrs, n_obs_ages, n_fleets = 1, level = 2e3) array(level, dim = c(1, n_yrs, 1, n_obs_ages, 1, n_fleets))
aa_use <- function(n_yrs, n_obs_ages, n_fleets = 1) array(1, dim = c(1, n_yrs, 1, n_obs_ages, 1, n_fleets))


test_that("each at-age data source is predicted through its fleet's ageing error", {

  arrays <- ae_test_arrays()
  ae <- ae_by_fleet(list(misread_matrix(4), misread_matrix(4, p = 0.2)))

  for(source in c("catch", "discard", "srv_index")) {
    out <- at_age_call(arrays, source, ae, n_obs_ages = 4)
    for(f in 1:2) for(r in 1:2) for(y in 1:3) {
      expect_equal(out$pred[r, y, 1, , 1, f], expected_at_age(arrays, source, ae, r, y, f), tolerance = 1e-12)
    }
  }

  # the likelihood is evaluated against those predictions: the same call on arrays mapped by hand
  mapped <- arrays
  for(f in 1:2) for(r in 1:2) for(y in 1:3) mapped$SrvIAA[1, r, y, 1, , 1, f] <- expected_at_age(arrays, "srv_index", ae, r, y, f)
  expect_equal(at_age_call(arrays, "srv_index", ae, 4)$nLL, at_age_call(mapped, "srv_index", NULL, 4)$nLL, tolerance = 1e-12)

  # the population-specific form reads the same matrix
  pop_out <- at_age_call(arrays, "srv_index", ae, n_obs_ages = 4, pop = TRUE)
  for(f in 1:2) for(r in 1:2) for(y in 1:3) {
    expect_equal(pop_out$pred[1, r, y, 1, , 1, f], expected_at_age(arrays, "srv_index", ae, r, y, f), tolerance = 1e-12)
  }
})

test_that("an identity ageing error leaves the predictions and likelihood exactly as they were", {

  arrays <- ae_test_arrays()
  identity_ae <- ae_by_fleet(list(diag(4), diag(4)))

  for(source in c("catch", "discard", "srv_index")) {
    expect_identical(at_age_call(arrays, source, identity_ae, 4), at_age_call(arrays, source, NULL, 4))
  }

  # a fleet reading without error keeps its unmapped predictions while the other fleet is mapped
  mixed <- at_age_call(arrays, "catch", ae_by_fleet(list(diag(4), misread_matrix(4))), 4)
  direct <- at_age_call(arrays, "catch", NULL, 4)
  expect_identical(mixed$pred[, , , , , 1], direct$pred[, , , , , 1])
  expect_false(isTRUE(all.equal(mixed$pred[, , , , , 2], direct$pred[, , , , , 2])))
})

test_that("a collapsed plus group puts the at-age data on the observed ages", {

  arrays <- ae_test_arrays()
  plus_3 <- cbind(diag(4)[, 1:2], rowSums(diag(4)[, 3:4])) # ages 3 and 4 read as 3+
  ae <- ae_by_fleet(list(plus_3, plus_3))

  out <- at_age_call(arrays, "catch", ae, n_obs_ages = 3)
  expect_equal(dim(out$pred), c(2, 3, 1, 3, 1, 2))
  expect_equal(out$pred[, , 1, 3, 1, 1], arrays$CAA[1, , , 1, 3, 1, 1] + arrays$CAA[1, , , 1, 4, 1, 1], tolerance = 1e-12)
  expect_equal(out$pred[, , 1, 1:2, 1, 1], arrays$CAA[1, , , 1, 1:2, 1, 1], tolerance = 1e-12)

  # observations left on the model ages no longer line up with the matrix
  expect_error(at_age_call(arrays, "catch", ae, n_obs_ages = 4), "observed ages")
})

test_that("a model with a collapsed plus group fits catch at age on the observed ages", {

  n_yrs <- 20
  plus_4 <- cbind(diag(6)[, 1:3], rowSums(diag(6)[, 4:6])) # ages 4 to 6 read as 4+

  il <- build_at_age(n_yrs = n_yrs, n_ages = 6, ObsCatchAA = aa_obs(n_yrs, 4), UseCatchAA = aa_use(n_yrs, 4),
                     AgeingError = plus_4, AgeObsCorr_catch = "us")
  expect_equal(dim(il$par$ln_sigmaCAA), c(4, 1, 1))
  expect_equal(dim(il$par$trans_rho_catch_us)[1], 6) # 4 observed ages give 6 pairs

  rep <- at_age_rep(il)
  expect_equal(dim(rep$PredCatchAA), c(1, n_yrs, 1, 4, 1, 1))
  expect_equal(rep$PredCatchAA[1, , 1, 4, 1, 1], rowSums(rep$CAA[1, 1, , 1, 4:6, 1, 1]), tolerance = 1e-12)
  expect_equal(rep$PredCatchAA[1, , 1, 1:3, 1, 1], rep$CAA[1, 1, , 1, 1:3, 1, 1], tolerance = 1e-12, ignore_attr = TRUE)

  # the same data left on the model ages are refused at setup
  expect_error(build_at_age(n_yrs = n_yrs, n_ages = 6, ObsCatchAA = aa_obs(n_yrs, 6), UseCatchAA = aa_use(n_yrs, 6),
                            AgeingError = plus_4), "n_obs_ages")
})

test_that("each fishery fleet and the survey read their own ageing error", {

  n_yrs <- 20
  n_ages <- 5
  ae_fish <- array(c(misread_matrix(n_ages), diag(n_ages)), dim = c(n_ages, n_ages, 2)) # fleet 2 reads without error
  ae_srv <- array(misread_matrix(n_ages, p = 0.2), dim = c(n_ages, n_ages, 1))

  il <- build_at_age(
    n_yrs = n_yrs,
    n_ages = n_ages,
    n_fleets = 2,
    ObsCatchAA = aa_obs(n_yrs, n_ages, n_fleets = 2),
    UseCatchAA = aa_use(n_yrs, n_ages, n_fleets = 2),
    AgeingError_fish = ae_fish,
    AgeingError_srv = ae_srv,
    srv_extra = list(ObsSrvIdxAA = aa_obs(n_yrs, n_ages, level = 1e5), UseSrvIdxAA = aa_use(n_yrs, n_ages))
  )
  rep <- at_age_rep(il)

  for(f in 1:2) {
    expect_equal(rep$PredCatchAA[1, , 1, , 1, f], rep$CAA[1, 1, , 1, , 1, f] %*% ae_fish[, , f],
                 tolerance = 1e-12, ignore_attr = TRUE)
  }
  expect_equal(rep$PredSrvIdxAA[1, , 1, , 1, 1], rep$SrvIAA[1, 1, , 1, , 1, 1] %*% ae_srv[, , 1],
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("the separable correlation reads its predictions through the ageing error", {

  n_yrs <- 20
  n_ages <- 5
  misread <- misread_matrix(n_ages)

  il <- build_at_age(n_yrs = n_yrs, n_ages = n_ages, ObsCatchAA = aa_obs(n_yrs, n_ages), UseCatchAA = aa_use(n_yrs, n_ages),
                     AgeingError = misread, AgeObsCorr_catch = "2dar1")
  rep <- at_age_rep(il)
  expect_equal(rep$PredCatchAA[1, , 1, , 1, 1], rep$CAA[1, 1, , 1, , 1, 1] %*% misread, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("an observed age that no model age is read as is refused at setup", {

  n_yrs <- 20
  gap <- cbind(diag(5)[, 1:3], 0) # model ages 4 and 5 dropped, so nothing is read as observed age 4
  expect_error(build_at_age(n_yrs = n_yrs, n_ages = 5, ObsCatchAA = aa_obs(n_yrs, 4), UseCatchAA = aa_use(n_yrs, 4),
                            AgeingError = gap), "no model age is read as")
})

test_that("the operating model draws catch and survey index at age through the ageing error", {

  n_yrs <- 15
  n_ages <- 6
  n_obs <- 4
  plus_4 <- cbind(diag(n_ages)[, 1:3], rowSums(diag(n_ages)[, 4:6])) # ages 4 to 6 read as 4+

  sim_list <- Setup_Sim_Dim(
    n_sims = 1,
    n_yrs = n_yrs,
    n_regions = 1,
    n_ages = n_ages,
    n_obs_ages = n_obs,
    n_lens = NULL,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1
  )
  sim_list <- Setup_Sim_Containers(sim_list)

  curve <- function(slope, infl) array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
                                       dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
  use_aa <- array(1, dim = c(1, n_yrs, 1, n_obs, 1, 1))

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, curve(3, 2)),
    ret_sel_input = replicate(1, curve(3, 2)),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(0.15, dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 1, 1)),
    ln_sigmaCAA = array(log(0.2), dim = c(n_obs, 1, 1)),
    UseCatchAA = use_aa,
    use_catch_aa = 1
  )
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(1, curve(1, 3)),
    ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 1, 1)),
    ln_sigmaSrvIdxAA = array(log(0.2), dim = c(n_obs, 1, 1)),
    UseSrvIdxAA = use_aa,
    use_srv_idx_aa = 1
  )

  biol <- function(val) array(val, dim = c(1, 1, n_yrs, 1, n_ages, 1))
  sim_list <- suppressWarnings(Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, biol(1)),
    WAA_fish_input = replicate(1, array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol(1)),
    AgeingError_input = array(rep(plus_4, each = n_yrs), dim = c(n_yrs, n_ages, n_obs, 1))
  ))

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(0.3), dim = c(2, 1, 1)),
    recruitment_opt = "mean_rec",
    init_age_strc = 1
  )

  set.seed(11)
  om <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

  expect_equal(dim(om$ObsCatchAA), c(1, n_yrs, 1, n_obs, 1, 1, 1))
  for(y in seq_len(n_yrs)) {
    expect_equal(om$TrueCatchAA[1, y, 1, , 1, 1, 1], as.vector(om$CAA[1, 1, y, 1, , 1, 1, 1] %*% plus_4), tolerance = 1e-12)
    expect_equal(om$TrueSrvIdxAA[1, y, 1, , 1, 1, 1], as.vector(om$SrvIAA[1, 1, y, 1, , 1, 1, 1] %*% plus_4), tolerance = 1e-12)
  }
})

test_that("an at-age draw refuses an ageing error that does not match its observed ages", {

  numbers <- array(1, dim = c(1, 1, 6, 1))
  use <- array(1, dim = c(1, 4, 1))
  expect_error(SPoRC:::sim_at_age_cell(numbers, numbers, use, array(0, dim = c(1, 4, 1)), array(log(0.2), dim = c(4, 1)),
                                       1, 0, 0, FALSE, 1, ageing_error = diag(6)), "observed ages")
})
