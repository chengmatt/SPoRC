# At-age aggregation that changes part way through a series.
#
# The composition streams have always taken a setting per year and fleet; the
# at-age streams took one per fleet. They now share the grammar, so a fleet that
# reported sexes combined early and split later can be stated once rather than
# split into two fleets to work around the API.
#
# A bare value still means what it always did, this setting for the whole series,
# so nothing written before this needs changing.

aa_fixture <- function(type, nx = 2, nr = 1, ny = 10, na = 6, flagged_sex = 1,
                       flagged_region = NULL, extra = list()) {
  aa <- c(nr, ny, 1, na, nx, 1)
  use <- array(0, dim = aa)
  if(is.null(flagged_region)) use[, , , , flagged_sex, ] <- 1 else use[flagged_region, , , , , ] <- 1
  obs <- array(0, dim = aa); obs[use == 1] <- 100
  sweep_input(
    dims = list(n_regions = nr, n_sexes = nx, n_fish_fleets = 1, n_srv_fleets = 1,
                n_yrs = ny, n_ages = na),
    catch = c(list(ObsCatchAA = obs, UseCatchAA = use, CatchAA_Type = type,
                   UseCatch = array(0, dim = c(nr, ny, 1, 1))), extra))
}

aa_nll <- function(type, ...) {
  il <- aa_fixture(type, ...)
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  as.numeric(obj$fn(obj$par))
}

# The four at-age streams, each with whatever else it needs to carry an
# observation. A discard fraction is a property of the catch rather than of an
# age, so a fleet reporting discards at age reports them in numbers and needs a
# retention curve for there to be any discards at all.
AA_STREAMS <- list(
  list(name = "CatchAA", stage = "catch", extra = function(ny, nr)
    list(UseCatch = array(0, dim = c(nr, ny, 1, 1)))),
  list(name = "DiscardAA", stage = "catch", extra = function(ny, nr)
    list(discard_units = array("abd", dim = 1),
         ObsDiscard = array(50, dim = c(nr, ny, 1, 1)),
         UseDiscard = array(0, dim = c(nr, ny, 1, 1)))),
  list(name = "SrvIdxAA", stage = "srvidx", extra = function(ny, nr)
    list(UseSrvIdx = array(0, dim = c(nr, ny, 1, 1)), srv_idx_type = "none")))

stream_nll <- function(st, type, nx = 2, nr = 1, ny = 10, na = 6) {
  aa <- c(nr, ny, 1, na, nx, 1)
  use <- array(0, dim = aa); use[, , , , 1, ] <- 1
  obs <- array(0, dim = aa); obs[use == 1] <- 10

  args <- st$extra(ny, nr)
  args[[paste0("Obs", st$name)]] <- obs
  args[[paste0("Use", st$name)]] <- use
  args[[paste0(st$name, "_Type")]] <- type

  call_args <- list(dims = list(n_regions = nr, n_sexes = nx, n_fish_fleets = 1,
                                n_srv_fleets = 1, n_yrs = ny, n_ages = na))
  call_args[[st$stage]] <- args
  if(identical(st$name, "DiscardAA"))
    call_args$fishsel <- list(use_fixed_ret_sel = 0, ret_sel_model = "logist1_Fleet_1",
                              ret_fixed_sel_pars_spec = "est_all")

  il <- do.call(sweep_input, call_args)
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  as.numeric(obj$fn(obj$par))
}


test_that("a bare value still sets the whole series", {
  # The form every model written before the grammar uses.
  il <- aa_fixture("spltRspltS")

  expect_equal(dim(il$data$CatchAA_Type), c(10L, 1L))
  expect_true(all(il$data$CatchAA_Type == 3))
})


test_that("one value per fleet still sets each fleet's whole series", {
  aa <- c(1, 10, 1, 6, 1, 2)
  use <- array(1, dim = aa); obs <- array(100, dim = aa)
  il <- sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 2, n_srv_fleets = 1,
                n_yrs = 10, n_ages = 6),
    catch = list(ObsCatchAA = obs, UseCatchAA = use,
                 CatchAA_Type = c("agg", "spltRspltS"),
                 UseCatch = array(0, dim = c(1, 10, 1, 2))))

  expect_equal(dim(il$data$CatchAA_Type), c(10L, 2L))
  expect_true(all(il$data$CatchAA_Type[, 1] == 0))
  expect_true(all(il$data$CatchAA_Type[, 2] == 3))
})


test_that("the year grammar sets each block", {
  il <- aa_fixture(c("spltRspltS_Year_1-5_Fleet_1", "spltRaggS_Year_6-terminal_Fleet_1"))

  expect_equal(il$data$CatchAA_Type[, 1], c(rep(3, 5), rep(1, 5)))
})


test_that("a setting that changes by year changes the objective by year", {
  # The sharp one, run on all four at-age streams. Summing over sexes and
  # splitting them are different observations, so a stream doing one for five
  # years and the other for five must land strictly between the two models that
  # do one throughout. Without this the grammar could parse correctly and be
  # ignored downstream, which is exactly what a per-fleet read of a year by fleet
  # matrix would do.
  #
  # Aggregating over sex rather than region on purpose: the fixture has F = 0
  # outside region 1, so region aggregation is a no-op there and the comparison
  # would pass while proving nothing.
  for(st in AA_STREAMS) {
    agg <- stream_nll(st, "spltRaggS")
    splt <- stream_nll(st, "spltRspltS")
    mixed <- stream_nll(st, c("spltRspltS_Year_1-5_Fleet_1",
                              "spltRaggS_Year_6-terminal_Fleet_1"))

    expect_false(isTRUE(all.equal(agg, splt)),
                 label = paste(st$name, "summing over sexes differs from splitting them"))
    expect_gt(mixed, min(agg, splt))
    expect_lt(mixed, max(agg, splt))
  } # end st loop
})


test_that("every at-age stream builds a tape and a finite gradient when its setting varies", {
  # Evaluating the objective on doubles is not enough: an operation that drops
  # the tape leaves the value right and the model unfittable.
  for(st in AA_STREAMS) {
    aa <- c(1, 10, 1, 6, 2, 1)
    use <- array(0, dim = aa); use[, , , , 1, ] <- 1
    obs <- array(0, dim = aa); obs[use == 1] <- 10
    args <- st$extra(10, 1)
    args[[paste0("Obs", st$name)]] <- obs
    args[[paste0("Use", st$name)]] <- use
    args[[paste0(st$name, "_Type")]] <- c("spltRspltS_Year_1-5_Fleet_1",
                                          "spltRaggS_Year_6-terminal_Fleet_1")
    call_args <- list(dims = list(n_regions = 1, n_sexes = 2, n_fish_fleets = 1,
                                  n_srv_fleets = 1, n_yrs = 10, n_ages = 6))
    call_args[[st$stage]] <- args
    if(identical(st$name, "DiscardAA"))
      call_args$fishsel <- list(use_fixed_ret_sel = 0, ret_sel_model = "logist1_Fleet_1",
                                ret_fixed_sel_pars_spec = "est_all")

    il <- do.call(sweep_input, call_args)
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)

    expect_true(is.finite(obj$fn(obj$par)), label = paste(st$name, "objective"))
    expect_true(all(is.finite(obj$gr(obj$par))), label = paste(st$name, "gradient"))
  } # end st loop
})


test_that("a fleet fitting 2dar1 cannot change its aggregation between years", {
  # A separable correlation is defined over the whole block of years by ages, so
  # a setting that changes inside the block has no meaning. The message says what
  # to do instead rather than only what is wrong.
  expect_error(
    aa_nll(c("spltRspltS_Year_1-5_Fleet_1", "spltRaggS_Year_6-terminal_Fleet_1"),
           extra = list(CatchAA_LikeType = "lognormal", rho_catch_spec = "est_2dar1")),
    "2dar1")
})


test_that("a summed margin is checked in the years it applies to", {
  # Region is summed only from year 6, so flagging both regions is a problem in
  # those years and not before. The check runs per year rather than per fleet.
  aa <- c(2, 10, 1, 6, 1, 1)
  use <- array(0, dim = aa); use[, , , , , ] <- 1
  obs <- array(100, dim = aa)
  build <- function(type) sweep_input(
    dims = list(n_regions = 2, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                n_yrs = 10, n_ages = 6),
    catch = list(ObsCatchAA = obs, UseCatchAA = use, CatchAA_Type = type,
                 UseCatch = array(0, dim = c(2, 10, 1, 1))))

  expect_error(build(c("spltRspltS_Year_1-5_Fleet_1", "agg_Year_6-terminal_Fleet_1")),
               "in year 6")
  expect_no_error(build("spltRspltS"))
})


test_that("every at-age stream carries the setting as year by fleet", {
  # Eight streams reach the same helper. A stream wired to the wrong fleet count
  # or the wrong year count would still build, and would then be read with the
  # wrong shape in the objective.
  il <- aa_fixture("spltRspltS")
  n_yrs <- length(il$data$years)

  for(nm in c("CatchAA_Type", "CatchAA_pop_Type", "DiscardAA_Type", "DiscardAA_pop_Type",
              "SrvIdxAA_Type", "SrvIdxAA_pop_Type")) {
    expect_equal(dim(il$data[[nm]])[1], n_yrs, label = paste(nm, "year margin"))
    expect_equal(dim(il$data[[nm]])[2],
                 if(grepl("^Srv", nm)) il$data$n_srv_fleets else il$data$n_fish_fleets,
                 label = paste(nm, "fleet margin"))
  } # end nm loop
})


test_that("the at-age vocabulary is refused on the composition streams and the reverse", {
  # The two families share the grammar and not the values, so a value from the
  # wrong family has to be named rather than silently accepted.
  expect_error(at_age_type_matrix("spltRjntS", 1, 10, "CatchAA_Type"), "spltRjntS")
  expect_error(parse_year_fleet_spec("spltRaggS_Year_1-terminal_Fleet_1", "FishAgeComps_Type",
                                     1, 10, c(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999)),
               "spltRaggS")
})


test_that("a retrospective peel trims the aggregation alongside its observations", {
  # The setting carries a year margin now, so a peel that trimmed only the
  # observations would leave the two out of step.
  il <- aa_fixture(c("spltRspltS_Year_1-5_Fleet_1", "spltRaggS_Year_6-terminal_Fleet_1"))

  expect_equal(nrow(il$data$CatchAA_Type), 10L)
  trimmed <- il$data$CatchAA_Type[1:7, , drop = FALSE]
  expect_equal(nrow(trimmed), 7L)
  expect_equal(trimmed[, 1], c(rep(3, 5), rep(1, 2)))
})


test_that("the fits plot names a summed margin instead of naming a slot", {
  # A margin the fleet sums over holds its observation in slot one, and that slot
  # is not the first sex. Labelling the facet "Sex 1" would name a sex the
  # observation is not about, which is what the plot did before the setting was
  # readable per row.
  d <- c(1, 10, 1, 6, 2, 1)
  use <- array(0, dim = d); use[, , , , 1, ] <- 1
  obs <- array(0, dim = d); obs[use == 1] <- 100

  build <- function(type) sweep_input(
    dims = list(n_regions = 1, n_sexes = 2, n_fish_fleets = 1, n_srv_fleets = 1,
                n_yrs = 10, n_ages = 6),
    catch = list(ObsCatchAA = obs, UseCatchAA = use, CatchAA_Type = type,
                 UseCatch = array(0, dim = c(1, 10, 1, 1))))

  summed <- build("spltRaggS")
  p <- get_at_age_fits_plot(list(summed$data), list(at_age_rep(summed)), "m1",
                            stream = "CatchAA")
  expect_setequal(unique(p$data$Sex), "All sexes")

  # and a setting that changes part way through says so for the years it applies
  mixed <- build(c("spltRspltS_Year_1-5_Fleet_1", "spltRaggS_Year_6-terminal_Fleet_1"))
  pm <- get_at_age_fits_plot(list(mixed$data), list(at_age_rep(mixed)), "m1",
                             stream = "CatchAA")
  expect_setequal(unique(pm$data$Sex), c("Sex 1", "All sexes"))
  expect_setequal(unique(pm$data$Year[pm$data$Sex == "All sexes"]),
                  mixed$data$years[6:10])
})
