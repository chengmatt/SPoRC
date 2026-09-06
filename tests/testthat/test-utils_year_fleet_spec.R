# The year-by-fleet specification grammar, which twelve composition data sources each parsed with their
# own copy of the same forty lines, already diverged in their error messages.
#
# They now share one parser, so this checks the grammar once and the integration test at the bottom checks
# every source still reaches it. The parser takes its vocabulary as an argument, so at-age sources fit too.

CODES <- c(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999)


test_that("a single block covering the series fills every year", {
  m <- parse_year_fleet_spec("agg_Year_1-terminal_Fleet_1", "X", 1, 10, CODES)

  expect_equal(dim(m), c(10L, 1L))
  expect_true(all(m == 0))
})


test_that("blocks partition the series and take effect at their own years", {
  m <- parse_year_fleet_spec(c("agg_Year_1-4_Fleet_1", "spltRspltS_Year_5-terminal_Fleet_1"),
                             "X", 1, 10, CODES)

  expect_equal(as.vector(m), c(rep(0, 4), rep(1, 6)))
})


test_that("a later block overwrites an earlier one where they overlap", {
  # This is what makes a general setting followed by a carve-out work, and it is
  # the behavior the composition data sources have always had.
  m <- parse_year_fleet_spec(c("agg_Year_1-terminal_Fleet_1", "spltRjntS_Year_3-4_Fleet_1"),
                             "X", 1, 10, CODES)

  expect_equal(as.vector(m), c(0, 0, 2, 2, 0, 0, 0, 0, 0, 0))
})


test_that("fleets are filled independently", {
  m <- parse_year_fleet_spec(c("agg_Year_1-terminal_Fleet_1",
                               "spltRspltS_Year_1-5_Fleet_2",
                               "none_Year_6-terminal_Fleet_2"),
                             "X", 2, 10, CODES)

  expect_true(all(m[, 1] == 0))
  expect_equal(m[, 2], c(rep(1, 5), rep(999, 5)))
})


test_that("the vocabulary is whatever the caller passes", {
  # The at-age data sources name their dims differently from the composition
  # data sources, and share the grammar rather than the values.
  aa <- c(agg = 0, spltRaggS = 1, aggRspltS = 2, spltRspltS = 3)
  m <- parse_year_fleet_spec(c("aggRspltS_Year_1-3_Fleet_1", "spltRspltS_Year_4-terminal_Fleet_1"),
                             "X", 1, 6, aa)

  expect_equal(as.vector(m), c(2, 2, 2, 3, 3, 3))
  expect_error(parse_year_fleet_spec("spltRjntS_Year_1-terminal_Fleet_1", "X", 1, 6, aa),
               "spltRjntS")
})


test_that("the caller's own constraint is applied per entry", {
  # Constraints that depend on other settings, such as a likelihood that cannot
  # take an aggregated observation, are passed in rather than built in.
  no_agg_on_fleet2 <- function(value, fleet)
    if(value == "agg" && fleet == 2) "fleet 2 cannot be aggregated" else NULL

  expect_error(
    parse_year_fleet_spec(c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2"),
                          "X", 2, 10, CODES, check = no_agg_on_fleet2),
    "fleet 2 cannot be aggregated")

  expect_no_error(
    parse_year_fleet_spec(c("agg_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2"),
                          "X", 2, 10, CODES, check = no_agg_on_fleet2))
})


test_that("a specification that does not cover every year and fleet is refused", {
  # The failure this replaces reported only that an NA had appeared. Saying which
  # year and fleet is uncovered is the difference between a hint and an answer.
  expect_error(parse_year_fleet_spec("agg_Year_1-5_Fleet_1", "X", 1, 10, CODES),
               "leaves 5 year and fleet combinations unset")
  expect_error(parse_year_fleet_spec("agg_Year_1-terminal_Fleet_1", "X", 2, 10, CODES),
               "fleet 2")
})


test_that("malformed specifications say what is wrong with them", {
  expect_error(parse_year_fleet_spec("bogus_Year_1-terminal_Fleet_1", "X", 1, 10, CODES),
               "Value is one of agg, spltRspltS, spltRjntS, none")
  expect_error(parse_year_fleet_spec("agg_Year_1-terminal_Fleet_9", "X", 1, 10, CODES),
               "naming fleet 9 of 1 fleets")
  expect_error(parse_year_fleet_spec("agg_Year_1-terminal", "X", 1, 10, CODES),
               "not a year and fleet specification")
  expect_error(parse_year_fleet_spec("agg_Year_5-2_Fleet_1", "X", 1, 10, CODES),
               "naming years 5 to 2")
  expect_error(parse_year_fleet_spec("agg_Year_1-999_Fleet_1", "X", 1, 10, CODES),
               "of a 10 year model")
  expect_error(parse_year_fleet_spec("agg_Year_x-y_Fleet_1", "X", 1, 10, CODES),
               "could not be read")
})


test_that("the argument name being parsed appears in every message", {
  # A model sets twelve of these and the failure has to say which one.
  for(msg in c("bogus_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_9",
               "agg_Year_1-5_Fleet_1")) {
    expect_error(parse_year_fleet_spec(msg, "SrvLenComps_pop_Type", 1, 10, CODES),
                 "SrvLenComps_pop_Type", label = msg)
  } # end msg loop
})


test_that("every composition data source still reaches the shared grammar", {
  # The twelve call sites, checked through the setup functions rather than
  # directly, so a data source wired to the wrong fleet count or year count shows up.
  skip_if_not(exists("sweep_input"), "the sweep test setup is not loaded")

  data_sources <- list(
    list(stage = "fishidx", arg = "FishAgeComps_Type"),
    list(stage = "fishidx", arg = "FishLenComps_Type"),
    list(stage = "srvidx",  arg = "SrvAgeComps_Type"),
    list(stage = "srvidx",  arg = "SrvLenComps_Type"))

  for(s in data_sources) {
    ov <- list(); ov[[s$arg]] <- "bogus_Year_1-terminal_Fleet_1"
    args <- list(dims = list(
      n_regions = 1,
      n_sexes = 1,
      n_fish_fleets = 1,
      n_srv_fleets = 1,
      n_yrs = 10,
      n_ages = 6
    ))
    args[[s$stage]] <- ov
    expect_error(do.call(sweep_input, args), s$arg, label = s$arg)
  } # end s loop
})
