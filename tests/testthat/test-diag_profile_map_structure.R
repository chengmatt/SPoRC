library(SPoRC)
library(testthat)

# do_likelihood_profile() rebuilds the map for the profiled parameter at every grid value.
# That rebuild has to start from the map the model was fitted with, because the map is the
# only record of which positions were held fixed and which were estimated as a single
# shared parameter. Renumbering every position uniquely instead frees parameters the fitted
# model never estimated, which lowers the likelihood at every grid value and leaves the
# difference from the MLE uncomparable.
#
# The map factor carries that structure in its level labels: a shared parameter shows up as
# one label on several positions, and a fixed position shows up as NA.

par_array <- array(0, dim = c(2, 3))

test_that("a position the fitted map fixed is still fixed through the profile", {
  # Positions 2 and 5 are fixed in the fitted model, position 3 is the profile target.
  fitted_map <- factor(c("1", NA, "2", "3", NA, "4"))
  prof_map <- SPoRC:::build_profile_map(par_array, fitted_map, 3)

  expect_true(all(is.na(prof_map[c(2, 5)])))
})

test_that("two positions sharing a level still share after the rebuild", {
  # Positions 1 and 4 are one mirrored parameter, position 6 is the profile target.
  fitted_map <- factor(c("1", "2", "3", "1", "4", "5"))
  prof_map <- SPoRC:::build_profile_map(par_array, fitted_map, 6)

  expect_false(is.na(prof_map[1]))
  expect_identical(as.character(prof_map[1]), as.character(prof_map[4]))
})

test_that("the profile does not free more parameters than the fitted model estimated", {
  # Three free parameters in the fitted map: one shared across positions 1 and 4, plus
  # positions 3 and 6. Fixing position 6 leaves two.
  fitted_map <- factor(c("1", NA, "2", "1", NA, "3"))
  prof_map <- SPoRC:::build_profile_map(par_array, fitted_map, 6)

  expect_identical(nlevels(droplevels(prof_map)), 2L)
})

test_that("the target positions come back NA", {
  fitted_map <- factor(as.character(1:6))
  prof_map <- SPoRC:::build_profile_map(par_array, fitted_map, c(2, 5))

  expect_true(all(is.na(prof_map[c(2, 5)])))
  expect_false(any(is.na(prof_map[c(1, 3, 4, 6)])))
})

test_that("several targets supplied one index at a time all come back NA", {
  # The profile hands `idx` to `[<-` one element at a time, so a list of index vectors
  # reaches the map the same way it reaches the parameter values.
  fitted_map <- factor(as.character(1:6))
  prof_map <- SPoRC:::build_profile_map(par_array, fitted_map, list(c(1, 2), 5))

  expect_true(all(is.na(prof_map[c(1, 2, 5)])))
  expect_false(any(is.na(prof_map[c(3, 4, 6)])))
})

test_that("the rebuild is stable when fed back its own output", {
  # The sequential branch rebuilds the map once per grid value, so a rebuild that read its
  # own previous result would have to leave the structure unchanged. It holds the fitted
  # map aside instead, but the helper is idempotent either way.
  fitted_map <- factor(c("1", NA, "2", "1", NA, "3"))
  once <- SPoRC:::build_profile_map(par_array, fitted_map, 6)
  twice <- SPoRC:::build_profile_map(par_array, once, 6)

  expect_identical(as.character(once), as.character(twice))
})

test_that("a parameter with no fitted map gets every position estimated on its own", {
  prof_map <- SPoRC:::build_profile_map(par_array, NULL, 3)

  expect_true(is.na(prof_map[3]))
  expect_identical(nlevels(droplevels(prof_map)), 5L)
})

test_that("a map that does not match the parameter length falls back", {
  # A stale map entry cannot describe the parameter, so the profile treats it as absent
  # rather than recycling it onto the wrong positions.
  prof_map <- SPoRC:::build_profile_map(par_array, factor(as.character(1:4)), 3)

  expect_identical(length(prof_map), 6L)
  expect_identical(nlevels(droplevels(prof_map)), 5L)
})

test_that("a fully fixed parameter stays fully fixed", {
  prof_map <- SPoRC:::build_profile_map(par_array, factor(rep(NA_character_, 6)), 3)

  expect_true(all(is.na(prof_map)))
})

# The mirror check warns when a profile target is tied by the fitted map to a position the
# profile does not fix. The shared parameter stays estimated in that case, so the grid value
# never holds and the profile comes back flat.

test_that("a target mirrored onto an untargeted position is flagged", {
  fitted_map <- factor(c("1", "2", "3", "1", "4", "5"))
  expect_warning(SPoRC:::check_profile_mirrors(par_array, fitted_map, "fish_fixed_sel_pars", 1),
                 "share a map level with positions outside `idx`")
})

test_that("profiling every position of a shared parameter is not flagged", {
  fitted_map <- factor(c("1", "2", "3", "1", "4", "5"))
  expect_silent(SPoRC:::check_profile_mirrors(par_array, fitted_map, "fish_fixed_sel_pars", c(1, 4)))
})

test_that("an unshared target is not flagged", {
  fitted_map <- factor(c("1", NA, "2", "3", NA, "4"))
  expect_silent(SPoRC:::check_profile_mirrors(par_array, fitted_map, "ln_M", 3))
})

# The helper above is only worth anything if both execution branches actually route through
# it, so the tests below drive do_likelihood_profile() end to end and read back the map and
# the parameter values RTMB was handed. RTMB::MakeADFun is stubbed out, so nothing is
# optimized and no model is built; the stub records its arguments and then fails, which the
# profile's own error handling absorbs.

profile_call_args <- function(do_par, idx = 3) {

  recorder <- new.env()
  recorder$calls <- list()
  stub <- function(func, parameters, map, random = NULL, silent = TRUE, ...) {
    recorder$calls <- c(recorder$calls, list(list(map = map$ln_M, values = as.vector(parameters$ln_M))))
    list(par = numeric(0), fn = function(p) stop("stub"), gr = function(p) stop("stub"),
         report = function(p) stop("stub"), env = list(last.par.best = numeric(0)))
  }

  parameters <- list(ln_M = array(log(0.2), dim = c(2, 3)))
  # positions 2 and 5 held fixed, positions 1 and 4 estimated as a single shared parameter
  mapping <- list(ln_M = factor(c("1", NA, "2", "1", NA, "3")))

  local_mocked_bindings(MakeADFun = stub, .package = "RTMB")
  # run the parallel branch in this process, so the stub and the mocks are visible to it
  local_mocked_bindings(plan = function(...) invisible(NULL), .package = "future")
  local_mocked_bindings(future_lapply = function(X, FUN, ...) lapply(X, FUN),
                        .package = "future.apply")

  utils::capture.output(suppressMessages(try(
    do_likelihood_profile(data = list(), parameters = parameters, mapping = mapping,
                          what = "ln_M", idx = idx, min_val = -2, max_val = -1.5, inc = 0.5,
                          do_par = do_par),
    silent = TRUE)))

  recorder$calls
}

for(branch in c("sequential", "parallel")) {

  do_par <- branch == "parallel"

  test_that(paste("the", branch, "branch hands RTMB the fitted map structure"), {
    args <- profile_call_args(do_par)
    expect_length(args, 2)

    for(call in args) {
      # positions the fitted model fixed are still fixed
      expect_true(all(is.na(call$map[c(2, 5)])))
      # positions the fitted model estimated as one parameter still share a level
      expect_identical(as.character(call$map[1]), as.character(call$map[4]))
      expect_false(is.na(call$map[1]))
      # the profile target is fixed
      expect_true(is.na(call$map[3]))
      # two free parameters, the same count the fitted model had once the target is fixed
      expect_identical(nlevels(droplevels(call$map)), 2L)
    }
  })

  test_that(paste("the", branch, "branch hands RTMB the profile value"), {
    args <- profile_call_args(do_par)

    expect_equal(args[[1]]$values[3], -2)
    expect_equal(args[[2]]$values[3], -1.5)
    # nothing outside the target moves off its fitted value
    expect_equal(args[[1]]$values[-3], rep(log(0.2), 5))
  })

  test_that(paste("the", branch, "branch rebuilds each grid value from the fitted map"), {
    # Rebuilding from the map the previous grid value left behind would drift, so the two
    # grid values have to produce the same map.
    args <- profile_call_args(do_par)

    expect_identical(as.character(args[[1]]$map), as.character(args[[2]]$map))
  })

} # end branch loop
