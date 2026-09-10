# The population-specific data sources against the objective. Before this, nothing in the suite fit
# population-specific catch, indices or compositions, so those blocks could break without a failure.

library(SPoRC)
library(testthat)

test_that("every population-specific data source contributes to the objective", {

  # every index likelihood, because the two index blocks split them differently: a
  # lognormal fleet is evaluated by eval_index_osa_nLL, a normal one by
  # get_index_pop_nLL regionally and get_index_regional_nLL by population. testing
  # one likelihood leaves the other route free to do nothing
  for(like in c("lognormal", "normal", "mvn")) {

    input_list <- suppressWarnings(suppressMessages(pop_sources_input(like = like)))

    # a source that is silently off reads as a passing test everywhere else, so check the
    # use flags before the likelihoods rather than trusting the setup
    for(flag in c("UseCatch_pop", "UseFishIdx_pop", "UseFishAgeComps_pop",
                  "UseSrvIdx_pop", "UseSrvAgeComps_pop"))
      expect_true(any(input_list$data[[flag]] == 1), info = paste(like, flag))

    nLL <- pop_sources_nLL(input_list)

    # each term is finite and actually moves the total. the non-zero half is what
    # catches a block that stopped evaluating rather than one that got the wrong answer
    for(term in names(nLL)) {
      expect_true(is.finite(nLL[[term]]), info = paste(like, term))
      expect_gt(abs(nLL[[term]]), 0)
    }

  } # end like loop
})


test_that("the objective and its gradient are finite under population-specific sources", {

  for(cfg in list(list(), list(n_seas = 2), list(n_pop = 3, n_regions = 3), list(like = "normal"))) {

    input_list <- suppressWarnings(suppressMessages(do.call(pop_sources_input, cfg)))
    obj <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)

    label <- paste(names(cfg), unlist(cfg), collapse = " ", sep = "=")
    expect_true(is.finite(obj$fn(obj$par)), info = label)
    expect_true(all(is.finite(obj$gr(obj$par))), info = label)

  } # end cfg loop
})


test_that("the regional index block evaluates the fleets the lognormal route skips", {

  # get_index_regional_nLL fits everything except lognormal, so normal and mvn are the
  # only likelihoods that reach it and the only ones that check it does anything
  for(like in c("normal", "mvn")) {
    nLL <- pop_sources_nLL(suppressWarnings(suppressMessages(pop_sources_input(like = like))))
    expect_gt(abs(nLL$FishIdx), 0)
    expect_gt(abs(nLL$SrvIdx), 0)
  } # end like loop
})


test_that("the two index routes give a population fleet different likelihoods", {

  # the regional block fits everything except lognormal, the population-specific block fits
  # only normal, so a fleet's likelihood decides which of the two evaluates it
  lognormal <- pop_sources_nLL(suppressWarnings(suppressMessages(pop_sources_input(like = "lognormal"))))
  normal <- pop_sources_nLL(suppressWarnings(suppressMessages(pop_sources_input(like = "normal"))))

  # both routes have to produce something, or "they differ" is satisfied by one of them being zero
  expect_gt(abs(lognormal$FishIdx_pop), 0)
  expect_gt(abs(normal$FishIdx_pop), 0)
  expect_gt(abs(lognormal$SrvIdx_pop), 0)
  expect_gt(abs(normal$SrvIdx_pop), 0)

  expect_false(isTRUE(all.equal(lognormal$FishIdx_pop, normal$FishIdx_pop)))
  expect_false(isTRUE(all.equal(lognormal$SrvIdx_pop, normal$SrvIdx_pop)))
})
