library(SPoRC)
library(testthat)
data("dusky_rtmb_model")

# The bundled dusky object carries a $rep snapshot taken when it was built, so the
# objective is re-evaluated here from its stored data and parameters to make sure
# every check below reflects the current SPoRC_rtmb rather than that snapshot.
mle_pars <- dusky_rtmb_model$env$parList()

evaluate_at_mle <- function(data) {
  fit_model(data, mle_pars, dusky_rtmb_model$mapping,
            random = NULL, silent = TRUE, do_optim = FALSE)
}


test_that("jnLL equals the weighted sum of its reported likelihood components", {
  baseline <- evaluate_at_mle(dusky_rtmb_model$data)
  expect_jnLL_decomposes(baseline, label = "dusky")

  # the fixture has to actually exercise several components, or the check above
  # passes for the uninteresting reason that everything is zero
  contributions <- jnLL_contributions(baseline)
  expect_gt(sum(contributions$contribution != 0), 5)
})


test_that("expect_jnLL_decomposes fails on a decomposition that does not hold", {
  baseline <- evaluate_at_mle(dusky_rtmb_model$data)

  # a component dropped from, or double counted in, the jnLL sum
  mismatched <- baseline
  mismatched$rep$jnLL <- baseline$rep$jnLL + baseline$rep$Catch_nLL[1]
  expect_failure(expect_jnLL_decomposes(mismatched, label = "mismatched"))

  # a new likelihood reported but missing from the list of jnLL terms
  unaccounted <- baseline
  unaccounted$rep$Newthing_nLL <- 1
  expect_failure(expect_jnLL_decomposes(unaccounted, label = "unaccounted"))

  # a likelihood still listed as a jnLL term but no longer reported
  stale <- baseline
  stale$rep$sel_nLL <- NULL
  expect_failure(expect_jnLL_decomposes(stale, label = "stale"))
})


test_that("a term that resolves to more than one contribution is an error, not extra rows", {
  baseline <- evaluate_at_mle(dusky_rtmb_model$data)

  # Wt_F is a scalar-mode term, so an array weight makes Wt_F * sum(Fmort_nLL)
  # return one value per weight element. Without the length check those become one
  # row each and the component is counted once per element, which reads as a
  # decomposition that misses by a multiple of itself rather than as a wrong mode.
  vector_weight <- baseline
  vector_weight$data$Wt_F <- c(1, 1)
  expect_error(jnLL_contributions(vector_weight), "Fmort_nLL")
  expect_error(jnLL_contributions(vector_weight), "scalar")
})


test_that("each likelihood weight scales its own component and leaves the rest untouched", {
  baseline <- evaluate_at_mle(dusky_rtmb_model$data)
  contributions <- jnLL_contributions(baseline)
  components <- grep("_nLL$", names(baseline$rep), value = TRUE)

  # only weights whose components are non-zero here, so that doubling a weight
  # produces a difference the test can actually detect
  active <- unique(stats::na.omit(contributions$weight[contributions$contribution != 0]))
  expect_gt(length(active), 2)

  for(wt_name in active) {
    perturbed_data <- dusky_rtmb_model$data
    perturbed_data[[wt_name]] <- perturbed_data[[wt_name]] * 2
    perturbed <- evaluate_at_mle(perturbed_data)

    # weights enter only the jnLL sum, so no likelihood component itself may move
    for(nm in components) {
      expect_equal(perturbed$rep[[nm]], baseline$rep[[nm]],
                   info = paste0("doubling ", wt_name, " changed the ", nm, " component"))
    }

    # doubling a weight adds exactly one more copy of what it already contributed,
    # which for Wt_Rec covers both Rec_nLL and Init_Rec_nLL
    expected_delta <- sum(contributions$contribution[which(contributions$weight == wt_name)])
    expect_equal(perturbed$rep$jnLL - baseline$rep$jnLL, expected_delta,
                 tolerance = 1e-8,
                 info = paste0(wt_name, " is not applied to the likelihoods jnLL_terms assigns it"))
  }
})
