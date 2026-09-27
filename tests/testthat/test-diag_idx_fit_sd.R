library(SPoRC)
library(testthat)

# get_idx_fits() intervals must use the standard deviation the index likelihood used, which
# holds the estimated component on top of the reported errors. Exact, no tolerance.

idx_inputs <- function(n_pop = 1, n_regions = 1, n_yrs = 3, n_seas = 1, n_fleets = 1) {
  arr4 <- function(x) array(x, dim = c(n_regions, n_yrs, n_seas, n_fleets))
  arr5 <- function(x) array(x, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fleets))
  list(
    data = list(
      ObsSrvIdx = arr4(10), ObsSrvIdx_SE = arr4(0.2),
      ObsFishIdx = arr4(5), ObsFishIdx_SE = arr4(0.1),
      UseSrvIdx = arr4(1), UseFishIdx = arr4(1),
      ObsSrvIdx_pop = arr5(10), ObsSrvIdx_pop_SE = arr5(0.2),
      ObsFishIdx_pop = arr5(5), ObsFishIdx_pop_SE = arr5(0.1),
      UseSrvIdx_pop = arr5(1), UseFishIdx_pop = arr5(1),
      Wt_SrvIdx = 1, Wt_FishIdx = 1, Wt_SrvIdx_pop = 1, Wt_FishIdx_pop = 1,
      srv_q_blocks = array(1, dim = c(n_regions, n_yrs, n_fleets)),
      fish_q_blocks = array(1, dim = c(n_regions, n_yrs, n_fleets)),
      years = seq_len(n_yrs) + 2000
    ),
    rep = list(PredSrvIdx = arr5(10), PredFishIdx = arr5(6)),
    arr4 = arr4, arr5 = arr5
  )
}

test_that("a report without the estimated component leaves the reported errors alone", {
  inp <- idx_inputs()
  out <- SPoRC::get_idx_fits(inp$data, inp$rep, inp$data$years)
  expect_equal(unique(out$se[out$Type == "Survey"]), 0.2)
  expect_equal(unique(out$se[out$Type == "Fishery"]), 0.1)
})

test_that("the estimated component widens the intervals on every index", {
  inp <- idx_inputs()
  inp$rep <- c(inp$rep, list(SrvIdx_SD = inp$arr4(0.35), FishIdx_SD = inp$arr4(0.15),
                             SrvIdx_pop_SD = inp$arr5(0.4), FishIdx_pop_SD = inp$arr5(0.25)))
  out <- SPoRC::get_idx_fits(inp$data, inp$rep, inp$data$years)

  expect_equal(unique(out$se[out$Type == "Survey"]), 0.35)
  expect_equal(unique(out$se[out$Type == "Fishery"]), 0.15)
  expect_equal(unique(out$se[out$Type == "Pop Survey"]), 0.4)
  expect_equal(unique(out$se[out$Type == "Pop Fishery"]), 0.25)
  expect_equal(out$lci, exp(log(out$obs) - 1.96 * out$se))
  expect_equal(out$uci, exp(log(out$obs) + 1.96 * out$se))
})

test_that("the likelihood weight still divides the total standard deviation", {
  inp <- idx_inputs()
  inp$data$Wt_SrvIdx <- 4
  inp$rep <- c(inp$rep, list(SrvIdx_SD = inp$arr4(0.35)))
  out <- SPoRC::get_idx_fits(inp$data, inp$rep, inp$data$years)
  expect_equal(unique(out$se[out$Type == "Survey"]), 0.35 / 2)
})
