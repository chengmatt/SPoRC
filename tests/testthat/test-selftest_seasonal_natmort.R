# Replicate self test of seasonal M recovery. The OM runs 0.45 and 0.20 over seasons of 0.4 and 0.6, and
# each replicate refits both rates free from one shared start, so a recovered split is not the start.
#
# Two things checked: at near-zero error both rates come back, and the annual total is identified much
# better than the split, since only the within-year change in numbers separates the seasons.

library(SPoRC)
library(testthat)

seasonal_M_replicates <- function(M_true, sigmaR, idx_se, iss, n_reps, seed0 = 400) {
  seasdur <- seasonal_M_cfg$seasdur
  out <- vapply(seq_len(n_reps), function(i) {
    sim <- seasonal_M_sim(
      M_by_seas = M_true,
      seed = seed0 + i,
      sigmaR = sigmaR,
      idx_se = idx_se,
      iss = iss
    )
    il <- suppressWarnings(seasonal_M_input(
      sim, list(M_spec = "est_ln_M", M_seasblk_spec = list(1, 2)), sigmaR = sigmaR))
    il$par$ln_M[] <- log(0.30) # one shared start, away from both
    fit <- fit_model(
      il$data,
      il$par,
      il$map,
      random = NULL,
      silent = TRUE,
      do_optim = TRUE,
      newton_loops = 1
    )
    c(unique(as.vector(fit$rep$natmort[,,,1,,])),
      unique(as.vector(fit$rep$natmort[,,,2,,])),
      max(abs(fit$gr(fit$optim$par))))
  }, numeric(3))

  list(
    M1 = out[1, ],
    M2 = out[2, ],
    grad = out[3, ],
    annual = out[1, ] * seasdur[1] + out[2, ] * seasdur[2]
  )
}

mare <- function(x, truth) stats::median(abs(100 * (x - truth) / truth))


test_that("both seasonal rates are recovered when the data are near noiseless", {

  M_true <- c(0.45, 0.20)
  r <- seasonal_M_replicates(M_true, sigmaR = 0.02, idx_se = 0.01, iss = 5000, n_reps = 10)

  expect_true(all(r$grad < 1e-4))

  # finds what it's given, both seasons
  expect_lt(mare(r$M1, M_true[1]), 5)
  expect_lt(mare(r$M2, M_true[2]), 5)

  # and they stay apart rather than collapsing onto the shared start
  expect_gt(stats::median(r$M1 - r$M2), 0.15)
})


test_that("at realistic noise the annual total is identified far better than the split", {

  M_true <- c(0.45, 0.20)
  seasdur <- seasonal_M_cfg$seasdur
  annual_true <- sum(M_true * seasdur)
  r <- seasonal_M_replicates(M_true, sigmaR = 0.3, idx_se = 0.1, iss = 300, n_reps = 10)

  expect_true(all(r$grad < 1e-4))

  # duration weighted total lands close and near unbiased
  expect_lt(mare(r$annual, annual_true), 6)
  expect_lt(abs(stats::median(100 * (r$annual - annual_true) / annual_true)), 3)

  # each seasonal rate is worse than the total it sums to. The gap is expected,
  # it's what makes an estimated split hard to defend
  expect_gt(mare(r$M1, M_true[1]), mare(r$annual, annual_true))
  expect_gt(mare(r$M2, M_true[2]), mare(r$annual, annual_true))

  # errors run opposite ways, i.e. trading between seasons while holding the
  # total, rather than moving as a pair
  expect_lt(stats::median(r$M1 - M_true[1]) * stats::median(r$M2 - M_true[2]), 0)
})
