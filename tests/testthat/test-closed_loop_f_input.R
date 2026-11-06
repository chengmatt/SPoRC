# fishing and discard mortality supplied to condition_closed_loop_simulations(), which takes both as
# overrides through ... but until 2026-11-05 built them from the fit with zeros after it whatever was passed

data("dusky_rtmb_model")

test_that("supplied fishing and discard mortality reach the returned sim_list", {
  m <- dusky_rtmb_model
  closed_loop_yrs <- 3
  n_sims <- 2
  f_dim <- c(m$data$n_regions, length(m$data$years) + closed_loop_yrs, m$data$n_seas,
             m$data$n_fish_fleets, n_sims) # region, year, season, fleet, sim
  Fmort <- array(0.001 * seq_len(prod(f_dim)), dim = f_dim) # a different value in every cell, so the fit's array cannot match
  dmr <- array(0.25 + 0.001 * seq_len(prod(f_dim)), dim = f_dim)

  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = n_sims,
    m$data,
    m$parameters,
    m$mapping,
    sd_rep = m$sdrep,
    rep = m$rep,
    random = NULL,
    Fmort_input = Fmort,
    dmr_input = dmr
  )

  expect_identical(cl$Fmort, Fmort)
  expect_identical(cl$dmr, dmr)
})


test_that("without them the projection years are unfished and the conditioning years are the fit's", {
  m <- dusky_rtmb_model
  closed_loop_yrs <- 3
  n_hist <- length(m$data$years)
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = 2,
    m$data,
    m$parameters,
    m$mapping,
    sd_rep = m$sdrep,
    rep = m$rep,
    random = NULL
  )

  expect_equal(dim(cl$Fmort), c(m$data$n_regions, n_hist + closed_loop_yrs, m$data$n_seas, m$data$n_fish_fleets, 2)) # region, year, season, fleet, sim
  expect_equal(cl$Fmort[, seq_len(n_hist), , , 2], m$rep$Fmort[, seq_len(n_hist), , ])
  expect_true(all(cl$Fmort[, n_hist + seq_len(closed_loop_yrs), , , ] == 0))
})


test_that("a supplied fishing mortality fishes the projection years in the operating model", {
  # with the override ignored the projection years ran at zero F, so their catch was zero
  m <- dusky_rtmb_model
  closed_loop_yrs <- 3
  n_hist <- length(m$data$years)
  proj <- n_hist + seq_len(closed_loop_yrs) # projection years
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = 1,
    m$data,
    m$parameters,
    m$mapping,
    sd_rep = m$sdrep,
    rep = m$rep,
    random = NULL
  )
  Fmort <- cl$Fmort
  Fmort[, proj, , , ] <- Fmort[, n_hist, , , ] # the last fitted year's F carried on
  cl_f <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs,
    n_sims = 1,
    m$data,
    m$parameters,
    m$mapping,
    sd_rep = m$sdrep,
    rep = m$rep,
    random = NULL,
    Fmort_input = Fmort
  )

  set.seed(1)
  out <- Simulate_Pop_Static(cl_f)
  expect_true(all(out$TrueCatch[, proj, , , ] > 0))
})
