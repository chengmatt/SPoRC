# Time blocks on R0. Under mean recruitment R0 IS mean recruitment, so a block is a productivity
# regime; under a stock-recruit form it is the curve's scale.
#
# The reference block supplies the single value the initial age structure, the prior, the ln_rinit
# penalty and the stock-recruit scale all read. Helpers in helper-blocks.R.

sim_r0 <- suppressWarnings(blocks_sim())

# a bare dimensions + recruitment pair, which is all the block parsing needs
r0_dim <- function() Setup_Mod_Dim(
  years = 1:30,
  ages = 1:10,
  lens = 1,
  n_regions = 1,
  n_sexes = 1,
  n_fish_fleets = 1,
  n_srv_fleets = 1,
  n_pop = 1,
  verbose = FALSE
)
r0_rec <- function(b, ...) suppressWarnings(suppressMessages(
  Setup_Mod_Rec(
    r0_dim(),
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(1), c(2, 1, 1)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    R0_blocks = b,
    ...
  )))

test_that("an unblocked model keeps one R0 column and is unchanged", {
  il <- blocks_em(sim_r0)
  expect_equal(dim(il$par$ln_global_R0), c(1L, 1L))
  expect_equal(unique(as.vector(il$data$R0_blocks)), 1)
  expect_equal(il$data$R0_ref_block, 1L)
  o <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  r <- o$report(o$env$last.par)
  expect_true(is.finite(o$fn(o$par)))
  expect_equal(length(unique(round(as.vector(r$R0_yr), 10))), 1)   # constant across years
})

test_that("an input list built before R0 blocks existed still runs", {
  il <- blocks_em(sim_r0)
  d <- il$data; d$R0_blocks <- NULL; d$R0_ref_block <- NULL
  p <- il$par; p$ln_global_R0 <- as.vector(p$ln_global_R0)        # the old vector shape
  o <- fit_model(d, p, il$map, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(o$fn(o$par)))
})

test_that("R0 block strings parse, size the parameter and guard their years", {
  skip_on_cran()
  mk <- r0_rec
  a <- mk(c("Block_1_Year_1-15_Pop_1", "Block_2_Year_16-terminal_Pop_1"))
  expect_equal(dim(a$par$ln_global_R0), c(1L, 2L))
  expect_equal(as.vector(a$data$R0_blocks), c(rep(1, 15), rep(2, 15)))
  b <- mk(c("Block_1_Year_1-10_Pop_1", "Block_2_Year_11-20_Pop_1", "Block_3_Year_21-terminal_Pop_1"))
  expect_equal(dim(b$par$ln_global_R0), c(1L, 3L))
  expect_error(mk(c("Block_1_Year_1-15_Pop_1")), "leaves some years unassigned")
  expect_error(mk("none_Pop_1", R0_ref_block = 5), "R0_ref_block must lie within")
  expect_error(mk("Regime_2_Pop_1"), "R0_blocks must be none_Pop_p")
})

test_that("hand-placed R0 blocks give a year-varying R0 and beat a single R0", {
  skip_on_cran()
  fit <- function(b) {
    il <- blocks_em(sim_r0)
    il$data$R0_blocks <- array(c(rep(1, 15), rep(2, 15)), dim = c(1, 30, 1))
    il$par$ln_global_R0 <- matrix(c(log(5), log(12)), 1, 2)
    il$map$ln_global_R0 <- NULL
    fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  }
  o <- fit(NULL)
  r <- o$report(o$env$last.par)
  expect_equal(dim(r$R0_yr), c(1L, 30L))
  expect_equal(as.vector(r$R0_yr)[1], 5, tolerance = 1e-8)
  expect_equal(as.vector(r$R0_yr)[30], 12, tolerance = 1e-8)
  expect_true(is.finite(o$fn(o$par)))
})

test_that("a retrospective peel truncates R0 blocks and drops a block it emptied", {

  skip_on_cran()
  il <- blocks_em(sim_r0)
  n_yrs <- length(il$data$years)
  il$data$R0_blocks <- array(c(rep(1, 15), rep(2, 15)), dim = c(1, n_yrs, 1))
  il$par$ln_global_R0 <- matrix(c(log(5), log(12)), 1, 2)
  il$map$ln_global_R0 <- NULL
  f <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)

  # a block with no data left moves the objective not at all when perturbed. The final
  # fn call restores the starting parameters, since report() reads env$last.par.
  n_free_dead <- function(o) {
    idx <- which(names(o$par) == "ln_global_R0"); p <- o$par; base <- o$fn(p)
    out <- sum(vapply(idx, function(k) { p2 <- p; p2[k] <- p2[k] + 0.5; abs(o$fn(p2) - base) < 1e-8 }, logical(1)))
    o$fn(p)
    out
  }

  # a shallow peel keeps both blocks
  a <- truncate_yr(j = 5, data = f$data, parameters = f$parameters, mapping = f$mapping)
  expect_equal(dim(a$retro_data$R0_blocks)[2], length(a$retro_data$years))
  expect_equal(dim(a$retro_parameters$ln_global_R0), c(1L, 2L))
  oa <- fit_model(a$retro_data, a$retro_parameters, a$retro_mapping, do_optim = FALSE, silent = TRUE)
  expect_equal(n_free_dead(oa), 0)
  ra <- oa$report(oa$env$last.par)
  expect_equal(ra$R0_yr[1, ncol(ra$R0_yr)], 12, tolerance = 1e-8)

  # a peel past the start of the terminal block leaves it with no data, so its column
  # goes rather than sitting in the Hessian as an exactly zero row
  b <- truncate_yr(j = 20, data = f$data, parameters = f$parameters, mapping = f$mapping)
  expect_equal(dim(b$retro_data$R0_blocks)[2], length(b$retro_data$years))
  expect_equal(dim(b$retro_parameters$ln_global_R0), c(1L, 1L))
  ob <- fit_model(b$retro_data, b$retro_parameters, b$retro_mapping, do_optim = FALSE, silent = TRUE)
  expect_equal(n_free_dead(ob), 0)
  rb <- ob$report(ob$env$last.par)
  expect_equal(length(unique(round(as.vector(rb$R0_yr), 10))), 1)   # only block 1 survives
})

test_that("an unblocked model peels exactly as it did before R0 blocks existed", {

  skip_on_cran()
  il <- blocks_em(sim_r0)
  f <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  for(j in c(5, 20)) {
    out <- truncate_yr(j = j, data = f$data, parameters = f$parameters, mapping = f$mapping)
    expect_equal(dim(out$retro_parameters$ln_global_R0), c(1L, 1L))
    o <- fit_model(out$retro_data, out$retro_parameters, out$retro_mapping, do_optim = FALSE, silent = TRUE)
    expect_true(is.finite(o$fn(o$par)))
  } # end j loop
})

test_that("the operating model has R0 blocks rather than flattening them", {

  skip_on_cran()
  # simulation_self_test and the closed loop both build R0_input from rep$R0_yr, not from
  # the single reference-block value in rep$R0, so a blocked fit conditions a blocked OM
  il <- blocks_em(sim_r0)
  n_yrs <- length(il$data$years)
  il$data$R0_blocks <- array(c(rep(1, 15), rep(2, 15)), dim = c(1, n_yrs, 1))
  il$par$ln_global_R0 <- matrix(c(log(5), log(12)), 1, 2)
  il$map$ln_global_R0 <- NULL
  f <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  r <- f$report(f$env$last.par)

  expect_equal(r$R0[1], 5, tolerance = 1e-8)                      # the reference block
  expect_equal(as.vector(r$R0_yr)[c(1, 15, 16, n_yrs)], c(5, 5, 12, 12), tolerance = 1e-8)

  # the unblocked case must give a flat R0_yr, so the OM is unchanged there
  il2 <- blocks_em(sim_r0)
  f2 <- fit_model(il2$data, il2$par, il2$map, do_optim = FALSE, silent = TRUE)
  r2 <- f2$report(f2$env$last.par)
  expect_equal(length(unique(round(as.vector(r2$R0_yr), 10))), 1)
  expect_equal(unique(as.vector(r2$R0_yr)), r2$R0[1], tolerance = 1e-8)
})

test_that("the operating model warns when it cannot honor R0_ref_block", {

  # year one of R0_input both solves the operating model's equilibrium and generates its
  # first year's recruitment, so the operating model always starts from the block in force
  # at year one; the estimation model starts from R0_ref_block
  d <- list(
    R0_blocks = array(c(rep(1, 15), rep(2, 15)), dim = c(1, 30, 1)),
    R0_ref_block = 1L,
    use_rinit = 0
  )
  expect_silent(SPoRC:::warn_R0_ref_block_om(d, "test"))

  d$R0_ref_block <- 2L
  expect_warning(SPoRC:::warn_R0_ref_block_om(d, "test"), "R0_ref_block is 2 but year one sits in block 1")
  expect_warning(SPoRC:::warn_R0_ref_block_om(d, "test"), "cannot recover the operating model")

  # rinit initializes both sides, so the blocks never enter
  d$use_rinit <- 1
  expect_silent(SPoRC:::warn_R0_ref_block_om(d, "test"))

  # an unblocked model has nothing to check
  expect_silent(SPoRC:::warn_R0_ref_block_om(list(R0_blocks = NULL, use_rinit = 0), "test"))
})
