library(SPoRC)
library(testthat)

# get_model_rep_from_mcmc() replays the model report at every posterior draw.
# A real ADFun is not needed to exercise it: the function only calls
# rtmb_obj$report(par = ), so a stand-in that echoes its input pins the warmup
# handling, the chain collapse, and the draw ordering exactly.

# The function sets a global future plan and never restores it, so tests wrap
# their calls to avoid leaking a multisession plan into the rest of the suite.
with_restored_plan <- function(code) {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  force(code)
}

# samples[iter, chain, param]; p1 holds 1:10, p2 holds 11:20.
make_mcmc <- function(n_iter = 5, n_chain = 2, warmup = 2) {
  n_param <- 2
  smp <- array(seq_len(n_iter * n_chain * n_param),
               dim = c(n_iter, n_chain, n_param),
               dimnames = list(NULL, NULL, c("p1", "p2")))
  list(samples = smp, warmup = warmup)
}

# Echoes the draw back so the values landing in the report can be checked.
echo_obj <- list(report = function(par) list(p1 = par[["p1"]], p2 = par[["p2"]]))

test_that("warmup draws are discarded before the report is replayed", {
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(n_iter = 5, n_chain = 2, warmup = 2),
                                   what = "p1", n_cores = 1))
  # 3 post-warmup iterations across 2 chains.
  expect_equal(nrow(out$p1), 3 * 2)
  # Warmup values 1, 2 (chain 1) and 6, 7 (chain 2) must not appear.
  expect_false(any(out$p1$value %in% c(1, 2, 6, 7)))
  expect_setequal(out$p1$value, c(3, 4, 5, 8, 9, 10))
})

test_that("a different warmup length changes the number of draws", {
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(n_iter = 5, n_chain = 2, warmup = 1),
                                   what = "p1", n_cores = 1))
  expect_equal(nrow(out$p1), 4 * 2)
  expect_setequal(out$p1$value, c(2, 3, 4, 5, 7, 8, 9, 10))
})

test_that("chains are collapsed with iterations varying fastest", {
  # Draw order matters for anyone joining these tables back to the samples.
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(), what = "p1", n_cores = 1))
  ordered <- out$p1[order(out$p1$posterior_sample), ]
  expect_equal(ordered$posterior_sample, 1:6)
  expect_equal(ordered$value, c(3, 4, 5, 8, 9, 10))
})

test_that("parameter names from the samples dimnames reach the report function", {
  # The report is indexed by name inside echo_obj, so this would error if the
  # column names were not kept onto the collapsed matrix.
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(), what = c("p1", "p2"), n_cores = 1))
  expect_setequal(out$p1$value, c(3, 4, 5, 8, 9, 10))
  expect_setequal(out$p2$value, c(13, 14, 15, 18, 19, 20))
})

test_that("every requested report component is returned and named", {
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(), what = c("p1", "p2"), n_cores = 1))
  expect_named(out, c("p1", "p2"))
  expect_length(out, 2)
})

test_that("components not requested are dropped", {
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(), what = "p2", n_cores = 1))
  expect_named(out, "p2")
})

test_that("an array-valued report component keeps its index columns", {
  # A 2x3 report quantity melts to Var1/Var2 index columns per draw.
  obj <- list(report = function(par) list(SSB = matrix(par[["p1"]], nrow = 2, ncol = 3)))
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(obj, make_mcmc(), what = "SSB", n_cores = 1))

  expect_true(all(c("Var1", "Var2", "value", "posterior_sample") %in% names(out$SSB)))
  expect_equal(nrow(out$SSB), 6 * 2 * 3)
  expect_setequal(out$SSB$Var1, 1:2)
  expect_setequal(out$SSB$Var2, 1:3)
  # Each draw contributes one full copy of the array.
  expect_equal(as.vector(table(out$SSB$posterior_sample)), rep(6, 6))
})

test_that("the result is a data.table per component", {
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(), what = "p1", n_cores = 1))
  expect_s3_class(out$p1, "data.table")
})

test_that("a single-chain posterior currently fails", {
  # Dropping warmup uses the default drop = TRUE, so a one-chain posterior
  # loses its chain dimension, n_param reads as NA, and the matrix() call
  # errors. Pinned as current behavior: single-chain runs are not supported.
  expect_error(with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(n_iter = 5, n_chain = 1, warmup = 2),
                                   what = "p1", n_cores = 1)),
    "ncol")
})

test_that("a single post-warmup iteration currently fails for the same reason", {
  expect_error(with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(n_iter = 3, n_chain = 2, warmup = 2),
                                   what = "p1", n_cores = 1)),
    "ncol")
})

test_that("two or more chains and iterations work", {
  out <- with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(echo_obj, make_mcmc(n_iter = 4, n_chain = 2, warmup = 2),
                                   what = "p1", n_cores = 1))
  expect_equal(nrow(out$p1), 2 * 2)
})

test_that("the report is evaluated once per draw", {
  # Counts calls via a file, since the workers run in separate processes.
  tally <- tempfile()
  obj <- list(report = function(par) {
    cat("1\n", file = tally, append = TRUE)
    list(p1 = par[["p1"]])
  })
  with_restored_plan(
    SPoRC::get_model_rep_from_mcmc(obj, make_mcmc(), what = "p1", n_cores = 1))
  expect_equal(length(readLines(tally)), 6)
})
