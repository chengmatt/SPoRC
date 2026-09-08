library(SPoRC)
library(testthat)

# A free initial age structure takes the numbers at age 2 and older from ln_InitDevs outright, so
# no equilibrium is projected and init_F_par never reaches the objective.

init_F_rec <- function(...) {
  sweep_input(rec = utils::modifyList(list(init_F_form = "abs", init_F_par = array(log(0.05), dim = c(3, 1, 5))),
                                      list(...)),
              stop_after = "rec")
}


test_that("estimating the initial F under a free initialization is refused", {

  msg <- tryCatch(init_F_rec(init_age_strc = "free", init_F_spec = "est"),
                  error = function(e) conditionMessage(e))

  expect_match(msg, "init_F_spec = 'est' with init_age_strc = 'free'")
  # the message must say which setting to change, not only that the pairing is wrong
  expect_match(msg, "init_F_spec = 'fix'")
})


test_that("the sensible pairings still build", {

  # fixed under a free initialization: the parameter is there and mapped off
  free_fix <- init_F_rec(init_age_strc = "free", init_F_spec = "fix")
  expect_true(all(is.na(free_fix$map$init_F_par)))

  # estimated under an equilibrium initialization, where it sets the age structure
  equil_est <- init_F_rec(init_age_strc = "scalar_no_move", init_F_spec = "est")
  expect_false(any(is.na(equil_est$map$init_F_par)))
})


test_that("a free initialization ignores the initial F, and an equilibrium one does not", {

  jnLL <- function(init_age_strc, init_F) {
    il <- sweep_input(rec = list(init_age_strc = init_age_strc, init_F_form = "abs",
                                 init_F_spec = "fix",
                                 init_F_par = array(log(init_F), dim = c(3, 1, 5))))
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    obj$fn(obj$par)
  }

  expect_equal(jnLL("free", 0.02), jnLL("free", 0.4), tolerance = 0)
  expect_false(isTRUE(all.equal(jnLL("scalar_no_move", 0.02), jnLL("scalar_no_move", 0.4))))
})
