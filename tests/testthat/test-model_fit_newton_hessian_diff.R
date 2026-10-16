# Checks the differenced-Hessian refinement: fit_model(newton_hessian = "diff") differences the exact gradient for the Hessian and
# reruns nlminb with each parameter in its own standard errors. It should reach the exact-Hessian Newton fit's optimum, and return a
# Hessian close to the tape's for sdreport.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

data <- dusky_rtmb_model$data
parameters <- dusky_rtmb_model$parameters
mapping <- dusky_rtmb_model$mapping

test_that("newton_hessian = 'diff' reaches the exact-Hessian optimum and returns its Hessian", {

  exact <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 2, silent = TRUE)
  diff <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 2, newton_hessian = "diff", silent = TRUE)

  expect_equal(diff$optim$objective, exact$optim$objective, tolerance = 1e-8)
  expect_lt(max(abs(diff$gr(diff$optim$par))), 1e-3)
  expect_equal(dim(diff$hessian), rep(length(diff$par), 2))
  expect_equal(unname(diff$hessian), unname(as.matrix(diff$he(diff$optim$par))), tolerance = 1e-3)
  expect_error(fit_model(data, parameters, mapping, newton_loops = 0, newton_hessian = "tape", silent = TRUE), "exact or diff")

})
