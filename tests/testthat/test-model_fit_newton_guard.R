# Checks that a Newton step is only kept when it lowers the largest absolute gradient: with an indefinite Hessian the step could
# move off the nlminb optimum to a worse point, so refinement stops instead and the gradient never ends above where it started.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

data <- dusky_rtmb_model$data
parameters <- dusky_rtmb_model$parameters
mapping <- dusky_rtmb_model$mapping

test_that("Newton refinement never leaves a larger gradient than nlminb's", {

  start <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 0, silent = TRUE)
  refined <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 3, silent = TRUE)

  expect_lte(max(abs(refined$gr(refined$optim$par))), max(abs(start$gr(start$optim$par))))
  expect_lte(refined$optim$objective, start$optim$objective + 1e-8)

})
