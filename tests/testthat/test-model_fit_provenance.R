# A saved fit used to carry no record of the package that built it, so an old results file
# could only be dated by which report names it lacked. fit_model now stamps every fit.

library(SPoRC)
library(testthat)

test_that("fit_provenance reports the running package and a usable commit field", {
  prov <- SPoRC:::fit_provenance()
  expect_named(prov, c("package_version", "commit_sha", "commit_source",
                       "r_version", "rtmb_version", "fit_time"))
  expect_equal(prov$package_version, as.character(packageVersion("SPoRC")))
  expect_true(is.na(prov$commit_sha) || grepl("^[0-9a-f]{40}$", prov$commit_sha))
  expect_true(prov$commit_source %in% c("github", "local_git", "unknown"))
  expect_true(xor(is.na(prov$commit_sha), prov$commit_source != "unknown")) # sha and source agree
  expect_s3_class(prov$fit_time, "POSIXct")
})

test_that("fit_model attaches the provenance record", {
  input <- suppressWarnings(suppressMessages(objective_setup_input()))
  obj <- fit_model(input$data, input$par, input$map, do_optim = FALSE, silent = TRUE)
  stamp <- obj$provenance[setdiff(names(obj$provenance), "fit_time")]
  fresh <- SPoRC:::fit_provenance()[names(stamp)]
  expect_equal(stamp, fresh) # same build, only the time stamp moves
})
