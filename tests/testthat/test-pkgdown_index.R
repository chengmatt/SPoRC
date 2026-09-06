# The pkgdown index against what the package actually contains. Every exported topic and vignette has
# to appear in _pkgdown.yml or the site build fails.
#
# Nothing else depends on that file, so it goes stale whenever a feature adds an export or a vignette
# and the failure only shows up in CI after a push. This moves the check to where the change is made.

test_that("every exported topic and vignette is in the pkgdown index", {
  skip_if_not_installed("pkgdown")
  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  skip_if_not(file.exists(file.path(pkg_root, "_pkgdown.yml")), "not running from a source tree")

  # check_pkgdown() aborts on the first category that is short, so its condition holds the
  # missing names: reporting the message is what makes the failure actionable
  result <- tryCatch({ pkgdown::check_pkgdown(pkg_root); "ok" },
                     error = function(e) conditionMessage(e))
  expect_equal(result, "ok")
})
