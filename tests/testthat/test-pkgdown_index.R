# The pkgdown index against what the package actually contains. Every exported topic and vignette has
# to appear in _pkgdown.yml or the site build fails.
#
# Nothing else depends on that file, so it goes stale whenever a feature adds an export or a vignette
# and the failure only shows up in CI after a push. This moves the check to where the change is made.

test_that("every exported topic and vignette is in the pkgdown index", {
  skip_if_not_installed("pkgdown")
  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  skip_if_not(file.exists(file.path(pkg_root, "_pkgdown.yml")), "not running from a source tree")

  pkg <- pkgdown::as_pkgdown(pkg_root)

  # a vignette still being written is counted here rather than in _pkgdown.yml, so the rest of the
  # check keeps running. take it off this list when it goes into the file
  unlisted_vignettes <- character(0)
  pkg$meta$articles[[1]]$contents <- c(pkg$meta$articles[[1]]$contents, unlisted_vignettes)

  # each index is checked on its own, since check_pkgdown() stops at the first one that is short and
  # the articles come before the reference topics. the condition holds the missing names
  for(index in c("data_articles_index", "data_reference_index")) {
    result <- tryCatch({
      utils::getFromNamespace(index, "pkgdown")(pkg)
      "ok"
    }, error = function(e) conditionMessage(e))
    expect_equal(result, "ok")
  } # end index loop
})
