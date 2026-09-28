# The generated half of the options documentation.
#
# ta_option_reference.Rmd is built by reading the package rather than describing it, which is what makes it
# complete. A change to how arguments or help pages are structured could leave it returning fewer rows.
#
# These check that every argument is present and that every argument says something.

test_that("the reference covers every argument of every setup stage", {
  ref <- option_reference()

  for(s in setup_stage_order()) {
    want <- setdiff(names(formals(getExportedValue("SPoRC", s))),
                    c("input_list", "...", "verbose"))
    got <- ref$argument[ref$stage == sub("^Setup_Mod_", "", s)]
    expect_setequal(got, want)
  }
})


test_that("every argument has its own documentation", {
  # the description comes from the help page, so an empty one is either an undocumented
  # argument or a change in how the help is read, and the vignette would render blank.
  #
  # an installed package without its help index has nothing to read and no source man/ to
  # fall back on, so every description is blank for a reason that says nothing
  skip_if(length(SPoRC:::rd_database()) == 0, "the package Rd database is not available")

  ref <- option_reference()
  blank <- ref[!nzchar(ref$description), c("stage", "argument")]

  expect_equal(nrow(blank), 0,
               label = if(nrow(blank)) paste("arguments with no description:",
                                             paste(blank$argument, collapse = ", ")) else "none")
})


test_that("the reference reports which settings the written guide also discusses", {
  # The guide is hand-written and will always trail the code; the point of
  # reporting the gap is to keep it visible rather than to close it here.
  guide <- testthat::test_path("..", "..", "vignettes", "t_model_options.Rmd")
  skip_if_not(file.exists(guide), "the options guide is not beside the tests")

  ref <- option_reference(guide = guide)
  expect_type(ref$in_guide, "logical")
  expect_false(any(is.na(ref$in_guide)))
  # some are discussed and some are not, or the check is not reading the guide
  expect_gt(sum(ref$in_guide), 0)
  expect_gt(sum(!ref$in_guide), 0)
})


test_that("the stages are listed in the order a model is built", {
  # the vignette groups by this order, so a reader going top to bottom follows the stages.
  # dimensions first and weighting last are the two ends that mean anything
  o <- setup_stage_order()

  expect_equal(o[1], "Setup_Mod_Dim")
  expect_equal(o[length(o)], "Setup_Mod_Weighting")
  expect_setequal(o, grep("^Setup_Mod_", getNamespaceExports("SPoRC"), value = TRUE))
})
