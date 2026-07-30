library(SPoRC)
library(testthat)

test_that("build_pe_map assigns unique ids by collapsing shared dims", {

  dims <- c(region = 2, season = 3, fleet = 2)

  test_that("est_all (no sharing) gives a unique id per cell", {
    map <- SPoRC:::build_pe_map(dims, share_over = character(0))
    expect_equal(unname(dim(map)), unname(dims))
    expect_equal(sort(unique(as.vector(map))), 1:prod(dims))
  })

  test_that("sharing over one dim collapses ids across that dim only", {
    map <- SPoRC:::build_pe_map(dims, share_over = "region")
    expect_equal(length(unique(as.vector(map))), unname(dims["season"] * dims["fleet"]))
    # cells differing only in region must share the same id
    expect_equal(map[1, 2, 1], map[2, 2, 1])
    # cells differing in season must NOT share an id
    expect_false(map[1, 1, 1] == map[1, 2, 1])
  })

  test_that("sharing over all dims collapses to a single id", {
    map <- SPoRC:::build_pe_map(dims, share_over = names(dims))
    expect_equal(length(unique(as.vector(map))), 1)
  })

  test_that("errors on unnamed dims or invalid share_over", {
    expect_error(SPoRC:::build_pe_map(c(2, 3)), "fully named")
    expect_error(SPoRC:::build_pe_map(dims, share_over = "bogus"), "subset of names")
  })
})

test_that("build_shared_spec_map parses est_all/fix/est_shared_* spec strings", {

  dims <- c(r = 2, seas = 3, f = 2)
  names(dims) <- c("region", "season", "fleet")
  abbrev <- c(r = "region", seas = "season", f = "fleet")

  test_that("est_all gives a unique id per cell", {
    m <- SPoRC:::build_shared_spec_map(dims, "est_all", abbrev)
    expect_equal(length(levels(m)), prod(dims))
    expect_true(all(!is.na(m)))
  })

  test_that("fix maps everything to NA", {
    m <- SPoRC:::build_shared_spec_map(dims, "fix", abbrev)
    expect_true(all(is.na(m)))
  })

  test_that("est_shared_<dims> matches build_pe_map with the same share_over", {
    m1 <- SPoRC:::build_shared_spec_map(dims, "est_shared_r_seas", abbrev)
    m2 <- SPoRC:::build_pe_map(dims, share_over = c("region", "season"))
    expect_equal(as.integer(m1), as.integer(factor(as.vector(m2))))
  })

  test_that("invalid spec errors with the list of valid options", {
    expect_error(SPoRC:::build_shared_spec_map(dims, "est_shared_bogus", abbrev), "not recognized")
  })
})
