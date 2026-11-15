# The operating model and the estimation model follow tagged fish through the same season steps
# (get_tag_cohorts_at_liberty, start_tag_cohort_season, move_tag_cohorts, end_tag_cohort_season), so their
# predicted recaptures agree under every move timing, for releases part way through a season, and with
# state-space numbers at age, which tagged fish take as the rest of the population does.

test_that("operating and estimation model tag recaptures agree under every move timing, with mid-season releases", {
  for(mt in 0:2) {
    il <- suppressWarnings(suppressMessages(build_em(suppressWarnings(suppressMessages(build_om(move_timing = mt))), move_timing = mt)))
    il$data$conv_tag_t_tagging[] <- 0.5 # released half way through the season
    diffs <- lockstep_diffs(lockstep_om(il))
    expect_lt(diffs[["TagRecap"]], 1e-10, label = paste("tag recaptures, move_timing", mt))
  } # end mt loop
})

test_that("operating and estimation model tag recaptures agree with state-space numbers at age", {
  for(mt in c(0, 2)) {
    il <- suppressWarnings(suppressMessages(build_em(suppressWarnings(suppressMessages(build_om(move_timing = mt))), move_timing = mt, NAA_re = "iid")))
    run <- lockstep_om(il, random = "ln_NAA")
    expect_gt(diff(range(run$rep$NAA_scalar)), 0.1) # the state moves numbers at age, so it would show in the tags
    expect_lt(lockstep_diffs(run)[["TagRecap"]], 1e-10, label = paste("tag recaptures, move_timing", mt))
  } # end mt loop
})
