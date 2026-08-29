# What the setup functions say when an argument is wrong.
#
# An error message is part of the interface. These two were each wrong in a way
# that sent the reader somewhere unhelpful: a per-fleet setting given as a scalar
# was reported as an invalid value, and the value named in that report was one the
# same message listed as valid; and a composition type given in the form its own
# error listed was then rejected for its fleet. Both cost more time than the
# mistakes deserved, so both are pinned here.
#
# These check the message a caller actually sees, so they are written against the
# exported functions rather than against the internal validators.

err_msg <- function(expr) {
  tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
}

fleet_spec_input <- function() sweep_input(stop_after = "srvidx")


test_that("a per-fleet setting given once says so, rather than blaming its value", {
  # sweep_dims carries five fishery fleets, so a scalar is short by four.
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets
  base <- list(input_list = il, fish_sel_model = paste0("logist1_Fleet_", seq_len(n)))

  msg <- err_msg(do.call(Setup_Mod_Fishsel_and_Q, c(base, list(
    fish_q_spec = "est_all", fish_fixed_sel_pars_spec = rep("est_all", n)))))

  expect_match(msg, "fish_q_spec has 1 entry for 5 fleets")
  # the message carries the call that fixes it, not just the diagnosis
  expect_match(msg, 'rep\\("est_all", 5\\)', fixed = FALSE)
  # and it must not repeat the old mistake of listing the value as unrecognized
  expect_false(grepl("not correctly specified", msg))
})


test_that("the length message covers the selectivity settings too", {
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets
  base <- list(input_list = il, fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
               fish_q_spec = rep("fix", n))

  for(arg in c("fish_fixed_sel_pars_spec", "fish_sel_devs_spec", "fishsel_pe_pars_spec")) {
    msg <- err_msg(do.call(Setup_Mod_Fishsel_and_Q,
                           c(base, stats::setNames(list("est_all"), arg))))
    expect_match(msg, sprintf("%s has 1 entry for %d fleets", arg, n),
                 label = sprintf("length message for %s", arg))
  }
})


test_that("a genuinely unrecognized setting is still reported as one", {
  # The length guard must not swallow the case it was added beside.
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets

  msg <- err_msg(Setup_Mod_Fishsel_and_Q(
    il, fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
    fish_q_spec = rep("fix", n), fish_fixed_sel_pars_spec = rep("not_a_setting", n)))

  expect_match(msg, "fish_fixed_sel_pars_spec not correctly specified")
  expect_match(msg, "est_shared_r_s")
})


test_that("a correctly specified per-fleet setting is accepted", {
  # A guard that rejected everything would pass every test above.
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets

  expect_no_error(Setup_Mod_Fishsel_and_Q(
    il, fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
    fish_q_spec = rep("fix", n), fish_fixed_sel_pars_spec = rep("est_all", n)))
})


test_that("a composition type names the whole form it has to be given in", {
  # The old message listed 'agg' among the valid settings, so a caller who passed
  # exactly that was then told their fleet was invalid. Either message now carries
  # the full form. The form is named generically because one parser serves both
  # the composition streams and the at-age ones, which have different vocabularies
  # and are each listed by the message that raises them.
  n <- sweep_dims$n_fish_fleets

  bad_value <- err_msg(sweep_input(fishidx = list(
    FishAgeComps_Type = paste0("nope_Year_1-terminal_Fleet_", seq_len(n)))))

  expect_match(bad_value, "Value_Year_x-y_Fleet_f")
  expect_match(bad_value, "agg_Year_1-terminal_Fleet_1", fixed = TRUE)
  # the offending value is quoted back rather than left for the caller to find
  expect_match(bad_value, "nope_Year_1-terminal_Fleet_1", fixed = TRUE)
})


test_that("a composition type given as a bare setting names the form too", {
  # This is the case the old pair of messages handled worst: 'agg' passes the
  # value check, because it is a valid setting, and is then rejected for a fleet
  # the caller never wrote.
  n <- sweep_dims$n_fish_fleets

  bare <- err_msg(sweep_input(fishidx = list(FishAgeComps_Type = rep("agg", n))))

  expect_match(bare, "Value_Year_x-y_Fleet_f")
})


test_that("no setup message still carries the old misspelling", {
  # 'specfied' appeared in five messages. It is the kind of thing a reader
  # searching the source for their error will not find.
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package source is not laid out beside the tests")

  hits <- unlist(lapply(list.files(r_dir, pattern = "\\.R$", full.names = TRUE), function(f) {
    grep("specfied", readLines(f, warn = FALSE), value = TRUE)
  }))

  expect_equal(hits, character(0))
})


test_that("a plot destination that resolves to nothing is refused", {
  # here::here() drops a NULL, which left grDevices::pdf() with a zero-length file
  # name and sent a 16-page PDF called 'NA' into whatever the working directory
  # happened to be. Running the test suite therefore wrote into the repository,
  # and the file was committed three separate times before anyone noticed.
  for(bad in list(NULL, NA, character(0), c("a", "b"))) {
    expect_error(
      plot_all_basic(data = list(), rep = list(), sd_rep = list(),
                     model_names = "x", out_path = bad),
      "out_path must be a single directory",
      label = sprintf("out_path = %s", paste(deparse(bad), collapse = "")))
  }
})


test_that("a starting value of the wrong shape is refused where it is given", {
  # Starting values arrive through ... and were substituted for the model's own
  # default without being measured against it. A value of the wrong shape is not
  # rejected by that: it is carried into the objective, read position by
  # position, and indexes past its own end somewhere else entirely. What comes
  # back is RTMB's "'*this' is not a valid 'advector'", which names nothing the
  # caller wrote.
  #
  # The default carries the shape the model expects, so it is what the supplied
  # value is now checked against.
  expect_error(sweep_input(biol = list(M_spec = "est_ln_M", ln_M = rep(log(0.2), 7))),
               "starting value for ln_M is length 7")
  expect_error(sweep_input(fishsel = list(ln_fish_q = rep(0, 99))),
               "starting value for ln_fish_q is length 99")

  # the message says what the model wanted, not merely that it was unhappy
  expect_error(sweep_input(fishsel = list(ln_fish_q = rep(0, 99))),
               "model expects 3 by 1 by 5")
})


test_that("a correctly shaped starting value is still substituted", {
  # A guard that rejected every starting value would satisfy the test above and
  # break the feature.
  shape <- dim(sweep_input()$par$ln_fish_q)
  il <- expect_no_error(sweep_input(fishsel = list(ln_fish_q = array(-0.5, dim = shape))))

  expect_equal(as.numeric(unique(as.vector(il$par$ln_fish_q))), -0.5)
})


test_that("a parameter and its map must agree on length at model build", {
  # The same failure reached from the other side: a map is paired with its
  # parameter by position, so a length disagreement is either a mis-shaped
  # starting value or a map built from the wrong dimensions. Caught at the model
  # build with the block named, rather than inside the AD pass.
  il <- sweep_input()
  broken <- il
  broken$par$ln_fish_q <- broken$par$ln_fish_q[1]

  expect_error(fit_model(broken$data, broken$par, broken$map, do_optim = FALSE, silent = TRUE),
               "parameters and their maps disagree on length")
  expect_error(fit_model(broken$data, broken$par, broken$map, do_optim = FALSE, silent = TRUE),
               "ln_fish_q")

  expect_no_error(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
})


test_that("a single starting value where many are wanted is told how to recycle", {
  # The common slip has an obvious repair, and it is the same repair for every
  # argument, so the hint belongs in the shared guard rather than in a check
  # written for one parameter. It appears only where recycling is actually the
  # fix: a wrong length that is not one has no such repair to suggest.
  msg <- err_msg(sweep_input(fishsel = list(ln_fish_q = 0)))
  expect_match(msg, "Recycle it with rep(ln_fish_q, 15)", fixed = TRUE)

  no_hint <- err_msg(sweep_input(biol = list(M_spec = "est_ln_M", ln_M = rep(log(0.2), 7))))
  expect_match(no_hint, "starting value for ln_M is length 7")
  expect_false(grepl("Recycle", no_hint))
})


test_that("the per-fleet and starting-value guards stay distinct", {
  # Two guards, two different kinds of argument: check_fleet_spec_length covers a
  # configuration string given once per fleet, use_starting_value covers the shape
  # of a parameter array. Collapsing them would be wrong in both directions, so
  # each keeps its own wording and this records which is which.
  il <- sweep_input(stop_after = "srvidx")
  n <- il$data$n_fish_fleets

  fleet_setting <- err_msg(Setup_Mod_Fishsel_and_Q(
    il, fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
    fish_q_spec = "est_all", fish_fixed_sel_pars_spec = rep("est_all", n)))
  expect_match(fleet_setting, "fish_q_spec has 1 entry for 5 fleets")

  starting_value <- err_msg(sweep_input(fishsel = list(ln_fish_q = 0)))
  expect_match(starting_value, "starting value for ln_fish_q")
})
