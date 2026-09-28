# What the setup functions say when an argument is wrong, since an error message is part of the interface.
#
# Two were wrong in a way that sent the reader somewhere unhelpful: a per-fleet setting given as a scalar was
# reported as invalid, and a composition type given in the form its own error listed was rejected anyway.
#
# Written against the exported functions rather than the internal validators, so they check the message
# a caller actually sees.

err_msg <- function(expr) {
  tryCatch({
    force(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
}

fleet_spec_input <- function() sweep_input(stop_after = "srvidx")


test_that("a per-fleet setting given once says so, rather than blaming its value", {
  # sweep_dims has five fishery fleets, so a scalar is short by four.
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets
  base <- list(input_list = il, fish_sel_model = paste0("logist1_Fleet_", seq_len(n)))

  msg <- err_msg(do.call(Setup_Mod_Fishsel_and_Q, c(base, list(
    fish_q_spec = "est_all", fish_fixed_sel_pars_spec = rep("est_all", n)))))

  expect_match(msg, "fish_q_spec has 1 entry for 5 fleets")
  # the message holds the call that fixes it, not just the diagnosis
  expect_match(msg, 'rep\\("est_all", 5\\)', fixed = FALSE)
  # and it must not repeat the old mistake of listing the value as unrecognized
  expect_false(grepl("not correctly specified", msg))
})


test_that("the length message covers the selectivity settings too", {
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets
  base <- list(
    input_list = il,
    fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
    fish_q_spec = rep("fix", n)
  )

  for(arg in c("fish_fixed_sel_pars_spec", "fish_sel_devs_spec", "fishsel_pe_pars_spec")) {
    msg <- err_msg(do.call(Setup_Mod_Fishsel_and_Q,
                           c(base, stats::setNames(list("est_all"), arg))))
    expect_match(msg, sprintf("%s has 1 entry for %d fleets", arg, n),
                 label = sprintf("length message for %s", arg))
  }
})


test_that("a genuinely unrecognized setting is still reported as one", {
  # the length check must not swallow the case it was added beside
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets

  msg <- err_msg(Setup_Mod_Fishsel_and_Q(
    il,
    fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
    fish_q_spec = rep("fix", n),
    fish_fixed_sel_pars_spec = rep("not_a_setting", n)
  ))

  expect_match(msg, "fish_fixed_sel_pars_spec not correctly specified")
  expect_match(msg, "est_shared_r_s")
})


test_that("a correctly specified per-fleet setting is accepted", {
  # a check that refused everything would pass every test above
  il <- fleet_spec_input()
  n <- il$data$n_fish_fleets

  expect_no_error(Setup_Mod_Fishsel_and_Q(
    il,
    fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
    fish_q_spec = rep("fix", n),
    fish_fixed_sel_pars_spec = rep("est_all", n)
  ))
})


test_that("a composition type names the whole form it has to be given in", {
  # both messages name the full form, since one listed 'agg' among the valid settings and
  # then told a caller who passed it that their fleet was invalid
  n <- sweep_dims$n_fish_fleets

  bad_value <- err_msg(sweep_input(fishidx = list(
    FishAgeComps_Type = paste0("nope_Year_1-terminal_Fleet_", seq_len(n)))))

  expect_match(bad_value, "Value_Year_x-y_Fleet_f")
  expect_match(bad_value, "agg_Year_1-terminal_Fleet_1", fixed = TRUE)
  # the offending value is quoted back rather than left for the caller to find
  expect_match(bad_value, "nope_Year_1-terminal_Fleet_1", fixed = TRUE)
})


test_that("a composition type given as a bare setting names the form too", {
  # 'agg' passes the value check, being a valid setting, so without this it was then
  # refused for a fleet the caller never wrote
  n <- sweep_dims$n_fish_fleets

  bare <- err_msg(sweep_input(fishidx = list(FishAgeComps_Type = rep("agg", n))))

  expect_match(bare, "Value_Year_x-y_Fleet_f")
})


test_that("no setup message still holds the old misspelling", {
  # 'specfied' appeared in five messages. It is the kind of thing a reader
  # searching the source for their error will not find.
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package source is not laid out beside the tests")

  # an installed package keeps a lazy-load database in R/ rather than sources, so this
  # skips rather than passing, since a scan of no files finds no misspellings either
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  skip_if(length(r_files) == 0, "package sources are not available to scan")

  # unlist() over a list of empty character vectors returns NULL rather than
  # character(0), so the comparison is made on a character vector either way
  hits <- as.character(unlist(lapply(r_files, function(f) {
    grep("specfied", readLines(f, warn = FALSE), value = TRUE)
  })))

  expect_equal(hits, character(0))
})


test_that("a plot destination that resolves to nothing is refused", {
  # here::here() drops a NULL, which leaves grDevices::pdf() with no file name and writes
  # a PDF called 'NA' into whatever the working directory happens to be
  for(bad in list(NULL, NA, character(0), c("a", "b"))) {
    expect_error(
      plot_all_basic(
        data = list(),
        rep = list(),
        sd_rep = list(),
        model_names = "x",
        out_path = bad
      ),
      "out_path must be a single directory",
      label = sprintf("out_path = %s", paste(deparse(bad), collapse = "")))
  }
})


test_that("a starting value of the wrong shape is refused where it is given", {
  # a starting value of the wrong shape otherwise reaches the objective, runs off its own end
  # and comes back as an RTMB advector complaint, which names nothing the caller wrote.
  #
  # the default holds the shape the model expects, so it is what the value is checked
  # against
  expect_error(sweep_input(biol = list(M_spec = "est_ln_M", ln_M = rep(log(0.2), 7))),
               "starting value for ln_M is length 7")
  expect_error(sweep_input(fishsel = list(ln_fish_q = rep(0, 99))),
               "starting value for ln_fish_q is length 99")

  # the message says what the model wanted, not merely that it was unhappy
  expect_error(sweep_input(fishsel = list(ln_fish_q = rep(0, 99))),
               "model expects 3 by 1 by 5")
})


test_that("a correctly shaped starting value is still substituted", {
  # a check that refused every starting value would satisfy the test above
  shape <- dim(sweep_input()$par$ln_fish_q)
  il <- expect_no_error(sweep_input(fishsel = list(ln_fish_q = array(-0.5, dim = shape))))

  expect_equal(as.numeric(unique(as.vector(il$par$ln_fish_q))), -0.5)
})


test_that("a parameter and its map must agree on length at model build", {
  # a map is paired with its parameter by position, so a length disagreement is either a
  # mis-shaped starting value or a map built from the wrong dimensions
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
  # the hint is the same for every argument, so it sits in the shared check, and appears
  # only where repeating the value is actually the fix
  msg <- err_msg(sweep_input(fishsel = list(ln_fish_q = 0)))
  expect_match(msg, "Recycle it with rep(ln_fish_q, 15)", fixed = TRUE)

  no_hint <- err_msg(sweep_input(biol = list(M_spec = "est_ln_M", ln_M = rep(log(0.2), 7))))
  expect_match(no_hint, "starting value for ln_M is length 7")
  expect_false(grepl("Recycle", no_hint))
})


test_that("the per-fleet and starting-value guards stay distinct", {
  # check_fleet_spec_length covers a setting given once per fleet and use_starting_value
  # the shape of a parameter array, so each keeps its own wording
  il <- sweep_input(stop_after = "srvidx")
  n <- il$data$n_fish_fleets

  fleet_setting <- err_msg(Setup_Mod_Fishsel_and_Q(
    il,
    fish_sel_model = paste0("logist1_Fleet_", seq_len(n)),
    fish_q_spec = "est_all",
    fish_fixed_sel_pars_spec = rep("est_all", n)
  ))
  expect_match(fleet_setting, "fish_q_spec has 1 entry for 5 fleets")

  starting_value <- err_msg(sweep_input(fishsel = list(ln_fish_q = 0)))
  expect_match(starting_value, "starting value for ln_fish_q")
})
