# Checks dsem_series: the names it writes are the names the setup offers, a filter picks by array index,
# a dim left out takes every level, and the estimated switch drops series the map holds fixed.

data("sgl_rg_dusky_data")

test_that("recruitment series come out named as the setup names them", {

  three <- suppressMessages(sweep_input())
  expect_equal(dsem_series(three, "rec"), paste0("rec_Pop_1_Region_", 1:3))
  expect_equal(dsem_series(three, "rec", region = 2), "rec_Pop_1_Region_2")
  expect_equal(dsem_series(three, "rec", Region = c(1, 3), pop = 1), paste0("rec_Pop_1_Region_", c(1, 3))) # case does not matter

  # one region, one population: the bare label, and no movement series at all
  dusky <- suppressMessages(build_goa_dusky_input(sgl_rg_dusky_data))
  expect_equal(dsem_series(dusky, "rec"), "rec")
  expect_equal(dsem_series(dusky, "move"), character(0))

  # the names are the ones the arrows accept
  n_yrs <- length(three$data$years)
  s <- dsem_series(three, "rec", region = 2)
  d <- suppressMessages(Setup_Mod_DSEM(three, sprintf("%s <-> %s, 0, NA, 0.5", s, s), dsem_data = NULL, dsem_processes = "rec"))
  expect_equal(d$data$dsem_var_names, s)

})

test_that("movement series follow the map: every free cell, filtered by index", {

  il <- suppressMessages(sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "iid_y_a_s", Movement_cont_pe_pars_spec = "fix")))
  dims <- dim(il$par$move_devs) # [pop, from, to, year, seas, age, sex]

  # what the map says, written out by hand
  m <- array(as.integer(il$map$move_devs), dim = dims)
  est <- which(apply(!is.na(m), c(1, 2, 3, 5, 6, 7), any), arr.ind = TRUE)
  by_hand <- sprintf("move_Pop_%d_From_%d_To_%d_Seas_%d_Age_%d_Sex_%d", est[,1], est[,2], est[,3], est[,4], est[,5], est[,6])
  expect_equal(dsem_series(il, "move"), by_hand)
  expect_equal(length(by_hand), dims[2] * dims[3] * (dims[6] - 1) * dims[7]) # age one is out, recruits do not move

  expect_equal(dsem_series(il, "move", from = 2, to = 1, age = 3, sex = 2), "move_Pop_1_From_2_To_1_Seas_1_Age_3_Sex_2")
  expect_equal(dsem_series(il, "move", age = 2:3), grep("_Age_[23]_", by_hand, value = TRUE))
  expect_equal(dsem_series(il, "move", age = 1), character(0)) # fixed by the map
  expect_length(dsem_series(il, "move", age = 1, estimated = FALSE), dims[2] * dims[3] * dims[7]) # held, but there

  # a fixed movement model estimates no deviation, so nothing is offered as linkable
  fixed <- suppressMessages(sweep_input())
  expect_equal(dsem_series(fixed, "move"), character(0))
  expect_gt(length(dsem_series(fixed, "move", estimated = FALSE)), 0)

})

test_that("numbers at age series are the state's ages, and bad filters are refused", {

  build_naa <- build_goa_dusky_input
  body(build_naa) <- do.call(substitute, list(body(build_naa), list(Setup_Mod_Biologicals = quote(function(...) Setup_Mod_Biologicals(..., NAA_re = "iid", NAA_re_ages = 5:9, NAA_sigma_spec = "fix")))))
  il <- suppressMessages(build_naa(sgl_rg_dusky_data))
  age_idx <- match(5:9, il$data$ages)
  expect_equal(dsem_series(il, "NAA"), sprintf("NAA_Pop_1_Region_1_Seas_1_Age_%d_Sex_1", age_idx))
  expect_equal(dsem_series(il, "NAA", age = age_idx[1]), "NAA_Pop_1_Region_1_Seas_1_Age_2_Sex_1") # by array index, not age label
  expect_equal(dsem_series(il, "NAA", age = 1), character(0)) # outside the state

  expect_error(dsem_series(il, "biomass"), "process should be one of")
  expect_error(dsem_series(il, "NAA", year = 3), "not a dim")
  expect_error(dsem_series(il, "NAA", 3), "Name every filter")
  expect_error(dsem_series(il, "NAA", age = 99), "from 1 to")
  expect_error(dsem_series(il, "NAA", age = 2, Age = 3), "filtered twice")

})
