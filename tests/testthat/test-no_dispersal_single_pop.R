# A single population under rec_region_prop_spec = "no_dispersal" stopped with "Comparison is generally
# unsafe for AD types", since the recruitment penalties tested its regional share against zero.
#
# Its share comes from parameters fixed at -20, so it is now the same model as that share fixed by hand.

test_that("one population no_dispersal runs and matches the share fixed by hand", {
  rec_args <- list(RecDevs_spec = "est_shared_r", InitDevs_spec = "est_shared_r")

  il_nd <- sweep_input(rec = c(rec_args, list(rec_region_prop_spec = "no_dispersal")))
  il_hand <- sweep_input(rec = rec_args)
  il_hand$par$rec_region_prop_pars[] <- -20 # what no_dispersal sets outside the natal region
  il_hand$map$rec_region_prop_pars <- factor(rep(NA, length(il_hand$par$rec_region_prop_pars)))

  expect_equal(il_nd$data$rec_region_prop_spec, 1)
  expect_equal(il_hand$data$rec_region_prop_spec, 0)

  obj_nd <- fit_model(il_nd$data, il_nd$par, il_nd$map, do_optim = FALSE, silent = TRUE)
  obj_hand <- fit_model(il_hand$data, il_hand$par, il_hand$map, do_optim = FALSE, silent = TRUE)
  rep_nd <- obj_nd$report(obj_nd$par)
  rep_hand <- obj_hand$report(obj_hand$par)

  expect_equal(length(obj_nd$par), length(obj_hand$par))
  expect_equal(obj_nd$fn(obj_nd$par), obj_hand$fn(obj_hand$par), tolerance = 1e-12)
  expect_equal(as.vector(obj_nd$gr(obj_nd$par)), as.vector(obj_hand$gr(obj_hand$par)), tolerance = 1e-12)

  # no region is skipped: the shared deviations split one penalty across every region
  expect_equal(sum(rep_nd$Rec_nLL), sum(rep_hand$Rec_nLL), tolerance = 1e-12)
  expect_equal(sum(rep_nd$Init_Rec_nLL), sum(rep_hand$Init_Rec_nLL), tolerance = 1e-12)
  expect_true(all(rep_nd$Rec_nLL[1, -1, ] != 0))
})
