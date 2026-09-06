# The apical sex offset against the scale sex offset. A scale multiplies the finished curve so the first
# and last bins move with it; an apical offset is the height the limbs build to, so those bins stay put.

library(SPoRC)
library(testthat)

bins <- 0:40
zero_devs <- array(0, dim = c(1, 1, 6, 2, 1))
curve <- function(pars, apical = 1) {
  as.vector(SPoRC:::Get_Selex(
    Selex_Model = 4,
    TimeVary_Model = 0,
    pars = pars,
    ln_seldevs = zero_devs,
    Region = 1,
    Year = 1,
    Bin = bins,
    Sex = 1,
    apical = apical
  ))
}
pars <- c(8, 0, log(6), log(60), -3, -1) # peak at bin 8, both ends well below one

test_that("an apical of one is the ordinary curve", {
  expect_equal(curve(pars, apical = 1), curve(pars), tolerance = 1e-12)
})

test_that("the apical is the height the plateau sits at, and the ends do not move with it", {
  A <- 0.3
  base <- curve(pars)
  off <- curve(pars, apical = A)

  # the plateau sits at the apical rather than at one
  expect_equal(max(off), A, tolerance = 1e-8)
  expect_equal(max(base), 1, tolerance = 1e-8)

  # The two ends stay where their own parameters put them. The joiners are what
  # hold the curve between the limbs and the plateau, and they are around 1e-8
  # rather than exactly zero this far from the peak, so that is the level these
  # hold to.
  expect_equal(off[1], base[1], tolerance = 1e-6)
  expect_equal(off[length(bins)], base[length(bins)], tolerance = 1e-6)
  expect_equal(off[1], stats::plogis(pars[5]), tolerance = 1e-6)
  expect_equal(off[length(bins)], stats::plogis(pars[6]), tolerance = 1e-6)

  # and the middle of the curve is lowered
  expect_lt(off[which.max(base)], base[which.max(base)])
})

test_that("a scale offset moves the ends and an apical offset does not", {
  A <- 0.3
  base <- curve(pars)
  scaled <- base * A          # what a scale offset does
  apic <- curve(pars, apical = A)

  # both put the plateau at the apical
  expect_equal(max(scaled), A, tolerance = 1e-8)
  expect_equal(max(apic), A, tolerance = 1e-8)

  # they part company at the two ends
  expect_equal(scaled[1], base[1] * A, tolerance = 1e-8)
  expect_equal(apic[1], base[1], tolerance = 1e-6)
  expect_gt(abs(apic[length(bins)] - scaled[length(bins)]), 1e-3)
})

test_that("setup accepts the apical levels and refuses them off the double normal", {
  build <- function(sel_model, offset) {
    il <- Setup_Mod_Dim(
      years = 1:5,
      ages = 1:10,
      lens = NA,
      n_regions = 1,
      n_sexes = 2,
      n_fish_fleets = 1,
      n_srv_fleets = 1,
      n_seas = 1,
      n_pop = 1,
      natal_region = 1,
      verbose = FALSE
    )
    suppressMessages(Setup_Mod_Fishsel_and_Q(
      input_list = il,
      fish_sel_model = paste0(sel_model, "_Fleet_1"),
      cont_tv_fish_sel = "none_Fleet_1",
      fish_sel_blocks = "none_Fleet_1",
      fish_q_blocks = "none_Fleet_1",
      fish_fixed_sel_pars_spec = "est_all",
      fish_sel_sex_offset = offset,
      fish_q_spec = "fix"))
  }
  ok <- build("dbnrml", "apical")
  expect_equal(as.numeric(ok$data$fishsel_sex_apical_offset), 1)
  expect_equal(as.numeric(ok$data$fishsel_sex_scale_offset), 0)

  both <- build("dbnrml", "par_apical")
  expect_equal(as.numeric(both$data$fishsel_sex_apical_offset), 1)
  expect_equal(as.numeric(both$data$fishsel_sex_par_offset), 1)

  # the offset the apical parameter rides on is created, as it is for a scale
  expect_true("ln_fishsel_sex_scale" %in% names(ok$par))
  expect_equal(sum(!is.na(as.numeric(as.character(ok$map$ln_fishsel_sex_scale)))), 1)

  # and it is refused where it would have no meaning
  expect_error(build("logist1", "apical"), "double normal")
  expect_error(build("dbnrml", "nonsense"), "must be one of")
})
