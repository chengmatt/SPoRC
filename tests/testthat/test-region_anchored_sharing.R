# The sharing specs anchored on region one and copied fleet maps region blind,
# which silently mapped selectivity off for any fleet whose data sit in a later
# region. These lock the fixed behavior: sharing anchors on a fleet's first
# region with data, and a fleet-sharing copy lands the reference's parameters in
# the sharing fleet's own regions.

# a fleet's catch confined to one region: fleet 1 in region 1, fleet 2 in region 2
region_split_catch <- function() {
  use <- array(1, dim = c(3, 13, 1, 5))
  use[2:3, , , 1] <- 0
  use[c(1, 3), , , 2] <- 0
  use
}

test_that("est_shared_f lands the reference's parameters in the sharing fleet's own region", {
  il <- sweep_input(
    catch = list(UseCatch = region_split_catch()),
    fishsel = list(fish_fixed_sel_pars_spec = c("est_all", "est_shared_f_1",
                                                "est_all", "est_all", "est_all"))
  )
  m <- array(as.integer(il$map$fish_fixed_sel_pars), dim = dim(il$par$fish_fixed_sel_pars))

  # the reference holds its parameters in region 1, its only region with data
  expect_true(all(!is.na(m[1, 1:2, 1, , 1])))
  # the sharing fleet holds the same parameters in region 2, its own region
  expect_identical(m[2, 1:2, 1, , 2], m[1, 1:2, 1, , 1])
})

test_that("est_shared_r anchors on the fleet's first region with data", {
  il <- sweep_input(
    catch = list(UseCatch = region_split_catch()),
    fishsel = list(fish_fixed_sel_pars_spec = c("est_all", "est_shared_r",
                                                "est_all", "est_all", "est_all"))
  )
  m <- array(as.integer(il$map$fish_fixed_sel_pars), dim = dim(il$par$fish_fixed_sel_pars))

  # fleet 2's data start in region 2, and its parameters exist rather than
  # being mapped off by an anchor it never reaches
  expect_true(all(!is.na(m[2, 1:2, 1, , 2])))
})

test_that("sigmaR specs resolve and share the region margin", {
  fake <- list(data = list(n_pop = 1, n_regions = 3, rec_dd = 1,
                           rec_region_prop_spec = 0, natal_region = 1),
               par = list(ln_sigmaR = array(0, dim = c(2, 1, 3))), map = list())

  n_est <- function(spec) nlevels(SPoRC:::do_sigmaR_mapping(fake, spec)$map$ln_sigmaR)
  expect_equal(n_est("est_all"), 6)            # every period and region
  expect_equal(n_est("est_shared_r"), 2)       # one per period, regions shared
  expect_equal(n_est("est_shared_all"), 1)
  expect_equal(n_est("fix_early_est_late"), 3) # late period, per region
  expect_equal(n_est("fix"), 0)
})

test_that("sigmaR slots without recruits stay off under natal homing", {
  fake <- list(data = list(n_pop = 2, n_regions = 3, rec_dd = 1,
                           rec_region_prop_spec = 1, natal_region = c(1, 2)),
               par = list(ln_sigmaR = array(0, dim = c(2, 2, 3))), map = list())

  m <- array(as.integer(SPoRC:::do_sigmaR_mapping(fake, "est_all")$map$ln_sigmaR),
             dim = c(2, 2, 3))
  # each population keeps its natal slot and nothing else
  expect_true(all(!is.na(m[, 1, 1])) && all(is.na(m[, 1, 2:3])))
  expect_true(all(!is.na(m[, 2, 2])) && all(is.na(m[, 2, c(1, 3)])))
})
