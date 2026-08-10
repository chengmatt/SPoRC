library(SPoRC)
library(testthat)

# ln_F_mean_spec = "fix" is the free log-F parameterization: the mean is mapped
# off at its starting value (defaulting to 0) so the deviations carry all of
# log F. Setup also warns when own-mean deviation centering is combined with an
# estimated mean in configurations where nothing else reads it, since the two
# then trade off along an exactly flat ridge.

mk_base <- function(init_F_form = "prop", init_age_strc = 1) {
  input_list <- Setup_Mod_Dim(years = 1:5, ages = 1:3, lens = NULL,
                              n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                              n_pop = 1, natal_region = 1, verbose = FALSE)
  suppressWarnings(Setup_Mod_Rec(input_list = input_list, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                                 ln_sigmaR = array(log(0.4), c(2, 1, 1)), rec_model = "mean_rec",
                                 sigmaR_spec = "fix", init_age_strc = init_age_strc,
                                 equil_init_age_strc = if(init_age_strc == 4) 0 else 2,
                                 init_F_form = init_F_form, ln_global_R0 = log(5)))
}

mk_catch <- function(input_list, Use_F_pen = 1, ...) {
  ObsCatch <- array(5, dim = c(1, 5, 1, 1))
  UseCatch <- array(1, dim = c(1, 5, 1, 1))
  suppressWarnings(Setup_Mod_Catch_and_F(input_list = input_list, ObsCatch = ObsCatch,
                                         UseCatch = UseCatch, Use_F_pen = Use_F_pen, ...))
}

test_that("the default keeps ln_F_mean estimated with its usual starting value", {
  out <- mk_catch(mk_base())
  expect_false(any(is.na(out$map$ln_F_mean)))
  expect_equal(as.numeric(out$par$ln_F_mean), log(0.1))
})

test_that("fix maps the mean off at zero so the deviations are log F outright", {
  out <- mk_catch(mk_base(), ln_F_mean_spec = "fix", Fdev_pen_center = "own_mean")
  expect_true(all(is.na(out$map$ln_F_mean)))
  expect_equal(as.numeric(out$par$ln_F_mean), 0)
})

test_that("a supplied starting value is respected under fix", {
  out <- mk_catch(mk_base(), ln_F_mean_spec = "fix", Fdev_pen_center = "own_mean",
                  ln_F_mean = array(log(0.2), dim = c(1, 1, 1)))
  expect_true(all(is.na(out$map$ln_F_mean)))
  expect_equal(as.numeric(out$par$ln_F_mean), log(0.2))
})

test_that("an invalid spec is rejected", {
  expect_error(mk_catch(mk_base(), ln_F_mean_spec = "nah"), "ln_F_mean_spec")
})

test_that("a fixed zero mean rejects penalties centred on it and accepts the rest", {

  # iid or ar1 centred on the fixed zero mean shrinks the deviations toward F = 1
  expect_error(mk_catch(mk_base(), ln_F_mean_spec = "fix"), "F = 1")
  expect_error(mk_catch(mk_base(), ln_F_mean_spec = "fix", Fdev_model = "ar1"), "F = 1")

  # own-mean centering, a random walk, or no penalty leave the level free
  expect_s3_class(mk_catch(mk_base(), ln_F_mean_spec = "fix", Fdev_pen_center = "own_mean")$map$ln_F_mean, "factor")
  expect_s3_class(mk_catch(mk_base(), ln_F_mean_spec = "fix", Fdev_model = "rw")$map$ln_F_mean, "factor")
  expect_s3_class(mk_catch(mk_base(), ln_F_mean_spec = "fix", Use_F_pen = 0)$map$ln_F_mean, "factor")

})

test_that("own-mean centering with an estimated mean warns only when nothing reads the mean", {

  ridge_warning <- function(input_list, ...) {
    ObsCatch <- array(5, dim = c(1, 5, 1, 1))
    UseCatch <- array(1, dim = c(1, 5, 1, 1))
    w <- character(0)
    withCallingHandlers(
      Setup_Mod_Catch_and_F(input_list = input_list, ObsCatch = ObsCatch,
                            UseCatch = UseCatch, Use_F_pen = 1, ...),
      warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
    )
    any(grepl("mutually unidentified", w))
  }

  # absolute-rate initialization F: nothing pins the level
  expect_true(ridge_warning(mk_base(init_F_form = "abs"), Fdev_pen_center = "own_mean"))
  # free initial age structure: likewise
  expect_true(ridge_warning(mk_base(init_age_strc = 4), Fdev_pen_center = "own_mean"))
  # proportional initialization F reads exp(ln_F_mean), so the mean is identified
  expect_false(ridge_warning(mk_base(init_F_form = "prop"), Fdev_pen_center = "own_mean"))
  # fixing the mean removes the ridge
  expect_false(ridge_warning(mk_base(init_F_form = "abs"), Fdev_pen_center = "own_mean",
                             ln_F_mean_spec = "fix"))
  # fixed centering never had one
  expect_false(ridge_warning(mk_base(init_F_form = "abs")))

})
