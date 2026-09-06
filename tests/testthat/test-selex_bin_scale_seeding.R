# Starting values stated on a bin scale have to be read on the bin scale the model is using.
#
# A double normal seeds its peak at the middle of the bin vector, which differs between ages and lengths:
# on ages 1-6 and lengths 20-100, a peak read off the wrong vector lands at 3.5, below every length bin.

selex_seed_input <- function(selex_type = "length", n_lens = 9) {
  NY <- 8; NAG <- 6
  lens <- seq(20, 100, length.out = n_lens)
  # a size-age key spreading each age across the length bins, so length
  # compositions are something the model can fit
  sizeage <- vapply(seq_len(NAG), function(a) {
    w <- stats::dnorm(seq_len(n_lens), 1 + (a - 1) * (n_lens - 1) / (NAG - 1), 1)
    w / sum(w)
  }, numeric(n_lens))

  il <- Setup_Mod_Dim(
    n_pop = 1,
    years = seq_len(NY),
    ages = seq_len(NAG),
    lens = lens,
    n_regions = 1,
    n_sexes = 1,
    n_seas = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    verbose = FALSE
  )
  il <- Setup_Mod_Rec(
    il,
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    do_rec_bias_ramp = 0,
    init_age_strc = 1,
    ln_global_R0 = log(1e6)
  )

  bd <- c(1, 1, NY, 1, NAG, 1)
  sa <- array(0, dim = c(1, 1, NY, 1, n_lens, NAG, 1))
  for(y in seq_len(NY)) sa[1, 1, y, 1, , , 1] <- sizeage

  il <- suppressWarnings(Setup_Mod_Biologicals(
    il,
    WAA = array(1, dim = bd),
    WAA_fish = array(1, dim = c(bd, 1)),
    WAA_srv = array(1, dim = c(bd, 1)),
    MatAA = array(1, dim = bd),
    fit_lengths = 1,
    SizeAgeTrans = sa,
    M_spec = "fix",
    Fixed_natmort = array(0.2, dim = c(1, 1, NY, NAG, 1))
  ))
  il <- Setup_Mod_Movement(il, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  il <- Setup_Mod_Tagging(il, use_conv_fish_tagging = 0)
  il <- suppressWarnings(Setup_Mod_Catch_and_F(
    il,
    ObsCatch = array(1e4, dim = c(1, NY, 1, 1)),
    UseCatch = array(1, dim = c(1, NY, 1, 1)),
    sigmaC_spec = "fix",
    sigmaF_spec = "fix"
  ))

  off <- array(0, dim = c(1, NY, 1, 1))
  il <- do.call(Setup_Mod_FishIdx_and_Comps, list(
    input_list = il,
    ObsFishIdx = array(1e5, dim = c(1, NY, 1, 1)),
    ObsFishIdx_SE = array(0.2, dim = c(1, NY, 1, 1)),
    UseFishIdx = off,
    fish_idx_type = "none",
    ObsFishAgeComps = array(1 / NAG, dim = c(1, NY, 1, NAG, 1, 1)),
    UseFishAgeComps = off,
    ISS_FishAgeComps = array(100, dim = c(1, NY, 1, 1, 1)),
    FishAgeComps_LikeType = "none",
    ObsFishLenComps = array(1 / n_lens, dim = c(1, NY, 1, n_lens, 1, 1)),
    UseFishLenComps = array(1, dim = c(1, NY, 1, 1)),
    ISS_FishLenComps = array(100, dim = c(1, NY, 1, 1, 1)),
    FishLenComps_LikeType = "Multinomial",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "agg_Year_1-terminal_Fleet_1"
  ))

  il <- suppressWarnings(Setup_Mod_SrvIdx_and_Comps(
    il,
    ObsSrvIdx = array(1e5, dim = c(1, NY, 1, 1)),
    ObsSrvIdx_SE = array(0.2, dim = c(1, NY, 1, 1)),
    UseSrvIdx = off,
    srv_idx_type = "none",
    ObsSrvAgeComps = array(1 / NAG, dim = c(1, NY, 1, NAG, 1, 1)),
    UseSrvAgeComps = off,
    ISS_SrvAgeComps = array(100, dim = c(1, NY, 1, 1, 1)),
    SrvAgeComps_LikeType = "none",
    ObsSrvLenComps = array(0, dim = c(1, NY, 1, n_lens, 1, 1)),
    UseSrvLenComps = off,
    ISS_SrvLenComps = array(0, dim = c(1, NY, 1, 1, 1)),
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "agg_Year_1-terminal_Fleet_1"
  ))

  list(il = il, lens = lens, ages = seq_len(NAG))
}


test_that("a length-based double normal seeds its peak among the length bins", {
  # Currently FAILS. fish_selex_type is converted from "age"/"length" to 0/1
  # before the seeding call reads it, so `fish_selex_type == 'length'` there is
  # comparing a number against a string and is never true. The peak is seeded at
  # the middle of the age range whatever the selectivity is over.
  fx <- selex_seed_input()

  z <- Setup_Mod_Fishsel_and_Q(
    fx$il,
    fish_sel_model = "dbnrml_Fleet_1",
    fish_selex_type = "length",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix"
  )

  peak <- z$par$fish_fixed_sel_pars[1, 1, 1, 1, 1]
  expect_gte(peak, min(fx$lens))
  expect_lte(peak, max(fx$lens))
  expect_equal(peak, min(fx$lens) + 0.5 * diff(range(fx$lens)))
})


test_that("a length-based survey double normal seeds its peak among the length bins", {
  # Currently FAILS, for the same reason in srv_selex_type.
  fx <- selex_seed_input()
  il <- Setup_Mod_Fishsel_and_Q(
    fx$il,
    fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix"
  )

  z <- Setup_Mod_Srvsel_and_Q(
    il,
    srv_sel_model = "dbnrml_Fleet_1",
    srv_selex_type = "length",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "fix"
  )

  peak <- z$par$srv_fixed_sel_pars[1, 1, 1, 1, 1]
  expect_gte(peak, min(fx$lens))
  expect_lte(peak, max(fx$lens))
})


test_that("an age-based double normal still seeds its peak among the ages", {
  # The age path is what the seeding has always done, and it stays correct.
  fx <- selex_seed_input()

  z <- Setup_Mod_Fishsel_and_Q(
    fx$il,
    fish_sel_model = "dbnrml_Fleet_1",
    fish_selex_type = "age",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix"
  )

  peak <- z$par$fish_fixed_sel_pars[1, 1, 1, 1, 1]
  expect_equal(peak, min(fx$ages) + 0.5 * diff(range(fx$ages)))
})


test_that("an unrecognized selectivity domain is rejected rather than left undefined", {
  # Currently FAILS. Neither the 'age' nor the 'length' branch assigns `bins` for
  # any other value and there is no else, so the setup dies later on
  # "object 'bins' not found" instead of naming the argument at fault.
  fx <- selex_seed_input()

  expect_error(
    Setup_Mod_Fishsel_and_Q(
      fx$il,
      fish_sel_model = "logist1_Fleet_1",
      fish_selex_type = "lenght",
      fish_fixed_sel_pars_spec = "est_all",
      fish_q_spec = "fix"
    ),
    "fish_selex_type")
})
