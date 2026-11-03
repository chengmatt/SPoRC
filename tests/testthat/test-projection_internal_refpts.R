# The reference points solved on the tape against Get_Reference_Points on the packaged sablefish reports: every
# option (1e-5), the SPR options' derivatives against central differences of the post hoc solve (1e-3), the local
# MSY corner, and the internal projection on five regions against the post hoc one under independent SPR (1e-8).

library(SPoRC)
library(testthat)
data("sgl_rg_sable_data")
data("sgl_rg_sable_rep")
data("mlt_rg_sable_data")
data("mlt_rg_sable_rep")

what_codes <- c(SPR = 1, MSY = 2, independent_SPR = 3, independent_MSY = 4, global_SPR = 5, global_MSY = 6, local_MSY = 7)

# the arguments get_refpts_tape reads, from a packaged data list and report. the single-region report has no
# movement, straying or spawning movement, which the builder still reads, so those get their single-region values
tape_args <- function(d, r, code, natmort = NULL) {
  n_yrs <- length(d$years)
  n_pop <- d$n_pop
  n_regions <- d$n_regions
  n_ages <- length(d$ages)
  n_seas <- d$n_seas
  n_sexes <- d$n_sexes
  n_fish_fleets <- d$n_fish_fleets
  if(is.null(natmort)) natmort <- expand_natmort_seasons(r$natmort, n_seas) # the packaged reports predate the season dim
  Movement <- if(is.null(r$Movement)) array(1, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes)) else r$Movement
  sgl_move <- if(is.null(r$sgl_seas_spawning_movement)) array(1, dim = c(n_pop, n_regions, n_regions, n_yrs, n_ages, n_sexes)) else r$sgl_seas_spawning_movement
  stray_rate <- if(is.null(r$stray_rate)) array(0, dim = c(n_pop, n_yrs)) else r$stray_rate
  list(refpt_type = if(code <= 2) 1 else 2, refpt_what = code, SPR_x = 0.4, sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)),
       is_discard_fleet = array(0, dim = n_fish_fleets), local_bh_msy_newton_steps = 6, newton_steps = 2, log_F_start = log(0.1),
       avg_yrs = n_yrs, rec_yrs = 1:(n_yrs - 1),
       n_yrs = n_yrs, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages, n_seas = n_seas, n_fish_fleets = n_fish_fleets,
       seasdur = if(is.null(d$seasdur)) 1 else d$seasdur, spawn_seas = if(is.null(d$spawn_seas)) 1 else d$spawn_seas, t_spawn = 0,
       natal_region = if(is.null(d$natal_region)) 1 else d$natal_region, do_recruits_move = if(is.null(d$do_recruits_move)) 0 else d$do_recruits_move,
       move_timing = 0, expm_nsub = 0, rec_model = 1,
       Fmort = r$Fmort, dmr = r$dmr, fish_sel = r$fish_sel, ret_sel = r$ret_sel, natmort = natmort, WAA = d$WAA, MatAA = d$MatAA,
       Movement = Movement, Mrate = NULL, sgl_seas_spawning_movement = sgl_move, stray_rate = stray_rate,
       rec_seas_prop = r$rec_seas_prop, rec_region_prop = r$rec_region_prop, Rec = r$Rec, R0 = r$R0, h_trans = r$h_trans)
}

post_hoc <- function(d, r, what) {
  type <- if(what %in% c("SPR", "MSY")) "single_region" else "multi_region"
  suppressWarnings(Get_Reference_Points(data = d, rep = r, SPR_x = 0.4, type = type, what = what))
}

# the tape's derivative of the reference points in a natural mortality scale, against central differences of the post hoc solve
derivative_check <- function(d, r, what, h = 1e-4) {
  natmort6 <- expand_natmort_seasons(r$natmort, d$n_seas)
  f <- function(eps) {
    out <- do.call(SPoRC:::get_refpts_tape, tape_args(d, r, what_codes[[what]], natmort = natmort6 * exp(eps[1])))
    c(out$f_ref_pt, as.numeric(out$b_ref_pt))
  }
  tape <- as.numeric(RTMB::MakeTape(f, 0)$jacobian(0))
  rr <- r
  rr$natmort <- r$natmort * exp(h)
  up <- post_hoc(d, rr, what)
  rr$natmort <- r$natmort * exp(-h)
  dn <- post_hoc(d, rr, what)
  list(tape = tape, fd = (c(up$f_ref_pt, as.numeric(up$b_ref_pt)) - c(dn$f_ref_pt, as.numeric(dn$b_ref_pt))) / (2 * h))
}

test_that("every reference point option solved on the tape matches Get_Reference_Points", {

  for(what in names(what_codes)) {
    single <- what %in% c("SPR", "MSY")
    d <- if(single) sgl_rg_sable_data else mlt_rg_sable_data
    r <- if(single) sgl_rg_sable_rep else mlt_rg_sable_rep
    post <- post_hoc(d, r, what)
    tape <- do.call(SPoRC:::get_refpts_tape, tape_args(d, r, what_codes[[what]]))

    # a region at zero fishing under local MSY is left at the value solve's answer, so compare where F is real
    fished <- as.numeric(post$f_ref_pt) > 1e-3
    expect_equal(as.numeric(tape$f_ref_pt)[fished], as.numeric(post$f_ref_pt)[fished], tolerance = 1e-5, info = what)
    expect_equal(as.numeric(tape$b_ref_pt), as.numeric(post$b_ref_pt), tolerance = 1e-5, info = what)
    expect_equal(as.numeric(tape$virgin_b_ref_pt), as.numeric(post$virgin_b_ref_pt), tolerance = 1e-5, info = what)
    expect_equal(as.numeric(tape$pop_b_ref_pt), as.numeric(post$pop_b_ref_pt), tolerance = 1e-5, info = what)
    expect_equal(as.numeric(tape$virgin_pop_b_ref_pt), as.numeric(post$virgin_pop_b_ref_pt), tolerance = 1e-5, info = what)
  } # end what loop

})

test_that("the SPR reference points have the derivative the implicit function theorem gives", {

  chk <- derivative_check(sgl_rg_sable_data, sgl_rg_sable_rep, "SPR")
  expect_equal(chk$tape, chk$fd, tolerance = 1e-3)
  for(what in c("global_SPR", "independent_SPR")) {
    chk <- derivative_check(mlt_rg_sable_data, mlt_rg_sable_rep, what)
    expect_true(all(is.finite(chk$tape)), info = what)
    expect_equal(chk$tape, chk$fd, tolerance = 1e-3, info = what)
  } # end what loop

})

test_that("local MSY at a corner solution keeps a derivative where fishing happens and none where it does not", {

  d <- mlt_rg_sable_data
  r <- mlt_rg_sable_rep
  post <- post_hoc(d, r, "local_MSY")
  unfished <- as.numeric(post$f_ref_pt) < 1e-6
  expect_true(any(unfished)) # the packaged stock's local MSY closes regions, which is what this checks

  chk <- derivative_check(d, r, "local_MSY")
  f_rows <- 1:d$n_regions
  b_rows <- d$n_regions + 1:d$n_regions
  expect_true(all(chk$tape[f_rows][unfished] == 0))
  expect_equal(chk$tape[f_rows][!unfished], chk$fd[f_rows][!unfished], tolerance = 0.1)
  expect_equal(chk$tape[b_rows], chk$fd[b_rows], tolerance = 1e-2)

})

test_that("the internal projection on five regions matches the post hoc one under independent SPR and the threshold rule", {

  d <- mlt_rg_sable_data
  r <- mlt_rg_sable_rep
  n_yrs <- length(d$years)
  n_pop <- d$n_pop
  n_regions <- d$n_regions
  n_ages <- length(d$ages)
  n_seas <- d$n_seas
  n_sexes <- d$n_sexes
  n_fish_fleets <- d$n_fish_fleets
  n_proj_yrs <- 6
  natmort6 <- expand_natmort_seasons(r$natmort, n_seas)
  WAA_fish <- array(0, dim = c(dim(d$WAA), n_fish_fleets))
  for(f in 1:n_fish_fleets) WAA_fish[,,,,,,f] <- d$WAA
  sexratio <- if(is.null(r$sexratio)) array(0.5, dim = c(n_pop, n_regions, n_yrs, n_sexes)) else r$sexratio
  ln_RecDevs <- array(0, dim = c(n_pop, n_regions, n_yrs + n_proj_yrs - 1)) # only the projected years are read, and the packaged report has none, so they are zero

  internal <- SPoRC:::get_internal_projection(
    proj_n_yrs = n_proj_yrs, proj_f_opt = 4, proj_f_input = array(0, dim = c(n_regions, n_proj_yrs)),
    proj_f_ref_pt = array(0, dim = c(n_regions, n_proj_yrs)), proj_b_ref_pt = array(1, dim = c(n_pop, n_regions, n_proj_yrs)),
    proj_hcr_function = NULL, proj_hcr_scope = 1, proj_F_constraint = NULL, proj_catch_input = NULL,
    proj_refpt_type = 2, proj_refpt_what = 3, proj_spr_x = 0.4, proj_sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)),
    proj_is_discard_fleet = array(0, dim = n_fish_fleets), proj_local_bh_msy_newton_steps = 6, proj_refpt_newton_steps = 2,
    proj_refpt_log_f_start = log(0.1), proj_n_avg_yrs = 1, proj_sel_opt = 1, proj_move_opt = 1, proj_calc_rec_yrs = 1:(n_yrs - 1),
    proj_rec_opt = 1, proj_naa_opt = 1, proj_growth_opt = 1,
    n_yrs = n_yrs, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages, n_sexes = n_sexes, n_seas = n_seas, n_fish_fleets = n_fish_fleets,
    seasdur = if(is.null(d$seasdur)) 1 else d$seasdur, spawn_seas = if(is.null(d$spawn_seas)) 1 else d$spawn_seas,
    natal_region = d$natal_region, do_recruits_move = if(is.null(d$do_recruits_move)) 0 else d$do_recruits_move,
    move_timing = 0, expm_nsub = 0, t_spawn = 0, rec_model = 1,
    NAA = r$NAA, NAA0 = r$NAA0, Fmort = r$Fmort, Rec = r$Rec, natmort = natmort6, WAA = d$WAA, WAA_fish = WAA_fish, MatAA = d$MatAA,
    fish_sel = r$fish_sel, ret_sel = r$ret_sel, dmr = r$dmr, Movement = r$Movement, Mrate = NULL,
    sgl_seas_spawning_movement = r$sgl_seas_spawning_movement, sexratio = sexratio, stray_rate = r$stray_rate,
    rec_seas_prop = r$rec_seas_prop, rec_region_prop = r$rec_region_prop, R0 = r$R0, h_trans = r$h_trans, ln_RecDevs = ln_RecDevs,
    ln_NAA = NULL, n_est_naa_re = 0, naa_re_ages = NULL, naa_re_seas = NULL, naa_re_where = NULL) # no state in the packaged report

  post <- post_hoc(d, r, "independent_SPR")
  r6 <- r
  r6$natmort <- natmort6
  inputs <- build_projection_inputs(data = d, rep = r6, refpts = post, n_proj_yrs = n_proj_yrs, t_spawn = 0)
  external <- do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "HCR", HCR_function = hcr_threshold)))

  expect_equal(as.numeric(internal$proj_SSB), as.numeric(external$proj_SSB), tolerance = 1e-8)
  expect_equal(as.numeric(internal$proj_Catch), as.numeric(external$proj_Catch), tolerance = 1e-8)
  expect_equal(as.numeric(internal$refpts$f_ref_pt), as.numeric(post$f_ref_pt), tolerance = 1e-6)
  expect_equal(dim(internal$proj_status), c(n_regions, n_proj_yrs))
  expect_equal(as.numeric(internal$proj_status[,1]), as.numeric(external$proj_SSB[1,,1]) / as.numeric(post$b_ref_pt[1,]), tolerance = 1e-8)

})
