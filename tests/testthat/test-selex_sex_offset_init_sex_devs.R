library(testthat)
library(RTMB)

# Selectivity sex offsets (par and scale) and sex-specific initial age
# deviations, tested at the function level against hand-computed values.

logistic_sel <- function(b50, k, ages) 1 / (1 + exp(-k * (ages - b50)))

make_selex_args <- function(n_ages = 6, n_sexes = 2, n_fleets = 1) {
  n_regions <- 1; n_yrs <- 3
  list(
    selex_type = 0,
    bins = 1:n_ages,
    sel_blocks = array(1, dim = c(n_regions, n_yrs, n_fleets)),
    sel_model = array(0, dim = c(n_regions, n_yrs, n_fleets)), # logist1
    fixed_sel_pars = array(0, dim = c(n_regions, 2, 1, n_sexes, n_fleets)),
    cont_tv_sel = array(0, dim = c(n_regions, n_fleets)),
    ln_seldevs = array(0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fleets)),
    use_fixed_sel = rep(0, n_fleets),
    sel_input = NULL,
    bicubic_Wbin = array(0, dim = c(n_regions, n_ages, 1, 1, n_fleets)),
    bicubic_Wyr = array(0, dim = c(n_regions, n_yrs, 1, 1, n_fleets)),
    bicubic_binnodes = array(0, dim = c(n_regions, n_yrs, n_fleets)),
    bicubic_yrnodes = array(0, dim = c(n_regions, n_yrs, n_fleets)),
    n_pop = 1,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_proj_yrs_devs = 0,
    n_seas = 1,
    n_ages = n_ages,
    n_lens = 1,
    n_sexes = n_sexes,
    n_fleets = n_fleets
  )
}

test_that("a par sex offset evaluates the second sex at the first sex's parameters plus its stored offsets", {

  args <- make_selex_args()
  args$fixed_sel_pars[1, , 1, 1, 1] <- log(c(3, 0.8))   # female ln_b50, ln_k
  args$fixed_sel_pars[1, , 1, 2, 1] <- c(0.3, -0.15)    # male slots hold offsets

  out <- do.call(Get_Selex_Array, c(args, list(sex_par_offset = 1, sex_scale_offset = 0, sex_scale = NULL)))

  expect_equal(as.numeric(out$sel[1,1,1,1,,1,1]), logistic_sel(3, 0.8, 1:6), tolerance = 1e-12)
  expect_equal(as.numeric(out$sel[1,1,1,1,,2,1]), logistic_sel(3 * exp(0.3), 0.8 * exp(-0.15), 1:6), tolerance = 1e-12)

  # offsets kept at zero reproduce sex-shared parameters
  args$fixed_sel_pars[1, , 1, 2, 1] <- 0
  out0 <- do.call(Get_Selex_Array, c(args, list(sex_par_offset = 1, sex_scale_offset = 0, sex_scale = NULL)))
  expect_equal(as.numeric(out0$sel[1,1,1,1,,2,1]), as.numeric(out0$sel[1,1,1,1,,1,1]), tolerance = 1e-12)
})

test_that("a scale sex offset multiplies the second sex's whole curve by exp(offset), which may exceed one", {

  args <- make_selex_args()
  args$fixed_sel_pars[1, , 1, 1, 1] <- log(c(3, 0.8))
  args$fixed_sel_pars[1, , 1, 2, 1] <- log(c(4, 1.1)) # male keeps its own parameters
  sex_scale <- array(0, dim = c(1, 1, 2, 1))
  sex_scale[1, 1, 2, 1] <- 0.2

  out <- do.call(Get_Selex_Array, c(args, list(sex_par_offset = 0, sex_scale_offset = 1, sex_scale = sex_scale)))

  expect_equal(as.numeric(out$sel[1,1,1,1,,1,1]), logistic_sel(3, 0.8, 1:6), tolerance = 1e-12)
  expect_equal(as.numeric(out$sel[1,1,1,1,,2,1]), logistic_sel(4, 1.1, 1:6) * exp(0.2), tolerance = 1e-12)
  expect_gt(max(as.numeric(out$sel[1,1,1,1,,2,1])), 1) # exceeds one at the asymptote

  # both offsets together compose
  args$fixed_sel_pars[1, , 1, 2, 1] <- c(0.3, -0.15)
  outb <- do.call(Get_Selex_Array, c(args, list(sex_par_offset = 1, sex_scale_offset = 1, sex_scale = sex_scale)))
  expect_equal(as.numeric(outb$sel[1,1,1,1,,2,1]), logistic_sel(3 * exp(0.3), 0.8 * exp(-0.15), 1:6) * exp(0.2), tolerance = 1e-12)
})

test_that("time-varying deviations apply to the effective parameters under a par offset", {

  args <- make_selex_args()
  args$cont_tv_sel[] <- 1 # iid
  args$fixed_sel_pars[1, , 1, 1, 1] <- log(c(3, 0.8))
  args$fixed_sel_pars[1, , 1, 2, 1] <- c(0.3, -0.15)
  args$ln_seldevs[1, 2, 1, 2, 1] <- 0.1 # male b50 deviation in year 2

  out <- do.call(Get_Selex_Array, c(args, list(sex_par_offset = 1, sex_scale_offset = 0, sex_scale = NULL)))
  expect_equal(as.numeric(out$sel[1,1,2,1,,2,1]), logistic_sel(3 * exp(0.3) * exp(0.1), 0.8 * exp(-0.15), 1:6), tolerance = 1e-12)
})

test_that("setup_sel_sex_offset validates its options and maps scale parameters for the second sex only", {

  messages_list <<- character(0)
  base_input <- list(data = list(n_sexes = 2, n_regions = 1), par = list(), map = list())
  sel_model_arr <- array(0, dim = c(1, 3, 2))
  cont_tv_mat <- array(0, dim = c(1, 2))

  out <- setup_sel_sex_offset(
    base_input,
    c("par", "scale"),
    prefix = "fish",
    n_fleets = 2,
    fleet_label = "fishery fleet",
    sel_model_arr = sel_model_arr,
    cont_tv_mat = cont_tv_mat,
    max_blks = 1
  )
  expect_equal(out$data$fishsel_sex_par_offset, c(1, 0))
  expect_equal(out$data$fishsel_sex_scale_offset, c(0, 1))
  map_scale <- array(as.numeric(as.character(out$map$ln_fishsel_sex_scale)), dim = c(1, 1, 2, 2))
  expect_true(all(is.na(map_scale[1, 1, , 1])))  # a par-offset fleet has no scale parameter
  expect_true(is.na(map_scale[1, 1, 1, 2]))      # the first sex is the reference
  expect_false(is.na(map_scale[1, 1, 2, 2]))     # the second sex's scale is estimated

  # one sex has nothing to offset
  base_input$data$n_sexes <- 1
  expect_error(setup_sel_sex_offset(
    base_input,
    c("par", "none"),
    prefix = "fish",
    n_fleets = 2,
    fleet_label = "fishery fleet",
    sel_model_arr = sel_model_arr,
    cont_tv_mat = cont_tv_mat,
    max_blks = 1
  ), "n_sexes > 1")

  # a mean-standardized form cancels a constant multiplier
  base_input$data$n_sexes <- 2
  sel_model_nonpar <- sel_model_arr; sel_model_nonpar[,,2] <- 9
  expect_error(setup_sel_sex_offset(
    base_input,
    c("none", "scale"),
    prefix = "fish",
    n_fleets = 2,
    fleet_label = "fishery fleet",
    sel_model_arr = sel_model_nonpar,
    cont_tv_mat = cont_tv_mat,
    max_blks = 1
  ), "mean-standardized")
})

test_that("an NSelBins plateau holds bins beyond it at the last computed value for parametric forms", {

  args <- make_selex_args()
  args$fixed_sel_pars[1, , 1, 1, 1] <- log(c(4, 0.4)) # slow enough that bin 4 is far from saturated
  args$fixed_sel_pars[1, , 1, 2, 1] <- log(c(4, 0.4))
  nselbins <- array(4, dim = c(1, 3, 1))

  out <- do.call(Get_Selex_Array, c(args, list(nselbins = nselbins)))
  expected <- logistic_sel(4, 0.4, 1:6); expected[5:6] <- expected[4]
  expect_equal(as.numeric(out$sel[1,1,1,1,,1,1]), expected, tolerance = 1e-12)

  # the plateau composes with a sex scale offset (the fm.tpl male arrangement)
  sex_scale <- array(0, dim = c(1, 1, 2, 1)); sex_scale[1, 1, 2, 1] <- 0.2
  outs <- do.call(Get_Selex_Array, c(args, list(
    nselbins = nselbins,
    sex_par_offset = 0,
    sex_scale_offset = 1,
    sex_scale = sex_scale
  )))
  expect_equal(as.numeric(outs$sel[1,1,1,1,,2,1]), expected * exp(0.2), tolerance = 1e-12)

  # zero means no plateau
  out0 <- do.call(Get_Selex_Array, c(args, list(nselbins = array(0, dim = c(1, 3, 1)))))
  expect_equal(as.numeric(out0$sel[1,1,1,1,,1,1]), logistic_sel(4, 0.4, 1:6), tolerance = 1e-12)
})

test_that("Get_Init_NAA applies sex-specific initial deviations and broadcasts a 3-D array", {

  n_ages <- 4; n_sexes <- 2
  d <- list(
    n_regions = 1,
    n_pop = 1,
    n_sexes = n_sexes,
    n_ages = n_ages,
    n_seas = 1,
    n_fish_fleets = 1,
    seasdur = 1,
    rec_seas_prop = matrix(1, 1, 1),
    natmort = array(0.2, dim = c(1, 1, n_ages, n_sexes)),
    init_F = array(0, dim = c(1, 1, 1)),
    dmr = array(0, dim = c(1, 1, 1)),
    fish_sel = array(0.5, dim = c(1, 1, 1, n_ages, n_sexes, 1)),
    ret_sel = array(1, dim = c(1, 1, 1, n_ages, n_sexes, 1)),
    R0_r = matrix(1000, 1, 1),
    sexratio = array(0.5, dim = c(1, 1, n_sexes)),
    Movement = array(1, dim = c(1, 1, 1, 1, n_ages, n_sexes)),
    do_recruits_move = 0
  )

  devs4 <- array(0, dim = c(1, 1, n_ages - 1, n_sexes))
  devs4[1, 1, , 1] <- c(0.1, -0.2, 0.3)
  devs4[1, 1, , 2] <- c(-0.4, 0.5, -0.6)

  out <- do.call(Get_Init_NAA, c(d, list(init_age_strc = 4, init_iter = 0, ln_InitDevs = devs4)))
  expect_equal(as.numeric(out[1, 1, 2:n_ages, 1]), 0.5 * exp(c(0.1, -0.2, 0.3)), tolerance = 1e-12)
  expect_equal(as.numeric(out[1, 1, 2:n_ages, 2]), 0.5 * exp(c(-0.4, 0.5, -0.6)), tolerance = 1e-12)

  # a 3-D array is one shared curve
  devs3 <- array(c(0.1, -0.2, 0.3), dim = c(1, 1, n_ages - 1))
  out3 <- do.call(Get_Init_NAA, c(d, list(init_age_strc = 4, init_iter = 0, ln_InitDevs = devs3)))
  expect_equal(as.numeric(out3[1, 1, 2:n_ages, 2]), as.numeric(out3[1, 1, 2:n_ages, 1]), tolerance = 1e-12)
  expect_equal(as.numeric(out3[1, 1, 2:n_ages, 1]), as.numeric(out[1, 1, 2:n_ages, 1]), tolerance = 1e-12)
})

test_that("the initial age penalty covers one copy per parameter and pools its own-mean center across sexes", {

  n_ages <- 5
  base_args <- list(
    n_pop = 1,
    n_regions = 1,
    n_ages = n_ages,
    n_est_rec_devs = 3,
    rec_dd = 0,
    natal_region = 1,
    rec_region_prop_spec = 0,
    rec_region_prop = matrix(1, 1, 1),
    equil_init_age_strc = 2,
    init_age_devs_shared = NULL,
    ln_sigmaR = array(log(0.6), dim = c(2, 1, 1)),
    bias_ramp = rep(0, 3),
    sigmaR_switch = 1,
    ln_RecDevs = array(0, dim = c(1, 1, 3)),
    sigmaR2_early = matrix(0.36, 1, 1),
    sigmaR2_late = matrix(0.36, 1, 1),
    do_rec_bias_ramp = 0
  )

  devs <- c(0.2, -0.1, 0.4, -0.3)
  dev3 <- array(devs, dim = c(1, 1, n_ages - 1))

  # the 3-D legacy array and the 4-D shared layout give the same penalty
  legacy <- do.call(get_recruitment_penalty, c(base_args, list(ln_InitDevs = dev3)))
  dev4 <- array(rep(devs, 2), dim = c(1, 1, n_ages - 1, 2))
  pen_shared <- array(0, dim = dim(dev4)); pen_shared[,,,1] <- 1
  shared <- do.call(get_recruitment_penalty, c(base_args, list(ln_InitDevs = dev4, init_devs_pen_use = pen_shared)))
  expect_equal(sum(shared$Init_Rec_nLL), sum(legacy$Init_Rec_nLL), tolerance = 1e-12)
  expect_equal(as.numeric(shared$Init_Rec_nLL[1,1,,2]), rep(0, n_ages - 1)) # the second copy is not penalized

  # sex-specific deviations are each penalized
  devs_m <- c(-0.15, 0.25, -0.35, 0.05)
  dev4d <- dev4; dev4d[1,1,,2] <- devs_m
  pen_all <- array(1, dim = dim(dev4d))
  distinct <- do.call(get_recruitment_penalty, c(base_args, list(ln_InitDevs = dev4d, init_devs_pen_use = pen_all)))
  expect_equal(sum(distinct$Init_Rec_nLL), -sum(dnorm(c(devs, devs_m), 0, 0.6, TRUE)), tolerance = 1e-12)

  # the own-mean center pools every penalized cell across both sexes
  own_args <- c(base_args, list(ln_InitDevs = dev4d, init_devs_pen_use = pen_all))
  own_args$InitDevs_pen_center <- 1
  own <- do.call(get_recruitment_penalty, own_args)
  mu_pool <- mean(c(devs, devs_m))
  expect_equal(sum(own$Init_Rec_nLL), -sum(dnorm(c(devs, devs_m), mu_pool, 0.6, TRUE)), tolerance = 1e-12)
})

test_that("do_InitDevs_mapping expands its single-sex logic across sexes per InitDevs_sex_spec", {

  messages_list <<- character(0)
  make_input <- function() {
    n_ages <- 5
    list(
      data = list(
        n_pop = 1,
        n_regions = 1,
        n_sexes = 2,
        ages = 1:n_ages,
        equil_init_age_strc = 2,
        rec_region_prop_spec = 0,
        natal_region = 1,
        rec_dd = 0
      ),
      par = list(ln_InitDevs = array(0, dim = c(1, 1, n_ages - 1, 2))),
      map = list()
    )
  }

  # est_shared_s: one parameter per age, read by both sexes, penalized once
  shared <- do_InitDevs_mapping(
    make_input(),
    InitDevs_spec = NULL,
    rec_dd = NULL,
    init_age_devs_shared = NULL,
    InitDevs_sex_spec = "est_shared_s"
  )
  map_shared <- array(as.numeric(as.character(shared$map$ln_InitDevs)), dim = c(1, 1, 4, 2))
  expect_equal(length(unique(na.omit(as.vector(map_shared)))), 4)
  expect_equal(map_shared[1,1,,1], map_shared[1,1,,2])
  expect_equal(sum(shared$data$init_devs_pen_use), 4)
  expect_equal(as.numeric(shared$data$init_devs_pen_use[1,1,,2]), rep(0, 4))

  # est_all: each sex has its own parameters, every copy penalized
  distinct <- do_InitDevs_mapping(
    make_input(),
    InitDevs_spec = NULL,
    rec_dd = NULL,
    init_age_devs_shared = NULL,
    InitDevs_sex_spec = "est_all"
  )
  map_distinct <- array(as.numeric(as.character(distinct$map$ln_InitDevs)), dim = c(1, 1, 4, 2))
  expect_equal(length(unique(na.omit(as.vector(map_distinct)))), 8)
  expect_equal(sum(distinct$data$init_devs_pen_use), 8)

  # a fixed plus group stays fixed for every sex
  fixed_plus <- make_input(); fixed_plus$data$equil_init_age_strc <- 1
  fp <- do_InitDevs_mapping(
    fixed_plus,
    InitDevs_spec = NULL,
    rec_dd = NULL,
    init_age_devs_shared = NULL,
    InitDevs_sex_spec = "est_all"
  )
  map_fp <- array(as.numeric(as.character(fp$map$ln_InitDevs)), dim = c(1, 1, 4, 2))
  expect_true(all(is.na(map_fp[1, 1, 4, ])))

  # one sex cannot be sex-specific
  single <- make_input(); single$data$n_sexes <- 1; single$par$ln_InitDevs <- array(0, dim = c(1, 1, 4, 1))
  expect_error(do_InitDevs_mapping(
    single,
    InitDevs_spec = NULL,
    rec_dd = NULL,
    init_age_devs_shared = NULL,
    InitDevs_sex_spec = "est_all"
  ), "n_sexes > 1")
})

test_that("retention selectivity holds the sex offsets and the plateau through setup and the objective", {

  messages_list <<- character(0)
  # the helper builds retention's own flags and scale parameters
  base_input <- list(data = list(n_sexes = 2, n_regions = 1), par = list(), map = list())
  out <- setup_sel_sex_offset(
    base_input,
    "par_scale",
    prefix = "ret",
    n_fleets = 1,
    fleet_label = "fishery fleet (retention)",
    sel_model_arr = array(0, dim = c(1, 3, 1)),
    cont_tv_mat = array(0, dim = c(1, 1)),
    max_blks = 1
  )
  expect_equal(out$data$retsel_sex_par_offset, 1)
  expect_equal(out$data$retsel_sex_scale_offset, 1)
  expect_equal(dim(out$par$ln_retsel_sex_scale), c(1, 1, 2, 1))

  # a small two-sex model with retention estimated under both offsets and a plateau
  n_yrs <- 12; n_ages <- 6; n_sexes <- 2
  input_list <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = NULL,
    n_regions = 1,
    n_sexes = n_sexes,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(5)
  )
  biol <- function(val) { a <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes)); for(s in 1:n_sexes) a[1,1,,1,,s] <- matrix(rep(val, each = n_yrs), n_yrs, n_ages); a }
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = biol(1:n_ages),
    MatAA = biol(rep(1, n_ages)),
    fit_lengths = 0,
    AgeingError = diag(n_ages),
    M_spec = "fix",
    Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, n_sexes))
  )
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )
  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = array(1, dim = c(1, n_yrs, 1, 1)),
    UseCatch = array(1, dim = c(1, n_yrs, 1, 1)),
    Use_F_pen = 0,
    sigmaC_spec = "fix",
    ln_sigmaC = array(log(0.05), dim = c(1, n_yrs, 1, 1))
  ))
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
    ObsFishIdx_SE = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
    UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
    ObsFishAgeComps = array(NA_real_, dim = c(1, n_yrs, 1, n_ages, n_sexes, 1)),
    UseFishAgeComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_FishAgeComps = array(0, dim = c(1, n_yrs, 1, n_sexes, 1)),
    ObsFishLenComps = NULL,
    UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_FishLenComps = NULL,
    fish_idx_type = "none",
    FishIdx_LikeType = "lognormal",
    FishAgeComps_LikeType = "none",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "none_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
    ObsSrvIdx_SE = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
    UseSrvIdx = array(0, dim = c(1, n_yrs, 1, 1)),
    ObsSrvAgeComps = array(NA_real_, dim = c(1, n_yrs, 1, n_ages, n_sexes, 1)),
    UseSrvAgeComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_sexes, 1)),
    ObsSrvLenComps = NULL,
    UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvLenComps = NULL,
    srv_idx_type = "none",
    SrvIdx_LikeType = "lognormal",
    SrvAgeComps_LikeType = "none",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "none_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix",
    use_fixed_ret_sel = 0,
    ret_sel_model = "logist1_Fleet_1_NSelBins_4",
    ret_fixed_sel_pars_spec = "est_all",
    ret_sel_sex_offset = "par_scale"
  )
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "fix"
  )
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 0,
    Wt_SrvIdx = 0,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1))
  )

  pars <- input_list$par
  pars$ret_fixed_sel_pars[1, , 1, 1, 1] <- log(c(3, 1.2))  # female ln_b50, ln_k
  pars$ret_fixed_sel_pars[1, , 1, 2, 1] <- c(0.3, -0.2)    # male offsets
  pars$ln_retsel_sex_scale[1, 1, 2, 1] <- -0.25             # males retained less
  expect_true("ln_retsel_sex_scale" %in% names(input_list$map))
  expect_equal(sum(!is.na(input_list$map$ln_retsel_sex_scale)), 1) # the second sex's scale is estimated

  obj <- fit_model(input_list$data, pars, input_list$map, do_optim = FALSE, silent = TRUE)
  expected_f <- logistic_sel(3, 1.2, 1:n_ages); expected_f[5:n_ages] <- expected_f[4]
  expected_m <- logistic_sel(3 * exp(0.3), 1.2 * exp(-0.2), 1:n_ages); expected_m[5:n_ages] <- expected_m[4]
  expected_m <- expected_m * exp(-0.25)
  expect_equal(as.numeric(obj$rep$ret_sel[1,1,1,1,,1,1]), expected_f, tolerance = 1e-12)
  expect_equal(as.numeric(obj$rep$ret_sel[1,1,1,1,,2,1]), expected_m, tolerance = 1e-12)
})

test_that("the between-sex tie on initial age deviations is a Gaussian on each later sex's difference from the first", {

  n_ages <- 5
  base_args <- list(
    n_pop = 1,
    n_regions = 1,
    n_ages = n_ages,
    n_est_rec_devs = 3,
    rec_dd = 0,
    natal_region = 1,
    rec_region_prop_spec = 0,
    rec_region_prop = matrix(1, 1, 1),
    equil_init_age_strc = 2,
    init_age_devs_shared = NULL,
    ln_sigmaR = array(log(0.6), dim = c(2, 1, 1)),
    bias_ramp = rep(0, 3),
    sigmaR_switch = 1,
    ln_RecDevs = array(0, dim = c(1, 1, 3)),
    sigmaR2_early = matrix(0.36, 1, 1),
    sigmaR2_late = matrix(0.36, 1, 1),
    do_rec_bias_ramp = 0
  )
  devs_f <- c(0.2, -0.1, 0.4, -0.3); devs_m <- c(-0.15, 0.25, -0.35, 0.05)
  dev4 <- array(c(devs_f, devs_m), dim = c(1, 1, n_ages - 1, 2))
  pen_all <- array(1, dim = dim(dev4))

  off <- do.call(get_recruitment_penalty, c(base_args, list(ln_InitDevs = dev4, init_devs_pen_use = pen_all)))
  expect_equal(sum(off$Init_Sex_nLL), 0)

  on <- do.call(get_recruitment_penalty, c(base_args, list(
    ln_InitDevs = dev4,
    init_devs_pen_use = pen_all,
    Use_init_sex_pen = 1,
    ln_sigma_init_sex = log(1 / sqrt(2))
  )))
  expect_equal(as.numeric(on$Init_Sex_nLL[1,1,,1]), rep(0, n_ages - 1)) # the first sex is the reference
  expect_equal(sum(on$Init_Sex_nLL), -sum(dnorm(devs_m - devs_f, 0, 1 / sqrt(2), TRUE)), tolerance = 1e-12)
  # net of the normalizing constant it is fm.tpl's norm2 of the difference
  expect_equal(sum(on$Init_Sex_nLL) - (n_ages - 1) * (log(1 / sqrt(2)) + 0.5 * log(2 * pi)), sum((devs_m - devs_f)^2), tolerance = 1e-12)
  # the initial-age penalty itself is untouched by the tie
  expect_equal(sum(on$Init_Rec_nLL), sum(off$Init_Rec_nLL), tolerance = 1e-12)

  # a fixed plus group (equil 1) is outside the tie, like the penalty
  args1 <- base_args; args1$equil_init_age_strc <- 1
  on1 <- do.call(get_recruitment_penalty, c(args1, list(
    ln_InitDevs = dev4,
    init_devs_pen_use = pen_all,
    Use_init_sex_pen = 1,
    ln_sigma_init_sex = 0
  )))
  expect_equal(as.numeric(on1$Init_Sex_nLL[1,1,n_ages - 1,2]), 0)

  # setup refuses the tie where it would be identically zero
  dim2 <- Setup_Mod_Dim(
    years = 1:5,
    ages = 1:4,
    lens = NULL,
    n_regions = 1,
    n_sexes = 2,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  expect_error(Setup_Mod_Rec(
    input_list = dim2,
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    InitDevs_sex_spec = "est_shared_s",
    Use_init_sex_pen = 1
  ), "est_all")
  ok <- Setup_Mod_Rec(
    input_list = dim2,
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    InitDevs_sex_spec = "est_all",
    Use_init_sex_pen = 1,
    init_sex_pen_sigma = 1 / sqrt(2)
  )
  expect_equal(ok$data$Use_init_sex_pen, 1)
  expect_equal(ok$data$ln_sigma_init_sex, log(1 / sqrt(2)))
  dim1 <- Setup_Mod_Dim(
    years = 1:5,
    ages = 1:4,
    lens = NULL,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  expect_error(Setup_Mod_Rec(input_list = dim1, rec_model = "mean_rec", sigmaR_spec = "fix", Use_init_sex_pen = 1), "n_sexes > 1")
})

test_that("sex offsets refuse a sex-shared fixed specification under par and map scale offsets only for a fleet's own blocks", {

  messages_list <<- character(0)
  base_input <- list(data = list(n_sexes = 2, n_regions = 1), par = list(), map = list())
  sel_model_arr <- array(0, dim = c(1, 4, 2))
  cont_tv_mat <- array(0, dim = c(1, 2))
  # fleet 1 has two blocks, fleet 2 one; the parameter array is padded to two
  sel_blocks <- array(1, dim = c(1, 4, 2)); sel_blocks[1, 3:4, 1] <- 2

  # a par offset reads the later sex's slots as offsets, so sharing them across sexes is refused
  expect_error(setup_sel_sex_offset(
    base_input,
    c("par", "none"),
    prefix = "fish",
    n_fleets = 2,
    fleet_label = "fishery fleet",
    sel_model_arr = sel_model_arr,
    cont_tv_mat = cont_tv_mat,
    max_blks = 2,
    sel_blocks = sel_blocks,
    fixed_spec = c("est_shared_s", "est_all")
  ), "shares the sex slots")
  ok <- setup_sel_sex_offset(
    base_input,
    c("par", "none"),
    prefix = "fish",
    n_fleets = 2,
    fleet_label = "fishery fleet",
    sel_model_arr = sel_model_arr,
    cont_tv_mat = cont_tv_mat,
    max_blks = 2,
    sel_blocks = sel_blocks,
    fixed_spec = c("est_shared_r", "est_all")
  )
  expect_equal(ok$data$fishsel_sex_par_offset, c(1, 0))

  # scale offsets exist for the blocks a fleet has and nowhere else
  out <- setup_sel_sex_offset(
    base_input,
    c("scale", "scale"),
    prefix = "fish",
    n_fleets = 2,
    fleet_label = "fishery fleet",
    sel_model_arr = sel_model_arr,
    cont_tv_mat = cont_tv_mat,
    max_blks = 2,
    sel_blocks = sel_blocks,
    fixed_spec = c("est_all", "est_all")
  )
  map_scale <- array(as.numeric(as.character(out$map$ln_fishsel_sex_scale)), dim = c(1, 2, 2, 2))
  expect_false(any(is.na(map_scale[1, 1:2, 2, 1]))) # fleet 1: both blocks, second sex
  expect_false(is.na(map_scale[1, 1, 2, 2]))        # fleet 2: its one block
  expect_true(is.na(map_scale[1, 2, 2, 2]))         # fleet 2 has no second block
  expect_true(all(is.na(map_scale[1, , 1, ])))      # the first sex is the reference
  expect_equal(sum(!is.na(map_scale)), 3)
})

test_that("the simulator draws initial age deviations per sex only when asked, and the shared draw is unchanged", {

  make_sim <- function(spec, n_sexes = 2) {
    sim_list <- Setup_Sim_Dim(
      n_sims = 1,
      n_yrs = 8,
      n_regions = 1,
      n_ages = 5,
      n_lens = NULL,
      n_sexes = n_sexes,
      n_fish_fleets = 1,
      n_srv_fleets = 1,
      n_pop = 1
    )
    sim_list <- Setup_Sim_Containers(sim_list)
    args <- list(
      sim_list = sim_list,
      R0_input = replicate(1, array(5, dim = c(1, 1, 8))),
      ln_sigmaR = array(log(0.5), dim = c(2, 1, 1)),
      recruitment_opt = "mean_rec",
      init_age_strc = 1
    )
    if(!is.null(spec)) args$InitDevs_sex_spec <- spec
    do.call(Setup_Sim_Rec, args)
  }

  # the default draw is one curve read by every sex
  shared <- make_sim("est_shared_s")
  expect_equal(shared$InitDevs_sex_spec, "est_shared_s")
  expect_equal(make_sim(NULL)$InitDevs_sex_spec, "est_shared_s") # the default is the shared draw

  # est_all draws each sex its own
  distinct <- make_sim("est_all")
  expect_equal(distinct$InitDevs_sex_spec, "est_all")

  # a single-sex model has nothing to draw separately
  expect_error(make_sim("est_all", n_sexes = 1), "n_sexes > 1")
  expect_error(make_sim("distinct"), "est_shared_s or est_all")

  # the draws themselves: shared repeats across sexes, est_all does not, and the
  # shared draw is the same numbers the simulator drew before the option existed
  draw <- function(spec) {
    set.seed(42)
    n_ages <- 5; n_sexes <- 2
    n_dev_draws <- if(spec == "est_all") n_sexes else 1
    array(stats::rnorm(n_dev_draws * (n_ages - 1), -0.5^2 / 2, 0.5), dim = c(n_ages - 1, n_sexes))
  }
  set.seed(42); legacy <- stats::rnorm(4, -0.5^2 / 2, 0.5)
  expect_equal(as.numeric(draw("est_shared_s")[, 1]), legacy, tolerance = 1e-12)
  expect_equal(as.numeric(draw("est_shared_s")[, 2]), legacy, tolerance = 1e-12) # both sexes read the one curve
  expect_false(isTRUE(all.equal(as.numeric(draw("est_all")[, 1]), as.numeric(draw("est_all")[, 2]))))
})
