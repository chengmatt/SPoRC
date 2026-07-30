library(SPoRC)
library(testthat)

test_that("Get_Movement works", {
  test_that("get_movement_dp_design_matrix returns correct dimensions", {
    dat <- data.frame(
      regions = rep(1:3, 2),
      years   = rep(1:2, each = 3),
      depth   = rnorm(6)
    )
    out <- get_movement_dp_design_matrix(
      data               = dat,
      preference_formula = ~ depth,
      diffusion_formula  = ~ depth
    )

    expect_equal(out$n_theta, 2L)           # intercept + depth
    expect_equal(out$n_gamma, 2L)
    expect_equal(nrow(out$X_zk), nrow(dat))
    expect_equal(nrow(out$W_zk), nrow(dat))
    expect_equal(ncol(out$X_zk), out$n_theta)
    expect_equal(ncol(out$W_zk), out$n_gamma)
  })

  # ── helpers ─────────────────────────────────────────────────────────────────

  make_fixed_move <- function(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes) {
    val <- 1 / n_regions
    array(val, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes))
  }

  call_fixed <- function(n_pop = 1, n_regions = 3, n_yrs = 2, n_seas = 1,
                         n_ages = 2, n_sexes = 1) {
    FM <- make_fixed_move(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)
    Get_Movement(
      move_type           = 0,
      do_recruits_move    = 1,
      n_pop               = n_pop,
      n_regions           = n_regions,
      n_yrs               = n_yrs,
      n_proj_yrs_devs     = 0,
      n_ages              = n_ages,
      n_sexes             = n_sexes,
      n_seas              = n_seas,
      move_pars           = array(0, c(n_pop, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
      move_devs           = array(0, c(n_pop, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
      use_fixed_movement  = 1,
      Fixed_Movement      = FM,
      log_move_diffusion_pars = numeric(0),
      move_preference_pars    = numeric(0),
      area_r                  = rep(1, n_regions),
      adjacency_mat           = matrix(1, n_regions, n_regions),
      ctmc_diffusion_bounds   = 0
    )
  }

  # ── fixed movement ───────────────────────────────────────────────────────────

  test_that("use_fixed_movement = 1 passes Fixed_Movement through unchanged", {
    res <- call_fixed()
    FM  <- make_fixed_move(1, 3, 2, 1, 2, 1)
    expect_equal(res$Movement, FM, ignore_attr = TRUE)
    expect_null(res$Mrate)
    expect_equal(res$move_pen, 0)
  })

  # ── unstructured multinomial logit movement (move_type = 0) ─────────────────

  make_unstructured_call <- function(move_pars_val = 0,
                                     n_pop = 1, n_regions = 3,
                                     n_yrs = 2, n_proj = 0,
                                     n_seas = 1, n_ages = 2, n_sexes = 1) {
    n_total <- n_yrs + n_proj
    mp <- array(move_pars_val, c(n_pop, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes))
    md <- array(0,             c(n_pop, n_regions, n_regions - 1, n_total, n_seas, n_ages, n_sexes))
    Get_Movement(
      move_type           = 0,
      do_recruits_move    = 1,
      n_pop               = n_pop,
      n_regions           = n_regions,
      n_yrs               = n_yrs,
      n_proj_yrs_devs     = n_proj,
      n_ages              = n_ages,
      n_sexes             = n_sexes,
      n_seas              = n_seas,
      move_pars           = mp,
      move_devs           = md,
      use_fixed_movement  = 0,
      Fixed_Movement      = NULL,
      log_move_diffusion_pars = numeric(0),
      move_preference_pars    = numeric(0),
      area_r                  = rep(1, n_regions),
      adjacency_mat           = matrix(1, n_regions, n_regions),
      ctmc_diffusion_bounds   = 0
    )
  }

  test_that("unstructured movement rows sum to 1 (equal logit values => uniform)", {
    res <- make_unstructured_call(move_pars_val = 0)
    for (p in 1:1) for (r in 1:3) for (y in 1:2) for (s in 1:1) for (a in 1:2) for (sx in 1:1) {
      expect_equal(sum(res$Movement[p, r, , y, s, a, sx]), 1,
                   tolerance = 1e-10, label = sprintf("row sum p%d r%d y%d", p, r, y))
    }
  })

  test_that("unstructured movement is uniform when all logit pars are zero", {
    res <- make_unstructured_call(move_pars_val = 0)
    expect_equal(unique(as.numeric(res$Movement)), 1/3, tolerance = 1e-10)
  })

  test_that("unstructured movement output dimensions are correct", {
    res <- make_unstructured_call(n_pop = 1, n_regions = 3,
                                  n_yrs = 3, n_proj = 1,
                                  n_ages = 4, n_sexes = 2)
    expect_equal(unname(dim(res$Movement)), c(1, 3, 3, 4, 1, 4, 2))
  })

  test_that("unstructured movement: projection years use last historical year's pars", {
    mp <- array(0, c(1, 3, 2, 2, 1, 2, 1))
    mp[1, , , 2, , , ] <- 5
    md <- array(0, c(1, 3, 2, 3, 1, 2, 1))

    res <- Get_Movement(
      move_type = 0, do_recruits_move = 1,
      n_pop = 1, n_regions = 3, n_yrs = 2, n_proj_yrs_devs = 1,
      n_ages = 2, n_sexes = 1, n_seas = 1,
      move_pars = mp, move_devs = md,
      use_fixed_movement = 0, Fixed_Movement = NULL,
      log_move_diffusion_pars = numeric(0), move_preference_pars = numeric(0),
      area_r = rep(1, 3), adjacency_mat = matrix(1, 3, 3), ctmc_diffusion_bounds = 0
    )
    expect_equal(res$Movement[, , , 3, , , ], res$Movement[, , , 2, , , ])
    expect_false(isTRUE(all.equal(res$Movement[, , , 3, , , ], res$Movement[, , , 1, , , ])))
  })

  # ── CTMC movement (move_type = 1) ────────────────────────────────────────────

  make_ctmc_dat <- function(n_pop = 1, n_regions = 3, n_yrs = 2,
                            n_seas = 1, n_ages = 2, n_sexes = 1) {
    expand.grid(
      pop     = seq_len(n_pop),
      regions = seq_len(n_regions),
      years   = seq_len(n_yrs),
      seas    = seq_len(n_seas),
      ages    = seq_len(n_ages),
      sexes   = seq_len(n_sexes)
    )
  }

  make_adjacency <- function(n_regions) {
    # Off-diagonal = 1 (adjacent); diagonal = 0 so D_ss diagonal correction is valid
    mat <- matrix(1L, n_regions, n_regions)
    diag(mat) <- 0L
    mat
  }

  make_ctmc_call <- function(n_pop = 1, n_regions = 3,
                             n_yrs = 2, n_proj = 0,
                             n_seas = 1, n_ages = 2, n_sexes = 1,
                             pref_pars = NULL,
                             diff_pars = NULL,
                             ctmc_diffusion_bounds = 0,
                             do_recruits_move = 1) {
    n_total <- n_yrs + n_proj
    dat <- make_ctmc_dat(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)

    # Intercept-only formulas (1 column) => pars must be length 1
    if (is.null(pref_pars)) pref_pars <- 0
    if (is.null(diff_pars)) diff_pars <- 0

    md <- array(0, c(n_pop, n_regions, n_regions - 1, n_total, n_seas, n_ages, n_sexes))

    Get_Movement(
      move_type               = 1,
      do_recruits_move        = do_recruits_move,
      n_pop                   = n_pop,
      n_regions               = n_regions,
      n_yrs                   = n_yrs,
      n_proj_yrs_devs         = n_proj,
      n_ages                  = n_ages,
      n_sexes                 = n_sexes,
      n_seas                  = n_seas,
      move_pars               = array(0, c(n_pop, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
      move_devs               = md,
      use_fixed_movement      = 0,
      Fixed_Movement          = NULL,
      ctmc_move_dat           = dat,
      preference_formula      = ~ 1,
      diffusion_formula       = ~ 1,
      log_move_diffusion_pars = diff_pars,
      move_preference_pars    = pref_pars,
      area_r                  = rep(1, n_regions),
      adjacency_mat           = make_adjacency(n_regions),
      ctmc_diffusion_bounds   = ctmc_diffusion_bounds
    )
  }

  test_that("CTMC movement: output dimensions are correct", {
    res <- make_ctmc_call(n_pop = 1, n_regions = 3, n_yrs = 2,
                          n_seas = 1, n_ages = 2, n_sexes = 1)
    expect_equal(unname(dim(res$Movement)), c(1, 3, 3, 2, 1, 2, 1))
    expect_equal(unname(dim(res$Mrate)),    c(1, 3, 3, 2, 1, 2, 1))
  })

  test_that("CTMC movement: Movement rows sum to 1", {
    res <- make_ctmc_call()
    for (p in 1:1) for (r in 1:3) for (y in 1:2) for (s in 1:1) for (a in 1:2) for (sx in 1:1) {
      expect_equal(sum(res$Movement[p, r, , y, s, a, sx]), 1,
                   tolerance = 1e-9,
                   label = sprintf("row sum p%d r%d y%d a%d", p, r, y, a))
    }
  })

  test_that("CTMC movement: all Movement entries are in [0, 1]", {
    res <- make_ctmc_call()
    expect_true(all(res$Movement >= 0 - 1e-10))
    expect_true(all(res$Movement <= 1 + 1e-10))
  })

  test_that("CTMC movement: zero preference gives move_pen = 0", {
    res_zero <- make_ctmc_call(pref_pars = 0)
    expect_equal(res_zero$move_pen, 0)
  })

  test_that("CTMC movement: non-zero preference increments move_pen", {
    res_nonzero <- make_ctmc_call(pref_pars = 1)
    expect_gt(res_nonzero$move_pen, 0)
  })

  test_that("CTMC movement: Mrate rows sum to ~ 0 (generator property)", {
    # Mrate stores t(Q_ss), so dim [pop, from, to, ...]: 'from' = rows of Q.
    # Rows of a valid generator Q sum to 0.
    res <- make_ctmc_call()
    for (p in 1:1) for (y in 1:2) for (s in 1:1) for (a in 1:2) for (sx in 1:1) {
      Q <- res$Mrate[p, , , y, s, a, sx]
      expect_equal(unname(rowSums(Q)), rep(0, nrow(Q)), tolerance = 1e-9,
                   label = sprintf("Q row sums y%d a%d", y, a))
    }
  })

  test_that("CTMC movement: projection year lookup is capped at n_yrs", {
    res <- make_ctmc_call(n_yrs = 2, n_proj = 1)
    expect_equal(unname(dim(res$Movement))[4], 3)
    expect_equal(res$Movement[, , , 3, , , ], res$Movement[, , , 2, , , ])
  })

  test_that("CTMC movement: do_recruits_move = 0 leaves age-1 movement as zeros", {
    res <- make_ctmc_call(do_recruits_move = 0, n_ages = 3)
    expect_equal(
      as.numeric(res$Movement[, , , , , 1, ]),
      rep(0, length(res$Movement[, , , , , 1, ]))
    )
    # age 2 rows should sum to 1
    expect_equal(sum(res$Movement[1, 1, , 1, 1, 2, 1]), 1, tolerance = 1e-9)
  })

  test_that("CTMC movement: ctmc_diffusion_bounds = 1 still yields valid movement", {
    # pref_pars length 1 to match intercept-only ~ 1 formula
    res <- make_ctmc_call(ctmc_diffusion_bounds = 1, pref_pars = 2)
    for (p in 1:1) for (r in 1:3) for (y in 1:2) for (s in 1:1) for (a in 1:2) for (sx in 1:1) {
      expect_equal(sum(res$Movement[p, r, , y, s, a, sx]), 1,
                   tolerance = 1e-9,
                   label = sprintf("bounds row sum r%d y%d a%d", r, y, a))
    }
  })

  test_that("CTMC movement: Mrate is NULL when use_fixed_movement = 1", {
    res <- call_fixed()
    expect_null(res$Mrate)
  })

  test_that("CTMC movement: Mrate is non-NULL for move_type = 1", {
    res <- make_ctmc_call()
    expect_false(is.null(res$Mrate))
  })

})
