# The operating model and the estimation model must agree on the variance each recruitment deviation
# has, and on the correction that variance implies. Deterministic, to 1e-8.

data("dusky_rtmb_model")

# the estimation model's correction, recovered by profiling the recruitment penalty. the penalty is
# quadratic in a common deviation value, so three evaluations locate its minimum outright, and the
# value every deviation would have to take to sit there is minus the correction the penalty centers on
em_correction <- function(model, sd_line = 0.4, rho = 0.5, switch_yr = 1, ramp_off = FALSE, dsem_variance = "conditional") {

  b <- list(data = dusky_rtmb_model$data, par = dusky_rtmb_model$parameters,
            map = dusky_rtmb_model$mapping, verbose = FALSE, store_config = FALSE)
  n_yrs <- length(b$data$years)

  args <- list(input_list = b, sigmaR_switch = switch_yr, rec_model = "mean_rec", init_age_strc = 1,
               ln_global_R0 = log(2.7), t_spawn = b$data$t_spawn, sigmaR_spec = "fix",
               ln_sigmaR = array(log(sd_line), dim = c(2, b$data$n_pop, b$data$n_regions)),
               RecDevs_model = if(model == "dsem") "dsem" else model)
  if(ramp_off) { args$do_rec_bias_ramp <- 1; args$bias_year <- rep(n_yrs, 4) } else args$do_rec_bias_ramp <- 0
  if(model == "ar1") {
    args$RecDevs_rho <- array(atanh(rho), dim = c(b$data$n_pop, b$data$n_regions))
    args$RecDevs_rho_spec <- "fix"
  }
  il <- suppressMessages(suppressWarnings(do.call(Setup_Mod_Rec, args)))

  if(model == "dsem") {
    il <- suppressMessages(Setup_Mod_DSEM(il, dsem_data = NULL, dsem_variance = dsem_variance,
          dsem_arrows = c(sprintf("rec -> rec, 1, rho, %.17g", rho), sprintf("rec <-> rec, 0, sd_rec, %.17g", sd_line))))
  }

  il$data$bias_correct_pe <- 1 # "rec", the packaged data list predates both switches
  il$data$bias_correct_oe <- 0
  dat <- sync_dev_map_data(il$data, il$map)

  penalty_at <- function(x) {
    par <- il$par
    par$ln_RecDevs[] <- x
    obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, dat), par, map = il$map, silent = TRUE)
    sum(obj$report(obj$par)$Rec_nLL)
  }

  h <- 0.05
  f_lo <- penalty_at(-h)
  f_mid <- penalty_at(0)
  f_hi <- penalty_at(h)
  h * (f_hi - f_lo) / (2 * (f_hi - 2 * f_mid + f_lo)) # minus the quadratic's minimizer, which is the correction
}

# the marginal variance the objective settled on for each linked recruitment cell. half of it is the
# line the objective moves the dsem's mean grid by, which is that branch's correction
em_dsem_margvar <- function(sd_line = 0.4, rho = 0.5, dsem_variance = "conditional") {

  b <- list(data = dusky_rtmb_model$data, par = dusky_rtmb_model$parameters,
            map = dusky_rtmb_model$mapping, verbose = FALSE, store_config = FALSE)
  il <- suppressMessages(suppressWarnings(Setup_Mod_Rec(
    input_list = b, sigmaR_switch = 1, rec_model = "mean_rec", init_age_strc = 1,
    ln_global_R0 = log(2.7), t_spawn = b$data$t_spawn, sigmaR_spec = "fix",
    ln_sigmaR = array(log(sd_line), dim = c(2, b$data$n_pop, b$data$n_regions)),
    do_rec_bias_ramp = 0, RecDevs_model = "dsem")))
  il <- suppressMessages(Setup_Mod_DSEM(il, dsem_data = NULL, dsem_variance = dsem_variance,
        dsem_arrows = c(sprintf("rec -> rec, 1, rho, %.17g", rho), sprintf("rec <-> rec, 0, sd_rec, %.17g", sd_line))))

  il$data$bias_correct_pe <- 1
  il$data$bias_correct_oe <- 0
  dat <- sync_dev_map_data(il$data, il$map)
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, dat), il$par, map = il$map, silent = TRUE)
  rec_col <- il$data$dsem_link_col[which(il$data$dsem_link_par == "ln_RecDevs")]
  as.numeric(obj$report(obj$par)$dsem_margvar_grid[, rec_col])
}

test_that("an iid recruitment deviation carries half its own variance", {

  expect_equal(em_correction("iid", sd_line = 0.4), 0.5 * 0.4^2, tolerance = 1e-8) # the sd squared

  # and none of it once the ramp is moved past the last year
  expect_equal(em_correction("iid", ramp_off = TRUE), 0, tolerance = 1e-10)

})

test_that("an ar1 recruitment deviation carries half its stationary variance", {

  for(rho in c(0, 0.5, 0.9)) {
    expect_equal(em_correction("ar1", sd_line = 0.4, rho = rho), 0.5 * 0.4^2 / (1 - rho^2), tolerance = 1e-8)
  }

  # at rho zero it is the independent case, and the ramp switches it off
  expect_equal(em_correction("ar1", rho = 0), em_correction("iid"), tolerance = 1e-10)
  expect_equal(em_correction("ar1", ramp_off = TRUE), 0, tolerance = 1e-10)

})

test_that("a dsem recruitment deviation carries half its marginal variance under the arrows", {

  sd_line <- 0.4
  rho <- 0.5

  # conditional: the sd line is the innovation sd, so the settled correction exceeds half its square
  cond <- 0.5 * em_dsem_margvar(sd_line = sd_line, rho = rho, dsem_variance = "conditional")
  expect_equal(cond[1], 0.5 * sd_line^2, tolerance = 1e-8)                       # year one holds only its own shock
  expect_equal(cond[length(cond)], 0.5 * sd_line^2 / (1 - rho^2), tolerance = 1e-4) # settled on the stationary value
  expect_gt(cond[length(cond)], cond[1])

  # marginal: every cell sits on the sd line, so the correction is flat and matches the iid form
  marg <- 0.5 * em_dsem_margvar(sd_line = sd_line, rho = rho, dsem_variance = "marginal")
  expect_equal(unique(round(marg, 8)), 0.5 * sd_line^2, tolerance = 1e-8)

  # a linked cell's penalty is the dsem density rather than Rec_nLL, so it cannot be profiled the way
  # the other models are. test-dsem_objective.R holds the equivalence instead: a bare sd line at sigmaR
  # gives the same objective as the iid penalty, center included
})

test_that("a random walk carries no correction, having no stationary variance", {

  expect_equal(em_correction("rw"), 0, tolerance = 1e-10)
  # and a walk refuses the ramp outright, since there is no correction to ramp
  expect_error(em_correction("rw", ramp_off = TRUE), "no stationary variance")

})

test_that("the operating model applies the same correction the estimation model adds back", {

  sd_line <- 0.4
  rho <- 0.5
  n_yrs <- 30
  n_ages <- 10
  n_sims <- 2

  om_devs <- function(model, bias_correct_pe) {
    logi <- function(k, a50) 1 / (1 + exp(-k * ((1:n_ages) - a50)))
    arr7 <- function(v, ny, nf) array(rep(v, each = ny), dim = c(1, 1, ny, 1, n_ages, 1, nf))
    arr6 <- function(v) array(rep(v, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
    sl <- Setup_Sim_Dim(n_sims = n_sims, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = 1,
                        n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 1, n_pop = 1,
                        bias_correct_pe = bias_correct_pe)
    sl <- Setup_Sim_Containers(sl)
    sl <- Setup_Sim_Fishing(sl, fish_sel_input = replicate(n_sims, arr7(logi(3, 5), n_yrs, 1)))
    sl <- Setup_Sim_Survey(sl, srv_sel_input = replicate(n_sims, arr7(logi(1, 3), n_yrs, 1)),
                           srv_q_input = array(1, dim = c(1, n_yrs, 1, n_sims)))
    sl <- Setup_Sim_Biologicals(sl, natmort_input = replicate(n_sims, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
                               WAA_input = replicate(n_sims, arr6(5 * logi(3, 3))),
                               WAA_fish_input = replicate(n_sims, arr7(5 * logi(3, 3), n_yrs, 1)),
                               WAA_srv_input = replicate(n_sims, arr7(5 * logi(3, 3), n_yrs, 1)),
                               MatAA_input = replicate(n_sims, arr6(logi(3, 3))))
    sl <- Setup_Sim_Tagging(sl, use_conv_fish_tagging = 0)
    sl$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, n_sims))
    args <- list(sim_list = sl, R0_input = replicate(n_sims, array(5, dim = c(1, 1, n_yrs))),
                 ln_sigmaR = array(log(sd_line), dim = c(2, 1, 1)), recruitment_opt = "mean_rec", init_age_strc = 1,
                 RecDevs_model = model)
    if(model == "ar1") args$RecDevs_rho <- array(rho, dim = c(1, 1))
    sl <- suppressMessages(suppressWarnings(do.call(Setup_Sim_Rec, args)))
    set.seed(11)
    om <- suppressWarnings(suppressMessages(Simulate_Pop_Static(sim_list = sl, output_path = NULL)))
    as.numeric(om$ln_RecDevs[1, 1, , ])
  }

  # the same replicate drawn twice off one seed, once centered and once not. rnorm shifts its mean and
  # keeps its deviate, so the gap between the two draws is the correction outright rather than a mean of it
  om_correction <- function(model) om_devs(model, "none") - om_devs(model, "rec")

  expect_equal(unique(round(om_correction("iid"), 10)), 0.5 * sd_line^2, tolerance = 1e-8)
  expect_equal(unique(round(om_correction("ar1"), 10)), 0.5 * sd_line^2 / (1 - rho^2), tolerance = 1e-8)
  expect_equal(unique(round(om_correction("rw"), 10)), 0, tolerance = 1e-12)

  # and both sides land on the same number, which is the point
  expect_equal(unique(round(om_correction("iid"), 8)), round(em_correction("iid"), 8), tolerance = 1e-8)
  expect_equal(unique(round(om_correction("ar1"), 8)), round(em_correction("ar1", rho = rho), 8), tolerance = 1e-8)

})
