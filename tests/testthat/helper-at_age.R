# A minimal at-age model with configurable regions, sexes and fishery fleets.
# The at-age streams are stored over regions and sexes whatever a fleet reports,
# so the tests that pin those margins need a model that has them.

build_at_age <- function(n_yrs = 20, n_ages = 5, n_regions = 1, n_sexes = 1, n_fleets = 1,
                         ObsCatchAA = NULL, UseCatchAA = NULL, srv_extra = list(), ...) {

  yrs <- seq_len(n_yrs); ages <- seq_len(n_ages)
  d1 <- c(1, n_regions, n_yrs, 1, n_ages, n_sexes)
  fl <- seq_len(n_fleets)

  il <- Setup_Mod_Dim(n_pop = 1, years = yrs, ages = ages, lens = NA,
                      n_regions = n_regions, n_sexes = n_sexes, n_seas = 1,
                      n_fish_fleets = n_fleets, n_srv_fleets = 1, verbose = FALSE)
  il <- Setup_Mod_Rec(il, rec_model = "mean_rec", sigmaR_spec = "fix",
                      do_rec_bias_ramp = 0, init_age_strc = 1, ln_global_R0 = log(1e6))
  il <- suppressWarnings(Setup_Mod_Biologicals(
    il, WAA = array(1, dim = d1), WAA_fish = array(1, dim = c(d1, n_fleets)),
    WAA_srv = array(1, dim = c(d1, 1)), MatAA = array(1, dim = d1),
    fit_lengths = 0, M_spec = "fix",
    Fixed_natmort = array(0.2, dim = c(1, n_regions, n_yrs, n_ages, n_sexes))))
  il <- Setup_Mod_Movement(il, use_fixed_movement = 1,
                           Fixed_Movement = if(n_regions == 1) NA else
                             array(1 / n_regions, dim = c(1, n_regions, n_regions, n_yrs, 1, n_ages, n_sexes)),
                           do_recruits_move = 0)
  il <- Setup_Mod_Tagging(il, use_conv_fish_tagging = 0)

  il <- suppressWarnings(Setup_Mod_Catch_and_F(
    il,
    ObsCatch = array(1e4, dim = c(n_regions, n_yrs, 1, n_fleets)),
    UseCatch = array(0, dim = c(n_regions, n_yrs, 1, n_fleets)),
    ObsCatchAA = ObsCatchAA, UseCatchAA = UseCatchAA,
    sigmaC_spec = "fix", sigmaF_spec = "fix", ...))

  il <- Setup_Mod_FishIdx_and_Comps(
    il, ObsFishIdx = array(NA, dim = c(n_regions, n_yrs, 1, n_fleets)),
    ObsFishIdx_SE = array(NA, dim = c(n_regions, n_yrs, 1, n_fleets)),
    UseFishIdx = array(0, dim = c(n_regions, n_yrs, 1, n_fleets)),
    ObsFishAgeComps = array(0, dim = c(n_regions, n_yrs, 1, n_ages, n_sexes, n_fleets)),
    UseFishAgeComps = array(0, dim = c(n_regions, n_yrs, 1, n_fleets)),
    ISS_FishAgeComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fleets)),
    ObsFishLenComps = array(0, dim = c(n_regions, n_yrs, 1, 1, n_sexes, n_fleets)),
    UseFishLenComps = array(0, dim = c(n_regions, n_yrs, 1, n_fleets)),
    ISS_FishLenComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fleets)),
    fish_idx_type = rep("none", n_fleets),
    FishAgeComps_LikeType = rep("none", n_fleets), FishLenComps_LikeType = rep("none", n_fleets),
    FishAgeComps_Type = paste0("agg_Year_1-terminal_Fleet_", fl),
    FishLenComps_Type = paste0("agg_Year_1-terminal_Fleet_", fl))

  srv_args <- list(
    input_list = il, ObsSrvIdx = array(1e5, dim = c(n_regions, n_yrs, 1, 1)),
    ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_yrs, 1, 1)),
    UseSrvIdx = array(1, dim = c(n_regions, n_yrs, 1, 1)),
    ObsSrvAgeComps = array(0, dim = c(n_regions, n_yrs, 1, n_ages, n_sexes, 1)),
    UseSrvAgeComps = array(0, dim = c(n_regions, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, 1)),
    ObsSrvLenComps = array(0, dim = c(n_regions, n_yrs, 1, 1, n_sexes, 1)),
    UseSrvLenComps = array(0, dim = c(n_regions, n_yrs, 1, 1)),
    ISS_SrvLenComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, 1)),
    srv_idx_type = "abd", SrvAgeComps_LikeType = "none", SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "agg_Year_1-terminal_Fleet_1")
  # the survey index at age is set here, so its arguments arrive separately
  if(any(names(srv_extra) == "UseSrvIdxAA")) srv_args$UseSrvIdx[] <- 0
  il <- do.call(Setup_Mod_SrvIdx_and_Comps, c(srv_args[setdiff(names(srv_args), names(srv_extra))], srv_extra))

  il <- Setup_Mod_Fishsel_and_Q(
    il, cont_tv_fish_sel = paste0("none_Fleet_", fl), fish_sel_blocks = paste0("none_Fleet_", fl),
    fish_sel_model = paste0("logist1_Fleet_", fl), fish_q_blocks = paste0("none_Fleet_", fl),
    fish_fixed_sel_pars_spec = rep("est_all", n_fleets), fish_q_spec = rep("fix", n_fleets))
  il <- Setup_Mod_Srvsel_and_Q(
    il, cont_tv_srv_sel = "none_Fleet_1", srv_sel_blocks = "none_Fleet_1",
    srv_sel_model = "logist1_Fleet_1", srv_q_blocks = "none_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "est_all")

  Setup_Mod_Weighting(il, Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1,
                      Wt_F = 1, Wt_Tagging = 0,
                      Wt_FishAgeComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fleets)),
                      Wt_FishLenComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, n_fleets)),
                      Wt_SrvAgeComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, 1)),
                      Wt_SrvLenComps = array(0, dim = c(n_regions, n_yrs, 1, n_sexes, 1)))
}

at_age_rep <- function(il) fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep

# Reference multivariate normal log density, so the correlation tests can state
# the covariance they expect without taking a dependency on another package.
dmvn_ref <- function(x, Sigma) {
  R <- chol(Sigma)
  z <- backsolve(R, x, transpose = TRUE)
  return(-0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(z^2)))
}
