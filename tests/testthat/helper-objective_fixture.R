# A small single-region operating model with retention and discarding, used to switch on
# the priors, penalties, and selectivity forms that no other test turns on. Discarding is
# needed because the discard mortality rate penalty only applies to cells with discard
# data, and retention is needed because the retained-selectivity prior acts on
# ret_fixed_sel_pars.
#
# Two versions exist: an age-only one, and one carrying length bins for the tests that
# estimate selectivity over lengths. Each is simulated once and reused, since none of these
# tests change the data, only which term the estimation model applies to it.

#' Size-age transition array for the length-structured fixture
#'
#' Spreads each age across the length bins with a normal kernel, normalized so the lengths
#' for a given age sum to one. Selectivity at age is then a weighted average of selectivity
#' at length.
#'
#' @keywords internal
objective_fixture_sizeage <- function(n_lens, n_ages) {
  vapply(seq_len(n_ages), function(a) {
    center <- 1 + (a - 1) * (n_lens - 1) / (n_ages - 1)
    w <- stats::dnorm(seq_len(n_lens), center, 1)
    w / sum(w)
  }, numeric(n_lens))
}


objective_fixture_sim <- local({
  cached <- list()

  function(n_lens = NULL) {
    key <- if(is.null(n_lens)) "age" else paste0("len", n_lens)
    if(!is.null(cached[[key]])) return(cached[[key]])

    set.seed(123)
    sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = 30, n_regions = 1, n_ages = 10,
                              n_lens = n_lens, n_sexes = 1, n_fish_fleets = 1,
                              n_srv_fleets = 1, n_pop = 1)
    sim_list <- Setup_Sim_Containers(sim_list)

    age_curve <- function(slope, infl, scale = 1) {
      array(rep(scale / (1 + exp(-slope * ((1:sim_list$n_ages) - infl))), each = sim_list$n_yrs),
            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                    sim_list$n_ages, sim_list$n_sexes, 1))
    }

    sim_list <- Setup_Sim_Fishing(
      sim_list = sim_list,
      fish_sel_input = replicate(sim_list$n_sims, age_curve(3, 2)),
      # retention rises with age and tops out at 0.5, so a real fraction is discarded
      ret_sel_input = replicate(sim_list$n_sims, age_curve(3, 5, scale = 0.5)),
      dmr_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                     sim_list$n_fish_fleets, sim_list$n_sims))
    )

    sim_list <- Setup_Sim_Survey(
      sim_list = sim_list,
      srv_sel_input = replicate(sim_list$n_sims, age_curve(1, 3))
    )

    biol_dim <- c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                  sim_list$n_ages, sim_list$n_sexes)
    waa <- array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs), dim = biol_dim)
    mat <- array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs), dim = biol_dim)

    size_age <- NULL
    if(!is.null(n_lens)) {
      one_year <- objective_fixture_sizeage(n_lens, sim_list$n_ages)
      size_age <- array(0, dim = c(biol_dim[1:4], n_lens, sim_list$n_ages, sim_list$n_sexes,
                                   sim_list$n_sims))
      for(y in seq_len(sim_list$n_yrs)) size_age[1, 1, y, 1, , , 1, 1] <- one_year
    }

    sim_list <- suppressWarnings(Setup_Sim_Biologicals(
      sim_list = sim_list,
      natmort_input = replicate(sim_list$n_sims, array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                   sim_list$n_yrs, sim_list$n_ages,
                                                                   sim_list$n_sexes))),
      WAA_input = replicate(sim_list$n_sims, waa),
      WAA_fish_input = replicate(sim_list$n_sims, array(waa, dim = c(biol_dim, sim_list$n_fish_fleets))),
      WAA_srv_input = replicate(sim_list$n_sims, array(waa, dim = c(biol_dim, sim_list$n_srv_fleets))),
      MatAA_input = replicate(sim_list$n_sims, mat),
      SizeAgeTrans_input = size_age
    ))

    sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
    sim_list$Movement <- array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
                                          sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                                          sim_list$n_sexes, sim_list$n_sims))

    sim_list <- Setup_Sim_Rec(
      sim_list = sim_list,
      R0_input = replicate(sim_list$n_sims, array(5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs))),
      ln_sigmaR = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_region)),
      recruitment_opt = "mean_rec", init_age_strc = 1
    )

    set.seed(777)
    cached[[key]] <<- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
    cached[[key]]
  }
})


#' Estimation model for the objective branch fixture
#'
#' Builds the input list with every prior and penalty switched off, then applies any
#' overrides so a single term can be turned on and its effect on jnLL isolated.
#'
#' @param rec,catch_f,fishsel,srvsel Named lists of arguments overriding the defaults passed
#'   to \code{Setup_Mod_Rec}, \code{Setup_Mod_Catch_and_F}, \code{Setup_Mod_Fishsel_and_Q},
#'   and \code{Setup_Mod_Srvsel_and_Q}.
#' @param n_lens Number of length bins. When given, the fixture carries length
#'   compositions and fits them, which selectivity estimated over lengths requires.
#'
#' @keywords internal
objective_fixture_input <- function(rec = list(), catch_f = list(), fishsel = list(),
                                    srvsel = list(), n_lens = NULL) {

  fit_lengths <- as.integer(!is.null(n_lens))
  sim_obj <- objective_fixture_sim(n_lens)
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)

  input_list <- Setup_Mod_Dim(
    years = 1:sim_obj$n_years, ages = 1:sim_obj$n_ages,
    # the simulation records the number of length bins; the model wants the bins themselves
    lens = if(is.null(n_lens)) sim_obj$n_lens else seq_len(n_lens),
    n_regions = sim_obj$n_regions, n_sexes = sim_obj$n_sexes,
    n_fish_fleets = sim_obj$n_fish_fleets, n_srv_fleets = sim_obj$n_srv_fleets,
    n_pop = sim_obj$n_pop, natal_region = sim_obj$natal_region, verbose = FALSE
  )

  n_yrs <- length(input_list$data$years)
  n_regions <- input_list$data$n_regions
  n_seas <- input_list$data$n_seas
  n_sexes <- input_list$data$n_sexes
  n_fish_fleets <- input_list$data$n_fish_fleets
  n_srv_fleets <- input_list$data$n_srv_fleets

  input_list <- do.call(Setup_Mod_Rec, modifyList(list(
    input_list = input_list, do_rec_bias_ramp = 0, sigmaR_switch = 1,
    ln_sigmaR = array(log(0.5), c(2, input_list$data$n_pop, n_regions)),
    rec_model = "mean_rec", sigmaR_spec = "fix",
    init_age_strc = 1, equil_init_age_strc = 2, ln_global_R0 = log(5)
  ), rec))

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list, WAA = sim_data$WAA, MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv,
    fit_lengths = fit_lengths,
    SizeAgeTrans = if(fit_lengths == 1) sim_data$SizeAgeTrans else NA,
    AgeingError = sim_data$AgeingError, M_spec = "fix",
    Fixed_natmort = array(0.3, dim = c(input_list$data$n_pop, n_regions, n_yrs,
                                       length(input_list$data$ages), n_sexes))
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)

  input_list <- do.call(Setup_Mod_Catch_and_F, modifyList(list(
    input_list = input_list, ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch,
    Use_F_pen = 1, sigmaC_spec = "fix", ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(n_regions, n_seas, n_fish_fleets)),
    ObsDiscard = sim_data$ObsDiscard, UseDiscard = sim_data$UseDiscard,
    sigma_dmr_spec = "fix", dmr_mean_spec = "est_all", ln_sigmaD = sim_data$ln_sigmaD
  ), catch_f))

  len_like <- if(fit_lengths == 1) "Multinomial" else "none"
  len_type <- if(fit_lengths == 1) "agg_Year_1-terminal_Fleet_1" else "none_Year_1-terminal_Fleet_1"

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = sim_data$UseFishIdx,
    ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = sim_data$UseFishAgeComps, UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = sim_data$ISS_FishLenComps,
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial", FishLenComps_LikeType = len_like,
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = len_type
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps, UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial", SrvLenComps_LikeType = len_like,
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = len_type
  )

  input_list <- do.call(Setup_Mod_Fishsel_and_Q, modifyList(list(
    input_list = input_list,
    fish_sel_model = "logist1_Fleet_1", fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "est_all",
    ret_sel_model = "asymplogist1_Fleet_1", ret_fixed_sel_pars_spec = "est_all",
    use_fixed_ret_sel = 0
  ), fishsel))

  input_list <- do.call(Setup_Mod_Srvsel_and_Q, modifyList(list(
    input_list = input_list, srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "est_all"
  ), srvsel))

  fish_wt <- array(1, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
  srv_wt <- array(1, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))

  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1,
    Wt_Discard = 1, Wt_D = 1,
    Wt_FishAgeComps = fish_wt, Wt_FishLenComps = fish_wt,
    Wt_SrvAgeComps = srv_wt, Wt_SrvLenComps = srv_wt
  )
}


#' Evaluate the objective without optimizing
#'
#' These tests compare jnLL at a common set of parameter values, so the model is only ever
#' evaluated, never fitted.
#'
#' @keywords internal
evaluate_input <- function(input_list) {
  fit_model(input_list$data, input_list$par, input_list$map,
            random = NULL, silent = TRUE, do_optim = FALSE)
}
