# Reference points for two populations that home and stray: global SPR against year classes followed for 800 years
# (1e-6), and global SPR and local MSY against where a 300 year projection settles (1e-5). Two regions, spawning part
# way through the year, and mortality, movement and homing that differ by population, region and age.

# Setup ----------------------------------------------------------------------

# two populations, one born in each region, fished at a different rate in each region
refpt_pop_model <- local({
  cache <- NULL
  function() {
    if(!is.null(cache)) return(cache)

    n_pop <- 2
    n_regions <- 2
    n_ages <- 6
    n_yrs <- 12

    natmort <- array(0, dim = c(n_pop, n_regions, n_yrs, n_ages, 1)) # per year
    spawn_move <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs, n_ages, 1)) # natal homing [pop, from, to, year, age, sex]
    Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs, 1, n_ages, 1)) # [pop, from, to, year, season, age, sex]

    for(p in 1:n_pop) {
      for(a in 1:n_ages) {

        for(r in 1:n_regions) natmort[p,r,,a,] <- c(0.2, 0.3)[r] * c(1, 1.15)[p] * c(1.5, 1.2, 1, 0.9, 0.8, 0.7)[a]

        # a share of the fish outside the natal region return to it to spawn, older fish more
        home <- p
        away <- 3 - p
        returning <- 0.3 + 0.08 * a
        spawn_move[p,home,home,,a,] <- 1
        spawn_move[p,away,home,,a,] <- returning
        spawn_move[p,away,away,,a,] <- 1 - returning

        # movement fractions, rows summing to one
        stay <- c(0.85 - 0.05 * a, 0.6 + 0.04 * a) * c(1, 0.9)[p]
        for(y in 1:n_yrs) {
          Movement[p,1,,y,1,a,1] <- c(stay[1], 1 - stay[1])
          Movement[p,2,,y,1,a,1] <- c(1 - stay[2], stay[2])
        } # end y loop

      } # end a loop
    } # end p loop

    il <- suppressWarnings(suppressMessages(sweep_input(
      dims = list(n_pop = n_pop, natal_region = c(1, 2), n_regions = n_regions, n_sexes = 1, n_fish_fleets = 1, n_yrs = n_yrs, n_ages = n_ages),
      rec = list(
        rec_model = "bh_rec",
        rec_dd = "local", # one curve per population
        h_spec = "fix",
        init_age_strc = 2,
        t_spawn = 0.4, # spawning 40 pct of the way into the year
        sgl_seas_spawning_movement = spawn_move
      ),
      biol = list(Fixed_natmort = natmort),
      move = list(Fixed_Movement = Movement, do_recruits_move = 0)
    )))

    # a tenth of population 1's spawning biomass and a quarter of population 2's count toward the other population
    il$data$use_fixed_stray_rate <- 1
    il$data$fixed_stray_rate[1,] <- 0.10
    il$data$fixed_stray_rate[2,] <- 0.25
    il$par$ln_F_mean[] <- log(c(0.15, 0.08)) # fishing mortality by region

    fit <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))
    cache <<- list(data = il$data, rep = fit$rep)
    cache
  }
})

# spawning biomass per recruit [population, region], one year class of each population fished at F in both regions
refpt_pop_per_recruit <- function(Fmort, data, rep, rec_region_prop) {

  n_pop <- data$n_pop
  n_regions <- data$n_regions
  n_ages <- length(data$ages)
  n_yrs <- length(data$years)
  sbpr <- matrix(0, n_pop, n_regions)

  for(p in 1:n_pop) {

    numbers <- 0.5 * rec_region_prop[p,] # female recruits by region, for one recruit of this population

    for(year in 1:800) {

      a <- min(year, n_ages) # the plus group keeps its own values every year after

      # total mortality by region, with dead discards
      retained <- rep$fish_sel[p,,n_yrs,1,a,1,1] * rep$ret_sel[p,,n_yrs,1,a,1,1]
      dead_discards <- rep$fish_sel[p,,n_yrs,1,a,1,1] * (1 - rep$ret_sel[p,,n_yrs,1,a,1,1]) * rep$dmr[,n_yrs,1,1]
      Z <- rep$natmort[p,,n_yrs,1,a,1] + Fmort * (retained + dead_discards)

      # recruits do not move in this model
      moved <- numbers
      if(a > 1) moved <- as.vector(numbers %*% rep$Movement[p,,,n_yrs,1,a,1])

      # fish home to spawn, then die where they spawn up to spawning
      homed <- as.vector(moved %*% rep$sgl_seas_spawning_movement[p,,,n_yrs,a,1])
      at_spawning <- homed * exp(-data$t_spawn * Z)
      sbpr[p,] <- sbpr[p,] + at_spawning * data$WAA[p,,n_yrs,1,a,1] * data$MatAA[p,,n_yrs,1,a,1]

      numbers <- moved * exp(-Z) # homing is for spawning only

    } # end year loop
  } # end p loop

  sbpr
}

# a population's own spawning biomass in its natal region plus the other populations' strays there
refpt_pop_effective <- function(biomass, data, stray_rate) {
  effective <- rep(0, data$n_pop)
  for(p in 1:data$n_pop) {
    home <- data$natal_region[p]
    n_born_here <- sum(data$natal_region == home)
    effective[p] <- biomass[p,home]
    for(p_other in 1:data$n_pop) {
      if(p_other != p) effective[p] <- effective[p] + stray_rate[p_other] / n_born_here * biomass[p_other,home]
    } # end p_other loop
  } # end p loop
  effective
}

# Tests ----------------------------------------------------------------------

test_that("the global F at a target SPR leaves that share of the populations' unfished effective spawning biomass", {

  model <- refpt_pop_model()
  data <- model$data
  rep <- model$rep
  n_yrs <- length(data$years)

  ref <- suppressWarnings(Get_Reference_Points(data = data, rep = rep, SPR_x = 0.4, type = "multi_region", what = "global_SPR"))

  total_rec <- rowSums(ref$mean_rec) # mean recruitment by population
  rec_region_prop <- ref$mean_rec / total_rec # share of each population's recruits by region
  stray_rate <- rep$stray_rate[,n_yrs]

  fished <- refpt_pop_per_recruit(as.numeric(ref$f_ref_pt[1]), data, rep, rec_region_prop)
  unfished <- refpt_pop_per_recruit(0, data, rep, rec_region_prop)

  # one F for both regions, at the target
  expect_equal(as.numeric(ref$f_ref_pt[1]), as.numeric(ref$f_ref_pt[2]), tolerance = 1e-12)
  expect_equal(sum(refpt_pop_effective(fished, data, stray_rate)) / sum(refpt_pop_effective(unfished, data, stray_rate)), 0.4, tolerance = 1e-6)

  # biomass reference points by population and region are per recruit values times the population's mean recruitment
  expect_equal(ref$b_ref_pt, fished * total_rec, tolerance = 1e-6)
  expect_equal(ref$virgin_b_ref_pt, unfished * total_rec, tolerance = 1e-6)

  # population reference points: effective spawning biomass in the natal region, zero elsewhere
  for(p in 1:data$n_pop) {
    home <- data$natal_region[p]
    expect_equal(ref$pop_b_ref_pt[p,home], refpt_pop_effective(fished * total_rec, data, stray_rate)[p], tolerance = 1e-6)
    expect_equal(ref$virgin_pop_b_ref_pt[p,home], refpt_pop_effective(unfished * total_rec, data, stray_rate)[p], tolerance = 1e-6)
    expect_equal(ref$pop_b_ref_pt[p,3 - home], 0)
  } # end p loop

})

test_that("homing, straying and the two populations all enter the answer in this setup", {

  # if any of these did nothing, the test above could not tell a calculation that left it out
  model <- refpt_pop_model()
  data <- model$data
  rep <- model$rep
  n_yrs <- length(data$years)
  ref <- suppressWarnings(Get_Reference_Points(data = data, rep = rep, SPR_x = 0.4, type = "multi_region", what = "global_SPR"))
  unfished <- refpt_pop_per_recruit(0, data, rep, ref$mean_rec / rowSums(ref$mean_rec))

  # the populations differ, and strays add to each population's effective spawning biomass
  expect_gt(abs(unfished[1,1] / unfished[2,2] - 1), 0.05)
  effective <- refpt_pop_effective(unfished, data, rep$stray_rate[,n_yrs])
  expect_gt(effective[1] / unfished[1,1], 1.01)
  expect_gt(effective[2] / unfished[2,2], 1.01)

  # fish are found outside their natal region, so homing moves spawners
  expect_gt(unfished[1,2] / unfished[1,1], 0.05)

})

# Long Projections -----------------------------------------------------------

# the model projected for 300 years at an annual F by region
refpt_pop_project <- function(ref, f_by_region, recruitment_opt) {
  model <- refpt_pop_model()
  ref$f_ref_pt <- f_by_region
  inputs <- suppressWarnings(suppressMessages(build_projection_inputs(data = model$data, rep = model$rep, refpts = ref, n_proj_yrs = 300,
                                                                      t_spawn = model$data$t_spawn, recruitment_opt = recruitment_opt)))
  suppressWarnings(suppressMessages(do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "Input")))))
}

test_that("a long projection at the global SPR F settles on the biomass reference points", {

  model <- refpt_pop_model()
  ref <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, SPR_x = 0.4, type = "multi_region", what = "global_SPR"))
  out <- refpt_pop_project(ref, ref$f_ref_pt, recruitment_opt = "mean_rec")

  # spawning biomass by population and region, and each population's effective spawning biomass in its natal region
  expect_equal(out$proj_SSB[,,300], ref$b_ref_pt, tolerance = 1e-5)
  for(p in 1:model$data$n_pop) {
    expect_equal(out$proj_eff_SSB[p,300], ref$pop_b_ref_pt[p,model$data$natal_region[p]], tolerance = 1e-5)
  } # end p loop

})

test_that("a long projection at the local Fmsy settles on Bmsy, and takes more catch than fishing harder or lighter", {

  model <- refpt_pop_model()
  ref <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, type = "multi_region", what = "local_MSY"))

  at_fmsy <- refpt_pop_project(ref, ref$f_ref_pt, recruitment_opt = "bh_rec")
  expect_equal(at_fmsy$proj_SSB[,,300], ref$b_ref_pt, tolerance = 1e-5)
  expect_equal(at_fmsy$proj_SSB[,,300], at_fmsy$proj_SSB[,,299], tolerance = 1e-10) # settled

  # long-run catch summed over populations, regions and fleets
  catch_at_fmsy <- sum(at_fmsy$proj_Catch[,,300,,])
  lighter <- refpt_pop_project(ref, 0.9 * ref$f_ref_pt, recruitment_opt = "bh_rec")
  harder <- refpt_pop_project(ref, 1.1 * ref$f_ref_pt, recruitment_opt = "bh_rec")
  expect_gt(catch_at_fmsy, sum(lighter$proj_Catch[,,300,,]))
  expect_gt(catch_at_fmsy, sum(harder$proj_Catch[,,300,,]))

  # fishing region 1 a little, where Fmsy is zero in this model, also takes less
  shifted <- refpt_pop_project(ref, ref$f_ref_pt + c(0.02, 0), recruitment_opt = "bh_rec")
  expect_gt(catch_at_fmsy, sum(shifted$proj_Catch[,,300,,]))

})

# Independent Reference Points -----------------------------------------------

# spawning biomass and retained yield per recruit by population in one region, with no movement in or out
refpt_pop_closed_region <- function(Fmort, r, data, rep) {

  n_pop <- data$n_pop
  n_ages <- length(data$ages)
  n_yrs <- length(data$years)
  sbpr <- ypr <- rep(0, n_pop)

  for(p in 1:n_pop) {

    numbers <- 0.5 # one recruit, half of them female

    for(year in 1:800) {

      a <- min(year, n_ages) # the plus group keeps its own values every year after
      ret_F <- Fmort * rep$fish_sel[p,r,n_yrs,1,a,1,1] * rep$ret_sel[p,r,n_yrs,1,a,1,1]
      disc_F <- Fmort * rep$fish_sel[p,r,n_yrs,1,a,1,1] * (1 - rep$ret_sel[p,r,n_yrs,1,a,1,1]) * rep$dmr[r,n_yrs,1,1]
      Z <- rep$natmort[p,r,n_yrs,1,a,1] + ret_F + disc_F

      sbpr[p] <- sbpr[p] + numbers * exp(-data$t_spawn * Z) * data$WAA[p,r,n_yrs,1,a,1] * data$MatAA[p,r,n_yrs,1,a,1]
      ypr[p] <- ypr[p] + ret_F / Z * numbers * (1 - exp(-Z)) * data$WAA_fish[p,r,n_yrs,1,a,1,1]
      numbers <- numbers * exp(-Z)

    } # end year loop
  } # end p loop

  list(sbpr = sbpr, ypr = ypr)
}

# each population's own value plus the other populations' strays, shared among those born in its natal region
refpt_pop_with_strays <- function(by_pop, data, stray_rate) {
  out <- rep(0, data$n_pop)
  for(p in 1:data$n_pop) {
    n_born_here <- sum(data$natal_region == data$natal_region[p])
    out[p] <- by_pop[p]
    for(p_other in 1:data$n_pop) {
      if(p_other != p) out[p] <- out[p] + stray_rate[p_other] / n_born_here * by_pop[p_other]
    } # end p_other loop
  } # end p loop
  out
}

test_that("independent SPR treats each region as closed and reaches the target there", {

  model <- refpt_pop_model()
  data <- model$data
  rep <- model$rep
  stray_rate <- rep$stray_rate[,length(data$years)]
  ref <- suppressWarnings(Get_Reference_Points(data = data, rep = rep, SPR_x = 0.4, type = "multi_region", what = "independent_SPR"))

  for(r in 1:data$n_regions) {

    fished <- refpt_pop_closed_region(as.numeric(ref$f_ref_pt[r]), r, data, rep)
    unfished <- refpt_pop_closed_region(0, r, data, rep)
    label <- sprintf("region %d", r)

    expect_equal(sum(refpt_pop_with_strays(fished$sbpr, data, stray_rate)) / sum(refpt_pop_with_strays(unfished$sbpr, data, stray_rate)), 0.4,
                 tolerance = 1e-6, label = paste("realized SPR,", label))

    # per recruit values times each population's mean recruitment in the region
    expect_equal(ref$b_ref_pt[,r], fished$sbpr * ref$mean_rec[,r], tolerance = 1e-6, label = paste("b_ref_pt,", label))
    expect_equal(ref$virgin_b_ref_pt[,r], unfished$sbpr * ref$mean_rec[,r], tolerance = 1e-6, label = paste("virgin_b_ref_pt,", label))
    expect_equal(ref$pop_b_ref_pt[,r], refpt_pop_with_strays(fished$sbpr * ref$mean_rec[,r], data, stray_rate), tolerance = 1e-6,
                 label = paste("pop_b_ref_pt,", label))
    expect_equal(ref$virgin_pop_b_ref_pt[,r], refpt_pop_with_strays(unfished$sbpr * ref$mean_rec[,r], data, stray_rate), tolerance = 1e-6,
                 label = paste("virgin_pop_b_ref_pt,", label))

  } # end r loop

})

test_that("independent MSY maximizes each closed region's equilibrium yield over the populations' curves", {

  model <- refpt_pop_model()
  data <- model$data
  rep <- model$rep
  stray_rate <- rep$stray_rate[,length(data$years)]
  ref <- suppressWarnings(Get_Reference_Points(data = data, rep = rep, type = "multi_region", what = "independent_MSY"))

  for(r in 1:data$n_regions) {

    R0_region <- as.numeric(rep$R0) * rep$rec_region_prop[,r] # each population's unfished recruitment in this region
    h <- rep$h_trans[,r] # steepness
    unfished <- refpt_pop_closed_region(0, r, data, rep)
    unfished_effective <- refpt_pop_with_strays(unfished$sbpr, data, stray_rate)

    # equilibrium recruitment by population under Beverton-Holt, and the yield it gives
    equilibrium <- function(Fmort) {
      fished <- refpt_pop_closed_region(Fmort, r, data, rep)
      fished_effective <- refpt_pop_with_strays(fished$sbpr, data, stray_rate)
      rec <- R0_region * (4 * h * fished_effective - (1 - h) * unfished_effective) / ((5 * h - 1) * fished_effective)
      list(rec = rec, sbpr = fished$sbpr, yield = sum(rec * fished$ypr))
    }

    peak <- stats::optimize(function(Fmort) equilibrium(Fmort)$yield, interval = c(0.001, 2), maximum = TRUE, tol = 1e-10)
    at_peak <- equilibrium(peak$maximum)
    label <- sprintf("region %d", r)

    expect_equal(as.numeric(ref$f_ref_pt[r]), peak$maximum, tolerance = 1e-5, label = paste("Fmsy,", label))
    expect_equal(ref$b_ref_pt[,r], at_peak$sbpr * at_peak$rec, tolerance = 1e-5, label = paste("b_ref_pt,", label))
    expect_equal(ref$virgin_b_ref_pt[,r], unfished$sbpr * R0_region, tolerance = 1e-6, label = paste("virgin_b_ref_pt,", label))
    expect_equal(ref$pop_b_ref_pt[,r], refpt_pop_with_strays(at_peak$sbpr * at_peak$rec, data, stray_rate), tolerance = 1e-5,
                 label = paste("pop_b_ref_pt,", label))
    expect_equal(ref$virgin_pop_b_ref_pt[,r], refpt_pop_with_strays(unfished$sbpr * R0_region, data, stray_rate), tolerance = 1e-6,
                 label = paste("virgin_pop_b_ref_pt,", label))

    # each population's recruits are split between the regions here, so its total R0 would give a different answer
    expect_lt(max(rep$rec_region_prop[,r]), 0.9)

  } # end r loop

})
