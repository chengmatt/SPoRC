# label_dims names the dims of the report, parameter and data lists without changing a value. Dims are
# checked on a model where no two dims are the same size, so a dim given the wrong name has the wrong
# length, and on rex for the length bins, conditional age-at-length and growth.

# 4 populations, 3 regions, 13 years, 6 seasons, 7 ages, 2 sexes, 5 fishery and 11 survey fleets
label_case <- local({
  cached <- list()
  function(which = "sweep") {
    if(!is.null(cached[[which]])) return(cached[[which]])
    il <- if(which == "sweep") {
      sweep_input(dims = list(n_pop = 4, n_seas = 6, n_srv_fleets = 11, natal_region = c(1, 2, 3, 1)),
                  rec = list(rec_dd = "local"))
    } else suppressWarnings(suppressMessages(build_goa_rex_input(mlt_rg_goa_rex_data)))
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    cached[[which]] <<- list(data = il$data, par = il$par, rep = obj$report(obj$env$last.par.best))
    cached[[which]]
  }
})

# vectors that list indices rather than holding one value per year, age or region
index_vectors <- c("years", "ages", "lens", "bias_year", "move_pairs", "move_re_pops", "move_re_sexes",
                   "move_re_years", "move_re_seas", "move_re_ages", "fishsel_devs_min_shared_bins",
                   "retsel_devs_min_shared_bins", "srvsel_devs_min_shared_bins", "growth_semipar_bins",
                   "dsem_var_names")

# how long each named dim should be in this model
expected_extents <- function(data) {
  n_yrs <- length(data$years)
  n_ages <- length(data$ages)
  n_lens <- length(data$lens)
  list(Pop = data$n_pop, Region = data$n_regions, Region_From = data$n_regions, Region_To = data$n_regions,
       Year = n_yrs + c(0, 1, data$n_proj_yrs_devs), Rec_Year = seq_len(n_yrs + data$n_proj_yrs_devs),
       Seas = data$n_seas, Age = n_ages, Init_Age = n_ages - 1, Sex = data$n_sexes,
       Fish_Fleet = data$n_fish_fleets, Srv_Fleet = data$n_srv_fleets, Len = n_lens,
       Obs_Age = c(dim(data$AgeingError)[3], n_ages), # data sources not fit keep the model ages
       Obs_Len = c(if(is.null(data$LenBinMap)) n_lens else ncol(data$LenBinMap), n_lens),
       CAAL_Len = c(if(is.null(data$CAAL_LenBinMap)) n_lens else ncol(data$CAAL_LenBinMap), n_lens))
}


test_that("every array the objective reports has its dims listed", {

  objective_text <- paste(deparse(body(SPoRC_rtmb)), collapse = "\n")
  reported <- unique(gsub("RTMB::REPORT\\(|\\)", "", regmatches(objective_text, gregexpr("RTMB::REPORT\\([A-Za-z0-9_.]+\\)", objective_text))[[1]]))
  expect_gt(length(reported), 100) # the search found the report section

  # likelihood totals and priors reported as one number
  scalars <- c("jnLL", "M_nLL", "h_nLL", "R0_nLL", "rinit_nLL", "fish_q_nLL", "srv_q_nLL", "sel_nLL",
               "Movement_nLL", "TagRep_nLL", "rec_prop_nLL", "growth_tv_nLL", "growth_semipar_nLL",
               "NAA_state_nLL", "dsem_nLL", "dsem_obs_nLL", "extra_nLL")

  # a new report needs its dims added to dim_roles
  expect_equal(setdiff(reported, c(names(dim_roles), scalars)), character(0))
  expect_equal(anyDuplicated(names(dim_roles)), 0)

})


test_that("every parameter and data array has its dims listed", {

  for(which in c("sweep", "rex")) {
    case <- label_case(which)
    for(list_name in c("par", "data")) {
      lst <- case[[list_name]]
      is_array <- vapply(lst, function(v) is.atomic(v) && !is.null(v) && (!is.null(dim(v)) || length(v) > 1), logical(1))
      unlisted <- setdiff(names(lst)[is_array], c(names(dim_roles), index_vectors))
      expect_equal(unlisted, character(0), info = paste(which, list_name))
    } # end list_name loop
  } # end which loop

})


test_that("each named dim has the length of what it is named for", {

  for(which in c("sweep", "rex")) {

    case <- label_case(which)
    expected <- expected_extents(case$data)

    n_checked <- 0
    for(list_name in c("rep", "par", "data")) {
      lst <- case[[list_name]]
      for(x_name in intersect(names(lst), names(dim_roles))) {

        x <- lst[[x_name]]
        if(!is.atomic(x) || is.null(x)) next
        roles <- dim_roles[[x_name]]
        extents <- if(is.null(dim(x))) length(x) else dim(x)
        if(is.null(dim(x)) && length(x) == 1) next # a scalar default, such as a weight of one, stands for any shape
        expect_equal(length(extents), length(roles), info = paste(which, list_name, x_name))

        for(i in seq_along(roles)) {
          if(roles[i] %in% names(expected)) expect_true(extents[i] %in% expected[[roles[i]]], info = paste(which, list_name, x_name, roles[i]))
        } # end i loop
        n_checked <- n_checked + 1

      } # end x_name loop
    } # end list_name loop
    expect_gt(n_checked, 400)

  } # end which loop

})


test_that("labeling changes names only, and leaves the lists it was given alone", {

  case <- label_case()
  for(list_name in c("rep", "par", "data")) {
    lst <- case[[list_name]]
    before <- lst
    lab <- label_dims(lst, case$data)
    expect_identical(lst, before)
    expect_identical(names(lab), names(lst))
    for(x_name in names(lst)) expect_identical(unname(lab[[x_name]]), unname(lst[[x_name]]), info = paste(list_name, x_name))
  } # end list_name loop

  # an array of a shape the table does not describe comes back untouched
  rep <- case$rep
  rep$SSB <- matrix(1, 2, 2)
  expect_null(dimnames(label_dims(rep, case$data)$SSB))

})


test_that("dims carry the model's years, ages and the names they are melted under", {

  case <- label_case()
  data <- case$data
  rep <- label_dims(case$rep, data)
  par <- label_dims(case$par, data)
  dat <- label_dims(data)

  expect_equal(names(dimnames(rep$NAA)), c("Pop", "Region", "Year", "Seas", "Age", "Sex"))
  expect_equal(dimnames(rep$NAA)$Year, as.character(c(data$years, max(data$years) + 1)))
  expect_equal(dimnames(rep$NAA)$Age, as.character(data$ages))
  expect_equal(names(dimnames(rep$Movement))[2:3], c("Region_From", "Region_To"))
  expect_equal(names(dimnames(rep$CAA))[7], "Fleet")
  expect_equal(dimnames(rep$Init_Rec_nLL)$Age, as.character(data$ages[-1])) # initial deviations start at age two
  expect_equal(names(rep$R0), as.character(1:4)) # vectors are named by population

  expect_equal(names(dimnames(par$ln_F_devs)), c("Region", "Year", "Seas", "Fleet"))
  expect_equal(names(dimnames(par$fish_fixed_sel_pars)), c("Region", "Sel_Par", "Block", "Sex", "Fleet"))
  expect_equal(dimnames(par$ln_RecDevs)$Year, as.character(data$years))

  expect_equal(names(dimnames(dat$ObsCatch)), c("Region", "Year", "Seas", "Fleet"))
  expect_equal(names(dimnames(dat$ObsFishAgeComps)), c("Region", "Year", "Seas", "Age", "Sex", "Fleet"))
  expect_equal(names(dimnames(dat$AgeingError)), c("Year", "Age", "Obs_Age")) # both age dims keep their own name

  melted <- reshape2::melt(rep$SSB)
  expect_equal(names(melted), c("Pop", "Region", "Year", "value"))
  expect_equal(sort(unique(melted$Year)), data$years)

})


test_that("length bins and conditional age-at-length are labeled on rex", {

  case <- label_case("rex")
  data <- case$data
  rep <- label_dims(case$rep, data)
  dat <- label_dims(data)

  expect_equal(names(dimnames(rep$SizeAgeTrans_srv))[5:6], c("Len", "Age"))
  expect_equal(dimnames(rep$SizeAgeTrans_srv)$Len, as.character(data$lens))
  expect_equal(names(dimnames(dat$ObsSrv_caal)), c("Region", "Year", "Seas", "Len", "Age", "Sex", "Fleet"))
  expect_equal(dimnames(dat$ObsSrvLenComps)$Len, as.character(data$lens))

})


test_that("user labels go on every dim of that kind, and labels already there are kept", {

  case <- label_case()
  lab <- label_dims(case$rep, case$data,
                    labels = list(Region = c("WGOA", "CGOA", "EGOA"),
                                  Sex = c("Female", "Male"),
                                  Srv_Fleet = paste0("srv", 1:11)
                                  ))

  expect_equal(dimnames(lab$NAA)$Sex, c("Female", "Male"))
  expect_equal(dimnames(lab$Movement)$Region_To, c("WGOA", "CGOA", "EGOA"))
  expect_equal(dimnames(lab$srv_sel)$Fleet, paste0("srv", 1:11))
  expect_equal(dimnames(lab$fish_sel)$Fleet, as.character(1:5)) # fishery fleets keep their numbers
  expect_equal(lab$SSB[, "EGOA", 1], case$rep$SSB[, 3, 1], ignore_attr = TRUE)

  # names on a dim that counts parameters are kept, numbered dims from the objective are not
  par <- case$par
  dimnames(par$fish_fixed_sel_pars) <- list(NULL, c("a50", "slope"), NULL, NULL, NULL)
  expect_equal(dimnames(label_dims(par, case$data)$fish_fixed_sel_pars)$Sel_Par, c("a50", "slope"))
  expect_equal(dimnames(lab$Movement)$Year, as.character(case$data$years))

})


test_that("labels that do not fit the model are refused", {

  case <- label_case()
  expect_error(label_dims(case$rep, case$data, labels = list(Sex = "Female")), "labels\\$Sex has 1 entries")
  expect_error(label_dims(case$rep, case$data, labels = list(Year = 1:13)), "Years, ages and length bins are read from data")
  expect_error(label_dims(case$rep, case$data, labels = c(Sex = "Female")), "named list")
  expect_error(label_dims(case$rep), "data must be the model's data list")

})
