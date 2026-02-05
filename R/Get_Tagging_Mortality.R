#' Get values for total tag mortality and tag fishing mortality
#'
#' @param tag_selex Tag selectivity options, == 0 (uniform, with F from fleet 1 (dominant fleet)), == 1 (sex-averaged selectivity, with F from fleet 1 (dominant fleet)), == 2 (sex-specific selectivity, with F from fleet 1 (dominant fleet)), 3 (uniform with F summed across fleets), 4 (sex averaged with F summed across fleets, weighted sum), 5 (sex-specific with F summed across fleets, weighted sum)
#' @param tag_natmort Tag natural mortality options == 0 (averaged across sexes and ages), == 1 (averaged across sexes, but unique for ages), == 2 (sex-specific, but averaged across ages), == 3 (sex-and age-specific)
#' @param Fmort Array of fishing mortality, dimensioned by n_region, n_years, n_fish_fleets
#' @param natmort Array of fishing mortality, dimensioned by n_region, n_years, n_ages, n_sexes
#' @param Tag_Shed Scalar chronic tag shedding rate (annual)
#' @param fish_sel Array of fishery selectivity, dimensioned by n_region, n_years, n_ages, n_sexes, n_fish_fleetss
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param y Year index
#' @param what Whether to return Z or F (total or fishing mortality)
#' @param n_fish_fleets Number of fishery fleets
#' @param seas Season index
#' @param seasdur Fraction of the year the season occurs in
#'
#' @returns Z or F values from tagging specifications
#' @keywords internal
#'
Get_Tagging_Mortality <- function(tag_selex,
                                  tag_natmort,
                                  Fmort,
                                  natmort,
                                  Tag_Shed,
                                  fish_sel,
                                  n_regions,
                                  n_ages,
                                  n_sexes,
                                  n_fish_fleets,
                                  y,
                                  seas,
                                  seasdur,
                                  what
                                  ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Parameterizations for Tag Selectivity
  # Uniform selectivity, with F from dominant fleet
  if(tag_selex == 0) tmp_F = array(Fmort[,y,seas,1], dim = c(n_regions, 1, n_ages, n_sexes))

  # Sex averaged selectivity from dominant fleet
  if(tag_selex == 1) {
    tmp_fish_sel = array(apply(fish_sel[,y,,,1], c(1,2), mean), dim = c(n_regions, n_ages, n_sexes))
    tmp_F = array(Fmort[,y,seas,1] * tmp_fish_sel, dim = c(n_regions, 1, n_ages, n_sexes))
  }

  # Sex-specific selectivity from dominant fleet
  if(tag_selex == 2) tmp_F = array(Fmort[,y,seas,1] * fish_sel[,y,,,1], dim = c(n_regions, 1, n_ages, n_sexes))

  # Uniform selectivity, with F from all fleets
  if(tag_selex == 3) tmp_F = array(rowSums(Fmort[,y,seas,]), dim = c(n_regions, 1, n_ages, n_sexes))

  # Sex averaged selectivity with weighted sum from all fleets
  if(tag_selex == 4) {
    tmp_F = array(0, dim = c(n_regions, 1, n_ages, n_sexes))
    tmp_fish_sel = apply(fish_sel[,y,,,,drop = FALSE], c(1,3,5), mean)
    # loop through to populate elements
    for(r in 1:n_regions) for(a in 1:n_ages) for(f in 1:n_fish_fleets) for(s in 1:n_sexes)
      tmp_F[r,1,a,s] = tmp_F[r,1,a,s] + Fmort[r,y,seas,f] * tmp_fish_sel[r,a,f]
  }

  # Sex-specific selectivity with weighted sum from all fleets
  if(tag_selex == 5) {
    tmp_F = array(0, dim = c(n_regions, 1, n_ages, n_sexes))
    # loop through to populate elements
    for(r in 1:n_regions) for(a in 1:n_ages) for(s in 1:n_sexes)
      tmp_F[r,1,a,s] = sum(Fmort[r,y,seas,] * fish_sel[r,y,a,s,])
  }

  # Parameterizations for natural mortality
  # Mean natural mortality across ages and sexes
  if(tag_natmort == 0) tmp_natmort = array(apply(natmort[,y,,,drop=FALSE], 1, mean), dim = c(n_regions, 1, n_ages, n_sexes))

  # Age-specific, sex-aggregated natural mortality
  if(tag_natmort == 1) {
    tmp_mean_natmort = apply(natmort[,y,,,drop=FALSE], c(1,3), mean)
    tmp_natmort = array(0, dim = c(n_regions, 1, n_ages, n_sexes))
    for(i in 1:n_regions) for(j in 1:n_ages) tmp_natmort[i,,j,] = tmp_mean_natmort[i,j] # fill in array
  }

  # Sex-specific, age-aggregated natural mortality
  if(tag_natmort == 2) {
    tmp_mean_natmort = apply(natmort[,y,,,drop = FALSE], c(1, 4), mean)
    tmp_natmort = array(0, dim = c(n_regions, 1, n_ages, n_sexes))
    for(i in 1:n_regions) for(j in 1:n_sexes) tmp_natmort[i,,,j] = tmp_mean_natmort[i,j] # fill in array
  }

  # Age and sex-specific natural mortality
  if(tag_natmort == 3) tmp_natmort = natmort[,y,,,drop = FALSE]

  if(what == "Z") val = tmp_F + ((Tag_Shed + tmp_natmort) * seasdur[seas])
  if(what == "F") val = tmp_F

  return(val)
}
