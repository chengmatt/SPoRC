# Shared routines for checking the internal one-step-ahead compositions against the ordinary likelihood: data
# made whole fish, so what the OSA path packs is the data itself, and the gaps between the two paths.

osa_comp_sources <- c("FishAgeComps", "FishLenComps", "SrvAgeComps", "SrvLenComps",
                      "FishAgeComps_discard", "FishLenComps_discard",
                      "FishAgeComps_pop", "FishLenComps_pop", "SrvAgeComps_pop", "SrvLenComps_pop",
                      "FishAgeComps_discard_pop", "FishLenComps_discard_pop",
                      "Fish_caal", "Srv_caal")

# one composition data source as whole fish with the sample size their total and a weight of one
as_whole_fish <- function(data, source, n_fish = 200) {

  use <- data[[paste0("Use", source)]]
  if(is.null(use) || !any(use == 1)) return(data)
  obs <- data[[paste0("Obs", source)]]
  iss <- data[[paste0("ISS_", source)]]
  type_mat <- data[[paste0(source, "_Type")]]
  is_caal <- grepl("_caal$", source) # a length row ahead of the ages
  is_pop <- grepl("_pop$", source) # a population ahead of the region
  whole <- function(block) round(block / sum(block) * n_fish)

  # every array read as [pop, region, year, season, row, bin, sex, fleet], with a pop or row of one when absent
  d <- dim(obs)
  n_pop <- if(is_pop) d[1] else 1
  n_rows <- if(is_caal) d[4] else 1
  obs_full <- array(obs, dim = c(n_pop, d[if(is_pop) 2 else 1], d[if(is_pop) 3 else 2], d[if(is_pop) 4 else 3], n_rows,
                                 d[length(d) - 2], d[length(d) - 1], d[length(d)]))
  iss_full <- array(iss, dim = c(dim(obs_full)[1:5], dim(obs_full)[7:8]))
  use_full <- array(use, dim = c(dim(obs_full)[1:5], dim(obs_full)[8]))
  fd <- dim(obs_full)

  for(p in 1:fd[1]) for(r in 1:fd[2]) for(y in 1:fd[3]) for(seas in 1:fd[4]) for(k in 1:fd[5]) for(f in 1:fd[8]) {
    if(use_full[p,r,y,seas,k,f] != 1) next
    ct <- type_mat[y,f]
    if(ct == 2) {
      # joint across sexes: one sample over the bin by sex stack
      block <- obs_full[p,r,y,seas,k,,,f]
      if(sum(block) == 0) next
      obs_full[p,r,y,seas,k,,,f] <- whole(block)
      iss_full[p,r,y,seas,k,1,f] <- sum(obs_full[p,r,y,seas,k,,,f])
    } else {
      # split by sex, or aggregated into the first sex: a sample per sex
      for(s in (if(ct == 1) 1:fd[7] else 1)) {
        block <- obs_full[p,r,y,seas,k,,s,f]
        if(sum(block) == 0) next
        obs_full[p,r,y,seas,k,,s,f] <- whole(block)
        iss_full[p,r,y,seas,k,s,f] <- sum(obs_full[p,r,y,seas,k,,s,f])
      } # end s loop
    }
  } # end p, r, y, seas, k, f loops

  data[[paste0("Obs", source)]] <- array(obs_full, dim = dim(obs))
  data[[paste0("ISS_", source)]] <- array(iss_full, dim = dim(iss))
  data[[paste0("Wt_", source)]][] <- 1
  data
}

# the ordinary and internal OSA paths at the start and at a nudged parameter vector: the objective gap, a
# constant when the two are the same likelihood, and the worst gradient gap relative to the largest gradient
osa_path_gaps <- function(data, par, map, random = NULL) {

  data$comp_const_obs <- 0
  data$addtocomp <- 1e-12
  for(source in osa_comp_sources) data <- as_whole_fish(data, source)
  data$do_internal_comp_osa <- FALSE
  ordinary <- fit_model(data, par, map, random = random, do_optim = FALSE, silent = TRUE)
  data$do_internal_comp_osa <- TRUE
  internal <- fit_model(data, par, map, random = random, do_optim = FALSE, silent = TRUE)

  set.seed(1)
  nudged <- ordinary$par + stats::rnorm(length(ordinary$par), 0, 0.02)
  gaps <- list()
  for(at in list(ordinary$par, nudged)) {
    gradient <- ordinary$gr(at)
    gaps[[length(gaps) + 1]] <- c(objective = internal$fn(at) - ordinary$fn(at),
                                  gradient = max(abs(internal$gr(at) - gradient)) / max(abs(gradient)))
  } # end at loop
  gaps
}
