# Shared routines for the conditional age-at-length row tests: rex sole at the assessment's estimate, its survey
# age-at-length summed onto coarser rows, and a map pairing adjacent model length bins.

rex_caal_input <- function() seed_goa_rex_mle(suppressWarnings(suppressMessages(build_goa_rex_input(mlt_rg_goa_rex_data))), mlt_rg_goa_rex_data)

# the rex survey age-at-length data summed onto the rows of a map, as data recorded on those rows would be
caal_on_rows <- function(data, caal_len_bin_map) {
  n_rows <- ncol(caal_len_bin_map)
  obs <- array(0, dim = replace(dim(data$ObsSrv_caal), 4, n_rows))
  iss <- array(0, dim = replace(dim(data$ISS_Srv_caal), 4, n_rows))
  use <- array(0, dim = replace(dim(data$UseSrv_caal), 4, n_rows))
  wt <- array(0, dim = replace(dim(data$Wt_Srv_caal), 4, n_rows))
  for(k in 1:n_rows) for(l in which(caal_len_bin_map[,k] == 1)) {
    obs[,,,k,,,] <- obs[,,,k,,,] + data$ObsSrv_caal[,,,l,,,]
    iss[,,,k,,] <- iss[,,,k,,] + data$ISS_Srv_caal[,,,l,,]
    use[,,,k,] <- pmax(use[,,,k,], data$UseSrv_caal[,,,l,])
    wt[,,,k,,] <- data$Wt_Srv_caal[,,,l,,] # one weight per fleet, the same in every bin
  } # end k, l loops
  data$CAAL_LenBinMap <- caal_len_bin_map
  data$ObsSrv_caal <- obs
  data$ISS_Srv_caal <- iss
  data$UseSrv_caal <- use
  data$Wt_Srv_caal <- wt
  data
}

# a 0/1 map of model length bins onto rows that each take two adjacent bins, the last alone when the count is odd
bin_pairs <- function(n_lens) {
  pairs <- matrix(0, n_lens, ceiling(n_lens / 2))
  for(l in 1:n_lens) pairs[l, ceiling(l / 2)] <- 1
  pairs
}
