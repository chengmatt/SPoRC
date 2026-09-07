# Purpose: compare the seeded SPoRC run against the wham fit_0 report, quantity by quantity
# Date Created: 9/7/26

library(here)

out_dir <- here("dev", "wham_bsb_bridge", "output")
dat <- readRDS(file.path(out_dir, "02_sporc_data.rds"))
sporc <- readRDS(file.path(out_dir, "04_seeded.rds"))$rep
wham <- dat$wham

n_pop <- dat$dims$n_pop
n_regions <- dat$dims$n_regions
n_yrs <- dat$dims$n_yrs
n_ages <- dat$dims$n_ages

# relative error over the cells wham holds away from zero
rel_err <- function(sporc_vals, wham_vals) {
  keep <- which(abs(wham_vals) > 1e-8)
  if(!length(keep)) return(NA_real_)
  max(abs(sporc_vals[keep] - wham_vals[keep]) / abs(wham_vals[keep]))
}

comparison <- data.frame(quantity = character(), n_cells = integer(), max_rel_err = numeric())

add_row <- function(comparison, label, sporc_vals, wham_vals) {
  rbind(comparison, data.frame(quantity = label,
                               n_cells = sum(abs(wham_vals) > 1e-8),
                               max_rel_err = rel_err(sporc_vals, wham_vals)))
}

# numbers at age on january 1, the state wham imposes and SPoRC is seeded at
comparison <- add_row(comparison, "NAA jan 1",
                      sporc$NAA[, , 1:n_yrs, 1, , 1], wham$NAA)

# one step ahead prediction from the seasonal dynamics, ages 2 and above
comparison <- add_row(comparison, "NAA predicted, ages 2+",
                      sporc$NAA_pred[, , 2:n_yrs, 1, 2:n_ages, 1], wham$pred_NAA[, , 2:n_yrs, 2:n_ages])

# spawning biomass, summed over the regions each population occupies
comparison <- add_row(comparison, "SSB",
                      t(apply(sporc$SSB, c(1, 3), sum)), wham$SSB)

# predicted catch, summed over populations and seasons to match wham's annual total
pred_catch <- apply(sporc$PredCatch, c(2, 3, 5), sum)
for(f in 1:length(dat$dims$fleet_region)) {
  comparison <- add_row(comparison, paste0("predicted catch, fleet ", f),
                        pred_catch[dat$dims$fleet_region[f], , dat$dims$fleet_gear[f]],
                        wham$pred_catch[, f])
}

# predicted catch at age, summed over populations and seasons
pred_caa <- apply(sporc$CAA, c(2, 3, 5, 7), sum)
for(f in 1:length(dat$dims$fleet_region)) {
  comparison <- add_row(comparison, paste0("predicted catch at age, fleet ", f),
                        pred_caa[dat$dims$fleet_region[f], , , dat$dims$fleet_gear[f]],
                        wham$pred_CAA[f, , ])
}

# predicted survey indices, read in the season each index sits in
pred_srv <- apply(sporc$PredSrvIdx, c(2, 3, 4, 5), sum)
for(i in 1:length(dat$dims$index_region)) {
  comparison <- add_row(comparison, paste0("predicted index ", i, ", season ", dat$dims$index_seas[i]),
                        pred_srv[dat$dims$index_region[i], , dat$dims$index_seas[i], dat$dims$index_gear[i]],
                        wham$pred_indices[, i])
}

# predicted survey index at age, in the index's own season
for(i in 1:length(dat$dims$index_region)) {
  seas <- dat$dims$index_seas[i]
  iaa <- apply(sporc$SrvIAA[, dat$dims$index_region[i], , seas, , 1, dat$dims$index_gear[i]], c(2, 3), sum)
  comparison <- add_row(comparison, paste0("predicted index at age ", i),
                        iaa * wham$q[1, i], wham$pred_IAA[i, , ])
}

# wham survives fish for the index timing and then applies the whole season's movement
# before reading the index; SPoRC moves at the end of the season, so rebuild wham's order
pred_srv_moved <- matrix(0, n_yrs, length(dat$dims$index_region))
for(i in 1:length(dat$dims$index_region)) {
  r <- dat$dims$index_region[i]
  seas <- dat$dims$index_seas[i]
  blk <- dat$dims$index_block[i]
  t_i <- dat$t_srv[r, seas, dat$dims$index_gear[i]]

  for(y in 1:n_yrs) {
    total <- 0
    for(p in 1:n_pop) for(a in 1:n_ages) {
      survived <- sporc$NAA[p, , y, seas, a, 1] * exp(-t_i * sporc$ZAA[p, , y, seas, a, 1])
      after_move <- as.vector(t(survived) %*% dat$Fixed_Movement[p, , , y, seas, a, 1])
      total <- total + after_move[r] * dat$selAA[blk, y, a]
    } # end a loop
    pred_srv_moved[y, i] <- wham$q[1, i] * total
  } # end y loop

  comparison <- add_row(comparison, paste0("predicted index ", i, ", movement applied first"),
                        pred_srv_moved[, i], wham$pred_indices[, i])
}

write.csv(comparison, file.path(out_dir, "05_bridge_comparison.csv"), row.names = FALSE)

# year by year spawning biomass, both models side by side
ssb_sporc <- t(apply(sporc$SSB, c(1, 3), sum))
ssb_table <- data.frame(year = dat$dims$years,
                        SSB_north_sporc = ssb_sporc[, 1],
                        SSB_north_wham = wham$SSB[, 1],
                        SSB_south_sporc = ssb_sporc[, 2],
                        SSB_south_wham = wham$SSB[, 2])

write.csv(ssb_table, file.path(out_dir, "05_ssb_by_year.csv"), row.names = FALSE)
